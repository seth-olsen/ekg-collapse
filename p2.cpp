/*

ekg scalar field collapse to black hole 

compile with:

   g++ -std=c++11 -g -Wall -O2 p1.cpp -o p1 -lbbhutil

input parameters in terminal as:
(do NOT include .sdf extension in outfile)

 ./p1 <outfile> <lastpt> <save_pt> <nsteps> <save_step>
      <lam> <r2m> <rmin> <rmax> <dspn> <tol> <maxit>
      <ic_Dsq> <ic_r0> <ic_Amp> <check_step> <zero_pi> 
      <somm_cond> <dspn_bound> <wr_ires> <wr_res>
      <wr_sol> <wr_itn> <hold_const>

where the value of any <parameter> can be set with the
following ordered pair of command line arguments:

 -parameter parameter_value

default values can be found at the start of main()
*/

#include <iostream>
#include <algorithm> // for max_element()
#include <fstream> // for mass/itn files
#include <ctime> // for quick time analysis
#include "fda-fns.h"
#include "fda-io.h"

using namespace std;

#define _USE_MATH_DEFINES
#ifndef M_PI
cerr << "NO PI\n";
const double M_PI = 4.0*atan(1.0);
#endif

// **********************************************************
// **********************************************************
//                       FUNCTIONS
// **********************************************************
// **********************************************************


//***********************NEED TO CHANGE THIS******************************
// compute and write mass
void mass_check(const vector<double>& f1, const vector<double>& f2,
		double t, double dr, double rmin, ofstream& out_stream,
		const vector<double>& f_vec, const vector<double>& beta_vec)
{
  double mass = 0.0, rval = rmin;
  int k = 0;
  for (auto val : f2) {
    mass += (0.5*f_vec[k]*(sq(f1[k]) + sq(val)) +
	     beta_vec[k]*f1[k]*val) * 4*M_PI*rval*rval*dr;
    ++k;
    rval += dr;
  }
  out_stream << t <<","<< mass << endl;
  return;
}

// get coarsened arrays from fields for writing
void get_wr_arr(const vector<double>& f1, const vector<double>& f2,
		    vector<double>& wr1, vector<double>& wr2,
		    int one_past_last, int savept)
{
  int k, s = 0;
  for (k = 0; k < one_past_last; ++k) {
    wr1[k] = f1[s];
    wr2[k] = f2[s];
    s += savept;
  }
  return;
}

void get_wr_arr_ires(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& oldxi, const vector<double>& oldpi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, vector<double>& wrxi, vector<double>& wrpi,
		     vector<double>& iresxi, vector<double>& irespi,
		     int lastwrite, int savept, double lam, double dr)
{
  double rval = rmin;
  wrxi[0] = xi[0];
  wrpi[0] = pi[0];
  iresxi[0] = xi[0] - iresxi_f(xi, pi, alpha, beta, psi, 0, lam)
    - oldxi[0] - iresxi_f(oldxi, oldpi, alpha, beta, psi, 0, lam);
  irespi[0] = pi[0] - irespi_f(xi, pi, alpha, beta, psi, 0, lam, dr, rval)
    - oldpi[0] - irespi_f(oldxi, oldpi, alpha, beta, psi, 0, lam, dr, rval);
  int k, s = savept;
  for (k = 1; k < lastwrite; ++k) {
    rval += savept*dr;
    wrxi[k] = xi[s];
    wrpi[k] = pi[s];
    iresxi[k] = xi[s] - iresxi_c(xi, pi, alpha, beta, psi, s, lam)
      - oldxi[s] - iresxi_c(oldxi, oldpi, alpha, beta, psi, s, lam);
    irespi[k] = pi[s] - irespi_c(xi, pi, alpha, beta, psi, s, lam, dr, rval)
      - oldpi[s] - irespi_c(oldxi, oldpi, alpha, beta, psi, s, lam, dr, rval);
    s += savept;
  }
  wrxi[lastwrite] = xi[s];
  wrpi[lastwrite] = pi[s];
  return;
  return;
}

// **********************************************************
// **********************************************************
//             INITIAL AND BOUNDARY CONDITIONS
// **********************************************************
// **********************************************************

inline double ic_alpha(double r, double r2m) { return 1.0; }

inline double ic_beta(double r, double r2m) { return 0; }

inline double ic_psi(double r, double r2m) { return 1.0; }

// for gaussian field or sin(coeff*r)*cos(coeff*t)/(coeff*r)
inline double ic_sol(double r, double amp, double dsq, double r0)
{ return amp * exp(-(r - r0)*(r - r0)/dsq); }

inline double ic_xi(double r, double amp, double dsq, double r0)
{ return -2 * (r - r0) * amp * exp(-(r - r0)*(r - r0)/dsq) / dsq; }

inline double ic_pi(double r, double amp, double dsq, double r0, double ctot, double cxi)
{ return ctot*(cxi*ic_xi(r, amp, dsq, r0) + ic_sol(r, amp, dsq, r0)/r); }

inline void update_sol(const vector<double>& xi, const vector<double>& pi,
		       const vector<double>& alpha, const vector<double>& beta,
		       const vector<double>& psi, vector<double>& sol,
		       double dt, int savept, int wr_shape) {
  int k, s;
  for (k = 0; k < wr_shape; ++k) {
    s = k * savept;
    sol[k] += dt * (alpha[s]*sqin(psi[s])*pi[s] + beta[s]*xi[s]);
  }
}

inline double fda_xi(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double lam)
{
  return ( 0.25*lam* dx_ekg_c(beta, xi, alpha, pi, psi, ind) );
}

inline double fda_pi(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double lam, double r, double dr)
{
  return ( 0.25*lam*p4in(psi[ind])* sqin(r)*dp_ekg_c(beta, pi, alpha, xi, psi, ind, dr, r)
	   - (lam/3.0)* pi[ind]*( d_c(beta, ind) + beta[ind]*(6*d_c(psi, ind)/psi[ind]
							     + 4*dr/r) ) );
  // or dpi_c(beta, pi, alpha, xi, psi, ind, dr, r)*sqin(r) replaced by
  // d3pi_c(beta, pi, alpha, xi, psi, ind, dr, r)
}

inline double fda0_xi(const vector<double>& xi, const vector<double>& pi,
		      const vector<double>& alpha, const vector<double>& beta,
		      const vector<double>& psi, double lam)
{
  return ( 0.25 * lam * dx_ekg_f(beta, xi, alpha, pi, psi, 0) );
}

inline double fda0_pi(const vector<double>& xi, const vector<double>& pi,
		      const vector<double>& alpha, const vector<double>& beta,
		      const vector<double>& psi, double lam, double r, double dr)
{
  return ( 0.25*lam*p4in(psi[0])* sqin(r)*dp_ekg_f(beta, pi, alpha, xi, psi, 0, dr, r)
	   - (lam/3.0)* pi[0]*( d_f(beta, 0) + beta[0]*(6*d_f(psi, 0)/psi[0]
							     + 4*dr/r) ) );
  // or dpi_f(beta, pi, alpha, xi, psi, ind, dr, r)*sqin(r) replaced by
  // d3pi_f(beta, pi, alpha, xi, psi, ind, dr, r)
}

// **********************************************************
// **********************************************************
//                         PROGRAM
// **********************************************************
// **********************************************************

int main(int argc, char **argv)
{
  // **********************************************************
  // ******************** PARAMETERS **************************
  // **********************************************************
  
  // user-set parameters
  string outfile = "p2";
  int lastpt = 1000; // grid size
  int save_pt = 1; // write only every (save_pt)th grid point
  int nsteps = 4000; // time steps
  int save_step = 8; // write only every (save_step)th time step
  double lam = 0.25; // dt/dr
  double r2m = 2.0;
  double rmin = 0.0;
  double rmax = 100.0;
  double dspn = 0.7; // dissipation coefficient
  double tol = 0.000000000001; // iterative method tolerance
  int maxit = 25; // max iterations for debugging
  double ic_Dsq = 4.0; // gaussian width
  double ic_r0 = 50.0; // gaussian center
  double ic_Amp = 1.0; // gaussian amplitude
  int check_step = 100; // for monitoring invariant mass
  // note: set bools in command line with integers 1=true or 0=false
  bool zero_pi = false; // zero initial time derivative?
  bool somm_cond = true; // sommerfeld condition at outer bound?
  bool dspn_bound = false; // dissipate boundary points?
  bool wr_ires = false; // write ires?
  bool wr_res = false; // write res?
  bool wr_sol = false; // write sol?
  bool wr_itn = false; // write itn counts?
  // variable to hold constant across resolutions
  string hold_const = "lambda"; // "lambda", "dt", or "dr"

  map<string, string *> p_str {{"-outfile",&outfile},
      {"-hold_const",&hold_const}};
  map<string, int *> p_int {{"-lastpt",&lastpt}, {"-save_pt", &save_pt},
      {"-nsteps", &nsteps}, {"-save_step",&save_step}, {"-maxit",&maxit},
      {"-check_step", &check_step}};
  map<string, double *> p_dbl {{"-lam",&lam}, {"-r2m",&r2m}, {"-rmin",&rmin},
      {"-rmax",&rmax}, {"-dspn",&dspn}, {"-tol",&tol}, {"-ic_Dsq",&ic_Dsq},
      {"-ic_r0",&ic_r0}, {"-ic_Amp",&ic_Amp}};
  map<string, bool *> p_bool {{"-zero_pi",&zero_pi},
      {"-somm_cond",&somm_cond}, {"-dspn_bound",&dspn_bound},
      {"-wr_ires",&wr_ires}, {"-wr_res",&wr_res},
      {"-wr_sol",&wr_sol}, {"-wr_itn",&wr_itn}};
  map<string, string> params;
  param_collect(argv, argc, params);
  param_set(params, p_str, p_int, p_dbl, p_bool);

  // check that grid size (lastpt = npts-1) is divisible by save_pt 
  if (lastpt % save_pt != 0) {
    cout << "ERROR: save_pt = " << save_pt << " entered for grid size " << lastpt << endl;
    save_pt -= lastpt % save_pt;
    cout << "--> corrected: save_pt = " << save_pt << endl;
  }
  
  // OBTAIN RESOLUTION FACTORS
  string resn_str;
  int num_resns;
  cout << "enter integer number of resolutions:" << endl;
  cin >> resn_str;
  num_resns = atoi(&resn_str[0]);
  vector<int> resolutions(num_resns);
  for (int k = 0; k < num_resns; ++k) {
    cout << "enter integer resolution factor at level " << k << ":" <<endl;
    cin >> resn_str;
    resolutions[k] = atoi(&resn_str[0]);
  }

  // bbhutil parameters for writing data to sdf
  int lastwr = lastpt/save_pt;
  int wr_shape = lastwr + 1;
  vector<double> wr_xi(wr_shape), wr_pi(wr_shape);
  double *field_arr[2] = {&wr_xi[0], &wr_pi[0]};
  int *bbh_shape = &wr_shape;
  int bbh_rank = 1;
  double coord_lims[2] = {rmin, rmax};
  double *coords = &coord_lims[0];
  
  int ires_size = ((wr_ires) ? wr_shape : 1);
  vector<double> iresxi(ires_size, 0.0), irespi(ires_size, 0.0);
  double *ires_arr[2] = {&iresxi[0], &irespi[0]};

  vector<double> sol(((wr_sol) ? wr_shape : 1), 0.0);	
  
  // **************************************************************
  // **************************************************************
  //                 LOOP PROGRAM OVER RESOLUTIONS
  // **************************************************************
  // **************************************************************

  int lastpt0 = lastpt; 
  int save_pt0 = save_pt;
  int nsteps0 = nsteps;
  int save_step0 = save_step;
  string outfile0 = outfile;
  double lam0 = lam;
  
  for (int factor : resolutions) {
    if (hold_const == "lambda") {
      lastpt = lastpt0 * factor;
      save_pt = save_pt0 * factor;
      nsteps = nsteps0 * factor;
      save_step = save_step0 * factor;
      outfile = to_string(factor) + "-" + outfile0;
    }
    else if (hold_const == "dt") {
      lastpt = lastpt0 * factor;
      save_pt = save_pt0 * factor;
      lam = lam0 * factor;
      outfile = to_string(factor) + "dr-" + outfile0;
    }
    else if (hold_const == "dr") {
      nsteps = nsteps0 * factor;
      save_step = save_step0 * factor;
      lam = lam0 / ((double) factor);
      outfile = to_string(factor) + "dt-" + outfile0;
    }
    else { cout << "ERROR: hold_const must be 'lambda' or 'dt' or 'dr'" << endl; }
  
    // derived parameters
    int npts = lastpt + 1;
    double dr = (rmax - rmin) / ((double) lastpt);
    double dt = lam * dr;
    double somm_coeff = 0.75*lam + 0.5*dt/rmax; // for outer bc
    
  // OUTPUT parameter data
  cout << "\noutfile name = " << outfile << "\ngrid size = " << lastpt << " (" << save_pt
       << "/write)\ntime steps = " << nsteps << " (" << save_step << "/write)\nlambda = "
       << lam << "\nr2m = " << r2m << "\nrmin = " << rmin << "\nrmax = " << rmax
       << "\ndissipation = " << dspn << "\niterative tolerance = " << tol << "\nmaximum iterations = "
       << maxit << "\nic_Dsq = " << ic_Dsq << "\nic_r0 = " << ic_r0 << "\nic_Amp = " << ic_Amp
       << "\nmass check step = " << check_step << "\nmaximum evolution time = " << nsteps*dt
       << "\ndr = " << dr << "\ndt = " << dt << endl;
  ofstream specs;
  string specs_name = outfile + ".txt";
  specs.open(specs_name, ofstream::out);
  specs << "\noutfile name = " << outfile << "\ngrid size = " << lastpt << " (" << save_pt
	<< "/write)\ntime steps = " << nsteps << " (" << save_step << "/write)\nlambda = "
	<< lam << "\nr2m = " << r2m << "\nrmin = " << rmin << "\nrmax = " << rmax
	<< "\ndissipation = " << dspn << "\niterative tolerance = " << tol << "\nmaximum iterations = "
	<< maxit << "\nic_Dsq = " << ic_Dsq << "\nic_r0 = " << ic_r0 << "\nic_Amp = " << ic_Amp
	<< "\nmass check step = " << check_step << "\nmaximum evolution time = " << nsteps*dt
	<< "\ndr = " << dr << "\ndt = " << dt << "\n\noptions:\nzero pi_0 = " << boolalpha << zero_pi
	<< "\nsommerfeld bc = " << somm_cond << "\ndissipation at bound = " << dspn_bound << endl;
  specs.close();
  
  // **********************************************************
  // ***************** OBJECT DECLARATIONS ********************
  // **********************************************************
  
  // outfiles
  string solname = "sol-" + outfile + ".sdf";
  string outfileXi_name = "Xi-" + outfile + ".sdf";
  string outfilePi_name = "Pi-" + outfile + ".sdf";
  char *name_arr[2] = {&outfileXi_name[0], &outfilePi_name[0]};
  string iresXi_name = "iresXi-" + outfile + ".sdf";
  string iresPi_name = "iresPi-" + outfile + ".sdf";
  char *iresname_arr[2] = {&iresXi_name[0], &iresPi_name[0]};  
  string mass_file = "mass-" + outfile + ".csv";
  ofstream ofs_mass;
  ofs_mass.open(mass_file, ofstream::out);
  ofs_mass << "save=" << save_step <<","<< "check=" << check_step << endl;
  string itn_file = "itns-" + outfile + ".csv";
  ofstream ofs_itn;
  if (wr_itn) { ofs_itn.open(itn_file, ofstream::out); }

  // fields and residuals
  vector<double> xi(npts, 0.0), pi(npts, 0.0);
  vector<double> old_xi(npts, 0.0), old_pi(npts, 0.0);
  vector<double> bxi(npts, 0.0), bpi(npts, 0.0);
  vector<double> resxi(npts, 0.0), respi(npts, 0.0);
  vector<double> alpha(npts, 0.0), beta(npts, 0.0), psi(npts, 0.0), ap2(npts, 0.0);
  ires_size = ((wr_ires) ? npts : 1);
  vector<double> xic1(ires_size, 0.0),  xic2(ires_size, 0.0),
    pic1(ires_size, 0.0), pic2(ires_size, 0.0),
    d1(ires_size, 0.0), d2(ires_size, 0.0);

  // *********************************************
  // **************** DEBUG **********************
  // *********************************************
  string resxi_fname = "resXi-" + outfile + ".sdf";
  string respi_fname = "resPi-" + outfile + ".sdf";
  char *resname_arr[2] = {&resxi_fname[0], &respi_fname[0]};
  int maxit_count = 0;

  time_t start_time = time(NULL); // time for rough performance measure
  
  // **********************************************************
  // ******************* INITIAL DATA ************************
  // **********************************************************

  int i, j, itn = 0; // declare loop integers
  double res = 1.0; // declare residual indicator
  double r = rmin, t = 0.0; // declare position and time variables
  for (j = 0; j < npts; ++j) {
    alpha[j] = ic_alpha(r, r2m);
    beta[j] = ic_beta(r, r2m);
    psi[j] = ic_psi(r, r2m);
    xi[j] = ic_xi(r, ic_Amp, ic_Dsq, ic_r0);
    if (!zero_pi) { pi[j] = ic_pi(r, ic_Amp, ic_Dsq, ic_r0, sq(psi[j])/alpha[j], 1 - beta[j]); }
    if ((wr_sol) && (j%save_pt == 0)) {
      sol[j/save_pt] = ic_sol(r, ic_Amp, ic_Dsq, ic_r0); }
    r += dr;
  }
  if (rmin == 0) {
    dirichlet0(xi);
    neumann0(pi);
  }

  // **********************************************************
  // ******************* TIME STEPPING ************************
  // *******************   & WRITING   ************************
  // **********************************************************
  
  gft_set_multi(); // start bbhutil file i/o
  for (i = 0; i < nsteps; ++i) {
    t = i * dt;
    // *************  WRITE the fields (& res, ires) at t(n)  **************
    if (i % save_step == 0) {
      if (wr_ires) {
	get_wr_arr_ires(xi, pi, old_xi, old_pi, alpha, beta, psi, wr_xi, wr_pi, 
			iresxi, irespi, lastwr, save_pt, lam, dr);
	wr_step(ires_arr, 2, iresname_arr, t, bbh_shape, bbh_rank, coords);
      }
      else { get_wr_arr(xi, pi, wr_xi, wr_pi, wr_shape, save_pt); }
      
      wr_step(field_arr, 2, name_arr, t, bbh_shape, bbh_rank, coords);
      
      if (wr_res) {
	get_wr_arr(resxi, respi, wr_xi, wr_pi, wr_shape, save_pt);
	wr_step(field_arr, 2, resname_arr, t, bbh_shape, bbh_rank, coords);
      }
      if (wr_sol) { gft_out_bbox(&solname[0], t, bbh_shape, bbh_rank,
				      coords, &sol[0]); }	
    }
    // now set old_xi/pi to t(n) so that xi/pi can be updated to t(n+1)
    old_xi = xi;
    old_pi = pi;

    // **********************************************************
    //         SOLVE EVOLUTION EQUATION Ax = b for x(n+1)
    // **********************************************************   
    // create rhs vec b in A_lhs.x(n+1) = A_rhs.x(n) := b
    r = rmin;
    if (rmin != 0) {
      bxi[0] = old_xi[0] + fda0_xi(old_xi, old_pi, alpha, beta, psi, lam);
      bpi[0] = old_pi[0] + fda0_pi(old_xi, old_pi, alpha, beta, psi, lam, rmin, dr);
    }
    for (j = 1; j < lastpt; ++j) {
      r += dr;
      bxi[j] = old_xi[j] + fda_xi(old_xi, old_pi, alpha, beta, psi, j, lam);
      bpi[j] = old_pi[j] + fda_pi(old_xi, old_pi, alpha, beta, psi, j, lam, r, dr);
    }
    // reset itn and set res > tol to enter GAUSS-SEIDEL ITERATIVE SOLVER
    itn = 0, res = 1.0;
    while (res > tol) {
      
      r = rmin;
      // UPDATE INTERIOR and collect residuals
      if (rmin == 0) { neumann0(pi); }
      else {
	xi[0] = bxi[0] + fda0_xi(xi, pi, alpha, beta, psi, lam);
	pi[0] = bpi[0] + fda0_pi(xi, pi, alpha, beta, psi, lam, rmin, dr);
      }
      r += dr;
      xi[1] = bxi[1] + fda_xi(xi, pi, alpha, beta, psi, 1, lam);
      pi[1] = bpi[1] + fda_pi(xi, pi, alpha, beta, psi, 1, lam, r, dr);
      for (j = 2; j < lastpt; ++j) {
        r += dr;
        xi[j] = bxi[j] + fda_xi(xi, pi, alpha, beta, psi, j, lam);
	pi[j] = bpi[j] + fda_pi(xi, pi, alpha, beta, psi, j, lam, r, dr);
	resxi[j-1] = abs(xi[j-1] - bxi[j-1] -
			 fda_xi(xi, pi, alpha, beta, psi, j-1, lam));
	respi[j-1] = abs(pi[j-1] - bpi[j-1] -
			 fda_pi(xi, pi, alpha, beta, psi, j-1, lam, r-dr, dr));
      }
      resxi[0] = abs(xi[0] - bxi[0] -
		     fda0_xi(xi, pi, alpha, beta, psi, lam));
      respi[0] = abs(pi[0] - bpi[0] -
		     fda0_pi(xi, pi, alpha, beta, psi, lam, rmin, dr));
      // UPDATE BOUNDARY
      if (somm_cond) {
	sommerfeld(xi, xi, lastpt, lam, somm_coeff);
	sommerfeld(pi, pi, lastpt, lam, somm_coeff);
      }
      
      resxi[lastpt-1] = abs(xi[lastpt-1] - bxi[lastpt-1] -
			    fda_xi(xi, pi, alpha, beta, psi, lastpt-1, lam));
      respi[lastpt-1] = abs(pi[lastpt-1] - bpi[lastpt-1] -
			    fda_pi(xi, pi, alpha, beta, psi, lastpt-1, lam, r, dr));
      
      // CHECK RESIDUAL
      res = max(*max_element(resxi.begin(), resxi.end()), *max_element(respi.begin(), respi.end())); // can also use 1-norm or 2-norm

      ++itn; 
      if (itn % maxit == 0) {
	res = 0.0;
	++maxit_count;
	if (i % 500*factor == 0) { cout << i << " res= " << res << " at " << itn << endl; }
      }
    }
    if (wr_itn) { ofs_itn << itn << endl; } // record itn count
    // ****************** ITERATIVE SOLUTION COMPLETE ******************
    
    // ****************** kreiss-oliger DISSIPATION ********************
    // at ind next to boundaries can call dissipate on ind+/-1 or ind+/-2
      if (rmin == 0) {
	xi[1] += antidiss1(dspn, old_xi);
	pi[1] += symdiss1(dspn, old_pi);
	if (dspn_bound) {
	  pi[0] += symdiss0(dspn, old_pi);
	}
      }
    for (j = 2; j < lastpt-1; ++j) {
      xi[j] += dissipate(dspn, old_xi, j);
      pi[j] += dissipate(dspn, old_pi, j);
    }

    // ****************** WRITE MASS & update field **********************
    //if (i % check_step*save_step == 0) { mass_check(xi, pi, t, dr, rmin, ofs_mass, alpha, beta); }
    if (wr_sol) { update_sol(xi, pi, alpha, beta, psi, sol, dt, save_pt, wr_shape); } 
  }
  // ******************** DONE TIME STEPPING *********************
  
  // write final time step
  if (nsteps % save_step == 0) {
    get_wr_arr(xi, pi, wr_xi, wr_pi, lastwr+1, save_pt);
    wr_step(field_arr, 2, name_arr, nsteps*dt, bbh_shape, bbh_rank, coords);
  }
  // close outfiles
  gft_close_all();
  ofs_mass.close();
  if (wr_itn) { ofs_itn.close(); }
  // print resolution runtime
  cout << difftime(time(NULL), start_time) << " seconds elapsed" << endl;
  //*************DEBUG************
  cout << maxit_count << " steps reached maxit=" << maxit << "\n" << endl;
  
  }
  // ******************** DONE LOOPING OVER RESOLUTIONS *********************
  
  return 0;
}
