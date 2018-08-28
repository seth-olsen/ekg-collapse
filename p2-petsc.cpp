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
#include <mpi.h>
#include <petscksp.h> // PETSc library

static char help[] = "EKG field collapse using PETSc\n";

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
// multiply by 4*M_PI*r^2 for dm/dr of scalar
inline double dmdr_scalar(double xival, double pival, double alphaval,
			  double betaval, double psival) {
  return pw4(psival)*betaval*xival*pival + 0.5*sq(psival)*alphaval*(sq(xival) + sq(pival)); }

inline double mass_aspect(const vector<double>& alpha, const vector<double>& beta,
			  const vector<double>& psi, int k, double dr, double r) {
  return r*pow(psi[k],6)*sqin(alpha[k])*sq(r*ddr_c(beta,k,dr) - beta[k])/18.0 -
    2*sq(r)*ddr_c(psi,k,dr)*(psi[k] + r*ddr_c(psi,k,dr)); }

inline double mass_aspect0(const vector<double>& alpha, const vector<double>& beta,
			  const vector<double>& psi, double dr, double r) {
  return r*pow(psi[0],6)*sqin(alpha[0])*sq(r*ddr_f(beta,0,dr) - beta[0])/18.0 -
    2*sq(r)*ddr_f(psi,0,dr)*(psi[0] + r*ddr_f(psi,0,dr)); }

inline double mass_aspectR(const vector<double>& alpha, const vector<double>& beta,
			  const vector<double>& psi, int k, double dr, double r) {
  return r*pow(psi[k],6)*sqin(alpha[k])*sq(r*ddr_b(beta,k,dr) - beta[k])/18.0 -
    2*sq(r)*ddr_b(psi,k,dr)*(psi[k] + r*ddr_b(psi,k,dr)); }
    
// compute and write mass
void mass_check(const vector<double>& xi, const vector<double>& pi,
		const vector<double>& alpha, const vector<double>& beta,
		const vector<double>& psi, double dr, double rmin,
		double t, ofstream& out_stream)
{
  double mass = 0.0, rval = rmin;
  int k = 0;
  for (auto xik : xi) {
    mass += 4*M_PI*rval*rval*dr*dmdr_scalar(xik, pi[k], alpha[k], beta[k], psi[k]);
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
		     int lastwrite, int savept, double lam, double dr, double rmin)
{
  double rval = rmin;
  wrxi[0] = xi[0];
  wrpi[0] = pi[0];
  if (rmin != 0) {
    iresxi[0] = xi[0] - iresxi_f(xi, pi, alpha, beta, psi, 0, lam)
      - oldxi[0] - iresxi_f(oldxi, oldpi, alpha, beta, psi, 0, lam);
    irespi[0] = pi[0] - irespi_f(xi, pi, alpha, beta, psi, 0, lam, dr, rval)
    - oldpi[0] - irespi_f(oldxi, oldpi, alpha, beta, psi, 0, lam, dr, rval);
  }
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

// **********************************************************
// **********************************************************
//                    FDA IMPLEMENTATIONS
// **********************************************************
// **********************************************************

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

// field = oldfield + fda(field) + fda(oldfield)
inline double fda_xi(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double lam)
{
  return ( 0.25*lam* dx_ekg_c(beta, xi, alpha, pi, psi, ind) );
}

inline double fda_pi(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double lam, double dr, double r)
{
  return ( 0.25*lam*p4in(psi[ind])* sqin(r)*dp_ekg_c(beta, pi, alpha, xi, psi, ind, dr, r)
	   - (lam/3.0)* pi[ind]*( d_c(beta, ind) + beta[ind]*(6*d_c(psi, ind)/psi[ind]
							     + 4*dr/r) ) );
  // or dpi_c(beta, pi, alpha, xi, psi, ind, dr, r)*sqin(r) replaced by
  // d3pi_c(beta, pi, alpha, xi, psi, ind, dr, r)
}
// only need these if no BC at rmin
inline double fda0_xi(const vector<double>& xi, const vector<double>& pi,
		      const vector<double>& alpha, const vector<double>& beta,
		      const vector<double>& psi, double lam)
{
  return ( 0.25 * lam * dx_ekg_f(beta, xi, alpha, pi, psi, 0) );
}

inline double fda0_pi(const vector<double>& xi, const vector<double>& pi,
		      const vector<double>& alpha, const vector<double>& beta,
		      const vector<double>& psi, double lam, double dr, double r)
{
  return ( 0.25*lam*p4in(psi[0])* sqin(r)*dp_ekg_f(beta, pi, alpha, xi, psi, 0, dr, r)
	   - (lam/3.0)* pi[0]*( d_f(beta, 0) + beta[0]*(6*d_f(psi, 0)/psi[0]
							     + 4*dr/r) ) );
  // or dpi_f(beta, pi, alpha, xi, psi, ind, dr, r)*sqin(r) replaced by
  // d3pi_f(beta, pi, alpha, xi, psi, ind, dr, r)
}

inline double fda_psi(const vector<double>& xi, const vector<double>& pi,
		      const vector<double>& alpha, const vector<double>& beta,
		      const vector<double>& psi, int ind, double lam, double dr, double r) {
  return 0.25*lam*(beta[ind]*(psi[ind+1] - psi[ind-1]) +
		   psi[ind]*(beta[ind+1] +4*beta[ind]/r - beta[ind-1])); }

inline double fda_respsi(const vector<double>& xi, const vector<double>& pi,
			 const vector<double>& alpha, const vector<double>& beta,
			 const vector<double>& psi, int ind, double dr, double r) {
  return ddr2_c(psi,ind,dr) + 2*ddr_c(psi,ind,dr)/r + sqin(alpha[ind])*
    (ddr_c(beta,ind,dr) - beta[ind]/r)*psi[ind]*pw4(psi[ind])/12.0 + 
    M_PI*(sq(xi[ind]) + sq(pi[ind]))*psi[ind] ; }

inline double fda_resbeta(const vector<double>& xi, const vector<double>& pi,
			  const vector<double>& alpha, const vector<double>& beta,
			  const vector<double>& psi, int ind, double dr, double r) {
  return  ddr2_c(beta,ind,dr) + 12*M_PI*sqin(psi[ind])*alpha[ind]*xi[ind]*pi[ind]
    + (2/r + 6*ddr_c(psi,ind,dr)/psi[ind] - ddr_c(alpha,ind,dr)/alpha[ind])*
    (ddr_c(beta,ind,dr) - beta[ind]/r); }

inline double fda_resalpha(const vector<double>& xi, const vector<double>& pi,
			   const vector<double>& alpha, const vector<double>& beta,
			   const vector<double>& psi, int ind, double dr, double r) {
  return ddr2_c(alpha,ind,dr) + 2*(1/r + ddr_c(psi,ind,dr)/psi[ind])*ddr_c(alpha,ind,dr)
    - 2*pw4(psi[ind])*sq(ddr_c(beta,ind,dr) - beta[ind]/r)/(3*alpha[ind])
    - 8*M_PI*alpha[ind]*sq(pi[ind]) ; }

inline double fdaR_respsi(const vector<double>& psi, int ind, double dr, double r) {
  return d_b(psi,ind) + 2*dr*(psi[ind] - 1)/r; }

inline double fdaR_resbeta(const vector<double>& beta, int ind, double dr, double r) {
  return d_b(beta,ind) + 2*dr*beta[ind]/r; }

inline double fdaR_resalpha(const vector<double>& alpha, int ind, double dr, double r) {
  return d_b(alpha,ind) + 2*dr*(alpha[ind] - 1)/r; }

// ***********************  JACOBIAN  ***********************

// ***********************  row alpha(ind)   ***********************

inline double jac_aa(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -2*sqin(dr) - 8*M_PI*sq(pi[ind]) +
    2*pw4(psi[ind])*sqin(alpha[ind])*sq(ddr_c(beta,ind,dr) - beta[ind]/r); }

inline double jac_aa_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return sqin(dr) + p_m*(ddr_c(psi,ind,dr)/psi[ind] + 1/r)/dr; }

inline double jac_ab(const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return 4*pw4(psi[ind])*(ddr_c(beta,ind,dr) - beta[ind]/r) / (3*alpha[ind]*r); }

inline double jac_ab_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return -p_m*2*pw4(psi[ind])*(ddr_c(beta,ind,dr) - beta[ind]/r) / (3*alpha[ind]*dr); }

inline double jac_ap(const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -2*(ddr_c(alpha,ind,dr)*ddr_c(psi,ind,dr)*sqin(psi[ind]) +
	     4*pow(psi[ind],3)*sq(ddr_c(beta,ind,dr) - beta[ind]/r) / (3*alpha[ind])); }

inline double jac_ap_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return p_m*ddr_c(alpha,ind, dr) / (dr*psi[ind]); }


// ***********************  row beta(ind)   ***********************

inline double jac_ba(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return 12*M_PI*xi[ind]*pi[ind]*sqin(psi[ind]) +
    ddr_c(alpha,ind,dr)*sqin(alpha[ind])*(ddr_c(beta,ind,dr) - beta[ind]/r); }

inline double jac_ba_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return p_m*0.5*(ddr_c(beta,ind,dr) - beta[ind]/r) / (dr*alpha[ind]); }

inline double jac_bb(const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -2*sqin(dr) - (2/r - ddr_c(alpha,ind,dr)/alpha[ind] + 6*ddr_c(psi,ind,dr))/dr; }

inline double jac_bb_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return sqin(dr) + p_m*0.5*(2/r - ddr_c(alpha,ind,dr)/alpha[ind] + 6*ddr_c(psi,ind,dr))/dr; }

inline double jac_bp(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -24*M_PI*xi[ind]*pi[ind]*alpha[ind]/pow(psi[ind],3) -
    6*ddr_c(psi,ind,dr)*sqin(psi[ind])*(ddr_c(beta,ind,dr) - beta[ind]/r); }

inline double jac_bp_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return p_m*3*(ddr_c(beta,ind,dr) - beta[ind]/r) / (dr*psi[ind]); }

// ***********************  row psi(ind)   ***********************

inline double jac_pa(const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -pow(psi[ind],5)*(ddr_c(beta,ind,dr) - beta[ind]/r) / (6*pow(alpha[ind],3)); }

inline double jac_pa_pm() { return 0; }

inline double jac_pb(const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -pow(psi[ind],5)*sqin(alpha[ind]) / (12*r); }

inline double jac_pb_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return p_m*pow(psi[ind],5)*sqin(alpha[ind]) / (24*dr); }

inline double jac_pp(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -2*sqin(dr) + M_PI*(sq(xi[ind]) + sq(pi[ind])) +
    5*pw4(psi[ind])*sqin(alpha[ind])*(ddr_c(beta,ind,dr) - beta[ind]/r) / 12.0; }

inline double jac_pp_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return sqin(dr) + p_m/(r*dr); }


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
  double r2m = 0.0;
  double rmin = 0.0;
  double rmax = 100.0;
  double dspn = 0.5; // dissipation coefficient
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
  bool wr_mass = false; // write mass aspect?
  // variable to hold constant across resolutions
  string hold_const = "lambda"; // "lambda", "dt", or "dr"
  int nresn = 1; // 1, 2, or 3
  int resn0 = 8, resn1 = 16, resn2 = 32; // in order of priority
  int *resns[3] = {&resn0, &resn1, &resn2};

  map<string, string *> p_str {{"-outfile",&outfile},
      {"-hold_const",&hold_const}};
  map<string, int *> p_int {{"-lastpt",&lastpt}, {"-save_pt", &save_pt},
      {"-nsteps", &nsteps}, {"-save_step",&save_step}, {"-maxit",&maxit},
      {"-check_step", &check_step}, {"-nresn", &nresn},
      {"-resn0", &resn0}, {"-resn1", &resn1}, {"-resn2", &resn2}};
  map<string, double *> p_dbl {{"-lam",&lam}, {"-r2m",&r2m}, {"-rmin",&rmin},
      {"-rmax",&rmax}, {"-dspn",&dspn}, {"-tol",&tol}, {"-ic_Dsq",&ic_Dsq},
      {"-ic_r0",&ic_r0}, {"-ic_Amp",&ic_Amp}};
  map<string, bool *> p_bool {{"-zero_pi",&zero_pi},
      {"-somm_cond",&somm_cond}, {"-dspn_bound",&dspn_bound},
      {"-wr_ires",&wr_ires}, {"-wr_res",&wr_res},
      {"-wr_sol",&wr_sol}, {"-wr_itn",&wr_itn}, {"-wr_mass",&wr_mass}};
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
  vector<int> resolutions(nresn);
  for (int k = 0; k < nresn; ++k) {
    resolutions[k] = *resns[k];
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
  vector<double> maspect(((wr_mass) ? wr_shape : 1), 0.0);
  
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
    cout << param_print(outfile,lastpt,save_pt,nsteps,save_step,lam,r2m,rmin,rmax,
			dspn,tol,maxit,ic_Dsq,ic_r0,ic_Amp,check_step,dr,dt,
			zero_pi,somm_cond,dspn_bound);
    // *** BACK TO MAIN() INDENT FOR SPACE ***   
  ofstream specs;
  string specs_name = outfile + ".txt";
  specs.open(specs_name, ofstream::out);
  specs << param_print(outfile,lastpt,save_pt,nsteps,save_step,lam,r2m,rmin,rmax,
			dspn,tol,maxit,ic_Dsq,ic_r0,ic_Amp,check_step,dr,dt,
			zero_pi,somm_cond,dspn_bound);
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
  string maspect_name = "maspect-" + outfile + ".sdf";

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
  double res = tol + 1.0; // declare residual indicator
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
			iresxi, irespi, lastwr, save_pt, lam, dr, rmin);
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
      bpi[0] = old_pi[0] + fda0_pi(old_xi, old_pi, alpha, beta, psi, lam, dr, rmin);
    }
    for (j = 1; j < lastpt; ++j) {
      r += dr;
      bxi[j] = old_xi[j] + fda_xi(old_xi, old_pi, alpha, beta, psi, j, lam);
      bpi[j] = old_pi[j] + fda_pi(old_xi, old_pi, alpha, beta, psi, j, lam, dr, r);
    }
    // reset itn and set res > tol to enter GAUSS-SEIDEL ITERATIVE SOLVER
    itn = 0, res = tol + 1.0;
    while (res > tol) {
      
      r = rmin;
      // UPDATE INTERIOR and collect residuals
      if (rmin == 0) { neumann0(pi); }
      else {
	xi[0] = bxi[0] + fda0_xi(xi, pi, alpha, beta, psi, lam);
	pi[0] = bpi[0] + fda0_pi(xi, pi, alpha, beta, psi, lam, dr, rmin);
      }
      r += dr;
      xi[1] = bxi[1] + fda_xi(xi, pi, alpha, beta, psi, 1, lam);
      pi[1] = bpi[1] + fda_pi(xi, pi, alpha, beta, psi, 1, lam, dr, r);
      for (j = 2; j < lastpt; ++j) {
        r += dr;
        xi[j] = bxi[j] + fda_xi(xi, pi, alpha, beta, psi, j, lam);
	pi[j] = bpi[j] + fda_pi(xi, pi, alpha, beta, psi, j, lam, dr, r);
	resxi[j-1] = abs(xi[j-1] - bxi[j-1] -
			 fda_xi(xi, pi, alpha, beta, psi, j-1, lam));
	respi[j-1] = abs(pi[j-1] - bpi[j-1] -
			 fda_pi(xi, pi, alpha, beta, psi, j-1, lam, dr, r-dr));
      }
      resxi[0] = abs(xi[0] - bxi[0] -
		     fda0_xi(xi, pi, alpha, beta, psi, lam));
      respi[0] = abs(pi[0] - bpi[0] -
		     fda0_pi(xi, pi, alpha, beta, psi, lam, dr, rmin));
      // UPDATE BOUNDARY
      if (somm_cond) {
	sommerfeld(xi, xi, lastpt, lam, somm_coeff);
	sommerfeld(pi, pi, lastpt, lam, somm_coeff);
      }
      
      resxi[lastpt-1] = abs(xi[lastpt-1] - bxi[lastpt-1] -
			    fda_xi(xi, pi, alpha, beta, psi, lastpt-1, lam));
      respi[lastpt-1] = abs(pi[lastpt-1] - bpi[lastpt-1] -
			    fda_pi(xi, pi, alpha, beta, psi, lastpt-1, lam, dr, r));
      
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
    if (wr_mass && i % check_step*save_step == 0) {
      //mass_check(xi, pi, alpha, beta, psi, dr, rmin, t, ofs_mass); }
      maspect[0] = mass_aspect0(alpha, beta, psi, dr, rmin);
      for (j = 1; j < lastwr; ++j) {
	maspect[j] = mass_aspect(alpha, beta, psi, j, dr, rmin + j*save_pt*dr); }
      maspect[lastwr] = mass_aspectR(alpha, beta, psi, lastwr, dr, rmax);
      gft_out_bbox(&maspect_name[0], t, bbh_shape, bbh_rank, coords, &maspect[0]);
    }
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

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ***************************************************************************

inline void set_inds(int *cols, int indphi) {
  cols[0] = indphi-2, cols[1] = indphi-1, cols[2] = indphi,
    cols[3] = indphi+1, cols[4] = indphi+2, cols[5] = indphi+3;
  return; }

inline void set_phivals(double *vals, int rind, double dr, double rmin,
			double r2m, double lam) {
  vals[0] = -cbeta(r(rind-1,dr,rmin), r2m, lam),
    vals[1] = -cf(r(rind-1,dr,rmin), r2m, lam),
    vals[2] = 0, vals[3] = 0,
    vals[4] = cbeta(r(rind+1,dr,rmin),r2m,lam),
    vals[5] = cf(r(rind+1,dr,rmin),r2m,lam);
  return; }

inline void set_pivals(double *vals, int rind, double dr, double rmin,
			double r2m, double lam) {
  vals[0] = -rr2(rind-1,rind,dr,rmin)*cf(r(rind-1,dr,rmin), r2m, lam),
    vals[1] = -rr2(rind-1,rind,dr,rmin)*cbeta(r(rind-1,dr,rmin), r2m, lam),
    vals[2] = 0, vals[3] = 0,
    vals[4] = rr2(rind+1,rind,dr,rmin)*cf(r(rind+1,dr,rmin),r2m,lam),
    vals[5] = rr2(rind+1,rind,dr,rmin)*cbeta(r(rind+1,dr,rmin),r2m,lam);
  return; }

inline void dspn_inds(int *cols, int ind) {
  cols[0] = ind-4, cols[1] = ind-2, cols[2] = ind,
    cols[3] = ind+2, cols[4] = ind+4;
  return; }

// ***************************************************************************
// PETSC STUFF
// ***************************************************************************
void petsc_part() {

   // derived parameters
  int npts = lastpt + 1;
  int lastwrite = lastpt / save_pt;
  int write_shape = lastwrite + 1;
  double dr = (rmax - rmin) / ((double) lastpt);
  double dt = lam * dr;
  double cN = 0.75*lam + 0.5*dt/rmax; // for outer bc
  double cdiss = -dspn * 0.0625;
  double drin = 1.0 / dr;

  // for mass monitoring
  vector<double> mass_data, mass_times;
  
  time_t start_time = time(NULL);

  // petsc object declaration
  PetscMPIInt rank;
  PetscMPIInt size;
  PetscErrorCode ierr;
  PetscInt N = 3 * npts; // size of vectors
  PetscInt last_row = N - 1;
  PetscInt diag_per_row = 1, offd_per_row = 8;
  PetscInt i, rstart, rend; 
  // matrices and vectors
  Mat A_lhs, A_rhs, dspn_mat;
  Vec abp, old_abp, abp_res;
  KSP ksp;
  PC pc;
  // checking matrix
  //PetscViewer lhs_viewer, rhs_viewer;

  // BOUNDARY CONDITIONS
  PetscInt col_inds0[3] = {0, 3, 6}, col_val1[1] = {1}, col_inds2[3] = {2, 5, 8};
  PetscInt col_indsNm3[3] = {N-9, N-6, N-3}, col_indsNm2[3] = {N-8, N-5, N-2}, col_indsNm1[3] = {N-7, N-4, N-1};
  PetscScalar mat_vals02[6] = {-3, 4, -1};
  PetscScalar mat_valsNm123[3] = {1, -4, 2*dr/rmax + 3};
  PetscScalar dspn_vals[5] = {cdiss, -4*cdiss, 6*cdiss, -4*cdiss, cdiss};
  
  //start petsc/mpi
  ierr = PetscInitialize(&argc, &argv, (char *)0, help);
  if (ierr) { return ierr; }
  MPI_Comm comm = PETSC_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // create vectors
  ierr = VecCreate(PETSC_DECIDE, N, &abp); CHKERRQ(ierr);
  ierr = VecDuplicate(abp, &old_abp); CHKERRQ(ierr);
  ierr = VecDuplicate(abp, &dspn_vec); CHKERRQ(ierr);
  // set initial values
  ierr = VecGetOwnershipRange(abp, &rstart, &rend); CHKERRQ(ierr);
  PetscInt n = rend - rstart;
  PetscInt indices[n];
  PetscScalar ic_abp[n];
  int j = 0;
  for (i = rstart; i < rend; ++i) {
    indices[j] = i;
    ic_abp[j] = (i % 2 == 0) ?
      ic_phi(r(i/2,dr,rmin), ic_Amp, ic_Dsq, ic_r0) :
      ic_pi(r((i-1)/2,dr,rmin), ic_Amp, ic_Dsq, ic_r0, r2m);
    ++j;
  }
  ierr = VecSetValues(abp, n, &indices[0], &ic_abp[0], INSERT_VALUES); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(abp); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(abp); CHKERRQ(ierr);
  //TEST ierr = VecView(uj, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  // create matrices
  int nlocal = n, nglobal;
  int *nrows = NULL;
  // get global min of rows owned by single process
  if (rank == 0) { nrows = (int *) malloc(sizeof(int) * size); }
  MPI_Gather(&n, 1, MPI_INT, nrows, 1, MPI_INT, 0, comm);
  if (rank == 0) { nlocal = *min_element(nrows, nrows + size); }
  MPI_Scatter(&nlocal, 1, MPI_INT, &nglobal, 1, MPI_INT, 0, comm);
  free(nrows);
  // set max number of nonzero off-diagonal matrix elements per row
  offd_per_row = max(3, 6-nglobal);

  ierr = MatCreateAIJ(comm, PETSC_DECIDE, PETSC_DECIDE, N, N, diag_per_row, NULL, offd_per_row, NULL, &A_rhs); CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(A_rhs, &rstart, &rend); CHKERRQ(ierr);
  n = rend - rstart;
  // boundary cases
  int on_start = 0, off_end = 0;
  if (rstart == 0) {
    ++on_start;
    ierr = MatSetValues(A_rhs, 1, &rstart, 6, col_inds01, mat_vals0, INSERT_VALUES); CHKERRQ(ierr);
  }
  if (rstart < 2 && rend > 1) {
    ++on_start;
    int pt = 1;
    ierr = MatSetValues(A_rhs, 1, &pt, 6, col_inds01, mat_vals1, INSERT_VALUES); CHKERRQ(ierr);
  }
  if (rend == N) {
    ++off_end;
    ierr = MatSetValues(A_rhs, 1, &last_row, 3, col_indsNm1, mat_valsNm12, INSERT_VALUES); CHKERRQ(ierr);
  }
  if (rstart < N-1 && rend > N-2) {
    ++off_end;
    int pt = N-2;
    ierr = MatSetValues(A_rhs, 1, &pt, 3, col_indsNm2, mat_valsNm12, INSERT_VALUES); CHKERRQ(ierr);
  }
  // interior ****MIGHT NEED TO INSERT ZEROS TO PREVENT COMPRESSION****
  PetscScalar mat_vals[6];
  PetscInt col_inds[6];
  int i_start = rstart + on_start, i_end = rend - off_end;
  for (i = i_start; i < i_end; ++i) {
    if (i % 2 == 0) {
      set_inds(col_inds, i);
      set_phivals(mat_vals, i/2, dr, rmin, r2m, lam);
    }
    else {
      set_inds(col_inds, i-1);
      set_phivals(mat_vals, (i-1)/2, dr, rmin, r2m, lam);
    }
    ierr = MatSetValues(A_rhs, 1, &i, 6, col_inds, mat_vals, INSERT_VALUES); CHKERRQ(ierr);
  }
  // matrix assembly
  ierr = MatAssemblyBegin(A_rhs, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A_rhs, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  
  // now have A, need A_rhs = I+A and A_lhs = I-A
  ierr = MatDuplicate(A_rhs, MAT_COPY_VALUES, &A_lhs); CHKERRQ(ierr);
  ierr = MatScale(A_lhs, -1.0); CHKERRQ(ierr);
  ierr = MatShift(A_lhs, 1.0); CHKERRQ(ierr);
  ierr = MatShift(A_rhs, 1.0); CHKERRQ(ierr);

  ierr = MatCreateAIJ(comm, PETSC_DECIDE, PETSC_DECIDE, N, N, 5, NULL, 4, NULL, &dspn_mat); CHKERRQ(ierr);
  PetscInt dspn_cols[5];
  for (i = 4; i < N-4; ++i) {
    dspn_inds(dspn_cols, i);
    ierr = MatSetValues(dspn_mat, 1, &i, 5, dspn_cols, dspn_vals, INSERT_VALUES); CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(dspn_mat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(dspn_mat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr); 

  // now check matrices
  /*
  ierr = PetscViewerASCIIOpen(comm, "matrix-lhs.m", &lhs_viewer); CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(lhs_viewer, PETSC_VIEWER_ASCII_MATHEMATICA); CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(comm, "matrix-rhs.m", &rhs_viewer); CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(rhs_viewer, PETSC_VIEWER_ASCII_MATHEMATICA); CHKERRQ(ierr);
  ierr = MatView(A_lhs, lhs_viewer); CHKERRQ(ierr);
  ierr = MatView(A_rhs, rhs_viewer); CHKERRQ(ierr);
  ierr = PetscViewerPopFormat(lhs_viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&lhs_viewer); CHKERRQ(ierr);
  ierr = PetscViewerPopFormat(rhs_viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&rhs_viewer); CHKERRQ(ierr);
  */

  // create linear solver context
  ierr = KSPCreate(comm, &ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, A_lhs, A_lhs); CHKERRQ(ierr);
  // block jacobi with ILU(0) is default
  //ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
  //ierr = PCSetType(pc,PCJACOBI); CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp, 1e-9, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(ierr);
  //ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERRQ(ierr);

  // bbhutil parameters for writing data to sdf
  int bbh_rank = 1;
  int *bbh_shape = &npts;
  double coord_lims[2] = {r2m, rmax};
  double *coords = &coord_lims[0];
  
  // creat context to gather abp to root for checking mass and writing
  ierr = VecScatterCreateToZero(abp, &ctx_root, &root_vec); CHKERRQ(ierr);
  vector<double> phi((rank == 0) ? write_shape : 1);
  vector<double> pi((rank == 0) ? write_shape : 1);
  double *field_arr[2] = {&phi[0], &pi[0]};
  string outfilePhi_name = "Phi-" + outfile + ".sdf";
  string outfilePi_name = "Pi-" + outfile + ".sdf";
  char *name_arr[2] = {&outfilePhi_name[0], &outfilePi_name[0]};

  ierr = VecScatterBegin(ctx_root, abp, root_vec, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(ctx_root, abp, root_vec, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
  // write t=0 data
  if (rank == 0) {
    ierr = VecGetArray(root_vec, &root_array); CHKERRQ(ierr);
    get_wr(root_array, phi, pi, write_shape, save_pt);
    ierr = VecRestoreArray(root_vec, &root_array); CHKERRQ(ierr);
    gft_set_multi();
    write_step(field_arr, 2, name_arr, 0, bbh_shape, bbh_rank, coords);
  }

  // evolve
  int step_lim = nsteps + 1;
  for (i = 1; i < step_lim; ++i) {
    // time stepping
    ierr = MatMult(A_rhs, abp, old_abp); CHKERRQ(ierr);
    ierr = KSPSolve(ksp, old_abp, abp); CHKERRQ(ierr);
    // dissipate
    ierr = MatMult(dspn_mat, old_abp, dspn_vec); CHKERRQ(ierr);
    ierr = VecAXPY(abp, 1, dspn_vec; CHKERRQ(ierr);
    // writing
    if (i % save_step == 0) {
      ierr = VecScatterBegin(ctx_root, abp, root_vec, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecScatterEnd(ctx_root, abp, root_vec, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
      if (rank == 0) {
	ierr = VecGetArray(root_vec, &root_array); CHKERRQ(ierr);
	get_wr(root_array, phi, pi, write_shape, save_pt);
	if (i % save_step*check_step == 0) {
	  double mass = 0;
	  for (j = 0; j < npts; ++j) {
	    mass += mass_integrand(r(j,dr,rmin), r2m, dr,
				   root_array[2*j], root_array[2*j+1]); }
	  mass_times.push_back(i*dt);
	  mass_data.push_back(mass);
	}
	ierr = VecRestoreArray(root_vec, &root_array); CHKERRQ(ierr);
	write_step(field_arr, 2, name_arr, i*dt, bbh_shape, bbh_rank, coords);
      }
    }
  }

  // close and destroy objects then finalize mpi/petsc
  ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
  //ierr = PCDestroy(&pc); CHKERRQ(ierr);
  ierr = MatDestroy(&A_rhs); CHKERRQ(ierr);
  ierr = MatDestroy(&A_lhs); CHKERRQ(ierr);
  //ierr = PetscViewerDestroy(&lhs_viewer); CHKERRQ(ierr);
  //ierr = PetscViewerDestroy(&rhs_viewer); CHKERRQ(ierr);
  ierr = VecDestroy(&old_abp); CHKERRQ(ierr);
  ierr = VecDestroy(&abp); CHKERRQ(ierr);
  ierr = VecDestroy(&root_vec); CHKERRQ(ierr);  
  ierr = VecScatterDestroy(&ctx_root); CHKERRQ(ierr);
  }
  return;
}
