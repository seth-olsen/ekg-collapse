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

#include "ekg-proc.h"

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
    if ((wr_mass) && (i % (check_step*save_step) == 0)) {
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
