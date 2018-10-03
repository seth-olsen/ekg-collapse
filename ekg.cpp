/*

ekg scalar field collapse to black hole 

compile with lapack-makefile

input parameters in terminal as:
(do NOT include .sdf extension in outfile)

 ./p1 <outfile> <lastpt> <save_pt> <nsteps> <save_step>
      <lam> <r2m> <rmin> <rmax> <dspn> <tol> <ell_tol> <maxit>
      <ic_Dsq> <ic_r0> <ic_Amp> <check_step> <zero_pi> 
      <somm_cond> <dspn_bound> <wr_ires> <wr_res> <wr_sol>
      <wr_itn> <wr_mtot> <wr_mass> <wr_xp> <wr_abp> 
      <wr_alpha> <wr_beta> <wr_psi> <hold_const>

where the value of any <parameter> can be set with the
following ordered pair of command line arguments:

 -parameter parameter_value

default values can be found at the start of main()
*/

#include "lapacke.h"
#include "ekg-proc.h"
#include "ekg-fns.h"
#include "fda-fns.h"
#include "fda-io.h"
#include <iostream>
#include <algorithm> // for max_element()
#include <fstream> // for mass/itn files
#include <ctime> // for quick time analysis
#include <cmath> // for ICs
#include <vector> // for everything
#include <string> // for parameter input
#include <map> // for parameter input
#include <cstdlib> // for atoi() and atof()
#include "bbhutil.h" // for output to .sdf

int main(int argc, char **argv)
{
  // **********************************************************
  // ******************** PARAMETERS **************************
  // **********************************************************
  
  // user-set parameters
  string outfile = "ekg";
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
  double ell_tol = tol; // will not auto change if tol changed
  int maxit = 25; // max iterations for debugging
  int ell_maxit = 2*maxit; // will not auto change if maxit set w/cmd line
  // update only accepted once it takes less that itn_tol after ell update
  int itn_tol = 2;
  int ell_itn_tol = 1; // max itns for ell soln to be accepted
  double ic_Dsq = 4.0; // gaussian width
  double ic_r0 = 50.0; // gaussian center
  double ic_Amp = 0.01; // gaussian amplitude
  int check_step = 1; // for monitoring invariant mass
  // note: set bools in command line with integers 1=true or 0=false
  bool psi_hyp = false;
  bool zero_pi = false; // zero initial time derivative?
  bool somm_cond = true; // sommerfeld condition at outer bound?
  bool dspn_bound = false; // dissipate boundary points?
  bool dr3up = false; // update pi with d/dr^3 scheme?
  bool elldiag = false; // solve elliptic equation w/o jac off-diag terms?
  bool static_metric = false; // ignore scalar field's effect on metric?
  bool wr_ires = false; // write ires? won't trigger without wr_xp
  bool wr_abpires = false; // write ires for metric vars? won't trigger without wr_abp
  bool wr_res = false; // write res? won't trigger without wr_xp
  bool wr_sol = false; // write sol?
  bool wr_itn = false; // write itn counts?
  bool wr_mtot = true; // write total mass?
  bool wr_mass = true; // write mass aspect?
  bool wr_xp = false; // write xi and pi?
  bool wr_abp = false; // write metric fields (alpha, beta, psi)?
  bool wr_alpha = true; // only triggered if wr_abp=true
  bool wr_beta = true; // only triggered if wr_abp=true
  bool wr_psi = true; // only triggered if wr_abp=true
  // variable to hold constant across resolutions
  string hold_const = "lambda"; // "lambda", "dt", or "dr"
  int nresn = 1; // 1, 2, or 3
  int resn0 = 4, resn1 = 8, resn2 = 16; // in order of priority
  int *resns[3] = {&resn0, &resn1, &resn2};

  map<string, string *> p_str {{"-outfile",&outfile},
      {"-hold_const",&hold_const}};
  map<string, int *> p_int {{"-lastpt",&lastpt}, {"-save_pt", &save_pt},
      {"-nsteps", &nsteps}, {"-save_step",&save_step}, {"-maxit",&maxit},
      {"-itn_tol",&itn_tol}, {"-ell_itn_tol",&ell_itn_tol},
      {"-check_step", &check_step}, {"-nresn", &nresn},
      {"-resn0", &resn0}, {"-resn1", &resn1}, {"-resn2", &resn2}};
  map<string, double *> p_dbl {{"-lam",&lam}, {"-r2m",&r2m}, {"-rmin",&rmin},
      {"-rmax",&rmax}, {"-dspn",&dspn}, {"-tol",&tol}, {"-ell_tol",&ell_tol},
      {"-ic_Dsq",&ic_Dsq}, {"-ic_r0",&ic_r0}, {"-ic_Amp",&ic_Amp}};
  map<string, bool *> p_bool {{"-psi_hyp",&psi_hyp}, {"-zero_pi",&zero_pi},
      {"-somm_cond",&somm_cond}, {"-dspn_bound",&dspn_bound}, {"-dr3up",&dr3up},
      {"-wr_ires",&wr_ires}, {"-wr_res",&wr_res}, {"-wr_sol",&wr_sol},
      {"-wr_itn",&wr_itn}, {"-wr_mtot",&wr_mtot}, {"-wr_mass",&wr_mass},
      {"-wr_xp",&wr_xp}, {"-wr_abp",&wr_abp}, {"-wr_abpires",&wr_abpires},
      {"-wr_alpha",&wr_alpha}, {"-wr_beta",&wr_beta}, {"-wr_psi",&wr_psi}};
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

  // set hyperbolic update scheme
  void (*hyp_update_fn)(const vector<double>&, const vector<double>&,
			vector<double>&, vector<double>&,
			vector<double>&, vector<double>&,
			const vector<double>&, const vector<double>&,
			const vector<double>&, double, double, double,
			double, int, double);
  if (dr3up) { hyp_update_fn = gs_dr3update; }
  else { hyp_update_fn = gs_update; }

  int (*ell_solver)(vector<double>&, vector<double>&,
		    const vector<double>&, const vector<double>&,
		    vector<double>&, vector<double>&,
		    vector<double>&, int, double, double,
		    int, int, int, int, int, vector<int>&,
		    int, int, double, double, int *);
  if (static_metric) { ell_solver = ell_solve_static_metric; }
  else if (elldiag) { ell_solver = ell_solve_abp_diag; }
  else { ell_solver = ell_solve_abp_full; }

  // bbhutil parameters for writing data to sdf
  int lastwr = lastpt/save_pt;
  int wr_shape = lastwr + 1;
  vector<double> wr_xi(wr_shape), wr_pi(wr_shape);
  double *field_arr[2] = {&wr_xi[0], &wr_pi[0]};
  int *bbh_shape = &wr_shape;
  int bbh_rank = 1;
  double coord_lims[2] = {rmin, rmax};
  double *coords = &coord_lims[0];

  int n_ell = ((psi_hyp) ? 2 : 3); // number of elliptic eqn fields
  int vlen = ((wr_ires) ? wr_shape : 1);
  vector<double> iresxi(vlen, 0.0), irespi(vlen, 0.0);
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
  string alphaname = "alpha-" + outfile + ".sdf";
  string betaname = "beta-" + outfile + ".sdf";
  string psiname = "psi-" + outfile + ".sdf";
  string outfileXi_name = "Xi-" + outfile + ".sdf";
  string outfilePi_name = "Pi-" + outfile + ".sdf";
  char *name_arr[2] = {&outfileXi_name[0], &outfilePi_name[0]};
  string iresXi_name = "iresXi-" + outfile + ".sdf";
  string iresPi_name = "iresPi-" + outfile + ".sdf";
  char *iresname_arr[2] = {&iresXi_name[0], &iresPi_name[0]};
  string iresalpha_name = "ires-" + alphaname;
  string iresbeta_name = "ires-" + betaname;
  string irespsi_name = "ires-" + psiname;
  char *abpiresname_arr[3] = {&iresalpha_name[0], &iresbeta_name[0], &irespsi_name[0]};
  string itn_file = "itns-" + outfile + ".csv";
  ofstream ofs_itn;
  if (wr_itn) { ofs_itn.open(itn_file, ofstream::out); }
  string maspect_name = "maspect-" + outfile + ".sdf";
  string mass_file = "mass-" + outfile + ".csv";
  ofstream ofs_mass;
  if (wr_mtot) {
    ofs_mass.open(mass_file, ofstream::out);
    ofs_mass << "save=" << save_step <<","<< "check=" << check_step << endl;
  }

  // fields and residuals
  vector<double> xi(npts, 0.0), pi(npts, 0.0);
  vector<double> old_xi(npts, 0.0), old_pi(npts, 0.0);
  vector<double> bxi(npts, 0.0), bpi(npts, 0.0);
  vector<double> resxi(npts, 0.0), respi(npts, 0.0);
  vector<double> alpha(npts, 0.0), beta(npts, 0.0), psi(npts, 0.0);
  vlen = ((psi_hyp) ? npts : 1);
  vector<double> old_psi(vlen, 0.0), bpsi(vlen, 0.0), respsi(vlen, 0.0);
  vlen = ((wr_ires) ? npts : 1);
  vector<double> xic1(vlen, 0.0),  xic2(vlen, 0.0),
    pic1(vlen, 0.0), pic2(vlen, 0.0),
    d1(vlen, 0.0), d2(vlen, 0.0);

  // *********************************************
  // **************** DEBUG **********************
  // *********************************************
  string resxi_fname = "resXi-" + outfile + ".sdf";
  string respi_fname = "resPi-" + outfile + ".sdf";
  char *resname_arr[2] = {&resxi_fname[0], &respi_fname[0]};
  int maxit_count = 0, ell_maxit_count = 0;

  // *********************************************
  // **************** LAPACK **********************
  // *********************************************
  
  // lapack object declaration
  lapack_int N = n_ell * npts;
  lapack_int kl = n_ell * 2;
  lapack_int ku = n_ell * 2;
  lapack_int nrhs = 1;
  lapack_int ldab = 2*kl + ku + 1;
  lapack_int ldb = N;
  // lapack_int info; --> now defined in ell_solver
  vector<lapack_int> ipiv(N);
  // matrices and vectors
  vector<double> jac(ldab*N, 0.0);
  vector<double> abpres(ldb, 0.0);

  time_t start_time = time(NULL); // time for rough performance measure
  
// **********************************************************
// ******************* INITIAL DATA ************************
// **********************************************************

  int i, j, itn = itn_tol + 1, ell_itn = ell_itn_tol + 1, hyp_ell_itn = 0; // declare loop integers
  double res = tol + 1; // declare residual indicators
  // double ell_res = ell_tol + 1; --> now defined in ell_solver
  double r = rmin, t = 0; // declare position and time variables
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

  // SOLVE ELLIPTIC EQUATIONS FOR t=0
  ell_itn = (*ell_solver)(jac, abpres, xi, pi, alpha, beta, psi, lastpt, dr, rmin,
			  N, kl, ku, nrhs, ldab, ipiv, ldb, ell_maxit, ell_tol, 0,
			  &ell_maxit_count);
  while (ell_itn > ell_itn_tol) {
    ell_itn = (*ell_solver)(jac, abpres, xi, pi, alpha, beta, psi, lastpt, dr, rmin,
			    N, kl, ku, nrhs, ldab, ipiv, ldb, ell_maxit, ell_tol, 0,
			    &ell_maxit_count);
  }
  
// **********************************************************
// ******************* TIME STEPPING ************************
// *******************   & WRITING   ************************
// **********************************************************
  
  gft_set_multi(); // start bbhutil file i/o
  for (i = 0; i < nsteps; ++i) {
    t = i * dt;
    // *************  WRITE the fields (& res, ires, sol) at t(n)  **************
    if (i % save_step == 0) {
      writing(lastwr, wr_shape, wr_xi, wr_pi, field_arr, bbh_shape, bbh_rank, coords,
	      wr_sol, solname, sol, wr_xp, name_arr, xi, pi, old_xi, old_pi, alpha, beta,
	      psi, iresxi, irespi, resxi, respi, wr_res, resname_arr, wr_ires, ires_arr,
	      iresname_arr, wr_abp, wr_alpha, alphaname, wr_beta, betaname,
	      wr_psi, psiname, wr_abpires, abpiresname_arr, save_pt, lam, dr, rmin, t);
    }
    // now set old_xi/pi to t(n) so that xi/pi can be updated to t(n+1)
    old_xi = xi;
    old_pi = pi;
// ******************************************************************
// ******************************************************************
//         SOLVE HYPERBOLIC & ELLIPTIC EQUATIONS ITERATIVELY
// ******************************************************************
// ******************************************************************
    
    r = rmin;
    set_rhs(bxi, bpi, old_xi, old_pi, alpha, beta, psi, lam, dr, r, lastpt);
    // reset res > tol to enter gauss-seidel solver for hyperbolic equations
    // reset ell_itn > 1 to enter direct linear solver for the elliptic equations
    // solution is accepted when elliptic equations take less than 2 updates
    hyp_ell_itn = 0;
    ell_itn = ell_itn_tol + 1; itn = itn_tol + 1;
    while (ell_itn > ell_itn_tol || itn > itn_tol) {
// ***********************************************************************
// ***************** START HYPERBOLIC ITERATIVE SOLUTION *****************
// ***********************************************************************
      itn = 0; res = tol + 1;
      while (res > tol) {
	(*hyp_update_fn)(bxi, bpi, resxi, respi, xi, pi, alpha, beta, psi, lam, dr, rmin, rmin,
			 lastpt, somm_coeff);
	// CHECK RESIDUAL
	res = max(*max_element(resxi.begin(), resxi.end()),
		  *max_element(respi.begin(), respi.end())); // can use 1-norm or 2-norm
	++itn; 
	if (itn % maxit == 0) {
	  if (i % 100*factor == 0) { cout << i << " res= " << res << " at " << itn << endl; }
	  res = 0.0;
	  ++maxit_count;
	}
      }
      if (wr_itn) { ofs_itn << itn << endl; } // record itn count
// ****************************************************************************
// ****************** HYPERBOLIC ITERATIVE SOLUTION COMPLETE ******************
// ****************************************************************************
    
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
      
// ***********************************************************************
// ****************** START ELLIPTIC ITERATIVE SOLUTION ******************
// ***********************************************************************
      ell_itn = (*ell_solver)(jac, abpres, xi, pi, alpha, beta, psi, lastpt, dr, rmin,
			      N, kl, ku, nrhs, ldab, ipiv, ldb, ell_maxit, ell_tol, t,
			      &ell_maxit_count);
      while (ell_itn > ell_itn_tol) {
	ell_itn = (*ell_solver)(jac, abpres, xi, pi, alpha, beta, psi, lastpt, dr, rmin,
				N, kl, ku, nrhs, ldab, ipiv, ldb, ell_maxit, ell_tol, t,
				&ell_maxit_count);
      }
      
// **************************************************************************
// ****************** ELLIPTIC ITERATIVE SOLUTION COMPLETE ******************
// **************************************************************************
    ++hyp_ell_itn;
    if (hyp_ell_itn > ell_maxit) {
      ell_itn = 0;
      itn = 0;
      cout << i << " hyp_ell_itn reached " << hyp_ell_itn << " at t = " << t << endl;
    }
    
    }   
// ***********************************************************************
// ***********************************************************************
// ****************** FULL ITERATIVE SOLUTION COMPLETE *******************
// ***********************************************************************
// ***********************************************************************
    
    // ****************** WRITE MASS & update field **********************
    if (i % (check_step*save_step) == 0) {
      if (wr_mass) {
	maspect[0] = mass_aspect0(alpha, beta, psi, dr, rmin);
	r = rmin;
	for (j = 1; j < lastwr; ++j) {
	  r += save_pt*dr;
	  maspect[j] = mass_aspect(alpha, beta, psi, j, dr, r);
	}
	maspect[lastwr] = mass_aspectR(alpha, beta, psi, lastwr, dr, rmax);
	gft_out_bbox(&maspect_name[0], t, bbh_shape, bbh_rank, coords, &maspect[0]);
	if (wr_mtot) { ofs_mass << i <<","<< t <<","<< maspect[lastwr] << endl; }
      }
      else if (wr_mtot) {
	ofs_mass << i <<","<< t <<","<< mass_aspectR(alpha, beta, psi, lastwr, dr, rmax) << endl;
      }
    }
    if (wr_sol) { update_sol(xi, pi, alpha, beta, psi, sol, dt, save_pt, wr_shape); } 
  }
  // ******************** DONE TIME STEPPING *********************
  
  // write final time step
  if (nsteps % save_step == 0) {
    writing(lastwr, wr_shape, wr_xi, wr_pi, field_arr, bbh_shape, bbh_rank, coords,
	    wr_sol, solname, sol, wr_xp, name_arr, xi, pi, old_xi, old_pi, alpha, beta,
	    psi, iresxi, irespi, resxi, respi, wr_res, resname_arr, wr_ires, ires_arr,
	    iresname_arr, wr_abp, wr_alpha, alphaname, wr_beta, betaname,
	    wr_psi, psiname, wr_abpires, abpiresname_arr,
	    save_pt, lam, dr, rmin, nsteps*dt);
  }
  // close outfiles
  gft_close_all();
  if (wr_itn) { ofs_itn.close(); }
  if (wr_mtot) { ofs_mass.close(); }
  // print resolution runtime
  cout << difftime(time(NULL), start_time) << " seconds elapsed" << endl;
  //*************DEBUG************
  cout << maxit_count << " hyperbolic steps reached maxit=" << maxit << endl;
  cout << ell_maxit_count << " elliptic steps reached maxit=" << ell_maxit << "\n" << endl;
  
  }
  // ******************** DONE LOOPING OVER RESOLUTIONS *********************
  
  return 0;
}
