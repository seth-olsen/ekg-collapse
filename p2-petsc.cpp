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

#include "ekg-fns.h"
#include <petscksp.h> // PETSc library

static char help[] = "EKG field collapse using PETSc\n";

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
  double ell_tol = tol; // will not auto change if tol changed
  double petsc_rtol = 0.000000001; // will not auto change if tol changed
  int maxit = 25; // max iterations for debugging
  int ell_maxit = 2*maxit; // will not auto change if maxit changed
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
  int resn0 = 4, resn1 = 8, resn2 = 16; // in order of priority
  int *resns[3] = {&resn0, &resn1, &resn2};

  map<string, string *> p_str {{"-outfile",&outfile},
      {"-hold_const",&hold_const}};
  map<string, int *> p_int {{"-lastpt",&lastpt}, {"-save_pt", &save_pt},
      {"-nsteps", &nsteps}, {"-save_step",&save_step}, {"-maxit",&maxit},
      {"-check_step", &check_step}, {"-nresn", &nresn},
      {"-resn0", &resn0}, {"-resn1", &resn1}, {"-resn2", &resn2}};
  map<string, double *> p_dbl {{"-lam",&lam}, {"-r2m",&r2m}, {"-rmin",&rmin},
      {"-rmax",&rmax}, {"-dspn",&dspn}, {"-tol",&tol}, {"-ell_tol",&ell_tol},
      {"-petsc_rtol",&petsc_rtol},
      {"-ic_Dsq",&ic_Dsq}, {"-ic_r0",&ic_r0}, {"-ic_Amp",&ic_Amp}};
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

  // initialize petsc/mpi
  PetscMPIInt rank;
  PetscMPIInt size;
  PetscErrorCode ierr;
  PetscInt nz = 9, petsc_itn;
  PetscScalar mat_vals[nz];
  PetscInt col_inds[nz];
  ierr = PetscInitialize(&argc, &argv, (char *)0, help);
  if (ierr) { return ierr; }
  MPI_Comm comm = PETSC_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  
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
  int maxit_count = 0, ell_maxit_count = 0;

  // *********************************************
  // **************** PETSC **********************
  // *********************************************

  // petsc object declaration  
  PetscInt N = 3 * npts; // size of vectors
  // matrices and vectors
  Mat jac;
  Vec abpres;
  KSP ksp;
  PC pc;
  // checking matrix
  //PetscViewer viewer;
  PetscInt indices[N];
  PetscScalar res_vals[N];
  PetscScalar *deltas;

  time_t start_time = time(NULL); // time for rough performance measure
  
// **********************************************************
// ******************* INITIAL DATA ************************
// **********************************************************

  int i, j, itn = 0, ell_itn = 2, hyp_ell_itn = 0; // declare loop integers
  double res = tol + 1.0, ell_res = ell_tol + 1;// declare residual indicators
  double r = rmin, t = 0.0; // declare position and time variables
  for (j = 0; j < npts; ++j) {
    alpha[j] = ic_alpha(r, r2m);
    beta[j] = ic_beta(r, r2m);
    psi[j] = ic_psi(r, r2m);
    xi[j] = ic_xi(r, ic_Amp, ic_Dsq, ic_r0);
    if (!zero_pi) { pi[j] = ic_pi(r, ic_Amp, ic_Dsq, ic_r0, sq(psi[j])/alpha[j], 1 - beta[j]); }
    if ((wr_sol) && (j%save_pt == 0)) {
      sol[j/save_pt] = ic_sol(r, ic_Amp, ic_Dsq, ic_r0); }
    indices[3*j] = 3*j, indices[3*j+1] = 3*j+1, indices[3*j+2] = 3*j+2;
    r += dr;
  }
  if (rmin == 0) {
    dirichlet0(xi);
    neumann0(pi);
  }

  // create petsc vector and matrix
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &abpres); CHKERRQ(ierr);
  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, N, N, nz, NULL, &jac); CHKERRQ(ierr);
  if (size == 1) {
    // BOUNDARY CONDITIONS
    PetscInt col_inds0[3] = {0, 3, 6}, col_inds1[1] = {1}, col_inds2[3] = {2, 5, 8};
    PetscInt col_indsNm3[3] = {N-9, N-6, N-3}, col_indsNm2[3] = {N-8, N-5, N-2}, col_indsNm1[3] = {N-7, N-4, N-1};
    PetscScalar mat_vals02[3] = {-3, 4, -1}, mat_vals1[1] = {1};
    PetscScalar mat_valsNm123[3] = {1, -4, 2*dr/rmax + 3};
    ierr = MatSetValues(jac, 1, &indices[0], 3, col_inds0, mat_vals02, INSERT_VALUES);
    CHKERRQ(ierr);
    ierr = MatSetValues(jac, 1, &indices[1], 1, col_inds1, mat_vals1, INSERT_VALUES);
    CHKERRQ(ierr);
    ierr = MatSetValues(jac, 1, &indices[2], 3, col_inds2, mat_vals02, INSERT_VALUES);
    CHKERRQ(ierr);
    ierr = MatSetValues(jac, 1, &indices[N-3], 3, col_indsNm3, mat_valsNm123, INSERT_VALUES);
    CHKERRQ(ierr);
    ierr = MatSetValues(jac, 1, &indices[N-2], 3, col_indsNm2, mat_valsNm123, INSERT_VALUES);
    CHKERRQ(ierr);
    ierr = MatSetValues(jac, 1, &indices[N-1], 3, col_indsNm1, mat_valsNm123, INSERT_VALUES);
    CHKERRQ(ierr);
  }
  else { cout << "\n\n\n\n****ERROR: NOT USING SINGLE PROCESSOR****\n\n\n\n" << endl; }

  // create linear solver context
  ierr = KSPCreate(comm, &ksp); CHKERRQ(ierr);
  ierr = KSPSetType(ksp, KSPPREONLY); CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCLU); CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp, petsc_rtol, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(ierr);
  //ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
  
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
    ell_itn = 2;
    while (ell_itn > 1) {
// ***********************************************************************
// ***************** START HYPERBOLIC ITERATIVE SOLUTION *****************
// ***********************************************************************
      itn = 0, res = tol + 1;
      while (res > tol) {
	
	r = rmin;
	gs_update(bxi, bpi, resxi, respi, xi, pi, alpha, beta, psi, lam, dr, r, rmin,
		  lastpt, somm_coeff);
	// CHECK RESIDUAL
	res = max(*max_element(resxi.begin(), resxi.end()),
		  *max_element(respi.begin(), respi.end())); // can use 1-norm or 2-norm
	++itn; 
	if (itn % maxit == 0) {
	  res = 0.0;
	  ++maxit_count;
	  if (i % 500*factor == 0) { cout << i << " res= " << res << " at " << itn << endl; }
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
      // get initial residual
      r = rmin;
      get_ell_res(res_vals, xi, pi, alpha, beta, psi, lastpt, N, dr, r);
      ell_res = max(*max_element(res_vals, res_vals+N),
		    abs(*min_element(res_vals, res_vals+N)));
      ell_itn = 0;
      // ************
      //ell_res = 0; // DEBUGGING
      // ************
      while (ell_res > ell_tol) {
	r = rmin;
	for (j = 1; j < lastpt; ++j) {
	  r += dr;
	  set_inds(col_inds, 3*j);
	  set_alphavals(mat_vals, xi, pi, alpha, beta, psi, j, dr, r);
	  ierr = MatSetValues(jac, 1, &indices[3*j], 9, col_inds, mat_vals, INSERT_VALUES);
	  CHKERRQ(ierr);
	  set_betavals(mat_vals, xi, pi, alpha, beta, psi, j, dr, r);
	  ierr = MatSetValues(jac, 1, &indices[3*j+1], 9, col_inds, mat_vals, INSERT_VALUES);
	  CHKERRQ(ierr);	
	  set_psivals(mat_vals, xi, pi, alpha, beta, psi, j, dr, r);
	  ierr = MatSetValues(jac, 1, &indices[3*j+2], 9, col_inds, mat_vals, INSERT_VALUES);
	  CHKERRQ(ierr);	
	}
      // matrix/vector assembly & checking
	ierr = MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	
	ierr = VecSetValues(abpres, N, &indices[0], &res_vals[0], INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(abpres); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(abpres); CHKERRQ(ierr);

	/*
	if (i == 0) {
	  PetscViewer viewer;
	  ierr = PetscViewerASCIIOpen(comm, "jac-matrix.m", &viewer); CHKERRQ(ierr);
	  ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATHEMATICA); CHKERRQ(ierr);
	  ierr = MatView(jac, viewer); CHKERRQ(ierr);
	  ierr = PetscViewerPopFormat(viewer); CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
	}
	*/
      
	// solve linear system
	ierr = KSPSetOperators(ksp, jac, jac); CHKERRQ(ierr);
	ierr = KSPSolve(ksp, abpres, abpres); CHKERRQ(ierr);
	/*
	ierr = KSPGetIterationNumber(ksp, &petsc_itn); CHKERRQ(ierr);
	if (petsc_itn < 2 || petsc_itn > ell_maxit) {
	  cout << i << " petsc_itn= " << petsc_itn << " at ell_itn " << ell_itn << endl; }
	*/
	// add displacements to metric functions
	ierr = VecGetArray(abpres, &deltas); CHKERRQ(ierr);
	for (j = 0; j < npts; ++j) {
	  alpha[j] -= deltas[3*j];
	  beta[j] -= deltas[3*j+1];
	  psi[j] -= deltas[3*j+2];
	}
	ierr = VecRestoreArray(abpres, &deltas); CHKERRQ(ierr);
	
	// get new residual
	get_ell_res(res_vals, xi, pi, alpha, beta, psi, lastpt, N, dr, rmin);
	ell_res = max(*max_element(res_vals, res_vals+N),
		      abs(*min_element(res_vals, res_vals+N)));
	++ell_itn; 
	if (ell_itn % ell_maxit == 0) {
	  ell_res = 0.0;
	  ++ell_maxit_count;
	  if (i % 500*factor == 0) { cout << i << " ell_res= " << res << " at " << ell_itn << endl; }
	}
      }
// **************************************************************************
// ****************** ELLIPTIC ITERATIVE SOLUTION COMPLETE ******************
// **************************************************************************
    ++hyp_ell_itn;
    }
    
// ***********************************************************************
// ***********************************************************************
// ****************** FULL ITERATIVE SOLUTION COMPLETE *******************
// ***********************************************************************
// ***********************************************************************
    
    // ****************** WRITE MASS & update field **********************
    if (wr_mass && i % check_step*save_step == 0) {
      //mass_check(xi, pi, alpha, beta, psi, dr, rmin, t, ofs_mass); }
      maspect[0] = mass_aspect0(alpha, beta, psi, dr, rmin);
      r = rmin;
      for (j = 1; j < lastwr; ++j) {
	r += save_pt*dr;
	maspect[j] = mass_aspect(alpha, beta, psi, j, dr, r);
      }
      maspect[lastwr] = mass_aspectR(alpha, beta, psi, lastwr, dr, rmax);
      gft_out_bbox(&maspect_name[0], t, bbh_shape, bbh_rank, coords, &maspect[0]);
    }
    if (wr_sol) { update_sol(xi, pi, alpha, beta, psi, sol, dt, save_pt, wr_shape); } 
  }
  // ******************** DONE TIME STEPPING *********************

  // close and destroy objects then finalize mpi/petsc
  ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
  ierr = MatDestroy(&jac); CHKERRQ(ierr);
  ierr = VecDestroy(&abpres); CHKERRQ(ierr);
  
  // write final time step
  if (nsteps % save_step == 0) {
    get_wr_arr(xi, pi, wr_xi, wr_pi, lastwr+1, save_pt);
    wr_step(field_arr, 2, name_arr, nsteps*dt, bbh_shape, bbh_rank, coords);
  }
  // close outfiles
  gft_close_all();
  if (wr_itn) { ofs_itn.close(); }
  // print resolution runtime
  cout << difftime(time(NULL), start_time) << " seconds elapsed" << endl;
  //*************DEBUG************
  cout << maxit_count << " hyperbolic steps reached maxit=" << maxit << endl;
  cout << ell_maxit_count << " elliptic steps reached maxit=" << ell_maxit << "\n" << endl;
  
  }
  // ******************** DONE LOOPING OVER RESOLUTIONS *********************
  ierr = PetscFinalize();
  if (ierr) { return ierr; }
  
  return 0;
}

  //TEST ierr = VecView(abp, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

    /*
  ierr = PetscViewerASCIIOpen(comm, "jac-matrix.m", &viewer); CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATHEMATICA); CHKERRQ(ierr);
  ierr = MatView(jac, viewer); CHKERRQ(ierr);
  ierr = PetscViewerPopFormat(viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  */
