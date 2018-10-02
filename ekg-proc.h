#ifndef EKG_PROC_H_INCLUDED
#define EKG_PROC_H_INCLUDED

#include <iostream>
#include <algorithm> // for max_element()
#include <fstream> // for mass/itn files
#include <ctime> // for quick time analysis
#include "lapacke.h"
#include "ekg-fns.h"
#include "fda-fns.h"
#include "fda-io.h"
#include <vector> // for everything
#include <cmath> // for ICs
#include <string> // for parameter input
#include <map> // for parameter input
#include <cstdlib> // for atoi() and atof()
#include "bbhutil.h" // for output to .sdf


void update_sol(const vector<double>& xi, const vector<double>& pi,
		const vector<double>& alpha, const vector<double>& beta,
		const vector<double>& psi, vector<double>& sol,
		double dt, int savept, int wr_shape);
// compute and write mass
void mass_check(const vector<double>& xi, const vector<double>& pi,
		const vector<double>& alpha, const vector<double>& beta,
		const vector<double>& psi, double dr, double rmin,
		double t, ofstream& out_stream);
// get coarsened arrays from fields for writing
void get_wr_arr(const vector<double>& f1, const vector<double>& f2,
		vector<double>& wr1, vector<double>& wr2,
		int one_past_last, int savept);
void get_wr_arr_ires(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& oldxi, const vector<double>& oldpi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, vector<double>& wrxi, vector<double>& wrpi,
		     vector<double>& iresxi, vector<double>& irespi,
		     int lastwrite, int savept, double lam, double dr, double rmin);
inline void writing(int lastwr, int wr_shape, vector<double>& wr_xi, vector<double>& wr_pi,
		    double *field_arr[2], int *bbh_shape, int bbh_rank, double *coords,
		    bool wr_sol, string solname, vector<double>& sol, bool wr_xp,
		    char *name_arr[2], const vector<double>& xi, const vector<double>& pi,
		    const vector<double>& old_xi, const vector<double>& old_pi,
		    const vector<double>& alpha, const vector<double>& beta,
		    const vector<double>& psi, vector<double>& iresxi, vector<double>& irespi,
		    const vector<double>& resxi, const vector<double>& respi, bool wr_res,
		    char *resname_arr[2], bool wr_ires, double *ires_arr[2], char *iresname_arr[2],
		    bool wr_abp, bool wr_alpha, string alphaname, bool wr_beta, string betaname,
		    bool wr_psi, string psiname, bool wr_abpires, char *abpiresname_arr[3],
		    int save_pt, double lam, double dr, double rmin, double t);
// set rhs of A.x(t_n+1) = b(t_n) for xi, pi (xpp version for including psi)
inline void set_rhs(vector<double>& bxi, vector<double>& bpi,
		    const vector<double>& old_xi, const vector<double>& old_pi,
		    const vector<double>& alpha, const vector<double>& beta,
		    const vector<double>& psi, double lam, double dr, double r,
		    int lastpt);
inline void set_rhs_xpp(vector<double>& bxi, vector<double>& bpi, vector<double> bpsi,
			const vector<double>& old_xi, const vector<double>& old_pi,
			const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, vector<double>& old_psi,
			double lam, double dr, double r, int lastpt);
// perform gauss-seidel update on xi, pi (xpp version for including psi)
void gs_update(const vector<double>& bxi, const vector<double>& bpi,
	       vector<double>& resxi, vector<double>& respi,
	       vector<double>& xi, vector<double>& pi,
	       const vector<double>& alpha, const vector<double>& beta,
	       const vector<double>& psi, double lam, double dr, double r,
	       double rmin, int lastpt, double somm_coeff);
void gs_dr3update(const vector<double>& bxi, const vector<double>& bpi,
		  vector<double>& resxi, vector<double>& respi,
		  vector<double>& xi, vector<double>& pi,
		  const vector<double>& alpha, const vector<double>& beta,
		  const vector<double>& psi, double lam, double dr, double r,
		  double rmin, int lastpt, double somm_coeff);
void gs_up_xpp(const vector<double>& bxi, const vector<double>& bpi,
	       const vector<double>& bpsi, vector<double>& resxi,
	       vector<double>& respi, vector<double>& respsi,
	       vector<double>& xi, vector<double>& pi,
	       const vector<double>& alpha, const vector<double>& beta,
	       vector<double>& psi, double lam, double dr, double r,
	       double rmin, int lastpt, double somm_coeff);
// petsc matrix setup
inline void set_inds(int *cols, int indalpha);
inline void set_alphavals(double *vals, const vector<double>& xi, const vector<double>& pi,
			  const vector<double>& alpha, const vector<double>& beta,
			  const vector<double>& psi, int ind, double dr, double r);
inline void set_betavals(double *vals, const vector<double>& xi, const vector<double>& pi,
			 const vector<double>& alpha, const vector<double>& beta,
			 const vector<double>& psi, int ind, double dr, double r);
inline void set_psivals(double *vals, const vector<double>& xi, const vector<double>& pi,
			const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, double dr, double r);
// set res_vals to be residual = L[f] for rhs of jacobian.(-delta) = residual
void get_ell_res(vector<double>& abpres, const vector<double>& xi, const vector<double>& pi,
		 const vector<double>& alpha, const vector<double>& beta,
		 const vector<double>& psi, int lastpt, int N, double dr, double r);
void get_pres(vector<double>& pres, const vector<double>& xi, const vector<double>& pi,
	      const vector<double>& alpha, const vector<double>& beta,
	      const vector<double>& psi, int lastpt, int N, double dr, double r);
void get_bres(vector<double>& bres, const vector<double>& xi, const vector<double>& pi,
	      const vector<double>& alpha, const vector<double>& beta,
	      const vector<double>& psi, int lastpt, int N, double dr, double r);
void get_ares(vector<double>& ares, const vector<double>& xi, const vector<double>& pi,
	      const vector<double>& alpha, const vector<double>& beta,
	      const vector<double>& psi, int lastpt, int N, double dr, double r);
//  for LAPACKE_dgbsv(): jac[ (kl + ku + 1) + (ldab - 1)*j + i ]  =  jac[ i, j ]
inline void set_jac_vecCM(vector<double>& jac, const vector<double>& xi, const vector<double>& pi,
			  const vector<double>& alpha, const vector<double>& beta, const vector<double>& psi,
			  int N, int ldab, int kl, int ku, int last, double dr, double r);
inline void set_jac_psiCM(vector<double>& jac, const vector<double>& xi, const vector<double>& pi,
			  const vector<double>& alpha, const vector<double>& beta, const vector<double>& psi,
			  int N, int ldab, int kl, int ku, int one_past_last, double dr, double r);
inline void set_jac_betaCM(vector<double>& jac, const vector<double>& xi, const vector<double>& pi,
			   const vector<double>& alpha, const vector<double>& beta, const vector<double>& psi,
			   int N, int ldab, int kl, int ku, int one_past_last, double dr, double r);
inline void set_jac_alphaCM(vector<double>& jac, const vector<double>& xi, const vector<double>& pi,
			    const vector<double>& alpha, const vector<double>& beta, const vector<double>& psi,
			    int N, int ldab, int kl, int ku, int one_past_last, double dr, double r);
int ell_solve_abp_diag(vector<double>& jac, vector<double>& abpres,
		       const vector<double>& xi, const vector<double>& pi,
		       const vector<double>& alpha, const vector<double>& beta,
		       const vector<double>& psi, int lastpt, double dr, double rmin,
		       int N, int kl, int ku, int nrhs, int ldab, vector<int>& ipiv,
		       int ldb, int ell_maxit, double ell_tol, double t, int *ell_maxit_count);
int ell_solve_abp_full(vector<double>& jac, vector<double>& abpres,
		       const vector<double>& xi, const vector<double>& pi,
		       vector<double>& alpha, vector<double>& beta,
		       vector<double>& psi, int lastpt, double dr, double rmin,
		       int N, int kl, int ku, int nrhs, int ldab, vector<int>& ipiv,
		       int ldb, int ell_maxit, double ell_tol, double t, int *ell_maxit_count);
int ell_solve_static_metric(vector<double>& jac, vector<double>& abpres,
			    const vector<double>& xi, const vector<double>& pi,
			    vector<double>& alpha, vector<double>& beta,
			    vector<double>& psi, int lastpt, double dr, double rmin,
			    int N, int kl, int ku, int nrhs, int ldab, vector<int>& ipiv,
			    int ldb, int ell_maxit, double ell_tol, double t,
			    int *ell_maxit_count) { return 0; }
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
void update_sol(const vector<double>& xi, const vector<double>& pi,
		const vector<double>& alpha, const vector<double>& beta,
		const vector<double>& psi, vector<double>& sol,
		double dt, int savept, int wr_shape) {
  int k, s;
  for (k = 0; k < wr_shape; ++k) {
    s = k * savept;
    sol[k] += dt * (alpha[s]*sqin(psi[s])*pi[s] + beta[s]*xi[s]);
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////
// WRITING
////////////////////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////////////////////
// GET COARSENED ARRAY FOR WRITING
void get_wr_arr(const vector<double>& f1, const vector<double>& f2, vector<double>& wr1,
		vector<double>& wr2, int one_past_last, int savept)
{
  int k, s = 0;
  for (k = 0; k < one_past_last; ++k) {
    wr1[k] = f1[s];
    wr2[k] = f2[s];
    s += savept;
  }
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////////////////////
inline void writing(int lastwr, int wr_shape, vector<double>& wr_xi, vector<double>& wr_pi,
		    double *field_arr[2], int *bbh_shape, int bbh_rank, double *coords,
		    bool wr_sol, string solname, vector<double>& sol, bool wr_xp,
		    char *name_arr[2], const vector<double>& xi, const vector<double>& pi,
		    const vector<double>& old_xi, const vector<double>& old_pi,
		    const vector<double>& alpha, const vector<double>& beta,
		    const vector<double>& psi, vector<double>& iresxi, vector<double>& irespi,
		    const vector<double>& resxi, const vector<double>& respi, bool wr_res,
		    char *resname_arr[2], bool wr_ires, double *ires_arr[2], char *iresname_arr[2],
		    bool wr_abp, bool wr_alpha, string alphaname, bool wr_beta, string betaname,
		    bool wr_psi, string psiname, bool wr_abpires, char *abpiresname_arr[3],
		    int save_pt, double lam, double dr, double rmin, double t)
{
  if (wr_sol) {
    gft_out_bbox(&solname[0], t, bbh_shape, bbh_rank, coords, &sol[0]); }
  if (wr_xp) {
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
  }
  if (wr_abp) {
    if (wr_alpha) {
      get_wr_f(alpha, wr_xi, wr_shape, save_pt);
      gft_out_bbox(&alphaname[0], t, bbh_shape, bbh_rank, coords, &wr_xi[0]);
    }
    if (wr_beta) {
      get_wr_f(beta, wr_xi, wr_shape, save_pt);
      gft_out_bbox(&betaname[0], t, bbh_shape, bbh_rank, coords, &wr_xi[0]);
    }
    if (wr_psi) {
      get_wr_f(psi, wr_xi, wr_shape, save_pt);
      gft_out_bbox(&psiname[0], t, bbh_shape, bbh_rank, coords, &wr_xi[0]);
    }
  }
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////
// set rhs of A.x(t_n+1) = b(t_n)
////////////////////////////////////////////////////////////////////////////////////////////////
inline void set_rhs(vector<double>& bxi, vector<double>& bpi,
		    const vector<double>& old_xi, const vector<double>& old_pi,
		    const vector<double>& alpha, const vector<double>& beta,
		    const vector<double>& psi, double lam, double dr, double r,
		    int lastpt) {
  if (r != 0) {
    bxi[0] = old_xi[0] + fda0_xi(old_xi, old_pi, alpha, beta, psi, lam);
    bpi[0] = old_pi[0] + fda0_pi(old_xi, old_pi, alpha, beta, psi, lam, dr, r);
  }
  for (int j = 1; j < lastpt; ++j) {
    r += dr;
    bxi[j] = old_xi[j] + fda_xi(old_xi, old_pi, alpha, beta, psi, j, lam);
    bpi[j] = old_pi[j] + fda_pi(old_xi, old_pi, alpha, beta, psi, j, lam, dr, r);
  }
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////
inline void set_rhs_xpp(vector<double>& bxi, vector<double>& bpi, vector<double>& bpsi,
			const vector<double>& old_xi, const vector<double>& old_pi,
			const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, vector<double>& old_psi,
			double lam, double dr, double r, int lastpt) {
  old_psi = psi;
  if (r != 0) {
    bxi[0] = old_xi[0] + fda0_xi(old_xi, old_pi, alpha, beta, old_psi, lam);
    bpi[0] = old_pi[0] + fda0_pi(old_xi, old_pi, alpha, beta, old_psi, lam, dr, r);
    bpsi[0] = old_psi[0] + fda0_psi(old_xi, old_pi, alpha, beta, old_psi, lam, dr, r);
  }
  for (int j = 1; j < lastpt; ++j) {
    r += dr;
    bxi[j] = old_xi[j] + fda_xi(old_xi, old_pi, alpha, beta, old_psi, j, lam);
    bpi[j] = old_pi[j] + fda_pi(old_xi, old_pi, alpha, beta, old_psi, j, lam, dr, r);
    bpsi[j] = old_psi[j] + fda_psi(old_xi, old_pi, alpha, beta, old_psi, j, lam, dr, r);
  }
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////
// perform gauss-seidel update (xi and pi only)
////////////////////////////////////////////////////////////////////////////////////////////////
void gs_update(const vector<double>& bxi, const vector<double>& bpi,
	       vector<double>& resxi, vector<double>& respi,
	       vector<double>& xi, vector<double>& pi,
	       const vector<double>& alpha, const vector<double>& beta,
	       const vector<double>& psi, double lam, double dr, double r,
	       double rmin, int lastpt, double somm_coeff) {
  if (r == 0) { neumann0(pi); }
  else {
    xi[0] = bxi[0] + fda0_xi(xi, pi, alpha, beta, psi, lam);
    pi[0] = bpi[0] + fda0_pi(xi, pi, alpha, beta, psi, lam, dr, r);
  }
  r += dr;
  xi[1] = bxi[1] + fda_xi(xi, pi, alpha, beta, psi, 1, lam);
  pi[1] = bpi[1] + fda_pi(xi, pi, alpha, beta, psi, 1, lam, dr, r);
  for (int j = 2; j < lastpt; ++j) {
    r += dr;
    xi[j] = bxi[j] + fda_xi(xi, pi, alpha, beta, psi, j, lam);
    pi[j] = bpi[j] + fda_pi(xi, pi, alpha, beta, psi, j, lam, dr, r);
    resxi[j-1] = abs(xi[j-1] - bxi[j-1] -
		     fda_xi(xi, pi, alpha, beta, psi, j-1, lam));
    respi[j-1] = abs(pi[j-1] - bpi[j-1] -
		     fda_pi(xi, pi, alpha, beta, psi, j-1, lam, dr, r-dr));
  }
  if (rmin == 0) { respi[0] = abs(neumann0res(pi)); }
  else {
    resxi[0] = abs(xi[0] - bxi[0] -
		   fda0_xi(xi, pi, alpha, beta, psi, lam));
    respi[0] = abs(pi[0] - bpi[0] -
		   fda0_pi(xi, pi, alpha, beta, psi, lam, dr, rmin));
  }
  // UPDATE BOUNDARY
  sommerfeld(xi, xi, lastpt, lam, somm_coeff);
  sommerfeld(pi, pi, lastpt, lam, somm_coeff);
  
  resxi[lastpt-1] = abs(xi[lastpt-1] - bxi[lastpt-1] -
			fda_xi(xi, pi, alpha, beta, psi, lastpt-1, lam));
  respi[lastpt-1] = abs(pi[lastpt-1] - bpi[lastpt-1] -
			fda_pi(xi, pi, alpha, beta, psi, lastpt-1, lam, dr, r));
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////
// perform gauss-seidel update on xi, pi, and psi
////////////////////////////////////////////////////////////////////////////////////////////////
void gs_up_xpp(const vector<double>& bxi, const vector<double>& bpi,
	       const vector<double>& bpsi, vector<double>& resxi,
	       vector<double>& respi, vector<double>& respsi,
	       vector<double>& xi, vector<double>& pi,
	       const vector<double>& alpha, const vector<double>& beta,
	       vector<double>& psi, double lam, double dr, double r,
	       double rmin, int lastpt, double somm_coeff) {
  //***********************************************
  //*********   UNFINISHED ************************
  //***********************************************
  if (r == 0) { neumann0(pi); }
  else {
    xi[0] = bxi[0] + fda0_xi(xi, pi, alpha, beta, psi, lam);
    pi[0] = bpi[0] + fda0_pi(xi, pi, alpha, beta, psi, lam, dr, r);
  }
  r += dr;
  xi[1] = bxi[1] + fda_xi(xi, pi, alpha, beta, psi, 1, lam);
  pi[1] = bpi[1] + fda_pi(xi, pi, alpha, beta, psi, 1, lam, dr, r);
  for (int j = 2; j < lastpt; ++j) {
    r += dr;
    xi[j] = bxi[j] + fda_xi(xi, pi, alpha, beta, psi, j, lam);
    pi[j] = bpi[j] + fda_pi(xi, pi, alpha, beta, psi, j, lam, dr, r);
    resxi[j-1] = abs(xi[j-1] - bxi[j-1] -
		     fda_xi(xi, pi, alpha, beta, psi, j-1, lam));
    respi[j-1] = abs(pi[j-1] - bpi[j-1] -
		     fda_pi(xi, pi, alpha, beta, psi, j-1, lam, dr, r-dr));
  }
  if (rmin == 0) { respi[0] = abs(neumann0res(pi)); }
  else {
    resxi[0] = abs(xi[0] - bxi[0] -
		   fda0_xi(xi, pi, alpha, beta, psi, lam));
    respi[0] = abs(pi[0] - bpi[0] -
		   fda0_pi(xi, pi, alpha, beta, psi, lam, dr, rmin));
  }
  // UPDATE BOUNDARY
  sommerfeld(xi, xi, lastpt, lam, somm_coeff);
  sommerfeld(pi, pi, lastpt, lam, somm_coeff);
  
  resxi[lastpt-1] = abs(xi[lastpt-1] - bxi[lastpt-1] -
			fda_xi(xi, pi, alpha, beta, psi, lastpt-1, lam));
  respi[lastpt-1] = abs(pi[lastpt-1] - bpi[lastpt-1] -
			fda_pi(xi, pi, alpha, beta, psi, lastpt-1, lam, dr, r));
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////
// perform gauss-seidel update with d/dr3 scheme (xi and pi only)
////////////////////////////////////////////////////////////////////////////////////////////////
void gs_dr3update(const vector<double>& bxi, const vector<double>& bpi,
		  vector<double>& resxi, vector<double>& respi,
		  vector<double>& xi, vector<double>& pi,
		  const vector<double>& alpha, const vector<double>& beta,
		  const vector<double>& psi, double lam, double dr, double r,
		  double rmin, int lastpt, double somm_coeff) {
  if (r == 0) { neumann0(pi); }
  else {
    xi[0] = bxi[0] + fda0_xi(xi, pi, alpha, beta, psi, lam);
    pi[0] = bpi[0] + fda0_pi(xi, pi, alpha, beta, psi, lam, dr, r);
  }
  r += dr;
  xi[1] = bxi[1] + fda_xi(xi, pi, alpha, beta, psi, 1, lam);
  pi[1] = bpi[1] + fda_dr3_pi(xi, pi, alpha, beta, psi, 1, lam, dr, r);
  for (int j = 2; j < lastpt; ++j) {
    r += dr;
    xi[j] = bxi[j] + fda_xi(xi, pi, alpha, beta, psi, j, lam);
    pi[j] = bpi[j] + fda_dr3_pi(xi, pi, alpha, beta, psi, j, lam, dr, r);
    resxi[j-1] = abs(xi[j-1] - bxi[j-1] -
		     fda_xi(xi, pi, alpha, beta, psi, j-1, lam));
    respi[j-1] = abs(pi[j-1] - bpi[j-1] -
		     fda_dr3_pi(xi, pi, alpha, beta, psi, j-1, lam, dr, r-dr));
  }
  if (rmin == 0) { respi[0] = abs(neumann0res(pi)); }
  else {
    resxi[0] = abs(xi[0] - bxi[0] -
		   fda0_xi(xi, pi, alpha, beta, psi, lam));
    respi[0] = abs(pi[0] - bpi[0] -
		   fda0_pi(xi, pi, alpha, beta, psi, lam, dr, rmin));
  }
  // UPDATE BOUNDARY
  sommerfeld(xi, xi, lastpt, lam, somm_coeff);
  sommerfeld(pi, pi, lastpt, lam, somm_coeff);
  
  resxi[lastpt-1] = abs(xi[lastpt-1] - bxi[lastpt-1] -
			fda_xi(xi, pi, alpha, beta, psi, lastpt-1, lam));
  respi[lastpt-1] = abs(pi[lastpt-1] - bpi[lastpt-1] -
			fda_dr3_pi(xi, pi, alpha, beta, psi, lastpt-1, lam, dr, r));
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////
// petsc set mat vals
////////////////////////////////////////////////////////////////////////////////////////////////
inline void set_inds(int *cols, int indalpha) {
  cols[0] = indalpha-3, cols[1] = indalpha-2, cols[2] = indalpha-1,
    cols[3] = indalpha, cols[4] = indalpha+1, cols[5] = indalpha+2,
    cols[6] = indalpha+3, cols[7] = indalpha+4, cols[8] = indalpha+5;
  return; }
////////////////////////////////////////////////////////////////////////////////////////////////
inline void set_alphavals(double *vals, const vector<double>& xi, const vector<double>& pi,
			  const vector<double>& alpha, const vector<double>& beta,
			  const vector<double>& psi, int ind, double dr, double r) {
  vals[0] = jac_aa_pm(alpha, beta, psi, ind, -1, dr, r),
    vals[1] = jac_ab_pm(alpha, beta, psi, ind, -1, dr, r),
    vals[2] = jac_ap_pm(alpha, beta, psi, ind, -1, dr, r),
    vals[3] = jac_aa(xi, pi, alpha, beta, psi, ind, dr, r),
    vals[4] = jac_ab(alpha, beta, psi, ind, dr, r),
    vals[5] = jac_ap(alpha, beta, psi, ind, dr, r),
    vals[6] = jac_aa_pm(alpha, beta, psi, ind, 1, dr, r),
    vals[7] = jac_ab_pm(alpha, beta, psi, ind, 1, dr, r),
    vals[8] = jac_ap_pm(alpha, beta, psi, ind, 1, dr, r);
  return; }
////////////////////////////////////////////////////////////////////////////////////////////////
inline void set_betavals(double *vals, const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  vals[0] = jac_ba_pm(alpha, beta, psi, ind, -1, dr, r),
    vals[1] = jac_bb_pm(alpha, beta, psi, ind, -1, dr, r),
    vals[2] = jac_bp_pm(alpha, beta, psi, ind, -1, dr, r),
    vals[3] = jac_ba(xi, pi, alpha, beta, psi, ind, dr, r),
    vals[4] = jac_bb(alpha, beta, psi, ind, dr, r),
    vals[5] = jac_bp(xi, pi, alpha, beta, psi, ind, dr, r),
    vals[6] = jac_ba_pm(alpha, beta, psi, ind, 1, dr, r),
    vals[7] = jac_bb_pm(alpha, beta, psi, ind, 1, dr, r),
    vals[8] = jac_bp_pm(alpha, beta, psi, ind, 1, dr, r);
  return; }
////////////////////////////////////////////////////////////////////////////////////////////////
inline void set_psivals(double *vals, const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  vals[0] = jac_pa_pm(),
    vals[1] = jac_pb_pm(alpha, beta, psi, ind, -1, dr, r),
    vals[2] = jac_pp_pm(alpha, beta, psi, ind, -1, dr, r),
    vals[3] = jac_pa(alpha, beta, psi, ind, dr, r),
    vals[4] = jac_pb(alpha, beta, psi, ind, dr, r),
    vals[5] = jac_pp(xi, pi, alpha, beta, psi, ind, dr, r),
    vals[6] = jac_pa_pm(),
    vals[7] = jac_pb_pm(alpha, beta, psi, ind, 1, dr, r),
    vals[8] = jac_pp_pm(alpha, beta, psi, ind, 1, dr, r);
  return; }
////////////////////////////////////////////////////////////////////////////////////////////////
// set res_vals to be residual = L[f] for rhs of jacobian.(-delta) = residual
////////////////////////////////////////////////////////////////////////////////////////////////
void get_ell_res(vector<double>& abpres, const vector<double>& xi, const vector<double>& pi,
		 const vector<double>& alpha, const vector<double>& beta,
		 const vector<double>& psi, int lastpt, int N, double dr, double r) {
  abpres[0] = neumann0res(alpha);
  abpres[1] = dirichlet0res(beta);
  abpres[2] = neumann0res(psi);
  for (int k = 1; k < lastpt; ++k) {
    r += dr;
    abpres[3*k] = fda_resalpha(xi, pi, alpha, beta, psi, k, dr, r);
    abpres[3*k+1] = fda_resbeta(xi, pi, alpha, beta, psi, k, dr, r);
    abpres[3*k+2] = fda_respsi(xi, pi, alpha, beta, psi, k, dr, r);
  }
  r += dr;
  abpres[N-3] = fdaR_resalpha(alpha, lastpt, dr, r);
  abpres[N-2] = fdaR_resbeta(beta, lastpt, dr, r);
  abpres[N-1] = fdaR_respsi(psi, lastpt, dr, r);
  return;
}

void get_pres(vector<double>& pres, const vector<double>& xi, const vector<double>& pi,
	      const vector<double>& alpha, const vector<double>& beta,
	      const vector<double>& psi, int lastpt, int N, double dr, double r) {
  pres[0] = -3*psi[0] + 4*psi[1] - psi[2];
  for (int k = 1; k < lastpt; ++k) {
    r += dr;
    pres[k] = fda_respsi(xi, pi, alpha, beta, psi, k, dr, r);
  }
  r += dr;
  pres[N-1] = fdaR_respsi(psi, lastpt, dr, r);
  return;
}

void get_bres(vector<double>& bres, const vector<double>& xi, const vector<double>& pi,
	      const vector<double>& alpha, const vector<double>& beta,
	      const vector<double>& psi, int lastpt, int N, double dr, double r) {
  bres[0] = beta[0];
  for (int k = 1; k < lastpt; ++k) {
    r += dr;
    bres[k] = fda_resbeta(xi, pi, alpha, beta, psi, k, dr, r);
  }
  r += dr;
  bres[N-1] = fdaR_resbeta(beta, lastpt, dr, r);
  return;
}

void get_ares(vector<double>& ares, const vector<double>& xi, const vector<double>& pi,
	      const vector<double>& alpha, const vector<double>& beta,
	      const vector<double>& psi, int lastpt, int N, double dr, double r) {
  ares[0] = -3*alpha[0] + 4*alpha[1] - alpha[2];
  for (int k = 1; k < lastpt; ++k) {
    r += dr;
    ares[k] = fda_resalpha(xi, pi, alpha, beta, psi, k, dr, r);
  }
  r += dr;
  ares[N-1] = fdaR_resalpha(alpha, lastpt, dr, r);
  return;
}


////////////////////////////////////////////////////////////////////////////////////////////////
//  for LAPACKE_dgbsv(): jac[ (kl + ku + 1) + (ldab - 1)*j + i ]  =  jac[ i, j ]
////////////////////////////////////////////////////////////////////////////////////////////////
inline void set_jac_vecCM(vector<double>& jac, const vector<double>& xi, const vector<double>& pi,
			  const vector<double>& alpha, const vector<double>& beta, const vector<double>& psi,
			  int N, int ldab, int kl, int ku, int one_past_last, double dr, double r) {
  int k = kl + ku;
  double rp1 = r + dr;
  // col 0
  jac[k] = -3; jac[++k] = 0; jac[++k] = 0;
  jac[++k] = jac_aa_pm(alpha, beta, psi, 1, -1, dr, rp1);
  jac[++k] = jac_ba_pm(alpha, beta, psi, 1, -1, dr, rp1);
  jac[++k] = jac_pa_pm(); jac[++k] = 0;
  // col 1
  k += kl + 5;
  jac[++k] = 0; jac[++k] = 1; jac[++k] = 0;
  jac[++k] = jac_ab_pm(alpha, beta, psi, 1, -1, dr, rp1);
  jac[++k] = jac_bb_pm(alpha, beta, psi, 1, -1, dr, rp1);
  jac[++k] = jac_pb_pm(alpha, beta, psi, 1, -1, dr, rp1); jac[++k] = 0; jac[++k] = 0;
  // col 2
  k += kl + 4;
  jac[++k] = 0; jac[++k] = 0; jac[++k] = -3;
  jac[++k] = jac_ap_pm(alpha, beta, psi, 1, -1, dr, rp1);
  jac[++k] = jac_bp_pm(alpha, beta, psi, 1, -1, dr, rp1);
  jac[++k] = jac_pp_pm(alpha, beta, psi, 1, -1, dr, rp1); jac[++k] = 0; jac[++k] = 0; jac[++k] = 0;
  // col 3
  r = rp1; rp1 += dr;
  k += kl + 3;
  jac[++k] = 4; jac[++k] = 0; jac[++k] = 0;
  jac[++k] = jac_aa(xi, pi, alpha, beta, psi, 1, dr, r);
  jac[++k] = jac_ba(xi, pi, alpha, beta, psi, 1, dr, r);
  jac[++k] = jac_pa(alpha, beta, psi, 1, dr, r);
  jac[++k] = jac_aa_pm(alpha, beta, psi, 2, -1, dr, rp1);
  jac[++k] = jac_ba_pm(alpha, beta, psi, 2, -1, dr, rp1);
  jac[++k] = jac_pa_pm(); jac[++k] = 0;
  // col 4
  k += kl + 2;
  jac[++k] = 0; jac[++k] = 0; jac[++k] = 0;
  jac[++k] = jac_ab(alpha, beta, psi, 1, dr, r);
  jac[++k] = jac_bb(alpha, beta, psi, 1, dr, r);
  jac[++k] = jac_pb(alpha, beta, psi, 1, dr, r);
  jac[++k] = jac_ab_pm(alpha, beta, psi, 2, -1, dr, rp1);
  jac[++k] = jac_bb_pm(alpha, beta, psi, 2, -1, dr, rp1);
  jac[++k] = jac_pb_pm(alpha, beta, psi, 2, -1, dr, rp1); jac[++k] = 0; jac[++k] = 0;
  // col 5
  k += kl + 1;
  jac[++k] = 0; jac[++k] = 0; jac[++k] = 4;
  jac[++k] = jac_ap(alpha, beta, psi, 1, dr, r);
  jac[++k] = jac_bp(xi, pi, alpha, beta, psi, 1, dr, r);
  jac[++k] = jac_pp(xi, pi, alpha, beta, psi, 1, dr, r);
  jac[++k] = jac_ap_pm(alpha, beta, psi, 2, -1, dr, rp1);
  jac[++k] = jac_bp_pm(alpha, beta, psi, 2, -1, dr, rp1);
  jac[++k] = jac_pp_pm(alpha, beta, psi, 2, -1, dr, rp1); jac[++k] = 0; jac[++k] = 0; jac[++k] = 0;
  // col 6
  double rm1 = r; r = rp1; rp1 += dr;
  k += kl;
  jac[++k] = -1; jac[++k] = 0; jac[++k] = 0;
  jac[++k] = jac_aa_pm(alpha, beta, psi, 1, 1, dr, rm1);
  jac[++k] = jac_ba_pm(alpha, beta, psi, 1, 1, dr, rm1);
  jac[++k] = jac_pa_pm();
  jac[++k] = jac_aa(xi, pi, alpha, beta, psi, 2, dr, r);
  jac[++k] = jac_ba(xi, pi, alpha, beta, psi, 2, dr, r);
  jac[++k] = jac_pa(alpha, beta, psi, 2, dr, r);
  jac[++k] = jac_aa_pm(alpha, beta, psi, 3, -1, dr, rp1);
  jac[++k] = jac_ba_pm(alpha, beta, psi, 3, -1, dr, rp1);
  jac[++k] = jac_pa_pm(); jac[++k] = 0;
  // col 7
  k += kl;
  jac[++k] = 0; jac[++k] = 0;
  jac[++k] = jac_ab_pm(alpha, beta, psi, 1, 1, dr, rm1);
  jac[++k] = jac_bb_pm(alpha, beta, psi, 1, 1, dr, rm1);
  jac[++k] = jac_pb_pm(alpha, beta, psi, 1, 1, dr, rm1);
  jac[++k] = jac_ab(alpha, beta, psi, 2, dr, r);
  jac[++k] = jac_bb(alpha, beta, psi, 2, dr, r);
  jac[++k] = jac_pb(alpha, beta, psi, 2, dr, r);
  jac[++k] = jac_ab_pm(alpha, beta, psi, 3, -1, dr, rp1);
  jac[++k] = jac_bb_pm(alpha, beta, psi, 3, -1, dr, rp1);
  jac[++k] = jac_pb_pm(alpha, beta, psi, 3, -1, dr, rp1); jac[++k] = 0; jac[++k] = 0;
  // col 8
  k += kl;
  jac[++k] = -1;
  jac[++k] = jac_ap_pm(alpha, beta, psi, 1, 1, dr, rm1);
  jac[++k] = jac_bp_pm(alpha, beta, psi, 1, 1, dr, rm1);
  jac[++k] = jac_pp_pm(alpha, beta, psi, 1, 1, dr, rm1);
  jac[++k] = jac_ap(alpha, beta, psi, 2, dr, r);
  jac[++k] = jac_bp(xi, pi, alpha, beta, psi, 2, dr, r);
  jac[++k] = jac_pp(xi, pi, alpha, beta, psi, 2, dr, r);
  jac[++k] = jac_ap_pm(alpha, beta, psi, 3, -1, dr, rp1);
  jac[++k] = jac_bp_pm(alpha, beta, psi, 3, -1, dr, rp1);
  jac[++k] = jac_pp_pm(alpha, beta, psi, 3, -1, dr, rp1); jac[++k] = 0; jac[++k] = 0; jac[++k] = 0;
  int j, jm1 = 2, jp1 = 4;
  for (j = 3; j < one_past_last; ++j) {
    // col d/d(alpha)
    rm1 = r; r = rp1; rp1 += dr;
    k += kl; jac[++k] = 0; jac[++k] = 0; jac[++k] = 0;
    jac[++k] = jac_aa_pm(alpha, beta, psi, jm1, 1, dr, rm1);
    jac[++k] = jac_ba_pm(alpha, beta, psi, jm1, 1, dr, rm1);
    jac[++k] = jac_pa_pm();
    jac[++k] = jac_aa(xi, pi, alpha, beta, psi, j, dr, r);
    jac[++k] = jac_ba(xi, pi, alpha, beta, psi, j, dr, r);
    jac[++k] = jac_pa(alpha, beta, psi, j, dr, r);
    jac[++k] = jac_aa_pm(alpha, beta, psi, jp1, -1, dr, rp1);
    jac[++k] = jac_ba_pm(alpha, beta, psi, jp1, -1, dr, rp1);
    jac[++k] = jac_pa_pm(); jac[++k] = 0;
    // col d/d(beta)
    k += kl; jac[++k] = 0; jac[++k] = 0;
    jac[++k] = jac_ab_pm(alpha, beta, psi, jm1, 1, dr, rm1);
    jac[++k] = jac_bb_pm(alpha, beta, psi, jm1, 1, dr, rm1);
    jac[++k] = jac_pb_pm(alpha, beta, psi, jm1, 1, dr, rm1);
    jac[++k] = jac_ab(alpha, beta, psi, j, dr, r);
    jac[++k] = jac_bb(alpha, beta, psi, j, dr, r);
    jac[++k] = jac_pb(alpha, beta, psi, j, dr, r);
    jac[++k] = jac_ab_pm(alpha, beta, psi, jp1, -1, dr, rp1);
    jac[++k] = jac_bb_pm(alpha, beta, psi, jp1, -1, dr, rp1);
    jac[++k] = jac_pb_pm(alpha, beta, psi, jp1, -1, dr, rp1); jac[++k] = 0; jac[++k] = 0;
    // col d/d(psi)
    k += kl; jac[++k] = 0;
    jac[++k] = jac_ap_pm(alpha, beta, psi, jm1, 1, dr, rm1);
    jac[++k] = jac_bp_pm(alpha, beta, psi, jm1, 1, dr, rm1);
    jac[++k] = jac_pp_pm(alpha, beta, psi, jm1, 1, dr, rm1);
    jac[++k] = jac_ap(alpha, beta, psi, j, dr, r);
    jac[++k] = jac_bp(xi, pi, alpha, beta, psi, j, dr, r);
    jac[++k] = jac_pp(xi, pi, alpha, beta, psi, j, dr, r);
    jac[++k] = jac_ap_pm(alpha, beta, psi, jp1, -1, dr, rp1);
    jac[++k] = jac_bp_pm(alpha, beta, psi, jp1, -1, dr, rp1);
    jac[++k] = jac_pp_pm(alpha, beta, psi, jp1, -1, dr, rp1); jac[++k] = 0; jac[++k] = 0; jac[++k] = 0;
    jm1 = j; ++jp1;
  }
  j = one_past_last;
  // col N-9
  rm1 = r; r = rp1; rp1 += dr;
  k += kl; jac[++k] = 0; jac[++k] = 0; jac[++k] = 0;
  jac[++k] = jac_aa_pm(alpha, beta, psi, jm1, 1, dr, rm1);
  jac[++k] = jac_ba_pm(alpha, beta, psi, jm1, 1, dr, rm1);
  jac[++k] = jac_pa_pm();
  jac[++k] = jac_aa(xi, pi, alpha, beta, psi, j, dr, r);
  jac[++k] = jac_ba(xi, pi, alpha, beta, psi, j, dr, r);
  jac[++k] = jac_pa(alpha, beta, psi, j, dr, r);
  jac[++k] = jac_aa_pm(alpha, beta, psi, jp1, -1, dr, rp1);
  jac[++k] = jac_ba_pm(alpha, beta, psi, jp1, -1, dr, rp1);
  jac[++k] = jac_pa_pm(); jac[++k] = 1;
  // col N-8
  k += kl; jac[++k] = 0; jac[++k] = 0;
  jac[++k] = jac_ab_pm(alpha, beta, psi, jm1, 1, dr, rm1);
  jac[++k] = jac_bb_pm(alpha, beta, psi, jm1, 1, dr, rm1);
  jac[++k] = jac_pb_pm(alpha, beta, psi, jm1, 1, dr, rm1);
  jac[++k] = jac_ab(alpha, beta, psi, j, dr, r);
  jac[++k] = jac_bb(alpha, beta, psi, j, dr, r);
  jac[++k] = jac_pb(alpha, beta, psi, j, dr, r);
  jac[++k] = jac_ab_pm(alpha, beta, psi, jp1, -1, dr, rp1);
  jac[++k] = jac_bb_pm(alpha, beta, psi, jp1, -1, dr, rp1);
  jac[++k] = jac_pb_pm(alpha, beta, psi, jp1, -1, dr, rp1); jac[++k] = 0; jac[++k] = 1;
  // col N-7
  k += kl; jac[++k] = 0;
  jac[++k] = jac_ap_pm(alpha, beta, psi, jm1, 1, dr, rm1);
  jac[++k] = jac_bp_pm(alpha, beta, psi, jm1, 1, dr, rm1);
  jac[++k] = jac_pp_pm(alpha, beta, psi, jm1, 1, dr, rm1);
  jac[++k] = jac_ap(alpha, beta, psi, j, dr, r);
  jac[++k] = jac_bp(xi, pi, alpha, beta, psi, j, dr, r);
  jac[++k] = jac_pp(xi, pi, alpha, beta, psi, j, dr, r);
  jac[++k] = jac_ap_pm(alpha, beta, psi, jp1, -1, dr, rp1);
  jac[++k] = jac_bp_pm(alpha, beta, psi, jp1, -1, dr, rp1);
  jac[++k] = jac_pp_pm(alpha, beta, psi, jp1, -1, dr, rp1); jac[++k] = 0; jac[++k] = 0; jac[++k] = 1;
  jm1 = j; ++j;
  // col N-6
  rm1 = r; r = rp1;
  k += kl; jac[++k] = 0; jac[++k] = 0; jac[++k] = 0;
  jac[++k] = jac_aa_pm(alpha, beta, psi, jm1, 1, dr, rm1);
  jac[++k] = jac_ba_pm(alpha, beta, psi, jm1, 1, dr, rm1);
  jac[++k] = jac_pa_pm();
  jac[++k] = jac_aa(xi, pi, alpha, beta, psi, j, dr, r);
  jac[++k] = jac_ba(xi, pi, alpha, beta, psi, j, dr, r);
  jac[++k] = jac_pa(alpha, beta, psi, j, dr, r);
  jac[++k] = -4; jac[++k] = 0; jac[++k] = 0;
  // col N-5
  k += 1 + kl; jac[++k] = 0; jac[++k] = 0;
  jac[++k] = jac_ab_pm(alpha, beta, psi, jm1, 1, dr, rm1);
  jac[++k] = jac_bb_pm(alpha, beta, psi, jm1, 1, dr, rm1);
  jac[++k] = jac_pb_pm(alpha, beta, psi, jm1, 1, dr, rm1);
  jac[++k] = jac_ab(alpha, beta, psi, j, dr, r);
  jac[++k] = jac_bb(alpha, beta, psi, j, dr, r);
  jac[++k] = jac_pb(alpha, beta, psi, j, dr, r);
  jac[++k] = 0; jac[++k] = -4; jac[++k] = 0;
  // col N-4
  k += 2 + kl; jac[++k] = 0;
  jac[++k] = jac_ap_pm(alpha, beta, psi, jm1, 1, dr, rm1);
  jac[++k] = jac_bp_pm(alpha, beta, psi, jm1, 1, dr, rm1);
  jac[++k] = jac_pp_pm(alpha, beta, psi, jm1, 1, dr, rm1);
  jac[++k] = jac_ap(alpha, beta, psi, j, dr, r);
  jac[++k] = jac_bp(xi, pi, alpha, beta, psi, j, dr, r);
  jac[++k] = jac_pp(xi, pi, alpha, beta, psi, j, dr, r);
  jac[++k] = 0; jac[++k] = 0; jac[++k] = -4;
  ++jm1;
  // col N-3
  double cR = 3 + (2*dr / (r + dr));
  k += 3 + kl; jac[++k] = 0; jac[++k] = 0; jac[++k] = 0;
  jac[++k] = jac_aa_pm(alpha, beta, psi, jm1, 1, dr, r);
  jac[++k] = jac_ba_pm(alpha, beta, psi, jm1, 1, dr, r);
  jac[++k] = jac_pa_pm();
  jac[++k] = cR; jac[++k] = 0; jac[++k] = 0;
  // col N-2
  k += 4 + kl; jac[++k] = 0; jac[++k] = 0;
  jac[++k] = jac_ab_pm(alpha, beta, psi, jm1, 1, dr, r);
  jac[++k] = jac_bb_pm(alpha, beta, psi, jm1, 1, dr, r);
  jac[++k] = jac_pb_pm(alpha, beta, psi, jm1, 1, dr, r);
  jac[++k] = 0; jac[++k] = cR; jac[++k] = 0;
  // col N-1
  k += 5 + kl; jac[++k] = 0;
  jac[++k] = jac_ap_pm(alpha, beta, psi, jm1, 1, dr, r);
  jac[++k] = jac_bp_pm(alpha, beta, psi, jm1, 1, dr, r);
  jac[++k] = jac_pp_pm(alpha, beta, psi, jm1, 1, dr, r);
  jac[++k] = 0; jac[++k] = 0; jac[++k] = cR;
  return;
}

inline void set_jac_alphaCM(vector<double>& jac, const vector<double>& xi, const vector<double>& pi,
			    const vector<double>& alpha, const vector<double>& beta, const vector<double>& psi,
			    int N, int ldab, int kl, int ku, int one_past_last, double dr, double r) {
  int k = kl + ku;
  double rp1 = r + dr;
  // col 0
  jac[k] = -3;
  jac[++k] = jac_aa_pm(alpha, beta, psi, 1, -1, dr, rp1);
  jac[++k] = 0;
  // col 1
  r = rp1; rp1 += dr;
  k += kl + 1;
  jac[++k] = 4;
  jac[++k] = jac_aa(xi, pi, alpha, beta, psi, 1, dr, r);
  jac[++k] = jac_aa_pm(alpha, beta, psi, 2, -1, dr, rp1);
  jac[++k] = 0;
  // col 2
  double rm1 = r; r = rp1; rp1 += dr;
  k += kl;
  jac[++k] = -1;
  jac[++k] = jac_aa_pm(alpha, beta, psi, 1, 1, dr, rm1);
  jac[++k] = jac_aa(xi, pi, alpha, beta, psi, 2, dr, r);
  jac[++k] = jac_aa_pm(alpha, beta, psi, 3, -1, dr, rp1);
  jac[++k] = 0;
  int j, jm1 = 2, jp1 = 4;
  for (j = 3; j < one_past_last; ++j) {
    // col d/d(alpha)
    rm1 = r; r = rp1; rp1 += dr;
    k += kl;
    jac[++k] = 0;
    jac[++k] = jac_aa_pm(alpha, beta, psi, jm1, 1, dr, rm1);
    jac[++k] = jac_aa(xi, pi, alpha, beta, psi, j, dr, r);
    jac[++k] = jac_aa_pm(alpha, beta, psi, jp1, -1, dr, rp1);
    jac[++k] = 0;
    jm1 = j; ++jp1;
  }
  j = one_past_last;
  // col N-3
  rm1 = r; r = rp1; rp1 += dr;
  k += kl;
  jac[++k] = 0;
  jac[++k] = jac_aa_pm(alpha, beta, psi, jm1, 1, dr, rm1);
  jac[++k] = jac_aa(xi, pi, alpha, beta, psi, j, dr, r);
  jac[++k] = jac_aa_pm(alpha, beta, psi, jp1, -1, dr, rp1);
  jac[++k] = 1;
  jm1 = j; ++j;
  // col N-2
  rm1 = r; r = rp1;
  k += kl;
  jac[++k] = 0;
  jac[++k] = jac_aa_pm(alpha, beta, psi, jm1, 1, dr, rm1);
  jac[++k] = jac_aa(xi, pi, alpha, beta, psi, j, dr, r);
  jac[++k] = -4;
  ++jm1;
  // col N-1
  double cR = 3 + (2*dr / (r + dr));
  k += 1 + kl;
  jac[++k] = 0;
  jac[++k] = jac_aa_pm(alpha, beta, psi, jm1, 1, dr, r);
  jac[++k] = cR;
  return;
}

inline void set_jac_betaCM(vector<double>& jac, const vector<double>& xi, const vector<double>& pi,
			   const vector<double>& alpha, const vector<double>& beta, const vector<double>& psi,
			   int N, int ldab, int kl, int ku, int one_past_last, double dr, double r) {
  int k = kl + ku;
  double rp1 = r + dr;
  // col 0
  jac[k] = 1;
  jac[++k] = jac_bb_pm(alpha, beta, psi, 1, -1, dr, rp1);
  jac[++k] = 0;
  // col 1
  r = rp1; rp1 += dr;
  k += kl + 1;
  jac[++k] = 0;
  jac[++k] = jac_bb(alpha, beta, psi, 1, dr, r);
  jac[++k] = jac_bb_pm(alpha, beta, psi, 2, -1, dr, rp1);
  jac[++k] = 0;
  // col 2
  double rm1 = r; r = rp1; rp1 += dr;
  k += kl;
  jac[++k] = 0;
  jac[++k] = jac_bb_pm(alpha, beta, psi, 1, 1, dr, rm1);
  jac[++k] = jac_bb(alpha, beta, psi, 2, dr, r);
  jac[++k] = jac_bb_pm(alpha, beta, psi, 3, -1, dr, rp1);
  jac[++k] = 0;
  int j, jm1 = 2, jp1 = 4;
  for (j = 3; j < one_past_last; ++j) {
    // col d/d(beta)
    rm1 = r; r = rp1; rp1 += dr;
    k += kl;
    jac[++k] = 0;
    jac[++k] = jac_bb_pm(alpha, beta, psi, jm1, 1, dr, rm1);
    jac[++k] = jac_bb(alpha, beta, psi, j, dr, r);
    jac[++k] = jac_bb_pm(alpha, beta, psi, jp1, -1, dr, rp1);
    jac[++k] = 0;
    jm1 = j; ++jp1;
  }
  j = one_past_last;
  // col N-3
  rm1 = r; r = rp1; rp1 += dr;
  k += kl;
  jac[++k] = 0;
  jac[++k] = jac_bb_pm(alpha, beta, psi, jm1, 1, dr, rm1);
  jac[++k] = jac_bb(alpha, beta, psi, j, dr, r);
  jac[++k] = jac_bb_pm(alpha, beta, psi, jp1, -1, dr, rp1);
  jac[++k] = 1;
  jm1 = j; ++j;
  // col N-2
  rm1 = r; r = rp1;
  k += kl;
  jac[++k] = 0;
  jac[++k] = jac_bb_pm(alpha, beta, psi, jm1, 1, dr, rm1);
  jac[++k] = jac_bb(alpha, beta, psi, j, dr, r);
  jac[++k] = -4;
  ++jm1;
  // col N-1
  double cR = 3 + (2*dr / (r + dr));
  k += 1 + kl;
  jac[++k] = 0;
  jac[++k] = jac_bb_pm(alpha, beta, psi, jm1, 1, dr, r);
  jac[++k] = cR;
  return;
}


inline void set_jac_psiCM(vector<double>& jac, const vector<double>& xi, const vector<double>& pi,
			  const vector<double>& alpha, const vector<double>& beta, const vector<double>& psi,
			  int N, int ldab, int kl, int ku, int one_past_last, double dr, double r) {
  int k = kl + ku;
  double rp1 = r + dr;
  // col 0
  jac[k] = -3;
  jac[++k] = jac_pp_pm(alpha, beta, psi, 1, -1, dr, rp1);
  jac[++k] = 0;
  // col 1
  r = rp1; rp1 += dr;
  k += kl + 1;
  jac[++k] = 4;
  jac[++k] = jac_pp(xi, pi, alpha, beta, psi, 1, dr, r);
  jac[++k] = jac_pp_pm(alpha, beta, psi, 2, -1, dr, rp1);
  jac[++k] = 0;
  // col 2
  double rm1 = r; r = rp1; rp1 += dr;
  k += kl;
  jac[++k] = -1;
  jac[++k] = jac_pp_pm(alpha, beta, psi, 1, 1, dr, rm1);
  jac[++k] = jac_pp(xi, pi, alpha, beta, psi, 2, dr, r);
  jac[++k] = jac_pp_pm(alpha, beta, psi, 3, -1, dr, rp1);
  jac[++k] = 0;
  int j, jm1 = 2, jp1 = 4;
  for (j = 3; j < one_past_last; ++j) {
    // col d/d(psi)
    rm1 = r; r = rp1; rp1 += dr;
    k += kl;
    jac[++k] = 0;
    jac[++k] = jac_pp_pm(alpha, beta, psi, jm1, 1, dr, rm1);
    jac[++k] = jac_pp(xi, pi, alpha, beta, psi, j, dr, r);
    jac[++k] = jac_pp_pm(alpha, beta, psi, jp1, -1, dr, rp1);
    jac[++k] = 0;
    jm1 = j; ++jp1;
  }
  j = one_past_last;
  // col N-3
  rm1 = r; r = rp1; rp1 += dr;
  k += kl;
  jac[++k] = 0;
  jac[++k] = jac_pp_pm(alpha, beta, psi, jm1, 1, dr, rm1);
  jac[++k] = jac_pp(xi, pi, alpha, beta, psi, j, dr, r);
  jac[++k] = jac_pp_pm(alpha, beta, psi, jp1, -1, dr, rp1);
  jac[++k] = 1;
  jm1 = j; ++j;
  // col N-2
  rm1 = r; r = rp1;
  k += kl;
  jac[++k] = 0;
  jac[++k] = jac_pp_pm(alpha, beta, psi, jm1, 1, dr, rm1);
  jac[++k] = jac_pp(xi, pi, alpha, beta, psi, j, dr, r);
  jac[++k] = -4;
  ++jm1;
  // col N-1
  double cR = 3 + (2*dr / (r + dr));
  k += 1 + kl;
  jac[++k] = 0;
  jac[++k] = jac_pp_pm(alpha, beta, psi, jm1, 1, dr, r);
  jac[++k] = cR;
  return;
}



int ell_solve_abp_diag(vector<double>& jac, vector<double>& abpres,
		       const vector<double>& xi, const vector<double>& pi,
		       vector<double>& alpha, vector<double>& beta,
		       vector<double>& psi, int lastpt, double dr, double rmin,
		       int N, int kl, int ku, int nrhs, int ldab, vector<int>& ipiv,
		       int ldb, int ell_maxit, double ell_tol, double t, int *ell_maxit_count) {
  int p_itn = 0, info = 0, j;
  // get initial residual
  get_pres(abpres, xi, pi, alpha, beta, psi, lastpt, N, dr, rmin);
  double ell_res = max(*max_element(abpres.begin(), abpres.end()),
		       abs(*min_element(abpres.begin(), abpres.end())));
  // if ell_res > ell_tol, solve jac.x = abpres and update abpres -= x
  while (ell_res > ell_tol && p_itn < ell_maxit) {
    set_jac_psiCM(jac, xi, pi, alpha, beta, psi, N, ldab, kl, ku, lastpt-2, dr, rmin);
    info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N, kl, ku, nrhs, &jac[0], ldab, &ipiv[0], &abpres[0], ldb);
    if (info != 0) { cout << "t = " << t << " p_info= " << info << " at p_itn = " << p_itn << endl; }
    for (j = 0; j < N; ++j) {
      psi[j] -= abpres[j];
    }	
    // get new residual
    get_pres(abpres, xi, pi, alpha, beta, psi, lastpt, N, dr, rmin);
    ell_res = max(*max_element(abpres.begin(), abpres.end()),
		  abs(*min_element(abpres.begin(), abpres.end())));
    ++p_itn; 
  }

  int b_itn = 0;
  get_bres(abpres, xi, pi, alpha, beta, psi, lastpt, N, dr, rmin);
  ell_res = max(*max_element(abpres.begin(), abpres.end()),
		abs(*min_element(abpres.begin(), abpres.end())));
  while (ell_res > ell_tol && b_itn < ell_maxit) {
    set_jac_betaCM(jac, xi, pi, alpha, beta, psi, N, ldab, kl, ku, lastpt-2, dr, rmin);
    info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N, kl, ku, nrhs, &jac[0], ldab, &ipiv[0], &abpres[0], ldb);
    if (info != 0) { cout << "t = " << t << " b_info= " << info << " at b_itn = " << b_itn << endl; }
    for (j = 0; j < N; ++j) {
      beta[j] -= abpres[j];
    }
    get_bres(abpres, xi, pi, alpha, beta, psi, lastpt, N, dr, rmin);
    ell_res = max(*max_element(abpres.begin(), abpres.end()),
		  abs(*min_element(abpres.begin(), abpres.end())));
    ++b_itn; 
  }

  int a_itn = 0;
  get_ares(abpres, xi, pi, alpha, beta, psi, lastpt, N, dr, rmin);
  ell_res = max(*max_element(abpres.begin(), abpres.end()),
		abs(*min_element(abpres.begin(), abpres.end())));
  while (ell_res > ell_tol && a_itn < ell_maxit) {
    set_jac_alphaCM(jac, xi, pi, alpha, beta, psi, N, ldab, kl, ku, lastpt-2, dr, rmin);
    info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N, kl, ku, nrhs, &jac[0], ldab, &ipiv[0], &abpres[0], ldb);
    if (info != 0) { cout << "t = " << t << " a_info= " << info << " at a_itn = " << a_itn << endl; }
    for (j = 0; j < N; ++j) {
      alpha[j] -= abpres[j];
    }
    get_bres(abpres, xi, pi, alpha, beta, psi, lastpt, N, dr, rmin);
    ell_res = max(*max_element(abpres.begin(), abpres.end()),
		  abs(*min_element(abpres.begin(), abpres.end())));
    ++a_itn; 
  }
  int ell_itn = max(p_itn, max(b_itn, a_itn));
  if (ell_itn == ell_maxit) {
      cout << "t = " << t << " ell_res= " << ell_res << " at " << ell_itn << endl;
      ell_res = 0.0;
      ell_itn = 0;
      ++(*ell_maxit_count);
    }
  return ;
}


int ell_solve_abp_full(vector<double>& jac, vector<double>& abpres,
		       const vector<double>& xi, const vector<double>& pi,
		       vector<double>& alpha, vector<double>& beta,
		       vector<double>& psi, int lastpt, double dr, double rmin,
		       int N, int kl, int ku, int nrhs, int ldab, vector<int>& ipiv,
		       int ldb, int ell_maxit, double ell_tol, double t, int *ell_maxit_count) {
  int ell_itn = 0, info = 0, one_past_last = lastpt + 1, j;
  // get initial residual
  get_ell_res(abpres, xi, pi, alpha, beta, psi, lastpt, N, dr, rmin);
  double ell_res = max(*max_element(abpres.begin(), abpres.end()),
		       abs(*min_element(abpres.begin(), abpres.end())));
  // if ell_res > ell_tol, solve jac.x = abpres and update abpres -= x
  while (ell_res > ell_tol) {
    set_jac_vecCM(jac, xi, pi, alpha, beta, psi, N, ldab, kl, ku, lastpt-2, dr, rmin);
    info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N, kl, ku, nrhs, &jac[0], ldab, &ipiv[0], &abpres[0], ldb);
    if (info != 0) { cout << "t = " << t << " info= " << info << " at ell_itn = " << ell_itn << endl; }
    for (j = 0; j < one_past_last; ++j) {
      alpha[j] -= abpres[3*j];
      beta[j] -= abpres[3*j+1];
      psi[j] -= abpres[3*j+2];
    }	
    // get new residual
    get_ell_res(abpres, xi, pi, alpha, beta, psi, lastpt, N, dr, rmin);
    ell_res = max(*max_element(abpres.begin(), abpres.end()),
		  abs(*min_element(abpres.begin(), abpres.end())));
    ++ell_itn;
    if (ell_itn == ell_maxit) {
      cout << "t = " << t << " ell_res= " << ell_res << " at " << ell_itn << endl;
      ell_res = 0.0;
      ell_itn = 0;
      ++(*ell_maxit_count);
    }
  }  
  return ell_itn;
}


#endif

