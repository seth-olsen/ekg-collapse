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
  return r*pw6(psi[k])*sq(r*ddr_c(beta,k,dr) - beta[k]) / (18*sq(alpha[k])) -
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

// set rhs of A.x(t_n+1) = b(t_n)
inline void set_rhs(vector<double>& bxi, vector<double>& bpi,
		    const vector<double>& old_xi, const vector<double>& old_pi,
		    const vector<double>& alpha, const vector<double>& beta,
		    const vector<double>& psi, double lam, double dr, double r,
		    int j, int lastpt) {
  if (r != 0) {
    bxi[0] = old_xi[0] + fda0_xi(old_xi, old_pi, alpha, beta, psi, lam);
    bpi[0] = old_pi[0] + fda0_pi(old_xi, old_pi, alpha, beta, psi, lam, dr, r);
  }
  for (j = 1; j < lastpt; ++j) {
    r += dr;
    bxi[j] = old_xi[j] + fda_xi(old_xi, old_pi, alpha, beta, psi, j, lam);
    bpi[j] = old_pi[j] + fda_pi(old_xi, old_pi, alpha, beta, psi, j, lam, dr, r);
  }
  return;
}

// perform gauss-seidel update
inline void gs_update(const vector<double>& bxi, const vector<double>& bpi,
		      vector<double>& resxi, vector<double>& respi,
		      vector<double>& xi, vector<double>& pi,
		      const vector<double>& alpha, const vector<double>& beta,
		      const vector<double>& psi, double lam, double dr, double r,
		      double rmin, int j, int lastpt, bool somm_cond, double somm_coeff) {
  if (r == 0) { neumann0(pi); }
  else {
    xi[0] = bxi[0] + fda0_xi(xi, pi, alpha, beta, psi, lam);
    pi[0] = bpi[0] + fda0_pi(xi, pi, alpha, beta, psi, lam, dr, r);
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
  return;
}
  

inline double fda_psi(const vector<double>& xi, const vector<double>& pi,
		      const vector<double>& alpha, const vector<double>& beta,
		      const vector<double>& psi, int ind, double lam, double dr, double r) {
  return 0.25*lam*(beta[ind]*(psi[ind+1] - psi[ind-1]) +
		   psi[ind]*(beta[ind+1] + 4*beta[ind]/r - beta[ind-1])); }

inline double fda_respsi(const vector<double>& xi, const vector<double>& pi,
			 const vector<double>& alpha, const vector<double>& beta,
			 const vector<double>& psi, int ind, double dr, double r) {
  return ddr2_c(psi,ind,dr) + 2*ddr_c(psi,ind,dr)/r +
    (0.5*d_c(beta,ind) - dr*beta[ind]/r)*pw5(psi[ind]) / (12*sq(alpha[ind])) + 
    M_PI*(sq(xi[ind]) + sq(pi[ind]))*psi[ind] ; }

inline double fda_resbeta(const vector<double>& xi, const vector<double>& pi,
			  const vector<double>& alpha, const vector<double>& beta,
			  const vector<double>& psi, int ind, double dr, double r) {
  return  ddr2_c(beta,ind,dr) + 12*M_PI*sqin(psi[ind])*alpha[ind]*xi[ind]*pi[ind]
    + (2/r + 6*ddr_c(psi,ind,dr)/psi[ind] - ddr_c(alpha,ind,dr)/alpha[ind])*
    (0.5*d_c(beta,ind) - dr*beta[ind]/r); }

inline double fda_resalpha(const vector<double>& xi, const vector<double>& pi,
			   const vector<double>& alpha, const vector<double>& beta,
			   const vector<double>& psi, int ind, double dr, double r) {
  return ddr2_c(alpha,ind,dr) + 2*(1/r + ddr_c(psi,ind,dr)/psi[ind])*ddr_c(alpha,ind,dr)
    - 2*pw4(psi[ind])*sq(0.5*d_c(beta,ind) - dr*beta[ind]/r)/(3*alpha[ind])
    - 8*M_PI*alpha[ind]*sq(pi[ind]) ; }

inline double fdaR_respsi(const vector<double>& psi, int ind, double dr, double r) {
  return d_b(psi,ind) + 2*dr*(psi[ind] - 1)/r; }

inline double fdaR_resbeta(const vector<double>& beta, int ind, double dr, double r) {
  return d_b(beta,ind) + 2*dr*beta[ind]/r; }

inline double fdaR_resalpha(const vector<double>& alpha, int ind, double dr, double r) {
  return d_b(alpha,ind) + 2*dr*(alpha[ind] - 1)/r; }

// set res_vals to be residual = L[f] for rhs of jacobian.(-delta) = residual
void get_ell_res(double *res_vals, const vector<double>& xi, const vector<double>& pi,
		 const vector<double>& alpha, const vector<double>& beta,
		 const vector<double>& psi, int lastpt, int N, double dr, double r) {
  res_vals[0] = neumann0res(alpha);
  res_vals[1] = dirichlet0res(beta);
  res_vals[2] = neumann0res(psi);
  for (int k = 1; k < lastpt; ++k) {
    r += dr;
    // MAYBE NEED TO CHANGE BACK TO NO sq(dr)
    res_vals[3*k] = fda_resalpha(xi, pi, alpha, beta, psi, k, dr, r);
    res_vals[3*k+1] = fda_resbeta(xi, pi, alpha, beta, psi, k, dr, r);
    res_vals[3*k+2] = fda_respsi(xi, pi, alpha, beta, psi, k, dr, r);
  }
  r += dr;
  res_vals[N-3] = fdaR_resalpha(alpha, lastpt, dr, r);
  res_vals[N-2] = fdaR_resbeta(beta, lastpt, dr, r);
  res_vals[N-1] = fdaR_respsi(psi, lastpt, dr, r);
  return;
}

// ***********************  JACOBIAN  ***********************

// ***********************  row alpha(ind)   ***********************

inline double jac_aa(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -2 - 8*M_PI*sq(dr)*sq(pi[ind]) +
    2*pw4(psi[ind])*sq(0.5*d_c(beta,ind) - dr*beta[ind]/r) / (3*sq(alpha[ind])); }

inline double jac_aa_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return 1 + p_m*(0.5*d_c(psi,ind)/psi[ind] + dr/r); }

inline double jac_ab(const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return 4*dr*pw4(psi[ind])*(0.5*d_c(beta,ind) - dr*beta[ind]/r) / (3*alpha[ind]*r); }

inline double jac_ab_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return -p_m*2*pw4(psi[ind])*(0.5*d_c(beta,ind) - dr*beta[ind]/r) / (3*alpha[ind]); }

inline double jac_ap(const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -0.5*d_c(alpha,ind)*d_c(psi,ind) / sq(psi[ind]) -
	     8*pw3(psi[ind])*sq(0.5*d_c(beta,ind) - dr*beta[ind]/r) / (3*alpha[ind]); }

inline double jac_ap_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return p_m*0.5*d_c(alpha,ind) / psi[ind]; }


// ***********************  row beta(ind)   ***********************

inline double jac_ba(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return 12*M_PI*sq(dr)*xi[ind]*pi[ind] / sq(psi[ind]) +
    0.5*d_c(alpha,ind)*(0.5*d_c(beta,ind) - dr*beta[ind]/r) / sq(alpha[ind]); }

inline double jac_ba_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return -p_m*0.5*(0.5*d_c(beta,ind) - dr*beta[ind]/r) / alpha[ind]; }

inline double jac_bb(const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -2*(1 + sq(dr)/sq(r)) - 0.5*dr*(d_c(alpha,ind)/alpha[ind] + 6*d_c(psi,ind)/psi[ind])/r; }

inline double jac_bb_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return 1 + p_m*0.25*(4*dr/r - d_c(alpha,ind)/alpha[ind] + 6*d_c(psi,ind)/psi[ind]); }

inline double jac_bp(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -6*(4*M_PI*sq(dr)*xi[ind]*pi[ind]*alpha[ind] / pw3(psi[ind]) +
	     0.5*d_c(psi,ind)*(0.5*d_c(beta,ind) - dr*beta[ind]/r) / sq(psi[ind])); }

inline double jac_bp_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return p_m*3*(0.5*d_c(beta,ind) - dr*beta[ind]/r) / psi[ind]; }

// ***********************  row psi(ind)   ***********************

inline double jac_pa(const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -dr*pw5(psi[ind])*(0.5*d_c(beta,ind) - dr*beta[ind]/r) / (6*pw3(alpha[ind])); }

inline double jac_pa_pm() { return 0; }

inline double jac_pb(const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -sq(dr)*pw5(psi[ind]) / (12*sq(alpha[ind])*r); }

inline double jac_pb_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return p_m*dr*pw5(psi[ind]) / (24*sq(alpha[ind])); }

inline double jac_pp(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -2 + M_PI*sq(dr)*(sq(xi[ind]) + sq(pi[ind])) +
    5*dr*pw4(psi[ind])*(0.5*d_c(beta,ind) - dr*beta[ind]/r) / (12*sq(alpha[ind])); }

inline double jac_pp_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return 1 + p_m*dr/r; }

inline void set_inds(int *cols, int indalpha) {
  cols[0] = indalpha-3, cols[1] = indalpha-2, cols[2] = indalpha-1,
    cols[3] = indalpha, cols[4] = indalpha+1, cols[5] = indalpha+2,
    cols[6] = indalpha+3, cols[7] = indalpha+4, cols[8] = indalpha+5;
  return; }

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
  double ell_tol = tol; // will not auto change if tol changed
  double petsc_tol = tol; // will not auto change if tol changed
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
      {"-petsc_tol",&petsc_tol},
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
  ierr = KSPSetTolerances(ksp, petsc_tol, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(ierr);
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
    set_rhs(bxi, bpi, old_xi, old_pi, alpha, beta, psi, lam, dr, r, j, lastpt);
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
		  j, lastpt, somm_cond, somm_coeff);
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
// ***********************************************************************
// ***********************************************************************
// *********************     ERROR BELOW HERE     ************************
// ***********************************************************************
// ***********************************************************************
// ***********************************************************************


      
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
// ***********************************************************************
// *********************     ERROR ABOVE HERE     ************************
// ***********************************************************************
// ***********************************************************************
// ***********************************************************************



    
// ***********************************************************************
// ***********************************************************************
// ****************** FULL ITERATIVE SOLUTION COMPLETE *******************
// ***********************************************************************
// ***********************************************************************
    
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
