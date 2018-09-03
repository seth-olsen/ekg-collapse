#include <iostream>
#include <algorithm> // for max_element()
#include <fstream> // for mass/itn files
#include <ctime> // for quick time analysis
#include "fda-fns.h"
#include "fda-io.h"

#define _USE_MATH_DEFINES
#ifndef M_PI
cerr << "NO PI\n";
const double M_PI = 4.0*atan(1.0);
#endif

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
  return r*pw6(psi[0])*sq(r*ddr_f(beta,0,dr) - beta[0]) / (18*sq(alpha[0])) -
    2*sq(r)*ddr_f(psi,0,dr)*(psi[0] + r*ddr_f(psi,0,dr)); }

inline double mass_aspectR(const vector<double>& alpha, const vector<double>& beta,
			  const vector<double>& psi, int k, double dr, double r) {
  return r*pw6(psi[k])*sq(r*ddr_b(beta,k,dr) - beta[k]) / (18*sq(alpha[0])) -
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

// perform gauss-seidel update
inline void gs_update(const vector<double>& bxi, const vector<double>& bpi,
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
  

inline double fda_psi(const vector<double>& xi, const vector<double>& pi,
		      const vector<double>& alpha, const vector<double>& beta,
		      const vector<double>& psi, int ind, double lam, double dr, double r) {
  return 0.25*lam*(beta[ind]*(psi[ind+1] - psi[ind-1]) +
		   psi[ind]*(beta[ind+1] + 4*beta[ind]/r - beta[ind-1])); }

inline double fda_respsi(const vector<double>& xi, const vector<double>& pi,
			 const vector<double>& alpha, const vector<double>& beta,
			 const vector<double>& psi, int ind, double dr, double r) {
  return d2_c(psi,ind) + dr*d_c(psi,ind)/r +
    pw5(psi[ind])*sq(0.5*d_c(beta,ind) - dr*beta[ind]/r) / (12*sq(alpha[ind])) + 
    M_PI*sq(dr)*psi[ind]*(sq(xi[ind]) + sq(pi[ind])); }

inline double fda_resbeta(const vector<double>& xi, const vector<double>& pi,
			  const vector<double>& alpha, const vector<double>& beta,
			  const vector<double>& psi, int ind, double dr, double r) {
  return  d2_c(beta,ind) + 12*M_PI*sq(dr)*alpha[ind]*xi[ind]*pi[ind] / sq(psi[ind])
    + (2*dr/r + 3*d_c(psi,ind)/psi[ind] - 0.5*d_c(alpha,ind)/alpha[ind])*
    (0.5*d_c(beta,ind) - dr*beta[ind]/r); }

inline double fda_resalpha(const vector<double>& xi, const vector<double>& pi,
			   const vector<double>& alpha, const vector<double>& beta,
			   const vector<double>& psi, int ind, double dr, double r) {
  return d2_c(alpha,ind) + d_c(alpha,ind)*(dr/r + 0.5*d_c(psi,ind)/psi[ind])
    - 2*pw4(psi[ind])*sq(0.5*d_c(beta,ind) - dr*beta[ind]/r) / (3*alpha[ind])
    - 8*M_PI*sq(dr)*alpha[ind]*sq(pi[ind]) ; }

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
  return -pw5(psi[ind])*sq(0.5*d_c(beta,ind) - dr*beta[ind]/r) / (6*pw3(alpha[ind])); }

inline double jac_pa_pm() { return 0; }

inline double jac_pb(const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -dr*pw5(psi[ind])*(0.5*d_c(beta,ind) - dr*beta[ind]/r) / (6*sq(alpha[ind])*r); }

inline double jac_pb_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return p_m*pw5(psi[ind])*(0.5*d_c(beta,ind) - dr*beta[ind]/r) / (12*sq(alpha[ind])); }

inline double jac_pp(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -2 + M_PI*sq(dr)*(sq(xi[ind]) + sq(pi[ind])) +
    5*pw4(psi[ind])*sq(0.5*d_c(beta,ind) - dr*beta[ind]/r) / (12*sq(alpha[ind])); }

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
