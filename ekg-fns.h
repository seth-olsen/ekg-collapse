#include "fda-fns.h"
#include "fda-io.h"

// multiply by 4*M_PI*r^2 for dm/dr of scalar
inline double dmdr_scalar(double xival, double pival, double alphaval,
			  double betaval, double psival) {
  return pw4(psival)*betaval*xival*pival + 0.5*sq(psival)*alphaval*(sq(xival) + sq(pival)); }
////////////////////////////////////////////////////////////////////////////////////////////////
// mass aspect function
////////////////////////////////////////////////////////////////////////////////////////////////
inline double mass_aspect(const vector<double>& alpha, const vector<double>& beta,
			  const vector<double>& psi, int k, double dr, double r) {
  return r*pw6(psi[k])*sq(r*ddr_c(beta,k,dr) - beta[k]) / (18*sq(alpha[k])) -
    2*sq(r)*ddr_c(psi,k,dr)*(psi[k] + r*ddr_c(psi,k,dr)); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double mass_aspect0(const vector<double>& alpha, const vector<double>& beta,
			  const vector<double>& psi, double dr, double r) {
  return r*pw6(psi[0])*sq(r*ddr_f(beta,0,dr) - beta[0]) / (18*sq(alpha[0])) -
    2*sq(r)*ddr_f(psi,0,dr)*(psi[0] + r*ddr_f(psi,0,dr)); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double mass_aspectR(const vector<double>& alpha, const vector<double>& beta,
			  const vector<double>& psi, int k, double dr, double r) {
  return r*pw6(psi[k])*sq(r*ddr_b(beta,k,dr) - beta[k]) / (18*sq(alpha[0])) -
    2*sq(r)*ddr_b(psi,k,dr)*(psi[k] + r*ddr_b(psi,k,dr)); }
////////////////////////////////////////////////////////////////////////////////////////////////
// field = oldfield + fda(field) + fda(oldfield)
////////////////////////////////////////////////////////////////////////////////////////////////
inline double fda_xi(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double lam)
{ return ( 0.25*lam* dx_ekg_c(beta, xi, alpha, pi, psi, ind) ); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double fda_pi(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double lam, double dr, double r) {
  return ( 0.25*lam*p4in(psi[ind])* sqin(r)*dp_ekg_c(beta, pi, alpha, xi, psi, ind, dr, r)
	   - (lam/3.0)* pi[ind]*( d_c(beta, ind) + beta[ind]*(6*d_c(psi, ind)/psi[ind] + 4*dr/r) ) ); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double fda_psi(const vector<double>& xi, const vector<double>& pi,
		      const vector<double>& alpha, const vector<double>& beta,
		      const vector<double>& psi, int ind, double lam, double dr, double r) {
  return 0.25*lam*(beta[ind]*d_c(psi,ind) + psi[ind]*(4*dr*beta[ind]/r + d_c(beta,ind))/6.0); }
////////////////////////////////////////////////////////////////////////////////////////////////
// only need these if no BC at rmin
////////////////////////////////////////////////////////////////////////////////////////////////
inline double fda0_xi(const vector<double>& xi, const vector<double>& pi,
		      const vector<double>& alpha, const vector<double>& beta,
		      const vector<double>& psi, double lam)
{ return ( 0.25 * lam * dx_ekg_f(beta, xi, alpha, pi, psi, 0) ); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double fda0_pi(const vector<double>& xi, const vector<double>& pi,
		      const vector<double>& alpha, const vector<double>& beta,
		      const vector<double>& psi, double lam, double dr, double r)
{ return ( 0.25*lam*p4in(psi[0])* sqin(r)*dp_ekg_f(beta, pi, alpha, xi, psi, 0, dr, r)
	   - (lam/3.0)* pi[0]*( d_f(beta, 0) + beta[0]*(6*d_f(psi, 0)/psi[0] + 4*dr/r) ) ); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double fda0_psi(const vector<double>& xi, const vector<double>& pi,
		       const vector<double>& alpha, const vector<double>& beta,
		       const vector<double>& psi, double lam, double dr, double r) {
  return 0.25*lam*(beta[0]*d_f(psi,0) + psi[0]*(4*dr*beta[0]/r + d_f(beta,0))/6.0); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double fda_respsi(const vector<double>& xi, const vector<double>& pi,
			 const vector<double>& alpha, const vector<double>& beta,
			 const vector<double>& psi, int ind, double dr, double r) {
  return d2_c(psi,ind) + dr*d_c(psi,ind)/r +
    pw5(psi[ind])*sq(0.5*d_c(beta,ind) - dr*beta[ind]/r) / (12*sq(alpha[ind])) + 
    M_PI*sq(dr)*psi[ind]*(sq(xi[ind]) + sq(pi[ind])); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double fda_resbeta(const vector<double>& xi, const vector<double>& pi,
			  const vector<double>& alpha, const vector<double>& beta,
			  const vector<double>& psi, int ind, double dr, double r) {
  return  d2_c(beta,ind) + 12*M_PI*sq(dr)*alpha[ind]*xi[ind]*pi[ind] / sq(psi[ind])
    + (2*dr/r + 3*d_c(psi,ind)/psi[ind] - 0.5*d_c(alpha,ind)/alpha[ind])*
    (0.5*d_c(beta,ind) - dr*beta[ind]/r); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double fda_resalpha(const vector<double>& xi, const vector<double>& pi,
			   const vector<double>& alpha, const vector<double>& beta,
			   const vector<double>& psi, int ind, double dr, double r) {
  return d2_c(alpha,ind) + d_c(alpha,ind)*(dr/r + 0.5*d_c(psi,ind)/psi[ind])
    - 2*pw4(psi[ind])*sq(0.5*d_c(beta,ind) - dr*beta[ind]/r) / (3*alpha[ind])
    - 8*M_PI*sq(dr)*alpha[ind]*sq(pi[ind]) ; }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double fdaR_respsi(const vector<double>& psi, int ind, double dr, double r) {
  return d_b(psi,ind) + 2*dr*(psi[ind] - 1)/r; }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double fdaR_resbeta(const vector<double>& beta, int ind, double dr, double r) {
  return d_b(beta,ind) + 2*dr*beta[ind]/r; }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double fdaR_resalpha(const vector<double>& alpha, int ind, double dr, double r) {
  return d_b(alpha,ind) + 2*dr*(alpha[ind] - 1)/r; }
////////////////////////////////////////////////////////////////////////////////////////////////
// **********************************************************
// **********************************************************
//             INITIAL AND BOUNDARY CONDITIONS
// **********************************************************
// **********************************************************
////////////////////////////////////////////////////////////////////////////////////////////////
inline double ic_alpha(double r, double r2m) { return 1.0; }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double ic_beta(double r, double r2m) { return 0; }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double ic_psi(double r, double r2m) { return 1.0; }
////////////////////////////////////////////////////////////////////////////////////////////////
// for gaussian field or sin(coeff*r)*cos(coeff*t)/(coeff*r)
inline double ic_sol(double r, double amp, double dsq, double r0)
{ return amp * exp(-(r - r0)*(r - r0)/dsq); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double ic_xi(double r, double amp, double dsq, double r0)
{ return -2 * (r - r0) * amp * exp(-(r - r0)*(r - r0)/dsq) / dsq; }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double ic_pi(double r, double amp, double dsq, double r0, double ctot, double cxi)
{ return ctot*(cxi*ic_xi(r, amp, dsq, r0) + ic_sol(r, amp, dsq, r0)/r); }
////////////////////////////////////////////////////////////////////////////////////////////////
// ***********************  JACOBIAN  ***********************
////////////////////////////////////////////////////////////////////////////////////////////////
// ***********************  row alpha(ind)   ***********************
////////////////////////////////////////////////////////////////////////////////////////////////
inline double jac_aa(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -2 - 8*M_PI*sq(dr)*sq(pi[ind]) +
    2*pw4(psi[ind])*sq(0.5*d_c(beta,ind) - dr*beta[ind]/r) / (3*sq(alpha[ind])); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double jac_aa_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return 1 + p_m*(0.5*d_c(psi,ind)/psi[ind] + dr/r); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double jac_ab(const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return 4*dr*pw4(psi[ind])*(0.5*d_c(beta,ind) - dr*beta[ind]/r) / (3*alpha[ind]*r); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double jac_ab_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return -p_m*2*pw4(psi[ind])*(0.5*d_c(beta,ind) - dr*beta[ind]/r) / (3*alpha[ind]); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double jac_ap(const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -0.5*d_c(alpha,ind)*d_c(psi,ind) / sq(psi[ind]) -
	     8*pw3(psi[ind])*sq(0.5*d_c(beta,ind) - dr*beta[ind]/r) / (3*alpha[ind]); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double jac_ap_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return p_m*0.5*d_c(alpha,ind) / psi[ind]; }
////////////////////////////////////////////////////////////////////////////////////////////////
// ***********************  row beta(ind)   ***********************
////////////////////////////////////////////////////////////////////////////////////////////////
inline double jac_ba(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return 12*M_PI*sq(dr)*xi[ind]*pi[ind] / sq(psi[ind]) +
    0.5*d_c(alpha,ind)*(0.5*d_c(beta,ind) - dr*beta[ind]/r) / sq(alpha[ind]); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double jac_ba_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return -p_m*0.5*(0.5*d_c(beta,ind) - dr*beta[ind]/r) / alpha[ind]; }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double jac_bb(const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -2*(1 + sq(dr)/sq(r)) - 0.5*dr*(d_c(alpha,ind)/alpha[ind] + 6*d_c(psi,ind)/psi[ind])/r; }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double jac_bb_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return 1 + p_m*0.25*(4*dr/r - d_c(alpha,ind)/alpha[ind] + 6*d_c(psi,ind)/psi[ind]); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double jac_bp(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -6*(4*M_PI*sq(dr)*xi[ind]*pi[ind]*alpha[ind] / pw3(psi[ind]) +
	     0.5*d_c(psi,ind)*(0.5*d_c(beta,ind) - dr*beta[ind]/r) / sq(psi[ind])); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double jac_bp_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return p_m*3*(0.5*d_c(beta,ind) - dr*beta[ind]/r) / psi[ind]; }
////////////////////////////////////////////////////////////////////////////////////////////////
// ***********************  row psi(ind)   ***********************
////////////////////////////////////////////////////////////////////////////////////////////////
inline double jac_pa(const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -pw5(psi[ind])*sq(0.5*d_c(beta,ind) - dr*beta[ind]/r) / (6*pw3(alpha[ind])); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double jac_pa_pm() { return 0; }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double jac_pb(const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -dr*pw5(psi[ind])*(0.5*d_c(beta,ind) - dr*beta[ind]/r) / (6*sq(alpha[ind])*r); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double jac_pb_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return p_m*pw5(psi[ind])*(0.5*d_c(beta,ind) - dr*beta[ind]/r) / (12*sq(alpha[ind])); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double jac_pp(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -2 + M_PI*sq(dr)*(sq(xi[ind]) + sq(pi[ind])) +
    5*pw4(psi[ind])*sq(0.5*d_c(beta,ind) - dr*beta[ind]/r) / (12*sq(alpha[ind])); }
////////////////////////////////////////////////////////////////////////////////////////////////
inline double jac_pp_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return 1 + p_m*dr/r; }
