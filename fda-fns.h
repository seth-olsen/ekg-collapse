#include <vector> // for everything
#include <cmath> // for ICs

using namespace std;

// x^2 and 1/x^2
inline double sq(double x) { return (x * x); }
inline double pw3(double x) { return x*sq(x); }
inline double pw4(double x) { return sq(sq(x)); }
inline double pw5(double x) { return x*pw4(x); }
inline double pw6(double x) { return sq(x)*pw4(x); }
inline double sqin(double x) { return (1.0 / sq(x)); }
inline double p3in(double x) { return (1.0 / pw3(x)); }
inline double p4in(double x) { return (1.0 / pw4(x)); }
inline double p5in(double x) { return (1.0 / pw5(x)); }
inline double p6in(double x) { return (1.0 / pw6(x)); }


// ******** DIFFERENCES AND DERIVATIVES ********

inline double d_c(const vector<double>& u, int ind)
{ return u[ind+1] - u[ind-1]; }

inline double d_f(const vector<double>& u, int ind)
{ return -3*u[ind] + 4*u[ind+1] - u[ind+2]; }

inline double d_b(const vector<double>& u, int ind)
{ return 3*u[ind] - 4*u[ind-1] + u[ind-2]; }

inline double ddr_c(const vector<double>& u, int ind, double dr)
{ return 0.5 * d_c(u, ind) / dr; }

inline double ddr_f(const vector<double>& u, int ind, double dr)
{ return 0.5 * d_f(u, ind) / dr; }

inline double ddr_b(const vector<double>& u, int ind, double dr)
{ return 0.5 * d_b(u, ind) / dr; }

inline double d2_c(const vector<double>& u, int ind)
{ return u[ind-1] - 2*u[ind] + u[ind+1]; }

inline double d2_f(const vector<double>& u, int ind)
{ return 2*u[ind] - 5*u[ind+1] + 4*u[ind+2] - u[ind+3]; }

inline double d2_b(const vector<double>& u, int ind)
{ return 2*u[ind] - 5*u[ind-1] + 4*u[ind-2] - u[ind-3]; }

inline double ddr2_c(const vector<double>& u, int ind, double dr)
{ return d2_c(u, ind) / sq(dr); }

inline double ddr2_f(const vector<double>& u, int ind, double dr)
{ return d2_f(u, ind) / sq(dr); }

inline double ddr2_b(const vector<double>& u, int ind, double dr)
{ return d2_b(u, ind) / sq(dr); }


// ******* DIFFERENCE OPERATORS FOR EKG-COLLAPSE

// centered differencing operator for fn1*field1 + fn2*field2/fn3^2
// --> mult by 1/(2*dr) to get d(fn1*field1 + fn2*field2/fn3^2)/dr at O(dr^2)
inline double dx_ekg_c(const vector<double>& fn1, const vector<double>& field1,
		    const vector<double>& fn2, const vector<double>& field2,
		    const vector<double>& fn3, int ind)
{ return fn1[ind+1]*field1[ind+1] + fn2[ind+1]*field2[ind+1]*sqin(fn3[ind+1])
    - fn1[ind-1]*field1[ind-1] - fn2[ind-1]*field2[ind-1]*sqin(fn3[ind-1]); }

// centered differencing operator for r^2*fn3^4*(fn1*field1 + fn2*field2/fn3^2)
// --> mult by 1/(2*dr) to get d(r^2*fn3^4*(fn1*field1 + fn2*field2/fn3^2))/dr at O(dr^2)
inline double dp_ekg_c(const vector<double>& fn1, const vector<double>& field1,
		    const vector<double>& fn2, const vector<double>& field2,
		    const vector<double>& fn3, int ind, double dr, double r)
{ return sq(r+dr)*pw4(fn3[ind+1])*(fn1[ind+1]*field1[ind+1] + fn2[ind+1]*field2[ind+1]*sqin(fn3[ind+1]))
    - sq(r-dr)*pw4(fn3[ind-1])*(fn1[ind-1]*field1[ind-1] + fn2[ind-1]*field2[ind-1]*sqin(fn3[ind-1])); }

// forward differencing operator for dn1*field1 + fn2*field2/fn3^2
// --> mult by 1/(2*dr) to get d(fn1*field1 + fn2*field2/fn3^2)/dr at O(dr^2)
inline double dx_ekg_f(const vector<double>& fn1, const vector<double>& field1,
		    const vector<double>& fn2, const vector<double>& field2,
		    const vector<double>& fn3, int ind)
{ return -3*(fn1[ind]*field1[ind] + fn2[ind]*field2[ind]*sqin(fn3[ind]))
    + 4*(fn1[ind+1]*field1[ind+1] + fn2[ind+1]*field2[ind+1]*sqin(fn3[ind+1]))
    - (fn1[ind+2]*field1[ind+2] + fn2[ind+2]*field2[ind+2]*sqin(fn3[ind+2])); }

// forward differencing operator for r^2*fn3^4*(fn1*field1 + fn2*field2/fn3^2)
// --> mult by 1/(2*dr) to get d(r^2*fn3^4*(fn1*field1 + fn2*field2/fn3^2))/dr at O(dr^2)
inline double dp_ekg_f(const vector<double>& fn1, const vector<double>& field1,
		    const vector<double>& fn2, const vector<double>& field2,
		    const vector<double>& fn3, int ind, double dr, double r)
{ return -3*sq(r)*pw4(fn3[ind])*(fn1[ind]*field1[ind] + fn2[ind]*field2[ind]*sqin(fn3[ind]))
    + 4*sq(r+dr)*pw4(fn3[ind+1])*(fn1[ind+1]*field1[ind+1] + fn2[ind+1]*field2[ind+1]*sqin(fn3[ind+1]))
    - sq(r+2*dr)*pw4(fn3[ind+2])*(fn1[ind+2]*field1[ind+2] + fn2[ind+2]*field2[ind+2]*sqin(fn3[ind+2])); }

// centered differencing operator for r^2*fn3^4*(fn1*field1 + fn2*field2) [using d/d(r^3)]
// --> mult by 1/(2*dr) to get (1/r^2)*d(r^2*fn3^4*(fn1*field1 + fn2*field2))/dr at O(dr^2)
inline double d3p_ekg_c(const vector<double>& fn1, const vector<double>& field1,
		     const vector<double>& fn2, const vector<double>& field2,
		     const vector<double>& fn3, int ind, double dr, double r)
{ return ( (sq(r+dr)*pw4(fn3[ind+1])*(fn1[ind+1]*field1[ind+1] + fn2[ind+1]*field2[ind+1]*sqin(fn3[ind+1]))
	  - sq(r-dr)*pw4(fn3[ind-1])*(fn1[ind-1]*field1[ind-1] + fn2[ind-1]*field2[ind-1]*sqin(fn3[ind-1])))
	  * 3.0 / (3*r*r + dr*dr) ); }

// forward differencing operator for r^2*fn3^4*(fn1*field1 + fn2*field2) [using d/d(r^3)]
// --> mult by 1/(2*dr) to get (1/r^2)*d(r^2*fn3^4*(fn1*field1 + fn2*field2))/dr at O(dr^2)
inline double d3p_ekg_f(const vector<double>& fn1, const vector<double>& field1,
		     const vector<double>& fn2, const vector<double>& field2,
		     const vector<double>& fn3, int ind, double dr, double r)
{ return ( (-3*sq(r)*pw4(fn3[ind])*(fn1[ind]*field1[ind] + fn2[ind]*field2[ind]*sqin(fn3[ind]))
	    + 4*sq(r+dr)*pw4(fn3[ind+1])*(fn1[ind+1]*field1[ind+1] + fn2[ind+1]*field2[ind+1]*sqin(fn3[ind+1]))
	    - sq(r+2*dr)*pw4(fn3[ind+2])*(fn1[ind+2]*field1[ind+2] + fn2[ind+2]*field2[ind+2]*sqin(fn3[ind+2])))
	   * 3.0 / (3*r*r - 2*dr*dr) ); }



// ********* DIFFERENCE OPERATORS FROM MASSLESS-SCALAR

// centered differencing operator for fn1*field1 + fn2*field2
// --> mult by 1/(2*dr) to get d(fn1*field1 + fn2*field2)/dr at O(dr^2)
inline double dif_c(const vector<double>& fn1, const vector<double>& field1,
		    const vector<double>& fn2, const vector<double>& field2,
		    int ind)
{ return fn1[ind+1]*field1[ind+1] + fn2[ind+1]*field2[ind+1]
    - fn1[ind-1]*field1[ind-1] - fn2[ind-1]*field2[ind-1]; }

// centered differencing operator for r^2*(fn1*field1 + fn2*field2)
// --> mult by 1/(2*dr) to get d(r^2*(fn1*field1 + fn2*field2))/dr at O(dr^2)
inline double r2dif_c(const vector<double>& fn1, const vector<double>& field1,
		      const vector<double>& fn2, const vector<double>& field2,
		      int ind, double dr, double r)
{ return sq(r+dr)*(fn1[ind+1]*field1[ind+1] + fn2[ind+1]*field2[ind+1])
    - sq(r-dr)*(fn1[ind-1]*field1[ind-1] + fn2[ind-1]*field2[ind-1]); }

// forward differencing operator for dn1*field1 + fn2*field2
// --> mult by 1/(2*dr) to get d(fn1*field1 + fn2*field2)/dr at O(dr^2)
inline double dif_f(const vector<double>& fn1, const vector<double>& field1,
		    const vector<double>& fn2, const vector<double>& field2,
		    int ind)
{ return -3*(fn1[ind]*field1[ind] + fn2[ind]*field2[ind])
    + 4*(fn1[ind+1]*field1[ind+1] + fn2[ind+1]*field2[ind+1])
    - (fn1[ind+2]*field1[ind+2] + fn2[ind+2]*field2[ind+2]); }

// forward differencing operator for r^2*(fn1*field1 + fn2*field2)
// --> mult by 1/(2*dr) to get d(r^2*(fn1*field1 + fn2*field2))/dr at O(dr^2)
inline double r2dif_f(const vector<double>& fn1, const vector<double>& field1,
		      const vector<double>& fn2, const vector<double>& field2,
		      int ind, double dr, double r)
{ return -3*sq(r)*(fn1[ind]*field1[ind] + fn2[ind]*field2[ind])
    + 4*sq(r+dr)*(fn1[ind+1]*field1[ind+1] + fn2[ind+1]*field2[ind+1])
    - sq(r+2*dr)*(fn1[ind+2]*field1[ind+2] + fn2[ind+2]*field2[ind+2]); }

// centered differencing operator for r^2*(fn1*field1 + fn2*field2) [using d/d(r^3)]
// --> mult by 1/(2*dr) to get (1/r^2)*d(r^2*(fn1*field1 + fn2*field2))/dr at O(dr^2)
inline double r2d3_c(const vector<double>& fn1, const vector<double>& field1,
		     const vector<double>& fn2, const vector<double>& field2,
		     int ind, double dr, double r)
{ return ( (sq(r+dr)*(fn1[ind+1]*field1[ind+1] + fn2[ind+1]*field2[ind+1])
	  - sq(r-dr)*(fn1[ind-1]*field1[ind-1] + fn2[ind-1]*field2[ind-1]))
	  * 3.0 / (3*r*r + dr*dr) ); }

// forward differencing operator for r^2*(fn1*field1 + fn2*field2) [using d/d(r^3)]
// --> mult by 1/(2*dr) to get (1/r^2)*d(r^2*(fn1*field1 + fn2*field2))/dr at O(dr^2)
inline double r2d3_f(const vector<double>& fn1, const vector<double>& field1,
		     const vector<double>& fn2, const vector<double>& field2,
		     int ind, double dr, double r)
{ return ( (-3*sq(r)*(fn1[ind]*field1[ind] + fn2[ind]*field2[ind])
	    + 4*sq(r+dr)*(fn1[ind+1]*field1[ind+1] + fn2[ind+1]*field2[ind+1])
	    - sq(r+2*dr)*(fn1[ind+2]*field1[ind+2] + fn2[ind+2]*field2[ind+2]))
	   * 3.0 / (3*r*r - 2*dr*dr) ); }


// d(ua/ub^2)/dr * 2*dr
inline double dab2_c(const vector<double>& ua, const vector<double>& ub, int ind) {
  return ( ua[ind+1]*sqin(ub[ind+1]) - ua[ind-1]*sqin(ub[ind-1]) ); }

inline double dab2_f(const vector<double>& ua, const vector<double>& ub, int ind) {
  return ( -3*ua[ind]*sqin(ub[ind]) + 4*ua[ind+1]*sqin(ub[ind+1]) - ua[ind+2]*sqin(ub[ind+2]) ); }

// d(r^2*ua^4)/dr * 2*dr
inline double dr2a4_c(const vector<double>& ua, int ind, double dr, double r) {
  return ( sq(r+dr)*pw4(ua[ind+1]) - sq(r-dr)*pw4(ua[ind-1]) ); }

inline double dr2a4_f(const vector<double>& ua, int ind, double dr, double r) {
  return ( -3*sq(r)*pw4(ua[ind]) + 4*sq(r+dr)*pw4(ua[ind+1]) - sq(r+2*dr)*pw4(ua[ind+2]) ); }



// *********** BOUNDARY CONDITIONS ****************


inline void dirichlet0(vector<double>& field) {
  field[0] = 0;
  return;
}

inline double dirichlet0res(const vector<double>& field) {;
  return field[0];
}

inline void neumann0(vector<double>& field) {
  field[0] = (4*field[1] - field[2]) / 3.0;
  return;
}

inline double neumann0res(const vector<double>& field) {
  return -3*field[0] + 4*field[1] - field[2];
}

inline void sommerfeld(vector<double>& field, const vector<double>& oldfield,
		   int ind, double lambda, double coeff) {
      field[ind] = (lambda*(field[ind-1] + oldfield[ind-1])
		 - 0.25*lambda*(field[ind-2] + oldfield[ind-2])
		 + (1 - coeff)*oldfield[ind]) / (1 + coeff);
  return;
}

inline double sommerfeldres(const vector<double>& field, const vector<double>& oldfield,
		   int ind, double lambda, double coeff) {
  return (coeff+1)*field[ind] + (coeff-1)*oldfield[ind]
    - lambda*(field[ind-1] + oldfield[ind-1])
    + 0.25*lambda*(field[ind-2] + oldfield[ind-2]);
}


// ********** DISSIPATION FUNCTIONS **************


// kreiss-oliger dissipation (p.23 choptuik notes)
inline double dissipate(double eps, const vector<double>& u, int ind)
{ return -eps * 0.0625 * ( u[ind-2] - 4*u[ind-1] + 6*u[ind]
			   - 4*u[ind+1] + u[ind+2] ); }

inline double symdiss1(double eps, const vector<double>& u)
{ return -eps * 0.0625 * ( u[3] - 4*u[2] + 7*u[1] - 4*u[0] ); }

inline double antidiss1(double eps, const vector<double>& u)
{ return -eps * 0.0625 * ( u[3] - 4*u[2] + 5*u[1] - 4*u[0] ); }

inline double symdiss0(double eps, const vector<double>& u)
{ return -eps * 0.0625 * ( 2*u[2] - 8*u[1] + 6*u[0] ); }

inline double antidiss0(double eps, const vector<double>& u)
{ return -eps * 0.0625 * ( 6*u[0] ); }



// ******** IRES FUNCTIONS *************************


// centered ires_f1 = f1 - ires_c(f1, f2, ind, c1, d1, c2, d2)
//                   - oldf1 - ires_c(oldf1, oldf2, ind, c1, d1, c2, d2)
inline double iresxi_c(const vector<double>& xi, const vector<double>& pi,
		       const vector<double>& alpha, const vector<double>& beta,
		       const vector<double>& psi, int ind, double lam) {
  return 0.25*lam*(xi[ind]*d_c(beta, ind) + beta[ind]*d_c(xi, ind) +
		   pi[ind]*dab2_c(alpha, psi, ind) + alpha[ind]*sqin(psi[ind])*d_c(pi, ind)); }

// forward ires_f1 = f1 - ires_c(f1, f2, ind, c1, d1, c2, d2)
//                   - oldf1 - ires_c(oldf1, oldf2, ind, c1, d1, c2, d2)
inline double iresxi_f(const vector<double>& xi, const vector<double>& pi,
		       const vector<double>& alpha, const vector<double>& beta,
		       const vector<double>& psi, int ind, double lam) {
  return 0.25*lam*(xi[ind]*d_f(beta, ind) + beta[ind]*d_f(xi, ind) +
		   pi[ind]*dab2_f(alpha, psi, ind) + alpha[ind]*sqin(psi[ind])*d_f(pi, ind)); }

inline double irespi_c(const vector<double>& xi, const vector<double>& pi,
		       const vector<double>& alpha, const vector<double>& beta,
		       const vector<double>& psi, int ind, double lam, double dr, double r) {
  return 0.25*lam*(pi[ind]*d_c(beta, ind) + beta[ind]*d_c(pi, ind) +
		   xi[ind]*dab2_c(alpha, psi, ind) + alpha[ind]*sqin(psi[ind])*d_c(xi, ind) +
		   (pi[ind]*beta[ind] + xi[ind]*alpha[ind]*sqin(psi[ind]))*sqin(r*sq(psi[ind]))
		   *dr2a4_c(psi, ind, dr, r))
    - (lam/3.0)* pi[ind]*( d_c(beta, ind) + beta[ind]*(6*d_c(psi, ind)/psi[ind] + 4*dr/r) ); }

inline double irespi_f(const vector<double>& xi, const vector<double>& pi,
		       const vector<double>& alpha, const vector<double>& beta,
		       const vector<double>& psi, int ind, double lam, double dr, double r) {
  return 0.25*lam*(pi[ind]*d_f(beta, ind) + beta[ind]*d_f(pi, ind) +
		   xi[ind]*dab2_f(alpha, psi, ind) + alpha[ind]*sqin(psi[ind])*d_f(xi, ind) +
		   (pi[ind]*beta[ind] + xi[ind]*alpha[ind]*sqin(psi[ind]))*sqin(r*sq(psi[ind]))
		   *dr2a4_f(psi, ind, dr, r))
    - (lam/3.0)* pi[ind]*( d_f(beta, ind) + beta[ind]*(6*d_f(psi, ind)/psi[ind] + 4*dr/r) ); }


// ***********************  JACOBIAN  ***********************

// ***********************  row alpha(ind)   ***********************

inline double oldjac_aa(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -2/sq(dr) - 8*M_PI*sq(pi[ind]) +
    2*pw4(psi[ind])*sq(ddr_c(beta,ind,dr) - beta[ind]/r) / (3*sq(alpha[ind])); }

inline double oldjac_aa_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return sqin(dr) + p_m*(ddr_c(psi,ind,dr)/psi[ind] + 1/r)/dr; }

inline double oldjac_ab(const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return 4*pw4(psi[ind])*(ddr_c(beta,ind,dr) - beta[ind]/r) / (3*alpha[ind]*r); }

inline double oldjac_ab_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return -p_m*2*pw4(psi[ind])*(ddr_c(beta,ind,dr) - beta[ind]/r) / (3*alpha[ind]*dr); }

inline double oldjac_ap(const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -2*(ddr_c(alpha,ind,dr)*ddr_c(psi,ind,dr)*sqin(psi[ind]) +
	     4*pw3(psi[ind])*sq(ddr_c(beta,ind,dr) - beta[ind]/r) / (3*alpha[ind])); }

inline double oldjac_ap_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return p_m*ddr_c(alpha,ind, dr) / (dr*psi[ind]); }


// ***********************  row beta(ind)   ***********************

inline double oldjac_ba(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return 12*M_PI*xi[ind]*pi[ind] / sq(psi[ind]) +
    ddr_c(alpha,ind,dr)*(ddr_c(beta,ind,dr) - beta[ind]/r) / sq(alpha[ind]); }

inline double oldjac_ba_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return -p_m*0.5*(ddr_c(beta,ind,dr) - beta[ind]/r) / (dr*alpha[ind]); }

inline double oldjac_bb(const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -2/sq(dr) - (2/r - ddr_c(alpha,ind,dr)/alpha[ind] + 6*ddr_c(psi,ind,dr)/psi[ind])/r; }

inline double oldjac_bb_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return sqin(dr) + p_m*0.5*(2/r - ddr_c(alpha,ind,dr)/alpha[ind] + 6*ddr_c(psi,ind,dr)/psi[ind])/dr; }

inline double oldjac_bp(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -6*(4*M_PI*xi[ind]*pi[ind]*alpha[ind]*p3in(psi[ind]) +
	     ddr_c(psi,ind,dr)*sqin(psi[ind])*(ddr_c(beta,ind,dr) - beta[ind]/r)); }

inline double oldjac_bp_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return p_m*3*(ddr_c(beta,ind,dr) - beta[ind]/r) / (dr*psi[ind]); }

// ***********************  row psi(ind)   ***********************

inline double oldjac_pa(const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -pw5(psi[ind])*(ddr_c(beta,ind,dr) - beta[ind]/r) / (6*pw3(alpha[ind])); }

inline double oldjac_pa_pm() { return 0; }

inline double oldjac_pb(const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -pw5(psi[ind]) / (sq(alpha[ind])*12*r); }

inline double oldjac_pb_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return p_m*pw5(psi[ind])*sqin(alpha[ind]) / (24*dr); }

inline double oldjac_pp(const vector<double>& xi, const vector<double>& pi,
		     const vector<double>& alpha, const vector<double>& beta,
		     const vector<double>& psi, int ind, double dr, double r) {
  return -2/sq(dr) + M_PI*(sq(xi[ind]) + sq(pi[ind])) +
    5*pw4(psi[ind])*(ddr_c(beta,ind,dr) - beta[ind]/r) / (12*sq(alpha[ind])); }

inline double oldjac_pp_pm(const vector<double>& alpha, const vector<double>& beta,
			const vector<double>& psi, int ind, int p_m, double dr, double r) {
  return sqin(dr) + p_m/(r*dr); }

// OLD VERSIONS DON'T HAVE 1/dr^2 FACTORED OUT

inline double oldfda_respsi(const vector<double>& xi, const vector<double>& pi,
			 const vector<double>& alpha, const vector<double>& beta,
			 const vector<double>& psi, int ind, double dr, double r) {
  return ddr2_c(psi,ind,dr) + 2*ddr_c(psi,ind,dr)/r +
    (ddr_c(beta,ind,dr) - beta[ind]/r)*pw5(psi[ind]) / (12*sq(alpha[ind])) + 
    M_PI*(sq(xi[ind]) + sq(pi[ind]))*psi[ind] ; }

inline double oldfda_resbeta(const vector<double>& xi, const vector<double>& pi,
			  const vector<double>& alpha, const vector<double>& beta,
			  const vector<double>& psi, int ind, double dr, double r) {
  return  ddr2_c(beta,ind,dr) + 12*M_PI*sqin(psi[ind])*alpha[ind]*xi[ind]*pi[ind]
    + (2/r + 6*ddr_c(psi,ind,dr)/psi[ind] - ddr_c(alpha,ind,dr)/alpha[ind])*
    (ddr_c(beta,ind,dr) - beta[ind]/r); }

inline double oldfda_resalpha(const vector<double>& xi, const vector<double>& pi,
			   const vector<double>& alpha, const vector<double>& beta,
			   const vector<double>& psi, int ind, double dr, double r) {
  return ddr2_c(alpha,ind,dr) + 2*(1/r + ddr_c(psi,ind,dr)/psi[ind])*ddr_c(alpha,ind,dr)
    - 2*pw4(psi[ind])*sq(ddr_c(beta,ind,dr) - beta[ind]/r)/(3*alpha[ind])
    - 8*M_PI*alpha[ind]*sq(pi[ind]) ; }

