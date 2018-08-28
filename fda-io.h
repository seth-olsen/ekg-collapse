#include <string> // for parameter input
#include <map> // for parameter input
#include <cstdlib> // for atoi() and atof()
#include "bbhutil.h" // for output to .sdf

using namespace std;

void param_collect(char **source, int num, map<string, string>& dest) {
  for (int arg = 1; arg < num; ++arg) {
    if (source[arg][0] == '-') {
      dest[source[arg]] = source[arg+1];
    }
  }
}

void param_set(map<string, string>& p_all, map<string, string *>& p_str,
	       map<string, int *>& p_int, map<string, double *>& p_dbl,
	       map<string, bool *>& p_bool) {
  for (pair<string, string> p : p_all) {
    if (p_str.count(p.first)) { *p_str[p.first] = p.second; }
    else if (p_int.count(p.first)) { *p_int[p.first] = atoi(&p.second[0]); }
    else if (p_dbl.count(p.first)) { *p_dbl[p.first] = atof(&p.second[0]); }
    else if (p_bool.count(p.first)) { *p_bool[p.first] = (bool) atoi(&p.second[0]); }
  }
}

// write fields using bbhutil
void wr_step(double* fields[], int nfields, char* files[],
		double time, int* shape, int rank, double* coordinates) {
  for (int k = 0; k < nfields; ++k) {
    gft_out_bbox(files[k], time, shape, rank, coordinates, fields[k]);
  }
  return;
}


// read fields using bbhutil
void read_step(char *files[], int times[], double *fields[], int nfields) {
  for (int k = 0; k < nfields; ++k) {
    gft_read_brief(files[k], times[k], fields[k]);
  }
  return;
}

string param_print(string outfile, int lastpt, int save_pt = 1, int nsteps,
		   int save_step = 8, double lam = 0.25, double r2m,
		   double rmin, double rmax, double dspn, double tol, int maxit,
		   double ic_Dsq, double ic_r0, double ic_Amp, int check_step,
		   bool zero_pi, bool somm_cond = true, bool dspn_bound) {
  string z_p = (zero_pi) ? "true" : "false";
  string s_c = (somm_cond) ? "true" : "false";
  string d_b = (dspn_bound) ? "true" : "false";
  return "\noutfile name = " + outfile + "\ngrid size = " +
    to_string(lastpt) + " (" + to_string(save_pt) + "/write)\ntime steps = " +
    to_string(nsteps) + " (" + to_string(save_step) + "/write)\nlambda = " +
    to_string(lam) + "\nr2m = " + to_string(r2m) + "\nrmin = " + to_string(rmin)
    + "\nrmax = " + to_string(rmax) + "\ndissipation = " + to_string(dspn) +
    "\niterative tolerance = " + to_string(tol) + "\nmaximum iterations = " +
    to_string(maxit) + "\nic_Dsq = " + to_string(ic_Dsq) + "\nic_r0 = " +
    to_string(ic_r0) + "\nic_Amp = " + to_string(ic_Amp) + "\nmass check step = "
    + to_string(check_step) + "\nmaximum evolution time = " + to_string(nsteps*dt)
    + "\ndr = " + to_string(dr) + "\ndt = " + to_string(dt) +
    "\n\noptions:\nzero pi_0 = " + z_p + "\nsommerfeld bc = " + s_c +
    "\ndissipation at bound = " + d_b + "\n";
}

