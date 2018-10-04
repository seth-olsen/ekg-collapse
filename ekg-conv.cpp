/*

this program will compute

q(t) = ||u_4h(x,t) - u_2h(x,t)||/||u_2h(x,t) - u_h(x,t)||

with (u_4h, u_2h, u_h) from resolution levels (0, 1, 2) 
and save to outfile (.csv) with the line format:
   step_number (= t/dt) , q_phi(t) , q_pi(t) , q_phi+pi(t)

input parameters in terminal as:
(do NOT include .sdf extension in outfile)

 ./p1-ctest <outfile> <lastpt> <save_pt> <nsteps> <save_step>
 <lam> <r2m> <rmin> <rmax> <dspn> <tol> <maxit> <ic_Dsq> 
 <ic_r0> <ic_Amp> <check_step> <zero_pi> <sommerfeld> 
 <dspn_bound> <write_ires> <write_itn> <hold_const>
 <same_times> <same_grids>

where the value of any <parameter> can be set with the
following ordered pair of command line arguments:

 -parameter parameter_value

default values can be found at the start of main()
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "fda-io.h"
using namespace std;

vector<string> get_unames(string pre, vector<string>& fnames) {
  vector<string> unames;
  for (string name : fnames) { unames.push_back(pre + "-" + name + ".sdf"); }
  return unames;
}

int main(int argc, char **argv)
{
  // coarse simulation parameters
  string outfile = "ekg";
  string outname = "0";
  string pre1 = "Xi", pre2 = "Pi", pre3 = "0", pre4 = "0", pre5 = "0",
    pre6 = "0", pre7 = "0", pre8 = "0", pre9 = "0", pre10 = "0";
  int lastpt = 1000; // grid size
  int save_pt = 1; // write only every (save_pt)th grid point
  int nsteps = 4000; // time steps
  int save_step = 8; // write only every (save_step)th time step
  double lam = 0.25; // dt/dr
  double r2m = 0;
  double rmin = 0;
  double rmax = 100.0;
  double dspn = 0.5; // dissipation coefficient
  double tol = 0.000000000001; // iterative method tolerance
  double ell_tol = tol;
  int maxit = 25; // max iterations for debugging
  double ic_Dsq = 4.0; // gaussian width
  double ic_r0 = 50.0; // gaussian center
  double ic_Amp = 1.0; // gaussian amplitude
  bool zero_pi = false; // zero initial time derivative?
  bool sommerfeld = true; // sommerfeld condition at outer bound?
  bool dspn_bound = false; // dissipate boundary points?
  bool psi_hyp = false; // psi evolved with hyperbolic eom?
  // variable to hold constant across resolutions
  string hold_const = "lambda"; // "lambda", "dt", or "dr"
  bool same_times = true;
  bool same_grids = true;
  // resolution factors
  int resn0 = 4;
  int resn1 = 2*resn0;
  int resn2 = 4*resn0; // 4h, 2h, and h

  // get parameters from command line
  map<string, string *> p_str {{"-outfile",&outfile}, {"-pre1",&pre1},
      {"-pre2",&pre2}, {"-pre3",&pre3}, {"-pre4",&pre4}, {"-pre5",&pre5},
      {"-pre6",&pre6}, {"-pre7",&pre7}, {"-pre8",&pre8}, {"-pre9",&pre9},
      {"-pre10",&pre10}, {"-hold_const",&hold_const}, {"-outname",&outname}};
  map<string, int *> p_int {{"-lastpt",&lastpt}, {"-save_pt", &save_pt},
      {"-nsteps", &nsteps}, {"-save_step",&save_step}, {"-maxit",&maxit},
      {"-resn0", &resn0}, {"-resn1", &resn1}, {"-resn2", &resn2}};
  map<string, double *> p_dbl {{"-lam",&lam}, {"-r2m",&r2m}, {"-rmin",&rmin},
      {"-rmax",&rmax}, {"-dspn",&dspn}, {"-tol",&tol}, {"-ic_Dsq",&ic_Dsq},
      {"-ic_r0",&ic_r0}, {"-ic_Amp",&ic_Amp}, {"-ell_tol",&ell_tol}};
  map<string, bool *> p_bool {
      {"-sommerfeld",&sommerfeld}, {"-dspn_bound",&dspn_bound},
      {"-same_times",&same_times}, {"-same_grids",&same_grids},
      {"-psi_hyp",&psi_hyp}};
  map<string, string> params;
  param_collect(argv, argc, params);
  param_set(params, p_str, p_int, p_dbl, p_bool);

  vector<string> prefixes{pre1};
  if (pre2 != "0") { prefixes.push_back(pre2); }
  if (pre3 != "0") { prefixes.push_back(pre3); }
  if (pre4 != "0") { prefixes.push_back(pre4); }
  if (pre5 != "0") { prefixes.push_back(pre5); }
  if (pre6 != "0") { prefixes.push_back(pre6); }
  if (pre7 != "0") { prefixes.push_back(pre7); }
  if (pre8 != "0") { prefixes.push_back(pre8); }
  if (pre9 != "0") { prefixes.push_back(pre9); }
  if (pre10 != "0") { prefixes.push_back(pre10); }
  int nwr = prefixes.size();

  // derived parameters from coarse file
  int gs = lastpt / save_pt;
  int num_steps = nsteps / save_step;
  double dr = (rmax - rmin) / ((double) lastpt);
  double dt = lam * dr;
  int npts0 = gs + 1;
  int npts1 = (same_grids) ? npts0 : 2*gs + 1;
  int npts2 = (same_grids) ? npts0 : 4*gs + 1;
  
  vector< vector<double> > norms(nwr, vector<double>(2, 0.0));
  vector<double> zeros(2, 0.0);
  vector< vector<double> > u_4h(nwr, vector<double>(npts0, 0.0)),
    u_2h(nwr, vector<double>(npts1, 0.0)), u_h(nwr, vector<double>(npts2, 0.0));
  vector< vector<double *> > field_arr;

  vector<string> fnames{to_string(resn0) + "-" + outfile,
	to_string(resn1) + "-" + outfile, to_string(resn2) + "-" + outfile};
  vector< vector<string> > unames;
  vector< vector<char *> > name_arr;
  // output file
  ofstream ofs;
  if (outname == "0") { outname = "conv-" + fnames[0] + ".csv"; }
  ofs.open(outname, ofstream::out);
  ofs <<  "coarse,mid,fine,constant,points,times,same_times,same_grids\n"
      <<  fnames[0]+","+fnames[1]+","+fnames[2]+","+hold_const+"," << npts0
      <<","<< num_steps <<","<< boolalpha << same_times <<","<< same_grids
      << "\n\ndspn,dspn_bound,zero_pi,sommerfeld,tol,maxit\n" << dspn <<","
      << dspn_bound <<","<< zero_pi <<","<< sommerfeld <<","<< tol <<","
      << maxit << "\n\nr2m,rmin,rmax,ic_Dsq,ic_r0,ic_Amp\n" << r2m <<","
      << rmin <<","<< rmax <<","<< ic_Dsq <<","<< ic_r0 <<","<< ic_Amp
      << "\n\ncoarse grid:\nlastpt,save_pt,nsteps,save_step,\n" << lastpt
      <<","<< save_pt <<","<< nsteps <<","<< save_step << "\n\nlam,dr,dt,"
      << "tmax\n" << lam <<","<< dr <<","<< dt <<","<< dt*nsteps
      << "\n\ntime";
  for (int k = 0; k < nwr; ++k) {
    ofs << ",Q" << prefixes[k] << "(t)";
    unames.push_back(get_unames(prefixes[k], fnames));
    name_arr.push_back({&unames[k][0][0], &unames[k][1][0], &unames[k][2][0]});
    field_arr.push_back({&u_4h[k][0], &u_2h[k][0], &u_h[k][0]});
  }
  ofs << endl;  
  
  // iterate through time steps
  int t1 = ((same_times) ? 1 : 2);
  int t2 = ((same_times) ? 1 : 4);
  int r1 = ((same_grids) ? 1 : 2);
  int r2 = ((same_grids) ? 1 : 4);
  int times[3];

  gft_set_multi();
  for (int t = 0; t < num_steps; ++t) {
    // compute ||u_4h(x,t) - u_2h(x,t)||/||u_2h(x,t) - u_h(x,t)||    
    times[0] = t+1, times[1] = t1*t+1, times[2] = t2*t+1;
    for (int k = 0; k < nwr; ++k) {
      read_step(name_arr[k], times, field_arr[k], 3);
      norms[k] = zeros;
    }
    // iterate through grid points
    for (int j = 0; j < npts0; ++j) {
      for (int k = 0; k < nwr; ++k) {
	norms[k][0] += abs(u_4h[k][j] - u_2h[k][r1*j]);
	norms[k][1] += abs(u_2h[k][r1*j] - u_h[k][r2*j]);
      }
    }
    // write time, q_u,
    ofs << t*save_step*dt;
    for (int k = 0; k < nwr; ++k) { ofs <<","<< norms[k][0]/norms[k][1]; }
    ofs << endl;
  }
  gft_close_all();
  
  cout << outname << "  written with:" << endl;
  cout << "grid points used = " << npts0 << "  " << ((same_grids) ?
   "(same grids)" : "(dif grids)") << "\ntime steps used = " << num_steps
       << "  " << ((same_times) ? "(same times)" : "(dif times)") << endl;
  ofs.close();

  return 0;
}
