#include "ekg-debug-proc.h"

int main(int argc, char**argv)
{
  int info = 0;
  int n = 6, ku = 2, kl = 2, nrhs = 1;
  int ldab = 2*kl + ku + 1;
  int ldb = n;
  vector<int> ipiv(n);
  vector<double> mat_vec(ldb*ldab), res(ldb);
  vector<int> vals{4,5,6,10,11,12,16,17,18,19,23,24,25,26,27,30,31,32,33,37,38,39};
  vector<int> zvals{13,20};
  for (int num : vals) { mat_vec[num] = num; }
  for (int num : zvals) { mat_vec[num] = 0; }
  for (int num = 0; num < ldb; ++num) { res[num] = num; }
  info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, n, kl, ku, nrhs, &mat_vec[0], ldab, &ipiv[0], &res[0], ldb);
  cout << "info = " << info << "\nres = " << endl;
  for (int num = 0; num < ldb; ++num) { cout << res[num] << endl; }


  return 0;
}
  
  




