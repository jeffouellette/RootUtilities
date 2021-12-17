#ifndef __LinAlg_cxx__
#define __LinAlg_cxx__

#include <LinAlg.h>


/**
 * Returns a TMatrix corresponding to the given TH2.
 */
TMatrixD GetTMatrixD (const TH2* h) {
  const int nx = h->GetNbinsX ();
  const int ny = h->GetNbinsY ();

  TMatrixD mat (nx, ny);

  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      mat[ix][iy] = h->GetBinContent (ix+1, iy+1);
    }
  }

  return mat;
}


#endif
