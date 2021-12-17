#ifndef __LinAlg_h__
#define __LinAlg_h__

#include <TH2.h>
#include <TMatrixD.h>


/**
 * Returns a TMatrix corresponding to the given TH2.
 */
TMatrixD GetTMatrixD (const TH2* h);


#endif
