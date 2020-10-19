#ifndef __Utilities_cxx__
#define __Utilities_cxx__

#include "Utilities.h"

#include <TSystemDirectory.h>
#include <TList.h>
#include <TSystemFile.h>
#include <TF1.h>
#include <TROOT.h>
#include <TLine.h>
#include <TLatex.h>
#include <TMarker.h>
#include <TPave.h>
#include <TH1.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TPad.h>

#include <math.h>

#include <iostream>
#include <cmath>

using namespace std;

/**
 * Returns a TString summarizing a measurement.
 * By default, there is only 1 significant digit in the error.
 * E.g., FormatMeasurement (40.58, 1.29) returns "40#pm1".
 * Or, FormatMeasurement (40.58, 1.29, 2) returns "40.6#pm1.3".
 */
const char* FormatMeasurement (double val, double err, const int n) {
  assert (n > 0); // sanity check

  TString out = "";

  string valStr = Form ("%g", val);
  string errStr = Form ("%g", err);

  if (err == 0) {
    short valDec = 0;
    while (valDec < (short)valStr.length() && valStr[valDec] != '.')
      valDec++;
    if (valDec == (short)valStr.length()) {
      valStr = valStr + ".";
    }

    if (valDec > n) { // if decimal is after least significant digit
      const double factorOfTen = pow (10, n-valDec); // e.g., for "1520" and n=2, get 0.01
      val = floor (factorOfTen * val + 0.5) / factorOfTen;
    }
    else { // if decimal is before least significant digit, e.g. for "15.24" and n=3.
      const double factorOfTen = pow (10, n-valDec+1); // e.g. for 15.24 and n=3 get 10
      val = floor (factorOfTen * val + 0.5) / factorOfTen;
    }
    return Form ("%g", val);
  }

  if (err < 1) {
    // find the first significant digit
    int errStart = 0;
    while (errStr[errStart] == '0' || errStr[errStart] == '.')
      errStart++;
    errStart = errStart - 2; // first 2 characters are necessarly "0."

    // round the value and error to the appropriate decimal place
    const double factorOfTen = pow (10, errStart+n);
    val = floor (factorOfTen * val + 0.5) / factorOfTen;
    err = floor (factorOfTen * err + 0.5) / factorOfTen;

    // recast to string
    valStr = Form ("%g", val);
    errStr = Form ("%g", err);

    // find where the decimal place is, append it if it's not present
    short valDec = 0;
    while (valDec < (short)valStr.length() && valStr[valDec] != '.')
      valDec++;
    if (valDec == (short)valStr.length()) {
      valStr = valStr + ".";
    }

    // pad with zeroes
    while ((short)valStr.length() < valDec + errStart + 1 + n)
      valStr = valStr + "0";
    while ((short)errStr.length() < errStart + 2 + n)
      errStr = errStr + "0";

    //// find where we are truncating the value
    //const short valCut = valDec + errStart + 2 + n;

    //// now truncate
    //valStr = valStr.substr (0, valCut);
    //errStr = errStr.substr (0, errStart + 2 + n);
  }
  else { // now if err>1
    // find the decimal place
    short errDec = 0;
    while (errDec < (short)errStr.length() && errStr[errDec] != '.')
      errDec++;

    // round the value and error to the appropriate decimal place
    if (errDec < (short)errStr.length() - 1) {
      const double factorOfTen = pow (10., n - errDec);
      if (n - errDec >= 0) {
        val = floor (factorOfTen * val + 0.5) / factorOfTen;
        err = floor (factorOfTen * err + 0.5) / factorOfTen;
      } else {
        val = floor (factorOfTen * val) / factorOfTen;
        err = floor (factorOfTen * err) / factorOfTen;
      }
    }

    // recast to string
    valStr = Form ("%g", val);
    errStr = Form ("%g", err);
  }

  // now save the value string to the output
  out = valStr + " #pm " + errStr + "";

  return out.Data();
}


/**
 * Truncates numbers to the desired precision based on the most precise uncertainty (stat. or syst.)
 */
void FormatMeasurement (string& s_val, string& s_stat, string& s_syst, const int n) {
  assert (n > 0); // sanity check

  double val = atof (s_val.c_str ());
  double stat = atof (s_stat.c_str ());
  double syst = atof (s_syst.c_str ());

  int i = 0;
  while (fabs (stat) < pow (10, n-1) || fabs (syst) < pow (10, n-1) && i < 10) { // i < 10 as an arbitrary cut (so this doesn't go forever)
    val *= 10;
    stat *= 10;
    syst *= 10;
    i++;
  }

  val = round (val);
  stat = round (stat);
  syst = round (syst);

  val = val / pow (10, i);
  stat = stat / pow (10, i);
  syst = syst / pow (10, i);

  s_val = to_string (val);
  s_stat = to_string (stat);
  s_syst = to_string (syst);

  while (s_val.back () == '0' && s_stat.back () == '0' && s_syst.back () == '0') {
    s_val.pop_back ();
    s_stat.pop_back ();
    s_syst.pop_back ();
  }

  return;
}


/**
 * Returns a linearly spaced array. The 0th element is lo, and the num-th element is hi.
 */
double* linspace (double lo, double hi, int num) {
  double* arr = new double[num+1];
  double delta = ((double)(hi)-(double)(lo))/(double)(num);
  for (int i = 0; i <= num; i++) {
    arr[i] = lo + i * delta;
  }
  return arr;
}


/**
 * Returns a logarithmically spaced array, where the 0th element is lo and the num-th element is hi.
 */
double* logspace (double lo, double hi, int num) {
  double loghi = TMath::Log2(hi);
  if (lo == 0) {
    double* arr = linspace(TMath::Log2(hi/(100*num)), loghi, num);
    for (int i = 0; i <= num; i++) {
      arr[i] = TMath::Power(2, arr[i]);
    }
    return arr;
  } else {
    double loglo = TMath::Log2(lo);
    double* arr = linspace(loglo, loghi, num);
    for (int i = 0; i <= num; i++) {
      arr[i] = TMath::Power(2, arr[i]);
    }
    return arr;
  }
}


/**
 * Returns the equivalent angle in the range 0 to 2pi.
 */
double InTwoPi (double phi) {
  while (phi < 0 || 2*M_PI <= phi) {
   if (phi < 0) phi += 2*M_PI;
   else phi -= 2*M_PI;
  }
  return phi;
}


/**
 * Returns the difference between two angles in 0 to pi.
 * If sign is true, will return a signed version such that phi2 = phi1 + dphi
 */
double DeltaPhi (double phi1, double phi2, const bool sign) {
  phi1 = InTwoPi(phi1);
  phi2 = InTwoPi(phi2);
  double dphi = abs(phi1 - phi2);
  while (dphi > M_PI) dphi = abs (dphi - 2*M_PI);

  if (sign && InTwoPi (phi2 + dphi) == phi1)
     dphi *= -1;

  return dphi;
}


/**
 * Returns dR between two eta, phi coordinates.
 */
double DeltaR (const double eta1, const double eta2, const double phi1, const double phi2 ) {
 const double deta = eta1 - eta2;
 const double dphi = DeltaPhi (phi1, phi2, false);
 return sqrt( pow( deta, 2 ) + pow( dphi, 2 ) );
}


/**
 * Sets all the errors in this histogram to 0.
 */
void ResetHistErrors (TH1D* h) {
  for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
    h->SetBinError (ix, 0);
  }
}


/**
 * Sets all the errors in this TGAE to 0.
 */
void ResetTGAEErrors (TGAE* g) {
  for (int ix = 0; ix < g->GetN (); ix++) {
    g->SetPointEYhigh (ix, 0);
    g->SetPointEYlow (ix, 0);
  }
}


/**
 * Sets all the x errors in this TGAE to 0.
 */
void ResetXErrors (TGAE* tg) {
  for (int ix = 0; ix < tg->GetN (); ix++) {
    tg->SetPointEXlow (ix, 0);
    tg->SetPointEXhigh (ix, 0);
  }
  return;
}


/**
 * Adds nSigma statistical error variations to this histogram
 */
void AddStatVar (TH1D* h, const bool upvar, const float nSigma) {
  if (!h)
    return;
  for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
    h->SetBinContent (ix, h->GetBinContent (ix) + (upvar ? 1.:-1.) * nSigma * h->GetBinError (ix));
  }
  return;
}



/**
 * Adds independent systematic errors in quadrature, storing the sum in master
 */
void AddErrorsInQuadrature (TH1D* master, TH1D* sys) {
  for (int ix = 1; ix <= master->GetNbinsX (); ix++) {
    if (master->GetBinContent (ix) != sys->GetBinContent (ix))
      cout << "Warning: In Utilities.cxx::AddErrorsInQuadrature: unequal nominal bin contents between systematics! Results will be skewed." << endl;

    const float newErr = sqrt (pow (master->GetBinError (ix), 2) + pow (sys->GetBinError (ix), 2));
    master->SetBinError (ix, newErr);
  }
}


/**
 * Adds independent systematic errors in quadrature, storing the sum in master
 */
void AddErrorsInQuadrature (TGAE* master, TGAE* sys, const bool doXErrs) {
  for (int ix = 0; ix < master->GetN (); ix++) {
    //double xm=0, ym=0, xs=0, ys=0;
    //master->GetPoint (ix, xm, ym);
    //sys->GetPoint (ix, xs, ys);

    //if (xm != xs || ym != ys)
    //  cout << "Warning: In Utilities.cxx::AddErrorsInQuadrature: Bins don't match!" << endl;

    float newErr = sqrt (pow (master->GetErrorYhigh (ix), 2) + pow (sys->GetErrorYhigh (ix), 2));
    master->SetPointEYhigh (ix, newErr);
    newErr = sqrt (pow (master->GetErrorYlow (ix), 2) + pow (sys->GetErrorYlow (ix), 2));
    master->SetPointEYlow (ix, newErr);

    if (doXErrs) {
      newErr = sqrt (pow (master->GetErrorXhigh (ix), 2) + pow (sys->GetErrorXhigh (ix), 2));
      master->SetPointEXhigh (ix, newErr);
      newErr = sqrt (pow (master->GetErrorXlow (ix), 2) + pow (sys->GetErrorXlow (ix), 2));
      master->SetPointEXlow (ix, newErr);
    }
  }
}


/**
 * Calculates simple systematics as maximum variations on the nominal.
 * Intended for combining up/down variations in an expandable way.
 */
void CalcSystematics (TH1D* sys, TH1D* var) {
  for (int ix = 1; ix <= sys->GetNbinsX (); ix++) {
    const float newErr = fabs (var->GetBinContent (ix) - sys->GetBinContent (ix));
    sys->SetBinError (ix, fmax (newErr, sys->GetBinError (ix)));
  }
}


/**
 * Calculates simple systematics as maximum variations on the nominal.
 * Intended for combining up/down variations in an expandable way.
 * dir represents the variation direction, -1 = down, 0 = either, 1 = up
 */
void CalcSystematics (TGAE* sys, TH1D* var, const bool applyBothWays) {
  for (int ix = 0; ix < sys->GetN (); ix++) {
    double x, y;
    sys->GetPoint (ix, x, y);
    const float newErr = var->GetBinContent (ix+1) - y;

    if (applyBothWays || newErr > 0) {
      sys->SetPointEYhigh (ix, fmax (fabs (newErr), sys->GetErrorYhigh (ix)));
    }
    if (applyBothWays || newErr < 0) {
      sys->SetPointEYlow (ix, fmax (fabs (newErr), sys->GetErrorYlow (ix)));
    }
  }
}


void CalcSystematics (TGAE* graph, TH1* optimal, const TH1* sys_hi, const TH1* sys_lo) {
  for (int ix = 1; ix <= optimal->GetNbinsX(); ix++) {
    const double content = optimal->GetBinContent (ix);
    const double lo = sys_lo->GetBinContent (ix);
    const double hi = sys_hi->GetBinContent (ix);

    const double err_lo = content - lo;
    const double err_hi = hi - content;

  

    graph->SetPoint (ix-1, optimal->GetBinCenter (ix), content);
    graph->SetPointEXlow (ix-1, optimal->GetBinCenter (ix) - optimal->GetBinLowEdge (ix));
    graph->SetPointEXhigh (ix-1, optimal->GetBinLowEdge (ix+1) - optimal->GetBinCenter (ix));

    if (err_lo < 0 && err_hi < 0) {
      graph->SetPointEYlow (ix-1, -err_hi);
      graph->SetPointEYhigh (ix-1, -err_lo);
    }
    else if (err_lo >= 0 && err_hi >= 0) {
      graph->SetPointEYlow (ix-1, err_lo);
      graph->SetPointEYhigh (ix-1, err_hi);
    }
    else {
      graph->SetPointEYlow (ix-1, sqrt (0.5*(err_lo*err_lo + err_hi*err_hi)));
      graph->SetPointEYhigh (ix-1, sqrt (0.5*(err_lo*err_lo + err_hi*err_hi)));
    }

  }
  return;
}


/**
 * Calculates the systematic errors on optimal, storing the results in graph.
 */
void CalcSystematics (TGAE* graph, TGAE* optimal, const TGraph* sys_hi, const TGraph* sys_lo, const bool doXErrs) {
  for (int ix = 0; ix < optimal->GetN(); ix++) {
    double x, y, xl, yl, xh, yh;
    optimal->GetPoint (ix, x, y);
    sys_lo->GetPoint (ix, xl, yl);
    sys_hi->GetPoint (ix, xh, yh);

    const double err_lo = y - yl;
    const double err_hi = yh - y;

    graph->SetPoint (ix, x, y);

    if (!doXErrs) {
      graph->SetPointEXlow (ix, optimal->GetErrorXlow (ix));
      graph->SetPointEXhigh (ix, optimal->GetErrorXhigh (ix));
    }
    else {
      const double xerr = fmax (fabs (xh - x), fabs (x - xl));
      graph->SetPointEXlow (ix, xerr);
      graph->SetPointEXhigh (ix, xerr);
    }

    if (err_lo < 0 && err_hi < 0) {
      graph->SetPointEYlow (ix, -err_hi);
      graph->SetPointEYhigh (ix, -err_lo);
    }
    else if (err_lo >= 0 && err_hi >= 0) {
      graph->SetPointEYlow (ix, err_lo);
      graph->SetPointEYhigh (ix, err_hi);
    }
    else {
      graph->SetPointEYlow (ix, sqrt (0.5*(err_lo*err_lo + err_hi*err_hi)));
      graph->SetPointEYhigh (ix, sqrt (0.5*(err_lo*err_lo + err_hi*err_hi)));
    }
  }
  return;
}


/**
 * Sets the bin contents in target as the error in errors / central values in centralValues
 */
void SaveRelativeErrors (TH1D* errors, TH1D* centralValues) {
  assert (errors->GetNbinsX () == centralValues->GetNbinsX ());
  for (int ix = 1; ix <= errors->GetNbinsX (); ix++) {
    if (centralValues->GetBinContent (ix) != 0)
      errors->SetBinContent (ix, errors->GetBinError (ix) / centralValues->GetBinContent (ix));
    errors->SetBinError (ix, 0);
  }
}


/**
 * Sets the bin contents in highs and lows as the respective errors in centralValues
 */
void SaveAbsoluteErrors (TGAE* errors, TGAE* centralValues, TH1D* highs, TH1D* lows) {
  assert (errors->GetN () == centralValues->GetN ());
  for (int ix = 0; ix < centralValues->GetN (); ix++) {
    double x = 0, y = 0, eyhi = 0, eylo = 0;
    eyhi = errors->GetErrorYhigh (ix);
    eylo = errors->GetErrorYlow (ix);

    highs->SetBinContent (ix+1, eyhi);
    lows->SetBinContent (ix+1, -eylo);

    highs->SetBinError (ix+1, 0);
    lows->SetBinError (ix+1, 0);
  }
}


/**
 * Sets the bin contents in highs and lows as the respective errors / central values in centralValues
 */
void SaveRelativeErrors (TGAE* errors, TGAE* centralValues, TH1D* highs, TH1D* lows) {
  assert (errors->GetN () == centralValues->GetN ());
  for (int ix = 0; ix < centralValues->GetN (); ix++) {
    double x = 0, y = 0, eyhi = 0, eylo = 0;
    centralValues->GetPoint (ix, x, y);
    eyhi = errors->GetErrorYhigh (ix);
    eylo = errors->GetErrorYlow (ix);

    if (y != 0) {
      highs->SetBinContent (ix+1, eyhi / y);
      lows->SetBinContent (ix+1, -eylo / y);
    }
    highs->SetBinError (ix+1, 0);
    lows->SetBinError (ix+1, 0);
  }
}



/**
 * Creates a projection of a TH3 with the specified axes by integrating between min and max on the 3rd axis of the TH3.
 */
TH2D* Project2D (TString name, TH3D* h3, const TString xaxis, const TString yaxis, const int min, const int max, const bool exclusive) {
  int nx, ny, nz;
  double* xbins = NULL;
  double* ybins = NULL;

  if (xaxis == yaxis || (!(xaxis == "x" || xaxis == "y" || xaxis == "z") || !(yaxis == "x" || yaxis == "y" || yaxis == "z")))
   return NULL;

  // determine which axis in the TH3 is the xaxis in the TH2
  if      (xaxis == "x") { nx = h3->GetNbinsX(); xbins = (double*)h3->GetXaxis()->GetXbins()->GetArray(); }
  else if (xaxis == "y") { nx = h3->GetNbinsY(); xbins = (double*)h3->GetYaxis()->GetXbins()->GetArray(); }
  else if (xaxis == "z") { nx = h3->GetNbinsZ(); xbins = (double*)h3->GetZaxis()->GetXbins()->GetArray(); }

  // determine which axis in the TH3 is the yaxis in the TH2
  if      (yaxis == "x") { ny = h3->GetNbinsX(); ybins = (double*)h3->GetXaxis()->GetXbins()->GetArray(); }
  else if (yaxis == "y") { ny = h3->GetNbinsY(); ybins = (double*)h3->GetYaxis()->GetXbins()->GetArray(); }
  else if (yaxis == "z") { ny = h3->GetNbinsZ(); ybins = (double*)h3->GetZaxis()->GetXbins()->GetArray(); }

  if (!xbins || !ybins) return NULL;

  // determine which axis in TH3 is the extra "z" axis
  if      (xaxis != "x" && yaxis != "x") nz = h3->GetNbinsX();
  else if (xaxis != "y" && yaxis != "y") nz = h3->GetNbinsY();
  else if (xaxis != "z" && yaxis != "z") nz = h3->GetNbinsZ();

  if (name == "") name = TString(h3->GetName()) + "_" + xaxis + yaxis;
  TH2D* h2 = new TH2D (name, "", nx, xbins, ny, ybins);

  for (int ix = 1; ix <= nx; ix++) {
   for (int iy = 1; iy <= ny; iy++) {
    double content = 0;
    double var = 0;

    int h3ix, h3iy, h3iz;
    for (int iz = 1; iz <= nz; iz++) { // loop over extra "z" axis (not in TH2)
     if (!exclusive && (iz < min || max < iz))
      continue;
     else if (exclusive && min <= iz && iz <= max)
      continue;

     if      (xaxis == "x" && yaxis == "y") { h3ix = ix; h3iy = iy; h3iz = iz; }
     else if (xaxis == "x" && yaxis == "z") { h3ix = ix; h3iz = iy; h3iy = iz; }
     else if (xaxis == "y" && yaxis == "x") { h3iy = ix; h3ix = iy; h3iz = iz; }
     else if (xaxis == "y" && yaxis == "z") { h3iy = ix; h3iz = iy; h3ix = iz; }
     else if (xaxis == "z" && yaxis == "x") { h3iz = ix; h3ix = iy; h3iy = iz; }
     else if (xaxis == "z" && yaxis == "y") { h3iz = ix; h3iy = iy; h3ix = iz; }

     content += h3->GetBinContent (h3ix, h3iy, h3iz);
     var += pow (h3->GetBinError (h3ix, h3iy, h3iz), 2);

    } // end loop over extra axis

    h2->SetBinContent (ix, iy, content);
    h2->SetBinError (ix, iy, sqrt (var));
   }
  }

  if      (xaxis == "x") h2->GetXaxis()->SetTitle (h3->GetXaxis()->GetTitle());
  else if (xaxis == "y") h2->GetXaxis()->SetTitle (h3->GetYaxis()->GetTitle());
  else if (xaxis == "z") h2->GetXaxis()->SetTitle (h3->GetZaxis()->GetTitle());

  if      (yaxis == "x") h2->GetYaxis()->SetTitle (h3->GetXaxis()->GetTitle());
  else if (yaxis == "y") h2->GetYaxis()->SetTitle (h3->GetYaxis()->GetTitle());
  else if (yaxis == "z") h2->GetYaxis()->SetTitle (h3->GetZaxis()->GetTitle());

  return h2;
}


/**
 * Separates each point on a TGAE by delta along the x axis, so that the errors don't overlap.
 */
void deltaize (TGAE* tg, const double delta, const bool logx, const double x1, const double x2) {
  double x, y, exh, exl;

  if (x1*x2 < 0) {
    cout << "Error in deltaize: x1 and x2 have different signs! Assuming both are negative." << endl;
  }
  const bool useX2OverX1Factor = (x1 > 0 && x2 > 0);

  for (int n = 0; n < tg->GetN (); n++) {
    tg->GetPoint (n, x, y);
    exh = x + tg->GetErrorXhigh (n);
    exl = x - tg->GetErrorXlow (n);

    if (logx && useX2OverX1Factor) tg->SetPoint (n, x * exp (0.5 * delta * log (x2/x1)), y);
    else if (logx) tg->SetPoint (n, x*delta, y);
    else tg->SetPoint (n, x + delta, y);

    tg->GetPoint (n, x, y);
    exh = fabs (exh - x);
    exl = fabs (x - exl);

    tg->SetPointEXhigh (n, exh);
    tg->SetPointEXlow (n, exl);
  }
}


/**
 * Offsets each point on a TGAE by delta along the y axis.
 */
void OffsetYAxis (TGAE* g, const double delta, const bool logx) {
  if (logx && delta < 0) {
    cout << "delta < 0 and logx set, these are incompatible options! Returning." << endl;
    return;
  }
  double x, y, yhe, yle;
  for (int n = 0; n < g->GetN (); n++) {
    g->GetPoint (n, x, y);
    yhe = g->GetErrorYhigh (n);
    yle = g->GetErrorYlow (n);

    if (logx) {
      g->SetPoint (n, x, y*delta);
      g->SetPointEYhigh (n, yhe*delta);
      g->SetPointEYlow (n, yle*delta);
    }
    else {
      g->SetPoint (n, x, y*(1+delta));
    }
  }
  return;
}


/**
 * Sets each point error to some constant value.
 * If mult is true, then will be multiplicative (intended for a log scale).
 */
void SetConstantXErrors (TGAE* tg, const double err, const bool mult, const double x1, const double x2) {
  double x, y, exh, exl;

  if (x1*x2 < 0) {
    cout << "Error in deltaize: x1 and x2 have different signs! Assuming both are negative." << endl;
  }
  const bool useX2OverX1Factor = (x1 > 0 && x2 > 0);

  for (int n = 0; n < tg->GetN (); n++) {
    tg->GetPoint (n, x, y);
    if (mult && useX2OverX1Factor) {
      exh = fabs (x * (exp (0.5*err * log (x2/x1)) - 1));
      exl = fabs (x * (1 - exp (-0.5*err * log(x2/x1))));
    }
    else if (mult) {
      exh = fabs (x * (exp (0.5*err) - 1));
      exl = fabs (x * (1 - exp (-0.5*err)));
    }
    else {
      exh = err;
      exl = err;
    }
    tg->SetPointEXhigh (n, exh);
    tg->SetPointEXlow (n, exl);
  }
  return;
}


/**
 * Makes a TGAE from the input histogram.
 */
TGAE* make_graph (TH1* h, const float cutoff) {
  TGAE* tg = new TGAE ();

  const float xlo = h->GetBinLowEdge (1);
  const float xhi = h->GetBinLowEdge (h->GetNbinsX ()) + h->GetBinWidth (h->GetNbinsX ());

  for (int n = 0; n < h->GetNbinsX (); n++) {
    if (cutoff >= 0 && h->GetBinContent (n+1) <= cutoff) {
      continue;
    }
    else {
      tg->SetPoint (tg->GetN (), h->GetBinCenter (n+1), h->GetBinContent (n+1));
      tg->SetPointError (tg->GetN ()-1, h->GetBinWidth (n+1) / 2, h->GetBinWidth (n+1) / 2, h->GetBinError (n+1), h->GetBinError (n+1));
    }
  }

  tg->GetXaxis ()->SetRangeUser (xlo, xhi);
  return tg;
}


/**
 * Converts a TEfficiency to a TGAE
 */
TGAE* TEff2TGAE (TEfficiency* e) {
  const int nx = e->GetTotalHistogram ()->GetNbinsX ();
  TGAE* g = new TGAE (nx);
  for (int ix = 0; ix < nx; ix++) {
    g->SetPoint (ix, e->GetTotalHistogram ()->GetBinCenter (ix+1), e->GetEfficiency (ix+1));
    g->SetPointEXhigh (ix, 0.5*e->GetTotalHistogram ()->GetBinWidth (ix+1));
    g->SetPointEXlow (ix, 0.5*e->GetTotalHistogram ()->GetBinWidth (ix+1));
    g->SetPointEYhigh (ix, e->GetEfficiencyErrorUp (ix+1));
    g->SetPointEYlow (ix, e->GetEfficiencyErrorLow (ix+1));
  }
  return g;
}


/**
 * Recenters a TGraph according to the midpoints of matched
 */
void RecenterGraph (TGraph* g, TGraph* matched) {
  double x, y, dummy;
  assert (g->GetN () == matched->GetN ());
  for (int n = 0; n < g->GetN (); n++) {
    g->GetPoint (n, x, y);
    matched->GetPoint (n, x, dummy);
    g->SetPoint (n, x, y);
  }
  return;
}


/**
 * Recenters a TGAE point for a log scale.
 */
void RecenterGraph (TGAE* g) {
  double x, y, newx, xlo, xhi, logDelta;
  for (int n = 0; n < g->GetN (); n++) {
    g->GetPoint (n, x, y);
    xlo = x - g->GetErrorXlow (n);
    xhi = x + g->GetErrorXhigh (n);
    logDelta = log10 (xhi) - log10 (xlo);
    newx = pow (10, log10 (xlo) + 0.5*logDelta);
    g->SetPoint (n, newx, y);
    g->SetPointEXlow (n, newx-xlo);
    g->SetPointEXhigh (n, xhi-newx);
  }
}


/**
 * Applies new binning to a histogram
 * BE CAREFUL: if bins edges don't overlap, this can lead to unexpected behavior!
 */ 
void RebinSomeBins (TH1D** _h, int nbins, double* bins) {
  TH1D* h = (*_h);

  if (nbins >= h->GetNbinsX ()) {
    cout << "More new bins than old bins, returning." << endl;
    return;
  }

  const TString name = h->GetName ();
  h->SetName ("temp");

  int noldbins = h->GetNbinsX ();
  double* oldbins = new double[noldbins+1];
  for (int ix = 1; ix <= noldbins; ix++) {
    oldbins[ix-1] = h->GetBinLowEdge (ix);
  }
  oldbins[noldbins] = h->GetBinLowEdge (noldbins) + h->GetBinWidth (noldbins);

  TH1D* hnew = new TH1D (name, "", nbins, bins);
  hnew->Sumw2 ();
  for (int ix = 0; ix < nbins; ix++) {
    for (int ixprime = 0; ixprime < noldbins; ixprime++) {
      if (bins[ix] <= oldbins[ixprime] && oldbins[ixprime+1] <= bins[ix+1]) {
        hnew->SetBinContent (ix+1, hnew->GetBinContent (ix+1) + h->GetBinContent (ixprime+1));
        hnew->SetBinError (ix+1, hnew->GetBinError (ix+1) + pow (h->GetBinError (ixprime+1), 2));
      }
      else if (oldbins[ixprime] <= bins[ix] && bins[ix+1] <= oldbins[ixprime+1]) {
        hnew->SetBinContent (ix+1, hnew->GetBinContent (ix+1) + h->GetBinContent (ixprime+1));
        hnew->SetBinError (ix+1, hnew->GetBinError (ix+1) + pow (h->GetBinError (ixprime+1), 2));
      }
    }
    hnew->SetBinError (ix+1, sqrt (hnew->GetBinError (ix+1)));
  }
  delete[] oldbins;

  delete h;
  *_h = hnew;
}


/**
 * Adds a to h without propagating errors (e.g. for subtracting a background)
 */
void AddNoErrors (TH1D* h, TH1D* a, const float sf) {
  assert (h->GetNbinsX () == a->GetNbinsX ());
  for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
    h->SetBinContent (ix, h->GetBinContent (ix) + sf * a->GetBinContent (ix));
  }
}


/**
 * Multiplies h by a raised to power without propagating errors (e.g. for subtracting a background)
 */
void MultiplyNoErrors (TH1D* h, TH1D* a, const float power) {
  for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
    h->SetBinContent (ix, h->GetBinContent (ix) * pow (a->GetBinContent (ix), power));
  }
}



/**
 * Returns true if marker is "open"
 */
bool IsOpenMarker (const Style_t ms) {
  return ms == kOpenCircle ||
         ms == kOpenSquare ||
         ms == kOpenTriangleUp ||
         ms == kOpenDiamond ||
         ms == kOpenCross ||
         ms == kOpenStar ||
         ms == kOpenTriangleDown ||
         ms == kOpenDiamondCross ||
         ms == kOpenSquareDiagonal ||
         ms == kOpenThreeTriangles ||
         ms == kOpenFourTrianglesX ||
         ms == kOpenDoubleDiamond ||
         ms == kOpenFourTrianglesPlus ||
         ms == kOpenCrossX;
}


/**
 * Returns true if marker is "full"
 */
bool IsFullMarker (const Style_t ms) {
  return ms == kFullCircle ||
         ms == kFullSquare ||
         ms == kFullTriangleUp ||
         ms == kFullDiamond ||
         ms == kFullCross ||
         ms == kFullStar ||
         ms == kFullTriangleDown ||
         ms == kFullThreeTriangles ||
         ms == kFullFourTrianglesX ||
         ms == kFullDoubleDiamond ||
         ms == kFullFourTrianglesPlus ||
         ms == kFullCrossX;
}


/**
 * Returns the "open" version of a "full" marker
 */
Style_t FullToOpenMarker (const Style_t ms) {
  switch (ms) {
    case kFullCircle: return kOpenCircle;
    case kFullSquare: return kOpenSquare;
    case kFullDiamond: return kOpenDiamond;
    case kFullCross: return kOpenCross;
    case kFullTriangleUp: return kOpenTriangleUp;
    case kFullTriangleDown: return kOpenTriangleDown;
    case kFullStar: return kOpenStar;
    case kFullCrossX: return kOpenCrossX;
    case kFullFourTrianglesPlus: return kOpenFourTrianglesPlus;
    case kFullFourTrianglesX: return kOpenFourTrianglesX;
    case kFullThreeTriangles: return kOpenThreeTriangles;
    case kFullDoubleDiamond: return kOpenDoubleDiamond;
    default: return kDot;
  }
}


/**
 * Returns the "full" version of an "open" marker
 */
Style_t OpenToFullMarker (const Style_t ms) {
  switch (ms) {
    case kOpenCircle: return kFullCircle;
    case kOpenSquare: return kFullSquare;
    case kOpenDiamond: return kFullDiamond;
    case kOpenCross: return kFullCross;
    case kOpenTriangleUp: return kFullTriangleUp;
    case kOpenTriangleDown: return kFullTriangleDown;
    case kOpenStar: return kFullStar;
    case kOpenCrossX: return kFullCrossX;
    case kOpenFourTrianglesPlus: return kFullFourTrianglesPlus;
    case kOpenFourTrianglesX: return kFullFourTrianglesX;
    case kOpenThreeTriangles: return kFullThreeTriangles;
    case kOpenDoubleDiamond: return kFullDoubleDiamond;
    default: return kDot;
  }
}


/**
 * Sets up a canvas with appropriate margins, useful for doing a colz plot
 */
void FormatTH2Canvas (TCanvas* c, const bool zAxisSpace) {
  if (zAxisSpace) c->SetRightMargin (0.18);
  else c->SetRightMargin (0.06);

  c->SetLeftMargin (0.15); 
}


/**
 * Returns a TGraphErrors with the same content as a TH1
 */
TGraphAsymmErrors* TH1TOTGraph (TH1 *h1) {
  if (!h1)
    std::cout << "TH1TOTGraph: histogram not found !" << std::endl;

  TGraphAsymmErrors* g1 = new TGraphAsymmErrors ();

  double x, y, ex, ey;
  for (int i = 1; i <= h1->GetNbinsX (); i++) {
    y = h1->GetBinContent (i);
    ey = h1->GetBinError (i);
    x = h1->GetBinCenter (i);
    ex = h1->GetBinWidth (i);
    //   cout << " x,y = " << x << " " << y << " ex,ey = " << ex << " " << ey << endl;

    g1->SetPoint (i-1, x, y);
    g1->SetPointError (i-1, ex, ex, ey, ey);
  }

  //g1->Print();
  return g1;
}



/**
 * Divides TH1 n over d assuming correlated binomial errors
 */
void BinomialDivide (TH1* n, TH1* d) {
  for (int ix = 1; ix <= n->GetNbinsX (); ix++) {
    float eff = n->GetBinContent (ix);
    if (d->GetBinContent (ix) != 0)
      eff = eff / d->GetBinContent (ix);

    float var = (eff - eff*eff);
    if (d->GetBinContent (ix) != 0)
      var = var / d->GetBinContent (ix);

    n->SetBinContent (ix, eff);
    n->SetBinError (ix, sqrt (fabs (var)));
  }
  return;
}



/**
 * The infamous myText draws a text label on your favorite plot
 */
void myText (double x, double y, Color_t color, const char* text, double tsize) {
  //double tsize=0.04;
  TLatex l;
  //l.SetTextAlign (12);
  l.SetTextSize (tsize); 
  l.SetNDC ();
  l.SetTextColor (color);
  l.DrawLatex (x, y, text);
}
 



void myBoxText (double x, double y, double boxsize, int mcolor, const char* text, const double tsize) {
  TLatex l;
  l.SetTextAlign (12);
  l.SetTextSize (tsize); 
  l.SetNDC ();
  l.DrawLatex (x, y, text);

  double y1 = y - 0.4*tsize;
  double y2 = y + 0.4*tsize;
  double x2 = x - 0.3*tsize;
  double x1 = x2 - boxsize;

  printf ("x1= %f x2= %f y1= %f y2= %f \n", x1, x2, y1, y2);

  TPave* mbox= new TPave(x1, y1, x2, y2, 0, "NDC");

  mbox->SetFillColor (mcolor);
  mbox->SetFillStyle (1001);
  mbox->Draw ();

  TLine mline;
  mline.SetLineWidth (4);
  mline.SetLineColor (1);
  mline.SetLineStyle (1);
  double y_new = (y1+y2)/2.;
  mline.DrawLineNDC (x1, y_new, x2, y_new);
}




void myBoxTextNoLine (double x, double y, double boxsize, int mcolor, const char* text, const double tsize) {
  TLatex l;
  l.SetTextAlign (12);
  l.SetTextSize (tsize); 
  l.SetNDC ();
  l.DrawLatex (x, y, text);

  double y1 = y - 0.4*tsize;
  double y2 = y + 0.4*tsize;
  double x2 = x - 0.3*tsize;
  double x1 = x2 - boxsize;

  printf ("x1= %f x2= %f y1= %f y2= %f \n", x1, x2, y1, y2);

  TPave* mbox= new TPave(x1, y1, x2, y2, 0, "NDC");

  mbox->SetFillColor (mcolor);
  mbox->SetFillStyle (1001);
  mbox->Draw ();
}




void myLineText (double x, double y, int color, int lstyle, const char* text, float lsize, double tsize) {
  TLine *markerLine = new TLine(x-(0.8*tsize)-0.02*lsize, y, x-(0.8*tsize)+0.02*lsize, y);
  markerLine->SetNDC ();
  markerLine->SetLineColor (color);
  markerLine->SetLineStyle (lstyle);
  markerLine->SetLineWidth (2);
  markerLine->Draw ();

  if (text[0] != '\0') {
    TLatex l;
    l.SetTextAlign (12);
    l.SetTextSize (tsize);
    l.SetNDC ();
    l.DrawLatex (x, y, text);
  }
}




void myLineColorText (double x, double y, int color, int lstyle, const char* text, float lsize, double tsize) {
  TLine *markerLine = new TLine(x-(0.8*tsize)-0.02*lsize, y, x-(0.8*tsize)+0.02*lsize, y);
  markerLine->SetNDC ();
  markerLine->SetLineColor (color);
  markerLine->SetLineStyle (lstyle);
  markerLine->SetLineWidth (2);
  markerLine->Draw ();

  if (text[0] != '\0') {
    TLatex l;
    l.SetTextAlign (12);
    l.SetTextSize (tsize);
    l.SetTextColor (color);
    l.SetNDC ();
    l.DrawLatex (x, y, text);
  }
}




void myMarkerText (double x, double y, int color, int mstyle, const char* text, float msize, double tsize, bool doOutline) {
//  double tsize=0.032;
  //TMarker* marker = new TMarker (x-(0.8*tsize), y, 8);
  TMarker* marker = new TMarker (x-(0.8*tsize), y+0.35*tsize, 8);
  marker->SetMarkerColor (color);
  marker->SetNDC ();
  marker->SetMarkerStyle (mstyle);
  marker->SetMarkerSize (msize);
  marker->Draw ();

  //TLine* markerLine = new TLine (x-(0.8*tsize)-0.02, y, x-(0.8*tsize)+0.02, y);
  TLine* markerLine = new TLine (x-(0.8*tsize)-0.02, y+0.35*tsize, x-(0.8*tsize)+0.02, y+0.35*tsize);
  markerLine->SetNDC ();
  markerLine->SetLineColor (color);
  markerLine->SetLineStyle (1);
  markerLine->SetLineWidth (2);
  markerLine->Draw ();

  if (doOutline && IsFullMarker (mstyle)) {
    //TMarker* marker2 = new TMarker (x-(0.8*tsize), y, 8);
    TMarker* marker2 = new TMarker (x-(0.8*tsize), y+0.35*tsize, 8);
    marker2->SetMarkerColor (kBlack);
    marker2->SetNDC ();
    marker2->SetMarkerStyle (FullToOpenMarker (mstyle));
    marker2->SetMarkerSize (msize);
    marker2->Draw ();
  }

  if (text[0] != '\0') {
    TLatex l;
    //l.SetTextAlign (12);
    l.SetTextSize (tsize); 
    l.SetNDC ();
    l.DrawLatex (x, y, text);
  }
}




void myMarkerTextNoLine (double x, double y, int color, int mstyle, const char* text, float msize, double tsize) {
//  double tsize=0.032;
  TMarker* marker = new TMarker(x-(0.8*tsize), y+0.35*tsize, 8);
  marker->SetMarkerColor (color);
  marker->SetNDC ();
  marker->SetMarkerStyle (mstyle);
  marker->SetMarkerSize (msize);
  marker->Draw ();

  if (text[0] != '\0') {
    TLatex l;
    //l.SetTextAlign (12);
    l.SetTextSize (tsize); 
    l.SetNDC ();
    l.DrawLatex (x,y,text);
  }
}




void myOnlyBoxText (double x, double y, double boxsize, int mcolor, int lcolor, int lstyle, const char* text, double tsize, int bstyle, double balpha) {
  double y1 = y-0.25*tsize;
  double y2 = y+0.25*tsize;
  double x2 = x-0.8*tsize+0.02*boxsize;
  double x1 = x-0.8*tsize-0.02*boxsize;
  //double x2 =  x-0.15*tsize;
  //double x1 = x-0.95*tsize;
  //printf("x1= %f x2= %f y1= %f y2= %f \n", x1, x2, y1, y2);
  TPave* mbox= new TPave (x1, y1, x2, y2, 0, "NDC");
  mbox->SetFillColorAlpha (mcolor, balpha);
  mbox->SetFillStyle (bstyle);
  mbox->Draw ();
  TLine mline;
  mline.SetLineWidth (1);
  mline.SetLineColor (lcolor);
  //mline.SetLineStyle (lstyle);
  mline.SetLineStyle (lstyle);
  //double y_new =  (y1+y2)/2.;
  //mline.DrawLineNDC (x1, y_new, x2, y_new);
  mline.DrawLineNDC (x1, y1, x2, y1);
  mline.DrawLineNDC (x1, y2, x2, y2);
  mline.DrawLineNDC (x1, y1, x1, y2);
  mline.DrawLineNDC (x2, y1, x2, y2);

  if (text[0] != '\0') {
    TLatex l;
    l.SetTextAlign (12);
    l.SetTextSize (tsize);
    l.SetNDC ();
    l.DrawLatex (x, y, text);
  }
}




TBox* TBoxNDC (const double x1, const double y1, const double x2, const double y2) {
  TPad* p;
  if (gDirectory->Get ("box_pad"))
    p = (TPad*) gDirectory->Get ("box_pad");
  else {
    p = new TPad ("box_pad", "box_pad", 0., 0., 1., 1.);
    p->SetFillStyle (0);
  } 
  p->Draw ();
  p->cd ();
  TBox* b = new TBox (x1, y1, x2, y2);
  return b;
}




void myMarkerAndBoxAndLineText (double x, double y, const double bsize, const int bstyle, const int bcolor, const double balpha, const int mcolor, const int mstyle, const double msize, const char* text, const double tsize) {
  const double y1 = y - (0.25*tsize) - (0.004*bsize) + 0.25*tsize;
  const double y2 = y + (0.25*tsize) + (0.004*bsize) + 0.25*tsize;
  const double x2 = x - (0.8*tsize) + 0.02;
  const double x1 = x - (0.8*tsize) + 0.02 - (0.04*bsize);

  TPave *mbox= new TPave (x1, y1, x2, y2, 0, "NDC");
  mbox->SetFillColorAlpha (bcolor, balpha);
  mbox->SetFillStyle (bstyle);
  mbox->Draw ();

  TLine mline;
  mline.SetLineWidth (1);
  mline.SetLineColor (mcolor);
  mline.SetLineStyle (1);
  mline.DrawLineNDC (x1, y1, x2, y1);
  mline.DrawLineNDC (x1, y2, x2, y2);
  mline.DrawLineNDC (x1, y1, x1, y2);
  mline.DrawLineNDC (x2, y1, x2, y2);

  if (mstyle != -1) {

    TMarker *marker = new TMarker(x - (0.8*tsize)+0.02-0.02*bsize, y+0.25*tsize, 8);
    marker->SetNDC();
    marker->SetMarkerColor (IsOpenMarker (mstyle) ? kBlack : mcolor);
    marker->SetMarkerStyle (mstyle);
    marker->SetMarkerSize (msize);

    //TLine *markerLine = new TLine(marker->GetX()-0.018*msize, marker->GetY(), marker->GetX()+0.018*msize, marker->GetY());
    TLine *markerLine = new TLine ();
    markerLine->SetNDC();
    markerLine->SetLineColor(mcolor);
    markerLine->SetLineStyle(1);
    markerLine->SetLineWidth(2);

    markerLine->DrawLineNDC (0.9*x1+0.1*x2, 0.5*(y1+y2), 0.1*x1+0.9*x2, 0.5*(y1+y2));
    markerLine->DrawLineNDC (0.5*(x1+x2), 0.9*y1+0.1*y2, 0.5*(x1+x2), 0.1*y1+0.9*y2);

    //if (IsOpenMarker (mstyle) {
    //  Rectangle_t bb = marker->GetBBox ();

    //  const double xbbl_user = gPad->PixelToUser (bb.fX);
    //  const short xbbl = bb.fX- 0.5*bb.fWidth ();
    //  const short ybbl = bb.fY + 0.5*bb.fWidth ();

    //  
    //}

    marker->Draw();

    if (FullToOpenMarker (mstyle) != kDot) {
      TMarker *marker2 = new TMarker(x - (0.8*tsize)+0.02-0.02*bsize, y+0.25*tsize, 8);
      marker2->SetNDC ();
      marker2->SetMarkerColor (kBlack);
      marker2->SetMarkerStyle (FullToOpenMarker (mstyle));
      marker2->SetMarkerSize (msize);
      marker2->Draw ();
    }
  }
  
  TLatex l;
  l.SetTextAlign (11);
  l.SetTextSize (tsize);
  l.SetNDC ();
  l.DrawLatex (x, y, text);
}

#endif
