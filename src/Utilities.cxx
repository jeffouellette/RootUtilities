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
#include <TGraphSmooth.h>

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
      const double factorOfTen = std::pow (10, n-valDec); // e.g., for "1520" and n=2, get 0.01
      val = floor (factorOfTen * val + 0.5) / factorOfTen;
    }
    else { // if decimal is before least significant digit, e.g. for "15.24" and n=3.
      const double factorOfTen = std::pow (10, n-valDec+1); // e.g. for 15.24 and n=3 get 10
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
    const double factorOfTen = std::pow (10, errStart+n);
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
      const double factorOfTen = std::pow (10., n - errDec);
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
  while (std::fabs (stat) < std::pow (10, n-1) || std::fabs (syst) < std::pow (10, n-1) && i < 10) { // i < 10 as an arbitrary cut (so this doesn't go forever)
    val *= 10;
    stat *= 10;
    syst *= 10;
    i++;
  }

  val = round (val);
  stat = round (stat);
  syst = round (syst);

  val = val / std::pow (10, i);
  stat = stat / std::pow (10, i);
  syst = syst / std::pow (10, i);

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
 * Sets the errors in h assuming the variances are stored in the entries of h2.
 */
void SetVariances (TH1D* h, TH2D* h2) {
  const int nb = h->GetNbinsX ();
  assert (nb == h2->GetNbinsX () && nb == h2->GetNbinsY ());

  for (int iX = 1; iX <= nb; iX++) {
    if (isnan (h2->GetBinContent (iX, iX)) || h2->GetBinContent (iX, iX) < 0)
      std::cout << "Problem detected with covariance matrix " << h2->GetName () << std::endl; 
    assert (h2->GetBinContent (iX, iX) >= 0);
    h->SetBinError (iX, std::sqrt (h2->GetBinContent (iX, iX)));
  }
}


/**
 * Divides two histograms assuming the entries in one are binomial samples of the other.
 */
void BinomialDivide (TH1* out, TH1* num, TH1* den, TH1* h_sumwgt2) {
  assert (out->GetNbinsX () == num->GetNbinsX ());
  assert (out->GetNbinsX () == den->GetNbinsX ());
  assert (h_sumwgt2 == nullptr || out->GetNbinsX () == h_sumwgt2->GetNbinsX ());
  assert (out->GetNbinsY () == num->GetNbinsY ());
  assert (out->GetNbinsY () == den->GetNbinsY ());
  assert (h_sumwgt2 == nullptr || out->GetNbinsY () == h_sumwgt2->GetNbinsY ());

  for (int iX = 1; iX <= out->GetNbinsX (); iX++) {
    for (int iY = 1; iY <= out->GetNbinsY (); iY++) {

      const double passes = num->GetBinContent (iX, iY);
      const double trials = den->GetBinContent (iX, iY);

      if (trials > 0) {
        const double sumwgt2 = h_sumwgt2->GetBinContent (iX, iY);
        const double rate = passes/trials;
        out->SetBinContent (iX, iY, rate);
        out->SetBinError (iX, iY, std::sqrt ((rate*(1.-rate))*sumwgt2/std::pow(trials,2))); // normal approximation to uncertainty on efficiency
      }
      else {
        out->SetBinContent (iX, iY, -1); // set to -1 so it's not plotted
        out->SetBinError (iX, iY, 1);    // set to 1 since we don't know what it is
      }
    }
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
  double dphi = std::abs(phi1 - phi2);
  while (dphi > M_PI) dphi = std::abs (dphi - 2*M_PI);

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
 return std::sqrt( std::pow( deta, 2 ) + std::pow( dphi, 2 ) );
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
 * Cuts off extra points in a TGraph beyond a specified range.
 */
void TrimGraph (TGraph* g, const double xmin, const double xmax) {
  double x, y;
  for (int i = 0; i< g->GetN (); i++) {
    g->GetPoint (i, x, y);
    if (x < xmin || x > xmax) {
      g->RemovePoint (i);
      i--;
    }
  }
  return;
}


/**
 * Calculates errors in a TH1D from a TH2D storing the quadrature sum of each entry.
 */
void CalcUncertainties (TH1D* h, TH2D* h2, const double n) {
  if (n <= 1)
    std::cout << "Warning: n <= 1 for histograms " << h->GetName () << " and " << h2->GetName () << std::endl;
  if (n == 0)
    return; // if no counts just return
  assert (n > 1); // check number of entries is at least 1 (otherwise the sample variance is undefined!)
  assert (h->GetNbinsX () == h2->GetNbinsX ());  // check h2 is a valid sum(x^2) of h.
  assert (h2->GetNbinsX () == h2->GetNbinsY ()); // check h2 is "square" -- same bins in X & Y.
  h->Scale (1./n);
  for (int iX = 1; iX <= h2->GetNbinsX (); iX++)
    for (int iY = 1; iY <= h2->GetNbinsY (); iY++)
      h2->SetBinContent (iX, iY, h2->GetBinContent (iX, iY) - n * (h->GetBinContent (iX)) * (h->GetBinContent (iY)));
  h->Scale (1.);
  h2->Scale (std::pow (n*(n-1.), -1));
  SetVariances (h, h2);
  return;
}


/**
 * Calculates errors in a TH1D from a TH2D storing the quadrature sum of each entry. Assumes weights are used, so a histogram with the number of events to the powers 0, 1, and 2 is needed.
 */
void CalcUncertainties (TH1D* h, TH2D* h2, TH1D* hn) {
  assert (h->GetNbinsX () == h2->GetNbinsX ());  // check h2 is a valid sum(x^2) of h.
  assert (h2->GetNbinsX () == h2->GetNbinsY ()); // check h2 is "square" -- same bins in X & Y.
  assert (hn->GetNbinsX () >= 3); // check counts histogram has at least 3 bins

  const double n0 = hn->GetBinContent (1);
  const double n1 = hn->GetBinContent (2);
  const double n2 = hn->GetBinContent (3);

  if (n0 <= 1) {
    if (n0 == 0)
      return; // if no counts just return
    else if (n0 == 1)
      std::cout << "Warning: n0 == 1 for histograms " << h->GetName () << " and " << h2->GetName () << std::endl;
    else
      std::cout << "Warning: n0 invalid for histograms " << h->GetName () << " and " << h2->GetName () << std::endl;
  }
  assert (n0 >= 1); // check number of entries is greater than 1; if equal to 1 then define all uncertainties to be 100%

  h->Scale (1./n1);
  for (int iX = 1; iX <= h2->GetNbinsX (); iX++)
    for (int iY = 1; iY <= h2->GetNbinsY (); iY++)
      h2->SetBinContent (iX, iY, (n0 == 1 ? h->GetBinContent (iX) * h->GetBinContent (iY) : h2->GetBinContent (iX, iY) - n1 * (h->GetBinContent (iX)) * (h->GetBinContent (iY))));
  //h->Scale (1., "width");
  if (n0 > 1)
    h2->Scale (std::pow (n0*(n1*n1-n2)/n1, -1));//, "width");
  //else
  //  h2->Scale (std::pow (n1, -1));
  SetVariances (h, h2);
  return;
}


/**
 * Adds the maximum systematic uncertainty from v1 and v2 to sys in quadrature.
 */
void AddMaxSystematic (TGAE* sys, TGAE* v1, TGAE* v2) {
  assert (sys->GetN () == v1->GetN ());
  assert (sys->GetN () == v2->GetN ());
  for (int ix = 0; ix < sys->GetN (); ix++) {
    double upErr = fmax (v1->GetErrorYhigh (ix), v2->GetErrorYhigh (ix));
    double downErr = fmax (v1->GetErrorYlow (ix), v2->GetErrorYlow (ix));
    sys->SetPointEYhigh (ix, std::sqrt (std::pow (sys->GetErrorYhigh (ix), 2) + std::pow (upErr, 2)));
    sys->SetPointEYlow (ix, std::sqrt (std::pow (sys->GetErrorYlow (ix), 2) + std::pow (downErr, 2)));
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

    const float newErr = std::sqrt (std::pow (master->GetBinError (ix), 2) + std::pow (sys->GetBinError (ix), 2));
    master->SetBinError (ix, newErr);
  }
}


/**
 * Adds independent systematic errors in quadrature, storing the sum in master
 */
void AddErrorsInQuadrature (TGAE* master, TH1D* nom, TH1D* sys) {
  assert (nom->GetNbinsX () == sys->GetNbinsX ());
  assert (nom->GetNbinsX () == master->GetN ());
  double eym, eys, newErr;
  for (int ix = 0; ix < master->GetN (); ix++) {
    eys = sys->GetBinContent (ix+1) - nom->GetBinContent (ix+1);

    {
      eym = master->GetErrorYhigh (ix);
      newErr = std::sqrt (eym*eym + eys*eys);
      master->SetPointEYhigh (ix, newErr);
    }
    {
      eym = master->GetErrorYlow (ix);
      newErr = std::sqrt (eym*eym + eys*eys);
      master->SetPointEYlow (ix, newErr);
    }
  }
  return;
}


/**
 * Adds independent systematic errors in quadrature, storing the sum in master
 */
void AddErrorsInQuadrature (TGAE* master, TGAE* sys, const bool doXErrs, const bool symmetrize) {
  for (int ix = 0; ix < master->GetN (); ix++) {
    //double xm=0, ym=0, xs=0, ys=0;
    //master->GetPoint (ix, xm, ym);
    //sys->GetPoint (ix, xs, ys);

    //if (xm != xs || ym != ys)
    //  cout << "Warning: In Utilities.cxx::AddErrorsInQuadrature: Bins don't match!" << endl;
    //
    float newErr = std::sqrt (std::pow (master->GetErrorYhigh (ix), 2) + std::pow (std::fmax (sys->GetErrorYhigh (ix), symmetrize ? sys->GetErrorYlow (ix) : 0), 2));
    master->SetPointEYhigh (ix, newErr);
    newErr = std::sqrt (std::pow (master->GetErrorYlow (ix), 2) + std::pow (std::fmax (sys->GetErrorYlow (ix), symmetrize ? sys->GetErrorYhigh (ix) : 0), 2));
    master->SetPointEYlow (ix, newErr);

    if (doXErrs) {
      newErr = std::sqrt (std::pow (master->GetErrorXhigh (ix), 2) + std::pow (sys->GetErrorXhigh (ix), 2));
      master->SetPointEXhigh (ix, newErr);
      newErr = std::sqrt (std::pow (master->GetErrorXlow (ix), 2) + std::pow (sys->GetErrorXlow (ix), 2));
      master->SetPointEXlow (ix, newErr);
    }
  }
}


/**
 * Adds independent systematic errors in quadrature, storing the sum in master
 * Allows user to send an array of systematics, not necessarily in order, and only adds the maximum uncertainty from these.
 */
void AddErrorsInQuadrature (TGAE* master, TGAE** sys, const std::vector <int>* indices, const bool symmetrize) {
  assert (indices->size () >= 1);
  
  std::vector <double> eyhi (master->GetN ());
  std::vector <double> eylo (master->GetN ());

  for (int index : *indices) {
    for (int ix = 0; ix < sys[index]->GetN (); ix++) {
      eyhi[ix] = std::fmax (eyhi[ix], std::fmax (sys[index]->GetErrorYhigh (ix), symmetrize ? sys[index]->GetErrorYlow (ix) : 0));
      eylo[ix] = std::fmax (eylo[ix], std::fmax (sys[index]->GetErrorYlow (ix), symmetrize ? sys[index]->GetErrorYhigh (ix) : 0));
    }
  }

  for (int ix = 0; ix < master->GetN (); ix++) {
    float newErr = std::sqrt (std::pow (master->GetErrorYhigh (ix), 2) + std::pow (eyhi[ix], 2));
    master->SetPointEYhigh (ix, newErr);
    newErr = std::sqrt (std::pow (master->GetErrorYlow (ix), 2) + std::pow (eylo[ix], 2));
    master->SetPointEYlow (ix, newErr);
  }
}


/**
 * Adds independent relative systematic errors in quadrature, storing the sum in master
 */
void AddRelErrorsInQuadrature (TGAE* master, TGAE* sys, const bool ignoreWarning, const bool symmetrize) {
  assert (master->GetN () == sys->GetN ());
  double xm, ym, xs, ys, eym, eys, newErr;
  for (int ix = 0; ix < master->GetN (); ix++) {

    master->GetPoint (ix, xm, ym);
    sys->GetPoint (ix, xs, ys);

    if (xm != xs && !ignoreWarning)
      cout << "Warning: In Utilities.cxx::AddRelErrorsInQuadrature: Bins don't match!" << endl;

    {
      eym = master->GetErrorYhigh (ix);
      eys = std::fmax (sys->GetErrorYhigh (ix), symmetrize ? sys->GetErrorYlow (ix) : 0);
      // convert to relative uncertainties
      if (ym != 0) eym = eym / ym;
      if (ys != 0) eys = eys / ys;
      newErr = std::sqrt (eym*eym + eys*eys);
      master->SetPointEYhigh (ix, newErr);
    }
    {
      eym = master->GetErrorYlow (ix);
      eys = std::fmax (sys->GetErrorYlow (ix), symmetrize ? sys->GetErrorYhigh (ix) : 0);
      // convert to relative uncertainties
      if (ym != 0) eym = eym / ym;
      if (ys != 0) eys = eys / ys;
      newErr = std::sqrt (eym*eym + eys*eys);
      master->SetPointEYlow (ix, newErr);
    }
  }
  return;
}


/**
 * Calculates simple systematics as maximum variations on the nominal.
 * Intended for combining up/down variations in an expandable way.
 */
void CalcSystematics (TH1D* sys, TH1D* var) {
  for (int ix = 1; ix <= sys->GetNbinsX (); ix++) {
    const float newErr = std::fabs (var->GetBinContent (ix) - sys->GetBinContent (ix));
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
      sys->SetPointEYhigh (ix, fmax (std::fabs (newErr), sys->GetErrorYhigh (ix)));
    }
    if (applyBothWays || newErr < 0) {
      sys->SetPointEYlow (ix, fmax (std::fabs (newErr), sys->GetErrorYlow (ix)));
    }
  }
}


/**
 * Calculates simple systematics as maximum variations on the nominal.
 * Intended for combining up/down variations in an expandable way.
 */
void CalcSystematics (TGAE* sys, TH1D* nom, TH1D* var) {
  assert (nom->GetNbinsX () == var->GetNbinsX ());
  for (int ix = 0; ix < nom->GetNbinsX (); ix++) {
    double x = nom->GetBinCenter (ix+1);
    double y = nom->GetBinContent (ix+1);

    sys->SetPoint (ix, x, y);

    const float newErr = var->GetBinContent (ix+1) - y;

    if (newErr > 0) sys->SetPointEYhigh (ix, newErr);
    else            sys->SetPointEYlow (ix, -newErr);
    sys->SetPointEXhigh (ix, nom->GetBinWidth (ix+1) * 0.5);
    sys->SetPointEXlow (ix, nom->GetBinWidth (ix+1) * 0.5);
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
      graph->SetPointEYlow (ix-1, std::sqrt (0.5*(err_lo*err_lo + err_hi*err_hi)));
      graph->SetPointEYhigh (ix-1, std::sqrt (0.5*(err_lo*err_lo + err_hi*err_hi)));
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
      const double xerr = fmax (std::fabs (xh - x), std::fabs (x - xl));
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
      graph->SetPointEYlow (ix, std::sqrt (0.5*(err_lo*err_lo + err_hi*err_hi)));
      graph->SetPointEYhigh (ix, std::sqrt (0.5*(err_lo*err_lo + err_hi*err_hi)));
    }
  }
  return;
}


/**
 * Extension of CalcSystematics (TGAE* sys, TH1D* nom, TH1D* var) for smoothing uncertainties.
 */
void SmoothSystematics (TGAE* sys, TH1D* nom, TH1D* var, const TString funcform) {
  TGraphErrors* tg = new TGraphErrors ();
  double x, y;
  double xlo = DBL_MAX, xhi = DBL_MIN;
  for (int i = 0; i < sys->GetN (); i++) {
    sys->GetPoint (i, x, y);
    xlo = std::fmin (x, xlo);
    xhi = std::fmax (x, xhi);
    if (y != 0) {
      tg->SetPoint (i, x, (sys->GetErrorYhigh (i) - sys->GetErrorYlow (i)) / y);
      double yv = var->GetBinContent (i+1);
      double yve = var->GetBinError (i+1);
      double yn = nom->GetBinContent (i+1);
      double yne = nom->GetBinError (i+1);
      tg->SetPointError (i, 0.5*nom->GetBinWidth (i+1), std::sqrt (std::pow (yv*yne/(yn*yn), 2) + std::pow (yve/yn, 2)));
    }
  }

  TF1* func = new TF1 ("functemp", funcform.Data (), xlo, xhi);
  tg->Fit (func, "RN0Q");

  for (int i = 0; i < sys->GetN (); i++) {
    sys->GetPoint (i, x, y);
    double newErr = func->Eval (x) * y;
    if (newErr > 0) {
      sys->SetPointEYhigh (i, std::fabs (newErr));
      sys->SetPointEYlow (i, 0);
    }
    else {
      sys->SetPointEYlow (i, std::fabs (newErr));
      sys->SetPointEYhigh (i, 0);
    }
  } 
  SaferDelete (&func);
  SaferDelete (&tg);
  return;
}


/**
 * Extension of CalcSystematics (TGAE* sys, TH1D* nom, TH1D* var) for smoothing uncertainties.
 * Uses LOWESS regression from TGraphSmooth.
 */
void SmoothSystematics (TGAE* sys) { //, TH1D* nom, TH1D* var) {
  TGraph* tg = new TGraph ();
  double x, y;
  double xlo = DBL_MAX, xhi = DBL_MIN;
  for (int i = 0; i < sys->GetN (); i++) {
    sys->GetPoint (i, x, y);
    xlo = std::fmin (x, xlo);
    xhi = std::fmax (x, xhi);
    if (y != 0) {
      tg->SetPoint (i, std::log (x), (sys->GetErrorYhigh (i) - sys->GetErrorYlow (i)) / y);
      //double yv = var->GetBinContent (i+1);
      //double yve = var->GetBinError (i+1);
      //double yn = nom->GetBinContent (i+1);
      //double yne = nom->GetBinError (i+1);
      //tg->SetPointError (i, 0.5*nom->GetBinWidth (i+1), std::sqrt (std::pow (yv*yne/(yn*yn), 2) + std::pow (yve/yn, 2)));
    }
  }
  TGraphSmooth* gsmooth = new TGraphSmooth ("normal");
  TGraph* gs = gsmooth->SmoothLowess (tg, "", 0.67, 5);

  for (int i = 0; i < sys->GetN (); i++) {
    sys->GetPoint (i, x, y);
    double newErr = gs->Eval (std::log (x)) * y;
    if (newErr > 0) {
      sys->SetPointEYhigh (i, std::fabs (newErr));
      sys->SetPointEYlow (i, 0);
    }
    else {
      sys->SetPointEYlow (i, std::fabs (newErr));
      sys->SetPointEYhigh (i, 0);
    }
  } 
  SaferDelete (&gsmooth);
  SaferDelete (&tg);
  return;
}


/**
 * Sets the bin contents in target as the errors in errors / central values in centralValues
 */
void SaveRelativeErrors (TGAE* target, TGAE* errors, TGAE* centralValues, const float sf) {
  assert (target->GetN () == errors->GetN ());
  assert (target->GetN () == centralValues->GetN ());

  double x, y, eyhi, eylo;
  for (int i = 0; i < target->GetN (); i++) {
    centralValues->GetPoint (i, x, y);
    if (y == 0) errors->GetPoint (i, x, y); // if 0 then take y from errors instead
    if (y == 0) y = 1; // if still 0 set to 1 to avoid divide by 0
    eyhi = errors->GetErrorYhigh (i);
    eylo = errors->GetErrorYlow (i);
    if (eyhi > eylo)      target->SetPoint (i, x, sf * eyhi / y);
    else if (eylo > eyhi) target->SetPoint (i, x, -sf * eylo / y);
    else if (eylo == 0)   target->SetPoint (i, x, 0);
    else std::cout << "not sure waht to do, please instruct!" << std::endl;
  }
  return;
}


/**
 * Sets the bin contents in target as the error in errors / central values in centralValues
 */
void SaveRelativeErrors (TH1D* errors, TH1D* centralValues, bool useCentVals) {
  assert (errors->GetNbinsX () == centralValues->GetNbinsX ());
  for (int ix = 1; ix <= errors->GetNbinsX (); ix++) {
    if (centralValues->GetBinContent (ix) != 0)
      errors->SetBinContent (ix, (useCentVals ? errors->GetBinContent (ix) : errors->GetBinError (ix)) / centralValues->GetBinContent (ix));
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
 * Returns a TGAE with central values equal to the uncertainties in g.
 */
TGAE* ConvertErrorsToCentralValues (TGAE* g, const bool doHighErrs, const float sf) {
  TGAE* gnew = (TGAE*) g->Clone ();

  ResetTGAEErrors (gnew);
  ResetXErrors (gnew);
  double x, y;
  for (int ix = 0; ix < gnew->GetN (); ix++) {
    gnew->GetPoint (ix, x, y);
    if (y != 0)
      gnew->SetPoint (ix, x, sf * (doHighErrs ? g->GetErrorYhigh (ix) : -g->GetErrorYlow (ix)) / y);
  }

  return gnew;
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
     var += std::pow (h3->GetBinError (h3ix, h3iy, h3iz), 2);

    } // end loop over extra axis

    h2->SetBinContent (ix, iy, content);
    h2->SetBinError (ix, iy, std::sqrt (var));
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

    if (logx && useX2OverX1Factor) tg->SetPoint (n, x * std::exp (0.5 * delta * std::log (x2/x1)), y);
    else if (logx) tg->SetPoint (n, x*delta, y);
    else tg->SetPoint (n, x + delta, y);

    tg->GetPoint (n, x, y);
    exh = std::fabs (exh - x);
    exl = std::fabs (x - exl);

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
      exh = std::fabs (x * (std::exp (0.5*err * std::log (x2/x1)) - 1));
      exl = std::fabs (x * (1 - std::exp (-0.5*err * std::log(x2/x1))));
    }
    else if (mult) {
      exh = std::fabs (x * (std::exp (0.5*err) - 1));
      exl = std::fabs (x * (1 - std::exp (-0.5*err)));
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
 * Makes a TGraphError from the input histogram.
 */
TGE* TH1ToTGE (TH1* h, const float cutoff) {
  TGE* tg = new TGE ();

  const float xlo = h->GetBinLowEdge (1);
  const float xhi = h->GetBinLowEdge (h->GetNbinsX ()) + h->GetBinWidth (h->GetNbinsX ());

  for (int n = 0; n < h->GetNbinsX (); n++) {
    if (cutoff >= 0 && h->GetBinContent (n+1) <= cutoff) {
      continue;
    }
    else {
      tg->SetPoint (tg->GetN (), h->GetBinCenter (n+1), h->GetBinContent (n+1));
      tg->SetPointError (tg->GetN ()-1, h->GetBinWidth (n+1) / 2, h->GetBinError (n+1));
    }
  }

  tg->GetXaxis ()->SetRangeUser (xlo, xhi);
  return tg;
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
    logDelta = std::log10 (xhi) - std::log10 (xlo);
    newx = std::pow (10, log10 (xlo) + 0.5*logDelta);
    g->SetPoint (n, newx, y);
    g->SetPointEXlow (n, newx-xlo);
    g->SetPointEXhigh (n, xhi-newx);
  }
}



/**
 * Scales a TGAE and its errors by the central values in a TH1D and/or a constant.
 */
void ScaleGraph (TGAE* g, const TH1D* h, const double sf) {
  if (h == nullptr) {
    double x, y;
    for (int i = 0; i < g->GetN (); i++) {
      g->GetPoint (i, x, y);
      g->SetPoint (i, x, y * sf);
      g->SetPointEYlow (i, g->GetErrorYlow (i) * sf);
      g->SetPointEYhigh (i, g->GetErrorYhigh (i) * sf);
    }
  }
  else {
    assert (g->GetN () == h->GetNbinsX ());
    double x, y, yhist;
    for (int i = 0; i < g->GetN (); i++) {
      yhist = h->GetBinContent (i+1);
      if (yhist == 0)
        continue;
      g->GetPoint (i, x, y);
      g->SetPoint (i, x, y * sf / yhist);
      g->SetPointEYlow (i, g->GetErrorYlow (i) * sf / yhist);
      g->SetPointEYhigh (i, g->GetErrorYhigh (i) * sf / yhist);
    }
  }
  return;
}



/**
 * Scales a TGAE and its errors by the central values in another TGraph and/or a constant.
 */
void ScaleByGraph (TGAE* g, const TGraph* den, const double sf) {
  assert (g->GetN () == den->GetN ());
  double x, y, yden;
  for (int i = 0; i < g->GetN (); i++) {
    den->GetPoint (i, x, yden);
    if (yden == 0)
      continue;
    g->GetPoint (i, x, y);
    g->SetPoint (i, x, y * sf / yden);
    g->SetPointEYlow (i, g->GetErrorYlow (i) * sf / yden);
    g->SetPointEYhigh (i, g->GetErrorYhigh (i) * sf / yden);
  }
  return;
}



/**
 * Applies new binning to a histogram
 * BE CAREFUL: if bins edges don't overlap, this can lead to unexpected behavior!
 */ 
void RebinSomeBins (TH1D** _h, int nbins, double* bins, const bool doWidths) {
  TH1D* h = (*_h);

  if (nbins == h->GetNbinsX ()) {
    if (nbins > h->GetNbinsX ())
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
        hnew->SetBinContent (ix+1, hnew->GetBinContent (ix+1) + h->GetBinContent (ixprime+1)*(doWidths ? h->GetBinWidth (ixprime+1) : 1));
        hnew->SetBinError (ix+1, hnew->GetBinError (ix+1) + std::pow (h->GetBinError (ixprime+1)*(doWidths ? h->GetBinWidth (ixprime+1) : 1), 2));
      }
      else if (oldbins[ixprime] <= bins[ix] && bins[ix+1] <= oldbins[ixprime+1]) {
        hnew->SetBinContent (ix+1, hnew->GetBinContent (ix+1) + h->GetBinContent (ixprime+1)*(doWidths ? h->GetBinWidth (ixprime+1) : 1));
        hnew->SetBinError (ix+1, hnew->GetBinError (ix+1) + std::pow (h->GetBinError (ixprime+1)*(doWidths ? h->GetBinWidth (ixprime+1) : 1), 2));
      }
    }
    hnew->SetBinError (ix+1, std::sqrt (hnew->GetBinError (ix+1)));
  }
  delete[] oldbins;

  delete h;

  if (doWidths)
    hnew->Scale (1, "width");

  *_h = hnew;
}



/**
 * Applies new binning to a 2D histogram
 * Assumes histogram has same bins along X & Y.
 * BE CAREFUL: if bins edges don't overlap, this can lead to unexpected behavior!
 */ 
void RebinSomeBins2D (TH2D** _h, int nbins, double* bins, const bool doWidths) {
  TH2D* h = (*_h);

  if (nbins >= h->GetNbinsX ()) {
    cout << "More new bins than old bins, returning." << endl;
    return;
  }

  const TString name = h->GetName ();
  h->SetName ("temp");

  int noldbins = h->GetNbinsX ();
  double* oldbins = new double[noldbins+1];
  for (int ix = 1; ix <= noldbins; ix++) {
    oldbins[ix-1] = h->GetXaxis ()->GetBinLowEdge (ix);
  }
  oldbins[noldbins] = h->GetXaxis ()->GetBinLowEdge (noldbins) + h->GetXaxis ()->GetBinWidth (noldbins);

  TH2D* hnew = new TH2D (name, "", nbins, bins, nbins, bins);
  hnew->Sumw2 ();
  for (int ix = 0; ix < nbins; ix++) {
    for (int iy = 0; iy < nbins; iy++) {
      for (int ixprime = 0; ixprime < noldbins; ixprime++) {
        for (int iyprime = 0; iyprime < noldbins; iyprime++) {
          if (bins[ix] <= oldbins[ixprime] && oldbins[ixprime+1] <= bins[ix+1] && bins[iy] <= oldbins[iyprime] && oldbins[iyprime+1] <= bins[iy+1]) {
            hnew->SetBinContent (ix+1, iy+1, hnew->GetBinContent (ix+1, iy+1) + h->GetBinContent (ixprime+1, iyprime+1)*(doWidths ? h->GetXaxis ()->GetBinWidth (ixprime+1) * h->GetYaxis ()->GetBinWidth (iyprime+1) : 1));
            hnew->SetBinError (ix+1, iy+1, hnew->GetBinError (ix+1, iy+1) + std::pow (h->GetBinError (ixprime+1, iyprime+1)*(doWidths ? h->GetXaxis ()->GetBinWidth (ixprime+1) * h->GetYaxis ()->GetBinWidth (iyprime+1) : 1), 2));
          }
          else if (oldbins[ixprime] <= bins[ix] && bins[ix+1] <= oldbins[ixprime+1] && oldbins[iyprime] <= bins[iy] && bins[iy+1] <= oldbins[iyprime+1]) {
            hnew->SetBinContent (ix+1, iy+1, hnew->GetBinContent (ix+1, iy+1) + h->GetBinContent (ixprime+1, iyprime+1)*(doWidths ? h->GetXaxis ()->GetBinWidth (ixprime+1) * h->GetYaxis ()->GetBinWidth (iyprime+1) : 1));
            hnew->SetBinError (ix+1, iy+1, hnew->GetBinError (ix+1, iy+1) + std::pow (h->GetBinError (ixprime+1, iyprime+1)*(doWidths ? h->GetXaxis ()->GetBinWidth (ixprime+1) * h->GetYaxis ()->GetBinWidth (iyprime+1) : 1), 2));
          }
        }
        hnew->SetBinError (ix+1, iy+1, std::sqrt (hnew->GetBinError (ix+1, iy+1)));
      }
    }
  }
  delete[] oldbins;

  delete h;

  if (doWidths)
    hnew->Scale (1, "width");

  *_h = hnew;
}



/**
 * Un-scales a histogram by the bin width and also applies an optional constant scaling factor (such as a number of events or a luminosity).
 */
void UnscaleWidth (TH1D* h, const float sf) {
  for (short iX = 1; iX <= h->GetNbinsX (); iX++) {
    h->SetBinContent (iX, h->GetBinContent (iX) * h->GetXaxis ()->GetBinWidth (iX) * sf);
    h->SetBinError   (iX, h->GetBinError   (iX) * h->GetXaxis ()->GetBinWidth (iX) * sf);
  }
  return;
}



/**
 * Scales a histogram by the bin width and also applies an optional constant scaling factor (such as a number of events or a luminosity).
 */
void ScaleWidth (TH1D* h, const float sf) {
  for (short iX = 1; iX <= h->GetNbinsX (); iX++) {
    h->SetBinContent (iX, h->GetBinContent (iX) / (h->GetXaxis ()->GetBinWidth (iX) * sf));
    h->SetBinError   (iX, h->GetBinError   (iX) / (h->GetXaxis ()->GetBinWidth (iX) * sf));
  }
  return;
}



/**
 * Draws histogram as a graph with some plotting settings.
 */
void myDraw (TH1D* h, const Color_t col, const Style_t mstyle, const float msize, const Style_t lstyle, const int lwidth, const bool doMOutline, const char* opt) {
  TGAE* g = make_graph (h);
  myDraw (g, col, mstyle, msize, lstyle, lwidth, opt, doMOutline);
  SaferDelete (&g);
  return;
} 



/**
 * Draws histogram as a graph with some plotting settings.
 */
void myDrawHist (TH1D* h, const Color_t col, const Style_t lstyle, const int lwidth, const char* opt) {
  TH1D* hdraw = (TH1D*) h->Clone ();
  hdraw->SetLineColor (col);
  hdraw->SetLineStyle (lstyle);
  hdraw->SetLineWidth (lwidth);
  hdraw->DrawCopy (opt);
  delete hdraw;
  return;
} 



/**
 * Draws a graph with some plotting settings.
 */
void myDraw (TGraph* g, const Color_t col, const Style_t mstyle, const float msize, const Style_t lstyle, const int lwidth, const char* opt, const bool doMOutline) {
  g->SetLineColor (col);
  g->SetMarkerColor (col);
  g->SetMarkerStyle (mstyle);
  g->SetMarkerSize (msize);
  g->SetLineStyle (lstyle);
  g->SetLineWidth (lwidth);
  ((TGAE*) g->Clone ())->Draw (opt);

  if (doMOutline && IsFullMarker (mstyle) && col != kBlack) {
    g->SetLineWidth (0);
    g->SetMarkerColor (kBlack);
    g->SetMarkerStyle (FullToOpenMarker (g->GetMarkerStyle ()));
    ((TGAE*) g->Clone ())->Draw (opt);
  }
  return;
} 



/**
 * Draws a function with some plotting settings.
 */
void myDraw (TF1* f, const Color_t col, const Style_t lstyle, const int lwidth, const char* opt) {
  f->SetLineColor (col);
  f->SetLineStyle (lstyle);
  f->SetLineWidth (lwidth);
  ((TF1*) f->Clone ())->Draw (opt);
  return;
} 



/**
 * Draws a filled area between gup and gdown with the given settings.
 */
void myDrawFill (TGraph* gup, TGraph* gdown, const Color_t col, const float falpha, const Style_t fstyle) {
  assert (gup->GetN () == gdown->GetN ());
  TGAE* g = new TGAE (2 * gup->GetN ());

  double x, y;
  for (int i = 0; i < gup->GetN (); i++) {
    gup->GetPoint (i, x, y);
    g->SetPoint (i, x, y);
  }
  for (int i = gdown->GetN ()-1; i >= 0; i--) {
    gdown->GetPoint (i, x, y);
    g->SetPoint (2 * gdown->GetN () - i - 1, x, y);
  }

  g->SetFillColorAlpha (col, falpha);
  g->SetFillStyle (fstyle);
  g->SetLineWidth (0);
  ((TGAE*) g->Clone ())->Draw ("f");
  delete g;
  return;
}



/**
 * Draws a graph as a systematic with some plotting settings.
 */
void myDrawSyst (TGAE* g, const Color_t col, const Style_t lstyle, const int lwidth, const float falpha, const char* opt) {
  g->SetLineColor (col);
//  g->SetMarkerStyle (0);
//  g->SetMarkerSize (0);
  g->SetLineStyle (lstyle);
  g->SetLineWidth (lwidth);
  g->SetFillColorAlpha (col, falpha);
  ((TGAE*) g->Clone ())->Draw (opt);
  return;
} 



/**
 * Draws a graph as a systematic with filled box errors.
 */
void myDrawSystFill (TGAE* g, const Color_t col, const float falpha, const Style_t fstyle) {
  TGAE* gp = (TGAE*) g->Clone ();
  gp->SetFillColorAlpha (col, falpha);
  gp->Draw ("2");
  return;
}



/**
 * Takes a TGAE and fills two other TGAEs offset by one up error and one down error.
 */
void MakeGupAndGdown (TGAE* g, TGAE* gu, TGAE* gd) {
  double x, y;
  for (int i = 0; i < g->GetN (); i++) {
    g->GetPoint (i, x, y);
    gu->SetPoint (gu->GetN (), x, y + g->GetErrorYhigh (i));
    gd->SetPoint (gd->GetN (), x, y - g->GetErrorYlow (i));
  }
  return;
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
    h->SetBinContent (ix, h->GetBinContent (ix) * std::pow (a->GetBinContent (ix), power));
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
    n->SetBinError (ix, std::sqrt (std::fabs (var)));
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




void myBoxText2 (double x, double y, int color, int mstyle, const char* text, float msize, double tsize, bool doOutline) {
  TLatex l;
  l.SetTextAlign (12);
  l.SetTextSize (tsize); 
  l.SetNDC ();
  l.DrawLatex (x, y, text);

  double y1 = y - 0.5*tsize;
  double y2 = y + 0.5*tsize;
  double x1 = x - 0.061;
  double x2 = x - 0.013;

  TPave* mbox= new TPave (x1, y1, x2, y2, 1, "NDC");
  mbox->SetLineColor (color);
  mbox->SetLineWidth (1);
  mbox->SetFillStyle (0);
  mbox->Draw ();

  TLine* markerLine = new TLine ();
  markerLine->SetLineColor (color);
  markerLine->SetLineStyle (1);
  markerLine->SetLineWidth (2);
  markerLine->DrawLineNDC (0.5*(x1+x2), y1, 0.5*(x1+x2), y2);

  TMarker* marker = new TMarker (0.5*(x1+x2), 0.5*(y1+y2), 8);
  marker->SetMarkerColor (color);
  marker->SetNDC ();
  marker->SetMarkerStyle (mstyle);
  marker->SetMarkerSize (msize);
  marker->Draw ();

  if (doOutline && IsFullMarker (mstyle)) {
    //TMarker* marker2 = new TMarker (x-(0.8*tsize), y, 8);
    TMarker* marker2 = new TMarker (0.5*(x1+x2), 0.5*(y1+y2), 8);
    marker2->SetMarkerColor (kBlack);
    marker2->SetNDC ();
    marker2->SetMarkerStyle (FullToOpenMarker (mstyle));
    marker2->SetMarkerSize (msize);
    marker2->Draw ();
  }
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




void myLineText2 (double x, double y, int color, int mstyle, const char* text, float msize, double tsize, bool doOutline) {
  TLatex l;
  l.SetTextAlign (12);
  l.SetTextSize (tsize); 
  l.SetNDC ();
  l.DrawLatex (x, y, text);

  double y1 = y - 0.5*tsize;
  double y2 = y + 0.5*tsize;
  double x1 = x - 0.061;
  double x2 = x - 0.013;

  TLine* markerLine = new TLine ();
  markerLine->SetLineColor (color);
  markerLine->SetLineStyle (1);
  markerLine->SetLineWidth (2);
  markerLine->DrawLineNDC (0.5*(x1+x2), y1, 0.5*(x1+x2), y2);

  TMarker* marker = new TMarker (0.5*(x1+x2), 0.5*(y1+y2), 8);
  marker->SetMarkerColor (color);
  marker->SetNDC ();
  marker->SetMarkerStyle (mstyle);
  marker->SetMarkerSize (msize);
  marker->Draw ();

  if (doOutline && IsFullMarker (mstyle)) {
    //TMarker* marker2 = new TMarker (x-(0.8*tsize), y, 8);
    TMarker* marker2 = new TMarker (0.5*(x1+x2), 0.5*(y1+y2), 8);
    marker2->SetMarkerColor (kBlack);
    marker2->SetNDC ();
    marker2->SetMarkerStyle (FullToOpenMarker (mstyle));
    marker2->SetMarkerSize (msize);
    marker2->Draw ();
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




void mySimpleMarkerAndBoxAndLineText (double x, double y, const double bsize, const int bstyle, const int bcolor, const double balpha, const int mcolor, const int mstyle, const double msize, const char* text, const double tsize, const int lstyle) {

  const double y1 = y - (0.25*tsize) - (0.004*bsize) + 0.25*tsize;
  const double y2 = y + (0.25*tsize) + (0.004*bsize) + 0.25*tsize;
  const double x2 = x - (0.8*tsize) + 0.01;
  const double x1 = x - (0.8*tsize) + 0.01 - (0.04*bsize);

  TPave *mbox= new TPave (x1, y1, x2, y2, 0, "NDC");
  mbox->SetFillColorAlpha (bcolor, balpha);
  mbox->SetFillStyle (bstyle);
  mbox->Draw ();

  TLine *markerLine = new TLine ();
  markerLine->SetNDC ();
  markerLine->SetLineColor (mcolor);
  markerLine->SetLineStyle (lstyle);
  markerLine->SetLineWidth (3);
  markerLine->DrawLineNDC (x1, 0.5*(y1+y2), x2, 0.5*(y1+y2));

  if (msize > 0) {
    TMarker *marker = new TMarker(x2-0.02*bsize, y+0.25*tsize, 8);
    marker->SetNDC();
    marker->SetMarkerColor (mcolor);
    marker->SetMarkerStyle (mstyle);
    marker->SetMarkerSize (msize);
    marker->Draw ();
  }

 
  TLatex l;
  l.SetTextAlign (11);
  l.SetTextSize (tsize);
  l.SetNDC ();
  l.DrawLatex (x, y, text);
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
