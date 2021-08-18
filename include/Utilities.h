#ifndef __Utilities_h__
#define __Utilities_h__

#include <TFile.h>
#include <TString.h>
#include <TProfile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>
#include <TCanvas.h>

#include <string>

using namespace std;

typedef TGraphAsymmErrors TGAE;

/**
 * Safely deletes an object so that its pointer is also deleted.
 */
template <typename T> inline void SaferDelete (T** t) {
  if (*t) { delete (*t); (*t) = nullptr; }
}

/**
 * Returns a TString summarizing a measurement.
 * By default, there is only 1 significant digit in the error.
 * E.g., FormatMeasurement (40.58, 1.29) returns "40#pm1".
 * Or, FormatMeasurement (40.58, 1.29, 2) returns "40.6#pm1.3".
 */
const char* FormatMeasurement (double val, double err, const int n=1);


/**
 * Truncates numbers to the desired precision based on the most precise uncertainty (stat. or syst.)
 */
void FormatMeasurement (string& s_val, string& s_stat, string& s_syst, const int n);


/**
 * Sets the errors in h assuming the variances are stored in the entries of h2.
 */
void SetVariances (TH1D* h, TH2D* h2);


/**
 * Divides two histograms assuming the entries in one are binomial samples of the other.
 */
void BinomialDivide (TH1* out, TH1* num, TH1* den, TH1* h_sumwgt2 = nullptr);


/**
 * Modifies the directory strings to point to the correct locations.
 */
void SetupDirectories (const TString dataSubDir, const TString thisWorkPath);


/**
 * Clears sub-directory information from the directory strings
 */
void ResetDirectories ();


/**
 * Returns a linearly spaced array. The 0th element is lo, and the num-th element is hi.
 */
double* linspace (double lo, double hi, int num);


/**
 * Returns a logarithmically spaced array, where the 0th element is lo and the num-th element is hi.
 */
double* logspace (double lo, double hi, int num);


/**
 * Returns the equivalent angle in the range 0 to 2pi.
 */
double InTwoPi (double phi);


/**
 * Returns the difference between two angles in 0 to pi.
 * If sign is true, will return a signed version such that phi2 = phi1 + dphi
 */
double DeltaPhi (double phi1, double phi2, const bool sign = false);


/**
 * Returns dR between two eta, phi coordinates.
 */
double DeltaR (const double eta1, const double eta2, const double phi1, const double phi2 );


/**
 * Sets all the errors in this histogram to 0.
 */
void ResetHistErrors (TH1D* h);


/**
 * Sets all the errors in this TGAE to 0.
 */
void ResetTGAEErrors (TGAE* g);


/**
 * Sets all the x errors in this TGAE to 0.
 */
void ResetXErrors (TGAE* tg);


/**
 * Calculates errors in a TH1D from a TH2D storing the quadrature sum of each entry.
 */
void CalcUncertainties (TH1D* h, TH2D* h2, const double nEvts);


/**
 * Calculates errors in a TH1D from a TH2D storing the quadrature sum of each entry. Assumes weights are used, so a histogram with the number of events to the powers 0, 1, and 2 is needed.
 */
void CalcUncertainties (TH1D* h, TH2D* h2, TH1D* hn);


/**
 * Adds the maximum systematic uncertainty from v1 and v2 to sys in quadrature.
 */
void AddMaxSystematic (TGAE* sys, TGAE* v1, TGAE* v2);


/**
 * Adds nSigma statistical error variations to this histogram
 */
void AddStatVar (TH1D* h, const bool upvar = true, const float nSigma = 1);


/**
 * Adds independent systematic errors in quadrature, storing the sum in master
 */
void AddErrorsInQuadrature (TH1D* master, TH1D* sys);


/**
 * Adds independent systematic errors in quadrature, storing the sum in master
 */
void AddErrorsInQuadrature (TGAE* master, TH1D* sys);


/**
 * Adds independent systematic errors in quadrature, storing the sum in master
 */
void AddErrorsInQuadrature (TGAE* master, TGAE* sys, const bool doXErrs = false, const bool symmetrize = false);


/**
 * Adds independent systematic errors in quadrature, storing the sum in master
 * Allows user to send an array of systematics, not necessarily in order, and only adds the maximum uncertainty from these.
 */
void AddErrorsInQuadrature (TGAE* master, TGAE** sys, const std::vector <int>* indices, const bool symmetrize = false);


/**
 * Adds independent relative systematic errors in quadrature, storing the sum in master
 */
void AddRelErrorsInQuadrature (TGAE* master, TGAE* sys, const bool ignoreWarning = false, const bool symmetrize = false);


/**
 * Calculates simple systematics as maximum variations on the nominal.
 * Intended for combining up/down variations in an expandable way.
 */
void CalcSystematics (TH1D* sys, TH1D* var);


/**
 * Calculates simple systematics as maximum variations on the nominal.
 * Intended for combining up/down variations in an expandable way.
 */
void CalcSystematics (TGAE* sys, TH1D* var, const bool applyBothWays = true);


/**
 * Calculates simple systematics as maximum variations on the nominal.
 * Intended for combining up/down variations in an expandable way.
 */
void CalcSystematics (TGAE* sys, TH1D* nom, TH1D* var);


/**
 * Calculates the systematic errors on optimal, storing the results in graph.
 */
void CalcSystematics (TGAE* graph, TH1* optimal, const TH1* sys_hi, const TH1* sys_lo);


/**
 * Calculates the systematic errors on optimal, storing the results in graph.
 */
void CalcSystematics (TGAE* graph, TGAE* optimal, const TGraph* sys_hi, const TGraph* sys_lo, const bool doXErrs = false);


/**
 * Sets the bin contents in target as the errors in errors / central values in centralValues
 */
void SaveRelativeErrors (TGAE* target, TGAE* errors, TGAE* centralValues, const float sf = 1);


/**
 * Sets the bin contents in target as the error / central values in centralValues
 */
void SaveRelativeErrors (TH1D* errors, TH1D* centralValues, bool useCentVals=false);


/**
 * Sets the bin contents in highs and lows as the respective errors in centralValues
 */
void SaveAbsoluteErrors (TGAE* errors, TGAE* centralValues, TH1D* highs, TH1D* lows);


/**
 * Sets the bin contents in highs and lows as the respective errors / central values in centralValues
 */
void SaveRelativeErrors (TGAE* errors, TGAE* centralValues, TH1D* highs, TH1D* lows);


/**
 * Returns a TGAE with central values equal to the uncertainties in g.
 */
TGAE* ConvertErrorsToCentralValues (TGAE* g, const bool doHighErrs = true, const float sf = 1);


/**
 * Creates a projection of a TH3 with the specified axes by integrating between min and max on the 3rd axis of the TH3.
 */
TH2D* Project2D (TString name, TH3D* h3, const TString xaxis, const TString yaxis, const int min, const int max, const bool exclusive = false);


/**
 * Separates each point on a TGAE by delta along the x axis, so that the errors don't overlap.
 */
void deltaize (TGAE* tg, const double delta = 0, const bool logx = false, const double x1 = 0, const double x2 = 0);


/**
 * Offsets each point on a TGAE by delta along the y axis.
 */
void OffsetYAxis (TGAE* g, const double delta, const bool logx);


/**
 * Sets each point error to some constant value.
 * If mult is true, then will be multiplicative (intended for a log scale).
 */
void SetConstantXErrors (TGAE* tg, const double err, const bool mult, const double x1 = 0, const double x2 = 0);


/**
 * Makes a TGAE from the input histogram.
 */
TGAE* make_graph (TH1* h, const float cutoff = -1);


/**
 * Converts a TEfficiency to a TGAE
 */
TGAE* TEff2TGAE (TEfficiency* e);


/**
 * Recenters a TGraph according to the midpoints of matched
 */
void RecenterGraph (TGraph* g, TGraph* matched);


/**
 * Recenters a TGAE point for a log scale.
 */
void RecenterGraph (TGAE* g);


/**
 * Scales a TGAE and its errors by the central values in a TH1D.
 */
void ScaleGraph (TGAE* g, const TH1D* h = nullptr, const double sf = 1.);


/**
 * Applies new binning to a histogram
 * BE CAREFUL: if bins edges don't overlap, this can lead to unexpected behavior!
 */
void RebinSomeBins (TH1D** _h, int nbins, double* bins, const bool doWidths = false);


/**
 * Draws histogram as a graph with some plotting settings.
 */
void myDraw (TH1D* h, const Color_t col, const Style_t mstyle, const float msize, const Style_t lstyle = 1, const int lwidth = 2, const bool doMOutline = true);


/**
 * Draws histogram as a graph with some plotting settings.
 */
void myDrawHist (TH1D* h, const Color_t col, const Style_t lstyle = 1, const int lwidth = 2);


/**
 * Draws a graph with some plotting settings.
 */
void myDraw (TGAE* g, const Color_t col, const Style_t mstyle, const float msize, const Style_t lstyle = 1, const int lwidth = 2, const char* opt = "P", const bool doMOutline = true);


/**
 * Draws a function with some plotting settings.
 */
void myDraw (TF1* f, const Color_t col, const Style_t lstyle = 1, const int lwidth = 2, const char* opt = "same");


/**
 * Draws a filled area between gup and gdown with the given settings.
 */
void myDrawFill (TGAE* gup, TGAE* gdown, const Color_t col, const float falpha, const Style_t fstyle = 1001);


/**
 * Draws a graph as a systematic with some plotting settings.
 */
void myDrawSyst (TGAE* g, const Color_t col, const Style_t lstyle = 1, const int lwidth = 1, const float falpha = 0., const char* opt = "5");


/**
 * Draws a graph as a systematic with filled box errors.
 */
void myDrawSystFill (TGAE* g, const Color_t col, const float falpha, const Style_t fstyle);


/**
 * Adds a to h without propagating errors (e.g. for subtracting a background)
 */
void AddNoErrors (TH1D* h, TH1D* a, const float sf = 1.);


/**
 * Multiplies h by a raised to power without propagating errors (e.g. for subtracting a background)
 */
void MultiplyNoErrors (TH1D* h, TH1D* a, const float power = 1.);

bool IsOpenMarker (const Style_t ms);

bool IsFullMarker (const Style_t ms);

Style_t FullToOpenMarker (const Style_t ms);

Style_t OpenToFullMarker (const Style_t ms);

void FormatTH2Canvas (TCanvas* c, const bool zAxisSpace=true);

TGraphAsymmErrors* TH1TOTGraph (TH1 *h1);

void BinomialDivide (TH1* n, TH1* d);

/**
 * The infamous myText draws a text label on your favorite plot
 */
void myText (double x, double y,  Color_t color, const char *text, double tsize=0.04);

void myBoxText (double x, double y, double boxsize, int mcolor, const char *text, const double tsize=0.06);

void myBoxText2 (double x, double y, int color, int mstyle, const char *text, float msize=1.25, double tsize=0.032, bool doOutline=false);

void myBoxTextNoLine (double x, double y, double boxsize, int mcolor, const char *text, const double tsize=0.06);

void myLineText (double x, double y, int color, int lstyle, const char *text, float lsize, double tsize);

void myLineText2 (double x, double y, int color, int mstyle, const char *text, float msize=1.25, double tsize=0.032, bool doOutline=false);

void myLineColorText (double x, double y, int color, int lstyle, const char *text, float lsize, double tsize);

void myMarkerText (double x, double y, int color, int mstyle, const char *text, float msize=1.25, double tsize=0.032, bool doOutline=false);

void myMarkerTextNoLine (double x, double y, int color, int mstyle, const char *text, float msize=1.25, double tsize=0.032);

void myOnlyBoxText (double x, double y, double boxsize, int mcolor, int lcolor, int lstyle, const char *text, double tsize, int bstyle, double balpha);

TBox* TBoxNDC (const double x1, const double y1, const double x2, const double y2);

void mySimpleMarkerAndBoxAndLineText (double x, double y, const double bsize, const int bstyle, const int bcolor, const double balpha, const int mcolor, const int mstyle, const double msize, const char* text, const double tsize=0.032);

void myMarkerAndBoxAndLineText (double x, double y, const double bsize, const int bstyle, const int bcolor, const double balpha, const int mcolor, const int mstyle, const double msize, const char* text, const double tsize=0.032);

#endif
