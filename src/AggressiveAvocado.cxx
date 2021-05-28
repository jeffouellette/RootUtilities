#ifndef __AggressiveAvocado_cxx__
#define __AggressiveAvocado_cxx__

#include <iostream>

#include "../include/AggressiveAvocado.h"

#include "TROOT.h"


void SetAggressiveAvocadoStyle () {
  static TStyle* aggAvocadoStyle = 0;
  std::cout << "\nReluctantly applying esoteric avocodo-derived plot style settings...\n" << std::endl ;
  if ( aggAvocadoStyle==0 ) aggAvocadoStyle = AggressiveAvocadoStyle();
  gROOT->SetStyle("AggressiveAvocado");
  gROOT->ForceStyle();
}




TStyle* AggressiveAvocadoStyle() {
  TStyle* aggAvocadoStyle = new TStyle("AggressiveAvocado","Aggressive avocado style");

  aggAvocadoStyle->SetPalette (kAvocado);

  // use plain black on white colors
  Int_t icol=TColor::GetColorPalette (0); // 
  aggAvocadoStyle->SetFrameBorderMode(kWhite);
  aggAvocadoStyle->SetFrameFillColor(icol);
  aggAvocadoStyle->SetCanvasBorderMode(kWhite);
  aggAvocadoStyle->SetCanvasColor(icol);
  aggAvocadoStyle->SetCanvasDefW (1000);
  aggAvocadoStyle->SetCanvasDefH (812);
  aggAvocadoStyle->SetPadBorderMode(kWhite);
  aggAvocadoStyle->SetPadColor(kWhite);
  aggAvocadoStyle->SetStatColor(kWhite);
  aggAvocadoStyle->SetAxisColor(kWhite,"x");
  aggAvocadoStyle->SetAxisColor(kWhite,"y");
  aggAvocadoStyle->SetAxisColor(kWhite,"z");
  //aggAvocadoStyle->SetFillColor(icol); // don't use: white fill color for *all* objects

  // set the paper & margin sizes
  aggAvocadoStyle->SetPaperSize(20,26);

  // set margin sizes
  aggAvocadoStyle->SetPadTopMargin(0.05);
  aggAvocadoStyle->SetPadRightMargin(0.05);
  aggAvocadoStyle->SetPadBottomMargin(0.16);
  aggAvocadoStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  aggAvocadoStyle->SetTitleXOffset(1.4);
  aggAvocadoStyle->SetTitleYOffset(1.4);

  // use large fonts
  Int_t tcol=kWhite;
  //Int_t font=72; // Helvetica italics
  Int_t font=42; // Helvetica
  Double_t tsize=0.05;
  aggAvocadoStyle->SetTextFont(font);

  aggAvocadoStyle->SetTextSize(tsize);
  aggAvocadoStyle->SetLabelColor(tcol,"x");
  aggAvocadoStyle->SetTitleColor(tcol,"x");
  aggAvocadoStyle->SetLabelColor(tcol,"y");
  aggAvocadoStyle->SetTitleColor(tcol,"y");
  aggAvocadoStyle->SetLabelColor(tcol,"z");
  aggAvocadoStyle->SetTitleColor(tcol,"z");

  aggAvocadoStyle->SetLabelFont(font,"x");
  aggAvocadoStyle->SetTitleFont(font,"x");
  aggAvocadoStyle->SetLabelFont(font,"y");
  aggAvocadoStyle->SetTitleFont(font,"y");
  aggAvocadoStyle->SetLabelFont(font,"z");
  aggAvocadoStyle->SetTitleFont(font,"z");
  
  aggAvocadoStyle->SetLabelSize(tsize,"x");
  aggAvocadoStyle->SetTitleSize(tsize,"x");
  aggAvocadoStyle->SetLabelSize(tsize,"y");
  aggAvocadoStyle->SetTitleSize(tsize,"y");
  aggAvocadoStyle->SetLabelSize(tsize,"z");
  aggAvocadoStyle->SetTitleSize(tsize,"z");

  // use bold lines and markers
  aggAvocadoStyle->SetMarkerStyle(20);
  aggAvocadoStyle->SetMarkerSize(1.2);
  aggAvocadoStyle->SetHistLineWidth(2.);
  aggAvocadoStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars 
  //aggAvocadoStyle->SetErrorX(0.001);
  // get rid of error bar caps
  aggAvocadoStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  aggAvocadoStyle->SetOptTitle(0);
  //aggAvocadoStyle->SetOptStat(1111);
  aggAvocadoStyle->SetOptStat(0);
  //aggAvocadoStyle->SetOptFit(1111);
  aggAvocadoStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  aggAvocadoStyle->SetPadTickX(1);
  aggAvocadoStyle->SetPadTickY(1);

  return aggAvocadoStyle;

}


#endif // __AggressiveAvocado_cxx__
