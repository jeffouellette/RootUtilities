#ifndef __MyColors_h__
#define __MyColors_h__

#include <TColor.h>

const Color_t myBlue = (Color_t) TColor::GetColor (45, 64, 245);
const Color_t myPurple = (Color_t) TColor::GetColor (130,  10, 130);
const Color_t myRed = (Color_t) TColor::GetColor (255,  12,  73);
const Color_t myGreen = (Color_t) TColor::GetColor ( 54, 167,  80);
const Color_t myOrange = (Color_t) TColor::GetColor (255,  68,   0);

const Color_t colors[] = {myBlue, myGreen, myOrange, myRed, myPurple};
const int nColors = sizeof (colors) / sizeof (colors[0]);
const Color_t systColors[] = {kRed+1, kAzure-2, kGreen+2, kViolet-3, kMagenta, kCyan+1, kOrange-3, kGreen-7, kBlue+1, kPink-5};
const int nSystColors = sizeof (systColors) / sizeof (systColors[0]);


#endif
