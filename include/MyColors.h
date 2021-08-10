#ifndef __MyColors_h__
#define __MyColors_h__

#include <TColor.h>

const Color_t myBlue = (Color_t) TColor::GetColor (45, 64, 245);
const Color_t myPurple = (Color_t) TColor::GetColor (130,  10, 130);
const Color_t myRed = (Color_t) TColor::GetColor (255,  12,  73);
const Color_t myGreen = (Color_t) TColor::GetColor ( 54, 167,  80);
const Color_t myOrange = (Color_t) TColor::GetColor (255,  68,   0);
const Color_t myLiteRed = (Color_t) TColor::GetColor (255, 153, 153);
const Color_t myCyan = (Color_t) TColor::GetColor (102, 204, 204);
const Color_t myLitePurple = (Color_t) TColor::GetColor (204, 102, 255);
const Color_t myPeaGreen = (Color_t) TColor::GetColor (153, 204, 51);
//const Color_t myLiteYellow = (Color_t) TColor::GetColor (255, 255, 102); // too light
//const Color_t myLiteYellow = (Color_t) TColor::GetColor (220, 220, 20);
const Color_t myLiteYellow = kOrange-2;
const Color_t myLiteGreen = (Color_t) TColor::GetColor (102, 255, 153);
const Color_t myLiteBlue = (Color_t) TColor::GetColor (0, 153, 255);
const Color_t myViolet = (Color_t) TColor::GetColor (102, 0, 153);
const Color_t myMaroon = (Color_t) TColor::GetColor (204, 51, 51);

const Color_t colors[] = {myBlue, myGreen, myOrange, myRed, myPurple};
const int nColors = sizeof (colors) / sizeof (colors[0]);
const Color_t pastels[] = {myLiteRed, myLitePurple, myCyan, myLiteGreen, myPeaGreen, myLiteYellow};
const int nPastels = sizeof (pastels) / sizeof (pastels[0]);
const Color_t systColors[] = {kRed+1, kAzure-2, kGreen+2, kViolet-3, kMagenta, kCyan+1, kOrange-3, kGreen-7, kBlue+1, kPink-5};
const int nSystColors = sizeof (systColors) / sizeof (systColors[0]);
const Color_t manyColors[] = {kRed+1, kPink-3, kPink+9, kPink+6, kMagenta-3, kViolet-3, kViolet+7, kBlue+1, kAzure-2, kAzure+7, kAzure+6, kCyan-7, kTeal, kTeal-5, kGreen-3, kSpring+7, kSpring+5, kOrange-5, kOrange-2, kOrange-3, kOrange+6, kOrange+10};
const int nManyColors = sizeof (manyColors) / sizeof (manyColors[0]);

const Color_t lauraSystColors[] = {(Color_t) TColor::GetColor (102, 102, 204), (Color_t) TColor::GetColor (255, 102, 255), (Color_t) TColor::GetColor (255, 102, 102), (Color_t) TColor::GetColor (255, 204, 0), (Color_t) TColor::GetColor (204, 204, 51), (Color_t) TColor::GetColor (51, 204, 51), (Color_t) TColor::GetColor (102, 255, 255), (Color_t) TColor::GetColor (255, 0, 102)};
const Color_t lauraColors[] = {(Color_t) TColor::GetColor (0, 0, 153), (Color_t) TColor::GetColor (204, 0, 204), (Color_t) TColor::GetColor (204, 0, 0), (Color_t) TColor::GetColor (255, 102, 0), (Color_t) TColor::GetColor (153, 153, 0), (Color_t) TColor::GetColor (0, 153, 0), (Color_t) TColor::GetColor (1, 204, 204), (Color_t) TColor::GetColor (255, 102, 153)};


#endif
