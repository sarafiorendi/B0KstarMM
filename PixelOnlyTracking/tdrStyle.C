// ##############
// # tdrStyle.C #
// ##############


#include "TStyle.h"


// Turns the grid lines on (true) or off (false)
void tdrGrid(bool gridOn)
{
  tdrStyle->SetPadGridX(gridOn);
  tdrStyle->SetPadGridY(gridOn);
}


void tdrHerrors(bool HerrorsOn)
{
  if (HerrorsOn == true)
    tdrStyle->SetErrorX(0.5);
  else if (HerrorsOn == false)
    tdrStyle->SetErrorX(0.);
}


void tdrMyStyle()
{
  tdrStyle->SetOptStat(111110);
  tdrStyle->SetStatX(0.9);
  tdrStyle->SetStatY(0.9);
  tdrStyle->SetStatW(0.2);
  tdrStyle->SetStatH(0.1);

  tdrStyle->SetTitleYOffset(1.3);
  tdrStyle->SetNdivisions(509, "XYZ");
  tdrStyle->SetCanvasDefW(700);
  tdrStyle->SetPadRightMargin(0.07);
}


void tdrMyAxesSize()
{
  tdrStyle->SetLabelSize(0.04, "XYZ");
}


void CMS(double intLumi)
{
  TLatex latex;

  latex.SetNDC();
  latex.SetTextSize(0.04);
  
  latex.SetTextAlign(31); // align right
  latex.DrawLatex(0.98,0.965,"#font[22]{CMS preliminary 2010}, #sqrt{s} = 7 TeV");

  // if (intLumi > 0.)
  //   {
  //     latex.SetTextAlign(31); // align right
  //     latex.DrawLatex(0.98,0.88,Form("#int #font[12]{L}dt = %.1f nb^{-1}",intLumi));
  //   }

  // latex.SetTextAlign(11); // align left
  // latex.DrawLatex(0.02,0.965,"#font[22]{CMS preliminary 2010}");
}


void setTDRStyle()
{
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
  
  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); // Height of canvas
  tdrStyle->SetCanvasDefW(600); // Width of canvas
  tdrStyle->SetCanvasDefX(0);   // Position on screen
  tdrStyle->SetCanvasDefY(0);

  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  tdrStyle->SetPadBorderSize(1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  // For the histo:
//   tdrStyle->SetHistFillColor(1);
//   tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(0.5);
  // tdrStyle->SetNumberContours(20);

  tdrStyle->SetEndErrorSize(2);
  tdrStyle->SetErrorX(0.);

  // tdrStyle->SetMarkerStyle(20);

  // For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  // For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(0.01);
  // tdrStyle->SetDateY(0.01);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS: SetOptStat("mr")
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  tdrStyle->SetStatStyle(1001);
  tdrStyle->SetStatX(0.);
  tdrStyle->SetStatY(0.);

  // Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

  // For the global title:
  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  tdrStyle->SetTitleH(0);     // Set the height of the title box
  tdrStyle->SetTitleW(0);     // Set the width of the title box
  tdrStyle->SetTitleX(0);     // Set the position of the title box
  tdrStyle->SetTitleY(0.985); // Set the position of the title box
  tdrStyle->SetTitleStyle(1001);
  tdrStyle->SetTitleBorderSize(2);

  // For the axis titles:
  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  tdrStyle->SetTitleOffset(0.9, "X");
  tdrStyle->SetTitleOffset(1.25, "Y");

  // For the axis labels:
  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:
  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  // Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(3.);
  // tdrStyle->SetLineStyleString(int i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(0.5);
  // tdrStyle->SetBarWidth(0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(int ncolors = 0, int* colors = 0);
  // tdrStyle->SetTimeOffset(double toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();
}
