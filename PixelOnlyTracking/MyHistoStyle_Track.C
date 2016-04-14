// #########################
// # MyHistoStyle_Tracks.C #
// #########################

{
// ##################
// # STYLE SETTINGS #
// ##################
  gROOT->Reset();

  setTDRStyle();
  gROOT->SetStyle("tdrStyle");
  tdrMyAxesSize();
  tdrMyStyle();

  // gROOT->SetStyle("Plain");
  // gStyle->SetOptStat(111110);
  // gStyle->SetStatX(0.9);
  // gStyle->SetStatY(0.9);
  // gStyle->SetStatW(0.2);
  // gStyle->SetStatH(0.1);
  // gStyle->SetOptFit(1);
  // gStyle->SetOptTitle(false);

  gROOT->ForceStyle();

// ################
// # FILE READING #
// ################
  TFile* file0 = TFile::Open("MyPixAnalysisData_Tk_Merged.root");
  TFile* file1 = TFile::Open("MyPixAnalysisMC_noTP_Tk_Merged.root");
  TFile* file2 = TFile::Open("MyPixAnalysisMC_Tk_Merged.root");

// ###################
// # Histo Variables #
// ###################
  TH1F* histo0;
  TH1F* histo1;
  TH1F* histo2;
  TH1F* histo3;
  TH1F* histo4;

  TH2F* histo2D0;
  TH2F* histo2D1;
  TH2F* histo2D2;

  TH1F* histo1D0;
  TH1F* histo1D1;
  TH1F* histo1D2;

  TPaveStats* PS;

  TCanvas* c[200];
  stringstream CanvasName;
  unsigned int itC = -1;

  double MaxEtaTkTrk = 2.3;
  double EtaThr = 1.5;

  double tmp0;
  double tmp1;
  double tmp2;

  bool WithMC = true;

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs eta");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs eta");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs eta");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());

  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.7,1.2);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Eff_eta.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs pt thr1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs pt thr1");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs pt thr1");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.7,1.2);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"|\\eta| < 1.5");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Eff_pt_etathr1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs pt thr2");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs pt thr2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs pt thr2");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.7,1.2);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"1.5 < |\\eta| < 2.3");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Eff_pt_etathr2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs pt");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs pt");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs pt");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.7,1.2);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Eff_pt.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs phi thr1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs phi thr1");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs phi thr1");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.5,1.4);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"-2.3 < \\eta < -1.5");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Eff_phi_etathr1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs phi thr2");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs phi thr2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs phi thr2");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.5,1.4);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"-1.5 < \\eta < 1.5");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Eff_phi_etathr2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs phi thr3");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs phi thr3");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/Efficiency pixeltracks vs phi thr3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.5,1.4);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"1.5 < \\eta < 2.3");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Eff_phi_etathr3.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs eta");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs eta");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs eta");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.7,1.2);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("AlgoEff_eta.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs pt thr1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs pt thr1");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs pt thr1");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.7,1.2);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"|\\eta| < 1.5");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("AlgoEff_pt_etathr1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs pt thr2");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs pt thr2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs pt thr2");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.7,1.2);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"1.5 < |\\eta| < 2.3");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("AlgoEff_pt_etathr2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs pt");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs pt");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs pt");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.7,1.2);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("AlgoEff_pt.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs phi thr1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs phi thr1");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs phi thr1");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.5,1.4);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"-2.3 < \\eta < -1.5");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("AlgoEff_phi_etathr1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs phi thr2");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs phi thr2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs phi thr2");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.5,1.4);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"-1.5 < \\eta < 1.5");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("AlgoEff_phi_etathr2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs phi thr3");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs phi thr3");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/Algoefficiency pixeltracks vs phi thr3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.5,1.4);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"1.5 < \\eta < 2.3");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("AlgoEff_phi_etathr3.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs eta");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs eta");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs eta");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.7,1.2);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Pur_eta.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs pt thr1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs pt thr1");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs pt thr1");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.7,1.2);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"|\\eta| < 1.5");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Pur_pt_etathr1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs pt thr2");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs pt thr2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs pt thr2");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.7,1.2);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"1.5 < |\\eta| < 2.3");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Pur_pt_etathr2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs pt");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs pt");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs pt");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.7,1.2);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Pur_pt.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs phi thr1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs phi thr1");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs phi thr1");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.5,1.4);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"-2.3 < \\eta < -1.5");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Pur_phi_etathr1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs phi thr2");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs phi thr2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs phi thr2");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.5,1.4);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"-1.5 < \\eta < 1.5");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Pur_phi_etathr2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs phi thr3");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs phi thr3");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity pixeltracks vs phi thr3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.5,1.4);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"1.5 < \\eta < 2.3");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Pur_phi_etathr3.png");

  tdrHerrors(true);

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt resolution vs pt thr1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt resolution vs pt thr1");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt resolution vs pt thr1");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.35);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"|\\eta| < 1.5");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("PtRes_etathr1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt resolution vs pt thr2");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt resolution vs pt thr2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt resolution vs pt thr2");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.35);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"1.5 < |\\eta| < 2.3");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("PtRes_etathr2.png");

  tdrHerrors(false);

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt resolution vs eta thr1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt resolution vs eta thr1");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt resolution vs eta thr1");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.2);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"0.75 < p_{t} < 1.25 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("PtRes_ptthr1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt resolution vs eta thr2");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt resolution vs eta thr2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt resolution vs eta thr2");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.2);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"1.25 < p_{t} < 1.75 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("PtRes_ptthr2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt resolution vs eta thr3");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt resolution vs eta thr3");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt resolution vs eta thr3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.2);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"1.75 < p_{t} < 2.25 G(eV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("PtRes_ptthr3.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt resolution vs eta thr4");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt resolution vs eta thr4");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt resolution vs eta thr4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.2);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"2.25 < p_{t} < 2.75 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("PtRes_ptthr4.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt resolution vs eta thr5");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt resolution vs eta thr5");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt resolution vs eta thr5");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.2);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"2.75 < p_{t} < 8.0 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("PtRes_ptthr5.png");

  tdrHerrors(true);

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt pulls vs pt thr1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt pulls vs pt thr1");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt pulls vs pt thr1");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,2.5);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"|\\eta| < 1.5");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("PtPulls_etathr1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt pulls vs pt thr2");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt pulls vs pt thr2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt pulls vs pt thr2");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,2.5);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"1.5 < |\\eta| < 2.3");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("PtPulls_etathr2.png");

  tdrHerrors(false);

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt pulls vs eta thr1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt pulls vs eta thr1");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt pulls vs eta thr1");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,2.5);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"0.75 < p_{t} < 1.25 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("PtPulls_ptthr1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt pulls vs eta thr2");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt pulls vs eta thr2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt pulls vs eta thr2");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,2.5);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"1.25 < p_{t} < 1.75 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("PtPulls_ptthr2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt pulls vs eta thr3");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt pulls vs eta thr3");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt pulls vs eta thr3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,2.5);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"1.75 < p_{t} < 2.25 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("PtPulls_ptthr3.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt pulls vs eta thr4");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt pulls vs eta thr4");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt pulls vs eta thr4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,2.5);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"2.25 < p_{t} < 2.75 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("PtPulls_ptthr4.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt pulls vs eta thr5");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt pulls vs eta thr5");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/pt pulls vs eta thr5");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,2.5);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"2.75 < p_{t} < 8.0 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("PtPulls_ptthr5.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/SkewPixTrkpt vs pt thr1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/SkewPixTrkpt vs pt thr1");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/SkewPixTrkpt vs pt thr1");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->GetYaxis()->SetRangeUser(-1.,25.);
  histo0->Draw();
  
  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetLineColor(kBlue);
      histo2->SetLineWidth(2.0);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"|\\eta| < 1.5");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw("");

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Skew_etathr1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/SkewPixTrkpt vs pt thr2");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/SkewPixTrkpt vs pt thr2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/SkewPixTrkpt vs pt thr2");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->GetYaxis()->SetRangeUser(-1.,25.);
  histo0->Draw();
  
  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetLineColor(kBlue);
      histo2->SetLineWidth(2.0);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"1.5 < |\\eta| < 2.3");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Skew_etathr2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 resolution vs pt thr1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 resolution vs pt thr1");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 resolution vs pt thr1");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.03);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"|\\eta| < 1.5");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("D0Res_etathr1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 resolution vs pt thr2");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 resolution vs pt thr2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 resolution vs pt thr2");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.03);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"1.5 < |\\eta| < 2.3");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("D0Res_etathr2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 resolution vs eta thr1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 resolution vs eta thr1");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 resolution vs eta thr1");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.04);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"0.75 < p_{t} < 1.25 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("D0Res_ptthr1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 resolution vs eta thr2");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 resolution vs eta thr2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 resolution vs eta thr2");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.04);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"1.25 < p_{t} < 1.75 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("D0Res_ptthr2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 resolution vs eta thr3");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 resolution vs eta thr3");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 resolution vs eta thr3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.04);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"1.75 < p_{t} < 2.25 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("D0Res_ptthr3.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 resolution vs eta thr4");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 resolution vs eta thr4");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 resolution vs eta thr4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.04);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"2.25 < p_{t} < 2.75 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("D0Res_ptthr4.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 resolution vs eta thr5");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 resolution vs eta thr5");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 resolution vs eta thr5");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.04);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"2.75 < p_{t} < 8.0 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("D0Res_ptthr5.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 pulls vs pt thr1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 pulls vs pt thr1");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 pulls vs pt thr1");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,2.5);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"|\\eta| < 1.5");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("D0Pulls_etathr1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 pulls vs pt thr2");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 pulls vs pt thr2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/d0 pulls vs pt thr2");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,2.5);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"1.5 < |\\eta| < 2.3");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("D0Pulls_etathr2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz resolution vs pt thr1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz resolution vs pt thr1");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz resolution vs pt thr1");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.03);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"|\\eta| < 1.5");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("DzRes_etathr1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz resolution vs pt thr2");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz resolution vs pt thr2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz resolution vs pt thr2");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.1);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"1.5 < |\\eta| < 2.3");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("DzRes_etathr2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz resolution vs eta thr1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz resolution vs eta thr1");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz resolution vs eta thr1");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.12);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"0.75 < p_{t} < 1.25 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("DzRes_ptthr1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz resolution vs eta thr2");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz resolution vs eta thr2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz resolution vs eta thr2");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.12);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"1.25 < p_{t} < 1.75 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("DzRes_ptthr2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz resolution vs eta thr3");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz resolution vs eta thr3");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz resolution vs eta thr3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.12);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"1.75 < p_{t} < 2.25 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("DzRes_ptthr3.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz resolution vs eta thr4");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz resolution vs eta thr4");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz resolution vs eta thr4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.12);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"2.25 < p_{t} < 2.75 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("DzRes_ptthr4.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz resolution vs eta thr5");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz resolution vs eta thr5");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz resolution vs eta thr5");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.12);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"2.75 < p_{t} < 8.0 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("DzRes_ptthr5.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz pulls vs eta thr1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz pulls vs eta thr1");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz pulls vs eta thr1");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,2.5);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"0.75 < p_{t} < 1.25 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("DzPulls_ptthr1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz pulls vs eta thr2");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz pulls vs eta thr2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz pulls vs eta thr2");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,2.5);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"1.25 < p_{t} < 1.75 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("DzPulls_ptthr2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz pulls vs eta thr3");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz pulls vs eta thr3");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz pulls vs eta thr3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,2.5);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"1.75 < p_{t} < 2.25 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("DzPulls_ptthr3.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz pulls vs eta thr4");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz pulls vs eta thr4");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz pulls vs eta thr4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,2.5);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"2.25 < p_{t} < 2.75 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("DzPulls_ptthr4.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz pulls vs eta thr5");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz pulls vs eta thr5");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz pulls vs eta thr5");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,2.5);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"2.75 < p_{t} < 8.0 (GeV/c)");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("DzPulls_ptthr5.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz pulls vs pt thr1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz pulls vs pt thr1");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz pulls vs pt thr1");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,2.5);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"|\\eta| < 1.5");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("DzPulls_etathr1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz pulls vs pt thr2");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz pulls vs pt thr2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/dz pulls vs pt thr2");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,2.5);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"1.5 < |\\eta| < 2.3");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("DzPulls_etathr2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity Efficiency Tracks");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity Efficiency Tracks");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/CommonPxTkPlots/Purity Efficiency Tracks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.7,1.2);
  histo0->Draw();
  
  cout << "@@@ Data @@@" << endl;
  cout << "Purity: " << histo0->GetBinContent(histo0->GetXaxis()->FindBin("Purity")) << "\t+/- " << histo0->GetBinError(histo0->GetXaxis()->FindBin("Purity")) << endl;
  cout << "Efficiency: " << histo0->GetBinContent(histo0->GetXaxis()->FindBin("Efficiency")) << "\t+/- " << histo0->GetBinError(histo0->GetXaxis()->FindBin("Efficiency")) << endl;
  cout << "Algo-Efficiency: " << histo0->GetBinContent(histo0->GetXaxis()->FindBin("Algo-Efficiency")) << "\t+/- " << histo0->GetBinError(histo0->GetXaxis()->FindBin("Algo-Efficiency")) << endl;

  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  cout << "@@@ Simulation @@@" << endl;
  cout << "Purity: " << histo1->GetBinContent(histo1->GetXaxis()->FindBin("Purity")) << " +/- " << histo1->GetBinError(histo1->GetXaxis()->FindBin("Purity")) << endl;
  cout << "Efficiency: " << histo1->GetBinContent(histo1->GetXaxis()->FindBin("Efficiency")) << " +/- " << histo1->GetBinError(histo1->GetXaxis()->FindBin("Efficiency")) << endl;
  cout << "Algo-Efficiency: " << histo1->GetBinContent(histo1->GetXaxis()->FindBin("Algo-Efficiency")) << " +/- " << histo1->GetBinError(histo1->GetXaxis()->FindBin("Algo-Efficiency")) << endl;

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  cout << "@@@ MC-Truth @@@" << endl;
  cout << "Purity: " << histo2->GetBinContent(histo2->GetXaxis()->FindBin("Purity")) << " +/- " << histo2->GetBinError(histo2->GetXaxis()->FindBin("Purity")) << endl;
  cout << "Efficiency: " << histo2->GetBinContent(histo2->GetXaxis()->FindBin("Efficiency")) << " +/- " << histo2->GetBinError(histo2->GetXaxis()->FindBin("Efficiency")) << endl;
  cout << "Algo-Efficiency: " << histo2->GetBinContent(histo2->GetXaxis()->FindBin("Algo-Efficiency")) << " +/- " << histo2->GetBinError(histo0->GetXaxis()->FindBin("Algo-Efficiency")) << endl;

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Pur_Eff_Trk.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/PixelTracks/Doublecounter pixeltracks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  c[itC]->SetLogy(1);

  histo0->SetFillColor(kBlue);
  histo0->Draw("bar");

  histo0->SetStats(true);

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_TrkPar_PixDouble_Count.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/PixelTracks/Doublecounter pixeltracks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  c[itC]->SetLogy(1);

  histo0->SetFillColor(kBlue);
  histo0->Draw("bar");
  
  histo0->SetStats(true);

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_Reco_PixDouble_Count.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/Doublecounter pixeltracks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  c[itC]->SetLogy(1);

  histo0->SetFillColor(kBlue);
  histo0->Draw("bar");
  
  histo0->SetStats(true);

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Data_PixDouble_Count.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/PixelTracks/Doublecounter pixeltracks eta distribution");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());

  histo0->SetFillColor(kBlue);
  histo0->Draw("bar");

  histo0->SetStats(true);

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_TrkPar_PixDouble_Eta.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/PixelTracks/Doublecounter pixeltracks eta distribution");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());

  histo0->SetFillColor(kBlue);
  histo0->Draw("bar");
  
  histo0->SetStats(true);

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_Reco_PixDouble_Eta.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/PixelTracks/Doublecounter pixeltracks eta distribution");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());

  histo0->SetFillColor(kBlue);
  histo0->Draw("bar");
  
  histo0->SetStats(true);

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Data_PixDouble_Eta.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/Doublecounter pixeltracks pt distribution");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());

  histo0->SetFillColor(kBlue);
  histo0->Draw("bar");
  
  histo0->SetStats(true);

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_TrkPar_PixDouble_Pt.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/PixelTracks/Doublecounter pixeltracks pt distribution");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());

  histo0->SetFillColor(kBlue);
  histo0->Draw("bar");
  
  histo0->SetStats(true);

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_Reco_PixDouble_Pt.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/PixelTracks/Doublecounter pixeltracks pt distribution");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());

  histo0->SetFillColor(kBlue);
  histo0->Draw("bar");
  
  histo0->SetStats(true);

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Data_PixDouble_Pt.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/Doublecounter pixeltracks eta relative distribution");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/Doublecounter pixeltracks eta relative distribution");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/Doublecounter pixeltracks eta relative distribution");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.08);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("PixDouble_EtaRel.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/Doublecounter pixeltracks pt relative distribution");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/Doublecounter pixeltracks pt relative distribution");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/Doublecounter pixeltracks pt relative distribution");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.025);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  histo0->SetStats(false);
  histo1->SetStats(false);
  if (WithMC == true) histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("PixDouble_PtRel.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/PixelTracks/Doublecounter trackertracks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  c[itC]->SetLogy(1);

  histo0->SetFillColor(kBlue);
  histo0->Draw("bar");
  
  histo0->SetStats(true);

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_Reco_TrkDouble_Count.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/Doublecounter trackertracks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  c[itC]->SetLogy(1);

  histo0->SetFillColor(kBlue);
  histo0->Draw("bar");
  
  histo0->SetStats(true);

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Data_TrkDouble_Count.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/PixelTracks/Doublecounter trackertracks eta distribution");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());

  histo0->SetFillColor(kBlue);
  histo0->Draw("bar");
  
  histo0->SetStats(true);

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_Reco_TrkDouble_Eta.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/Doublecounter trackertracks eta distribution");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());

  histo0->SetFillColor(kBlue);
  histo0->Draw("bar");
  
  histo0->SetStats(true);

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Data_TrkDouble_Eta.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/PixelTracks/Doublecounter trackertracks pt distribution");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());

  histo0->SetFillColor(kBlue);
  histo0->Draw("bar");
  
  histo0->SetStats(true);

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_Reco_TrkDouble_Pt.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/Doublecounter trackertracks pt distribution");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());

  histo0->SetFillColor(kBlue);
  histo0->Draw("bar");
  
  histo0->SetStats(true);

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Data_TrkDouble_Pt.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/Doublecounter trackertracks eta relative distribution");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/Doublecounter trackertracks eta relative distribution");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,0.0008);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("TrkDouble_EtaRel.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/Doublecounter trackertracks pt relative distribution");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/Doublecounter trackertracks pt relative distribution");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYax is()->SetRangeUser(0.,0.0007);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Data");
  leg->AddEntry(histo1,"Simulation");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("TrkDouble_PtRel.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/TrackFit/dz_pt_1_eta_0.1");
  histo1 = (TH1F*)file2->Get("MyProcess/TrackFit/dz_pt_1_eta_1.5");
  histo2 = (TH1F*)file2->Get("MyProcess/TrackFit/dz_pt_1_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"p_{t} = 1 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("MC_TrkPar_dzRes1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/TrackFit/dz_pt_6_eta_0.1");
  histo1 = (TH1F*)file2->Get("MyProcess/TrackFit/dz_pt_6_eta_1.5");
  histo2 = (TH1F*)file2->Get("MyProcess/TrackFit/dz_pt_6_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->GetYaxis()->SetRangeUser(0.,350.);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"p_{t} = 6 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("MC_TrkPar_dzRes2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/TrackFit/dz_pt_1_eta_0.1");
  histo1 = (TH1F*)file1->Get("MyProcess/TrackFit/dz_pt_1_eta_1.5");
  histo2 = (TH1F*)file1->Get("MyProcess/TrackFit/dz_pt_1_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"p_{t} = 1 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("MC_Reco_dzRes1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/TrackFit/dz_pt_6_eta_0.1");
  histo1 = (TH1F*)file1->Get("MyProcess/TrackFit/dz_pt_6_eta_1.5");
  histo2 = (TH1F*)file1->Get("MyProcess/TrackFit/dz_pt_6_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->GetYaxis()->SetRangeUser(0.,350.);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"p_{t} = 6 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("MC_Reco_dzRes2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/TrackFit/dz_pt_1_eta_0.1");
  histo1 = (TH1F*)file0->Get("MyProcess/TrackFit/dz_pt_1_eta_1.5");
  histo2 = (TH1F*)file0->Get("MyProcess/TrackFit/dz_pt_1_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"p_{t} = 1 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("Data_dzRes1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/TrackFit/dz_pt_6_eta_0.1");
  histo1 = (TH1F*)file0->Get("MyProcess/TrackFit/dz_pt_6_eta_1.5");
  histo2 = (TH1F*)file0->Get("MyProcess/TrackFit/dz_pt_6_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->GetYaxis()->SetRangeUser(0.,400.);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"p_{t} = 6 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("Data_dzRes2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/TrackFit/dzerr_pt_1_eta_0.1");
  histo1 = (TH1F*)file1->Get("MyProcess/TrackFit/dzerr_pt_1_eta_1.5");
  histo2 = (TH1F*)file1->Get("MyProcess/TrackFit/dzerr_pt_1_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.25,0.7,0.55,0.9,"p_{t} = 1 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("MC_Reco_dzErr1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/TrackFit/dzerr_pt_6_eta_0.1");
  histo1 = (TH1F*)file1->Get("MyProcess/TrackFit/dzerr_pt_6_eta_1.5");
  histo2 = (TH1F*)file1->Get("MyProcess/TrackFit/dzerr_pt_6_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.25,0.7,0.55,0.9,"p_{t} = 6 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("MC_Reco_dzErr2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/TrackFit/dzerr_pt_1_eta_0.1");
  histo1 = (TH1F*)file0->Get("MyProcess/TrackFit/dzerr_pt_1_eta_1.5");
  histo2 = (TH1F*)file0->Get("MyProcess/TrackFit/dzerr_pt_1_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.25,0.7,0.55,0.9,"p_{t} = 1 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("Data_dzErr1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/TrackFit/dzerr_pt_6_eta_0.1");
  histo1 = (TH1F*)file0->Get("MyProcess/TrackFit/dzerr_pt_6_eta_1.5");
  histo2 = (TH1F*)file0->Get("MyProcess/TrackFit/dzerr_pt_6_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.25,0.7,0.55,0.9,"p_{t} = 6 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("Data_dzErr2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/TrackFit/d0_pt_1_eta_0.1");
  histo1 = (TH1F*)file2->Get("MyProcess/TrackFit/d0_pt_1_eta_1.5");
  histo2 = (TH1F*)file2->Get("MyProcess/TrackFit/d0_pt_1_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"p_{t} = 1 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("MC_TrkPar_d0Res1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/TrackFit/d0_pt_6_eta_0.1");
  histo1 = (TH1F*)file2->Get("MyProcess/TrackFit/d0_pt_6_eta_1.5");
  histo2 = (TH1F*)file2->Get("MyProcess/TrackFit/d0_pt_6_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"p_{t} = 6 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("MC_TrkPar_d0Res2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/TrackFit/d0_pt_1_eta_0.1");
  histo1 = (TH1F*)file1->Get("MyProcess/TrackFit/d0_pt_1_eta_1.5");
  histo2 = (TH1F*)file1->Get("MyProcess/TrackFit/d0_pt_1_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"p_{t} = 1 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("MC_Reco_d0Res1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/TrackFit/d0_pt_6_eta_0.1");
  histo1 = (TH1F*)file1->Get("MyProcess/TrackFit/d0_pt_6_eta_1.5");
  histo2 = (TH1F*)file1->Get("MyProcess/TrackFit/d0_pt_6_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->GetYaxis()->SetRangeUser(0.,300.);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"p_{t} = 6 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("MC_Reco_d0Res2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/TrackFit/d0_pt_1_eta_0.1");
  histo1 = (TH1F*)file0->Get("MyProcess/TrackFit/d0_pt_1_eta_1.5");
  histo2 = (TH1F*)file0->Get("MyProcess/TrackFit/d0_pt_1_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"p_{t} = 1 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("Data_d0Res1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/TrackFit/d0_pt_6_eta_0.1");
  histo1 = (TH1F*)file0->Get("MyProcess/TrackFit/d0_pt_6_eta_1.5");
  histo2 = (TH1F*)file0->Get("MyProcess/TrackFit/d0_pt_6_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->GetYaxis()->SetRangeUser(0.,450.);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"p_{t} = 6 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("Data_d0Res2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/TrackFit/d0err_pt_1_eta_0.1");
  histo1 = (TH1F*)file1->Get("MyProcess/TrackFit/d0err_pt_1_eta_1.5");
  histo2 = (TH1F*)file1->Get("MyProcess/TrackFit/d0err_pt_1_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.25,0.7,0.55,0.9,"p_{t} = 1 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("MC_Reco_d0Err1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/TrackFit/d0err_pt_6_eta_0.1");
  histo1 = (TH1F*)file1->Get("MyProcess/TrackFit/d0err_pt_6_eta_1.5");
  histo2 = (TH1F*)file1->Get("MyProcess/TrackFit/d0err_pt_6_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.25,0.7,0.55,0.9,"p_{t} = 6 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("MC_Reco_d0Err2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/TrackFit/d0err_pt_1_eta_0.1");
  histo1 = (TH1F*)file0->Get("MyProcess/TrackFit/d0err_pt_1_eta_1.5");
  histo2 = (TH1F*)file0->Get("MyProcess/TrackFit/d0err_pt_1_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.25,0.7,0.55,0.9,"p_{t} = 1 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("Data_d0Err1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/TrackFit/d0err_pt_6_eta_0.1");
  histo1 = (TH1F*)file0->Get("MyProcess/TrackFit/d0err_pt_6_eta_1.5");
  histo2 = (TH1F*)file0->Get("MyProcess/TrackFit/d0err_pt_6_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.25,0.7,0.55,0.9,"p_{t} = 6 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("Data_d0Err2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/TrackFit/pt_pt_1_eta_0.1");
  histo1 = (TH1F*)file2->Get("MyProcess/TrackFit/pt_pt_1_eta_1.5");
  histo2 = (TH1F*)file2->Get("MyProcess/TrackFit/pt_pt_1_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"p_{t} = 1 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("MC_TrkPar_ptRes1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/TrackFit/pt_pt_6_eta_0.1");
  histo1 = (TH1F*)file2->Get("MyProcess/TrackFit/pt_pt_6_eta_1.5");
  histo2 = (TH1F*)file2->Get("MyProcess/TrackFit/pt_pt_6_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->GetYaxis()->SetRangeUser(0.,20.);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"p_{t} = 6 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("MC_TrkPar_ptRes2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/TrackFit/pt_pt_1_eta_0.1");
  histo1 = (TH1F*)file1->Get("MyProcess/TrackFit/pt_pt_1_eta_1.5");
  histo2 = (TH1F*)file1->Get("MyProcess/TrackFit/pt_pt_1_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"p_{t} = 1 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("MC_Reco_ptRes1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/TrackFit/pt_pt_6_eta_0.1");
  histo1 = (TH1F*)file1->Get("MyProcess/TrackFit/pt_pt_6_eta_1.5");
  histo2 = (TH1F*)file1->Get("MyProcess/TrackFit/pt_pt_6_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->GetYaxis()->SetRangeUser(0.,20.);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"p_{t} = 6 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("MC_Reco_ptRes2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/TrackFit/pt_pt_1_eta_0.1");
  histo1 = (TH1F*)file0->Get("MyProcess/TrackFit/pt_pt_1_eta_1.5");
  histo2 = (TH1F*)file0->Get("MyProcess/TrackFit/pt_pt_1_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
   histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"p_{t} = 1 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("Data_ptRes1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/TrackFit/pt_pt_6_eta_0.1");
  histo1 = (TH1F*)file0->Get("MyProcess/TrackFit/pt_pt_6_eta_1.5");
  histo2 = (TH1F*)file0->Get("MyProcess/TrackFit/pt_pt_6_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->GetYaxis()->SetRangeUser(0.,30.);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
   histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"p_{t} = 6 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("Data_ptRes2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/TrackFit/pterr_pt_1_eta_0.1");
  histo1 = (TH1F*)file1->Get("MyProcess/TrackFit/pterr_pt_1_eta_1.5");
  histo2 = (TH1F*)file1->Get("MyProcess/TrackFit/pterr_pt_1_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->GetXaxis()->SetRangeUser(0.,2000.);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.25,0.7,0.55,0.9,"p_{t} = 1 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("MC_Reco_ptErr1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/TrackFit/pterr_pt_6_eta_0.1");
  histo1 = (TH1F*)file1->Get("MyProcess/TrackFit/pterr_pt_6_eta_1.5");
  histo2 = (TH1F*)file1->Get("MyProcess/TrackFit/pterr_pt_6_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.25,0.7,0.55,0.9,"p_{t} = 6 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("MC_Reco_ptErr2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/TrackFit/pterr_pt_1_eta_0.1");
  histo1 = (TH1F*)file0->Get("MyProcess/TrackFit/pterr_pt_1_eta_1.5");
  histo2 = (TH1F*)file0->Get("MyProcess/TrackFit/pterr_pt_1_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->GetXaxis()->SetRangeUser(0.,2000.);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.25,0.7,0.55,0.9,"p_{t} = 1 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("Data_ptErr1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/TrackFit/pterr_pt_6_eta_0.1");
  histo1 = (TH1F*)file0->Get("MyProcess/TrackFit/pterr_pt_6_eta_1.5");
  histo2 = (TH1F*)file0->Get("MyProcess/TrackFit/pterr_pt_6_eta_2.3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");
    
  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.25,0.7,0.55,0.9,"p_{t} = 6 (GeV/c)");
  leg->AddEntry(histo0,"\\eta = 0.1");
  leg->AddEntry(histo1,"\\eta = 1.5");
  leg->AddEntry(histo2,"\\eta = 2.3");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kRed);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kBlue);
  c[itC]->SaveAs("Data_ptErr2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo2D0 = (TH2F*)file0->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi");
  histo2D1 = (TH2F*)file1->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi");
  histo2D2 = (TH2F*)file2->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo1D0 = new TH1F("H1D0","H1D0", histo2D0->GetNbinsX(), histo2D0->GetBinLowEdge(1), histo2D0->GetBinLowEdge(histo2D0->GetNbinsX()) + histo2D0->GetBinWidth(histo2D0->GetNbinsX()));
  histo1D1 = new TH1F("H1D1","H1D1", histo2D1->GetNbinsX(), histo2D1->GetBinLowEdge(1), histo2D1->GetBinLowEdge(histo2D1->GetNbinsX()) + histo2D1->GetBinWidth(histo2D1->GetNbinsX()));
  histo1D2 = new TH1F("H1D2","H1D2", histo2D2->GetNbinsX(), histo2D2->GetBinLowEdge(1), histo2D2->GetBinLowEdge(histo2D2->GetNbinsX()) + histo2D2->GetBinWidth(histo2D2->GetNbinsX()));
  
  for (int i = 1; i <= histo2D0->GetNbinsX(); i++)
    {
      for (int j = 1; j <= histo2D0->GetNbinsY(); j++)
	{
	  histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) + histo2D0->GetBinContent(i,j));
	  histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) + histo2D1->GetBinContent(i,j));
	  histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) + histo2D2->GetBinContent(i,j));
	}
    }

  tmp0 = histo1D0->GetMaximum();
  tmp1 = histo1D1->GetMaximum();
  tmp2 = histo1D2->GetMaximum();
  for (int i = 1; i <= histo1D0->GetNbinsX(); i++)
    {
      histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) / tmp0); 
      histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) / tmp1); 
      histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) / tmp2); 
    }

  histo1D0->SetLineColor(kBlack);
  histo1D0->SetLineWidth(2.0);
  histo1D0->SetXTitle("Track \\eta");
  histo1D0->SetYTitle("Norm. entries");
  histo1D0->GetYaxis()->SetRangeUser(0.,1.5);
  histo1D0->Draw();

  histo1D1->SetLineColor(kRed);
  histo1D1->SetLineWidth(2.0);
  histo1D1->Draw("same");

  if (WithMC == true)
    {
      histo1D2->SetLineColor(kBlue);
      histo1D2->SetLineWidth(2.0);
      histo1D2->Draw("same");
    }

  histo1D0->SetStats(false);
  histo1D1->SetStats(false);
  if (WithMC == true) histo1D2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo1D0,"Data");
  leg->AddEntry(histo1D1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo1D2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("GenTracks_Eta_Eff.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo2D0 = (TH2F*)file0->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi algo");
  histo2D1 = (TH2F*)file1->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi algo");
  histo2D2 = (TH2F*)file2->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi algo");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo1D0 = new TH1F("H1D0","H1D0", histo2D0->GetNbinsX(), histo2D0->GetBinLowEdge(1), histo2D0->GetBinLowEdge(histo2D0->GetNbinsX()) + histo2D0->GetBinWidth(histo2D0->GetNbinsX()));
  histo1D1 = new TH1F("H1D1","H1D1", histo2D1->GetNbinsX(), histo2D1->GetBinLowEdge(1), histo2D1->GetBinLowEdge(histo2D1->GetNbinsX()) + histo2D1->GetBinWidth(histo2D1->GetNbinsX()));
  histo1D2 = new TH1F("H1D2","H1D2", histo2D2->GetNbinsX(), histo2D2->GetBinLowEdge(1), histo2D2->GetBinLowEdge(histo2D2->GetNbinsX()) + histo2D2->GetBinWidth(histo2D2->GetNbinsX()));
  
  for (int i = 1; i <= histo2D0->GetNbinsX(); i++)
    {
      for (int j = 1; j <= histo2D0->GetNbinsY(); j++)
	{
	  histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) + histo2D0->GetBinContent(i,j));
	  histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) + histo2D1->GetBinContent(i,j));
	  histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) + histo2D2->GetBinContent(i,j));
	}
    }

  tmp0 = histo1D0->GetMaximum();
  tmp1 = histo1D1->GetMaximum();
  tmp2 = histo1D2->GetMaximum();
  for (int i = 1; i <= histo1D0->GetNbinsX(); i++)
    {
      histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) / tmp0); 
      histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) / tmp1); 
      histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) / tmp2); 
    }

  histo1D0->SetLineColor(kBlack);
  histo1D0->SetLineWidth(2.0);
  histo1D0->SetXTitle("Track \\eta");
  histo1D0->SetYTitle("Norm. entries");
  histo1D0->GetYaxis()->SetRangeUser(0.,1.5);
  histo1D0->Draw();

  histo1D1->SetLineColor(kRed);
  histo1D1->SetLineWidth(2.0);
  histo1D1->Draw("same");

  if (WithMC == true)
    {
      histo1D2->SetLineColor(kBlue);
      histo1D2->SetLineWidth(2.0);
      histo1D2->Draw("same");
    }

  histo1D0->SetStats(false);
  histo1D1->SetStats(false);
  if (WithMC == true) histo1D2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo1D0,"Data");
  leg->AddEntry(histo1D1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo1D2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("GenTracks_Eta_AlgoEff.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo2D0 = (TH2F*)file0->Get("MyProcess/PixelTracks/Pixeltracks distribution in eta and phi");
  histo2D1 = (TH2F*)file1->Get("MyProcess/PixelTracks/Pixeltracks distribution in eta and phi");
  histo2D2 = (TH2F*)file2->Get("MyProcess/PixelTracks/Pixeltracks distribution in eta and phi");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo1D0 = new TH1F("H1D0","H1D0", histo2D0->GetNbinsX(), histo2D0->GetBinLowEdge(1), histo2D0->GetBinLowEdge(histo2D0->GetNbinsX()) + histo2D0->GetBinWidth(histo2D0->GetNbinsX()));
  histo1D1 = new TH1F("H1D1","H1D1", histo2D1->GetNbinsX(), histo2D1->GetBinLowEdge(1), histo2D1->GetBinLowEdge(histo2D1->GetNbinsX()) + histo2D1->GetBinWidth(histo2D1->GetNbinsX()));
  histo1D2 = new TH1F("H1D2","H1D2", histo2D2->GetNbinsX(), histo2D2->GetBinLowEdge(1), histo2D2->GetBinLowEdge(histo2D2->GetNbinsX()) + histo2D2->GetBinWidth(histo2D2->GetNbinsX()));
  
  for (int i = 1; i <= histo2D0->GetNbinsX(); i++)
    {
      for (int j = 1; j <= histo2D0->GetNbinsY(); j++)
	{
	  histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) + histo2D0->GetBinContent(i,j));
	  histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) + histo2D1->GetBinContent(i,j));
	  histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) + histo2D2->GetBinContent(i,j));
	}
    }

  tmp0 = histo1D0->GetMaximum();
  tmp1 = histo1D1->GetMaximum();
  tmp2 = histo1D2->GetMaximum();
  for (int i = 1; i <= histo1D0->GetNbinsX(); i++)
    {
      histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) / tmp0); 
      histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) / tmp1); 
      histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) / tmp2); 
    }

  histo1D0->SetLineColor(kBlack);
  histo1D0->SetLineWidth(2.0);
  histo1D0->SetXTitle("Track \\eta");
  histo1D0->SetYTitle("Norm. entries");
  histo1D0->GetYaxis()->SetRangeUser(0.,1.5);
  histo1D0->Draw();

  histo1D1->SetLineColor(kRed);
  histo1D1->SetLineWidth(2.0);
  histo1D1->Draw("same");

  if (WithMC == true)
    {
      histo1D2->SetLineColor(kBlue);
      histo1D2->SetLineWidth(2.0);
      histo1D2->Draw("same");
    }

  histo1D0->SetStats(false);
  histo1D1->SetStats(false);
  if (WithMC == true) histo1D2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo1D0,"Data");
  leg->AddEntry(histo1D1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo1D2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("PixTracks_Eta.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo2D0 = (TH2F*)file0->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and pt");
  histo2D1 = (TH2F*)file1->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and pt");
  histo2D2 = (TH2F*)file2->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and pt");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo1D0 = new TH1F("H1D0","H1D0", histo2D0->GetNbinsY(), histo2D0->GetBinLowEdge(1), histo2D0->GetBinLowEdge(histo2D0->GetNbinsX()) + histo2D0->GetBinWidth(histo2D0->GetNbinsY()));
  histo1D1 = new TH1F("H1D1","H1D1", histo2D1->GetNbinsY(), histo2D1->GetBinLowEdge(1), histo2D1->GetBinLowEdge(histo2D1->GetNbinsX()) + histo2D1->GetBinWidth(histo2D1->GetNbinsY()));
  histo1D2 = new TH1F("H1D2","H1D2", histo2D2->GetNbinsY(), histo2D2->GetBinLowEdge(1), histo2D2->GetBinLowEdge(histo2D2->GetNbinsX()) + histo2D2->GetBinWidth(histo2D2->GetNbinsY()));
  
  for (int i = 1; i <= histo2D0->GetNbinsY(); i++)
    {
      for (int j = 1; j <= histo2D0->GetNbinsX(); j++)
	{
	  histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) + histo2D0->GetBinContent(j,i));
	  histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) + histo2D1->GetBinContent(j,i));
	  histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) + histo2D2->GetBinContent(j,i));
	}
    }

  tmp0 = histo1D0->GetMaximum();
  tmp1 = histo1D1->GetMaximum();
  tmp2 = histo1D2->GetMaximum();
  for (int i = 1; i <= histo1D0->GetNbinsX(); i++)
    {
      histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) / tmp0); 
      histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) / tmp1); 
      histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) / tmp2); 
    }

  histo1D0->SetLineColor(kBlack);
  histo1D0->SetLineWidth(2.0);
  histo1D0->SetXTitle("Track p_{t} (GeV/c)");
  histo1D0->SetYTitle("Norm. entries");
  histo1D0->GetYaxis()->SetRangeUser(0.,1.5);
  histo1D0->Draw();

  histo1D1->SetLineColor(kRed);
  histo1D1->SetLineWidth(2.0);
  histo1D1->Draw("same");

  if (WithMC == true)
    {
      histo1D2->SetLineColor(kBlue);
      histo1D2->SetLineWidth(2.0);
      histo1D2->Draw("same");
    }

  histo1D0->SetStats(false);
  histo1D1->SetStats(false);
  if (WithMC == true) histo1D2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo1D0,"Data");
  leg->AddEntry(histo1D1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo1D2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("GenTracks_Pt_Eff.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo2D0 = (TH2F*)file0->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and pt algo");
  histo2D1 = (TH2F*)file1->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and pt algo");
  histo2D2 = (TH2F*)file2->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and pt algo");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo1D0 = new TH1F("H1D0","H1D0", histo2D0->GetNbinsY(), histo2D0->GetBinLowEdge(1), histo2D0->GetBinLowEdge(histo2D0->GetNbinsX()) + histo2D0->GetBinWidth(histo2D0->GetNbinsY()));
  histo1D1 = new TH1F("H1D1","H1D1", histo2D1->GetNbinsY(), histo2D1->GetBinLowEdge(1), histo2D1->GetBinLowEdge(histo2D1->GetNbinsX()) + histo2D1->GetBinWidth(histo2D1->GetNbinsY()));
  histo1D2 = new TH1F("H1D2","H1D2", histo2D2->GetNbinsY(), histo2D2->GetBinLowEdge(1), histo2D2->GetBinLowEdge(histo2D2->GetNbinsX()) + histo2D2->GetBinWidth(histo2D2->GetNbinsY()));
  
  for (int i = 1; i <= histo2D0->GetNbinsY(); i++)
    {
      for (int j = 1; j <= histo2D0->GetNbinsX(); j++)
	{
	  histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) + histo2D0->GetBinContent(j,i));
	  histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) + histo2D1->GetBinContent(j,i));
	  histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) + histo2D2->GetBinContent(j,i));
	}
    }

  tmp0 = histo1D0->GetMaximum();
  tmp1 = histo1D1->GetMaximum();
  tmp2 = histo1D2->GetMaximum();
  for (int i = 1; i <= histo1D0->GetNbinsX(); i++)
    {
      histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) / tmp0); 
      histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) / tmp1); 
      histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) / tmp2); 
    }

  histo1D0->SetLineColor(kBlack);
  histo1D0->SetLineWidth(2.0);
  histo1D0->SetXTitle("Track p_{t} (GeV/c)");
  histo1D0->SetYTitle("Norm. entries");
  histo1D0->GetYaxis()->SetRangeUser(0.,1.5);
  histo1D0->Draw();

  histo1D1->SetLineColor(kRed);
  histo1D1->SetLineWidth(2.0);
  histo1D1->Draw("same");

  if (WithMC == true)
    {
      histo1D2->SetLineColor(kBlue);
      histo1D2->SetLineWidth(2.0);
      histo1D2->Draw("same");
    }

  histo1D0->SetStats(false);
  histo1D1->SetStats(false);
  if (WithMC == true) histo1D2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo1D0,"Data");
  leg->AddEntry(histo1D1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo1D2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("GenTracks_Pt_AlgoEff.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo2D0 = (TH2F*)file0->Get("MyProcess/PixelTracks/Pixeltracks distribution in eta and pt");
  histo2D1 = (TH2F*)file1->Get("MyProcess/PixelTracks/Pixeltracks distribution in eta and pt");
  histo2D2 = (TH2F*)file2->Get("MyProcess/PixelTracks/Pixeltracks distribution in eta and pt");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo1D0 = new TH1F("H1D0","H1D0", histo2D0->GetNbinsY(), histo2D0->GetBinLowEdge(1), histo2D0->GetBinLowEdge(histo2D0->GetNbinsX()) + histo2D0->GetBinWidth(histo2D0->GetNbinsY()));
  histo1D1 = new TH1F("H1D1","H1D1", histo2D1->GetNbinsY(), histo2D1->GetBinLowEdge(1), histo2D1->GetBinLowEdge(histo2D1->GetNbinsX()) + histo2D1->GetBinWidth(histo2D1->GetNbinsY()));
  histo1D2 = new TH1F("H1D2","H1D2", histo2D2->GetNbinsY(), histo2D2->GetBinLowEdge(1), histo2D2->GetBinLowEdge(histo2D2->GetNbinsX()) + histo2D2->GetBinWidth(histo2D2->GetNbinsY()));
  
  for (int i = 1; i <= histo2D0->GetNbinsY(); i++)
    {
      for (int j = 1; j <= histo2D0->GetNbinsX(); j++)
	{
	  histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) + histo2D0->GetBinContent(j,i));
	  histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) + histo2D1->GetBinContent(j,i));
	  histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) + histo2D2->GetBinContent(j,i));
	}
    }

  tmp0 = histo1D0->GetMaximum();
  tmp1 = histo1D1->GetMaximum();
  tmp2 = histo1D2->GetMaximum();
  for (int i = 1; i <= histo1D0->GetNbinsX(); i++)
    {
      histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) / tmp0); 
      histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) / tmp1); 
      histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) / tmp2); 
    }

  histo1D0->SetLineColor(kBlack);
  histo1D0->SetLineWidth(2.0);
  histo1D0->SetXTitle("Track p_{t} (GeV/c)");
  histo1D0->SetYTitle("Norm. entries");
  histo1D0->GetYaxis()->SetRangeUser(0.,1.5);
  histo1D0->Draw();

  histo1D1->SetLineColor(kRed);
  histo1D1->SetLineWidth(2.0);
  histo1D1->Draw("same");

  if (WithMC == true)
    {
      histo1D2->SetLineColor(kBlue);
      histo1D2->SetLineWidth(2.0);
      histo1D2->Draw("same");
    }

  histo1D0->SetStats(false);
  histo1D1->SetStats(false);
  if (WithMC == true) histo1D2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo1D0,"Data");
  leg->AddEntry(histo1D1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo1D2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("PixTracks_Pt.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo2D0 = (TH2F*)file0->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi");
  histo2D1 = (TH2F*)file1->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi");
  histo2D2 = (TH2F*)file2->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo1D0 = new TH1F("H1D0","H1D0", histo2D0->GetNbinsY(), histo2D0->ProjectionY()->GetBinLowEdge(1), histo2D0->ProjectionY()->GetBinLowEdge(histo2D0->GetNbinsY()) + histo2D0->ProjectionY()->GetBinWidth(histo2D0->GetNbinsY()));
  histo1D1 = new TH1F("H1D1","H1D1", histo2D1->GetNbinsY(), histo2D1->ProjectionY()->GetBinLowEdge(1), histo2D1->ProjectionY()->GetBinLowEdge(histo2D1->GetNbinsY()) + histo2D1->ProjectionY()->GetBinWidth(histo2D1->GetNbinsY()));
  histo1D2 = new TH1F("H1D2","H1D2", histo2D2->GetNbinsY(), histo2D2->ProjectionY()->GetBinLowEdge(1), histo2D2->ProjectionY()->GetBinLowEdge(histo2D2->GetNbinsY()) + histo2D2->ProjectionY()->GetBinWidth(histo2D2->GetNbinsY()));
  
  for (int i = 1; i <= histo2D0->GetNbinsY(); i++)
    {
      for (int j = histo2D0->GetXaxis()->FindBin(-MaxEtaTkTrk); j <= histo2D0->GetXaxis()->FindBin(-EtaThr); j++)
	{
	  histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) + histo2D0->GetBinContent(j,i));
	  histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) + histo2D1->GetBinContent(j,i));
	  histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) + histo2D2->GetBinContent(j,i));
	}
    }

  tmp0 = histo1D0->GetMaximum();
  tmp1 = histo1D1->GetMaximum();
  tmp2 = histo1D2->GetMaximum();
  for (int i = 1; i <= histo1D0->GetNbinsX(); i++)
    {
      histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) / tmp0); 
      histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) / tmp1); 
      histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) / tmp2); 
    }
  
  histo1D0->SetLineColor(kBlack);
  histo1D0->SetLineWidth(2.0);
  histo1D0->SetXTitle("Track \\phi (deg)");
  histo1D0->SetYTitle("Norm. entries");
  histo1D0->GetYaxis()->SetRangeUser(0.,1.5);
  histo1D0->Draw();

  histo1D1->SetLineColor(kRed);
  histo1D1->SetLineWidth(2.0);
  histo1D1->Draw("same");

  if (WithMC == true)
    {
      histo1D2->SetLineColor(kBlue);
      histo1D2->SetLineWidth(2.0);
      histo1D2->Draw("same");
    }

  histo1D0->SetStats(false);
  histo1D1->SetStats(false);
  if (WithMC == true) histo1D2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo1D0,"Data");
  leg->AddEntry(histo1D1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo1D2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("GenTracks_Phi1_Eff.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo2D0 = (TH2F*)file0->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi");
  histo2D1 = (TH2F*)file1->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi");
  histo2D2 = (TH2F*)file2->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo1D0 = new TH1F("H1D0","H1D0", histo2D0->GetNbinsY(), histo2D0->ProjectionY()->GetBinLowEdge(1), histo2D0->ProjectionY()->GetBinLowEdge(histo2D0->GetNbinsY()) + histo2D0->ProjectionY()->GetBinWidth(histo2D0->GetNbinsY()));
  histo1D1 = new TH1F("H1D1","H1D1", histo2D1->GetNbinsY(), histo2D1->ProjectionY()->GetBinLowEdge(1), histo2D1->ProjectionY()->GetBinLowEdge(histo2D1->GetNbinsY()) + histo2D1->ProjectionY()->GetBinWidth(histo2D1->GetNbinsY()));
  histo1D2 = new TH1F("H1D2","H1D2", histo2D2->GetNbinsY(), histo2D2->ProjectionY()->GetBinLowEdge(1), histo2D2->ProjectionY()->GetBinLowEdge(histo2D2->GetNbinsY()) + histo2D2->ProjectionY()->GetBinWidth(histo2D2->GetNbinsY()));
  
  for (int i = 1; i <= histo2D0->GetNbinsY(); i++)
    {
      for (int j = histo2D0->GetXaxis()->FindBin(-EtaThr); j <= histo2D0->GetXaxis()->FindBin(EtaThr); j++)
	{
	  histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) + histo2D0->GetBinContent(j,i));
	  histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) + histo2D1->GetBinContent(j,i));
	  histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) + histo2D2->GetBinContent(j,i));
	}
    }

  tmp0 = histo1D0->GetMaximum();
  tmp1 = histo1D1->GetMaximum();
  tmp2 = histo1D2->GetMaximum();
  for (int i = 1; i <= histo1D0->GetNbinsX(); i++)
    {
      histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) / tmp0); 
      histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) / tmp1); 
      histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) / tmp2); 
    }

  histo1D0->SetLineColor(kBlack);
  histo1D0->SetLineWidth(2.0);
  histo1D0->SetXTitle("Track \\phi (deg)");
  histo1D0->SetYTitle("Norm. entries");
  histo1D0->GetYaxis()->SetRangeUser(0.,1.5);
  histo1D0->Draw();

  histo1D1->SetLineColor(kRed);
  histo1D1->SetLineWidth(2.0);
  histo1D1->Draw("same");

  if (WithMC == true)
    {
      histo1D2->SetLineColor(kBlue);
      histo1D2->SetLineWidth(2.0);
      histo1D2->Draw("same");
    }

  histo1D0->SetStats(false);
  histo1D1->SetStats(false);
  if (WithMC == true) histo1D2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo1D0,"Data");
  leg->AddEntry(histo1D1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo1D2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("GenTracks_Phi2_Eff.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo2D0 = (TH2F*)file0->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi");
  histo2D1 = (TH2F*)file1->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi");
  histo2D2 = (TH2F*)file2->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo1D0 = new TH1F("H1D0","H1D0", histo2D0->GetNbinsY(), histo2D0->ProjectionY()->GetBinLowEdge(1), histo2D0->ProjectionY()->GetBinLowEdge(histo2D0->GetNbinsY()) + histo2D0->ProjectionY()->GetBinWidth(histo2D0->GetNbinsY()));
  histo1D1 = new TH1F("H1D1","H1D1", histo2D1->GetNbinsY(), histo2D1->ProjectionY()->GetBinLowEdge(1), histo2D1->ProjectionY()->GetBinLowEdge(histo2D1->GetNbinsY()) + histo2D1->ProjectionY()->GetBinWidth(histo2D1->GetNbinsY()));
  histo1D2 = new TH1F("H1D2","H1D2", histo2D2->GetNbinsY(), histo2D2->ProjectionY()->GetBinLowEdge(1), histo2D2->ProjectionY()->GetBinLowEdge(histo2D2->GetNbinsY()) + histo2D2->ProjectionY()->GetBinWidth(histo2D2->GetNbinsY()));
  
  for (int i = 1; i <= histo2D0->GetNbinsY(); i++)
    {
      for (int j = histo2D0->GetXaxis()->FindBin(EtaThr); j <= histo2D0->GetXaxis()->FindBin(MaxEtaTkTrk); j++)
	{
	  histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) + histo2D0->GetBinContent(j,i));
	  histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) + histo2D1->GetBinContent(j,i));
	  histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) + histo2D2->GetBinContent(j,i));
	}
    }

  tmp0 = histo1D0->GetMaximum();
  tmp1 = histo1D1->GetMaximum();
  tmp2 = histo1D2->GetMaximum();
  for (int i = 1; i <= histo1D0->GetNbinsX(); i++)
    {
      histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) / tmp0); 
      histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) / tmp1); 
      histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) / tmp2); 
    }

  histo1D0->SetLineColor(kBlack);
  histo1D0->SetLineWidth(2.0);
  histo1D0->SetXTitle("Track \\phi (deg)");
  histo1D0->SetYTitle("Norm. entries");
  histo1D0->GetYaxis()->SetRangeUser(0.,1.5);
  histo1D0->Draw();

  histo1D1->SetLineColor(kRed);
  histo1D1->SetLineWidth(2.0);
  histo1D1->Draw("same");

  if (WithMC == true)
    {
      histo1D2->SetLineColor(kBlue);
      histo1D2->SetLineWidth(2.0);
      histo1D2->Draw("same");
    }

  histo1D0->SetStats(false);
  histo1D1->SetStats(false);
  if (WithMC == true) histo1D2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo1D0,"Data");
  leg->AddEntry(histo1D1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo1D2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("GenTracks_Phi3_Eff.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo2D0 = (TH2F*)file0->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi algo");
  histo2D1 = (TH2F*)file1->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi algo");
  histo2D2 = (TH2F*)file2->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi algo");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo1D0 = new TH1F("H1D0","H1D0", histo2D0->GetNbinsY(), histo2D0->ProjectionY()->GetBinLowEdge(1), histo2D0->ProjectionY()->GetBinLowEdge(histo2D0->GetNbinsY()) + histo2D0->ProjectionY()->GetBinWidth(histo2D0->GetNbinsY()));
  histo1D1 = new TH1F("H1D1","H1D1", histo2D1->GetNbinsY(), histo2D1->ProjectionY()->GetBinLowEdge(1), histo2D1->ProjectionY()->GetBinLowEdge(histo2D1->GetNbinsY()) + histo2D1->ProjectionY()->GetBinWidth(histo2D1->GetNbinsY()));
  histo1D2 = new TH1F("H1D2","H1D2", histo2D2->GetNbinsY(), histo2D2->ProjectionY()->GetBinLowEdge(1), histo2D2->ProjectionY()->GetBinLowEdge(histo2D2->GetNbinsY()) + histo2D2->ProjectionY()->GetBinWidth(histo2D2->GetNbinsY()));
  
  for (int i = 1; i <= histo2D0->GetNbinsY(); i++)
    {
      for (int j = histo2D0->GetXaxis()->FindBin(-MaxEtaTkTrk); j <= histo2D0->GetXaxis()->FindBin(-EtaThr); j++)
	{
	  histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) + histo2D0->GetBinContent(j,i));
	  histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) + histo2D1->GetBinContent(j,i));
	  histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) + histo2D2->GetBinContent(j,i));
	}
    }

  tmp0 = histo1D0->GetMaximum();
  tmp1 = histo1D1->GetMaximum();
  tmp2 = histo1D2->GetMaximum();
  for (int i = 1; i <= histo1D0->GetNbinsX(); i++)
    {
      histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) / tmp0); 
      histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) / tmp1); 
      histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) / tmp2); 
    }

  histo1D0->SetLineColor(kBlack);
  histo1D0->SetLineWidth(2.0);
  histo1D0->SetXTitle("Track \\phi (deg)");
  histo1D0->SetYTitle("Norm. entries");
  histo1D0->GetYaxis()->SetRangeUser(0.,1.5);
  histo1D0->Draw();

  histo1D1->SetLineColor(kRed);
  histo1D1->SetLineWidth(2.0);
  histo1D1->Draw("same");

  if (WithMC == true)
    {
      histo1D2->SetLineColor(kBlue);
      histo1D2->SetLineWidth(2.0);
      histo1D2->Draw("same");
    }

  histo1D0->SetStats(false);
  histo1D1->SetStats(false);
  if (WithMC == true) histo1D2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo1D0,"Data");
  leg->AddEntry(histo1D1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo1D2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("GenTracks_Phi1_AlgoEff.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo2D0 = (TH2F*)file0->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi algo");
  histo2D1 = (TH2F*)file1->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi algo");
  histo2D2 = (TH2F*)file2->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi algo");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo1D0 = new TH1F("H1D0","H1D0", histo2D0->GetNbinsY(), histo2D0->ProjectionY()->GetBinLowEdge(1), histo2D0->ProjectionY()->GetBinLowEdge(histo2D0->GetNbinsY()) + histo2D0->ProjectionY()->GetBinWidth(histo2D0->GetNbinsY()));
  histo1D1 = new TH1F("H1D1","H1D1", histo2D1->GetNbinsY(), histo2D1->ProjectionY()->GetBinLowEdge(1), histo2D1->ProjectionY()->GetBinLowEdge(histo2D1->GetNbinsY()) + histo2D1->ProjectionY()->GetBinWidth(histo2D1->GetNbinsY()));
  histo1D2 = new TH1F("H1D2","H1D2", histo2D2->GetNbinsY(), histo2D2->ProjectionY()->GetBinLowEdge(1), histo2D2->ProjectionY()->GetBinLowEdge(histo2D2->GetNbinsY()) + histo2D2->ProjectionY()->GetBinWidth(histo2D2->GetNbinsY()));
  
  for (int i = 1; i <= histo2D0->GetNbinsY(); i++)
    {
      for (int j = histo2D0->GetXaxis()->FindBin(-EtaThr); j <= histo2D0->GetXaxis()->FindBin(EtaThr); j++)
	{
	  histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) + histo2D0->GetBinContent(j,i));
	  histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) + histo2D1->GetBinContent(j,i));
	  histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) + histo2D2->GetBinContent(j,i));
	}
    }

  tmp0 = histo1D0->GetMaximum();
  tmp1 = histo1D1->GetMaximum();
  tmp2 = histo1D2->GetMaximum();
  for (int i = 1; i <= histo1D0->GetNbinsX(); i++)
    {
      histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) / tmp0); 
      histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) / tmp1); 
      histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) / tmp2); 
    }

  histo1D0->SetLineColor(kBlack);
  histo1D0->SetLineWidth(2.0);
  histo1D0->SetXTitle("Track \\phi (deg)");
  histo1D0->SetYTitle("Norm. entries");
  histo1D0->GetYaxis()->SetRangeUser(0.,1.5);
  histo1D0->Draw();

  histo1D1->SetLineColor(kRed);
  histo1D1->SetLineWidth(2.0);
  histo1D1->Draw("same");

  if (WithMC == true)
    {
      histo1D2->SetLineColor(kBlue);
      histo1D2->SetLineWidth(2.0);
      histo1D2->Draw("same");
    }

  histo1D0->SetStats(false);
  histo1D1->SetStats(false);
  if (WithMC == true) histo1D2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo1D0,"Data");
  leg->AddEntry(histo1D1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo1D2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("GenTracks_Phi2_AlgoEff.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo2D0 = (TH2F*)file0->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi algo");
  histo2D1 = (TH2F*)file1->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi algo");
  histo2D2 = (TH2F*)file2->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi algo");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo1D0 = new TH1F("H1D0","H1D0", histo2D0->GetNbinsY(), histo2D0->ProjectionY()->GetBinLowEdge(1), histo2D0->ProjectionY()->GetBinLowEdge(histo2D0->GetNbinsY()) + histo2D0->ProjectionY()->GetBinWidth(histo2D0->GetNbinsY()));
  histo1D1 = new TH1F("H1D1","H1D1", histo2D1->GetNbinsY(), histo2D1->ProjectionY()->GetBinLowEdge(1), histo2D1->ProjectionY()->GetBinLowEdge(histo2D1->GetNbinsY()) + histo2D1->ProjectionY()->GetBinWidth(histo2D1->GetNbinsY()));
  histo1D2 = new TH1F("H1D2","H1D2", histo2D2->GetNbinsY(), histo2D2->ProjectionY()->GetBinLowEdge(1), histo2D2->ProjectionY()->GetBinLowEdge(histo2D2->GetNbinsY()) + histo2D2->ProjectionY()->GetBinWidth(histo2D2->GetNbinsY()));
  
  for (int i = 1; i <= histo2D0->GetNbinsY(); i++)
    {
      for (int j = histo2D0->GetXaxis()->FindBin(EtaThr); j <= histo2D0->GetXaxis()->FindBin(MaxEtaTkTrk); j++)
	{
	  histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) + histo2D0->GetBinContent(j,i));
	  histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) + histo2D1->GetBinContent(j,i));
	  histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) + histo2D2->GetBinContent(j,i));
	}
    }

  tmp0 = histo1D0->GetMaximum();
  tmp1 = histo1D1->GetMaximum();
  tmp2 = histo1D2->GetMaximum();
  for (int i = 1; i <= histo1D0->GetNbinsX(); i++)
    {
      histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) / tmp0); 
      histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) / tmp1); 
      histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) / tmp2); 
    }

  histo1D0->SetLineColor(kBlack);
  histo1D0->SetLineWidth(2.0);
  histo1D0->SetXTitle("Track \\phi (deg)");
  histo1D0->SetYTitle("Norm. entries");
  histo1D0->GetYaxis()->SetRangeUser(0.,1.5);
  histo1D0->Draw();

  histo1D1->SetLineColor(kRed);
  histo1D1->SetLineWidth(2.0);
  histo1D1->Draw("same");

  if (WithMC == true)
    {
      histo1D2->SetLineColor(kBlue);
      histo1D2->SetLineWidth(2.0);
      histo1D2->Draw("same");
    }

  histo1D0->SetStats(false);
  histo1D1->SetStats(false);
  if (WithMC == true) histo1D2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo1D0,"Data");
  leg->AddEntry(histo1D1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo1D2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("GenTracks_Phi3_AlgoEff.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo2D0 = (TH2F*)file0->Get("MyProcess/PixelTracks/Pixeltracks distribution in eta and phi");
  histo2D1 = (TH2F*)file1->Get("MyProcess/PixelTracks/Pixeltracks distribution in eta and phi");
  histo2D2 = (TH2F*)file2->Get("MyProcess/PixelTracks/Pixeltracks distribution in eta and phi");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo1D0 = new TH1F("H1D0","H1D0", histo2D0->GetNbinsY(), histo2D0->ProjectionY()->GetBinLowEdge(1), histo2D0->ProjectionY()->GetBinLowEdge(histo2D0->GetNbinsY()) + histo2D0->ProjectionY()->GetBinWidth(histo2D0->GetNbinsY()));
  histo1D1 = new TH1F("H1D1","H1D1", histo2D1->GetNbinsY(), histo2D1->ProjectionY()->GetBinLowEdge(1), histo2D1->ProjectionY()->GetBinLowEdge(histo2D1->GetNbinsY()) + histo2D1->ProjectionY()->GetBinWidth(histo2D1->GetNbinsY()));
  histo1D2 = new TH1F("H1D2","H1D2", histo2D2->GetNbinsY(), histo2D2->ProjectionY()->GetBinLowEdge(1), histo2D2->ProjectionY()->GetBinLowEdge(histo2D2->GetNbinsY()) + histo2D2->ProjectionY()->GetBinWidth(histo2D2->GetNbinsY()));
  
  for (int i = 1; i <= histo2D0->GetNbinsY(); i++)
    {
      for (int j = histo2D0->GetXaxis()->FindBin(-MaxEtaTkTrk); j <= histo2D0->GetXaxis()->FindBin(-EtaThr); j++)
	{
	  histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) + histo2D0->GetBinContent(j,i));
	  histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) + histo2D1->GetBinContent(j,i));
	  histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) + histo2D2->GetBinContent(j,i));
	}
    }

  tmp0 = histo1D0->GetMaximum();
  tmp1 = histo1D1->GetMaximum();
  tmp2 = histo1D2->GetMaximum();
  for (int i = 1; i <= histo1D0->GetNbinsX(); i++)
    {
      histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) / tmp0); 
      histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) / tmp1); 
      histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) / tmp2); 
    }

  histo1D0->SetLineColor(kBlack);
  histo1D0->SetLineWidth(2.0);
  histo1D0->SetXTitle("Track \\phi (deg)");
  histo1D0->SetYTitle("Norm. entries");
  histo1D0->GetYaxis()->SetRangeUser(0.,1.5);
  histo1D0->Draw();

  histo1D1->SetLineColor(kRed);
  histo1D1->SetLineWidth(2.0);
  histo1D1->Draw("same");

  if (WithMC == true)
    {
      histo1D2->SetLineColor(kBlue);
      histo1D2->SetLineWidth(2.0);
      histo1D2->Draw("same");
    }

  histo1D0->SetStats(false);
  histo1D1->SetStats(false);
  if (WithMC == true) histo1D2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo1D0,"Data");
  leg->AddEntry(histo1D1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo1D2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("PixTracks_Phi1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo2D0 = (TH2F*)file0->Get("MyProcess/PixelTracks/Pixeltracks distribution in eta and phi");
  histo2D1 = (TH2F*)file1->Get("MyProcess/PixelTracks/Pixeltracks distribution in eta and phi");
  histo2D2 = (TH2F*)file2->Get("MyProcess/PixelTracks/Pixeltracks distribution in eta and phi");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo1D0 = new TH1F("H1D0","H1D0", histo2D0->GetNbinsY(), histo2D0->ProjectionY()->GetBinLowEdge(1), histo2D0->ProjectionY()->GetBinLowEdge(histo2D0->GetNbinsY()) + histo2D0->ProjectionY()->GetBinWidth(histo2D0->GetNbinsY()));
  histo1D1 = new TH1F("H1D1","H1D1", histo2D1->GetNbinsY(), histo2D1->ProjectionY()->GetBinLowEdge(1), histo2D1->ProjectionY()->GetBinLowEdge(histo2D1->GetNbinsY()) + histo2D1->ProjectionY()->GetBinWidth(histo2D1->GetNbinsY()));
  histo1D2 = new TH1F("H1D2","H1D2", histo2D2->GetNbinsY(), histo2D2->ProjectionY()->GetBinLowEdge(1), histo2D2->ProjectionY()->GetBinLowEdge(histo2D2->GetNbinsY()) + histo2D2->ProjectionY()->GetBinWidth(histo2D2->GetNbinsY()));
  
  for (int i = 1; i <= histo2D0->GetNbinsY(); i++)
    {
      for (int j = histo2D0->GetXaxis()->FindBin(-EtaThr); j <= histo2D0->GetXaxis()->FindBin(EtaThr); j++)
	{
	  histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) + histo2D0->GetBinContent(j,i));
	  histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) + histo2D1->GetBinContent(j,i));
	  histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) + histo2D2->GetBinContent(j,i));
	}
    }

  tmp0 = histo1D0->GetMaximum();
  tmp1 = histo1D1->GetMaximum();
  tmp2 = histo1D2->GetMaximum();
  for (int i = 1; i <= histo1D0->GetNbinsX(); i++)
    {
      histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) / tmp0); 
      histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) / tmp1); 
      histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) / tmp2); 
    }

  histo1D0->SetLineColor(kBlack);
  histo1D0->SetLineWidth(2.0);
  histo1D0->SetXTitle("Track \\phi (deg)");
  histo1D0->SetYTitle("Norm. entries");
  histo1D0->GetYaxis()->SetRangeUser(0.,1.5);
  histo1D0->Draw();

  histo1D1->SetLineColor(kRed);
  histo1D1->SetLineWidth(2.0);
  histo1D1->Draw("same");

  if (WithMC == true)
    {
      histo1D2->SetLineColor(kBlue);
      histo1D2->SetLineWidth(2.0);
      histo1D2->Draw("same");
    }

  histo1D0->SetStats(false);
  histo1D1->SetStats(false);
  if (WithMC == true) histo1D2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo1D0,"Data");
  leg->AddEntry(histo1D1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo1D2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("PixTracks_Phi2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo2D0 = (TH2F*)file0->Get("MyProcess/PixelTracks/Pixeltracks distribution in eta and phi");
  histo2D1 = (TH2F*)file1->Get("MyProcess/PixelTracks/Pixeltracks distribution in eta and phi");
  histo2D2 = (TH2F*)file2->Get("MyProcess/PixelTracks/Pixeltracks distribution in eta and phi");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo1D0 = new TH1F("H1D0","H1D0", histo2D0->GetNbinsY(), histo2D0->ProjectionY()->GetBinLowEdge(1), histo2D0->ProjectionY()->GetBinLowEdge(histo2D0->GetNbinsY()) + histo2D0->ProjectionY()->GetBinWidth(histo2D0->GetNbinsY()));
  histo1D1 = new TH1F("H1D1","H1D1", histo2D1->GetNbinsY(), histo2D1->ProjectionY()->GetBinLowEdge(1), histo2D1->ProjectionY()->GetBinLowEdge(histo2D1->GetNbinsY()) + histo2D1->ProjectionY()->GetBinWidth(histo2D1->GetNbinsY()));
  histo1D2 = new TH1F("H1D2","H1D2", histo2D2->GetNbinsY(), histo2D2->ProjectionY()->GetBinLowEdge(1), histo2D2->ProjectionY()->GetBinLowEdge(histo2D2->GetNbinsY()) + histo2D2->ProjectionY()->GetBinWidth(histo2D2->GetNbinsY()));
  
  for (int i = 1; i <= histo2D0->GetNbinsY(); i++)
    {
      for (int j = histo2D0->GetXaxis()->FindBin(EtaThr); j <= histo2D0->GetXaxis()->FindBin(MaxEtaTkTrk); j++)
	{
	  histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) + histo2D0->GetBinContent(j,i));
	  histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) + histo2D1->GetBinContent(j,i));
	  histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) + histo2D2->GetBinContent(j,i));
	}
    }

  tmp0 = histo1D0->GetMaximum();
  tmp1 = histo1D1->GetMaximum();
  tmp2 = histo1D2->GetMaximum();
  for (int i = 1; i <= histo1D0->GetNbinsX(); i++)
    {
      histo1D0->SetBinContent(i, histo1D0->GetBinContent(i) / tmp0); 
      histo1D1->SetBinContent(i, histo1D1->GetBinContent(i) / tmp1); 
      histo1D2->SetBinContent(i, histo1D2->GetBinContent(i) / tmp2); 
    }

  histo1D0->SetLineColor(kBlack);
  histo1D0->SetLineWidth(2.0);
  histo1D0->SetXTitle("Track \\phi (deg)");
  histo1D0->SetYTitle("Norm. entries");
  histo1D0->GetYaxis()->SetRangeUser(0.,1.5);
  histo1D0->Draw();

  histo1D1->SetLineColor(kRed);
  histo1D1->SetLineWidth(2.0);
  histo1D1->Draw("same");

  if (WithMC == true)
    {
      histo1D2->SetLineColor(kBlue);
      histo1D2->SetLineWidth(2.0);
      histo1D2->Draw("same");
    }

  histo1D0->SetStats(false);
  histo1D1->SetStats(false);
  if (WithMC == true) histo1D2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo1D0,"Data");
  leg->AddEntry(histo1D1,"Simulation");
  if (WithMC == true) leg->AddEntry(histo1D2,"MC-Truth");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("PixTracks_Phi3.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/PixelTracks/Trackertracks efficiency vs eta");
  histo1 = (TH1F*)file2->Get("MyProcess/PixelTracks/Total trackertracks efficiency vs eta");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetXaxis()->SetRangeUser(-2.3,2.3);
  histo0->GetYaxis()->SetRangeUser(0.7,1.2);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Gen.track eff. for 3 pix.hit tracks");
  leg->AddEntry(histo1,"Total gen.track eff.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("TrkEff_Eta.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/PixelTracks/Trackertracks efficiency vs pt");
  histo1 = (TH1F*)file2->Get("MyProcess/PixelTracks/Total trackertracks efficiency vs pt");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.7,1.2);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Gen.track eff. for 3 pix.hit tracks");
  leg->AddEntry(histo1,"Total gen.track eff.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("TrkEff_Pt.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/PixelTracks/Trackertracks efficiency vs phi thr1");
  histo1 = (TH1F*)file2->Get("MyProcess/PixelTracks/Total trackertracks efficiency vs phi thr1");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.5,1.4);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"-2.3 < \\eta < -1.5");
  leg->AddEntry(histo0,"Gen.track eff. for 3 pix.hit tracks");
  leg->AddEntry(histo1,"Total gen.track eff.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("TrkEff_Phi_etathr1.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/PixelTracks/Trackertracks efficiency vs phi thr2");
  histo1 = (TH1F*)file2->Get("MyProcess/PixelTracks/Total trackertracks efficiency vs phi thr2");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.5,1.4);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"-1.5 < \\eta < 1.5");
  leg->AddEntry(histo0,"Gen.track eff. for 3 pix.hit tracks");
  leg->AddEntry(histo1,"Total gen.track eff.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("TrkEff_Phi_etathr2.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/PixelTracks/Trackertracks efficiency vs phi thr3");
  histo1 = (TH1F*)file2->Get("MyProcess/PixelTracks/Total trackertracks efficiency vs phi thr3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.5,1.4);
  histo0->Draw();
  
  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"1.5 < \\eta < 2.3");
  leg->AddEntry(histo0,"Gen.track eff. for 3 pix.hit tracks");
  leg->AddEntry(histo1,"Total gen.track eff.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("TrkEff_Phi_etathr3.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/Pixeltracks pt distribution ver1");
  histo1 = (TH1F*)file0->Get("MyProcess/PixelTracks/Pixeltracks pt distribution ver2");
  histo2 = (TH1F*)file0->Get("MyProcess/PixelTracks/Pixeltracks pt distribution ver3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);
  histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Num. Eff. and Algo.Eff.");
  leg->AddEntry(histo1,"Num. Pur.");
  leg->AddEntry(histo2,"Den. Pur.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Data_PixTracks_Pt.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/Pixeltracks eta distribution ver1");
  histo1 = (TH1F*)file0->Get("MyProcess/PixelTracks/Pixeltracks eta distribution ver2");
  histo2 = (TH1F*)file0->Get("MyProcess/PixelTracks/Pixeltracks eta distribution ver3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);
  histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Num. Eff. and Algo.Eff.");
  leg->AddEntry(histo1,"Num. Pur.");
  leg->AddEntry(histo2,"Den. Pur.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Data_PixTracks_Eta.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/Pixeltracks chi2DoF distribution ver1");
  histo1 = (TH1F*)file0->Get("MyProcess/PixelTracks/Pixeltracks chi2DoF distribution ver2");
  histo2 = (TH1F*)file0->Get("MyProcess/PixelTracks/Pixeltracks chi2DoF distribution ver3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);
  histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Num. Eff. and Algo.Eff.");
  leg->AddEntry(histo1,"Num. Pur.");
  leg->AddEntry(histo2,"Den. Pur.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Data_PixTracks_chi2DoF.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/Pixeltracks phi distribution ver1");
  histo1 = (TH1F*)file0->Get("MyProcess/PixelTracks/Pixeltracks phi distribution ver2");
  histo2 = (TH1F*)file0->Get("MyProcess/PixelTracks/Pixeltracks phi distribution ver3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);
  histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Num. Eff. and Algo.Eff.");
  leg->AddEntry(histo1,"Num. Pur.");
  leg->AddEntry(histo2,"Den. Pur.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Data_PixTracks_Phi.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/PixelTracks/Pixeltracks pt distribution ver1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/Pixeltracks pt distribution ver2");
  histo2 = (TH1F*)file1->Get("MyProcess/PixelTracks/Pixeltracks pt distribution ver3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);
  histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Num. Eff. and Algo.Eff.");
  leg->AddEntry(histo1,"Num. Pur.");
  leg->AddEntry(histo2,"Den. Pur.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_Reco_PixTracks_Pt.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/PixelTracks/Pixeltracks eta distribution ver1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/Pixeltracks eta distribution ver2");
  histo2 = (TH1F*)file1->Get("MyProcess/PixelTracks/Pixeltracks eta distribution ver3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);
  histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Num. Eff. and Algo.Eff.");
  leg->AddEntry(histo1,"Num. Pur.");
  leg->AddEntry(histo2,"Den. Pur.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_Reco_PixTracks_Eta.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/PixelTracks/Pixeltracks chi2DoF distribution ver1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/Pixeltracks chi2DoF distribution ver2");
  histo2 = (TH1F*)file1->Get("MyProcess/PixelTracks/Pixeltracks chi2DoF distribution ver3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);
  histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Num. Eff. and Algo.Eff.");
  leg->AddEntry(histo1,"Num. Pur.");
  leg->AddEntry(histo2,"Den. Pur.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_Reco_PixTracks_chi2DoF.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/PixelTracks/Pixeltracks phi distribution ver1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/Pixeltracks phi distribution ver2");
  histo2 = (TH1F*)file1->Get("MyProcess/PixelTracks/Pixeltracks phi distribution ver3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);
  histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Num. Eff. and Algo.Eff.");
  leg->AddEntry(histo1,"Num. Pur.");
  leg->AddEntry(histo2,"Den. Pur.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_Reco_PixTracks_Phi.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/PixelTracks/Pixeltracks pt distribution ver1");
  histo1 = (TH1F*)file2->Get("MyProcess/PixelTracks/Pixeltracks pt distribution ver2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/Pixeltracks pt distribution ver3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);
  histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Num. Eff. and Algo.Eff.");
  leg->AddEntry(histo1,"Num. Pur.");
  leg->AddEntry(histo2,"Den. Pur.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_TrkPar_PixTracks_Pt.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/PixelTracks/Pixeltracks eta distribution ver1");
  histo1 = (TH1F*)file2->Get("MyProcess/PixelTracks/Pixeltracks eta distribution ver2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/Pixeltracks eta distribution ver3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);
  histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Num. Eff. and Algo.Eff.");
  leg->AddEntry(histo1,"Num. Pur.");
  leg->AddEntry(histo2,"Den. Pur.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_TrkPar_PixTracks_Eta.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/PixelTracks/Pixeltracks chi2DoF distribution ver1");
  histo1 = (TH1F*)file2->Get("MyProcess/PixelTracks/Pixeltracks chi2DoF distribution ver2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/Pixeltracks chi2DoF distribution ver3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);
  histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Num. Eff. and Algo.Eff.");
  leg->AddEntry(histo1,"Num. Pur.");
  leg->AddEntry(histo2,"Den. Pur.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_TrkPar_PixTracks_chi2DoF.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/PixelTracks/Pixeltracks phi distribution ver1");
  histo1 = (TH1F*)file2->Get("MyProcess/PixelTracks/Pixeltracks phi distribution ver2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/Pixeltracks phi distribution ver3");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);
  histo2->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Num. Eff. and Algo.Eff.");
  leg->AddEntry(histo1,"Num. Pur.");
  leg->AddEntry(histo2,"Den. Pur.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_TrkPar_PixTracks_Phi.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/Trackertracks pt distribution ver1");
  histo1 = (TH1F*)file0->Get("MyProcess/PixelTracks/Trackertracks pt distribution ver2");
  histo2 = (TH1F*)file0->Get("MyProcess/PixelTracks/Trackertracks pt distribution ver3");
  histo3 = (TH1F*)file0->Get("MyProcess/PixelTracks/Trackertracks pt distribution ver4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("same");

  histo3->SetLineColor(kGreen);
  histo3->SetLineWidth(2.0);
  histo3->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);
  histo2->SetStats(false);
  histo3->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Num. Eff. and Algo.Eff.");
  leg->AddEntry(histo1,"Num. Pur.");
  leg->AddEntry(histo2,"Den. Eff.");
  leg->AddEntry(histo3,"Den. Algo.Eff.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Data_TrkTracks_Pt.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/Trackertracks eta distribution ver1");
  histo1 = (TH1F*)file0->Get("MyProcess/PixelTracks/Trackertracks eta distribution ver2");
  histo2 = (TH1F*)file0->Get("MyProcess/PixelTracks/Trackertracks eta distribution ver3");
  histo3 = (TH1F*)file0->Get("MyProcess/PixelTracks/Trackertracks eta distribution ver4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("same");

  histo3->SetLineColor(kGreen);
  histo3->SetLineWidth(2.0);
  histo3->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);
  histo2->SetStats(false);
  histo3->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Num. Eff. and Algo.Eff.");
  leg->AddEntry(histo1,"Num. Pur.");
  leg->AddEntry(histo2,"Den. Eff.");
  leg->AddEntry(histo3,"Den. Algo.Eff.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Data_TrkTracks_Eta.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/Trackertracks chi2DoF distribution ver1");
  histo1 = (TH1F*)file0->Get("MyProcess/PixelTracks/Trackertracks chi2DoF distribution ver2");
  histo2 = (TH1F*)file0->Get("MyProcess/PixelTracks/Trackertracks chi2DoF distribution ver3");
  histo3 = (TH1F*)file0->Get("MyProcess/PixelTracks/Trackertracks chi2DoF distribution ver4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("same");

  histo3->SetLineColor(kGreen);
  histo3->SetLineWidth(2.0);
  histo3->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);
  histo2->SetStats(false);
  histo3->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Num. Eff. and Algo.Eff.");
  leg->AddEntry(histo1,"Num. Pur.");
  leg->AddEntry(histo2,"Den. Eff.");
  leg->AddEntry(histo3,"Den. Algo.Eff.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Data_TrkTracks_chi2DoF.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/Trackertracks phi distribution ver1");
  histo1 = (TH1F*)file0->Get("MyProcess/PixelTracks/Trackertracks phi distribution ver2");
  histo2 = (TH1F*)file0->Get("MyProcess/PixelTracks/Trackertracks phi distribution ver3");
  histo3 = (TH1F*)file0->Get("MyProcess/PixelTracks/Trackertracks phi distribution ver4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("same");

  histo3->SetLineColor(kGreen);
  histo3->SetLineWidth(2.0);
  histo3->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);
  histo2->SetStats(false);
  histo3->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Num. Eff. and Algo.Eff.");
  leg->AddEntry(histo1,"Num. Pur.");
  leg->AddEntry(histo2,"Den. Eff.");
  leg->AddEntry(histo3,"Den. Algo.Eff.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Data_TrkTracks_Phi.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/PixelTracks/Trackertracks pt distribution ver1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/Trackertracks pt distribution ver2");
  histo2 = (TH1F*)file1->Get("MyProcess/PixelTracks/Trackertracks pt distribution ver3");
  histo3 = (TH1F*)file1->Get("MyProcess/PixelTracks/Trackertracks pt distribution ver4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("same");

  histo3->SetLineColor(kGreen);
  histo3->SetLineWidth(2.0);
  histo3->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);
  histo2->SetStats(false);
  histo3->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Num. Eff. and Algo.Eff.");
  leg->AddEntry(histo1,"Num. Pur.");
  leg->AddEntry(histo2,"Den. Eff.");
  leg->AddEntry(histo3,"Den. Algo.Eff.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_Reco_TrkTracks_Pt.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/PixelTracks/Trackertracks eta distribution ver1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/Trackertracks eta distribution ver2");
  histo2 = (TH1F*)file1->Get("MyProcess/PixelTracks/Trackertracks eta distribution ver3");
  histo3 = (TH1F*)file1->Get("MyProcess/PixelTracks/Trackertracks eta distribution ver4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("same");

  histo3->SetLineColor(kGreen);
  histo3->SetLineWidth(2.0);
  histo3->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);
  histo2->SetStats(false);
  histo3->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Num. Eff. and Algo.Eff.");
  leg->AddEntry(histo1,"Num. Pur.");
  leg->AddEntry(histo2,"Den. Eff.");
  leg->AddEntry(histo3,"Den. Algo.Eff.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_Reco_TrkTracks_Eta.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/PixelTracks/Trackertracks chi2DoF distribution ver1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/Trackertracks chi2DoF distribution ver2");
  histo2 = (TH1F*)file1->Get("MyProcess/PixelTracks/Trackertracks chi2DoF distribution ver3");
  histo3 = (TH1F*)file1->Get("MyProcess/PixelTracks/Trackertracks chi2DoF distribution ver4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("same");

  histo3->SetLineColor(kGreen);
  histo3->SetLineWidth(2.0);
  histo3->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);
  histo2->SetStats(false);
  histo3->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Num. Eff. and Algo.Eff.");
  leg->AddEntry(histo1,"Num. Pur.");
  leg->AddEntry(histo2,"Den. Eff.");
  leg->AddEntry(histo3,"Den. Algo.Eff.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_Reco_TrkTracks_chi2DoF.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/PixelTracks/Trackertracks phi distribution ver1");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/Trackertracks phi distribution ver2");
  histo2 = (TH1F*)file1->Get("MyProcess/PixelTracks/Trackertracks phi distribution ver3");
  histo3 = (TH1F*)file1->Get("MyProcess/PixelTracks/Trackertracks phi distribution ver4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("same");

  histo3->SetLineColor(kGreen);
  histo3->SetLineWidth(2.0);
  histo3->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);
  histo2->SetStats(false);
  histo3->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Num. Eff. and Algo.Eff.");
  leg->AddEntry(histo1,"Num. Pur.");
  leg->AddEntry(histo2,"Den. Eff.");
  leg->AddEntry(histo3,"Den. Algo.Eff.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_Reco_TrkTracks_Phi.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/PixelTracks/Trackertracks pt distribution ver1");
  histo1 = (TH1F*)file2->Get("MyProcess/PixelTracks/Trackertracks pt distribution ver2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/Trackertracks pt distribution ver3");
  histo3 = (TH1F*)file2->Get("MyProcess/PixelTracks/Trackertracks pt distribution ver4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("same");

  histo3->SetLineColor(kGreen);
  histo3->SetLineWidth(2.0);
  histo3->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);
  histo2->SetStats(false);
  histo3->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Num. Eff. and Algo.Eff.");
  leg->AddEntry(histo1,"Num. Pur.");
  leg->AddEntry(histo2,"Den. Eff.");
  leg->AddEntry(histo3,"Den. Algo.Eff.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_TrkPar_TrkTracks_Pt.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/PixelTracks/Trackertracks eta distribution ver1");
  histo1 = (TH1F*)file2->Get("MyProcess/PixelTracks/Trackertracks eta distribution ver2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/Trackertracks eta distribution ver3");
  histo3 = (TH1F*)file2->Get("MyProcess/PixelTracks/Trackertracks eta distribution ver4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("same");

  histo3->SetLineColor(kGreen);
  histo3->SetLineWidth(2.0);
  histo3->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);
  histo2->SetStats(false);
  histo3->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Num. Eff. and Algo.Eff.");
  leg->AddEntry(histo1,"Num. Pur.");
  leg->AddEntry(histo2,"Den. Eff.");
  leg->AddEntry(histo3,"Den. Algo.Eff.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_TrkPar_TrkTracks_Eta.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/PixelTracks/Trackertracks phi distribution ver1");
  histo1 = (TH1F*)file2->Get("MyProcess/PixelTracks/Trackertracks phi distribution ver2");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelTracks/Trackertracks phi distribution ver3");
  histo3 = (TH1F*)file2->Get("MyProcess/PixelTracks/Trackertracks phi distribution ver4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlack);
  histo0->SetLineWidth(2.0);
  histo0->Draw();

  histo1->SetLineColor(kRed);
  histo1->SetLineWidth(2.0);
  histo1->Draw("same");

  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2.0);
  histo2->Draw("same");

  histo3->SetLineColor(kGreen);
  histo3->SetLineWidth(2.0);
  histo3->Draw("same");

  histo0->SetStats(false);
  histo1->SetStats(false);
  histo2->SetStats(false);
  histo3->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"Num. Eff. and Algo.Eff.");
  leg->AddEntry(histo1,"Num. Pur.");
  leg->AddEntry(histo2,"Den. Eff.");
  leg->AddEntry(histo3,"Den. Algo.Eff.");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_TrkPar_TrkTracks_Phi.png");

// // ########################################
// // # CHANGE CHARACTERISTICS OF ALL HISTOs #
// // ########################################
//   histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi algo");
//   histo1 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Matched pixeltracks distribution in eta and phi eff");
//   histo2 = (TH1F*)histo1->Clone("hc1");

//   itC++;
//   CanvasName.str("");
//   CanvasName << "c" << itC;
//   TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());

//   histo2->Divide(histo0);
//   histo2->GetZaxis()->SetRangeUser(0.,1.);

//   histo2->Draw("gcolz");

//   histo2->SetStats(false);

// // ###################
// // # UPDATE AND SAVE #
// // ###################
//   CMS(0.);
//   c[itC]->Update();
//   c[itC]->SaveAs("Data_AlgoEff_eta_phi.png");

// // ########################################
// // # CHANGE CHARACTERISTICS OF ALL HISTOs #
// // ########################################
//   histo0 = (TH1F*)file1->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi algo");
//   histo1 = (TH1F*)file1->Get("MyProcess/PixelTracks/CommonPxTkPlots/Matched pixeltracks distribution in eta and phi eff");
//   histo2 = (TH1F*)histo1->Clone("hc2");

//   itC++;
//   CanvasName.str("");
//   CanvasName << "c" << itC;
//   TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());

//   histo2->Divide(histo0);
//   histo2->GetZaxis()->SetRangeUser(0.,1.);

//   histo2->Draw("gcolz");

//   histo2->SetStats(false);

// // ###################
// // # UPDATE AND SAVE #
// // ###################
//   CMS(0.);
//   c[itC]->Update();
//   c[itC]->SaveAs("MC_Reco_AlgoEff_eta_phi.png");

// // ########################################
// // # CHANGE CHARACTERISTICS OF ALL HISTOs #
// // ########################################
//   histo0 = (TH1F*)file0->Get("MyProcess/PixelTracks/Trackertracks distribution in eta and phi");
//   histo1 = (TH1F*)file0->Get("MyProcess/PixelTracks/CommonPxTkPlots/Matched pixeltracks distribution in eta and phi eff");
//   histo2 = (TH1F*)histo1->Clone("hc3");

//   itC++;
//   CanvasName.str("");
//   CanvasName << "c" << itC;
//   TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());

//   histo2->Divide(histo0);
//   histo2->GetZaxis()->SetRangeUser(0.,1.);

//   histo2->Draw("gcolz");

//   histo2->SetStats(false);

// // ###################
// // # UPDATE AND SAVE #
// // ###################
//   CMS(0.);
//   c[itC]->Update();
//   c[itC]->SaveAs("Data_Eff_eta_phi.png");
}
