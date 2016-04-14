// #########################
// # MyHistoStyle_Vertex.C #
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

//   gROOT->SetStyle("Plain");
//   gStyle->SetOptStat(111110);
//   gStyle->SetStatX(0.9);
//   gStyle->SetStatY(0.9);
//   gStyle->SetStatW(0.2);
//   gStyle->SetStatH(0.1);
//   gStyle->SetOptFit(1);
//   gStyle->SetOptTitle(false);

  gROOT->ForceStyle();

// ################
// # FILE READING #
// ################
  TFile* file0 = TFile::Open("MyPixAnalysisData_Vx_Merged.root");
  TFile* file1 = TFile::Open("MyPixAnalysisMC_noTP_Vx_Merged.root");
  TFile* file2 = TFile::Open("MyPixAnalysisMC_Vx_Merged.root");

// ###################
// # Histo Variables #
// ###################
  TH1F* histo0;
  TH1F* histo1;
  TH1F* histo2;

  TH2F* histo02;
  TH2F* histo12;
  TH2F* histo22;

  TPaveStats* PS;

  TCanvas* c[200];
  stringstream CanvasName;
  unsigned int itC = -1;

  bool WithMC = true;

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position difference along X vs Ntrks");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position difference along X vs Ntrks");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position difference along X vs Ntrks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,300.);
  histo0->GetXaxis()->SetRangeUser(0.,30.);
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
  c[itC]->SaveAs("Dif_X.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position difference along Y vs Ntrks");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position difference along Y vs Ntrks");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position difference along Y vs Ntrks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,300.);
  histo0->GetXaxis()->SetRangeUser(0.,30.);
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
  c[itC]->SaveAs("Dif_Y.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position difference along Z vs Ntrks");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position difference along Z vs Ntrks");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position difference along Z vs Ntrks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,300.);
  histo0->GetXaxis()->SetRangeUser(0.,30.);
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
  c[itC]->SaveAs("Dif_Z.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position pull along X vs Ntrks");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position pull along X vs Ntrks");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position pull along X vs Ntrks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,2.);
  histo0->GetXaxis()->SetRangeUser(0.,30.);
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
  c[itC]->SaveAs("Pull_X.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position pull along Y vs Ntrks");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position pull along Y vs Ntrks");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position pull along Y vs Ntrks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,2.);
  histo0->GetXaxis()->SetRangeUser(0.,30.);
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
  c[itC]->SaveAs("Pull_Y.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position pull along Z vs Ntrks");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position pull along Z vs Ntrks");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex position pull along Z vs Ntrks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,2.);
  histo0->GetXaxis()->SetRangeUser(0.,30.);
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
  c[itC]->SaveAs("Pull_Z.png");

  tdrStyle->SetStatW(0.12);
  tdrStyle->SetStatH(0.07);

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/VxNtrks/x_diff_Ntrk_20");
  histo1 = (TH1F*)file2->Get("MyProcess/VxNtrks/x_diff_Ntrk_10");
  histo2 = (TH1F*)file2->Get("MyProcess/VxNtrks/x_diff_Ntrk_4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlue);
  histo0->SetLineWidth(2.0);
  histo0->SetXTitle("Res X (\\mum)");
  histo0->SetYTitle("Entries (#)");
  histo0->GetYaxis()->SetRangeUser(0.,28000.);
  histo0->Draw();
  
  histo1->SetLineColor(kBlack);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kRed);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");

  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"");
  leg->AddEntry(histo0,"Ntrks 20");
  leg->AddEntry(histo1,"Ntrks 10");
  leg->AddEntry(histo2,"Ntrks  4");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlue);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.45);
  PS->SetY2NDC(0.65);
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.25);
  PS->SetY2NDC(0.45);
  PS->SetTextColor(kRed);
  c[itC]->SaveAs("MC_TrkPar_X_Ntrks.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/VxNtrks/y_diff_Ntrk_20");
  histo1 = (TH1F*)file2->Get("MyProcess/VxNtrks/y_diff_Ntrk_10");
  histo2 = (TH1F*)file2->Get("MyProcess/VxNtrks/y_diff_Ntrk_4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlue);
  histo0->SetLineWidth(2.0);
  histo0->SetXTitle("Res Y (\\mum)");
  histo0->SetYTitle("Entries (#)");
  histo0->GetYaxis()->SetRangeUser(0.,28000.);
  histo0->Draw();
  
  histo1->SetLineColor(kBlack);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kRed);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");

  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"");
  leg->AddEntry(histo0,"Ntrks 20");
  leg->AddEntry(histo1,"Ntrks 10");
  leg->AddEntry(histo2,"Ntrks  4");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlue);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.45);
  PS->SetY2NDC(0.65);
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.25);
  PS->SetY2NDC(0.45);
  PS->SetTextColor(kRed);
  c[itC]->SaveAs("MC_TrkPar_Y_Ntrks.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/VxNtrks/z_diff_Ntrk_20");
  histo1 = (TH1F*)file2->Get("MyProcess/VxNtrks/z_diff_Ntrk_10");
  histo2 = (TH1F*)file2->Get("MyProcess/VxNtrks/z_diff_Ntrk_4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlue);
  histo0->SetLineWidth(2.0);
  histo0->SetXTitle("Res Z (cm)");
  histo0->SetYTitle("Entries (#)");
  histo0->GetYaxis()->SetRangeUser(0.,34000.);
  histo0->Draw();
  
  histo1->SetLineColor(kBlack);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kRed);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");

  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"");
  leg->AddEntry(histo0,"Ntrks 20");
  leg->AddEntry(histo1,"Ntrks 10");
  leg->AddEntry(histo2,"Ntrks  4");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlue);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.45);
  PS->SetY2NDC(0.65);
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.25);
  PS->SetY2NDC(0.45);
  PS->SetTextColor(kRed);
  c[itC]->SaveAs("MC_TrkPar_Z_Ntrks.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/VxNtrks/x_diff_Ntrk_20");
  histo1 = (TH1F*)file1->Get("MyProcess/VxNtrks/x_diff_Ntrk_10");
  histo2 = (TH1F*)file1->Get("MyProcess/VxNtrks/x_diff_Ntrk_4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlue);
  histo0->SetLineWidth(2.0);
  histo0->SetXTitle("Res X (\\mum)");
  histo0->SetYTitle("Entries (#)");
  histo0->GetYaxis()->SetRangeUser(0.,25000.);
  histo0->Draw();
  
  histo1->SetLineColor(kBlack);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kRed);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");

  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"");
  leg->AddEntry(histo0,"Ntrks 20");
  leg->AddEntry(histo1,"Ntrks 10");
  leg->AddEntry(histo2,"Ntrks  4");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlue);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.45);
  PS->SetY2NDC(0.65);
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.25);
  PS->SetY2NDC(0.45);
  PS->SetTextColor(kRed);
  c[itC]->SaveAs("MC_Reco_X_Ntrks.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/VxNtrks/y_diff_Ntrk_20");
  histo1 = (TH1F*)file1->Get("MyProcess/VxNtrks/y_diff_Ntrk_10");
  histo2 = (TH1F*)file1->Get("MyProcess/VxNtrks/y_diff_Ntrk_4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());

  histo0->SetLineColor(kBlue);
  histo0->SetLineWidth(2.0);
  histo0->SetXTitle("Res Y (\\mum)");
  histo0->SetYTitle("Entries (#)");
  histo0->GetYaxis()->SetRangeUser(0.,20000.);
  histo0->Draw();
  
  histo1->SetLineColor(kBlack);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kRed);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");

  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"");
  leg->AddEntry(histo0,"Ntrks 20");
  leg->AddEntry(histo1,"Ntrks 10");
  leg->AddEntry(histo2,"Ntrks  4");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlue);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.45);
  PS->SetY2NDC(0.65);
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.25);
  PS->SetY2NDC(0.45);
  PS->SetTextColor(kRed);
  c[itC]->SaveAs("MC_Reco_Y_Ntrks.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/VxNtrks/z_diff_Ntrk_20");
  histo1 = (TH1F*)file1->Get("MyProcess/VxNtrks/z_diff_Ntrk_10");
  histo2 = (TH1F*)file1->Get("MyProcess/VxNtrks/z_diff_Ntrk_4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlue);
  histo0->SetLineWidth(2.0);
  histo0->SetXTitle("Res Z (cm)");
  histo0->SetYTitle("Entries (#)");
  histo0->GetYaxis()->SetRangeUser(0.,20000.);
  histo0->Draw();
  
  histo1->SetLineColor(kBlack);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kRed);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");

  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"");
  leg->AddEntry(histo0,"Ntrks 20");
  leg->AddEntry(histo1,"Ntrks 10");
  leg->AddEntry(histo2,"Ntrks  4");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlue);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.45);
  PS->SetY2NDC(0.65);
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.25);
  PS->SetY2NDC(0.45);
  PS->SetTextColor(kRed);
  c[itC]->SaveAs("MC_Reco_Z_Ntrks.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/VxNtrks/x_diff_Ntrk_20");
  histo1 = (TH1F*)file0->Get("MyProcess/VxNtrks/x_diff_Ntrk_10");
  histo2 = (TH1F*)file0->Get("MyProcess/VxNtrks/x_diff_Ntrk_4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlue);
  histo0->SetLineWidth(2.0);
  histo0->SetXTitle("Res X (\\mum)");
  histo0->SetYTitle("Entries (#)");
  histo0->GetYaxis()->SetRangeUser(0.,20000.);
  histo0->Draw();
  
  histo1->SetLineColor(kBlack);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kRed);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");

  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"");
  leg->AddEntry(histo0,"Ntrks 20");
  leg->AddEntry(histo1,"Ntrks 10");
  leg->AddEntry(histo2,"Ntrks  4");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlue);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.45);
  PS->SetY2NDC(0.65);
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.25);
  PS->SetY2NDC(0.45);
  PS->SetTextColor(kRed);
  c[itC]->SaveAs("Data_X_Ntrks.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/VxNtrks/y_diff_Ntrk_20");
  histo1 = (TH1F*)file0->Get("MyProcess/VxNtrks/y_diff_Ntrk_10");
  histo2 = (TH1F*)file0->Get("MyProcess/VxNtrks/y_diff_Ntrk_4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlue);
  histo0->SetLineWidth(2.0);
  histo0->SetXTitle("Res Y (\\mum)");
  histo0->SetYTitle("Entries (#)");
  histo0->GetYaxis()->SetRangeUser(0.,20000.);
  histo0->Draw();
  
  histo1->SetLineColor(kBlack);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kRed);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");

  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"");
  leg->AddEntry(histo0,"Ntrks 20");
  leg->AddEntry(histo1,"Ntrks 10");
  leg->AddEntry(histo2,"Ntrks  4");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlue);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.45);
  PS->SetY2NDC(0.65);
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.25);
  PS->SetY2NDC(0.45);
  PS->SetTextColor(kRed);
  c[itC]->SaveAs("Data_Y_Ntrks.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/VxNtrks/z_diff_Ntrk_20");
  histo1 = (TH1F*)file0->Get("MyProcess/VxNtrks/z_diff_Ntrk_10");
  histo2 = (TH1F*)file0->Get("MyProcess/VxNtrks/z_diff_Ntrk_4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlue);
  histo0->SetLineWidth(2.0);
  histo0->SetXTitle("Res Z (cm)");
  histo0->SetYTitle("Entries (#)");
  histo0->GetYaxis()->SetRangeUser(0.,35000.);
  histo0->Draw();
  
  histo1->SetLineColor(kBlack);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kRed);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");

  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"");
  leg->AddEntry(histo0,"Ntrks 20");
  leg->AddEntry(histo1,"Ntrks 10");
  leg->AddEntry(histo2,"Ntrks  4");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlue);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.45);
  PS->SetY2NDC(0.65);
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.25);
  PS->SetY2NDC(0.45);
  PS->SetTextColor(kRed);
  c[itC]->SaveAs("Data_Z_Ntrks.png");

  tdrStyle->SetStatW(0.2);
  tdrStyle->SetStatH(0.1);

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelVertices/CommonPxTkPlots/Purity Efficiency Vertices");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelVertices/CommonPxTkPlots/Purity Efficiency Vertices");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelVertices/CommonPxTkPlots/Purity Efficiency Vertices");

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

  histo1->SetMarkerStyle(21);
  histo1->SetMarkerColor(kRed);
  histo1->Draw("same");

  cout << "@@@ Simulation @@@" << endl;
  cout << "Purity: " << histo1->GetBinContent(histo1->GetXaxis()->FindBin("Purity")) << " +/- " << histo1->GetBinError(histo1->GetXaxis()->FindBin("Purity")) << endl;
  cout << "Efficiency: " << histo1->GetBinContent(histo1->GetXaxis()->FindBin("Efficiency")) << " +/- " << histo1->GetBinError(histo1->GetXaxis()->FindBin("Efficiency")) << endl;

  if (WithMC == true)
    {
      histo2->SetMarkerStyle(22);
      histo2->SetMarkerColor(kBlue);
      histo2->Draw("same");
    }

  cout << "@@@ MC-Truth @@@" << endl;
  cout << "Purity: " << histo2->GetBinContent(histo2->GetXaxis()->FindBin("Purity")) << " +/- " << histo2->GetBinError(histo2->GetXaxis()->FindBin("Purity")) << endl;
  cout << "Efficiency: " << histo2->GetBinContent(histo2->GetXaxis()->FindBin("Efficiency")) << " +/- " << histo2->GetBinError(histo2->GetXaxis()->FindBin("Efficiency")) << endl;

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
  c[itC]->SaveAs("Pur_Eff_Vx.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/PixelVertices/CommonPxTkPlots/pt distribution of pixelvertices tracks");
  histo1 = (TH1F*)file2->Get("MyProcess/PixelVertices/CommonPxTkPlots/pt distribution of pixelvertices nomatch tracks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetFillColor(kBlue);
  histo0->Draw("bar");
  
  histo1->SetFillColor(kRed);
  histo1->Draw("bar,same");

  histo0->SetStats(false);
  histo1->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"All Tracks");
  leg->AddEntry(histo1,"Non Matched Tracks");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_TrkPar_Pix_pt.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/PixelVertices/CommonPxTkPlots/pt distribution of pixelvertices tracks");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelVertices/CommonPxTkPlots/pt distribution of pixelvertices nomatch tracks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetFillColor(kBlue);
  histo0->Draw("bar");
  
  histo1->SetFillColor(kRed);
  histo1->Draw("bar,same");

  histo0->SetStats(false);
  histo1->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"All Tracks");
  leg->AddEntry(histo1,"Non Matched Tracks");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_Reco_Pix_pt.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelVertices/CommonPxTkPlots/pt distribution of pixelvertices tracks");
  histo1 = (TH1F*)file0->Get("MyProcess/PixelVertices/CommonPxTkPlots/pt distribution of pixelvertices nomatch tracks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetFillColor(kBlue);
  histo0->Draw("bar");
  
  histo1->SetFillColor(kRed);
  histo1->Draw("bar,same");

  histo0->SetStats(false);
  histo1->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"All Tracks");
  leg->AddEntry(histo1,"Non Matched Tracks");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Data_Pix_pt.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/PixelVertices/CommonPxTkPlots/eta distribution of pixelvertices nomatch tracks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetFillColor(kRed);
  histo0->GetYaxis()->SetRangeUser(0.,20000.);
  histo0->Draw("bar");
  
  histo0->SetStats(true);

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_TrkPar_Pix_eta.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/PixelVertices/CommonPxTkPlots/eta distribution of pixelvertices nomatch tracks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetFillColor(kRed);
  histo0->GetYaxis()->SetRangeUser(0.,120000.);
  histo0->Draw("bar");
  
  histo0->SetStats(true);

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_Reco_Pix_eta.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelVertices/CommonPxTkPlots/eta distribution of pixelvertices nomatch tracks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetFillColor(kRed);
  histo0->GetYaxis()->SetRangeUser(0.,400000.);
  histo0->Draw("bar");
  
  histo0->SetStats(true);

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Data_Pix_eta.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/PixelVertices/CommonPxTkPlots/pt distribution of trackervertices tracks");
  histo1 = (TH1F*)file2->Get("MyProcess/PixelVertices/CommonPxTkPlots/pt distribution of trackervertices nomatch tracks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetFillColor(kBlue);
  histo0->Draw("bar");
  
  histo1->SetFillColor(kRed);
  histo1->Draw("bar,same");

  histo0->SetStats(false);
  histo1->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"All Tracks");
  leg->AddEntry(histo1,"Non Matched Tracks");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_TrkPar_Trk_pt.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/PixelVertices/CommonPxTkPlots/pt distribution of trackervertices tracks");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelVertices/CommonPxTkPlots/pt distribution of trackervertices nomatch tracks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetFillColor(kBlue);
  histo0->Draw("bar");
  
  histo1->SetFillColor(kRed);
  histo1->Draw("bar,same");

  histo0->SetStats(false);
  histo1->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"All Tracks");
  leg->AddEntry(histo1,"Non Matched Tracks");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_Reco_Trk_pt.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelVertices/CommonPxTkPlots/pt distribution of trackervertices tracks");
  histo1 = (TH1F*)file0->Get("MyProcess/PixelVertices/CommonPxTkPlots/pt distribution of trackervertices nomatch tracks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetFillColor(kBlue);
  histo0->Draw("bar");
  
  histo1->SetFillColor(kRed);
  histo1->Draw("bar,same");

  histo0->SetStats(false);
  histo1->SetStats(false);

  TLegend* leg = new TLegend(0.55,0.7,0.9,0.9,"");
  leg->AddEntry(histo0,"All Tracks");
  leg->AddEntry(histo1,"Non Matched Tracks");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Data_Trk_pt.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file2->Get("MyProcess/PixelVertices/CommonPxTkPlots/eta distribution of trackervertices nomatch tracks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetFillColor(kRed);
  histo0->GetYaxis()->SetRangeUser(0.,80000.);
  histo0->Draw("bar");
  
  histo0->SetStats(true);

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_TrkPar_Trk_eta.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/PixelVertices/CommonPxTkPlots/eta distribution of trackervertices nomatch tracks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetFillColor(kRed);
  histo0->GetYaxis()->SetRangeUser(0.,60000.);
  histo0->Draw("bar");
  
  histo0->SetStats(true);

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("MC_Reco_Trk_eta.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelVertices/CommonPxTkPlots/eta distribution of trackervertices nomatch tracks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetFillColor(kRed);
  histo0->GetYaxis()->SetRangeUser(0.,200000.);
  histo0->Draw("bar");
  
  histo0->SetStats(true);

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  c[itC]->SaveAs("Data_Trk_eta.png");

  tdrStyle->SetTitleYOffset(1.18);

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex resolution along X vs Ntrks");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex resolution along X vs Ntrks");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex resolution along X vs Ntrks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,300.);
  histo0->GetXaxis()->SetRangeUser(0.,30.);
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
  c[itC]->SaveAs("Res_X.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex resolution along Y vs Ntrks");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex resolution along Y vs Ntrks");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex resolution along Y vs Ntrks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,300.);
  histo0->GetXaxis()->SetRangeUser(0.,30.);
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
  c[itC]->SaveAs("Res_Y.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex resolution along Z vs Ntrks");
  histo1 = (TH1F*)file1->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex resolution along Z vs Ntrks");
  histo2 = (TH1F*)file2->Get("MyProcess/PixelVertices/CommonPxTkPlots/Vertex resolution along Z vs Ntrks");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetMarkerStyle(20);
  histo0->SetMarkerColor(kBlack);
  histo0->GetYaxis()->SetRangeUser(0.,300.);
  histo0->GetXaxis()->SetRangeUser(0.,30.);
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
  c[itC]->SaveAs("Res_Z.png");

  tdrStyle->SetTitleYOffset(1.3);

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/VxNtrks/x_err2_Ntrk_20");
  histo1 = (TH1F*)file1->Get("MyProcess/VxNtrks/x_err2_Ntrk_10");
  histo2 = (TH1F*)file1->Get("MyProcess/VxNtrks/x_err2_Ntrk_4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlue);
  histo0->SetLineWidth(2.0);
  histo0->GetYaxis()->SetRangeUser(0.,6000.);
  histo0->SetXTitle("\\sigma^{2}_{TrkX} (\\mum)^{2}");
  histo0->SetYTitle("Entries (#)");
  histo0->Draw();
  
  histo1->SetLineColor(kBlack);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kRed);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");

  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"");
  leg->AddEntry(histo0,"Ntrks 20");
  leg->AddEntry(histo1,"Ntrks 10");
  leg->AddEntry(histo2,"Ntrks  4");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlue);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kRed);
  c[itC]->SaveAs("MC_Reco_err2X_Ntrks.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/VxNtrks/y_err2_Ntrk_20");
  histo1 = (TH1F*)file1->Get("MyProcess/VxNtrks/y_err2_Ntrk_10");
  histo2 = (TH1F*)file1->Get("MyProcess/VxNtrks/y_err2_Ntrk_4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlue);
  histo0->SetLineWidth(2.0);
  histo0->GetYaxis()->SetRangeUser(0.,6000.);
  histo0->SetXTitle("\\sigma^{2}_{TrkY} (\\mum)^{2}");
  histo0->SetYTitle("Entries (#)");
  histo0->Draw();
  
  histo1->SetLineColor(kBlack);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kRed);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");

  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"");
  leg->AddEntry(histo0,"Ntrks 20");
  leg->AddEntry(histo1,"Ntrks 10");
  leg->AddEntry(histo2,"Ntrks  4");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlue);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kRed);
  c[itC]->SaveAs("MC_Reco_err2Y_Ntrks.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file1->Get("MyProcess/VxNtrks/z_err2_Ntrk_20");
  histo1 = (TH1F*)file1->Get("MyProcess/VxNtrks/z_err2_Ntrk_10");
  histo2 = (TH1F*)file1->Get("MyProcess/VxNtrks/z_err2_Ntrk_4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlue);
  histo0->SetLineWidth(2.0);
  histo0->GetYaxis()->SetRangeUser(0.,6000.);
  histo0->SetXTitle("\\sigma^{2}_{TrkZ} (\\mum)^{2}");
  histo0->SetYTitle("Entries (#)");
  histo0->Draw();
  
  histo1->SetLineColor(kBlack);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kRed);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");

  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"");
  leg->AddEntry(histo0,"Ntrks 20");
  leg->AddEntry(histo1,"Ntrks 10");
  leg->AddEntry(histo2,"Ntrks  4");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlue);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kRed);
  c[itC]->SaveAs("MC_Reco_err2Z_Ntrks.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/VxNtrks/x_err2_Ntrk_20");
  histo1 = (TH1F*)file0->Get("MyProcess/VxNtrks/x_err2_Ntrk_10");
  histo2 = (TH1F*)file0->Get("MyProcess/VxNtrks/x_err2_Ntrk_4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlue);
  histo0->SetLineWidth(2.0);
  histo0->GetYaxis()->SetRangeUser(0.,14000.);
  histo0->SetXTitle("\\sigma^{2}_{TrkX} (\\mum)^{2}");
  histo0->SetYTitle("Entries (#)");
  histo0->Draw();
  
  histo1->SetLineColor(kBlack);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kRed);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");

  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"");
  leg->AddEntry(histo0,"Ntrks 20");
  leg->AddEntry(histo1,"Ntrks 10");
  leg->AddEntry(histo2,"Ntrks  4");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlue);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kRed);
  c[itC]->SaveAs("Data_err2X_Ntrks.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/VxNtrks/y_err2_Ntrk_20");
  histo1 = (TH1F*)file0->Get("MyProcess/VxNtrks/y_err2_Ntrk_10");
  histo2 = (TH1F*)file0->Get("MyProcess/VxNtrks/y_err2_Ntrk_4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlue);
  histo0->SetLineWidth(2.0);
  histo0->GetYaxis()->SetRangeUser(0.,14000.);
  histo0->SetXTitle("\\sigma^{2}_{TrkY} (\\mum)^{2}");
  histo0->SetYTitle("Entries (#)");
  histo0->Draw();
  
  histo1->SetLineColor(kBlack);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kRed);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");

  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"");
  leg->AddEntry(histo0,"Ntrks 20");
  leg->AddEntry(histo1,"Ntrks 10");
  leg->AddEntry(histo2,"Ntrks  4");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlue);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kRed);
  c[itC]->SaveAs("Data_err2Y_Ntrks.png");

// ########################################
// # CHANGE CHARACTERISTICS OF ALL HISTOs #
// ########################################
  histo0 = (TH1F*)file0->Get("MyProcess/VxNtrks/z_err2_Ntrk_20");
  histo1 = (TH1F*)file0->Get("MyProcess/VxNtrks/z_err2_Ntrk_10");
  histo2 = (TH1F*)file0->Get("MyProcess/VxNtrks/z_err2_Ntrk_4");

  itC++;
  CanvasName.str("");
  CanvasName << "c" << itC;
  TCanvas* c[itC] = new TCanvas(CanvasName.str().c_str(), CanvasName.str().c_str());
  
  histo0->SetLineColor(kBlue);
  histo0->SetLineWidth(2.0);
  histo0->GetYaxis()->SetRangeUser(0.,14000.);
  histo0->SetXTitle("\\sigma^{2}_{TrkZ} (\\mum)^{2}");
  histo0->SetYTitle("Entries (#)");
  histo0->Draw();
  
  histo1->SetLineColor(kBlack);
  histo1->SetLineWidth(2.0);
  histo1->Draw("sames");

  histo2->SetLineColor(kRed);
  histo2->SetLineWidth(2.0);
  histo2->Draw("sames");

  histo0->SetStats(true);
  histo1->SetStats(true);
  histo2->SetStats(true);

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.9,"");
  leg->AddEntry(histo0,"Ntrks 20");
  leg->AddEntry(histo1,"Ntrks 10");
  leg->AddEntry(histo2,"Ntrks  4");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

// ###################
// # UPDATE AND SAVE #
// ###################
  CMS(0.);
  c[itC]->Update();
  PS = (TPaveStats*)histo0->GetListOfFunctions()->FindObject("stats");
  PS->SetTextColor(kBlue);
  PS = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.55);
  PS->SetY2NDC(0.725);
  PS->SetTextColor(kBlack);
  PS = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  PS->SetY1NDC(0.375);
  PS->SetY2NDC(0.55);
  PS->SetTextColor(kRed);
  c[itC]->SaveAs("Data_err2Z_Ntrks.png");
}
