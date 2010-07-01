#include "/afs/cern.ch/user/s/syu/scripts/chi2Nbins.h"
#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"

void displayThreefiles()

{
  setTDRStyle();
  gStyle->SetOptStat(0);

  std::string var1 = "h_EB_comb3IsoSB_EGdata_pt15";
  std::string var2 = "h_EB_comb3IsoSB_EGdata_pt15";
  std::string var3 = "h_EB_comb3Iso_bkg_pt15";

  std::string xtitle = "Iso_{ECAL}+Iso_{HCAL}+Iso_{TRK} (GeV)";
  std::string ytitle = "A.U.";

  TH1F* h1;
  TH1F* h2;
  TH1F* h3;

  char tempName[300];

  TCanvas* c1 = new TCanvas("c1","",500,500);

  TFile *f1 = TFile::Open("SBDataTemplate_proj_comb.root");
  TFile *f2 = TFile::Open("SBMCTemplate_proj_comb.root");
  TFile *f3 = TFile::Open("SBMCTemplate_proj_comb.root");

  h1  = (TH1F*)(f1->Get(var1.data()));
  h1->SetXTitle(xtitle.data());  
  h1->Rebin(12);
  h1->SetYTitle(ytitle.data());  

  h2  = (TH1F*)(f2->Get(var2.data()));
  h2->SetXTitle(xtitle.data());
  h2->Rebin(12);
  h2->SetYTitle(ytitle.data());  

  h3  = (TH1F*)(f3->Get(var3.data()));
  h3->SetXTitle(xtitle.data());
  h3->Rebin(12);
  h3->SetYTitle(ytitle.data());  

  h1->GetYaxis()->SetDecimals();
  h2->GetYaxis()->SetDecimals();
  h3->GetYaxis()->SetDecimals();

  h1->SetLineColor(1);
  h1->SetMarkerColor(1);
  h1->SetMarkerSize(1);
  h1->SetMarkerStyle(21);


  h2->SetLineColor(2);
  h2->SetMarkerColor(2);
  h2->SetMarkerSize(1);
  h2->SetMarkerStyle(23);

  h3->SetLineColor(kGreen-3);
  h3->SetMarkerColor(kGreen-3);
  h3->SetMarkerSize(1);
  h3->SetMarkerStyle(26);



  float scale1 = 1.0/(float)h1->Integral();
  float scale2 = 1.0/(float)h2->Integral();
  float scale3 = 1.0/(float)h3->Integral();

   h1->SetTitle("");
   h2->SetTitle("");
   h1->Sumw2();
   h1->Scale(scale1);
   h2->Sumw2();
   h2->Scale(scale2);
   h3->Sumw2();
   h3->Scale(scale3);

   float max1   = h1->GetBinError(h1->GetMaximumBin()) + h1->GetMaximum();
   float max2   = h2->GetBinError(h2->GetMaximumBin()) + h2->GetMaximum();


   h1->Draw("histe");
   h2->Draw("histe,same");
   h3->Draw("histe,same");


   TLegend* leg = new TLegend(0.40,0.27,0.60,0.42);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->SetTextSize(0.04);
   leg->SetBorderSize(0);
   leg->AddEntry(h1,"Data sideband");
   leg->AddEntry(h2,"MC sideband");   
   leg->AddEntry(h3,"MC background truth");   
   leg->Draw("same");

   std::string filename;
   std::string psname = "Iso_datadriven_background";
   filename = psname + ".eps";
   c1->Print(filename.data());
   filename = psname + ".gif";
   c1->Print(filename.data());
}
		     
