#include "/afs/cern.ch/user/s/syu/scripts/chi2Nbins.h"
#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"

void displayThreefiles()

{
  setTDRStyle();
  gStyle->SetOptStat(0);

  std::string var1 = "h_EB_comb3IsoSB_EGdata_pt20";
  std::string var2 = "h_EB_comb3IsoSB_EGdata_pt20";
  std::string var3 = "h_EB_comb3Iso_bkg_pt20";

//   std::string xtitle = "Iso_{ECAL}+Iso_{HCAL}+Iso_{TRK} (GeV)";
  std::string xtitle = "Iso";
  std::string ytitle = "A.U.";

  TH1F* h1;
  TH1F* h2;
  TH1F* h3;

  char tempName[300];

  TCanvas* c1 = new TCanvas("c1","",500,500);

  TFile *f1 = TFile::Open("Integral_SBDataTemplate_131511_139239_362_PMMC_proj_comb.root");
  TFile *f2 = TFile::Open("Intergral_SBMCTemplate_131511_139239_362_PMMC_proj_comb.root");
  TFile *f3 = TFile::Open("Intergral_SBMCTemplate_131511_139239_362_PMMC_proj_comb.root");

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

  h1->GetYaxis()->SetNdivisions(5);
  h2->GetYaxis()->SetNdivisions(5);
  h3->GetYaxis()->SetNdivisions(5);

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
  h2->SetLineWidth(3);

  h3->SetLineColor(kGreen-3);
  h3->SetMarkerColor(kGreen-3);
  h3->SetMarkerSize(1);
  h3->SetMarkerStyle(26);
  h3->SetLineStyle(2);
  h3->SetLineWidth(3);


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


//   h2->Draw("histe");
//   h1->Draw("histe,same");
//   h3->Draw("histe,same");

  h1->Draw("histe");
  h2->Draw("hist,same");
  h3->Draw("hist,same");


  TLegend* leg = new TLegend(0.40,0.27,0.60,0.42);
  //    leg->SetHeader("p_{T}[50,80]");
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->AddEntry(h1,"Data sideband");
  leg->AddEntry(h2,"MC sideband","l");   
  leg->AddEntry(h3,"MC background truth","l");   
  leg->Draw("same");

  std::string filename;
  std::string psname = "Iso_datadriven_background";
  filename = psname + ".eps";
  c1->Print(filename.data());
  filename = psname + ".gif";
  c1->Print(filename.data());

  double chi2;
  int nbins;

  chi2NbinsCompare(h1,h2,chi2,nbins);

  cout << "chi2 test for sideband data and sideband MC = " << 
    TMath::Prob(chi2,nbins) << endl;
  cout << "KS test for sideband data and sideband MC = " << 
    h2->KolmogorovTest(h1,"X") << endl;
  cout << "KS test 2 for sideband data and sideband MC = " << 
    h2->KolmogorovTest(h1) << endl;


  chi2NbinsCompare(h2,h3,chi2,nbins);

  cout << "chi2 test for sideband MC and background MC = " << 
    TMath::Prob(chi2,nbins) << endl;

  cout << "KS test for sideband MC and background MC = " << 
    h2->KolmogorovTest(h3,"X") << endl;

  cout << "KS test2 for sideband MC and background MC = " << 
    h2->KolmogorovTest(h3) << endl;


}
		     
