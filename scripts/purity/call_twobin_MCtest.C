#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TH1.h>
#include <string>
#include "purity_twobin.C"
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMath.h> 
#include <TCanvas.h>
#include <TLegend.h>
#include "setTDRStyle.C"

using namespace std;

// the pt and eta binning
const double fBinsEta[]={0,1.45,1.7,2.5};
const double fBinsPt[]={15,20,30,50,80,120};
const int nPtBin = sizeof(fBinsPt)/sizeof(fBinsPt[0])-1;
const int nEtaBin = (sizeof(fBinsEta)/sizeof(fBinsEta[0]))/2;

void call_twobin_MCtest(){
  setTDRStyle();
  
  double mcpurity_EB[8];
  double mcpurity_EE[8];
  double purity_EB[8];
  double purity_EE[8];
  double purity_err_EB[8];
  double purity_err_EE[8];

  double purity_template_EB[8];
  double purity_template_EE[8];
  double purity_err_template_EB[8];
  double purity_err_template_EE[8];
 
  double purity_twobin_EB[8]={0};
  double purity_err_twobin_EB[8]={0};

  double purity_twobin_EE[8]={0};
  double purity_err_twobin_EE[8]={0};

  double nsig_twobin_EB[8]={0};
  double nsig_err_twobin_EB[8]={0};

  double nsig_twobin_EE[8]={0};
  double nsig_err_twobin_EE[8]={0};

 

  int nbin=5;
  double mcsig_EB[8] = {4782.0, 4140.2, 1202.3, 180.9, 32.37, 5.80, 0.502, 0.053};
  double mcsig_EE[8] = {2960.7, 2356.1,  727.8, 108.6, 15.80, 2.75, 0.188, 0.028};
  
  double mcbkg_EB[8] = {32211, 18360, 3121.4, 339.2, 24.51, 2.48, 0.415, 0.0518};
  double mcbkg_EE[8] = {48219, 23910, 4464.0, 492.1, 36.12, 3.28, 0.385, 0.0466};

  double ptbin[20] = {17.5, 25, 40, 65, 100, 145, 200, 265};
  double pterr[20] = {2.5, 5., 10, 15, 20, 25, 30, 35};
  int pt[10] = {15,20,30,50,80,120,170,230};  

  //unbinned fit
  double sig_EB[8] =     { 5225.2,  4090.8, 1066.8, 212.7, 35.1, 10.6};
  double sig_err_EB[8] = {  122.1,    99.5,   45.4,  16.3,  6.5,  3.5};
  double bkg_EB[8] =     {31777.8, 18419.2, 3267.2, 317.3, 31.0, 12.4};
  double bkg_err_EB[8] = {  203.6,   155.6,   65.3,  19.3,  6.2,  3.8};

  double sig_EE[8] =     { 3049.8,  2208.0,  686.4, 136.1, 15.55,  6.12};
  double sig_err_EE[8] = {  127.0,    91.2,   40.1,  14.2,  4.89,  2.59};
  double bkg_EE[8] =     {48139.5, 24068.0, 4517.6, 473.9, 49.45, 13.88};
  double bkg_err_EE[8] = {  247.4,   173.7,   74.0,  23.2,  7.59,  3.79};

  //template fit
  double sig_template_EB[8] =     {5218.761,4185.771,1109.753,202.912,32.384, 3.8};
  double sig_err_template_EB[8] = {125.577,102.830, 46.883,16.390,5.990, 1.6};
  double bkg_template_EB[8] =     {31774.322,18314.933,3214.038,317.172,24.508, 4.5};
  double bkg_err_template_EB[8] = {205.7, 157.1, 65.6, 19.6, 5.3, 1.9};
			  
  double sig_template_EE[8] =     {2626.388,2269.472,795.176,142.129,11.441, 1.8};
  double sig_err_template_EE[8] = {131.006,94.176,44.260,14.566,4.062, 1.0};
  double bkg_template_EE[8] =     {48553.637,23996.410,4396.580, 457.519,40.528,4.2};
  double bkg_err_template_EE[8] = {251.2, 174.9, 74.6, 23.0, 6.7, 1.8};
  


  for(int i=0; i<nbin; i++) {
    mcpurity_EB[i] = mcsig_EB[i]/(mcsig_EB[i]+mcbkg_EB[i]);
    mcpurity_EE[i] = mcsig_EE[i]/(mcsig_EE[i]+mcbkg_EE[i]);

    purity_EB[i] = sig_EB[i]/(sig_EB[i]+bkg_EB[i]);
    purity_EE[i] = sig_EE[i]/(sig_EE[i]+bkg_EE[i]);

    purity_err_EB[i] = purity_EB[i]*TMath::Sqrt(1/sig_EB[i] + 1/bkg_EB[i]);
    purity_err_EE[i] = purity_EE[i]*TMath::Sqrt(1/sig_EE[i] + 1/bkg_EE[i]);

    purity_template_EB[i] = sig_template_EB[i]/(sig_template_EB[i]+bkg_template_EB[i]);
    purity_template_EE[i] = sig_template_EE[i]/(sig_template_EE[i]+bkg_template_EE[i]);

    purity_err_template_EB[i] = purity_template_EB[i]*TMath::Sqrt(1/sig_template_EB[i] + 1/bkg_template_EB[i]);
    purity_err_template_EE[i] = purity_template_EE[i]*TMath::Sqrt(1/sig_template_EE[i] + 1/bkg_template_EE[i]);
    cout << "MC truth purity[" << i << "] = "  <<  
      mcpurity_EB[i] << "\t" << mcpurity_EE[i] << endl;
  }    

  // first try two bin results
  TH1F* hTemplate_S[2][5];
  TH1F* hTemplate_B[2][5];
  TH1F* hTemplate_data[2][5];
  TH1F* hTemplate_data_template[2][5];


  TFile *inf_data = new TFile("/mc/QCD_mess/pcncu18_figures/fakedata_sumIso_histo.root");

  TFile *inf_template = new TFile("/mc/QCD_mess/pcncu18_figures/template_sumIso_histo.root");

//   TFile *inf_template = new TFile("/mc/QCD_mess/pcncu18_figures/fakedata_sumIso_histo.root");
  
  char tmp[1000];
  for(int ieta = 0; ieta < nEtaBin; ieta++){
     for(int ipt=0; ipt < nPtBin; ipt++){
       sprintf(tmp,"Template_Et_%d_%d_Eta_%.2f_%.2f",
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
       cout << "looking for histogram " << tmp << " in file " << 
	 inf_data->GetName() << endl;
       hTemplate_data[ieta][ipt] = (TH1F*)inf_data->FindObjectAny(tmp);

       sprintf(tmp,"Template_Et_%d_%d_Eta_%.2f_%.2f",
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
       cout << "looking for histogram " << tmp << " in file " << 
	 inf_template->GetName() << endl;
       hTemplate_data_template[ieta][ipt] = (TH1F*)inf_template->FindObjectAny(tmp);

       sprintf(tmp,"TemplateS_Et_%d_%d_Eta_%.2f_%.2f",
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
       cout << "looking for histogram " << tmp << " in file " << 
	 inf_template->GetName() << endl;
       hTemplate_S[ieta][ipt] = (TH1F*)inf_template->FindObjectAny(tmp);

       sprintf(tmp,"TemplateB_Et_%d_%d_Eta_%.2f_%.2f",
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
       cout << "looking for histogram " << tmp << " in file " << 
	 inf_template->GetName() << endl;
       hTemplate_B[ieta][ipt] = (TH1F*)inf_template->FindObjectAny(tmp);

     }
  }


  for(int ieta = 0; ieta < nEtaBin; ieta++){
     for(int ipt=0; ipt < nPtBin; ipt++){
       
       
       double tempsig=0, temperr=0;
       purity_twobin(hTemplate_data[ieta][ipt],
		     hTemplate_S[ieta][ipt],
		     hTemplate_B[ieta][ipt],"test",
		     tempsig,temperr);
       
       if(ieta==0){
	 purity_twobin_EB[ipt]=tempsig;
	 purity_err_twobin_EB[ipt]=temperr; 
	 nsig_twobin_EB[ipt]=tempsig*(hTemplate_data[ieta][ipt]->Integral()+
				      hTemplate_data_template[ieta][ipt]->Integral());
	 nsig_err_twobin_EB[ipt]=temperr*(hTemplate_data[ieta][ipt]->Integral()+
					  hTemplate_data_template[ieta][ipt]->Integral());
       }
	 
       else if(ieta==1){
	 purity_twobin_EE[ipt] = tempsig;
	 purity_err_twobin_EE[ipt] = temperr;
	 nsig_twobin_EE[ipt]=tempsig*(hTemplate_data[ieta][ipt]->Integral()+
				      hTemplate_data_template[ieta][ipt]->Integral());
	 nsig_err_twobin_EE[ipt]=temperr*(hTemplate_data[ieta][ipt]->Integral()+
					  hTemplate_data_template[ieta][ipt]->Integral());
       }

     }
  }

  

//   TH1F* hpurity_EB = (TH1F*)inf_data->FindObjectAny("hTruthPurity_Eta_0.00_1.45");
//   TH1F* hpurity_EE = (TH1F*)inf_data->FindObjectAny("hTruthPurity_Eta_1.70_2.50");
//   for(int i=1;i<=hpurity_EB->GetNbinsX(); i++)
//     {

//        mcpurity_EB[i-1] = hpurity_EB->GetBinContent(i);
//        mcpurity_EE[i-1] = hpurity_EE->GetBinContent(i);

//      }
  
  
  TCanvas *c1 = new TCanvas("c1","",600,600);
  c1->Draw();
  TH2F *h2 = new TH2F("h2","",100, 0., 140, 100, 0, 1);
  h2->SetNdivisions(505,"XY");
  h2->SetXTitle("p_{T} (GeV)");
  h2->SetYTitle("purity");
  h2->Draw();
  TGraph *tg_mc_purity_EB = new TGraph(nbin, ptbin, mcpurity_EB);
  tg_mc_purity_EB->SetMarkerSize(1.);
  tg_mc_purity_EB->Draw("p same");
  TGraphErrors *tgrs_purity_EB = new TGraphErrors(nbin, ptbin, purity_EB, pterr, purity_err_EB);
  tgrs_purity_EB->SetMarkerColor(4);
  tgrs_purity_EB->SetMarkerStyle(25);
  tgrs_purity_EB->SetLineColor(4);
  tgrs_purity_EB->SetLineWidth(2);
  tgrs_purity_EB->Draw("p e same");

  TGraphErrors *tgrs_purity_template_EB = new TGraphErrors(nbin, ptbin, purity_template_EB, pterr, purity_err_template_EB);
  tgrs_purity_template_EB->SetMarkerColor(6);
  tgrs_purity_template_EB->SetMarkerStyle(26);
  tgrs_purity_template_EB->SetLineColor(6);
  tgrs_purity_template_EB->SetLineWidth(2);
  tgrs_purity_template_EB->Draw("p e same");

  TGraphErrors *tgrs_purity_twobin_EB = new TGraphErrors(nbin, ptbin, purity_twobin_EB, pterr, purity_err_twobin_EB);
  tgrs_purity_twobin_EB->SetMarkerColor(kGreen);
  tgrs_purity_twobin_EB->SetMarkerStyle(3 );
  tgrs_purity_twobin_EB->SetLineColor(kGreen);
  tgrs_purity_twobin_EB->SetLineWidth(2);
  tgrs_purity_twobin_EB->Draw("p e same");

  tgrs_purity_template_EB->Draw("p e same");
  tg_mc_purity_EB->Draw("p same");

  TLegend *tleg = new TLegend(0.15, 0.65, 0.55, 0.92);
  char text[50];  
  tleg->SetHeader("EB Photons");
  tleg->SetFillColor(0);
  tleg->SetShadowColor(0);
  tleg->SetBorderSize(0);
  tleg->AddEntry(tg_mc_purity_EB,"MC truth","p");
  tleg->AddEntry(tgrs_purity_EB,"unbinned fit","pl");
  tleg->AddEntry(tgrs_purity_template_EB,"template fit","pl");
  tleg->AddEntry(tgrs_purity_twobin_EB,"two bin","pl");
  tleg->Draw();

  //  c1->SaveAs("twobin_test/allfitter_purity_EB.pdf");
  //  c1->SaveAs("twobin_test/allfitter_purity_EB.C");
  

  TCanvas *c2 = new TCanvas("c2","",600,600);
  c2->Draw();
//   TH2F *h2 = new TH2F("h2","",100, 0., 140, 100, 0, 1);
//   h2->SetNdivisions(505,"XY");
//   h2->SetXTitle("p_{T} (GeV)");
//   h2->SetYTitle("purity");

  h2->Draw();
  TGraph *tg_mc_purity_EE = new TGraph(nbin, ptbin, mcpurity_EE);
  tg_mc_purity_EE->SetMarkerSize(1.);
  tg_mc_purity_EE->Draw("p same");
  TGraphErrors *tgrs_purity_EE = new TGraphErrors(nbin, ptbin, purity_EE, pterr, purity_err_EE);
  tgrs_purity_EE->SetMarkerColor(4);
  tgrs_purity_EE->SetMarkerStyle(25);
  tgrs_purity_EE->SetLineColor(4);
  tgrs_purity_EE->SetLineWidth(2);
  tgrs_purity_EE->Draw("p e same");

  TGraphErrors *tgrs_purity_template_EE = new TGraphErrors(nbin, ptbin, purity_template_EE, pterr, purity_err_template_EE);
  tgrs_purity_template_EE->SetMarkerColor(6);
  tgrs_purity_template_EE->SetMarkerStyle(26);
  tgrs_purity_template_EE->SetLineColor(6);
  tgrs_purity_template_EE->SetLineWidth(2);
  tgrs_purity_template_EE->Draw("p e same");

  TGraphErrors *tgrs_purity_twobin_EE = new TGraphErrors(nbin, ptbin, purity_twobin_EE, pterr, purity_err_twobin_EE);
  tgrs_purity_twobin_EE->SetMarkerColor(kGreen);
  tgrs_purity_twobin_EE->SetMarkerStyle(3);
  tgrs_purity_twobin_EE->SetLineColor(kGreen);
  tgrs_purity_twobin_EE->SetLineWidth(2);
  tgrs_purity_twobin_EE->Draw("p e same");

  tgrs_purity_template_EE->Draw("p e same");
  tg_mc_purity_EE->Draw("p same");

  TLegend *tleg = new TLegend(0.15, 0.65, 0.55, 0.92);
  char text[50];  
  tleg->SetHeader("EE Photons");
  tleg->SetFillColor(0);
  tleg->SetShadowColor(0);
  tleg->SetBorderSize(0);
  tleg->AddEntry(tg_mc_purity_EE,"MC truth","p");
  tleg->AddEntry(tgrs_purity_EE,"unbinned fit","pl");
  tleg->AddEntry(tgrs_purity_template_EE,"template fit","pl");
  tleg->AddEntry(tgrs_purity_twobin_EE,"two bin","pl");
  tleg->Draw();

  c2->SaveAs("twobin_test/allfitter_purity_EE.pdf");
  c2->SaveAs("twobin_test/allfitter_purity_EE.C");


  TCanvas *c3 = new TCanvas("c3","",600,600);
  c3->Draw();
  gPad->SetLogy();
  //TH2F *h3 = new TH2F("h3","",100, 0., 140, 100, 0., 6000);  
  TH2F *h3 = new TH2F("h3","",100, 0., 140, 100, 0.1, 10000);  
  h3->SetNdivisions(505,"XY");
  h3->SetXTitle("p_{T} (GeV)");
  h3->SetYTitle("Signal yield / bin");
  h3->Draw();
  TGraph *tg_mc_sig_EB = new TGraph(nbin, ptbin, mcsig_EB);
  tg_mc_sig_EB->SetMarkerSize(1);
  tg_mc_sig_EB->Draw("p same");
  TGraphErrors *tgrs_sig_EB = new TGraphErrors(nbin, ptbin, sig_EB, pterr, sig_err_EB);
  tgrs_sig_EB->SetMarkerColor(4);
  tgrs_sig_EB->SetMarkerStyle(25);
  tgrs_sig_EB->SetLineColor(4);
  tgrs_sig_EB->SetLineWidth(2);
  tgrs_sig_EB->Draw("p e same");

  TGraphErrors *tgrs_sig_template_EB = new TGraphErrors(nbin, ptbin, sig_template_EB, pterr, sig_err_template_EB);
  tgrs_sig_template_EB->SetMarkerColor(6);
  tgrs_sig_template_EB->SetMarkerStyle(26);
  tgrs_sig_template_EB->SetLineColor(6);
  tgrs_sig_template_EB->SetLineWidth(2);
  tgrs_sig_template_EB->Draw("p e same");

  TGraphErrors *tgrs_sig_twobin_EB = new TGraphErrors(nbin, ptbin, nsig_twobin_EB, pterr, nsig_err_twobin_EB);
  tgrs_sig_twobin_EB->SetMarkerColor(kGreen);
  tgrs_sig_twobin_EB->SetMarkerStyle(3);
  tgrs_sig_twobin_EB->SetLineColor(kGreen);
  tgrs_sig_twobin_EB->SetLineWidth(2);
  tgrs_sig_twobin_EB->Draw("p e same");



  TLegend *tleg = new TLegend(0.193, 0.168, 0.642, 0.437);
  char text[50];  
  tleg->SetHeader("EB Photons");
  tleg->SetFillColor(0);
  tleg->SetShadowColor(0);
  tleg->SetBorderSize(0);
  tleg->AddEntry(tg_mc_sig_EB,"MC truth","p");
  tleg->AddEntry(tgrs_sig_EB,"unbinned fit","pl");
  tleg->AddEntry(tgrs_sig_template_EB,"template fit","pl");
  tleg->AddEntry(tgrs_sig_twobin_EB,"two bin","pl");
  tleg->Draw();

  c3->SaveAs("twobin_test/allfitter_yield_EB.pdf");

  TCanvas *c4 = new TCanvas("c4","",600,600);
  c4->Draw();
  gPad->SetLogy();

  h3->Draw();
  TGraph *tg_mc_sig_EE = new TGraph(nbin, ptbin, mcsig_EE);
  tg_mc_sig_EE->SetMarkerSize(1.);
  tg_mc_sig_EE->Draw("p same");
  TGraphErrors *tgrs_sig_EE = new TGraphErrors(nbin, ptbin, sig_EE, pterr, sig_err_EE);
  tgrs_sig_EE->SetMarkerColor(4);
  tgrs_sig_EE->SetMarkerStyle(25);
  tgrs_sig_EE->SetLineColor(4);
  tgrs_sig_EE->SetLineWidth(2);
  tgrs_sig_EE->Draw("p e same");
  TGraphErrors *tgrs_sig_template_EE = new TGraphErrors(nbin, ptbin, sig_template_EE, pterr, sig_err_template_EE);
  tgrs_sig_template_EE->SetMarkerColor(6);
  tgrs_sig_template_EE->SetMarkerStyle(26);
  tgrs_sig_template_EE->SetLineColor(6);
  tgrs_sig_template_EE->SetLineWidth(2);
  tgrs_sig_template_EE->Draw("p e same");

  TGraphErrors *tgrs_sig_twobin_EE = new TGraphErrors(nbin, ptbin, nsig_twobin_EE, pterr, nsig_err_twobin_EE);
  tgrs_sig_twobin_EE->SetMarkerColor(kGreen);
  tgrs_sig_twobin_EE->SetMarkerStyle(3);
  tgrs_sig_twobin_EE->SetLineColor(kGreen);
  tgrs_sig_twobin_EE->SetLineWidth(2);
  tgrs_sig_twobin_EE->Draw("p e same");



  TLegend *tleg = new TLegend(0.193, 0.168, 0.642, 0.437);
  char text[50];  
  tleg->SetHeader("EE Photons");
  tleg->SetFillColor(0);
  tleg->SetShadowColor(0);
  tleg->SetBorderSize(0);
  tleg->AddEntry(tg_mc_sig_EE,"MC truth","p");
  tleg->AddEntry(tgrs_sig_EE,"unbinned fit","pl");
  tleg->AddEntry(tgrs_sig_template_EE,"template fit","pl");
  tleg->AddEntry(tgrs_sig_twobin_EE,"two bin","pl");
  tleg->Draw();


  c4->SaveAs("twobin_test/allfitter_yield_EE.pdf");


}
