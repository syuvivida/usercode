#include <iostream>
#include <TFile.h>
#include <TH2.h>
#include <TProfile.h>
#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"

void plotCorr(char* infname)
{
   TFile *inf = new TFile(infname);
   TH2F* data = (TH2F*)inf->FindObjectAny("h_EB_comb3Iso_SIEIE_EGdata_SIG");
   TH2F* mc = (TH2F*)inf->FindObjectAny("h_EB_comb3Iso_SIEIE_bkg_sum_BKG");

   setTDRStyle();
   const float ymin = -1.0;
   const float ymax = 11.0;

   TProfile* pf_data = (TProfile*)data->ProfileY("pf_data",1,120);
   pf_data->SetTitle("");
   pf_data->SetMinimum(ymin);
   pf_data->SetMaximum(ymax);
   pf_data->GetXaxis()->SetRangeUser(0.01,0.02);
   pf_data->SetLineColor(1);
   pf_data->SetMarkerColor(1);
   pf_data->SetMarkerSize(1);
   pf_data->SetMarkerStyle(8);
   pf_data->GetXaxis()->SetNdivisions(5);
   pf_data->GetXaxis()->SetDecimals();

   pf_data->GetYaxis()->SetTitle("<Iso>");
   pf_data->GetXaxis()->SetTitle("#sigma_{i#eta i#eta}");

   TProfile* pf_mc = (TProfile*)mc->ProfileY("pf_mc",1,120);
   pf_mc->SetTitle("");
   pf_mc->SetMinimum(ymin);
   pf_mc->SetMaximum(ymax);
   pf_mc->GetXaxis()->SetRangeUser(0.01,0.02);
   pf_mc->SetLineColor(2);
   pf_mc->SetMarkerColor(2);
   pf_mc->SetMarkerSize(1);
   pf_mc->SetMarkerStyle(23);

   pf_mc->GetYaxis()->SetTitle("<Iso>");
   pf_mc->GetXaxis()->SetTitle("#sigma_{i#eta i#eta}");
   pf_mc->GetXaxis()->SetNdivisions(5);
   pf_mc->GetXaxis()->SetDecimals();

   TCanvas* c1 = new TCanvas("c1","",500,500);
   pf_data->Draw("");
   pf_mc->Draw("same");
   TLegend* leg = new TLegend(0.40,0.27,0.60,0.42);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->SetTextSize(0.04);
   leg->SetBorderSize(0);
   leg->AddEntry(pf_data,"Data sideband");
   leg->AddEntry(pf_mc,"MC sideband");   
   leg->Draw("same");


   c1->Print("Iso_SIEIE_corr.eps");
   c1->Print("Iso_SIEIE_corr.gif");
   c1->Print("Iso_SIEIE_corr.C");

   
}
