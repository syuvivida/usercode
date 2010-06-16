
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <vector>
#include <string>
#include <TCanvas.h>
#include <TLegend.h>
#include <string>
#include "TCut.h"
#include "tightSelection.h"

using namespace std;


void ratioErr(float n1, float n1err, float n2, float n2err,
	      float& ratio, float& err)
{
  if(fabs(n1+n2)<1e-6){ratio=0;err=0;return;}
  ratio = (n1)/(n1+n2);
  
  err= pow(1/(n1+n2) - n1/(n1+n2)/(n1+n2),2)*n1err*n1err+
    pow(n1/(n1+n2)/(n1+n2),2)*n2err*n2err; 

  err = sqrt(err);

}

void histoErr(TH1F* hpass, TH1F* hfail, TH1F* heff)
{
  for(int i=1; i<=hpass->GetNbinsX(); i++){
   
    float n1 = hpass->GetBinContent(i);
    float n2 = hfail->GetBinContent(i);
    float n1err = hpass->GetBinError(i);
    float n2err = hfail->GetBinError(i);

    if(fabs(n1)<1e-6 && fabs(n2)<1e-6)continue;
    
    float ratio=0;
    float err=0;

    ratioErr(n1,n1err,n2,n2err,ratio,err);

    heff->SetBinContent(i,ratio);
    heff->SetBinError(i,err);

  }

  return;
}


TH1F* getPlot(TH1F *hTemplate, char *var, TCut cut, char *filename, char* treeName, float scale, int pthat_min, int pthat_max)
{
   TFile *inf = new TFile(filename);
   TTree *t = (TTree*)inf->FindObjectAny(treeName);
   TH1F * tmp = (TH1F*)hTemplate->Clone();
   tmp->SetName(Form("tmp%d",scale));
   char pthat_cut[1000];
   sprintf (pthat_cut,"(ptHat>%d && ptHat<%d)",pthat_min, pthat_max);
   TCut cut_tot;
   cut_tot = cut&&pthat_cut;
   cout << "current cut = " << cut_tot.GetTitle() << endl;
   t->Draw(Form("%s>>tmp%d",var,scale),cut_tot);
   tmp->Sumw2();
   cout << "scale = " << scale << endl;
   tmp->Scale(scale);
   cout<<tmp->GetEntries()<<endl;

   t->Delete();
   return tmp;
}

void MC_eff_Plot()
{
//    gROOT->LoadMacro("tdrstyle.C");
//    setTDRStyle();

   // Input files
   FILE *sTable = fopen("inputFile.txt","r");
   // photon Pt distribution
   char *var1 = "pt";
   // photon Eta distribution
   char *var2 = "eta";


   TCut kineCut[2];
   kineCut[0] = "pt > 15.0 && abs(eta) < 1.45";
   kineCut[1] = "pt > 15.0 && abs(eta) > 1.7 && abs(eta) < 2.5";


   TCut IDCut[2];
   
   IDCut[0] = rsCutEB;
   IDCut[1] = rsCutEE;

   TCut preSelection[2];
   TCut phoSelection_pass[2];
   TCut phoSelection_fail[2];

   
   for(int idec=0; idec<2; idec++){
     
     // photon selection for Denominator
     preSelection[idec] = kineCut[idec] + sigCut;
     // photon selections for Numerator
     phoSelection_pass[idec] = preSelection[idec] +  (IDCut[idec] + eventCut);
     phoSelection_fail[idec] = preSelection[idec] + !(IDCut[idec] + eventCut);

   }

   // Histograms for photon Pt distribution
   const float fBinsPt[]={15,20,30,50,80,120};
   const int nPtBin = sizeof(fBinsPt)/sizeof(fBinsPt[0])-1;
   
   TH1F *hTemplate_Pt = new TH1F("hTemplate_Pt", "", nPtBin, fBinsPt);
   hTemplate_Pt->SetXTitle("p_{T}(#gamma) [GeV/c]");
   TH1F *hs_Pt = (TH1F*)hTemplate_Pt->Clone();
   hs_Pt->SetName("hs_Pt");

   TH1F *denominator_Pt[2];
   TH1F *numerator_Pt_pass[2];
   TH1F *numerator_Pt_fail[2];
   TH1F *phoEff_Pt[2];

   for (int a=0; a<2; a++) {
     denominator_Pt[a] = (TH1F*)hTemplate_Pt->Clone(); 
     numerator_Pt_pass[a] = (TH1F*)hTemplate_Pt->Clone();
     numerator_Pt_fail[a] = (TH1F*)hTemplate_Pt->Clone();
     phoEff_Pt[a] = (TH1F*)hTemplate_Pt->Clone();

     denominator_Pt[a]->SetName(Form("preSelection_Pt%d", a));
     numerator_Pt_pass[a]->SetName(Form("phoSelection_pass_Pt%d", a));
     numerator_Pt_fail[a]->SetName(Form("phoSelection_fail_Pt%d", a));
     phoEff_Pt[a] ->SetName(Form("phoEff_Pt%d",a));
 
     denominator_Pt[a]->Reset();
     numerator_Pt_pass[a]->Reset();
     numerator_Pt_fail[a]->Reset();
     phoEff_Pt[a]->Reset();

     denominator_Pt[a]->Sumw2();
     numerator_Pt_pass[a]->Sumw2();
     numerator_Pt_fail[a]->Sumw2();

   }

   // Histograms for photon Eta distribution
   TH1F *hTemplate_Eta = new TH1F("hTemplate_Eta", "", 60, -3, 3);
   hTemplate_Eta->SetXTitle("#eta(#gamma)");
   TH1F *hs_Eta = (TH1F*) hTemplate_Eta->Clone();
   hs_Eta->SetName("hs_Eta");
   TH1F *denominator_Eta = (TH1F*) hTemplate_Eta->Clone();
   denominator_Eta->SetName("preSelection_Eta");
   TH1F *numerator_Eta_pass = (TH1F*) hTemplate_Eta->Clone();
   numerator_Eta_pass->SetName("phoSelection_pass_Eta");
   TH1F *numerator_Eta_fail = (TH1F*) hTemplate_Eta->Clone();
   numerator_Eta_fail->SetName("phoSelection_fail_Eta");
   TH1F *phoEff_Eta = (TH1F*)hTemplate_Eta->Clone();
   phoEff_Eta->SetName("phoEff_Eta");

   TH1F *phoEff_Eta_debug = (TH1F*)hTemplate_Eta->Clone();
   phoEff_Eta_debug->SetName("phoEff_Eta_debug");

   denominator_Eta->Reset();
   denominator_Eta->Sumw2();


   numerator_Eta_pass->Reset();
   numerator_Eta_pass->Sumw2();

   numerator_Eta_fail->Reset();
   numerator_Eta_fail->Sumw2();

   phoEff_Eta->Reset();
   
   phoEff_Eta_debug->Reset();
   phoEff_Eta_debug->Sumw2();

   int flag=1;
   int nfile=0;
   while (flag!=-1){
      // Input ROOT files
      char filename[200];
      flag=fscanf(sTable,"%s",filename);
      // Cross section
      char tmp[100];
      flag=fscanf(sTable,"%s",tmp);
      float cross=atof(tmp);
      // Number of event
      flag=fscanf(sTable,"%s",tmp);
      float nevt=atof(tmp);

      // ptHat min
      flag=fscanf(sTable,"%s",tmp);
      int ptmin=atof(tmp);

      // ptHat max
      flag=fscanf(sTable,"%s",tmp);
      int ptmax=atof(tmp);

      if (flag!=-1) {
         cout <<filename<<" "<<cross<<" "<<nevt<<endl;

         for (int a=0; a<2; a++) {
	   // Fill the histograms for photon Pt distribution
           hs_Pt->Reset();
           hs_Pt = getPlot(hTemplate_Pt, var1, preSelection[a], filename, "Analysis", cross/nevt, ptmin, ptmax);
           denominator_Pt[a]->Add(hs_Pt);

           hs_Pt->Reset();
           hs_Pt = getPlot(hTemplate_Pt, var1, phoSelection_pass[a], filename, "Analysis", cross/nevt, ptmin, ptmax);
           numerator_Pt_pass[a]->Add(hs_Pt);

           hs_Pt->Reset();
           hs_Pt = getPlot(hTemplate_Pt, var1, phoSelection_fail[a], filename, "Analysis", cross/nevt, ptmin, ptmax);
           numerator_Pt_fail[a]->Add(hs_Pt);


	   // Fill the histograms for photon Eta distribution

	   hs_Eta->Reset();
	   hs_Eta = getPlot(hTemplate_Eta, var2, preSelection[a], filename, "Analysis", cross/nevt, ptmin, ptmax);
	   denominator_Eta->Add(hs_Eta);

	   hs_Eta->Reset();
	   hs_Eta = getPlot(hTemplate_Eta, var2, phoSelection_pass[a], filename, "Analysis", cross/nevt, ptmin, ptmax);
	   numerator_Eta_pass->Add(hs_Eta);
	   
	   hs_Eta->Reset();
	   hs_Eta = getPlot(hTemplate_Eta, var2, phoSelection_fail[a], filename, "Analysis", cross/nevt, ptmin, ptmax);
	   numerator_Eta_fail->Add(hs_Eta);
         } // loop over detectors

      }
   } // finish reading files


   for (int a=0; a<2; a++) {
     
     histoErr(numerator_Pt_pass[a],numerator_Pt_fail[a],phoEff_Pt[a]);
     phoEff_Pt[a]->GetYaxis()->SetDecimals();
     phoEff_Pt[a]->SetTitleOffset(1.4,"Y");

   }
   phoEff_Pt[0]->SetTitle("|#eta(#gamma)|<1.45");
   phoEff_Pt[1]->SetTitle("1.7 < |#eta(#gamma)|<2.5");

   histoErr(numerator_Eta_pass,numerator_Eta_fail,phoEff_Eta);
   phoEff_Eta->GetYaxis()->SetDecimals();
   phoEff_Eta->SetTitleOffset(1.4,"Y");
   phoEff_Eta->SetTitle("Efficiency as a function of eta");

   phoEff_Eta_debug->Divide(numerator_Eta_pass,denominator_Eta,1.0,1.0,"B");


   for(int i=1;i<= 60;i++)
     {
       cout << "Numerator = " << numerator_Eta_pass->GetBinContent(i) << endl;
       cout << "Denominator = " << denominator_Eta->GetBinContent(i) << endl;
       cout << "efficiency = " << phoEff_Eta->GetBinContent(i) << " \t" << phoEff_Eta_debug->GetBinContent(i) << endl;
       cout << "efficiency error = " << phoEff_Eta->GetBinError(i) << " \t" << phoEff_Eta_debug->GetBinError(i) << endl;

     }


// barrel and endcap on separate canvases
   TCanvas *c1 = new TCanvas("c1", "", 500, 500);
   phoEff_Pt[0]->GetYaxis()->SetTitle("ID Efficiency");
   phoEff_Pt[0]->SetMarkerColor(2);
   phoEff_Pt[0]->SetMarkerStyle(21);
   phoEff_Pt[0]->SetMarkerSize(0.7);
   phoEff_Pt[0]->Draw("e");
   
   TCanvas *c2 = new TCanvas("c2", "", 500, 500);
   phoEff_Pt[1]->GetYaxis()->SetTitle("ID Efficiency");
   phoEff_Pt[1]->SetMarkerColor(2);
   phoEff_Pt[1]->SetMarkerStyle(21);
   phoEff_Pt[1]->SetMarkerSize(0.7);
   phoEff_Pt[1]->Draw("e");

// eta canvas 
   TCanvas *c3 = new TCanvas("c3", "", 500, 500);
   c3->cd();
   phoEff_Eta->GetYaxis()->SetTitle("ID Efficiency");
   phoEff_Eta->SetMarkerColor(2);
   phoEff_Eta->SetMarkerStyle(21);
   phoEff_Eta->SetMarkerSize(0.7);
   phoEff_Eta->Draw("e");


   TFile *file = new TFile("Efficiency.root", "UPDATE");
   for (int a=0; a<2; a++) {
     phoEff_Pt[a]->Write();
     denominator_Pt[a]->Write();
     numerator_Pt_pass[a]->Write();
     numerator_Pt_fail[a]->Write();
     
   }

   denominator_Eta->Write();
   numerator_Eta_pass->Write();
   numerator_Eta_fail->Write();
   phoEff_Eta->Write();
   phoEff_Eta_debug->Write();
   file->Close();
}

