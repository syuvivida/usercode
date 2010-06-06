#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <string>
#include <TCut.h>
#include <TCanvas.h>
#include <TLegend.h>
#include "tightSelection.h"
#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"

using namespace std;

void compareNoUE(const char* var, int ieta,bool logy=false)
{  
  setTDRStyle();
  gStyle->SetOptStat(0);
  TH1F* hTemplate = new TH1F("hTemplate","",30,-1.0,5.0);
  TH1F* hnoUE = (TH1F*)hTemplate->Clone();
  TH1F* hD6T = (TH1F*)hTemplate->Clone();
  TH1F* hATLAS = (TH1F*)hTemplate->Clone();

  TTree* tnoUE;
  TTree* tD6T;
  TTree* tATLAS;

  TCut allCut = ieta==0? sigCut + rsCutEB: sigCut + rsCutEE;
  cout << "cutting on " << allCut.GetTitle() << endl;

  TFile *inf0 = new TFile("MPA-noue_pthat15.root");
  tnoUE = (TTree*)inf0->FindObjectAny("Analysis"); 
  TH1F *htmp = (TH1F*)hTemplate->Clone();
  htmp->SetName("htmp");
  htmp->Reset();
  tnoUE->Draw(Form("%s>>htmp",var),allCut);
  hnoUE = (TH1F*)htmp->Clone();
  cout << hnoUE->Integral() << endl;

  TFile *inf1 = new TFile("MPA-d6t_pthat15.root");
  tD6T = (TTree*)inf1->FindObjectAny("Analysis");
  TH1F *htmp1 = (TH1F*)hTemplate->Clone();
  htmp1->SetName("htmp1");
  htmp1->Reset();
  tD6T->Draw(Form("%s>>htmp1",var),allCut);
  hD6T = (TH1F*)htmp1->Clone();
  cout << hD6T->Integral() << endl;


  TFile *inf2 = new TFile("MPA-atlas_pthat15.root");
  tATLAS   = (TTree*)inf2->FindObjectAny("Analysis");
  TH1F *htmp2 = (TH1F*)hTemplate->Clone();
  htmp2->SetName("htmp2");
  htmp2->Reset();
  tATLAS->Draw(Form("%s>>htmp2",var),allCut);
  hATLAS = (TH1F*)htmp2->Clone();
  cout << hATLAS->Integral() << endl;
  

  hnoUE->SetName("hnoUE");
  hnoUE->SetLineColor(2);
  hnoUE->SetMarkerColor(2);
  hnoUE->SetMarkerSize(1);
  hnoUE->SetMarkerStyle(24);
  hnoUE->SetXTitle(var);
  hD6T->SetName("hD6T");
  hD6T->SetLineColor(4);
  hD6T->SetMarkerColor(4);
  hD6T->SetMarkerSize(1);
  hD6T->SetMarkerStyle(21);
  hD6T->SetXTitle(var);
  hATLAS->SetName("hATLAS");
  hATLAS->SetLineColor(3);
  hATLAS->SetMarkerColor(3);
  hATLAS->SetMarkerSize(1);
  hATLAS->SetMarkerStyle(26);
  hATLAS->SetXTitle(var);

  hnoUE->Sumw2();
  hnoUE->Scale(1.0/(float)hnoUE->Integral());
  hD6T->Sumw2();
  hD6T->Scale(1.0/(float)hD6T->Integral());
  hATLAS->Sumw2();
  hATLAS->Scale(1.0/(float)hATLAS->Integral());

  TCanvas* c1 = new TCanvas("c1","",500,500);
  if(logy)c1->SetLogy(1);
  else c1->SetLogy(0);
  hnoUE->Draw("he");
  hD6T->Draw("hesame");
  hATLAS->Draw("hesame");

  std::string decName;
  decName = ieta==0? "EB": "EE";

  TLegend* leg = new TLegend(0.50,0.62,0.70,0.87);
  leg->SetHeader(Form("PhotonJet_Pt15 (%s)",decName.data()));
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->AddEntry(hnoUE,"No UE");
  leg->AddEntry(hD6T,"D6T");
  leg->AddEntry(hATLAS,"ATLAS");
  leg->Draw("same");

  std::string logName;
  logName = logy? "log":"";
  c1->Print(Form("%s_UECompPhoJetPt15_pt20_%s%s.eps",var,decName.data(),
		 logName.data()));
  c1->Print(Form("%s_UECompPhoJetPt15_pt20_%s%s.gif",var,decName.data(),
		 logName.data()));

}
