#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TProfile.h>
#include "tightSelection.h"
#include <TCut.h>

TProfile* getPlot(TProfile *h1, char *var, TCut cut, char *filename, char* treeName)
{
  TFile *inf = new TFile(filename);
  TTree *t = (TTree*)inf->FindObjectAny(treeName);
  TProfile * tmp = (TProfile*)h1->Clone();
  tmp->SetName("tmp");
  t->Draw(Form("%s>>tmp",var),cut);  
  return tmp;
}


TH2F* getPlot(TH2F *h1, char *var, TCut cut, char *filename, char* treeName)
{
  TFile *inf = new TFile(filename);
  TTree *t = (TTree*)inf->FindObjectAny(treeName);
  TH2F * tmp = (TH2F*)h1->Clone();
  tmp->SetName("tmp");
  t->Draw(Form("%s>>tmp",var),cut);  
  return tmp;
}

void plotPhotonType(char* var="(ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04):genCalIsoDR04")
{
  
  TCut allCut = simpleCut + hardScatterCut;
   
  TProfile *hTemplate;
  TProfile *htmp= new TProfile("htmp","",100,0,20);
  htmp->SetXTitle("genCalIso [GeV]");
  htmp->SetYTitle("recCalIso [GeV]");

  hTemplate= getPlot(htmp,var,allCut,"MPA_PhotonJetPt15_31X.root","Analysis");
  hTemplate->SetName("hTemplate");
  hTemplate->Draw();
  gStyle->SetOptFit(11111);
  hTemplate->Fit("pol1","","",15,20);
  TF1* f1 = hTemplate->GetFunction("pol1");

  float p0 = f1->GetParameter(0);
  float p1 = f1->GetParameter(1);

  char tempName[1000];
  sprintf(tempName,
 	  "((ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04)-genCalIsoDR04*%f):((ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04)-genCalIsoDR04*(%f))",p1,-1.0/p1);
  cout << tempName << endl;

  TProfile *hTemplate_decompos;
  const int nbin=50;
  const float min = 0.0;
  const float max = 10.0;
  TProfile *htmp2= new TProfile("htmp2","",nbin,min,max);
  hTemplate_decompos= getPlot(htmp2,tempName,allCut,"MPA_PhotonJetPt15_31X.root","Analysis");
  hTemplate_decompos->SetName("hTemplate_decompos");
  hTemplate_decompos->Draw();  
  hTemplate_decompos->SetYTitle(Form("recCalIso-genCalIso*%.2f",p1));
  hTemplate_decompos->SetXTitle(Form("recCalIso+genCalIso*%.2f",1.0/p1));
  gStyle->SetOptFit(11111);
//   hTemplate_decompos->Fit("pol1");
  
  
  TCanvas* c1 = new TCanvas("c1","",500,1000);
  c1->Divide(1,2);
  c1->cd(1);
  hTemplate->Draw("");
  c1->cd(2);
  hTemplate_decompos->SetErrorOption("S");
  hTemplate_decompos->Draw("");
  c1->Print("bestGenCalIsoDR04.eps");
  c1->Print("bestGenCalIsoDR04.gif");
 
  TCanvas* c2 = new TCanvas("c2","",500,1000);
  c2->Divide(1,2);
  c2->cd(1);
  TH1F* hMean = new TH1F("hMean","",nbin,min,max);
  hMean->SetXTitle(Form("recCalIso+genCalIso*%.2f",1.0/p1));
  hMean->SetTitleSize(0.06,"Y");
  hMean->SetTitleOffset(1.2,"Y");
  hMean->SetYTitle(Form("Mean of recCalIso-genCalIso*%.2f",p1));
  for(int i=1; i <= nbin; i++)
    hMean->SetBinContent(i,hTemplate_decompos->GetBinContent(i));
  hMean->Draw();
  c2->cd(2);
  TH1F* hRMS = new TH1F("hRMS","",nbin,min,max);
  hRMS->SetXTitle(Form("recCalIso+genCalIso*%.2f",1.0/p1));
  hRMS->SetTitleSize(0.06,"Y");
  hRMS->SetTitleOffset(1.2,"Y");
  hRMS->SetYTitle(Form("RMS of recCalIso-genCalIso*%.2f",p1));
  for(int i=1; i <= nbin; i++)
    hRMS->SetBinContent(i,hTemplate_decompos->GetBinError(i));
  hRMS->Draw();
  c2->Print("bestGenCalIsoDR04_sup.eps");
  c2->Print("bestGenCalIsoDR04_sup.gif");


  int bestDeComposXBin = 11;
  
  float bestDeComposX = hMean->GetBinCenter(bestDeComposXBin);
  float bestDeComposY = hMean->GetBinContent(bestDeComposXBin);

  cout << "bestDeComposX = " << bestDeComposX << endl;
  cout << "bestDeComposY = " << bestDeComposY << endl;
  float bestGenIso = (bestDeComposX - bestDeComposY)/((1.0)/p1 + p1);
  float bestRecIso =  bestDeComposX - (1.0/p1) * bestGenIso;

  cout << "bestGenIso = " << bestGenIso << endl;
  cout << "bestRecIso = " << bestRecIso << endl;


}


