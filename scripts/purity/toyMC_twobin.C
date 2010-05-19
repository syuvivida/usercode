#include <TRandom.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <string>
#include "purity_twobin.C"
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

using namespace std;

void toyMC_twobin()
{

  TFile *inf = new TFile("template_geniso5.root");
  double fBinsPt[]={15,20,30,50};
  const int nPtBin = sizeof(fBinsPt)/sizeof(fBinsPt[0])-1;
  const int nEtaBin = 2; // 0~1.45, 1.55~2.5


  TH1D* htmp_data;
  TH1D* htmp_sig;
  TH1D* htmp_bkg;

  TH1D* hrdm_data;
  TH1D* hrdm_sig;
  TH1D* hrdm_bkg;

  const int nData = 300;
  const int nTemplate = 1000;
  const double purityInput = 0.2;
  const int nCount = 1000;

  std::string histoName;
  std::string decName = "EB";
  std::string tmpName;
  char tmp[1000];

  TH1D* hPull    = new TH1D("hPull","",100,-5,5);
  TH1D* hCentral = new TH1D("hCentral","",100,0,1);
  TH1D* hError   = new TH1D("hError","",100,0,1);

  for(int ic=0; ic < nCount; ic++){
  for(int ipt=0; ipt < 1; ipt++){
 
    sprintf(tmp,"ecalhcalIso_EGdata_pt%d",(int)fBinsPt[ipt]);
    tmpName = tmp; 
    histoName = "h_" + decName + "_" + tmpName; 
    htmp_data = (TH1D*)inf->FindObjectAny(histoName.data());
    hrdm_data = (TH1D*)htmp_data->Clone();
    hrdm_data->Reset();
    
    sprintf(tmp,"ecalhcalIso_sig_pt%d",(int)fBinsPt[ipt]);
    tmpName = tmp;
    histoName = "h_" + decName + "_" + tmpName; 
    htmp_sig = (TH1D*)inf->FindObjectAny(histoName.data());
    hrdm_sig = (TH1D*)htmp_sig ->Clone();
    hrdm_sig ->Reset();
    

    sprintf(tmp,"ecalhcalIso_bkg_pt%d",(int)fBinsPt[ipt]);
    tmpName = tmp;
    histoName = "h_" + decName + "_" + tmpName; 
    htmp_bkg = (TH1D*)inf->FindObjectAny(histoName.data());
    hrdm_bkg = (TH1D*)htmp_bkg ->Clone();
    hrdm_bkg ->Reset();
     

    int ndata_rdm_sig = gRandom->Poisson(nData*purityInput);      
    int ndata_rdm_bkg = gRandom->Poisson(nData*(1-purityInput));
  
    hrdm_data->FillRandom(htmp_sig, ndata_rdm_sig);
    hrdm_data->FillRandom(htmp_bkg, ndata_rdm_bkg);

    cout << "hrdm_data->GetEntries() = " <<hrdm_data->GetEntries() << endl;

    int ntemplate_rdm = gRandom->Poisson(nTemplate);

    hrdm_sig ->FillRandom(htmp_sig, ntemplate_rdm);
    hrdm_bkg ->FillRandom(htmp_bkg, ntemplate_rdm);


    cout << "hrdm_sig->GetEntries() = " <<hrdm_sig->GetEntries() << endl;
    cout << "hrdm_bkg->GetEntries() = " <<hrdm_bkg->GetEntries() << endl;

    double sigfrac;
    double sigfrac_err;
    sprintf(tmp,"ecalhcalIso_EGdata_pt%d",(int)fBinsPt[ipt]);
    tmpName = tmp;
    purity_twobin(hrdm_data,hrdm_sig,hrdm_bkg,decName + "_" + tmpName,sigfrac,sigfrac_err);
    if(sigfrac<1e-6)continue;
    hPull->Fill((sigfrac-purityInput)/sigfrac_err);
    hCentral->Fill(sigfrac);
    hError->Fill(sigfrac_err);

  }

  }
  TCanvas* c1 = new TCanvas("c1","",1500,500);
  c1->Divide(3,1);
  c1->cd(1);
  gStyle->SetOptStat(0);
  hPull->Draw();
  hPull->SetXTitle("(central-input)/error");
  hPull->Fit("gaus");
  c1->cd(2);
  gStyle->SetOptStat(1);
  hCentral->Draw();
  hCentral->SetXTitle("central value");
  c1->cd(3);
  gStyle->SetOptStat(1);
  hError->Draw();
  hError->SetXTitle("error");

}
