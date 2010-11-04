
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <string>
#include <TCanvas.h>
#include "fitLifeTime.C"
#include <string>
#include <iostream>

using namespace std;

const int REBINNINGS=3;

const double fBinsPt[]={21.,23.,26.,30.,35.,40.,45.,50.,60.,85.,120.,300,10000};
const int nPtBin = sizeof(fBinsPt)/sizeof(fBinsPt[0])-1;
const int nEtaBin = 2;

void call_fitLifeTime(){

  TFile* finFile = new TFile("template_comb3Iso_template.root");
  
  TH1D* hTemplate = (TH1D*)finFile->FindObjectAny("h_EB_comb3Iso_EGdata_pt21");
  hTemplate->Reset();

  char* dec[2] = {"EB","EE"};
  char tmp[1000];

  TH1D* hTemplate_S[nEtaBin][nPtBin];

  for(int ieta = 0; ieta < nEtaBin; ieta++){
     for(int ipt=0; ipt < nPtBin; ipt++){

       sprintf(tmp,"h_%s_comb3Iso_sig_pt%d",dec[ieta],(int)fBinsPt[ipt]);
       cout << "looking for histogram " << tmp << " in file " << 
	 finFile->GetName() << endl;
       hTemplate_S[ieta][ipt] = (TH1D*)finFile->FindObjectAny(tmp);
       hTemplate_S[ieta][ipt]->Rebin(REBINNINGS);

     }
  }


  for(int ieta = 0; ieta < nEtaBin; ieta++){
     for(int ipt=0; ipt < nPtBin; ipt++){

       fitLifeTime(hTemplate_S[ieta][ipt]);


     }
  }





}
