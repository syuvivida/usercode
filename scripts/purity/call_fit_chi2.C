
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include "fit_chi2.C"

void call_fit_chi2(std::string filename)
{

  TFile *inf = new TFile(filename.data());
//   TH1F *h_data = (TH1F*)inf->FindObjectAny("testfitterMixData");
  TH1F *h_sig  = (TH1F*)inf->FindObjectAny("testfitterMixSig");
  TH1F *h_bkg  = (TH1F*)inf->FindObjectAny("testfitterMixBkg");

  TH1F* h_data = h_sig->Clone();
  h_data->Reset();
  h_data->Sumw2();
  h_data->Add(h_sig,h_bkg,1.0,1.0);
 
  fit_chi2(h_data,h_sig,h_bkg);


}
