
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <string>
#include "purity_twobin.C"

using namespace std;

void call_twobin(std::string filename)
{

  TFile *inf = new TFile(filename.data());
  std::string remword="_histo.root";
  size_t pos = filename.find(remword);
  std::string forOutput = filename;  
  if(pos!= std::string::npos)
    forOutput.swap(forOutput.erase(pos,remword.length()));   

  std::string histo = forOutput + "MixData";
  cout << "Finding " << histo << endl;
  TH1F *h_data = (TH1F*)inf->FindObjectAny(histo.data());

  histo = forOutput + "MixSig";
  cout << "Finding " << histo << endl;
  TH1F *h_sig  = (TH1F*)inf->FindObjectAny(histo.data());

  histo = forOutput + "MixBkg";
  cout << "Finding " << histo << endl;
  TH1F *h_bkg  = (TH1F*)inf->FindObjectAny(histo.data());

//   TH1F* h_data = (TH1F*)h_sig->Clone();
//   h_data->Reset();
//   h_data->Sumw2();
//   h_data->Add(h_sig,h_bkg,1.0,1.0);
 
  purity_twobin(h_data,h_sig,h_bkg,forOutput);


}
