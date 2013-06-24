#include <string>
#include <iostream>
#include <fstream>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TStyle.h>
#include <TAxis.h>

using namespace std;

const int nRadius=7;

void readPhiR(string inFileName){
  
  double r[nRadius];
  for(int i=0;i<nRadius; i++)r[i]=0.1*(i+1);
  double er[nRadius]={0};
  double phir[nRadius]={0};
  double ephir[nRadius]={0};
  
  ifstream fin;
  fin.open(inFileName.data());
  for(int i=0; i < nRadius; i++)
    {
      fin >> phir[i];
    }

  fin.close();

  TGraphErrors* h_out = new TGraphErrors(
					 nRadius, r, phir,er,ephir);
  h_out->SetName("h_out");
  h_out->SetMarkerStyle(8);
  h_out->SetMarkerColor(4);
  h_out->SetMarkerSize(1);
  h_out->GetXaxis()->SetTitle("jet raidus: r");
  h_out->GetYaxis()->SetTitle("#Phi(r)");
  // h_out->GetYaxis()->SetTitleOffset(1.5);
  h_out->SetTitle("Predictions of CMS 7 TeV Inclusive Jet");
  h_out->SetLineWidth(2);


  string remword  =".out";
  string outFileName = inFileName;
  size_t pos  = inFileName.find(remword);

  if(pos!= std::string::npos)
    outFileName.swap(outFileName.erase(pos,remword.length()));
  
  TFile* outFile = new TFile(Form("histo_%s.root",outFileName.data()),
			     "recreate");       
  h_out->Write();
  outFile->Close();


}
