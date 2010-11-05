#include <TMath.h>
#include <TF1.h>
#include <TH1.h>
#include <iostream>
#include <fstream>
#include <TMinuit.h>

using namespace std;

Double_t exp_conv (Double_t *v, Double_t *par)
{
  Double_t ctau = par[1];
  Double_t sigma = par[3]; //using narrow prompt width
  Double_t x = v[0]-par[2];

  Double_t arg1 = TMath::Exp( 0.5*sigma*sigma/ctau/ctau - x/ctau );
  Double_t arg2 = 1.0 - TMath::Freq( (sigma/ctau - x/sigma) );
  Double_t func = par[0]/ctau * arg1 * arg2;

  if (func<=0) func=1e-10;
  return func;
}


void fitLifeTime(TH1D* hfit)
{
  const float fit_lo_edge =  -1.0;
  const float fit_hi_edge =   5.0;
  const int BINNINGS = 3;

//   TH1D* hfit = (TH1D*)gROOT->FindObject(histoName);
//   hfit->Rebin(BINNINGS);

  TF1* f1 = new TF1("f1",exp_conv, -1., 20.,11);
  f1->SetParName(0, "normalization");
  f1->SetParName(1, "exp ctau");
  f1->SetParName(2, "offset");
  f1->SetParName(3, "exp sigma");
  
  f1->SetParameters(hfit->Integral(), 1., 0.6, 0.3);
  f1->SetNpx(2500);
  f1->SetLineColor(2);

  hfit->Fit("f1","","",fit_lo_edge,fit_hi_edge);
  gMinuit->mnmatu(1);

  hfit->Fit("f1","","",fit_lo_edge,fit_hi_edge);

  hfit->Fit("f1","","",fit_lo_edge,fit_hi_edge);
  int fit_status=hfit->Fit("f1","","",fit_lo_edge,fit_hi_edge);
  if(fit_status!=0)
  {	
     f1->SetParameters(1000., 1., 0.6, 0.3);
     hfit->Fit("f1","","",fit_lo_edge,fit_hi_edge);
   }


  ofstream fout;
  fout.open("sigGausMean.dat",ios::app | ios::out);
  fout << hfit->GetName() << "\t" << f1->GetParameter(2) << " " << 
    f1->GetParError(2) << endl;
  fout.close();
  
  cout << " ==================== Fit result ==================" << endl;
  for(int i=0; i< f1-> GetNpar(); i++)
    cout << " Parameter " << i << " = " << f1->GetParameter(i) << endl;
  
  cout << " ==================================================" << endl;




}

