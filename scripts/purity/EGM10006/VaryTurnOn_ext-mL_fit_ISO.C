#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TF1.h"
#include <fstream>
#include <iostream>
#include <TMinuit.h>
//#include "vector.h"
#include <TMath.h>
#include "TVirtualFitter.h"
#include "TFile.h"
#include "/afs/cern.ch/user/s/syu/scratch0/VaryShape/chi2Nbins.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TLine.h"
#include "TRandom2.h"

using namespace std;

#define NPAR 11

TF1* fchi2;
//TH1F* h1;
TF1* bsFin;

// Remember to change the values of xh and xl in all_tree.C too!!!
//Double_t xh = 6.5;
//Double_t xl = 4.5;
Double_t xl = -1.;
Double_t xh = 1.;
Double_t bin_size = 0.05;

const double _two_pi = 2.0 * TMath::Pi();
Double_t fit_lo_edge = -1.;
Double_t fit_hi_edge = 1.;


vector<Double_t> dataColl;
vector<Double_t> sigColl;
vector<Double_t> bkgColl;

vector<Double_t> totalColl;
vector<Double_t> ctauColl;

vector<Double_t> info;
vector<Double_t> info_err;

vector<Double_t> Para;
vector<Double_t> Para_err;

Double_t SigPDFnorm = 0.;
Double_t BkgPDFnorm = 0.;

Int_t para_index = 0;
Double_t para_sigma = 1.;

Int_t Opt_SavePDF = 1;

//par[0] = fs Jpsi signal fraction;
//-----mass part-----------------------------------------------------
//par[1] = g norm; par[2] g1 mean; par[3] g1 width; 
//par[4] g2 ratio;  par[5] g2 mean; par[6] g2 width;
//par[7] bg norm; par[8] bg slope;


Double_t g(Double_t *v, Double_t *par)
{
  Double_t arg = 0;
  arg = (v[0] - par[2]) / par[3];
  
  //Double_t area = par[1]/(par[3]*sqrt(2*3.1415926));
  Double_t norm_area = 1./(par[3]*sqrt(2*3.1415926));
  Double_t gaus = norm_area*TMath::Exp(-0.5*arg*arg);

  if (gaus<=0) gaus=1e-10;
  return gaus;
}

Double_t bifg(Double_t *v, Double_t *par)
{
  Double_t arg = 0;
  if ( v[0] > par[1] ) arg = (v[0] - par[1]) / par[2];
  else arg = (v[0] - par[1]) / par[3];
  
  //Double_t area = par[1]/(par[3]*sqrt(2*3.1415926));
  Double_t norm_area = 1./((par[2]+par[3])*0.5*sqrt(2*3.1415926));
  Double_t gaus = norm_area*TMath::Exp(-0.5*arg*arg);

  if (gaus<=0) gaus=1e-10;
  return gaus;
}

Double_t exp_conv (Double_t *v, Double_t *par)
{
  Double_t ctau = par[1];
  //Double_t sigma = par[14];
  Double_t sigma = par[3]; //using narrow prompt width
  Double_t x = v[0]-par[2];
  
  Double_t arg1 = TMath::Exp( 0.5*sigma*sigma/ctau/ctau - x/ctau );
  Double_t arg2 = 1.0 - TMath::Freq( (sigma/ctau - x/sigma) );
  //Double_t func = 1.0/ctau * arg1 * arg2;
  Double_t func = par[0]/ctau * arg1 * arg2;
  //Double_t func = bifg(v,par);
  if (func<=0) func=1e-10;
  return func;
}

Double_t exp_conv_norm(Double_t *v, Double_t *par)
{
  Double_t func = exp_conv(v,par) / SigPDFnorm;
  if (func<=0) func=1e-10;
  return func;
}

Double_t expinv_power(Double_t *v, Double_t *par){

  Double_t x = v[0]-par[6];
  Double_t func=0.;
  if (x>0.) {
    Double_t fitval = 1.- TMath::Exp(par[5]*x);
    //func = par[4] * (1-par[9]*x) *TMath::Power(1-par[7]*x,par[8]) * fitval;
    func = par[4] * TMath::Power(1-par[7]*x,par[8]) * fitval;
    //func = par[4] * fitval;
  }
  return func;
}

Double_t expinv_power_norm (Double_t *v, Double_t *par)
{
  Double_t func = expinv_power(v,par)/BkgPDFnorm;
  if (func<=0) func=1e-10;
  return func;
}

Double_t sum_norm(Double_t *v, Double_t *par)
{
  Double_t func1 = exp_conv(v,par) / SigPDFnorm;
  Double_t func2 = expinv_power(v,par)/BkgPDFnorm;
  Double_t func = func1+func2;
  if (func<=0) func=1e-10;
  return func;
}

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

  Double_t Lsum=0.;
  Double_t Nevt=dataColl.size();
  Double_t fs = par[0];
  Double_t fb = par[1];

  double tmppar[12];
  for(int i=0; i<12; i++){
    tmppar[i] = par[i+2];
  }
  TF1 *funcS = new TF1("funcS", exp_conv, -1., 20., 12);
  funcS->SetParameters(tmppar);
  Double_t Signorm = funcS->Integral(-1., 20.);
  TF1 *funcB = new TF1("funcB", expinv_power, -1., 20., 12);
  funcB->SetParameters(tmppar);
  Double_t Bkgnorm = funcB->Integral(-1., 20.);

  SigPDFnorm = Signorm;
  BkgPDFnorm = Bkgnorm;

  for ( int i=0; i<dataColl.size(); i++ ) {
    //PDF for signal and background
    Double_t x = dataColl[i];
    Double_t Ls = exp_conv(&x,tmppar)/Signorm;
    Double_t Lb = expinv_power(&x,tmppar)/Bkgnorm;
    //Get Log Likelihood		  
    Lsum += TMath::Log( (fs*Ls + fb*Lb) / (fs+fb) );
  }
  //f=2*( -1*Lsum + (fs+fb) - Nevt*TMath::Log(fs+fb) );
  f=2*( -1*Lsum + (fs+fb) - Nevt*TMath::Log(fs+fb) );
}


//___________________________________________________________________________
Double_t* Ifit(int shift, Double_t& dataResult, Double_t& dataErr, std::string dataFile, 
	       TH1D* hsig, TH1D* hbkg, TH1D* hEGdata, Double_t* FitPar,
	       int ptbin=30, char EBEE[10]="EB", int fit_data=2)
{
  
  printf(" *** calling Ifit for %s , ptbin %d *** \n\n", EBEE,ptbin);

  cout << "The number of bins are: " << endl;
  cout << "hdata nbins = " << hEGdata->GetNbinsX() << endl;
  cout << "hsig nbins = " << hsig->GetNbinsX() << endl;
  cout << "hbkg nbins = " << hbkg->GetNbinsX() << endl;
  

  TCanvas *c1 = new TCanvas("HF1", "Histos1", 0, 0, 600, 600);
  gStyle->SetOptFit(0);

  if(fit_data != 3) dataColl.clear();
  sigColl.clear();
  bkgColl.clear();

  totalColl.clear();
  ctauColl.clear();
  Para.clear();
  Para_err.clear();

  info.clear();
  info_err.clear();

  float ptmax=0.;  
  if(ptbin== 21) ptmax= 23;
  if(ptbin== 23) ptmax= 26;
  if(ptbin== 26) ptmax= 30;
  if(ptbin== 30) ptmax= 35;
  if(ptbin== 35) ptmax= 40;
  if(ptbin== 40) ptmax= 45;
  if(ptbin== 45) ptmax= 50;
  if(ptbin== 50) ptmax= 60;
  if(ptbin== 60) ptmax= 85;
  if(ptbin== 85) ptmax= 120;
  if(ptbin== 120) ptmax= 300;
  if(ptbin== 300) ptmax= 500;



  Double_t* fitted = new Double_t[6];
  fitted[0] = 0.;    fitted[1] = 0.;    fitted[2] = 0.;    fitted[3] = 0.;
  fitted[4] = 0.;    fitted[5] = 0.;

  char hname[30];


  hsig->SetLineColor(1);
  hbkg->SetLineColor(1);
  hsig->SetNdivisions(505,"XY");
  hbkg->SetNdivisions(505,"XY");
  hsig->SetTitle("");
  hbkg->SetTitle("");
  hsig->SetXTitle("combined ISO (GeV)");
  hbkg->SetXTitle("combined ISO (GeV)");

  TH1F *hsum = (TH1F*)hsig->Clone();  
  hsum->Add(hbkg,1);
  float ntemplate = 1.;
  if (hsum->Integral()>1.) ntemplate = hsum->Integral();
  float sigfrac = hsig->Integral()/ntemplate*0.8;

  TH1F *hsum_norm = (TH1F*)hsum->Clone();  
  hsum_norm->Scale(1./hsum->Integral());

  TH1F *hdata = new TH1F();
  int ndata=0;
  if ( fit_data==1 ) {
    hdata = (TH1F*)hEGdata->Clone();
    ndata = (int)hdata->Integral();    
    for(int ibin=1; ibin<=hdata->GetNbinsX(); ibin++){
      for(int ipoint=0; ipoint<hdata->GetBinContent(ibin); ipoint++) {
	dataColl.push_back(hdata->GetBinCenter(ibin));
      }
    }
    ndata = dataColl.size();
  }else if (fit_data==2 ){
      hdata = (TH1F*)hEGdata->Clone();
    hdata -> Reset();
    dataColl.clear();
    FILE *infile =  fopen(dataFile.data(),"r");  
    float xdata, xdata1, xdata2; // combined isolation, pt, eta

    int flag = 1;
    while (flag!=-1){
      flag =fscanf(infile,"%f %f %f",&xdata, &xdata1, &xdata2);
      if( xdata1 >= ptbin && xdata1 < ptmax && xdata<20.) {
	if((strcmp(EBEE,"EB")==0 && TMath::Abs(xdata2)<1.45) ||
	   (strcmp(EBEE,"EE")==0 && TMath::Abs(xdata2)<2.5 && TMath::Abs(xdata2)>1.7) ) {
 	  dataColl.push_back(xdata);
	  hdata->Fill(xdata);
 	}
      } 
    }// keep reading files as long as text exists
    ndata = dataColl.size();

    printf("test print data 2 %2.3f \n", dataColl[2]);    
//     cout << "ndata in dataColl = " << ndata << endl;
    if ( ndata == 0 ) {
      printf(" no data to fit \n");
      return fitted;
    }
  }

  if(ndata==0) {
    printf(" ---  no events in the fit \n");
    return fitted;
  }
    
  //test fit the template and get PDFs
  TCanvas *c10 = new TCanvas("c10","c10",1000,500);
  c10->Divide(2,1);
  c10->cd(1);

  double par[20] = {hsig->GetMaximum(), 1., 0.6, 0.3,
		    hbkg->GetMaximum(), -.45, -0.05, 0.03, 1., 1., 1., 1.};

  if(strcmp(EBEE,"EE")==0) { par[2]=-0.1, par[3]=0.2; par[6]=-0.15; par[7]=0.02; };
  int fit_status;

  TF1 *f1 = new TF1("f1", exp_conv, -1., 20., 11);
  TF1 *fmcsigfit = new TF1("fmcsigfit", exp_conv, -1., 20., 11);
  fmcsigfit->SetLineColor(4);
  fmcsigfit->SetLineWidth(2);

  f1->SetNpx(10000);
  f1->SetParameters(par);
  f1->SetLineWidth(2);
  c10->cd(1);
  fit_status = hsig->Fit(f1,"","",-1., 5.);

  hsig->Draw();
  f1->Draw("same");
  if ( fit_status > 0 ) {
     printf("fit signal template failed. QUIT \n");
     return fitted;
  }
  if(para_index>0 && para_index<4){
    double tmppar = f1->GetParameter(para_index);
    f1->SetParameter(para_index, tmppar+para_sigma*f1->GetParError(para_index));
  }

  TF1 *fmcsig = (TF1*)f1->Clone();
  TF1 *fmcsigcorr = (TF1*)f1->Clone();
  fmcsig->SetNpx(10000);
  fmcsigcorr->SetNpx(10000);
  fmcsigfit->SetNpx(10000);
  
  TCanvas *c101 = new TCanvas("c101","c101",1000,500);
  c101->Divide(2,1);
  c101->cd(1);

  fmcsig->SetLineColor(1);
//   fmcsig->Draw();
//   f1->Draw("same");
  TH1F *htmp1 = (TH1F*)fmcsig->GetHistogram();
//   TH1F *htmp2 = (TH1F*)fmcsigcorr->GetHistogram();
  
  TH2F *htmp2 = new TH2F("htmp2","",210, -1., 20., 100, 0., htmp1->GetMaximum()*1.25);

  htmp2->SetNdivisions(505,"XY");
  htmp2->SetXTitle("Iso");
  htmp2->SetYTitle("A.U.");
  htmp2->SetLineColor(1);

//   htmp2->Draw();
//   htmp1->Draw("same");
//   htmp2->Add(htmp1,-1);
//   htmp2->Divide(htmp1);
  htmp2->GetXaxis()->SetRangeUser(-1., 10.);
  htmp2->SetMinimum(-1.);
  //htmp2->SetMaximum(1.5);
  htmp2->Draw();
  fmcsig->Draw("same");
//   fmcsigcorr->Draw("same");
  
  TLegend *tleg1 = new TLegend(0.5, 0.7, 0.93, 0.92);
  tleg1->SetHeader("");
  tleg1->SetFillColor(0);
  tleg1->SetShadowColor(0);
  tleg1->SetBorderSize(0);
  tleg1->AddEntry(fmcsig,"Zee data","l");
  //tleg1->AddEntry(fmcsigcorr,"corrected shape","l");
  tleg1->AddEntry(fmcsigfit,"shape from data","l");
  tleg1->Draw();

  //return fitted;
       
  SigPDFnorm = f1->Integral(-1.,20.);
  printf("status %d, sig area %3.3f \n", fit_status,f1->Integral(-1., 20.));


  f1->SetParameter(2,f1->GetParameter(2)+0.2);
  f1->SetParameter(3,f1->GetParameter(3)+0.1);

  Para.push_back(f1->GetParameter(0));
  Para.push_back(f1->GetParameter(1));
  Para.push_back(f1->GetParameter(2));
  Para.push_back(f1->GetParameter(3));

  Para_err.push_back(f1->GetParError(0));
  Para_err.push_back(f1->GetParError(1));
  Para_err.push_back(f1->GetParError(2));
  Para_err.push_back(f1->GetParError(3));

  c10->cd(2);
  TF1 *fbkgfit = new TF1("fbkgfit", expinv_power, -1., 20., 11);  

  TF1 *f3 = new TF1("f3", expinv_power, -1., 20., 11);
  fbkgfit->SetNpx(10000);  
  fbkgfit->SetLineColor(4);
  fbkgfit->SetLineWidth(2);

  f3->SetNpx(10000);
  f3->SetLineWidth(2);
  f3->SetParameters(f1->GetParameters());
    
  f3->SetParLimits(5,-5.,0.);
  f3->SetParLimits(6,-0.5,0.);
  f3->SetParLimits(7,0.001,0.2);
  f3->SetParLimits(8,0.5,5.);
  if ( strcmp(EBEE,"EB")==0 ){  
//     f3->FixParameter(8,1.);
//     f3->FixParameter(6,-0.1);
    f3->SetParLimits(8,1.,1.5);
  }

  float bkg_bend_power = 1.;
  if ( ptbin==21 ) bkg_bend_power = 4.5;
  if ( ptbin==23 ) bkg_bend_power = 4.;
  if ( ptbin==26 ) bkg_bend_power = 3.5;
  if ( ptbin==30 ) bkg_bend_power = 2.6;
  if ( ptbin==35 ) bkg_bend_power = 2.2;
  if ( ptbin==40 ) bkg_bend_power = 2.;
  if ( ptbin==45 ) bkg_bend_power = 2.;
  if ( ptbin==50 ) bkg_bend_power = 1.8;
  if ( ptbin==60 ) bkg_bend_power = 1.5;
  if ( ptbin==85 ) bkg_bend_power = 1.;
  if ( ptbin==120 ) bkg_bend_power = 1.;


  if ( strcmp(EBEE,"EE")==0 ){  
    f3->SetParameter(8,bkg_bend_power);
    f3->SetParLimits(8,bkg_bend_power-1., bkg_bend_power+1.);
  }

  f3->FixParameter(0,f3->GetParameter(0));
  f3->FixParameter(1,f3->GetParameter(1));
  f3->FixParameter(2,f3->GetParameter(2));
  f3->FixParameter(3,f3->GetParameter(3));

  hbkg->SetMaximum(hbkg->GetMaximum()*3.);
  fit_status = hbkg->Fit(f3,"b","",-1., 20.);
  hbkg->Draw();
  if ( fit_status > 0 ) {
    printf("fit background template failed. QUIT \n");    
    return fitted;
  }else {
    f3->Draw("same");
  }

  TF1 *fmcbkg = (TF1*)f3->Clone();
  fmcbkg->SetLineColor(1);
  c101->cd(2);

  htmp1 = (TH1F*)fmcbkg->GetHistogram();
  htmp2 = new TH2F("htmp2","",210, -1., 20., 100, 0., htmp1->GetMaximum()*1.25);

  htmp2->SetNdivisions(505,"XY");
  htmp2->SetXTitle("Iso");
  htmp2->SetYTitle("A.U.");
  htmp2->SetLineColor(1);

  htmp2->GetXaxis()->SetRangeUser(-1., 20.);
  htmp2->SetMinimum(-1.);
  htmp2->SetMaximum(1.5);
  htmp2->Draw();
  fmcbkg->Draw("same");

  TLegend *tleg2 = new TLegend(0.25, 0.2, 0.6, 0.42);
  tleg2->SetHeader("");
  tleg2->SetFillColor(0);
  tleg2->SetShadowColor(0);
  tleg2->SetBorderSize(0);
  if ( strcmp(EBEE,"EB")==0 ){  
    tleg2->AddEntry(fmcbkg,"MC shape","l");
  }else {
    tleg2->AddEntry(fmcbkg,"Data SB shape","l");
  }
  tleg2->AddEntry(fbkgfit,"shape from data","l");
  tleg2->Draw();
  
  if(para_index>4){
    double tmppar = f3->GetParameter(para_index);
    f3->SetParameter(para_index, tmppar+para_sigma*f3->GetParError(para_index));
  }

//   f3->SetParameter(5,-0.5);
//   f3->SetParameter(6,-0.05);
//   f3->SetParameter(7,0.02);
//   f3->SetParameter(8,1.);

  Para.push_back(f3->GetParameter(4));
  Para.push_back(f3->GetParameter(5));
  Para.push_back(f3->GetParameter(6));
  Para.push_back(f3->GetParameter(7)); 
  Para.push_back(f3->GetParameter(8)); 

  Para_err.push_back(f3->GetParError(4));
  Para_err.push_back(f3->GetParError(5));
  Para_err.push_back(f3->GetParError(6));
  Para_err.push_back(f3->GetParError(7));
  Para_err.push_back(f3->GetParError(8));

  BkgPDFnorm = f3->Integral(-1., 20.);
  printf("status %d, bkg area %3.3f \n", fit_status,f3->Integral(-1., 20.)/hdata->GetBinWidth(2));

  //test PDFs
  TCanvas *c11 = new TCanvas("c11","c11",1000,500);
  c11->Divide(2,1);
  c11->cd(1);
  TF1 *f11 = new TF1("f11",exp_conv_norm, -1., 20., 11);
  f11->SetNpx(10000);
  f11->SetParameters(f3->GetParameters());
  f11->Draw();
  printf(" SIG PDF area %2.3f \n", f11->Integral(-1., 20.));

  c11->cd(2);
  TF1 *f12 = new TF1("f12",expinv_power_norm, -1., 20., 11);
  f12->SetNpx(10000);
  f12->SetParameters(f3->GetParameters());
  f12->Draw();
  printf(" BKG PDF area %2.3f \n", f12->Integral(-1., 20.));

  //c1->cd();

  printf(" --------- before the fit ------------- \n");
  printf("Nsig %2.3f, Nbg %2.3f, Ntemplate %3.3f \n", hsig->Integral(), hbkg->Integral(), ntemplate);
  printf("Purity %2.3f, init size %4.3f,  fit sample size %4d\n", hsig->Integral()/hsum->Integral(), hsum->Integral(), ndata);
  printf(" -------------------------------------- \n");



  printf( " -----  Got %d, %d, %d events for fit ----- \n ", dataColl.size(),
	  sigColl.size(), bkgColl.size() );  

  //--------------------------------------------------
  //init parameters for fit
  Double_t vstart[11] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
  vstart[0] = sigfrac*ndata;
  vstart[1] = (1-sigfrac)*ndata;
  for (int ii=0; ii<9; ii++) {    
    vstart[ii+2] = Para[ii]; //8 shape parameters
  }
  TMinuit *gMinuit = new TMinuit(NPAR);  
  gMinuit->Command("SET STR 1");
  gMinuit->SetFCN(fcn);
  Double_t arglist[11];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  arglist[0] = 1;
  gMinuit->mnexcm("SET PRINT", arglist ,1,ierflg);

  Double_t step[] = { 1.,1.,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,};

  for ( int ii=0; ii<9; ii++){
    printf(" para %d, %.5f, err %.5f \n", ii, Para[ii], Para_err[ii]);
  }

  float sigma = 3.;
  gMinuit->mnparm(0,  "Signal yield"  , vstart[0],  step[0], 0., ndata*2.  , ierflg);
  gMinuit->mnparm(1,  "background yield"  , vstart[1],  step[1], 0., ndata*2. , ierflg);

//   gMinuit->mnparm(2,  "constant"     , Para[0],  0.00,  Para[0], Para[0] , ierflg);
//   gMinuit->mnparm(3,  "exp tail"     , Para[1],  0.01,  Para[1]-sigma*Para_err[1], Para[1]+sigma*Para_err[1], ierflg);
//   gMinuit->mnparm(4,  "exg mean"     , Para[2],  0.01,  Para[2]-sigma*Para_err[2], Para[2]+sigma*Para_err[2], ierflg);
//   gMinuit->mnparm(5,  "exg width"    , Para[3],  0.01,  Para[3]-sigma*Para_err[3], Para[3]+sigma*Para_err[3], ierflg);
//   gMinuit->mnparm(6,  "constant"     , Para[4],  0.00,  Para[4]                  , Para[4]                  , ierflg);
//   gMinuit->mnparm(7,  "bg exp turnon", Para[5],  0.01,  Para[5]-sigma*Para_err[5], Para[5]+sigma*Para_err[5], ierflg);
//   gMinuit->mnparm(8,  "bg x offset  ", Para[6],  0.01,  Para[6]-sigma*Para_err[6], Para[6]+sigma*Para_err[6], ierflg);
//   gMinuit->mnparm(9,  "bg bend slope", Para[7],  0.01,  0.001                    , 0.1                      , ierflg);
// //   gMinuit->mnparm(10, "bg bend power", Para[8],  0.01,  Para[8]-sigma*Para_err[8], Para[8]+sigma*Para_err[8], ierflg);
//   gMinuit->mnparm(10, "bg bend power", Para[8],  0.01,  0.5                       , 5.                       , ierflg);

//   gMinuit->mnparm(2,  "constant"     , Para[0], TMath::Abs(Para[0]*0.0) ,  Para[0], Para[0], ierflg);
//   gMinuit->mnparm(3,  "exp tail"     , Para[1], TMath::Abs(Para[1]*0.01) ,  Para[1]-sigma*Para_err[1], Para[1]+sigma*Para_err[1], ierflg); 
// //   gMinuit->mnparm(3,  "exp tail"     , Para[1], TMath::Abs(Para[1]*0.1) ,  0.8    , 1.3    , ierflg);
//   gMinuit->mnparm(4,  "exg mean"     , Para[2], TMath::Abs(Para[2]*0.1) ,  0.5    , 1.0    , ierflg);
//   gMinuit->mnparm(5,  "exg width"    , Para[3], TMath::Abs(Para[3]*0.1) ,  0.25   , 0.5    , ierflg);
//   gMinuit->mnparm(6,  "constant"     , Para[4], TMath::Abs(Para[4]*0.0) ,  Para[4], Para[4], ierflg);
//   gMinuit->mnparm(7,  "bg exp turnon", Para[5], TMath::Abs(Para[5]*0.1) ,  -0.7   , -0.3   , ierflg);
//   gMinuit->mnparm(8,  "bg x offset  ", Para[6], TMath::Abs(Para[6]*0.0) ,  -0.15  , -0.05  , ierflg);
//   gMinuit->mnparm(9,  "bg bend slope", Para[7], TMath::Abs(Para[7]*0.1) ,  0.01   , 0.05   , ierflg);
//   gMinuit->mnparm(10, "bg bend power", Para[8], TMath::Abs(Para[8]*0.1) ,  0.5    , 1.5    , ierflg);

  gMinuit->mnparm(2,  "constant"     , Para[0],  0.00,  Para[0], Para[0] , ierflg);
  gMinuit->mnparm(3,  "exp tail"     , Para[1],  0.00,  Para[1]-sigma*Para_err[1], Para[1]+sigma*Para_err[1], ierflg);
  gMinuit->mnparm(4,  "exg mean"     , Para[2],  0.00,  Para[2]-sigma*Para_err[2], Para[2]+sigma*Para_err[2], ierflg);
  gMinuit->mnparm(5,  "exg width"    , Para[3],  0.00,  Para[3]-sigma*Para_err[3], Para[3]+sigma*Para_err[3], ierflg);
  gMinuit->mnparm(6,  "constant"     , Para[4],  0.00,  Para[4]                  , Para[4]                  , ierflg);
  gMinuit->mnparm(7,  "bg exp turnon", Para[5],  0.00,  Para[5]-sigma*Para_err[5], Para[5]+sigma*Para_err[5], ierflg);
  gMinuit->mnparm(8,  "bg x offset  ", Para[6],  0.00,  Para[6]-sigma*Para_err[6], Para[6]+sigma*Para_err[6], ierflg);
  gMinuit->mnparm(9,  "bg bend slope", Para[7],  0.00,  0.001                    , 0.1                      , ierflg);
  gMinuit->mnparm(10, "bg bend power", Para[8],  0.00,  Para[8]-sigma*Para_err[8], Para[8]+sigma*Para_err[8], ierflg);
  
  printf(" --------------------------------------------------------- \n");
  printf(" Now ready for minimization step \n --------------------------------------------------------- \n");
  
  arglist[0] = 500; // number of iteration
  gMinuit->mnexcm("MIGRAD", arglist,1,ierflg);
  //can do scan
//   arglist[0] = 0;
//   gMinuit->mnexcm("SCAN", arglist,1,ierflg);

  printf (" -------------------------------------------- \n");
  printf("Finished.  ierr = %d \n", ierflg);

  double para[NPAR+1],errpara[NPAR+1];

  double tmp_errpara[NPAR+1];

  for(int j=0; j<=NPAR-1;j++) { tmp_errpara[j]=0.1; }
  for(int j=2; j<=NPAR-1;j++) { 
    if(Para_err[j-2]!=0.) tmp_errpara[j]=TMath::Abs(Para_err[j-2]); 
  }
  
  int ni=6;       if ( strcmp(EBEE,"EE")==0 ) { ni=6; }//if(ptbin==21) ni=0;}
  
  if ( ierflg == 0 ) {
    for(int i=0; i<ni; i++) {
      float istep[10] = {0.,0.,0.,0.,0.,0.,0.};
      if (i<(ni-1)) {
	istep[i] = 0.001;
      }else {
	for (int j=0; j<ni-1; j++) {istep[j] = 0.001;}
      }

      for(int j=0; j<=NPAR-1;j++) {
	gMinuit->GetParameter(j, para[j], errpara[j]);
	if (errpara[j] != 0. ) {
	  tmp_errpara[j] = TMath::Abs(errpara[j]);
	}
      }

      if ( strcmp(EBEE,"EB")==0 ) {

	sigma = 10.;

 	if ( i==(ni-1) ) { sigma=5.;istep[1]=istep[4]=0.; }
	if ( ptbin==21 && i==1 ){ sigma=3.; }
 	if ( ptbin==21 && i==(ni-1) ){ sigma=20.; }
	if ( ptbin==23 && i==0 ){ para[7]=-0.5; }
	if ( ptbin==23 && i==1 ){ istep[1]=0.; istep[3]=0.01; }
 	if ( ptbin==23 && i==3 ){ istep[1]=0.01; istep[3]=0.0; }
	if ( ptbin==23 && i==(ni-1) ){ sigma=20.; }
 	if ( ptbin==26 && i==1 ){ sigma=5.; }	
	if ( ptbin==26 && i==(ni-1) ){ sigma=20.; }
	if ( ptbin==30 && i==(ni-1) ){ sigma=3.; }
 	if ( ptbin==35 && i==(ni-1) ) { sigma=10.; }
 	if ( ptbin==40 && i==(ni-1) ) { sigma=5.; istep[4]=0.01; }
 	if ( ptbin==45 && i==(ni-1) ) { sigma=10.; }
	if ( ptbin==60 && i==0 ) { para[3]=1.; para[4]=0.6; para[5]=0.32; para[7]=-0.45; para[9]=0.025; para[10] = 1.;}
 	if ( ptbin==60 && i==(ni-1) ) { sigma=5.; istep[4]=0.01;}
	if ( ptbin>=85 && i==(ni-1) ){ sigma=3.; }
	if ( ptbin==300 ) { istep[2]=istep[3]=istep[4]=0.; }// para[7] = -5.11907e-02; istep[1]=0.; }
	float tmp8=0.;
	
// 	if( i!= (ni-1) ) {
	  gMinuit->mnparm(0,  "Signal yield"  ,   para[0],  1., para[0]-100.*tmp_errpara[0], para[0]+100.*tmp_errpara[0], ierflg);
	  gMinuit->mnparm(1,  "background yield", para[1],  1., para[1]-100.*tmp_errpara[1], para[1]+100.*tmp_errpara[1], ierflg);
	  gMinuit->mnparm(2,  "constant"     , para[2],  0., para[2]-100.*tmp_errpara[2], para[2]+100.*tmp_errpara[2], ierflg);
	  gMinuit->mnparm(6,  "constant"     , para[6],  0., para[6]-100.*tmp_errpara[6], para[6]+100.*tmp_errpara[6], ierflg);
	  gMinuit->mnparm(3,  "exp tail"     , para[3],  istep[4],  para[3]-sigma*tmp_errpara[3], para[3]+sigma*tmp_errpara[3], ierflg);
	  gMinuit->mnparm(4,  "exg mean"     , para[4],  istep[3],  para[4]-sigma*tmp_errpara[4], para[4]+sigma*tmp_errpara[4], ierflg);
	  gMinuit->mnparm(5,  "exg width"    , para[5],  istep[2],  para[5]-sigma*tmp_errpara[5], para[5]+sigma*tmp_errpara[5], ierflg);
	  gMinuit->mnparm(7,  "bg exp turnon", para[7],  istep[1],  para[7]-sigma*tmp_errpara[7], para[7]+sigma*tmp_errpara[7], ierflg);
	  gMinuit->mnparm(8,  "bg x offset  ", para[8],  tmp8    ,  para[8]-sigma*tmp_errpara[8], para[8]+sigma*tmp_errpara[8], ierflg);
	  gMinuit->mnparm(9,  "bg bend slope", para[9],  istep[0],  para[9]-sigma*tmp_errpara[9], para[9]+sigma*tmp_errpara[9], ierflg);      
	  float sigma10=5.;
	  if ( para[10]-sigma10*tmp_errpara[10] < 1. )// && i!=(ni-1))
	    gMinuit->mnparm(10, "bg bend power", para[10],  istep[0], 1.,  para[10]+sigma10*tmp_errpara[10], ierflg);      
	  else
	    gMinuit->mnparm(10, "bg bend power", para[10],  istep[0], para[10]-sigma10*tmp_errpara[10],  para[10]+sigma10*tmp_errpara[10], ierflg);      
// 	}else {
// 	  gMinuit->mnparm(2,  "constant"     , Para[0], TMath::Abs(Para[0]*0.0) ,  Para[0], Para[0], ierflg);
// 	  //gMinuit->mnparm(3,  "exp tail"     , Para[1], TMath::Abs(Para[1]*0.01) ,  Para[1]-sigma*Para_err[1], Para[1]+sigma*Para_err[1], ierflg); 
// 	  gMinuit->mnparm(3,  "exp tail"     , Para[1], TMath::Abs(Para[1]*0.0) ,  0.8    , 1.3    , ierflg);
// 	  gMinuit->mnparm(4,  "exg mean"     , Para[2], TMath::Abs(Para[2]*0.1) ,  0.5    , 1.0    , ierflg);
// 	  gMinuit->mnparm(5,  "exg width"    , Para[3], TMath::Abs(Para[3]*0.1) ,  0.25   , 0.5    , ierflg);
// 	  gMinuit->mnparm(6,  "constant"     , Para[4], TMath::Abs(Para[4]*0.0) ,  Para[4], Para[4], ierflg);
// 	  gMinuit->mnparm(7,  "bg exp turnon", Para[5], TMath::Abs(Para[5]*0.0) ,  -0.7   , -0.3   , ierflg);
// 	  gMinuit->mnparm(8,  "bg x offset  ", Para[6], TMath::Abs(Para[6]*0.0) ,  -0.15  , -0.05  , ierflg);
// 	  gMinuit->mnparm(9,  "bg bend slope", Para[7], TMath::Abs(Para[7]*0.1) ,  0.01   , 0.05   , ierflg);
// 	  gMinuit->mnparm(10, "bg bend power", Para[8], TMath::Abs(Para[8]*0.1) ,  0.5    , 1.5    , ierflg);
// 	}


	if( ptbin >=300 ) { 
	  gMinuit->mnparm(3,  "exp tail"  , 1.257281,  0.0,  para[1]-3.*tmp_errpara[1], para[1]+3.*tmp_errpara[1], ierflg);
	  gMinuit->mnparm(4,  "exg mean"  , 0.856906,  0.0,  para[2]-3.*tmp_errpara[2], para[2]+3.*tmp_errpara[2], ierflg);
	  gMinuit->mnparm(5,  "exg width" , 0.320847,  0.0,  para[3]-3.*tmp_errpara[3], para[3]+3.*tmp_errpara[3], ierflg);
	}      

    }else{	

	sigma=10.;
	if ( i==0 ) { para[10] = bkg_bend_power; tmp_errpara[10] = 0.3; }
 	if ( i==(ni-1) ) { sigma=3.;istep[1]=istep[4]=0.; } //test of not changing signal template
     	if ( i==(ni-1) ) { istep[4]=0.;}

   	if ( ptbin==21 && i==(ni-1) ) { sigma=20.;}
  	if ( ptbin==23 && i==0 ) { sigma=5.;}
  	if ( ptbin==23 && i==(ni-1) ) { sigma=10.;}
	if ( ptbin<30 && ptbin>21 && i==1 ){ istep[1]=0.; istep[3]=0.01; }
 	if ( ptbin<30 && ptbin>21 && i==3 ){ istep[1]=0.01; istep[3]=0.0; }
	if ( ptbin==26 && i==1 ) { para[7] = -0.8; }
	if ( ptbin==26 && i==(ni-1) ) { sigma=10.; }
  	if ( ptbin==30 && i==(ni-1) ) { sigma=10.; }
 	if ( ptbin==35) {para[7] = -0.75;}
 	if ( ptbin==40 && i==0) {para[7] = -0.65; para[10] = 2.;}
	if ( ptbin==45 && i==(ni-1) ) {sigma=5.;}
	if ( ptbin==85 && i==(ni-1) ) {sigma=10.; istep[4]=0.01;}
	if (ptbin >= 85 ) { para[10] = bkg_bend_power; tmp_errpara[10] = 1.; }

	if ( ptbin==120 ) { para[7] = -0.615255; istep[1]=0.;}

	
//     	if ( ptbin==120 && i==0 ) { 
// 	  para[3] = 1.446454; para[4]=-0.016373; para[5]=0.163238;
// 	  istep[2]=istep[3]=istep[4]=0.; sigma=5.; tmp_errpara[10]=0.2;
// 	}
//     	if ( ptbin==120 && i==(ni-1) ) { istep[2]=istep[3]=istep[4]=0.; sigma=5.;}

	gMinuit->mnparm(0,  "Signal yield"  ,   para[0],  1., para[0]-100.*tmp_errpara[0], para[0]+100.*tmp_errpara[0], ierflg);
	gMinuit->mnparm(1,  "background yield", para[1],  1., para[1]-100.*tmp_errpara[1], para[1]+100.*tmp_errpara[1], ierflg);
	gMinuit->mnparm(2,  "constant"     , para[2],  0.,  para[2], para[2] , ierflg);
	gMinuit->mnparm(6,  "constant"     , para[6],  0.,  para[6], para[6], ierflg);	
	gMinuit->mnparm(3,  "exp tail"     , para[3],  istep[4],  para[3]-sigma*tmp_errpara[3], para[3]+sigma*tmp_errpara[3], ierflg);
	gMinuit->mnparm(4,  "exg mean"     , para[4],  istep[3],  para[4]-sigma*tmp_errpara[4], para[4]+sigma*tmp_errpara[4], ierflg);
	gMinuit->mnparm(5,  "exg width"    , para[5],  istep[2],  para[5]-sigma*tmp_errpara[5], para[5]+sigma*tmp_errpara[5], ierflg);
	gMinuit->mnparm(7,  "bg exp turnon", para[7],  istep[1],  para[7]-sigma*tmp_errpara[7], para[7]+sigma*tmp_errpara[7], ierflg);
	gMinuit->mnparm(8,  "bg x offset  ", para[8],  0.00,      para[8]-sigma*tmp_errpara[8], para[8]+sigma*tmp_errpara[8], ierflg);
	gMinuit->mnparm(9,  "bg bend slope", para[9],  istep[0],  para[9]-sigma*tmp_errpara[9], para[9]+sigma*tmp_errpara[9], ierflg);	
  
	float minerr=1.;
	//if ( tmp_errpara[10] > 0.5) tmp_errpara[10] = 0.5;
	float sigma10=5.;
	if ( para[10]-sigma10*tmp_errpara[10] < 1. ) 
	  gMinuit->mnparm(10, "bg bend power", para[10],  istep[0], minerr,  para[10]+sigma10*tmp_errpara[10], ierflg);
	else 
	  gMinuit->mnparm(10, "bg bend power", para[10],  istep[0], para[10]-sigma10*tmp_errpara[10],  para[10]+sigma10*tmp_errpara[10], ierflg);

      }
      printf(" ************ \n");
      printf("  do %d th fit  \n", i);
      if(i==5 && dataFile.find("toy")    != std::string::npos)
	{
	  cout << "dataResult = " << dataResult << "\t dataErr = " << dataErr << endl;
	  // fixed turn on at +- 1 sigma
	  gMinuit->mnparm(7,  "bg exp turnon", dataResult-(float)shift*dataErr,  0.00,  para[7]-sigma*tmp_errpara[7], para[7]+sigma*tmp_errpara[7], ierflg);

	}
      else if(dataFile.find("toy")    == std::string::npos)
	{
	  dataResult = para[7];
	  dataErr    = tmp_errpara[7];
	}
      arglist[0] = 500; // number of iteration
      gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);      
      if ( ierflg != 0 ) {
  	printf("fit failed at %d iteration \n", i);
  	c1->cd();	c1->Draw();  	hdata->Draw("phe");
  	return fitted;
      }
    }
  }
 
  Double_t amin,edm,errdef; 
  if ( ierflg == 0 ) {
    for(int j=0; j<=NPAR-1;j++) {
      gMinuit->GetParameter(j, para[j],errpara[j]);
      info.push_back(para[j]);
      info_err.push_back(errpara[j]);
      printf("Parameter  %d = %f +- %f\n",j,para[j],errpara[j]);	
    }
    para[NPAR] = dataColl.size();
    printf(" fitted yield %2.3f \n", (para[0]+para[1])/ndata );
    
    info.push_back(sigColl.size());
    
    for(int j=0; j<=NPAR-1;j++) {
      tmp_errpara[j] = errpara[j];
      if( tmp_errpara[j] == 0. ) tmp_errpara[j] = par[j]*.1;      
    }
    //do minos if fit sucessed.

  }
  if (ierflg != 0 )  {
    printf(" *********** Fit failed! ************\n");
    gMinuit->GetParameter(0, para[0],errpara[0]);
    gMinuit->GetParameter(1, para[1],errpara[1]);
    para[0]=0.; errpara[0]=0.;

    c1->cd();
    c1->Draw();  
    //gPad->SetLogy();
    hdata->SetNdivisions(505,"XY");
    hdata->SetXTitle("comb. ISO (GeV)");
    hdata->SetYTitle("Entries");
    hdata->SetTitle("");
    hdata->SetMarkerStyle(8);
    hdata->SetMinimum(0.);
    if ( hdata->GetMaximum()<10.) hdata->SetMaximum(15.);
    else hdata->SetMaximum(hdata->GetMaximum()*1.25);
    if ( strcmp(EBEE,"EE")==0 &&ptbin == 15 ) hdata->SetMaximum(hdata->GetMaximum()*1.25);
   
    hdata->Draw("phe");  

    return fitted;    
  }

  
  // Print results
//   Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  gMinuit->mnprin(1,amin);  
  gMinuit->mnmatu(1);
  printf(" ========= happy ending !? =========================== \n");
  
  printf("FCN =  %3.3f \n", amin);

  //use new PDF form
  double tmppar[12];
  for(int ii=0; ii<9; ii++){
    tmppar[ii] = para[ii+2];
    fmcsigfit->SetParameter(ii,tmppar[ii]);
    fbkgfit->SetParameter(ii,tmppar[ii]);
  }

  c101->cd(1);
  
  //fmcsigfit->SetParameters(tmppar);
  //fmcsigfit->SetParameter(2,0.1);
  //fmcsigfit->SetLineStyle(2);

  fmcsigfit->Draw("same");
  c101->cd(2);

  fbkgfit->SetParameter(4,fbkgfit->GetParameter(4)*fmcbkg->Integral(-1., 20.)/fbkgfit->Integral(-1., 20.));
  fbkgfit->Draw("same");

  char fname[100];
  sprintf(fname,"plots/template_Ifit%s_%d.pdf",EBEE,ptbin);
  c101->SaveAs(fname);


  f11->SetParameters(tmppar);
  SigPDFnorm = f11->Integral(-1., 20.);
  f12->SetParameters(tmppar);
  BkgPDFnorm = f12->Integral(-1., 20.);


  // plot
  c1->cd();
  c1->Draw();  
  //gPad->SetLogy();
   hdata->SetNdivisions(505,"XY");
   hdata->SetXTitle("comb. ISO (GeV)");
   hdata->SetYTitle("Entries");
   hdata->SetTitle("");
   hdata->SetMarkerStyle(8);
   hdata->SetMinimum(0.);
   if ( hdata->GetMaximum()<10.) hdata->SetMaximum(15.);
   else hdata->SetMaximum(hdata->GetMaximum()*1.5);
   if ( strcmp(EBEE,"EE")==0 &&ptbin == 15 ) hdata->SetMaximum(hdata->GetMaximum()*1.2);

   hdata->Draw("p e ");

  f11->SetParameter(0, para[0]*f11->GetParameter(0)/f11->Integral(-1., 20.)*hdata->GetBinWidth(2));
//   f11->SetFillColor(5);
  f11->SetLineColor(4);
  //f11->SetFillColor(603);
  f11->SetLineWidth(2);
//   f11->SetFillStyle(3001);
  f11->Draw("same");

  f12->SetParameter(4, para[1]*f12->GetParameter(4)/f12->Integral(-1., 20.)*hdata->GetBinWidth(2));
//   f12->SetFillColor(8);
  f12->SetLineColor(2);
  //f12->SetFillColor(603);
  f12->SetLineWidth(2);
//   f12->SetFillStyle(3013);
  f12->Draw("same");

  TF1 *f13 = new TF1("f13",sum_norm, -1., 20 ,11);
  f13->SetNpx(10000);
  f13->SetParameters(f12->GetParameters());
  f13->SetParameter(0, para[0]*f11->GetParameter(0)/f11->Integral(-1., 20.)*hdata->GetBinWidth(2));
  f13->SetParameter(4, para[1]*f12->GetParameter(4)/f12->Integral(-1., 20.)*hdata->GetBinWidth(2));  
  f13->SetLineWidth(2);
  f13->SetLineColor(1);
  f13->Draw("same");
  f11->Draw("same");
  hdata->Draw("pe same");

//   cout << "The number of bins are: " << endl;
//   cout << "hdata nbins = " << hdata->GetNbinsX() << endl;
//   cout << "hsig nbins = " << hsig->GetNbinsX() << endl;
//   cout << "hbkg nbins = " << hbkg->GetNbinsX() << endl;

  // get chi2/NDF
  double chi2ForThisBin=0;
  int nbinForThisBin=0;
  chi2Nbins(f13, hdata, chi2ForThisBin, nbinForThisBin);
  for(int epar=0; epar < 11; epar++)
    {
//       cout << "f11 parameter " << epar << " = " << 
// 	f11->GetParameter(epar) << endl;
      FitPar[epar] = f11->GetParameter(epar);
    }

  for(int epar=0; epar < 11; epar++)
    {
//       cout << "f12 parameter " << epar << " = " << 
// 	f12->GetParameter(epar) << endl;
      FitPar[epar+11] = f12->GetParameter(epar);
    }

  for(int epar=0; epar < 11; epar++)
    {
//       cout << "f13 parameter " << epar << " = " << 
// 	f13->GetParameter(epar) << endl;
      FitPar[epar+22] = f13->GetParameter(epar);

    }

//   cout << "hdata integral = " << hdata->Integral() << endl;
//   cout << endl;

//   printf("fit area %3.2f; sig area %3.2f; bg area %3.2f\n", f13->Integral(-1., 20.)/hdata->GetBinWidth(2),  f11->Integral(-1., 20.)/hdata->GetBinWidth(2),f12->Integral(-1., 20.)/hdata->GetBinWidth(2));

//   for(int i=0; i<12; i++){
//     printf(" fit para %d = %4.3f \n", i, f13->GetParameter(i));
//   }

   TLegend *tleg = new TLegend(0.5, 0.7, 0.93, 0.92);
   char text[50];
   sprintf(text,"%s Pt %d ~ %.0f GeV",EBEE, ptbin, ptmax);
   tleg->SetHeader(text);
   tleg->SetFillColor(0);
   tleg->SetShadowColor(0);
   tleg->SetBorderSize(0);
   sprintf(text,"#chi^{2}/NDF = %.1f/%d",chi2ForThisBin,nbinForThisBin);
   tleg->AddEntry(hdata,text,"");
   sprintf(text,"Data %.1f events",hdata->Integral());
   tleg->AddEntry(hdata,text,"pl");
   sprintf(text,"Fitted %.1f events",para[0]+para[1]);//f13->Integral(-1., 20.)/hdata->GetBinWidth(2));
   tleg->AddEntry(f13,text,"l");
   sprintf(text,"SIG %.1f #pm %.1f events",para[0], errpara[0]);
   tleg->AddEntry(f11,text,"f");
   sprintf(text,"BKG %.1f #pm %.1f events",para[1], errpara[1]);
   tleg->AddEntry(f12,text,"f");
   tleg->Draw();


   gPad->RedrawAxis();

   printf("%s, ptbin %d, Data %.1f events \n",EBEE, ptbin, hdata->Integral());
   printf("Fitted %.1f (in 5GeV) %.1f events \n",para[0]+para[1],f13->Integral(-1.,5.));
   printf("SIG %.1f #pm %.1f events \n",para[0], errpara[0]);
   printf("SIG (in 5GeV) %.1f #pm %.1f events \n",f11->Integral(-1.,5.)/hdata->GetBinWidth(2), f11->Integral(-1.,5.)*errpara[0]/para[0]/hdata->GetBinWidth(2));
   printf("BKG %.1f #pm %.1f events \n",para[1], errpara[1]);
   printf("BKG (in 5GeV) %.1f #pm %.1f events \n",f12->Integral(-1.,5.)/hdata->GetBinWidth(2), f12->Integral(-1.,5.)*errpara[1]/para[1]/hdata->GetBinWidth(2));
   
   float purity = f11->Integral(-1.,5.)/hdata->GetBinWidth(2)/(f11->Integral(-1.,5.)/hdata->GetBinWidth(2)+f12->Integral(-1.,5.)/hdata->GetBinWidth(2));
   float purity_err = purity*errpara[0]/para[0];
   printf("Purity (in 5GeV) %.3f #pm %.3f  \n", purity, purity_err);


//   hsig->Scale(para[0]/hsig->Integral());
//   hbkg->Scale(para[1]/hbkg->Integral());
//   hbkg->Add(hsig);

//   hsig->SetLineColor(1);
//   hsig->SetFillColor(5);
//   hsig->SetFillStyle(3001);

//   hbkg->SetLineWidth(2);


//   hsig->Draw("same");
//   hbkg->Draw("same");


  sprintf(fname,"plots/unbinned_free_Ifit%s_%d.pdf",EBEE,ptbin);
  if (para_index>0) sprintf(fname,"plots/unbinned_Ifit%s_%d_para%d_sigma%1.0f.pdf",EBEE,ptbin,para_index,para_sigma);
  if(Opt_SavePDF == 1) {
    c1->SaveAs(fname);


  } else {

   c1->Close();
   c10->Close();
   c101->Close();
   c11->Close();

  }

  printf("----- fit results with signal projection   ----------- \n");

  fitted[0] = para[0];
  fitted[1] = errpara[0];
  fitted[2] = para[1];
  fitted[3] = errpara[1];
  fitted[4] = f11->Integral(-1.,5.)/hdata->GetBinWidth(2);
  fitted[5] = f11->Integral(-1.,5.)*errpara[0]/para[0]/hdata->GetBinWidth(2);

  return fitted;
}

