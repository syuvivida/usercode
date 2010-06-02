#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TF1.h"
#include <fstream>
#include <iostream>
#include <TMinuit.h>
#include <vector>
#include <TMath.h>
#include "TVirtualFitter.h"
#include <TPaveText.h>
#include "TFile.h"
#include "chi2Nbins.h"


using namespace std;

#define NPAR 2


vector<Double_t> dataColl;
vector<Double_t> sigColl;
vector<Double_t> bkgColl;

vector<Double_t> info;
vector<Double_t> info_err;

vector<Double_t> Para;

Double_t SigPDFnorm = 0.;
Double_t BkgPDFnorm = 0.;


Double_t g(Double_t *v, Double_t *par)
{
  Double_t arg = 0;
  arg = (v[0] - par[2]) / par[3];
  
  //Double_t area = par[1]/(par[3]*sqrt(2*3.1415926));
  Double_t norm_area = 1/(par[3]*sqrt(2*3.1415926));
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
    tmppar[i] = Para[i];
  }
  for ( int i=0; i<dataColl.size(); i++ ) {
    //PDF for signal and background
    Double_t x = dataColl[i];
    Double_t Ls = exp_conv_norm(&x,tmppar);
    Double_t Lb = expinv_power_norm(&x,tmppar);
    //Get Log Likelihood		  
    Lsum += TMath::Log( (fs*Ls + fb*Lb) / (fs+fb) );
  }
  f=2*( -1*Lsum + (fs+fb) - Nevt*TMath::Log(fs+fb) );
}



//___________________________________________________________________________
Double_t* Ifit(TH1F* dataInput, TH1F* sigTemplate, TH1F* bkgTemplate, 
	       int fit_data=1, std::string dataText="EGdata_comb3Iso_et_0531.dat",
	       double etamin=-1., double etamax=-1.,
	       double ptmin=-1., double ptmax=-1.)
{

  std::string histoName = dataInput->GetName();

  TCanvas *c1 = new TCanvas("HF1", "Histos1", 0, 0, 600, 600);
  double count=0;
  dataColl.clear();
  sigColl.clear();
  bkgColl.clear();

  Para.clear();

  Double_t* fitted = new Double_t[8];
  fitted[0] = fitted[1] = fitted[2] = fitted[3] = 0.0;
  fitted[4] = fitted[5] = fitted[6] = fitted[7] = 0.0;

  TH1F* hsig = (TH1F*)sigTemplate->Clone();
  hsig->SetName("hsig");

  TH1F* hbkg = (TH1F*)bkgTemplate->Clone();
  hbkg->SetName("hbkg");

  hsig->SetLineColor(1);
  hbkg->SetLineColor(1);

  TH1F *hsum = (TH1F*)hsig->Clone();
  hsum->Add(hbkg,1);
  float ntemplate = 1.;
  if (hsum->Integral()>1.) ntemplate = hsum->Integral();
  float sigfrac = hsig->Integral()/ntemplate;

  TH1F *hsum_norm = (TH1F*)hsum->Clone();  
  hsum_norm->Scale(1./hsum->Integral());

  TH1F *hdata;
  int ndata=0;
  hdata = (TH1F*)dataInput->Clone();
  hdata -> SetName("hdata");

  // for weighted MC
  if(fit_data==0){
    ndata = hdata->Integral();    
    for(int ibin=1; ibin<=hdata->GetNbinsX(); ibin++){
      for(int ipoint=0; ipoint<hdata->GetBinContent(ibin); ipoint++) {
	dataColl.push_back(hdata->GetBinCenter(ibin));
      }
    }
  }
  // for real data
  else if(fit_data==1){
    FILE *infile =  fopen(dataText.data(),"r");  
    float xdata, xdata1, xdata2; // combined isolation, pt, eta

    cout << "Requiring data within pt = " << ptmin << "--- " << ptmax << endl;
    cout << "Requiring data within eta = " << etamin << "--- " << etamax << endl;
    int flag = 1;
    while (flag!=-1){
      flag =fscanf(infile,"%f %f %f",&xdata, &xdata1, &xdata2);
      if( xdata1 >= ptmin && xdata1 < ptmax && xdata<11. && 
	  fabs(xdata2) > etamin && fabs(xdata2) < etamax)
	{
 	  dataColl.push_back(xdata);
 	}
    } // keep reading files as long as text exists
    ndata = dataColl.size();
  }
  cout << "There are " << ndata << " data points " << endl;
  
  if(ndata==0) {
    printf(" ---  no evetns in the fit \n");
    return fitted;
  }
    
  //test fit the template and get PDFs
  TCanvas *c10 = new TCanvas("c10","",1000,500);
  c10->Divide(2,1);
  c10->cd(1);
  
  hsig->Scale(1./hsig->Integral()); 
  hbkg->Scale(1./hbkg->Integral());  
  
  double par[20] = {1., 1., 0.5, 0.3,
		    hbkg->GetMaximum(),-.3,-1., 0.01, 0.5, 0.01, 1., 1.};
  int fit_status;

  TF1 *f1 = new TF1("f1", exp_conv, -1., 11., 12);
  f1->SetParameters(par);

  c10->cd(1);
  fit_status = hsig->Fit(f1,"");
  hsig->Draw();
  SigPDFnorm = f1->Integral(-1.,11.);
  printf("status %d, sig area %3.3f \n", fit_status,f1->Integral(-1.,11.));
//   if ( fit_status > 0 ) {
//      printf("fit signal template failed. QUIT \n");
//      return fitted;
//   }

  Para.push_back(f1->GetParameter(0));
  Para.push_back(f1->GetParameter(1));
  Para.push_back(f1->GetParameter(2));
  Para.push_back(f1->GetParameter(3)); 

  c10->cd(2);
  
  TF1 *f3 = new TF1("f3", expinv_power, -1., 11., 12);
  f3->SetParameters(f1->GetParameters());
  f3->SetParLimits(5,-10.,1.);
  f3->SetParLimits(6,-1.,2.);
  f3->SetParLimits(7,0.,0.09);
//   f3->SetParLimits(8,0.4,0.6);

//   if(etamin > 1.55 && fabs(ptmin-15.)<1e-6 && 
//      fabs(ptmax-20.)<1e-6)
//     {
//       cout << "find EE in pt bin 15--20" << endl;
//       f3->FixParameter(5,-0.1);
//     }
  
  if(etamin > 1.55 && fabs(ptmin-50.)<1e-6 && 
     fabs(ptmax-80.)<1e-6)
    {
      cout << "find EE in pt bin 50--80" << endl;
      f3->FixParameter(5,-0.1);
    }
  if(etamax < 1.55 && fabs(ptmin-80.)<1e-6 && 
     fabs(ptmax-120.)<1e-6)
    {
      cout << "find EB in pt bin 80--120" << endl;
      f3->FixParameter(5,-0.1);
    }

//    f3->FixParameter(6,-0.05);
  f3->FixParameter(8,0.5);
  f3->FixParameter(0,f3->GetParameter(0));
  f3->FixParameter(1,f3->GetParameter(1));
  f3->FixParameter(2,f3->GetParameter(2));
  f3->FixParameter(3,f3->GetParameter(3));

  hbkg->SetMaximum(hbkg->GetMaximum()*3.);
  fit_status = hbkg->Fit(f3,"b");
  hbkg->Draw();
  printf("status %d, bkg area %3.3f \n", fit_status,f3->Integral(-1.,11.)/hdata->GetBinWidth(2));
//   if ( fit_status > 0 ) {
//     printf("fit background template failed. QUIT \n");
//     return fitted;
//   }

  Para.push_back(f3->GetParameter(4));
  Para.push_back(f3->GetParameter(5));
  Para.push_back(f3->GetParameter(6));
  Para.push_back(f3->GetParameter(7)); 
  Para.push_back(f3->GetParameter(8)); 
  BkgPDFnorm = f3->Integral(-1.,11);
  c10->SaveAs(Form("plots/PDFFit_%s.gif",dataInput->GetName()));
  

  //test PDFs
  TCanvas *c11 = new TCanvas("c11","",1000,500);
  c11->Divide(2,1);
  c11->cd(1);
  TF1 *f11 = new TF1("f11",exp_conv_norm, -1., 11., 12);
  f11->SetParameters(f3->GetParameters());
  f11->Draw();
  printf(" SIG PDF area %2.3f \n", f11->Integral(-1.,11.));

  c11->cd(2);
  TF1 *f12 = new TF1("f12",expinv_power_norm, -1., 11., 12);
  f12->SetParameters(f3->GetParameters());
  f12->Draw();
  printf(" BKG PDF area %2.3f \n", f12->Integral(-1.,11.));


  //======= Determine the ratio of isolation < 5 GeV to isolation < 11 GeV
  double purityMaxReach   = 5.0;
  double scale_signal     = f11->Integral(-1,purityMaxReach);
  double scale_background = f12->Integral(-1,purityMaxReach);

  cout << "scale_signal = " << scale_signal << endl;
  cout << "scale_background =" << scale_background << endl;


  c1->cd();
  printf(" --------- before the fit ------------- \n");
  printf("Nsig %2.3f, Nbg %2.3f, Ntemplate %3.3f \n", hsig->Integral(), hbkg->Integral(), ntemplate);
  printf("Purity %2.3f, init size %4.3f,  fit sample size %4d\n", hsig->Integral()/hsum->Integral(), hsum->Integral(), ndata);
  printf(" -------------------------------------- \n");


  printf( " -----  Got %d, %d, %d events for fit ----- \n ", dataColl.size(),
	  sigColl.size(), bkgColl.size() );  

  //========== Dump fit parameters
  for(unsigned int ipara = 0; ipara < Para.size(); ipara++)
    {
      cout << "Parameter " << ipara << " = " << Para[ipara] << endl;

    }


  //--------------------------------------------------
  //init parameters for fit
  Double_t vstart[10] = {1., 1.};
  vstart[0] = sigfrac*ndata;
  vstart[1] = (1-sigfrac)*ndata;
 
  TMinuit *gMinuit = new TMinuit(NPAR);  
  gMinuit->Command("SET STR 1");
  gMinuit->SetFCN(fcn);
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  arglist[0] = 1;
  gMinuit->mnexcm("SET PRINT", arglist ,1,ierflg);

  Double_t step[] = { 0.1, 0.1,};

  gMinuit->mnparm(0,  "Signal yield"  , vstart[0],  step[0], 0., ndata*2.  , ierflg);
  gMinuit->mnparm(1,  "background yield"  , vstart[1],  step[1], 0., ndata*2. , ierflg);
  

  printf(" --------------------------------------------------------- \n");
  printf(" Setting stragety = 2 \n ----------------------\n");
  
  arglist[0] = 2;
  gMinuit->mnexcm("SET STRAT", arglist ,1,ierflg);

  printf(" --------------------------------------------------------- \n");
  printf(" Now ready for minimization step \n ----------------------\n");

  arglist[0] = 2000; // number of iteration
  arglist[1] = 1.;
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
  printf (" -------------------------------------------- \n");
  printf("Finished.  ierr = %d \n", ierflg);

  info.clear();
  info_err.clear();

  double para[NPAR+1],errpara[NPAR+1];
  if ( ierflg == 0 ) 
    {
      for(int j=0; j<=NPAR-1;j++) {
        gMinuit->GetParameter(j, para[j],errpara[j]);
        para[NPAR] = dataColl.size();
        info.push_back(para[j]);
        info_err.push_back(errpara[j]);
        printf("Parameter (yeild) %d = %f +- %f\n",j,para[j],errpara[j]);
	
      }
      printf(" fitted yield %2.3f \n", (para[0]+para[1])/ndata );

      info.push_back(sigColl.size());

    }
  else {
    printf(" *********** Fit failed! ************\n");
    gMinuit->GetParameter(0, para[0],errpara[0]);
    gMinuit->GetParameter(1, para[1],errpara[1]);
    para[0]=0.; errpara[0]=0.;
  }

  
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  gMinuit->mnprin(1,amin);  
  gMinuit->mnmatu(1);
  printf(" ========= happy ending !? =========================== \n");
  
  printf("FCN =  %3.3f \n", amin);


  // plot
  c1->Draw();  
  //gPad->SetLogy();
   hdata->SetNdivisions(505,"XY");
   hdata->SetXTitle("comb. ISO (GeV)");
   hdata->SetYTitle("Entries");
   hdata->SetTitle("");
   hdata->SetMarkerStyle(8);
   hdata->SetMinimum(0.);
   hdata->SetMaximum(hdata->GetMaximum()*1.5);

   hdata->Draw("p e");

  f11->SetParameter(0, para[0]*f11->GetParameter(0)/f11->Integral(-1,11)*hdata->GetBinWidth(2));
  f11->SetFillColor(5);
  f11->SetLineColor(1);
  //f11->SetFillColor(603);
  f11->SetLineWidth(1);
  f11->SetFillStyle(3001);
  f11->Draw("same");

  f12->SetParameter(4, para[1]*f12->GetParameter(4)/f12->Integral(-1,11)*hdata->GetBinWidth(2));
  f12->SetFillColor(8);
  f12->SetLineColor(1);
  //f12->SetFillColor(603);
  f12->SetLineWidth(1);
  f12->SetFillStyle(3013);
  f12->Draw("same");

  TF1 *f13 = new TF1("f13",sum_norm, -1, 11,12);
  f13->SetParameters(f12->GetParameters());
  f13->SetParameter(0, para[0]*f11->GetParameter(0)/f11->Integral(-1,11)*hdata->GetBinWidth(2));
  f13->SetParameter(4, para[1]*f12->GetParameter(4)/f12->Integral(-1,11)*hdata->GetBinWidth(2));  
  f13->SetLineWidth(2);
  f13->SetLineColor(1);
  f13->Draw("same");
  

  printf("fit area %3.2f; sig area %3.2f; bg area %3.2f\n", f13->Integral(-1,11)/hdata->GetBinWidth(2),  f11->Integral(-1,11)/hdata->GetBinWidth(2),f12->Integral(-1,11)/hdata->GetBinWidth(2));

  char text[1000];
  // get chi2/NDF  
  double chi2ForThisBin=0;
  int nbinForThisBin=0;
  chi2Nbins(f13, hdata, chi2ForThisBin, nbinForThisBin);

  TPaveText *pavetex = new TPaveText(0.43, 0.87, 0.90, 0.92,"NDCBR");
  pavetex->SetBorderSize(0);
  pavetex->SetFillColor(0);
  pavetex->SetFillStyle(0);
  pavetex->SetLineWidth(3);
  pavetex->SetTextAlign(12);
  pavetex->SetTextSize(0.03);
  pavetex->AddText(Form("#chi^{2}/NDF=%.1f/%d",chi2ForThisBin, nbinForThisBin));
  pavetex->Draw();

  TLegend *tleg = new TLegend(0.43, 0.60, 0.90, 0.87);
  tleg->SetHeader(dataInput->GetTitle());

  tleg->SetTextSize(0.03);
  tleg->SetFillColor(0);
  tleg->SetShadowColor(0);
  tleg->SetBorderSize(0);
  sprintf(text,"Data %5.1f events",hdata->Integral());
  tleg->AddEntry(hdata,text,"pl");
  sprintf(text,"Fitted %5.1f events",f13->Integral(-1,11)/hdata->GetBinWidth(2));
  tleg->AddEntry(f13,text,"l");
  sprintf(text,"SIG %5.1f #pm %5.1f events",para[0], errpara[0]);
  tleg->AddEntry(f11,text,"f");
  sprintf(text,"BKG %5.1f #pm %5.1f events",para[1], errpara[1]);
  tleg->AddEntry(f12,text,"f");
  tleg->Draw();


  gPad->RedrawAxis();

  char fname[1000];
  sprintf(fname,"plots/unbinned_Ifit_%s.eps",dataInput->GetName());
  c1->SaveAs(fname);
  sprintf(fname,"plots/unbinned_Ifit_%s.gif",dataInput->GetName());
  c1->SaveAs(fname);

  printf("----- fit results with signal projection   ----------- \n");


  fitted[0] = para[0];
  fitted[1] = errpara[0];
  fitted[2] = para[1];
  fitted[3] = errpara[1];

  fitted[4] = para[0]*scale_signal;
  fitted[5] = errpara[0]*scale_signal;
  fitted[6] = para[1]*scale_background;
  fitted[7] = errpara[1]*scale_background;


  return fitted;
}


