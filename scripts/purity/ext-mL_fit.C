#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TF1.h"
#include <fstream>
#include <iostream>
#include <TMinuit.h>
#include "vector.h"
#include <TMath.h>
#include "TVirtualFitter.h"
#include "TFile.h"


using namespace std;

#define NPAR 2

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

//par[0] = fs Jpsi signal fraction;
//-----mass part-----------------------------------------------------
//par[1] = g norm; par[2] g1 mean; par[3] g1 width; 
//par[4] g2 ratio;  par[5] g2 mean; par[6] g2 width;
//par[7] bg norm; par[8] bg slope;

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

  Double_t Lsum=0.;
  Double_t Nevt=0.;
  Double_t fs = par[0];
  Double_t fb = par[1];

  for ( int i=0; i<dataColl.size(); i++ ) {
    Nevt += dataColl[i];
    //PDF for signal and background
    Double_t Ls = sigColl[i];
    Double_t Lb = bkgColl[i];	
    for (int data=0; data<dataColl[i]; data++) {	
      //Get Log Likelihood
      Lsum += TMath::Log( (fs*Ls + fb*Lb) / (fs+fb) );
    }
  }
  f=2*( -1*Lsum + (fs+fb)  - Nevt*TMath::Log(fs+fb) );
  
}

void pulltest(int ptbin=15, char EBEE[10]="EB", float input=0.5){


  TH1F *h1 = new TH1F("h1","",100,-10., 10.);
  TH1F *h2 = new TH1F("h2","",3000, 0., 3000);

  h1->SetNdivisions(505,"XY");
  h2->SetNdivisions(505,"XY");

  int nexp=1000;
  Double_t Nevt=0.;

  for (int i=0; i<nexp; i++) {
    Ifit(ptbin,EBEE);
    Nevt=0.;
    for ( int ii=0; ii<dataColl.size(); ii++ ) {
      Nevt += dataColl[ii];
    }
    printf("fit purity %2.2f +- %2.2f err with %d events. \n", info[0], info_err[0], Nevt);
    h1->Fill((info[0]/Nevt-input)/(info_err[0]/Nevt));
    h2->Fill(info[0]);
  }    

  TCanvas *c2 = new TCanvas("c2","",1000,500);
  c2->Divide(2,1);
  c2->cd(1);
  char txt[100];
  sprintf(txt, "(purity-input)/error");
  h1->SetXTitle(txt);
  h1->Fit("gaus");
  h1->Draw();
  c2->cd(2);
  sprintf(txt, "fitted signal (input %d)", input*Nevt);
  h2->SetXTitle(txt);
  h2->Fit("gaus");
  h2->GetXaxis()->SetRangeUser(0., Nevt*1.2);
  if ( input >0.8 )  h2->GetXaxis()->SetRangeUser(0., Nevt*1.4);
  h2->Draw();  
  sprintf(txt, "plots/extmLfit_pull_%s_pt%d.pdf", EBEE, ptbin);
  c2->SaveAs(txt);

  
}



//___________________________________________________________________________
Double_t* Ifit(TH1F* dataInput, TH1F* sigTemplate, TH1F* bkgTemplate, 
	       int fit_data=1)
{

  TCanvas *c1 = new TCanvas("HF1", "Histos1", 0, 0, 600, 600);
  double count=0;
  dataColl.clear();
  sigColl.clear();
  bkgColl.clear();

  totalColl.clear();
  ctauColl.clear();

  TH1F *hsum = new TH1F();
  float ntemplate = 1.;
  float sigfrac = 0.1;
  TH1F *hsum_norm = new TH1F();
  TH1F *hdata = new TH1F();

  float ndata=0;
  if ( fit_data>0 ) {
    hdata = (TH1F*)dataInput->Clone();
    ndata = hdata->Integral();
  }else {
    hsum = (TH1F*)sigTemplate->Clone();
    hsum->Add(bkgTemplate,1);
    cout << "For histogram " << sigTemplate->GetName() << " and " << bkgTemplate->GetName() << " sum = " << 
      hsum->Integral() << endl;

    if (hsum->Integral()>1.) ntemplate = hsum->Integral();
    sigfrac = sigTemplate->Integral()/ntemplate;

    hsum_norm = (TH1F*)hsum->Clone();  
    hsum_norm->Scale(1./hsum->Integral());

    hdata = (TH1F*)hsum_norm->Clone();
    //ndata = (int) gRandom->Poisson(hsum->Integral());
    ndata=ntemplate;
    hdata->FillRandom(hsum_norm, ndata);
  }
  if(ndata==0) {
    printf(" ---  no events in the fit \n");
    Double_t* fitted = new Double_t[4];
    fitted[0] = 0.;
    fitted[1] = 0.;
    fitted[2] = 0.;
    fitted[3] = 0.;
    return fitted;
  }
    
  printf(" --------- before the fit ------------- \n");
  printf("Nsig %2.3f, Nbg %2.3f, Ntemplate %3.3f \n", sigTemplate->Integral(), bkgTemplate->Integral(), ntemplate);
//   printf("Purity %2.3f, init size %4.3f,  test sample size %4d\n", sigTemplate->Integral()/hsum->Integral(), hsum->Integral(), ndata);
  printf(" -------------------------------------- \n");

  int nbins = hdata->GetNbinsX();

  sigTemplate->Scale(1./sigTemplate->Integral());
  bkgTemplate->Scale(1./bkgTemplate->Integral());  
  for (int ibin=1; ibin<=nbins; ibin++) {
    dataColl.push_back(hdata->GetBinContent(ibin));
    sigColl.push_back(sigTemplate->GetBinContent(ibin));
    bkgColl.push_back(bkgTemplate->GetBinContent(ibin));    
  }
  printf( " -----  Got %d, %d, %d events for fit ----- \n ", dataColl.size(),
	  sigColl.size(), bkgColl.size() );  
  if ( dataColl.size() != sigColl.size() || sigColl.size()!=bkgColl.size() ) {
    printf(" error ...  inconsistent hit collection size \n");
    return;
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
  printf(" Now ready for minimization step \n --------------------------------------------------------- \n");
  
  arglist[0] = 2000; // number of iteration
  arglist[1] = 1.;
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
  printf (" -------------------------------------------- \n");
  printf("Finished.  ierr = %2.2f \n", ierflg);

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

      //do minos if fit sucessed.
//       printf("         ---------------------------------------------------------\n");
//       printf("          Now call for minos step \n");
//       printf("         ---------------------------------------------------------\n");
      
//       arglist[0] = 200; // number of iteration
//       arglist[1] = 1;
//       gMinuit->mnexcm("MINOS", arglist ,2,ierflg);
//       printf("         --------------------------------------------------------- \n");
//       printf("         Done Minos.  ierr = %d \n", ierflg);
//       Double_t amin;
//       gMinuit->mnprin(1,amin);
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

  double yerr[20];
  for(int i=0;i<20;i++){
    yerr[i] = 0.;
  }

  sigTemplate->Scale(para[0]);
  bkgTemplate->Scale(para[1]);
  TH1F *hfit = (TH1F*)sigTemplate->Clone();
  hfit->Add(bkgTemplate);


  sigTemplate->SetLineColor(1);
  sigTemplate->SetFillColor(5);
  sigTemplate->SetFillStyle(3001);

  bkgTemplate->SetLineWidth(2);
  // plot
  c1->Draw();  
  //gPad->SetLogy();
  hdata->SetLineColor(1);
  hdata->SetNdivisions(505,"XY");
  hdata->SetXTitle("comb. ISO (GeV)");
  hdata->SetYTitle("Entries");
  hdata->SetTitleOffset(1.4,"Y");
  hdata->SetTitle();
  hdata->SetMarkerStyle(8);
  hdata->SetMinimum(0.);
  hdata->SetMaximum(hdata->GetMaximum()*1.4);
  hdata->Draw("p e");
  sigTemplate->Draw("hist same");
  bkgTemplate->SetMarkerStyle(0);
  bkgTemplate->SetFillColor(8);
  bkgTemplate->SetLineWidth(1);
  bkgTemplate->SetFillStyle(3013);
  bkgTemplate->SetError(yerr);
  bkgTemplate->Draw("hist same");
  hfit->SetMarkerStyle(0);
  hfit->SetLineColor(1);
  hfit->SetLineWidth(2);
  hfit->SetError(yerr);
  hfit->Draw("hist same");


  TLegend *tleg = new TLegend(0.4, 0.65, 0.95, 0.9);
  char text[50];
  tleg->SetHeader(dataInput->GetTitle());
  tleg->SetFillColor(0);
  tleg->SetShadowColor(0);
  tleg->SetBorderSize(0);
  sprintf(text,"Data %5.1f events",hdata->Integral());
  tleg->AddEntry(hdata,text,"pl");
  sprintf(text,"Fitted %5.1f events",hfit->Integral());
  tleg->AddEntry(hfit,text,"l");
  sprintf(text,"SIG %5.1f #pm %5.1f events",para[0], errpara[0]);
  tleg->AddEntry(sigTemplate,text,"f");
  sprintf(text,"BKG %5.1f #pm %5.1f events",para[1], errpara[1]);
  tleg->AddEntry(bkgTemplate,text,"f");
  tleg->Draw();

  gPad->RedrawAxis();


  cout << hdata->GetName() << endl;
  char fname[300];
  sprintf(fname,"plots/SBfit_Ifit_%s.pdf",hdata->GetName());
  c1->SaveAs(fname);

  printf("----- fit results with signal projection   ----------- \n");

  //   ftemplate->Close();

  Double_t* fitted = new Double_t[4];
  fitted[0] = para[0];
  fitted[1] = errpara[0];
  fitted[2] = para[1];
  fitted[3] = errpara[1];
  return fitted;
}

