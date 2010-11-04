#include <iostream>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TH1.h>
#include <TRandom2.h>
#include <TSystem.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h> 
#include <TCanvas.h>
#include <TLegend.h>
#include "Nov03_ext-mL_fit_ISO.C"

using namespace std;

// the pt and eta binning
const double fBinsPt[]={21.,23.,26.,30.,35.,40.,45.,50.,60.,85.,120.,300,10000};
const double fBinsEta[2] = {0.5,2.0};
const int nPtBin = sizeof(fBinsPt)/sizeof(fBinsPt[0])-1;
const int nEtaBin = 2;

const int REBINNINGS_TEMP=3;
const int REBINNINGS_DATA=5;
const int NEXP = 200;

const float fit_lo = -1.;
const float fit_hi = 20.;

Double_t mySigPDFnorm = 1.0;
Double_t myBkgPDFnorm = 1.0;

Double_t mysum_norm(Double_t *v, Double_t *par)
{
  Double_t func1 = exp_conv(v,par)/mySigPDFnorm;
  Double_t func2 = expinv_power(v,par)/myBkgPDFnorm;
  Double_t func = func1+func2;
  if (func<=0) func=1e-10;
  return func;
}


void toyMC_testFitter(){
	
  gSystem->mkdir("toysPlot");
  char tmp[1000];


  TH1D* htoyResult_pull[nEtaBin][nPtBin];
  TH1D* htoyResult_bias[nEtaBin][nPtBin];

  TFile* finFile = new TFile("template_comb3Iso_template.root");  
  TH1D* hTemplate = (TH1D*)finFile->FindObjectAny("h_EB_comb3Iso_EGdata_pt21");
  hTemplate->Reset();

  
  char* dec[2] = {"EB","EE"};
  for(int ieta=0; ieta<nEtaBin; ieta++){
    for(int ipt=0; ipt < nPtBin; ipt++){

      htoyResult_pull[ieta][ipt] = new TH1D(Form("hpull_%s_pt_%d",		
						 dec[ieta],
						 (int)fBinsPt[ipt]),
 					    "",50,-5.0,5.0);


      htoyResult_bias[ieta][ipt] = new TH1D(Form("hbias_%s_pt_%d",		
						 dec[ieta],
						 (int)fBinsPt[ipt]),
 					    "",100,-0.5,0.5);


    }
  }



  TH1D* hfit_sig;
  TH1D* hfit_bkg;

  TH1D* hTemplate_S[nEtaBin][nPtBin];
  TH1D* hTemplate_B[nEtaBin][nPtBin];
  TH1D* hdata_data[nEtaBin][nPtBin];
  TH1D* htemp;


  for(int ieta=0; ieta< nEtaBin; ieta++){
    for(int ipt=0; ipt < nPtBin; ipt++){

      if(ieta!=0 || (ipt!=7 && ipt!=4))continue;
//       if(ieta!=0 || (ipt!=4))continue;
      
       // getting histograms from data root file
       sprintf(tmp,"h_%s_comb3Iso_EGdata_pt%d",dec[ieta],(int)fBinsPt[ipt]);
       cout << "looking for histogram " << tmp << " in file " << 
	 finFile->GetName() << endl;
       hdata_data[ieta][ipt] = (TH1D*)finFile->FindObjectAny(tmp);
       hdata_data[ieta][ipt]->Rebin(REBINNINGS_DATA);

       // filling unbinned data
       htemp = (TH1D*)hdata_data[ieta][ipt]->Clone("htemp");
       htemp->Reset();


       sprintf(tmp,"h_%s_comb3Iso_sig_pt%d",dec[ieta],(int)fBinsPt[ipt]);
       cout << "looking for histogram " << tmp << " in file " << 
	 finFile->GetName() << endl;
       hTemplate_S[ieta][ipt] = (TH1D*)finFile->FindObjectAny(tmp);
       hTemplate_S[ieta][ipt]->Rebin(REBINNINGS_TEMP);

       if(ieta==0)
       sprintf(tmp,"h_%s_comb3Iso_bkg_pt%d",dec[ieta],(int)fBinsPt[ipt]);
       else if(ieta==1)
	 sprintf(tmp,"h_%s_comb3IsoSB_EGdata_pt%d",dec[ieta],(int)fBinsPt[ipt]);
       cout << "looking for histogram " << tmp << " in file " << 
	 finFile->GetName() << endl;
       hTemplate_B[ieta][ipt] = (TH1D*)finFile->FindObjectAny(tmp);
       hTemplate_B[ieta][ipt]->Rebin(REBINNINGS_TEMP);




       const int NRETURN = 3*NPAR;
       Double_t myFitPar[NRETURN]={0};
       Double_t* FuncFitResult;
       FuncFitResult = Ifit("EGdata_comb3Iso_et.dat",
			    hTemplate_S[ieta][ipt],hTemplate_B[ieta][ipt],
			    hdata_data[ieta][ipt], myFitPar,
			    (int)fBinsPt[ipt], dec[ieta],2);

       Double_t nsig_input    = FuncFitResult[0];
       Double_t nsigerr_input = FuncFitResult[1];
       Double_t nbkg_input    = FuncFitResult[2];
       Double_t nbkgerr_input = FuncFitResult[3];

	    
       Double_t sigFitPar[NPAR]={0};
       for(int ipar=0; ipar<NPAR; ipar++)
	 sigFitPar[ipar] = myFitPar[ipar];

       Double_t bkgFitPar[NPAR]={0};
       for(int ipar=0; ipar<NPAR; ipar++)
	 bkgFitPar[ipar] = myFitPar[ipar+NPAR];


       Double_t sumFitPar[NPAR]={0};
       for(int ipar=0; ipar<NPAR; ipar++)
	 sumFitPar[ipar] = myFitPar[ipar+NPAR*2];


       TF1* fsig = new TF1("fsig", exp_conv, fit_lo, fit_hi, 11);
       fsig->SetParameters(sigFitPar);       
       fsig->SetParameter(0,1.0);

       mySigPDFnorm = fsig->Integral(fit_lo,fit_hi);
       cout << "mySigPDFnorm = " << mySigPDFnorm << endl;
      
       TF1* fbkg = new TF1("fbkg", expinv_power, fit_lo, fit_hi, 11);
       fbkg->SetParameters(bkgFitPar);
       fbkg->SetParameter(4,1.0);

       myBkgPDFnorm = fbkg->Integral(fit_lo, fit_hi);
       cout << "myBkgPDFnorm = " << myBkgPDFnorm << endl;
 
       TF1* fsum = new TF1("fsum",mysum_norm, fit_lo, fit_hi,11);
       fsum->SetParameters(sumFitPar);

       cout << "Using nsig_input = " << nsig_input << endl;
       cout << "Using nbkg_input = " << nbkg_input << endl;

       fsum->SetParameter(0, nsig_input*hdata_data[ieta][ipt]->GetBinWidth(2));
       fsum->SetParameter(4, nbkg_input*hdata_data[ieta][ipt]->GetBinWidth(2));
       fsum->SetLineColor(2);
       fsum->SetNpx(2500);

       cout << "fsum integral = " << fsum->Integral(fit_lo,fit_hi) << endl;
       for(int ipar=0; ipar<NPAR; ipar++)cout << "fsum par " << ipar << " = " << fsum->GetParameter(ipar) << endl;


//        FILE *infile =  fopen("EGdata_comb3Iso_et.dat","r");  
//        float xdata, xdata1, xdata2;
//        for (int i=0;i<datapoint;i++){
// 	 fscanf(infile,"%f %f %f",&xdata, &xdata1, &xdata2);
// 	 if( xdata1 >= fBinsPt[ipt] && xdata1 < fBinsPt[ipt+1] && xdata<20.) {
// 	   if(TMath::Abs(xdata2)<1.45) {
// 	     htemp->Fill(xdata);
// 	   }
// 	 }
//        }


       
//        TCanvas* myCanvas = new TCanvas("myCanvas","myCanvas");
//        htemp->Draw();
//        fsum->Draw("same");

       
//        // loops over toys
       for(int iexp=0; iexp<NEXP; iexp++){

	 TH1D* htoyMC_data = (TH1D*)hdata_data[ieta][ipt]->Clone("htoyMC_data");
	 htoyMC_data->Reset();
	 
	 TH1D* htoyMC_sig  = (TH1D*)hTemplate_S[ieta][ipt]->Clone("htoyMC_sig");
// 	 htoyMC_sig->Reset();
	 
	 TH1D* htoyMC_bkg  = (TH1D*)hTemplate_B[ieta][ipt]->Clone("htoyMC_bkg");
// 	 htoyMC_bkg->Reset();




	 UInt_t nowSeed = (unsigned long)gSystem->Now();
	 gRandom->SetSeed(nowSeed);
	 int nsiggen  = gRandom->Poisson(nsig_input);
	 int nbkggen  = gRandom->Poisson(nbkg_input);
	 int ndata = nsiggen + nbkggen;

	 // reset toy MC data
	 htoyMC_data->Reset();
	 ofstream fout;
	 std::string toyData = "toy.dat";
	 fout.open(toyData.data());
	 for(int ieve=0; ieve < nsiggen; ieve++)
	   {
	     Double_t xvalue = fsig->GetRandom(fit_lo,fit_hi);
	     fout << xvalue << " " << 
	       0.5*(fBinsPt[ipt]+fBinsPt[ipt+1]) << " " << fBinsEta[ieta] << endl;
	     htoyMC_data->Fill(xvalue);
	   }
	     
	 for(int ieve=0; ieve < nbkggen; ieve++)
	   {
	     Double_t xvalue = fbkg->GetRandom(fit_lo,fit_hi);
	     fout << xvalue << " " << 
	       0.5*(fBinsPt[ipt]+fBinsPt[ipt+1]) << " " << fBinsEta[ieta] << endl;
	     htoyMC_data->Fill(xvalue);
	   }
	     
	 fout.close();

// 	 htoyMC_data->FillRandom("fsig",nsiggen);
// 	 htoyMC_data->FillRandom("fbkg",nbkggen);


// 	 htoyMC_data->Draw();
	 
	 cout << "Generated ndata = " << ndata << "\t nsiggen = " << nsiggen << " \t nbkggen = " << 
	   nbkggen << endl;


	 Double_t* toyFitResult;
	 Double_t toyMyFitPar[NRETURN]={0};

	 toyFitResult =  Ifit(toyData.data(), 
			      htoyMC_sig, htoyMC_bkg,
			      htoyMC_data, toyMyFitPar,
			      (int)fBinsPt[ipt], dec[ieta],2);

  	Double_t nsigtoyfit    = toyFitResult[0];
  	Double_t errnsigtoyfit = toyFitResult[1];
  	Double_t nbkgtoyfit    = toyFitResult[2];
  	Double_t errnbkgtoyfit = toyFitResult[3];

	 Double_t toySumFitPar[NPAR]={0};
	 for(int ipar=0; ipar<NPAR; ipar++)
	 toySumFitPar[ipar] = toyMyFitPar[ipar+NPAR*2];


 	 fsum->SetParameters(toySumFitPar);
   	 fsum->SetParameter(0, nsigtoyfit*hdata_data[ieta][ipt]->GetBinWidth(2));
  	 fsum->SetParameter(4, nbkgtoyfit*hdata_data[ieta][ipt]->GetBinWidth(2));
//    	 fsum->SetParameter(0, (float)nsiggen*hdata_data[ieta][ipt]->GetBinWidth(2));
//   	 fsum->SetParameter(4, (float)nbkggen*hdata_data[ieta][ipt]->GetBinWidth(2));

	 if(iexp%4==0){
 	 TCanvas* myCanvas = new TCanvas("myCanvas","SHIT");
	 htoyMC_data->SetMaximum(htoyMC_data->GetMaximum()*1.5);
 	 htoyMC_data->Draw();
 	 fsum->Draw("same");
	 myCanvas->Print(Form("toysPlot/fit_%03i.gif",iexp));
	 delete myCanvas;
	 }

	 cout << "fsum integral = " << fsum->Integral(-1,20) << endl;
	 for(int ipar=0; ipar<NPAR; ipar++)cout << "fsum par " << ipar << " = " << fsum->GetParameter(ipar) << endl;


	
	cout << "Total data = " << nsiggen + nbkggen << endl;
	cout << "toyMC_data->Integral() = " << htoyMC_data->Integral() << endl;
	cout << "Expected nsig = " << nsiggen << " and Fitted nsig = " << nsigtoyfit << endl;
	cout << "Expected nbkg = " << nbkggen << " and Fitted nbkg = " << nbkgtoyfit << endl;
	cout << "Input nsig = " << nsig_input << endl;

    	Double_t pull = (nsigtoyfit - nsig_input)/errnsigtoyfit;
    
   	htoyResult_pull[ieta][ipt]->Fill(pull);

	Double_t bias = (nsigtoyfit - nsig_input)/nsig_input;
	htoyResult_bias[ieta][ipt]->Fill(bias);


       } // end loops of toys

       TCanvas* myToyCanvas = new TCanvas("myToyCanvas","",1000,500);
       myToyCanvas->Divide(2,1);
       myToyCanvas->cd(1);

       htoyResult_pull[ieta][ipt]->SetFillColor(kAzure-4);
       htoyResult_pull[ieta][ipt]->SetFillStyle(1001);
       htoyResult_pull[ieta][ipt]->Draw();

       myToyCanvas->cd(2);

       htoyResult_bias[ieta][ipt]->SetFillColor(kViolet-9);
       htoyResult_bias[ieta][ipt]->SetFillStyle(1001);
       htoyResult_bias[ieta][ipt]->Draw();
    
       delete fsig;
       delete fbkg;


    } // end of loop over pt bins


  } 


  TFile* outFile = new TFile("fittertest.root",
			     "recreate");

  for(int ieta=0; ieta < nEtaBin; ieta++){
    for(int ipt=0; ipt < nPtBin; ipt++){      
 
      if(htoyResult_pull[ieta][ipt]->GetEntries()>0)
       htoyResult_pull[ieta][ipt]->Write();

      if(htoyResult_bias[ieta][ipt]->GetEntries()>0)
	htoyResult_bias[ieta][ipt]->Write();      
    }
  }
   
  outFile->Close();

}
