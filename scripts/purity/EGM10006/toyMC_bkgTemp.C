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
#include <TROOT.h>
#include "Nov09_fit.C"

using namespace std;

// the pt and eta binning
const Double_t fBinsPt[]={21.,23.,26.,30.,35.,40.,45.,50.,60.,85.,120.,300};
const Double_t fBinsEta[2] = {0.5,2.0};

const int nPtBin = sizeof(fBinsPt)/sizeof(fBinsPt[0])-1;
const int nEtaBin = 2;

const Double_t fNsigPt[2][nPtBin]=
  {
    {
      26304.067089,
      22608.144062,
      21972.048313,
      14271.476950,
      7053.907180,
      4792.810030,
      2845.355508,
      2906.144032,
      2391.605544,
      661.728019,
      229.530958
    },
    {
      18264.588466,
      16262.348287,
      16284.902091,
      10317.982854,
      5798.697786,
      3171.995241,
      1718.942019,
      1840.277296,
      1366.803312,
      333.989036,
      106.529129
    }
  };


const int REBINNINGS_TEMP=3;
const int REBINNINGS_DATA=5;
const int NEXP = 200;

const Double_t fit_lo = -1.;
const Double_t fit_hi = 20.;

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


void toyMC_bkgTemp(int runIeta=-1, int runIpt=-1){
	
  cout << "ROOT version = " << gROOT->GetVersion() << endl;
  //   gSystem->mkdir("toysPlot");
  char tmp[1000];


  TH1D* htoyResult_pull[nEtaBin][nPtBin];
  TH1D* htoyResult_bias[nEtaBin][nPtBin];

  TFile *fsumFile = new TFile("/afs/cern.ch/user/s/syu/scratch0/LxplusArea/proj_comb_comb3Iso_template.root");
  TFile* finFile = new TFile("/afs/cern.ch/user/s/syu/scratch0/LxplusArea/template_comb3Iso_template.root");  
  TFile* zeeFile = new TFile("/afs/cern.ch/user/s/syu/scratch0/LxplusArea/anadipho_Zee_Vg_3pb.root");
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
  TH1D* hZeeTemplate_S[nEtaBin];
  TH1D* hdata_data[nEtaBin][nPtBin];
  TH1D* htemp;
  TH1D* hTemplate_ShiftB[nEtaBin][nPtBin];

  for(int ieta=0; ieta< nEtaBin; ieta++){

    // getting a different signal template
    if(ieta==0){
      sprintf(tmp,"h_%s_combIso",dec[ieta]);
      cout << "looking for histogram " << tmp << " in file " << 
	zeeFile->GetName() << endl;
      hZeeTemplate_S[ieta]= (TH1D*)zeeFile->FindObjectAny(tmp);	
      hZeeTemplate_S[ieta]->Rebin(REBINNINGS_TEMP);
    }
    else {
      sprintf(tmp,"h_%s_comb3Iso_sig_sum_SIG",dec[ieta]); //no pt dep.
      cout << "looking for histogram " << tmp << " in file " << 
	fsumFile->GetName() << endl;
      hZeeTemplate_S[ieta]= (TH1D*)fsumFile->FindObjectAny(tmp);	
      hZeeTemplate_S[ieta]->Rebin(REBINNINGS_TEMP);
    }

    for(int ipt=0; ipt < nPtBin; ipt++){

      if(runIeta>=0 && ieta!=runIeta)continue;
      if(runIpt>=0 && ipt!=runIpt)continue;
      
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

      if(ieta==0 && fBinsPt[ipt]>=50)
	sprintf(tmp,"h_%s_comb3Iso_bkg_pt%d",dec[ieta],50);
      else if(ieta==0 )
	sprintf(tmp,"h_%s_comb3Iso_bkg_pt%d",dec[ieta],(int)fBinsPt[ipt]);  
      else if(ieta==1 && fBinsPt[ipt]>=60)
	sprintf(tmp,"h_%s_comb3IsoSB_EGdata_pt%d",dec[ieta],60);       
      else if(ieta==1)
	sprintf(tmp,"h_%s_comb3IsoSB_EGdata_pt%d",dec[ieta],(int)fBinsPt[ipt]);
      cout << "looking for histogram " << tmp << " in file " << 
	finFile->GetName() << endl;
      hTemplate_B[ieta][ipt] = (TH1D*)finFile->FindObjectAny(tmp);
      hTemplate_B[ieta][ipt]->Rebin(REBINNINGS_TEMP);

      // getting a different background template
      if(ieta==0 && fBinsPt[ipt]>=50)
	sprintf(tmp,"h_%s_comb3IsoSB_EGdata_pt%d",dec[ieta],50);
      else if(ieta==0 )
	sprintf(tmp,"h_%s_comb3IsoSB_EGdata_pt%d",dec[ieta],(int)fBinsPt[ipt]);  
      else if(ieta==1 && fBinsPt[ipt]>=60)
	sprintf(tmp,"h_%s_comb3Iso_bkg_pt%d",dec[ieta],60);       
      else if(ieta==1)
	sprintf(tmp,"h_%s_comb3Iso_bkg_pt%d",dec[ieta],(int)fBinsPt[ipt]);
      cout << "looking for histogram " << tmp << " in file " << 
	finFile->GetName() << endl;
      hTemplate_ShiftB[ieta][ipt] = (TH1D*)finFile->FindObjectAny(tmp);

      hTemplate_ShiftB[ieta][ipt]->Rebin(REBINNINGS_TEMP);



      const int NRETURN = 3*NPAR;
      Double_t myFitPar[NRETURN]={0};
      Double_t* FuncFitResult;
      FuncFitResult = Ifit("/afs/cern.ch/user/s/syu/scratch0/LxplusArea/EGdata_comb3Iso_et.dat",
			   hZeeTemplate_S[ieta],hTemplate_ShiftB[ieta][ipt],
			   hdata_data[ieta][ipt], myFitPar,
			   (int)fBinsPt[ipt], dec[ieta],2);

      Double_t nsig_input    = FuncFitResult[0];
      Double_t nsigerr_input = FuncFitResult[1];
      Double_t nbkg_input    = FuncFitResult[2];
      Double_t nbkgerr_input = FuncFitResult[3];
      Double_t nsig_input5GeV = FuncFitResult[4];
      Double_t nsigerr_input5GeV = FuncFitResult[5];

       // force the parameters since EE pt=40 fails
       if(ieta==1 && ipt==5)
	 {
	   nsig_input = 3172.0;
	   nbkg_input = 10031.0;
	   nsig_input5GeV = 3158.7;
	   
	   Double_t tempPar[NRETURN]={
	     22517.049862,
	     0.900766,
	     0.044772,
	     0.191920,
	     313.878244,
	     -0.545069,
	     -0.281830,
	     0.026143,
	     2.549494,
	     0.,
	     0.,
	     22517.049862,
	     0.900766,
	     0.044772,
	     0.191920,
	     313.878244,
	     -0.545069,
	     -0.281830,
	     0.026143,
	     2.549494,
	     0.,
	     0.,
	     22517.049862,
	     0.900766,
	     0.044772,
	     0.191920,
	     313.878244,
	     -0.545069,
	     -0.281830,
	     0.026143,
	     2.549494,
	     0.,
	     0.
	   };

	   for(int ipar=0; ipar <NRETURN; ipar++)
	     myFitPar[ipar] = tempPar[ipar];
	 } // end if EE, Et=40 GeV


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
      cout << "Current parameter = " << fbkg->GetParameter(5) << endl;

      myBkgPDFnorm = fbkg->Integral(fit_lo, fit_hi);
      cout << "myBkgPDFnorm = " << myBkgPDFnorm << endl;
 
      TF1* fsum = new TF1("fsum",mysum_norm, fit_lo, fit_hi,11);
      fsum->SetParameters(sumFitPar);


      nsig_input = fNsigPt[ieta][ipt];

      cout << "Using nsig_input = " << nsig_input << endl;
      cout << "Using nbkg_input = " << nbkg_input << endl;
      cout << "Using nsig_input5GeV = " << nsig_input5GeV << endl;

      fsum->SetParameter(0, nsig_input*hdata_data[ieta][ipt]->GetBinWidth(2));
      fsum->SetParameter(4, nbkg_input*hdata_data[ieta][ipt]->GetBinWidth(2));
      fsum->SetLineColor(2);
      fsum->SetNpx(2500);

      //        cout << "fsum integral = " << fsum->Integral(fit_lo,fit_hi) << endl;
      //        for(int ipar=0; ipar<NPAR; ipar++)cout << "fsum par " << ipar << " = " << fsum->GetParameter(ipar) << endl;

//       FILE *infile =  fopen("/afs/cern.ch/user/s/syu/scratch0/LxplusArea/EGdata_comb3Iso_et.dat","r");  
//       Double_t xdata, xdata1, xdata2; // combined isolation, pt, eta
//       int flag = 1;
//       while (flag!=-1){
// 	flag =fscanf(infile,"%f %f %f",&xdata, &xdata1, &xdata2);
// 	if( xdata1 >= fBinsPt[ipt] && xdata1 < fBinsPt[ipt+1] && xdata<20.) {
// 	  if((ieta==0 && TMath::Abs(xdata2)<1.5) ||
// 	     (ieta==1 && TMath::Abs(xdata2)>1.5) ) {
// 	    htemp->Fill(xdata);
// 	  }
// 	} 
//       }// keep reading files as long as text exists

    
//       TCanvas* myCanvas = new TCanvas("myCanvas","myCanvas");
//       htemp->Draw();
//       fsum->Draw("same");

       
      //        // loops over toys
      for(int iexp=0; iexp<NEXP; iexp++){

	TH1D* htoyMC_data = (TH1D*)hdata_data[ieta][ipt]->Clone("htoyMC_data");
	htoyMC_data->Reset();
	 
	TH1D* htoyMC_sig  = (TH1D*)hZeeTemplate_S[ieta]->Clone("htoyMC_sig");
	 
	TH1D* htoyMC_bkg  = (TH1D*)hTemplate_B[ieta][ipt]->Clone("htoyMC_bkg");




	UInt_t nowSeed = (unsigned long)gSystem->Now();
	gRandom->SetSeed(nowSeed);
	int nsiggen  = gRandom->Poisson(nsig_input);
	int nbkggen  = gRandom->Poisson(nbkg_input);
	int ndata = nsiggen + nbkggen;

	// reset toy MC data
	htoyMC_data->Reset();
	ofstream fout;
	std::string toyData = Form("/tmp/syu/changeBkgTemptoy_%d_%d.dat",ieta,ipt);

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
	Double_t nsigtoyfit_5GeV = toyFitResult[4];
	Double_t errnsigtoyfit_5GeV = toyFitResult[5];


	if(errnsigtoyfit < 1e-6 || errnbkgtoyfit < 1e-6 || 
	   nsigtoyfit < 1e-6 || nbkgtoyfit < 1e-6 || 
	   nsig_input < 1e-6 || fabs(nsigtoyfit - nsig_input)<1e-4)continue;
	   

	Double_t toySumFitPar[NPAR]={0};
	for(int ipar=0; ipar<NPAR; ipar++)
	  toySumFitPar[ipar] = toyMyFitPar[ipar+NPAR*2];


//  	fsum->SetParameters(toySumFitPar);
//  	fsum->SetParameter(0, nsigtoyfit*hdata_data[ieta][ipt]->GetBinWidth(2));
//  	fsum->SetParameter(4, nbkgtoyfit*hdata_data[ieta][ipt]->GetBinWidth(2));
	fsum->SetParameter(0, (Double_t)nsiggen*hdata_data[ieta][ipt]->GetBinWidth(2));
	fsum->SetParameter(4, (Double_t)nbkggen*hdata_data[ieta][ipt]->GetBinWidth(2));

	// 	if(iexp%20==0){
 	 	  TCanvas* myCanvas = new TCanvas("myCanvas","SHIT");
 	 	  htoyMC_data->SetMaximum(htoyMC_data->GetMaximum()*1.5);
 	 	  htoyMC_data->Draw();
 	 	  fsum->Draw("same");
// 	 	  myCanvas->Print(Form("toysPlot/fit_%03i.gif",iexp));
// 	 	  delete myCanvas;
	// 	}

	// 	 cout << "fsum integral = " << fsum->Integral(-1,20) << endl;
	// 	 for(int ipar=0; ipar<NPAR; ipar++)cout << "fsum par " << ipar << " = " << fsum->GetParameter(ipar) << endl;

	
	Double_t pull = (nsigtoyfit - nsig_input)/errnsigtoyfit;
	Double_t bias = (nsigtoyfit - nsig_input)/nsig_input;

    
   	htoyResult_pull[ieta][ipt]->Fill(pull);
	htoyResult_bias[ieta][ipt]->Fill(bias);


      } // end loops of toys

      delete fsig;
      delete fbkg;


    } // end of loop over pt bins


  }


  TFile* outFile = new TFile(Form("changBkgTemp_%s_pt%d.root",
				  dec[runIeta], (int)fBinsPt[runIpt]),
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
