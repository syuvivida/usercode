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
#include "ext-mL_fit.C"
#include "ext-mL_fit_unbinned.C"
#include "setTDRStyle.C"

using namespace std;

// the pt and eta binning
const Double_t fBinsEta[]={0,1.45,1.7,2.5};
const Double_t fBinsPt[]={15,20,30,50,80,120};
const int nPtBin = sizeof(fBinsPt)/sizeof(fBinsPt[0])-1;
const int nEtaBin = (sizeof(fBinsEta)/sizeof(fBinsEta[0]))/2;

const int REBINNINGS=1;

const int NEXP = 1000;

Double_t nsig_data[nEtaBin][nPtBin]=
  {
     {
       1050.8,
       317.8,
       92.2,
       18.7,
       3.0
     },
     {
       930.7,
       258.7,
       60.1,
       15.4,
       0.4
     }

  };


void ratioErr(Double_t n1, Double_t n1err, Double_t n2, Double_t n2err,
	      Double_t& ratio, Double_t& err)
{
  if(fabs(n1+n2)<1e-6){ratio=0;err=0;return;}
  ratio = (n1)/(n1+n2);
  
  err= pow(1/(n1+n2) - n1/(n1+n2)/(n1+n2),2)*n1err*n1err+
    pow(n1/(n1+n2)/(n1+n2),2)*n2err*n2err; 

  err = sqrt(err);

}



void toyMC_background(){
	

  char tmp[1000];
  setTDRStyle();
  // settings for purity TGraphAsymmetryErrors
  Double_t fBinsPtMidPoint[nPtBin]={0};
  Double_t fBinsPtError[nPtBin]={0};
  
  for(int ipt=0; ipt < nPtBin; ipt++)
    {
      fBinsPtMidPoint[ipt] = 0.5*(fBinsPt[ipt+1]+fBinsPt[ipt]);
      fBinsPtError[ipt] = 0.5*(fBinsPt[ipt+1]-fBinsPt[ipt]);
    }


  Double_t nsig_func[nEtaBin][nPtBin]={{0}};
  Double_t nbkg_func[nEtaBin][nPtBin]={{0}};
  Double_t nsig_err_func[nEtaBin][nPtBin]={{0}};
  Double_t nbkg_err_func[nEtaBin][nPtBin]={{0}};


  TH1F* hTemplate_S[nEtaBin][nPtBin];
  TH1F* hTemplate_B[nEtaBin][nPtBin];
  TH1F* hTemplate_data[nEtaBin][nPtBin];
  TH1F* hdata_S[nEtaBin][nPtBin];
  TH1F* hdata_B[nEtaBin][nPtBin];
  TH1F* hdata_data[nEtaBin][nPtBin];


  
  std::string dataFile     = "rawYield_SBMC.root";
  std::string templateFile = "rawYield_SBdata.root";

  TFile* inf_data = new TFile(dataFile.data());
  TFile* inf_template = new TFile(templateFile.data());

  TH1F* htoyResult_pull[nEtaBin][nPtBin];
  TH1F* htoyResult_diff[nEtaBin][nPtBin];

  
  char* dec[2] = {"EB","EE"};
  for(int ieta=0; ieta<nEtaBin; ieta++){
    for(int ipt=0; ipt < nPtBin; ipt++){

      htoyResult_pull[ieta][ipt] = new TH1F(Form("hpull_Eta_%.2f_%.2f_Et_%d_%d",		
						 fBinsEta[ieta*2],fBinsEta[ieta*2+1],
						 (int)fBinsPt[ipt], (int)fBinsPt[ipt+1]),
 					    "",200,-20.0,20.0);


      htoyResult_diff[ieta][ipt] = new TH1F(Form("hdiff_Eta_%.2f_%.2f_Et_%d_%d",		
						 fBinsEta[ieta*2],fBinsEta[ieta*2+1],
						 (int)fBinsPt[ipt], (int)fBinsPt[ipt+1]),
 					    "",50,-1.0,1.0);


      // getting histograms from data root file
      sprintf(tmp,"hOutputData_%s_pt%d",dec[ieta],(int)fBinsPt[ipt]);
      cout << "looking for histogram " << tmp << " in file " << 
	inf_data->GetName() << endl;
      hdata_data[ieta][ipt] = (TH1F*)inf_data->FindObjectAny(tmp);
      hdata_data[ieta][ipt]->Rebin(REBINNINGS);
      // setting titles
      sprintf(tmp,"%.2f < |#eta(#gamma)| < %.2f, %d < p_{T}(#gamma) < %d GeV",
	      fBinsEta[ieta*2],fBinsEta[ieta*2+1],
	      (int)fBinsPt[ipt], (int)fBinsPt[ipt+1]);
      hdata_data[ieta][ipt]->SetTitle(tmp);


      sprintf(tmp,"hOutputSig_%s_pt%d",dec[ieta],(int)fBinsPt[ipt]);
      cout << "looking for histogram " << tmp << " in file " << 
	inf_data->GetName() << endl;
      hdata_S[ieta][ipt] = (TH1F*)inf_data->FindObjectAny(tmp);
      hdata_S[ieta][ipt]->Rebin(REBINNINGS);

      sprintf(tmp,"hOutputBkg_%s_pt%d",dec[ieta],(int)fBinsPt[ipt]);
      cout << "looking for histogram " << tmp << " in file " << 
	inf_data->GetName() << endl;
      hdata_B[ieta][ipt] = (TH1F*)inf_data->FindObjectAny(tmp);
      hdata_B[ieta][ipt]->Rebin(REBINNINGS);


      // getting histogram for template root file
      sprintf(tmp,"hOutputData_%s_pt%d",dec[ieta],(int)fBinsPt[ipt]);
      cout << "looking for histogram " << tmp << " in file " << 
	inf_template->GetName() << endl;
      hTemplate_data[ieta][ipt] = (TH1F*)inf_template->FindObjectAny(tmp);
      hTemplate_data[ieta][ipt]->Rebin(REBINNINGS);
      sprintf(tmp,"%.2f < |#eta(#gamma)| < %.2f, %d < p_{T}(#gamma) < %d GeV",
	      fBinsEta[ieta*2],fBinsEta[ieta*2+1],
	      (int)fBinsPt[ipt], (int)fBinsPt[ipt+1]);
      hTemplate_data[ieta][ipt]->SetTitle(tmp);

      sprintf(tmp,"hOutputSig_%s_pt%d",dec[ieta],(int)fBinsPt[ipt]);
      cout << "looking for histogram " << tmp << " in file " << 
	inf_template->GetName() << endl;
      hTemplate_S[ieta][ipt] = (TH1F*)inf_template->FindObjectAny(tmp);
      hTemplate_S[ieta][ipt]->Rebin(REBINNINGS);

      sprintf(tmp,"hOutputBkg_%s_pt%d",dec[ieta],(int)fBinsPt[ipt]);
      cout << "looking for histogram " << tmp << " in file " << 
	inf_template->GetName() << endl;
      hTemplate_B[ieta][ipt] = (TH1F*)inf_template->FindObjectAny(tmp);
      hTemplate_B[ieta][ipt]->Rebin(REBINNINGS);

    }
  }


  TH1F* htoyMC_data = (TH1F*)hdata_data[0][0]->Clone();
//   htoyMC_data->SetName("htoyMC_data");
  htoyMC_data->Reset();
  
  TH1F* htoyMC_sig  = (TH1F*)hTemplate_S[0][0]->Clone();
//   htoyMC_sig->SetName("htoyMC_sig");
  htoyMC_sig->Reset();

  TH1F* htoyMC_bkg  = (TH1F*)hTemplate_B[0][0]->Clone();
//   htoyMC_bkg->SetName("htoyMC_bkg");
  htoyMC_bkg->Reset();

  TH1F* hfit_sig;
  TH1F* hfit_bkg;


  for(int ieta=0; ieta< nEtaBin; ieta++){
    for(int ipt=0; ipt < nPtBin; ipt++){
      
      // first obtain fit for the template
      hfit_sig  = (TH1F*)hdata_S[ieta][ipt]->Clone();
      hfit_sig  -> SetName("hfit_sig");
      if(ieta==0)
	hfit_sig  -> Rebin(2);      
//       hfit_sig->Scale(1.0/(float)hfit_sig->Integral());

      hfit_bkg  = (TH1F*)hdata_B[ieta][ipt]->Clone();
      hfit_bkg  -> SetName("hfit_bkg");
      if(ieta==0)
	hfit_bkg  -> Rebin(2);
//       hfit_bkg->Scale(1.0/(float)hfit_bkg->Integral());

      double sigPar[20] = {hfit_sig->GetMaximum(), 1., 0.5, 0.3,
			hfit_bkg->GetMaximum(),-.3,-1., 0.01, 0.5, 0.01, 1., 1.};

      TF1 *fsig = new TF1("fsig", exp_conv, -1., 11., 12);
      fsig->SetParameters(sigPar);
      
      hfit_sig->Fit(fsig,"","",-1,5.0);
      
      TF1 *fbkg = new TF1("fbkg", expinv_power, -1., 11., 12);

      fbkg->SetParameters(fsig->GetParameters());
      fbkg->SetParLimits(5,-10.,1.);
      fbkg->SetParLimits(6,-1.,2.);
      fbkg->SetParLimits(7,0.,0.09);

      if(ieta==1 && ipt==3)
	{
// 	  cout << "find EE in pt bin 50--80" << endl;
	  fbkg->FixParameter(5,-0.1);
	}
      if(ieta==0 && ipt==4)
	{
// 	  cout << "find EB in pt bin 80--120" << endl;
	  fbkg->FixParameter(5,-0.1);
	}

      fbkg->FixParameter(8,0.5);
      fbkg->FixParameter(0,fbkg->GetParameter(0)); 
      fbkg->FixParameter(1,fbkg->GetParameter(1));
      fbkg->FixParameter(2,fbkg->GetParameter(2));
      fbkg->FixParameter(3,fbkg->GetParameter(3));

      hfit_bkg->Fit(fbkg,"b");

//       for(int i=0;i<12;i++)cout << "parameter " << i << " = " << 
// 	fbkg->GetParameter(i) << endl;
 
       htoyMC_sig  = (TH1F*)hTemplate_S[ieta][ipt]->Clone();
       htoyMC_bkg  = (TH1F*)hTemplate_B[ieta][ipt]->Clone();


       // loops over toys
       for(int iexp=0; iexp<NEXP; iexp++){
 	htoyMC_data->Reset();

 	UInt_t nowSeed = (unsigned long)gSystem->Now();
 	gRandom->SetSeed(nowSeed);
  	int ndata = hdata_data[ieta][ipt]->Integral();
 	int nsig  = gRandom->Poisson(nsig_data[ieta][ipt]);
 	int nbkg  = gRandom->Poisson(ndata - nsig_data[ieta][ipt]);

 // 	htoyMC_data->FillRandom(hTemplate_S[ieta][ipt],nsig);
 // 	htoyMC_data->FillRandom(hTemplate_B[ieta][ipt],nbkg);
 	htoyMC_data->FillRandom("fsig",nsig);
 	htoyMC_data->FillRandom("fbkg",nbkg);


	htoyMC_data->Draw();

 	cout << "ndata = " << ndata << "\t nsig = " << nsig << " \t nbkg = " << 
 	  nbkg << endl;

 	Double_t scaleNumber[20];
 	for(int i=0;i<20;i++)scaleNumber[i]=1.;


  	// 2nd, get unbinned fit
  	Double_t* FuncFitResult;

  	FuncFitResult = Ifit(htoyMC_data,
  			     htoyMC_sig,
  			     htoyMC_bkg,0,"EGdata_100628.dat",
  			     fBinsEta[ieta*2],fBinsEta[ieta*2+1],
  			     fBinsPt[ipt],fBinsPt[ipt+1],
  			     NULL
  			     );	



  	Double_t nsigfunc    = FuncFitResult[0];
  	Double_t errnsigfunc = FuncFitResult[1];
  	Double_t nbkgfunc    = FuncFitResult[2];
  	Double_t errnbkgfunc = FuncFitResult[3];


   	Double_t pull = (nsigfunc - nsig_data[ieta][ipt])/errnsigfunc;
   	Double_t diff = (nsigfunc - nsig_data[ieta][ipt])/nsig_data[ieta][ipt];
    
//  //  	Double_t* TemplateFitResult;
//  //  	TemplateFitResult = IfitBin(htoyMC_data,
//  //  				    htoyMC_sig,
//  //  				    htoyMC_bkg);


//  //  	Double_t nsigtemplate    = TemplateFitResult[0];
//  //  	Double_t errnsigtemplate = TemplateFitResult[1];
//  //  	Double_t nbkgtemplate    = TemplateFitResult[2];
//  //  	Double_t errnbkgtemplate = TemplateFitResult[3];

//  //   	Double_t pull = (nsigtemplate - nsig_data[ieta][ipt])/errnsigtemplate;
    
	htoyResult_diff[ieta][ipt]->Fill(diff);
  	htoyResult_pull[ieta][ipt]->Fill(pull);

       } // end loops of toys
    
      delete fsig;
      delete fbkg;

    } // end of loop over pt bins

  } 

  for(int ieta=0; ieta < nEtaBin; ieta++){
    TFile* outFile = new TFile(Form("fitterTest_%s.root",dec[ieta]),
			       "recreate");
   
    for(int ipt=0; ipt < nPtBin; ipt++){
     
      htoyResult_pull[ieta][ipt]->Write();
      htoyResult_diff[ieta][ipt]->Write();
    }
    outFile->Close();
  }
    
   

}
