

#include <iostream>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TH1.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h> 
#include <TCanvas.h>
#include <TLegend.h>
#include "purity_twobin.C"
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
    //     {
    //       933.7,
    //       265.6,
    //       109.7,
    //       21.6,
    //       3.0
    //     },
    //     {
    //       712.5,
    //       210.1,
    //       74.2,
    //       16.1,
    //       0.8
    //     }
    {
      964.3,
      274.1,
      112.0,
      22.8,
      3.0
    },
    {
      548.7,
      210.8,
      70.7,
      14.9,
      0.8
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

// calculate errors from weighted histograms
void histoErrorEiko(TH1F* h, int maxbin, Double_t& total_value, Double_t& total_err)
{

  total_value = total_err = 0;
  for(int i=1; i<= maxbin; i++)
    {
      total_value += h->GetBinContent(i);
      total_err   += pow(h->GetBinError(i),2);
      
    }

  total_err = sqrt(total_err);
  return;

}


Double_t sigPar[12]={0.58544,1.39853,0.482806,0.263328,64.9728,-1.28153,-0.108433,0.0446752,0.5,0.01,1,1};


void toyMC_data(bool dataDriven=true){
	

  TF1 *fsig = new TF1("fsig", exp_conv, -1., 11., 12);
  for(int ipar=0; ipar<4; ipar++)
    fsig->SetParameter(ipar,sigPar[ipar]);


  TF1 *fbkg = new TF1("fbkg", expinv_power, -1., 11., 12);
  for(int ipar=0;ipar<12;ipar++)
    {
      fbkg->SetParameter(ipar,sigPar[ipar]);
    }


  cout << "use dataDriven = " << dataDriven << endl;
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


  
  std::string dataFile     = "rawYield.root";
  std::string templateFile = "rawYield.root";

  TFile* inf_data = new TFile(dataFile.data());
  TFile* inf_template = new TFile(templateFile.data());

  TH1F* htoyResult_pull[nEtaBin][nPtBin];

  
  char* dec[2] = {"EB","EE"};
  for(int ieta=0; ieta<nEtaBin; ieta++){
    for(int ipt=0; ipt < nPtBin; ipt++){

      htoyResult_pull[ieta][ipt] = new TH1F(Form("hpull_Eta_%.2f_%.2f_Et_%d_%d",		
						 fBinsEta[ieta*2],fBinsEta[ieta*2+1],
						 (int)fBinsPt[ipt], (int)fBinsPt[ipt+1]),
					    "",50,-5.0,5.0);


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
  htoyMC_data->SetName("htoyMC_data");
  htoyMC_data->Reset();

  for(int ieta=0; ieta< nEtaBin; ieta++){
    for(int ipt=0; ipt < nPtBin; ipt++){

      // loops over toys
      for(int iexp=0; iexp<NEXP; iexp++){
	htoyMC_data->Reset();
	UInt_t nowSeed = (unsigned long)gSystem->Now();
	gRandom->SetSeed(nowSeed);
	int ndata = gRandom->Poisson(hdata_data[ieta][ipt]->Integral());      
	int nsig  = gRandom->Poisson(nsig_data[ieta][ipt]);
	int nbkg  = ndata - nsig;

	// first fill signal
	htoyMC_data->FillRandom(hTemplate_S[ieta][ipt], nsig);
	//       htoyMC_data->FillRandom("fsig", nsig);
	// second fill background
	htoyMC_data->FillRandom(hTemplate_B[ieta][ipt], nbkg);
	//       htoyMC_data->FillRandom("fbkg", nbkg);

	cout << "ndata = " << ndata << "\t nsig = " << nsig << " \t nbkg = " << 
	  nbkg << endl;

	Double_t scaleEff = 1.0;
	Double_t scaleNumber[20];
	for(int i=0;i<20;i++)scaleNumber[i]=1.;


	// 2nd, get unbinned fit
	//       Double_t* FuncFitResult;

	//       FuncFitResult = Ifit(htoyMC_data,
	// 			   hTemplate_S[ieta][ipt],
	// 			   hTemplate_B[ieta][ipt],0,"RS_100604_fix.dat",
	// 			   fBinsEta[ieta*2],fBinsEta[ieta*2+1],
	// 			   fBinsPt[ipt],fBinsPt[ipt+1],
	// 			   NULL,
	// 			   scaleNumber);	

 

	//       Double_t nsigfunc    = FuncFitResult[0]/scaleEff;
	//       Double_t errnsigfunc = FuncFitResult[1]/scaleEff;
	//       Double_t nbkgfunc    = FuncFitResult[2]/scaleEff;
	//       Double_t errnbkgfunc = FuncFitResult[3]/scaleEff;
  
      
	//       Double_t pull = (nsigfunc - nsig_data[ieta][ipt])/errnsigfunc;
      
	Double_t* TemplateFitResult;
	TemplateFitResult = IfitBin(htoyMC_data,
				    hTemplate_S[ieta][ipt],
				    hTemplate_B[ieta][ipt]);

 
	Double_t nsigtemplate    = TemplateFitResult[0]/scaleEff;
	Double_t errnsigtemplate = TemplateFitResult[1]/scaleEff;
	Double_t nbkgtemplate    = TemplateFitResult[2]/scaleEff;
	Double_t errnbkgtemplate = TemplateFitResult[3]/scaleEff;

	Double_t pull = (nsigtemplate - nsig_data[ieta][ipt])/errnsigtemplate;
      
      
	htoyResult_pull[ieta][ipt]->Fill(pull);
      } // end loops of toys
    
    } // end of loop over pt bins

  } 

  for(int ieta=0; ieta < nEtaBin; ieta++){
    TFile* outFile = new TFile(Form("pull_%s.root",dec[ieta]),
			       "recreate");
   
    for(int ipt=0; ipt < nPtBin; ipt++){
     
      htoyResult_pull[ieta][ipt]->Write();
    }
    outFile->Close();
  }
    
   

}
