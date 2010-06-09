/*==============================================================================================
  Use Rong-Shyang's style histograms to call template binned likelihood fitter, unbined fitter, 
 two-bin fitter and compare results with MC truth. In addition, signal yields are saved in 
 txt and tex formats.

  .L call_allfitters_rs.C++

  (1) if fit data, and get corrected yield:

  call_allfitters_rs(true, true, luminosity in nb)

  // template MC samples were produced normalized to 100/nb

  (2) if fit fake data from MC and get corrected yield:
 
  call_allfitters_rs(false, true)


 Shin-Shan Eiko Yu, 2010.06.01

  ==============================================================================================*/



#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TH1.h>
#include <string>
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
const double fBinsEta[]={0,1.45,1.7,2.5};
const double fBinsPt[]={15,20,30,50,80,120};
const int nPtBin = sizeof(fBinsPt)/sizeof(fBinsPt[0])-1;
const int nEtaBin = (sizeof(fBinsEta)/sizeof(fBinsEta[0]))/2;

const int REBINNINGS=6;

double tempNumbers[10][2]={
  {7158.5 , 147.3},
  {4788.7 , 109.5} ,
  {896.7 , 40.8},
  {222.8 , 17.0},
  {27.1 , 6.1},
  {6712.3 , 127.6},
  {3619.0 , 92.8},
  {559.7 , 33.9},
  {161.9 , 15.5},
  {17.1 , 4.8}
};

double tempPurityNumbers[10][2]={
  {0.296 +- 0.004},
  {0.344 +- 0.006},
  {0.445 +- 0.012},
  {0.684 +- 0.021},
  {0.732 +- 0.061},
  {0.188 +- 0.003},
  {0.226 +- 0.005},
  {0.227 +- 0.011},
  {0.510 +- 0.027},
  {0.548 +- 0.081}

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
void histoErrorEiko(TH1F* h, int maxbin, double& total_value, double& total_err)
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

void call_allfitters_rs(bool fitData=false, bool dataDriven=false,bool doEffCorr=false, double lumi=8.02){
	
  cout << "do efficiency correction = " << doEffCorr << endl;
  cout << "fit data = " << fitData << endl;
  cout << "use dataDriven = " << dataDriven << endl;
  std::string isDataName = fitData? "data": "MC";

  setTDRStyle();
  // settings for purity TGraphAsymmetryErrors
  double fBinsPtMidPoint[nPtBin]={0};
  double fBinsPtError[nPtBin]={0};
  
  for(int ipt=0; ipt < nPtBin; ipt++)
    {
      fBinsPtMidPoint[ipt] = 0.5*(fBinsPt[ipt+1]+fBinsPt[ipt]);
      fBinsPtError[ipt] = 0.5*(fBinsPt[ipt+1]-fBinsPt[ipt]);
    }

  // first get efficiency
  double eff[nEtaBin][nPtBin]={{0}};
  double efferr[nEtaBin][nPtBin]={{0}};

  TFile* inf_eff = new TFile("Efficiency_allMC.root");
  TH1F* phoEff[nEtaBin];
  char tmp[1000];
  
  for(int i=0;i<nEtaBin;i++)
    {
      sprintf(tmp,"phoEff_Pt%d",i);
      phoEff[i] = (TH1F*)inf_eff->FindObjectAny(tmp);
      for(int ipt=0; ipt<nPtBin;ipt++)
	{
	  eff[i][ipt]    =phoEff[i]->GetBinContent(ipt+1);
	  efferr[i][ipt] =phoEff[i]->GetBinError(ipt+1); 
	  cout << "Efficiency in ieta=" << i << ", ipt = " << ipt << " bin = " << eff[i][ipt] << 
	    " +- " << efferr[i][ipt] << endl;
	}
   
    }

  double mcpurity[nEtaBin][nPtBin]={{0}};
  double mcpurity_err[nEtaBin][nPtBin]={{0}};
  double nsig_mc[nEtaBin][nPtBin]={{0}};
  double nbkg_mc[nEtaBin][nPtBin]={{0}};
  double nsig_err_mc[nEtaBin][nPtBin]={{0}};
  double nbkg_err_mc[nEtaBin][nPtBin]={{0}};

  double purity_func[nEtaBin][nPtBin]={{0}};
  double purity_err_func[nEtaBin][nPtBin]={{0}};
  double nsig_func[nEtaBin][nPtBin]={{0}};
  double nbkg_func[nEtaBin][nPtBin]={{0}};
  double nsig_err_func[nEtaBin][nPtBin]={{0}};
  double nbkg_err_func[nEtaBin][nPtBin]={{0}};

  double purity_template[nEtaBin][nPtBin]={{0}};
  double purity_err_template[nEtaBin][nPtBin]={{0}};
  double nsig_template[nEtaBin][nPtBin]={{0}};
  double nbkg_template[nEtaBin][nPtBin]={{0}};
  double nsig_err_template[nEtaBin][nPtBin]={{0}};
  double nbkg_err_template[nEtaBin][nPtBin]={{0}};
 
  double purity_2bin[nEtaBin][nPtBin]={{0}};
  double purity_err_2bin[nEtaBin][nPtBin]={{0}};
  double nsig_2bin[nEtaBin][nPtBin]={{0}};
  double nsig_err_2bin[nEtaBin][nPtBin]={{0}};


  TH1F* hTemplate_S[nEtaBin][nPtBin];
  TH1F* hTemplate_B[nEtaBin][nPtBin];
  TH1F* hTemplate_data[nEtaBin][nPtBin];
  TH1F* hdata_S[nEtaBin][nPtBin];
  TH1F* hdata_B[nEtaBin][nPtBin];
  TH1F* hdata_data[nEtaBin][nPtBin];


  
  std::string dataFile     = fitData? "TightEESBData_template.root":"fakedata_template.root";
  std::string templateFile = fitData? "TightEESBData_template.root":"template_template.root";
//   std::string dataFile     = fitData? "template_comb3Iso.root":"template_comb3Iso_test.root";
//   std::string templateFile = fitData? "template_comb3Iso.root":"template_comb3Iso_template.root";


  TFile* inf_data = new TFile(dataFile.data());
  TFile* inf_template = new TFile(templateFile.data());
  
  char* dec[2] = {"EB","EE"};
  for(int ieta = 0; ieta < nEtaBin; ieta++){
     for(int ipt=0; ipt < nPtBin; ipt++){
       // getting histograms from data root file
       sprintf(tmp,"h_%s_comb3Iso_EGdata_pt%d",dec[ieta],(int)fBinsPt[ipt]);
       cout << "looking for histogram " << tmp << " in file " << 
	 inf_data->GetName() << endl;
       hdata_data[ieta][ipt] = (TH1F*)inf_data->FindObjectAny(tmp);
       hdata_data[ieta][ipt]->Rebin(REBINNINGS);
       // setting titles
       sprintf(tmp,"%.2f < |#eta(#gamma)| < %.2f, %d < p_{T}(#gamma) < %d GeV",
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1],
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1]);
       hdata_data[ieta][ipt]->SetTitle(tmp);


       sprintf(tmp,"h_%s_comb3Iso_sig_pt%d",dec[ieta],(int)fBinsPt[ipt]);
       cout << "looking for histogram " << tmp << " in file " << 
	 inf_data->GetName() << endl;
       hdata_S[ieta][ipt] = (TH1F*)inf_data->FindObjectAny(tmp);
       hdata_S[ieta][ipt]->Rebin(REBINNINGS);

       sprintf(tmp,"h_%s_comb3Iso_bkg_pt%d",dec[ieta],(int)fBinsPt[ipt]);
       cout << "looking for histogram " << tmp << " in file " << 
	 inf_data->GetName() << endl;
       hdata_B[ieta][ipt] = (TH1F*)inf_data->FindObjectAny(tmp);
       hdata_B[ieta][ipt]->Rebin(REBINNINGS);


       // getting histogram for template root file
       sprintf(tmp,"h_%s_comb3Iso_EGdata_pt%d",dec[ieta],(int)fBinsPt[ipt]);
       cout << "looking for histogram " << tmp << " in file " << 
	 inf_template->GetName() << endl;
       hTemplate_data[ieta][ipt] = (TH1F*)inf_template->FindObjectAny(tmp);
       hTemplate_data[ieta][ipt]->Rebin(REBINNINGS);
       sprintf(tmp,"%.2f < |#eta(#gamma)| < %.2f, %d < p_{T}(#gamma) < %d GeV",
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1],
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1]);
       hTemplate_data[ieta][ipt]->SetTitle(tmp);

       sprintf(tmp,"h_%s_comb3Iso_sig_pt%d",dec[ieta],(int)fBinsPt[ipt]);
       cout << "looking for histogram " << tmp << " in file " << 
	 inf_template->GetName() << endl;
       hTemplate_S[ieta][ipt] = (TH1F*)inf_template->FindObjectAny(tmp);
       hTemplate_S[ieta][ipt]->Rebin(REBINNINGS);

       sprintf(tmp,"h_%s_comb3Iso_bkg_pt%d",dec[ieta],(int)fBinsPt[ipt]);
       cout << "looking for histogram " << tmp << " in file " << 
	 inf_template->GetName() << endl;
       hTemplate_B[ieta][ipt] = (TH1F*)inf_template->FindObjectAny(tmp);
       hTemplate_B[ieta][ipt]->Rebin(REBINNINGS);

     }
  }

  // files for data driven methods
  std::string SpikeFile       = "spike_withID.root";
  TFile* inf_dataSpike = new TFile(SpikeFile.data());
  TH1F* hdata_Spike  = (TH1F*)inf_dataSpike->FindObjectAny("h_EB_comb3Iso_EGdata_SIG");
  hdata_Spike->Rebin(REBINNINGS);
  cout << "looking for histogram " << hdata_Spike->GetName() << " in file " << 
	 inf_dataSpike->GetName() << endl;


//   std::string SBFile       = "template_comb3Iso_SBMC.root";
  std::string SBFile       = "TightEESBData_template.root";
  TFile* inf_dataSB = new TFile(SBFile.data());
  TH1F* hdata_SB[nEtaBin];
  for(int ieta=0; ieta < 2; ieta++)
    {
      sprintf(tmp,"h_%s_comb3IsoSB_EGdata_pt15",dec[ieta]);
       cout << "looking for histogram " << tmp << " in file " << 
	 inf_dataSB->GetName() << endl;
       hdata_SB[ieta]= (TH1F*)inf_dataSB->FindObjectAny(tmp);
       hdata_SB[ieta]->Rebin(REBINNINGS);
    }



  int IntegralBin= hdata_data[0][0]->GetNbinsX()*0.5;

  for(int ieta = 0; ieta < nEtaBin; ieta++){
    for(int ipt=0; ipt < nPtBin; ipt++){
      

      double scaleEff = doEffCorr? eff[ieta][ipt]: 1.0;

      // 1st, get MC truth
      double nsigmc=0, errnsigmc=0;
      double nbkgmc=0, errnbkgmc=0;

      double For5GeVnsigmc=0, For5GeVerrnsigmc=0;
      double For5GeVnbkgmc=0, For5GeVerrnbkgmc=0;
      
      if(fitData)
	{
	  histoErrorEiko(hTemplate_S[ieta][ipt], hTemplate_S[ieta][ipt]->GetNbinsX(), 
			 nsigmc, errnsigmc);
	  histoErrorEiko(hTemplate_B[ieta][ipt], hTemplate_B[ieta][ipt]->GetNbinsX(), 
			 nbkgmc, errnbkgmc);	  
 	  nsigmc  *= (1.0/scaleEff) * (lumi/100.0);
 	  errnsigmc  *= (1.0/scaleEff) * (lumi/100.0);

 	  nbkgmc  *= (1.0/scaleEff) * (lumi/100.0);
 	  errnbkgmc  *= (1.0/scaleEff) * (lumi/100.0);

	  histoErrorEiko(hTemplate_S[ieta][ipt], IntegralBin, 
			 For5GeVnsigmc, For5GeVerrnsigmc);
	  histoErrorEiko(hTemplate_B[ieta][ipt], IntegralBin, 
			 For5GeVnbkgmc, For5GeVerrnbkgmc);	  
 	  For5GeVnsigmc  *= (1.0/scaleEff) * (lumi/100.0);
 	  For5GeVerrnsigmc  *= (1.0/scaleEff) * (lumi/100.0);

 	  For5GeVnbkgmc  *= (1.0/scaleEff) * (lumi/100.0);
 	  For5GeVerrnbkgmc  *= (1.0/scaleEff) * (lumi/100.0);

	}
      else
	{
	  histoErrorEiko(hdata_S[ieta][ipt], hdata_S[ieta][ipt]->GetNbinsX(), 
			 nsigmc, errnsigmc);
	  histoErrorEiko(hdata_B[ieta][ipt], hdata_B[ieta][ipt]->GetNbinsX(), 
			 nbkgmc, errnbkgmc);	  
 	  nsigmc  *= (1.0/scaleEff);
 	  errnsigmc  *= (1.0/scaleEff);

 	  nbkgmc  *= (1.0/scaleEff);
 	  errnbkgmc  *= (1.0/scaleEff);

	  histoErrorEiko(hdata_S[ieta][ipt], IntegralBin, 
			 For5GeVnsigmc, For5GeVerrnsigmc);
	  histoErrorEiko(hdata_B[ieta][ipt], IntegralBin, 
			 For5GeVnbkgmc, For5GeVerrnbkgmc);	  
 	  For5GeVnsigmc  *= (1.0/scaleEff);
 	  For5GeVerrnsigmc  *= (1.0/scaleEff);

 	  For5GeVnbkgmc  *= (1.0/scaleEff);
 	  For5GeVerrnbkgmc  *= (1.0/scaleEff);
	}



      if(nsigmc+nbkgmc < 1e-6)continue;

      Double_t puritymc = 0;
      Double_t errpuritymc = 0;
 
      ratioErr(For5GeVnsigmc, For5GeVerrnsigmc, For5GeVnbkgmc, For5GeVerrnbkgmc,
	       puritymc, errpuritymc);

      nsig_mc[ieta][ipt] = nsigmc;
      nsig_err_mc[ieta][ipt] = errnsigmc;
      nbkg_mc[ieta][ipt] = nbkgmc;
      nbkg_err_mc[ieta][ipt] = errnbkgmc;
      
      mcpurity[ieta][ipt] = puritymc;
      mcpurity_err[ieta][ipt] = errpuritymc;

      

      // 2nd, get unbinned fit
      Double_t* FuncFitResult;
      if(fitData && !dataDriven)
 	FuncFitResult = Ifit(hdata_data[ieta][ipt],
 			     hTemplate_S[ieta][ipt],
 			     hTemplate_B[ieta][ipt],1,"RS_100604_fix.dat",
 			     fBinsEta[ieta*2],fBinsEta[ieta*2+1],
 			     fBinsPt[ipt],fBinsPt[ipt+1]);	

      else if(fitData && dataDriven && ieta==0)
  	FuncFitResult = Ifit(hdata_data[ieta][ipt], 
     			     hdata_Spike,
   			     hdata_SB[ieta],1,"RS_100604_fix.dat",
  			     fBinsEta[ieta*2],fBinsEta[ieta*2+1],
  			     fBinsPt[ipt],fBinsPt[ipt+1],
 			     hTemplate_B[ieta][ipt]
			     );	

      else if(fitData && dataDriven && ieta==1)
   	FuncFitResult = Ifit(hdata_data[ieta][ipt],
   			     hTemplate_S[ieta][ipt],
    			     hdata_SB[ieta],1,"RS_100604_fix.dat",
   			     fBinsEta[ieta*2],fBinsEta[ieta*2+1],
   			     fBinsPt[ipt],fBinsPt[ipt+1],
 			     hTemplate_B[ieta][ipt]);	
      else if(!fitData && dataDriven)
   	FuncFitResult = Ifit(hdata_data[ieta][ipt],
   			     hTemplate_S[ieta][ipt],
    			     hdata_SB[ieta],0,"RS_100604_fix.dat",
   			     fBinsEta[ieta*2],fBinsEta[ieta*2+1],
   			     fBinsPt[ipt],fBinsPt[ipt+1],
  			     hTemplate_B[ieta][ipt]
			     );	
      else
	FuncFitResult = Ifit(hdata_data[ieta][ipt],
			     hTemplate_S[ieta][ipt],
			     hTemplate_B[ieta][ipt],0,"RS_100604_fix.dat",
			     fBinsEta[ieta*2],fBinsEta[ieta*2+1],
			     fBinsPt[ipt],fBinsPt[ipt+1]);	
 
      Double_t nsigfunc    = FuncFitResult[0]/scaleEff;
      Double_t errnsigfunc = FuncFitResult[1]/scaleEff;
      Double_t nbkgfunc    = FuncFitResult[2]/scaleEff;
      Double_t errnbkgfunc = FuncFitResult[3]/scaleEff;
  
      Double_t For5GeVnsigfunc    = FuncFitResult[4]/scaleEff;
      Double_t For5GeVerrnsigfunc = FuncFitResult[5]/scaleEff;
      Double_t For5GeVnbkgfunc    = FuncFitResult[6]/scaleEff;
      Double_t For5GeVerrnbkgfunc = FuncFitResult[7]/scaleEff;

      Double_t purityfunc = 0;
      Double_t errpurityfunc = 0;
 
      ratioErr(For5GeVnsigfunc, For5GeVerrnsigfunc, For5GeVnbkgfunc, For5GeVerrnbkgfunc,
	       purityfunc, errpurityfunc);

      purity_func[ieta][ipt]    = purityfunc;
      purity_err_func[ieta][ipt]= errpurityfunc;
      nsig_func[ieta][ipt]      = nsigfunc;
      nsig_err_func[ieta][ipt]  = errnsigfunc;
      nbkg_func[ieta][ipt]      = nbkgfunc;
      nbkg_err_func[ieta][ipt]  = errnbkgfunc;


      // 3rd, get Template Fit result
      Double_t* TemplateFitResult;
      if(!dataDriven)
       TemplateFitResult = IfitBin(hdata_data[ieta][ipt],
 				  hTemplate_S[ieta][ipt],
 				  hTemplate_B[ieta][ipt]);
      else if(dataDriven && !fitData)
	TemplateFitResult = IfitBin(hdata_data[ieta][ipt],
				    hTemplate_S[ieta][ipt],
				    hdata_SB[ieta]);

      else if(dataDriven && ieta==0)
	TemplateFitResult = IfitBin(hdata_data[ieta][ipt],
 				    hdata_Spike,
 				    hdata_SB[ieta]);
      else if(dataDriven && ieta==1)
  	TemplateFitResult = IfitBin(hdata_data[ieta][ipt],
  				    hTemplate_S[ieta][ipt],
//   				    hTemplate_B[ieta][ipt]);
 				    hdata_SB[ieta]);

 
      Double_t nsigtemplate    = TemplateFitResult[0]/scaleEff;
      Double_t errnsigtemplate = TemplateFitResult[1]/scaleEff;
      Double_t nbkgtemplate    = TemplateFitResult[2]/scaleEff;
      Double_t errnbkgtemplate = TemplateFitResult[3]/scaleEff;

      Double_t For5GeVnsigtemplate    = TemplateFitResult[4]/scaleEff;
      Double_t For5GeVerrnsigtemplate = TemplateFitResult[5]/scaleEff;
      Double_t For5GeVnbkgtemplate    = TemplateFitResult[6]/scaleEff;
      Double_t For5GeVerrnbkgtemplate = TemplateFitResult[7]/scaleEff;

      Double_t puritytemplate = 0;
      Double_t errpuritytemplate = 0;
 
      ratioErr(For5GeVnsigtemplate, For5GeVerrnsigtemplate, For5GeVnbkgtemplate, For5GeVerrnbkgtemplate,
	       puritytemplate, errpuritytemplate);

//       purity_template[ieta][ipt]    = tempPurityNumbers[ipt+ieta*5][0];
//       purity_err_template[ieta][ipt]= tempPurityNumbers[ipt+ieta*5][1];
//       nsig_template[ieta][ipt]      = tempNumbers[ipt+ieta*5][0];
//       nsig_err_template[ieta][ipt]  = tempNumbers[ipt+ieta*5][1];
      purity_template[ieta][ipt]    = puritytemplate;
      purity_err_template[ieta][ipt]= errpuritytemplate;
      nsig_template[ieta][ipt]      = nsigtemplate;
      nsig_err_template[ieta][ipt]  = errnsigtemplate;
      nbkg_template[ieta][ipt]      = nbkgtemplate;
      nbkg_err_template[ieta][ipt]  = errnbkgtemplate;



      // 4th, get two bin result
      Double_t* TwoBinResult;
      if(!dataDriven)
       TwoBinResult = purity_twobin(hdata_data[ieta][ipt],
				    hTemplate_S[ieta][ipt],
				    hTemplate_B[ieta][ipt],IntegralBin);
      else if(dataDriven && ieta==0)
	TwoBinResult = purity_twobin(hdata_data[ieta][ipt],
				     hTemplate_S[ieta][ipt],
 				     hdata_SB[ieta],IntegralBin);
      else if(dataDriven && ieta==1)
	TwoBinResult = purity_twobin(hdata_data[ieta][ipt],
				     hTemplate_S[ieta][ipt],
				     hdata_SB[ieta],IntegralBin);

      Double_t ndata = hdata_data[ieta][ipt]->Integral(1,IntegralBin);
      purity_2bin[ieta][ipt]     = TwoBinResult[2];
      purity_err_2bin[ieta][ipt] = TwoBinResult[3];
      nsig_2bin[ieta][ipt]       = purity_2bin[ieta][ipt]*ndata/scaleEff;
      nsig_err_2bin[ieta][ipt]   = purity_err_2bin[ieta][ipt]*ndata/scaleEff;

    } // end of loop over pt bins
  } // end of loop over eta bins


  ofstream fout;
  std::string corr= doEffCorr? "_corr":"";
  fout.open(Form("yield%s.txt",corr.data()));

  ofstream fout2;
  fout2.open(Form("purity%s.txt",corr.data()));
 
  ofstream texout;
  texout.open(Form("yield%s.tex",corr.data()));

  ofstream texout2;
  texout2.open(Form("purity%s.tex",corr.data()));
 
  // printing out the results
  for(int ieta=0; ieta < nEtaBin; ieta++){
    for(int ipt=0; ipt < nPtBin; ipt++){
      cout << endl;
      cout << " ===== for ieta=" << ieta << ", ipt=" << ipt << " ======== " << endl;
      cout << " MC truth : " << endl;
      cout << " purity = " << mcpurity[ieta][ipt] << " +- " << mcpurity_err[ieta][ipt] << endl;
      cout << " nsig = " << nsig_mc[ieta][ipt] << " +- " << nsig_err_mc[ieta][ipt] << endl;
      cout << " nbkg = " << nbkg_mc[ieta][ipt] << " +- " << nbkg_err_mc[ieta][ipt] << endl;
      cout << endl;

      fout << Form("%d-%d",(int)fBinsPt[ipt], (int)fBinsPt[ipt+1]) << " ";
      fout << Form("& %.1f +- %.1f", nsig_mc[ieta][ipt],nsig_err_mc[ieta][ipt]) << " ";

      fout2 << Form("%d-%d",(int)fBinsPt[ipt], (int)fBinsPt[ipt+1]) << " ";
      fout2 << Form("& %.3f +- %.3f", mcpurity[ieta][ipt],mcpurity_err[ieta][ipt]) << " ";

      texout << Form("%d-%d",(int)fBinsPt[ipt], (int)fBinsPt[ipt+1]) << " ";
      texout << Form("& %.1f $\\pm$ %.1f", nsig_mc[ieta][ipt],nsig_err_mc[ieta][ipt]) << " ";

      texout2 << Form("%d-%d",(int)fBinsPt[ipt], (int)fBinsPt[ipt+1]) << " ";
      texout2 << Form("& %.3f $\\pm$ %.3f", mcpurity[ieta][ipt],mcpurity_err[ieta][ipt]) << " ";

      cout << " Unbinned fit result : " << endl;
      cout << " purity = " << purity_func[ieta][ipt] << " +- " << purity_err_func[ieta][ipt] << endl;
      cout << " nsig = " <<   nsig_func[ieta][ipt] << " +- " << nsig_err_func[ieta][ipt] << endl;
      cout << " nbkg = " <<   nbkg_func[ieta][ipt] << " +- " << nbkg_err_func[ieta][ipt] << endl;
      cout << endl;

      fout  << Form("& %.1f +- %.1f", nsig_func[ieta][ipt], nsig_err_func[ieta][ipt]) << " ";
      fout2 << Form("& %.3f +- %.3f", purity_func[ieta][ipt], purity_err_func[ieta][ipt]) << " ";

      texout << Form("& %.1f $\\pm$ %.1f", nsig_func[ieta][ipt], nsig_err_func[ieta][ipt]) << " ";
      texout2 << Form("& %.3f $\\pm$ %.3f", purity_func[ieta][ipt], purity_err_func[ieta][ipt]) << " ";

      cout << " Template Fit result : " << endl;
      cout << " purity = " << purity_template[ieta][ipt] << " +- " << purity_err_template[ieta][ipt] << endl;
      cout << " nsig = " <<   nsig_template[ieta][ipt] << " +- " << nsig_err_template[ieta][ipt] << endl;
      cout << " nbkg = " <<   nbkg_template[ieta][ipt] << " +- " << nbkg_err_template[ieta][ipt] << endl;
      cout << endl;

      fout << Form("& %.1f +- %.1f", nsig_template[ieta][ipt], nsig_err_template[ieta][ipt]) << " ";
      fout2 << Form("& %.3f +- %.3f", purity_template[ieta][ipt], purity_err_template[ieta][ipt]) << " ";

      texout << Form("& %.1f $\\pm$ %.1f", nsig_template[ieta][ipt], nsig_err_template[ieta][ipt]) << " ";
      texout2 << Form("& %.3f $\\pm$ %.3f", purity_template[ieta][ipt], purity_err_template[ieta][ipt]) << " ";

      cout << " Two-bin result : " << endl;
      cout << " purity = " << purity_2bin[ieta][ipt] << " +- " << purity_err_2bin[ieta][ipt] << endl;
      cout << " nsig = " <<   nsig_2bin[ieta][ipt] << " +- " << nsig_err_2bin[ieta][ipt] << endl;

      fout << Form("& %.1f +- %.1f", nsig_2bin[ieta][ipt], nsig_err_2bin[ieta][ipt]) << " ";
      fout2 << Form("& %.3f +- %.3f", purity_2bin[ieta][ipt], purity_err_2bin[ieta][ipt]) << " ";

      texout << "\\\\" << endl;
      fout << endl;

      texout2 << "\\\\" << endl;
      fout2 << endl;

    }
  }

  texout.close();
  fout.close();

  texout2.close();
  fout2.close();
  
  TCanvas *c1 = new TCanvas("c1","",600,600);
  c1->Draw();
  TH2F *h2 = new TH2F("h2","",100, 0., 140, 100, 0, 1.2);
  h2->SetNdivisions(505,"XY");
  h2->SetXTitle("p_{T} (GeV)");
  h2->SetYTitle("purity");
  h2->Draw();

  TGraphAsymmErrors *tg_mc_purity_EB = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, 
							     mcpurity[0],
							     fBinsPtError,fBinsPtError,
							     mcpurity_err[0],mcpurity_err[0]);
  tg_mc_purity_EB->SetName("tg_mc_purity_EB");
  tg_mc_purity_EB->SetMarkerSize(1.);
  tg_mc_purity_EB->SetLineWidth(2);
  tg_mc_purity_EB->Draw("p e same");


  TGraphAsymmErrors *tgrs_purity_EB = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, 
							    purity_func[0], 
							    fBinsPtError, fBinsPtError, 
							    purity_err_func[0], purity_err_func[0]);
  tgrs_purity_EB->SetName("tgrs_purity_EB");
  tgrs_purity_EB->SetMarkerColor(4);
  tgrs_purity_EB->SetMarkerStyle(25);
  tgrs_purity_EB->SetLineColor(4);
  tgrs_purity_EB->SetLineWidth(2);
  tgrs_purity_EB->Draw("p e same");

  TGraphAsymmErrors *tgrs_purity_template_EB = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, 
								     purity_template[0], 
								     fBinsPtError,fBinsPtError, 
								     purity_err_template[0],
								     purity_err_template[0]);
  tgrs_purity_template_EB->SetName("tgrs_purity_template_EB");
  tgrs_purity_template_EB->SetMarkerColor(6);
  tgrs_purity_template_EB->SetMarkerStyle(26);
  tgrs_purity_template_EB->SetLineColor(6);
  tgrs_purity_template_EB->SetLineWidth(2);
  tgrs_purity_template_EB->Draw("p e same");

  TGraphAsymmErrors *tgrs_purity_twobin_EB = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, 
								   purity_2bin[0], 
								   fBinsPtError,fBinsPtError, 
								   purity_err_2bin[0],
								   purity_err_2bin[0]);
  tgrs_purity_twobin_EB->SetName("tgrs_purity_twobin_EB");
  tgrs_purity_twobin_EB->SetMarkerColor(kGreen);
  tgrs_purity_twobin_EB->SetMarkerStyle(3 );
  tgrs_purity_twobin_EB->SetLineColor(kGreen);
  tgrs_purity_twobin_EB->SetLineWidth(2);
//   tgrs_purity_twobin_EB->Draw("p e same");

  tgrs_purity_template_EB->Draw("p e same");
  tg_mc_purity_EB->Draw("p same");

  TLegend *tleg = new TLegend(0.15, 0.65, 0.55, 0.92);
  tleg->SetHeader("EB Photons");
  tleg->SetFillColor(0);
  tleg->SetFillStyle(0);
  tleg->SetShadowColor(0);
  tleg->SetBorderSize(0);
  tleg->AddEntry(tg_mc_purity_EB,"MC Prediction","pl");
  tleg->AddEntry(tgrs_purity_EB,"Function Fit","pl");
//   tleg->AddEntry(tgrs_purity_template_EB,"Function Fit with #sigma_{i#etai#eta} SB","pl");
  tleg->AddEntry(tgrs_purity_template_EB,"Template Fit","pl");
//   tleg->AddEntry(tgrs_purity_twobin_EB,"two bin","pl");
  tleg->Draw();


  c1->SaveAs(Form("plots/%s_allfitter_purity_EB.eps",isDataName.data()));
  c1->SaveAs(Form("plots/%s_allfitter_purity_EB.gif",isDataName.data()));
  c1->SaveAs(Form("plots/%s_allfitter_purity_EB.C",isDataName.data()));
  

  TCanvas *c2 = new TCanvas("c2","",600,600);
  c2->Draw();
  h2->Draw();

  TGraphAsymmErrors *tg_mc_purity_EE = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, 
							     mcpurity[1],
							     fBinsPtError,fBinsPtError,
							     mcpurity_err[1],mcpurity_err[1]);
  tg_mc_purity_EE->SetName("tg_mc_purity_EE");
  tg_mc_purity_EE->SetMarkerSize(1.);
  tg_mc_purity_EE->SetLineWidth(2);
  tg_mc_purity_EE->Draw("p e same");


  TGraphAsymmErrors *tgrs_purity_EE = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, 
							    purity_func[1], 
							    fBinsPtError, fBinsPtError, 
							    purity_err_func[1],purity_err_func[1]);
  tgrs_purity_EE->SetName("tgrs_purity_EE");
  tgrs_purity_EE->SetMarkerColor(4);
  tgrs_purity_EE->SetMarkerStyle(25);
  tgrs_purity_EE->SetLineColor(4);
  tgrs_purity_EE->SetLineWidth(2);
  tgrs_purity_EE->Draw("p e same");

  TGraphAsymmErrors *tgrs_purity_template_EE = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, 
								     purity_template[1], 
								     fBinsPtError,fBinsPtError, 
								     purity_err_template[1],
								     purity_err_template[1]);
  tgrs_purity_template_EE->SetName("tgrs_purity_template_EE");
  tgrs_purity_template_EE->SetMarkerColor(6);
  tgrs_purity_template_EE->SetMarkerStyle(26);
  tgrs_purity_template_EE->SetLineColor(6);
  tgrs_purity_template_EE->SetLineWidth(2);
  tgrs_purity_template_EE->Draw("p e same");

  TGraphAsymmErrors *tgrs_purity_twobin_EE = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, 
								   purity_2bin[1], 
								   fBinsPtError, fBinsPtError,
								   purity_err_2bin[1], 
								   purity_err_2bin[1]);
  tgrs_purity_twobin_EE->SetName("tgrs_purity_twobin_EE");
  tgrs_purity_twobin_EE->SetMarkerColor(kGreen);
  tgrs_purity_twobin_EE->SetMarkerStyle(3);
  tgrs_purity_twobin_EE->SetLineColor(kGreen);
  tgrs_purity_twobin_EE->SetLineWidth(2);
//   tgrs_purity_twobin_EE->Draw("p e same");

  tgrs_purity_template_EE->Draw("p e same");
  tg_mc_purity_EE->Draw("p same");

  TLegend *tleg2 = new TLegend(0.15, 0.65, 0.55, 0.92);
  tleg2->SetHeader("EE Photons");
  tleg2->SetFillColor(0);
  tleg2->SetFillStyle(0);
  tleg2->SetShadowColor(0);
  tleg2->SetBorderSize(0);
  tleg2->AddEntry(tg_mc_purity_EE,"MC Prediction","pl");
  tleg2->AddEntry(tgrs_purity_EE,"Function Fit","pl");
//   tleg2->AddEntry(tgrs_purity_template_EE,"Function Fit with #sigma_{i#etai#eta} SB","pl");
  tleg2->AddEntry(tgrs_purity_template_EE,"Template Fit","pl");
//   tleg2->AddEntry(tgrs_purity_twobin_EE,"two bin","pl");
  tleg2->Draw();

  c2->SaveAs(Form("plots/%s_allfitter_purity_EE.eps",isDataName.data()));
  c2->SaveAs(Form("plots/%s_allfitter_purity_EE.gif",isDataName.data()));
  c2->SaveAs(Form("plots/%s_allfitter_purity_EE.C",isDataName.data()));


  corr = doEffCorr?"Corrected":"";
  TCanvas *c3 = new TCanvas("c3","",600,600);
  c3->Draw();
  gPad->SetLogy();
  TH2F *h3 = new TH2F("h3","",100, 0., 140, 100, 0.1, 10000);  
  h3->SetNdivisions(505,"XY");
  h3->SetXTitle("p_{T} (GeV)");
  h3->SetYTitle(Form("%s Signal yield / bin",corr.data()));
  h3->Draw();

  TGraphAsymmErrors *tg_mc_sig_EB = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, 
							  nsig_mc[0],
							  fBinsPtError,fBinsPtError,
							  nsig_err_mc[0],nsig_err_mc[0]);
  tg_mc_sig_EB->SetName("tg_mc_sig_EB");
  tg_mc_sig_EB->SetMarkerSize(1.);
  tg_mc_sig_EB->SetLineWidth(2);
  tg_mc_sig_EB->Draw("p e same");


  TGraphAsymmErrors *tgrs_sig_EB = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, nsig_func[0], 
							 fBinsPtError, fBinsPtError, 
							 nsig_err_func[0], nsig_err_func[0]);
  tgrs_sig_EB->SetName("tgrs_sig_EB");
  tgrs_sig_EB->SetMarkerColor(4);
  tgrs_sig_EB->SetMarkerStyle(25);
  tgrs_sig_EB->SetLineColor(4);
  tgrs_sig_EB->SetLineWidth(2);
  tgrs_sig_EB->Draw("p e same");

  TGraphAsymmErrors *tgrs_sig_template_EB = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, 
								  nsig_template[0], 
								  fBinsPtError,fBinsPtError, 
								  nsig_err_template[0],
								  nsig_err_template[0]);
  tgrs_sig_template_EB->SetName("tgrs_sig_template_EB");
  tgrs_sig_template_EB->SetMarkerColor(6);
  tgrs_sig_template_EB->SetMarkerStyle(26);
  tgrs_sig_template_EB->SetLineColor(6);
  tgrs_sig_template_EB->SetLineWidth(2);
  tgrs_sig_template_EB->Draw("p e same");

  TGraphAsymmErrors *tgrs_sig_twobin_EB = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, nsig_2bin[0], 
								fBinsPtError, fBinsPtError, 
								nsig_err_2bin[0],nsig_err_2bin[0]);
  tgrs_sig_twobin_EB->SetName("tgrs_sig_twobin_EB");
  tgrs_sig_twobin_EB->SetMarkerColor(kGreen);
  tgrs_sig_twobin_EB->SetMarkerStyle(3);
  tgrs_sig_twobin_EB->SetLineColor(kGreen);
  tgrs_sig_twobin_EB->SetLineWidth(2);
//   tgrs_sig_twobin_EB->Draw("p e same");

  TLegend *tleg3 = new TLegend(0.193, 0.168, 0.642, 0.437);
  tleg3->SetHeader("EB Photons");
  tleg3->SetFillColor(0);
  tleg3->SetFillStyle(0);
  tleg3->SetShadowColor(0);
  tleg3->SetBorderSize(0);
  tleg3->AddEntry(tg_mc_sig_EB,"MC Prediction","pl");
  tleg3->AddEntry(tgrs_sig_EB,"Function Fit","pl");
  tleg3->AddEntry(tgrs_sig_template_EB,"Template Fit","pl");
//   tleg3->AddEntry(tgrs_sig_template_EB,"Function Fit with #sigma_{i#etai#eta} SB","pl");
//   tleg3->AddEntry(tgrs_sig_twobin_EB,"two bin","pl");
  tleg3->Draw();

  c3->SaveAs(Form("plots/%s_allfitter_yield_EB.eps",isDataName.data()));
  c3->SaveAs(Form("plots/%s_allfitter_yield_EB.gif",isDataName.data()));
  c3->SaveAs(Form("plots/%s_allfitter_yield_EB.C",isDataName.data()));

  TCanvas *c4 = new TCanvas("c4","",600,600);
  c4->Draw();
  gPad->SetLogy();

  h3->Draw();

  TGraphAsymmErrors *tg_mc_sig_EE = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, 
							  nsig_mc[1],
							  fBinsPtError,fBinsPtError,
							  nsig_err_mc[1],nsig_err_mc[1]);
  tg_mc_sig_EE->SetName("tg_mc_sig_EE");
  tg_mc_sig_EE->SetMarkerSize(1.);
  tg_mc_sig_EE->SetLineWidth(2);
  tg_mc_sig_EE->Draw("p e same");



  TGraphAsymmErrors *tgrs_sig_EE = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, nsig_func[1], 
							 fBinsPtError,fBinsPtError, 
							 nsig_err_func[1],nsig_err_func[1]);
  tgrs_sig_EE->SetName("tgrs_sig_EE");
  tgrs_sig_EE->SetMarkerColor(4);
  tgrs_sig_EE->SetMarkerStyle(25);
  tgrs_sig_EE->SetLineColor(4);
  tgrs_sig_EE->SetLineWidth(2);
  tgrs_sig_EE->Draw("p e same");

  TGraphAsymmErrors *tgrs_sig_template_EE = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, 
								  nsig_template[1], 
								  fBinsPtError,fBinsPtError, 
								  nsig_err_template[1],
								  nsig_err_template[1]);
  tgrs_sig_template_EE->SetName("tgrs_sig_template_EE");
  tgrs_sig_template_EE->SetMarkerColor(6);
  tgrs_sig_template_EE->SetMarkerStyle(26);
  tgrs_sig_template_EE->SetLineColor(6);
  tgrs_sig_template_EE->SetLineWidth(2);
  tgrs_sig_template_EE->Draw("p e same");

  TGraphAsymmErrors *tgrs_sig_twobin_EE = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, nsig_2bin[1], 
							   fBinsPtError, fBinsPtError, 
							   nsig_err_2bin[1],nsig_err_2bin[1]);
  tgrs_sig_twobin_EE->SetName("tgrs_sig_twobin_EE");
  tgrs_sig_twobin_EE->SetMarkerColor(kGreen);
  tgrs_sig_twobin_EE->SetMarkerStyle(3);
  tgrs_sig_twobin_EE->SetLineColor(kGreen);
  tgrs_sig_twobin_EE->SetLineWidth(2);
//   tgrs_sig_twobin_EE->Draw("p e same");


  TLegend *tleg4 = new TLegend(0.193, 0.168, 0.642, 0.437);
  tleg4->SetHeader("EE Photons");
  tleg4->SetFillColor(0);
  tleg4->SetFillStyle(0);
  tleg4->SetShadowColor(0);
  tleg4->SetBorderSize(0);
  tleg4->AddEntry(tg_mc_sig_EE,"MC Prediction","pl");
  tleg4->AddEntry(tgrs_sig_EE,"Function Fit","pl");
  tleg4->AddEntry(tgrs_sig_template_EE,"Template Fit","pl");
//   tleg4->AddEntry(tgrs_sig_template_EE,"Function Fit with #sigma_{i#etai#eta} SB","pl");
//   tleg4->AddEntry(tgrs_sig_twobin_EE,"two bin","pl");
  tleg4->Draw();

  c4->SaveAs(Form("plots/%s_allfitter_yield_EE.eps",isDataName.data()));
  c4->SaveAs(Form("plots/%s_allfitter_yield_EE.gif",isDataName.data()));
  c4->SaveAs(Form("plots/%s_allfitter_yield_EE.C",isDataName.data()));

  TFile* outFile = new TFile("isotemplate_yield_rongshyang.root","recreate"); 
  tg_mc_purity_EB->Write();
  tgrs_purity_EB->Write();
  tgrs_purity_template_EB->Write();
  tgrs_purity_twobin_EB->Write();

  tg_mc_purity_EE->Write();
  tgrs_purity_EE->Write();
  tgrs_purity_template_EE->Write();
  tgrs_purity_twobin_EE->Write();

  tg_mc_sig_EB->Write();
  tgrs_sig_EB->Write();
  tgrs_sig_template_EB->Write();

  tg_mc_sig_EE->Write();
  tgrs_sig_EE->Write();
  tgrs_sig_template_EE->Write();

  outFile->Close();



}
