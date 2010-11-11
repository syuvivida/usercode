#include <fstream>
#include <iostream>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TSystem.h>

using namespace std;

void readPullBias(char*prefix="pull_bias"){

  gSystem->mkdir(Form("figures_%s",prefix));
  char* dec[2] = {"EB","EE"};
  const int fBinsPt[]={21,23,26,30,35,40,45,50,60,85,120,300};
  const int nPtBin = sizeof(fBinsPt)/sizeof(fBinsPt[0])-1;
  const int NRUNS=8;

  ofstream foutm;
  foutm.open("mean.dat");

  ofstream fouts;
  fouts.open("sigma.dat");

  TH1D* hpull = new TH1D("hpull","",50,-5.0,5.0);
  hpull->SetXTitle("(Fit-Input)/Fit_err");
  TH1D* hbias = new TH1D("hbias","",100,-0.5,0.5);
  hbias->SetXTitle("(Fit-Input)/Input");

  TCanvas* c1 = new TCanvas("c1","",500,500);

  char tmp[1000];
  for(int ieta=0; ieta < 2; ieta++){
    for(int ipt=0; ipt < nPtBin; ipt++){
      hpull->Reset();
      hbias->Reset();
      for(int irun=0; irun < NRUNS; irun++){
	sprintf(tmp,"%s%d_%s_pt%d.dat",prefix,irun,dec[ieta],fBinsPt[ipt]);		cout << "Opening file " << tmp << endl;
	FILE *infile =  fopen(tmp,"r");  
	cout << "Opened file " << tmp << endl;

	Double_t pull=-100., bias=-100.; // combined isolation, pt, eta
	int flag = 1;
	if(infile==NULL){
	  cout << "File does not exist" << endl;
	  continue;
	}
	  
	while (flag!=-1){
	  flag =fscanf(infile,"%lf %lf",&pull, &bias);
	  if(flag!=-1){
	    hpull->Fill(pull);
	    hbias->Fill(bias);
	  }
 
	}// keep reading files as long as text exists

      } // end of iterations
      
      cout << hpull->GetEntries() << endl;
      cout << hbias->GetEntries() << endl;
 
      hpull->SetTitle(Form("%s, pt = %d-%d GeV",dec[ieta],
			   fBinsPt[ipt],fBinsPt[ipt+1]));
      
      hbias->SetTitle(Form("%s, pt = %d-%d GeV",dec[ieta],
			   fBinsPt[ipt],fBinsPt[ipt+1]));
      if(hpull->GetEntries()<1)continue;
      hpull->Fit("gaus");
      TF1* f1 = hpull->GetFunction("gaus");
      foutm << dec[ieta] << " pt" << fBinsPt[ipt] << " pull: " << f1->GetParameter(1) << " +- " << f1->GetParError(1) << endl;
      fouts << dec[ieta] << " pt" << fBinsPt[ipt] << " pull: " << f1->GetParameter(2) << " +- " << f1->GetParError(2) << endl;

      c1->Print(Form("figures_%s/pull_%s_%s_pt%d.png",
		     prefix,prefix,dec[ieta],fBinsPt[ipt]));

      if(hbias->GetEntries()<1)continue;
      hbias->Fit("gaus");
      TF1* f2 = hbias->GetFunction("gaus");
      foutm << dec[ieta] << " pt" << fBinsPt[ipt] << " bias: " << f2->GetParameter(1) << " +- " << f2->GetParError(1) << endl;
      fouts << dec[ieta] << " pt" << fBinsPt[ipt] << " bias: " << f2->GetParameter(2) << " +- " << f2->GetParError(2) << endl;
      c1->Print(Form("figures_%s/bias_%s_%s_pt%d.png",
		     prefix,prefix,dec[ieta],fBinsPt[ipt]));

    } // end of loop over pt bins
  } // end of loop over eta bins



}
