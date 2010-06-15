
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1.h>
#include <string>
#include <TCanvas.h>
#include <TLegend.h>
#include "chi2Nbins.h"
using namespace std;


const float yield[2][5]=
  {{0}};
const float yieldErr[2][5]=
  {{0}};

// the pt and eta binning
const double fBinsEta[]={0,1.45,1.7,2.5};
const double fBinsPt[]={15,20,30,50,80,120};
const int nPtBin = 5;
const int nEtaBin = 2;

void compareFittedDataMC(const int firstPtBin=1)
{
  ifstream fin;
  fin.open("rawYield.dat");
  for(int j=0; j<2; j++)
    for(int i=0; i<5; i++)
      fin >> yield[j][i] >> yieldErr[j][i];
  
  TFile* inf_data = new TFile("rawYield.root");

// 						   Form("h_EB_%s_EGdata_pt15",var));
  TH1F* hTemplate = (TH1F*)inf_data->FindObjectAny("hOutputData_EB_pt15");
  hTemplate->SetName("hTemplate");
  hTemplate->Reset();

  TH1F* hAllData[nEtaBin];
  TH1F* hAllSig[nEtaBin];
  TH1F* hAllBkg[nEtaBin];
  
  char* dec[2] = {"EB","EE"};

  for(int i=0; i<nEtaBin; i++)
    {
      hAllData[i] = (TH1F*)hTemplate->Clone();
      hAllData[i]->SetName(Form("hAllData_%s",dec[i]));
      hAllData[i]->Reset();

      hAllSig[i] = (TH1F*)hTemplate->Clone();
      hAllSig[i]->SetName(Form("hAllSig_%s",dec[i]));
      hAllSig[i]->Reset();

      hAllBkg[i] = (TH1F*)hTemplate->Clone();
      hAllBkg[i]->SetName(Form("hAllBkg_%s",dec[i]));
      hAllBkg[i]->Reset();
    }
  char tmp[1000];
  for(int ieta=0; ieta< nEtaBin; ieta++)
    {
      for(int ipt=firstPtBin; ipt< nPtBin; ipt++)
	{
// 	  TH1F* hnow = (TH1F*)inf_data->FindObjectAny(
// 						      Form("h_%s_%s_EGdata_pt%d",dec[ieta],var,(int)fBinsPt[ipt]));
	  TH1F* hnow = (TH1F*)inf_data->FindObjectAny(
						      Form("hOutputData_%s_pt%d",dec[ieta],(int)fBinsPt[ipt]));

	  hnow->SetName("hnow");
	  cout << "Here 2" << endl;
	  
	  float ndata = hnow->Integral();

	  hAllData[ieta]->Add(hnow);


	  TH1F* hnow2 = (TH1F*)inf_data->FindObjectAny(
						       Form("hOutputSig_%s_pt%d",dec[ieta],(int)fBinsPt[ipt]));
	  hnow2->SetName("hnow2");
	  
	  
	  hnow2->Sumw2();
	  hnow2->Scale(yield[ieta][ipt]/(float)hnow2->Integral());
	  hAllSig[ieta]->Add(hnow2);

	  TH1F* hnow3 = (TH1F*)inf_data->FindObjectAny(
						       Form("hOutputBkg_%s_pt%d",dec[ieta],(int)fBinsPt[ipt]));

// 	  TH1F* hnow3 = (TH1F*)inf_data->FindObjectAny(
// 						       Form("h_%s_%s_bkg_pt%d",dec[ieta],var,(int)fBinsPt[ipt]));
	  
	  hnow3->SetName("hnow3");
	  
	  hnow3->Sumw2();
	  float nbkg = ndata -yield[ieta][ipt];
	  hnow3->Scale(nbkg/(float)hnow3->Integral());
	  hnow3->Add(hnow2);
	  hAllBkg[ieta]->Add(hnow3);

	}
    }
  
  TCanvas* c1 = new TCanvas("c1","",500,500);


   for(int ieta=0; ieta<2; ieta++)
    {
      cout << "Total data yields = " << hAllData[ieta]->Integral() << endl;
      cout << "Total signal yields = " << hAllSig[ieta]->Integral() << endl;
      cout << "Total background yields = " << hAllBkg[ieta]->Integral() << endl;


      float sum   = 0;
      float sumErr= 0;
      for(int ipt=firstPtBin;ipt<5;ipt++)
	{
	  sum += yield[ieta][ipt];
	  sumErr += yieldErr[ieta][ipt]*yieldErr[ieta][ipt];
	}
      sumErr = sqrt(sumErr);
      cout << "Total signal = " << sum << " +- " << sumErr << endl;

      hAllSig[ieta]->SetLineColor(5);
      hAllSig[ieta]->SetFillColor(5);
      hAllSig[ieta]->SetFillStyle(3001);

      hAllBkg[ieta]->SetLineColor(4);
      hAllBkg[ieta]->SetLineStyle(2);
      hAllBkg[ieta]->SetLineWidth(2);

      hAllData[ieta]->SetXTitle("Iso_{ECAL}+Iso_{HCAL}+Iso_{TRK} (GeV)");
      hAllData[ieta]->SetYTitle("Entries");
      hAllData[ieta]->SetTitleOffset(1.4,"Y");
      hAllData[ieta]->SetTitle("");

      hAllBkg[ieta]->SetXTitle("Iso_{ECAL}+Iso_{HCAL}+Iso_{TRK} (GeV)");
      hAllBkg[ieta]->SetYTitle("Entries");
      hAllBkg[ieta]->SetTitleOffset(1.4,"Y");
      hAllBkg[ieta]->SetTitle("");

//        if(ieta==0)
	 hAllData[ieta]->Draw("e");
//        else
// 	 hAllBkg[ieta]->Draw("hist");
       hAllData[ieta]->Draw("e,same");
      hAllSig[ieta]->Draw("hist,same");
      hAllBkg[ieta]->Draw("hist,same");

      TLegend* leg = new TLegend(0.561,0.661,0.760,0.911);
      leg->SetHeader(Form("%s photons",dec[ieta]));
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.04);
      leg->SetBorderSize(0);
      leg->AddEntry(hAllData[ieta],"Data");
      leg->AddEntry(hAllBkg[ieta],"Fit Result");
      leg->AddEntry(hAllSig[ieta],"Estimated Signal");
      leg->Draw("same");

      double chi2ForThisBin=0;
      int nbinForThisBin=0;
      chi2NbinsHisto(hAllBkg[ieta], hAllData[ieta], 
		     chi2ForThisBin, nbinForThisBin);
      
      c1->Print(Form("PtAbove20_yield_%s.eps",dec[ieta]));
      c1->Print(Form("PtAbove20_yield_%s.pdf",dec[ieta]));
    }
 

}
