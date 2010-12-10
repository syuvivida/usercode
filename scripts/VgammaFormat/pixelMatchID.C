#define pixelMatchID_cxx
#include "pixelMatchID.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>

void pixelMatchID::Loop()
{
  if (fChain == 0) return;


  TGraphAsymmErrors *phoEff_Pt[2];
  char* dec[2] = {"EB","EE"};

  TH1F* hpt_denominator[2];
  TH1F* hpt_numerator[2];

  for (int a=0; a<2; a++) 
    {
      phoEff_Pt[a] = new TGraphAsymmErrors();
      phoEff_Pt[a] -> SetNameTitle(Form("Efficiency_%s",dec[a]),Form("%s",dec[a]));
      hpt_denominator[a]       = new TH1F(Form("hpt_denominator%d",a),"",60,0,300);
      hpt_numerator[a]       = new TH1F(Form("hpt_numerator%d",a),"",60,0,300);
    }

  TGraphAsymmErrors *phoEff_Eta;
  phoEff_Eta = new TGraphAsymmErrors();
  phoEff_Eta       -> SetNameTitle("Efficiency_Eta","");
  TH1F* heta_denominator    = new TH1F("heta_denominator","",50,-2.5,2.5);
  TH1F* heta_numerator      = new TH1F("heta_numerator","",50,-2.5,2.5);

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    if (jentry%100==0)
      printf("%4.1f%% started.\r",(float)jentry/(float)nentries*100.);                                                                       
     
    for(int iph=0; iph < nPho; iph++){
      int decIndex=-1;
      if(fabs(phoSCEta[iph])<1.4442)decIndex=0;
      else if(fabs(phoSCEta[iph])>1.566 && fabs(phoSCEta[iph])< 2.5)decIndex=1;
      else continue;

//       if(phoEt[iph]<21)continue;
      if(phoGenIndex[iph]<0)continue;
      if(phoGenMomPID[iph]!=22)continue;
      if(phoHoverE[iph]>0.05)continue;
      if(phoSigmaIEtaIEta[iph]>0.01 && decIndex==0)continue;
      if(phoSigmaIEtaIEta[iph]>0.028 && decIndex==1)continue;

      hpt_denominator[decIndex]->Fill(phoEt[iph]);
      heta_denominator->Fill(phoEta[iph]);
	
      if(!phohasPixelSeed[iph])
	{
	  hpt_numerator[decIndex]->Fill(phoEt[iph]);
	  heta_numerator->Fill(phoEta[iph]);
	}

    } // end of looping over photons

  } // end of looping over entries


   
  for(int a=0; a<2; a++)
    phoEff_Pt[a]->BayesDivide(hpt_numerator[a],hpt_denominator[a]);
  phoEff_Eta->BayesDivide(heta_numerator, heta_denominator);
 
  std::string remword=".root";
   size_t pos = inputFile_.find(remword);
   std::string forOutput = inputFile_;  
   if(pos!= std::string::npos)
     forOutput.swap(forOutput.erase(pos,remword.length()));   
   std::string endfix = "_efficiency.root";
   std::string outputFile = forOutput + endfix;

   TFile *file = new TFile(outputFile.data(), "UPDATE");
  for (int a=0; a<2; a++) {
    phoEff_Pt[a]->Write();
  }
  phoEff_Eta->Write();
  file->Close();



}
