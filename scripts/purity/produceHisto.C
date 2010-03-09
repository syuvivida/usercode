// Author: Shin-Shan Eiko Yu, National Central University, Taiwan

#define produceHisto_cxx
#include "produceHisto.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TMath.h>
#include <TVector3.h>

using namespace std;

void produceHisto::Loop()
{

  cout << "The settings are: " << endl;
  cout << "useLeadingPhotonOnly_ = " << useLeadingPhotonOnly_ << endl;
  cout << "matching_ = " << matching_ << endl;

  // check if this is a QCD or a photon+jet samples
  // should check the particle list in GenpBlock when they are available
  int isQCD = 0;
  if(inputFile_.find("QCD") != std::string::npos)isQCD=1;

  // check the pthat sample
  const int pthat[NMC+1]={15, 30, 80, 170, 300, 470, 800, 10000000};
  char name[300];
  Float_t PTHAT_MIN = -999.;
  Float_t PTHAT_MAX = -999.;
  int thisMC = -1;
  for(int i=0; i < NMC; i++)
    {
      sprintf(name,"%d",pthat[i]);
      if(inputFile_.find(name) != std::string::npos){
	PTHAT_MIN = pthat[i];
	PTHAT_MAX = pthat[i+1];
	thisMC = i;
      }
    }
 

  cout << "Restricted to pthat = " << PTHAT_MIN << " to " << PTHAT_MAX 
       << " GeV" << endl;


   if (fChain == 0) return;

   // First, define histograms
   // pt binning
   Float_t fBinsY[]={15,20,27,35,45,57,72,90,120,150,200,300,400,550};
   const int nYbin = sizeof(fBinsY)/sizeof(fBinsY[0])-1;
   Float_t fBinsZ[]={0,0.9,1.45,1.55,2.5};
   const int nZbin = sizeof(fBinsZ)/sizeof(fBinsZ[0])-1;
   TH1F* hecaliso[nYbin][nZbin]; 
   TH1F* hhcaliso[nYbin][nZbin]; 
   TH1F* hcaliso[nYbin][nZbin];
   TH1F* htrkiso[nYbin][nZbin];  

   TH1F* h_pthat = new TH1F("h_pthat","",1000,0,1000);

   char title[300];
   
   for(int i=0; i < nYbin; i++){     
     for(int j=0; j < nZbin; j++){


       sprintf(name,"hcaliso%02i_%02i",i,j);
       sprintf(title,"Calorimeter isolation pt = %3.0f--%3.0f GeV,"
	       " %1.2f < |eta| < %1.2f",
	       fBinsY[i],fBinsY[i+1], fBinsZ[j], fBinsZ[j+1]);
       hcaliso[i][j] = new TH1F(name,title,280,-2.,12.0);
       hcaliso[i][j] -> Sumw2();


       sprintf(name,"hecaliso%02i_%02i",i,j);
       sprintf(title,"ECAL isolation pt = %3.0f--%3.0f GeV,"
	       " %1.2f < |eta| < %1.2f",
	       fBinsY[i],fBinsY[i+1], fBinsZ[j], fBinsZ[j+1]);
       hecaliso[i][j] = new TH1F(name,title,180,-2.,7.0);
       hecaliso[i][j] -> Sumw2();


       sprintf(name,"hhcaliso%02i_%02i",i,j);
       sprintf(title,"HCAL isolation pt = %3.0f--%3.0f GeV,"
	       " %1.2f < |eta| < %1.2f",
	       fBinsY[i],fBinsY[i+1], fBinsZ[j], fBinsZ[j+1]);
       hhcaliso[i][j] = new TH1F(name,title,100,0,5.0);
       hhcaliso[i][j] -> Sumw2();


       sprintf(name,"htrkiso%02i_%02i",i,j);
       sprintf(title,"Track isolation pt = %3.0f--%3.0f GeV,"
	       " %1.2f < |eta| < %1.2f",
	       fBinsY[i],fBinsY[i+1], fBinsZ[j], fBinsZ[j+1]);
       htrkiso[i][j] = new TH1F(name,title,90,0,9.0);
       htrkiso[i][j] -> Sumw2();
     }

   }


   // Now loop over entries
   Long64_t ncount =0;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      // cut away overlapping pthat

      if(ptHat < PTHAT_MIN || ptHat > PTHAT_MAX)continue;
      h_pthat->Fill(ptHat);

      // first loose photon ID cuts
      for (int ph=0; ph < nPhotons; ph++){

      // There is a bad generator-level matching requirements now
      // This has to be updated

	bool hasGenPho=false;
	bool hasHardPho=false;
	for(int nmatch=0; nmatch < nGenMatch[ph] && nmatch < 9; nmatch++)
	  {
	    if(genMatchPdgId[ph][nmatch]==22 
	       && genMatchStatus[ph][nmatch]==1)
	      hasGenPho=true;
	    if(genMatchPdgId[ph][nmatch]==22 
	       && genMatchStatus[ph][nmatch]==3)
	      hasHardPho=true;
	  }

  	if( ( ((!hasHardPho || !hasGenPho) && isQCD==0) || 
	     (!hasGenPho && isQCD==1) ) && matching_ )continue;

	bool isBarrel = (isEB[ph]==1);
	bool isEndCap = (isEE[ph]==1);
	bool inAnyGap = (isEBEEGap[ph]==1) || 
	  (isEB[ph]==1 && isEBGap[ph]==1) || 
	  (isEE[ph]==1 && isEEGap[ph]==1);
	Float_t thiset = et[ph];
	TVector3 temp(momentumX[ph],
		      momentumY[ph],
		      momentumZ[ph]
		      );
	//	Float_t thiseta = temp.Eta();
	Float_t thiseta     = caloPositionEta[ph];	
	Float_t ecalIso     = ecalRecHitSumEtConeDR04[ph];
	Float_t hcalIso     = hcalTowerSumEtConeDR04[ph];
	Float_t calIso      = hcalIso + ecalIso;
	Float_t trkIso      = trkSumPtHollowConeDR04[ph];
	Float_t trkIsoDR03  = trkSumPtHollowConeDR03[ph];
	Float_t hadem       = hadronicOverEm[ph];

	
// 	cout << jentry << " : " << ph << "\t" 
// 	     << inAnyGap << "\t" << thiset << "\t" << thiseta << "\t" 
// 	     << ecalIso  << "\t" << hcalIso << "\t" << trkIso << "\t" 
// 	     << hadem << endl;

	if(inAnyGap)continue;

	// basic Et and supercluster eta requirements
 	if(thiset        < fBinsY[0])continue;
 	if(thiset        >= fBinsY[nYbin])continue;
 	if(fabs(thiseta) < fBinsZ[0])continue;
 	if(fabs(thiseta) >= fBinsZ[nZbin])continue;
 	if(fabs(thiseta) >= 2.5)continue;

	// Barrel cuts
 	if(isBarrel && ecalIso  > 5.0 + 0.004*thiset)continue;
 	if(isBarrel && hcalIso  > 5.0)continue;
 	if(isBarrel && trkIso > 9.0 )continue;
	if(isBarrel && hadem > 0.15)continue;

	// Endcap cuts
 	if(isEndCap && ecalIso  > 5.0 + 0.0021*thiset)continue;
 	if(isEndCap && hcalIso  > 5.0)continue;
 	if(isEndCap && trkIso > 9.0 )continue;
	if(isEndCap && hadem > 0.15)continue;

	int EtBin = -1;
	EtBin = TMath::BinarySearch(nYbin+1, fBinsY,thiset);
	
	int EtaBin = -1;
	EtaBin = TMath::BinarySearch(nZbin+1,fBinsZ,fabs(thiseta));

//  	cout << "EtBin = " << EtBin << " and EtaBin = " << EtaBin << endl;

	hecaliso[EtBin][EtaBin]->Fill(ecalIso);
	hhcaliso[EtBin][EtaBin]->Fill(hcalIso);
	htrkiso[EtBin][EtaBin]->Fill(trkIso);

	hcaliso[EtBin][EtaBin]->Fill(calIso);

	if(useLeadingPhotonOnly_)break;

      } // end of loop over photons

   } // end of loop over entries

   std::string remword=".root";
   size_t pos = inputFile_.find(remword);
  
   if(pos!= std::string::npos)
     inputFile_.swap(inputFile_.erase(pos,remword.length()));
   

   // dump histogram to a root file
   std::string histoFile = inputFile_ + "_histo.root";

   TFile* outFile = new TFile(histoFile.data(),"recreate");
   for(int i=0; i < nYbin; i++){     
     for(int j=0; j < nZbin; j++){

        hecaliso[i][j] ->Write();
        hhcaliso[i][j] ->Write();
        htrkiso[i][j]  ->Write();
	hcaliso[i][j]  ->Write();
      }

    }
   
   h_pthat->Write();
   outFile->Close();
   cout << "ncount = " << ncount << endl;

 }
