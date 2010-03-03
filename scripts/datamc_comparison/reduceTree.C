#define reduceTree_cxx
#include "reduceTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TMath.h>
#include <TVector3.h>

using namespace std;

void reduceTree::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   std::string remword=".root";
   size_t pos = inputFile_.find(remword);
  
   if(pos!= std::string::npos)
     inputFile_.swap(inputFile_.erase(pos,remword.length()));
   
   // dump histogram to a root file
   std::string outputFile_ = inputFile_ + "_small.root";

 
   TFile* newfile = new TFile(outputFile_.data(),"recreate");
   TTree* newtree = fChain->CloneTree(0);

   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      bool pass = true;
      // first loose photon ID cuts
      for (int ph=0; ph < nPhotons; ph++){

	bool isBarrel = (isEB[ph]==1);
	bool isEndCap = (isEE[ph]==1);
	bool inAnyGap = (isEBEEGap[ph]==1) || 
	  (isEB[ph]==1 && isEBGap[ph]==1) || 
	  (isEE[ph]==1 && isEEGap[ph]==1);

	if(inAnyGap)continue;

	float thiseta = caloPositionEta[ph];	
	if(fabs(thiseta) > 2.5)continue;


	float ecalIso = ecalRecHitSumEtConeDR04[ph];
	float hcalIso = hcalTowerSumEtConeDR04[ph];
	float trkIso  = trkSumPtHollowConeDR04[ph];
	float hadem   = hadronicOverEm[ph];
	float thiset = et[ph];

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
	
	pass = true;
	break;
      } // end of loop over photons
      if(pass)newtree->Fill();

   } // end of loop over entries

   cout << "Total saved" << endl;
   newtree->Print();
   newtree->AutoSave();
   delete newfile;
}
