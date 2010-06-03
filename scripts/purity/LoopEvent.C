#define LoopEvent_cxx
#include "LoopEvent.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>
#include <iostream>

using namespace std;

void LoopEvent::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   ofstream fout;
   fout.open("EGdata_100603.dat");
   
   
   Long64_t ntotal = 0;
   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry%1000==0)cout << jentry << endl;
      // if (Cut(ientry) < 0) continue;

      if( TTBit[36] || TTBit[37] || TTBit[38] || TTBit[39])continue;
      if(vtxIsFake)continue;
      if(!(vtxNdof > 4 && fabs(vtxZ) <= 15))continue;
      if(!HLT_Photon15_L1R)continue;

      for(int i=0; i < nPhotons; i++){
	
	if(isEB[i] && !((seedSeverity[i]!=3 && seedSeverity[i]!=4 ) && (seedRecoFlag[i] != 2)))
	  continue;

	Float_t sumIso = ecalRecHitSumEtConeDR04[i]+hcalTowerSumEtConeDR04[i]+trkSumPtHollowConeDR04[i];
	
	if(!(pt[i] > 15.0 && fabs(eta[i])<1.45 && hadronicOverEm[i] <  0.05 
	     && sigmaIetaIeta[i] < 0.01 && sumIso < 11.0)
	   && !(pt[i] > 15.0 && fabs(eta[i]) > 1.7 && fabs(eta[i]) < 2.5 
		&& hadronicOverEm[i] <  0.05 && fabs(ESRatio[i]) > 0.1 && sumIso < 11.0))continue;
	
	fout << sumIso << " " << pt[i] << " " << eta[i] << endl;
	ntotal++;
      }


   }

   cout << "ntotal = " << ntotal << endl;
   fout.close();


}
