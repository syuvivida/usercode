#define mass_phodijet_cxx
#include "mass_phodijet.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void mass_phodijet::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   Long64_t nGoodPho=0;
   Long64_t nEventWithPho=0;

   Long64_t nGoodJet=0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      Int_t maxPhoIndex  = -1;
      Float_t maxPhoEt   = -999;

      // first checking if there's a good photon
      for(int ipho=0; ipho< nPho; ipho++){
	if(!isGoodPho(ientry, ipho))continue;
	nGoodPho++;
	if(phoEt[ipho] > maxPhoEt)
	  {
	    maxPhoEt    = phoEt[ipho];
	    maxPhoIndex = ipho;
	  }
      } // end of loop over photons
      
      if(maxPhoIndex > -1)
	nEventWithPho++;


      // second checking if there's a good jet
      for(int ijet=0; ijet < nJet; ijet++){
	if(!isGoodJet(ientry, ijet))continue;
	nGoodJet++;
      } // end of loop over jets


   } // end of loop over entries

   cout << "Total number of good photons = " << nGoodPho << endl;
   cout << "Total number of events with good photons = " << nEventWithPho << endl;
   cout << "Total number of good jets = " << nGoodJet << endl;

} // end of loop function


// check if this reco-photon is a good photon
Bool_t mass_phodijet::isGoodPho(Long64_t entry, Int_t ipho)
{


  bool isEB=false;
  bool isEE=false;
  if(fabs(phoSCEta[ipho]) < 1.4442)isEB=true;
  if(fabs(phoSCEta[ipho]) > 1.566 && 
     fabs(phoSCEta[ipho]) < 2.5) isEE=true;

  if(!isEB && !isEE)return false;
  if(phoEt[ipho] < 20)return false;
  if(phoHoverE[ipho] > 0.05)return false;
  if(phohasPixelSeed[ipho]==1)return false;
  if(phoEcalIsoDR04[ipho] > 4.2 +0.006 * phoEt[ipho])return false;
  if(phoHcalIsoDR04[ipho] > 2.2 +0.0025* phoEt[ipho])return false;
  if(phoTrkIsoHollowDR04[ipho] > 2.0 +0.001* phoEt[ipho])return false;
  if(isEB && phoSigmaIEtaIEta[ipho] > 0.01)return false;
  if(isEE && phoSigmaIEtaIEta[ipho] > 0.028)return false;

  return true;

}

// check if this reco-jet is a good jet
Bool_t mass_phodijet::isGoodJet(Long64_t entry, Int_t ijet)
{
  if(jetPt[ijet] < 20.0)return false;
  if(fabs(jetEta[ijet]) > 3.0)return false;
  if(jetNConstituents[ijet] <= 1)return false;
  if(jetCHF[ijet] <= 0.0)return false;
  if(jetCEF[ijet] >= 0.99)return false;
  if(jetNCH[ijet] <= 0)return false;
  if(jetNHF[ijet] >= 0.99)return false;
  if(jetNEF[ijet] >= 0.99)return false;
  
  return true;

}

