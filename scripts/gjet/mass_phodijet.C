#define mass_phodijet_cxx
#include "mass_phodijet.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>


void mass_phodijet::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   Long64_t nTotalGoodPho=0;
   Long64_t nTotalEventWithPho=0;   
   Long64_t nTotalGoodJet=0;
   Long64_t nTotalEventWithPho2Jet=0;

   TH1F* h_npho = new TH1F("h_npho","Number of good photons",10,-0.5,9.5);
   TH1F* h_njet = new TH1F("h_njet","Number of good jets",50,-0.5,49.5);

   TH1F* h_dR   = new TH1F("h_dR","#Delta R between leading photon and any jet", 50,0,1);
   TH1F* h_mjj  = new TH1F("h_mjj","Dijet mass", 500,0,1000);
   TH1F* h_mgjj  = new TH1F("h_mgjj","photon+Dijet mass", 500,0,1000);
   TH2F* h_mjj_mgjj  = new TH2F("h_mjj_mgjj","Dijet mass vs 3-body mass", 500,0,1000,500,0,1000);



   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      Int_t   maxPhoIndex  = -1;
      Float_t maxPhoEt     = -999.;
      Int_t   nGoodPho     = 0;

      // first checking if there's a good photon
      for(int ipho=0; ipho< nPho; ipho++){
	if(!isGoodPho(ientry, ipho))continue;
	nTotalGoodPho++;
	nGoodPho++;
	if(phoEt[ipho] > maxPhoEt)
	  {
	    maxPhoIndex = ipho;
	    maxPhoEt    = phoEt[ipho];
	  }
      } // end of loop over photons
      
      if(nGoodPho>0)
	nTotalEventWithPho++;
      h_npho->Fill(nGoodPho);

      TLorentzVector pho_p4(0,0,0,0);
      if(nGoodPho>0){
	pho_p4.SetPtEtaPhiE(phoEt[maxPhoIndex],
			    phoEta[maxPhoIndex],
			    phoPhi[maxPhoIndex],
			    phoE[maxPhoIndex]);
      }

      // second checking if there're two good jets
      
      Int_t   maxJetIndex        = -1;
      Float_t maxJetPt           = -999.;
      Int_t   secondmaxJetIndex  = -1;
      Float_t secondmaxJetPt     = -999.; 
      Int_t   nGoodJet           = 0;

      for(int ijet=0; ijet < nJet; ijet++){
	if(!isGoodJet(ientry, ijet))continue;

	TLorentzVector tempjet_p4(0,0,0,0);
	tempjet_p4.SetPtEtaPhiE(jetPt[ijet],
				jetEta[ijet],
				jetPhi[ijet],
				jetEn[ijet]);

	Float_t dR = maxPhoIndex>=0? pho_p4.DeltaR(tempjet_p4): 999.0;
	h_dR->Fill(dR);
	// remove overlap between leading photon and jets
	if(dR < 0.5)continue;

	nTotalGoodJet++;
	nGoodJet++;

	if(jetPt[ijet] > maxJetPt)
	  {
	    maxJetIndex = ijet;
	    maxJetPt    = jetPt[ijet];
	  }
	else if(jetPt[ijet] > secondmaxJetPt)
	  {
	    secondmaxJetIndex = ijet;
	    secondmaxJetPt    = jetPt[ijet];
	  }
	
      } // end of loop over jets

     
      h_njet->Fill(nGoodJet);
 
      if(maxPhoIndex < 0 || maxJetIndex < 0 || secondmaxJetIndex < 0)continue;

      nTotalEventWithPho2Jet++;

      TLorentzVector jet1_p4(0,0,0,0);
      TLorentzVector jet2_p4(0,0,0,0);


      jet1_p4.SetPtEtaPhiE(jetPt[maxJetIndex],
			   jetEta[maxJetIndex],
			   jetPhi[maxJetIndex],
			   jetEn[maxJetIndex]);


      jet2_p4.SetPtEtaPhiE(jetPt[secondmaxJetIndex],
			   jetEta[secondmaxJetIndex],
			   jetPhi[secondmaxJetIndex],
			   jetEn[secondmaxJetIndex]);

     
      Float_t mjj  = (jet1_p4+jet2_p4).M();
      Float_t mgjj = (pho_p4+jet1_p4+jet2_p4).M();

      h_mjj->Fill(mjj);
      h_mgjj->Fill(mgjj);
      h_mjj_mgjj->Fill(mjj,mgjj);
 
      
   } // end of loop over entries

   cout << "Total number of good photons = " << nTotalGoodPho << endl;
   cout << "Total number of events with good photons = " << nTotalEventWithPho << endl;
   cout << "Total number of good jets = " << nTotalGoodJet << endl;
   cout << "Total number of events with a photon and 2 jets = " << nTotalEventWithPho2Jet << endl;

  TFile* outFile = new TFile("output_histo.root","recreate");               

  h_npho->Write();
  h_njet->Write();
  h_dR->Write();
  h_mjj->Write();
  h_mgjj->Write();
  h_mjj_mgjj->Write();
  outFile->Close();      


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
  if(jetPt[ijet] < 10.0)return false;
  if(fabs(jetEta[ijet]) > 3.0)return false;
  if(jetNConstituents[ijet] <= 1)return false;
  if(jetNHF[ijet] >= 0.99)return false;
  if(jetNEF[ijet] >= 0.99)return false;

  // for the tracker region
  if(fabs(jetEta[ijet])<2.4 && jetCHF[ijet] <= 0.0)return false;
  if(fabs(jetEta[ijet])<2.4 && jetCEF[ijet] >= 0.99)return false;
  if(fabs(jetEta[ijet])<2.4 && jetNCH[ijet] <= 0)return false;

  return true;

}

