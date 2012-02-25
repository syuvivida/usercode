#define debugCount_cxx
#include "debugCount.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TProfile.h>
#include <TLorentzVector.h>


const Float_t BARREL_MAXETA=1.4442;
const Float_t ENDCAP_MINETA=1.566;
const Float_t ENDCAP_MAXETA=2.5;

const Int_t   NETA=2;
const Float_t fEtaBin[]={0,ENDCAP_MINETA,ENDCAP_MAXETA};

void debugCount::Loop() 
{

  cout << "here version 2" << endl;
   if (fChain == 0) return;

   cout << "there" << endl;
   Long64_t nentries = fChain->GetEntriesFast();

   cout << "There are " << nentries << " entries" << endl; 
   // declare and define histograms
   TH1F* h_dec = new TH1F("h_dec","",2,fEtaBin);
   TH1F* h_dR  = new TH1F("h_dR","",100,0,5);

   // pt distribution
   TH1F* h_pt_template = new TH1F("h_pt_template","",500,0,500);
   TH1F* h_pt_dirgamma = (TH1F*) h_pt_template->Clone("h_pt_dirgamma");
   TH1F* h_pt_1stjet = (TH1F*) h_pt_template->Clone("h_pt_1stjet");
   TH1F* h_pt_2ndjet = (TH1F*) h_pt_template->Clone("h_pt_2ndjet");
   
   // rapidity distribution
   TH1F* h_y_template = new TH1F("h_y_template","",100,-5,5);
   TH1F* h_y_dirgamma = (TH1F*) h_y_template->Clone("h_y_dirgamma");
   TH1F* h_y_1stjet = (TH1F*) h_y_template->Clone("h_y_1stjet");
   TH1F* h_y_2ndjet = (TH1F*) h_y_template->Clone("h_y_2ndjet");


   Long64_t nPass[30]={0};

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //      if(jentry>50000)break;
      nPass[0] ++;

      Int_t TriggerHLT75 = HLTIndex[30];

      if(TriggerHLT75    > 0 && HLT[TriggerHLT75]>0)      
      nPass[1]++;
      else continue;
      


      Int_t ngood_vtx=IsVtxGood;
      if(ngood_vtx==0)continue;

      nPass[2] ++;

      bool findLeadingPhoton = false;
      int leadingPhotonIndex = -1;
      Float_t phoMaxPt = -9999.;

     
      // now find a good leading photon
      for(int ipho=0; ipho < nPho; ipho++){
	
	if(!isGoodPho(ientry,ipho))continue;


	h_pt_dirgamma->Fill(phoEt[ipho]);
	h_y_dirgamma ->Fill(phoSCEta[ipho]);

	if(phoEt[ipho] > phoMaxPt)
	  {
	    phoMaxPt = phoEt[ipho];
	    leadingPhotonIndex= ipho;
	    findLeadingPhoton = true;
	  }

	
      } // end of leading photon search

      if(!findLeadingPhoton)continue;
      nPass[3]++;

      TLorentzVector l4_pho(0,0,0,0);
      l4_pho.SetPtEtaPhiE(
			  phoEt[leadingPhotonIndex],
			  phoEta[leadingPhotonIndex],
			  phoPhi[leadingPhotonIndex],
			  phoE[leadingPhotonIndex]
			  );

      Float_t leadingPhotonEt = phoEt[leadingPhotonIndex];

      // first check which reco jet is the one from the highest and 
      // second et gen jet
      bool findLeadingJet = false;
      Int_t leadingJetIndex=-1;
      Float_t jetMaxPt=-9999.;

      bool findSecondLeadingJet = false;
      Int_t secondLeadingJetIndex=-1;
      Float_t secondJetMaxPt=-9999.;

      for(int ijet=0; ijet < nJet; ijet++){

	TLorentzVector l4_tempjet(0,0,0,0);
	l4_tempjet.SetPtEtaPhiE(
				jetPt[ijet],
				jetEta[ijet],
				jetPhi[ijet],
				jetEn[ijet]
				);

	Float_t dR = l4_tempjet.DeltaR(l4_pho);
	h_dR->Fill(dR);

	if(dR<0.5)continue;
      
	if(!isGoodLooseJet(ientry,ijet))continue;

	if(jetPt[ijet] > jetMaxPt)
	  {
	    secondJetMaxPt = jetMaxPt;
	    secondLeadingJetIndex = leadingJetIndex;
	    jetMaxPt = jetPt[ijet];
	    leadingJetIndex= ijet;
	  }
	else if(jetPt[ijet] > secondJetMaxPt)
	  {
	    secondJetMaxPt = jetPt[ijet];
	    secondLeadingJetIndex = ijet;
	  }


      } // end of loop over jets


      // only study
      if(leadingJetIndex >=0 
	 )
	findLeadingJet = true;
	 
      if(secondLeadingJetIndex >=0 
	 )
	findSecondLeadingJet = true;


      TLorentzVector l4_1stjet(0,0,0,0);
      if(findLeadingJet)
	{
	  nPass[4]++;
	  l4_1stjet.SetPtEtaPhiE(
				  jetPt[leadingJetIndex],
				  jetEta[leadingJetIndex],
				  jetPhi[leadingJetIndex],
				  jetEn[leadingJetIndex]
				  );
	  h_y_1stjet->Fill(l4_1stjet.Rapidity());
 	  h_pt_1stjet->Fill(l4_1stjet.Pt());
	  
	}
      

      

      TLorentzVector l4_2ndjet(0,0,0,0);           
      if(findSecondLeadingJet)
	{
	  l4_2ndjet.SetPtEtaPhiE(
				  jetPt[secondLeadingJetIndex],
				  jetEta[secondLeadingJetIndex],
				  jetPhi[secondLeadingJetIndex],
				  jetEn[secondLeadingJetIndex]
				  );

	   h_y_2ndjet->Fill(l4_2ndjet.Rapidity());
	   h_pt_2ndjet->Fill(l4_2ndjet.Pt());
	}


      if(!findLeadingPhoton || !findLeadingJet)continue;
      nPass[5]++;


      if(leadingPhotonEt >= 85 && leadingPhotonEt < 95)nPass[7]++;

      if(leadingPhotonEt >= 95 && leadingPhotonEt < 110)nPass[8]++;
      
      if(leadingPhotonEt >= 110 && leadingPhotonEt < 130)nPass[9]++;

      if(leadingPhotonEt >= 130 && leadingPhotonEt < 160)nPass[10]++;

      if(leadingPhotonEt >= 160 && leadingPhotonEt < 200)nPass[11]++;
      

  
   } // end of loop over entries


   for(int i=0; i<30; i++)
     if(nPass[i]>0)
       cout << "nPass["<< i << "]=" << nPass[i] << endl;
   
   TFile* outFile = new TFile("/home/syu/CVSCode/test.root",
			      "recreate");
   h_y_dirgamma -> Write();
   h_y_1stjet -> Write();
   h_y_2ndjet -> Write();

   h_pt_dirgamma -> Write();
   h_pt_1stjet -> Write();
   h_pt_2ndjet -> Write();

   h_dR -> Write();


   outFile->Close();     

}


Bool_t debugCount::isFidJet (Long64_t entry, Int_t ijet)
{
  if(jetPt[ijet] < 10.0)return false;
  if(fabs(jetEta[ijet]) > 3.0)return false;
  return true;

}


// check if this reco-jet is a good loose jet
Bool_t debugCount::isGoodLooseJet(Long64_t entry, Int_t ijet)
{
  if(fabs(jetEta[ijet]) > 2.4)return false;
  if(jetPt[ijet] < 30.0)return false;
  if(jetNConstituents[ijet] <= 1)return false;
  if(jetNHF[ijet] >= 0.99)return false;
  if(jetNEF[ijet] >= 0.99)return false;

  // for the tracker region
  if(fabs(jetEta[ijet])<2.4 && jetCHF[ijet] <= 0.0)return false;
  if(fabs(jetEta[ijet])<2.4 && jetCEF[ijet] >= 0.99)return false;
  if(fabs(jetEta[ijet])<2.4 && jetNCH[ijet] <= 0)return false;

  return true;

}


// check if this reco-jet is a good medium jet
Bool_t debugCount::isGoodMediumJet(Long64_t entry, Int_t ijet)
{
  if(jetPt[ijet] < 10.0)return false;
  if(fabs(jetEta[ijet]) > 3.0)return false;
  if(jetNConstituents[ijet] <= 1)return false;
  if(jetNHF[ijet] >= 0.95)return false;
  if(jetNEF[ijet] >= 0.95)return false;

  // for the tracker region
  if(fabs(jetEta[ijet])<2.4 && jetCHF[ijet] <= 0.0)return false;
  if(fabs(jetEta[ijet])<2.4 && jetCEF[ijet] >= 0.99)return false;
  if(fabs(jetEta[ijet])<2.4 && jetNCH[ijet] <= 0)return false;

  return true;

}

// check if this reco-jet is a good medium jet
Bool_t debugCount::isGoodTightJet(Long64_t entry, Int_t ijet)
{
  if(jetPt[ijet] < 10.0)return false;
  if(fabs(jetEta[ijet]) > 3.0)return false;
  if(jetNConstituents[ijet] <= 1)return false;
  if(jetNHF[ijet] >= 0.90)return false;
  if(jetNEF[ijet] >= 0.90)return false;

  // for the tracker region
  if(fabs(jetEta[ijet])<2.4 && jetCHF[ijet] <= 0.0)return false;
  if(fabs(jetEta[ijet])<2.4 && jetCEF[ijet] >= 0.99)return false;
  if(fabs(jetEta[ijet])<2.4 && jetNCH[ijet] <= 0)return false;

  return true;

}

Bool_t debugCount::isGoodPho(Long64_t entry, Int_t ipho)
{


  if(fabs(phoSCEta[ipho]) > 1.44)return false;
  if(phoHoverE[ipho] > 0.05)return false;
  if(phohasPixelSeed[ipho]==1)return false;

  if(phoEcalIsoDR04[ipho] > 4.2 +0.003 * phoEt[ipho])return false;
  if(phoHcalIsoDR04[ipho] > 2.2 +0.001* phoEt[ipho])return false;
  if(phoTrkIsoHollowDR04[ipho] > 2.0 +0.001* phoEt[ipho])return false;
 

  return true;

}

