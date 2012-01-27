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

      Int_t HLT_Photon30_CaloIdVL = HLTIndex[63];
      
      Int_t HLT_Photon30_CaloIdVL_IsoL = HLTIndex[64];

      Int_t HLT_Photon50_CaloIdVL = HLTIndex[65];

      Int_t HLT_Photon50_CaloIdVL_IsoL = HLTIndex[66];

      Int_t HLT_Photon75_CaloIdVL = HLTIndex[30];

      Int_t HLT_Photon75_CaloIdVL_IsoL = HLTIndex[29];

      Int_t HLT_Photon90_CaloIdVL = HLTIndex[67];

      Int_t HLT_Photon90_CaloIdVL_IsoL = HLTIndex[15];

      Int_t HLT_Photon135 = HLTIndex[31];
      
      if(HLT_Photon75_CaloIdVL      > 0 && HLT[HLT_Photon75_CaloIdVL]>0)
	nPass[1]++;
      if(HLT_Photon75_CaloIdVL_IsoL > 0 && HLT[HLT_Photon75_CaloIdVL_IsoL]>0)
        nPass[15]++;

      if(HLT_Photon30_CaloIdVL      > 0 && HLT[HLT_Photon30_CaloIdVL]>0)
	nPass[16]++;
      if(HLT_Photon30_CaloIdVL_IsoL > 0 && HLT[HLT_Photon30_CaloIdVL_IsoL]>0)
	nPass[17]++;

      if(HLT_Photon50_CaloIdVL      > 0 && HLT[HLT_Photon50_CaloIdVL]>0)
        nPass[18]++;
      if(HLT_Photon50_CaloIdVL_IsoL > 0 && HLT[HLT_Photon50_CaloIdVL_IsoL]>0)
	nPass[19]++;

      if(HLT_Photon90_CaloIdVL      > 0 && HLT[HLT_Photon90_CaloIdVL]>0)
        nPass[20]++;
      if(HLT_Photon90_CaloIdVL_IsoL > 0 && HLT[HLT_Photon90_CaloIdVL_IsoL]>0)
	nPass[21]++;

      if(HLT_Photon135              > 0 && HLT[HLT_Photon135]>0)
        nPass[22]++;


      if(HLT_Photon75_CaloIdVL      > 0 && HLT[HLT_Photon75_CaloIdVL]>1)
	nPass[4]++; // prescaled

      if(HLT_Photon75_CaloIdVL      > 0 && HLT[HLT_Photon75_CaloIdVL]!=1)continue;
      nPass[3]++;
      



      Int_t ngood_vtx=IsVtxGood;
      if(ngood_vtx==0)continue;

      nPass[5] ++;



      bool findLeadingPhoton = false;
      int leadingPhotonIndex = -1;
      Float_t phoMaxPt = -9999.;

      // now find a good leading photon
      for(int ipho=0; ipho < nPho; ipho++){
	
	if(!isGoodPho(ientry,ipho))continue;

	nPass[6]++;

	h_pt_dirgamma->Fill(phoEt[ipho]);
	h_y_dirgamma ->Fill(phoSCEta[ipho]);

	if( fabs(phoSCEta[ipho])<0.9 )nPass[7]++;
	else if( fabs(phoSCEta[ipho])>0.9 && 
		 fabs(phoSCEta[ipho])<BARREL_MAXETA )
	  nPass[8]++;
	else if( fabs(phoSCEta[ipho])>ENDCAP_MINETA &&
		 fabs(phoSCEta[ipho])< 2.1)
	  nPass[9]++;

	else if( fabs(phoSCEta[ipho])> 2.1 && 
		 fabs(phoSCEta[ipho])< ENDCAP_MAXETA)
	  nPass[10]++;
	
	if(fabs(phoSCEta[ipho]) < BARREL_MAXETA && phoSigmaIEtaIEta[ipho] < 0.010)
	  nPass[13]++;

	if(fabs(phoSCEta[ipho]) < ENDCAP_MAXETA && 
	   fabs(phoSCEta[ipho]) > ENDCAP_MINETA 
	   && phoSigmaIEtaIEta[ipho] < 0.030)
	  nPass[14]++;

	


	if(phoEt[ipho] > phoMaxPt)
	  {
	    phoMaxPt = phoEt[ipho];
	    leadingPhotonIndex= ipho;
	    findLeadingPhoton = true;
	  }
	
      } // end of leading photon search

      if(!findLeadingPhoton)continue;
	
      nPass[11]++;

      // first check which reco jet is the one from the highest and 
      // second et gen jet
      bool findLeadingJet = false;
      Int_t leadingJetIndex=-1;
      Float_t jetMaxPt=-9999.;

      bool findSecondLeadingJet = false;
      Int_t secondLeadingJetIndex=-1;
      Float_t secondJetMaxPt=-9999.;

      for(int ijet=0; ijet < nJet; ijet++){

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


      TLorentzVector gen_1stjet(0,0,0,0);
      if(findLeadingJet)
	{
	  nPass[12]++;
	  gen_1stjet.SetPtEtaPhiE(
				  jetPt[leadingJetIndex],
				  jetEta[leadingJetIndex],
				  jetPhi[leadingJetIndex],
				  jetEn[leadingJetIndex]
				  );
	  h_y_1stjet->Fill(gen_1stjet.Rapidity());
 	  h_pt_1stjet->Fill(gen_1stjet.Pt());
	  
	}
      

      

      TLorentzVector gen_2ndjet(0,0,0,0);           
      if(findSecondLeadingJet)
	{
	  gen_2ndjet.SetPtEtaPhiE(
				  jetPt[secondLeadingJetIndex],
				  jetEta[secondLeadingJetIndex],
				  jetPhi[secondLeadingJetIndex],
				  jetEn[secondLeadingJetIndex]
				  );

	   h_y_2ndjet->Fill(gen_2ndjet.Rapidity());
	   h_pt_2ndjet->Fill(gen_2ndjet.Pt());
	}

      
  
   } // end of loop over entries


   for(int i=0; i<30; i++)
     if(nPass[i]>0)
       cout << "nPass["<< i << "]=" << nPass[i] << endl;
   
   /*
   std::string remword="/data2/syu/ggNTuple/";
   size_t pos = inputFile_.find(remword);
  
   if(pos!= std::string::npos)
     inputFile_.swap(inputFile_.erase(pos,remword.length()));
     


   TFile* outFile = new TFile(Form("debugCount_pstar%dto%d_"
				   "yb%.1lf""to""%.1lf_%s_%s",
				   (Int_t)pstarmin, (Int_t)pstarmax,
				   ybmin, ybmax,
				   MODE[mode],inputFile_.data()),"recreate");               
				   
   */

   TFile* outFile = new TFile("/home/syu/ggNtuple_scripts/CMSSW_4_2_8_patch7/src/runJob/debug.root",
			      "recreate");
   h_y_dirgamma -> Write();
   h_y_1stjet -> Write();
   h_y_2ndjet -> Write();

   h_pt_dirgamma -> Write();
   h_pt_1stjet -> Write();
   h_pt_2ndjet -> Write();


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


  Bool_t isEB=false;
  Bool_t isEE=false;
  if(fabs(phoSCEta[ipho]) < BARREL_MAXETA)isEB=true;
  if(fabs(phoSCEta[ipho]) > ENDCAP_MINETA && 
     fabs(phoSCEta[ipho]) < ENDCAP_MAXETA) isEE=true;

  if(!isEB && !isEE)return false;
  if(phoEt[ipho] < 85.0)return false;
  if(phoEt[ipho] > 95.0)return false;
  if(phoHoverE[ipho] > 0.05)return false;
  if(phohasPixelSeed[ipho]==1)return false;
//   if(phoEcalIsoDR04[ipho] > 4.2 +0.006 * phoEt[ipho])return false;
//   if(phoHcalIsoDR04[ipho] > 2.2 +0.0025* phoEt[ipho])return false;
//   if(phoTrkIsoHollowDR04[ipho] > 2.0 +0.001* phoEt[ipho])return false;

  if(phoEcalIsoDR04[ipho] > 4.2 +0.003 * phoEt[ipho])return false;
  if(phoHcalIsoDR04[ipho] > 2.2 +0.001* phoEt[ipho])return false;
  if(phoTrkIsoHollowDR04[ipho] > 2.0 +0.001* phoEt[ipho])return false;
 
  //  Float_t sumIso = phoEcalIsoDR04[ipho] + phoHcalIsoDR04[ipho] + phoTrkIsoHollowDR04[ipho];
  
  //  if(sumIso < 5.0)return false;
 
  //  if(isEB && phoSigmaIEtaIEta[ipho] > 0.010)return false;
  //  if(isEE && phoSigmaIEtaIEta[ipho] > 0.030)return false;

  return true;

}

