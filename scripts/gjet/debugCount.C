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

   TH1F* h_pt_trig30 = (TH1F*) h_pt_template->Clone("h_pt_trig30");
   TH1F* h_pt_trig50 = (TH1F*) h_pt_template->Clone("h_pt_trig50");
   TH1F* h_pt_trig75 = (TH1F*) h_pt_template->Clone("h_pt_trig75");
   TH1F* h_pt_trig90 = (TH1F*) h_pt_template->Clone("h_pt_trig90");
   TH1F* h_pt_trig135 = (TH1F*) h_pt_template->Clone("h_pt_trig135");
   
   // rapidity distribution
   TH1F* h_y_template = new TH1F("h_y_template","",100,-5,5);
   TH1F* h_y_dirgamma = (TH1F*) h_y_template->Clone("h_y_dirgamma");
   TH1F* h_y_1stjet = (TH1F*) h_y_template->Clone("h_y_1stjet");
   TH1F* h_y_2ndjet = (TH1F*) h_y_template->Clone("h_y_2ndjet");


   const int NCOUNTS=40;
   TH1I* h_npass = new TH1I("h_npass","",NCOUNTS,-0.5,NCOUNTS-0.5);
   h_npass->SetXTitle("nPass");

   Long64_t nPass[NCOUNTS]={0};

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //if(jentry>50000)break;

      nPass[0] ++;
      h_npass->Fill(0);

      Int_t HLT_Photon30_CaloIdVL = HLTIndex[63];
      
      Int_t HLT_Photon30_CaloIdVL_IsoL = HLTIndex[64];

      Int_t HLT_Photon50_CaloIdVL = HLTIndex[65];

      Int_t HLT_Photon50_CaloIdVL_IsoL = HLTIndex[66];

      Int_t HLT_Photon75_CaloIdVL = HLTIndex[30];

      Int_t HLT_Photon75_CaloIdVL_IsoL = HLTIndex[29];

      Int_t HLT_Photon90_CaloIdVL = HLTIndex[67];

      Int_t HLT_Photon90_CaloIdVL_IsoL = HLTIndex[15];

      Int_t HLT_Photon135 = HLTIndex[31];
      

      Bool_t fire_HLT_Photon30_CaloIdVL = HLT_Photon30_CaloIdVL>0 && HLT[HLT_Photon30_CaloIdVL]>0;
      Bool_t fire_HLT_Photon30_CaloIdVL_IsoL = HLT_Photon30_CaloIdVL_IsoL>0 && HLT[HLT_Photon30_CaloIdVL_IsoL]>0;
      Bool_t fire_HLT_Photon50_CaloIdVL = HLT_Photon50_CaloIdVL>0 && HLT[HLT_Photon50_CaloIdVL]>0;

      Bool_t fire_HLT_Photon50_CaloIdVL_IsoL = HLT_Photon50_CaloIdVL_IsoL>0 && HLT[HLT_Photon50_CaloIdVL_IsoL]>0;
      Bool_t fire_HLT_Photon75_CaloIdVL = HLT_Photon75_CaloIdVL>0 && HLT[HLT_Photon75_CaloIdVL]>0;
      Bool_t fire_HLT_Photon75_CaloIdVL_IsoL = HLT_Photon75_CaloIdVL_IsoL>0 && HLT[HLT_Photon75_CaloIdVL_IsoL]>0;
      Bool_t fire_HLT_Photon90_CaloIdVL = HLT_Photon90_CaloIdVL>0 && HLT[HLT_Photon90_CaloIdVL]>0;
      Bool_t fire_HLT_Photon90_CaloIdVL_IsoL = HLT_Photon90_CaloIdVL_IsoL>0 && HLT[HLT_Photon90_CaloIdVL_IsoL]>0;
      Bool_t fire_HLT_Photon135 = HLT_Photon135>0 && HLT[HLT_Photon135]>0;

      if(fire_HLT_Photon30_CaloIdVL)
	{
	  nPass[1]++;
	  h_npass->Fill(1);
	}
      if(fire_HLT_Photon30_CaloIdVL_IsoL)
	{
	  nPass[2]++;
	  h_npass->Fill(2);
	}

      if(fire_HLT_Photon50_CaloIdVL)
	{
	  nPass[3]++;
	  h_npass->Fill(3);
	}
      if(fire_HLT_Photon50_CaloIdVL_IsoL)
	{
	  nPass[4]++;
	  h_npass->Fill(4);
	}

      if(fire_HLT_Photon75_CaloIdVL)
	{ 
	  nPass[5]++;
	  h_npass->Fill(5);
	}

      if(fire_HLT_Photon75_CaloIdVL_IsoL)
	{
	  nPass[6]++;
	  h_npass->Fill(6);
	}

      if(fire_HLT_Photon90_CaloIdVL)
	{
	  nPass[7]++;
	  h_npass->Fill(7);
	}
      if(fire_HLT_Photon90_CaloIdVL_IsoL)
	{
	  nPass[8]++;
	  h_npass->Fill(8);
	}
      if(fire_HLT_Photon135)
	{
	  nPass[9]++;
	  h_npass->Fill(9);
	}

      /*
      if(HLT_Photon75_CaloIdVL      > 0 && HLT[HLT_Photon75_CaloIdVL]>1)
	{
	  nPass[10]++; // prescaled
	  h_npass->Fill(10);
	}

      if(HLT_Photon75_CaloIdVL      > 0 && HLT[HLT_Photon75_CaloIdVL]!=1)
	{
	  continue;
	  nPass[11]++;
	  h_npass->Fill(11);
	}
      */


      Int_t ngood_vtx=IsVtxGood;
      if(ngood_vtx==0)continue;

      nPass[12] ++;
      h_npass->Fill(12);


      bool findLeadingPhoton = false;
      int leadingPhotonIndex = -1;
      Float_t phoMaxPt = -9999.;

      // now find a good leading photon
      for(int ipho=0; ipho < nPho; ipho++){
	
	
	// before applying photon cuts
	if(fire_HLT_Photon30_CaloIdVL && phoEt[ipho]>=40 && phoEt[ipho] < 60 && fabs(phoSCEta[ipho])<2.5)
	  {
	    nPass[30]++;
	    h_npass->Fill(30);
	  }
 	if(fire_HLT_Photon50_CaloIdVL && phoEt[ipho]>=60 && phoEt[ipho] < 85 && fabs(phoSCEta[ipho])<2.5)
	  {
	    nPass[31]++;
	    h_npass->Fill(31);
	  }
   
  	if(fire_HLT_Photon75_CaloIdVL && phoEt[ipho]>=85 && phoEt[ipho] < 100 && fabs(phoSCEta[ipho])<2.5)
	  {
	    nPass[32]++;
	    h_npass->Fill(32);
	  }
   
	if(fire_HLT_Photon90_CaloIdVL && phoEt[ipho]>=100 && phoEt[ipho] < 145 && fabs(phoSCEta[ipho])<2.5)
	  {
	    nPass[33]++;
	    h_npass->Fill(33);
	  }
   
	if(fire_HLT_Photon135 && phoEt[ipho]>=145 && phoEt[ipho] < 300 && fabs(phoSCEta[ipho])<2.5)
	  {
	    nPass[34]++;
	    h_npass->Fill(34);
	  }
   


	// after applying photon IDs
	if(!isGoodPho(ientry,ipho))continue;

	nPass[13]++;
	h_npass->Fill(13);

	h_pt_dirgamma->Fill(phoEt[ipho]);
	h_y_dirgamma ->Fill(phoSCEta[ipho]);

	if( fabs(phoSCEta[ipho])<0.9 )
	  {
	    nPass[14]++;
	    h_npass->Fill(14);
	  }
	else if( fabs(phoSCEta[ipho])>0.9 && 
		 fabs(phoSCEta[ipho])<BARREL_MAXETA )
	  {
	    nPass[15]++;
	    h_npass->Fill(15);
	  }

	else if( fabs(phoSCEta[ipho])>ENDCAP_MINETA &&
		 fabs(phoSCEta[ipho])< 2.1)
	  {
	    nPass[16]++;
	    h_npass->Fill(16);	    
	  }
	else if( fabs(phoSCEta[ipho])> 2.1 && 
		 fabs(phoSCEta[ipho])< ENDCAP_MAXETA)
	  {
	    nPass[17]++;
	    h_npass->Fill(17);
	  }
	
	if(fabs(phoSCEta[ipho]) < BARREL_MAXETA && phoSigmaIEtaIEta[ipho] < 0.010)
	  {
	    nPass[18]++;
	    h_npass->Fill(18);
	  }
	if(fabs(phoSCEta[ipho]) < ENDCAP_MAXETA && 
	   fabs(phoSCEta[ipho]) > ENDCAP_MINETA 
	   && phoSigmaIEtaIEta[ipho] < 0.028)
	  {
	    nPass[19]++;
	    h_npass->Fill(19);
	  }
	    
	if(phoEt[ipho] > phoMaxPt)
	  {
	    phoMaxPt = phoEt[ipho];
	    leadingPhotonIndex= ipho;
	    findLeadingPhoton = true;
	  }
	

	// now use only one trigger path for each pt bin
	if(fire_HLT_Photon30_CaloIdVL)
	  {
	    h_pt_trig30->Fill(phoEt[ipho]);
	    if(phoEt[ipho]>=40 && phoEt[ipho] < 60)
	      {
		nPass[20]++;
		h_npass->Fill(20);
	      }
	  }
   

	if(fire_HLT_Photon50_CaloIdVL)
	  {
	    h_pt_trig50->Fill(phoEt[ipho]);
	    if(phoEt[ipho]>=60 && phoEt[ipho] < 85)
	      {
		nPass[21]++;
		h_npass->Fill(21);
	      }
	  }
	

	if(fire_HLT_Photon75_CaloIdVL)
	  {
	    h_pt_trig75->Fill(phoEt[ipho]);
	    if(phoEt[ipho]>=85 && phoEt[ipho] < 100)
	      {
		nPass[22]++;
		h_npass->Fill(22);
	      }
	  }


	if(fire_HLT_Photon90_CaloIdVL)
	  {
	    h_pt_trig90->Fill(phoEt[ipho]);
	    if(phoEt[ipho]>=100 && phoEt[ipho] < 145)
	      {
		nPass[23]++;
		h_npass->Fill(23);
	      }
	  }

	if(fire_HLT_Photon135)
	  {
	    h_pt_trig135->Fill(phoEt[ipho]);
	    if(phoEt[ipho]>=145 && phoEt[ipho] < 300)
	      {
		nPass[24]++;
		h_npass->Fill(24);
	      }
	  }
 
	

      } // end of leading photon search

      if(!findLeadingPhoton)continue;


	
      nPass[25]++;
      h_npass->Fill(25);

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
	  nPass[26]++;
	  h_npass->Fill(26);
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


   for(int i=0; i<NCOUNTS; i++)
       if(nPass[i]>0)
	 cout << "nPass["<< i << "]=" << nPass[i] << endl;

   std::string remword="/data4/syu/";

   if(_inputDirName.find("v4") != std::string::npos)
     remword="/data2/syu/";

   size_t pos = _inputDirName.find(remword);
  
   if(pos!= std::string::npos)
     _inputDirName.swap(_inputDirName.erase(pos,remword.length()));
     


   TFile* outFile = new TFile(Form("/home/syu/ggNtuple_scripts/CMSSW_4_2_8_patch7/src/runJob/11005_%s.root",_inputDirName.data()),"recreate");               
				   
   h_y_dirgamma -> Write();
   h_y_1stjet -> Write();
   h_y_2ndjet -> Write();

   h_pt_dirgamma -> Write();
   h_pt_1stjet -> Write();
   h_pt_2ndjet -> Write();
   
   h_npass->Write();

   h_pt_trig30->Write();
   h_pt_trig50->Write();
   h_pt_trig75->Write();
   h_pt_trig90->Write();
   h_pt_trig135->Write();
 
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
  //  if(phoEt[ipho] < 85.0)return false;
  //  if(phoEt[ipho] > 95.0)return false;
  if(phoHoverE[ipho] > 0.05)return false;
  if(phohasPixelSeed[ipho]==1)return false;
//   if(phoEcalIsoDR04[ipho] > 4.2 +0.006 * phoEt[ipho])return false;
//   if(phoHcalIsoDR04[ipho] > 2.2 +0.0025* phoEt[ipho])return false;
//   if(phoTrkIsoHollowDR04[ipho] > 2.0 +0.001* phoEt[ipho])return false;

//  if(phoEcalIsoDR04[ipho] > 4.2 +0.003 * phoEt[ipho])return false;
//  if(phoHcalIsoDR04[ipho] > 2.2 +0.001* phoEt[ipho])return false;
//  if(phoTrkIsoHollowDR04[ipho] > 2.0 +0.001* phoEt[ipho])return false;
 
  Float_t sumIso = phoEcalIsoDR04[ipho] + phoHcalIsoDR04[ipho] + phoTrkIsoHollowDR04[ipho];
  
  //  if(sumIso < 5.0)return false;
 
  if(isEB && phoSigmaIEtaIEta[ipho] > 0.010)return false;
  if(isEE && phoSigmaIEtaIEta[ipho] > 0.028)return false;

  return true;

}

