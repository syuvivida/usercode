#include <TH1.h>
#include <TFile.h>
#include <TChain.h>
#include <TBranch.h>
#include <iostream>
#include <fstream>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystemDirectory.h>
#include <TProfile.h>					
#include <TList.h>
#include <TF1.h>
#include <vector>
#include <TLorentzVector.h>
#include "format.h"

using namespace std;

void TrigEff_withL3Match(TChain* pho, 
	  std::string filestring, bool isBarrel, bool applyR9)
{
//   int nbin=200;
  int nbin=100;
  float xmin=0.0;
  float xmax=200.0;

  int ncount[10]={0};
  TH1F* h_etadeno = new TH1F("h_etadeno","Reconstructed and matched photon "
			     "#eta before photon trigger cuts", 120,-3.0,3.0);
  TH1F* h_etanumr = new TH1F("h_etanumr","Reconstructed and matched photon "
			     "#eta after photon trigger cuts", 120,-3.0,3.0);


  TH1F* h_getdeno = new TH1F("h_getdeno","Reconstructed and matched photon "
			     "Et before photon trigger cuts", nbin,xmin,xmax);
  TH1F* h_getnumr = new TH1F("h_getnumr","Reconstructed and matched photon "
			     "Et after photon trigger cuts", nbin,xmin,xmax);
  TH1F* h_recgetdeno = new TH1F("h_recgetdeno","Reconstructed photon Et "
				"before photon trigger cuts", nbin,xmin,xmax);
  TH1F* h_recgetnumr = new TH1F("h_recgetnumr","Reconstructed photon Et "
				"after photon trigger cuts", nbin,xmin,xmax);
  TH1F* h_gengetdeno = new TH1F("h_gengetdeno","Generated photon Et before "
				"photon trigger cuts", nbin,xmin,xmax);
  TH1F* h_gengetnumr = new TH1F("h_gengetnumr","Generated photon Et after "
				"photon trigger cuts", nbin,xmin,xmax);
  TH1F* h_jetgetdeno = new TH1F("h_jetgetdeno","Reconstructed and jet photon"
				" Et before photon trigger cuts", 
				nbin,xmin,xmax);
  TH1F* h_jetgetnumr = new TH1F("h_jetgetnumr","Reconstructed and jet photon"
				" Et after photon trigger cuts", 
				nbin,xmin,xmax);
  TH1F* h_qgetdeno   = new TH1F("h_qgetdeno","Reconstructed and quark photon"
				" Et before photon trigger cuts", 
				nbin,xmin,xmax);
  TH1F* h_qgetnumr   = new TH1F("h_qgetnumr","Reconstructed and quark photon"
				" Et after photon trigger cuts", 
				nbin,xmin,xmax);
  TH1F* h_ggetdeno   = new TH1F("h_ggetdeno","Reconstructed and gluon photon"
				" Et before photon trigger cuts", 
				nbin,xmin,xmax);
  TH1F* h_ggetnumr   = new TH1F("h_ggetnumr","Reconstructed and gluon photon"
				" Et after photon trigger cuts", 
				nbin,xmin,xmax);

  TH1F* h_mugetdeno = new TH1F("h_mugetdeno","Reconstructed photon Et "
			       "before photon trigger cuts", nbin,xmin,xmax);
  TH1F* h_mugetnumr = new TH1F("h_mugetnumr","Reconstructed photon Et "
			       "after photon trigger cuts", nbin,xmin,xmax);

  TH1F* h_25allgetdeno = new TH1F("h_25allgetdeno","Reconstructed photon Et "
			       "before photon trigger cuts", nbin,xmin,xmax);
  TH1F* h_25allgetnumr = new TH1F("h_25allgetnumr","Reconstructed photon Et "
			       "after photon trigger cuts", nbin,xmin,xmax);


  TH1F*  h_isoallgetdeno =  new TH1F("h_isoallgetdeno","Reconstructed and matched photon "
			       "Et before photon trigger cuts", nbin,xmin,xmax);
  TH1F* h_isoallgetnumr =  new TH1F("h_isoallgetnumr","Reconstructed and matched photon "
			       "Et after photon trigger cuts", nbin,xmin,xmax);


  TH1F* h_l1mugetdeno = new TH1F("h_l1mugetdeno","Reconstructed photon Et "
			       "before photon trigger cuts", nbin,xmin,xmax);
  TH1F* h_l1mugetnumr = new TH1F("h_l1mugetnumr","Reconstructed photon Et "
			       "after photon trigger cuts", nbin,xmin,xmax);

  // HLT_Photon10_L1R / L1EG5
  TH1F* h_get10deno = new TH1F("h_get10deno","",nbin,xmin,xmax);
  TH1F* h_get10numr = new TH1F("h_get10numr","",nbin,xmin,xmax);

  // L1EG5 / reconstructed photon
  TH1F* h_get5deno = new TH1F("h_get5deno","",nbin,xmin,xmax);
  TH1F* h_get5numr = new TH1F("h_get5numr","",nbin,xmin,xmax);

  // L1EG5 / reconstructed jet photon
  TH1F* h_jetget5deno = new TH1F("h_jetget5deno","",nbin,xmin,xmax);
  TH1F* h_jetget5numr = new TH1F("h_jetget5numr","",nbin,xmin,xmax);


  PhoInfoBranches PhoInfo;
  PhoInfo.Register(pho);
  EvtInfoBranches EvtInfo;
  EvtInfo.Register(pho);
  GenInfoBranches GenInfo;
  GenInfo.Register(pho);
  
  bool isQCD = true;
  
  // histogram definition
  Long64_t nentries = (Long64_t)pho->GetEntries();
  
  // loop over event entries
  for( Long64_t j=0; j< nentries; j++){

    // this line is very important!!
    pho->GetEntry(j);

    if(j<5)
      {
	for(int ab=0; ab<GenInfo.Size; ab++)
	  {
	    if(GenInfo.PID[ab]==22 && GenInfo.MPID[ab]==22)
	      {
		isQCD=false; break;
	      }
	    else if(GenInfo.PID[ab]==22 && GenInfo.MPID[ab]==0)
	      {
		isQCD=false; break;
	      }
	  }

	cout << "I am QCD " << isQCD << endl;

      }

    // loop over photons
    for(int i=0; i<PhoInfo.Size; i++)
      {
 	float thisEt = PhoInfo.Et[i];
	float thisEta = PhoInfo.Eta[i];
	bool hasLoosePhoton = false;
  	int HLT = PhoInfo.Trig[i];

     	bool passL1EG5Trigger = HLT & TRIGGER::HLT_L1SingleEG5;
     	bool passL1EG8Trigger = HLT & TRIGGER::HLT_L1SingleEG8;
	bool passPhoton10Trigger = HLT & TRIGGER::HLT_Photon10_L1R;
  	bool passPhoton15Trigger = HLT & TRIGGER::HLT_Photon15_L1R;
	bool passPhoton25Trigger = HLT & TRIGGER::HLT_Photon25_L1R;
	bool passPhoton20IsoTrigger = HLT & TRIGGER::HLT_Photon20_LooseEcalIso_TrackIso_L1R;

	ncount[0]++;
	// has to contain a loose photon
     	if(PhoInfo.IsLoose[i]!=1)continue;
	ncount[1]++;

	if(applyR9 && PhoInfo.R9[i] < 0.93)continue;
	ncount[2]++;
	
	if(isBarrel && PhoInfo.Location[i]!=1)continue;
	ncount[3]++;

	if(!isBarrel && PhoInfo.Location[i]!=2)continue;
	ncount[4]++;
	
// 	// requires only barrel

	hasLoosePhoton = true;

	// matched to a genp photon
	int pid = PhoInfo.GenPID[i];
	int mpid = PhoInfo.GenMomPID[i]; 
	bool isFromHard = pid==22
	  && mpid==22;
	
	bool isFromJet = pid==22
	  && abs(mpid)> 50;
	
	bool isFromQuark = pid==22
	  && abs(mpid)<10;

	bool isFromGluon = pid==22
	  && abs(mpid)==21;
	
	
 	// denomintor for HLT_15_L1R using L1_SingleEG5
    	if(hasLoosePhoton && passL1EG5Trigger)
 	  {
	    ncount[5]++;
	    h_recgetdeno->Fill(thisEt);
	    if(isFromHard)
	      {
		ncount[6]++;
		h_getdeno->Fill(thisEt);
		h_etadeno->Fill(thisEta);
	      }
	    else if(isFromJet)
	      h_jetgetdeno->Fill(thisEt);
	    else if(isFromQuark)
	      h_qgetdeno->Fill(thisEt);
	    else if(isFromGluon)
	      h_ggetdeno->Fill(thisEt);

	    // numerator

  	    if(passPhoton15Trigger)
	      {
		h_recgetnumr->Fill(thisEt);
		if(isFromHard)
		  {
		    h_getnumr->Fill(thisEt);
		    h_etanumr->Fill(thisEta);
		  }
		else if(isFromJet)
		  h_jetgetnumr->Fill(thisEt);
		else if(isFromQuark)
		  h_qgetnumr->Fill(thisEt);
		else if(isFromGluon)
		  h_ggetnumr->Fill(thisEt);
	      }
	  
	  } // if there is a loose photon that passes base trigger	    

	bool isEitherHardOrJet = (isQCD && isFromJet) || (!isQCD && isFromHard);
// 	// denomintor for HLT_15_L1R and L1_EG5 using L1_SingleEG8
    	if(hasLoosePhoton && isEitherHardOrJet)
 	  {
	    h_l1mugetdeno->Fill(thisEt);
	    if(passL1EG8Trigger)
	      {
		h_mugetdeno->Fill(thisEt);
		h_l1mugetnumr->Fill(thisEt);
	      }
	    // numerator
 	    if(passL1EG8Trigger && passPhoton15Trigger)
	      h_mugetnumr->Fill(thisEt);
	  } // if there is a loose photon that passes base trigger	    


// 	// Photon_25_L1R
    	if(hasLoosePhoton && passL1EG5Trigger && isEitherHardOrJet)
 	  {
	    h_25allgetdeno->Fill(thisEt);

	    // numerator
 	    if(passPhoton25Trigger)
	      h_25allgetnumr->Fill(thisEt);
	  
	  } // if there is a loose photon that passes base trigger	    



// 	// Photon_20_Iso
    	if(hasLoosePhoton && passL1EG5Trigger && isEitherHardOrJet)
 	  {
	    h_isoallgetdeno->Fill(thisEt);

	    // numerator
 	    if(passPhoton20IsoTrigger)
	      h_isoallgetnumr->Fill(thisEt);
	  
	  } // if there is a loose photon that passes base trigger	    


	// Photon_10_R relative to L1EG5
	if(hasLoosePhoton && passL1EG5Trigger && isEitherHardOrJet)
	  {
	    h_get10deno->Fill(thisEt);	    
	    if(passPhoton10Trigger)
	      h_get10numr->Fill(thisEt);
	  }

	
	// L1EG5 relative to reconstructed photon
  	if(hasLoosePhoton)
  	  {
  	    if(isFromHard)
	      {
		h_get5deno->Fill(thisEt);
		h_gengetdeno->Fill(PhoInfo.GenPt[i]);
	      }

  	    if(isFromHard && passL1EG5Trigger)
	      {
		h_get5numr->Fill(thisEt);
		h_gengetnumr->Fill(PhoInfo.GenPt[i]);
	      }


  	    if(isFromJet)
  	      h_jetget5deno->Fill(thisEt);
  	    if(isFromJet && passL1EG5Trigger)
  	      h_jetget5numr->Fill(thisEt);

  	  }





      } // end of loop over photon objects


    
    
  } // end of loop over event entries
  
  std::string detector  = isBarrel? "barrel" : "endcap";
  std::string R9name    = applyR9?  "R9"     : "default";
  
  std::string histoFile =  "/home/syu/Photon_Code/CMSSW_3_1_2/src/syu/MyObjectCounter/test/script_output/test_" + filestring + "_" + detector +  "_" + R9name
    + ".root";

  TFile* outFile = new TFile(histoFile.data(),"recreate");

  h_etadeno->Write();
  h_etanumr->Write();

  h_mugetdeno->Write();
  h_mugetnumr->Write();


  h_l1mugetdeno->Write();
  h_l1mugetnumr->Write();

  h_25allgetdeno->Write();
  h_25allgetnumr->Write();

  h_isoallgetdeno->Write();
  h_isoallgetnumr->Write();


  h_getdeno->Write();
  h_getnumr->Write();
  h_recgetdeno->Write();
  h_recgetnumr->Write();
  h_gengetdeno->Write();
  h_gengetnumr->Write();


  h_jetgetdeno->Write();
  h_jetgetnumr->Write();
  h_qgetdeno  ->Write();
  h_qgetnumr  ->Write();
  h_ggetdeno  ->Write();
  h_ggetnumr  ->Write();


  h_get10deno->Write();
  h_get10numr->Write();

  h_get5deno->Write();
  h_get5numr->Write();

  h_jetget5deno->Write();
  h_jetget5numr->Write();

  outFile->Close();
  for(int i=0; i<10;i++)
    if(ncount[i]>0)cout << "ncount[" << i << "]=" << ncount[i] << endl;
  

}


