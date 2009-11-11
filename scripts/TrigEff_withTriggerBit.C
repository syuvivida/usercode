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

void TrigEff_withTriggerBit(std::string dirname, bool applyR9)
{
  int nbin=200;
  float xmin=0.0;
  float xmax=200.0;


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


  // HLT_Photon15_L1R / HLT_Photon10_L1R
  TH1F* h_get15to10deno = new TH1F("h_get15to10deno","",nbin,xmin,xmax);
  TH1F* h_get15to10numr = new TH1F("h_get15to10numr","",nbin,xmin,xmax);


  // chain in all the ncu ntuples in the same directory
  TChain* pho = new TChain("RECOTrigger/root");

  TSystemDirectory *base = new TSystemDirectory("root","root");
  std::string outerDir = "/home/syu/Photon_Code/CMSSW_3_1_2/src/syu/MyObjectCounter/test/triggerdata/";
  std::string filename = outerDir+ dirname + "/res";

  base->SetDirectory(filename.data());
  TList *listOfFiles = base->GetListOfFiles();
  TIter fileIt(listOfFiles);
  TFile *fileH = new TFile();
  int nfile=0;
  while(fileH = (TFile*)fileIt()) {
    std::string fileN = fileH->GetName();
//     std::string baseString = "triggerbit";
    std::string baseString = "root";
    if( fileH->IsFolder())  continue;
    if(fileN.find(baseString) == std::string::npos)continue;
    cout << fileN.data() << endl;
    nfile++;
    pho->Add(fileN.data());
  }

  std::cout << "Opened " << nfile << " files" << std::endl;

  PhoInfoBranches PhoInfo;
  PhoInfo.Register(pho);
  EvtInfoBranches EvtInfo;
  EvtInfo.Register(pho);
  
  // histogram definition
  Long64_t nentries = (Long64_t)pho->GetEntries();
  
  // loop over event entries
  for( Long64_t j=0; j< nentries; j++){

    // this line is very important!!
    pho->GetEntry(j);



    bool hasLoosePhoton = false;
    bool isFilled       = false;
    bool isFilled_Iso   = false;
    bool isFilled_mu    =false;
    bool isFilled_25    = false;
    bool isFilled_10    = false;
    bool isFilled_15to10= false;
    
    for(int i=0; i<PhoInfo.Size; i++)
      {
 	float thisEt = PhoInfo.Et[i];
	float thisEta = PhoInfo.Eta[i];

//  	int HLT = PhoInfo.Trig[i];
  	int HLT = EvtInfo.HLT;

	bool passMuonTrigger = HLT & TRIGGER::HLT_Mu5;
   	bool passL1EG5Trigger = HLT & TRIGGER::HLT_L1SingleEG5;
   	bool passL1EG8Trigger = HLT & TRIGGER::HLT_L1SingleEG8;
	bool passPhoton10Trigger = HLT & TRIGGER::HLT_Photon10_L1R;
  	bool passPhoton15Trigger = HLT & TRIGGER::HLT_Photon15_L1R;
	bool passPhoton25Trigger = HLT & TRIGGER::HLT_Photon25_L1R;
	bool passPhoton20IsoTrigger = HLT & TRIGGER::HLT_Photon20_LooseEcalIso_TrackIso_L1R;

	
	// has to contain a loose photon
     	if(!PhoInfo.IsLoose[i])continue;
	if(applyR9 && PhoInfo.R9[i] < 0.93)continue;
	
//  	if(applyR9 && PhoInfo.SCNCrystal[i] > 50)continue;

// 	// requires only barrel
//  	if(fabs(thisEta) > 1.4)continue;

	hasLoosePhoton = true;

	TLorentzVector rec_pho(0,0,0,0);
	rec_pho.SetPtEtaPhiE(
			     PhoInfo.Et[i],
			     PhoInfo.Eta[i],
			     PhoInfo.Phi[i],
			     PhoInfo.E[i]);
	

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
	
	// plotted vs genp photon Et
	
	
// 	// denomintor for HLT_15_L1R using L1_SingleEG5
    	if(hasLoosePhoton && passL1EG5Trigger && !isFilled)
 	  {
	    isFilled=true;
	    h_recgetdeno->Fill(thisEt);
	    if(isFromHard)
	      {
		h_getdeno->Fill(thisEt);
		h_etadeno->Fill(thisEta);
	      }
	    else if(isFromJet)
	      h_jetgetdeno->Fill(thisEt);
	    else if(isFromQuark)
	      h_qgetdeno->Fill(thisEt);
	    else if(isFromGluon)
	      h_ggetdeno->Fill(thisEt);

	    if(PhoInfo.GenPt[i]>0)
	      h_gengetdeno->Fill(PhoInfo.GenPt[i]);

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
		if(PhoInfo.GenPt[i]>0)
		  h_gengetnumr->Fill(PhoInfo.GenPt[i]);
	      }

	  
	  } // if there is a loose photon that passes base trigger	    


// 	// denomintor for HLT_15_L1R and L1_EG5 using HLT_Mu5
    	if(hasLoosePhoton && passMuonTrigger && !isFilled_mu)
 	  {
	    isFilled_mu=true;
	    h_mugetdeno->Fill(thisEt);
	    h_l1mugetdeno->Fill(thisEt);

	    // numerator
 	    if(passPhoton15Trigger)
	      h_mugetnumr->Fill(thisEt);

	    if(passL1EG5Trigger)
	      h_l1mugetnumr->Fill(thisEt);
	  
	  } // if there is a loose photon that passes base trigger	    


// 	// Photon_25_L1R
    	if(hasLoosePhoton && passL1EG5Trigger && !isFilled_25)
 	  {
	    isFilled_25=true;
	    h_25allgetdeno->Fill(thisEt);

	    // numerator
 	    if(passPhoton25Trigger)
	      h_25allgetnumr->Fill(thisEt);
	  
	  } // if there is a loose photon that passes base trigger	    



// 	// Photon_20_Iso
    	if(hasLoosePhoton && passL1EG5Trigger && !isFilled_Iso)
 	  {
	    isFilled_Iso=true;
	    h_isoallgetdeno->Fill(thisEt);

	    // numerator
 	    if(passPhoton20IsoTrigger)
	      h_isoallgetnumr->Fill(thisEt);
	  
	  } // if there is a loose photon that passes base trigger	    


	// Photon_10_R relative to L1EG5
	if(hasLoosePhoton && passL1EG5Trigger && !isFilled_10 && isFromHard)
	  {
	    isFilled_10=true;
	    h_get10deno->Fill(thisEt);
	    
	    if(passPhoton10Trigger)
	      h_get10numr->Fill(thisEt);
	  }


	// Photon_15_R relative to Photon_10_R
	if(hasLoosePhoton && passPhoton10Trigger && !isFilled_15to10 
	   && isFromHard)
	  {
	    isFilled_15to10=true;
	    h_get15to10deno->Fill(thisEt);
	    
	    if(passPhoton15Trigger)
	      h_get15to10numr->Fill(thisEt);
	  }




      } // end of loop over photon objects


    
    
  } // end of loop over event entries

  std::string histoFile = outerDir + dirname + "/" + dirname + "_debug2_default.root";
  if(applyR9) histoFile = outerDir + dirname + "/" + dirname + "_debug2_R9.root";

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

  h_get15to10deno ->Write();
  h_get15to10numr ->Write();

  outFile->Close();
  
  

}


