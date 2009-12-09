#ifndef FORMAT_HH
#define FORMAT_HH

// here, I define the format of root trees, for event tree, 
// photon tree, lepton tree, generator tree

#include <TROOT.h>
#include <TSystem.h>
#include <TObject.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>

#define MAX_LEPTONS 	64
#define MAX_PHOTONS     64
#define MAX_TRIGGER    120
#define MAX_GENS        32
#define dummy          -9999

namespace syu{
class EvtInfoBranches {
public:
  int   RunNo;
  int   EvtNo;
  int   HLT;
  
  void Initialize(){
    
    RunNo = EvtNo = HLT = dummy;

  }

  void Register(TTree *root) {
	root->Branch("EvtInfo.RunNo"	 , &RunNo     , "EvtInfo.RunNo/I"     );
	root->Branch("EvtInfo.EvtNo"	 , &EvtNo     , "EvtInfo.EvtNo/I"     );   
	root->Branch("EvtInfo.HLT"	 , &HLT       , "EvtInfo.HLT/I"     );   
  }										    
}; // end of EvtInfoBranches


class PhoInfoBranches {
public:
  int	Size; 
  int   Index[MAX_PHOTONS];
  int   Location[MAX_PHOTONS];
  int   IsLoose[MAX_PHOTONS];
  int   Trig[MAX_PHOTONS];
  float L1Pt[MAX_PHOTONS];
  float L3Pt[MAX_PHOTONS];
  float E[MAX_PHOTONS];
  float Et[MAX_PHOTONS];
  float Pz[MAX_PHOTONS];
  float Eta[MAX_PHOTONS];
  float Phi[MAX_PHOTONS];
  float R9[MAX_PHOTONS];
  float TrkPtSum[MAX_PHOTONS];
  float EcalRecHitEtSum[MAX_PHOTONS];
  float HcalTowerEtSum[MAX_PHOTONS];
  float HoverE[MAX_PHOTONS];
  float SCE[MAX_PHOTONS];
  float SCEt[MAX_PHOTONS];
  float SCEta[MAX_PHOTONS];
  float SCPhi[MAX_PHOTONS];
  float SCEtaWidth[MAX_PHOTONS];
  float SCPhiWidth[MAX_PHOTONS];
  int   SCNCrystal[MAX_PHOTONS];
  float SigEta[MAX_PHOTONS];
  float SigIEta[MAX_PHOTONS];
  int   GenPID[MAX_PHOTONS];
  int   GenGMomPID[MAX_PHOTONS];
  int   GenMomPID[MAX_PHOTONS];
  float GenPt[MAX_PHOTONS];
  float GenMomPt[MAX_PHOTONS];


  void Initialize(){

    Size = 0;
    for(int i=0; i< MAX_PHOTONS; i++){

      Index[i]     =dummy;
      Location[i]  =dummy;
      IsLoose[i]   =dummy;
      Trig[i]      =0;
      L1Pt[i]      =dummy;
      L3Pt[i]      =dummy;
      E[i]    	   =dummy;
      Et[i]   	   =dummy;
      Pz[i]   	   =dummy;
      Eta[i]  	   =dummy;
      Phi[i]  	   =dummy;
      R9[i]   	   =dummy;
      TrkPtSum[i]  =dummy;
      EcalRecHitEtSum[i]=dummy;
      HcalTowerEtSum[i]=dummy;
      HoverE[i]	   =dummy;
      SCE[i]	   =dummy;
      SCEt[i]	   =dummy;
      SCEta[i]	   =dummy;
      SCPhi[i]	   =dummy;
      SCEtaWidth[i]=dummy;
      SCPhiWidth[i]=dummy;
      SCNCrystal[i]=dummy;
      SigEta[i]    =dummy;
      SigIEta[i]   =dummy;
      GenPID[i]    =dummy;
      GenGMomPID[i]=dummy;
      GenMomPID[i] =dummy;
      GenPt[i]     =dummy;
      GenMomPt[i]  =dummy;

    }
  }
  
  //Candidate* CandRef[MAX_PHOTONS]; // backward pointer to the PAT objects
  
  void Register(TTree *root) {
	root->Branch("PhoInfo.Size"	  , &Size	   , "PhoInfo.Size/I"			  );
	root->Branch("PhoInfo.Index"	  , &Index[0]	   , "PhoInfo.Index[PhoInfo.Size]/I"	  );
	root->Branch("PhoInfo.Location"	  , &Location[0]   , "PhoInfo.Location[PhoInfo.Size]/I"	  );
	root->Branch("PhoInfo.IsLoose"	  , &IsLoose[0]	   , "PhoInfo.IsLoose[PhoInfo.Size]/I"	  );
	root->Branch("PhoInfo.Trig"	  , &Trig[0]	   , "PhoInfo.Trig[PhoInfo.Size]/I"	  );
	root->Branch("PhoInfo.L1Pt"	  , &L1Pt[0]	   , "PhoInfo.L1Pt[PhoInfo.Size]/F"	  );
	root->Branch("PhoInfo.L3Pt"	  , &L3Pt[0]	   , "PhoInfo.L3Pt[PhoInfo.Size]/F"	  );
	root->Branch("PhoInfo.E"          , &E[0]          , "PhoInfo.E[PhoInfo.Size]/F"          );
	root->Branch("PhoInfo.Et"         , &Et[0]         , "PhoInfo.Et[PhoInfo.Size]/F"         );
	root->Branch("PhoInfo.Pz"         , &Pz[0]         , "PhoInfo.Pz[PhoInfo.Size]/F"         );
	root->Branch("PhoInfo.Eta"        , &Eta[0]        , "PhoInfo.Eta[PhoInfo.Size]/F"        );
	root->Branch("PhoInfo.Phi"        , &Phi[0]        , "PhoInfo.Phi[PhoInfo.Size]/F"        );
	root->Branch("PhoInfo.R9"         , &R9[0]         , "PhoInfo.R9[PhoInfo.Size]/F"         );
	root->Branch("PhoInfo.TrkPtSum"   , &TrkPtSum[0]   , "PhoInfo.TrkPtSum[PhoInfo.Size]/F"   );
	root->Branch("PhoInfo.EcalRecHitEtSum", &EcalRecHitEtSum[0], "PhoInfo.EcalRecHitEtSum[PhoInfo.Size]/F");
	root->Branch("PhoInfo.HcalTowerEtSum",  &HcalTowerEtSum[0], "PhoInfo.HcalTowerEtSum[PhoInfo.Size]/F");
	root->Branch("PhoInfo.HoverE"     , &HoverE[0]     , "PhoInfo.HoverE[PhoInfo.Size]/F"     );
	root->Branch("PhoInfo.SCE"        , &SCE[0]        , "PhoInfo.SCE[PhoInfo.Size]/F");
	root->Branch("PhoInfo.SCEt"       , &SCEt[0]       , "PhoInfo.SCEt[PhoInfo.Size]/F");	
	root->Branch("PhoInfo.SCEta"      , &SCEta[0]      , "PhoInfo.SCEta[PhoInfo.Size]/F");	
	root->Branch("PhoInfo.SCPhi"      , &SCPhi[0]      , "PhoInfo.SCPhi[PhoInfo.Size]/F");	
	root->Branch("PhoInfo.SCEtaWidth" , &SCEtaWidth[0] , "PhoInfo.SCEtaWidth[PhoInfo.Size]/F");	
	root->Branch("PhoInfo.SCPhiWidth" , &SCPhiWidth[0] , "PhoInfo.SCPhiWidth[PhoInfo.Size]/F");	
	root->Branch("PhoInfo.SCNCrystal",  &SCNCrystal[0] , "PhoInfo.SCNCrystal[PhoInfo.Size]/I");
	root->Branch("PhoInfo.SigEta",      &SigEta[0]     , "PhoInfo.SigEta[PhoInfo.Size]/F");
	root->Branch("PhoInfo.SigIEta",     &SigIEta[0]    , "PhoInfo.SigIEta[PhoInfo.Size]/F");
        root->Branch("PhoInfo.GenPID",      &GenPID[0]     , "PhoInfo.GenPID[PhoInfo.Size]/I");
	root->Branch("PhoInfo.GenGMomPID" , &GenGMomPID[0] , "PhoInfo.GenGMomPID[PhoInfo.Size]/I" );
	root->Branch("PhoInfo.GenMomPID"  , &GenMomPID[0]  , "PhoInfo.GenMomPID[PhoInfo.Size]/I" );
	root->Branch("PhoInfo.GenPt"   , &GenPt[0]   , "PhoInfo.GenPt[PhoInfo.Size]/F" );
	root->Branch("PhoInfo.GenMomPt"   , &GenMomPt[0]   , "PhoInfo.GenMomPt[PhoInfo.Size]/F" );

  }
};  // end of PhoInfoBranches

class LepInfoBranches {
public:
  int	Size; 
  int	Index[MAX_LEPTONS];
  int	LeptonType[MAX_LEPTONS];
  int	Charge[MAX_LEPTONS];
  float Pt[MAX_LEPTONS];
  float Eta[MAX_LEPTONS];
  float Phi[MAX_LEPTONS];
  float TrackIso[MAX_LEPTONS];
  float EcalIso[MAX_LEPTONS];
  float HcalIso[MAX_LEPTONS];
  
  //Candidate* CandRef[MAX_LEPTONS]; // backward pointer to the PAT objects
  void Initialize(){

    Size = 0;
    for(int i=0; i < MAX_LEPTONS; i++){

      Index[i]      = dummy;
      LeptonType[i] = dummy;
      Charge[i]     =  0;
      Pt[i]         = dummy;
      Eta[i]        = dummy;
      Phi[i]        = dummy;
      TrackIso[i]   = dummy;
      EcalIso[i]    = dummy;
      HcalIso[i]    = dummy;
  
    }
  }
  
  void Register(TTree *root) {
	root->Branch("LepInfo.Size"	  , &Size	   , "LepInfo.Size/I"			  );
	root->Branch("LepInfo.Index"	  , &Index[0]	   , "LepInfo.Index[LepInfo.Size]/I"	  );
	root->Branch("LepInfo.LeptonType" , &LeptonType[0] , "LepInfo.LeptonType[LepInfo.Size]/I" );
	root->Branch("LepInfo.Charge"	  , &Charge[0]     , "LepInfo.Charge[LepInfo.Size]/I"	  );
	root->Branch("LepInfo.Pt"	  , &Pt[0]	   , "LepInfo.Pt[LepInfo.Size]/F"	  );
	root->Branch("LepInfo.Eta"	  , &Eta[0]	   , "LepInfo.Eta[LepInfo.Size]/F"	  );
	root->Branch("LepInfo.Phi"	  , &Phi[0]	   , "LepInfo.Phi[LepInfo.Size]/F"	  );
	root->Branch("LepInfo.TrackIso"   , &TrackIso[0]   , "LepInfo.TrackIso[LepInfo.Size]/F"   );
	root->Branch("LepInfo.EcalIso"    , &EcalIso[0]    , "LepInfo.EcalIso[LepInfo.Size]/F"    );   
	root->Branch("LepInfo.HcalIso"    , &HcalIso[0]    , "LepInfo.HcalIso[LepInfo.Size]/F"    ); 
  }  
}; // end of LepInfoBranches

   


class GenInfoBranches {
public:
  int	Size; 
  int   PID[MAX_GENS];
  int   MPID[MAX_GENS];
  float Mass[MAX_GENS];
  float Pt[MAX_GENS];
  float Eta[MAX_GENS];
  float Phi[MAX_GENS];
  int   Trig[MAX_GENS];
  float L1Pt[MAX_GENS];
  float L3Pt[MAX_GENS];

  
  //Candidate* CandRef[MAX_GENS]; // backward pointer to the PAT objects
  void Initialize(){

    Size = 0;
    for(int i=0; i < MAX_GENS; i++){
      
      PID[i] = dummy;
      MPID[i] = dummy;
      Mass[i] = dummy;
      Pt[i]   = dummy;
      Eta[i]  = dummy;
      Phi[i]  = dummy;
      Trig[i]      =0;
      L1Pt[i]      =dummy;
      L3Pt[i]      =dummy;
    
  
    }
  }
  
  void Register(TTree *root) {
	root->Branch("GenInfo.Size"	  , &Size	   , "GenInfo.Size/I"			  );
	root->Branch("GenInfo.PID"	  , &PID[0]	   , "GenInfo.PID[GenInfo.Size]/I"	  );
	root->Branch("GenInfo.MPID"	  , &MPID[0]	   , "GenInfo.MPID[GenInfo.Size]/I"	  );
	root->Branch("GenInfo.Mass"	  , &Mass[0]	   , "GenInfo.Mass[GenInfo.Size]/F"	  );
	root->Branch("GenInfo.Pt"	  , &Pt[0]	   , "GenInfo.Pt[GenInfo.Size]/F"	  );
	root->Branch("GenInfo.Eta"	  , &Eta[0]	   , "GenInfo.Eta[GenInfo.Size]/F"	  );
	root->Branch("GenInfo.Phi"	  , &Phi[0]	   , "GenInfo.Phi[GenInfo.Size]/F"	  );
	root->Branch("GenInfo.Trig"	  , &Trig[0]	   , "GenInfo.Trig[GenInfo.Size]/I"	  );
	root->Branch("GenInfo.L1Pt"	  , &L1Pt[0]	   , "GenInfo.L1Pt[GenInfo.Size]/F"	  );
	root->Branch("GenInfo.L3Pt"	  , &L3Pt[0]	   , "GenInfo.L3Pt[GenInfo.Size]/F"	  );


  }  
}; // end of GenInfoBranches

   

}

#endif
