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

#define MAX_VTXS        20
#define MAX_LEPTONS 	64
#define MAX_PHOTONS     64
#define MAX_TRIGGER    120
#define MAX_GENS        32
#define dummy          -9999

namespace syu{


class EvtInfoBranches {
public:
  int   isData;
  int   RunNo;
  int   EvtNo;
  int   HLT;
  
  void Initialize(){
    
    isData= RunNo = EvtNo = HLT = dummy;

  }

  void Register(TTree *root) {
	root->Branch("EvtInfo.isData"	 , &isData    , "EvtInfo.isData/I"    );
	root->Branch("EvtInfo.RunNo"	 , &RunNo     , "EvtInfo.RunNo/I"     );
	root->Branch("EvtInfo.EvtNo"	 , &EvtNo     , "EvtInfo.EvtNo/I"     );   
	root->Branch("EvtInfo.HLT"	 , &HLT       , "EvtInfo.HLT/I"       );   

  }										    
}; // end of EvtInfoBranches


class VtxInfoBranches {

public:
  int Size;
  int Index[MAX_VTXS];
  float pos[MAX_VTXS][3];
  float poserr[MAX_VTXS][3];
  float rho[MAX_VTXS];
  float chi2[MAX_VTXS];
  int   ndof[MAX_VTXS];
  int   ntrks[MAX_VTXS];


  void Initialize(){

    Size = 0;
    for(int i=0; i < MAX_VTXS; i++){
      
      Index[i] = dummy;

      for(int j=0; j < 3; j++){
	pos[i][3]    = dummy;
	poserr[i][3] = dummy;
      }  // loop over x,y,z

      rho[i]       = dummy;
      chi2[i]      = dummy;
      ndof[i]      = dummy;
      ntrks[i]     = dummy;
      
    } // loop over number of vertices
    
  }

  void Register(TTree *root) {
    root->Branch("VtxInfo.Size"	        , &Size	         , "VtxInfo.Size/I"			  );
    root->Branch("VtxInfo.Index"	, Index 	 , "VtxInfo.Index[VtxInfo.Size]/I"	  );
    root->Branch("VtxInfo.pos"          , pos            , "VtxInfo.pos[VtxInfo.Size][3]/F"     );
    root->Branch("VtxInfo.poserr"       , poserr         , "VtxInfo.poserr[VtxInfo.Size][3]/F"  );
    root->Branch("VtxInfo.rho"          , rho            , "VtxInfo.rho[VtxInfo.Size]/F"        );
    root->Branch("VtxInfo.chi2"         , chi2           , "VtxInfo.chi2[VtxInfo.Size]/F"       );
    root->Branch("VtxInfo.ndof"         , ndof           , "VtxInfo.ndof[VtxInfo.Size]/I"       );
    root->Branch("VtxInfo.ntrks"        , ntrks          , "VtxInfo.ntrks[VtxInfo.Size]/I"      );
    
  }  
  
};



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
	root->Branch("PhoInfo.Index"	  , Index   	   , "PhoInfo.Index[PhoInfo.Size]/I"	  );
	root->Branch("PhoInfo.Location"	  , Location       , "PhoInfo.Location[PhoInfo.Size]/I"	  );
	root->Branch("PhoInfo.IsLoose"	  , IsLoose   	   , "PhoInfo.IsLoose[PhoInfo.Size]/I"	  );
	root->Branch("PhoInfo.Trig"	  , Trig   	   , "PhoInfo.Trig[PhoInfo.Size]/I"	  );
	root->Branch("PhoInfo.L1Pt"	  , L1Pt   	   , "PhoInfo.L1Pt[PhoInfo.Size]/F"	  );
	root->Branch("PhoInfo.L3Pt"	  , L3Pt   	   , "PhoInfo.L3Pt[PhoInfo.Size]/F"	  );
	root->Branch("PhoInfo.E"          , E              , "PhoInfo.E[PhoInfo.Size]/F"          );
	root->Branch("PhoInfo.Et"         , Et             , "PhoInfo.Et[PhoInfo.Size]/F"         );
	root->Branch("PhoInfo.Pz"         , Pz             , "PhoInfo.Pz[PhoInfo.Size]/F"         );
	root->Branch("PhoInfo.Eta"        , Eta            , "PhoInfo.Eta[PhoInfo.Size]/F"        );
	root->Branch("PhoInfo.Phi"        , Phi            , "PhoInfo.Phi[PhoInfo.Size]/F"        );
	root->Branch("PhoInfo.R9"         , R9             , "PhoInfo.R9[PhoInfo.Size]/F"         );
	root->Branch("PhoInfo.TrkPtSum"   , TrkPtSum       , "PhoInfo.TrkPtSum[PhoInfo.Size]/F"   );
	root->Branch("PhoInfo.EcalRecHitEtSum", EcalRecHitEtSum, "PhoInfo.EcalRecHitEtSum[PhoInfo.Size]/F");
	root->Branch("PhoInfo.HcalTowerEtSum",  HcalTowerEtSum, "PhoInfo.HcalTowerEtSum[PhoInfo.Size]/F");
	root->Branch("PhoInfo.HoverE"     , HoverE         , "PhoInfo.HoverE[PhoInfo.Size]/F"     );
	root->Branch("PhoInfo.SCE"        , SCE            , "PhoInfo.SCE[PhoInfo.Size]/F"        );
	root->Branch("PhoInfo.SCEt"       , SCEt           , "PhoInfo.SCEt[PhoInfo.Size]/F"       );	
	root->Branch("PhoInfo.SCEta"      , SCEta          , "PhoInfo.SCEta[PhoInfo.Size]/F"      );	
	root->Branch("PhoInfo.SCPhi"      , SCPhi          , "PhoInfo.SCPhi[PhoInfo.Size]/F"      );	
	root->Branch("PhoInfo.SCEtaWidth" , SCEtaWidth     , "PhoInfo.SCEtaWidth[PhoInfo.Size]/F" );	
	root->Branch("PhoInfo.SCPhiWidth" , SCPhiWidth     , "PhoInfo.SCPhiWidth[PhoInfo.Size]/F" );	
	root->Branch("PhoInfo.SCNCrystal" , SCNCrystal     , "PhoInfo.SCNCrystal[PhoInfo.Size]/I" );
	root->Branch("PhoInfo.SigEta"     , SigEta         , "PhoInfo.SigEta[PhoInfo.Size]/F"     );
	root->Branch("PhoInfo.SigIEta"    , SigIEta        , "PhoInfo.SigIEta[PhoInfo.Size]/F"    );
        root->Branch("PhoInfo.GenPID"     , GenPID         , "PhoInfo.GenPID[PhoInfo.Size]/I"     );
	root->Branch("PhoInfo.GenGMomPID" , GenGMomPID     , "PhoInfo.GenGMomPID[PhoInfo.Size]/I" );
	root->Branch("PhoInfo.GenMomPID"  , GenMomPID      , "PhoInfo.GenMomPID[PhoInfo.Size]/I"  );
	root->Branch("PhoInfo.GenPt"      , GenPt          , "PhoInfo.GenPt[PhoInfo.Size]/F"      );
	root->Branch("PhoInfo.GenMomPt"   , GenMomPt       , "PhoInfo.GenMomPt[PhoInfo.Size]/F"   );

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
	root->Branch("LepInfo.Index"	  , Index          , "LepInfo.Index[LepInfo.Size]/I"	  );
	root->Branch("LepInfo.LeptonType" , LeptonType     , "LepInfo.LeptonType[LepInfo.Size]/I" );
	root->Branch("LepInfo.Charge"	  , Charge         , "LepInfo.Charge[LepInfo.Size]/I"	  );
	root->Branch("LepInfo.Pt"	  , Pt   	   , "LepInfo.Pt[LepInfo.Size]/F"	  );
	root->Branch("LepInfo.Eta"	  , Eta  	   , "LepInfo.Eta[LepInfo.Size]/F"	  );
	root->Branch("LepInfo.Phi"	  , Phi 	   , "LepInfo.Phi[LepInfo.Size]/F"	  );
	root->Branch("LepInfo.TrackIso"   , TrackIso       , "LepInfo.TrackIso[LepInfo.Size]/F"   );
	root->Branch("LepInfo.EcalIso"    , EcalIso        , "LepInfo.EcalIso[LepInfo.Size]/F"    );   
	root->Branch("LepInfo.HcalIso"    , HcalIso        , "LepInfo.HcalIso[LepInfo.Size]/F"    ); 
  }  
}; // end of LepInfoBranches

   


class GenInfoBranches {
public:
  float ptHat;
  int	Size; 
  int   Index[MAX_GENS];
  int   GenIndex[MAX_GENS];
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

    ptHat = dummy;
    Size = 0;
    for(int i=0; i < MAX_GENS; i++){
      
      Index[i]     = dummy;
      GenIndex[i]  = dummy;
      PID[i]       = dummy;
      MPID[i]      = dummy;
      Mass[i]      = dummy;
      Pt[i]        = dummy;
      Eta[i]       = dummy;
      Phi[i]       = dummy;
      Trig[i]      = 0;
      L1Pt[i]      = dummy;
      L3Pt[i]      = dummy;
    
  
    }
  }
  
  void Register(TTree *root) {
	root->Branch("GenInfo.ptHat"      , &ptHat         , "GenInfo.ptHat/F"     );
	root->Branch("GenInfo.Size"	  , &Size	   , "GenInfo.Size/I"			  );
	root->Branch("GenInfo.Index"	  , Index 	   , "GenInfo.Index[GenInfo.Size]/I"	  );
	root->Branch("GenInfo.GenIndex"	  , GenIndex 	   , "GenInfo.GenIndex[GenInfo.Size]/I"	  );
	root->Branch("GenInfo.PID"	  , PID 	   , "GenInfo.PID[GenInfo.Size]/I"	  );
	root->Branch("GenInfo.MPID"	  , MPID	   , "GenInfo.MPID[GenInfo.Size]/I"	  );
	root->Branch("GenInfo.Mass"	  , Mass	   , "GenInfo.Mass[GenInfo.Size]/F"	  );
	root->Branch("GenInfo.Pt"	  , Pt  	   , "GenInfo.Pt[GenInfo.Size]/F"	  );
	root->Branch("GenInfo.Eta"	  , Eta 	   , "GenInfo.Eta[GenInfo.Size]/F"	  );
	root->Branch("GenInfo.Phi"	  , Phi 	   , "GenInfo.Phi[GenInfo.Size]/F"	  );
	root->Branch("GenInfo.Trig"	  , Trig 	   , "GenInfo.Trig[GenInfo.Size]/I"	  );
	root->Branch("GenInfo.L1Pt"	  , L1Pt	   , "GenInfo.L1Pt[GenInfo.Size]/F"	  );
	root->Branch("GenInfo.L3Pt"	  , L3Pt	   , "GenInfo.L3Pt[GenInfo.Size]/F"	  );


  }  
}; // end of GenInfoBranches

   

}

#endif
