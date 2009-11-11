#ifndef FORMAT_HH
#define FORMAT_HH

#define MAX_LEPTONS 	64
#define MAX_PHOTONS     64
#define MAX_TRIGGER    120
#define dummy          -9999

namespace TRIGGER{
  enum {    
    HLT_L1SingleEG5                                = 0x1 << 1,                                               
    HLT_L1SingleEG8                                = 0x1 << 2,                                               
    HLT_Photon10_L1R                               = 0x1 << 3,                                              
    HLT_Photon10_LooseEcalIso_TrackIso_L1R         = 0x1 << 4,                                               
    HLT_Photon15_L1R                               = 0x1 << 5,            
    
    HLT_Photon15_TrackIso_L1R                      = 0x1 << 6, 

    HLT_Photon15_LooseEcalIso_L1R                  = 0x1 << 7,

    HLT_Photon20_LooseEcalIso_TrackIso_L1R         = 0x1 << 8, 

    HLT_Photon25_L1R                               = 0x1 << 9, 

    HLT_Mu5                                        = 0x1 << 10                   
  };

}



class EvtInfoBranches {
public:
  int   RunNo;
  int   EvtNo;
  int   HLT;
  
  void Initialize(){
    
    RunNo = EvtNo = HLT=dummy;

  }

  void Register(TTree *root) {
	root->Branch("EvtInfo.RunNo"	 , &RunNo     , "EvtInfo.RunNo/I"     );
	root->Branch("EvtInfo.EvtNo"	 , &EvtNo     , "EvtInfo.EvtNo/I"     );   
	root->Branch("EvtInfo.HLT"	 , &HLT    , "EvtInfo.HLT/I"     );   
  }										    
}; // end of EvtInfoBranches


class PhoInfoBranches {
public:
  int	Size; 
  int   Index[MAX_PHOTONS];
  int   Location[MAX_PHOTONS];
  int   IsLoose[MAX_PHOTONS];
  int   Trig[MAX_PHOTONS];
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
  int   GenPID[MAX_PHOTONS];
  int   GenGMomPID[MAX_PHOTONS];
  int   GenMomPID[MAX_PHOTONS];
  float GenPt[MAX_PHOTONS];
  float GenMomPt[MAX_PHOTONS];


  void Initialize(){

    Size = 0;
    for(int i=0; i< MAX_PHOTONS; i++){

      Index[i]     =dummy;
      IsLoose[i]   =dummy;
      Location[i]  =dummy;
      Trig[i]      =0;
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
      SCE[i]	   =dummy;
      HoverE[i]	   =dummy;
      SCEt[i]	   =dummy;
      SCEta[i]	   =dummy;
      SCPhi[i]	   =dummy;
      SCEtaWidth[i]=dummy;
      SCPhiWidth[i]=dummy;
      SCNCrystal[i]=dummy;
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

   


#endif
