#ifndef ZEEFORMAT_HH
#define ZEEFORMAT_HH

#define MAX_GENS        64
#define dummy          -9999

namespace zee{
class EvtInfoBranches {
public:
  int   RunNo;
  int   EvtNo;
  int   HLT;
  float RecZMass;
  float GenZMass;
  
  void Initialize(){
    
    RunNo = EvtNo = HLT = dummy;
    RecZMass = GenZMass = dummy;

  }

  void Register(TTree *root) {
	root->Branch("EvtInfo.RunNo"	 , &RunNo     , "EvtInfo.RunNo/I"     );
	root->Branch("EvtInfo.EvtNo"	 , &EvtNo     , "EvtInfo.EvtNo/I"     );   
	root->Branch("EvtInfo.HLT"	 , &HLT    , "EvtInfo.HLT/I"     ); 
	root->Branch("EvtInfo.RecZMass"  , &RecZMass, "EvtInfo.RecZMass/F");
	root->Branch("EvtInfo.GenZMass"  , &GenZMass, "EvtInfo.GenZMass/F");
  }										    
}; // end of EvtInfoBranches


class PhoInfoBranches {
public:
  int   Location[2];
  int   IsLoose[2];
  int   Trig[2];
  float L1Pt[2];
  float L3Pt[2];
  float E[2];
  float Et[2];
  float Pz[2];
  float Eta[2];
  float Phi[2];
  float R9[2];
  float TrkPtSum[2];
  float EcalRecHitEtSum[2];
  float HcalTowerEtSum[2];
  float HoverE[2];
  float SCE[2];
  float SCEt[2];
  float SCEta[2];
  float SCPhi[2];
  float SCEtaWidth[2];
  float SCPhiWidth[2];
  int   SCNCrystal[2];
  float SigEta[2];
  float SigIEta[2];
  int   EleClass[2];
  int   GenPID[2];
  int   GenGMomPID[2];
  int   GenMomPID[2];
  float GenPt[2];
  float GenMomPt[2];
  

  void Initialize(){

    for(int i=0; i< 2; i++){

      IsLoose[i]   =dummy;
      Location[i]  =dummy;
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
      SCE[i]	   =dummy;
      HoverE[i]	   =dummy;
      SCEt[i]	   =dummy;
      SCEta[i]	   =dummy;
      SCPhi[i]	   =dummy;
      SCEtaWidth[i]=dummy;
      SCPhiWidth[i]=dummy;
      SCNCrystal[i]=dummy;
      SigEta[i]    =dummy;
      SigIEta[i]   =dummy;
      EleClass[i]  =dummy;
      GenPID[i]    =dummy;
      GenGMomPID[i]=dummy;
      GenMomPID[i] =dummy;
      GenPt[i]     =dummy;
      GenMomPt[i]  =dummy;

    }
  }
  
  //Candidate* CandRef[2]; // backward pointer to the PAT objects
  
  void Register(TTree *root) {
	root->Branch("PhoInfo.Location"	  , &Location[0]   , "PhoInfo.Location[2]/I"	  );
	root->Branch("PhoInfo.IsLoose"	  , &IsLoose[0]	   , "PhoInfo.IsLoose[2]/I"	  );
	root->Branch("PhoInfo.Trig"	  , &Trig[0]	   , "PhoInfo.Trig[2]/I"	  );
	root->Branch("PhoInfo.L1Pt"	  , &L1Pt[0]	   , "PhoInfo.L1Pt[2]/F"	  );
	root->Branch("PhoInfo.L3Pt"	  , &L3Pt[0]	   , "PhoInfo.L3Pt[2]/F"	  );
	root->Branch("PhoInfo.E"          , &E[0]          , "PhoInfo.E[2]/F"          );
	root->Branch("PhoInfo.Et"         , &Et[0]         , "PhoInfo.Et[2]/F"         );
	root->Branch("PhoInfo.Pz"         , &Pz[0]         , "PhoInfo.Pz[2]/F"         );
	root->Branch("PhoInfo.Eta"        , &Eta[0]        , "PhoInfo.Eta[2]/F"        );
	root->Branch("PhoInfo.Phi"        , &Phi[0]        , "PhoInfo.Phi[2]/F"        );
	root->Branch("PhoInfo.R9"         , &R9[0]         , "PhoInfo.R9[2]/F"         );
	root->Branch("PhoInfo.TrkPtSum"   , &TrkPtSum[0]   , "PhoInfo.TrkPtSum[2]/F"   );
	root->Branch("PhoInfo.EcalRecHitEtSum", &EcalRecHitEtSum[0], "PhoInfo.EcalRecHitEtSum[2]/F");
	root->Branch("PhoInfo.HcalTowerEtSum",  &HcalTowerEtSum[0], "PhoInfo.HcalTowerEtSum[2]/F");
	root->Branch("PhoInfo.HoverE"     , &HoverE[0]     , "PhoInfo.HoverE[2]/F"     );
	root->Branch("PhoInfo.SCE"        , &SCE[0]        , "PhoInfo.SCE[2]/F");
	root->Branch("PhoInfo.SCEt"       , &SCEt[0]       , "PhoInfo.SCEt[2]/F");	
	root->Branch("PhoInfo.SCEta"      , &SCEta[0]      , "PhoInfo.SCEta[2]/F");	
	root->Branch("PhoInfo.SCPhi"      , &SCPhi[0]      , "PhoInfo.SCPhi[2]/F");	
	root->Branch("PhoInfo.SCEtaWidth" , &SCEtaWidth[0] , "PhoInfo.SCEtaWidth[2]/F");	
	root->Branch("PhoInfo.SCPhiWidth" , &SCPhiWidth[0] , "PhoInfo.SCPhiWidth[2]/F");	
	root->Branch("PhoInfo.SCNCrystal",  &SCNCrystal[0] , "PhoInfo.SCNCrystal[2]/I");
	root->Branch("PhoInfo.SigEta",      &SigEta[0]     , "PhoInfo.SigEta[2]/F");
	root->Branch("PhoInfo.SigIEta",     &SigIEta[0]    , "PhoInfo.SigIEta[2]/F");
	root->Branch("PhoInfo.EleClass",    &EleClass[0]   , "PhoInfo.EleClass[2]/I");
        root->Branch("PhoInfo.GenPID",      &GenPID[0]     , "PhoInfo.GenPID[2]/I");
	root->Branch("PhoInfo.GenGMomPID" , &GenGMomPID[0] , "PhoInfo.GenGMomPID[2]/I" );
	root->Branch("PhoInfo.GenMomPID"  , &GenMomPID[0]  , "PhoInfo.GenMomPID[2]/I" );
	root->Branch("PhoInfo.GenPt"   , &GenPt[0]   , "PhoInfo.GenPt[2]/F" );
	root->Branch("PhoInfo.GenMomPt"   , &GenMomPt[0]   , "PhoInfo.GenMomPt[2]/F" );

  }
};  // end of PhoInfoBranches


class GenInfoBranches {
public:
  int	Size; 
  int   PID[MAX_GENS];
  int   MPID[MAX_GENS];
  float Mass[MAX_GENS];
  float Pt[MAX_GENS];
  float Eta[MAX_GENS];
  float Phi[MAX_GENS];

  
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
  }  
}; // end of GenInfoBranches

   

}

#endif
