#ifndef ZEEFORMAT_HH
#define ZEEFORMAT_HH

#include <TROOT.h>
#include <TSystem.h>
#include <TObject.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>


#define MAX_GENS        32
#define MAX_ZPAIR       10
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
  int   ZPairSize;
  int   Location[MAX_ZPAIR][2];
  int   IsLoose[MAX_ZPAIR][2];
  int   Trig[MAX_ZPAIR][2];
  float L1Pt[MAX_ZPAIR][2];
  float L3Pt[MAX_ZPAIR][2];
  float E[MAX_ZPAIR][2];
  float Et[MAX_ZPAIR][2];
  float Pz[MAX_ZPAIR][2];
  float Eta[MAX_ZPAIR][2];
  float Phi[MAX_ZPAIR][2];
  float R9[MAX_ZPAIR][2];
  float TrkPtSum[MAX_ZPAIR][2];
  float EcalRecHitEtSum[MAX_ZPAIR][2];
  float HcalTowerEtSum[MAX_ZPAIR][2];
  float HoverE[MAX_ZPAIR][2];
  float SCE[MAX_ZPAIR][2];
  float SCEt[MAX_ZPAIR][2];
  float SCEta[MAX_ZPAIR][2];
  float SCPhi[MAX_ZPAIR][2];
  float SCEtaWidth[MAX_ZPAIR][2];
  float SCPhiWidth[MAX_ZPAIR][2];
  int   SCNCrystal[MAX_ZPAIR][2];
  float SigEta[MAX_ZPAIR][2];
  float SigIEta[MAX_ZPAIR][2];

  float eE[MAX_ZPAIR][2];
  float eEt[MAX_ZPAIR][2];
  float ePz[MAX_ZPAIR][2];
  float eEta[MAX_ZPAIR][2];
  float ePhi[MAX_ZPAIR][2];
  float eR9[MAX_ZPAIR][2];
  float eTrkPtSum[MAX_ZPAIR][2];
  float eEcalRecHitEtSum[MAX_ZPAIR][2];
  float eHcalTowerEtSum[MAX_ZPAIR][2];
  float eHoverE[MAX_ZPAIR][2];
  float eSCE[MAX_ZPAIR][2];
  float eSCEt[MAX_ZPAIR][2];
  float eSCEta[MAX_ZPAIR][2];
  float eSCPhi[MAX_ZPAIR][2];
  float eSCEtaWidth[MAX_ZPAIR][2];
  float eSCPhiWidth[MAX_ZPAIR][2];
  int   eSCNCrystal[MAX_ZPAIR][2];
  float eSigEta[MAX_ZPAIR][2];
  float eSigIEta[MAX_ZPAIR][2];
  int   eCharge[MAX_ZPAIR][2];
  int   eClass[MAX_ZPAIR][2];

  float genPt[MAX_ZPAIR][2];
  float genEta[MAX_ZPAIR][2];
  float genPhi[MAX_ZPAIR][2];
  

  void Initialize(){

    ZPairSize = 0;
    for(int i=0; i<MAX_ZPAIR;i++)
      for(int j=0; j< 2; j++){

      IsLoose[i][j]        =dummy;
      Location[i][j]       =dummy;
      Trig[i][j]           =0;
      L1Pt[i][j]           =dummy;
      L3Pt[i][j]           =dummy;
      E[i][j]    	   =dummy;
      Et[i][j]   	   =dummy;
      Pz[i][j]   	   =dummy;
      Eta[i][j]  	   =dummy;
      Phi[i][j]  	   =dummy;
      R9[i][j]   	   =dummy;
      TrkPtSum[i][j]       =dummy;
      EcalRecHitEtSum[i][j]=dummy;
      HcalTowerEtSum[i][j] =dummy;
      SCE[i][j]	           =dummy;
      HoverE[i][j]	   =dummy;
      SCEt[i][j]	   =dummy;
      SCEta[i][j]	   =dummy;
      SCPhi[i][j]	   =dummy;
      SCEtaWidth[i][j]     =dummy;
      SCPhiWidth[i][j]     =dummy;
      SCNCrystal[i][j]     =dummy;
      SigEta[i][j]         =dummy;
      SigIEta[i][j]        =dummy;

      eE[i][j]             =dummy;
      eEt[i][j]            =dummy;
      ePz[i][j]            =dummy;
      eEta[i][j]           =dummy;
      ePhi[i][j]           =dummy;
      eR9[i][j]            =dummy;
      eTrkPtSum[i][j]      =dummy;
      eEcalRecHitEtSum[i][j] =dummy;
      eHcalTowerEtSum[i][j]=dummy;
      eHoverE[i][j]        =dummy;
      eSCE[i][j]           =dummy;
      eSCEt[i][j]          =dummy;
      eSCEta[i][j]         =dummy;
      eSCPhi[i][j]         =dummy;
      eSCEtaWidth[i][j]    =dummy;
      eSCPhiWidth[i][j]    =dummy;
      eSCNCrystal[i][j]    =dummy;
      eSigEta[i][j]        =dummy;
      eSigIEta[i][j]       =dummy;
      eCharge[i][j]        =dummy;
      eClass[i][j]         =dummy;

      genPt[i][j]          =dummy;
      genEta[i][j]         =dummy;
      genPhi[i][j]         =dummy;

    }
  }
  
  
  void Register(TTree *root) {
    root->Branch("PhoInfo.ZPairSize"    , &ZPairSize	   , "PhoInfo.ZPairSize/I");
    root->Branch("PhoInfo.Location"	  , Location   , "PhoInfo.Location[PhoInfo.ZPairSize][2]/I"	  );
    root->Branch("PhoInfo.IsLoose"	  , IsLoose	   , "PhoInfo.IsLoose[PhoInfo.ZPairSize][2]/I"	  );
    root->Branch("PhoInfo.Trig"	  , Trig	   , "PhoInfo.Trig[PhoInfo.ZPairSize][2]/I"	  );
    root->Branch("PhoInfo.L1Pt"	  , L1Pt	   , "PhoInfo.L1Pt[PhoInfo.ZPairSize][2]/F"	  );
    root->Branch("PhoInfo.L3Pt"	  , L3Pt	   , "PhoInfo.L3Pt[PhoInfo.ZPairSize][2]/F"	  );
    root->Branch("PhoInfo.E"          , E          , "PhoInfo.E[PhoInfo.ZPairSize][2]/F"          );
    root->Branch("PhoInfo.Et"         , Et         , "PhoInfo.Et[PhoInfo.ZPairSize][2]/F"         );
    root->Branch("PhoInfo.Pz"         , Pz         , "PhoInfo.Pz[PhoInfo.ZPairSize][2]/F"         );
    root->Branch("PhoInfo.Eta"        , Eta        , "PhoInfo.Eta[PhoInfo.ZPairSize][2]/F"        );
    root->Branch("PhoInfo.Phi"        , Phi        , "PhoInfo.Phi[PhoInfo.ZPairSize][2]/F"        );
    root->Branch("PhoInfo.R9"         , R9         , "PhoInfo.R9[PhoInfo.ZPairSize][2]/F"         );
    root->Branch("PhoInfo.TrkPtSum"   , TrkPtSum   , "PhoInfo.TrkPtSum[PhoInfo.ZPairSize][2]/F"   );
    root->Branch("PhoInfo.EcalRecHitEtSum", EcalRecHitEtSum, "PhoInfo.EcalRecHitEtSum[PhoInfo.ZPairSize][2]/F");
    root->Branch("PhoInfo.HcalTowerEtSum",  HcalTowerEtSum, "PhoInfo.HcalTowerEtSum[PhoInfo.ZPairSize][2]/F");
    root->Branch("PhoInfo.HoverE"     , HoverE     , "PhoInfo.HoverE[PhoInfo.ZPairSize][2]/F"     );
    root->Branch("PhoInfo.SCE"        , SCE        , "PhoInfo.SCE[PhoInfo.ZPairSize][2]/F");
    root->Branch("PhoInfo.SCEt"       , SCEt       , "PhoInfo.SCEt[PhoInfo.ZPairSize][2]/F");	
    root->Branch("PhoInfo.SCEta"      , SCEta      , "PhoInfo.SCEta[PhoInfo.ZPairSize][2]/F");	
    root->Branch("PhoInfo.SCPhi"      , SCPhi      , "PhoInfo.SCPhi[PhoInfo.ZPairSize][2]/F");	
    root->Branch("PhoInfo.SCEtaWidth" , SCEtaWidth , "PhoInfo.SCEtaWidth[PhoInfo.ZPairSize][2]/F");	
    root->Branch("PhoInfo.SCPhiWidth" , SCPhiWidth , "PhoInfo.SCPhiWidth[PhoInfo.ZPairSize][2]/F");	
    root->Branch("PhoInfo.SCNCrystal",  SCNCrystal , "PhoInfo.SCNCrystal[PhoInfo.ZPairSize][2]/I");
    root->Branch("PhoInfo.SigEta",      SigEta     , "PhoInfo.SigEta[PhoInfo.ZPairSize][2]/F");
    root->Branch("PhoInfo.SigIEta",     SigIEta    , "PhoInfo.SigIEta[PhoInfo.ZPairSize][2]/F");

    root->Branch("PhoInfo.eE"          , eE          , "PhoInfo.eE[PhoInfo.ZPairSize][2]/F"          );
    root->Branch("PhoInfo.eEt"         , eEt         , "PhoInfo.eEt[PhoInfo.ZPairSize][2]/F"         );
    root->Branch("PhoInfo.ePz"         , ePz         , "PhoInfo.ePz[PhoInfo.ZPairSize][2]/F"         );
    root->Branch("PhoInfo.eEta"        , eEta        , "PhoInfo.eEta[PhoInfo.ZPairSize][2]/F"        );
    root->Branch("PhoInfo.ePhi"        , ePhi        , "PhoInfo.ePhi[PhoInfo.ZPairSize][2]/F"        );
    root->Branch("PhoInfo.eR9"         , eR9         , "PhoInfo.eR9[PhoInfo.ZPairSize][2]/F"         );
    root->Branch("PhoInfo.eTrkPtSum"   , eTrkPtSum   , "PhoInfo.eTrkPtSum[PhoInfo.ZPairSize][2]/F"   );
    root->Branch("PhoInfo.eEcalRecHitEtSum", eEcalRecHitEtSum, "PhoInfo.eEcalRecHitEtSum[PhoInfo.ZPairSize][2]/F");
    root->Branch("PhoInfo.eHcalTowerEtSum",  eHcalTowerEtSum, "PhoInfo.eHcalTowerEtSum[PhoInfo.ZPairSize][2]/F");
    root->Branch("PhoInfo.eHoverE"     , eHoverE     , "PhoInfo.eHoverE[PhoInfo.ZPairSize][2]/F"     );
    root->Branch("PhoInfo.eSCE"        , eSCE        , "PhoInfo.eSCE[PhoInfo.ZPairSize][2]/F");
    root->Branch("PhoInfo.eSCEt"       , eSCEt       , "PhoInfo.eSCEt[PhoInfo.ZPairSize][2]/F");	
    root->Branch("PhoInfo.eSCEta"      , eSCEta      , "PhoInfo.eSCEta[PhoInfo.ZPairSize][2]/F");	
    root->Branch("PhoInfo.eSCPhi"      , eSCPhi      , "PhoInfo.eSCPhi[PhoInfo.ZPairSize][2]/F");	
    root->Branch("PhoInfo.eSCEtaWidth" , eSCEtaWidth , "PhoInfo.eSCEtaWidth[PhoInfo.ZPairSize][2]/F");	
    root->Branch("PhoInfo.eSCPhiWidth" , eSCPhiWidth , "PhoInfo.eSCPhiWidth[PhoInfo.ZPairSize][2]/F");	
    root->Branch("PhoInfo.eSCNCrystal",  eSCNCrystal , "PhoInfo.eSCNCrystal[PhoInfo.ZPairSize][2]/I");
    root->Branch("PhoInfo.eSigEta",      eSigEta     , "PhoInfo.eSigEta[PhoInfo.ZPairSize][2]/F");
    root->Branch("PhoInfo.eSigIEta",     eSigIEta    , "PhoInfo.eSigIEta[PhoInfo.ZPairSize][2]/F");
    root->Branch("PhoInfo.eClass",    eClass   , "PhoInfo.eClass[PhoInfo.ZPairSize][2]/I");

    root->Branch("PhoInfo.genPt",      genPt    , "PhoInfo.genPt[PhoInfo.ZPairSize][2]/F");
    root->Branch("PhoInfo.genEta",     genEta    , "PhoInfo.genEta[PhoInfo.ZPairSize][2]/F");
    root->Branch("PhoInfo.genPhi",     genPhi    , "PhoInfo.genPhi[PhoInfo.ZPairSize][2]/F");
   

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
      
      PID[i]  =dummy;
      MPID[i] =dummy;
      Mass[i] =dummy;
      Pt[i]   =dummy;
      Eta[i]  =dummy;
      Phi[i]  =dummy;
       
  
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
