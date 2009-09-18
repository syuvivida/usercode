#define MAX_LEPTONS 	64
#define MAX_PHOTONS     64
#define dummy          -9999


class EvtInfoBranches {
public:
  int   RunNo;
  int   EvtNo;
  
  void Register(TTree *root) {
	root->SetBranchAddress("EvtInfo.RunNo"     , &RunNo     );
	root->SetBranchAddress("EvtInfo.EvtNo"     , &EvtNo     );     
  }  										    
};

class PhoInfoBranches {
public:
  int	Size; 
  int   Index[MAX_PHOTONS];
  float E[MAX_PHOTONS];
  float Et[MAX_PHOTONS];
  float Pz[MAX_PHOTONS];
  float Eta[MAX_PHOTONS];
  float Phi[MAX_PHOTONS];
  float R9[MAX_PHOTONS];
  float TrkIso[MAX_PHOTONS];
  float EcalIso[MAX_PHOTONS];
  float HcalIso[MAX_PHOTONS];
  float HoverE[MAX_PHOTONS];
  float SCE[MAX_PHOTONS];
  float SCEt[MAX_PHOTONS];
  float SCEta[MAX_PHOTONS];
  float SCPhi[MAX_PHOTONS];
  float SCEtaWidth[MAX_PHOTONS];
  float SCPhiWidth[MAX_PHOTONS];
  int   GenGMomPID[MAX_PHOTONS];
  int   GenMomPID[MAX_PHOTONS];
  float GenMomPt[MAX_PHOTONS];

  
  //Candidate* CandRef[MAX_PHOTONS]; // backward pointer to the PAT objects
  
  void Register(TTree *root) {
	root->SetBranchAddress("PhoInfo.Size"	  , &Size	     );
	root->SetBranchAddress("PhoInfo.Index"	  , &Index[0]	     );
	root->SetBranchAddress("PhoInfo.E"          , &E[0]            );
	root->SetBranchAddress("PhoInfo.Et"         , &Et[0]           );
	root->SetBranchAddress("PhoInfo.Pz"         , &Pz[0]           );
	root->SetBranchAddress("PhoInfo.Eta"        , &Eta[0]          );
	root->SetBranchAddress("PhoInfo.Phi"        , &Phi[0]          );
	root->SetBranchAddress("PhoInfo.R9"         , &R9[0]           );
	root->SetBranchAddress("PhoInfo.TrkIso"     , &TrkIso[0]       );
	root->SetBranchAddress("PhoInfo.EcalIso"    , &EcalIso[0]      );
	root->SetBranchAddress("PhoInfo.HcalIso"    , &HcalIso[0]      );
	root->SetBranchAddress("PhoInfo.HoverE"     , &HoverE[0]       );
	root->SetBranchAddress("PhoInfo.SCE"        , &SCE[0]          );
	root->SetBranchAddress("PhoInfo.SCEt"       , &SCEt[0]         );
	root->SetBranchAddress("PhoInfo.SCEta"      , &SCEta[0]        );
	root->SetBranchAddress("PhoInfo.SCPhi"      , &SCPhi[0]        );
	root->SetBranchAddress("PhoInfo.SCEtaWidth" , &SCEtaWidth[0]   );	
	root->SetBranchAddress("PhoInfo.SCPhiWidth" , &SCPhiWidth[0]   );	
	root->SetBranchAddress("PhoInfo.GenGMomPID" , &GenGMomPID[0]   );
	root->SetBranchAddress("PhoInfo.GenMomPID"  , &GenMomPID[0]    );
	root->SetBranchAddress("PhoInfo.GenMomPt"   , &GenMomPt[0]     );

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


  void Register(TTree *root) {
	root->SetBranchAddress("LepInfo.Size"	    , &Size	     );
	root->SetBranchAddress("LepInfo.Index"      , &Index[0]      );
	root->SetBranchAddress("LepInfo.LeptonType" , &LeptonType[0] );
	root->SetBranchAddress("LepInfo.Charge"     , &Charge[0]     );
	root->SetBranchAddress("LepInfo.Pt"	    , &Pt[0]	     );
	root->SetBranchAddress("LepInfo.Eta"	    , &Eta[0]	     );
	root->SetBranchAddress("LepInfo.Phi"	    , &Phi[0]	     );
	root->SetBranchAddress("LepInfo.TrackIso"   , &TrackIso[0]   );
	root->SetBranchAddress("LepInfo.EcalIso"    , &EcalIso[0]    );   
	root->SetBranchAddress("LepInfo.HcalIso"    , &HcalIso[0]    ); 

  }    
};


