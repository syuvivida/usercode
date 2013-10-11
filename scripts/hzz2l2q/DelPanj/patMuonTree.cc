#include "DelPanj/TreeMaker/interface/patMuonTree.h"

patMuonTree::patMuonTree(std::string name, TTree* tree, const edm::ParameterSet& iConfig):
  mu2012ID_ ( iConfig.getParameter<edm::ParameterSet>("mu2012IDSet")),
  muoRhoIsoInputTag_(iConfig.getParameter<edm::InputTag>("muoRhoIso"))
{
  tree_=tree; 
  patMuonLabel_ = iConfig.getParameter<edm::InputTag>("patMuons");
  SetBranches();
}


patMuonTree::~patMuonTree(){
  delete tree_;

} 


void
patMuonTree::Fill(const edm::Event& iEvent){
  Clear();

  bool isData = iEvent.isRealData();

  // rho for muon
  edm::Handle<double> muo_rho_event;
  iEvent.getByLabel(muoRhoIsoInputTag_,muo_rho_event);
  muoRho_ = *(muo_rho_event.product());

  mu2012ID_.SetData(isData);
  mu2012ID_.SetRho(muoRho_);


  edm::Handle<pat::MuonCollection> patMuonHandle;
  if(not iEvent.getByLabel(patMuonLabel_,patMuonHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
	     <<patMuonLabel_<<std::endl; exit(0);}

  pat::MuonCollection muColl(*(patMuonHandle.product()));
  std::sort(muColl.begin(),muColl.end(),PtGreater());

  nmu=muColl.size();
  pat::MuonCollection::const_iterator mu;
  
  for(mu=muColl.begin(); mu!=muColl.end(); mu++){

    patMuonPt_.push_back(mu->pt());
    patMuonEta_.push_back(mu->eta());
    patMuonPhi_.push_back(mu->phi());
    patMuonM_.push_back(mu->mass());
    patMuonTrkIso_.push_back(mu->trackIso());
    patMuonHcalIso_.push_back(mu->hcalIso());
    patMuonEcalIso_.push_back(mu->ecalIso());
    patMuonCharge_.push_back(mu->charge());
    patMuonNumChambers_.push_back(mu->numberOfChambers());
    patMuonNumMatches_.push_back(mu->numberOfMatches());
    patMuonStationMask_.push_back(mu->stationMask());
    
    bool isTrackMuon =mu->isTrackerMuon();
    bool isGlobalMuon =mu->isGlobalMuon();
    bool isStandAloneMuon =mu->isStandAloneMuon();

    reco::TrackRef muTrk;
    if(isStandAloneMuon)
      muTrk = mu->standAloneMuon();
    else if(isTrackMuon)
      muTrk = mu->track();
    else if(isGlobalMuon)
      muTrk = mu->globalTrack();
      
    patMuonChi2Ndoff_.push_back(muTrk->normalizedChi2());
    patMuonNhits_.push_back(muTrk->numberOfValidHits());

    patMuonDxy_.push_back(mu->dB());
    patMuonDz_.push_back(mu->userFloat("dzVtx"));

    double iso1 = 999.;
    double iso2 = 999.;
    double iso3 = 999.;
    double isoPU = -999.; 
  
    iso1  =  mu->pfIsolationR04().sumChargedHadronPt;
    iso2  =  mu->pfIsolationR04().sumNeutralHadronEt;
    iso3  =  mu->pfIsolationR04().sumPhotonEt;
    isoPU =  mu->pfIsolationR04().sumPUPt;    

    patMuonChHadIso_.push_back(iso1);
    patMuonNeHadIso_.push_back(iso2);
    patMuonGamIso_.push_back(iso3);
    patMuonPUPt_.push_back(isoPU);

    std::map<std::string, bool> Pass = mu2012ID_.CutRecord(*mu);
    int passOrNot = PassAll(Pass);
    patMuonPassID_.push_back(passOrNot);

  
  }
}


void
patMuonTree::AddBranch(std::vector<double>* vec, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);
}


void
patMuonTree::AddBranch(std::vector<int>* vec, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);
}

void 
patMuonTree::AddBranch(double* x, std::string name){
  std::string brName="Electron"+name;
  tree_->Branch(name.c_str(),x,(brName+"/D").c_str());
}

void 
patMuonTree::AddBranch(std::vector<std::string>* vec, std::string name){
  std::string brName="patMuon"+name;
  tree_->Branch(name.c_str(),vec);
}





void 
patMuonTree::SetBranches(){
  AddBranch(&nmu,"NumMu");
  AddBranch(&patMuonPt_, "patMuonPt_");
  AddBranch(&patMuonEta_, "patMuonEta_");
  AddBranch(&patMuonPhi_, "patMuonPhi_");
  AddBranch(&patMuonM_, "patMuonM_");
  AddBranch(&patMuonTrkIso_, "patMuonTrkIso_");
  AddBranch(&patMuonHcalIso_, "patMuonHcalIso_");
  AddBranch(&patMuonEcalIso_, "patMuonEcalIso_");
  AddBranch(&patMuonCharge_, "patMuonCharge_");
  AddBranch(&patMuonNumChambers_, "patMuonNumChambers_");
  AddBranch(&patMuonNumMatches_, "patMuonNumMatches_");
  AddBranch(&patMuonStationMask_, "patMuonStationMask_");
  AddBranch(&patMuonNumSegments_, "patMuonNumSegments_");
  AddBranch(&patMuonChi2Ndoff_, "patMuonChi2Ndoff_");
  AddBranch(&patMuonNhits_, "patMuonNhits_");
  AddBranch(&patMuonDxy_, "patMuonDxy_");
  AddBranch(&patMuonDz_, "patMuonDz_");
  
  AddBranch(&patMuonPassID_, "patMuonPassID_");

  AddBranch(&patMuonChHadIso_, "patMuonChHadIso_");
  AddBranch(&patMuonNeHadIso_, "patMuonNeHadIso_");
  AddBranch(&patMuonGamIso_, "patMuonGamIso_");
  AddBranch(&patMuonPUPt_, "patMuonPUPt_");


}


void
patMuonTree::Clear(){
  patMuonPt_.clear();
  patMuonEta_.clear();
  patMuonPhi_.clear();
  patMuonM_.clear();
  patMuonTrkIso_.clear();
  patMuonHcalIso_.clear();
  patMuonEcalIso_.clear();
  patMuonCharge_.clear();
  patMuonNumChambers_.clear();
  patMuonNumMatches_.clear();
  patMuonStationMask_.clear();
  patMuonNumSegments_.clear();
  patMuonChi2Ndoff_.clear();
  patMuonNhits_.clear();
  patMuonDxy_.clear();
  patMuonDz_.clear();

  patMuonPassID_.clear();

  patMuonChHadIso_.clear();
  patMuonNeHadIso_.clear();
  patMuonGamIso_.clear();
  patMuonPUPt_.clear();


  nmu=0;
}
