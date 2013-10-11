#include "DelPanj/TreeMaker/interface/patElecTree.h"

#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
// for conversion finder
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "RecoEgamma/EgammaTools/interface/ConversionInfo.h"




/*
  NOTE: For the users of particle flow.

  Really a large number of issues have to be 
  explained for the case of the particle flow
  electrons. This module at the moment treats
  the electrons as if they were the non particle flow
  objects. We will make all our selections similar to 
  the regular objects. Particle flow electrons need 
  deeper study.
*/



patElecTree::patElecTree(std::string name, TTree* tree, const edm::ParameterSet& iConfig):
  e2012ID_ ( iConfig.getParameter<edm::ParameterSet>("e2012IDSet")),
  eleRhoIsoInputTag_(iConfig.getParameter<edm::InputTag>("eleRhoIso"))
{
  tree_=tree; 
  patElecLabel_ = iConfig.getParameter<edm::InputTag>("patElectrons");
  SetBranches();
}


patElecTree::~patElecTree(){
  delete tree_;

} 



void
patElecTree::Fill(const edm::Event& iEvent){
  Clear();
  edm::Handle<pat::ElectronCollection> patElecHandle;
  if(not iEvent.getByLabel(patElecLabel_,patElecHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
	     <<patElecLabel_<<std::endl; exit(0);}


  // rho for electron
  edm::Handle<double> ele_rho_event;
  iEvent.getByLabel(eleRhoIsoInputTag_,ele_rho_event);
  eleRho_ = *(ele_rho_event.product());

  bool isData = iEvent.isRealData();

  e2012ID_.SetData(isData);
  e2012ID_.SetRho(eleRho_);

  pat::ElectronCollection eleColl(*(patElecHandle.product()));
  std::sort(eleColl.begin(),eleColl.end(),PtGreater());

  nele=eleColl.size();
  pat::ElectronCollection::const_iterator ele;
  
  for(ele=eleColl.begin(); ele!=eleColl.end(); ele++){


    patElecEt_.push_back(ele->et());
    patElecEnergy_.push_back(ele->energy());
    patElecPt_.push_back(ele->p4().pt());
    patElecEta_.push_back(ele->eta());
    patElecPhi_.push_back(ele->phi());
    patElecM_.push_back(ele->mass());

    double supercluster_eta =-999;    
    supercluster_eta = ele->superCluster()->eta();
    patElecScEta_.push_back(supercluster_eta);

     
     
    //     double supercluster_e =-999;
    double sigihih = -999;
    double deletain = -999;
    double delphiin = -999;
    double hoe = -999;


    sigihih          = ele->sigmaIetaIeta();
    delphiin         = ele->deltaPhiSuperClusterTrackAtVtx();
    deletain         = ele->deltaEtaSuperClusterTrackAtVtx();
    hoe              = ele->hadronicOverEm();
    
    //     //Don't know how this gonna fare in case of track-only PFElectrons.
    patElecSigIhIh_.push_back(sigihih);
    patElecDelEtaIn_.push_back(deletain);
    patElecDelPhiIn_.push_back(delphiin);
    patElecHoE_.push_back(hoe);

    patElecTrkIso_.push_back(ele->dr03TkSumPt());
    patElecHcalIso_.push_back(ele->dr03HcalTowerSumEt());
    patElecEcalIso_.push_back(ele->dr03EcalRecHitSumEt());
    double ooemoop = fabs(
			  1.0/std::max((double)1e-3,(double)ele->ecalEnergy()) - 
			  1.0/std::max((double)1e-3,(double)sqrt(ele->trackMomentumAtVtx().mag2()))
			  );

    patElecEoverP_.push_back(ooemoop);
    patElecDxy_.push_back(ele->userFloat("dxy"));
    patElecDz_.push_back(ele->userFloat("dz"));
    patElecMissingHits_.push_back( ele->gsfTrack().get()->trackerExpectedHitsInner().numberOfHits() );
    patElecDist_.push_back(ele->convDist() );
    patElecDeltaCotTheta_.push_back(ele->convDcot() ) ;

    patElecChHadIso_.push_back(ele->chargedHadronIso());
    patElecNeHadIso_.push_back(ele->neutralHadronIso());
    patElecGamIso_.push_back(ele->photonIso());
    patElecInBarrel_.push_back(ele->isEB());
    patElecInEndcap_.push_back(ele->isEE());
    

    std::map<std::string, bool> Pass    = e2012ID_.CutRecord(*ele); 
    int passOrNot = PassAll(Pass);
    patElecPassID_.push_back(passOrNot);
    patElecPassVBTFID_.push_back(ele->electronID("eidVBTFCom95"));
    patElecHasConv_.push_back(ele->userFloat("hasMatchConv"));
  }
}


void
patElecTree::AddBranch(std::vector<double>* vec, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);
}

void
patElecTree::AddBranch(std::vector<int>* vec, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);
}

void 
patElecTree::AddBranch(double* x, std::string name){
  std::string brName="Electron"+name;
  tree_->Branch(name.c_str(),x,(brName+"/D").c_str());
}

void 
patElecTree::AddBranch(std::vector<std::string>* vec, std::string name){
  std::string brName="patElec"+name;
  tree_->Branch(name.c_str(),vec);
}





void 
patElecTree::SetBranches(){
  AddBranch(&patElecEt_, "patElecEt_");
  AddBranch(&patElecEnergy_, "patElecEnergy_");
  AddBranch(&patElecPt_, "patElecPt_");
  AddBranch(&patElecEta_, "patElecEta_");
  AddBranch(&patElecPhi_, "patElecPhi_");
  AddBranch(&patElecM_, "patElecM_");
  AddBranch(&patElecScEta_, "patElecScEta_");
  AddBranch(&patElecSigIhIh_, "patElecSigIhIh_");
  AddBranch(&patElecDelEtaIn_, "patElecDelEtaIn_");
  AddBranch(&patElecDelPhiIn_, "patElecDelPhiIn_");
  AddBranch(&patElecHoE_, "patElecHoE_");
  AddBranch(&patElecTrkIso_, "patElecTrkIso_");
  AddBranch(&patElecHcalIso_, "patElecHcalIso_");
  AddBranch(&patElecEcalIso_, "patElecEcalIso_");
  AddBranch(&patElecEoverP_, "patElecEoverP_");

  AddBranch(&patElecDxy_, "patElecDxy_");
  AddBranch(&patElecDz_, "patElecDz_");

  AddBranch(&patElecChHadIso_, "patElecChHadIso_");
  AddBranch(&patElecNeHadIso_, "patElecNeHadIso_");
  AddBranch(&patElecGamIso_, "patElecGamIso_");
  
  // conversion rejection
  AddBranch(&patElecMissingHits_, "patElecMissingHits_");
  AddBranch(&patElecDist_, "patElecDist_" );
  AddBranch(&patElecDeltaCotTheta_, "patElecDeltaCotTheta_");
  AddBranch(&patElecInBarrel_,"patElecInBarrel_");
  AddBranch(&patElecInEndcap_,"patElecInEndcap_");

  AddBranch(&patElecHasConv_,"patElecHasConv_");
  AddBranch(&patElecPassID_, "patElecPassID_");

  AddBranch(&patElecPassVBTFID_, "patElecPassVBTFID_");


}


void
patElecTree::Clear(){
  patElecEt_.clear();
  patElecEnergy_.clear();
  patElecPt_.clear();
  patElecEta_.clear();
  patElecPhi_.clear();
  patElecM_.clear();
  patElecScEta_.clear();
  patElecSigIhIh_.clear();
  patElecDelEtaIn_.clear();
  patElecDelPhiIn_.clear();
  patElecHoE_.clear();
  patElecTrkIso_.clear();
  patElecHcalIso_.clear();
  patElecEcalIso_.clear();
  patElecEoverP_.clear();
  patElecDxy_.clear();
  patElecDz_.clear();
  patElecChHadIso_.clear();
  patElecNeHadIso_.clear();
  patElecGamIso_.clear();
  patElecMissingHits_.clear();
  patElecDist_.clear();
  patElecDeltaCotTheta_.clear();
  
  patElecHasConv_.clear();
  patElecInBarrel_.clear();
  patElecInEndcap_.clear();

  patElecPassID_.clear();
  patElecPassVBTFID_.clear();

}

















