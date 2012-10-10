#include "DelPanj/TreeMaker/interface/jetTree.h"
#include <CLHEP/Vector/LorentzVector.h>
#include "FWCore/Framework/interface/ESHandle.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "TMath.h"
typedef math::XYZTLorentzVector LorentzVector;

jetTree::jetTree(std::string desc, TTree* tree, const edm::ParameterSet& iConfig):
  baseTree(desc, tree),
  parSet_(iConfig.getParameter<edm::ParameterSet>(desc.c_str()))
{
  JetLabel_ = parSet_.getParameter<edm::InputTag>("Jets");
  SetBranches();
}


jetTree::~jetTree(){
  delete tree_;
}


void
jetTree::Fill(const edm::Event& iEvent, edm::EventSetup const& iSetup){
  Clear();
  
  edm::Handle<std::vector<pat::Jet> > JetHandle;
  if(not iEvent.getByLabel(JetLabel_,JetHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
	     <<JetLabel_<<std::endl; exit(0);}


  const std::vector<pat::Jet>* jets = JetHandle.product();
  std::vector<pat::Jet>::const_iterator jet =jets->begin();

  for(;jet!=jets->end();jet++){

    //Stuff common for all jets.
    JetPt_.push_back(jet->pt());
    JetEta_.push_back(jet->eta());
    JetPhi_.push_back(jet->phi());
    JetM_.push_back(jet->mass());
    JetEn_.push_back(jet->energy());

    if(jet->isPFJet()){      
      JetCharMulti_.push_back(jet->chargedMultiplicity());
      JetNeutEmEFr_.push_back(jet->neutralEmEnergyFraction ());
      JetCharHadEFr_.push_back(jet->chargedHadronEnergyFraction ());
      JetNeutHadEFr_.push_back(jet->neutralHadronEnergyFraction ());
      JetCharEmEFr_.push_back(jet->chargedEmEnergyFraction ());
    }
    else{
      double dummy = -99999.0;
      JetCharMulti_.push_back(dummy);
      JetNeutEmEFr_.push_back(dummy);
      JetCharHadEFr_.push_back(dummy);
      JetNeutHadEFr_.push_back(dummy);
      JetCharEmEFr_.push_back(dummy);
    }
    

  }
}

void
jetTree::SetBranches(){
  
  AddBranch(&JetPt_, "Pt_");
  AddBranch(&JetEta_, "Eta_");
  AddBranch(&JetPhi_, "Phi_");
  AddBranch(&JetM_, "M_");
  AddBranch(&JetEn_, "En_");
  
  AddBranch(&JetCharMulti_, "CharMulti_");
  AddBranch(&JetNeutEmEFr_, "NeutEmEFr_");
  AddBranch(&JetCharHadEFr_, "CharHadEFr_");
  AddBranch(&JetNeutHadEFr_, "NeutHadEFr_");
  AddBranch(&JetCharEmEFr_, "CharEmEFr_");

}


void
jetTree::Clear(){

  JetPt_.clear();
  JetEta_.clear();
  JetPhi_.clear();
  JetM_.clear();
  JetEn_.clear();

  JetCharMulti_.clear();
  JetNeutEmEFr_.clear();
  JetCharHadEFr_.clear();
  JetNeutHadEFr_.clear();
  JetCharEmEFr_.clear();


}
