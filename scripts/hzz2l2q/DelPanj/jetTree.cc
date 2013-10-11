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
  int pass=0;
  int njets=0;
  for(;jet!=jets->end();jet++){
    njets++;
    //Stuff common for all jets.
    JetPt_.push_back(jet->pt());
    JetEta_.push_back(jet->eta());
    JetPhi_.push_back(jet->phi());
    JetM_.push_back(jet->mass());
    JetEn_.push_back(jet->energy());
    JetBeta_.push_back(jet->userFloat("puBeta"));

    
    bool passID = passLooseJetID(jet);
    if(passID) pass = 1;

    if(jet->isPFJet()){      
      JetPassID_.push_back(pass);
      JetCharMulti_.push_back(jet->chargedMultiplicity());
      JetNeutEmEFr_.push_back(jet->neutralEmEnergyFraction ());
      JetCharHadEFr_.push_back(jet->chargedHadronEnergyFraction ());
      JetNeutHadEFr_.push_back(jet->neutralHadronEnergyFraction ());
      JetCharEmEFr_.push_back(jet->chargedEmEnergyFraction ());
    }
    else{
      double dummy = -99999.0;
      JetPassID_.push_back(0);
      JetCharMulti_.push_back(dummy);
      JetNeutEmEFr_.push_back(dummy);
      JetCharHadEFr_.push_back(dummy);
      JetNeutHadEFr_.push_back(dummy);
      JetCharEmEFr_.push_back(dummy);
    }
    

  } // end of loop over jets

  if(iEvent.id().event()==150075191 || iEvent.id().event()==778689)
    cout << "Number of jets in event " << iEvent.id().event() << "= " << njets << endl;

}


bool jetTree::passLooseJetID(std::vector<pat::Jet>::const_iterator recjet)
{
  double eta = recjet->eta();
  if(recjet->getPFConstituents().size() <= 1)return false;                                                                               
  if(recjet->neutralHadronEnergyFraction() >= 0.99)return false;
  if(recjet->neutralEmEnergyFraction() >= 0.99)return false;
  //   // for the tracker region
  if(fabs(eta)<2.4 && recjet->chargedHadronEnergyFraction()<= 0.0)return false;
  if(fabs(eta)<2.4 && recjet->chargedEmEnergyFraction() >= 0.99)return false;
  if(fabs(eta)<2.4 && recjet->chargedMultiplicity() <= 0)return false;
  return true;

}



void
jetTree::SetBranches(){
  
  AddBranch(&JetPt_, "Pt_");
  AddBranch(&JetEta_, "Eta_");
  AddBranch(&JetPhi_, "Phi_");
  AddBranch(&JetM_, "M_");
  AddBranch(&JetEn_, "En_");
  AddBranch(&JetBeta_, "Beta_");
  AddBranch(&JetPassID_, "PassID_");
  
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
  JetBeta_.clear();
  JetPassID_.clear();

  JetCharMulti_.clear();
  JetNeutEmEFr_.clear();
  JetCharHadEFr_.clear();
  JetNeutHadEFr_.clear();
  JetCharEmEFr_.clear();


}
