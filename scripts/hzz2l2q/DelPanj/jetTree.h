#ifndef __JET_TREE_H_
#define __JET_TREE_H_





#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DelPanj/TreeMaker/interface/utils.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DelPanj/TreeMaker/interface/baseTree.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
using namespace std;
using namespace edm;



class jetTree  : public baseTree{

 public:
  jetTree(std::string name, TTree* tree, const edm::ParameterSet& cfg);//name=patJetAk05
  ~jetTree();
  void Fill(const edm::Event& iEvent, edm::EventSetup const& iSetup) ; 
  bool passLooseJetID(std::vector<pat::Jet>::const_iterator recJet);
  void SetBranches();
  void Clear();
 
 private:
  jetTree();
  edm::InputTag JetLabel_;
  bool usePFObjects_;
  edm::ParameterSet parSet_;	 

  bool NJet_;
  //Branches common to all the jets.
  std::vector<double> JetPt_;
  std::vector<double> JetEta_;
  std::vector<double> JetPhi_;
  std::vector<double> JetM_;
  std::vector<double> JetEn_;
  std::vector<double> JetBeta_; // only works for H->ZZ->2l2q
  std::vector<double> JetCharMulti_;
  std::vector<double> JetNeutEmEFr_;
  std::vector<double> JetCharHadEFr_;
  std::vector<double> JetNeutHadEFr_;
  std::vector<double> JetCharEmEFr_;
  std::vector<int>    JetPassID_;

};

#endif
