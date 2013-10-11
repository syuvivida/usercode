#ifndef __ZZ_TREE_HH_
#define __ZZ_TREE_HH_

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include "TTree.h"
#include "TLorentzVector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DelPanj/TreeMaker/interface/utils.h"
#include "DelPanj/TreeMaker/interface/muSelector.h"
#include "DelPanj/TreeMaker/interface/eSelector.h"
#include "pandolf/BTagSFUtil/src/BTagSFUtil.h" 

#include "DataFormats/PatCandidates/interface/Jet.h"


using namespace std;
using namespace edm;
using namespace reco;

class ZZTree{


 public:
  ZZTree(std::string name, TTree* tree, const edm::ParameterSet& cfg);
  ~ZZTree();
  void Fill(const edm::Event& iEvent, edm::EventSetup const& iSetup);
  void SetBranches();
  void Clear();

 private:
  ZZTree();

  BTagSFUtil* btsfutiljp;
  void AddBranch(double* x, std::string name);
  void AddBranch(int* x, std::string name);
  void AddBranch(std::vector<double>*, std::string name);
  void AddBranch(std::vector<int>*, std::string name);
  void AddBranchArray(const int arraySize, double* x, std::string name);

  bool passLooseJetID(const pat::Jet* recJet);

  edm::InputTag hzzeejj_;
  edm::InputTag hzzmmjj_;
  edm::InputTag eleRhoIsoInputTag_;
  edm::InputTag muoRhoIsoInputTag_;
  bool merged_;

  eSelector e2012ID_;
  muSelector mu2012ID_;

  TTree* tree_;

  double eleRho_;
  double muoRho_;
  double metSig_;
  unsigned int nJets_;
 
  std::vector<double> higgsPt_;
  std::vector<double> higgsEta_;
  std::vector<double> higgsPhi_;
  std::vector<double> higgsM_;
  std::vector<double> higgsMRefit_;

  std::vector<double> zllPt_;
  std::vector<double> zllEta_;
  std::vector<double> zllPhi_;
  std::vector<double> zllM_;
  std::vector<double> zlldR_; // deltaR between two leptons

  std::vector<double> zjjPt_;
  std::vector<double> zjjEta_;
  std::vector<double> zjjPhi_;
  std::vector<double> zjjM_;
  std::vector<double> zjjMRefit_;
  std::vector<double> zjjdR_; // deltaR between two jets   

  std::vector<int>    jetIndex_;
  std::vector<int>    jetHiggsIndex_;
  std::vector<double> jetE_;
  std::vector<double> jetPt_;
  std::vector<double> jetEta_;
  std::vector<double> jetPhi_;
  std::vector<double> jetTau1_;
  std::vector<double> jetTau2_;
  std::vector<double> jetQV_;
  
  std::vector<double> heliLD_;
  std::vector<double> heliLDRefit_;

  std::vector<int>    nBTags_;
  std::vector<int>    lepType_;
  std::vector<int>    passBit_;

  /*
  std::vector<double> eID01;
  std::vector<double> eID02;
  std::vector<double> eID03;
  std::vector<double> eID04;
  std::vector<double> eID05;
  std::vector<double> eID06;
  std::vector<double> eID07;
  std::vector<double> eID08;
  std::vector<double> eID09;
  std::vector<double> eID10;
  std::vector<double> eID11;
  std::vector<double> eID12;
  std::vector<double> eID13;
  std::vector<double> eID14;
  std::vector<double> eID15;
  std::vector<double> eID16;
  std::vector<double> eID17;

  std::vector<double> muID01;
  std::vector<double> muID02;
  std::vector<double> muID03;
  std::vector<double> muID04;
  std::vector<double> muID05;
  std::vector<double> muID06;
  std::vector<double> muID07;
  std::vector<double> muID08;
  std::vector<double> muID09;
  std::vector<double> muID10;
  std::vector<double> muID11;
  std::vector<double> muID12;
  std::vector<double> muID13;
  std::vector<double> muID14;
  std::vector<double> muID15;
  */

};

#endif

