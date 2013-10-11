#ifndef __MUON_TREE_H_
#define __MUON_TREE_H_

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DelPanj/TreeMaker/interface/utils.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DelPanj/TreeMaker/interface/muSelector.h"

using namespace std;
using namespace edm;

class patMuonTree{


 public:
  patMuonTree(std::string name, TTree* tree, const edm::ParameterSet& cfg);
  ~patMuonTree();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void Clear();

 private:
//TTree* tree_;
  //Dont Allow User to Call the Default Constructor.
  patMuonTree();


  void AddBranch(double* x, std::string name);
  void AddBranch(std::vector<int>*, std::string name);
  void AddBranch(std::vector<double>*, std::string name);
  void AddBranch(std::vector<std::string>*, std::string name);

  edm::InputTag patMuonLabel_;
  edm::InputTag muoRhoIsoInputTag_;
  muSelector mu2012ID_;
  double muoRho_;

  int Nmuons;
  std::vector<int> patMuonType;
 //pt, eta, phi, M : Enough to caluclate
  //px,py,pz etc.==>No need to store later
  TTree* tree_;
  double nmu; 
  std::vector<std::string> patMuonType_; 
  std::vector<double> patMuonPt_;
  std::vector<double> patMuonEta_;
  std::vector<double> patMuonPhi_;
  std::vector<double> patMuonM_;
  std::vector<double> patMuonSumTrkPt_;
  std::vector<double> patMuonTrkIso_;
  std::vector<double> patMuonHcalIso_;
  std::vector<double> patMuonEcalIso_;
  std::vector<double> patMuonCharge_;
  std::vector<double> patMuonNumChambers_;
  std::vector<double> patMuonNumMatches_;
  std::vector<double> patMuonStationMask_;
  std::vector<double> patMuonNumSegments_;
  std::vector<double> patMuonChi2Ndoff_;
  std::vector<double> patMuonNhits_;
  std::vector<double> patMuonDxy_;
  std::vector<double> patMuonDz_;


  std::vector<double> patMuonChHadIso_;
  std::vector<double> patMuonNeHadIso_;
  std::vector<double> patMuonGamIso_;
  std::vector<double> patMuonPUPt_;

  std::vector<int> patMuonPassID_;

};

#endif

