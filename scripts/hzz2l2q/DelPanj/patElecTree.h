#ifndef __ELEC_TREE_H_
#define __ELEC_TREE_H_


/*

Author: Anil P Singh
 Panjab University

*/

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DelPanj/TreeMaker/interface/utils.h"
#include "DelPanj/TreeMaker/interface/eSelector.h"

using namespace std;
using namespace edm;

class patElecTree{


 public:
  patElecTree(std::string name, TTree* tree, const edm::ParameterSet& cfg);
  ~patElecTree();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void Clear();

 private:
//TTree* tree_;
  //Dont Allow User to Call the Default Constructor.
  patElecTree();


  void AddBranch(double* x, std::string name);
  void AddBranch(std::vector<int>*, std::string name);
  void AddBranch(std::vector<double>*, std::string name);
  void AddBranch(std::vector<std::string>*, std::string name);
  edm::InputTag dcsTag_;
  bool isData_;

  edm::InputTag patElecLabel_;

  edm::InputTag eleRhoIsoInputTag_;
  eSelector e2012ID_;
  double eleRho_;

  int Nelecs;
  std::vector<int> patElecType;
 //pt, eta, phi, M : Enough to caluclate
  //px,py,pz etc.==>No need to store later
  TTree* tree_;
  double nele; 
  std::vector<double> patElecEt_;
  std::vector<double> patElecEnergy_;
  std::vector<double> patElecPt_;
  std::vector<double> patElecEta_;
  std::vector<double> patElecPhi_;
  std::vector<double> patElecM_;
  std::vector<double> patElecScEta_;
  std::vector<double> patElecSigIhIh_;
  std::vector<double> patElecDelEtaIn_;
  std::vector<double> patElecDelPhiIn_;
  std::vector<double> patElecHoE_;
  std::vector<double> patElecTrkIso_;
  std::vector<double> patElecHcalIso_;
  std::vector<double> patElecEcalIso_;
  std::vector<double> patElecEoverP_;
  std::vector<double> patElecDxy_;
  std::vector<double> patElecDz_;
  std::vector<double> patElecChHadIso_;
  std::vector<double> patElecNeHadIso_;
  std::vector<double> patElecGamIso_;

  // conversion rejection
  std::vector<double> patElecMissingHits_ ;
  std::vector<double> patElecDist_ ;
  std::vector<double> patElecDeltaCotTheta_ ;
  std::vector<int>    patElecHasConv_;
  std::vector<double> patElecInBarrel_;
  std::vector<double> patElecInEndcap_;

  std::vector<int>    patElecPassID_;
  std::vector<int>    patElecPassVBTFID_;

};

#endif

