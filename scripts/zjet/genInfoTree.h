#ifndef _GEN_INFO_TREE_H_
#define _GEN_INFO_TREE_H_

#include <memory>
#include <vector>
#include <string>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"

//
// class declaration
//

class genInfoTree{

 public:
  genInfoTree(std::string name, TTree* tree, const edm::ParameterSet& cfg);
  ~genInfoTree();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void Clear();
  
  
  // ----------member data ---------------------------
  void AddBranch(double* x, std::string name);
  void AddBranch(std::vector<double>*, std::string name);
  void AddBranch(std::vector<int>*, std::string name);
  
  edm::InputTag genPartLabel_;
  
  
  TTree* tree_;

  double ptHat_;      // added by Eiko
  double mcWeight_;   // added by Eiko

  std::vector<double> genParE_;
  std::vector<double> genParPt_;
  std::vector<double> genParEta_;
  std::vector<double> genParPhi_;
  std::vector<int>    genParSt_;
  std::vector<int>    genParQ_;
  std::vector<int>    genParId_;
  std::vector<int>    genMomParId_; // added by Eiko
  std::vector<int>    genParIndex_;    // added by Eiko

  std::vector<double> genLepE_;
  std::vector<double> genLepPt_;
  std::vector<double> genLepEta_;
  std::vector<double> genLepPhi_;
  std::vector<double> genLepPhoE_;
  std::vector<double> genLepPhoPt_;
  std::vector<double> genLepPhoEta_;
  std::vector<double> genLepPhoPhi_;
  std::vector<int>    genLepId_;

  std::vector<double> genJetE_;
  std::vector<double> genJetPt_;
  std::vector<double> genJetEta_; 
  std::vector<double> genJetPhi_; 
};

#endif
