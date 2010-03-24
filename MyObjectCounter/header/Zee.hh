#ifndef ZEE_HH
#define ZEE_HH


#include <TROOT.h>
#include <TSystem.h>
#include <TObject.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"

// the most important code
#include "syu/MyObjectCounter/header/MyAlg.hh"
#include "syu/MyObjectCounter/header/zeeformat.hh" 


using namespace edm;
using namespace std;
using namespace reco;
using namespace trigger;
using namespace math;
using namespace ROOT;
using namespace zee;


class Zee : public edm::EDAnalyzer {
public:
  explicit Zee(const edm::ParameterSet&) ;
  ~Zee();  

    
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) ;
  virtual void endJob() ;

  TTree *root;
  EvtInfoBranches  EvtInfo;
  PhoInfoBranches  PhoInfo;
  GenInfoBranches  GenInfo;

  partEtMap<reco::GsfElectron>::Type myEleEtMap;
  partEtMap<reco::Photon>::Type      myPhoEtMap;
  elePhoMap<reco::GsfElectron, reco::Photon>::Type 
                                 myElePhoMap;
  partGenMap<reco::GsfElectron>::Type myEleGenMap;
  partL1Map<reco::Photon>::Type  myPhoL1Map;
  partL3Map<reco::Photon>::Type  myPhoL3Map_HLTL1EG5;
  partL3Map<reco::Photon>::Type  myPhoL3Map_HLTL1EG8;
  partL3Map<reco::Photon>::Type  myPhoL3Map_HLT10;
  partL3Map<reco::Photon>::Type  myPhoL3Map_HLT15;
  partL3Map<reco::Photon>::Type  myPhoL3Map_HLT20Iso;
  partL3Map<reco::Photon>::Type  myPhoL3Map_HLT25;
  partL3Map<reco::Photon>::Type  myPhoL3Map_HLTMu5;

  int              _nIn;
  int              _nOut;
  MyAlg            _alg;


  // histogram
  TH1F* h_trkisodeno;
  TH1F* h_trkisonumr;
  TH1F* h_sigietadeno;
  TH1F* h_sigietanumr;
  TH2F* h_xy;
  TH2F* h_rz;
  TH1F* h_dE;



};


#endif

