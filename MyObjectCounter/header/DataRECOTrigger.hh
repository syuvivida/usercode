#ifndef DATARECOTRIGGER_HH
#define DATARECOTRIGGER_HH

#include <TROOT.h>
#include <TSystem.h>
#include <TObject.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TProfile.h>
#include <TLorentzVector.h>


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"


#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

// the most important code
#include "syu/MyObjectCounter/header/MyAlg.hh"
#include "syu/MyObjectCounter/header/format.hh" 

using namespace edm;
using namespace reco;
using namespace std;
using namespace trigger;
using namespace math;
using namespace ROOT;
using namespace syu;


class DataRECOTrigger : public edm::EDAnalyzer {
public:
  explicit DataRECOTrigger(const edm::ParameterSet&) ;
  ~DataRECOTrigger();  

    
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) ;
  virtual void endJob() ;
  
  TTree *root;
  syu::EvtInfoBranches  EvtInfo;
  syu::PhoInfoBranches  PhoInfo;
  syu::VtxInfoBranches  VtxInfo;

  partEtMap<reco::Photon>::Type  myPhoEtMap;
  partGenMap<reco::Photon>::Type myPhoGenMap;
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
  int              _pdgCode;


  // histograms


  TH1F* h_getdeno;
  TH1F* h_getnumr;
  TH1F* h_recgetdeno;
  TH1F* h_recgetnumr;
  TH1F* h_gengetdeno;
  TH1F* h_gengetnumr;
  TH1F* h_jetgetdeno;
  TH1F* h_jetgetnumr;
  TH1F* h_qgetdeno;
  TH1F* h_qgetnumr;
  TH1F* h_ggetdeno;
  TH1F* h_ggetnumr;
  TH1F* h_debug1;
  TH1F* h_debug2;
  TH1F* h_eta1;
  TH1F* h_eta2;
  TH1F* h_mugetdeno;
  TH1F* h_mugetnumr;


};



#endif

