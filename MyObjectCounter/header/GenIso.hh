#ifndef GENISO_HH
#define GENISO_HH
// system include files
#include <memory>

#include <TROOT.h>
#include <TSystem.h>
#include <TObject.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
// #include <TGraphAsymmErrors.h>
#include <TLorentzVector.h>


#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/UtilAlgos/interface/TFileDirectory.h"


#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

// the most important code
#include "syu/MyObjectCounter/header/MyAlg.hh"

using namespace edm;
using namespace reco;
using namespace std;
using namespace math;
using namespace ROOT;


class GenIso : public edm::EDFilter {
public:
  explicit GenIso(const edm::ParameterSet&) ;
  ~GenIso();  

    
private:
  virtual void beginJob(const edm::EventSetup&) ;
//   virtual void analyze(const edm::Event&, const edm::EventSetup&) ;
  virtual bool filter(edm::Event&, const edm::EventSetup&) ;
  virtual void endJob() ;
  
 
  int              _nIn;
  int              _nOut;
  MyAlg            _alg;
  int              _pdgCode;
  double           _etaMax;
  double           _ptMin;
  double           _isoMaxCut;

  edm::Service<TFileService> fs;

  // sanity check
  TH1F* fPtHat;

  // for isolation
  TH1F* fIsoPt00;
  TH1F* fIsoPt01;
  TH1F* fEffIso;

  TH1F* fTrkIsoPt00;
  TH1F* fTrkIsoPt01;
  TH1F* fEffTrkIso;

  TH1F* fIsoEta00;
  TH1F* fIsoEta01;
  TH1F* fEffIsoEta;

  TH1F* fIsoCone0300;
  TH1F* fIsoCone0301;
  TH1F* fEffIsoCone03;

  // separate in 3 cases
  // for photons with pt that match to pthat
  TH1F* fIsoPt00_match;
  TH1F* fIsoPt01_match;
  TH1F* fEffIso_match;
  TH1F* fIsoEta00_match;
  TH1F* fIsoEta01_match;
  TH1F* fEffIsoEta_match;

  // for photons with pt that are 50% lower than the pthat
  TH1F* fIsoPt00_lopt;
  TH1F* fIsoPt01_lopt;
  TH1F* fEffIso_lopt;
  TH1F* fIsoEta00_lopt;
  TH1F* fIsoEta01_lopt;
  TH1F* fEffIsoEta_lopt;

  // for photons with pt that are 50% higher than the pthat
  TH1F* fIsoPt00_hipt;
  TH1F* fIsoPt01_hipt;
  TH1F* fEffIso_hipt;
  TH1F* fIsoEta00_hipt;
  TH1F* fIsoEta01_hipt;
  TH1F* fEffIsoEta_hipt;


  TH1F* fIsoPt00_mismatch;
  TH1F* fIsoPt01_mismatch;
  TH1F* fEffIso_mismatch;
  TH1F* fIsoEta00_mismatch;
  TH1F* fIsoEta01_mismatch;
  TH1F* fEffIsoEta_mismatch;


  // check angular separation between photon and parton
  
  TH1F* dphigj_match;
  TH1F* dphigj_lopt;    // photon
  TH1F* dphigj_hipt;    // photon
  TH1F* dphigj_mismatch; // photon with status = 1 has lower pt than status = 3

  TH1F* detagj_match;
  TH1F* detagj_lopt;    // photon
  TH1F* detagj_hipt;    // photon
  TH1F* detagj_mismatch; // photon with status = 1 has lower pt than status = 3

  TH1F* drgj_match;
  TH1F* drgj_lopt;    // photon                            
  TH1F* drgj_hipt;    // photon                                                                      
  TH1F* drgj_mismatch; // photon with status = 1 has lower pt than status = 3                                   


//   TGraphAsymmErrors* fEffIsoCone03;


};



#endif

