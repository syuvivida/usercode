#ifndef GENISO_HH
#define GENISO_HH

#include <TROOT.h>
#include <TSystem.h>
#include <TObject.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TLorentzVector.h>


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
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


class GenIso : public edm::EDAnalyzer {
public:
  explicit GenIso(const edm::ParameterSet&) ;
  ~GenIso();  

    
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) ;
  virtual void endJob() ;
  
 
  int              _nIn;
  int              _nOut;
  MyAlg            _alg;
  int              _pdgCode;
  double           _etaMax;
  double           _ptMin;

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


};



#endif

