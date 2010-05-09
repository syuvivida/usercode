#ifndef QCDGENISO_HH
#define QCDGENISO_HH
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

#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

// the most important code
#include "syu/MyObjectCounter/header/MyAlg.hh"

using namespace edm;
using namespace reco;
using namespace std;
using namespace math;
using namespace ROOT;


class QCDGenIso : public edm::EDFilter {
public:
  explicit QCDGenIso(const edm::ParameterSet&) ;
  ~QCDGenIso();  

    
private:
  virtual void beginJob() ;
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
  TH1F* fPtHard;
  TH1F* fPtStablePho;
  TH1F* fPtBremPho;
  TH1F* fPtFragPho;
  TH1F* fPtDecayPho;

  // isolation distribution
  TH1F* histIsoDR04_all;
  TH1F* histIsoDR04_brem;
  TH1F* histIsoDR04_frag;
  TH1F* histIsoDR04_decay;

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





//   TGraphAsymmErrors* fEffIsoCone03;


};



#endif

