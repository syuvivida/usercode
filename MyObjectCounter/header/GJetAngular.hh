#ifndef GJETANGULAR_HH
#define GJETANGULAR_HH
// system include files
#include <memory>

#include <TROOT.h>
#include <TSystem.h>
#include <TObject.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>


#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "Math/LorentzVector.h"
#include "Math/GenVector/Boost.h"
// the most important code
#include "syu/MyObjectCounter/header/MyAlg.hh"

using namespace edm;
using namespace reco;
using namespace std;
using namespace math;
using namespace ROOT;


class GJetAngular : public edm::EDFilter {
public:
  explicit GJetAngular(const edm::ParameterSet&) ;
  ~GJetAngular();  

    
private:
  virtual void beginJob() ;
//   virtual void analyze(const edm::Event&, const edm::EventSetup&) ;
  virtual bool filter(edm::Event&, const edm::EventSetup&) ;
  virtual void endJob() ;
  
 
  int              _nIn;
  int              _nOut;
  MyAlg            _alg;

  double           _etaMax;
  double           _ptMin;


  edm::Service<TFileService> fs;

  // sanity check
  TH1F* fPtHat;
  TH1F* fPStar;
  TH1F* fPtHard;
  TH1F* fPtPho;
  TH1F* fPtStablePho;
  
  TH1F* fCosThetaStar;
  TH1F* fTanH;
  TH2F* fPtHatVsCosThetaStar;
  TH2F* fPStarVsCosThetaStar;

  TProfile* pfPStarVsCosThetaStar;
  TProfile* pfPtHatVsCosThetaStar;

};



#endif

