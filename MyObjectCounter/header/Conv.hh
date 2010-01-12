#ifndef CONV_HH
#define CONV_HH

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

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/UtilAlgos/interface/TFileDirectory.h"

#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

// the most important code
#include "syu/MyObjectCounter/header/MyAlg.hh"
#include "syu/MyObjectCounter/header/format.hh" 
#include "RecoEgamma/EgammaMCTools/interface/ElectronMCTruth.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruth.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruthFinder.h"

using namespace edm;
using namespace reco;
using namespace std;
using namespace trigger;
using namespace math;
using namespace ROOT;
using namespace syu;



class Conv : public edm::EDAnalyzer {
public:
  explicit Conv(const edm::ParameterSet&) ;
  ~Conv();  

    
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) ;
  virtual void endJob() ;
  
  TTree *root;
  syu::EvtInfoBranches  EvtInfo;
  syu::PhoInfoBranches  PhoInfo;
  syu::GenInfoBranches  GenInfo;

  partEtMap<reco::Photon>::Type  myPhoEtMap;
  partGenMap<reco::Photon>::Type myPhoGenMap;
  partEtMap<reco::GsfElectron>::Type myEleEtMap;
  std::vector<reco::GenParticleCollection::const_iterator> myHardGen;
  std::vector<reco::GenParticleCollection::const_iterator> myConvPho;

  PhotonMCTruthFinder*  thePhotonMCTruthFinder_;

 
  int              _nIn;
  int              _nOut;
  MyAlg            _alg;
  int              _pdgCode;
  
  TH1F*            heta_all;
  TH1F*            heta_original;
  TProfile*        heta_counte;
  TProfile*        heta_conv;

  TH1F*            hpt_all;
  TH1F*            hpt_original;
  TProfile*        hpt_counte;
  TProfile*        hpt_conv;

  TH1F*            he_original;
  TProfile*        he_counte;
  TProfile*        he_conv;


};



#endif

