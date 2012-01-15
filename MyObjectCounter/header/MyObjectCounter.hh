#ifndef MYOBJECTCOUNTER_HH
#define MYOBJECTCOUNTER_HH


#include <TROOT.h>
#include <TSystem.h>
#include <TObject.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TLorentzVector.h>


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "syu/MyObjectCounter/header/MyAlg.hh"
#include "syu/MyObjectCounter/header/format.hh" 


using namespace edm;
using namespace reco;
using namespace std;
using namespace pat;
using namespace math;
using namespace ROOT;
using namespace syu;


class MyObjectCounter : public edm::EDAnalyzer {
public:
  explicit MyObjectCounter(const edm::ParameterSet&) ;
  ~MyObjectCounter();  

    
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) ;
  virtual void endJob() ;
  
  TTree *root;
  EvtInfoBranches  EvtInfo;
  PhoInfoBranches  PhoInfo;
  LepInfoBranches  LepInfo;

  int              _nIn;
  int              _nOut;
  MyAlg            _alg; 

  // histograms
  TH1F* h_leppt;
  TH1F* h_lepeta;
  TH1F* h_get;


};



#endif


