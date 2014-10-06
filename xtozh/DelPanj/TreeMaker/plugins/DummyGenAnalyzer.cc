#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <iostream>
#include <string>
#include <TH1F.h>
#include <TFile.h>
using namespace std;
using namespace edm;



class DummyGenAnalyzer : public EDAnalyzer {
private: 

  TFile * output;
  TH1F* h_zpt;
  TH1F* h_higgspt;
  TH1F* h_xpt;
  TH1F* h_zy;
  TH1F* h_higgsy;
  TH1F* h_xy;


public:
  explicit DummyGenAnalyzer( const edm::ParameterSet & cfg ) : 
    fileName_(cfg.getUntrackedParameter<std::string>("histoutputFile"))
  {
  }
  void beginJob(){
    output = new TFile(fileName_.data(), "RECREATE");
    
    h_zpt   = new TH1F("h_zpt","",1000,0,1000);
    h_zpt->SetXTitle("p_{T}(Z) [GeV]");
    h_zpt->Sumw2();

    h_zy    = new TH1F("h_zy","",60,-3.0,3.0); 
    h_zy->SetXTitle("y_{Z}");
    h_zy->Sumw2();

    h_higgspt   = new TH1F("h_higgspt","",1000,0,1000);
    h_higgspt->SetXTitle("p_{T}(h) [GeV]");
    h_higgspt->Sumw2();

    h_higgsy    = new TH1F("h_higgsy","",60,-3.0,3.0); 
    h_higgsy->SetXTitle("y_{h}");
    h_higgsy->Sumw2();


    h_xpt   = new TH1F("h_xpt","",1000,0,1000);
    h_xpt->SetXTitle("p_{T}(X) [GeV]");
    h_xpt->Sumw2();

    h_xy    = new TH1F("h_xy","",60,-3.0,3.0); 
    h_xy->SetXTitle("y_{X}");
    h_xy->Sumw2();


  }

  ~DummyGenAnalyzer(){
    delete output;
  }

private:
  void analyze( const Event & iEvent, const EventSetup & iSetup ) {


    using namespace edm;
    edm::Handle<reco::GenParticleCollection> genParticleHandle;
    if(not iEvent.getByLabel("genParticles", genParticleHandle))
    {
      std::cout<<
	"GenAnalyzer: Generator Level Information not found\n"
	       <<std::endl;
    }
    unsigned int genIndex=0;
    const reco::GenParticleCollection* genColl= &(*genParticleHandle);

    // now loop
    reco::GenParticleCollection::const_iterator geni = genColl->begin();
    for(; geni!=genColl->end();geni++){ 

      reco::GenParticle gen = *geni;
      int PID = gen.pdgId();
      int status= gen.status();
    if( gen.pdgId()!=9000001 && gen.pdgId()!=1023 &&
	(gen.status()!=3) && 
     	(gen.status()<21 || gen.status()>29))continue;

    genIndex++;

    int momPID = gen.numberOfMothers() ==1?
      gen.mother()->pdgId():-1;

    if(PID==23 && (momPID==1023 || momPID==9000001) && ( (status>=21 && status<=29) || status==3))
      {
	h_zpt->Fill(gen.pt());
	h_zy->Fill(gen.rapidity());
      }
    else if(PID==25 && (momPID==1023 || momPID==9000001) && ((status>=21 && status<=29) || status==3))
      {
	h_higgspt->Fill(gen.pt());
	h_higgsy->Fill(gen.rapidity());

      }

    else if((PID==9000001  || PID==1023) && (status==62 || status==3))
      {
	h_xpt->Fill(gen.pt());
	h_xy->Fill(gen.rapidity());
      }  
    } // end of loop over particles

    

  }    
      
  void endJob(){
    h_zpt->Write();
    h_higgspt->Write();
    h_xpt->Write();
    h_zy->Write();
    h_higgsy->Write();
    h_xy->Write();
    output->Write();
    output->Close();
  }   

  void beginRun(edm::Run const& iRun, edm::EventSetup const& es){


  }


  std::string fileName_;
};

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( DummyGenAnalyzer );


