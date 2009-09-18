// -----------------------------------------------------
// MyObjectCounter.cc 
// -- a very simple example analyzer for PAT tutorial
// ----------------------------------------------------- 
// Shin-Shan Yu


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

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/UtilAlgos/interface/TFileDirectory.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"

#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"


// this file contains the format of lepton, photon, and event structures
#include "format.hh" 
                     

using namespace edm;
using namespace reco;
using namespace std;
using namespace pat;
using namespace math;
using namespace ROOT;


class Histo_struct {

public:
  TH1F* h_leppt;
  TH1F* h_lepeta;
  TH1F* h_get;
  void BookHistograms(edm::Service<TFileService> fo)
  {
    h_leppt  = fo->make<TH1F>("h_leppt","",50,0,200);
    h_lepeta = fo->make<TH1F>("h_lepeta","",24,-6.0,6.0);
    h_get    = fo->make<TH1F>("h_get", "", 50,0,200);
  }

};


class MyObjectCounter : public edm::EDAnalyzer {
public:
  explicit MyObjectCounter(const edm::ParameterSet&) ;
  ~MyObjectCounter();  

    
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void dumpGenInfo(const edm::Event&); 
  virtual void analyze(const edm::Event&, const edm::EventSetup&) ;
  virtual void endJob() ;
  
  TTree *root;
  EvtInfoBranches  EvtInfo;
  PhoInfoBranches  PhoInfo;
  LepInfoBranches  LepInfo;
  Histo_struct     HistoInfo;

  edm::InputTag    _phoLabel;
  edm::InputTag    _muoLabel;
  edm::InputTag    _eleLabel;
  edm::InputTag    _trigger;
  edm::InputTag    _triggerEvent;
  std::string      _matcherName;
  bool             _dumpHEP;

};


MyObjectCounter::MyObjectCounter(const edm::ParameterSet& iConfig):
  _phoLabel(iConfig.getParameter< edm::InputTag >( "phoLabel" ) ),
  _muoLabel(iConfig.getParameter< edm::InputTag >( "muoLabel" ) ),
  _eleLabel(iConfig.getParameter< edm::InputTag >( "eleLabel" ) ),
  _trigger( iConfig.getParameter< edm::InputTag >( "trigger" ) ),
  _triggerEvent( iConfig.getParameter< edm::InputTag >( "triggerEvent" ) ),
  _matcherName( iConfig.getParameter< std::string >( "matcherName" ) )

{  
  _dumpHEP = iConfig.getUntrackedParameter<bool>("dumpHEP", true);
 
}


MyObjectCounter::~MyObjectCounter()
{
}

void MyObjectCounter::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  TFileDirectory results = TFileDirectory( fs->mkdir("MyObjectCounter") );
  HistoInfo.BookHistograms(fs);
  root = new TTree("root","root");
  EvtInfo.Register(root);  
  PhoInfo.Register(root);
  LepInfo.Register(root);


  EvtInfo.Initialize();  
  PhoInfo.Initialize();
  LepInfo.Initialize();
}

void MyObjectCounter::endJob() 
{
}

// dump generator-level information

void MyObjectCounter::dumpGenInfo(const edm::Event& iEvent)
{

  edm::Handle<reco::GenParticleCollection> GenHandle;  
  bool hasGenParticle = iEvent.getByLabel("genParticles", GenHandle);
  if(!hasGenParticle)return;

  std::cout << "============================ " << 
    "Run " <<  iEvent.id().run() << " Evt " <<  iEvent.id().event() <<
    " ============================= " << std::endl << std::endl;

  int genIndex = 0;
  for( std::vector<GenParticle>::const_iterator it_gen = GenHandle->begin(); 
       it_gen != GenHandle->end(); it_gen++ ) {

    printf("%4i", genIndex);
    printf("%7i", it_gen->pdgId());
    printf("%6i", it_gen->status());
    printf("   ");
    printf("%9.3f",it_gen->p4().mass());
    printf("%9.3f",it_gen->energy());
    printf("%9.3f",it_gen->pt());
    printf("\n");
    genIndex ++;
  }

  std::cout << std::endl 
	    << "==========================================================="
            << "===============" 
	    << std::endl;

}



// analyzing reconstructed electrons, muons, and photons

void MyObjectCounter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  if(_dumpHEP)
    dumpGenInfo(iEvent);

  // filter trigger path
  edm::Handle<TriggerResults> TrgResultsHandle;
  bool with_TriggerResults = iEvent.getByLabel(InputTag("TriggerResults::HLT"),TrgResultsHandle);
  if (with_TriggerResults) {
  
    TriggerNames TrgNames( *TrgResultsHandle );   
    unsigned int TrgIdx1 = TrgNames.triggerIndex( "HLT_Photon10_L1R" );
    unsigned int TrgIdx2 = TrgNames.triggerIndex( "HLT_Photon15_L1R" );
	
    if (!TrgResultsHandle->accept(TrgIdx1) &&
	!TrgResultsHandle->accept(TrgIdx2)) return; 	
  }
  else
    std::cout << "TriggerResults::HLT not found" << std::endl;
   
 // PAT trigger information
  edm::Handle< TriggerEvent > triggerEvent;
  bool with_TriggerEvent = iEvent.getByLabel( _triggerEvent, triggerEvent );

  edm::Handle< TriggerObjectCollection > triggerObjects;
  bool with_TriggerObjects = iEvent.getByLabel( _trigger, triggerObjects );


  // Start to fill the main root branches
  // initialize the variables of ntuples
  EvtInfo.Initialize();  

  
  EvtInfo.RunNo  = iEvent.id().run();
  EvtInfo.EvtNo  = iEvent.id().event();
  
  
  // look for photons
  edm::Handle<std::vector<pat::Photon>   > PhotonHandle;
  bool hasPho = iEvent.getByLabel(_phoLabel,   PhotonHandle);

  if(hasPho){
    
    // first do photon trigger matching if information 
    // is available
    if(with_TriggerEvent && with_TriggerObjects)
      {
	// PAT trigger helper for trigger matching information
	const pat::helper::TriggerMatchHelper matchHelper;
	const TriggerObjectMatch * triggerMatch( triggerEvent->triggerObjectMatchResult( _matcherName ) );

	for ( size_t iPhoton = 0; iPhoton < PhotonHandle->size(); ++iPhoton ){	 
	  const reco::CandidateBaseRef candBaseRef( PhotonRef( PhotonHandle, iPhoton ) );
	  const TriggerObjectRef trigRef( matchHelper.triggerMatchObject( candBaseRef, triggerMatch, iEvent, *triggerEvent ) );
	  // fill histograms
	  if ( trigRef.isAvailable() ) { 
	    std::cout << candBaseRef->pt() << "\t" << trigRef->pt() << std::endl;
	  }
	}// loop over pat::photons

      }
    else
      {
	std::cout << "I can not do trigger matching without"
	  "trigger objects or trigger events" << std::endl;
      }
  


    // initialize the variables of ntuples
    PhoInfo.Initialize();


  
    for( std::vector<pat::Photon>::const_iterator it_ph = PhotonHandle->begin(); 
	 it_ph != PhotonHandle->end(); it_ph++ ) {
       
      if (PhoInfo.Size>=MAX_PHOTONS) {
	fprintf(stderr,"ERROR: number of photons exceeds the size of array.\n");
	exit(1);
      }


	
      HistoInfo.h_get->Fill(it_ph->et());
	
      PhoInfo.Index[PhoInfo.Size]        = PhoInfo.Size;  
      PhoInfo.E    [PhoInfo.Size] 	 = it_ph->energy(); 
      PhoInfo.Et   [PhoInfo.Size] 	 = it_ph->et();     
      PhoInfo.Pz   [PhoInfo.Size] 	 = it_ph->pz();     
      PhoInfo.Eta  [PhoInfo.Size] 	 = it_ph->eta();    
      PhoInfo.Phi  [PhoInfo.Size] 	 = it_ph->phi();    
      PhoInfo.R9   [PhoInfo.Size] 	 = it_ph->r9();     
      PhoInfo.TrkIso[PhoInfo.Size]       = it_ph->trackIso();
      PhoInfo.EcalIso[PhoInfo.Size]      = it_ph->ecalIso();
      PhoInfo.HcalIso[PhoInfo.Size]      = it_ph->hcalIso();
      PhoInfo.HoverE[PhoInfo.Size]       = it_ph->hadronicOverEm();
	
      PhoInfo.SCE[PhoInfo.Size]          = it_ph->superCluster()->energy();
      PhoInfo.SCEta[PhoInfo.Size]        = it_ph->superCluster()->eta();
      PhoInfo.SCPhi[PhoInfo.Size]        = it_ph->superCluster()->phi();
      PhoInfo.SCEtaWidth[PhoInfo.Size]   = it_ph->superCluster()->etaWidth();
      PhoInfo.SCPhiWidth[PhoInfo.Size]   = it_ph->superCluster()->phiWidth();
	
      float et = it_ph->superCluster()->energy()/cosh(it_ph->superCluster()->eta());
      PhoInfo.SCEt[PhoInfo.Size]       = et;

      if(it_ph->genPhoton()){
	  
	if(it_ph->genPhoton()->mother()){
	  PhoInfo.GenMomPID[PhoInfo.Size]  = it_ph->genPhoton()->mother()->pdgId();
	  PhoInfo.GenMomPt[PhoInfo.Size]   = it_ph->genPhoton()->mother()->pt();
	  if(it_ph->genPhoton()->mother()->mother())
	    PhoInfo.GenGMomPID[PhoInfo.Size] = it_ph->genPhoton()->mother()->mother()->pdgId();
	}
     
      } // check generator-level matches
	

      PhoInfo.Size++;
    }
  } // if there are photons
  

  
  // look for muons, electrons
  edm::Handle<std::vector<pat::Muon> >     MuonHandle;
  edm::Handle<std::vector<pat::Electron> > ElectronHandle;

  bool hasMuo = iEvent.getByLabel(_muoLabel,     MuonHandle);
  bool hasEle = iEvent.getByLabel(_eleLabel, ElectronHandle);
  
  
  // initialize the variables of ntuples
  if(hasMuo || hasEle){
    LepInfo.Initialize();

  
    for( std::vector<pat::Muon>::const_iterator it_mu = MuonHandle->begin(); 
	 it_mu != MuonHandle->end(); it_mu++ ) {
       
      if (LepInfo.Size>=MAX_LEPTONS) {
	fprintf(stderr,"ERROR: number of leptons exceeds the size of array.\n");
	exit(1);
      }
	

      if (it_mu->pt()<10.) continue;
      if (it_mu->isGlobalMuon()) continue;


      LepInfo.Index      [LepInfo.Size] = LepInfo.Size;
      LepInfo.LeptonType [LepInfo.Size] = 13;
      LepInfo.Charge     [LepInfo.Size] = it_mu->charge();
      LepInfo.Pt         [LepInfo.Size] = it_mu->pt();
      LepInfo.Eta        [LepInfo.Size] = it_mu->eta();
      LepInfo.Phi        [LepInfo.Size] = it_mu->phi();
      LepInfo.TrackIso   [LepInfo.Size] = it_mu->trackIso();
      LepInfo.EcalIso    [LepInfo.Size] = it_mu->ecalIso();
      LepInfo.HcalIso    [LepInfo.Size] = it_mu->hcalIso();	  
      LepInfo.Size++;
    }

  
    for( std::vector<pat::Electron>::const_iterator it_el = ElectronHandle->begin(); 
	 it_el != ElectronHandle->end(); it_el++ ) {
       
      if (LepInfo.Size>=MAX_LEPTONS) {
	fprintf(stderr,"ERROR: number of leptons exceeds the size of array.\n");
	exit(1);
      }
       

      if (it_el->pt()<10.) continue; 
 
      LepInfo.Index      [LepInfo.Size] = LepInfo.Size;
      LepInfo.LeptonType [LepInfo.Size] = 11;
      LepInfo.Charge     [LepInfo.Size] = it_el->charge();
      LepInfo.Pt         [LepInfo.Size] = it_el->pt();
      LepInfo.Eta        [LepInfo.Size] = it_el->eta();
      LepInfo.Phi        [LepInfo.Size] = it_el->phi();
      LepInfo.TrackIso   [LepInfo.Size] = it_el->trackIso();
      LepInfo.EcalIso    [LepInfo.Size] = it_el->ecalIso();
      LepInfo.HcalIso    [LepInfo.Size] = it_el->hcalIso();				
      LepInfo.Size++;

      HistoInfo.h_leppt -> Fill(it_el->pt());
      HistoInfo.h_lepeta-> Fill(it_el->eta());

    }
  } // if there are muons or electrons


  
  root->Fill();
  
}


//define this as a plug-in
DEFINE_FWK_MODULE(MyObjectCounter);


