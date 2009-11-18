// -----------------------------------------------------
// TriggerBitCounter.cc 
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
using namespace syu;


class TriggerBitCounter : public edm::EDAnalyzer {
public:
  explicit TriggerBitCounter(const edm::ParameterSet&) ;
  ~TriggerBitCounter();  

    
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void dumpGenInfo(const edm::Event&); 
  virtual void analyze(const edm::Event&, const edm::EventSetup&) ;
  virtual void endJob() ;
  bool IsLoosePhoton(reco::Photon it_ph);
  
  TTree *root;
  EvtInfoBranches  EvtInfo;
  PhoInfoBranches  PhoInfo;
  LepInfoBranches  LepInfo;

  edm::InputTag    _phoLabel;
  edm::InputTag    _muoLabel;
  edm::InputTag    _eleLabel;
  edm::InputTag    _trigger;
  edm::InputTag    _triggerEvent;
  std::string      _matcherName;
  bool             _dumpHEP;
  int              _nIn;
  int              _nOut;

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


TriggerBitCounter::TriggerBitCounter(const edm::ParameterSet& iConfig):
  _phoLabel(iConfig.getParameter< edm::InputTag >( "phoLabel" ) ),
  _muoLabel(iConfig.getParameter< edm::InputTag >( "muoLabel" ) ),
  _eleLabel(iConfig.getParameter< edm::InputTag >( "eleLabel" ) ),
  _trigger( iConfig.getParameter< edm::InputTag >( "trigger" ) ),
  _triggerEvent( iConfig.getParameter< edm::InputTag >( "triggerEvent" ) ),
  _matcherName( iConfig.getParameter< std::string >( "matcherName" ) ),
  _nIn(0), _nOut(0)

{  
  _dumpHEP = iConfig.getUntrackedParameter<bool>("dumpHEP", true);
 
}


TriggerBitCounter::~TriggerBitCounter()
{
}

void TriggerBitCounter::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  TFileDirectory results = TFileDirectory( fs->mkdir("TriggerBitCounter") );

  int nbin=100;
  double xmin=0.0;
  double xmax=200.0;

  h_eta1  = fs->make<TH1F>("h_eta1","#eta of photons "
			   "before loose photon cuts", 60,-3.0,3.0);
  h_eta2  = fs->make<TH1F>("h_eta2","#eta of photons "
			   "after loose photon cuts", 60,-3.0,3.0);

  h_debug1  = fs->make<TH1F>("h_debug1","Reconstructed photon "
			     "Et before loose photon cuts", nbin,xmin,xmax);
  h_debug2 = fs->make<TH1F>("h_debug2","Reconstructed photon "
			    "Et before loose photon cuts", nbin,xmin,xmax);
  h_getdeno = fs->make<TH1F>("h_getdeno","Reconstructed and matched photon "
			     "Et before photon trigger cuts", nbin,xmin,xmax);
  h_getnumr = fs->make<TH1F>("h_getnumr","Reconstructed and matched photon "
			     "Et after photon trigger cuts", nbin,xmin,xmax);
  h_recgetdeno = fs->make<TH1F>("h_recgetdeno","Reconstructed photon Et "
				"before photon trigger cuts", nbin,xmin,xmax);
  h_recgetnumr = fs->make<TH1F>("h_recgetnumr","Reconstructed photon Et "
				"after photon trigger cuts", nbin,xmin,xmax);
  h_gengetdeno = fs->make<TH1F>("h_gengetdeno","Generated photon Et before "
				"photon trigger cuts", nbin,xmin,xmax);
  h_gengetnumr = fs->make<TH1F>("h_gengetnumr","Generated photon Et after "
				"photon trigger cuts", nbin,xmin,xmax);
  h_jetgetdeno = fs->make<TH1F>("h_jetgetdeno","Reconstructed and jet photon"
				" Et before photon trigger cuts", 
				nbin,xmin,xmax);
  h_jetgetnumr = fs->make<TH1F>("h_jetgetnumr","Reconstructed and jet photon"
				" Et after photon trigger cuts", 
				nbin,xmin,xmax);
  h_qgetdeno   = fs->make<TH1F>("h_qgetdeno","Reconstructed and quark photon"
				" Et before photon trigger cuts", 
				nbin,xmin,xmax);
  h_qgetnumr   = fs->make<TH1F>("h_qgetnumr","Reconstructed and quark photon"
				" Et after photon trigger cuts", 
				nbin,xmin,xmax);
  h_ggetdeno   = fs->make<TH1F>("h_ggetdeno","Reconstructed and gluon photon"
				" Et before photon trigger cuts", 
				nbin,xmin,xmax);
  h_ggetnumr   = fs->make<TH1F>("h_ggetnumr","Reconstructed and gluon photon"
				" Et after photon trigger cuts", 
				nbin,xmin,xmax);

  h_mugetdeno = fs->make<TH1F>("h_mugetdeno","Reconstructed photon Et "
			       "before photon trigger cuts", nbin,xmin,xmax);
  h_mugetnumr = fs->make<TH1F>("h_mugetnumr","Reconstructed photon Et "
			       "after photon trigger cuts", nbin,xmin,xmax);

  root = new TTree("root","root");
  EvtInfo.Register(root);  
  PhoInfo.Register(root);
  LepInfo.Register(root);
  EvtInfo.Initialize();  
  PhoInfo.Initialize();
  LepInfo.Initialize();
  _nIn= _nOut = 0;
}

void TriggerBitCounter::endJob() 
{
  std::cout << "TriggerBitCounter has " << _nIn << " input events and " << _nOut  << " events" << endl;
}

// dump generator-level information

void TriggerBitCounter::dumpGenInfo(const edm::Event& iEvent)
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

void TriggerBitCounter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  _nIn++;
  if(_dumpHEP)
    dumpGenInfo(iEvent);

  int thisEvent_trigger =0;

  // filter trigger path
  edm::Handle<TriggerResults> TrgResultsHandle;
  bool with_TriggerResults = iEvent.getByLabel(InputTag("TriggerResults::HLT"),TrgResultsHandle);
  if (with_TriggerResults) {
  
    TriggerNames TrgNames( *TrgResultsHandle );   

    vector<string> hlNames_ = TrgNames.triggerNames();
    if(_nIn < 3)
      for (size_t i=0; i<TrgNames.size(); ++i) {
	cout<<"HLT bit = "<<i<<"   "<<hlNames_[i]<<endl;
      }
 
    int NTrigger = TrgNames.size();

    int tempIndex = (unsigned int)TrgNames.triggerIndex( "HLT_L1SingleEG5"); 
    if(tempIndex < NTrigger && TrgResultsHandle->accept(tempIndex))
      thisEvent_trigger |= TRIGGER::HLT_L1SingleEG5;
    
    tempIndex = (unsigned int)TrgNames.triggerIndex( "HLT_L1SingleEG8");
    if(tempIndex < NTrigger && TrgResultsHandle->accept(tempIndex))
      thisEvent_trigger |= TRIGGER::HLT_L1SingleEG8;

    tempIndex = (unsigned int)TrgNames.triggerIndex( "HLT_Photon10_L1R");
    if(tempIndex < NTrigger && TrgResultsHandle->accept(tempIndex))
      thisEvent_trigger |= TRIGGER::HLT_Photon10_L1R;

    tempIndex = (unsigned int)TrgNames.triggerIndex( "HLT_Photon10_LooseEcalIso_TrackIso_L1R");
    if(tempIndex < NTrigger && TrgResultsHandle->accept(tempIndex))
      thisEvent_trigger |= TRIGGER::HLT_Photon10_LooseEcalIso_TrackIso_L1R;
    
    tempIndex = (unsigned int)TrgNames.triggerIndex( "HLT_Photon15_L1R");
    if(tempIndex < NTrigger && TrgResultsHandle->accept(tempIndex))
      thisEvent_trigger |= TRIGGER::HLT_Photon15_L1R;

    tempIndex = (unsigned int)TrgNames.triggerIndex( "HLT_Photon15_TrackIso_L1R");
    if(tempIndex < NTrigger && TrgResultsHandle->accept(tempIndex))
      thisEvent_trigger |= TRIGGER::HLT_Photon15_TrackIso_L1R;

    tempIndex = (unsigned int)TrgNames.triggerIndex( "HLT_Photon15_LooseEcalIso_L1R");
    if(tempIndex < NTrigger && TrgResultsHandle->accept(tempIndex))
      thisEvent_trigger |= TRIGGER::HLT_Photon15_LooseEcalIso_L1R;
    
    tempIndex = (unsigned int)TrgNames.triggerIndex( "HLT_Photon20_LooseEcalIso_TrackIso_L1R");
    if(tempIndex < NTrigger && TrgResultsHandle->accept(tempIndex))
      thisEvent_trigger |= TRIGGER::HLT_Photon20_LooseEcalIso_TrackIso_L1R;

    tempIndex = (unsigned int)TrgNames.triggerIndex( "HLT_Photon25_L1R");
    if(tempIndex < NTrigger && TrgResultsHandle->accept(tempIndex))
      thisEvent_trigger |= TRIGGER::HLT_Photon25_L1R;

    tempIndex = (unsigned int)TrgNames.triggerIndex( "HLT_Mu5");
    if(tempIndex < NTrigger && TrgResultsHandle->accept(tempIndex))
      thisEvent_trigger |= TRIGGER::HLT_Mu5;

  }

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
  EvtInfo.HLT    = thisEvent_trigger;
  
 
  
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
  


    // initialize the variables of ntuples
    PhoInfo.Initialize();

    bool _isFilled = false;
    bool _isFilled_Mu = false;

    
    for( std::vector<pat::Photon>::const_iterator it_ph = PhotonHandle->begin(); 
	 it_ph != PhotonHandle->end(); it_ph++ ) {
       
      if (PhoInfo.Size>=MAX_PHOTONS) {
	fprintf(stderr,"ERROR: number of photons exceeds the size of array.\n");
	exit(1);
      }

      float ecalIso = it_ph->ecalRecHitSumEtConeDR04();
      float hcalIso = it_ph->hcalTowerSumEtConeDR04();
      float trkIso  = it_ph->trkSumPtHollowConeDR04();

      PhoInfo.Index[PhoInfo.Size]        = PhoInfo.Size;  
      PhoInfo.E    [PhoInfo.Size] 	 = it_ph->energy(); 
      PhoInfo.Et   [PhoInfo.Size] 	 = it_ph->et();     
      PhoInfo.Pz   [PhoInfo.Size] 	 = it_ph->pz();     
      PhoInfo.Eta  [PhoInfo.Size] 	 = it_ph->eta();    
      PhoInfo.Phi  [PhoInfo.Size] 	 = it_ph->phi();    
      PhoInfo.R9   [PhoInfo.Size] 	 = it_ph->r9();     
      PhoInfo.TrkPtSum[PhoInfo.Size]     = trkIso;
      PhoInfo.EcalRecHitEtSum[PhoInfo.Size]= ecalIso;
      PhoInfo.HcalTowerEtSum[PhoInfo.Size] = hcalIso;
      PhoInfo.HoverE[PhoInfo.Size]       = it_ph->hadronicOverEm();

      PhoInfo.SCE[PhoInfo.Size]          = it_ph->superCluster()->energy();
      PhoInfo.SCEta[PhoInfo.Size]        = it_ph->superCluster()->eta();
      PhoInfo.SCPhi[PhoInfo.Size]        = it_ph->superCluster()->phi();
      PhoInfo.SCEtaWidth[PhoInfo.Size]   = it_ph->superCluster()->etaWidth();
      PhoInfo.SCPhiWidth[PhoInfo.Size]   = it_ph->superCluster()->phiWidth();
	
      float scet = it_ph->superCluster()->energy()/cosh(it_ph->superCluster()->eta());
      PhoInfo.SCEt[PhoInfo.Size]       = scet;

      PhoInfo.SCNCrystal[PhoInfo.Size] = it_ph->superCluster()->size();

      if(it_ph->genPhoton()){
	  
	PhoInfo.GenPID[PhoInfo.Size]       = 22;
	PhoInfo.GenPt[PhoInfo.Size]        = it_ph->genPhoton()->pt();
	if(it_ph->genPhoton()->mother()){
	  PhoInfo.GenMomPID[PhoInfo.Size]  = it_ph->genPhoton()->mother()->pdgId();
	  PhoInfo.GenMomPt[PhoInfo.Size]   = it_ph->genPhoton()->mother()->pt();
	  if(it_ph->genPhoton()->mother()->mother())
	    PhoInfo.GenGMomPID[PhoInfo.Size] = it_ph->genPhoton()->mother()->mother()->pdgId();
	}
     
      } // check generator-level matches
      else
	
	PhoInfo.GenPID[PhoInfo.Size]       = -1;
      
      reco::Photon rPho = reco::Photon(*it_ph);
      bool passOffline = IsLoosePhoton(rPho);

      PhoInfo.IsLoose[PhoInfo.Size] = passOffline? 1: 0;
   
      PhoInfo.Size++;

      // need to pass the following cuts to be loose photons
      float et = it_ph->et();
      float eta = it_ph->eta();
      float genpt = it_ph->genPhoton()? it_ph->genPhoton()->pt():-999;

      h_debug1->Fill(et);
      h_eta1  ->Fill(eta);
 
      // loose photon EB

      if(!passOffline)continue;
      
      h_debug2->Fill(et);
      h_eta2  ->Fill(eta);

      int genMomPID = it_ph->genPhoton() && it_ph->genPhoton()->mother()? 
	it_ph->genPhoton()->mother()->pdgId(): 0;

      bool isFromHardScattering = (genMomPID ==22);

      bool isFromJet   = (abs(genMomPID)> 50);
      
      bool isFromQuark = (abs(genMomPID)< 10) && genMomPID!=0;
      
      bool isFromGluon = (abs(genMomPID)==21);

      // Use SingleEG5 as the denominator

      if((thisEvent_trigger & TRIGGER::HLT_L1SingleEG5) && !_isFilled)
	{
	  _isFilled=true;
	  h_recgetdeno->Fill(et);
	  if(isFromHardScattering)
	    {
	      h_getdeno->Fill(et);
	      h_gengetdeno->Fill(genpt);
	    }
	  else if(isFromJet)
	    h_jetgetdeno->Fill(et);
	  else if(isFromQuark)
	    h_qgetdeno->Fill(et);
	  else if(isFromGluon)
	    h_ggetdeno->Fill(et);
	  
	  if(thisEvent_trigger & TRIGGER::HLT_Photon15_L1R)
	    {
	      h_recgetnumr->Fill(et);
	      if(isFromHardScattering)
		{
		  h_getnumr->Fill(et);
		  h_gengetnumr->Fill(genpt);
		}
	      else if(isFromJet)
		h_jetgetnumr->Fill(et);
	      else if(isFromQuark)
		h_qgetnumr->Fill(et);
	      else if(isFromGluon)
		h_ggetnumr->Fill(et);
	    } // if pass 15 trigger
	} // if pass L1EG5 trigger


      // Use Muon trigger as the denominator

      if((thisEvent_trigger & TRIGGER::HLT_Mu5) && !_isFilled_Mu)
	{
	  _isFilled_Mu=true;
  	  h_mugetdeno->Fill(et);

 	  if(thisEvent_trigger & TRIGGER::HLT_Photon15_L1R)	    
	    {h_mugetnumr->Fill(et);}

	} // if passing muon trigger

    } // if number of photons > 0


  } // if there are photon info
  
  
  // look for muons, electrons
  edm::Handle<std::vector<pat::Muon> >     MuonHandle;
  edm::Handle<std::vector<pat::Electron> > ElectronHandle;

  bool hasMuo = iEvent.getByLabel(_muoLabel,     MuonHandle);
  bool hasEle = iEvent.getByLabel(_eleLabel, ElectronHandle);
  
  
  // initialize the variables of ntuples
  if(hasMuo || hasEle){
    LepInfo.Initialize();

    if(hasMuo){
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
    }
  
    if(hasEle){
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

//       h_leppt -> Fill(it_el->pt());
//       h_lepeta-> Fill(it_el->eta());

    }
    }
  } // if there are muons or electrons
  
  
  root->Fill();
  _nOut++;
  
}


bool TriggerBitCounter::IsLoosePhoton(reco::Photon it_ph)
{
  float et = it_ph.et();
  float eta = it_ph.eta();

  float ecalIso = it_ph.ecalRecHitSumEtConeDR04();
  float hcalIso = it_ph.hcalTowerSumEtConeDR04();
  float trkIso  = it_ph.trkSumPtHollowConeDR04();

  bool isBarrel = it_ph.isEB();
  bool isEndCap = it_ph.isEE();
  bool inAnyGap = it_ph.isEBEEGap() || 
    (it_ph.isEB() && it_ph.isEBGap()) || 
    (it_ph.isEE() && it_ph.isEEGap());

  if(inAnyGap)return false;

  // Barrel cuts
  if(isBarrel && ecalIso  > 5.0 + 0.0045*et)return false;
  if(isBarrel && hcalIso  > 5.0)return false;
  if(isBarrel && trkIso > 9.0 )return false;
  if(isBarrel && it_ph.hadronicOverEm() > 0.15)return false;

  // Endcap cuts
  if(isEndCap && ecalIso  > 5.0 + 0.02*et)return false;
  if(isEndCap && hcalIso  > 7.0)return false;
  if(isEndCap && trkIso > 9.0 )return false;
  if(isEndCap && it_ph.hadronicOverEm() > 0.15)return false;

  if(et        < 10.0)return false;
  if(fabs(eta) > 1.44 && fabs(eta)<1.56)return false;
  if(fabs(eta) > 2.5)return false;

  return true;

}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerBitCounter);


