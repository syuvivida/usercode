// -----------------------------------------------------
// EikoTrigMatch.cc 
// -- a very simple example analyzer for RECO data
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
#include "CommonTools/Utils/interface/PtComparator.h"


#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"


#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

// this file contains the format of lepton, photon, and event structures
#include "format.hh" 
                     

using namespace edm;
using namespace reco;
using namespace std;
using namespace math;
using namespace ROOT;


class Histo_struct {

public:

  TH1F* h_dRL15;
  TH1F* h_dRL18;
  TH1F* h_dRHLT15;

  TH1F* h_getdeno;
  TH1F* h_getnumr;
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


  void BookHistograms(edm::Service<TFileService> fo)
  {
    
    int nbin=100;
    double xmin=0.0;
    double xmax=200.0;

    h_dRL15  = fo->make<TH1F>("h_dRL15","#Delta R between photons and"
			      "trigger L1EG5", 100,0,10);
    h_dRL18  = fo->make<TH1F>("h_dRL18","#Delta R between photons and"
			      "trigger L1EG8", 100,0,10);
    h_dRHLT15= fo->make<TH1F>("h_dRHLT15","#Delta R between photons and"
			      "trigger HLT15", 100,0,10);

    h_eta1  = fo->make<TH1F>("h_eta1","#eta of photons "
				"before loose photon cuts", 60,-3.0,3.0);
    h_eta2  = fo->make<TH1F>("h_eta2","#eta of photons "
				"after loose photon cuts", 60,-3.0,3.0);

    h_debug1  = fo->make<TH1F>("h_debug1","Reconstructed photon "
			       "Et before loose photon cuts", nbin,xmin,xmax);
    h_debug2 = fo->make<TH1F>("h_debug2","Reconstructed photon "
			       "Et before loose photon cuts", nbin,xmin,xmax);
    h_getdeno = fo->make<TH1F>("h_getdeno","Reconstructed and matched photon "
			       "Et before photon trigger cuts", nbin,xmin,xmax);
    h_getnumr = fo->make<TH1F>("h_getnumr","Reconstructed and matched photon "
			       "Et after photon trigger cuts", nbin,xmin,xmax);

    h_jetgetdeno = fo->make<TH1F>("h_jetgetdeno","Reconstructed and jet photon"
				  " Et before photon trigger cuts", 
				  nbin,xmin,xmax);
    h_jetgetnumr = fo->make<TH1F>("h_jetgetnumr","Reconstructed and jet photon"
				  " Et after photon trigger cuts", 
				  nbin,xmin,xmax);
    h_qgetdeno   = fo->make<TH1F>("h_qgetdeno","Reconstructed and quark photon"
				  " Et before photon trigger cuts", 
				  nbin,xmin,xmax);
    h_qgetnumr   = fo->make<TH1F>("h_qgetnumr","Reconstructed and quark photon"
				  " Et after photon trigger cuts", 
				  nbin,xmin,xmax);
    h_ggetdeno   = fo->make<TH1F>("h_ggetdeno","Reconstructed and gluon photon"
				  " Et before photon trigger cuts", 
				  nbin,xmin,xmax);
    h_ggetnumr   = fo->make<TH1F>("h_ggetnumr","Reconstructed and gluon photon"
				  " Et after photon trigger cuts", 
				  nbin,xmin,xmax);

  }

};


class EikoTrigMatch : public edm::EDAnalyzer {
public:
  explicit EikoTrigMatch(const edm::ParameterSet&) ;
  ~EikoTrigMatch();  

    
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void dumpGenInfo(const edm::Event&); 
  virtual void analyze(const edm::Event&, const edm::EventSetup&) ;
  virtual void endJob() ;
  
  TTree *root;
  EvtInfoBranches  EvtInfo;
  PhoInfoBranches  PhoInfo;
  Histo_struct     HistoInfo;

  bool             _dumpHEP;
  int              _nIn;
  int              _nOut;

};


EikoTrigMatch::EikoTrigMatch(const edm::ParameterSet& iConfig):
  _nIn(0), _nOut(0)

{  
  _dumpHEP = iConfig.getUntrackedParameter<bool>("dumpHEP", true);
 
}


EikoTrigMatch::~EikoTrigMatch()
{
}

void EikoTrigMatch::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  TFileDirectory results = TFileDirectory( fs->mkdir("EikoTrigMatch") );
  HistoInfo.BookHistograms(fs);
  root = new TTree("root","root");
  EvtInfo.Register(root);  
  PhoInfo.Register(root);
  EvtInfo.Initialize();  
  PhoInfo.Initialize();
  _nIn= _nOut = 0;
}

void EikoTrigMatch::endJob() 
{
  std::cout << "EikoTrigMatch has " << _nIn << " input events and " << _nOut  << " output events" << endl;
}

// dump generator-level information

void EikoTrigMatch::dumpGenInfo(const edm::Event& iEvent)
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

void EikoTrigMatch::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  // fill L1 and HLT info
  // get objects possed by each filter
  edm::Handle<trigger::TriggerEventWithRefs> triggerObj;
  iEvent.getByLabel("hltTriggerSummaryRAW",triggerObj); 
  if(!triggerObj.isValid()) { 
    cout << "Can't find hltTriggerSummaryRAW" << endl;
    return;
  }

  std::vector< edm::Ref<l1extra::L1EmParticleCollection> > L1Trig5Cands;
  std::vector< edm::Ref<l1extra::L1EmParticleCollection> > L1Trig8Cands;
  std::vector< edm::Ref<l1extra::L1EmParticleCollection> > tempL1Trig5Cands;
  std::vector< edm::Ref<l1extra::L1EmParticleCollection> > tempL1Trig8Cands;
  std::vector< edm::Ref<reco::RecoEcalCandidateCollection> > HLTTrigCands;

//   const edm::InputTag trig1("hltL1NonIsoHLTLEITISinglePhotonEt25EtFilter"
  const edm::InputTag trigL15("hltL1sRelaxedSingleEgammaEt5", "", "HLT");
  const edm::InputTag trigL18("hltL1sRelaxedSingleEgammaEt8", "", "HLT");
  const edm::InputTag trigHLT15("hltL1NonIsoSinglePhotonEt15LEIEtFilter","", "HLT");
  
  bool hasL15   = false;
  bool hasL18   = false;
  bool hasHLT15 = false;
    
  if ( triggerObj->filterIndex(trigL15)<triggerObj->size())hasL15=true;
  if ( triggerObj->filterIndex(trigL18)<triggerObj->size())hasL18=true;
  if ( triggerObj->filterIndex(trigHLT15)<triggerObj->size())hasHLT15=true;


  int typeL15= trigger::TriggerL1NoIsoEG;
  int temptypeL15= trigger::TriggerL1IsoEG;
  int typeL18= trigger::TriggerL1NoIsoEG;
  int temptypeL18= trigger::TriggerL1IsoEG;
  int typeHLT15 = trigger::TriggerCluster;


  ////////////////////////////////////////////////////////////
  //      Retrieve saved filter objects                     //
  ////////////////////////////////////////////////////////////
  if(hasL15)
    {
      triggerObj->getObjects(triggerObj->filterIndex(trigL15),typeL15,
			     L1Trig5Cands);

      triggerObj->getObjects(triggerObj->filterIndex(trigL15),temptypeL15,
			     tempL1Trig5Cands);      
      
      for(unsigned int itemp=0; itemp < tempL1Trig5Cands.size(); itemp++)
	L1Trig5Cands.push_back(tempL1Trig5Cands[itemp]);
      
    }
  if(hasL18)
    {
      triggerObj->getObjects(triggerObj->filterIndex(trigL18),typeL18,
			     L1Trig8Cands);

      triggerObj->getObjects(triggerObj->filterIndex(trigL18),temptypeL18,
			     tempL1Trig8Cands);      
      
      for(unsigned int itemp=0; itemp < tempL1Trig8Cands.size(); itemp++)
	L1Trig8Cands.push_back(tempL1Trig8Cands[itemp]);
     
    }

  if(hasHLT15)
    triggerObj->getObjects(triggerObj->filterIndex(trigHLT15),typeHLT15,
			   HLTTrigCands);


//   cout << hasL15 << "\t" << hasL18 << "\t" << hasHLT15 << endl;
//   cout << L1Trig5Cands.size() << "\t" << L1Trig8Cands.size() << "\t"
//        << HLTTrigCands.size() << endl;

  // Start to fill the main root branches
  // initialize the variables of ntuples
  EvtInfo.Initialize();  

  
  EvtInfo.RunNo  = iEvent.id().run();
  EvtInfo.EvtNo  = iEvent.id().event();
  EvtInfo.HLT    = thisEvent_trigger;
  
  // look for generator photons
  edm::Handle<reco::GenParticleCollection> GenHandle;  
  bool hasGenParticle = iEvent.getByLabel("genParticles", GenHandle);

  if(!hasGenParticle)return;

  // initialize the variables of ntuples
  PhoInfo.Initialize();
  bool _isFilled = false;


  for( std::vector<GenParticle>::const_iterator it_gen = 
	 GenHandle->begin(); 
       it_gen != GenHandle->end(); it_gen++ ) {

    float et = it_gen->et();
    float eta = it_gen->eta();
    
    if(it_gen->pdgId()!=22 || it_gen->status()!=1)continue;
    if(et < 2)continue;
    if(fabs(eta)>2.5)continue;
    if(fabs(eta)>1.44 && fabs(eta)<1.56)continue;

    if (PhoInfo.Size>=MAX_PHOTONS) {
      fprintf(stderr,"ERROR: number of photons exceeds the size of array.\n");
      exit(1);
    }
    
    PhoInfo.Index[PhoInfo.Size]          = PhoInfo.Size;  
    PhoInfo.E    [PhoInfo.Size] 	 = it_gen->energy(); 
    PhoInfo.Et   [PhoInfo.Size] 	 = it_gen->et();     
    PhoInfo.Pz   [PhoInfo.Size] 	 = it_gen->pz();     
    PhoInfo.Eta  [PhoInfo.Size] 	 = it_gen->eta();    
    PhoInfo.Phi  [PhoInfo.Size] 	 = it_gen->phi();    

    int trigMatchCode=0;

    // start doing L1 trigger matching
    float closestDeltaR = 0.5;
    float closestEcalCandIndex = -1;
    for (unsigned int j=0; j< L1Trig5Cands.size(); j++) {
      float deltaR = reco::deltaR(L1Trig5Cands[j]->momentum(),
				  it_gen->momentum());
//       cout << "deltaR L1EG5 : " << deltaR << endl;
      HistoInfo.h_dRL15->Fill(deltaR);

      if (deltaR < closestDeltaR) {
	closestDeltaR = deltaR;
	closestEcalCandIndex = j;
	break;
      }
    }
    if(closestEcalCandIndex >=0)trigMatchCode |= TRIGGER::HLT_L1SingleEG5;

    closestDeltaR = 0.5;
    closestEcalCandIndex = -1;
    for (unsigned int j=0; j< L1Trig8Cands.size(); j++) {
      float deltaR = reco::deltaR(L1Trig8Cands[j]->momentum(),
				  it_gen->momentum());
      
      HistoInfo.h_dRL18->Fill(deltaR);
//       cout << "deltaR L1EG8 : " << deltaR << endl;
      if (deltaR < closestDeltaR) {
	closestDeltaR = deltaR;
	closestEcalCandIndex = j;
	break;
      }
    }
    if(closestEcalCandIndex >=0)trigMatchCode |= TRIGGER::HLT_L1SingleEG8;

      // start doing HLT trigger matching
    closestDeltaR = 0.5;
    closestEcalCandIndex = -1;
    for (unsigned int j=0; j< HLTTrigCands.size(); j++) {
      float deltaR = reco::deltaR(HLTTrigCands[j]->momentum(),
				  it_gen->momentum());
      
      HistoInfo.h_dRHLT15->Fill(deltaR);
//       cout << "deltaR HLT : " << deltaR << endl;
      if (deltaR < closestDeltaR) {
	closestDeltaR = deltaR;
	closestEcalCandIndex = j;
	break;
      }
    }
    if(closestEcalCandIndex >=0)trigMatchCode |= TRIGGER::HLT_Photon15_L1R;
    
      

    PhoInfo.Trig[PhoInfo.Size] = trigMatchCode;

    // started doing MC-reconstruction matching
    
    float genpt   = -999;
    int genMomPID = 0;
    
    genpt = it_gen->pt();
    PhoInfo.GenPID[PhoInfo.Size]       = it_gen->pdgId();
    PhoInfo.GenPt[PhoInfo.Size]        = genpt;
    if(it_gen->mother()){
      genMomPID = it_gen->mother()->pdgId();
      PhoInfo.GenMomPID[PhoInfo.Size]  = genMomPID;
      PhoInfo.GenMomPt[PhoInfo.Size]   = it_gen->mother()->pt();
      if(it_gen->mother()->mother())
	PhoInfo.GenGMomPID[PhoInfo.Size] = it_gen->mother()->mother()->pdgId();
    }
  
//     cout << "genpt = " << genpt << " mother = " << genMomPID << endl;

    PhoInfo.Size++;

    
    // need to pass the following cuts to be loose photons
    HistoInfo.h_debug1->Fill(et);
    HistoInfo.h_eta1  ->Fill(eta);
    

    bool isFromHardScattering = (genMomPID ==22);
    
    bool isFromJet   = (abs(genMomPID)> 50);
    
    bool isFromQuark = (abs(genMomPID)< 10) && genMomPID!=0;
    
    bool isFromGluon = (abs(genMomPID)==21);

    // Use SingleEG5 as the denominator


    if((trigMatchCode & TRIGGER::HLT_L1SingleEG5) && !_isFilled)
      {
	_isFilled=true;
	if(isFromHardScattering)
	  HistoInfo.h_getdeno->Fill(et);
	else if(isFromJet)
	  HistoInfo.h_jetgetdeno->Fill(et);
	else if(isFromQuark)
	  HistoInfo.h_qgetdeno->Fill(et);
	else if(isFromGluon)
	  HistoInfo.h_ggetdeno->Fill(et);
	  
	if(trigMatchCode & TRIGGER::HLT_Photon15_L1R)
	  {
	    if(isFromHardScattering)	      
	      HistoInfo.h_getnumr->Fill(et);
	    else if(isFromJet)
	      HistoInfo.h_jetgetnumr->Fill(et);
	    else if(isFromQuark)
	      HistoInfo.h_qgetnumr->Fill(et);
	    else if(isFromGluon)
	      HistoInfo.h_ggetnumr->Fill(et);
	  } // if pass 15 trigger
      } // if pass L1EG5 trigger

  } // it_gen
  root->Fill();
  
  _nOut++;
  
}



//define this as a plug-in
DEFINE_FWK_MODULE(EikoTrigMatch);


