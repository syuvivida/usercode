// -----------------------------------------------------
// Zee.cc 
// -- a very simple example analyzer for RECO data
// ----------------------------------------------------- 
// Shin-Shan Yu


#include <TROOT.h>
#include <TSystem.h>
#include <TObject.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
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

// #include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

// this file contains the format of lepton, photon, and event structures
#include "zeeformat.hh" 
#include "trigformat.hh"
#include <map>                     

using namespace edm;
using namespace std;
using namespace reco;
using namespace trigger;
using namespace math;
using namespace ROOT;
using namespace zee;

typedef std::map<reco::GsfElectronCollection::const_iterator,reco::PhotonCollection::const_iterator> elePhoMap;
typedef std::map<reco::GsfElectronCollection::const_iterator,reco::GenParticleCollection::const_iterator> eleGenMap;
typedef std::map< reco::PhotonCollection::const_iterator,std::vector<l1extra::L1EmParticle>::const_iterator > phoL1Map;
typedef std::map< reco::PhotonCollection::const_iterator,trigger::TriggerObject > phoL3Map;



class Zee : public edm::EDAnalyzer {
public:
  explicit Zee(const edm::ParameterSet&) ;
  ~Zee();  

    
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void dumpGenInfo(const edm::Event&); 
  virtual void analyze(const edm::Event&, const edm::EventSetup&) ;
  virtual void endJob() ;
  bool isLoosePhoton(const reco::Photon& it_ph);
  void MatchElectronToPhoton(const edm::Event&);
  void MatchElectronToGenp(const edm::Event&);
  void MatchPhoToL1(const edm::Event&);
  void MatchPhoToL3(const edm::Event&, 
		    std::string trgPath,
		    std::string tag,
		    phoL3Map& mymap);
  void TurnOnHLTBit(std::string trgPath, 
		    int         trgCode);


  edm::InputTag    _eleLabel;
  edm::InputTag    _phoLabel;
  edm::InputTag    _genLabel;

  TTree *root;
  EvtInfoBranches  EvtInfo;
  PhoInfoBranches  PhoInfo;
  GenInfoBranches  GenInfo;
  elePhoMap        myElePhoMap;
  eleGenMap        myEleGenMap;
  phoL1Map         myPhoL1Map;
  phoL3Map         myPhoL3Map_HLTL1EG5;
  phoL3Map         myPhoL3Map_HLTL1EG8;
  phoL3Map         myPhoL3Map_HLT10;
  phoL3Map         myPhoL3Map_HLT15;
  phoL3Map         myPhoL3Map_HLT20Iso;
  phoL3Map         myPhoL3Map_HLT25;
  phoL3Map         myPhoL3Map_HLTMu5;


  int              _thisEvent_trigger;


  bool             _dumpHEP;
  int              _nIn;
  int              _nOut;

  edm::Handle<trigger::TriggerEvent> trgEvent;
  edm::Handle<TriggerResults> TrgResultsHandle;


  // histogram
  TH1F* h_trkisodeno;
  TH1F* h_trkisonumr;
  TH1F* h_sigietadeno;
  TH1F* h_sigietanumr;
  TH2F* h_xy;
  TH2F* h_rz;
  TH1F* h_dE;



};


Zee::Zee(const edm::ParameterSet& iConfig):
  _nIn(0), _nOut(0), _thisEvent_trigger(0),
  _phoLabel(iConfig.getParameter< edm::InputTag >( "phoLabel" ) ),
  _eleLabel(iConfig.getParameter< edm::InputTag >( "eleLabel" ) ),
  _genLabel(iConfig.getParameter< edm::InputTag >( "genLabel" ) )

{  
  _dumpHEP = iConfig.getUntrackedParameter<bool>("dumpHEP", false);
 
}


Zee::~Zee()
{
}

void Zee::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  TFileDirectory results = TFileDirectory( fs->mkdir("Zee") );


  h_trkisodeno = fs->make<TH1F>("h_trkisodeno","Track isolation denominator",
				200,0,100);
  h_trkisonumr = fs->make<TH1F>("h_trkisonumr","Track isolation numerator",
				200,0,100);
  h_sigietadeno = fs->make<TH1F>("h_sigietadeno","Sigma ieta ieta denominator",
				 200,0,1);
  h_sigietanumr = fs->make<TH1F>("h_sigietanumr","Sigma ieta ieta numerator",
				 200,0,1);
  
  h_xy          = fs->make<TH2F>("h_xy","",100,-10.0,10.0,100,-10.0,10.0);
  h_rz          = fs->make<TH2F>("h_rz","",100,0,10.0,100,-10.0,10.0);
  h_dE          = fs->make<TH1F>("h_dE","",100,0,2);


  root = new TTree("root","root");
  EvtInfo.Register(root);  
  PhoInfo.Register(root);
  GenInfo.Register(root);
  EvtInfo.Initialize();  
  PhoInfo.Initialize();
  GenInfo.Initialize();
}

void Zee::endJob() 
{
  std::cout << "Zee has " << _nIn << " input events and " << _nOut  << " events" << endl;
}

// dump generator-level information

void Zee::dumpGenInfo(const edm::Event& iEvent)
{
  edm::Handle<reco::GenParticleCollection> GenHandle;  
  bool hasGenParticle = iEvent.getByLabel(_genLabel, GenHandle);
  if(!hasGenParticle)return;

  std::cout << "============================ " << 
    "Run " <<  iEvent.id().run() << " Evt " <<  iEvent.id().event() <<
    " ============================= " << std::endl << std::endl;

  int genIndex = 0;

  printf("GenIndex  ");
  printf("PDG  ");
  printf("Status ");
  printf("Mother PID");
  printf("   ");
  printf("Mass ");
  printf("Energy ");
  printf("Pt ");
  printf("Eta ");
  printf("\n");


  for( GenParticleCollection::const_iterator it_gen = GenHandle->begin(); 
       it_gen != GenHandle->end(); it_gen++ ) {

    printf("%4i", genIndex);
    printf("%7i", it_gen->pdgId());
    printf("%6i", it_gen->status());
    if(it_gen->mother())
      printf("%7i", it_gen->mother()->pdgId());
    else
      printf("%7i", -1);
    printf("   ");
    printf("%9.3f",it_gen->p4().mass());
    printf("%9.3f",it_gen->energy());
    printf("%9.3f",it_gen->pt());
    printf("%9.3f",it_gen->eta());
    printf("\n");
    genIndex ++;
  }

  std::cout << std::endl 
	    << "==========================================================="
            << "===============" 
	    << std::endl;

}


// analyzing reconstructed electrons, muons, and photons

void Zee::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  _nIn++;
  _thisEvent_trigger = 0;
  
  if(_dumpHEP)
    dumpGenInfo(iEvent);

  myElePhoMap.clear();
  myEleGenMap.clear();
  myPhoL1Map.clear();
  myPhoL3Map_HLTL1EG5.clear();
  myPhoL3Map_HLTL1EG8.clear();
  myPhoL3Map_HLT10.clear();
  myPhoL3Map_HLT15.clear();
  myPhoL3Map_HLT20Iso.clear();
  myPhoL3Map_HLT25.clear();
  myPhoL3Map_HLTMu5.clear();

  MatchElectronToPhoton(iEvent);
  MatchElectronToGenp(iEvent);

 // photon - to L1 matching
  MatchPhoToL1(iEvent);


  // HLT_L1EG5
  MatchPhoToL3(iEvent,"hltL1sRelaxedSingleEgammaEt5",
	       "HLT", myPhoL3Map_HLTL1EG5);

  // HLT_L1EG8
  MatchPhoToL3(iEvent,"hltL1sRelaxedSingleEgammaEt8",
	       "HLT",myPhoL3Map_HLTL1EG8);


  // HLT_10L1R
  MatchPhoToL3(iEvent,"hltL1NonIsoHLTNonIsoSinglePhotonEt10HcalIsolFilter",
	       "HLT",myPhoL3Map_HLT10);
   
  // HLT_15L1R
  MatchPhoToL3(iEvent,"hltL1NonIsoHLTNonIsoSinglePhotonEt15HcalIsolFilter",
	       "HLT",myPhoL3Map_HLT15);


  // HLT_20IsoL1R
  MatchPhoToL3(iEvent,"hltL1NonIsoHLTLEITISinglePhotonEt20TrackIsolFilter",
	       "HLT",myPhoL3Map_HLT20Iso);

  // HLT_25L1R
  MatchPhoToL3(iEvent,"hltL1NonIsoHLTNonIsoSinglePhotonEt25HcalIsolFilter",
	       "HLT",myPhoL3Map_HLT25);

  // HLT_MU5
  MatchPhoToL3(iEvent,"hltSingleMu5L3Filtered5",
	       "HLT",myPhoL3Map_HLTMu5);


  if(_nIn < 3)
    {
      cout << "myPhoL3Map_HLTL1EG5_L3Obj.size() ="  << myPhoL3Map_HLTL1EG5.size() << endl;
      cout << "myPhoL3Map_HLTL1EG8_L3Obj.size() ="  << myPhoL3Map_HLTL1EG8.size() << endl;
      cout << "myPhoL3Map_HLT10_L3Obj.size() = " << myPhoL3Map_HLT10.size() << endl;
      cout << "myPhoL3Map_HLT15_L3Obj.size() =" << myPhoL3Map_HLT15.size() << endl;
      cout << "myPhoL3Map_HLTMu5_L3Obj.size() =" << myPhoL3Map_HLTMu5.size() << endl;
      cout << "myPhoL3Map_HLT20Iso_L3Obj.size() =" << myPhoL3Map_HLT20Iso.size() << endl;
      cout << "myPhoL3Map_HLT25_L3Obj.size() = " << myPhoL3Map_HLT25.size() << endl;
    }


  
  // Start to fill the main root branches
  // initialize the variables of ntuples
  // filter trigger path
  //1E31
  bool with_TriggerResults = iEvent.getByLabel(InputTag("TriggerResults::HLT"),TrgResultsHandle);
  // 8E29
//   bool with_TriggerResults = iEvent.getByLabel(InputTag("TriggerResults::HLT8E29"),TrgResultsHandle);
  
  if (with_TriggerResults) {
  
    TriggerNames TrgNames( *TrgResultsHandle );   

    vector<string> hlNames_ = TrgNames.triggerNames();
    if(_nIn < 3)
      for (size_t i=0; i<TrgNames.size(); ++i) {
	cout<<"HLT bit = "<<i<<"   "<<hlNames_[i]<<endl;
      }

    TurnOnHLTBit("HLT_L1SingleEG5",TRIGGER::HLT_L1SingleEG5);

    TurnOnHLTBit("HLT_L1SingleEG8",TRIGGER::HLT_L1SingleEG8);

    TurnOnHLTBit("HLT_Photon10_L1R",TRIGGER::HLT_Photon10_L1R);

    TurnOnHLTBit("HLT_Photon10_LooseEcalIso_TrackIso_L1R",
		 TRIGGER::HLT_Photon10_LooseEcalIso_TrackIso_L1R);
    
    TurnOnHLTBit("HLT_Photon15_L1R",TRIGGER::HLT_Photon15_L1R);
    
    TurnOnHLTBit("HLT_Photon15_TrackIso_L1R",TRIGGER::HLT_Photon15_TrackIso_L1R);

    TurnOnHLTBit("HLT_Photon15_LooseEcalIso_L1R",TRIGGER::HLT_Photon15_LooseEcalIso_L1R);

    TurnOnHLTBit("HLT_Photon20_LooseEcalIso_TrackIso_L1R",TRIGGER::HLT_Photon20_LooseEcalIso_TrackIso_L1R);

    TurnOnHLTBit("HLT_Photon25_L1R",TRIGGER::HLT_Photon25_L1R);

    TurnOnHLTBit("HLT_Mu5",TRIGGER::HLT_Mu5);
//     TurnOnHLTBit("HLT_L1MuOpen",TRIGGER::HLT_Mu5);

  } // there's trigger results

 

  EvtInfo.Initialize();  
  EvtInfo.RunNo  = iEvent.id().run();
  EvtInfo.EvtNo  = iEvent.id().event();
  EvtInfo.HLT    = _thisEvent_trigger;


  // fill generator information
  edm::Handle<reco::GenParticleCollection> GenHandle;  
  bool hasGenParticle = iEvent.getByLabel(_genLabel, GenHandle); 
  GenInfo.Initialize();
  if(hasGenParticle){
    for( std::vector<GenParticle>::const_iterator it_gen = 
	   GenHandle->begin(); 
	 it_gen != GenHandle->end() && GenInfo.Size < MAX_GENS; it_gen++ ) {
   
      if(it_gen->pdgId()==23 && GenInfo.Size < 10)EvtInfo.GenZMass = it_gen->mass();
      GenInfo.PID[GenInfo.Size] = it_gen->pdgId();
      GenInfo.MPID[GenInfo.Size] = it_gen->mother()?
	it_gen->mother()->pdgId():0;
      GenInfo.Mass[GenInfo.Size] = it_gen->mass();
      GenInfo.Pt[GenInfo.Size] = it_gen->pt();
      GenInfo.Eta[GenInfo.Size] = it_gen->eta();
      GenInfo.Phi[GenInfo.Size] = it_gen->phi();
      GenInfo.Size ++;

    } // check generator-level
  } // end of if hasGenParticle



//   const edm::InputTag ebColName = edm::InputTag("reducedEcalRecHitsEB");
//   const edm::InputTag esColName = edm::InputTag("reducedEcalRecHitsEE");

//   Handle<EcalRecHitCollection> EBReducedRecHits;
//   iEvent.getByLabel(ebColName, EBReducedRecHits);
//   Handle<EcalRecHitCollection> EEReducedRecHits;
//   iEvent.getByLabel(esColName, EEReducedRecHits);
//   EcalClusterLazyTools lazyTool(iEvent, iSetup, ebColName, esColName);

  // now loop over Z electron pairs
  PhoInfo.Initialize();
  int it_count=-1;
  for(elePhoMap::iterator it_e=myElePhoMap.begin(); 
      it_e!= myElePhoMap.end(); ++it_e)
    {
      it_count ++;

      int jt_count=-1;
      for(elePhoMap::const_iterator jt_e = myElePhoMap.begin();
	  jt_e != myElePhoMap.end(); ++jt_e){

	jt_count ++;
	if(it_count >= jt_count)continue; // each pair is used once only
	if(jt_e->first == it_e->first)continue;


	EvtInfo.RecZMass = (it_e->first->p4() 
			    + jt_e->first->p4()).M();

	reco::GsfElectronCollection::const_iterator eleMatch[2];
	eleMatch[0] = it_e->first;
	eleMatch[1] = jt_e->first;

	// initialize the variables of ntuples

	for (unsigned int ey=0; ey < 2; ey++){

	  // first fill electron variables
	  reco::GsfElectronCollection::const_iterator thisEle = 
	    eleMatch[ey];

	  // 	const reco::CaloClusterPtr eleSeed = (*thisEle).superCluster()->seed();
	  // 	float eR9 = thisEle->superCluster()->rawEnergy() > 1e-6? 
	  // 	  lazyTool.e3x3(*eleSeed)/thisEle->superCluster()->rawEnergy():-999;
	  float escet = thisEle->superCluster()->energy()/cosh(thisEle->superCluster()->eta());
	  PhoInfo.eE[PhoInfo.ZPairSize][ey]   = thisEle->energy();
	  PhoInfo.eEt[PhoInfo.ZPairSize][ey]  = thisEle->pt();
	  PhoInfo.ePz[PhoInfo.ZPairSize][ey]  = thisEle->pz();
	  PhoInfo.eEta[PhoInfo.ZPairSize][ey] = thisEle->eta();
	  PhoInfo.ePhi[PhoInfo.ZPairSize][ey] = thisEle->phi();
	  // 	PhoInfo.eR9[PhoInfo.ZPairSize][ey]  = eR9;
	  PhoInfo.eTrkPtSum[PhoInfo.ZPairSize][ey] = 
	    thisEle-> dr04TkSumPt();
	  PhoInfo.eEcalRecHitEtSum[PhoInfo.ZPairSize][ey] = 
	    thisEle->dr04EcalRecHitSumEt();
	  PhoInfo.eHcalTowerEtSum[PhoInfo.ZPairSize][ey] = 
	    thisEle->dr04HcalTowerSumEt();
	  PhoInfo.eHoverE[PhoInfo.ZPairSize][ey] = thisEle->hadronicOverEm();
	  PhoInfo.eSCE[PhoInfo.ZPairSize][ey] = thisEle->superCluster()->energy();	

	  PhoInfo.eSCEt[PhoInfo.ZPairSize][ey] = escet;
	  PhoInfo.eSCEta[PhoInfo.ZPairSize][ey] = thisEle->superCluster()->eta();
	  PhoInfo.eSCPhi[PhoInfo.ZPairSize][ey] = thisEle->superCluster()->phi();
	  PhoInfo.eSCEtaWidth[PhoInfo.ZPairSize][ey] = thisEle->superCluster()->etaWidth();
	  PhoInfo.eSCPhiWidth[PhoInfo.ZPairSize][ey] = thisEle->superCluster()->phiWidth();
	  PhoInfo.eSCNCrystal[PhoInfo.ZPairSize][ey] = thisEle->superCluster()->size();
	  PhoInfo.eSigEta[PhoInfo.ZPairSize][ey] = thisEle->sigmaEtaEta(); 
	  PhoInfo.eSigIEta[PhoInfo.ZPairSize][ey] = thisEle->sigmaIetaIeta();
	  PhoInfo.eCharge[PhoInfo.ZPairSize][ey]  = thisEle->charge();
	  PhoInfo.eClass[PhoInfo.ZPairSize][ey]  = thisEle->classification();

	  // second fill photon variables
	  reco::PhotonCollection::const_iterator it_ph = 
	    myElePhoMap[eleMatch[ey]];

	  float ecalIso = it_ph->ecalRecHitSumEtConeDR04();
	  float hcalIso = it_ph->hcalTowerSumEtConeDR04();
	  float trkIso  = it_ph->trkSumPtHollowConeDR04();
	
	  int location = -1; 
	  bool isBarrel = it_ph->isEB();
	  bool isEndCap = it_ph->isEE();
	  bool inAnyGap = it_ph->isEBEEGap() || 
	    (it_ph->isEB() && it_ph->isEBGap()) || 
	    (it_ph->isEE() && it_ph->isEEGap());

	  if(inAnyGap)location=0;
	  else if(isBarrel)location=1;
	  else if(isEndCap)location=2;

	  PhoInfo.Location[PhoInfo.ZPairSize][ey]     = location;
	  PhoInfo.E    [PhoInfo.ZPairSize][ey] 	 = it_ph->energy(); 
	  PhoInfo.Et   [PhoInfo.ZPairSize][ey] 	 = it_ph->et();     
	  PhoInfo.Pz   [PhoInfo.ZPairSize][ey] 	 = it_ph->pz();     
	  PhoInfo.Eta  [PhoInfo.ZPairSize][ey] 	 = it_ph->eta();    
	  PhoInfo.Phi  [PhoInfo.ZPairSize][ey] 	 = it_ph->phi();    
	  PhoInfo.R9   [PhoInfo.ZPairSize][ey] 	 = it_ph->r9();     
	  PhoInfo.TrkPtSum[PhoInfo.ZPairSize][ey]     = trkIso;
	  PhoInfo.EcalRecHitEtSum[PhoInfo.ZPairSize][ey]= ecalIso;
	  PhoInfo.HcalTowerEtSum[PhoInfo.ZPairSize][ey] = hcalIso;
	  PhoInfo.HoverE[PhoInfo.ZPairSize][ey]       = it_ph->hadronicOverEm();

	  PhoInfo.SCE[PhoInfo.ZPairSize][ey]          = it_ph->superCluster()->energy();
	  PhoInfo.SCEta[PhoInfo.ZPairSize][ey]        = it_ph->superCluster()->eta();
	  PhoInfo.SCPhi[PhoInfo.ZPairSize][ey]        = it_ph->superCluster()->phi();
	  PhoInfo.SCEtaWidth[PhoInfo.ZPairSize][ey]   = it_ph->superCluster()->etaWidth();
	  PhoInfo.SCPhiWidth[PhoInfo.ZPairSize][ey]   = it_ph->superCluster()->phiWidth();
	
	  float scet = it_ph->superCluster()->energy()/cosh(it_ph->superCluster()->eta());
	  PhoInfo.SCEt[PhoInfo.ZPairSize][ey]       = scet;

	  PhoInfo.SCNCrystal[PhoInfo.ZPairSize][ey] = it_ph->superCluster()->size();

	  PhoInfo.SigEta[PhoInfo.ZPairSize][ey]       = it_ph->sigmaEtaEta();
	  PhoInfo.SigIEta[PhoInfo.ZPairSize][ey]      = it_ph->sigmaIetaIeta();
      
	  bool passOffline = isLoosePhoton(*it_ph);
	  PhoInfo.IsLoose[PhoInfo.ZPairSize][ey] = passOffline? 1: 0;
 
	  // match to generator level

	  if(myEleGenMap.find(eleMatch[ey])!=myEleGenMap.end())
	    {
	      reco::GenParticleCollection::const_iterator thisGen = myEleGenMap[eleMatch[ey]];
	      PhoInfo.genPt[PhoInfo.ZPairSize][ey] = thisGen->pt();
	      PhoInfo.genEta[PhoInfo.ZPairSize][ey] = thisGen->eta();
	      PhoInfo.genPhi[PhoInfo.ZPairSize][ey] = thisGen->phi();
	    }


	  // ==================================================
	  // check trigger match
	  // ==================================================

	  // check L1 first
	  if(myPhoL1Map.find(it_ph)!=myPhoL1Map.end())
	    {
	      PhoInfo.L1Pt[PhoInfo.ZPairSize][ey] = myPhoL1Map[it_ph]->pt();
	    }

 	  // now check L3
	  int trigMatchCode=0;
	  float l3pt=-999;

	  if(myPhoL3Map_HLTMu5.find(it_ph)!=myPhoL3Map_HLTMu5.end())
	    {
	      trigMatchCode |= TRIGGER::HLT_Mu5;
	      l3pt = myPhoL3Map_HLTMu5[it_ph].pt();
	    }

	  // L1EG5
	  if( myPhoL3Map_HLTL1EG5.find(it_ph)!= myPhoL3Map_HLTL1EG5.end())
	    {
	      trigMatchCode |= TRIGGER::HLT_L1SingleEG5;
	      l3pt = myPhoL3Map_HLTL1EG5[it_ph].pt();
	    }

	  // L1EG8
	  if( myPhoL3Map_HLTL1EG8.find(it_ph)!= myPhoL3Map_HLTL1EG8.end())
	    {
	      trigMatchCode |= TRIGGER::HLT_L1SingleEG8;
	      l3pt = myPhoL3Map_HLTL1EG8[it_ph].pt();
	    }

	  // HLT_10_L1R
	  if(myPhoL3Map_HLT10.find(it_ph)!=myPhoL3Map_HLT10.end())
	    {
	      trigMatchCode |= TRIGGER::HLT_Photon10_L1R;
	      l3pt =myPhoL3Map_HLT10[it_ph].pt();
	    }


	  // HLT_15_L1R
	  if(myPhoL3Map_HLT15.find(it_ph)!=myPhoL3Map_HLT15.end())
	    {
	      trigMatchCode |= TRIGGER::HLT_Photon15_L1R;
	      l3pt =myPhoL3Map_HLT15[it_ph].pt();
	    }

	  // HLT_20Iso
	  if(myPhoL3Map_HLT20Iso.find(it_ph)!=myPhoL3Map_HLT20Iso.end())
	    {
	      trigMatchCode |= TRIGGER::HLT_Photon20_LooseEcalIso_TrackIso_L1R;
	      l3pt =myPhoL3Map_HLT20Iso[it_ph].pt();
	    }
     
	  // HLT_25_L1R
	  if(myPhoL3Map_HLT25.find(it_ph)!=myPhoL3Map_HLT25.end())
	    {
	      trigMatchCode |= TRIGGER::HLT_Photon25_L1R;
	      l3pt =myPhoL3Map_HLT25[it_ph].pt();
	    }


	  PhoInfo.Trig[PhoInfo.ZPairSize][ey] = trigMatchCode;
	  PhoInfo.L3Pt[PhoInfo.ZPairSize][ey] = l3pt;     


	} // end of ntuple ey loop
	
	PhoInfo.ZPairSize ++;
      
      
      } // end of loop over jt_e
    } // end of loop over it_e
    
  root->Fill();
   
  _nOut++;
  
}


void Zee::MatchElectronToPhoton(const edm::Event& iEvent)
{
  myElePhoMap.clear();

  // look for electron collection
  Handle<reco::GsfElectronCollection> electronColl;
  iEvent.getByLabel(_eleLabel, electronColl);
  bool eleIsValid = electronColl.isValid();
  if(!eleIsValid){cout << "Electron is not valid" << endl; return;}

  // look for photon collection
  Handle<reco::PhotonCollection> photonColl;
  iEvent.getByLabel(_phoLabel, photonColl);
  bool phoIsValid = photonColl.isValid();
  if(!phoIsValid){cout << "Photon is not valid" << endl; return;}
  

  // first loop over electrons
  unsigned int count=0;
  for(reco::GsfElectronCollection::const_iterator it_e = electronColl->begin();
      it_e != electronColl->end(); it_e++){

    if(it_e->pt()>10)count ++ ;
    double ptMax   = 0;
    bool findMatch = false;

    // now loop over photons
    reco::PhotonCollection::const_iterator thisPho;
    for (reco::PhotonCollection::const_iterator it_ph = photonColl->begin(); 
	 it_ph!=photonColl->end(); it_ph++){
  
      // find the highest pt photon that has the same supercluster seed as 
      // electron
      double thisPhoEnergy = it_ph->superCluster()->energy();
      double distx = fabs(it_e->superCluster()->seed()->x() -
			  it_ph->superCluster()->seed()->x());
      double disty = fabs(it_e->superCluster()->seed()->y() -
			  it_ph->superCluster()->seed()->y());
      double distrphi = sqrt(distx*distx + disty*disty);

      double distz = fabs(it_e->superCluster()->seed()->z() -
			  it_ph->superCluster()->seed()->z());

      h_xy->Fill(distx,disty);
      h_rz->Fill(distrphi,distz);
   
      double relE  = it_e->superCluster()->seed()->energy()>1e-3? 
	fabs(it_ph->superCluster()->seed()->energy() - 
	     it_e->superCluster()->seed()->energy())
	/it_e->superCluster()->seed()->energy(): -999;
      
      h_dE->Fill(relE);

      if( distx < 3.0 && disty < 3.0 && distz < 3.0 && relE < 0.2	 
	 && thisPhoEnergy > ptMax)
	{
	  ptMax   = thisPhoEnergy;
	  thisPho = it_ph;
	  findMatch = true;
	} // if matched
      
    } // end of loop over photon

    if(findMatch)
      myElePhoMap.insert(std::pair<reco::GsfElectronCollection::const_iterator,
			 reco::PhotonCollection::const_iterator>(it_e,thisPho));
    
  
  } // end of loop over electron

  if(myElePhoMap.size()!=count)
    cout << "There are " << count  << " electrons" << endl;

  return;

}


void Zee::MatchElectronToGenp(const edm::Event& iEvent)
{

  myEleGenMap.clear();

  // look for electron collection
  Handle<reco::GsfElectronCollection> electronColl;
  iEvent.getByLabel(_eleLabel, electronColl);
  bool eleIsValid = electronColl.isValid();
  if(!eleIsValid){cout << "Electron is not valid" << endl; return;}

  // look for Gen particle collection
  edm::Handle<reco::GenParticleCollection> GenHandle;  
  bool hasGenParticle = iEvent.getByLabel("genParticles", GenHandle);
  if(!hasGenParticle)return;

  unsigned int count=0;

  // loop over electron first
  for(reco::GsfElectronCollection::const_iterator it_e = electronColl->begin();
      it_e != electronColl->end(); it_e++){

    double ptMax   = 0;
    bool findMatch = false;

    count ++;

    // find the highest et gen particle from Z matched to electron 
    GenParticleCollection::const_iterator thisGen;
    for(GenParticleCollection::const_iterator it_gen = GenHandle->begin(); 
	it_gen != GenHandle->end(); it_gen++ ) {

      double thisGenEnergy = it_gen->pt();

      if(abs(it_gen->pdgId())!=11 || it_gen->status()!=1)continue;
      if(it_gen->mother() && abs(it_gen->mother()->pdgId())!=11)continue;
      if(it_gen->mother()->mother() 
	 && abs(it_gen->mother()->mother()->pdgId()!=23))continue;
      
      float dR = reco::deltaR(it_e->momentum(), it_gen->momentum());
      float relPt = fabs(it_e->pt()-it_gen->pt())/it_gen->pt();

      if(dR<0.5 && relPt < 1.0 && thisGenEnergy > ptMax)
	{
	  ptMax = thisGenEnergy;
	  thisGen = it_gen;
	  findMatch = true;
	} // if find a match
    } // end of loop over gen particle

    if(findMatch)
      myEleGenMap.insert(std::pair<reco::GsfElectronCollection::const_iterator,
			 reco::GenParticleCollection::const_iterator>
			 (it_e,thisGen));

  } // end of loop over electron

//   if(myEleGenMap.size()!=count)
//     cout << "There are " << count  << " electrons and " << myEleGenMap.size() << " generated." << endl;
  
  return;

}


void Zee::MatchPhoToL1(const edm::Event& iEvent)
{
 // look for photons
  if(myElePhoMap.size()==0)return;

  // check for L1 extra particles
  edm::Handle<l1extra::L1EmParticleCollection> l1EmNonIso;
  edm::Handle<l1extra::L1EmParticleCollection> l1EmIso;

  const edm::InputTag l1EmNonIsoTag = edm::InputTag("hltL1extraParticles","NonIsolated");
  const edm::InputTag l1EmIsoTag = edm::InputTag("hltL1extraParticles","Isolated");


  bool hasL1NoIso=false;
  bool hasL1Iso = false;

  hasL1NoIso = iEvent.getByLabel(l1EmNonIsoTag,l1EmNonIso);
  hasL1Iso = iEvent.getByLabel(l1EmIsoTag,l1EmIso);

  if(!hasL1NoIso && !hasL1Iso)return;


  for(elePhoMap::iterator it_e=myElePhoMap.begin(); 
      it_e!= myElePhoMap.end(); ++it_e)
    {

      reco::PhotonCollection::const_iterator it_ph = it_e->second;

      float ptMax = 0;
      std::vector<l1extra::L1EmParticle>::const_iterator maxPtIter;
      bool hasMatch = false;


      for( std::vector<l1extra::L1EmParticle>::const_iterator it_l1 = 
	     l1EmIso->begin(); 
	   it_l1 != l1EmIso->end(); it_l1++ ) {
    
	float deltaR = reco::deltaR(it_l1->eta(), it_l1->phi(),
				    it_ph->eta(), it_ph->phi());

	float thisL1Pt = it_l1->pt();

	if (deltaR < 0.5 && thisL1Pt > ptMax ) {
	
	  ptMax = thisL1Pt;
	  maxPtIter = it_l1;
	  hasMatch = true;
	}

      } // end of L1 isolated trigger objects


      for( std::vector<l1extra::L1EmParticle>::const_iterator it_l1 = 
	     l1EmNonIso->begin(); 
	   it_l1 != l1EmNonIso->end(); it_l1++ ) {
    
	float deltaR = reco::deltaR(it_l1->eta(), it_l1->phi(),
				    it_ph->eta(), it_ph->phi());

	float thisL1Pt = it_l1->pt();

	if (deltaR < 0.5 && thisL1Pt > ptMax ) {
	
	  ptMax = thisL1Pt;
	  maxPtIter = it_l1;
	  hasMatch = true;
	}

      } // end of L1 non-isolated trigger objects

    
      if(hasMatch)
	myPhoL1Map.insert(std::pair<reco::PhotonCollection::const_iterator,
			  std::vector<l1extra::L1EmParticle>::const_iterator>(it_ph,maxPtIter));    
    } // end of loop over photons

  return;

}

void Zee::MatchPhoToL3(const edm::Event& iEvent, 
		       std::string trgPath,
		       std::string tag,
		       phoL3Map& thisMap)
{
  thisMap.clear();

 // look for photons
  if(myElePhoMap.size()==0)return;

  // check for HLT objects
  const edm::InputTag hltsummaryTag = edm::InputTag("hltTriggerSummaryAOD","","HLT");
  bool has1E31TrigInfo=
    iEvent.getByLabel(hltsummaryTag, trgEvent);    
  if(!has1E31TrigInfo)return;

  const trigger::TriggerObjectCollection& TOC(trgEvent->getObjects());

  const edm::InputTag myLastFilter = edm::InputTag(trgPath,"",tag);
  if ( trgEvent->filterIndex(myLastFilter) >= trgEvent->sizeFilters())return;
  const trigger::Keys& keys( trgEvent->filterKeys( trgEvent->filterIndex(myLastFilter) ));


  for(elePhoMap::iterator it_e=myElePhoMap.begin(); 
      it_e!= myElePhoMap.end(); ++it_e)
    {

    reco::PhotonCollection::const_iterator it_ph = it_e->second;

    float ptMax = 0;
    trigger::TriggerObject maxPtIter;
    bool hasMatch = false;

    for ( unsigned int hlto = 0; hlto < keys.size(); hlto++ ) {
      size_type hltf = keys[hlto];
      const trigger::TriggerObject L3obj(TOC[hltf]);
      
      float deltaR = reco::deltaR(L3obj.eta(), L3obj.phi(),
				  it_ph->eta(), it_ph->phi());
      float thisL3Pt = L3obj.pt();

      if (deltaR < 0.5 && thisL3Pt > ptMax) {
	maxPtIter = L3obj;
	ptMax = thisL3Pt;
	hasMatch = true;
      }
    } // end of loop over trigger objects

    if(hasMatch)
      thisMap.insert(std::pair<reco::PhotonCollection::const_iterator,
			trigger::TriggerObject>(it_ph,maxPtIter));    
    

    }// end of loop over photons


   return;
}

void Zee::TurnOnHLTBit(std::string trgPath, 
		       int         trgCode)
{
    TriggerNames trgName( *TrgResultsHandle);   
    int NTrigger = trgName.size();
    int tempIndex = (unsigned int)trgName.triggerIndex( trgPath); 
    if(tempIndex < NTrigger && TrgResultsHandle->accept(tempIndex))
      _thisEvent_trigger |= trgCode;
    return;
}




bool Zee::isLoosePhoton(const reco::Photon& it_ph)
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
  if(isBarrel && ecalIso  > 5.0 + 0.004*et)return false;
  if(isBarrel && hcalIso  > 5.0)return false;
  if(isBarrel && trkIso > 9.0 )return false;
  if(isBarrel && it_ph.hadronicOverEm() > 0.15)return false;

  // Endcap cuts
  if(isEndCap && ecalIso  > 5.0 + 0.0021*et)return false;
  if(isEndCap && hcalIso  > 5.0)return false;
  if(isEndCap && trkIso > 9.0 )return false;
  if(isEndCap && it_ph.hadronicOverEm() > 0.15)return false;

  if(et        < 10.0)return false;
  if(fabs(eta) > 1.44 && fabs(eta)<1.56)return false;
  if(fabs(eta) > 2.5)return false;

  return true;

}

//define this as a plug-in
DEFINE_FWK_MODULE(Zee);


