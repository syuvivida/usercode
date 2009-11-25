// -----------------------------------------------------
// GenTrig.cc 
// -- a very simple example analyzer for Gen data
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


#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

// this file contains the format of lepton, photon, and event structures
#include "format.hh" 
#include "trigformat.hh"
                     
using namespace edm;
using namespace reco;
using namespace std;
using namespace trigger;
using namespace math;
using namespace ROOT;
using namespace syu;



typedef std::map< reco::GenParticleCollection::const_iterator,std::vector<l1extra::L1EmParticle>::const_iterator > phoL1Map;
typedef std::map< reco::GenParticleCollection::const_iterator,trigger::TriggerObject > phoL3Map;
typedef std::map< float,reco::GenParticleCollection::const_iterator, std::greater<float> > phoEtMap;


class GenTrig : public edm::EDAnalyzer {
public:
  explicit GenTrig(const edm::ParameterSet&) ;
  ~GenTrig();  

    
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void dumpGenInfo(const edm::Event&); 
  virtual void analyze(const edm::Event&, const edm::EventSetup&) ;
  virtual void endJob() ;
  
  TTree *root;
  EvtInfoBranches  EvtInfo;
  GenInfoBranches  GenInfo;


  edm::Handle<trigger::TriggerEvent> trgEvent;
  edm::Handle<TriggerResults> TrgResultsHandle;


  phoEtMap  myPhoEtMap;
  phoL1Map  myPhoL1Map;
  phoL3Map  myPhoL3Map_HLTL1EG5;
  phoL3Map  myPhoL3Map_HLTL1EG8;
  phoL3Map  myPhoL3Map_HLT10;
  phoL3Map  myPhoL3Map_HLT15;
  phoL3Map  myPhoL3Map_HLT20Iso;
  phoL3Map  myPhoL3Map_HLT25;
  phoL3Map  myPhoL3Map_HLTMu5;


  void MatchGenPhoToL1(const edm::Event&);

  void MatchGenPhoToL3(const edm::Event&, 
		       std::string trgPath,
		       std::string tag,
		       phoL3Map& mymap);

  void TurnOnHLTBit(std::string trgPath, 
		    int         trgCode);


  int              _thisEvent_trigger;
  bool             _dumpHEP;
  int              _pdgCode;
  int              _nIn;
  int              _nOut;

  edm::InputTag    _genLabel;


  // histograms
  TH1F* h_dRL15;
  TH1F* h_dRL18;
  TH1F* h_dRHLT15;

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
  TH1F* h_debug3;
  TH1F* h_debug4;
  TH1F* h_eta1;
  TH1F* h_eta2;
  TH1F* h_mugetdeno;
  TH1F* h_mugetnumr;


};


GenTrig::GenTrig(const edm::ParameterSet& iConfig):
  _nIn(0), _nOut(0), _thisEvent_trigger(0),
  _genLabel(iConfig.getParameter< edm::InputTag >( "genLabel" ) )


{  
  _dumpHEP = iConfig.getUntrackedParameter<bool>("dumpHEP", false);
  _pdgCode = iConfig.getUntrackedParameter<int>("pdgCode",  22);
}


GenTrig::~GenTrig()
{
}

void GenTrig::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  TFileDirectory results = TFileDirectory( fs->mkdir("GenTrig") );


  h_dRL15  = fs->make<TH1F>("h_dRL15","#Delta R between photons and"
			    "trigger L1EG5", 100,0,10);
  h_dRL18  = fs->make<TH1F>("h_dRL18","#Delta R between photons and"
			    "trigger L1EG8", 100,0,10);
  h_dRHLT15= fs->make<TH1F>("h_dRHLT15","#Delta R between photons and"
			    "trigger HLT15", 100,0,10);

    
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
  h_debug3 = fs->make<TH1F>("h_debug3","Reconstructed photon "
			    "Et before loose photon cuts", nbin,-10,10);
  h_debug4 = fs->make<TH1F>("h_debug4","Reconstructed photon "
			    "Et before loose photon cuts", nbin,-0.5,99.5);
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
  GenInfo.Register(root);
  EvtInfo.Initialize();  
  GenInfo.Initialize();
  //   _nIn= _nOut = _thisEvent_trigger = 0;

  std::cout << "Focusing on PDG code = " << _pdgCode << std::endl;
}

void GenTrig::endJob() 
{
  std::cout << "GenTrig has " << _nIn << " input events and " << _nOut  << " events" << endl;
}

// dump generator-level information

void GenTrig::dumpGenInfo(const edm::Event& iEvent)
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
  printf("Pt");
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
    printf("%9.3f",it_gen->vx());
    printf("%9.3f",it_gen->vy());
    printf("%9.3f",it_gen->vz());
    printf("\n");
    genIndex ++;
  }

  std::cout << std::endl 
	    << "==========================================================="
            << "===============" 
	    << std::endl;

}


// analyzing reconstructed electrons, muons, and photons

void GenTrig::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  _nIn++;
  _thisEvent_trigger = 0;
  
  if(_dumpHEP)
    dumpGenInfo(iEvent);

  myPhoEtMap.clear();
  myPhoL1Map.clear();
  myPhoL3Map_HLTL1EG5.clear();
  myPhoL3Map_HLTL1EG8.clear();
  myPhoL3Map_HLT10.clear();
  myPhoL3Map_HLT15.clear();
  myPhoL3Map_HLT20Iso.clear();
  myPhoL3Map_HLT25.clear();
  myPhoL3Map_HLTMu5.clear();

  // look for Gen particle collection
  edm::Handle<reco::GenParticleCollection> GenHandle;  
  bool hasGenParticle = iEvent.getByLabel(_genLabel, GenHandle);
  if(!hasGenParticle)return;

  // photon - to L1 matching
  MatchGenPhoToL1(iEvent);

  // HLT_L1EG5
  MatchGenPhoToL3(iEvent,"hltL1sRelaxedSingleEgammaEt5",
		  "HLT", myPhoL3Map_HLTL1EG5);

  // HLT_L1EG8
  MatchGenPhoToL3(iEvent,"hltL1sRelaxedSingleEgammaEt8",
		  "HLT",myPhoL3Map_HLTL1EG8);


  // HLT_10L1R
  MatchGenPhoToL3(iEvent,"hltL1NonIsoHLTNonIsoSinglePhotonEt10HcalIsolFilter",
		  "HLT",myPhoL3Map_HLT10);
   
  // HLT_15L1R
  MatchGenPhoToL3(iEvent,"hltL1NonIsoHLTNonIsoSinglePhotonEt15HcalIsolFilter",
		  "HLT",myPhoL3Map_HLT15);


  // HLT_20IsoL1R
  MatchGenPhoToL3(iEvent,"hltL1NonIsoHLTLEITISinglePhotonEt20TrackIsolFilter",
		  "HLT",myPhoL3Map_HLT20Iso);

  // HLT_25L1R
  MatchGenPhoToL3(iEvent,"hltL1NonIsoHLTNonIsoSinglePhotonEt25HcalIsolFilter",
		  "HLT",myPhoL3Map_HLT25);

  // HLT_MU5
  MatchGenPhoToL3(iEvent,"hltSingleMu5L3Filtered5",
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

 

  // Start to fill the main root branches
  // initialize the variables of ntuples
  EvtInfo.Initialize();  

  
  EvtInfo.RunNo  = iEvent.id().run();
  EvtInfo.EvtNo  = iEvent.id().event();
  EvtInfo.HLT    = _thisEvent_trigger;


  // sort photons by Et

  for (reco::GenParticleCollection::const_iterator it_gen = GenHandle->begin(); 
       it_gen!=GenHandle->end(); it_gen++){
    if(abs(it_gen->pdgId())!=_pdgCode || it_gen->status()!=1)continue;
    if(it_gen->pt() < 2.)continue;
    float et = it_gen->pt();
    myPhoEtMap.insert(std::pair<float,reco::GenParticleCollection::const_iterator>(et,it_gen));
   
  }

  
  // initialize the variables of ntuples
  GenInfo.Initialize();
  bool _isFilled = false;
  bool _isFilled_Mu = false;
  

  for(phoEtMap::iterator i=myPhoEtMap.begin(); 
      i!= myPhoEtMap.end(); ++i)
    {
      
      if (GenInfo.Size>=MAX_GENS) {
	fprintf(stderr,"ERROR: number of generator photons exceeds the size of array.\n");
	exit(1);
      }

      reco::GenParticleCollection::const_iterator it_gen = i->second;
            

      float genpt   = it_gen->pt();
      int genMomPID = it_gen->mother()? it_gen->mother()->pdgId():_pdgCode;
	
    
      GenInfo.PID[GenInfo.Size] = it_gen->pdgId();
      GenInfo.MPID[GenInfo.Size] = genMomPID;
      GenInfo.Mass[GenInfo.Size] = it_gen->mass();
      GenInfo.Pt[GenInfo.Size] = genpt;
      GenInfo.Eta[GenInfo.Size] = it_gen->eta();
      GenInfo.Phi[GenInfo.Size] = it_gen->phi();
    
      bool isFromHardScattering = (abs(genMomPID) ==_pdgCode);
	
      bool isFromJet   = (abs(genMomPID)> 50);
	
      bool isFromQuark = (abs(genMomPID)< 10) && genMomPID!=0;
	
      bool isFromGluon = (abs(genMomPID)==21);

      // ==================================================
      // check trigger match
      // ==================================================

      // check L1 first
      if(myPhoL1Map.find(it_gen)!=myPhoL1Map.end())
	{
	  GenInfo.L1Pt[GenInfo.Size] = myPhoL1Map[it_gen]->pt();
	}

      // now check L3
      int trigMatchCode=0;
      float l3pt=-999;

      if(myPhoL3Map_HLTMu5.find(it_gen)!=myPhoL3Map_HLTMu5.end())
	{
	  trigMatchCode |= TRIGGER::HLT_Mu5;
	  l3pt = myPhoL3Map_HLTMu5[it_gen].pt();
	}
	
      // L1EG5
      if( myPhoL3Map_HLTL1EG5.find(it_gen)!= myPhoL3Map_HLTL1EG5.end())
	{
	  trigMatchCode |= TRIGGER::HLT_L1SingleEG5,
	    l3pt = myPhoL3Map_HLTL1EG5[it_gen].pt();
	}

      // L1EG8
      if( myPhoL3Map_HLTL1EG8.find(it_gen)!= myPhoL3Map_HLTL1EG8.end())
	{
	  trigMatchCode |= TRIGGER::HLT_L1SingleEG8,
	    l3pt = myPhoL3Map_HLTL1EG8[it_gen].pt();
	}

      // HLT_10_L1R
      if(myPhoL3Map_HLT10.find(it_gen)!=myPhoL3Map_HLT10.end())
	{
	  trigMatchCode |= TRIGGER::HLT_Photon10_L1R;
	  l3pt =myPhoL3Map_HLT10[it_gen].pt();
	}


      // HLT_15_L1R
      if(myPhoL3Map_HLT15.find(it_gen)!=myPhoL3Map_HLT15.end())
	{
	  trigMatchCode |= TRIGGER::HLT_Photon15_L1R;
	  l3pt =myPhoL3Map_HLT15[it_gen].pt();
	}

      // HLT_20Iso
      if(myPhoL3Map_HLT20Iso.find(it_gen)!=myPhoL3Map_HLT20Iso.end())
	{
	  trigMatchCode |= TRIGGER::HLT_Photon20_LooseEcalIso_TrackIso_L1R;
	  l3pt =myPhoL3Map_HLT20Iso[it_gen].pt();
	}
     
      // HLT_25_L1R
      if(myPhoL3Map_HLT25.find(it_gen)!=myPhoL3Map_HLT25.end())
	{
	  trigMatchCode |= TRIGGER::HLT_Photon25_L1R;
	  l3pt =myPhoL3Map_HLT25[it_gen].pt();
	}



      GenInfo.Trig[GenInfo.Size] = trigMatchCode;
      GenInfo.L3Pt[GenInfo.Size] = l3pt;     

      GenInfo.Size++;

      float et = it_gen->et();
      float eta = it_gen->eta();

      h_debug1->Fill(et);
      h_eta1  ->Fill(eta);
 
      
      h_debug2->Fill(et);
      h_eta2  ->Fill(eta);

      // Use SingleEG5 as the denominator

      if((_thisEvent_trigger & TRIGGER::HLT_L1SingleEG5) && !_isFilled)
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
	  
	  if(_thisEvent_trigger & TRIGGER::HLT_Photon15_L1R)
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

      if((_thisEvent_trigger & TRIGGER::HLT_Mu5) && !_isFilled_Mu)
	{
	  _isFilled_Mu=true;
	  h_mugetdeno->Fill(et);

	  if(_thisEvent_trigger & TRIGGER::HLT_Photon15_L1R)	    
	    h_mugetnumr->Fill(et);

	} // if passing muon trigger


    } // if number of photons > 0
      
  
  root->Fill();
 
  _nOut++;
  
}



void GenTrig::MatchGenPhoToL1(const edm::Event& iEvent)
{

  // look for Gen particle collection
  edm::Handle<reco::GenParticleCollection> GenHandle;  
  bool hasGenParticle = iEvent.getByLabel(_genLabel, GenHandle);
  if(!hasGenParticle)return;

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


  for (reco::GenParticleCollection::const_iterator it_gen = GenHandle->begin(); 
       it_gen!=GenHandle->end(); it_gen++){

    if(abs(it_gen->pdgId())!=_pdgCode || it_gen->status()!=1)continue;
    if(it_gen->pt() < 2.)continue;

    float ptMax = 0;
    std::vector<l1extra::L1EmParticle>::const_iterator maxPtIter;
    bool hasMatch = false;


    for( std::vector<l1extra::L1EmParticle>::const_iterator it_l1 = 
	   l1EmIso->begin(); 
	 it_l1 != l1EmIso->end(); it_l1++ ) {
    
      float deltaR = reco::deltaR(it_l1->eta(), it_l1->phi(),
				  it_gen->eta(), it_gen->phi());

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
				  it_gen->eta(), it_gen->phi());

      float thisL1Pt = it_l1->pt();

      if (deltaR < 0.5 && thisL1Pt > ptMax ) {
	
	ptMax = thisL1Pt;
	maxPtIter = it_l1;
	hasMatch = true;
      }

    } // end of L1 non-isolated trigger objects

    
    if(hasMatch)
      myPhoL1Map.insert(std::pair<reco::GenParticleCollection::const_iterator,
			std::vector<l1extra::L1EmParticle>::const_iterator>(it_gen,maxPtIter));    
  } // end of loop over photons

  return;

}

void GenTrig::MatchGenPhoToL3(const edm::Event& iEvent, 
			      std::string trgPath,
			      std::string tag,
			      phoL3Map& thisMap)
{
  thisMap.clear();

  // look for Gen particle collection
  edm::Handle<reco::GenParticleCollection> GenHandle;  
  bool hasGenParticle = iEvent.getByLabel(_genLabel, GenHandle);
  if(!hasGenParticle)return;

  // check for HLT objects
  const edm::InputTag hltsummaryTag = edm::InputTag("hltTriggerSummaryAOD","","HLT");
  bool has1E31TrigInfo=
    iEvent.getByLabel(hltsummaryTag, trgEvent);    
  if(!has1E31TrigInfo)return;

  const trigger::TriggerObjectCollection& TOC(trgEvent->getObjects());

  const edm::InputTag myLastFilter = edm::InputTag(trgPath,"",tag);
  if ( trgEvent->filterIndex(myLastFilter) >= trgEvent->sizeFilters())return;
  const trigger::Keys& keys( trgEvent->filterKeys( trgEvent->filterIndex(myLastFilter) ));

  for (reco::GenParticleCollection::const_iterator it_gen = GenHandle->begin(); 
       it_gen!=GenHandle->end(); it_gen++){

    if(abs(it_gen->pdgId())!=_pdgCode || it_gen->status()!=1)continue;
    if(it_gen->pt() < 2.)continue;

    float ptMax = 0;
    trigger::TriggerObject maxPtIter;
    bool hasMatch = false;

    for ( unsigned int hlto = 0; hlto < keys.size(); hlto++ ) {
      size_type hltf = keys[hlto];
      const trigger::TriggerObject L3obj(TOC[hltf]);
      
      float deltaR = reco::deltaR(L3obj.eta(), L3obj.phi(),
				  it_gen->eta(), it_gen->phi());
      float thisL3Pt = L3obj.pt();

      if (deltaR < 0.5 && thisL3Pt > ptMax) {
	maxPtIter = L3obj;
	ptMax = thisL3Pt;
	hasMatch = true;
      }
    } // end of loop over trigger objects

    if(hasMatch)
      thisMap.insert(std::pair<reco::GenParticleCollection::const_iterator,
		     trigger::TriggerObject>(it_gen,maxPtIter));    


  }// end of loop over photons

  return;
}

void GenTrig::TurnOnHLTBit(std::string trgPath, 
			   int         trgCode)
{
  TriggerNames trgName( *TrgResultsHandle);   
  int NTrigger = trgName.size();
  int tempIndex = (unsigned int)trgName.triggerIndex( trgPath); 
  if(tempIndex < NTrigger && TrgResultsHandle->accept(tempIndex))
    _thisEvent_trigger |= trgCode;
  return;
}




//define this as a plug-in
DEFINE_FWK_MODULE(GenTrig);


