// -----------------------------------------------------
// MyObjectCounter.cc 
// -- a very simple example analyzer for PAT tutorial
// ----------------------------------------------------- 
// Shin-Shan Yu


#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "syu/MyObjectCounter/header/MyObjectCounter.hh"


MyObjectCounter::MyObjectCounter(const edm::ParameterSet& iConfig):
  _nIn(0), _nOut(0),_alg(iConfig)

{  
 
}


MyObjectCounter::~MyObjectCounter()
{
}

void MyObjectCounter::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  TFileDirectory results = TFileDirectory( fs->mkdir("MyObjectCounter") );

  h_leppt  = fs->make<TH1F>("h_leppt","",50,0,200);
  h_lepeta = fs->make<TH1F>("h_lepeta","",24,-6.0,6.0);
  h_get    = fs->make<TH1F>("h_get", "", 50,0,200);


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
  std::cout << "MyObjectCounter has " << _nIn << " input events and " << _nOut  << " events" << endl;
}


// analyzing reconstructed electrons, muons, and photons

void MyObjectCounter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  _nIn++;
  _alg.init(iEvent, false, false, true, true, true); 

  // filter trigger path
  edm::Handle<TriggerResults> TrgResultsHandle;
  bool with_TriggerResults = iEvent.getByLabel(InputTag("TriggerResults::HLT"),TrgResultsHandle);
  if (with_TriggerResults) {
  
    TriggerNames TrgNames( *TrgResultsHandle );   
    unsigned int TrgIdx1 = TrgNames.triggerIndex( "HLT_Photon10_L1R" );
    unsigned int TrgIdx2 = TrgNames.triggerIndex( "HLT_Photon15_L1R" );

    vector<string> hlNames_ = TrgNames.triggerNames();
    if(_nIn < 3){
      for (size_t i=0; i<TrgNames.size(); ++i) {
	cout<<"HLT bit = "<<i<<"   "<<hlNames_[i]<<endl;
      }
    }
	
    if (!TrgResultsHandle->accept(TrgIdx1) &&
	!TrgResultsHandle->accept(TrgIdx2)) return; 	
  }

   
  // Start to fill the main root branches
  // initialize the variables of ntuples
  EvtInfo.Initialize();  

  
  EvtInfo.RunNo  = iEvent.id().run();
  EvtInfo.EvtNo  = iEvent.id().event();
  
  
  // look for photons

  // initialize the variables of ntuples
  PhoInfo.Initialize();
 
  
  partEtMap<pat::Photon>::Type thisPhoMap = _alg.getPatPhoEtMap();

  for( partEtMap<pat::Photon>::Type::iterator mapIter = thisPhoMap.begin(); 
       mapIter != thisPhoMap.end(); mapIter++ ) {
   
    if (PhoInfo.Size>=MAX_PHOTONS) {
      fprintf(stderr,"ERROR: number of photons exceeds the size of array.\n");
      exit(1);
    }
	
    myContainer<pat::Photon>::myIter it_ph= mapIter->second;

    h_get->Fill(it_ph->et());
	
    PhoInfo.Index[PhoInfo.Size]        = PhoInfo.Size;  
    PhoInfo.E    [PhoInfo.Size]        = it_ph->energy(); 
    PhoInfo.Et   [PhoInfo.Size]        = it_ph->et();     
    PhoInfo.Pz   [PhoInfo.Size]        = it_ph->pz();     
    PhoInfo.Eta  [PhoInfo.Size]        = it_ph->eta();    
    PhoInfo.Phi  [PhoInfo.Size]        = it_ph->phi();    
    PhoInfo.R9   [PhoInfo.Size]        = it_ph->r9();     
    PhoInfo.HoverE[PhoInfo.Size]       = it_ph->hadronicOverEm();
	
    PhoInfo.SCE[PhoInfo.Size]          = it_ph->superCluster()->energy();
    PhoInfo.SCEta[PhoInfo.Size]        = it_ph->superCluster()->eta();
    PhoInfo.SCPhi[PhoInfo.Size]        = it_ph->superCluster()->phi();
    PhoInfo.SCEtaWidth[PhoInfo.Size]   = it_ph->superCluster()->etaWidth();
    PhoInfo.SCPhiWidth[PhoInfo.Size]   = it_ph->superCluster()->phiWidth();
	
    float et = it_ph->superCluster()->energy()/cosh(it_ph->superCluster()->eta());
    PhoInfo.SCEt[PhoInfo.Size]       = et;

    PhoInfo.SCNCrystal[PhoInfo.Size] = it_ph->superCluster()->size();

    if(it_ph->genPhoton()){
	  
      PhoInfo.GenPID[PhoInfo.Size]       = 22;
      if(it_ph->genPhoton()->mother()){
	PhoInfo.GenMomPID[PhoInfo.Size]  = it_ph->genPhoton()->mother()->pdgId();
	PhoInfo.GenMomPt[PhoInfo.Size]   = it_ph->genPhoton()->mother()->pt();
	if(it_ph->genPhoton()->mother()->mother())
	  PhoInfo.GenGMomPID[PhoInfo.Size] = it_ph->genPhoton()->mother()->mother()->pdgId();
      }
	
    } // check generator-level matches
    else
	
      PhoInfo.GenPID[PhoInfo.Size]       = -1;
    
    PhoInfo.Size++;
  }
  
  
  // initialize the variables of ntuples
  LepInfo.Initialize();
  partEtMap<pat::Muon>::Type thisMuoMap = _alg.getPatMuoEtMap();
  
  for( partEtMap<pat::Muon>::Type::iterator mapIter = thisMuoMap.begin(); 
       mapIter != thisMuoMap.end(); mapIter++ ) {
          
    if (LepInfo.Size>=MAX_LEPTONS) {
      fprintf(stderr,"ERROR: number of leptons exceeds the size of array.\n");
      exit(1);
    }
	
    myContainer<pat::Muon>::myIter it_mu = mapIter->second;

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


  partEtMap<pat::Electron>::Type thisEleMap = _alg.getPatEleEtMap();

  for( partEtMap<pat::Electron>::Type::iterator mapIter = thisEleMap.begin(); 
       mapIter != thisEleMap.end(); mapIter++ ) {          
  
    if (LepInfo.Size>=MAX_LEPTONS) {
      fprintf(stderr,"ERROR: number of leptons exceeds the size of array.\n");
      exit(1);
    }
       
    myContainer<pat::Electron>::myIter it_el = mapIter->second;

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

    h_leppt -> Fill(it_el->pt());
    h_lepeta-> Fill(it_el->eta());

  } // if there are muons or electrons


  
  root->Fill();
  _nOut++;
  
}


//define this as a plug-in
DEFINE_FWK_MODULE(MyObjectCounter);


