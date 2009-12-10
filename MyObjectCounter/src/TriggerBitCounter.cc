// -----------------------------------------------------
// TriggerBitCounter.cc 
// -- a very simple example analyzer for PAT tutorial
// ----------------------------------------------------- 
// Shin-Shan Yu

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "syu/MyObjectCounter/header/TriggerBitCounter.hh"

TriggerBitCounter::TriggerBitCounter(const edm::ParameterSet& iConfig):
  _nIn(0), _nOut(0), _alg(iConfig)

{  
 
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
}

void TriggerBitCounter::endJob() 
{
  std::cout << "TriggerBitCounter has " << _nIn << " input events and " << _nOut  << " output events" << endl;
}


// analyzing reconstructed electrons, muons, and photons

void TriggerBitCounter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  _nIn++;
  _alg.init(iEvent, false, false, true, true, true); 

  int thisEvent_trigger =0;


  _alg.turnOnHLTBit("HLT_L1SingleEG5",TRIGGER::HLT_L1SingleEG5);

  _alg.turnOnHLTBit("HLT_L1SingleEG8",TRIGGER::HLT_L1SingleEG8);
  
  _alg.turnOnHLTBit("HLT_Photon10_L1R",TRIGGER::HLT_Photon10_L1R);
  
  _alg.turnOnHLTBit("HLT_Photon10_LooseEcalIso_TrackIso_L1R",
			TRIGGER::HLT_Photon10_LooseEcalIso_TrackIso_L1R);
    
  _alg.turnOnHLTBit("HLT_Photon15_L1R",TRIGGER::HLT_Photon15_L1R);
    
  _alg.turnOnHLTBit("HLT_Photon15_TrackIso_L1R",TRIGGER::HLT_Photon15_TrackIso_L1R);
  
  _alg.turnOnHLTBit("HLT_Photon15_LooseEcalIso_L1R",TRIGGER::HLT_Photon15_LooseEcalIso_L1R);
    
  _alg.turnOnHLTBit("HLT_Photon20_LooseEcalIso_TrackIso_L1R",TRIGGER::HLT_Photon20_LooseEcalIso_TrackIso_L1R);

  _alg.turnOnHLTBit("HLT_Photon25_L1R",TRIGGER::HLT_Photon25_L1R);

  //  _alg.turnOnHLTBit("HLT_Mu5",TRIGGER::HLT_Mu5);
  _alg.turnOnHLTBit("HLT_L1MuOpen",TRIGGER::HLT_Mu5);


  
  // Start to fill the main root branches
  // initialize the variables of ntuples
  EvtInfo.Initialize();  

  
  EvtInfo.RunNo  = iEvent.id().run();
  EvtInfo.EvtNo  = iEvent.id().event();
  thisEvent_trigger = _alg.getThisEventTriggerBit();
  EvtInfo.HLT    = thisEvent_trigger;

  
  // initialize the variables of ntuples
  PhoInfo.Initialize();
  
  bool _isFilled = false;
  bool _isFilled_Mu = false;

  partEtMap<pat::Photon>::Type thisPhoMap = _alg.getPatPhoEtMap();

  for( partEtMap<pat::Photon>::Type::iterator mapIter = thisPhoMap.begin(); 
       mapIter != thisPhoMap.end(); mapIter++ ) {
    
    if (PhoInfo.Size>=MAX_PHOTONS) {
      fprintf(stderr,"ERROR: number of photons exceeds the size of array.\n");
      exit(1);
    }
	
    myContainer<pat::Photon>::myIter it_ph= mapIter->second;
    
    float ecalIso = it_ph->ecalRecHitSumEtConeDR04();
    float hcalIso = it_ph->hcalTowerSumEtConeDR04();
    float trkIso  = it_ph->trkSumPtHollowConeDR04();

    PhoInfo.Index[PhoInfo.Size]          = PhoInfo.Size;  
    PhoInfo.E    [PhoInfo.Size] 	 = it_ph->energy(); 
    PhoInfo.Et   [PhoInfo.Size] 	 = it_ph->et();     
    PhoInfo.Pz   [PhoInfo.Size] 	 = it_ph->pz();     
    PhoInfo.Eta  [PhoInfo.Size] 	 = it_ph->eta();    
    PhoInfo.Phi  [PhoInfo.Size] 	 = it_ph->phi();    
    PhoInfo.R9   [PhoInfo.Size] 	 = it_ph->r9();     
    PhoInfo.TrkPtSum[PhoInfo.Size]       = trkIso;
    PhoInfo.EcalRecHitEtSum[PhoInfo.Size]= ecalIso;
    PhoInfo.HcalTowerEtSum[PhoInfo.Size] = hcalIso;
    PhoInfo.HoverE[PhoInfo.Size]         = it_ph->hadronicOverEm();

    PhoInfo.SCE[PhoInfo.Size]            = it_ph->superCluster()->energy();
    PhoInfo.SCEta[PhoInfo.Size]          = it_ph->superCluster()->eta();
    PhoInfo.SCPhi[PhoInfo.Size]          = it_ph->superCluster()->phi();
    PhoInfo.SCEtaWidth[PhoInfo.Size]     = it_ph->superCluster()->etaWidth();
    PhoInfo.SCPhiWidth[PhoInfo.Size]     = it_ph->superCluster()->phiWidth();
	
    float scet = it_ph->superCluster()->energy()/cosh(it_ph->superCluster()->eta());
    PhoInfo.SCEt[PhoInfo.Size]           = scet;

    PhoInfo.SCNCrystal[PhoInfo.Size]     = it_ph->superCluster()->size();

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
      
    bool passOffline = _alg.isMyLoosePhoton<pat::Photon>(it_ph);

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

  root->Fill();
  _nOut++;
  
}


//define this as a plug-in
DEFINE_FWK_MODULE(TriggerBitCounter);


