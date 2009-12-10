// -----------------------------------------------------
// Zee.cc 
// -- a very simple example analyzer for RECO data
// ----------------------------------------------------- 
// Shin-Shan Yu

#include "syu/MyObjectCounter/header/Zee.hh" 


Zee::Zee(const edm::ParameterSet& iConfig):
  _nIn(0), _nOut(0), _alg(iConfig)
{  

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
  std::cout << "Zee has " << _nIn << " input events and " << _nOut  << " output events" << endl;
}


// analyzing reconstructed electrons, muons, and photons

void Zee::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  _nIn++;
  int thisEvent_trigger = 0;
  _alg.init(iEvent, true, true, true, true); 

  myEleEtMap.clear();   // sorted electron Et map
  myPhoEtMap.clear();   // sorted photon Et map
  myElePhoMap.clear();  // sorted electron photon matching
  myEleGenMap.clear();  // generator-level electron and reconstructed electron mapping
  myPhoL1Map.clear(); 
  myPhoL3Map_HLTL1EG5.clear();
  myPhoL3Map_HLTL1EG8.clear();
  myPhoL3Map_HLT10.clear();
  myPhoL3Map_HLT15.clear();
  myPhoL3Map_HLT20Iso.clear();
  myPhoL3Map_HLT25.clear();
  myPhoL3Map_HLTMu5.clear();
  
  myEleEtMap = _alg.getEleEtMap();
  myPhoEtMap = _alg.getPhoEtMap();

//   if(myEleEtMap.size()==0)return;
//   if(myPhoEtMap.size()==0)return;


  // matching between electron and photon
  _alg.matchEleToPho<reco::GsfElectron, reco::Photon>( myEleEtMap,
						       myPhoEtMap,
						       myElePhoMap);

  // matching between electron and generator-level particles
  _alg.matchPartToGen<reco::GsfElectron>(0.5,1.0,myEleEtMap,myEleGenMap);
  
  // photon - to L1 matching
  _alg.matchPartToL1<reco::Photon>(myPhoEtMap, myPhoL1Map);

// HLT_L1EG5
  _alg.matchPartToL3<reco::Photon>("hltL1sRelaxedSingleEgammaEt5", "HLT",
				   myPhoEtMap, myPhoL3Map_HLTL1EG5);

  // HLT_L1EG8
  _alg.matchPartToL3<reco::Photon>("hltL1sRelaxedSingleEgammaEt8", "HLT",
				   myPhoEtMap, myPhoL3Map_HLTL1EG8);

  // HLT_10L1R
  _alg.matchPartToL3<reco::Photon>("hltL1NonIsoHLTNonIsoSinglePhotonEt10HcalIsolFilter",
				   "HLT", myPhoEtMap, myPhoL3Map_HLT10);
   
  // HLT_15L1R
  _alg.matchPartToL3<reco::Photon>("hltL1NonIsoHLTNonIsoSinglePhotonEt15HcalIsolFilter",
				   "HLT", myPhoEtMap, myPhoL3Map_HLT15);

  // HLT_20IsoL1R
  _alg.matchPartToL3<reco::Photon>("hltL1NonIsoHLTLEITISinglePhotonEt20TrackIsolFilter",
				   "HLT", myPhoEtMap, myPhoL3Map_HLT20Iso);

  // HLT_25L1R
  _alg.matchPartToL3<reco::Photon>("hltL1NonIsoHLTNonIsoSinglePhotonEt25HcalIsolFilter",
				   "HLT", myPhoEtMap, myPhoL3Map_HLT25);

  // HLT_MU5
  _alg.matchPartToL3<reco::Photon>("hltSingleMu5L3Filtered5", "HLT",
				   myPhoEtMap, myPhoL3Map_HLTMu5);


 
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

  _alg.turnOnHLTBit("HLT_Mu5",TRIGGER::HLT_Mu5);
  //  _alg.turnOnHLTBit("HLT_L1MuOpen",TRIGGER::HLT_Mu5);
  


  EvtInfo.Initialize();  
  EvtInfo.RunNo  = iEvent.id().run();
  EvtInfo.EvtNo  = iEvent.id().event();
  thisEvent_trigger = _alg.getThisEventTriggerBit();
  EvtInfo.HLT    = thisEvent_trigger;


  // fill generator information
  edm::Handle<reco::GenParticleCollection> GenHandle;  
  bool hasGenParticle = iEvent.getByLabel("genParticles", GenHandle); 
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
  typedef elePhoMap<reco::GsfElectron, reco::Photon>::Type::iterator mapIter;

  for(mapIter it_e=myElePhoMap.begin(); it_e!= myElePhoMap.end(); ++it_e)
    {
      it_count ++;

      int jt_count=-1;
      for(mapIter jt_e = myElePhoMap.begin(); jt_e != myElePhoMap.end(); 
	  ++jt_e){
	
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
      
	  bool passOffline = _alg.isMyLoosePhoton<reco::Photon>(it_ph);
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


//define this as a plug-in
DEFINE_FWK_MODULE(Zee);


