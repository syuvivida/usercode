// -----------------------------------------------------
// RECOTrigger.cc 
// -- a very simple example analyzer for RECO data
// ----------------------------------------------------- 
// Shin-Shan Yu

#include "syu/MyObjectCounter/header/RECOTrigger.hh" 


RECOTrigger::RECOTrigger(const edm::ParameterSet& iConfig):
  _nIn(0), _nOut(0), _alg(iConfig)
{  
  _pdgCode = iConfig.getUntrackedParameter<int>("pdgCode",22);
}


RECOTrigger::~RECOTrigger()
{
}

void RECOTrigger::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  TFileDirectory results = TFileDirectory( fs->mkdir("RECOTrigger") );


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
  VtxInfo.Register(root);
  PhoInfo.Register(root);
  GenInfo.Register(root); 
  EvtInfo.Initialize();  
  VtxInfo.Initialize();
  PhoInfo.Initialize();
  GenInfo.Initialize();
  std::cout << "Focusing on PDG code = " << _pdgCode << std::endl;
}

void RECOTrigger::endJob() 
{
  std::cout << "RECOTrigger has " << _nIn << " input events and " << _nOut  << " output events" << endl;
}


// analyzing reconstructed electrons, muons, and photons

void RECOTrigger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  _nIn++;
  int thisEvent_trigger = 0;
  _alg.init(iEvent, true, false, true); 
  
  myPhoEtMap.clear();
  myPhoGenMap.clear();
  myPhoL1Map.clear();
  myPhoL3Map_HLTL1EG5.clear();
  myPhoL3Map_HLTL1EG8.clear();
  myPhoL3Map_HLT10.clear();
  myPhoL3Map_HLT15.clear();
  myPhoL3Map_HLT20Iso.clear();
  myPhoL3Map_HLT25.clear();
  myPhoL3Map_HLTMu5.clear();

  // sort photon particle collection
  myPhoEtMap = _alg.getPhoEtMap();

  if(myPhoEtMap.size()==0)return;

  // photon - to generator-level photon matching
  
  _alg.matchPartToGen<reco::Photon>(0.2, 1.0, myPhoEtMap, myPhoGenMap);

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

  //  _alg.turnOnHLTBit("HLT_Mu5",TRIGGER::HLT_Mu5);
  _alg.turnOnHLTBit("HLT_L1MuOpen",TRIGGER::HLT_Mu5);
  

 
  // Start to fill the main root branches
  // initialize the variables of ntuples
  EvtInfo.Initialize();  

  EvtInfo.isData = _alg.isData()? 1: 0;
  EvtInfo.RunNo  = iEvent.id().run();
  EvtInfo.EvtNo  = iEvent.id().event();
  thisEvent_trigger = _alg.getThisEventTriggerBit();
  EvtInfo.HLT    = thisEvent_trigger;

  // fill vertex information
  VtxInfo.Initialize(); 
  edm::Handle<reco::VertexCollection> vertexHandle = _alg.getVtxHandle();
  if(vertexHandle.isValid()){
    for (reco::VertexCollection::const_iterator vit=vertexHandle->begin(); 
	 vit!=vertexHandle->end(); vit++){
      if(vit->isFake())continue;
      if (VtxInfo.Size>= MAX_VTXS) {
	fprintf(stderr,"ERROR: number of vertices exceeds the size of array.\n");
	exit(1);
      }            

      VtxInfo.Index[VtxInfo.Size] = VtxInfo.Size;
      VtxInfo.pos[VtxInfo.Size][0] = vit->position().X();
      VtxInfo.pos[VtxInfo.Size][1] = vit->position().Y();
      VtxInfo.pos[VtxInfo.Size][2] = vit->position().Z();
      
      VtxInfo.poserr[VtxInfo.Size][0] = vit->xError();
      VtxInfo.poserr[VtxInfo.Size][1] = vit->yError();
      VtxInfo.poserr[VtxInfo.Size][2] = vit->zError();
      
      VtxInfo.rho[VtxInfo.Size]   = vit->position().Rho();
      VtxInfo.chi2[VtxInfo.Size]  = vit->chi2();
      VtxInfo.ndof[VtxInfo.Size]  = vit->ndof();
      VtxInfo.ntrks[VtxInfo.Size] = vit->tracksSize();

      VtxInfo.Size ++ ;
    
    } // end of loop over vertices

  } // if vertexHandle is valid

  // initialize the variables of ntuples

  // first fill generator information, not associated with reconstructed 
  // particle 
  GenInfo.Initialize();
  GenInfo.ptHat  = _alg.ptHat();
  edm::Handle<reco::GenParticleCollection> GenHandle = _alg.getGenHandle();
  if(GenHandle.isValid()){
    for( GenParticleCollection::const_iterator it_gen = 
	   GenHandle->begin(); 
	 it_gen != GenHandle->end() && GenInfo.Size < MAX_GENS; it_gen++ ) {

      if (GenInfo.Size>= MAX_GENS) {
	fprintf(stderr,"ERROR: number of generator photons exceeds the size of array.\n");
	exit(1);
      }            

      float genpt   = it_gen->pt();
      int genMomPID = it_gen->mother()? it_gen->mother()->pdgId():_pdgCode;
      if(it_gen->pdgId()!=_pdgCode || it_gen->pt() < 2.)continue;
   
      GenInfo.Index[GenInfo.Size] = GenInfo.Size;
      GenInfo.PID[GenInfo.Size]   = it_gen->pdgId();
      GenInfo.MPID[GenInfo.Size]  = genMomPID;
      GenInfo.Mass[GenInfo.Size]  = it_gen->mass();
      GenInfo.Pt[GenInfo.Size]    = genpt;
      GenInfo.Eta[GenInfo.Size]   = it_gen->eta();
      GenInfo.Phi[GenInfo.Size]   = it_gen->phi();
      GenInfo.Size ++;

    } // end of filling generator-level information    
  } // if have generator-level information

  PhoInfo.Initialize();
  bool _isFilled = false;
  bool _isFilled_Mu = false;
  
  for(partEtMap<reco::Photon>::Type::iterator i=myPhoEtMap.begin(); 
      i!= myPhoEtMap.end() && PhoInfo.Size < MAX_PHOTONS; ++i)
    {
      
      if (PhoInfo.Size >= MAX_PHOTONS) {
	fprintf(stderr,"ERROR: number of photons exceeds the size of array.\n");
	exit(1);
      }

      reco::PhotonCollection::const_iterator it_ph = i->second;
      
      int location = -1; 
      bool isBarrel = it_ph->isEB();
      bool isEndCap = it_ph->isEE();
      bool inAnyGap = it_ph->isEBEEGap() || 
	(it_ph->isEB() && it_ph->isEBGap()) || 
	(it_ph->isEE() && it_ph->isEEGap());

      if(inAnyGap)location=0;
      else if(isBarrel)location=1;
      else if(isEndCap)location=2;

      PhoInfo.Index[PhoInfo.Size]        = PhoInfo.Size;  
      PhoInfo.Location[PhoInfo.Size]     = location;
      PhoInfo.E    [PhoInfo.Size] 	 = it_ph->energy(); 
      PhoInfo.Et   [PhoInfo.Size] 	 = it_ph->et();     
      PhoInfo.Pz   [PhoInfo.Size] 	 = it_ph->pz();     
      PhoInfo.Eta  [PhoInfo.Size] 	 = it_ph->eta();    
      PhoInfo.Phi  [PhoInfo.Size] 	 = it_ph->phi();    
      PhoInfo.R9   [PhoInfo.Size] 	 = it_ph->r9();     

      float ecalIso = it_ph->ecalRecHitSumEtConeDR04();
      float hcalIso = it_ph->hcalTowerSumEtConeDR04();
      float trkIso  = it_ph->trkSumPtHollowConeDR04();

      PhoInfo.TrkPtSum[PhoInfo.Size]     = trkIso;
      PhoInfo.EcalRecHitEtSum[PhoInfo.Size]= ecalIso;
      PhoInfo.HcalTowerEtSum[PhoInfo.Size] = hcalIso;
      PhoInfo.HoverE[PhoInfo.Size]       = it_ph->hadronicOverEm();

      PhoInfo.SCE[PhoInfo.Size]          = it_ph->superCluster()->energy();
      PhoInfo.SCEta[PhoInfo.Size]        = it_ph->superCluster()->eta();
      PhoInfo.SCPhi[PhoInfo.Size]        = it_ph->superCluster()->phi();
      PhoInfo.SCEtaWidth[PhoInfo.Size]   = it_ph->superCluster()->etaWidth();
      PhoInfo.SCPhiWidth[PhoInfo.Size]   = it_ph->superCluster()->phiWidth();
	
      PhoInfo.SigEta[PhoInfo.Size]   = it_ph->sigmaEtaEta();
      PhoInfo.SigIEta[PhoInfo.Size]  = it_ph->sigmaIetaIeta();
 	
      float scet = it_ph->superCluster()->energy()/cosh(it_ph->superCluster()->eta());
      PhoInfo.SCEt[PhoInfo.Size]       = scet;

      PhoInfo.SCNCrystal[PhoInfo.Size] = it_ph->superCluster()->size();

      // ==================================================
      // started doing MC-reconstruction matching
      // ==================================================

      float genpt   = -999;
      float genmompt = -999;
      int genMomPID = 0;
	
      if(myPhoGenMap.find(it_ph)!=myPhoGenMap.end())
	{
	  GenParticleCollection::const_iterator it_gen = myPhoGenMap[it_ph];
	  PhoInfo.GenPID[PhoInfo.Size]       = _pdgCode;
	  genpt = it_gen->pt();
	  PhoInfo.GenPt[PhoInfo.Size]        = genpt;
	  if(it_gen->mother()){
	    genmompt  = it_gen->mother()->pt();
	    genMomPID = it_gen->mother()->pdgId();
	    PhoInfo.GenMomPID[PhoInfo.Size]  = genMomPID;
	    PhoInfo.GenMomPt[PhoInfo.Size]   = it_gen->mother()->pt();
	    if(it_gen->mother()->mother())
	      PhoInfo.GenGMomPID[PhoInfo.Size] = it_gen->mother()->mother()->pdgId();
	  }
	  else {
	    genMomPID = _pdgCode;
	    PhoInfo.GenMomPID[PhoInfo.Size]  = genMomPID;
	  }
	  
	} // if there is a matching

      bool isFromHardScattering = (abs(genMomPID) ==_pdgCode);

      bool isFromJet   = (abs(genMomPID)> 50);
   
      bool isFromQuark = (abs(genMomPID)< 10) && genMomPID!=0;
   
      bool isFromGluon = (abs(genMomPID)==21);

      // ==================================================
      // check trigger match
      // ==================================================

      // check L1 first
      if(myPhoL1Map.find(it_ph)!=myPhoL1Map.end())
	{
	  PhoInfo.L1Pt[PhoInfo.Size] = myPhoL1Map[it_ph]->pt();
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
	  trigMatchCode |= TRIGGER::HLT_L1SingleEG5,
	    l3pt = myPhoL3Map_HLTL1EG5[it_ph].pt();
	}

      // L1EG8
      if( myPhoL3Map_HLTL1EG8.find(it_ph)!= myPhoL3Map_HLTL1EG8.end())
	{
	  trigMatchCode |= TRIGGER::HLT_L1SingleEG8,
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


      PhoInfo.Trig[PhoInfo.Size] = trigMatchCode;
      PhoInfo.L3Pt[PhoInfo.Size] = l3pt;     

      bool passOffline = _alg.isMyLoosePhoton<reco::Photon>(it_ph);
      PhoInfo.IsLoose[PhoInfo.Size] = passOffline? 1: 0;

      PhoInfo.Size++;

      // need to pass the following cuts to be loose photons
      float et = it_ph->et();
      float eta = it_ph->eta();

      h_debug1->Fill(et);
      h_eta1  ->Fill(eta);
 
      // loose photon EB

      if(!passOffline)continue;
      
      h_debug2->Fill(et);
      h_eta2  ->Fill(eta);

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
	    h_mugetnumr->Fill(et);

	} // if passing muon trigger


    } // if number of photons > 0
  
  
  root->Fill();
 
  _nOut++;
  
}



//define this as a plug-in
DEFINE_FWK_MODULE(RECOTrigger);


