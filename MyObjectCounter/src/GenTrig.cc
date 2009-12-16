// -----------------------------------------------------
// GenTrig.cc 
// -- a very simple example analyzer for Gen data
// ----------------------------------------------------- 
// Shin-Shan Yu

#include "syu/MyObjectCounter/header/GenTrig.hh"




GenTrig::GenTrig(const edm::ParameterSet& iConfig):
  _nIn(0), _nOut(0), _alg(iConfig)
{  
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

  std::cout << "Focusing on PDG code = " << _pdgCode << std::endl;
}

void GenTrig::endJob() 
{
  std::cout << "GenTrig has " << _nIn << " input events and " << _nOut  << " output events" << endl;
}


// analyzing reconstructed electrons, muons, and photons

void GenTrig::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  _nIn++;
  int thisEvent_trigger = 0;
  _alg.init(iEvent, false, false, true); 
  

  myPhoEtMap.clear();
  myPhoL1Map.clear();
  myPhoL3Map_HLTL1EG5.clear();
  myPhoL3Map_HLTL1EG8.clear();
  myPhoL3Map_HLT10.clear();
  myPhoL3Map_HLT15.clear();
  myPhoL3Map_HLT20Iso.clear();
  myPhoL3Map_HLT25.clear();
  myPhoL3Map_HLTMu5.clear();

  // sort Gen particle collection
  myPhoEtMap = _alg.getGenEtMap();

  if(myPhoEtMap.size()==0)return;

  // photon - to L1 matching
  _alg.matchPartToL1<reco::GenParticle>(myPhoEtMap, myPhoL1Map);

// HLT_L1EG5
  _alg.matchPartToL3<reco::GenParticle>("hltL1sRelaxedSingleEgammaEt5", "HLT",
		     myPhoEtMap, myPhoL3Map_HLTL1EG5);

  // HLT_L1EG8
  _alg.matchPartToL3<reco::GenParticle>("hltL1sRelaxedSingleEgammaEt8", "HLT",
		     myPhoEtMap, myPhoL3Map_HLTL1EG8);


  // HLT_10L1R
  _alg.matchPartToL3<reco::GenParticle>("hltL1NonIsoHLTNonIsoSinglePhotonEt10HcalIsolFilter",
		     "HLT", myPhoEtMap, myPhoL3Map_HLT10);
   
  // HLT_15L1R
  _alg.matchPartToL3<reco::GenParticle>("hltL1NonIsoHLTNonIsoSinglePhotonEt15HcalIsolFilter",
		     "HLT", myPhoEtMap, myPhoL3Map_HLT15);


  // HLT_20IsoL1R
  _alg.matchPartToL3<reco::GenParticle>("hltL1NonIsoHLTLEITISinglePhotonEt20TrackIsolFilter",
		     "HLT", myPhoEtMap, myPhoL3Map_HLT20Iso);

  // HLT_25L1R
  _alg.matchPartToL3<reco::GenParticle>("hltL1NonIsoHLTNonIsoSinglePhotonEt25HcalIsolFilter",
		     "HLT", myPhoEtMap, myPhoL3Map_HLT25);

  // HLT_MU5
  _alg.matchPartToL3<reco::GenParticle>("hltSingleMu5L3Filtered5", "HLT",
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

  
  EvtInfo.RunNo  = iEvent.id().run();
  EvtInfo.EvtNo  = iEvent.id().event();
  thisEvent_trigger = _alg.getThisEventTriggerBit();
  EvtInfo.HLT    = thisEvent_trigger;


  
  // initialize the variables of ntuples
  GenInfo.Initialize();
  bool _isFilled = false;
  bool _isFilled_Mu = false;
  

  for(partEtMap<reco::GenParticle>::Type::iterator 
	i=myPhoEtMap.begin(); 
      i!= myPhoEtMap.end() && GenInfo.Size < MAX_GENS; ++i)
    {
      
      if (GenInfo.Size >= MAX_GENS) {
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
DEFINE_FWK_MODULE(GenTrig);


