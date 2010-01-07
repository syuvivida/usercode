// -----------------------------------------------------
// Conv.cc 
// -- a very simple example analyzer for RECO data
// ----------------------------------------------------- 
// Shin-Shan Yu

#include "syu/MyObjectCounter/header/Conv.hh" 

Conv::Conv(const edm::ParameterSet& iConfig):
  _nIn(0), _nOut(0), _alg(iConfig)
{  
  _pdgCode = iConfig.getUntrackedParameter<int>("pdgCode",22);
}


Conv::~Conv()
{
}

void Conv::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  TFileDirectory results = TFileDirectory( fs->mkdir("Conv") );

  heta_all  = fs->make<TH1F>("heta_all","#eta of all photons",
			     200,-4.0,4.0);
  heta_original = fs->make<TH1F>("heta_original","#eta of all photons from hard scattering",
			     200,-4.0,4.0);
  heta_counte  = fs->make<TProfile>("heta_counte","Conversion rate (count electrons)",
				   200,-4.0,4.0);
  heta_conv = fs->make<TProfile>("heta_conv","Conversion rate (use PhotonMCTruth)", 200,-4.0,4.0);


  hpt_all  = fs->make<TH1F>("hpt_all","p_{T} of all photons", 50,0,200);
  hpt_original = fs->make<TH1F>("hpt_original","p_{T} of all photons from hard scattering", 50, 0, 200);
  hpt_counte  = fs->make<TProfile>("hpt_counte","Conversion rate (count electrons)", 50, 0, 200);
  hpt_conv = fs->make<TProfile>("hpt_conv","Conversion rate (use PhotonMCTruth)", 50, 0, 200);

  root = new TTree("root","root");
  EvtInfo.Register(root);  
  PhoInfo.Register(root);
  GenInfo.Register(root);
  EvtInfo.Initialize();  
  PhoInfo.Initialize();
  GenInfo.Initialize();

  thePhotonMCTruthFinder_ = new PhotonMCTruthFinder();

  std::cout << "Focusing on PDG code = " << _pdgCode << std::endl;
}

void Conv::endJob() 
{
  std::cout << "Conv has " << _nIn << " input events and " << _nOut  << " output events" << endl;
}


// analyzing reconstructed electrons, muons, and photons

void Conv::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  _nIn++;
  _alg.init(iEvent, true, true, false); 
  
  myPhoEtMap.clear();
  myPhoGenMap.clear();
  myEleEtMap.clear();
  myHardGen.clear();

  // sort photon particle collection
  myPhoEtMap = _alg.getPhoEtMap();
  myEleEtMap = _alg.getEleEtMap();
  myHardGen  = _alg.getHardGenVec();

  // Start to fill the main root branches
  // initialize the variables of ntuples
  EvtInfo.Initialize();  

  
  EvtInfo.RunNo  = iEvent.id().run();
  EvtInfo.EvtNo  = iEvent.id().event();
  EvtInfo.HLT    = myEleEtMap.size();

 
  // initialize the variables of ntuples

  // first fill generator information, not associated with reconstructed 
  // particle 
  GenInfo.Initialize();
  edm::Handle<reco::GenParticleCollection> GenHandle;  
  bool hasGenParticle = iEvent.getByLabel("genParticles", GenHandle);
  if(hasGenParticle){
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
    
      GenInfo.PID[GenInfo.Size]  = it_gen->pdgId();
      GenInfo.MPID[GenInfo.Size] = genMomPID;
      GenInfo.Mass[GenInfo.Size] = it_gen->mass();
      GenInfo.Pt[GenInfo.Size]   = genpt;
      GenInfo.Eta[GenInfo.Size]  = it_gen->eta();
      GenInfo.Phi[GenInfo.Size]  = it_gen->phi();
      GenInfo.Size ++;

    } // end of filling generator-level information    
  } // if have generator-level information

  PhoInfo.Initialize();

//   std::cout  << " MCPhotonAnalyzer Looking for MC truth " << "\n";
 
  //get simtrack info
  std::vector<SimTrack> theSimTracks;
  std::vector<SimVertex> theSimVertices;
 
  edm::Handle<SimTrackContainer> SimTk;
  edm::Handle<SimVertexContainer> SimVtx;
  iEvent.getByLabel("g4SimHits",SimTk);
  iEvent.getByLabel("g4SimHits",SimVtx);
 
  theSimTracks.insert(theSimTracks.end(),SimTk->begin(),SimTk->end());
  theSimVertices.insert(theSimVertices.end(),SimVtx->begin(),SimVtx->end());
  /*
  for(int i=0; i< myHardGen.size(); i++){

    reco::GenParticleCollection::const_iterator thisGen = myHardGen[i];
    cout << "Hard scattering photon momentum " << i << " = " << 
      thisGen->p4().x() << ", " << 
      thisGen->p4().y() << ", " << 
      thisGen->p4().z() << endl;
 
  }
  */
 
  // check if there is any electrons in simulated track bank

  bool hasElectron = false;

  for (std::vector<SimTrack>::iterator iSimTk = theSimTracks.begin(); iSimTk != theSimTracks.end(); ++iSimTk){
    if (  (*iSimTk).noVertex() ) continue;
    if ( abs((*iSimTk).type()) == 11){       
      hasElectron = true;
    }
    
  } // loop over simulated tracks


  std::vector<PhotonMCTruth> mcPhotons=thePhotonMCTruthFinder_->find (theSimTracks,  theSimVertices);
//   std::cout << " MCPhotonAnalyzer mcPhotons size " <<  mcPhotons.size() << std::endl;
   
  for ( std::vector<PhotonMCTruth>::const_iterator iPho=mcPhotons.begin(); iPho !=mcPhotons.end(); ++iPho ){

    float eta = iPho->fourMomentum().eta();
    heta_all->Fill(eta);
    
    float et = iPho->fourMomentum().et();
    hpt_all->Fill(et);

    reco::GenParticleCollection::const_iterator matchedPart;
    bool findOne = 
      _alg.getMatchGen( matchedPart, 22, 1, 
 			iPho->fourMomentum().x(),
 			iPho->fourMomentum().y(),
 			iPho->fourMomentum().z(),
 			iPho->primaryVertex().x(),
 			iPho->primaryVertex().y(),
 			iPho->primaryVertex().z()
			);
    
    if(!findOne)continue;

    if(matchedPart->mother() && 
       (matchedPart->mother()->pdgId()!=22 ||
	matchedPart->mother()->status()!=3))
      continue;

    heta_original->Fill(eta);
    hpt_original->Fill(et);

    float rate = iPho->isAConversion()? 1:0;
    heta_conv->Fill(eta,rate);
    hpt_conv->Fill(et,rate);
    
    float rate_e = hasElectron? 1:0;
    heta_counte->Fill(eta,rate_e);
    hpt_counte->Fill(et,rate_e);
    
  } // end of loop over MC truth

  root->Fill();
  _nOut++;

}



//define this as a plug-in
DEFINE_FWK_MODULE(Conv);


