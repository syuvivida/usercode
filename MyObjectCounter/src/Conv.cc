// -----------------------------------------------------
// Conv.cc 
// -- a very simple example analyzer for RECO data
// ----------------------------------------------------- 
// Shin-Shan Yu

#include "syu/MyObjectCounter/header/Conv.hh" 
#include <algorithm>

Conv::Conv(const edm::ParameterSet& iConfig):
  _nIn(0), _nOut(0), _alg(iConfig)
{  
  _pdgCode = iConfig.getUntrackedParameter<int>("pdgCode",22);
}


Conv::~Conv()
{
}

void Conv::beginJob()
{
  edm::Service<TFileService> fs;
  TFileDirectory results = TFileDirectory( fs->mkdir("Conv") );

  heta_all  = fs->make<TH1F>("heta_all","#eta of all photons",
			     80,-4.0,4.0);
  heta_original = fs->make<TH1F>("heta_original","#eta of all photons from hard scattering",
			     80,-4.0,4.0);
  heta_counte  = fs->make<TProfile>("heta_counte","Conversion rate (count electrons)",
				   80,-4.0,4.0);
  heta_conv = fs->make<TProfile>("heta_conv","Conversion rate (use PhotonMCTruth)", 80,-4.0,4.0);


  he_original = fs->make<TH1F>("he_original","E of all photons from hard scattering", 100, 0, 500);
  he_counte  = fs->make<TProfile>("he_counte","Conversion rate (count electrons)", 100, 0, 500);
  he_conv = fs->make<TProfile>("he_conv","Conversion rate (use PhotonMCTruth)", 100, 0, 500);

  hpt_all  = fs->make<TH1F>("hpt_all","p_{T} of all photons", 50,0,250);
  hpt_original = fs->make<TH1F>("hpt_original","p_{T} of all photons from hard scattering", 50, 0, 250);
  hpt_counte  = fs->make<TProfile>("hpt_counte","Conversion rate (count electrons)", 50, 0, 250);
  hpt_conv = fs->make<TProfile>("hpt_conv","Conversion rate (use PhotonMCTruth)", 50, 0, 250);

  root = new TTree("root","root");
  EvtInfo.Register(root);  
  GenInfo.Register(root);
  EvtInfo.Initialize();  
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
  myConvPho.clear();

  // sort photon particle collection
  myPhoEtMap = _alg.getPhoEtMap();
  myEleEtMap = _alg.getEleEtMap();
  myHardGen  = _alg.getHardGenVec();
  myConvPho  = _alg.getConvPhoton();


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
  GenInfo.ptHat  = _alg.ptHat();
  int countIndex=0;
  edm::Handle<reco::GenParticleCollection> GenHandle = _alg.getGenHandle();  
  if(GenHandle.isValid()){
    for( GenParticleCollection::const_iterator it_gen = 
	   GenHandle->begin(); 
	 it_gen != GenHandle->end() && GenInfo.Size < MAX_GENS; it_gen++ ) {

      if (GenInfo.Size>= MAX_GENS) {
	fprintf(stderr,"ERROR: number of generator photons exceeds the size of array.\n");
	exit(1);
      }            
      countIndex++;
      float genpt   = it_gen->pt();
      int genMomPID = it_gen->mother()? it_gen->mother()->pdgId():_pdgCode;
      if(it_gen->pdgId()!=_pdgCode || it_gen->pt() < 2.)continue;
      
      float convRate = 0;
      if(std::find(myConvPho.begin(),myConvPho.end(),it_gen)!= myConvPho.end())
	convRate = 1.0;
      else convRate=0.0;
      //      if(genMomPID==22 && it_gen->mother()->status()==3)
      //	heta_counte->Fill(it_gen->eta(),convRate);

      GenInfo.Index[GenInfo.Size]    = GenInfo.Size;
      GenInfo.GenIndex[GenInfo.Size] = countIndex;
      GenInfo.PID[GenInfo.Size]  = it_gen->pdgId();
      GenInfo.MPID[GenInfo.Size] = genMomPID;
      GenInfo.Mass[GenInfo.Size] = it_gen->mass();
      GenInfo.Pt[GenInfo.Size]   = genpt;
      GenInfo.Eta[GenInfo.Size]  = it_gen->eta();
      GenInfo.Phi[GenInfo.Size]  = it_gen->phi();
      GenInfo.Size ++;

    } // end of filling generator-level information    
  } // if have generator-level information



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
  if(myHardGen.size()>1.0){
  for(int i=0; i< myHardGen.size(); i++){

    reco::GenParticleCollection::const_iterator thisGen = myHardGen[i];
    cout << "Hard scattering photon momentum " << i << " = " << 
      thisGen->p4().x() << ", " << 
      thisGen->p4().y() << ", " << 
      thisGen->p4().z() << endl;
 
  }
  
  _alg.dumpGenInfo(iEvent);
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

    float energy = iPho->fourMomentum().e();

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

    //    if(matchedPart->mother() && 
    //       (matchedPart->mother()->pdgId()!=22 ||
    //	matchedPart->mother()->status()!=3))
    //      continue;
    if(matchedPart->mother() && 
       matchedPart->mother()->pdgId()!=111)continue;

    heta_original->Fill(eta);
    hpt_original->Fill(et);
    he_original->Fill(energy);

    float rate = iPho->isAConversion()? 1:0;
    heta_conv->Fill(eta,rate);
    hpt_conv->Fill(et,rate);
    he_conv->Fill(energy,rate);
    
    float rate_e = hasElectron? 1:0;
    heta_counte->Fill(eta,rate_e);
    hpt_counte->Fill(et,rate_e);
    he_counte->Fill(energy,rate_e);
    
  } // end of loop over MC truth

  root->Fill();
  _nOut++;

}



//define this as a plug-in
DEFINE_FWK_MODULE(Conv);


