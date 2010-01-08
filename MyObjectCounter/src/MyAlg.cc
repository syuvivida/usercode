
#include "syu/MyObjectCounter/header/MyAlg.hh" 


MyAlg::MyAlg( const edm::ParameterSet & iConfig )  :
  _parameters( iConfig ), _event_trigger(0),_isData(true),_ptHat(-999.)
{
  _dumpHEP = iConfig.getUntrackedParameter<bool>("dumpHEP", false);
  _pdgCode = iConfig.getUntrackedParameter<int>("pdgCode",  22);
  _trigDeltaRMax  = iConfig.getUntrackedParameter<double>("trigDR",  0.5);
}

//_____________________________________________________________________________
MyAlg::~MyAlg() {
}

MyAlg& MyAlg::operator=( const MyAlg& original)
{
  if ( this == 0 ) {
    std::cerr << "EXITING: trying to use 'operator=' with an MatrixMeschach for which memory space is not reserved!!!!" << std::endl;
    std::exception();
  }
  
  if ( this != &original ) {
    _phoHandle        = original._phoHandle;
    _eleHandle        = original._eleHandle;
    _genHandle        = original._genHandle;
    _hepMCHandle      = original._hepMCHandle;
    _genEventHandle   = original._genEventHandle;
    _l1EmNonIsoHandle = original._l1EmNonIsoHandle;
    _l1EmIsoHandle    = original._l1EmIsoHandle;
    _trgEventHandle   = original._trgEventHandle;
    _trgResultsHandle = original._trgResultsHandle;
    _vertexHandle     = original._vertexHandle;
    _hardGenParticle  = original._hardGenParticle;
    _phoEtMap         = original._phoEtMap;
    _eleEtMap         = original._eleEtMap;
    _genEtMap         = original._genEtMap;

    _patPhoHandle     = original._patPhoHandle;
    _patEleHandle     = original._patEleHandle;
    _patMuoHandle     = original._patMuoHandle;

    _patPhoEtMap      = original._patPhoEtMap;
    _patEleEtMap      = original._patEleEtMap;
    _patMuoEtMap      = original._patMuoEtMap;

    _parameters       = original._parameters;
    _event_trigger    = original._event_trigger;
    _isData           = original._isData;
    _ptHat            = original._ptHat;
    _dumpHEP          = original._dumpHEP;
    _pdgCode          = original._pdgCode;
    _trigDeltaRMax    = original._trigDeltaRMax;
   
  }
   return *this;
}



void MyAlg::init(const edm::Event &event,
		 bool doPho,
		 bool doEle,
		 bool doHLT,
		 bool doPAT)
{
  _event_trigger = 0;
  _genEtMap.clear();
  _phoEtMap.clear();
  _eleEtMap.clear();
  _patPhoEtMap.clear();
  _patEleEtMap.clear();
  _patMuoEtMap.clear();
  _hardGenParticle.clear();
  _isData = event.isRealData();

  getHandles(event, doPho, doEle, doHLT, doPAT);
  if(_dumpHEP)dumpGenInfo(event);

}

void MyAlg::getHandles(const edm::Event  & event,
		       bool doPho,
		       bool doEle,
		       bool doHLT,
		       bool doPAT)

{
 
  if(doPho){
    edm::InputTag photonName       = _parameters.getParameter<edm::InputTag>("phoLabel"  );
    event.getByLabel(photonName, _phoHandle);
    sortParticles<reco::Photon>(_phoHandle, _phoEtMap);
 
  }

  if(doEle){
    edm::InputTag electronName     = _parameters.getParameter<edm::InputTag>("eleLabel");
    event.getByLabel(electronName, _eleHandle);
    sortParticles<reco::GsfElectron>(_eleHandle, _eleEtMap);
  }

  if(isMC()){
    
    event.getByLabel("generator",_hepMCHandle);
    event.getByLabel("generator",_genEventHandle);
    _ptHat = _genEventHandle.isValid() && _genEventHandle->hasBinningValues() ? 
      _genEventHandle->binningValues()[0]: -999.;

    edm::InputTag genParticlesName = _parameters.getParameter<edm::InputTag>("genLabel");
    event.getByLabel(genParticlesName, _genHandle);
    sortGenParticles();
  }

  if(doHLT){
    edm::InputTag HLTName          = _parameters.getParameter<edm::InputTag>("HLTLabel");
    event.getByLabel(HLTName, _trgResultsHandle);
  }

  const edm::InputTag l1EmNonIsoTag = edm::InputTag("hltL1extraParticles","NonIsolated");
  const edm::InputTag l1EmIsoTag = edm::InputTag("hltL1extraParticles","Isolated");
  const edm::InputTag hltsummaryTag = edm::InputTag("hltTriggerSummaryAOD","","HLT");

  event.getByLabel(l1EmNonIsoTag, _l1EmNonIsoHandle);
  event.getByLabel(l1EmIsoTag, _l1EmIsoHandle);
  event.getByLabel(hltsummaryTag, _trgEventHandle);

  const edm::InputTag vtxCollTag = edm::InputTag("offlinePrimaryVertices");
  event.getByLabel(vtxCollTag, _vertexHandle);

  if(doPAT)
    {
      edm::InputTag photonName  = 
	_parameters.getParameter<edm::InputTag>("patPho");
      event.getByLabel(photonName, _patPhoHandle);
      sortParticles<pat::Photon>(_patPhoHandle, _patPhoEtMap);

      edm::InputTag electronName  = 
	_parameters.getParameter<edm::InputTag>("patEle");
      event.getByLabel(electronName, _patEleHandle);
      sortParticles<pat::Electron>(_patEleHandle, _patEleEtMap);

      edm::InputTag muonName = 
	_parameters.getParameter<edm::InputTag>("patMuo");
      event.getByLabel(muonName, _patMuoHandle);
      sortParticles<pat::Muon>(_patMuoHandle, _patMuoEtMap);
      
    }

}


//_____________________________________________________________________________
// dump generator-level information
//_____________________________________________________________________________


void MyAlg::dumpGenInfo(const edm::Event& iEvent)
{

  if(!_genHandle.isValid())return;

  std::cout << "============================ " << 
    "Run " <<  iEvent.id().run() << " Evt " <<  iEvent.id().event() <<
    " ============================= " << std::endl << std::endl;

  int genIndex = 0;

  printf("GenIndex  ");
  printf("PDG  ");
  printf("Status ");
  printf("Mother");
  printf("   ");
  printf("Mass ");
  printf("Energy ");
  printf("Pt ");
  printf("Px ");
  printf("Py ");
  printf("Pz ");
  printf("Vx ");
  printf("Vy ");
  printf("Vz ");
  printf("\n");


  for( reco::GenParticleCollection::const_iterator it_gen = _genHandle->begin(); 
       it_gen != _genHandle->end(); it_gen++ ) {

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
    printf("%9.3f",it_gen->p4().x());
    printf("%9.3f",it_gen->p4().y());
    printf("%9.3f",it_gen->p4().z());
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
  return;

}


//_____________________________________________________________________________
void MyAlg::turnOnHLTBit(std::string trgPath, int trgCode)
{
  if(!_trgResultsHandle.isValid())return;
  edm::TriggerNames trgName( *_trgResultsHandle);   
  int NTrigger = trgName.size();
  int tempIndex = (unsigned int)trgName.triggerIndex( trgPath); 
  if(tempIndex < NTrigger && _trgResultsHandle->accept(tempIndex))
    _event_trigger |= trgCode;
  return;
}


//_____________________________________________________________________________


void MyAlg::print()
{
  std::cout << "thisEvent_trigger = " << _event_trigger << std::endl;
  std::cout << "isData = " << _isData << std::endl;
  std::cout << "ptHat = " << _ptHat << std::endl;
  std::cout << "dumpHEP = " << _dumpHEP << std::endl;
  std::cout << "pdgCode = " << _pdgCode << std::endl;
  std::cout << "trigDeltaRMax = " << _trigDeltaRMax << std::endl;
}


//_____________________________________________________________________________


void MyAlg::sortGenParticles(){

  _genEtMap.clear();
  _hardGenParticle.clear();
  if(isData())return;
  if(!_genHandle.isValid())return;
  int countIndex = -1;
  for (reco::GenParticleCollection::const_iterator it_gen = 
	 _genHandle->begin(); it_gen!=_genHandle->end(); it_gen++){
    countIndex ++;
    if(abs(it_gen->pdgId())!=_pdgCode || it_gen->status()!=1)continue;
    int mPID = -1;
    int mStatus = -1;
    if(it_gen->mother())
      { 
	mPID    = it_gen->mother()->pdgId();
	mStatus = it_gen->mother()->status();
      }
    if(countIndex < 50 && (mPID==-1 || 
			   (abs(mPID)==_pdgCode && mStatus==3)))
      _hardGenParticle.push_back(it_gen);
    

    if(it_gen->pt() < 2.)continue;
    float et = it_gen->pt();
    _genEtMap.insert(std::pair<float,reco::GenParticleCollection::const_iterator>(et,it_gen));
  }
  return;
}


bool MyAlg::getMatchGen(reco::GenParticleCollection::const_iterator& tempIter,
			int pdgCode, int status, 
			float px1, float py1, float pz1, 
			float vx1, float vy1, float vz1)
{
  if(isData())return false;
  if(!_genHandle.isValid())return false;
  int countIndex = -1;

  bool findOne = false;
  for (reco::GenParticleCollection::const_iterator it_gen = 
	 _genHandle->begin(); it_gen!=_genHandle->end(); it_gen++){
    countIndex ++;
    if(abs(it_gen->pdgId())!= pdgCode || it_gen->status()!= status)continue;
    
    float dpx = fabs(px1 - it_gen->p4().x());
    float dpy = fabs(py1 - it_gen->p4().y());
    float dpz = fabs(pz1 - it_gen->p4().z());
    float dvx = fabs(vx1 - it_gen->vx());
    float dvy = fabs(vy1 - it_gen->vy());
    float dvz = fabs(vz1 - it_gen->vz());

    
    if(dpx < 0.001 && dpy < 0.001 && dpz < 0.001 && 
       (vx1 < -9999 || vy1 < -9999 || vz1 < -9999 || 
	(dvx < 0.001 && dvy < 0.001 && dvz < 0.001)))
      {
	tempIter = it_gen;
	findOne = true;
	break;
      }

  }

  return findOne;

}

