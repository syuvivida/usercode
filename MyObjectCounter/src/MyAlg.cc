
#include "syu/MyObjectCounter/header/MyAlg.hh" 


MyAlg::MyAlg( const edm::ParameterSet & iConfig )  :
  _parameters( iConfig ), _event_trigger(0)
{
  _dumpHEP = iConfig.getUntrackedParameter<bool>("dumpHEP", false);
  _pdgCode = iConfig.getUntrackedParameter<int>("pdgCode",  22);
  _deltaRMax  = iConfig.getUntrackedParameter<double>("dR",  0.5);
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
    _parameters       = original._parameters;
    _phoHandle        = original._phoHandle;
    _eleHandle        = original._eleHandle;
    _genHandle        = original._genHandle;
    _l1EmNonIsoHandle = original._l1EmNonIsoHandle;
    _l1EmIsoHandle    = original._l1EmIsoHandle;
    _trgEventHandle   = original._trgEventHandle;
    _trgResultsHandle = original._trgResultsHandle;
    _phoEtMap         = original._phoEtMap;
    _eleEtMap         = original._eleEtMap;
    _genEtMap         = original._genEtMap;

    _patPhoHandle     = original._patPhoHandle;
    _patEleHandle     = original._patEleHandle;
    _patMuoHandle     = original._patMuoHandle;

    _patPhoEtMap      = original._patPhoEtMap;
    _patEleEtMap      = original._patEleEtMap;
    _patMuoEtMap      = original._patMuoEtMap;

    _dumpHEP          = original._dumpHEP;
    _pdgCode          = original._pdgCode;
    _deltaRMax        = original._deltaRMax;
    _event_trigger    = original._event_trigger;
   
  }
   return *this;
}

void MyAlg::print()
{
  std::cout << "dumpHEP = " << _dumpHEP << std::endl;
  std::cout << "pdgCode = " << _pdgCode << std::endl;
  std::cout << "deltaRMax = " << _deltaRMax << std::endl;
  std::cout << "thisEvent_trigger = " << _event_trigger << std::endl;
}


void MyAlg::init(const edm::Event &event,
		 bool doPho,
		 bool doEle,
		 bool doGen,
		 bool doHLT,
		 bool doPAT)
{
  _event_trigger = 0;
  _genEtMap.clear();
  _phoEtMap.clear();
  _eleEtMap.clear();
  getHandles(event, doPho, doEle, doGen, doHLT, doPAT);
  if(_dumpHEP)dumpGenInfo(event);

}

void MyAlg::getHandles(const edm::Event  & event,
		       bool doPho,
		       bool doEle,
		       bool doGen,
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

  if(doGen){
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


void MyAlg::sortGenParticles(){

  _genEtMap.clear();
  if(!_genHandle.isValid())return;
  int count = 0;
  for (reco::GenParticleCollection::const_iterator it_gen = 
	 _genHandle->begin(); 
       it_gen!=_genHandle->end() && count < MAX_GENS; it_gen++){
    if(abs(it_gen->pdgId())!=_pdgCode || it_gen->status()!=1)continue;
    if(it_gen->pt() < 2.)continue;
    count ++ ;
    float et = it_gen->pt();
    _genEtMap.insert(std::pair<float,reco::GenParticleCollection::const_iterator>(et,it_gen));
  }
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
// dump generator-level information


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
  printf("Mother PID");
  printf("   ");
  printf("Mass ");
  printf("Energy ");
  printf("Pt");
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
