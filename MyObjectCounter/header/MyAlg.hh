#ifndef MYALG_HH
#define MYALG_HH

#include <iostream>
#include <map>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "syu/MyObjectCounter/header/format.hh" 
#include "syu/MyObjectCounter/header/trigformat.hh"


using namespace std;
using namespace edm;
using namespace reco;
using namespace trigger;
using namespace math;
using namespace ROOT;


template <class T>
struct myContainer
{
  typedef typename std::vector<T> container;
  typedef typename container::const_iterator myIter;
};

template <class T>
struct partEtMap
{
  typedef typename std::map< float, typename myContainer<T>::myIter, 
			     std::greater<float> > Type;
};

template <class T>
struct partL1Map
{
  typedef typename std::map< typename myContainer<T>::myIter, 
			     l1extra::L1EmParticleCollection::const_iterator > 
  Type;
};


template <class T>
struct partL3Map
{
  typedef typename std::map< typename myContainer<T>::myIter, 
			     trigger::TriggerObject > Type;
  
};

class MyAlg {
//-----------------------------------------------------------------------------
//  methods
//-----------------------------------------------------------------------------
public:
  MyAlg(edm::ParameterSet const & parameters);
  MyAlg(){}
  ~MyAlg();

  MyAlg& operator=( const MyAlg& original);

  // initialization, calls getHandles and dumpGenInfo

  void init(const edm::Event & event,
	    bool doPho,
	    bool doEle,
	    bool doGen,
	    bool doHLT,
	    bool doPAT=false);

  void getHandles(const edm::Event  & event, 
		  bool doPho,
		  bool doEle,
		  bool doGen,
		  bool doHLT,
		  bool doPAT=false);

  void dumpGenInfo(const edm::Event&); 

  // turn on HLT bit information
  void turnOnHLTBit(std::string trgPath, int trgCode);
  int  getThisEventTriggerBit(){return _event_trigger;}

  // print out the private data members
  void print();

  void       sortGenParticles(); // with a PDG code cut
  template<class T> void sortParticles(
	       const edm::Handle<typename myContainer<T>::container >& handle,
	       typename partEtMap<T>::Type& etmap);

  template <class T> void matchPartToL1(typename partEtMap<T>::Type& etmap, 
					typename partL1Map<T>::Type& mymap );

  template <class T> void matchPartToL3(std::string trgPath, std::string tag,
					typename partEtMap<T>::Type& etmap,
					typename partL3Map<T>::Type& mymap);

  partEtMap<reco::Photon>::Type      getPhoEtMap(){return _phoEtMap;}
  partEtMap<reco::GsfElectron>::Type getEleEtMap(){return _eleEtMap;}
  partEtMap<reco::GenParticle>::Type getGenEtMap(){return _genEtMap;}

  partEtMap<pat::Photon>::Type       getPatPhoEtMap(){return _patPhoEtMap;}
  partEtMap<pat::Electron>::Type     getPatEleEtMap(){return _patEleEtMap;}
  partEtMap<pat::Muon>::Type         getPatMuoEtMap(){return _patMuoEtMap;}

  
private: 

   // Keep a version of the parameter set in question
  edm::Handle<reco::PhotonCollection>          _phoHandle;
  edm::Handle<reco::GsfElectronCollection>     _eleHandle;
  edm::Handle<reco::GenParticleCollection>     _genHandle;
  edm::Handle<l1extra::L1EmParticleCollection> _l1EmNonIsoHandle;
  edm::Handle<l1extra::L1EmParticleCollection> _l1EmIsoHandle;
  edm::Handle<trigger::TriggerEvent>           _trgEventHandle;
  edm::Handle<edm::TriggerResults>             _trgResultsHandle;


  partEtMap<reco::Photon>::Type                   _phoEtMap;
  partEtMap<reco::GsfElectron>::Type              _eleEtMap;
  partEtMap<reco::GenParticle>::Type              _genEtMap;


  // handles to pat 
  edm::Handle<std::vector<pat::Photon> >          _patPhoHandle;
  edm::Handle<std::vector<pat::Electron> >        _patEleHandle;
  edm::Handle<std::vector<pat::Muon> >            _patMuoHandle;

  partEtMap<pat::Photon>::Type                    _patPhoEtMap;
  partEtMap<pat::Electron>::Type                  _patEleEtMap;
  partEtMap<pat::Muon>::Type                      _patMuoEtMap;


  edm::ParameterSet         _parameters;
  bool _dumpHEP;
  int  _pdgCode;
  double _deltaRMax;
  int _event_trigger;

};

// === inline definition of template functions (has to be in the header files)
//

template < class T> void MyAlg::sortParticles(
const edm::Handle< typename myContainer<T>::container >& handle,  
typename partEtMap<T>::Type& sortedEtMap)
{ 
  sortedEtMap.clear(); 
  if(!handle.isValid())return;
  typename myContainer<T>::myIter it_part; 
  for ( it_part = handle->begin(); it_part != handle->end(); it_part++){
    float et = it_part->pt(); 
    sortedEtMap.insert(std::pair<float,
		       typename myContainer<T>::myIter>(et,it_part));  
  }
  return;
}


template <class T> void MyAlg::matchPartToL1(
				 typename partEtMap<T>::Type& etmap,
				 typename partL1Map<T>::Type& mymap)
{
  mymap.clear();

  if(!_l1EmIsoHandle.isValid() 
     && !_l1EmNonIsoHandle.isValid())return;

  typedef typename partEtMap<T>::Type::iterator mapIter;
 
  for(mapIter it_part=etmap.begin(); 
      it_part!= etmap.end(); ++it_part)
    {
      float ptMax = 0;
      std::vector<l1extra::L1EmParticle>::const_iterator maxPtIter;
      bool hasMatch = false;

      for( std::vector<l1extra::L1EmParticle>::const_iterator it_l1 = 
	   _l1EmIsoHandle->begin(); 
	   it_l1 != _l1EmIsoHandle->end(); it_l1++ ) {
    
	float deltaR = reco::deltaR(it_l1->eta(), it_l1->phi(),
				    it_part->second->eta(), 
				    it_part->second->phi());
	
	float thisL1Pt = it_l1->pt();
	
	if (deltaR < _deltaRMax && thisL1Pt > ptMax ) {
	  
	  ptMax = thisL1Pt;
	  maxPtIter = it_l1;
	  hasMatch = true;
      }

    } // end of L1 isolated trigger objects


    for( std::vector<l1extra::L1EmParticle>::const_iterator it_l1 = 
	  _l1EmNonIsoHandle->begin(); 
	 it_l1 != _l1EmNonIsoHandle->end(); it_l1++ ) {
    
      float deltaR = reco::deltaR(it_l1->eta(), it_l1->phi(),
				  it_part->second->eta(), 
				  it_part->second->phi());

      float thisL1Pt = it_l1->pt();

      if (deltaR < _deltaRMax && thisL1Pt > ptMax ) {
	
	ptMax = thisL1Pt;
	maxPtIter = it_l1;
	hasMatch = true;
      }

    } // end of L1 non-isolated trigger objects

    
    if(hasMatch)
      mymap.insert(std::pair<typename myContainer<T>::myIter,
		   l1extra::L1EmParticleCollection::const_iterator>(it_part->second,maxPtIter));    
  } // end of loop over photons

  return;




}


//_____________________________________________________________________________
//
template <class T> void MyAlg::matchPartToL3(std::string trgPath,
					     std::string tag,
					    typename partEtMap<T>::Type& etmap,
					    typename partL3Map<T>::Type& mymap)

{
  mymap.clear();

  // check for HLT objects
  if(!_trgEventHandle.isValid())return;


  const trigger::TriggerObjectCollection& 
    TOC(_trgEventHandle->getObjects());
  const edm::InputTag myLastFilter = edm::InputTag(trgPath,"",tag);
  trigger::size_type thisFilterIndex = 
    _trgEventHandle->filterIndex(myLastFilter);
  if ( thisFilterIndex >= _trgEventHandle->sizeFilters())return;
  const trigger::Keys& keys( _trgEventHandle->filterKeys(thisFilterIndex));

  typedef typename partEtMap<T>::Type::iterator mapIter;
  for(mapIter it_part=etmap.begin(); 
      it_part!= etmap.end(); ++it_part)
    {
 
      float ptMax = 0;
      trigger::TriggerObject maxPtIter;
      bool hasMatch = false;

      for ( unsigned int hlto = 0; hlto < keys.size(); hlto++ ) {
	size_type hltf = keys[hlto];
	const trigger::TriggerObject L3obj(TOC[hltf]);
      
	float deltaR = reco::deltaR(L3obj.eta(), L3obj.phi(),
				    it_part->second->eta(), it_part->second->phi());
	float thisL3Pt = L3obj.pt();

	if (deltaR < _deltaRMax && thisL3Pt > ptMax) {
	  maxPtIter = L3obj;
	  ptMax = thisL3Pt;
	  hasMatch = true;
	}
      } // end of loop over trigger objects

      if(hasMatch)
	mymap.insert(std::pair<typename myContainer<T>::myIter,
		     trigger::TriggerObject>(it_part->second,maxPtIter));    


    }// end of loop over photons

  return;
}






#endif
