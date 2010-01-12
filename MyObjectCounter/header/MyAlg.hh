#ifndef MYALG_HH
#define MYALG_HH

#include <iostream>
#include <map>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "RecoEgamma/EgammaMCTools/interface/ElectronMCTruth.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruth.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruthFinder.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
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
struct partGenMap
{
  typedef typename std::map<typename myContainer<T>::myIter,
			    reco::GenParticleCollection::const_iterator> Type;
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


// although this is a template class, but should be used only to electrons 
// and photons (reco or pat)
template <class T, class A>
struct elePhoMap
{
  typedef typename std::map< typename myContainer<T>::myIter,
			     typename myContainer<A>::myIter > Type;
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
	    bool doHLT,
	    bool doPAT=false);

  void getHandles(const edm::Event  & event, 
		  bool doPho,
		  bool doEle,
		  bool doHLT,
		  bool doPAT=false);

  void dumpGenInfo(const edm::Event&); 

  // turn on HLT bit information
  void turnOnHLTBit(std::string trgPath, int trgCode);
  int  getThisEventTriggerBit(){return _event_trigger;}

  // print out the private data members
  void print();

  bool isMC(){return (!_isData);}
  bool isData(){return _isData;}
  double ptHat(){return _ptHat;}

  // sort generator particles by Pt (with a PDG and pt > 2 GeV cut), 
  // then insert a pair of et value and const_iterator of GenParticles into a 
  // map (also a private member, genEtmap)

  void       sortGenParticles(); 

  // sort any class of particles by Et, then insert a pair of et value and
  // const_iterator of class T into a map (etmap)
  template<class T> void sortParticles(
	       const edm::Handle<typename myContainer<T>::container >& handle,
	       typename partEtMap<T>::Type& etmap);

  
  // match any class of particles (et sorted, etmap) to a generator-level 
  // particls as long as deltaR, relPt, status criteria are met, then insert a 
  // pair of const_iterators of class T and generator-level particles (genmap)

  template <class T> void matchPartToGen(float deltaRCut,
					 float relPtCut,
					 typename partEtMap<T>::Type  &etmap,
					 typename partGenMap<T>::Type &genmap,
					 int status=1
					 );

  // match any class of particles (et sorted, etmap) to a highest Et L1 object
  // as long as deltaR criteria are met, then insert 
  // pair of const_iterators of class T and L1 object (l1map)

  template <class T> void matchPartToL1(typename partEtMap<T>::Type& etmap, 
					typename partL1Map<T>::Type& l1map );


  // match any class of particles (et sorted, etmap) to a highest Et L3 object
  // as long as deltaR criteria are met, then insert 
  // pair of const_iterators of class T and L3 object (l3map)

  template <class T> void matchPartToL3(std::string trgPath, std::string tag,
					typename partEtMap<T>::Type& etmap,
					typename partL3Map<T>::Type& l3map);  

  // this mapping can be used for either from electron to photon, 
  // or from photon to electron, the first one is the mother
 
  template <class T, class A> void matchEleToPho(
		typename partEtMap<T>::Type& elemap,
		typename partEtMap<A>::Type& phomap,
		typename elePhoMap<T,A>::Type &mymap);


  // although this is a template, should only apply on photons (reco or pat)
  template <class T> bool isMyLoosePhoton(
			typename myContainer<T>::myIter it_ph);
  


  // get the presorted map of photons, electrons, muons, generator-level 
  // particles

  std::vector<reco::GenParticleCollection::const_iterator>     getHardGenVec(){return _hardGenParticle;}
  std::vector<reco::GenParticleCollection::const_iterator>     getConvPhoton(){return _myConvPhotons;}
  std::vector<PhotonMCTruth>         getPhotonMCTruth(){return _myPhotonMCTruth;}

  bool getMatchGen(reco::GenParticleCollection::const_iterator& tempIter,
		   int pdgCode, int status, float px1, float py1, float pz1, 
		   float vx1=-99999, float vy1=-99999, float vz1=-99999);

  partEtMap<reco::Photon>::Type      getPhoEtMap(){return _phoEtMap;}
  partEtMap<reco::GsfElectron>::Type getEleEtMap(){return _eleEtMap;}
  partEtMap<reco::GenParticle>::Type getGenEtMap(){return _genEtMap;}

  partEtMap<pat::Photon>::Type       getPatPhoEtMap(){return _patPhoEtMap;}
  partEtMap<pat::Electron>::Type     getPatEleEtMap(){return _patEleEtMap;}
  partEtMap<pat::Muon>::Type         getPatMuoEtMap(){return _patMuoEtMap;}

  // get each individual handle

  edm::Handle<reco::GenParticleCollection> getGenHandle(){return _genHandle;}
  edm::Handle<edm::TriggerResults> getTrgResultsHandle(){return _trgResultsHandle;}
  edm::Handle<reco::VertexCollection> getVtxHandle(){return _vertexHandle;}

private: 

  // handles to do sorting and mapping for reconstructed particles
  edm::Handle<reco::PhotonCollection>          _phoHandle;
  edm::Handle<reco::GsfElectronCollection>     _eleHandle;
  edm::Handle<reco::GenParticleCollection>     _genHandle;
  edm::Handle<HepMCProduct>                    _hepMCHandle;
  edm::Handle<GenEventInfoProduct>             _genEventHandle;
  edm::Handle<l1extra::L1EmParticleCollection> _l1EmNonIsoHandle;
  edm::Handle<l1extra::L1EmParticleCollection> _l1EmIsoHandle;
  edm::Handle<trigger::TriggerEvent>           _trgEventHandle;
  edm::Handle<edm::TriggerResults>             _trgResultsHandle;
  edm::Handle<reco::VertexCollection>          _vertexHandle;


  std::vector<reco::GenParticleCollection::const_iterator> _hardGenParticle;
  std::vector<reco::GenParticleCollection::const_iterator> _myConvPhotons;
  std::vector<PhotonMCTruth>                      _myPhotonMCTruth;
  PhotonMCTruthFinder*                            _thePhotonMCTruthFinder;

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
  // HLT bit of an event
  int _event_trigger;

  bool _isData;
  double _ptHat;

  // should dump HEPG information or not
  bool _dumpHEP;
  // pdgcode for matchPartToGen and sortGenParticles
  int  _pdgCode;
  // deltaR for matching between particles and L1/L3 trigger objects
  double _trigDeltaRMax;


};

// === inline definition of template functions (has to be in the header files)
//
// sort any class of particles by Et, then insert a pair of et value and
// const_iterator of class T into a map (sortedEtMap)

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


// match any class of particles (et sorted, etmap) to a generator-level 
// particls as long as deltaR, relPt, status criteria are met, then insert a 
// pair of const_iterators of class T and generator-level particles (genmap)
//_____________________________________________________________________________
//
template <class T> void MyAlg::matchPartToGen(
	     float deltaRCut, float relPtCut,
	     typename partEtMap<T>::Type& etmap,
	     typename partGenMap<T>::Type& genmap, int status)
{
  genmap.clear();
  if(!_genHandle.isValid())return;

  typedef typename partEtMap<T>::Type::iterator mapIter;
  for (mapIter it_part= etmap.begin();
       it_part != etmap.end(); ++it_part)
    {

      bool hasMatch = false;

      for( GenParticleCollection::const_iterator it_gen = 
	   _genHandle->begin(); it_gen != _genHandle->end(); it_gen++ ) {
  
	if(abs(it_gen->pdgId())!=_pdgCode)continue;
	if(it_gen->status()!=status)continue;
    
	float dR = reco::deltaR(it_part->second->momentum(), 
				it_gen->momentum());
	float relPt = fabs(it_part->second->pt()-it_gen->pt())/it_gen->pt();

	if(dR< deltaRCut && relPt < relPtCut)
	  {
	    hasMatch = true;
	    genmap.insert(std::pair< typename myContainer<T>::myIter,
			  reco::GenParticleCollection::const_iterator>
			 (it_part->second,it_gen));

	  } // if find a match

	if(hasMatch)break;
    
      } // end of loop of generator-level
    } // end of loop of class T loop

}


// match any class of particles (et sorted, etmap) to a highest Et L1 object
// as long as deltaR criteria are met, then insert 
// pair of const_iterators of class T and L1 object (l1map)

template <class T> void MyAlg::matchPartToL1(
				 typename partEtMap<T>::Type& etmap,
				 typename partL1Map<T>::Type& l1map)
{
  l1map.clear();

  if(!_l1EmIsoHandle.isValid() 
     && !_l1EmNonIsoHandle.isValid())return;

  typedef typename partEtMap<T>::Type::iterator mapIter;
  for(mapIter it_part=etmap.begin(); 
      it_part!= etmap.end(); ++it_part)
    {
      float ptMax = 0;
      std::vector<l1extra::L1EmParticle>::const_iterator maxL1PtIter;
      bool hasMatch = false;

      for( std::vector<l1extra::L1EmParticle>::const_iterator it_l1 = 
	   _l1EmIsoHandle->begin(); 
	   it_l1 != _l1EmIsoHandle->end(); it_l1++ ) {
    
	float deltaR = reco::deltaR(it_l1->eta(), it_l1->phi(),
				    it_part->second->eta(), 
				    it_part->second->phi());
	
	float thisL1Pt = it_l1->pt();
	
	if (deltaR < _trigDeltaRMax && thisL1Pt > ptMax ) {
	  
	  ptMax = thisL1Pt;
	  maxL1PtIter = it_l1;
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

      if (deltaR < _trigDeltaRMax && thisL1Pt > ptMax ) {
	
	ptMax = thisL1Pt;
	maxL1PtIter = it_l1;
	hasMatch = true;
      }

    } // end of L1 non-isolated trigger objects

    
    if(hasMatch)
      l1map.insert(std::pair<typename myContainer<T>::myIter,
		   l1extra::L1EmParticleCollection::const_iterator>(it_part->second,maxL1PtIter));    
  } // end of loop over photons

  return;




}


// match any class of particles (et sorted, etmap) to a highest Et L3 object
// as long as deltaR criteria are met, then insert 
// pair of const_iterators of class T and L3 object (l3map)

//_____________________________________________________________________________
//
template <class T> void MyAlg::matchPartToL3(std::string trgPath,
					     std::string tag,
					    typename partEtMap<T>::Type& etmap,
					    typename partL3Map<T>::Type& l3map)

{
  l3map.clear();

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
      trigger::TriggerObject maxL3PtIter;
      bool hasMatch = false;

      for ( unsigned int hlto = 0; hlto < keys.size(); hlto++ ) {
	size_type hltf = keys[hlto];
	const trigger::TriggerObject L3obj(TOC[hltf]);
      
	float deltaR = reco::deltaR(L3obj.eta(), L3obj.phi(),
				    it_part->second->eta(), it_part->second->phi());
	float thisL3Pt = L3obj.pt();

	if (deltaR < _trigDeltaRMax && thisL3Pt > ptMax) {
	  maxL3PtIter = L3obj;
	  ptMax = thisL3Pt;
	  hasMatch = true;
	}
      } // end of loop over trigger objects

      if(hasMatch)
	l3map.insert(std::pair<typename myContainer<T>::myIter,
		     trigger::TriggerObject>(it_part->second,maxL3PtIter));    


    }// end of loop over photons

  return;
}


//_______________________________________________________________
// this mapping can be used for either from electron to photon, 
// or from photon to electron, the first one is mother
//______________________________________________________________

template <class T, class A> void MyAlg::matchEleToPho(
	      typename partEtMap<T>::Type& elemap,
	      typename partEtMap<A>::Type& phomap,
	      typename elePhoMap<T,A>::Type& mymap)
{
  
 // first loop over electrons
  unsigned int count=0;
  typedef typename partEtMap<T>::Type::iterator mapIter_ab;
  typedef typename partEtMap<A>::Type::iterator mapIter_ey;
 
  for(mapIter_ab eIter = elemap.begin(); eIter!=elemap.end(); ++eIter){

    typename myContainer<T>::myIter it_e = eIter->second;
    if(it_e->pt()>10)count ++ ;
    double ptMax   = 0;
    bool findMatch = false;
    typename myContainer<A>::myIter bestMatchIter;
    // now loop over photons
    for(mapIter_ey pIter = phomap.begin(); pIter!=phomap.end(); ++pIter){
      
      typename myContainer<A>::myIter it_ph = pIter->second;
      // find the highest pt photon that has the same supercluster seed as 
      // electron
      double thisPhoEnergy = it_ph->superCluster()->energy();
      double distx = fabs(it_e->superCluster()->seed()->x() -
			  it_ph->superCluster()->seed()->x());
      double disty = fabs(it_e->superCluster()->seed()->y() -
			  it_ph->superCluster()->seed()->y());
//       double distrphi = sqrt(distx*distx + disty*disty);

      double distz = fabs(it_e->superCluster()->seed()->z() -
			  it_ph->superCluster()->seed()->z());
   
      double relE  = it_e->superCluster()->seed()->energy()>1e-3? 
	fabs(it_ph->superCluster()->seed()->energy() - 
	     it_e->superCluster()->seed()->energy())
	/it_e->superCluster()->seed()->energy(): -999;
      
      if( distx < 3.0 && disty < 3.0 && distz < 3.0 && relE < 0.2	 
	 && thisPhoEnergy > ptMax)
	{
	  ptMax   = thisPhoEnergy;
	  bestMatchIter = it_ph;
	  findMatch = true;
	} // if matched
      
    } // end of loop over photon

    if(findMatch)
      mymap.insert(std::pair<typename myContainer<T>::myIter,
		   typename myContainer<A>::myIter>(it_e,bestMatchIter));
    
  } // end of loop over electron

  if(mymap.size()!=count)
    cout << "There are " << count  << " electrons" << endl;
  

}



//------------------------------------------------------------------------
// although this is a template, should only apply on photons (reco or pat)
//_________________________________________________________________________

template <class T> bool MyAlg::isMyLoosePhoton(
typename myContainer<T>::myIter it_ph)
{
  float et = it_ph->et();
  float eta = it_ph->eta();

  float ecalIso = it_ph->ecalRecHitSumEtConeDR04();
  float hcalIso = it_ph->hcalTowerSumEtConeDR04();
  float trkIso  = it_ph->trkSumPtHollowConeDR04();

  bool isBarrel = it_ph->isEB();
  bool isEndCap = it_ph->isEE();
  bool inAnyGap = it_ph->isEBEEGap() || 
    (it_ph->isEB() && it_ph->isEBGap()) || 
    (it_ph->isEE() && it_ph->isEEGap());

  if(inAnyGap)return false;

  /*
  // Barrel cuts
  if(isBarrel && ecalIso  > 5.0 + 0.0045*et)return false;
  if(isBarrel && hcalIso  > 5.0)return false;
  if(isBarrel && trkIso > 9.0 )return false;
  if(isBarrel && it_ph->hadronicOverEm() > 0.15)return false;

  // Endcap cuts
  if(isEndCap && ecalIso  > 5.0 + 0.02*et)return false;
  if(isEndCap && hcalIso  > 7.0)return false;
  if(isEndCap && trkIso > 9.0 )return false;
  if(isEndCap && it_ph->hadronicOverEm() > 0.15)return false;
  */

  // Barrel cuts
  if(isBarrel && ecalIso  > 5.0 + 0.004*et)return false;
  if(isBarrel && hcalIso  > 5.0)return false;
  if(isBarrel && trkIso > 9.0 )return false;
  if(isBarrel && it_ph->hadronicOverEm() > 0.15)return false;

  // Endcap cuts
  if(isEndCap && ecalIso  > 5.0 + 0.0021*et)return false;
  if(isEndCap && hcalIso  > 5.0)return false;
  if(isEndCap && trkIso > 9.0 )return false;
  if(isEndCap && it_ph->hadronicOverEm() > 0.15)return false;


  if(et        < 10.0)return false;
  if(fabs(eta) > 1.44 && fabs(eta)<1.56)return false;
  if(fabs(eta) > 2.5)return false;

  return true;

}

#endif
