#include "NtupleClasses.hh"

ClassImp(SY_NT::PhoObj);

#define dummy -9999

namespace SY_NT{

  PhoObj::PhoObj():
    _isElectron(dummy),
    _isConverted(dummy),
    _R9(dummy),
    _ECalIso(dummy),
    _ECalIsoRecHit(dummy),
    _HCalIso(dummy),
    _HCalIsoRecHit(dummy),
    _TIso(dummy),
    _SCE(dummy),
    _SCEt(dummy),
    _SCEta(dummy),
    _SCPhi(dummy),
    _SCEtaWidth(dummy),
    _SCPhiWidth(dummy),
    _E(dummy),
    _Et(dummy),
    _Eta(dummy),
    _Phi(dummy),
    _HadEm(dummy),
    _CovEtaEta(dummy)
  {;}

  bool PhoObj::MatchedToEarlyPho()
  {
    if(_GenE <0)return false;
    return true;
  }

  bool PhoObj::IsLoosePhoton()
  {
    
//     if(_TIso > 9.0)return false;
//     if(_ECalIso > 5.0 + 0.015 * _Et)return false;
//     if(_HCalIso > 7.0) return false;
    if(_HadEm > 0.1)return false;
    return true;

  }


  GenObj::GenObj():
    _E(dummy), _Et(dummy), _Eta(dummy), _Phi(dummy)
  {;}

}
