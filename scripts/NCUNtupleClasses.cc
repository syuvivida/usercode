#include "NCUNtupleClasses.hh"

ClassImp(SY_NT::EvtObj);
ClassImp(SY_NT::PhoObj);
ClassImp(SY_NT::GenObj);

#define dummy -9999

namespace SY_NT{

  EvtObj::EvtObj():
    _Run(dummy),
    _Evt(dummy),
    _MET(dummy),
    _METPhi(dummy)
  { memset(_HLT,0x00,sizeof(_HLT));}

  PhoObj::PhoObj():
    _E(dummy),
    _Et(dummy),
    _Eta(dummy),
    _Phi(dummy),
    _R9(dummy),
    _TIso(dummy),
    _ECalIso(dummy),
    _HCalIso(dummy),
    _HadEm(dummy),
    _SigEtaEta(dummy),
    _SigIEtaIEta(dummy),
    _Location(dummy),    
    _GenIndex(dummy),
    _GenMomPID(dummy),
    _GenMomPt(dummy),
    _GenGMomPID(dummy),    
    _SCE(dummy),
    _SCEt(dummy),
    _SCEta(dummy),
    _SCPhi(dummy),
    _SCEtaWidth(dummy),
    _SCPhiWidth(dummy)
  {;}


  bool PhoObj::IsLoosePhoton()
  {
    // default
    if(_TIso > 9.0)return false;
    if(_ECalIso > 5.0 + 0.015 * _Et)return false;
    if(_HCalIso > 7.0) return false;
    if(_R9 < 0.93)return false;
    if(_HadEm > 0.1)return false;
    if(_Et < 10.0)return false;
    if(fabs(_Eta) > 1.44 && fabs(_Eta)<1.56)return false;
    if(fabs(_Eta) > 2.5)return false;
    return true;

  }


  GenObj::GenObj():
    _Index(dummy),
    _PID(dummy),
    _Pt(dummy),
    _Eta(dummy),
    _Phi(dummy),
    _Mass(dummy),
    _MomPID(dummy),
    _MomPt(dummy),
    _MomEta(dummy),
    _MomPhi(dummy),
    _MomMass(dummy),
    _GMomPID(dummy)
  {;}

}
