#ifndef NTUPLECLASSES_HH
#define NTUPLECLASSES_HH

//////////////////////////////////////////////////////
//
// Classes to fill my ntuples with
//
////////////////////////////////////////////////////////

#include "TObject.h"
#include "TChain.h"
#include <TBranch.h>
#include <TLeaf.h>
#include <vector>

namespace SY_NT{
  class PhoObj: public TObject{
  public:
    PhoObj();
    ~PhoObj() {;}


    // Setters
    void  SetIsElectron(bool value)          {_isElectron = value;}
    void  SetIsConverted(bool value)         {_isConverted = value;}
    void  SetR9(float value)                 {_R9 = value;}
    void  SetECalIso(float value)            {_ECalIso = value;}
    void  SetECalIsoRecHit(float value)      {_ECalIsoRecHit = value;}
    void  SetHCalIso(float value)            {_HCalIso = value;}
    void  SetHCalIsoRecHit(float value)      {_HCalIsoRecHit = value;}
    void  SetTIso(float value)               {_TIso = value;}
    void  SetSCE(float value)                {_SCE = value;}
    void  SetSCEt(float value)               {_SCEt  = value;}
    void  SetSCEta(float value)              {_SCEta = value;}
    void  SetSCPhi(float value)              {_SCPhi = value;}
    void  SetSCEtaWidth(float value)         {_SCEtaWidth = value;}
    void  SetSCPhiWidth(float value)         {_SCPhiWidth = value;}
    void  SetE(float value)                  {_E   = value;}
    void  SetEt(float value)                 {_Et  = value;}
    void  SetEta(float value)                {_Eta = value;}
    void  SetPhi(float value)                {_Phi = value;}
    void  SetHadEm(float value)              {_HadEm = value;}
    void  SetCovEtaEta(float value)          {_CovEtaEta = value;}
    
    void  SetGenE(float value)               {_GenE = value;}
    void  SetGenEt(float value)              {_GenEt = value;}
    void  SetGenEta(float value)             {_GenEta = value;}
    void  SetGenPhi(float value)             {_GenPhi = value;}
 
    // Getters
    bool  IsElectron()                       {return _isElectron;}
    bool  IsConverted()                      {return _isConverted;}
    float R9()                               {return _R9;}
    float ECalIso()                          {return _ECalIso;}
    float ECalIsoRecHit()                    {return _ECalIsoRecHit;}
    float HCalIso()                          {return _HCalIso;}
    float HCalIsoRecHit()                    {return _HCalIsoRecHit;}
    float CalIso()                           {return _ECalIso + _HCalIso;}
    float TIso()                             {return _TIso;}
    float SCE()                              {return _SCE;}
    float SCEt()                             {return _SCEt;}
    float SCEta()                            {return _SCEta;}
    float SCPhi()                            {return _SCPhi;}
    float SCEtaWidth()                       {return _SCEtaWidth;}
    float SCPhiWidth()                       {return _SCPhiWidth;}
    float E()                                {return _E;}
    float Et()                               {return _Et;}
    float Eta()                              {return _Eta;}
    float Phi()                              {return _Phi;}
    float HadEm()                            {return _HadEm;}
    float CovEtaEta()                        {return _CovEtaEta;}

    float GenE()                             {return _GenE;}
    float GenEt()                            {return _GenEt;}
    float GenEta()                           {return _GenEta;}
    float GenPhi()                           {return _GenPhi;}

    bool  MatchedToEarlyPho();    
    bool  IsLoosePhoton();
    
    
  private:
    bool  _isElectron;
    bool  _isConverted;
    float _R9;
    float _ECalIso;
    float _ECalIsoRecHit;
    float _HCalIso;
    float _HCalIsoRecHit;
    float _TIso;
    float _SCE;
    float _SCEt;
    float _SCEta;
    float _SCPhi;
    float _SCEtaWidth;
    float _SCPhiWidth;
    float _E;
    float _Et;
    float _Eta;
    float _Phi;
    float _HadEm;
    float _CovEtaEta;

    float _GenE;
    float _GenEt;
    float _GenEta;
    float _GenPhi;


    ClassDef(PhoObj,1)
 };
  

  class GenObj: public TObject{
  public:
    GenObj();
    ~GenObj() {;}

    // Setters
    void SetE(float value)                   {_E   = value;}
    void SetEt(float value)                  {_Et  = value;}
    void SetEta(float value)                 {_Eta = value;}
    void SetPhi(float value)                 {_Phi = value;}

    // Getters
    float E()                                { return _E;}
    float Et()                               { return _Et;}
    float Eta()                              { return _Eta;}
    float Phi()                              { return _Phi;}

  private:
    
    float _E;
    float _Et;
    float _Eta;
    float _Phi;

    ClassDef(GenObj,1)
 };
   

} // end of SY_NT name space

// load ntuples

void LoadPhos(TChain *ch, std::vector<SY_NT::PhoObj> & PhoVect)
{
  int npho=(int)ch->GetBranch("Photonindex")->GetLeaf("Photonindex")->GetValue();

  for(int i=0; i < npho; i++)
    {
      SY_NT::PhoObj thisObj;
      thisObj.SetIsElectron((int)ch->GetBranch("isAlsoElectron")->GetLeaf("isAlsoElectron")->GetValue(i));
      thisObj.SetIsConverted((int)ch->GetBranch("isConverted")->GetLeaf("isConverted")->GetValue(i));
      thisObj.SetR9(ch->GetBranch("r9")->GetLeaf("r9")->GetValue(i));
      thisObj.SetECalIso(ch->GetBranch("ecalIso")->GetLeaf("ecalIso")->GetValue(i));
      thisObj.SetECalIsoRecHit(ch->GetBranch("isolationEcalRecHit")->GetLeaf("isolationEcalRecHit")->GetValue(i));
      thisObj.SetHCalIso(ch->GetBranch("hcalIso")->GetLeaf("hcalIso")->GetValue(i));
      thisObj.SetHCalIsoRecHit(ch->GetBranch("isolationHcalRecHit")->GetLeaf("isolationHcalRecHit")->GetValue(i));
      thisObj.SetTIso(ch->GetBranch("trackIso")->GetLeaf("trackIso")->GetValue(i));
      thisObj.SetSCE(ch->GetBranch("scE")->GetLeaf("scE")->GetValue(i));
      thisObj.SetSCEt(ch->GetBranch("scEt")->GetLeaf("scEt")->GetValue(i));
      thisObj.SetSCEta(ch->GetBranch("scEta")->GetLeaf("scEta")->GetValue(i));
      thisObj.SetSCPhi(ch->GetBranch("scPhi")->GetLeaf("scPhi")->GetValue(i));
      thisObj.SetSCEtaWidth(ch->GetBranch("scEtaWidth")->GetLeaf("scEtaWidth")->GetValue(i));
      thisObj.SetSCPhiWidth(ch->GetBranch("scPhiWidth")->GetLeaf("scPhiWidth")->GetValue(i));
      thisObj.SetE(ch->GetBranch("PhotonE")->GetLeaf("PhotonE")->GetValue(i));
      thisObj.SetEt(ch->GetBranch("PhotonEt")->GetLeaf("PhotonEt")->GetValue(i));
      thisObj.SetEta(ch->GetBranch("PhotonEta")->GetLeaf("PhotonEta")->GetValue(i));
      thisObj.SetPhi(ch->GetBranch("PhotonPhi")->GetLeaf("PhotonPhi")->GetValue(i));
      thisObj.SetHadEm(ch->GetBranch("HoverE")->GetLeaf("HoverE")->GetValue(i));
      thisObj.SetCovEtaEta(ch->GetBranch("CovEtaEta")->GetLeaf("CovEtaEta")->GetValue(i));
      
      thisObj.SetGenE(ch->GetBranch("MatchedGenPE")->GetLeaf("MatchedGenPE")->GetValue(i));
      thisObj.SetGenEt(ch->GetBranch("MatchedGenPEt")->GetLeaf("MatchedGenPEt")->GetValue(i));      
      thisObj.SetGenEta(ch->GetBranch("MatchedGenPEta")->GetLeaf("MatchedGenPEta")->GetValue(i));
      thisObj.SetGenPhi(ch->GetBranch("MatchedGenPPhi")->GetLeaf("MatchedGenPPhi")->GetValue(i));


      PhoVect.push_back(thisObj);

    } // loop over nphotons

}  // end of LoadPhos


void LoadGens(TChain *ch, std::vector<SY_NT::GenObj> & GenVect)
{
  int ngen=(int)ch->GetBranch("GenPhotonindex")->GetLeaf("GenPhotonindex")->GetValue();
  
  for(int i=0; i < ngen; i++)
    {
      SY_NT::GenObj thisObj;
      thisObj.SetE(ch->GetBranch("GenPhotonE")->GetLeaf("GenPhotonE")->GetValue(i));
      thisObj.SetEt(ch->GetBranch("GenPhotonEt")->GetLeaf("GenPhotonEt")->GetValue(i));
      thisObj.SetEta(ch->GetBranch("GenPhotonEta")->GetLeaf("GenPhotonEta")->GetValue(i));
      thisObj.SetPhi(ch->GetBranch("GenPhotonPhi")->GetLeaf("GenPhotonPhi")->GetValue(i));
      
      GenVect.push_back(thisObj);


    }

}

#endif
