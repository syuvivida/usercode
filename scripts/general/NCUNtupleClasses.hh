#ifndef NCUNTUPLECLASSES_HH
#define NCUNTUPLECLASSES_HH

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
#include <iostream>

#define NTrigPath 102

namespace SY_NT{

  // Event Objects
  
  class EvtObj: public TObject {
  public:
    EvtObj();
    ~EvtObj() {;}

    // Setters 
    void SetRun(int value)                    {_Run = value;}
    void SetEvt(int value)                    {_Evt = value;}
    void SetMET(double value)                 {_MET = value;}
    void SetMETPhi(double value)              {_METPhi = value;}
    void SetHLT(int index, int value)             
    {if(index<0 || index> NTrigPath-1)
	std::cout << "Index is out of range 0~101!" << std::endl;
      else _HLT[index] = value;}

    // Getters
    
    int Run()                                 {return _Run;}
    int Evt()                                 {return _Evt;}
    double MET()                              {return _MET;}
    double METPhi()                           {return _METPhi;}
    int* HLT()                                {return _HLT;}

  private:
    int   _Run;
    int   _Evt;
    double _MET;
    double _METPhi;
    int   _HLT[NTrigPath];    

   ClassDef(EvtObj,1)

  };

  
  // photon objects
  class PhoObj: public TObject{
  public:
    PhoObj();
    ~PhoObj() {;}


    // Setters
    void  SetE(double value)                  {_E   = value;}
    void  SetEt(double value)                 {_Et  = value;}
    void  SetEta(double value)                {_Eta = value;}
    void  SetPhi(double value)                {_Phi = value;}
    void  SetR9(double value)                 {_R9 = value;}
    void  SetTIso(double value)               {_TIso = value;}
    void  SetECalIso(double value)            {_ECalIso = value;}
    void  SetHCalIso(double value)            {_HCalIso = value;}
    void  SetHadEm(double value)              {_HadEm = value;}
    void  SetSigEtaEta(double value)          {_SigEtaEta = value;}
    void  SetSigIEtaIEta(double value)        {_SigIEtaIEta = value;}
    void  SetLocation(int value)              {_Location = value;}

    void  SetGenIndex(int value)              {_GenIndex = value;}
    void  SetGenMomPID(int value)             {_GenMomPID = value;}
    void  SetGenMomPt(double value)           {_GenMomPt = value;}
    void  SetGenGMomPID(int value)            {_GenGMomPID = value;}

    void  SetSCE(double value)                {_SCE = value;}
    void  SetSCEt(double value)               {_SCEt  = value;}
    void  SetSCEta(double value)              {_SCEta = value;}
    void  SetSCPhi(double value)              {_SCPhi = value;}
    void  SetSCEtaWidth(double value)         {_SCEtaWidth = value;}
    void  SetSCPhiWidth(double value)         {_SCPhiWidth = value;}
    
 
    // Getters
    
    double  E()                               {return _E;}
    double  Et()                              {return _Et;}
    double  Eta()                             {return _Eta;}
    double  Phi()                             {return _Phi;}
    double  R9()                              {return _R9;}
    double  TIso()                            {return _TIso;}
    double  ECalIso()                         {return _ECalIso;}
    double  HCalIso()                         {return _HCalIso;}
    double  CalIso()                          {return _ECalIso + _HCalIso;}
    double  HadEm()                           {return _HadEm;}
    double  SigEtaEta()                       {return _SigEtaEta;}
    double  SigIEtaIEta()                     {return _SigIEtaIEta;}
    int     Location()                        {return _Location;}

    int     GenIndex()                        {return _GenIndex;}
    int     GenMomPID()                       {return _GenMomPID;}
    double  GenMomPt()                        {return _GenMomPt;}
    int     GenGMomPID()                      {return _GenGMomPID;}

    double  SCE()                             {return _SCE;}
    double  SCEt()                            {return _SCEt;}
    double  SCEta()                           {return _SCEta;}
    double  SCPhi()                           {return _SCPhi;}
    double  SCEtaWidth()                      {return _SCEtaWidth;}
    double  SCPhiWidth()                      {return _SCPhiWidth;}
    
    bool    IsLoosePhoton();
    
    
  private:



    double   _E;
    double   _Et;
    double   _Eta;
    double   _Phi;
    double   _R9;
    double   _TIso;
    double   _ECalIso;
    double   _HCalIso;
    double   _HadEm;
    double   _SigEtaEta;
    double   _SigIEtaIEta;
    int      _Location;

    int      _GenIndex;
    int      _GenMomPID;
    double   _GenMomPt;
    int      _GenGMomPID;

    double   _SCE;
    double   _SCEt;
    double   _SCEta;
    double   _SCPhi;
    double   _SCEtaWidth;
    double   _SCPhiWidth;


    ClassDef(PhoObj,1)
 };
  

  // generator-level objects

  class GenObj: public TObject{
  public:
    GenObj();
    ~GenObj() {;}


    // Setters
    void SetIndex(int value)                  {_Index = value;}
    void SetPID(int value)                    {_PID = value;}
    void SetPt(double value)                  {_Pt = value;}
    void SetEta(double value)                 {_Eta = value;}
    void SetPhi(double value)                 {_Phi = value;}
    void SetMass(double value)                {_Mass = value;}

    void SetMomPID(int value)                 {_MomPID = value;}
    void SetMomPt(double value)               {_MomPt = value;}
    void SetMomEta(double value)              {_MomEta = value;}
    void SetMomPhi(double value)              {_MomPhi = value;}
    void SetMomMass(double value)             {_MomMass = value;}
    void SetGMomPID(int value)                {_GMomPID = value;}


    // Getters
    int Index()                               {return _Index;}
    int PID()                                 {return _PID;}
    double Pt()                               {return _Pt;}
    double Eta()                              {return _Eta;}
    double Phi()                              {return _Phi;}
    double Mass()                             {return _Mass;}

    int MomPID()                              {return _MomPID;}
    double MomPt()                            {return _MomPt;}
    double MomEta()                           {return _MomEta;}
    double MomPhi()                           {return _MomPhi;}
    double MomMass()                          {return _MomMass;}
    int GMomPID()                             {return _GMomPID;}

  private:

    int    _Index;
    int    _PID;
    double _Pt;
    double _Eta;
    double _Phi;
    double _Mass;

    int    _MomPID;
    double _MomPt;
    double _MomEta;
    double _MomPhi;
    double _MomMass;
    int    _GMomPID;

    ClassDef(GenObj,1)
 };
   

} // end of SY_NT name space

// load ntuples

void LoadEvt(TChain* ch, SY_NT::EvtObj & thisEvt)
{
  thisEvt.SetRun((int)ch->GetBranch("run")->GetLeaf("run")->GetValue());
  thisEvt.SetEvt((int)ch->GetBranch("event")->GetLeaf("event")->GetValue());
  thisEvt.SetMET(ch->GetBranch("MET")->GetLeaf("MET")->GetValue());
  thisEvt.SetMETPhi(ch->GetBranch("METPhi")->GetLeaf("METPhi")->GetValue());
  for(int i=0; i<NTrigPath; i++)
    thisEvt.SetHLT(i,
		   (int)ch->GetBranch("HLT")->GetLeaf("HLT")->GetValue(i));

}

void LoadPhos(TChain *ch, std::vector<SY_NT::PhoObj> & PhoVect)
{
  int npho=(int)ch->GetBranch("nPho")->GetLeaf("nPho")->GetValue();


  for(int i=0; i < npho; i++)
    {
      SY_NT::PhoObj thisObj;

      thisObj.SetE(ch->GetBranch("phoE")->GetLeaf("phoE")->GetValue(i));
      thisObj.SetEt(ch->GetBranch("phoEt")->GetLeaf("phoEt")->GetValue(i));
      thisObj.SetEta(ch->GetBranch("phoEta")->GetLeaf("phoEta")->GetValue(i));
      thisObj.SetPhi(ch->GetBranch("phoPhi")->GetLeaf("phoPhi")->GetValue(i));
      thisObj.SetR9(ch->GetBranch("phoR9")->GetLeaf("phoR9")->GetValue(i));
      thisObj.SetTIso(ch->GetBranch("phoTrkIso")->GetLeaf("phoTrkIso")->GetValue(i));
      thisObj.SetECalIso(ch->GetBranch("phoEcalIso")->GetLeaf("phoEcalIso")->GetValue(i));
      thisObj.SetHCalIso(ch->GetBranch("phoHcalIso")->GetLeaf("phoHcalIso")->GetValue(i));
      thisObj.SetHadEm(ch->GetBranch("phoHoverE")->GetLeaf("phoHoverE")->GetValue(i));
      thisObj.SetSigEtaEta(ch->GetBranch("phoSigmaEtaEta")->GetLeaf("phoSigmaEtaEta")->GetValue(i));
      thisObj.SetSigIEtaIEta(ch->GetBranch("phoSigmaIEtaIEta")->GetLeaf("phoSigmaIEtaIEta")->GetValue(i));
      thisObj.SetLocation((int)ch->GetBranch("phoPos")->GetLeaf("phoPos")->GetValue(i));
      thisObj.SetGenIndex((int)ch->GetBranch("phoGenIndex")->GetLeaf("phoGenIndex")->GetValue(i));
      thisObj.SetGenMomPID((int)ch->GetBranch("phoGenMomPID")->GetLeaf("phoGenMomPID")->GetValue(i));
      thisObj.SetGenMomPt(ch->GetBranch("phoGenMomPt")->GetLeaf("phoGenMomPt")->GetValue(i));
      thisObj.SetGenGMomPID((int)ch->GetBranch("phoGenGMomPID")->GetLeaf("phoGenGMomPID")->GetValue(i));
      thisObj.SetSCE(ch->GetBranch("phoSCE")->GetLeaf("phoSCE")->GetValue(i));
      thisObj.SetSCEt(ch->GetBranch("phoSCEt")->GetLeaf("phoSCEt")->GetValue(i));
      thisObj.SetSCEta(ch->GetBranch("phoSCEta")->GetLeaf("phoSCEta")->GetValue(i));
      thisObj.SetSCPhi(ch->GetBranch("phoSCPhi")->GetLeaf("phoSCPhi")->GetValue(i));
      thisObj.SetSCEtaWidth(ch->GetBranch("phoSCEtaWidth")->GetLeaf("phoSCEtaWidth")->GetValue(i));
      thisObj.SetSCPhiWidth(ch->GetBranch("phoSCPhiWidth")->GetLeaf("phoSCPhiWidth")->GetValue(i));
      
      PhoVect.push_back(thisObj);
      
    } // loop over nphotons

}  // end of LoadPhos


void LoadGens(TChain *ch, std::vector<SY_NT::GenObj> & GenVect)
{
  int ngen=(int)ch->GetBranch("nMC")->GetLeaf("nMC")->GetValue();
  
  for(int i=0; i < ngen; i++)
    {
      SY_NT::GenObj thisObj;
      
      thisObj.SetIndex((int)ch->GetBranch("mcIndex")->GetLeaf("mcIndex")->GetValue(i));
      thisObj.SetPID((int)ch->GetBranch("mcPID")->GetLeaf("mcPID")->GetValue(i));
      thisObj.SetPt(ch->GetBranch("mcPt")->GetLeaf("mcPt")->GetValue(i));
      thisObj.SetEta(ch->GetBranch("mcEta")->GetLeaf("mcEta")->GetValue(i));
      thisObj.SetPhi(ch->GetBranch("mcPhi")->GetLeaf("mcPhi")->GetValue(i));
      thisObj.SetMass(ch->GetBranch("mcMass")->GetLeaf("mcMass")->GetValue(i));
      thisObj.SetMomPID((int)ch->GetBranch("mcMomPID")->GetLeaf("mcMomPID")->GetValue(i));
      thisObj.SetMomPt(ch->GetBranch("mcMomPt")->GetLeaf("mcMomPt")->GetValue(i));
      thisObj.SetMomEta(ch->GetBranch("mcMomEta")->GetLeaf("mcMomEta")->GetValue(i));
      thisObj.SetMomPhi(ch->GetBranch("mcMomPhi")->GetLeaf("mcMomPhi")->GetValue(i));
      thisObj.SetMomMass(ch->GetBranch("mcMomMass")->GetLeaf("mcMomMass")->GetValue(i));
      thisObj.SetGMomPID((int)ch->GetBranch("mcGMomPID")->GetLeaf("mcGMomPID")->GetValue(i)); 

      GenVect.push_back(thisObj);


    }

}

#endif
