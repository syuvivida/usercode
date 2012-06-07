//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jun  4 12:01:45 2012 by ROOT version 5.27/06b
// from TTree tree/tree
// found on file: /data2/syu/zjet_vectorNtuple/anne-marie_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola.root
//////////////////////////////////////////////////////////

#ifndef puweight_sys_h
#define puweight_sys_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string>
#include <iostream>

using namespace std;



   const Int_t kMaxEvtInfo_VertexX = 1;
   const Int_t kMaxEvtInfo_VertexY = 1;
   const Int_t kMaxEvtInfo_VertexZ = 1;
   const Int_t kMaxptHat = 1;
   const Int_t kMaxmcWeight = 1;
   const Int_t kMaxgenParE = 1;
   const Int_t kMaxgenParPt = 1;
   const Int_t kMaxgenParEta = 1;
   const Int_t kMaxgenParPhi = 1;
   const Int_t kMaxgenParQ = 1;
   const Int_t kMaxgenParId = 1;
   const Int_t kMaxgenParSt = 1;
   const Int_t kMaxgenMomParId = 1;
   const Int_t kMaxgenParIndex = 1;
   const Int_t kMaxgenJetE = 1;
   const Int_t kMaxgenJetPt = 1;
   const Int_t kMaxgenJetEta = 1;
   const Int_t kMaxgenJetPhi = 1;
   const Int_t kMaxfastJetRho = 1;
   const Int_t kMaxPFchsJetRho = 1;
   const Int_t kMaxlepIsoRho = 1;
   const Int_t kMaxpatElecEt = 1;
   const Int_t kMaxpatElecEnergy = 1;
   const Int_t kMaxpatElecPt = 1;
   const Int_t kMaxpatElecEta = 1;
   const Int_t kMaxpatElecPhi = 1;
   const Int_t kMaxpatElecM = 1;
   const Int_t kMaxpatElecScEn = 1;
   const Int_t kMaxpatElecScEt = 1;
   const Int_t kMaxpatElecScEta = 1;
   const Int_t kMaxpatElecScPhi = 1;
   const Int_t kMaxpatElecisEcalDriven = 1;
   const Int_t kMaxpatElecisTrackerDriven = 1;
   const Int_t kMaxpatElecSigIhIh = 1;
   const Int_t kMaxpatElecDelEtaIn = 1;
   const Int_t kMaxpatElecDelPhiIn = 1;
   const Int_t kMaxpatElecHoE = 1;
   const Int_t kMaxpatElecTrkIso = 1;
   const Int_t kMaxpatElecHcalIso = 1;
   const Int_t kMaxpatElecEcalIso = 1;
   const Int_t kMaxpatElecCharge = 1;
   const Int_t kMaxpatElecRelIsoComb = 1;
   const Int_t kMaxpatElecDxy = 1;
   const Int_t kMaxpatElecD0 = 1;
   const Int_t kMaxpatElecDsz = 1;
   const Int_t kMaxpatElecDz = 1;
   const Int_t kMaxpatElecDxyBS = 1;
   const Int_t kMaxpatElecDszBS = 1;
   const Int_t kMaxpatElecDzBS = 1;
   const Int_t kMaxpatElecMva = 1;
   const Int_t kMaxpatElecChHadSumPt03 = 1;
   const Int_t kMaxpatElecNeHadSumPt03 = 1;
   const Int_t kMaxpatElecGamSumPt03 = 1;
   const Int_t kMaxpatElecChHadSumPt04 = 1;
   const Int_t kMaxpatElecChHadIso = 1;
   const Int_t kMaxpatElecNeHadIso = 1;
   const Int_t kMaxpatElecGamIso = 1;
   const Int_t kMaxpatElecNeHadSumPt04 = 1;
   const Int_t kMaxpatElecGamSumPt04 = 1;
   const Int_t kMaxpatElecChHadSumPt05 = 1;
   const Int_t kMaxpatElecNeHadSumPt05 = 1;
   const Int_t kMaxpatElecGamSumPt05 = 1;
   const Int_t kMaxpatElecMissingHits = 1;
   const Int_t kMaxpatElecDist = 1;
   const Int_t kMaxpatElecDeltaCotTheta = 1;
   const Int_t kMaxpatElecConvRadius = 1;
   const Int_t kMaxpatElecInBarrel = 1;
   const Int_t kMaxpatElecInEndcap = 1;
   const Int_t kMaxpatJetPfAk05Pt = 1;
   const Int_t kMaxpatJetPfAk05Eta = 1;
   const Int_t kMaxpatJetPfAk05Phi = 1;
   const Int_t kMaxpatJetPfAk05M = 1;
   const Int_t kMaxpatJetPfAk05Rapidity = 1;
   const Int_t kMaxpatJetPfAk05Px = 1;
   const Int_t kMaxpatJetPfAk05Py = 1;
   const Int_t kMaxpatJetPfAk05Pz = 1;
   const Int_t kMaxpatJetPfAk05En = 1;
   const Int_t kMaxpatJetPfAk05UnCorrPt = 1;
   const Int_t kMaxpatJetPfAk05UnCorrPx = 1;
   const Int_t kMaxpatJetPfAk05UnCorrPy = 1;
   const Int_t kMaxpatJetPfAk05UnCorrPz = 1;
   const Int_t kMaxpatJetPfAk05UnCorrEnt = 1;
   const Int_t kMaxpatJetPfAk05PhotEn = 1;
   const Int_t kMaxpatJetPfAk05ElecEn = 1;
   const Int_t kMaxpatJetPfAk05MuonEn = 1;
   const Int_t kMaxpatJetPfAk05HfHadEn = 1;
   const Int_t kMaxpatJetPfAk05HfEmEn = 1;
   const Int_t kMaxpatJetPfAk05CharHadE = 1;
   const Int_t kMaxpatJetPfAk05NeutHadE = 1;
   const Int_t kMaxpatJetPfAk05CharEmE = 1;
   const Int_t kMaxpatJetPfAk05CharMuE = 1;
   const Int_t kMaxpatJetPfAk05NeutEmE = 1;
   const Int_t kMaxpatJetPfAk05MuonMulti = 1;
   const Int_t kMaxpatJetPfAk05NeutMulti = 1;
   const Int_t kMaxpatJetPfAk05CharMulti = 1;
   const Int_t kMaxpatJetPfAk05CharHadMulti = 1;
   const Int_t kMaxpatJetPfAk05NeutHadMulti = 1;
   const Int_t kMaxpatJetPfAk05PhotMulti = 1;
   const Int_t kMaxpatJetPfAk05ElecMulti = 1;
   const Int_t kMaxpatJetPfAk05HfHadMulti = 1;
   const Int_t kMaxpatJetPfAk05HfEmMulti = 1;
   const Int_t kMaxpatJetPfAk05PhotEnFr = 1;
   const Int_t kMaxpatJetPfAk05MuonEnFr = 1;
   const Int_t kMaxpatJetPfAk05HfHadEnFr = 1;
   const Int_t kMaxpatJetPfAk05HfEmEnFr = 1;
   const Int_t kMaxpatJetPfAk05NeutEmEFr = 1;
   const Int_t kMaxpatJetPfAk05CharHadEFr = 1;
   const Int_t kMaxpatJetPfAk05NeutHadEFr = 1;
   const Int_t kMaxpatJetPfAk05CharEmEFr = 1;
   const Int_t kMaxpatJetPfAk05CharMuEFr = 1;
   const Int_t kMaxpatJetPfAk05N60 = 1;
   const Int_t kMaxpatJetPfAk05N90 = 1;
   const Int_t kMaxpatJetPfAk05EmFr = 1;
   const Int_t kMaxpatJetPfAk05HadFr = 1;
   const Int_t kMaxpatJetPfAk05EmEbEn = 1;
   const Int_t kMaxpatJetPfAk05EmEeEn = 1;
   const Int_t kMaxpatJetPfAk05EmHfEn = 1;
   const Int_t kMaxpatJetPfAk05HadHbEn = 1;
   const Int_t kMaxpatJetPfAk05HadHeEn = 1;
   const Int_t kMaxpatJetPfAk05HadHfEn = 1;
   const Int_t kMaxpatJetPfAk05HadHoEn = 1;
   const Int_t kMaxpatJetPfAk05TotalEm = 1;
   const Int_t kMaxpatJetPfAk05TotalHad = 1;
   const Int_t kMaxpatJetPfAk05EmEbFr = 1;
   const Int_t kMaxpatJetPfAk05EmEeFr = 1;
   const Int_t kMaxpatJetPfAk05EmHfFr = 1;
   const Int_t kMaxpatJetPfAk05HadHbFr = 1;
   const Int_t kMaxpatJetPfAk05HadHeFr = 1;
   const Int_t kMaxpatJetPfAk05HadHfFr = 1;
   const Int_t kMaxpatJetPfAk05HadHoFr = 1;
   const Int_t kMaxpatJetPfAk05JesUncert = 1;
   const Int_t kMaxpatJetPfAk05GenPartonID = 1;
   const Int_t kMaxpatJetPfAk05GenPartonPt = 1;
   const Int_t kMaxpatJetPfAk05GenPartonEta = 1;
   const Int_t kMaxpatJetPfAk05GenPartonPhi = 1;
   const Int_t kMaxpatJetPfAk05GenPartonE = 1;

class puweight_sys {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        PU_weight;
   Double_t        PU_nTrueInt;
   Int_t           PU_nPUVert;
   Int_t           EvtInfo_EventNum;
   Int_t           EvtInfo_RunNum;
   Int_t           EvtInfo_LumiSection;
   Int_t           EvtInfo_BunchXing;
   Int_t           EvtInfo_NumVtx;
   vector<double>  *EvtInfo_VertexX_;
   vector<double>  *EvtInfo_VertexY_;
   vector<double>  *EvtInfo_VertexZ_;
   Double_t        ptHat_;
   Double_t        mcWeight_;
   vector<double>  *genParE_;
   vector<double>  *genParPt_;
   vector<double>  *genParEta_;
   vector<double>  *genParPhi_;
   vector<int>     *genParQ_;
   vector<int>     *genParId_;
   vector<int>     *genParSt_;
   vector<int>     *genMomParId_;
   vector<int>     *genParIndex_;
   vector<double>  *genJetE_;
   vector<double>  *genJetPt_;
   vector<double>  *genJetEta_;
   vector<double>  *genJetPhi_;
   vector<double>  *fastJetRho_;
   vector<double>  *PFchsJetRho_;
   vector<double>  *lepIsoRho_;
   vector<double>  *patElecEt_;
   vector<double>  *patElecEnergy_;
   vector<double>  *patElecPt_;
   vector<double>  *patElecEta_;
   vector<double>  *patElecPhi_;
   vector<double>  *patElecM_;
   vector<double>  *patElecScEn_;
   vector<double>  *patElecScEt_;
   vector<double>  *patElecScEta_;
   vector<double>  *patElecScPhi_;
   vector<double>  *patElecisEcalDriven_;
   vector<double>  *patElecisTrackerDriven_;
   vector<double>  *patElecSigIhIh_;
   vector<double>  *patElecDelEtaIn_;
   vector<double>  *patElecDelPhiIn_;
   vector<double>  *patElecHoE_;
   vector<double>  *patElecTrkIso_;
   vector<double>  *patElecHcalIso_;
   vector<double>  *patElecEcalIso_;
   vector<double>  *patElecCharge_;
   vector<double>  *patElecRelIsoComb_;
   vector<double>  *patElecDxy_;
   vector<double>  *patElecD0_;
   vector<double>  *patElecDsz_;
   vector<double>  *patElecDz_;
   vector<double>  *patElecDxyBS_;
   vector<double>  *patElecDszBS_;
   vector<double>  *patElecDzBS_;
   vector<double>  *patElecMva_;
   vector<double>  *patElecChHadSumPt03_;
   vector<double>  *patElecNeHadSumPt03_;
   vector<double>  *patElecGamSumPt03_;
   vector<double>  *patElecChHadSumPt04_;
   vector<double>  *patElecChHadIso_;
   vector<double>  *patElecNeHadIso_;
   vector<double>  *patElecGamIso_;
   vector<double>  *patElecNeHadSumPt04_;
   vector<double>  *patElecGamSumPt04_;
   vector<double>  *patElecChHadSumPt05_;
   vector<double>  *patElecNeHadSumPt05_;
   vector<double>  *patElecGamSumPt05_;
   vector<double>  *patElecMissingHits_;
   vector<double>  *patElecDist_;
   vector<double>  *patElecDeltaCotTheta_;
   vector<double>  *patElecConvRadius_;
   vector<double>  *patElecInBarrel_;
   vector<double>  *patElecInEndcap_;
   vector<double>  *patJetPfAk05Pt_;
   vector<double>  *patJetPfAk05Eta_;
   vector<double>  *patJetPfAk05Phi_;
   vector<double>  *patJetPfAk05M_;
   vector<double>  *patJetPfAk05Rapidity_;
   vector<double>  *patJetPfAk05Px_;
   vector<double>  *patJetPfAk05Py_;
   vector<double>  *patJetPfAk05Pz_;
   vector<double>  *patJetPfAk05En_;
   vector<double>  *patJetPfAk05UnCorrPt_;
   vector<double>  *patJetPfAk05UnCorrPx_;
   vector<double>  *patJetPfAk05UnCorrPy_;
   vector<double>  *patJetPfAk05UnCorrPz_;
   vector<double>  *patJetPfAk05UnCorrEnt_;
   vector<double>  *patJetPfAk05PhotEn_;
   vector<double>  *patJetPfAk05ElecEn_;
   vector<double>  *patJetPfAk05MuonEn_;
   vector<double>  *patJetPfAk05HfHadEn_;
   vector<double>  *patJetPfAk05HfEmEn_;
   vector<double>  *patJetPfAk05CharHadE_;
   vector<double>  *patJetPfAk05NeutHadE_;
   vector<double>  *patJetPfAk05CharEmE_;
   vector<double>  *patJetPfAk05CharMuE_;
   vector<double>  *patJetPfAk05NeutEmE_;
   vector<double>  *patJetPfAk05MuonMulti_;
   vector<double>  *patJetPfAk05NeutMulti_;
   vector<double>  *patJetPfAk05CharMulti_;
   vector<double>  *patJetPfAk05CharHadMulti_;
   vector<double>  *patJetPfAk05NeutHadMulti_;
   vector<double>  *patJetPfAk05PhotMulti_;
   vector<double>  *patJetPfAk05ElecMulti_;
   vector<double>  *patJetPfAk05HfHadMulti_;
   vector<double>  *patJetPfAk05HfEmMulti_;
   vector<double>  *patJetPfAk05PhotEnFr_;
   vector<double>  *patJetPfAk05MuonEnFr_;
   vector<double>  *patJetPfAk05HfHadEnFr_;
   vector<double>  *patJetPfAk05HfEmEnFr_;
   vector<double>  *patJetPfAk05NeutEmEFr_;
   vector<double>  *patJetPfAk05CharHadEFr_;
   vector<double>  *patJetPfAk05NeutHadEFr_;
   vector<double>  *patJetPfAk05CharEmEFr_;
   vector<double>  *patJetPfAk05CharMuEFr_;
   vector<double>  *patJetPfAk05N60_;
   vector<double>  *patJetPfAk05N90_;
   vector<double>  *patJetPfAk05EmFr_;
   vector<double>  *patJetPfAk05HadFr_;
   vector<double>  *patJetPfAk05EmEbEn_;
   vector<double>  *patJetPfAk05EmEeEn_;
   vector<double>  *patJetPfAk05EmHfEn_;
   vector<double>  *patJetPfAk05HadHbEn_;
   vector<double>  *patJetPfAk05HadHeEn_;
   vector<double>  *patJetPfAk05HadHfEn_;
   vector<double>  *patJetPfAk05HadHoEn_;
   vector<double>  *patJetPfAk05TotalEm_;
   vector<double>  *patJetPfAk05TotalHad_;
   vector<double>  *patJetPfAk05EmEbFr_;
   vector<double>  *patJetPfAk05EmEeFr_;
   vector<double>  *patJetPfAk05EmHfFr_;
   vector<double>  *patJetPfAk05HadHbFr_;
   vector<double>  *patJetPfAk05HadHeFr_;
   vector<double>  *patJetPfAk05HadHfFr_;
   vector<double>  *patJetPfAk05HadHoFr_;
   vector<double>  *patJetPfAk05JesUncert_;
   vector<int>     *patJetPfAk05GenPartonID_;
   vector<double>  *patJetPfAk05GenPartonPt_;
   vector<double>  *patJetPfAk05GenPartonEta_;
   vector<double>  *patJetPfAk05GenPartonPhi_;
   vector<double>  *patJetPfAk05GenPartonE_;
   vector<int>     *trigResults;
   vector<string>  *trigName;

   // List of branches
   TBranch        *b_PU_weight;   //!
   TBranch        *b_PU_nTrueInt;   //!
   TBranch        *b_PU_nPUVert;   //!
   TBranch        *b_EvtInfo_EventNum;   //!
   TBranch        *b_EvtInfo_RunNum;   //!
   TBranch        *b_EvtInfo_LumiSection;   //!
   TBranch        *b_EvtInfo_BunchXing;   //!
   TBranch        *b_EvtInfo_NumVtx;   //!
   TBranch        *b_EvtInfo_VertexX_;   //!
   TBranch        *b_EvtInfo_VertexY_;   //!
   TBranch        *b_EvtInfo_VertexZ_;   //!
   TBranch        *b_ptHat_;   //!
   TBranch        *b_mcWeight_;   //!
   TBranch        *b_genParE_;   //!
   TBranch        *b_genParPt_;   //!
   TBranch        *b_genParEta_;   //!
   TBranch        *b_genParPhi_;   //!
   TBranch        *b_genParQ_;   //!
   TBranch        *b_genParId_;   //!
   TBranch        *b_genParSt_;   //!
   TBranch        *b_genMomParId_;   //!
   TBranch        *b_genParIndex_;   //!
   TBranch        *b_genJetE_;   //!
   TBranch        *b_genJetPt_;   //!
   TBranch        *b_genJetEta_;   //!
   TBranch        *b_genJetPhi_;   //!
   TBranch        *b_fastJetRho_;   //!
   TBranch        *b_PFchsJetRho_;   //!
   TBranch        *b_lepIsoRho_;   //!
   TBranch        *b_patElecEt_;   //!
   TBranch        *b_patElecEnergy_;   //!
   TBranch        *b_patElecPt_;   //!
   TBranch        *b_patElecEta_;   //!
   TBranch        *b_patElecPhi_;   //!
   TBranch        *b_patElecM_;   //!
   TBranch        *b_patElecScEn_;   //!
   TBranch        *b_patElecScEt_;   //!
   TBranch        *b_patElecScEta_;   //!
   TBranch        *b_patElecScPhi_;   //!
   TBranch        *b_patElecisEcalDriven_;   //!
   TBranch        *b_patElecisTrackerDriven_;   //!
   TBranch        *b_patElecSigIhIh_;   //!
   TBranch        *b_patElecDelEtaIn_;   //!
   TBranch        *b_patElecDelPhiIn_;   //!
   TBranch        *b_patElecHoE_;   //!
   TBranch        *b_patElecTrkIso_;   //!
   TBranch        *b_patElecHcalIso_;   //!
   TBranch        *b_patElecEcalIso_;   //!
   TBranch        *b_patElecCharge_;   //!
   TBranch        *b_patElecRelIsoComb_;   //!
   TBranch        *b_patElecDxy_;   //!
   TBranch        *b_patElecD0_;   //!
   TBranch        *b_patElecDsz_;   //!
   TBranch        *b_patElecDz_;   //!
   TBranch        *b_patElecDxyBS_;   //!
   TBranch        *b_patElecDszBS_;   //!
   TBranch        *b_patElecDzBS_;   //!
   TBranch        *b_patElecMva_;   //!
   TBranch        *b_patElecChHadSumPt03_;   //!
   TBranch        *b_patElecNeHadSumPt03_;   //!
   TBranch        *b_patElecGamSumPt03_;   //!
   TBranch        *b_patElecChHadSumPt04_;   //!
   TBranch        *b_patElecChHadIso_;   //!
   TBranch        *b_patElecNeHadIso_;   //!
   TBranch        *b_patElecGamIso_;   //!
   TBranch        *b_patElecNeHadSumPt04_;   //!
   TBranch        *b_patElecGamSumPt04_;   //!
   TBranch        *b_patElecChHadSumPt05_;   //!
   TBranch        *b_patElecNeHadSumPt05_;   //!
   TBranch        *b_patElecGamSumPt05_;   //!
   TBranch        *b_patElecMissingHits_;   //!
   TBranch        *b_patElecDist_;   //!
   TBranch        *b_patElecDeltaCotTheta_;   //!
   TBranch        *b_patElecConvRadius_;   //!
   TBranch        *b_patElecInBarrel_;   //!
   TBranch        *b_patElecInEndcap_;   //!
   TBranch        *b_patJetPfAk05Pt_;   //!
   TBranch        *b_patJetPfAk05Eta_;   //!
   TBranch        *b_patJetPfAk05Phi_;   //!
   TBranch        *b_patJetPfAk05M_;   //!
   TBranch        *b_patJetPfAk05Rapidity_;   //!
   TBranch        *b_patJetPfAk05Px_;   //!
   TBranch        *b_patJetPfAk05Py_;   //!
   TBranch        *b_patJetPfAk05Pz_;   //!
   TBranch        *b_patJetPfAk05En_;   //!
   TBranch        *b_patJetPfAk05UnCorrPt_;   //!
   TBranch        *b_patJetPfAk05UnCorrPx_;   //!
   TBranch        *b_patJetPfAk05UnCorrPy_;   //!
   TBranch        *b_patJetPfAk05UnCorrPz_;   //!
   TBranch        *b_patJetPfAk05UnCorrEnt_;   //!
   TBranch        *b_patJetPfAk05PhotEn_;   //!
   TBranch        *b_patJetPfAk05ElecEn_;   //!
   TBranch        *b_patJetPfAk05MuonEn_;   //!
   TBranch        *b_patJetPfAk05HfHadEn_;   //!
   TBranch        *b_patJetPfAk05HfEmEn_;   //!
   TBranch        *b_patJetPfAk05CharHadE_;   //!
   TBranch        *b_patJetPfAk05NeutHadE_;   //!
   TBranch        *b_patJetPfAk05CharEmE_;   //!
   TBranch        *b_patJetPfAk05CharMuE_;   //!
   TBranch        *b_patJetPfAk05NeutEmE_;   //!
   TBranch        *b_patJetPfAk05MuonMulti_;   //!
   TBranch        *b_patJetPfAk05NeutMulti_;   //!
   TBranch        *b_patJetPfAk05CharMulti_;   //!
   TBranch        *b_patJetPfAk05CharHadMulti_;   //!
   TBranch        *b_patJetPfAk05NeutHadMulti_;   //!
   TBranch        *b_patJetPfAk05PhotMulti_;   //!
   TBranch        *b_patJetPfAk05ElecMulti_;   //!
   TBranch        *b_patJetPfAk05HfHadMulti_;   //!
   TBranch        *b_patJetPfAk05HfEmMulti_;   //!
   TBranch        *b_patJetPfAk05PhotEnFr_;   //!
   TBranch        *b_patJetPfAk05MuonEnFr_;   //!
   TBranch        *b_patJetPfAk05HfHadEnFr_;   //!
   TBranch        *b_patJetPfAk05HfEmEnFr_;   //!
   TBranch        *b_patJetPfAk05NeutEmEFr_;   //!
   TBranch        *b_patJetPfAk05CharHadEFr_;   //!
   TBranch        *b_patJetPfAk05NeutHadEFr_;   //!
   TBranch        *b_patJetPfAk05CharEmEFr_;   //!
   TBranch        *b_patJetPfAk05CharMuEFr_;   //!
   TBranch        *b_patJetPfAk05N60_;   //!
   TBranch        *b_patJetPfAk05N90_;   //!
   TBranch        *b_patJetPfAk05EmFr_;   //!
   TBranch        *b_patJetPfAk05HadFr_;   //!
   TBranch        *b_patJetPfAk05EmEbEn_;   //!
   TBranch        *b_patJetPfAk05EmEeEn_;   //!
   TBranch        *b_patJetPfAk05EmHfEn_;   //!
   TBranch        *b_patJetPfAk05HadHbEn_;   //!
   TBranch        *b_patJetPfAk05HadHeEn_;   //!
   TBranch        *b_patJetPfAk05HadHfEn_;   //!
   TBranch        *b_patJetPfAk05HadHoEn_;   //!
   TBranch        *b_patJetPfAk05TotalEm_;   //!
   TBranch        *b_patJetPfAk05TotalHad_;   //!
   TBranch        *b_patJetPfAk05EmEbFr_;   //!
   TBranch        *b_patJetPfAk05EmEeFr_;   //!
   TBranch        *b_patJetPfAk05EmHfFr_;   //!
   TBranch        *b_patJetPfAk05HadHbFr_;   //!
   TBranch        *b_patJetPfAk05HadHeFr_;   //!
   TBranch        *b_patJetPfAk05HadHfFr_;   //!
   TBranch        *b_patJetPfAk05HadHoFr_;   //!
   TBranch        *b_patJetPfAk05JesUncert_;   //!
   TBranch        *b_patJetPfAk05GenPartonID_;   //!
   TBranch        *b_patJetPfAk05GenPartonPt_;   //!
   TBranch        *b_patJetPfAk05GenPartonEta_;   //!
   TBranch        *b_patJetPfAk05GenPartonPhi_;   //!
   TBranch        *b_patJetPfAk05GenPartonE_;   //!
   TBranch        *b_trigResults;   //!
   TBranch        *b_trigName;   //!

   puweight_sys(std::string inputfile, TTree *tree=0);
   virtual ~puweight_sys();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(bool match=false);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Bool_t   isFidEle (Int_t iele);
   virtual Bool_t   isWP80VBTF11Ele(Int_t iele);

   virtual Bool_t   isFidJet (Int_t ijet);
   virtual Bool_t   isGoodLooseJet(Int_t ijet);
   virtual Int_t    matchRecoToParton(Int_t ijet); // return matched index in genInfo
   virtual Int_t    matchRecoToGenJet(Int_t ijet); // return matched index in genjet
   std::string _inputFile;
};

#endif

#ifdef puweight_sys_cxx
puweight_sys::puweight_sys(std::string filename, TTree *tree)
{
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename.data());
    if (!f) {
      f = new TFile(filename.data());
      f->cd("tree");
    }
    tree = (TTree*)gDirectory->Get("tree");

  }
  Init(tree);
  _inputFile = filename;
}

puweight_sys::~puweight_sys()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t puweight_sys::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t puweight_sys::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void puweight_sys::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   EvtInfo_VertexX_ = 0;
   EvtInfo_VertexY_ = 0;
   EvtInfo_VertexZ_ = 0;
   genParE_ = 0;
   genParPt_ = 0;
   genParEta_ = 0;
   genParPhi_ = 0;
   genParQ_ = 0;
   genParId_ = 0;
   genParSt_ = 0;
   genMomParId_ = 0;
   genParIndex_ = 0;
   genJetE_ = 0;
   genJetPt_ = 0;
   genJetEta_ = 0;
   genJetPhi_ = 0;
   fastJetRho_ = 0;
   PFchsJetRho_ = 0;
   lepIsoRho_ = 0;
   patElecEt_ = 0;
   patElecEnergy_ = 0;
   patElecPt_ = 0;
   patElecEta_ = 0;
   patElecPhi_ = 0;
   patElecM_ = 0;
   patElecScEn_ = 0;
   patElecScEt_ = 0;
   patElecScEta_ = 0;
   patElecScPhi_ = 0;
   patElecisEcalDriven_ = 0;
   patElecisTrackerDriven_ = 0;
   patElecSigIhIh_ = 0;
   patElecDelEtaIn_ = 0;
   patElecDelPhiIn_ = 0;
   patElecHoE_ = 0;
   patElecTrkIso_ = 0;
   patElecHcalIso_ = 0;
   patElecEcalIso_ = 0;
   patElecCharge_ = 0;
   patElecRelIsoComb_ = 0;
   patElecDxy_ = 0;
   patElecD0_ = 0;
   patElecDsz_ = 0;
   patElecDz_ = 0;
   patElecDxyBS_ = 0;
   patElecDszBS_ = 0;
   patElecDzBS_ = 0;
   patElecMva_ = 0;
   patElecChHadSumPt03_ = 0;
   patElecNeHadSumPt03_ = 0;
   patElecGamSumPt03_ = 0;
   patElecChHadSumPt04_ = 0;
   patElecChHadIso_ = 0;
   patElecNeHadIso_ = 0;
   patElecGamIso_ = 0;
   patElecNeHadSumPt04_ = 0;
   patElecGamSumPt04_ = 0;
   patElecChHadSumPt05_ = 0;
   patElecNeHadSumPt05_ = 0;
   patElecGamSumPt05_ = 0;
   patElecMissingHits_ = 0;
   patElecDist_ = 0;
   patElecDeltaCotTheta_ = 0;
   patElecConvRadius_ = 0;
   patElecInBarrel_ = 0;
   patElecInEndcap_ = 0;
   patJetPfAk05Pt_ = 0;
   patJetPfAk05Eta_ = 0;
   patJetPfAk05Phi_ = 0;
   patJetPfAk05M_ = 0;
   patJetPfAk05Rapidity_ = 0;
   patJetPfAk05Px_ = 0;
   patJetPfAk05Py_ = 0;
   patJetPfAk05Pz_ = 0;
   patJetPfAk05En_ = 0;
   patJetPfAk05UnCorrPt_ = 0;
   patJetPfAk05UnCorrPx_ = 0;
   patJetPfAk05UnCorrPy_ = 0;
   patJetPfAk05UnCorrPz_ = 0;
   patJetPfAk05UnCorrEnt_ = 0;
   patJetPfAk05PhotEn_ = 0;
   patJetPfAk05ElecEn_ = 0;
   patJetPfAk05MuonEn_ = 0;
   patJetPfAk05HfHadEn_ = 0;
   patJetPfAk05HfEmEn_ = 0;
   patJetPfAk05CharHadE_ = 0;
   patJetPfAk05NeutHadE_ = 0;
   patJetPfAk05CharEmE_ = 0;
   patJetPfAk05CharMuE_ = 0;
   patJetPfAk05NeutEmE_ = 0;
   patJetPfAk05MuonMulti_ = 0;
   patJetPfAk05NeutMulti_ = 0;
   patJetPfAk05CharMulti_ = 0;
   patJetPfAk05CharHadMulti_ = 0;
   patJetPfAk05NeutHadMulti_ = 0;
   patJetPfAk05PhotMulti_ = 0;
   patJetPfAk05ElecMulti_ = 0;
   patJetPfAk05HfHadMulti_ = 0;
   patJetPfAk05HfEmMulti_ = 0;
   patJetPfAk05PhotEnFr_ = 0;
   patJetPfAk05MuonEnFr_ = 0;
   patJetPfAk05HfHadEnFr_ = 0;
   patJetPfAk05HfEmEnFr_ = 0;
   patJetPfAk05NeutEmEFr_ = 0;
   patJetPfAk05CharHadEFr_ = 0;
   patJetPfAk05NeutHadEFr_ = 0;
   patJetPfAk05CharEmEFr_ = 0;
   patJetPfAk05CharMuEFr_ = 0;
   patJetPfAk05N60_ = 0;
   patJetPfAk05N90_ = 0;
   patJetPfAk05EmFr_ = 0;
   patJetPfAk05HadFr_ = 0;
   patJetPfAk05EmEbEn_ = 0;
   patJetPfAk05EmEeEn_ = 0;
   patJetPfAk05EmHfEn_ = 0;
   patJetPfAk05HadHbEn_ = 0;
   patJetPfAk05HadHeEn_ = 0;
   patJetPfAk05HadHfEn_ = 0;
   patJetPfAk05HadHoEn_ = 0;
   patJetPfAk05TotalEm_ = 0;
   patJetPfAk05TotalHad_ = 0;
   patJetPfAk05EmEbFr_ = 0;
   patJetPfAk05EmEeFr_ = 0;
   patJetPfAk05EmHfFr_ = 0;
   patJetPfAk05HadHbFr_ = 0;
   patJetPfAk05HadHeFr_ = 0;
   patJetPfAk05HadHfFr_ = 0;
   patJetPfAk05HadHoFr_ = 0;
   patJetPfAk05JesUncert_ = 0;
   patJetPfAk05GenPartonID_ = 0;
   patJetPfAk05GenPartonPt_ = 0;
   patJetPfAk05GenPartonEta_ = 0;
   patJetPfAk05GenPartonPhi_ = 0;
   patJetPfAk05GenPartonE_ = 0;
   trigResults = 0;
   trigName = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("PU_weight", &PU_weight, &b_PU_weight);
   fChain->SetBranchAddress("PU_nTrueInt", &PU_nTrueInt, &b_PU_nTrueInt);
   fChain->SetBranchAddress("PU_nPUVert", &PU_nPUVert, &b_PU_nPUVert);
   fChain->SetBranchAddress("EvtInfo_EventNum", &EvtInfo_EventNum, &b_EvtInfo_EventNum);
   fChain->SetBranchAddress("EvtInfo_RunNum", &EvtInfo_RunNum, &b_EvtInfo_RunNum);
   fChain->SetBranchAddress("EvtInfo_LumiSection", &EvtInfo_LumiSection, &b_EvtInfo_LumiSection);
   fChain->SetBranchAddress("EvtInfo_BunchXing", &EvtInfo_BunchXing, &b_EvtInfo_BunchXing);
   fChain->SetBranchAddress("EvtInfo_NumVtx", &EvtInfo_NumVtx, &b_EvtInfo_NumVtx);
   fChain->SetBranchAddress("EvtInfo_VertexX_", &EvtInfo_VertexX_, &b_EvtInfo_VertexX_);
   fChain->SetBranchAddress("EvtInfo_VertexY_", &EvtInfo_VertexY_, &b_EvtInfo_VertexY_);
   fChain->SetBranchAddress("EvtInfo_VertexZ_", &EvtInfo_VertexZ_, &b_EvtInfo_VertexZ_);
   fChain->SetBranchAddress("ptHat_", &ptHat_, &b_ptHat_);
   fChain->SetBranchAddress("mcWeight_", &mcWeight_, &b_mcWeight_);
   fChain->SetBranchAddress("genParE_", &genParE_, &b_genParE_);
   fChain->SetBranchAddress("genParPt_", &genParPt_, &b_genParPt_);
   fChain->SetBranchAddress("genParEta_", &genParEta_, &b_genParEta_);
   fChain->SetBranchAddress("genParPhi_", &genParPhi_, &b_genParPhi_);
   fChain->SetBranchAddress("genParQ_", &genParQ_, &b_genParQ_);
   fChain->SetBranchAddress("genParId_", &genParId_, &b_genParId_);
   fChain->SetBranchAddress("genParSt_", &genParSt_, &b_genParSt_);
   fChain->SetBranchAddress("genMomParId_", &genMomParId_, &b_genMomParId_);
   fChain->SetBranchAddress("genParIndex_", &genParIndex_, &b_genParIndex_);
   fChain->SetBranchAddress("genJetE_", &genJetE_, &b_genJetE_);
   fChain->SetBranchAddress("genJetPt_", &genJetPt_, &b_genJetPt_);
   fChain->SetBranchAddress("genJetEta_", &genJetEta_, &b_genJetEta_);
   fChain->SetBranchAddress("genJetPhi_", &genJetPhi_, &b_genJetPhi_);
   fChain->SetBranchAddress("fastJetRho_", &fastJetRho_, &b_fastJetRho_);
   fChain->SetBranchAddress("PFchsJetRho_", &PFchsJetRho_, &b_PFchsJetRho_);
   fChain->SetBranchAddress("lepIsoRho_", &lepIsoRho_, &b_lepIsoRho_);
   fChain->SetBranchAddress("patElecEt_", &patElecEt_, &b_patElecEt_);
   fChain->SetBranchAddress("patElecEnergy_", &patElecEnergy_, &b_patElecEnergy_);
   fChain->SetBranchAddress("patElecPt_", &patElecPt_, &b_patElecPt_);
   fChain->SetBranchAddress("patElecEta_", &patElecEta_, &b_patElecEta_);
   fChain->SetBranchAddress("patElecPhi_", &patElecPhi_, &b_patElecPhi_);
   fChain->SetBranchAddress("patElecM_", &patElecM_, &b_patElecM_);
   fChain->SetBranchAddress("patElecScEn_", &patElecScEn_, &b_patElecScEn_);
   fChain->SetBranchAddress("patElecScEt_", &patElecScEt_, &b_patElecScEt_);
   fChain->SetBranchAddress("patElecScEta_", &patElecScEta_, &b_patElecScEta_);
   fChain->SetBranchAddress("patElecScPhi_", &patElecScPhi_, &b_patElecScPhi_);
   fChain->SetBranchAddress("patElecisEcalDriven_", &patElecisEcalDriven_, &b_patElecisEcalDriven_);
   fChain->SetBranchAddress("patElecisTrackerDriven_", &patElecisTrackerDriven_, &b_patElecisTrackerDriven_);
   fChain->SetBranchAddress("patElecSigIhIh_", &patElecSigIhIh_, &b_patElecSigIhIh_);
   fChain->SetBranchAddress("patElecDelEtaIn_", &patElecDelEtaIn_, &b_patElecDelEtaIn_);
   fChain->SetBranchAddress("patElecDelPhiIn_", &patElecDelPhiIn_, &b_patElecDelPhiIn_);
   fChain->SetBranchAddress("patElecHoE_", &patElecHoE_, &b_patElecHoE_);
   fChain->SetBranchAddress("patElecTrkIso_", &patElecTrkIso_, &b_patElecTrkIso_);
   fChain->SetBranchAddress("patElecHcalIso_", &patElecHcalIso_, &b_patElecHcalIso_);
   fChain->SetBranchAddress("patElecEcalIso_", &patElecEcalIso_, &b_patElecEcalIso_);
   fChain->SetBranchAddress("patElecCharge_", &patElecCharge_, &b_patElecCharge_);
   fChain->SetBranchAddress("patElecRelIsoComb_", &patElecRelIsoComb_, &b_patElecRelIsoComb_);
   fChain->SetBranchAddress("patElecDxy_", &patElecDxy_, &b_patElecDxy_);
   fChain->SetBranchAddress("patElecD0_", &patElecD0_, &b_patElecD0_);
   fChain->SetBranchAddress("patElecDsz_", &patElecDsz_, &b_patElecDsz_);
   fChain->SetBranchAddress("patElecDz_", &patElecDz_, &b_patElecDz_);
   fChain->SetBranchAddress("patElecDxyBS_", &patElecDxyBS_, &b_patElecDxyBS_);
   fChain->SetBranchAddress("patElecDszBS_", &patElecDszBS_, &b_patElecDszBS_);
   fChain->SetBranchAddress("patElecDzBS_", &patElecDzBS_, &b_patElecDzBS_);
   fChain->SetBranchAddress("patElecMva_", &patElecMva_, &b_patElecMva_);
   fChain->SetBranchAddress("patElecChHadSumPt03_", &patElecChHadSumPt03_, &b_patElecChHadSumPt03_);
   fChain->SetBranchAddress("patElecNeHadSumPt03_", &patElecNeHadSumPt03_, &b_patElecNeHadSumPt03_);
   fChain->SetBranchAddress("patElecGamSumPt03_", &patElecGamSumPt03_, &b_patElecGamSumPt03_);
   fChain->SetBranchAddress("patElecChHadSumPt04_", &patElecChHadSumPt04_, &b_patElecChHadSumPt04_);
   fChain->SetBranchAddress("patElecChHadIso_", &patElecChHadIso_, &b_patElecChHadIso_);
   fChain->SetBranchAddress("patElecNeHadIso_", &patElecNeHadIso_, &b_patElecNeHadIso_);
   fChain->SetBranchAddress("patElecGamIso_", &patElecGamIso_, &b_patElecGamIso_);
   fChain->SetBranchAddress("patElecNeHadSumPt04_", &patElecNeHadSumPt04_, &b_patElecNeHadSumPt04_);
   fChain->SetBranchAddress("patElecGamSumPt04_", &patElecGamSumPt04_, &b_patElecGamSumPt04_);
   fChain->SetBranchAddress("patElecChHadSumPt05_", &patElecChHadSumPt05_, &b_patElecChHadSumPt05_);
   fChain->SetBranchAddress("patElecNeHadSumPt05_", &patElecNeHadSumPt05_, &b_patElecNeHadSumPt05_);
   fChain->SetBranchAddress("patElecGamSumPt05_", &patElecGamSumPt05_, &b_patElecGamSumPt05_);
   fChain->SetBranchAddress("patElecMissingHits_", &patElecMissingHits_, &b_patElecMissingHits_);
   fChain->SetBranchAddress("patElecDist_", &patElecDist_, &b_patElecDist_);
   fChain->SetBranchAddress("patElecDeltaCotTheta_", &patElecDeltaCotTheta_, &b_patElecDeltaCotTheta_);
   fChain->SetBranchAddress("patElecConvRadius_", &patElecConvRadius_, &b_patElecConvRadius_);
   fChain->SetBranchAddress("patElecInBarrel_", &patElecInBarrel_, &b_patElecInBarrel_);
   fChain->SetBranchAddress("patElecInEndcap_", &patElecInEndcap_, &b_patElecInEndcap_);
   fChain->SetBranchAddress("patJetPfAk05Pt_", &patJetPfAk05Pt_, &b_patJetPfAk05Pt_);
   fChain->SetBranchAddress("patJetPfAk05Eta_", &patJetPfAk05Eta_, &b_patJetPfAk05Eta_);
   fChain->SetBranchAddress("patJetPfAk05Phi_", &patJetPfAk05Phi_, &b_patJetPfAk05Phi_);
   fChain->SetBranchAddress("patJetPfAk05M_", &patJetPfAk05M_, &b_patJetPfAk05M_);
   fChain->SetBranchAddress("patJetPfAk05Rapidity_", &patJetPfAk05Rapidity_, &b_patJetPfAk05Rapidity_);
   fChain->SetBranchAddress("patJetPfAk05Px_", &patJetPfAk05Px_, &b_patJetPfAk05Px_);
   fChain->SetBranchAddress("patJetPfAk05Py_", &patJetPfAk05Py_, &b_patJetPfAk05Py_);
   fChain->SetBranchAddress("patJetPfAk05Pz_", &patJetPfAk05Pz_, &b_patJetPfAk05Pz_);
   fChain->SetBranchAddress("patJetPfAk05En_", &patJetPfAk05En_, &b_patJetPfAk05En_);
   fChain->SetBranchAddress("patJetPfAk05UnCorrPt_", &patJetPfAk05UnCorrPt_, &b_patJetPfAk05UnCorrPt_);
   fChain->SetBranchAddress("patJetPfAk05UnCorrPx_", &patJetPfAk05UnCorrPx_, &b_patJetPfAk05UnCorrPx_);
   fChain->SetBranchAddress("patJetPfAk05UnCorrPy_", &patJetPfAk05UnCorrPy_, &b_patJetPfAk05UnCorrPy_);
   fChain->SetBranchAddress("patJetPfAk05UnCorrPz_", &patJetPfAk05UnCorrPz_, &b_patJetPfAk05UnCorrPz_);
   fChain->SetBranchAddress("patJetPfAk05UnCorrEnt_", &patJetPfAk05UnCorrEnt_, &b_patJetPfAk05UnCorrEnt_);
   fChain->SetBranchAddress("patJetPfAk05PhotEn_", &patJetPfAk05PhotEn_, &b_patJetPfAk05PhotEn_);
   fChain->SetBranchAddress("patJetPfAk05ElecEn_", &patJetPfAk05ElecEn_, &b_patJetPfAk05ElecEn_);
   fChain->SetBranchAddress("patJetPfAk05MuonEn_", &patJetPfAk05MuonEn_, &b_patJetPfAk05MuonEn_);
   fChain->SetBranchAddress("patJetPfAk05HfHadEn_", &patJetPfAk05HfHadEn_, &b_patJetPfAk05HfHadEn_);
   fChain->SetBranchAddress("patJetPfAk05HfEmEn_", &patJetPfAk05HfEmEn_, &b_patJetPfAk05HfEmEn_);
   fChain->SetBranchAddress("patJetPfAk05CharHadE_", &patJetPfAk05CharHadE_, &b_patJetPfAk05CharHadE_);
   fChain->SetBranchAddress("patJetPfAk05NeutHadE_", &patJetPfAk05NeutHadE_, &b_patJetPfAk05NeutHadE_);
   fChain->SetBranchAddress("patJetPfAk05CharEmE_", &patJetPfAk05CharEmE_, &b_patJetPfAk05CharEmE_);
   fChain->SetBranchAddress("patJetPfAk05CharMuE_", &patJetPfAk05CharMuE_, &b_patJetPfAk05CharMuE_);
   fChain->SetBranchAddress("patJetPfAk05NeutEmE_", &patJetPfAk05NeutEmE_, &b_patJetPfAk05NeutEmE_);
   fChain->SetBranchAddress("patJetPfAk05MuonMulti_", &patJetPfAk05MuonMulti_, &b_patJetPfAk05MuonMulti_);
   fChain->SetBranchAddress("patJetPfAk05NeutMulti_", &patJetPfAk05NeutMulti_, &b_patJetPfAk05NeutMulti_);
   fChain->SetBranchAddress("patJetPfAk05CharMulti_", &patJetPfAk05CharMulti_, &b_patJetPfAk05CharMulti_);
   fChain->SetBranchAddress("patJetPfAk05CharHadMulti_", &patJetPfAk05CharHadMulti_, &b_patJetPfAk05CharHadMulti_);
   fChain->SetBranchAddress("patJetPfAk05NeutHadMulti_", &patJetPfAk05NeutHadMulti_, &b_patJetPfAk05NeutHadMulti_);
   fChain->SetBranchAddress("patJetPfAk05PhotMulti_", &patJetPfAk05PhotMulti_, &b_patJetPfAk05PhotMulti_);
   fChain->SetBranchAddress("patJetPfAk05ElecMulti_", &patJetPfAk05ElecMulti_, &b_patJetPfAk05ElecMulti_);
   fChain->SetBranchAddress("patJetPfAk05HfHadMulti_", &patJetPfAk05HfHadMulti_, &b_patJetPfAk05HfHadMulti_);
   fChain->SetBranchAddress("patJetPfAk05HfEmMulti_", &patJetPfAk05HfEmMulti_, &b_patJetPfAk05HfEmMulti_);
   fChain->SetBranchAddress("patJetPfAk05PhotEnFr_", &patJetPfAk05PhotEnFr_, &b_patJetPfAk05PhotEnFr_);
   fChain->SetBranchAddress("patJetPfAk05MuonEnFr_", &patJetPfAk05MuonEnFr_, &b_patJetPfAk05MuonEnFr_);
   fChain->SetBranchAddress("patJetPfAk05HfHadEnFr_", &patJetPfAk05HfHadEnFr_, &b_patJetPfAk05HfHadEnFr_);
   fChain->SetBranchAddress("patJetPfAk05HfEmEnFr_", &patJetPfAk05HfEmEnFr_, &b_patJetPfAk05HfEmEnFr_);
   fChain->SetBranchAddress("patJetPfAk05NeutEmEFr_", &patJetPfAk05NeutEmEFr_, &b_patJetPfAk05NeutEmEFr_);
   fChain->SetBranchAddress("patJetPfAk05CharHadEFr_", &patJetPfAk05CharHadEFr_, &b_patJetPfAk05CharHadEFr_);
   fChain->SetBranchAddress("patJetPfAk05NeutHadEFr_", &patJetPfAk05NeutHadEFr_, &b_patJetPfAk05NeutHadEFr_);
   fChain->SetBranchAddress("patJetPfAk05CharEmEFr_", &patJetPfAk05CharEmEFr_, &b_patJetPfAk05CharEmEFr_);
   fChain->SetBranchAddress("patJetPfAk05CharMuEFr_", &patJetPfAk05CharMuEFr_, &b_patJetPfAk05CharMuEFr_);
   fChain->SetBranchAddress("patJetPfAk05N60_", &patJetPfAk05N60_, &b_patJetPfAk05N60_);
   fChain->SetBranchAddress("patJetPfAk05N90_", &patJetPfAk05N90_, &b_patJetPfAk05N90_);
   fChain->SetBranchAddress("patJetPfAk05EmFr_", &patJetPfAk05EmFr_, &b_patJetPfAk05EmFr_);
   fChain->SetBranchAddress("patJetPfAk05HadFr_", &patJetPfAk05HadFr_, &b_patJetPfAk05HadFr_);
   fChain->SetBranchAddress("patJetPfAk05EmEbEn_", &patJetPfAk05EmEbEn_, &b_patJetPfAk05EmEbEn_);
   fChain->SetBranchAddress("patJetPfAk05EmEeEn_", &patJetPfAk05EmEeEn_, &b_patJetPfAk05EmEeEn_);
   fChain->SetBranchAddress("patJetPfAk05EmHfEn_", &patJetPfAk05EmHfEn_, &b_patJetPfAk05EmHfEn_);
   fChain->SetBranchAddress("patJetPfAk05HadHbEn_", &patJetPfAk05HadHbEn_, &b_patJetPfAk05HadHbEn_);
   fChain->SetBranchAddress("patJetPfAk05HadHeEn_", &patJetPfAk05HadHeEn_, &b_patJetPfAk05HadHeEn_);
   fChain->SetBranchAddress("patJetPfAk05HadHfEn_", &patJetPfAk05HadHfEn_, &b_patJetPfAk05HadHfEn_);
   fChain->SetBranchAddress("patJetPfAk05HadHoEn_", &patJetPfAk05HadHoEn_, &b_patJetPfAk05HadHoEn_);
   fChain->SetBranchAddress("patJetPfAk05TotalEm_", &patJetPfAk05TotalEm_, &b_patJetPfAk05TotalEm_);
   fChain->SetBranchAddress("patJetPfAk05TotalHad_", &patJetPfAk05TotalHad_, &b_patJetPfAk05TotalHad_);
   fChain->SetBranchAddress("patJetPfAk05EmEbFr_", &patJetPfAk05EmEbFr_, &b_patJetPfAk05EmEbFr_);
   fChain->SetBranchAddress("patJetPfAk05EmEeFr_", &patJetPfAk05EmEeFr_, &b_patJetPfAk05EmEeFr_);
   fChain->SetBranchAddress("patJetPfAk05EmHfFr_", &patJetPfAk05EmHfFr_, &b_patJetPfAk05EmHfFr_);
   fChain->SetBranchAddress("patJetPfAk05HadHbFr_", &patJetPfAk05HadHbFr_, &b_patJetPfAk05HadHbFr_);
   fChain->SetBranchAddress("patJetPfAk05HadHeFr_", &patJetPfAk05HadHeFr_, &b_patJetPfAk05HadHeFr_);
   fChain->SetBranchAddress("patJetPfAk05HadHfFr_", &patJetPfAk05HadHfFr_, &b_patJetPfAk05HadHfFr_);
   fChain->SetBranchAddress("patJetPfAk05HadHoFr_", &patJetPfAk05HadHoFr_, &b_patJetPfAk05HadHoFr_);
   fChain->SetBranchAddress("patJetPfAk05JesUncert_", &patJetPfAk05JesUncert_, &b_patJetPfAk05JesUncert_);
   fChain->SetBranchAddress("patJetPfAk05GenPartonID_", &patJetPfAk05GenPartonID_, &b_patJetPfAk05GenPartonID_);
   fChain->SetBranchAddress("patJetPfAk05GenPartonPt_", &patJetPfAk05GenPartonPt_, &b_patJetPfAk05GenPartonPt_);
   fChain->SetBranchAddress("patJetPfAk05GenPartonEta_", &patJetPfAk05GenPartonEta_, &b_patJetPfAk05GenPartonEta_);
   fChain->SetBranchAddress("patJetPfAk05GenPartonPhi_", &patJetPfAk05GenPartonPhi_, &b_patJetPfAk05GenPartonPhi_);
   fChain->SetBranchAddress("patJetPfAk05GenPartonE_", &patJetPfAk05GenPartonE_, &b_patJetPfAk05GenPartonE_);
   fChain->SetBranchAddress("trigResults", &trigResults, &b_trigResults);
   fChain->SetBranchAddress("trigName", &trigName, &b_trigName);
   Notify();
}

Bool_t puweight_sys::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void puweight_sys::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t puweight_sys::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef puweight_sys_cxx
