//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb 17 19:15:28 2012 by ROOT version 5.27/06b
// from TTree EventTree/Event data
// found on file: /data4/syu/TreeMC_G_Pt-80to120/ggtree_mc_1_1_Asg.root
//////////////////////////////////////////////////////////

#ifndef gg_angularmc_eff_h
#define gg_angularmc_eff_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream>
#include <string>

using namespace std;

   const Int_t MAXN = 1000;
   const Int_t kMaxeleESEffSigmaRR = 1;
   const Int_t kMaxeleESE1 = 1;
   const Int_t kMaxeleESE3 = 1;
   const Int_t kMaxeleESE5 = 1;
   const Int_t kMaxeleESE7 = 1;
   const Int_t kMaxeleESE11 = 1;
   const Int_t kMaxeleESE21 = 1;
   const Int_t kMaxphoESEffSigmaRR = 1;
   const Int_t kMaxphoESE1 = 1;
   const Int_t kMaxphoESE3 = 1;
   const Int_t kMaxphoESE5 = 1;
   const Int_t kMaxphoESE7 = 1;
   const Int_t kMaxphoESE11 = 1;
   const Int_t kMaxphoESE21 = 1;
   const Int_t kMaxnPFPho = 388;
   const Int_t kMaxPFPhoE = 1;
   const Int_t kMaxPFPhoEt = 1;
   const Int_t kMaxPFPhoEta = 1;
   const Int_t kMaxPFPhoPhi = 1;
   const Int_t kMaxPFPhoHoverE = 1;
   const Int_t kMaxnPFEle = 2;
   const Int_t kMaxPFElePt = 1;
   const Int_t kMaxPFEleEta = 1;
   const Int_t kMaxPFElePhi = 1;
   const Int_t kMaxPFEleEn = 1;
   const Int_t kMaxPFEleCharge = 1;
   const Int_t kMaxnPFMu = 4;
   const Int_t kMaxPFmuCharge = 1;
   const Int_t kMaxPFmuPhi = 1;
   const Int_t kMaxPFmuEta = 1;
   const Int_t kMaxPFmuPt = 1;
   const Int_t kMaxPFmuPz = 1;
   const Int_t kMaxPFmuChambers = 1;
   const Int_t kMaxPFmuD0 = 1;
   const Int_t kMaxPFmuDz = 1;
   const Int_t kMaxPFmuChi2NDF = 1;
   const Int_t kMaxPFmuNumberOfValidTrkHits = 1;
   const Int_t kMaxPFmuNumberOfValidPixelHits = 1;
   const Int_t kMaxPFmuNumberOfValidMuonHits = 1;

class gg_angularmc_eff {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           lumis;
   Bool_t          isData;
   Int_t           nHLT;
   Int_t           HLT[239];   //[nHLT]
   Int_t           HLTIndex[70];
   Float_t         bspotPos[3];
   Int_t           nVtx;
   Float_t         vtx[MAXN][3];   //[nVtx]
   Int_t           IsVtxGood;
   Int_t           nVtxBS;
   Float_t         vtxbs[MAXN][3];   //[nVtxBS]
   Float_t         vtxbsPtMod[MAXN];   //[nVtxBS]
   Float_t         vtxbsSumPt2[MAXN];   //[nVtxBS]
   vector<vector<int> > *vtxbsTkIndex;
   vector<vector<float> > *vtxbsTkWeight;
   Int_t           nTrk;
   Float_t         trkP[MAXN][3];   //[nTrk]
   Float_t         trkVtx[MAXN][3];   //[nTrk]
   Float_t         trkd0[MAXN];   //[nTrk]
   Float_t         trkd0Err[MAXN];   //[nTrk]
   Float_t         trkdz[MAXN];   //[nTrk]
   Float_t         trkdzErr[MAXN];   //[nTrk]
   Float_t         trkPtErr[MAXN];   //[nTrk]
   Int_t           trkQuality[MAXN];   //[nTrk]
   Int_t           nGoodTrk;
   Int_t           IsTracksGood;
   Float_t         pdf[7];
   Float_t         pthat;
   Float_t         processID;
   Int_t           nMC;
   Int_t           mcPID[MAXN];   //[nMC]
   Float_t         mcVtx[MAXN][3];   //[nMC]
   Float_t         mcPt[MAXN];   //[nMC]
   Float_t         mcMass[MAXN];   //[nMC]
   Float_t         mcEta[MAXN];   //[nMC]
   Float_t         mcPhi[MAXN];   //[nMC]
   Float_t         mcE[MAXN];   //[nMC]
   Float_t         mcEt[MAXN];   //[nMC]
   Int_t           mcGMomPID[MAXN];   //[nMC]
   Int_t           mcMomPID[MAXN];   //[nMC]
   Float_t         mcMomPt[MAXN];   //[nMC]
   Float_t         mcMomMass[MAXN];   //[nMC]
   Float_t         mcMomEta[MAXN];   //[nMC]
   Float_t         mcMomPhi[MAXN];   //[nMC]
   Int_t           mcIndex[MAXN];   //[nMC]
   Int_t           mcDecayType[MAXN];   //[nMC]
   Float_t         mcCalIsoDR03[MAXN];   //[nMC]
   Float_t         mcTrkIsoDR03[MAXN];   //[nMC]
   Float_t         mcCalIsoDR04[MAXN];   //[nMC]
   Float_t         mcTrkIsoDR04[MAXN];   //[nMC]
   Float_t         genMET;
   Float_t         genMETPhi;
   Int_t           nPUInfo;
   Int_t           nPU[3];   //[nPUInfo]
   Int_t           puBX[3];   //[nPUInfo]
   Float_t         MET;
   Float_t         METPhi;
   Float_t         METsumEt;
   Float_t         tcMET;
   Float_t         tcMETPhi;
   Float_t         tcMETsumEt;
   Float_t         tcMETmEtSig;
   Float_t         tcMETSig;
   Float_t         pfMET;
   Float_t         pfMETPhi;
   Float_t         pfMETsumEt;
   Float_t         pfMETmEtSig;
   Float_t         pfMETSig;
   Int_t           nEle;
   Int_t           eleTrg[MAXN][13];   //[nEle]
   Int_t           eleClass[MAXN];   //[nEle]
   Int_t           eleCharge[MAXN];   //[nEle]
   Float_t         eleEn[MAXN];   //[nEle]
   Float_t         eleSCRawEn[MAXN];   //[nEle]
   Float_t         eleSCEn[MAXN];   //[nEle]
   Float_t         eleESEn[MAXN];   //[nEle]
   Float_t         elePt[MAXN];   //[nEle]
   Float_t         eleEta[MAXN];   //[nEle]
   Float_t         elePhi[MAXN];   //[nEle]
   Float_t         eleEtaVtx[MAXN][50];   //[nEle]
   Float_t         elePhiVtx[MAXN][50];   //[nEle]
   Float_t         eleEtVtx[MAXN][50];   //[nEle]
   Float_t         eleSCEta[MAXN];   //[nEle]
   Float_t         eleSCPhi[MAXN];   //[nEle]
   Float_t         eleSCEtaWidth[MAXN];   //[nEle]
   Float_t         eleSCPhiWidth[MAXN];   //[nEle]
   Float_t         eleVtx[MAXN][3];   //[nEle]
   Float_t         eleD0[MAXN];   //[nEle]
   Float_t         eleDz[MAXN];   //[nEle]
   Float_t         eleD0Vtx[MAXN][50];   //[nEle]
   Float_t         eleDzVtx[MAXN][50];   //[nEle]
   Float_t         eleHoverE[MAXN];   //[nEle]
   Float_t         eleEoverP[MAXN];   //[nEle]
   Float_t         elePin[MAXN];   //[nEle]
   Float_t         elePout[MAXN];   //[nEle]
   Float_t         eleBrem[MAXN];   //[nEle]
   Float_t         eledEtaAtVtx[MAXN];   //[nEle]
   Float_t         eledPhiAtVtx[MAXN];   //[nEle]
   Float_t         eleSigmaIEtaIEta[MAXN];   //[nEle]
   Float_t         eleSigmaIEtaIPhi[MAXN];   //[nEle]
   Float_t         eleSigmaIPhiIPhi[MAXN];   //[nEle]
   Float_t         eleEmax[MAXN];   //[nEle]
   Float_t         eleE3x3[MAXN];   //[nEle]
   Float_t         eleE5x5[MAXN];   //[nEle]
   Float_t         eleE2x5Right[MAXN];   //[nEle]
   Float_t         eleE2x5Left[MAXN];   //[nEle]
   Float_t         eleE2x5Top[MAXN];   //[nEle]
   Float_t         eleE2x5Bottom[MAXN];   //[nEle]
   Float_t         eleRegrE[MAXN];   //[nEle]
   Float_t         eleRegrEerr[MAXN];   //[nEle]
   Float_t         eleSeedTime[MAXN];   //[nEle]
   Int_t           eleRecoFlag[MAXN];   //[nEle]
   Int_t           eleGenIndex[MAXN];   //[nEle]
   Int_t           eleGenGMomPID[MAXN];   //[nEle]
   Int_t           eleGenMomPID[MAXN];   //[nEle]
   Float_t         eleGenMomPt[MAXN];   //[nEle]
   Float_t         eleIsoTrkDR03[MAXN];   //[nEle]
   Float_t         eleIsoEcalDR03[MAXN];   //[nEle]
   Float_t         eleIsoHcalDR03[MAXN];   //[nEle]
   Float_t         eleIsoHcalSolidDR03[MAXN];   //[nEle]
   Float_t         eleIsoTrkDR04[MAXN];   //[nEle]
   Float_t         eleIsoEcalDR04[MAXN];   //[nEle]
   Float_t         eleIsoHcalDR04[MAXN];   //[nEle]
   Float_t         eleIsoHcalSolidDR04[MAXN];   //[nEle]
   Int_t           eleMissHits[MAXN];   //[nEle]
   Float_t         eleConvDist[MAXN];   //[nEle]
   Float_t         eleConvDcot[MAXN];   //[nEle]
   Float_t         eleESEffSigmaRR[MAXN][3];   //[nEle]
   Float_t         eleESE1[MAXN][2];   //[nEle]
   Float_t         eleESE3[MAXN][2];   //[nEle]
   Float_t         eleESE5[MAXN][2];   //[nEle]
   Float_t         eleESE7[MAXN][2];   //[nEle]
   Float_t         eleESE11[MAXN][2];   //[nEle]
   Float_t         eleESE21[MAXN][2];   //[nEle]
   Int_t           eleNBC[MAXN];   //[nEle]
   Float_t         eleBrLinear[MAXN];   //[nEle]
   Float_t         eleCetaCorrE[MAXN];   //[nEle]
   Float_t         eleCetaCorrEt[MAXN];   //[nEle]
   Float_t         eleBremCorrE[MAXN];   //[nEle]
   Float_t         eleBremCorrEt[MAXN];   //[nEle]
   Float_t         eleFullCorrE[MAXN];   //[nEle]
   Float_t         eleFullCorrEt[MAXN];   //[nEle]
   Int_t           nPho;
   Int_t           phoTrg[MAXN][8];   //[nPho]
   Int_t           phoTrgFilter[MAXN][30];   //[nPho]
   Bool_t          phoIsPhoton[MAXN];   //[nPho]
   Float_t         phoCaloPos[MAXN][3];   //[nPho]
   Float_t         phoE[MAXN];   //[nPho]
   Float_t         phoEt[MAXN];   //[nPho]
   Float_t         phoEta[MAXN];   //[nPho]
   Float_t         phoVtx[MAXN][3];   //[nPho]
   Float_t         phoPhi[MAXN];   //[nPho]
   Float_t         phoEtVtx[MAXN][50];   //[nPho]
   Float_t         phoEtaVtx[MAXN][50];   //[nPho]
   Float_t         phoPhiVtx[MAXN][50];   //[nPho]
   Float_t         phoR9[MAXN];   //[nPho]
   Float_t         phoCetaCorrE[MAXN];   //[nPho]
   Float_t         phoCetaCorrEt[MAXN];   //[nPho]
   Float_t         phoBremCorrE[MAXN];   //[nPho]
   Float_t         phoBremCorrEt[MAXN];   //[nPho]
   Float_t         phoFullCorrE[MAXN];   //[nPho]
   Float_t         phoFullCorrEt[MAXN];   //[nPho]
   Float_t         phoTrkIsoSolidDR03[MAXN];   //[nPho]
   Float_t         phoTrkIsoHollowDR03[MAXN];   //[nPho]
   Float_t         phoEcalIsoDR03[MAXN];   //[nPho]
   Float_t         phoHcalIsoDR03[MAXN];   //[nPho]
   Float_t         phoHcalIsoSolidDR03[MAXN];   //[nPho]
   Float_t         phoTrkIsoSolidDR04[MAXN];   //[nPho]
   Float_t         phoTrkIsoHollowDR04[MAXN];   //[nPho]
   Float_t         phoTrkIsoHollowNoDzDR03[MAXN];   //[nPho]
   Float_t         phoTrkIsoHollowNoDzDR04[MAXN];   //[nPho]
   Float_t         phoCiCTrkIsoDR03[MAXN][50];   //[nPho]
   Float_t         phoCiCTrkIsoDR04[MAXN][50];   //[nPho]
   Float_t         phoCiCdRtoTrk[MAXN];   //[nPho]
   Int_t           phoCiC4Id[MAXN][50];   //[nPho]
   Int_t           phoCiC6Id[MAXN][50];   //[nPho]
   Int_t           phoCiC6PFId[MAXN][50];   //[nPho]
   Float_t         phoEcalIsoDR04[MAXN];   //[nPho]
   Float_t         phoHcalIsoDR04[MAXN];   //[nPho]
   Float_t         phoHcalIsoSolidDR04[MAXN];   //[nPho]
   Float_t         phoHoverE[MAXN];   //[nPho]
   Float_t         phoSigmaIEtaIEta[MAXN];   //[nPho]
   Float_t         phoSigmaIEtaIPhi[MAXN];   //[nPho]
   Float_t         phoSigmaIPhiIPhi[MAXN];   //[nPho]
   Float_t         phoEmax[MAXN];   //[nPho]
   Float_t         phoEtop[MAXN];   //[nPho]
   Float_t         phoEbottom[MAXN];   //[nPho]
   Float_t         phoEleft[MAXN];   //[nPho]
   Float_t         phoEright[MAXN];   //[nPho]
   Float_t         phoE3x3[MAXN];   //[nPho]
   Float_t         phoE3x1[MAXN];   //[nPho]
   Float_t         phoE1x3[MAXN];   //[nPho]
   Float_t         phoE5x5[MAXN];   //[nPho]
   Float_t         phoE1x5[MAXN];   //[nPho]
   Float_t         phoE2x5Max[MAXN];   //[nPho]
   Float_t         phoE2x5Right[MAXN];   //[nPho]
   Float_t         phoE2x5Left[MAXN];   //[nPho]
   Float_t         phoE2x5Top[MAXN];   //[nPho]
   Float_t         phoE2x5Bottom[MAXN];   //[nPho]
   Float_t         phoPFIsoNeutral[MAXN][8];   //[nPho]
   Float_t         phoPFIsoPhoton[MAXN][8];   //[nPho]
   Float_t         phoPFIsoCharged[MAXN][20][8];   //[nVtx]
   Float_t         phoResiCorrE[MAXN];   //[nPho]
   Float_t         phoResiCorrEsigma[MAXN];   //[nPho]
   Float_t         phoRegrE[MAXN];   //[nPho]
   Float_t         phoRegrEerr[MAXN];   //[nPho]
   Float_t         phoetaC[MAXN];   //[nPho]
   Float_t         phophiC[MAXN];   //[nPho]
   Float_t         phoetaS[MAXN];   //[nPho]
   Float_t         phophiS[MAXN];   //[nPho]
   Float_t         phoetaM[MAXN];   //[nPho]
   Float_t         phophiM[MAXN];   //[nPho]
   Float_t         phoSeedTime[MAXN];   //[nPho]
   Int_t           phoSeedDetId1[MAXN];   //[nPho]
   Int_t           phoSeedDetId2[MAXN];   //[nPho]
   Int_t           phoRecoFlag[MAXN];   //[nPho]
   Int_t           phoPos[MAXN];   //[nPho]
   Int_t           phoGenIndex[MAXN];   //[nPho]
   Int_t           phoGenGMomPID[MAXN];   //[nPho]
   Int_t           phoGenMomPID[MAXN];   //[nPho]
   Float_t         phoGenMomPt[MAXN];   //[nPho]
   Float_t         phoSCE[MAXN];   //[nPho]
   Float_t         phoSCRawE[MAXN];   //[nPho]
   Float_t         phoESEn[MAXN];   //[nPho]
   Float_t         phoSCEt[MAXN];   //[nPho]
   Float_t         phoSCEta[MAXN];   //[nPho]
   Float_t         phoSCPhi[MAXN];   //[nPho]
   Float_t         phoSCEtaWidth[MAXN];   //[nPho]
   Float_t         phoSCPhiWidth[MAXN];   //[nPho]
   Float_t         phoSCBrem[MAXN];   //[nPho]
   Int_t           phoOverlap[MAXN];   //[nPho]
   Int_t           phohasPixelSeed[MAXN];   //[nPho]
   Int_t           phoIsConv[MAXN];   //[nPho]
   Int_t           phoNConv[MAXN];   //[nPho]
   Float_t         phoConvInvMass[MAXN];   //[nPho]
   Float_t         phoConvCotTheta[MAXN];   //[nPho]
   Float_t         phoConvEoverP[MAXN];   //[nPho]
   Float_t         phoConvZofPVfromTrks[MAXN];   //[nPho]
   Float_t         phoConvMinDist[MAXN];   //[nPho]
   Float_t         phoConvdPhiAtVtx[MAXN];   //[nPho]
   Float_t         phoConvdPhiAtCalo[MAXN];   //[nPho]
   Float_t         phoConvdEtaAtCalo[MAXN];   //[nPho]
   Float_t         phoConvTrkd0[MAXN][2];   //[nPho]
   Float_t         phoConvTrkPin[MAXN][2];   //[nPho]
   Float_t         phoConvTrkPout[MAXN][2];   //[nPho]
   Float_t         phoConvTrkdz[MAXN][2];   //[nPho]
   Float_t         phoConvTrkdzErr[MAXN][2];   //[nPho]
   Float_t         phoConvChi2[MAXN];   //[nPho]
   Float_t         phoConvChi2Prob[MAXN];   //[nPho]
   Float_t         phoConvCharge[MAXN][2];   //[nPho]
   Float_t         phoConvValidVtx[MAXN];   //[nPho]
   Float_t         phoConvLikeLihood[MAXN];   //[nPho]
   Float_t         phoConvP4[MAXN][4];   //[nPho]
   Float_t         phoConvVtx[MAXN][3];   //[nPho]
   Float_t         phoConvVtxErr[MAXN][3];   //[nPho]
   Float_t         phoConvPairMomentum[MAXN][3];   //[nPho]
   Float_t         phoConvRefittedMomentum[MAXN][3];   //[nPho]
   Float_t         phoESEffSigmaRR[MAXN][3];   //[nPho]
   Float_t         phoESE1[MAXN][2];   //[nPho]
   Float_t         phoESE3[MAXN][2];   //[nPho]
   Float_t         phoESE5[MAXN][2];   //[nPho]
   Float_t         phoESE7[MAXN][2];   //[nPho]
   Float_t         phoESE11[MAXN][2];   //[nPho]
   Float_t         phoESE21[MAXN][2];   //[nPho]
   Int_t           nMu;
   Int_t           muTrg[MAXN][6];   //[nMu]
   Float_t         muEta[MAXN];   //[nMu]
   Float_t         muPhi[MAXN];   //[nMu]
   Int_t           muCharge[MAXN];   //[nMu]
   Float_t         muPt[MAXN];   //[nMu]
   Float_t         muPz[MAXN];   //[nMu]
   Int_t           muGenIndex[MAXN];   //[nMu]
   Float_t         muIsoTrk[MAXN];   //[nMu]
   Float_t         muIsoCalo[MAXN];   //[nMu]
   Float_t         muIsoEcal[MAXN];   //[nMu]
   Float_t         muIsoHcal[MAXN];   //[nMu]
   Float_t         muChi2NDF[MAXN];   //[nMu]
   Float_t         muEmVeto[MAXN];   //[nMu]
   Float_t         muHadVeto[MAXN];   //[nMu]
   Int_t           muType[MAXN];   //[nMu]
   Bool_t          muID[MAXN][6];   //[nMu]
   Float_t         muD0[MAXN];   //[nMu]
   Float_t         muDz[MAXN];   //[nMu]
   Float_t         muD0Vtx[MAXN][50];   //[nMu]
   Float_t         muDzVtx[MAXN][50];   //[nMu]
   Int_t           muNumberOfValidTrkHits[MAXN];   //[nMu]
   Int_t           muNumberOfValidPixelHits[MAXN];   //[nMu]
   Int_t           muNumberOfValidMuonHits[MAXN];   //[nMu]
   Int_t           muStations[MAXN];   //[nMu]
   Int_t           muChambers[MAXN];   //[nMu]
   Int_t           nPFPho_;
   Float_t         PFPhoE_[kMaxnPFPho];   //[nPFPho_]
   Float_t         PFPhoEt_[kMaxnPFPho];   //[nPFPho_]
   Float_t         PFPhoEta_[kMaxnPFPho];   //[nPFPho_]
   Float_t         PFPhoPhi_[kMaxnPFPho];   //[nPFPho_]
   Float_t         PFPhoHoverE_[kMaxnPFPho];   //[nPFPho_]
   Int_t           nPFEle_;
   Float_t         PFElePt_[kMaxnPFEle];   //[nPFEle_]
   Float_t         PFEleEta_[kMaxnPFEle];   //[nPFEle_]
   Float_t         PFElePhi_[kMaxnPFEle];   //[nPFEle_]
   Float_t         PFEleEn_[kMaxnPFEle];   //[nPFEle_]
   Int_t           PFEleCharge[kMaxnPFEle];   //[nPFEle_]
   Int_t           nPFMu_;
   Int_t           PFmuCharge_[kMaxnPFMu];   //[nPFMu_]
   Float_t         PFmuPhi_[kMaxnPFMu];   //[nPFMu_]
   Float_t         PFmuEta_[kMaxnPFMu];   //[nPFMu_]
   Float_t         PFmuPt_[kMaxnPFMu];   //[nPFMu_]
   Float_t         PFmuPz_[kMaxnPFMu];   //[nPFMu_]
   Int_t           PFmuChambers_[kMaxnPFMu];   //[nPFMu_]
   Float_t         PFmuD0_[kMaxnPFMu];   //[nPFMu_]
   Float_t         PFmuDz_[kMaxnPFMu];   //[nPFMu_]
   Float_t         PFmuChi2NDF_[kMaxnPFMu];   //[nPFMu_]
   Int_t           PFmuNumberOfValidTrkHits_[kMaxnPFMu];   //[nPFMu_]
   Int_t           PFmuNumberOfValidPixelHits_[kMaxnPFMu];   //[nPFMu_]
   Int_t           PFmuNumberOfValidMuonHits_[kMaxnPFMu];   //[nPFMu_]
   Float_t         rho25;
   Int_t           nJet;
   Int_t           jetTrg[MAXN][14];   //[nJet]
   Float_t         jetEn[MAXN];   //[nJet]
   Float_t         jetPt[MAXN];   //[nJet]
   Float_t         jetEta[MAXN];   //[nJet]
   Float_t         jetPhi[MAXN];   //[nJet]
   Float_t         jetEt[MAXN];   //[nJet]
   Float_t         jetRawPt[MAXN];   //[nJet]
   Float_t         jetRawEn[MAXN];   //[nJet]
   Float_t         jetCHF[MAXN];   //[nJet]
   Float_t         jetNHF[MAXN];   //[nJet]
   Float_t         jetCEF[MAXN];   //[nJet]
   Float_t         jetNEF[MAXN];   //[nJet]
   Int_t           jetNCH[MAXN];   //[nJet]
   Float_t         jetHFHAE[MAXN];   //[nJet]
   Float_t         jetHFEME[MAXN];   //[nJet]
   Int_t           jetNConstituents[MAXN];   //[nJet]
   Float_t         jetTrackCountHiEffBJetTags[MAXN];   //[nJet]
   Float_t         jetTrackCountHiPurBJetTags[MAXN];   //[nJet]
   Float_t         jetSimpleSVHiEffBJetTags[MAXN];   //[nJet]
   Float_t         jetSimpleSVHiPurBJetTags[MAXN];   //[nJet]
   Int_t           jetPartonID[MAXN];   //[nJet]
   Int_t           jetGenJetIndex[MAXN];   //[nJet]
   Float_t         jetGenJetEn[MAXN];   //[nJet]
   Float_t         jetGenJetPt[MAXN];   //[nJet]
   Float_t         jetGenJetEta[MAXN];   //[nJet]
   Float_t         jetGenJetPhi[MAXN];   //[nJet]
   Int_t           jetGenPartonID[MAXN];   //[nJet]
   Float_t         jetGenEn[MAXN];   //[nJet]
   Float_t         jetGenPt[MAXN];   //[nJet]
   Float_t         jetGenEta[MAXN];   //[nJet]
   Float_t         jetGenPhi[MAXN];   //[nJet]
   Int_t           nZee;
   Float_t         ZeeMass[MAXN];   //[nZee]
   Float_t         ZeePt[MAXN];   //[nZee]
   Float_t         ZeeEta[MAXN];   //[nZee]
   Float_t         ZeePhi[MAXN];   //[nZee]
   Int_t           ZeeLeg1Index[MAXN];   //[nZee]
   Int_t           ZeeLeg2Index[MAXN];   //[nZee]
   Int_t           nZmumu;
   Float_t         ZmumuMass[MAXN];   //[nZmumu]
   Float_t         ZmumuPt[MAXN];   //[nZmumu]
   Float_t         ZmumuEta[MAXN];   //[nZmumu]
   Float_t         ZmumuPhi[MAXN];   //[nZmumu]
   Int_t           ZmumuLeg1Index[MAXN];   //[nZmumu]
   Int_t           ZmumuLeg2Index[MAXN];   //[nZmumu]
   Int_t           nWenu;
   Float_t         WenuMassTPfMET[MAXN];   //[nWenu]
   Float_t         WenuEtPfMET[MAXN];   //[nWenu]
   Float_t         WenuACopPfMET[MAXN];   //[nWenu]
   Int_t           WenuEleIndex[MAXN];   //[nWenu]
   Int_t           nWmunu;
   Float_t         WmunuMassTPfMET[MAXN];   //[nWmunu]
   Float_t         WmunuEtPfMET[MAXN];   //[nWmunu]
   Float_t         WmunuACopPfMET[MAXN];   //[nWmunu]
   Int_t           WmunuMuIndex[MAXN];   //[nWmunu]
   Int_t           nConv;
   Float_t         convP4[MAXN][4];   //[nConv]
   Float_t         convVtx[MAXN][3];   //[nConv]
   Float_t         convVtxErr[MAXN][3];   //[nConv]
   Float_t         convPairMomentum[MAXN][3];   //[nConv]
   Float_t         convRefittedMomentum[MAXN][3];   //[nConv]
   Int_t           convNTracks[MAXN];   //[nConv]
   Float_t         convPairInvMass[MAXN];   //[nConv]
   Float_t         convPairCotThetaSep[MAXN];   //[nConv]
   Float_t         convEoverP[MAXN];   //[nConv]
   Float_t         convDistOfMinApproach[MAXN];   //[nConv]
   Float_t         convDPhiTrksAtVtx[MAXN];   //[nConv]
   Float_t         convDPhiTrksAtEcal[MAXN];   //[nConv]
   Float_t         convDEtaTrksAtEcal[MAXN];   //[nConv]
   Float_t         convDxy[MAXN];   //[nConv]
   Float_t         convDz[MAXN];   //[nConv]
   Float_t         convLxy[MAXN];   //[nConv]
   Float_t         convLz[MAXN];   //[nConv]
   Float_t         convZofPrimVtxFromTrks[MAXN];   //[nConv]
   Int_t           convNHitsBeforeVtx[MAXN][2];   //[nConv]
   Int_t           convNSharedHits[MAXN];   //[nConv]
   Int_t           convValidVtx[MAXN];   //[nConv]
   Float_t         convMVALikelihood[MAXN];   //[nConv]
   Float_t         convChi2[MAXN];   //[nConv]
   Float_t         convChi2Probability[MAXN];   //[nConv]
   Float_t         convTk1Dz[MAXN];   //[nConv]
   Float_t         convTk2Dz[MAXN];   //[nConv]
   Float_t         convTk1DzErr[MAXN];   //[nConv]
   Float_t         convTk2DzErr[MAXN];   //[nConv]
   Int_t           convCh1Ch2[MAXN];   //[nConv]
   Float_t         convTk1D0[MAXN];   //[nConv]
   Float_t         convTk1Pout[MAXN];   //[nConv]
   Float_t         convTk1Pin[MAXN];   //[nConv]
   Float_t         convTk2D0[MAXN];   //[nConv]
   Float_t         convTk2Pout[MAXN];   //[nConv]
   Float_t         convTk2Pin[MAXN];   //[nConv]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_nHLT;   //!
   TBranch        *b_HLT;   //!
   TBranch        *b_HLTIndex;   //!
   TBranch        *b_bspotPos;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_vtx;   //!
   TBranch        *b_IsVtxGood;   //!
   TBranch        *b_nVtxBS;   //!
   TBranch        *b_vtxbs;   //!
   TBranch        *b_vtxbsPtMod;   //!
   TBranch        *b_vtxbsSumPt2;   //!
   TBranch        *b_vtxbsTkIndex;   //!
   TBranch        *b_vtxbsTkWeight;   //!
   TBranch        *b_nTrk;   //!
   TBranch        *b_trkP;   //!
   TBranch        *b_trkVtx;   //!
   TBranch        *b_trkd0;   //!
   TBranch        *b_trkd0Err;   //!
   TBranch        *b_trkdz;   //!
   TBranch        *b_trkdzErr;   //!
   TBranch        *b_trkPtErr;   //!
   TBranch        *b_trkQuality;   //!
   TBranch        *b_nGoodTrk;   //!
   TBranch        *b_IsTracksGood;   //!
   TBranch        *b_pdf;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_processID;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_mcPID;   //!
   TBranch        *b_mcVtx;   //!
   TBranch        *b_mcPt;   //!
   TBranch        *b_mcMass;   //!
   TBranch        *b_mcEta;   //!
   TBranch        *b_mcPhi;   //!
   TBranch        *b_mcE;   //!
   TBranch        *b_mcEt;   //!
   TBranch        *b_mcGMomPID;   //!
   TBranch        *b_mcMomPID;   //!
   TBranch        *b_mcMomPt;   //!
   TBranch        *b_mcMomMass;   //!
   TBranch        *b_mcMomEta;   //!
   TBranch        *b_mcMomPhi;   //!
   TBranch        *b_mcIndex;   //!
   TBranch        *b_mcDecayType;   //!
   TBranch        *b_mcCalIsoDR03;   //!
   TBranch        *b_mcTrkIsoDR03;   //!
   TBranch        *b_mcCalIsoDR04;   //!
   TBranch        *b_mcTrkIsoDR04;   //!
   TBranch        *b_genMET;   //!
   TBranch        *b_genMETPhi;   //!
   TBranch        *b_nPUInfo;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_puBX;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_METPhi;   //!
   TBranch        *b_METsumEt;   //!
   TBranch        *b_tcMET;   //!
   TBranch        *b_tcMETPhi;   //!
   TBranch        *b_tcMETsumEt;   //!
   TBranch        *b_tcMETmEtSig;   //!
   TBranch        *b_tcMETSig;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_pfMETsumEt;   //!
   TBranch        *b_pfMETmEtSig;   //!
   TBranch        *b_pfMETSig;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_eleTrg;   //!
   TBranch        *b_eleClass;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleEn;   //!
   TBranch        *b_eleSCRawEn;   //!
   TBranch        *b_eleSCEn;   //!
   TBranch        *b_eleESEn;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleEtaVtx;   //!
   TBranch        *b_elePhiVtx;   //!
   TBranch        *b_eleEtVtx;   //!
   TBranch        *b_eleSCEta;   //!
   TBranch        *b_eleSCPhi;   //!
   TBranch        *b_eleSCEtaWidth;   //!
   TBranch        *b_eleSCPhiWidth;   //!
   TBranch        *b_eleVtx;   //!
   TBranch        *b_eleD0;   //!
   TBranch        *b_eleDz;   //!
   TBranch        *b_eleD0Vtx;   //!
   TBranch        *b_eleDzVtx;   //!
   TBranch        *b_eleHoverE;   //!
   TBranch        *b_eleEoverP;   //!
   TBranch        *b_elePin;   //!
   TBranch        *b_elePout;   //!
   TBranch        *b_eleBrem;   //!
   TBranch        *b_eledEtaAtVtx;   //!
   TBranch        *b_eledPhiAtVtx;   //!
   TBranch        *b_eleSigmaIEtaIEta;   //!
   TBranch        *b_eleSigmaIEtaIPhi;   //!
   TBranch        *b_eleSigmaIPhiIPhi;   //!
   TBranch        *b_eleEmax;   //!
   TBranch        *b_eleE3x3;   //!
   TBranch        *b_eleE5x5;   //!
   TBranch        *b_eleE2x5Right;   //!
   TBranch        *b_eleE2x5Left;   //!
   TBranch        *b_eleE2x5Top;   //!
   TBranch        *b_eleE2x5Bottom;   //!
   TBranch        *b_eleRegrE;   //!
   TBranch        *b_eleRegrEerr;   //!
   TBranch        *b_eleSeedTime;   //!
   TBranch        *b_eleRecoFlag;   //!
   TBranch        *b_eleGenIndex;   //!
   TBranch        *b_eleGenGMomPID;   //!
   TBranch        *b_eleGenMomPID;   //!
   TBranch        *b_eleGenMomPt;   //!
   TBranch        *b_eleIsoTrkDR03;   //!
   TBranch        *b_eleIsoEcalDR03;   //!
   TBranch        *b_eleIsoHcalDR03;   //!
   TBranch        *b_eleIsoHcalSolidDR03;   //!
   TBranch        *b_eleIsoTrkDR04;   //!
   TBranch        *b_eleIsoEcalDR04;   //!
   TBranch        *b_eleIsoHcalDR04;   //!
   TBranch        *b_eleIsoHcalSolidDR04;   //!
   TBranch        *b_eleMissHits;   //!
   TBranch        *b_eleConvDist;   //!
   TBranch        *b_eleConvDcot;   //!
   TBranch        *b_eleESEffSigmaRR;   //!
   TBranch        *b_eleESE1;   //!
   TBranch        *b_eleESE3;   //!
   TBranch        *b_eleESE5;   //!
   TBranch        *b_eleESE7;   //!
   TBranch        *b_eleESE11;   //!
   TBranch        *b_eleESE21;   //!
   TBranch        *b_eleNBC;   //!
   TBranch        *b_eleBrLinear;   //!
   TBranch        *b_eleCetaCorrE;   //!
   TBranch        *b_eleCetaCorrEt;   //!
   TBranch        *b_eleBremCorrE;   //!
   TBranch        *b_eleBremCorrEt;   //!
   TBranch        *b_eleFullCorrE;   //!
   TBranch        *b_eleFullCorrEt;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoTrg;   //!
   TBranch        *b_phoTrgFilter;   //!
   TBranch        *b_phoIsPhoton;   //!
   TBranch        *b_phoCaloPos;   //!
   TBranch        *b_phoE;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoVtx;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoEtVtx;   //!
   TBranch        *b_phoEtaVtx;   //!
   TBranch        *b_phoPhiVtx;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoCetaCorrE;   //!
   TBranch        *b_phoCetaCorrEt;   //!
   TBranch        *b_phoBremCorrE;   //!
   TBranch        *b_phoBremCorrEt;   //!
   TBranch        *b_phoFullCorrE;   //!
   TBranch        *b_phoFullCorrEt;   //!
   TBranch        *b_phoTrkIsoSolidDR03;   //!
   TBranch        *b_phoTrkIsoHollowDR03;   //!
   TBranch        *b_phoEcalIsoDR03;   //!
   TBranch        *b_phoHcalIsoDR03;   //!
   TBranch        *b_phoHcalIsoSolidDR03;   //!
   TBranch        *b_phoTrkIsoSolidDR04;   //!
   TBranch        *b_phoTrkIsoHollowDR04;   //!
   TBranch        *b_phoTrkIsoHollowNoDzDR03;   //!
   TBranch        *b_phoTrkIsoHollowNoDzDR04;   //!
   TBranch        *b_phoCiCTrkIsoDR03;   //!
   TBranch        *b_phoCiCTrkIsoDR04;   //!
   TBranch        *b_phoCiCdRtoTrk;   //!
   TBranch        *b_phoCiC4Id;   //!
   TBranch        *b_phoCiC6Id;   //!
   TBranch        *b_phoCiC6PFId;   //!
   TBranch        *b_phoEcalIsoDR04;   //!
   TBranch        *b_phoHcalIsoDR04;   //!
   TBranch        *b_phoHcalIsoSolidDR04;   //!
   TBranch        *b_phoHoverE;   //!
   TBranch        *b_phoSigmaIEtaIEta;   //!
   TBranch        *b_phoSigmaIEtaIPhi;   //!
   TBranch        *b_phoSigmaIPhiIPhi;   //!
   TBranch        *b_phoEmax;   //!
   TBranch        *b_phoEtop;   //!
   TBranch        *b_phoEbottom;   //!
   TBranch        *b_phoEleft;   //!
   TBranch        *b_phoEright;   //!
   TBranch        *b_phoE3x3;   //!
   TBranch        *b_phoE3x1;   //!
   TBranch        *b_phoE1x3;   //!
   TBranch        *b_phoE5x5;   //!
   TBranch        *b_phoE1x5;   //!
   TBranch        *b_phoE2x5Max;   //!
   TBranch        *b_phoE2x5Right;   //!
   TBranch        *b_phoE2x5Left;   //!
   TBranch        *b_phoE2x5Top;   //!
   TBranch        *b_phoE2x5Bottom;   //!
   TBranch        *b_phoPFIsoNeutral;   //!
   TBranch        *b_phoPFIsoPhoton;   //!
   TBranch        *b_phoPFIsoCharged;   //!
   TBranch        *b_phoResiCorrE;   //!
   TBranch        *b_phoResiCorrEsigma;   //!
   TBranch        *b_phoRegrE;   //!
   TBranch        *b_phoRegrEerr;   //!
   TBranch        *b_phoetaC;   //!
   TBranch        *b_phophiC;   //!
   TBranch        *b_phoetaS;   //!
   TBranch        *b_phophiS;   //!
   TBranch        *b_phoetaM;   //!
   TBranch        *b_phophiM;   //!
   TBranch        *b_phoSeedTime;   //!
   TBranch        *b_phoSeedDetId1;   //!
   TBranch        *b_phoSeedDetId2;   //!
   TBranch        *b_phoRecoFlag;   //!
   TBranch        *b_phoPos;   //!
   TBranch        *b_phoGenIndex;   //!
   TBranch        *b_phoGenGMomPID;   //!
   TBranch        *b_phoGenMomPID;   //!
   TBranch        *b_phoGenMomPt;   //!
   TBranch        *b_phoSCE;   //!
   TBranch        *b_phoSCRawE;   //!
   TBranch        *b_phoESEn;   //!
   TBranch        *b_phoSCEt;   //!
   TBranch        *b_phoSCEta;   //!
   TBranch        *b_phoSCPhi;   //!
   TBranch        *b_phoSCEtaWidth;   //!
   TBranch        *b_phoSCPhiWidth;   //!
   TBranch        *b_phoSCBrem;   //!
   TBranch        *b_phoOverlap;   //!
   TBranch        *b_phohasPixelSeed;   //!
   TBranch        *b_phoIsConv;   //!
   TBranch        *b_phoNConv;   //!
   TBranch        *b_phoConvInvMass;   //!
   TBranch        *b_phoConvCotTheta;   //!
   TBranch        *b_phoConvEoverP;   //!
   TBranch        *b_phoConvZofPVfromTrks;   //!
   TBranch        *b_phoConvMinDist;   //!
   TBranch        *b_phoConvdPhiAtVtx;   //!
   TBranch        *b_phoConvdPhiAtCalo;   //!
   TBranch        *b_phoConvdEtaAtCalo;   //!
   TBranch        *b_phoConvTrkd0;   //!
   TBranch        *b_phoConvTrkPin;   //!
   TBranch        *b_phoConvTrkPout;   //!
   TBranch        *b_phoConvTrkdz;   //!
   TBranch        *b_phoConvTrkdzErr;   //!
   TBranch        *b_phoConvChi2;   //!
   TBranch        *b_phoConvChi2Prob;   //!
   TBranch        *b_phoConvCharge;   //!
   TBranch        *b_phoConvValidVtx;   //!
   TBranch        *b_phoConvLikeLihood;   //!
   TBranch        *b_phoConvP4;   //!
   TBranch        *b_phoConvVtx;   //!
   TBranch        *b_phoConvVtxErr;   //!
   TBranch        *b_phoConvPairMomentum;   //!
   TBranch        *b_phoConvRefittedMomentum;   //!
   TBranch        *b_phoESEffSigmaRR;   //!
   TBranch        *b_phoESE1;   //!
   TBranch        *b_phoESE3;   //!
   TBranch        *b_phoESE5;   //!
   TBranch        *b_phoESE7;   //!
   TBranch        *b_phoESE11;   //!
   TBranch        *b_phoESE21;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_muTrg;   //!
   TBranch        *b_muEta;   //!
   TBranch        *b_muPhi;   //!
   TBranch        *b_muCharge;   //!
   TBranch        *b_muPt;   //!
   TBranch        *b_muPz;   //!
   TBranch        *b_muGenIndex;   //!
   TBranch        *b_muIsoTrk;   //!
   TBranch        *b_muIsoCalo;   //!
   TBranch        *b_muIsoEcal;   //!
   TBranch        *b_muIsoHcal;   //!
   TBranch        *b_muChi2NDF;   //!
   TBranch        *b_muEmVeto;   //!
   TBranch        *b_muHadVeto;   //!
   TBranch        *b_muType;   //!
   TBranch        *b_muID;   //!
   TBranch        *b_muD0;   //!
   TBranch        *b_muDz;   //!
   TBranch        *b_muD0Vtx;   //!
   TBranch        *b_muDzVtx;   //!
   TBranch        *b_muNumberOfValidTrkHits;   //!
   TBranch        *b_muNumberOfValidPixelHits;   //!
   TBranch        *b_muNumberOfValidMuonHits;   //!
   TBranch        *b_muStations;   //!
   TBranch        *b_muChambers;   //!
   TBranch        *b_nPFPho_;   //!
   TBranch        *b_PFPhoE_;   //!
   TBranch        *b_PFPhoEt_;   //!
   TBranch        *b_PFPhoEta_;   //!
   TBranch        *b_PFPhoPhi_;   //!
   TBranch        *b_PFPhoHoverE_;   //!
   TBranch        *b_nPFEle_;   //!
   TBranch        *b_PFElePt_;   //!
   TBranch        *b_PFEleEta_;   //!
   TBranch        *b_PFElePhi_;   //!
   TBranch        *b_PFEleEn_;   //!
   TBranch        *b_PFEleCharge;   //!
   TBranch        *b_nPFMu_;   //!
   TBranch        *b_PFmuCharge_;   //!
   TBranch        *b_PFmuPhi_;   //!
   TBranch        *b_PFmuEta_;   //!
   TBranch        *b_PFmuPt_;   //!
   TBranch        *b_PFmuPz_;   //!
   TBranch        *b_PFmuChambers_;   //!
   TBranch        *b_PFmuD0_;   //!
   TBranch        *b_PFmuDz_;   //!
   TBranch        *b_PFmuChi2NDF_;   //!
   TBranch        *b_PFmuNumberOfValidTrkHits_;   //!
   TBranch        *b_PFmuNumberOfValidPixelHits_;   //!
   TBranch        *b_PFmuNumberOfValidMuonHits_;   //!
   TBranch        *b_rho25;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_jetTrg;   //!
   TBranch        *b_jetEn;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetEt;   //!
   TBranch        *b_jetRawPt;   //!
   TBranch        *b_jetRawEn;   //!
   TBranch        *b_jetCHF;   //!
   TBranch        *b_jetNHF;   //!
   TBranch        *b_jetCEF;   //!
   TBranch        *b_jetNEF;   //!
   TBranch        *b_jetNCH;   //!
   TBranch        *b_jetHFHAE;   //!
   TBranch        *b_jetHFEME;   //!
   TBranch        *b_jetNConstituents;   //!
   TBranch        *b_jetTrackCountHiEffBJetTags;   //!
   TBranch        *b_jetTrackCountHiPurBJetTags;   //!
   TBranch        *b_jetSimpleSVHiEffBJetTags;   //!
   TBranch        *b_jetSimpleSVHiPurBJetTags;   //!
   TBranch        *b_jetPartonID;   //!
   TBranch        *b_jetGenJetIndex;   //!
   TBranch        *b_jetGenJetEn;   //!
   TBranch        *b_jetGenJetPt;   //!
   TBranch        *b_jetGenJetEta;   //!
   TBranch        *b_jetGenJetPhi;   //!
   TBranch        *b_jetGenPartonID;   //!
   TBranch        *b_jetGenEn;   //!
   TBranch        *b_jetGenPt;   //!
   TBranch        *b_jetGenEta;   //!
   TBranch        *b_jetGenPhi;   //!
   TBranch        *b_nZee;   //!
   TBranch        *b_ZeeMass;   //!
   TBranch        *b_ZeePt;   //!
   TBranch        *b_ZeeEta;   //!
   TBranch        *b_ZeePhi;   //!
   TBranch        *b_ZeeLeg1Index;   //!
   TBranch        *b_ZeeLeg2Index;   //!
   TBranch        *b_nZmumu;   //!
   TBranch        *b_ZmumuMass;   //!
   TBranch        *b_ZmumuPt;   //!
   TBranch        *b_ZmumuEta;   //!
   TBranch        *b_ZmumuPhi;   //!
   TBranch        *b_ZmumuLeg1Index;   //!
   TBranch        *b_ZmumuLeg2Index;   //!
   TBranch        *b_nWenu;   //!
   TBranch        *b_WenuMassTPfMET;   //!
   TBranch        *b_WenuEtPfMET;   //!
   TBranch        *b_WenuACopPfMET;   //!
   TBranch        *b_WenuEleIndex;   //!
   TBranch        *b_nWmunu;   //!
   TBranch        *b_WmunuMassTPfMET;   //!
   TBranch        *b_WmunuEtPfMET;   //!
   TBranch        *b_WmunuACopPfMET;   //!
   TBranch        *b_WmunuMuIndex;   //!
   TBranch        *b_nConv;   //!
   TBranch        *b_convP4;   //!
   TBranch        *b_convVtx;   //!
   TBranch        *b_convVtxErr;   //!
   TBranch        *b_convPairMomentum;   //!
   TBranch        *b_convRefittedMomentum;   //!
   TBranch        *b_convNTracks;   //!
   TBranch        *b_convPairInvMass;   //!
   TBranch        *b_convPairCotThetaSep;   //!
   TBranch        *b_convEoverP;   //!
   TBranch        *b_convDistOfMinApproach;   //!
   TBranch        *b_convDPhiTrksAtVtx;   //!
   TBranch        *b_convDPhiTrksAtEcal;   //!
   TBranch        *b_convDEtaTrksAtEcal;   //!
   TBranch        *b_convDxy;   //!
   TBranch        *b_convDz;   //!
   TBranch        *b_convLxy;   //!
   TBranch        *b_convLz;   //!
   TBranch        *b_convZofPrimVtxFromTrks;   //!
   TBranch        *b_convNHitsBeforeVtx;   //!
   TBranch        *b_convNSharedHits;   //!
   TBranch        *b_convValidVtx;   //!
   TBranch        *b_convMVALikelihood;   //!
   TBranch        *b_convChi2;   //!
   TBranch        *b_convChi2Probability;   //!
   TBranch        *b_convTk1Dz;   //!
   TBranch        *b_convTk2Dz;   //!
   TBranch        *b_convTk1DzErr;   //!
   TBranch        *b_convTk2DzErr;   //!
   TBranch        *b_convCh1Ch2;   //!
   TBranch        *b_convTk1D0;   //!
   TBranch        *b_convTk1Pout;   //!
   TBranch        *b_convTk1Pin;   //!
   TBranch        *b_convTk2D0;   //!
   TBranch        *b_convTk2Pout;   //!
   TBranch        *b_convTk2Pin;   //!

   gg_angularmc_eff(std::string filename,TTree *tree=0);
   virtual ~gg_angularmc_eff();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(bool applyCOMCut=false, bool applyPileupCorr=true);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Int_t    phoDecCode(Long64_t entry, Int_t ipho);
   virtual Bool_t   isGoodPho(Long64_t entry, Int_t ipho, bool applyPileupCorr=false);
   virtual Double_t phoEcalIso(Long64_t entry, Int_t ipho, bool applyPileupCorr=false);
   virtual Double_t phoHcalIso(Long64_t entry, Int_t ipho, bool applyPileupCorr=false);
   virtual Double_t phoTrkIso(Long64_t entry, Int_t ipho, bool applyPileupCorr=false);
   virtual Bool_t   isGoodLooseJet(Long64_t entry, Int_t ijet);
   virtual Bool_t   isGoodMediumJet(Long64_t entry, Int_t ijet);
   virtual Bool_t   isGoodTightJet(Long64_t entry, Int_t ijet);
   virtual Bool_t   isFidJet (Long64_t entry, Int_t ijet);
   virtual Bool_t   isFidPho (Long64_t entry, Int_t ipho);
   std::string _inputFileName;

};

#endif

#ifdef gg_angularmc_eff_cxx
gg_angularmc_eff::gg_angularmc_eff(std::string filename, TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename.data());
     if (!f) {
       f = new TFile(filename.data());
       f->cd(Form("%s:/ggNtuplizer",filename.data()));
      }
      tree = (TTree*)gDirectory->Get("EventTree");

   }
   Init(tree);
   _inputFileName = filename;
}

gg_angularmc_eff::~gg_angularmc_eff()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gg_angularmc_eff::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gg_angularmc_eff::LoadTree(Long64_t entry)
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

void gg_angularmc_eff::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   vtxbsTkIndex = 0;
   vtxbsTkWeight = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

//    fChain->SetBranchAddress("run", &run, &b_run);
//    fChain->SetBranchAddress("event", &event, &b_event);
//    fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
//    fChain->SetBranchAddress("isData", &isData, &b_isData);
//    fChain->SetBranchAddress("nHLT", &nHLT, &b_nHLT);
//    fChain->SetBranchAddress("HLT", HLT, &b_HLT);
//    fChain->SetBranchAddress("HLTIndex", HLTIndex, &b_HLTIndex);
//    fChain->SetBranchAddress("bspotPos", bspotPos, &b_bspotPos);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("vtx", vtx, &b_vtx);
   fChain->SetBranchAddress("IsVtxGood", &IsVtxGood, &b_IsVtxGood);
//    fChain->SetBranchAddress("nVtxBS", &nVtxBS, &b_nVtxBS);
//    fChain->SetBranchAddress("vtxbs", vtxbs, &b_vtxbs);
//    fChain->SetBranchAddress("vtxbsPtMod", vtxbsPtMod, &b_vtxbsPtMod);
//    fChain->SetBranchAddress("vtxbsSumPt2", vtxbsSumPt2, &b_vtxbsSumPt2);
//    fChain->SetBranchAddress("vtxbsTkIndex", &vtxbsTkIndex, &b_vtxbsTkIndex);
//    fChain->SetBranchAddress("vtxbsTkWeight", &vtxbsTkWeight, &b_vtxbsTkWeight);
//    fChain->SetBranchAddress("nTrk", &nTrk, &b_nTrk);
//    fChain->SetBranchAddress("trkP", trkP, &b_trkP);
//    fChain->SetBranchAddress("trkVtx", trkVtx, &b_trkVtx);
//    fChain->SetBranchAddress("trkd0", trkd0, &b_trkd0);
//    fChain->SetBranchAddress("trkd0Err", trkd0Err, &b_trkd0Err);
//    fChain->SetBranchAddress("trkdz", trkdz, &b_trkdz);
//    fChain->SetBranchAddress("trkdzErr", trkdzErr, &b_trkdzErr);
//    fChain->SetBranchAddress("trkPtErr", trkPtErr, &b_trkPtErr);
//    fChain->SetBranchAddress("trkQuality", trkQuality, &b_trkQuality);
//    fChain->SetBranchAddress("nGoodTrk", &nGoodTrk, &b_nGoodTrk);
//    fChain->SetBranchAddress("IsTracksGood", &IsTracksGood, &b_IsTracksGood);
//    fChain->SetBranchAddress("pdf", pdf, &b_pdf);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
//    fChain->SetBranchAddress("processID", &processID, &b_processID);
//    fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
//    fChain->SetBranchAddress("mcPID", mcPID, &b_mcPID);
//    fChain->SetBranchAddress("mcVtx", mcVtx, &b_mcVtx);
//    fChain->SetBranchAddress("mcPt", mcPt, &b_mcPt);
//    fChain->SetBranchAddress("mcMass", mcMass, &b_mcMass);
//    fChain->SetBranchAddress("mcEta", mcEta, &b_mcEta);
//    fChain->SetBranchAddress("mcPhi", mcPhi, &b_mcPhi);
//    fChain->SetBranchAddress("mcE", mcE, &b_mcE);
//    fChain->SetBranchAddress("mcEt", mcEt, &b_mcEt);
//    fChain->SetBranchAddress("mcGMomPID", mcGMomPID, &b_mcGMomPID);
//    fChain->SetBranchAddress("mcMomPID", mcMomPID, &b_mcMomPID);
//    fChain->SetBranchAddress("mcMomPt", mcMomPt, &b_mcMomPt);
//    fChain->SetBranchAddress("mcMomMass", mcMomMass, &b_mcMomMass);
//    fChain->SetBranchAddress("mcMomEta", mcMomEta, &b_mcMomEta);
//    fChain->SetBranchAddress("mcMomPhi", mcMomPhi, &b_mcMomPhi);
//    fChain->SetBranchAddress("mcIndex", mcIndex, &b_mcIndex);
//    fChain->SetBranchAddress("mcDecayType", mcDecayType, &b_mcDecayType);
//    fChain->SetBranchAddress("mcCalIsoDR03", mcCalIsoDR03, &b_mcCalIsoDR03);
//    fChain->SetBranchAddress("mcTrkIsoDR03", mcTrkIsoDR03, &b_mcTrkIsoDR03);
//    fChain->SetBranchAddress("mcCalIsoDR04", mcCalIsoDR04, &b_mcCalIsoDR04);
//    fChain->SetBranchAddress("mcTrkIsoDR04", mcTrkIsoDR04, &b_mcTrkIsoDR04);
//    fChain->SetBranchAddress("genMET", &genMET, &b_genMET);
//    fChain->SetBranchAddress("genMETPhi", &genMETPhi, &b_genMETPhi);
//    fChain->SetBranchAddress("nPUInfo", &nPUInfo, &b_nPUInfo);
//    fChain->SetBranchAddress("nPU", nPU, &b_nPU);
//    fChain->SetBranchAddress("puBX", puBX, &b_puBX);
//    fChain->SetBranchAddress("MET", &MET, &b_MET);
//    fChain->SetBranchAddress("METPhi", &METPhi, &b_METPhi);
//    fChain->SetBranchAddress("METsumEt", &METsumEt, &b_METsumEt);
//    fChain->SetBranchAddress("tcMET", &tcMET, &b_tcMET);
//    fChain->SetBranchAddress("tcMETPhi", &tcMETPhi, &b_tcMETPhi);
//    fChain->SetBranchAddress("tcMETsumEt", &tcMETsumEt, &b_tcMETsumEt);
//    fChain->SetBranchAddress("tcMETmEtSig", &tcMETmEtSig, &b_tcMETmEtSig);
//    fChain->SetBranchAddress("tcMETSig", &tcMETSig, &b_tcMETSig);
//    fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
//    fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
//    fChain->SetBranchAddress("pfMETsumEt", &pfMETsumEt, &b_pfMETsumEt);
//    fChain->SetBranchAddress("pfMETmEtSig", &pfMETmEtSig, &b_pfMETmEtSig);
//    fChain->SetBranchAddress("pfMETSig", &pfMETSig, &b_pfMETSig);
//    fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
//    fChain->SetBranchAddress("eleTrg", eleTrg, &b_eleTrg);
//    fChain->SetBranchAddress("eleClass", eleClass, &b_eleClass);
//    fChain->SetBranchAddress("eleCharge", eleCharge, &b_eleCharge);
//    fChain->SetBranchAddress("eleEn", eleEn, &b_eleEn);
//    fChain->SetBranchAddress("eleSCRawEn", eleSCRawEn, &b_eleSCRawEn);
//    fChain->SetBranchAddress("eleSCEn", eleSCEn, &b_eleSCEn);
//    fChain->SetBranchAddress("eleESEn", eleESEn, &b_eleESEn);
//    fChain->SetBranchAddress("elePt", elePt, &b_elePt);
//    fChain->SetBranchAddress("eleEta", eleEta, &b_eleEta);
//    fChain->SetBranchAddress("elePhi", elePhi, &b_elePhi);
//    fChain->SetBranchAddress("eleEtaVtx", eleEtaVtx, &b_eleEtaVtx);
//    fChain->SetBranchAddress("elePhiVtx", elePhiVtx, &b_elePhiVtx);
//    fChain->SetBranchAddress("eleEtVtx", eleEtVtx, &b_eleEtVtx);
//    fChain->SetBranchAddress("eleSCEta", eleSCEta, &b_eleSCEta);
//    fChain->SetBranchAddress("eleSCPhi", eleSCPhi, &b_eleSCPhi);
//    fChain->SetBranchAddress("eleSCEtaWidth", eleSCEtaWidth, &b_eleSCEtaWidth);
//    fChain->SetBranchAddress("eleSCPhiWidth", eleSCPhiWidth, &b_eleSCPhiWidth);
//    fChain->SetBranchAddress("eleVtx", eleVtx, &b_eleVtx);
//    fChain->SetBranchAddress("eleD0", eleD0, &b_eleD0);
//    fChain->SetBranchAddress("eleDz", eleDz, &b_eleDz);
//    fChain->SetBranchAddress("eleD0Vtx", eleD0Vtx, &b_eleD0Vtx);
//    fChain->SetBranchAddress("eleDzVtx", eleDzVtx, &b_eleDzVtx);
//    fChain->SetBranchAddress("eleHoverE", eleHoverE, &b_eleHoverE);
//    fChain->SetBranchAddress("eleEoverP", eleEoverP, &b_eleEoverP);
//    fChain->SetBranchAddress("elePin", elePin, &b_elePin);
//    fChain->SetBranchAddress("elePout", elePout, &b_elePout);
//    fChain->SetBranchAddress("eleBrem", eleBrem, &b_eleBrem);
//    fChain->SetBranchAddress("eledEtaAtVtx", eledEtaAtVtx, &b_eledEtaAtVtx);
//    fChain->SetBranchAddress("eledPhiAtVtx", eledPhiAtVtx, &b_eledPhiAtVtx);
//    fChain->SetBranchAddress("eleSigmaIEtaIEta", eleSigmaIEtaIEta, &b_eleSigmaIEtaIEta);
//    fChain->SetBranchAddress("eleSigmaIEtaIPhi", eleSigmaIEtaIPhi, &b_eleSigmaIEtaIPhi);
//    fChain->SetBranchAddress("eleSigmaIPhiIPhi", eleSigmaIPhiIPhi, &b_eleSigmaIPhiIPhi);
//    fChain->SetBranchAddress("eleEmax", eleEmax, &b_eleEmax);
//    fChain->SetBranchAddress("eleE3x3", eleE3x3, &b_eleE3x3);
//    fChain->SetBranchAddress("eleE5x5", eleE5x5, &b_eleE5x5);
//    fChain->SetBranchAddress("eleE2x5Right", eleE2x5Right, &b_eleE2x5Right);
//    fChain->SetBranchAddress("eleE2x5Left", eleE2x5Left, &b_eleE2x5Left);
//    fChain->SetBranchAddress("eleE2x5Top", eleE2x5Top, &b_eleE2x5Top);
//    fChain->SetBranchAddress("eleE2x5Bottom", eleE2x5Bottom, &b_eleE2x5Bottom);
//    fChain->SetBranchAddress("eleRegrE", eleRegrE, &b_eleRegrE);
//    fChain->SetBranchAddress("eleRegrEerr", eleRegrEerr, &b_eleRegrEerr);
//    fChain->SetBranchAddress("eleSeedTime", eleSeedTime, &b_eleSeedTime);
//    fChain->SetBranchAddress("eleRecoFlag", eleRecoFlag, &b_eleRecoFlag);
//    fChain->SetBranchAddress("eleGenIndex", eleGenIndex, &b_eleGenIndex);
//    fChain->SetBranchAddress("eleGenGMomPID", eleGenGMomPID, &b_eleGenGMomPID);
//    fChain->SetBranchAddress("eleGenMomPID", eleGenMomPID, &b_eleGenMomPID);
//    fChain->SetBranchAddress("eleGenMomPt", eleGenMomPt, &b_eleGenMomPt);
//    fChain->SetBranchAddress("eleIsoTrkDR03", eleIsoTrkDR03, &b_eleIsoTrkDR03);
//    fChain->SetBranchAddress("eleIsoEcalDR03", eleIsoEcalDR03, &b_eleIsoEcalDR03);
//    fChain->SetBranchAddress("eleIsoHcalDR03", eleIsoHcalDR03, &b_eleIsoHcalDR03);
//    fChain->SetBranchAddress("eleIsoHcalSolidDR03", eleIsoHcalSolidDR03, &b_eleIsoHcalSolidDR03);
//    fChain->SetBranchAddress("eleIsoTrkDR04", eleIsoTrkDR04, &b_eleIsoTrkDR04);
//    fChain->SetBranchAddress("eleIsoEcalDR04", eleIsoEcalDR04, &b_eleIsoEcalDR04);
//    fChain->SetBranchAddress("eleIsoHcalDR04", eleIsoHcalDR04, &b_eleIsoHcalDR04);
//    fChain->SetBranchAddress("eleIsoHcalSolidDR04", eleIsoHcalSolidDR04, &b_eleIsoHcalSolidDR04);
//    fChain->SetBranchAddress("eleMissHits", eleMissHits, &b_eleMissHits);
//    fChain->SetBranchAddress("eleConvDist", eleConvDist, &b_eleConvDist);
//    fChain->SetBranchAddress("eleConvDcot", eleConvDcot, &b_eleConvDcot);
//    fChain->SetBranchAddress("eleESEffSigmaRR", eleESEffSigmaRR, &b_eleESEffSigmaRR);
//    fChain->SetBranchAddress("eleESE1", eleESE1, &b_eleESE1);
//    fChain->SetBranchAddress("eleESE3", eleESE3, &b_eleESE3);
//    fChain->SetBranchAddress("eleESE5", eleESE5, &b_eleESE5);
//    fChain->SetBranchAddress("eleESE7", eleESE7, &b_eleESE7);
//    fChain->SetBranchAddress("eleESE11", eleESE11, &b_eleESE11);
//    fChain->SetBranchAddress("eleESE21", eleESE21, &b_eleESE21);
//    fChain->SetBranchAddress("eleNBC", eleNBC, &b_eleNBC);
//    fChain->SetBranchAddress("eleBrLinear", eleBrLinear, &b_eleBrLinear);
//    fChain->SetBranchAddress("eleCetaCorrE", eleCetaCorrE, &b_eleCetaCorrE);
//    fChain->SetBranchAddress("eleCetaCorrEt", eleCetaCorrEt, &b_eleCetaCorrEt);
//    fChain->SetBranchAddress("eleBremCorrE", eleBremCorrE, &b_eleBremCorrE);
//    fChain->SetBranchAddress("eleBremCorrEt", eleBremCorrEt, &b_eleBremCorrEt);
//    fChain->SetBranchAddress("eleFullCorrE", eleFullCorrE, &b_eleFullCorrE);
//    fChain->SetBranchAddress("eleFullCorrEt", eleFullCorrEt, &b_eleFullCorrEt);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
//    fChain->SetBranchAddress("phoTrg", phoTrg, &b_phoTrg);
//    fChain->SetBranchAddress("phoTrgFilter", phoTrgFilter, &b_phoTrgFilter);
//    fChain->SetBranchAddress("phoIsPhoton", phoIsPhoton, &b_phoIsPhoton);
//    fChain->SetBranchAddress("phoCaloPos", phoCaloPos, &b_phoCaloPos);
   fChain->SetBranchAddress("phoE", phoE, &b_phoE);
   fChain->SetBranchAddress("phoEt", phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoEta", phoEta, &b_phoEta);
//    fChain->SetBranchAddress("phoVtx", phoVtx, &b_phoVtx);
   fChain->SetBranchAddress("phoPhi", phoPhi, &b_phoPhi);
//    fChain->SetBranchAddress("phoEtVtx", phoEtVtx, &b_phoEtVtx);
//    fChain->SetBranchAddress("phoEtaVtx", phoEtaVtx, &b_phoEtaVtx);
//    fChain->SetBranchAddress("phoPhiVtx", phoPhiVtx, &b_phoPhiVtx);
//    fChain->SetBranchAddress("phoR9", phoR9, &b_phoR9);
//    fChain->SetBranchAddress("phoCetaCorrE", phoCetaCorrE, &b_phoCetaCorrE);
//    fChain->SetBranchAddress("phoCetaCorrEt", phoCetaCorrEt, &b_phoCetaCorrEt);
//    fChain->SetBranchAddress("phoBremCorrE", phoBremCorrE, &b_phoBremCorrE);
//    fChain->SetBranchAddress("phoBremCorrEt", phoBremCorrEt, &b_phoBremCorrEt);
//    fChain->SetBranchAddress("phoFullCorrE", phoFullCorrE, &b_phoFullCorrE);
//    fChain->SetBranchAddress("phoFullCorrEt", phoFullCorrEt, &b_phoFullCorrEt);
//    fChain->SetBranchAddress("phoTrkIsoSolidDR03", phoTrkIsoSolidDR03, &b_phoTrkIsoSolidDR03);
//    fChain->SetBranchAddress("phoTrkIsoHollowDR03", phoTrkIsoHollowDR03, &b_phoTrkIsoHollowDR03);
//    fChain->SetBranchAddress("phoEcalIsoDR03", phoEcalIsoDR03, &b_phoEcalIsoDR03);
//    fChain->SetBranchAddress("phoHcalIsoDR03", phoHcalIsoDR03, &b_phoHcalIsoDR03);
//    fChain->SetBranchAddress("phoHcalIsoSolidDR03", phoHcalIsoSolidDR03, &b_phoHcalIsoSolidDR03);
//    fChain->SetBranchAddress("phoTrkIsoSolidDR04", phoTrkIsoSolidDR04, &b_phoTrkIsoSolidDR04);
   fChain->SetBranchAddress("phoTrkIsoHollowDR04", phoTrkIsoHollowDR04, &b_phoTrkIsoHollowDR04);
//    fChain->SetBranchAddress("phoTrkIsoHollowNoDzDR03", phoTrkIsoHollowNoDzDR03, &b_phoTrkIsoHollowNoDzDR03);
//    fChain->SetBranchAddress("phoTrkIsoHollowNoDzDR04", phoTrkIsoHollowNoDzDR04, &b_phoTrkIsoHollowNoDzDR04);
//    fChain->SetBranchAddress("phoCiCTrkIsoDR03", phoCiCTrkIsoDR03, &b_phoCiCTrkIsoDR03);
//    fChain->SetBranchAddress("phoCiCTrkIsoDR04", phoCiCTrkIsoDR04, &b_phoCiCTrkIsoDR04);
//    fChain->SetBranchAddress("phoCiCdRtoTrk", phoCiCdRtoTrk, &b_phoCiCdRtoTrk);
//    fChain->SetBranchAddress("phoCiC4Id", phoCiC4Id, &b_phoCiC4Id);
//    fChain->SetBranchAddress("phoCiC6Id", phoCiC6Id, &b_phoCiC6Id);
//    fChain->SetBranchAddress("phoCiC6PFId", phoCiC6PFId, &b_phoCiC6PFId);
   fChain->SetBranchAddress("phoEcalIsoDR04", phoEcalIsoDR04, &b_phoEcalIsoDR04);
   fChain->SetBranchAddress("phoHcalIsoDR04", phoHcalIsoDR04, &b_phoHcalIsoDR04);
//    fChain->SetBranchAddress("phoHcalIsoSolidDR04", phoHcalIsoSolidDR04, &b_phoHcalIsoSolidDR04);
   fChain->SetBranchAddress("phoHoverE", phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phoSigmaIEtaIEta", phoSigmaIEtaIEta, &b_phoSigmaIEtaIEta);
//    fChain->SetBranchAddress("phoSigmaIEtaIPhi", phoSigmaIEtaIPhi, &b_phoSigmaIEtaIPhi);
//    fChain->SetBranchAddress("phoSigmaIPhiIPhi", phoSigmaIPhiIPhi, &b_phoSigmaIPhiIPhi);
// //    fChain->SetBranchAddress("phoEmax", phoEmax, &b_phoEmax);
//    fChain->SetBranchAddress("phoEtop", phoEtop, &b_phoEtop);
//    fChain->SetBranchAddress("phoEbottom", phoEbottom, &b_phoEbottom);
//    fChain->SetBranchAddress("phoEleft", phoEleft, &b_phoEleft);
//    fChain->SetBranchAddress("phoEright", phoEright, &b_phoEright);
//    fChain->SetBranchAddress("phoE3x3", phoE3x3, &b_phoE3x3);
//    fChain->SetBranchAddress("phoE3x1", phoE3x1, &b_phoE3x1);
//    fChain->SetBranchAddress("phoE1x3", phoE1x3, &b_phoE1x3);
//    fChain->SetBranchAddress("phoE5x5", phoE5x5, &b_phoE5x5);
//    fChain->SetBranchAddress("phoE1x5", phoE1x5, &b_phoE1x5);
//    fChain->SetBranchAddress("phoE2x5Max", phoE2x5Max, &b_phoE2x5Max);
//    fChain->SetBranchAddress("phoE2x5Right", phoE2x5Right, &b_phoE2x5Right);
//    fChain->SetBranchAddress("phoE2x5Left", phoE2x5Left, &b_phoE2x5Left);
//    fChain->SetBranchAddress("phoE2x5Top", phoE2x5Top, &b_phoE2x5Top);
//    fChain->SetBranchAddress("phoE2x5Bottom", phoE2x5Bottom, &b_phoE2x5Bottom);
//    fChain->SetBranchAddress("phoPFIsoNeutral", phoPFIsoNeutral, &b_phoPFIsoNeutral);
//    fChain->SetBranchAddress("phoPFIsoPhoton", phoPFIsoPhoton, &b_phoPFIsoPhoton);
//    fChain->SetBranchAddress("phoPFIsoCharged", phoPFIsoCharged, &b_phoPFIsoCharged);
//    fChain->SetBranchAddress("phoResiCorrE", phoResiCorrE, &b_phoResiCorrE);
//    fChain->SetBranchAddress("phoResiCorrEsigma", phoResiCorrEsigma, &b_phoResiCorrEsigma);
//    fChain->SetBranchAddress("phoRegrE", phoRegrE, &b_phoRegrE);
//    fChain->SetBranchAddress("phoRegrEerr", phoRegrEerr, &b_phoRegrEerr);
//    fChain->SetBranchAddress("phoetaC", phoetaC, &b_phoetaC);
//    fChain->SetBranchAddress("phophiC", phophiC, &b_phophiC);
//    fChain->SetBranchAddress("phoetaS", phoetaS, &b_phoetaS);
//    fChain->SetBranchAddress("phophiS", phophiS, &b_phophiS);
//    fChain->SetBranchAddress("phoetaM", phoetaM, &b_phoetaM);
//    fChain->SetBranchAddress("phophiM", phophiM, &b_phophiM);
//    fChain->SetBranchAddress("phoSeedTime", phoSeedTime, &b_phoSeedTime);
//    fChain->SetBranchAddress("phoSeedDetId1", phoSeedDetId1, &b_phoSeedDetId1);
//    fChain->SetBranchAddress("phoSeedDetId2", phoSeedDetId2, &b_phoSeedDetId2);
//    fChain->SetBranchAddress("phoRecoFlag", phoRecoFlag, &b_phoRecoFlag);
//    fChain->SetBranchAddress("phoPos", phoPos, &b_phoPos);
   fChain->SetBranchAddress("phoGenIndex", phoGenIndex, &b_phoGenIndex);
   fChain->SetBranchAddress("phoGenGMomPID", phoGenGMomPID, &b_phoGenGMomPID);
   fChain->SetBranchAddress("phoGenMomPID", phoGenMomPID, &b_phoGenMomPID);
//    fChain->SetBranchAddress("phoGenMomPt", phoGenMomPt, &b_phoGenMomPt);
//    fChain->SetBranchAddress("phoSCE", phoSCE, &b_phoSCE);
//    fChain->SetBranchAddress("phoSCRawE", phoSCRawE, &b_phoSCRawE);
//    fChain->SetBranchAddress("phoESEn", phoESEn, &b_phoESEn);
//    fChain->SetBranchAddress("phoSCEt", phoSCEt, &b_phoSCEt);
   fChain->SetBranchAddress("phoSCEta", phoSCEta, &b_phoSCEta);
//    fChain->SetBranchAddress("phoSCPhi", phoSCPhi, &b_phoSCPhi);
//    fChain->SetBranchAddress("phoSCEtaWidth", phoSCEtaWidth, &b_phoSCEtaWidth);
//    fChain->SetBranchAddress("phoSCPhiWidth", phoSCPhiWidth, &b_phoSCPhiWidth);
//    fChain->SetBranchAddress("phoSCBrem", phoSCBrem, &b_phoSCBrem);
//    fChain->SetBranchAddress("phoOverlap", phoOverlap, &b_phoOverlap);
   fChain->SetBranchAddress("phohasPixelSeed", phohasPixelSeed, &b_phohasPixelSeed);
//    fChain->SetBranchAddress("phoIsConv", phoIsConv, &b_phoIsConv);
//    fChain->SetBranchAddress("phoNConv", phoNConv, &b_phoNConv);
//    fChain->SetBranchAddress("phoConvInvMass", phoConvInvMass, &b_phoConvInvMass);
//    fChain->SetBranchAddress("phoConvCotTheta", phoConvCotTheta, &b_phoConvCotTheta);
//    fChain->SetBranchAddress("phoConvEoverP", phoConvEoverP, &b_phoConvEoverP);
//    fChain->SetBranchAddress("phoConvZofPVfromTrks", phoConvZofPVfromTrks, &b_phoConvZofPVfromTrks);
//    fChain->SetBranchAddress("phoConvMinDist", phoConvMinDist, &b_phoConvMinDist);
//    fChain->SetBranchAddress("phoConvdPhiAtVtx", phoConvdPhiAtVtx, &b_phoConvdPhiAtVtx);
//    fChain->SetBranchAddress("phoConvdPhiAtCalo", phoConvdPhiAtCalo, &b_phoConvdPhiAtCalo);
//    fChain->SetBranchAddress("phoConvdEtaAtCalo", phoConvdEtaAtCalo, &b_phoConvdEtaAtCalo);
//    fChain->SetBranchAddress("phoConvTrkd0", phoConvTrkd0, &b_phoConvTrkd0);
//    fChain->SetBranchAddress("phoConvTrkPin", phoConvTrkPin, &b_phoConvTrkPin);
//    fChain->SetBranchAddress("phoConvTrkPout", phoConvTrkPout, &b_phoConvTrkPout);
//    fChain->SetBranchAddress("phoConvTrkdz", phoConvTrkdz, &b_phoConvTrkdz);
//    fChain->SetBranchAddress("phoConvTrkdzErr", phoConvTrkdzErr, &b_phoConvTrkdzErr);
//    fChain->SetBranchAddress("phoConvChi2", phoConvChi2, &b_phoConvChi2);
//    fChain->SetBranchAddress("phoConvChi2Prob", phoConvChi2Prob, &b_phoConvChi2Prob);
//    fChain->SetBranchAddress("phoConvCharge", phoConvCharge, &b_phoConvCharge);
//    fChain->SetBranchAddress("phoConvValidVtx", phoConvValidVtx, &b_phoConvValidVtx);
//    fChain->SetBranchAddress("phoConvLikeLihood", phoConvLikeLihood, &b_phoConvLikeLihood);
//    fChain->SetBranchAddress("phoConvP4", phoConvP4, &b_phoConvP4);
//    fChain->SetBranchAddress("phoConvVtx", phoConvVtx, &b_phoConvVtx);
//    fChain->SetBranchAddress("phoConvVtxErr", phoConvVtxErr, &b_phoConvVtxErr);
//    fChain->SetBranchAddress("phoConvPairMomentum", phoConvPairMomentum, &b_phoConvPairMomentum);
//    fChain->SetBranchAddress("phoConvRefittedMomentum", phoConvRefittedMomentum, &b_phoConvRefittedMomentum);
//    fChain->SetBranchAddress("phoESEffSigmaRR", phoESEffSigmaRR, &b_phoESEffSigmaRR);
//    fChain->SetBranchAddress("phoESE1", phoESE1, &b_phoESE1);
//    fChain->SetBranchAddress("phoESE3", phoESE3, &b_phoESE3);
//    fChain->SetBranchAddress("phoESE5", phoESE5, &b_phoESE5);
//    fChain->SetBranchAddress("phoESE7", phoESE7, &b_phoESE7);
//    fChain->SetBranchAddress("phoESE11", phoESE11, &b_phoESE11);
//    fChain->SetBranchAddress("phoESE21", phoESE21, &b_phoESE21);
//    fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
//    fChain->SetBranchAddress("muTrg", muTrg, &b_muTrg);
//    fChain->SetBranchAddress("muEta", muEta, &b_muEta);
//    fChain->SetBranchAddress("muPhi", muPhi, &b_muPhi);
//    fChain->SetBranchAddress("muCharge", muCharge, &b_muCharge);
//    fChain->SetBranchAddress("muPt", muPt, &b_muPt);
//    fChain->SetBranchAddress("muPz", muPz, &b_muPz);
//    fChain->SetBranchAddress("muGenIndex", muGenIndex, &b_muGenIndex);
//    fChain->SetBranchAddress("muIsoTrk", muIsoTrk, &b_muIsoTrk);
//    fChain->SetBranchAddress("muIsoCalo", muIsoCalo, &b_muIsoCalo);
//    fChain->SetBranchAddress("muIsoEcal", muIsoEcal, &b_muIsoEcal);
//    fChain->SetBranchAddress("muIsoHcal", muIsoHcal, &b_muIsoHcal);
//    fChain->SetBranchAddress("muChi2NDF", muChi2NDF, &b_muChi2NDF);
//    fChain->SetBranchAddress("muEmVeto", muEmVeto, &b_muEmVeto);
//    fChain->SetBranchAddress("muHadVeto", muHadVeto, &b_muHadVeto);
//    fChain->SetBranchAddress("muType", muType, &b_muType);
//    fChain->SetBranchAddress("muID", muID, &b_muID);
//    fChain->SetBranchAddress("muD0", muD0, &b_muD0);
//    fChain->SetBranchAddress("muDz", muDz, &b_muDz);
//    fChain->SetBranchAddress("muD0Vtx", muD0Vtx, &b_muD0Vtx);
//    fChain->SetBranchAddress("muDzVtx", muDzVtx, &b_muDzVtx);
//    fChain->SetBranchAddress("muNumberOfValidTrkHits", muNumberOfValidTrkHits, &b_muNumberOfValidTrkHits);
//    fChain->SetBranchAddress("muNumberOfValidPixelHits", muNumberOfValidPixelHits, &b_muNumberOfValidPixelHits);
//    fChain->SetBranchAddress("muNumberOfValidMuonHits", muNumberOfValidMuonHits, &b_muNumberOfValidMuonHits);
//    fChain->SetBranchAddress("muStations", muStations, &b_muStations);
//    fChain->SetBranchAddress("muChambers", muChambers, &b_muChambers);
//    fChain->SetBranchAddress("nPFPho_", &nPFPho_, &b_nPFPho_);
//    fChain->SetBranchAddress("PFPhoE_", PFPhoE_, &b_PFPhoE_);
//    fChain->SetBranchAddress("PFPhoEt_", PFPhoEt_, &b_PFPhoEt_);
//    fChain->SetBranchAddress("PFPhoEta_", PFPhoEta_, &b_PFPhoEta_);
//    fChain->SetBranchAddress("PFPhoPhi_", PFPhoPhi_, &b_PFPhoPhi_);
//    fChain->SetBranchAddress("PFPhoHoverE_", PFPhoHoverE_, &b_PFPhoHoverE_);
//    fChain->SetBranchAddress("nPFEle_", &nPFEle_, &b_nPFEle_);
//    fChain->SetBranchAddress("PFElePt_", PFElePt_, &b_PFElePt_);
//    fChain->SetBranchAddress("PFEleEta_", PFEleEta_, &b_PFEleEta_);
//    fChain->SetBranchAddress("PFElePhi_", PFElePhi_, &b_PFElePhi_);
//    fChain->SetBranchAddress("PFEleEn_", PFEleEn_, &b_PFEleEn_);
//    fChain->SetBranchAddress("PFEleCharge", PFEleCharge, &b_PFEleCharge);
//    fChain->SetBranchAddress("nPFMu_", &nPFMu_, &b_nPFMu_);
//    fChain->SetBranchAddress("PFmuCharge_", PFmuCharge_, &b_PFmuCharge_);
//    fChain->SetBranchAddress("PFmuPhi_", PFmuPhi_, &b_PFmuPhi_);
//    fChain->SetBranchAddress("PFmuEta_", PFmuEta_, &b_PFmuEta_);
//    fChain->SetBranchAddress("PFmuPt_", PFmuPt_, &b_PFmuPt_);
//    fChain->SetBranchAddress("PFmuPz_", PFmuPz_, &b_PFmuPz_);
//    fChain->SetBranchAddress("PFmuChambers_", PFmuChambers_, &b_PFmuChambers_);
//    fChain->SetBranchAddress("PFmuD0_", PFmuD0_, &b_PFmuD0_);
//    fChain->SetBranchAddress("PFmuDz_", PFmuDz_, &b_PFmuDz_);
//    fChain->SetBranchAddress("PFmuChi2NDF_", PFmuChi2NDF_, &b_PFmuChi2NDF_);
//    fChain->SetBranchAddress("PFmuNumberOfValidTrkHits_", PFmuNumberOfValidTrkHits_, &b_PFmuNumberOfValidTrkHits_);
//    fChain->SetBranchAddress("PFmuNumberOfValidPixelHits_", PFmuNumberOfValidPixelHits_, &b_PFmuNumberOfValidPixelHits_);
//    fChain->SetBranchAddress("PFmuNumberOfValidMuonHits_", PFmuNumberOfValidMuonHits_, &b_PFmuNumberOfValidMuonHits_);
   fChain->SetBranchAddress("rho25", &rho25, &b_rho25);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
//    fChain->SetBranchAddress("jetTrg", jetTrg, &b_jetTrg);
   fChain->SetBranchAddress("jetEn", jetEn, &b_jetEn);
   fChain->SetBranchAddress("jetPt", jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetEt", jetEt, &b_jetEt);
//    fChain->SetBranchAddress("jetRawPt", jetRawPt, &b_jetRawPt);
//    fChain->SetBranchAddress("jetRawEn", jetRawEn, &b_jetRawEn);
   fChain->SetBranchAddress("jetCHF", jetCHF, &b_jetCHF);
   fChain->SetBranchAddress("jetNHF", jetNHF, &b_jetNHF);
   fChain->SetBranchAddress("jetCEF", jetCEF, &b_jetCEF);
   fChain->SetBranchAddress("jetNEF", jetNEF, &b_jetNEF);
   fChain->SetBranchAddress("jetNCH", jetNCH, &b_jetNCH);
//    fChain->SetBranchAddress("jetHFHAE", jetHFHAE, &b_jetHFHAE);
//    fChain->SetBranchAddress("jetHFEME", jetHFEME, &b_jetHFEME);
   fChain->SetBranchAddress("jetNConstituents", jetNConstituents, &b_jetNConstituents);
//    fChain->SetBranchAddress("jetTrackCountHiEffBJetTags", jetTrackCountHiEffBJetTags, &b_jetTrackCountHiEffBJetTags);
//    fChain->SetBranchAddress("jetTrackCountHiPurBJetTags", jetTrackCountHiPurBJetTags, &b_jetTrackCountHiPurBJetTags);
//    fChain->SetBranchAddress("jetSimpleSVHiEffBJetTags", jetSimpleSVHiEffBJetTags, &b_jetSimpleSVHiEffBJetTags);
//    fChain->SetBranchAddress("jetSimpleSVHiPurBJetTags", jetSimpleSVHiPurBJetTags, &b_jetSimpleSVHiPurBJetTags);
   fChain->SetBranchAddress("jetPartonID", jetPartonID, &b_jetPartonID);
   fChain->SetBranchAddress("jetGenJetIndex", jetGenJetIndex, &b_jetGenJetIndex);
   fChain->SetBranchAddress("jetGenJetEn", jetGenJetEn, &b_jetGenJetEn);
   fChain->SetBranchAddress("jetGenJetPt", jetGenJetPt, &b_jetGenJetPt);
   fChain->SetBranchAddress("jetGenJetEta", jetGenJetEta, &b_jetGenJetEta);
   fChain->SetBranchAddress("jetGenJetPhi", jetGenJetPhi, &b_jetGenJetPhi);
   fChain->SetBranchAddress("jetGenPartonID", jetGenPartonID, &b_jetGenPartonID);
   fChain->SetBranchAddress("jetGenEn", jetGenEn, &b_jetGenEn);
   fChain->SetBranchAddress("jetGenPt", jetGenPt, &b_jetGenPt);
   fChain->SetBranchAddress("jetGenEta", jetGenEta, &b_jetGenEta);
   fChain->SetBranchAddress("jetGenPhi", jetGenPhi, &b_jetGenPhi);
//    fChain->SetBranchAddress("nZee", &nZee, &b_nZee);
//    fChain->SetBranchAddress("ZeeMass", ZeeMass, &b_ZeeMass);
//    fChain->SetBranchAddress("ZeePt", ZeePt, &b_ZeePt);
//    fChain->SetBranchAddress("ZeeEta", ZeeEta, &b_ZeeEta);
//    fChain->SetBranchAddress("ZeePhi", ZeePhi, &b_ZeePhi);
//    fChain->SetBranchAddress("ZeeLeg1Index", ZeeLeg1Index, &b_ZeeLeg1Index);
//    fChain->SetBranchAddress("ZeeLeg2Index", ZeeLeg2Index, &b_ZeeLeg2Index);
//    fChain->SetBranchAddress("nZmumu", &nZmumu, &b_nZmumu);
//    fChain->SetBranchAddress("ZmumuMass", ZmumuMass, &b_ZmumuMass);
//    fChain->SetBranchAddress("ZmumuPt", ZmumuPt, &b_ZmumuPt);
//    fChain->SetBranchAddress("ZmumuEta", ZmumuEta, &b_ZmumuEta);
//    fChain->SetBranchAddress("ZmumuPhi", ZmumuPhi, &b_ZmumuPhi);
//    fChain->SetBranchAddress("ZmumuLeg1Index", ZmumuLeg1Index, &b_ZmumuLeg1Index);
//    fChain->SetBranchAddress("ZmumuLeg2Index", ZmumuLeg2Index, &b_ZmumuLeg2Index);
//    fChain->SetBranchAddress("nWenu", &nWenu, &b_nWenu);
//    fChain->SetBranchAddress("WenuMassTPfMET", WenuMassTPfMET, &b_WenuMassTPfMET);
//    fChain->SetBranchAddress("WenuEtPfMET", WenuEtPfMET, &b_WenuEtPfMET);
//    fChain->SetBranchAddress("WenuACopPfMET", WenuACopPfMET, &b_WenuACopPfMET);
//    fChain->SetBranchAddress("WenuEleIndex", WenuEleIndex, &b_WenuEleIndex);
//    fChain->SetBranchAddress("nWmunu", &nWmunu, &b_nWmunu);
//    fChain->SetBranchAddress("WmunuMassTPfMET", WmunuMassTPfMET, &b_WmunuMassTPfMET);
//    fChain->SetBranchAddress("WmunuEtPfMET", WmunuEtPfMET, &b_WmunuEtPfMET);
//    fChain->SetBranchAddress("WmunuACopPfMET", WmunuACopPfMET, &b_WmunuACopPfMET);
//    fChain->SetBranchAddress("WmunuMuIndex", WmunuMuIndex, &b_WmunuMuIndex);
//    fChain->SetBranchAddress("nConv", &nConv, &b_nConv);
//    fChain->SetBranchAddress("convP4", convP4, &b_convP4);
//    fChain->SetBranchAddress("convVtx", convVtx, &b_convVtx);
//    fChain->SetBranchAddress("convVtxErr", convVtxErr, &b_convVtxErr);
//    fChain->SetBranchAddress("convPairMomentum", convPairMomentum, &b_convPairMomentum);
//    fChain->SetBranchAddress("convRefittedMomentum", convRefittedMomentum, &b_convRefittedMomentum);
//    fChain->SetBranchAddress("convNTracks", convNTracks, &b_convNTracks);
//    fChain->SetBranchAddress("convPairInvMass", convPairInvMass, &b_convPairInvMass);
//    fChain->SetBranchAddress("convPairCotThetaSep", convPairCotThetaSep, &b_convPairCotThetaSep);
//    fChain->SetBranchAddress("convEoverP", convEoverP, &b_convEoverP);
//    fChain->SetBranchAddress("convDistOfMinApproach", convDistOfMinApproach, &b_convDistOfMinApproach);
//    fChain->SetBranchAddress("convDPhiTrksAtVtx", convDPhiTrksAtVtx, &b_convDPhiTrksAtVtx);
//    fChain->SetBranchAddress("convDPhiTrksAtEcal", convDPhiTrksAtEcal, &b_convDPhiTrksAtEcal);
//    fChain->SetBranchAddress("convDEtaTrksAtEcal", convDEtaTrksAtEcal, &b_convDEtaTrksAtEcal);
//    fChain->SetBranchAddress("convDxy", convDxy, &b_convDxy);
//    fChain->SetBranchAddress("convDz", convDz, &b_convDz);
//    fChain->SetBranchAddress("convLxy", convLxy, &b_convLxy);
//    fChain->SetBranchAddress("convLz", convLz, &b_convLz);
//    fChain->SetBranchAddress("convZofPrimVtxFromTrks", convZofPrimVtxFromTrks, &b_convZofPrimVtxFromTrks);
//    fChain->SetBranchAddress("convNHitsBeforeVtx", convNHitsBeforeVtx, &b_convNHitsBeforeVtx);
//    fChain->SetBranchAddress("convNSharedHits", convNSharedHits, &b_convNSharedHits);
//    fChain->SetBranchAddress("convValidVtx", convValidVtx, &b_convValidVtx);
//    fChain->SetBranchAddress("convMVALikelihood", convMVALikelihood, &b_convMVALikelihood);
//    fChain->SetBranchAddress("convChi2", convChi2, &b_convChi2);
//    fChain->SetBranchAddress("convChi2Probability", convChi2Probability, &b_convChi2Probability);
//    fChain->SetBranchAddress("convTk1Dz", convTk1Dz, &b_convTk1Dz);
//    fChain->SetBranchAddress("convTk2Dz", convTk2Dz, &b_convTk2Dz);
//    fChain->SetBranchAddress("convTk1DzErr", convTk1DzErr, &b_convTk1DzErr);
//    fChain->SetBranchAddress("convTk2DzErr", convTk2DzErr, &b_convTk2DzErr);
//    fChain->SetBranchAddress("convCh1Ch2", convCh1Ch2, &b_convCh1Ch2);
//    fChain->SetBranchAddress("convTk1D0", convTk1D0, &b_convTk1D0);
//    fChain->SetBranchAddress("convTk1Pout", convTk1Pout, &b_convTk1Pout);
//    fChain->SetBranchAddress("convTk1Pin", convTk1Pin, &b_convTk1Pin);
//    fChain->SetBranchAddress("convTk2D0", convTk2D0, &b_convTk2D0);
//    fChain->SetBranchAddress("convTk2Pout", convTk2Pout, &b_convTk2Pout);
//    fChain->SetBranchAddress("convTk2Pin", convTk2Pin, &b_convTk2Pin);
   Notify();
}

Bool_t gg_angularmc_eff::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gg_angularmc_eff::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gg_angularmc_eff::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gg_angularmc_eff_cxx
