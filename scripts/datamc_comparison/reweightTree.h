//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb 11 22:27:19 2010 by ROOT version 5.22/00d
// from TTree Analysis/Analysis
// found on file: jan29mc_900GeV_patch4_purity.root
//////////////////////////////////////////////////////////

#ifndef reweightTree_h
#define reweightTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class reweightTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           orbit;
   Int_t           bunchCrossing;
   Int_t           luminosityBlock;
   Int_t           timesec;
   Int_t           nTTBits;
   Int_t           kMaxTTBits;
   UChar_t         TTBit[64];   //[nTTBits]
   UChar_t         TTBit34;
   UChar_t         TTBit40;
   UChar_t         TTBit41;
   UChar_t         TTBit0;
   Int_t           kMaxTrigFlag;
   Int_t           nHLTBits;
   UChar_t         HLTBit[123];   //[nHLTBits]
   Int_t           nHfTowersP;
   Int_t           nHfTowersN;
   UChar_t         clusVtxCut;
   Double_t        clusVtxQual;
   Int_t           nHits;
   Float_t         simVertexX;
   Float_t         simVertexY;
   Float_t         simVertexZ;
   Float_t         ptHat;
   Double_t        beamSpotX;
   Double_t        beamSpotY;
   Double_t        beamSpotZ;
   UChar_t         vtxIsFake;
   Double_t        vtxX;
   Double_t        vtxY;
   Double_t        vtxZ;
   Double_t        vtxXError;
   Double_t        vtxYError;
   Double_t        vtxZError;
   Int_t           vtxNTrk;
   Int_t           vtxNTrkWeight50;
   Double_t        vtxsumw;
   Double_t        vtxChi2;
   Double_t        vtxNdof;
   Double_t        vtxNormChi2;
   Int_t           nGdTrks;
   Int_t           isMETEmpty;
   Double_t        metEt;
   Double_t        metPx;
   Double_t        metPy;
   Double_t        metPz;
   Double_t        e_longitudinal;
   Int_t           nCorrections;
   Float_t         corEx;
   Float_t         corEy;
   Float_t         corSumEt;
   Float_t         uncorrectedPt;
   Float_t         uncorrectedPhi;
   UChar_t         isCaloMET;
   UChar_t         isRecoMET;
   Double_t        maxEtInEmTowers;
   Double_t        maxEtInHadTowers;
   Double_t        etFractionHadronic;
   Double_t        emEtFraction;
   Double_t        hadEtInHB;
   Double_t        hadEtInHO;
   Double_t        hadEtInHE;
   Double_t        hadEtInHF;
   Double_t        metSignificance;
   Double_t        CaloSETInpHF;
   Double_t        CaloSETInmHF;
   Double_t        CaloMETInpHF;
   Double_t        CaloMETInmHF;
   Double_t        CaloMETPhiInpHF;
   Double_t        CaloMETPhiInmHF;
   Int_t           nJets;
   Int_t           kMaxJets;
   Float_t         jetPt[20];   //[kMaxJets]
   Float_t         jetE[20];   //[kMaxJets]
   Float_t         jetP[20];   //[kMaxJets]
   Float_t         jetEta[20];   //[kMaxJets]
   Float_t         jetPhi[20];   //[kMaxJets]
   Float_t         jetCharge[20];   //[kMaxJets]
   Int_t           jetNtrk[20];   //[kMaxJets]
   UChar_t         isCaloJet[20];   //[kMaxJets]
   UChar_t         isPFJet[20];   //[kMaxJets]
   UChar_t         isBasicJet[20];   //[kMaxJets]
   Float_t         maxEInEmTowers[20];   //[kMaxJets]
   Float_t         maxEInHadTowers[20];   //[kMaxJets]
   Float_t         energyFractionHadronic[20];   //[kMaxJets]
   Float_t         emEnergyFraction[20];   //[kMaxJets]
   Float_t         hadEnergyInHB[20];   //[kMaxJets]
   Float_t         hadEnergyInHO[20];   //[kMaxJets]
   Float_t         hadEnergyInHE[20];   //[kMaxJets]
   Float_t         hadEnergyInHF[20];   //[kMaxJets]
   Float_t         emEnergyInEB[20];   //[kMaxJets]
   Float_t         emEnergyInEE[20];   //[kMaxJets]
   Float_t         emEnergyInHF[20];   //[kMaxJets]
   Float_t         towersArea[20];   //[kMaxJets]
   Int_t           n90[20];   //[kMaxJets]
   Int_t           n60[20];   //[kMaxJets]
   Float_t         fHPD[20];   //[kMaxJets]
   Double_t        PHOLEAD_energy;
   Double_t        PHOLEAD_pt;
   Double_t        PHOLEAD_eta;
   Double_t        PHOLEAD_phi;
   Double_t        PHOLEAD_p;
   Double_t        PHOLEAD_et;
   Double_t        PHOLEAD_momentumX;
   Double_t        PHOLEAD_momentumY;
   Double_t        PHOLEAD_momentumZ;
   Float_t         PHOLEAD_r9;
   Int_t           PHOLEAD_isEBGap;
   Int_t           PHOLEAD_isEEGap;
   Int_t           PHOLEAD_isEBEEGap;
   Int_t           PHOLEAD_isTransGap;
   Int_t           PHOLEAD_isEB;
   Int_t           PHOLEAD_isEE;
   Double_t        PHOLEAD_rawEnergy;
   Double_t        PHOLEAD_preshowerEnergy;
   Int_t           PHOLEAD_numOfPreshClusters;
   Int_t           PHOLEAD_clustersSize;
   Double_t        PHOLEAD_phiWidth;
   Double_t        PHOLEAD_etaWidth;
   Double_t        PHOLEAD_scEta;
   Double_t        PHOLEAD_scPhi;
   Float_t         PHOLEAD_eMax;
   Float_t         PHOLEAD_e2nd;
   Float_t         PHOLEAD_e2x2;
   Float_t         PHOLEAD_e3x2;
   Float_t         PHOLEAD_e3x3;
   Float_t         PHOLEAD_e4x4;
   Float_t         PHOLEAD_e5x5;
   Float_t         PHOLEAD_e2x5Right;
   Float_t         PHOLEAD_e2x5Left;
   Float_t         PHOLEAD_e2x5Top;
   Float_t         PHOLEAD_e2x5Bottom;
   Float_t         PHOLEAD_eRight;
   Float_t         PHOLEAD_eLeft;
   Float_t         PHOLEAD_eTop;
   Float_t         PHOLEAD_eBottom;
   Float_t         PHOLEAD_covPhiPhi;
   Float_t         PHOLEAD_covEtaPhi;
   Float_t         PHOLEAD_covEtaEta;
   Float_t         PHOLEAD_maxEnergyXtal;
   Float_t         PHOLEAD_sigmaEtaEta;
   Float_t         PHOLEAD_sigmaIetaIeta;
   Float_t         PHOLEAD_r1x5;
   Float_t         PHOLEAD_r2x5;
   Float_t         PHOLEAD_e1x5;
   Float_t         PHOLEAD_e2x5;
   Float_t         PHOLEAD_hadronicOverEm;
   Float_t         PHOLEAD_hadronicDepth1OverEm;
   Float_t         PHOLEAD_hadronicDepth2OverEm;
   Float_t         PHOLEAD_trackIso;
   Float_t         PHOLEAD_caloIso;
   Float_t         PHOLEAD_ecalIso;
   Float_t         PHOLEAD_hcalIso;
   Float_t         PHOLEAD_ecalRecHitSumEtConeDR04;
   Float_t         PHOLEAD_hcalTowerSumEtConeDR04;
   Float_t         PHOLEAD_hcalDepth1TowerSumEtConeDR04;
   Float_t         PHOLEAD_hcalDepth2TowerSumEtConeDR04;
   Float_t         PHOLEAD_trkSumPtSolidConeDR04;
   Float_t         PHOLEAD_trkSumPtHollowConeDR04;
   Int_t           PHOLEAD_nTrkSolidConeDR04;
   Int_t           PHOLEAD_nTrkHollowConeDR04;
   Float_t         PHOLEAD_ecalRecHitSumEtConeDR03;
   Float_t         PHOLEAD_hcalTowerSumEtConeDR03;
   Float_t         PHOLEAD_hcalDepth1TowerSumEtConeDR03;
   Float_t         PHOLEAD_hcalDepth2TowerSumEtConeDR03;
   Float_t         PHOLEAD_trkSumPtSolidConeDR03;
   Float_t         PHOLEAD_trkSumPtHollowConeDR03;
   Int_t           PHOLEAD_nTrkSolidConeDR03;
   Int_t           PHOLEAD_nTrkHollowConeDR03;
   UChar_t         PHOLEAD_hasConversionTracks;
   UChar_t         PHOLEAD_hasPixelSeed;
   Int_t           PHOLEAD_nTracks;
   UChar_t         PHOLEAD_isConverted;
   Float_t         PHOLEAD_convPairInvariantMass;
   Float_t         PHOLEAD_convpairCotThetaSeparation;
   Float_t         PHOLEAD_convPairMomentumMag;
   Float_t         PHOLEAD_convPairMomentumPerp;
   Float_t         PHOLEAD_convPairMomentumPhi;
   Float_t         PHOLEAD_convPairMomentumEta;
   Float_t         PHOLEAD_convPairMomentumX;
   Float_t         PHOLEAD_convPairMomentumY;
   Float_t         PHOLEAD_convPairMomentumZ;
   Float_t         PHOLEAD_convDistOfMinimumApproach;
   Float_t         PHOLEAD_convDPhiTracksAtVtx;
   Float_t         PHOLEAD_convDPhiTracksAtEcal;
   Float_t         PHOLEAD_convDEtaTracksAtEcal;
   UChar_t         PHOLEAD_convVtxValid;
   Float_t         PHOLEAD_convVtxEta;
   Float_t         PHOLEAD_convVtxPhi;
   Float_t         PHOLEAD_convVtxR;
   Float_t         PHOLEAD_convVtxX;
   Float_t         PHOLEAD_convVtxY;
   Float_t         PHOLEAD_convVtxZ;
   Float_t         PHOLEAD_convVtxChi2;
   Float_t         PHOLEAD_convVtxNdof;
   Float_t         PHOLEAD_convMVALikelihood;
   Float_t         PHOLEAD_convVtxChi2Prob;
   Float_t         PHOLEAD_convEoverP;
   Float_t         PHOLEAD_convzOfPrimaryVertexFromTracks;
   UChar_t         isGenMatched;
   Float_t         genMatchedPt;
   Float_t         genMatchedEta;
   Float_t         genMatchedPhi;
   Int_t           motherID;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_orbit;   //!
   TBranch        *b_bunchCrossing;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_timesec;   //!
   TBranch        *b_nTTBits;   //!
   TBranch        *b_kMaxTTBits;   //!
   TBranch        *b_TTBit;   //!
   TBranch        *b_TTBit34;   //!
   TBranch        *b_TTBit40;   //!
   TBranch        *b_TTBit41;   //!
   TBranch        *b_TTBit0;   //!
   TBranch        *b_kMaxTrigFlag;   //!
   TBranch        *b_nHLTBits;   //!
   TBranch        *b_HLTBit;   //!
   TBranch        *b_nHfTowersP;   //!
   TBranch        *b_nHfTowersN;   //!
   TBranch        *b_clusVtxCut;   //!
   TBranch        *b_clusVtxQual;   //!
   TBranch        *b_nHits;   //!
   TBranch        *b_simVertexX;   //!
   TBranch        *b_simVertexY;   //!
   TBranch        *b_simVertexZ;   //!
   TBranch        *b_ptHat;   //!
   TBranch        *b_beamSpotX;   //!
   TBranch        *b_beamSpotY;   //!
   TBranch        *b_beamSpotZ;   //!
   TBranch        *b_vtxIsFake;   //!
   TBranch        *b_vtxX;   //!
   TBranch        *b_vtxY;   //!
   TBranch        *b_vtxZ;   //!
   TBranch        *b_vtxXError;   //!
   TBranch        *b_vtxYError;   //!
   TBranch        *b_vtxZError;   //!
   TBranch        *b_vtxNTrk;   //!
   TBranch        *b_vtxNTrkWeight50;   //!
   TBranch        *b_vtxsumw;   //!
   TBranch        *b_vtxChi2;   //!
   TBranch        *b_vtxNdof;   //!
   TBranch        *b_vtxNormChi2;   //!
   TBranch        *b_nGdTrks;   //!
   TBranch        *b_isMETEmpty;   //!
   TBranch        *b_metEt;   //!
   TBranch        *b_metPx;   //!
   TBranch        *b_metPy;   //!
   TBranch        *b_metPz;   //!
   TBranch        *b_e_longitudinal;   //!
   TBranch        *b_nCorrections;   //!
   TBranch        *b_corEx;   //!
   TBranch        *b_corEy;   //!
   TBranch        *b_corSumEt;   //!
   TBranch        *b_uncorrectedPt;   //!
   TBranch        *b_uncorrectedPhi;   //!
   TBranch        *b_isCaloMET;   //!
   TBranch        *b_isRecoMET;   //!
   TBranch        *b_maxEtInEmTowers;   //!
   TBranch        *b_maxEtInHadTowers;   //!
   TBranch        *b_etFractionHadronic;   //!
   TBranch        *b_emEtFraction;   //!
   TBranch        *b_hadEtInHB;   //!
   TBranch        *b_hadEtInHO;   //!
   TBranch        *b_hadEtInHE;   //!
   TBranch        *b_hadEtInHF;   //!
   TBranch        *b_metSignificance;   //!
   TBranch        *b_CaloSETInpHF;   //!
   TBranch        *b_CaloSETInmHF;   //!
   TBranch        *b_CaloMETInpHF;   //!
   TBranch        *b_CaloMETInmHF;   //!
   TBranch        *b_CaloMETPhiInpHF;   //!
   TBranch        *b_CaloMETPhiInmHF;   //!
   TBranch        *b_nJets;   //!
   TBranch        *b_kMaxJets;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetE;   //!
   TBranch        *b_jetP;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetCharge;   //!
   TBranch        *b_jetNtrk;   //!
   TBranch        *b_isCaloJet;   //!
   TBranch        *b_isPFJet;   //!
   TBranch        *b_isBasicJet;   //!
   TBranch        *b_maxEInEmTowers;   //!
   TBranch        *b_maxEInHadTowers;   //!
   TBranch        *b_energyFractionHadronic;   //!
   TBranch        *b_emEnergyFraction;   //!
   TBranch        *b_hadEnergyInHB;   //!
   TBranch        *b_hadEnergyInHO;   //!
   TBranch        *b_hadEnergyInHE;   //!
   TBranch        *b_hadEnergyInHF;   //!
   TBranch        *b_emEnergyInEB;   //!
   TBranch        *b_emEnergyInEE;   //!
   TBranch        *b_emEnergyInHF;   //!
   TBranch        *b_towersArea;   //!
   TBranch        *b_n90;   //!
   TBranch        *b_n60;   //!
   TBranch        *b_fHPD;   //!
   TBranch        *b_PHOLEAD_energy;   //!
   TBranch        *b_PHOLEAD_pt;   //!
   TBranch        *b_PHOLEAD_eta;   //!
   TBranch        *b_PHOLEAD_phi;   //!
   TBranch        *b_PHOLEAD_p;   //!
   TBranch        *b_PHOLEAD_et;   //!
   TBranch        *b_PHOLEAD_momentumX;   //!
   TBranch        *b_PHOLEAD_momentumY;   //!
   TBranch        *b_PHOLEAD_momentumZ;   //!
   TBranch        *b_PHOLEAD_r9;   //!
   TBranch        *b_PHOLEAD_isEBGap;   //!
   TBranch        *b_PHOLEAD_isEEGap;   //!
   TBranch        *b_PHOLEAD_isEBEEGap;   //!
   TBranch        *b_PHOLEAD_isTransGap;   //!
   TBranch        *b_PHOLEAD_isEB;   //!
   TBranch        *b_PHOLEAD_isEE;   //!
   TBranch        *b_PHOLEAD_rawEnergy;   //!
   TBranch        *b_PHOLEAD_preshowerEnergy;   //!
   TBranch        *b_PHOLEAD_numOfPreshClusters;   //!
   TBranch        *b_PHOLEAD_clustersSize;   //!
   TBranch        *b_PHOLEAD_phiWidth;   //!
   TBranch        *b_PHOLEAD_etaWidth;   //!
   TBranch        *b_PHOLEAD_scEta;   //!
   TBranch        *b_PHOLEAD_scPhi;   //!
   TBranch        *b_PHOLEAD_eMax;   //!
   TBranch        *b_PHOLEAD_e2nd;   //!
   TBranch        *b_PHOLEAD_e2x2;   //!
   TBranch        *b_PHOLEAD_e3x2;   //!
   TBranch        *b_PHOLEAD_e3x3;   //!
   TBranch        *b_PHOLEAD_e4x4;   //!
   TBranch        *b_PHOLEAD_e5x5;   //!
   TBranch        *b_PHOLEAD_e2x5Right;   //!
   TBranch        *b_PHOLEAD_e2x5Left;   //!
   TBranch        *b_PHOLEAD_e2x5Top;   //!
   TBranch        *b_PHOLEAD_e2x5Bottom;   //!
   TBranch        *b_PHOLEAD_eRight;   //!
   TBranch        *b_PHOLEAD_eLeft;   //!
   TBranch        *b_PHOLEAD_eTop;   //!
   TBranch        *b_PHOLEAD_eBottom;   //!
   TBranch        *b_PHOLEAD_covPhiPhi;   //!
   TBranch        *b_PHOLEAD_covEtaPhi;   //!
   TBranch        *b_PHOLEAD_covEtaEta;   //!
   TBranch        *b_PHOLEAD_maxEnergyXtal;   //!
   TBranch        *b_PHOLEAD_sigmaEtaEta;   //!
   TBranch        *b_PHOLEAD_sigmaIetaIeta;   //!
   TBranch        *b_PHOLEAD_r1x5;   //!
   TBranch        *b_PHOLEAD_r2x5;   //!
   TBranch        *b_PHOLEAD_e1x5;   //!
   TBranch        *b_PHOLEAD_e2x5;   //!
   TBranch        *b_PHOLEAD_hadronicOverEm;   //!
   TBranch        *b_PHOLEAD_hadronicDepth1OverEm;   //!
   TBranch        *b_PHOLEAD_hadronicDepth2OverEm;   //!
   TBranch        *b_PHOLEAD_trackIso;   //!
   TBranch        *b_PHOLEAD_caloIso;   //!
   TBranch        *b_PHOLEAD_ecalIso;   //!
   TBranch        *b_PHOLEAD_hcalIso;   //!
   TBranch        *b_PHOLEAD_ecalRecHitSumEtConeDR04;   //!
   TBranch        *b_PHOLEAD_hcalTowerSumEtConeDR04;   //!
   TBranch        *b_PHOLEAD_hcalDepth1TowerSumEtConeDR04;   //!
   TBranch        *b_PHOLEAD_hcalDepth2TowerSumEtConeDR04;   //!
   TBranch        *b_PHOLEAD_trkSumPtSolidConeDR04;   //!
   TBranch        *b_PHOLEAD_trkSumPtHollowConeDR04;   //!
   TBranch        *b_PHOLEAD_nTrkSolidConeDR04;   //!
   TBranch        *b_PHOLEAD_nTrkHollowConeDR04;   //!
   TBranch        *b_PHOLEAD_ecalRecHitSumEtConeDR03;   //!
   TBranch        *b_PHOLEAD_hcalTowerSumEtConeDR03;   //!
   TBranch        *b_PHOLEAD_hcalDepth1TowerSumEtConeDR03;   //!
   TBranch        *b_PHOLEAD_hcalDepth2TowerSumEtConeDR03;   //!
   TBranch        *b_PHOLEAD_trkSumPtSolidConeDR03;   //!
   TBranch        *b_PHOLEAD_trkSumPtHollowConeDR03;   //!
   TBranch        *b_PHOLEAD_nTrkSolidConeDR03;   //!
   TBranch        *b_PHOLEAD_nTrkHollowConeDR03;   //!
   TBranch        *b_PHOLEAD_hasConversionTracks;   //!
   TBranch        *b_PHOLEAD_hasPixelSeed;   //!
   TBranch        *b_PHOLEAD_nTracks;   //!
   TBranch        *b_PHOLEAD_isConverted;   //!
   TBranch        *b_PHOLEAD_convPairInvariantMass;   //!
   TBranch        *b_PHOLEAD_convpairCotThetaSeparation;   //!
   TBranch        *b_PHOLEAD_convPairMomentumMag;   //!
   TBranch        *b_PHOLEAD_convPairMomentumPerp;   //!
   TBranch        *b_PHOLEAD_convPairMomentumPhi;   //!
   TBranch        *b_PHOLEAD_convPairMomentumEta;   //!
   TBranch        *b_PHOLEAD_convPairMomentumX;   //!
   TBranch        *b_PHOLEAD_convPairMomentumY;   //!
   TBranch        *b_PHOLEAD_convPairMomentumZ;   //!
   TBranch        *b_PHOLEAD_convDistOfMinimumApproach;   //!
   TBranch        *b_PHOLEAD_convDPhiTracksAtVtx;   //!
   TBranch        *b_PHOLEAD_convDPhiTracksAtEcal;   //!
   TBranch        *b_PHOLEAD_convDEtaTracksAtEcal;   //!
   TBranch        *b_PHOLEAD_convVtxValid;   //!
   TBranch        *b_PHOLEAD_convVtxEta;   //!
   TBranch        *b_PHOLEAD_convVtxPhi;   //!
   TBranch        *b_PHOLEAD_convVtxR;   //!
   TBranch        *b_PHOLEAD_convVtxX;   //!
   TBranch        *b_PHOLEAD_convVtxY;   //!
   TBranch        *b_PHOLEAD_convVtxZ;   //!
   TBranch        *b_PHOLEAD_convVtxChi2;   //!
   TBranch        *b_PHOLEAD_convVtxNdof;   //!
   TBranch        *b_PHOLEAD_convMVALikelihood;   //!
   TBranch        *b_PHOLEAD_convVtxChi2Prob;   //!
   TBranch        *b_PHOLEAD_convEoverP;   //!
   TBranch        *b_PHOLEAD_convzOfPrimaryVertexFromTracks;   //!
   TBranch        *b_isGenMatched;   //!
   TBranch        *b_genMatchedPt;   //!
   TBranch        *b_genMatchedEta;   //!
   TBranch        *b_genMatchedPhi;   //!
   TBranch        *b_motherID;   //!

   reweightTree(TTree *tree=0);
   virtual ~reweightTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef reweightTree_cxx
reweightTree::reweightTree(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("lsf_paoloskim_900_mc.root");
      if (!f) {
         f = new TFile("lsf_paoloskim_900_mc.root");
         f->cd("lsf_paoloskim_900_mc.root:/NTuples");
      }
      tree = (TTree*)gDirectory->Get("Analysis");

   }
   Init(tree);
}

reweightTree::~reweightTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t reweightTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t reweightTree::LoadTree(Long64_t entry)
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

void reweightTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("orbit", &orbit, &b_orbit);
   fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_bunchCrossing);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("timesec", &timesec, &b_timesec);
   fChain->SetBranchAddress("nTTBits", &nTTBits, &b_nTTBits);
   fChain->SetBranchAddress("kMaxTTBits", &kMaxTTBits, &b_kMaxTTBits);
   fChain->SetBranchAddress("TTBit", TTBit, &b_TTBit);
   fChain->SetBranchAddress("TTBit34", &TTBit34, &b_TTBit34);
   fChain->SetBranchAddress("TTBit40", &TTBit40, &b_TTBit40);
   fChain->SetBranchAddress("TTBit41", &TTBit41, &b_TTBit41);
   fChain->SetBranchAddress("TTBit0", &TTBit0, &b_TTBit0);
   fChain->SetBranchAddress("kMaxTrigFlag", &kMaxTrigFlag, &b_kMaxTrigFlag);
   fChain->SetBranchAddress("nHLTBits", &nHLTBits, &b_nHLTBits);
   fChain->SetBranchAddress("HLTBit", HLTBit, &b_HLTBit);
   fChain->SetBranchAddress("nHfTowersP", &nHfTowersP, &b_nHfTowersP);
   fChain->SetBranchAddress("nHfTowersN", &nHfTowersN, &b_nHfTowersN);
   fChain->SetBranchAddress("clusVtxCut", &clusVtxCut, &b_clusVtxCut);
   fChain->SetBranchAddress("clusVtxQual", &clusVtxQual, &b_clusVtxQual);
   fChain->SetBranchAddress("nHits", &nHits, &b_nHits);
   fChain->SetBranchAddress("simVertexX", &simVertexX, &b_simVertexX);
   fChain->SetBranchAddress("simVertexY", &simVertexY, &b_simVertexY);
   fChain->SetBranchAddress("simVertexZ", &simVertexZ, &b_simVertexZ);
   fChain->SetBranchAddress("ptHat", &ptHat, &b_ptHat);
   fChain->SetBranchAddress("beamSpotX", &beamSpotX, &b_beamSpotX);
   fChain->SetBranchAddress("beamSpotY", &beamSpotY, &b_beamSpotY);
   fChain->SetBranchAddress("beamSpotZ", &beamSpotZ, &b_beamSpotZ);
   fChain->SetBranchAddress("vtxIsFake", &vtxIsFake, &b_vtxIsFake);
   fChain->SetBranchAddress("vtxX", &vtxX, &b_vtxX);
   fChain->SetBranchAddress("vtxY", &vtxY, &b_vtxY);
   fChain->SetBranchAddress("vtxZ", &vtxZ, &b_vtxZ);
   fChain->SetBranchAddress("vtxXError", &vtxXError, &b_vtxXError);
   fChain->SetBranchAddress("vtxYError", &vtxYError, &b_vtxYError);
   fChain->SetBranchAddress("vtxZError", &vtxZError, &b_vtxZError);
   fChain->SetBranchAddress("vtxNTrk", &vtxNTrk, &b_vtxNTrk);
   fChain->SetBranchAddress("vtxNTrkWeight50", &vtxNTrkWeight50, &b_vtxNTrkWeight50);
   fChain->SetBranchAddress("vtxsumw", &vtxsumw, &b_vtxsumw);
   fChain->SetBranchAddress("vtxChi2", &vtxChi2, &b_vtxChi2);
   fChain->SetBranchAddress("vtxNdof", &vtxNdof, &b_vtxNdof);
   fChain->SetBranchAddress("vtxNormChi2", &vtxNormChi2, &b_vtxNormChi2);
   fChain->SetBranchAddress("nGdTrks", &nGdTrks, &b_nGdTrks);
   fChain->SetBranchAddress("isMETEmpty", &isMETEmpty, &b_isMETEmpty);
   fChain->SetBranchAddress("metEt", &metEt, &b_metEt);
   fChain->SetBranchAddress("metPx", &metPx, &b_metPx);
   fChain->SetBranchAddress("metPy", &metPy, &b_metPy);
   fChain->SetBranchAddress("metPz", &metPz, &b_metPz);
   fChain->SetBranchAddress("e_longitudinal", &e_longitudinal, &b_e_longitudinal);
   fChain->SetBranchAddress("nCorrections", &nCorrections, &b_nCorrections);
   fChain->SetBranchAddress("corEx", &corEx, &b_corEx);
   fChain->SetBranchAddress("corEy", &corEy, &b_corEy);
   fChain->SetBranchAddress("corSumEt", &corSumEt, &b_corSumEt);
   fChain->SetBranchAddress("uncorrectedPt", &uncorrectedPt, &b_uncorrectedPt);
   fChain->SetBranchAddress("uncorrectedPhi", &uncorrectedPhi, &b_uncorrectedPhi);
   fChain->SetBranchAddress("isCaloMET", &isCaloMET, &b_isCaloMET);
   fChain->SetBranchAddress("isRecoMET", &isRecoMET, &b_isRecoMET);
   fChain->SetBranchAddress("maxEtInEmTowers", &maxEtInEmTowers, &b_maxEtInEmTowers);
   fChain->SetBranchAddress("maxEtInHadTowers", &maxEtInHadTowers, &b_maxEtInHadTowers);
   fChain->SetBranchAddress("etFractionHadronic", &etFractionHadronic, &b_etFractionHadronic);
   fChain->SetBranchAddress("emEtFraction", &emEtFraction, &b_emEtFraction);
   fChain->SetBranchAddress("hadEtInHB", &hadEtInHB, &b_hadEtInHB);
   fChain->SetBranchAddress("hadEtInHO", &hadEtInHO, &b_hadEtInHO);
   fChain->SetBranchAddress("hadEtInHE", &hadEtInHE, &b_hadEtInHE);
   fChain->SetBranchAddress("hadEtInHF", &hadEtInHF, &b_hadEtInHF);
   fChain->SetBranchAddress("metSignificance", &metSignificance, &b_metSignificance);
   fChain->SetBranchAddress("CaloSETInpHF", &CaloSETInpHF, &b_CaloSETInpHF);
   fChain->SetBranchAddress("CaloSETInmHF", &CaloSETInmHF, &b_CaloSETInmHF);
   fChain->SetBranchAddress("CaloMETInpHF", &CaloMETInpHF, &b_CaloMETInpHF);
   fChain->SetBranchAddress("CaloMETInmHF", &CaloMETInmHF, &b_CaloMETInmHF);
   fChain->SetBranchAddress("CaloMETPhiInpHF", &CaloMETPhiInpHF, &b_CaloMETPhiInpHF);
   fChain->SetBranchAddress("CaloMETPhiInmHF", &CaloMETPhiInmHF, &b_CaloMETPhiInmHF);
   fChain->SetBranchAddress("nJets", &nJets, &b_nJets);
   fChain->SetBranchAddress("kMaxJets", &kMaxJets, &b_kMaxJets);
   fChain->SetBranchAddress("jetPt", jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetE", jetE, &b_jetE);
   fChain->SetBranchAddress("jetP", jetP, &b_jetP);
   fChain->SetBranchAddress("jetEta", jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetCharge", jetCharge, &b_jetCharge);
   fChain->SetBranchAddress("jetNtrk", jetNtrk, &b_jetNtrk);
   fChain->SetBranchAddress("isCaloJet", isCaloJet, &b_isCaloJet);
   fChain->SetBranchAddress("isPFJet", isPFJet, &b_isPFJet);
   fChain->SetBranchAddress("isBasicJet", isBasicJet, &b_isBasicJet);
   fChain->SetBranchAddress("maxEInEmTowers", maxEInEmTowers, &b_maxEInEmTowers);
   fChain->SetBranchAddress("maxEInHadTowers", maxEInHadTowers, &b_maxEInHadTowers);
   fChain->SetBranchAddress("energyFractionHadronic", energyFractionHadronic, &b_energyFractionHadronic);
   fChain->SetBranchAddress("emEnergyFraction", emEnergyFraction, &b_emEnergyFraction);
   fChain->SetBranchAddress("hadEnergyInHB", hadEnergyInHB, &b_hadEnergyInHB);
   fChain->SetBranchAddress("hadEnergyInHO", hadEnergyInHO, &b_hadEnergyInHO);
   fChain->SetBranchAddress("hadEnergyInHE", hadEnergyInHE, &b_hadEnergyInHE);
   fChain->SetBranchAddress("hadEnergyInHF", hadEnergyInHF, &b_hadEnergyInHF);
   fChain->SetBranchAddress("emEnergyInEB", emEnergyInEB, &b_emEnergyInEB);
   fChain->SetBranchAddress("emEnergyInEE", emEnergyInEE, &b_emEnergyInEE);
   fChain->SetBranchAddress("emEnergyInHF", emEnergyInHF, &b_emEnergyInHF);
   fChain->SetBranchAddress("towersArea", towersArea, &b_towersArea);
   fChain->SetBranchAddress("n90", n90, &b_n90);
   fChain->SetBranchAddress("n60", n60, &b_n60);
   fChain->SetBranchAddress("fHPD", fHPD, &b_fHPD);
   fChain->SetBranchAddress("PHOLEAD_energy", &PHOLEAD_energy, &b_PHOLEAD_energy);
   fChain->SetBranchAddress("PHOLEAD_pt", &PHOLEAD_pt, &b_PHOLEAD_pt);
   fChain->SetBranchAddress("PHOLEAD_eta", &PHOLEAD_eta, &b_PHOLEAD_eta);
   fChain->SetBranchAddress("PHOLEAD_phi", &PHOLEAD_phi, &b_PHOLEAD_phi);
   fChain->SetBranchAddress("PHOLEAD_p", &PHOLEAD_p, &b_PHOLEAD_p);
   fChain->SetBranchAddress("PHOLEAD_et", &PHOLEAD_et, &b_PHOLEAD_et);
   fChain->SetBranchAddress("PHOLEAD_momentumX", &PHOLEAD_momentumX, &b_PHOLEAD_momentumX);
   fChain->SetBranchAddress("PHOLEAD_momentumY", &PHOLEAD_momentumY, &b_PHOLEAD_momentumY);
   fChain->SetBranchAddress("PHOLEAD_momentumZ", &PHOLEAD_momentumZ, &b_PHOLEAD_momentumZ);
   fChain->SetBranchAddress("PHOLEAD_r9", &PHOLEAD_r9, &b_PHOLEAD_r9);
   fChain->SetBranchAddress("PHOLEAD_isEBGap", &PHOLEAD_isEBGap, &b_PHOLEAD_isEBGap);
   fChain->SetBranchAddress("PHOLEAD_isEEGap", &PHOLEAD_isEEGap, &b_PHOLEAD_isEEGap);
   fChain->SetBranchAddress("PHOLEAD_isEBEEGap", &PHOLEAD_isEBEEGap, &b_PHOLEAD_isEBEEGap);
   fChain->SetBranchAddress("PHOLEAD_isTransGap", &PHOLEAD_isTransGap, &b_PHOLEAD_isTransGap);
   fChain->SetBranchAddress("PHOLEAD_isEB", &PHOLEAD_isEB, &b_PHOLEAD_isEB);
   fChain->SetBranchAddress("PHOLEAD_isEE", &PHOLEAD_isEE, &b_PHOLEAD_isEE);
   fChain->SetBranchAddress("PHOLEAD_rawEnergy", &PHOLEAD_rawEnergy, &b_PHOLEAD_rawEnergy);
   fChain->SetBranchAddress("PHOLEAD_preshowerEnergy", &PHOLEAD_preshowerEnergy, &b_PHOLEAD_preshowerEnergy);
   fChain->SetBranchAddress("PHOLEAD_numOfPreshClusters", &PHOLEAD_numOfPreshClusters, &b_PHOLEAD_numOfPreshClusters);
   fChain->SetBranchAddress("PHOLEAD_clustersSize", &PHOLEAD_clustersSize, &b_PHOLEAD_clustersSize);
   fChain->SetBranchAddress("PHOLEAD_phiWidth", &PHOLEAD_phiWidth, &b_PHOLEAD_phiWidth);
   fChain->SetBranchAddress("PHOLEAD_etaWidth", &PHOLEAD_etaWidth, &b_PHOLEAD_etaWidth);
   fChain->SetBranchAddress("PHOLEAD_scEta", &PHOLEAD_scEta, &b_PHOLEAD_scEta);
   fChain->SetBranchAddress("PHOLEAD_scPhi", &PHOLEAD_scPhi, &b_PHOLEAD_scPhi);
   fChain->SetBranchAddress("PHOLEAD_eMax", &PHOLEAD_eMax, &b_PHOLEAD_eMax);
   fChain->SetBranchAddress("PHOLEAD_e2nd", &PHOLEAD_e2nd, &b_PHOLEAD_e2nd);
   fChain->SetBranchAddress("PHOLEAD_e2x2", &PHOLEAD_e2x2, &b_PHOLEAD_e2x2);
   fChain->SetBranchAddress("PHOLEAD_e3x2", &PHOLEAD_e3x2, &b_PHOLEAD_e3x2);
   fChain->SetBranchAddress("PHOLEAD_e3x3", &PHOLEAD_e3x3, &b_PHOLEAD_e3x3);
   fChain->SetBranchAddress("PHOLEAD_e4x4", &PHOLEAD_e4x4, &b_PHOLEAD_e4x4);
   fChain->SetBranchAddress("PHOLEAD_e5x5", &PHOLEAD_e5x5, &b_PHOLEAD_e5x5);
   fChain->SetBranchAddress("PHOLEAD_e2x5Right", &PHOLEAD_e2x5Right, &b_PHOLEAD_e2x5Right);
   fChain->SetBranchAddress("PHOLEAD_e2x5Left", &PHOLEAD_e2x5Left, &b_PHOLEAD_e2x5Left);
   fChain->SetBranchAddress("PHOLEAD_e2x5Top", &PHOLEAD_e2x5Top, &b_PHOLEAD_e2x5Top);
   fChain->SetBranchAddress("PHOLEAD_e2x5Bottom", &PHOLEAD_e2x5Bottom, &b_PHOLEAD_e2x5Bottom);
   fChain->SetBranchAddress("PHOLEAD_eRight", &PHOLEAD_eRight, &b_PHOLEAD_eRight);
   fChain->SetBranchAddress("PHOLEAD_eLeft", &PHOLEAD_eLeft, &b_PHOLEAD_eLeft);
   fChain->SetBranchAddress("PHOLEAD_eTop", &PHOLEAD_eTop, &b_PHOLEAD_eTop);
   fChain->SetBranchAddress("PHOLEAD_eBottom", &PHOLEAD_eBottom, &b_PHOLEAD_eBottom);
   fChain->SetBranchAddress("PHOLEAD_covPhiPhi", &PHOLEAD_covPhiPhi, &b_PHOLEAD_covPhiPhi);
   fChain->SetBranchAddress("PHOLEAD_covEtaPhi", &PHOLEAD_covEtaPhi, &b_PHOLEAD_covEtaPhi);
   fChain->SetBranchAddress("PHOLEAD_covEtaEta", &PHOLEAD_covEtaEta, &b_PHOLEAD_covEtaEta);
   fChain->SetBranchAddress("PHOLEAD_maxEnergyXtal", &PHOLEAD_maxEnergyXtal, &b_PHOLEAD_maxEnergyXtal);
   fChain->SetBranchAddress("PHOLEAD_sigmaEtaEta", &PHOLEAD_sigmaEtaEta, &b_PHOLEAD_sigmaEtaEta);
   fChain->SetBranchAddress("PHOLEAD_sigmaIetaIeta", &PHOLEAD_sigmaIetaIeta, &b_PHOLEAD_sigmaIetaIeta);
   fChain->SetBranchAddress("PHOLEAD_r1x5", &PHOLEAD_r1x5, &b_PHOLEAD_r1x5);
   fChain->SetBranchAddress("PHOLEAD_r2x5", &PHOLEAD_r2x5, &b_PHOLEAD_r2x5);
   fChain->SetBranchAddress("PHOLEAD_e1x5", &PHOLEAD_e1x5, &b_PHOLEAD_e1x5);
   fChain->SetBranchAddress("PHOLEAD_e2x5", &PHOLEAD_e2x5, &b_PHOLEAD_e2x5);
   fChain->SetBranchAddress("PHOLEAD_hadronicOverEm", &PHOLEAD_hadronicOverEm, &b_PHOLEAD_hadronicOverEm);
   fChain->SetBranchAddress("PHOLEAD_hadronicDepth1OverEm", &PHOLEAD_hadronicDepth1OverEm, &b_PHOLEAD_hadronicDepth1OverEm);
   fChain->SetBranchAddress("PHOLEAD_hadronicDepth2OverEm", &PHOLEAD_hadronicDepth2OverEm, &b_PHOLEAD_hadronicDepth2OverEm);
   fChain->SetBranchAddress("PHOLEAD_trackIso", &PHOLEAD_trackIso, &b_PHOLEAD_trackIso);
   fChain->SetBranchAddress("PHOLEAD_caloIso", &PHOLEAD_caloIso, &b_PHOLEAD_caloIso);
   fChain->SetBranchAddress("PHOLEAD_ecalIso", &PHOLEAD_ecalIso, &b_PHOLEAD_ecalIso);
   fChain->SetBranchAddress("PHOLEAD_hcalIso", &PHOLEAD_hcalIso, &b_PHOLEAD_hcalIso);
   fChain->SetBranchAddress("PHOLEAD_ecalRecHitSumEtConeDR04", &PHOLEAD_ecalRecHitSumEtConeDR04, &b_PHOLEAD_ecalRecHitSumEtConeDR04);
   fChain->SetBranchAddress("PHOLEAD_hcalTowerSumEtConeDR04", &PHOLEAD_hcalTowerSumEtConeDR04, &b_PHOLEAD_hcalTowerSumEtConeDR04);
   fChain->SetBranchAddress("PHOLEAD_hcalDepth1TowerSumEtConeDR04", &PHOLEAD_hcalDepth1TowerSumEtConeDR04, &b_PHOLEAD_hcalDepth1TowerSumEtConeDR04);
   fChain->SetBranchAddress("PHOLEAD_hcalDepth2TowerSumEtConeDR04", &PHOLEAD_hcalDepth2TowerSumEtConeDR04, &b_PHOLEAD_hcalDepth2TowerSumEtConeDR04);
   fChain->SetBranchAddress("PHOLEAD_trkSumPtSolidConeDR04", &PHOLEAD_trkSumPtSolidConeDR04, &b_PHOLEAD_trkSumPtSolidConeDR04);
   fChain->SetBranchAddress("PHOLEAD_trkSumPtHollowConeDR04", &PHOLEAD_trkSumPtHollowConeDR04, &b_PHOLEAD_trkSumPtHollowConeDR04);
   fChain->SetBranchAddress("PHOLEAD_nTrkSolidConeDR04", &PHOLEAD_nTrkSolidConeDR04, &b_PHOLEAD_nTrkSolidConeDR04);
   fChain->SetBranchAddress("PHOLEAD_nTrkHollowConeDR04", &PHOLEAD_nTrkHollowConeDR04, &b_PHOLEAD_nTrkHollowConeDR04);
   fChain->SetBranchAddress("PHOLEAD_ecalRecHitSumEtConeDR03", &PHOLEAD_ecalRecHitSumEtConeDR03, &b_PHOLEAD_ecalRecHitSumEtConeDR03);
   fChain->SetBranchAddress("PHOLEAD_hcalTowerSumEtConeDR03", &PHOLEAD_hcalTowerSumEtConeDR03, &b_PHOLEAD_hcalTowerSumEtConeDR03);
   fChain->SetBranchAddress("PHOLEAD_hcalDepth1TowerSumEtConeDR03", &PHOLEAD_hcalDepth1TowerSumEtConeDR03, &b_PHOLEAD_hcalDepth1TowerSumEtConeDR03);
   fChain->SetBranchAddress("PHOLEAD_hcalDepth2TowerSumEtConeDR03", &PHOLEAD_hcalDepth2TowerSumEtConeDR03, &b_PHOLEAD_hcalDepth2TowerSumEtConeDR03);
   fChain->SetBranchAddress("PHOLEAD_trkSumPtSolidConeDR03", &PHOLEAD_trkSumPtSolidConeDR03, &b_PHOLEAD_trkSumPtSolidConeDR03);
   fChain->SetBranchAddress("PHOLEAD_trkSumPtHollowConeDR03", &PHOLEAD_trkSumPtHollowConeDR03, &b_PHOLEAD_trkSumPtHollowConeDR03);
   fChain->SetBranchAddress("PHOLEAD_nTrkSolidConeDR03", &PHOLEAD_nTrkSolidConeDR03, &b_PHOLEAD_nTrkSolidConeDR03);
   fChain->SetBranchAddress("PHOLEAD_nTrkHollowConeDR03", &PHOLEAD_nTrkHollowConeDR03, &b_PHOLEAD_nTrkHollowConeDR03);
   fChain->SetBranchAddress("PHOLEAD_hasConversionTracks", &PHOLEAD_hasConversionTracks, &b_PHOLEAD_hasConversionTracks);
   fChain->SetBranchAddress("PHOLEAD_hasPixelSeed", &PHOLEAD_hasPixelSeed, &b_PHOLEAD_hasPixelSeed);
   fChain->SetBranchAddress("PHOLEAD_nTracks", &PHOLEAD_nTracks, &b_PHOLEAD_nTracks);
   fChain->SetBranchAddress("PHOLEAD_isConverted", &PHOLEAD_isConverted, &b_PHOLEAD_isConverted);
   fChain->SetBranchAddress("PHOLEAD_convPairInvariantMass", &PHOLEAD_convPairInvariantMass, &b_PHOLEAD_convPairInvariantMass);
   fChain->SetBranchAddress("PHOLEAD_convpairCotThetaSeparation", &PHOLEAD_convpairCotThetaSeparation, &b_PHOLEAD_convpairCotThetaSeparation);
   fChain->SetBranchAddress("PHOLEAD_convPairMomentumMag", &PHOLEAD_convPairMomentumMag, &b_PHOLEAD_convPairMomentumMag);
   fChain->SetBranchAddress("PHOLEAD_convPairMomentumPerp", &PHOLEAD_convPairMomentumPerp, &b_PHOLEAD_convPairMomentumPerp);
   fChain->SetBranchAddress("PHOLEAD_convPairMomentumPhi", &PHOLEAD_convPairMomentumPhi, &b_PHOLEAD_convPairMomentumPhi);
   fChain->SetBranchAddress("PHOLEAD_convPairMomentumEta", &PHOLEAD_convPairMomentumEta, &b_PHOLEAD_convPairMomentumEta);
   fChain->SetBranchAddress("PHOLEAD_convPairMomentumX", &PHOLEAD_convPairMomentumX, &b_PHOLEAD_convPairMomentumX);
   fChain->SetBranchAddress("PHOLEAD_convPairMomentumY", &PHOLEAD_convPairMomentumY, &b_PHOLEAD_convPairMomentumY);
   fChain->SetBranchAddress("PHOLEAD_convPairMomentumZ", &PHOLEAD_convPairMomentumZ, &b_PHOLEAD_convPairMomentumZ);
   fChain->SetBranchAddress("PHOLEAD_convDistOfMinimumApproach", &PHOLEAD_convDistOfMinimumApproach, &b_PHOLEAD_convDistOfMinimumApproach);
   fChain->SetBranchAddress("PHOLEAD_convDPhiTracksAtVtx", &PHOLEAD_convDPhiTracksAtVtx, &b_PHOLEAD_convDPhiTracksAtVtx);
   fChain->SetBranchAddress("PHOLEAD_convDPhiTracksAtEcal", &PHOLEAD_convDPhiTracksAtEcal, &b_PHOLEAD_convDPhiTracksAtEcal);
   fChain->SetBranchAddress("PHOLEAD_convDEtaTracksAtEcal", &PHOLEAD_convDEtaTracksAtEcal, &b_PHOLEAD_convDEtaTracksAtEcal);
   fChain->SetBranchAddress("PHOLEAD_convVtxValid", &PHOLEAD_convVtxValid, &b_PHOLEAD_convVtxValid);
   fChain->SetBranchAddress("PHOLEAD_convVtxEta", &PHOLEAD_convVtxEta, &b_PHOLEAD_convVtxEta);
   fChain->SetBranchAddress("PHOLEAD_convVtxPhi", &PHOLEAD_convVtxPhi, &b_PHOLEAD_convVtxPhi);
   fChain->SetBranchAddress("PHOLEAD_convVtxR", &PHOLEAD_convVtxR, &b_PHOLEAD_convVtxR);
   fChain->SetBranchAddress("PHOLEAD_convVtxX", &PHOLEAD_convVtxX, &b_PHOLEAD_convVtxX);
   fChain->SetBranchAddress("PHOLEAD_convVtxY", &PHOLEAD_convVtxY, &b_PHOLEAD_convVtxY);
   fChain->SetBranchAddress("PHOLEAD_convVtxZ", &PHOLEAD_convVtxZ, &b_PHOLEAD_convVtxZ);
   fChain->SetBranchAddress("PHOLEAD_convVtxChi2", &PHOLEAD_convVtxChi2, &b_PHOLEAD_convVtxChi2);
   fChain->SetBranchAddress("PHOLEAD_convVtxNdof", &PHOLEAD_convVtxNdof, &b_PHOLEAD_convVtxNdof);
   fChain->SetBranchAddress("PHOLEAD_convMVALikelihood", &PHOLEAD_convMVALikelihood, &b_PHOLEAD_convMVALikelihood);
   fChain->SetBranchAddress("PHOLEAD_convVtxChi2Prob", &PHOLEAD_convVtxChi2Prob, &b_PHOLEAD_convVtxChi2Prob);
   fChain->SetBranchAddress("PHOLEAD_convEoverP", &PHOLEAD_convEoverP, &b_PHOLEAD_convEoverP);
   fChain->SetBranchAddress("PHOLEAD_convzOfPrimaryVertexFromTracks", &PHOLEAD_convzOfPrimaryVertexFromTracks, &b_PHOLEAD_convzOfPrimaryVertexFromTracks);
   fChain->SetBranchAddress("isGenMatched", &isGenMatched, &b_isGenMatched);
   fChain->SetBranchAddress("genMatchedPt", &genMatchedPt, &b_genMatchedPt);
   fChain->SetBranchAddress("genMatchedEta", &genMatchedEta, &b_genMatchedEta);
   fChain->SetBranchAddress("genMatchedPhi", &genMatchedPhi, &b_genMatchedPhi);
   fChain->SetBranchAddress("motherID", &motherID, &b_motherID);
   Notify();
}

Bool_t reweightTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void reweightTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t reweightTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef reweightTree_cxx
