//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan 11 14:41:34 2011 by ROOT version 5.22/00d
// from TTree Analysis/Analysis
// found on file: EGData_run2010A.root
//////////////////////////////////////////////////////////

#ifndef LoopEntries_h
#define LoopEntries_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>

class LoopEntries {
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
   Float_t         instDelLumiLS;
   Int_t           nTTBits;
   Int_t           kMaxTTBits;
   UChar_t         TTBit[64];   //[nTTBits]
   UChar_t         TTBit34;
   UChar_t         TTBit40;
   UChar_t         TTBit41;
   UChar_t         TTBit0;
   Int_t           kMaxL1Bits;
   Int_t           nL1Bits;
   UChar_t         L1Bit[128];   //[nL1Bits]
   Int_t           kMaxL1Obj;
   Int_t           nL1EMIso;
   Float_t         l1EMIsoEnergy[4];   //[nL1EMIso]
   Float_t         l1EMIsoEt[4];   //[nL1EMIso]
   Float_t         l1EMIsoEta[4];   //[nL1EMIso]
   Float_t         l1EMIsoEPhi[4];   //[nL1EMIso]
   Int_t           nL1EMnonIso;
   Float_t         l1EMnonIsoEnergy[4];   //[nL1EMnonIso]
   Float_t         l1EMnonIsoEt[4];   //[nL1EMnonIso]
   Float_t         l1EMnonIsoEta[4];   //[nL1EMnonIso]
   Float_t         l1EMnonIsoEPhi[4];   //[nL1EMnonIso]
   UChar_t         HLT_Photon10_L1R;
   UChar_t         HLT_Photon10_Cleaned_L1R;
   UChar_t         HLT_Photon15_L1R;
   UChar_t         HLT_Photon15_Cleaned_L1R;
   UChar_t         HLT_Photon15_LooseEcalIso_L1R;
   UChar_t         HLT_Photon15_LooseEcalIso_Cleaned_L1R;
   UChar_t         HLT_Photon15_TrackIso_L1R;
   UChar_t         HLT_Photon15_TrackIso_Cleaned_L1R;
   UChar_t         HLT_Photon17_Isol_SC17HE_L1R_v1;
   UChar_t         HLT_Photon20_L1R;
   UChar_t         HLT_Photon20_Cleaned_L1R;
   UChar_t         HLT_Photon25_Cleaned_L1R;
   UChar_t         HLT_Photon30_L1R;
   UChar_t         HLT_Photon30_Cleaned_L1R;
   UChar_t         HLT_Photon30_L1R_8E29;
   UChar_t         HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1;
   UChar_t         HLT_Photon35_Isol_Cleaned_L1R_v1;
   UChar_t         HLT_Photon40_CaloId_Cleaned_L1R_v1;
   UChar_t         HLT_Photon40_Isol_Cleaned_L1R_v1;
   UChar_t         HLT_Photon50_L1R;
   UChar_t         HLT_Photon50_Cleaned_L1R;
   UChar_t         HLT_Photon50_Cleaned_L1R_v1;
   UChar_t         HLT_Photon50_NoHE_Cleaned_L1R;
   UChar_t         HLT_Photon70_Cleaned_L1R_v1;
   UChar_t         HLT_Photon70_NoHE_Cleaned_L1R_v1;
   UChar_t         HLT_Photon110_NoHE_Cleaned_L1R_v1;
   UChar_t         HLT_DoublePhoton5_L1R;
   UChar_t         HLT_DoublePhoton10_L1R;
   UChar_t         HLT_DoublePhoton15_L1R;
   UChar_t         HLT_DoublePhoton17_L1R;
   UChar_t         HLT_DoublePhoton20_L1R;
   UChar_t         HLT_DoublePhoton22_L1R_v1;
   Int_t           kMaxTrigFlag;
   Int_t           nHLTBits;
   UChar_t         HLTBit[156];   //[nHLTBits]
   Int_t           HLT_Photon20_Cleaned_L1R_prescale;
   Int_t           HLT_Photon30_Cleaned_L1R_prescale;
   Int_t           HLT_Photon50_Cleaned_L1R_v1_prescale;
   Int_t           nHfTowersP;
   Int_t           nHfTowersN;
   Double_t        beamSpotX;
   Double_t        beamSpotY;
   Double_t        beamSpotZ;
   Int_t           nVtxAll;
   Int_t           nVtxGood;
   Int_t           nVtxNotFake;
   Double_t        vtxX[13];   //[nVtxNotFake]
   Double_t        vtxY[13];   //[nVtxNotFake]
   Double_t        vtxZ[13];   //[nVtxNotFake]
   Double_t        vtxXError[13];   //[nVtxNotFake]
   Double_t        vtxYError[13];   //[nVtxNotFake]
   Double_t        vtxZError[13];   //[nVtxNotFake]
   Int_t           vtxNTrk[13];   //[nVtxNotFake]
   Int_t           vtxNTrkWeight05[13];   //[nVtxNotFake]
   Double_t        vtxChi2[13];   //[nVtxNotFake]
   Double_t        vtxNdof[13];   //[nVtxNotFake]
   Double_t        vtxNormChi2[13];   //[nVtxNotFake]
   Int_t           isMETEmpty;
   Double_t        metEt;
   Double_t        metPx;
   Double_t        metPy;
   Double_t        metPz;
   Double_t        metE_longitudinal;
   Int_t           metNCorrections;
   Float_t         metCorEx;
   Float_t         metCorEy;
   Float_t         metCorSumEt;
   Float_t         metUncorrectedPt;
   Float_t         metUncorrectedPhi;
   UChar_t         metIsCaloMET;
   UChar_t         metIsRecoMET;
   Int_t           nJets;
   Int_t           kMaxJets;
   Float_t         jetPt[20];   //[nJets]
   Float_t         jetE[20];   //[nJets]
   Float_t         jetP[20];   //[nJets]
   Float_t         jetEta[20];   //[nJets]
   Float_t         jetPhi[20];   //[nJets]
   Float_t         bDiscriminatorHighEff[20];   //[nJets]
   Float_t         bDiscriminatorHighPur[20];   //[nJets]
   Float_t         jetCharge[20];   //[nJets]
   Int_t           jetNtrk[20];   //[nJets]
   UChar_t         jetIsCaloJet[20];   //[nJets]
   UChar_t         jetIsPFJet[20];   //[nJets]
   UChar_t         jetIsBasicJet[20];   //[nJets]
   Float_t         jetMaxEInEmTowers[20];   //[nJets]
   Float_t         jetMaxEInHadTowers[20];   //[nJets]
   Float_t         jetEnergyFractionHadronic[20];   //[nJets]
   Float_t         jetEmEnergyFraction[20];   //[nJets]
   Float_t         jetHadEnergyInHB[20];   //[nJets]
   Float_t         jetHadEnergyInHO[20];   //[nJets]
   Float_t         jetHadEnergyInHE[20];   //[nJets]
   Float_t         jetHadEnergyInHF[20];   //[nJets]
   Float_t         jetEmEnergyInEB[20];   //[nJets]
   Float_t         jetEmEnergyInEE[20];   //[nJets]
   Float_t         jetEmEnergyInHF[20];   //[nJets]
   Float_t         jetTowersArea[20];   //[nJets]
   Int_t           jetN90[20];   //[nJets]
   Int_t           jetN60[20];   //[nJets]
   Float_t         jetFHPD[20];   //[nJets]
   Float_t         jetFRBX[20];   //[nJets]
   Float_t         jetN90Hits[20];   //[nJets]
   Int_t           nPhotons;
   Int_t           kMaxPhotons;
   Float_t         p[10];   //[nPhotons]
   Float_t         et[10];   //[nPhotons]
   Float_t         energy[10];   //[nPhotons]
   Float_t         momentumX[10];   //[nPhotons]
   Float_t         momentumY[10];   //[nPhotons]
   Float_t         momentumZ[10];   //[nPhotons]
   Float_t         pt[10];   //[nPhotons]
   Float_t         eta[10];   //[nPhotons]
   Float_t         phi[10];   //[nPhotons]
   Float_t         r9[10];   //[nPhotons]
   UChar_t         isEBGap[10];   //[nPhotons]
   UChar_t         isEEGap[10];   //[nPhotons]
   UChar_t         isEBEEGap[10];   //[nPhotons]
   UChar_t         isTransGap[10];   //[nPhotons]
   UChar_t         isEB[10];   //[nPhotons]
   UChar_t         isEE[10];   //[nPhotons]
   Float_t         rawEnergy[10];   //[nPhotons]
   Float_t         preshowerEnergy[10];   //[nPhotons]
   Int_t           numOfPreshClusters[10];   //[nPhotons]
   Float_t         ESRatio[10];   //[nPhotons]
   Int_t           clustersSize[10];   //[nPhotons]
   Int_t           scSize[10];   //[nPhotons]
   Float_t         phiWidth[10];   //[nPhotons]
   Float_t         etaWidth[10];   //[nPhotons]
   Float_t         scEta[10];   //[nPhotons]
   Float_t         scPhi[10];   //[nPhotons]
   Float_t         maxEnergyXtal[10];   //[nPhotons]
   Float_t         sigmaEtaEta[10];   //[nPhotons]
   Float_t         sigmaIetaIeta[10];   //[nPhotons]
   Float_t         sigmaIetaIphi[10];   //[nPhotons]
   Float_t         sigmaIphiIphi[10];   //[nPhotons]
   Float_t         r1x5[10];   //[nPhotons]
   Float_t         r2x5[10];   //[nPhotons]
   Float_t         e1x5[10];   //[nPhotons]
   Float_t         e2x5[10];   //[nPhotons]
   Float_t         seedTime[10];   //[nPhotons]
   Float_t         seedChi2[10];   //[nPhotons]
   Float_t         seedOutOfTimeChi2[10];   //[nPhotons]
   Int_t           seedRecoFlag[10];   //[nPhotons]
   Int_t           seedSeverity[10];   //[nPhotons]
   Float_t         tRight[10];   //[nPhotons]
   Float_t         tLeft[10];   //[nPhotons]
   Float_t         tTop[10];   //[nPhotons]
   Float_t         tBottom[10];   //[nPhotons]
   Float_t         eMax[10];   //[nPhotons]
   Float_t         e2nd[10];   //[nPhotons]
   Float_t         e2x2[10];   //[nPhotons]
   Float_t         e3x2[10];   //[nPhotons]
   Float_t         e3x3[10];   //[nPhotons]
   Float_t         e4x4[10];   //[nPhotons]
   Float_t         e5x5[10];   //[nPhotons]
   Float_t         e2overe8[10];   //[nPhotons]
   Float_t         e2x5Right[10];   //[nPhotons]
   Float_t         e2x5Left[10];   //[nPhotons]
   Float_t         e2x5Top[10];   //[nPhotons]
   Float_t         e2x5Bottom[10];   //[nPhotons]
   Float_t         eRight[10];   //[nPhotons]
   Float_t         eLeft[10];   //[nPhotons]
   Float_t         eTop[10];   //[nPhotons]
   Float_t         eBottom[10];   //[nPhotons]
   Float_t         covPhiPhi[10];   //[nPhotons]
   Float_t         covEtaPhi[10];   //[nPhotons]
   Float_t         covEtaEta[10];   //[nPhotons]
   Float_t         hadronicOverEm[10];   //[nPhotons]
   Float_t         hadronicDepth1OverEm[10];   //[nPhotons]
   Float_t         hadronicDepth2OverEm[10];   //[nPhotons]
   Float_t         trackIso[10];   //[nPhotons]
   Float_t         caloIso[10];   //[nPhotons]
   Float_t         ecalIso[10];   //[nPhotons]
   Float_t         hcalIso[10];   //[nPhotons]
   Float_t         compTrackIso[10];   //[nPhotons]
   Float_t         compEcalIso[10];   //[nPhotons]
   Float_t         compHcalIso[10];   //[nPhotons]
   Int_t           nPixelHits[10];   //[nPhotons]
   Float_t         c1[10];   //[nPhotons]
   Float_t         c2[10];   //[nPhotons]
   Float_t         c3[10];   //[nPhotons]
   Float_t         c4[10];   //[nPhotons]
   Float_t         c5[10];   //[nPhotons]
   Float_t         r1[10];   //[nPhotons]
   Float_t         r2[10];   //[nPhotons]
   Float_t         r3[10];   //[nPhotons]
   Float_t         r4[10];   //[nPhotons]
   Float_t         r5[10];   //[nPhotons]
   Float_t         cc1[10];   //[nPhotons]
   Float_t         cc2[10];   //[nPhotons]
   Float_t         cc3[10];   //[nPhotons]
   Float_t         cc4[10];   //[nPhotons]
   Float_t         cc5[10];   //[nPhotons]
   Float_t         cr1[10];   //[nPhotons]
   Float_t         cr2[10];   //[nPhotons]
   Float_t         cr3[10];   //[nPhotons]
   Float_t         cr4[10];   //[nPhotons]
   Float_t         cr5[10];   //[nPhotons]
   Float_t         dr11[10];   //[nPhotons]
   Float_t         dr12[10];   //[nPhotons]
   Float_t         dr13[10];   //[nPhotons]
   Float_t         dr14[10];   //[nPhotons]
   Float_t         dr21[10];   //[nPhotons]
   Float_t         dr22[10];   //[nPhotons]
   Float_t         dr23[10];   //[nPhotons]
   Float_t         dr24[10];   //[nPhotons]
   Float_t         dr31[10];   //[nPhotons]
   Float_t         dr32[10];   //[nPhotons]
   Float_t         dr33[10];   //[nPhotons]
   Float_t         dr34[10];   //[nPhotons]
   Float_t         dr41[10];   //[nPhotons]
   Float_t         dr42[10];   //[nPhotons]
   Float_t         dr43[10];   //[nPhotons]
   Float_t         dr44[10];   //[nPhotons]
   Float_t         t11[10];   //[nPhotons]
   Float_t         t12[10];   //[nPhotons]
   Float_t         t13[10];   //[nPhotons]
   Float_t         t14[10];   //[nPhotons]
   Float_t         t21[10];   //[nPhotons]
   Float_t         t22[10];   //[nPhotons]
   Float_t         t23[10];   //[nPhotons]
   Float_t         t24[10];   //[nPhotons]
   Float_t         t31[10];   //[nPhotons]
   Float_t         t32[10];   //[nPhotons]
   Float_t         t33[10];   //[nPhotons]
   Float_t         t34[10];   //[nPhotons]
   Float_t         t41[10];   //[nPhotons]
   Float_t         t42[10];   //[nPhotons]
   Float_t         t43[10];   //[nPhotons]
   Float_t         t44[10];   //[nPhotons]
   Float_t         ecalRecHitSumEtConeDR04[10];   //[nPhotons]
   Float_t         hcalTowerSumEtConeDR04[10];   //[nPhotons]
   Float_t         hcalDepth1TowerSumEtConeDR04[10];   //[nPhotons]
   Float_t         hcalDepth2TowerSumEtConeDR04[10];   //[nPhotons]
   Float_t         trkSumPtSolidConeDR04[10];   //[nPhotons]
   Float_t         trkSumPtHollowConeDR04[10];   //[nPhotons]
   Int_t           nTrkSolidConeDR04[10];   //[nPhotons]
   Int_t           nTrkHollowConeDR04[10];   //[nPhotons]
   Float_t         ecalRecHitSumEtConeDR03[10];   //[nPhotons]
   Float_t         hcalTowerSumEtConeDR03[10];   //[nPhotons]
   Float_t         hcalDepth1TowerSumEtConeDR03[10];   //[nPhotons]
   Float_t         hcalDepth2TowerSumEtConeDR03[10];   //[nPhotons]
   Float_t         trkSumPtSolidConeDR03[10];   //[nPhotons]
   Float_t         trkSumPtHollowConeDR03[10];   //[nPhotons]
   Int_t           nTrkSolidConeDR03[10];   //[nPhotons]
   Int_t           nTrkHollowConeDR03[10];   //[nPhotons]
   UChar_t         hasConversionTracks[10];   //[nPhotons]
   UChar_t         hasPixelSeed[10];   //[nPhotons]
   UChar_t         isLoose[10];   //[nPhotons]
   UChar_t         isTight[10];   //[nPhotons]
   Int_t           nTracks[10];   //[nPhotons]
   UChar_t         isConverted[10];   //[nPhotons]
   Float_t         convPairInvariantMass[10];   //[nPhotons]
   Float_t         convTrack2InnerMomentumPerp[10];   //[nPhotons]
   Float_t         convTrack2InnerMomentumEta[10];   //[nPhotons]
   Float_t         convTrack2InnerMomentumPhi[10];   //[nPhotons]
   Float_t         convTrack2OuterMomentumPerp[10];   //[nPhotons]
   Float_t         convTrack2OuterMomentumEta[10];   //[nPhotons]
   Float_t         convTrack2OuterMomentumPhi[10];   //[nPhotons]
   Float_t         convTrack1InnerMomentumPerp[10];   //[nPhotons]
   Float_t         convTrack1InnerMomentumEta[10];   //[nPhotons]
   Float_t         convTrack1InnerMomentumPhi[10];   //[nPhotons]
   Float_t         convTrack1OuterMomentumPerp[10];   //[nPhotons]
   Float_t         convTrack1OuterMomentumEta[10];   //[nPhotons]
   Float_t         convTrack1OuterMomentumPhi[10];   //[nPhotons]
   Float_t         convpairCotThetaSeparation[10];   //[nPhotons]
   Float_t         convPairMomentumMag[10];   //[nPhotons]
   Float_t         convPairMomentumPerp[10];   //[nPhotons]
   Float_t         convPairMomentumPhi[10];   //[nPhotons]
   Float_t         convPairMomentumEta[10];   //[nPhotons]
   Float_t         convPairMomentumX[10];   //[nPhotons]
   Float_t         convPairMomentumY[10];   //[nPhotons]
   Float_t         convPairMomentumZ[10];   //[nPhotons]
   Float_t         convPairRefittedMomentumMag[10];   //[nPhotons]
   Float_t         convPairRefittedMomentumPerp[10];   //[nPhotons]
   Float_t         convPairRefittedMomentumPhi[10];   //[nPhotons]
   Float_t         convPairRefittedMomentumEta[10];   //[nPhotons]
   Float_t         convPairRefittedMomentumX[10];   //[nPhotons]
   Float_t         convPairRefittedMomentumY[10];   //[nPhotons]
   Float_t         convPairRefittedMomentumZ[10];   //[nPhotons]
   Float_t         convDistOfMinimumApproach[10];   //[nPhotons]
   Float_t         convDPhiTracksAtVtx[10];   //[nPhotons]
   Float_t         convDPhiTracksAtEcal[10];   //[nPhotons]
   Float_t         convDEtaTracksAtEcal[10];   //[nPhotons]
   UChar_t         convVtxValid[10];   //[nPhotons]
   Float_t         convVtxEta[10];   //[nPhotons]
   Float_t         convVtxPhi[10];   //[nPhotons]
   Float_t         convVtxR[10];   //[nPhotons]
   Float_t         convVtxX[10];   //[nPhotons]
   Float_t         convVtxY[10];   //[nPhotons]
   Float_t         convVtxZ[10];   //[nPhotons]
   Float_t         convMVALikelihood[10];   //[nPhotons]
   Float_t         convVtxChi2Prob[10];   //[nPhotons]
   Float_t         convEoverP[10];   //[nPhotons]
   Float_t         convzOfPrimaryVertexFromTracks[10];   //[nPhotons]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_orbit;   //!
   TBranch        *b_bunchCrossing;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_timesec;   //!
   TBranch        *b_instDelLumiLS;   //!
   TBranch        *b_nTTBits;   //!
   TBranch        *b_kMaxTTBits;   //!
   TBranch        *b_TTBit;   //!
   TBranch        *b_TTBit34;   //!
   TBranch        *b_TTBit40;   //!
   TBranch        *b_TTBit41;   //!
   TBranch        *b_TTBit0;   //!
   TBranch        *b_kMaxL1Bits;   //!
   TBranch        *b_nL1Bits;   //!
   TBranch        *b_L1Bit;   //!
   TBranch        *b_kMaxL1Obj;   //!
   TBranch        *b_nL1EMIso;   //!
   TBranch        *b_l1EMIsoEnergy;   //!
   TBranch        *b_l1EMIsoEt;   //!
   TBranch        *b_l1EMIsoEta;   //!
   TBranch        *b_l1EMIsoEPhi;   //!
   TBranch        *b_nL1EMnonIso;   //!
   TBranch        *b_l1EMnonIsoEnergy;   //!
   TBranch        *b_l1EMnonIsoEt;   //!
   TBranch        *b_l1EMnonIsoEta;   //!
   TBranch        *b_l1EMnonIsoEPhi;   //!
   TBranch        *b_HLT_Photon10_L1R;   //!
   TBranch        *b_HLT_Photon10_Cleaned_L1R;   //!
   TBranch        *b_HLT_Photon15_L1R;   //!
   TBranch        *b_HLT_Photon15_Cleaned_L1R;   //!
   TBranch        *b_HLT_Photon15_LooseEcalIso_L1R;   //!
   TBranch        *b_HLT_Photon15_LooseEcalIso_Cleaned_L1R;   //!
   TBranch        *b_HLT_Photon15_TrackIso_L1R;   //!
   TBranch        *b_HLT_Photon15_TrackIso_Cleaned_L1R;   //!
   TBranch        *b_HLT_Photon17_Isol_SC17HE_L1R_v1;   //!
   TBranch        *b_HLT_Photon20_L1R;   //!
   TBranch        *b_HLT_Photon20_Cleaned_L1R;   //!
   TBranch        *b_HLT_Photon25_Cleaned_L1R;   //!
   TBranch        *b_HLT_Photon30_L1R;   //!
   TBranch        *b_HLT_Photon30_Cleaned_L1R;   //!
   TBranch        *b_HLT_Photon30_L1R_8E29;   //!
   TBranch        *b_HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1;   //!
   TBranch        *b_HLT_Photon35_Isol_Cleaned_L1R_v1;   //!
   TBranch        *b_HLT_Photon40_CaloId_Cleaned_L1R_v1;   //!
   TBranch        *b_HLT_Photon40_Isol_Cleaned_L1R_v1;   //!
   TBranch        *b_HLT_Photon50_L1R;   //!
   TBranch        *b_HLT_Photon50_Cleaned_L1R;   //!
   TBranch        *b_HLT_Photon50_Cleaned_L1R_v1;   //!
   TBranch        *b_HLT_Photon50_NoHE_Cleaned_L1R;   //!
   TBranch        *b_HLT_Photon70_Cleaned_L1R_v1;   //!
   TBranch        *b_HLT_Photon70_NoHE_Cleaned_L1R_v1;   //!
   TBranch        *b_HLT_Photon110_NoHE_Cleaned_L1R_v1;   //!
   TBranch        *b_HLT_DoublePhoton5_L1R;   //!
   TBranch        *b_HLT_DoublePhoton10_L1R;   //!
   TBranch        *b_HLT_DoublePhoton15_L1R;   //!
   TBranch        *b_HLT_DoublePhoton17_L1R;   //!
   TBranch        *b_HLT_DoublePhoton20_L1R;   //!
   TBranch        *b_HLT_DoublePhoton22_L1R_v1;   //!
   TBranch        *b_kMaxTrigFlag;   //!
   TBranch        *b_nHLTBits;   //!
   TBranch        *b_HLTBit;   //!
   TBranch        *b_HLT_Photon20_Cleaned_L1R_prescale;   //!
   TBranch        *b_HLT_Photon30_Cleaned_L1R_prescale;   //!
   TBranch        *b_HLT_Photon50_Cleaned_L1R_v1_prescale;   //!
   TBranch        *b_nHfTowersP;   //!
   TBranch        *b_nHfTowersN;   //!
   TBranch        *b_beamSpotX;   //!
   TBranch        *b_beamSpotY;   //!
   TBranch        *b_beamSpotZ;   //!
   TBranch        *b_nVtxAll;   //!
   TBranch        *b_nVtxGood;   //!
   TBranch        *b_nVtxNotFake;   //!
   TBranch        *b_vtxX;   //!
   TBranch        *b_vtxY;   //!
   TBranch        *b_vtxZ;   //!
   TBranch        *b_vtxXError;   //!
   TBranch        *b_vtxYError;   //!
   TBranch        *b_vtxZError;   //!
   TBranch        *b_vtxNTrk;   //!
   TBranch        *b_vtxNTrkWeight05;   //!
   TBranch        *b_vtxChi2;   //!
   TBranch        *b_vtxNdof;   //!
   TBranch        *b_vtxNormChi2;   //!
   TBranch        *b_isMETEmpty;   //!
   TBranch        *b_metEt;   //!
   TBranch        *b_metPx;   //!
   TBranch        *b_metPy;   //!
   TBranch        *b_metPz;   //!
   TBranch        *b_metE_longitudinal;   //!
   TBranch        *b_metNCorrections;   //!
   TBranch        *b_metCorEx;   //!
   TBranch        *b_metCorEy;   //!
   TBranch        *b_metCorSumEt;   //!
   TBranch        *b_metUncorrectedPt;   //!
   TBranch        *b_metUncorrectedPhi;   //!
   TBranch        *b_metIsCaloMET;   //!
   TBranch        *b_metIsRecoMET;   //!
   TBranch        *b_nJets;   //!
   TBranch        *b_kMaxJets;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetE;   //!
   TBranch        *b_jetP;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_bDiscriminatorHighEff;   //!
   TBranch        *b_bDiscriminatorHighPur;   //!
   TBranch        *b_jetCharge;   //!
   TBranch        *b_jetNtrk;   //!
   TBranch        *b_jetIsCaloJet;   //!
   TBranch        *b_jetIsPFJet;   //!
   TBranch        *b_jetIsBasicJet;   //!
   TBranch        *b_jetMaxEInEmTowers;   //!
   TBranch        *b_jetMaxEInHadTowers;   //!
   TBranch        *b_jetEnergyFractionHadronic;   //!
   TBranch        *b_jetEmEnergyFraction;   //!
   TBranch        *b_jetHadEnergyInHB;   //!
   TBranch        *b_jetHadEnergyInHO;   //!
   TBranch        *b_jetHadEnergyInHE;   //!
   TBranch        *b_jetHadEnergyInHF;   //!
   TBranch        *b_jetEmEnergyInEB;   //!
   TBranch        *b_jetEmEnergyInEE;   //!
   TBranch        *b_jetEmEnergyInHF;   //!
   TBranch        *b_jetTowersArea;   //!
   TBranch        *b_jetN90;   //!
   TBranch        *b_jetN60;   //!
   TBranch        *b_jetFHPD;   //!
   TBranch        *b_jetFRBX;   //!
   TBranch        *b_jetN90Hits;   //!
   TBranch        *b_nPhotons;   //!
   TBranch        *b_kMaxPhotons;   //!
   TBranch        *b_p;   //!
   TBranch        *b_et;   //!
   TBranch        *b_energy;   //!
   TBranch        *b_momentumX;   //!
   TBranch        *b_momentumY;   //!
   TBranch        *b_momentumZ;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_r9;   //!
   TBranch        *b_isEBGap;   //!
   TBranch        *b_isEEGap;   //!
   TBranch        *b_isEBEEGap;   //!
   TBranch        *b_isTransGap;   //!
   TBranch        *b_isEB;   //!
   TBranch        *b_isEE;   //!
   TBranch        *b_rawEnergy;   //!
   TBranch        *b_preshowerEnergy;   //!
   TBranch        *b_numOfPreshClusters;   //!
   TBranch        *b_ESRatio;   //!
   TBranch        *b_clustersSize;   //!
   TBranch        *b_scSize;   //!
   TBranch        *b_phiWidth;   //!
   TBranch        *b_etaWidth;   //!
   TBranch        *b_scEta;   //!
   TBranch        *b_scPhi;   //!
   TBranch        *b_maxEnergyXtal;   //!
   TBranch        *b_sigmaEtaEta;   //!
   TBranch        *b_sigmaIetaIeta;   //!
   TBranch        *b_sigmaIetaIphi;   //!
   TBranch        *b_sigmaIphiIphi;   //!
   TBranch        *b_r1x5;   //!
   TBranch        *b_r2x5;   //!
   TBranch        *b_e1x5;   //!
   TBranch        *b_e2x5;   //!
   TBranch        *b_seedTime;   //!
   TBranch        *b_seedChi2;   //!
   TBranch        *b_seedOutOfTimeChi2;   //!
   TBranch        *b_seedRecoFlag;   //!
   TBranch        *b_seedSeverity;   //!
   TBranch        *b_tRight;   //!
   TBranch        *b_tLeft;   //!
   TBranch        *b_tTop;   //!
   TBranch        *b_tBottom;   //!
   TBranch        *b_eMax;   //!
   TBranch        *b_e2nd;   //!
   TBranch        *b_e2x2;   //!
   TBranch        *b_e3x2;   //!
   TBranch        *b_e3x3;   //!
   TBranch        *b_e4x4;   //!
   TBranch        *b_e5x5;   //!
   TBranch        *b_e2overe8;   //!
   TBranch        *b_e2x5Right;   //!
   TBranch        *b_e2x5Left;   //!
   TBranch        *b_e2x5Top;   //!
   TBranch        *b_e2x5Bottom;   //!
   TBranch        *b_eRight;   //!
   TBranch        *b_eLeft;   //!
   TBranch        *b_eTop;   //!
   TBranch        *b_eBottom;   //!
   TBranch        *b_covPhiPhi;   //!
   TBranch        *b_covEtaPhi;   //!
   TBranch        *b_covEtaEta;   //!
   TBranch        *b_hadronicOverEm;   //!
   TBranch        *b_hadronicDepth1OverEm;   //!
   TBranch        *b_hadronicDepth2OverEm;   //!
   TBranch        *b_trackIso;   //!
   TBranch        *b_caloIso;   //!
   TBranch        *b_ecalIso;   //!
   TBranch        *b_hcalIso;   //!
   TBranch        *b_compTrackIso;   //!
   TBranch        *b_compEcalIso;   //!
   TBranch        *b_compHcalIso;   //!
   TBranch        *b_nPixelHits;   //!
   TBranch        *b_c1;   //!
   TBranch        *b_c2;   //!
   TBranch        *b_c3;   //!
   TBranch        *b_c4;   //!
   TBranch        *b_c5;   //!
   TBranch        *b_r1;   //!
   TBranch        *b_r2;   //!
   TBranch        *b_r3;   //!
   TBranch        *b_r4;   //!
   TBranch        *b_r5;   //!
   TBranch        *b_cc1;   //!
   TBranch        *b_cc2;   //!
   TBranch        *b_cc3;   //!
   TBranch        *b_cc4;   //!
   TBranch        *b_cc5;   //!
   TBranch        *b_cr1;   //!
   TBranch        *b_cr2;   //!
   TBranch        *b_cr3;   //!
   TBranch        *b_cr4;   //!
   TBranch        *b_cr5;   //!
   TBranch        *b_dr11;   //!
   TBranch        *b_dr12;   //!
   TBranch        *b_dr13;   //!
   TBranch        *b_dr14;   //!
   TBranch        *b_dr21;   //!
   TBranch        *b_dr22;   //!
   TBranch        *b_dr23;   //!
   TBranch        *b_dr24;   //!
   TBranch        *b_dr31;   //!
   TBranch        *b_dr32;   //!
   TBranch        *b_dr33;   //!
   TBranch        *b_dr34;   //!
   TBranch        *b_dr41;   //!
   TBranch        *b_dr42;   //!
   TBranch        *b_dr43;   //!
   TBranch        *b_dr44;   //!
   TBranch        *b_t11;   //!
   TBranch        *b_t12;   //!
   TBranch        *b_t13;   //!
   TBranch        *b_t14;   //!
   TBranch        *b_t21;   //!
   TBranch        *b_t22;   //!
   TBranch        *b_t23;   //!
   TBranch        *b_t24;   //!
   TBranch        *b_t31;   //!
   TBranch        *b_t32;   //!
   TBranch        *b_t33;   //!
   TBranch        *b_t34;   //!
   TBranch        *b_t41;   //!
   TBranch        *b_t42;   //!
   TBranch        *b_t43;   //!
   TBranch        *b_t44;   //!
   TBranch        *b_ecalRecHitSumEtConeDR04;   //!
   TBranch        *b_hcalTowerSumEtConeDR04;   //!
   TBranch        *b_hcalDepth1TowerSumEtConeDR04;   //!
   TBranch        *b_hcalDepth2TowerSumEtConeDR04;   //!
   TBranch        *b_trkSumPtSolidConeDR04;   //!
   TBranch        *b_trkSumPtHollowConeDR04;   //!
   TBranch        *b_nTrkSolidConeDR04;   //!
   TBranch        *b_nTrkHollowConeDR04;   //!
   TBranch        *b_ecalRecHitSumEtConeDR03;   //!
   TBranch        *b_hcalTowerSumEtConeDR03;   //!
   TBranch        *b_hcalDepth1TowerSumEtConeDR03;   //!
   TBranch        *b_hcalDepth2TowerSumEtConeDR03;   //!
   TBranch        *b_trkSumPtSolidConeDR03;   //!
   TBranch        *b_trkSumPtHollowConeDR03;   //!
   TBranch        *b_nTrkSolidConeDR03;   //!
   TBranch        *b_nTrkHollowConeDR03;   //!
   TBranch        *b_hasConversionTracks;   //!
   TBranch        *b_hasPixelSeed;   //!
   TBranch        *b_isLoose;   //!
   TBranch        *b_isTight;   //!
   TBranch        *b_nTracks;   //!
   TBranch        *b_isConverted;   //!
   TBranch        *b_convPairInvariantMass;   //!
   TBranch        *b_convTrack2InnerMomentumPerp;   //!
   TBranch        *b_convTrack2InnerMomentumEta;   //!
   TBranch        *b_convTrack2InnerMomentumPhi;   //!
   TBranch        *b_convTrack2OuterMomentumPerp;   //!
   TBranch        *b_convTrack2OuterMomentumEta;   //!
   TBranch        *b_convTrack2OuterMomentumPhi;   //!
   TBranch        *b_convTrack1InnerMomentumPerp;   //!
   TBranch        *b_convTrack1InnerMomentumEta;   //!
   TBranch        *b_convTrack1InnerMomentumPhi;   //!
   TBranch        *b_convTrack1OuterMomentumPerp;   //!
   TBranch        *b_convTrack1OuterMomentumEta;   //!
   TBranch        *b_convTrack1OuterMomentumPhi;   //!
   TBranch        *b_convpairCotThetaSeparation;   //!
   TBranch        *b_convPairMomentumMag;   //!
   TBranch        *b_convPairMomentumPerp;   //!
   TBranch        *b_convPairMomentumPhi;   //!
   TBranch        *b_convPairMomentumEta;   //!
   TBranch        *b_convPairMomentumX;   //!
   TBranch        *b_convPairMomentumY;   //!
   TBranch        *b_convPairMomentumZ;   //!
   TBranch        *b_convPairRefittedMomentumMag;   //!
   TBranch        *b_convPairRefittedMomentumPerp;   //!
   TBranch        *b_convPairRefittedMomentumPhi;   //!
   TBranch        *b_convPairRefittedMomentumEta;   //!
   TBranch        *b_convPairRefittedMomentumX;   //!
   TBranch        *b_convPairRefittedMomentumY;   //!
   TBranch        *b_convPairRefittedMomentumZ;   //!
   TBranch        *b_convDistOfMinimumApproach;   //!
   TBranch        *b_convDPhiTracksAtVtx;   //!
   TBranch        *b_convDPhiTracksAtEcal;   //!
   TBranch        *b_convDEtaTracksAtEcal;   //!
   TBranch        *b_convVtxValid;   //!
   TBranch        *b_convVtxEta;   //!
   TBranch        *b_convVtxPhi;   //!
   TBranch        *b_convVtxR;   //!
   TBranch        *b_convVtxX;   //!
   TBranch        *b_convVtxY;   //!
   TBranch        *b_convVtxZ;   //!
   TBranch        *b_convMVALikelihood;   //!
   TBranch        *b_convVtxChi2Prob;   //!
   TBranch        *b_convEoverP;   //!
   TBranch        *b_convzOfPrimaryVertexFromTracks;   //!

   LoopEntries(TTree *tree=0);
   virtual ~LoopEntries();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(std::string outputFile="/home/syu/386_QCD_NewData/output.dat");
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef LoopEntries_cxx
LoopEntries::LoopEntries(TTree *tree)
{
  FILE *fTable = fopen("inputFile.txt","r");
  TChain* tempChain = new TChain("NTuples/Analysis");
  int flag=1;   
  int nfile=0;
  char filename[100];
  char tmp[1000];
  while (flag!=-1){
    // first reading input file
    flag=fscanf(fTable,"%s",filename);
    std::string tempFile = filename;

    bool isData      = false;

    if(tempFile.find("EGData")    != std::string::npos)isData     =true;

    // read in x-section
    flag=fscanf(fTable,"%s",tmp);
    double cross=atof(tmp);
    // read in number of events
    flag=fscanf(fTable,"%s",tmp);
    double nevt=atof(tmp);

    flag=fscanf(fTable,"%s",tmp);
    int ptHatLo=atof(tmp);
     
    flag=fscanf(fTable,"%s",tmp);
    int ptHatHi=atof(tmp);

    if (flag!=-1) {
      if(isData) 
 	{ 
 	  cout << "filling data" << endl;
	  tempChain->Add(tempFile.data());
	  cout << "adding " << tempFile << endl;
 	} 
      nfile++; 
    }
	
  }

  tree = tempChain;
  Init(tree);


}

LoopEntries::~LoopEntries()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t LoopEntries::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t LoopEntries::LoadTree(Long64_t entry)
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

void LoopEntries::Init(TTree *tree)
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
/*    fChain->SetBranchAddress("orbit", &orbit, &b_orbit); */
/*    fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_bunchCrossing); */
/*    fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock); */
/*    fChain->SetBranchAddress("timesec", &timesec, &b_timesec); */
/*    fChain->SetBranchAddress("instDelLumiLS", &instDelLumiLS, &b_instDelLumiLS); */
/*    fChain->SetBranchAddress("nTTBits", &nTTBits, &b_nTTBits); */
/*    fChain->SetBranchAddress("kMaxTTBits", &kMaxTTBits, &b_kMaxTTBits); */
/*    fChain->SetBranchAddress("TTBit", TTBit, &b_TTBit); */
/*    fChain->SetBranchAddress("TTBit34", &TTBit34, &b_TTBit34); */
/*    fChain->SetBranchAddress("TTBit40", &TTBit40, &b_TTBit40); */
/*    fChain->SetBranchAddress("TTBit41", &TTBit41, &b_TTBit41); */
/*    fChain->SetBranchAddress("TTBit0", &TTBit0, &b_TTBit0); */
/*    fChain->SetBranchAddress("kMaxL1Bits", &kMaxL1Bits, &b_kMaxL1Bits); */
/*    fChain->SetBranchAddress("nL1Bits", &nL1Bits, &b_nL1Bits); */
/*    fChain->SetBranchAddress("L1Bit", L1Bit, &b_L1Bit); */
/*    fChain->SetBranchAddress("kMaxL1Obj", &kMaxL1Obj, &b_kMaxL1Obj); */
/*    fChain->SetBranchAddress("nL1EMIso", &nL1EMIso, &b_nL1EMIso); */
/*    fChain->SetBranchAddress("l1EMIsoEnergy", l1EMIsoEnergy, &b_l1EMIsoEnergy); */
/*    fChain->SetBranchAddress("l1EMIsoEt", l1EMIsoEt, &b_l1EMIsoEt); */
/*    fChain->SetBranchAddress("l1EMIsoEta", l1EMIsoEta, &b_l1EMIsoEta); */
/*    fChain->SetBranchAddress("l1EMIsoEPhi", l1EMIsoEPhi, &b_l1EMIsoEPhi); */
/*    fChain->SetBranchAddress("nL1EMnonIso", &nL1EMnonIso, &b_nL1EMnonIso); */
/*    fChain->SetBranchAddress("l1EMnonIsoEnergy", l1EMnonIsoEnergy, &b_l1EMnonIsoEnergy); */
/*    fChain->SetBranchAddress("l1EMnonIsoEt", l1EMnonIsoEt, &b_l1EMnonIsoEt); */
/*    fChain->SetBranchAddress("l1EMnonIsoEta", l1EMnonIsoEta, &b_l1EMnonIsoEta); */
/*    fChain->SetBranchAddress("l1EMnonIsoEPhi", l1EMnonIsoEPhi, &b_l1EMnonIsoEPhi); */
   fChain->SetBranchAddress("HLT_Photon10_L1R", &HLT_Photon10_L1R, &b_HLT_Photon10_L1R);
   fChain->SetBranchAddress("HLT_Photon10_Cleaned_L1R", &HLT_Photon10_Cleaned_L1R, &b_HLT_Photon10_Cleaned_L1R);
   fChain->SetBranchAddress("HLT_Photon15_L1R", &HLT_Photon15_L1R, &b_HLT_Photon15_L1R);
   fChain->SetBranchAddress("HLT_Photon15_Cleaned_L1R", &HLT_Photon15_Cleaned_L1R, &b_HLT_Photon15_Cleaned_L1R);
   fChain->SetBranchAddress("HLT_Photon15_LooseEcalIso_L1R", &HLT_Photon15_LooseEcalIso_L1R, &b_HLT_Photon15_LooseEcalIso_L1R);
   fChain->SetBranchAddress("HLT_Photon15_LooseEcalIso_Cleaned_L1R", &HLT_Photon15_LooseEcalIso_Cleaned_L1R, &b_HLT_Photon15_LooseEcalIso_Cleaned_L1R);
   fChain->SetBranchAddress("HLT_Photon15_TrackIso_L1R", &HLT_Photon15_TrackIso_L1R, &b_HLT_Photon15_TrackIso_L1R);
   fChain->SetBranchAddress("HLT_Photon15_TrackIso_Cleaned_L1R", &HLT_Photon15_TrackIso_Cleaned_L1R, &b_HLT_Photon15_TrackIso_Cleaned_L1R);
   fChain->SetBranchAddress("HLT_Photon17_Isol_SC17HE_L1R_v1", &HLT_Photon17_Isol_SC17HE_L1R_v1, &b_HLT_Photon17_Isol_SC17HE_L1R_v1);
   fChain->SetBranchAddress("HLT_Photon20_L1R", &HLT_Photon20_L1R, &b_HLT_Photon20_L1R);
   fChain->SetBranchAddress("HLT_Photon20_Cleaned_L1R", &HLT_Photon20_Cleaned_L1R, &b_HLT_Photon20_Cleaned_L1R);
   fChain->SetBranchAddress("HLT_Photon25_Cleaned_L1R", &HLT_Photon25_Cleaned_L1R, &b_HLT_Photon25_Cleaned_L1R);
   fChain->SetBranchAddress("HLT_Photon30_L1R", &HLT_Photon30_L1R, &b_HLT_Photon30_L1R);
   fChain->SetBranchAddress("HLT_Photon30_Cleaned_L1R", &HLT_Photon30_Cleaned_L1R, &b_HLT_Photon30_Cleaned_L1R);
   fChain->SetBranchAddress("HLT_Photon30_L1R_8E29", &HLT_Photon30_L1R_8E29, &b_HLT_Photon30_L1R_8E29);
   fChain->SetBranchAddress("HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1", &HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1, &b_HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1);
   fChain->SetBranchAddress("HLT_Photon35_Isol_Cleaned_L1R_v1", &HLT_Photon35_Isol_Cleaned_L1R_v1, &b_HLT_Photon35_Isol_Cleaned_L1R_v1);
   fChain->SetBranchAddress("HLT_Photon40_CaloId_Cleaned_L1R_v1", &HLT_Photon40_CaloId_Cleaned_L1R_v1, &b_HLT_Photon40_CaloId_Cleaned_L1R_v1);
   fChain->SetBranchAddress("HLT_Photon40_Isol_Cleaned_L1R_v1", &HLT_Photon40_Isol_Cleaned_L1R_v1, &b_HLT_Photon40_Isol_Cleaned_L1R_v1);
   fChain->SetBranchAddress("HLT_Photon50_L1R", &HLT_Photon50_L1R, &b_HLT_Photon50_L1R);
   fChain->SetBranchAddress("HLT_Photon50_Cleaned_L1R", &HLT_Photon50_Cleaned_L1R, &b_HLT_Photon50_Cleaned_L1R);
   fChain->SetBranchAddress("HLT_Photon50_Cleaned_L1R_v1", &HLT_Photon50_Cleaned_L1R_v1, &b_HLT_Photon50_Cleaned_L1R_v1);
   fChain->SetBranchAddress("HLT_Photon50_NoHE_Cleaned_L1R", &HLT_Photon50_NoHE_Cleaned_L1R, &b_HLT_Photon50_NoHE_Cleaned_L1R);
   fChain->SetBranchAddress("HLT_Photon70_Cleaned_L1R_v1", &HLT_Photon70_Cleaned_L1R_v1, &b_HLT_Photon70_Cleaned_L1R_v1);
   fChain->SetBranchAddress("HLT_Photon70_NoHE_Cleaned_L1R_v1", &HLT_Photon70_NoHE_Cleaned_L1R_v1, &b_HLT_Photon70_NoHE_Cleaned_L1R_v1);
   fChain->SetBranchAddress("HLT_Photon110_NoHE_Cleaned_L1R_v1", &HLT_Photon110_NoHE_Cleaned_L1R_v1, &b_HLT_Photon110_NoHE_Cleaned_L1R_v1);
   fChain->SetBranchAddress("HLT_DoublePhoton5_L1R", &HLT_DoublePhoton5_L1R, &b_HLT_DoublePhoton5_L1R);
   fChain->SetBranchAddress("HLT_DoublePhoton10_L1R", &HLT_DoublePhoton10_L1R, &b_HLT_DoublePhoton10_L1R);
   fChain->SetBranchAddress("HLT_DoublePhoton15_L1R", &HLT_DoublePhoton15_L1R, &b_HLT_DoublePhoton15_L1R);
   fChain->SetBranchAddress("HLT_DoublePhoton17_L1R", &HLT_DoublePhoton17_L1R, &b_HLT_DoublePhoton17_L1R);
   fChain->SetBranchAddress("HLT_DoublePhoton20_L1R", &HLT_DoublePhoton20_L1R, &b_HLT_DoublePhoton20_L1R);
   fChain->SetBranchAddress("HLT_DoublePhoton22_L1R_v1", &HLT_DoublePhoton22_L1R_v1, &b_HLT_DoublePhoton22_L1R_v1);
   fChain->SetBranchAddress("kMaxTrigFlag", &kMaxTrigFlag, &b_kMaxTrigFlag);
   fChain->SetBranchAddress("nHLTBits", &nHLTBits, &b_nHLTBits);
   fChain->SetBranchAddress("HLTBit", HLTBit, &b_HLTBit);
   fChain->SetBranchAddress("HLT_Photon20_Cleaned_L1R_prescale", &HLT_Photon20_Cleaned_L1R_prescale, &b_HLT_Photon20_Cleaned_L1R_prescale);
   fChain->SetBranchAddress("HLT_Photon30_Cleaned_L1R_prescale", &HLT_Photon30_Cleaned_L1R_prescale, &b_HLT_Photon30_Cleaned_L1R_prescale);
   fChain->SetBranchAddress("HLT_Photon50_Cleaned_L1R_v1_prescale", &HLT_Photon50_Cleaned_L1R_v1_prescale, &b_HLT_Photon50_Cleaned_L1R_v1_prescale);
   fChain->SetBranchAddress("nHfTowersP", &nHfTowersP, &b_nHfTowersP);
   fChain->SetBranchAddress("nHfTowersN", &nHfTowersN, &b_nHfTowersN);
   fChain->SetBranchAddress("beamSpotX", &beamSpotX, &b_beamSpotX);
   fChain->SetBranchAddress("beamSpotY", &beamSpotY, &b_beamSpotY);
   fChain->SetBranchAddress("beamSpotZ", &beamSpotZ, &b_beamSpotZ);
   fChain->SetBranchAddress("nVtxAll", &nVtxAll, &b_nVtxAll);
   fChain->SetBranchAddress("nVtxGood", &nVtxGood, &b_nVtxGood);
   fChain->SetBranchAddress("nVtxNotFake", &nVtxNotFake, &b_nVtxNotFake);
   fChain->SetBranchAddress("vtxX", vtxX, &b_vtxX);
   fChain->SetBranchAddress("vtxY", vtxY, &b_vtxY);
   fChain->SetBranchAddress("vtxZ", vtxZ, &b_vtxZ);
   fChain->SetBranchAddress("vtxXError", vtxXError, &b_vtxXError);
   fChain->SetBranchAddress("vtxYError", vtxYError, &b_vtxYError);
   fChain->SetBranchAddress("vtxZError", vtxZError, &b_vtxZError);
   fChain->SetBranchAddress("vtxNTrk", vtxNTrk, &b_vtxNTrk);
   fChain->SetBranchAddress("vtxNTrkWeight05", vtxNTrkWeight05, &b_vtxNTrkWeight05);
   fChain->SetBranchAddress("vtxChi2", vtxChi2, &b_vtxChi2);
   fChain->SetBranchAddress("vtxNdof", vtxNdof, &b_vtxNdof);
   fChain->SetBranchAddress("vtxNormChi2", vtxNormChi2, &b_vtxNormChi2);
   fChain->SetBranchAddress("isMETEmpty", &isMETEmpty, &b_isMETEmpty);
/*    fChain->SetBranchAddress("metEt", &metEt, &b_metEt); */
/*    fChain->SetBranchAddress("metPx", &metPx, &b_metPx); */
/*    fChain->SetBranchAddress("metPy", &metPy, &b_metPy); */
/*    fChain->SetBranchAddress("metPz", &metPz, &b_metPz); */
/*    fChain->SetBranchAddress("metE_longitudinal", &metE_longitudinal, &b_metE_longitudinal); */
/*    fChain->SetBranchAddress("metNCorrections", &metNCorrections, &b_metNCorrections); */
/*    fChain->SetBranchAddress("metCorEx", &metCorEx, &b_metCorEx); */
/*    fChain->SetBranchAddress("metCorEy", &metCorEy, &b_metCorEy); */
/*    fChain->SetBranchAddress("metCorSumEt", &metCorSumEt, &b_metCorSumEt); */
/*    fChain->SetBranchAddress("metUncorrectedPt", &metUncorrectedPt, &b_metUncorrectedPt); */
/*    fChain->SetBranchAddress("metUncorrectedPhi", &metUncorrectedPhi, &b_metUncorrectedPhi); */
/*    fChain->SetBranchAddress("metIsCaloMET", &metIsCaloMET, &b_metIsCaloMET); */
/*    fChain->SetBranchAddress("metIsRecoMET", &metIsRecoMET, &b_metIsRecoMET); */
/*    fChain->SetBranchAddress("nJets", &nJets, &b_nJets); */
/*    fChain->SetBranchAddress("kMaxJets", &kMaxJets, &b_kMaxJets); */
/*    fChain->SetBranchAddress("jetPt", jetPt, &b_jetPt); */
/*    fChain->SetBranchAddress("jetE", jetE, &b_jetE); */
/*    fChain->SetBranchAddress("jetP", jetP, &b_jetP); */
/*    fChain->SetBranchAddress("jetEta", jetEta, &b_jetEta); */
/*    fChain->SetBranchAddress("jetPhi", jetPhi, &b_jetPhi); */
/*    fChain->SetBranchAddress("bDiscriminatorHighEff", bDiscriminatorHighEff, &b_bDiscriminatorHighEff); */
/*    fChain->SetBranchAddress("bDiscriminatorHighPur", bDiscriminatorHighPur, &b_bDiscriminatorHighPur); */
/*    fChain->SetBranchAddress("jetCharge", jetCharge, &b_jetCharge); */
/*    fChain->SetBranchAddress("jetNtrk", jetNtrk, &b_jetNtrk); */
/*    fChain->SetBranchAddress("jetIsCaloJet", jetIsCaloJet, &b_jetIsCaloJet); */
/*    fChain->SetBranchAddress("jetIsPFJet", jetIsPFJet, &b_jetIsPFJet); */
/*    fChain->SetBranchAddress("jetIsBasicJet", jetIsBasicJet, &b_jetIsBasicJet); */
/*    fChain->SetBranchAddress("jetMaxEInEmTowers", jetMaxEInEmTowers, &b_jetMaxEInEmTowers); */
/*    fChain->SetBranchAddress("jetMaxEInHadTowers", jetMaxEInHadTowers, &b_jetMaxEInHadTowers); */
/*    fChain->SetBranchAddress("jetEnergyFractionHadronic", jetEnergyFractionHadronic, &b_jetEnergyFractionHadronic); */
/*    fChain->SetBranchAddress("jetEmEnergyFraction", jetEmEnergyFraction, &b_jetEmEnergyFraction); */
/*    fChain->SetBranchAddress("jetHadEnergyInHB", jetHadEnergyInHB, &b_jetHadEnergyInHB); */
/*    fChain->SetBranchAddress("jetHadEnergyInHO", jetHadEnergyInHO, &b_jetHadEnergyInHO); */
/*    fChain->SetBranchAddress("jetHadEnergyInHE", jetHadEnergyInHE, &b_jetHadEnergyInHE); */
/*    fChain->SetBranchAddress("jetHadEnergyInHF", jetHadEnergyInHF, &b_jetHadEnergyInHF); */
/*    fChain->SetBranchAddress("jetEmEnergyInEB", jetEmEnergyInEB, &b_jetEmEnergyInEB); */
/*    fChain->SetBranchAddress("jetEmEnergyInEE", jetEmEnergyInEE, &b_jetEmEnergyInEE); */
/*    fChain->SetBranchAddress("jetEmEnergyInHF", jetEmEnergyInHF, &b_jetEmEnergyInHF); */
/*    fChain->SetBranchAddress("jetTowersArea", jetTowersArea, &b_jetTowersArea); */
/*    fChain->SetBranchAddress("jetN90", jetN90, &b_jetN90); */
/*    fChain->SetBranchAddress("jetN60", jetN60, &b_jetN60); */
/*    fChain->SetBranchAddress("jetFHPD", jetFHPD, &b_jetFHPD); */
/*    fChain->SetBranchAddress("jetFRBX", jetFRBX, &b_jetFRBX); */
/*    fChain->SetBranchAddress("jetN90Hits", jetN90Hits, &b_jetN90Hits); */
   fChain->SetBranchAddress("nPhotons", &nPhotons, &b_nPhotons);
   fChain->SetBranchAddress("kMaxPhotons", &kMaxPhotons, &b_kMaxPhotons);
   fChain->SetBranchAddress("p", p, &b_p);
   fChain->SetBranchAddress("et", et, &b_et);
   fChain->SetBranchAddress("energy", energy, &b_energy);
   fChain->SetBranchAddress("momentumX", momentumX, &b_momentumX);
   fChain->SetBranchAddress("momentumY", momentumY, &b_momentumY);
   fChain->SetBranchAddress("momentumZ", momentumZ, &b_momentumZ);
   fChain->SetBranchAddress("pt", pt, &b_pt);
   fChain->SetBranchAddress("eta", eta, &b_eta);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("r9", r9, &b_r9);
   fChain->SetBranchAddress("isEBGap", isEBGap, &b_isEBGap);
   fChain->SetBranchAddress("isEEGap", isEEGap, &b_isEEGap);
   fChain->SetBranchAddress("isEBEEGap", isEBEEGap, &b_isEBEEGap);
   fChain->SetBranchAddress("isTransGap", isTransGap, &b_isTransGap);
   fChain->SetBranchAddress("isEB", isEB, &b_isEB);
   fChain->SetBranchAddress("isEE", isEE, &b_isEE);
   fChain->SetBranchAddress("rawEnergy", rawEnergy, &b_rawEnergy);
   fChain->SetBranchAddress("preshowerEnergy", preshowerEnergy, &b_preshowerEnergy);
   fChain->SetBranchAddress("numOfPreshClusters", numOfPreshClusters, &b_numOfPreshClusters);
   fChain->SetBranchAddress("ESRatio", ESRatio, &b_ESRatio);
   fChain->SetBranchAddress("clustersSize", clustersSize, &b_clustersSize);
   fChain->SetBranchAddress("scSize", scSize, &b_scSize);
   fChain->SetBranchAddress("phiWidth", phiWidth, &b_phiWidth);
   fChain->SetBranchAddress("etaWidth", etaWidth, &b_etaWidth);
   fChain->SetBranchAddress("scEta", scEta, &b_scEta);
   fChain->SetBranchAddress("scPhi", scPhi, &b_scPhi);
   fChain->SetBranchAddress("maxEnergyXtal", maxEnergyXtal, &b_maxEnergyXtal);
   fChain->SetBranchAddress("sigmaEtaEta", sigmaEtaEta, &b_sigmaEtaEta);
   fChain->SetBranchAddress("sigmaIetaIeta", sigmaIetaIeta, &b_sigmaIetaIeta);
   fChain->SetBranchAddress("sigmaIetaIphi", sigmaIetaIphi, &b_sigmaIetaIphi);
   fChain->SetBranchAddress("sigmaIphiIphi", sigmaIphiIphi, &b_sigmaIphiIphi);
   fChain->SetBranchAddress("r1x5", r1x5, &b_r1x5);
   fChain->SetBranchAddress("r2x5", r2x5, &b_r2x5);
   fChain->SetBranchAddress("e1x5", e1x5, &b_e1x5);
   fChain->SetBranchAddress("e2x5", e2x5, &b_e2x5);
   fChain->SetBranchAddress("seedTime", seedTime, &b_seedTime);
   fChain->SetBranchAddress("seedChi2", seedChi2, &b_seedChi2);
   fChain->SetBranchAddress("seedOutOfTimeChi2", seedOutOfTimeChi2, &b_seedOutOfTimeChi2);
   fChain->SetBranchAddress("seedRecoFlag", seedRecoFlag, &b_seedRecoFlag);
   fChain->SetBranchAddress("seedSeverity", seedSeverity, &b_seedSeverity);
   fChain->SetBranchAddress("tRight", tRight, &b_tRight);
   fChain->SetBranchAddress("tLeft", tLeft, &b_tLeft);
   fChain->SetBranchAddress("tTop", tTop, &b_tTop);
   fChain->SetBranchAddress("tBottom", tBottom, &b_tBottom);
   fChain->SetBranchAddress("eMax", eMax, &b_eMax);
   fChain->SetBranchAddress("e2nd", e2nd, &b_e2nd);
   fChain->SetBranchAddress("e2x2", e2x2, &b_e2x2);
   fChain->SetBranchAddress("e3x2", e3x2, &b_e3x2);
   fChain->SetBranchAddress("e3x3", e3x3, &b_e3x3);
   fChain->SetBranchAddress("e4x4", e4x4, &b_e4x4);
   fChain->SetBranchAddress("e5x5", e5x5, &b_e5x5);
   fChain->SetBranchAddress("e2overe8", e2overe8, &b_e2overe8);
   fChain->SetBranchAddress("e2x5Right", e2x5Right, &b_e2x5Right);
   fChain->SetBranchAddress("e2x5Left", e2x5Left, &b_e2x5Left);
   fChain->SetBranchAddress("e2x5Top", e2x5Top, &b_e2x5Top);
   fChain->SetBranchAddress("e2x5Bottom", e2x5Bottom, &b_e2x5Bottom);
   fChain->SetBranchAddress("eRight", eRight, &b_eRight);
   fChain->SetBranchAddress("eLeft", eLeft, &b_eLeft);
   fChain->SetBranchAddress("eTop", eTop, &b_eTop);
   fChain->SetBranchAddress("eBottom", eBottom, &b_eBottom);
   fChain->SetBranchAddress("covPhiPhi", covPhiPhi, &b_covPhiPhi);
   fChain->SetBranchAddress("covEtaPhi", covEtaPhi, &b_covEtaPhi);
   fChain->SetBranchAddress("covEtaEta", covEtaEta, &b_covEtaEta);
   fChain->SetBranchAddress("hadronicOverEm", hadronicOverEm, &b_hadronicOverEm);
   fChain->SetBranchAddress("hadronicDepth1OverEm", hadronicDepth1OverEm, &b_hadronicDepth1OverEm);
   fChain->SetBranchAddress("hadronicDepth2OverEm", hadronicDepth2OverEm, &b_hadronicDepth2OverEm);
   fChain->SetBranchAddress("trackIso", trackIso, &b_trackIso);
   fChain->SetBranchAddress("caloIso", caloIso, &b_caloIso);
   fChain->SetBranchAddress("ecalIso", ecalIso, &b_ecalIso);
   fChain->SetBranchAddress("hcalIso", hcalIso, &b_hcalIso);
   fChain->SetBranchAddress("compTrackIso", compTrackIso, &b_compTrackIso);
   fChain->SetBranchAddress("compEcalIso", compEcalIso, &b_compEcalIso);
   fChain->SetBranchAddress("compHcalIso", compHcalIso, &b_compHcalIso);
   fChain->SetBranchAddress("nPixelHits", nPixelHits, &b_nPixelHits);
/*    fChain->SetBranchAddress("c1", c1, &b_c1); */
/*    fChain->SetBranchAddress("c2", c2, &b_c2); */
/*    fChain->SetBranchAddress("c3", c3, &b_c3); */
/*    fChain->SetBranchAddress("c4", c4, &b_c4); */
/*    fChain->SetBranchAddress("c5", c5, &b_c5); */
/*    fChain->SetBranchAddress("r1", r1, &b_r1); */
/*    fChain->SetBranchAddress("r2", r2, &b_r2); */
/*    fChain->SetBranchAddress("r3", r3, &b_r3); */
/*    fChain->SetBranchAddress("r4", r4, &b_r4); */
/*    fChain->SetBranchAddress("r5", r5, &b_r5); */
/*    fChain->SetBranchAddress("cc1", cc1, &b_cc1); */
/*    fChain->SetBranchAddress("cc2", cc2, &b_cc2); */
/*    fChain->SetBranchAddress("cc3", cc3, &b_cc3); */
/*    fChain->SetBranchAddress("cc4", cc4, &b_cc4); */
/*    fChain->SetBranchAddress("cc5", cc5, &b_cc5); */
/*    fChain->SetBranchAddress("cr1", cr1, &b_cr1); */
/*    fChain->SetBranchAddress("cr2", cr2, &b_cr2); */
/*    fChain->SetBranchAddress("cr3", cr3, &b_cr3); */
/*    fChain->SetBranchAddress("cr4", cr4, &b_cr4); */
/*    fChain->SetBranchAddress("cr5", cr5, &b_cr5); */
/*    fChain->SetBranchAddress("dr11", dr11, &b_dr11); */
/*    fChain->SetBranchAddress("dr12", dr12, &b_dr12); */
/*    fChain->SetBranchAddress("dr13", dr13, &b_dr13); */
/*    fChain->SetBranchAddress("dr14", dr14, &b_dr14); */
/*    fChain->SetBranchAddress("dr21", dr21, &b_dr21); */
/*    fChain->SetBranchAddress("dr22", dr22, &b_dr22); */
/*    fChain->SetBranchAddress("dr23", dr23, &b_dr23); */
/*    fChain->SetBranchAddress("dr24", dr24, &b_dr24); */
/*    fChain->SetBranchAddress("dr31", dr31, &b_dr31); */
/*    fChain->SetBranchAddress("dr32", dr32, &b_dr32); */
/*    fChain->SetBranchAddress("dr33", dr33, &b_dr33); */
/*    fChain->SetBranchAddress("dr34", dr34, &b_dr34); */
/*    fChain->SetBranchAddress("dr41", dr41, &b_dr41); */
/*    fChain->SetBranchAddress("dr42", dr42, &b_dr42); */
/*    fChain->SetBranchAddress("dr43", dr43, &b_dr43); */
/*    fChain->SetBranchAddress("dr44", dr44, &b_dr44); */
/*    fChain->SetBranchAddress("t11", t11, &b_t11); */
/*    fChain->SetBranchAddress("t12", t12, &b_t12); */
/*    fChain->SetBranchAddress("t13", t13, &b_t13); */
/*    fChain->SetBranchAddress("t14", t14, &b_t14); */
/*    fChain->SetBranchAddress("t21", t21, &b_t21); */
/*    fChain->SetBranchAddress("t22", t22, &b_t22); */
/*    fChain->SetBranchAddress("t23", t23, &b_t23); */
/*    fChain->SetBranchAddress("t24", t24, &b_t24); */
/*    fChain->SetBranchAddress("t31", t31, &b_t31); */
/*    fChain->SetBranchAddress("t32", t32, &b_t32); */
/*    fChain->SetBranchAddress("t33", t33, &b_t33); */
/*    fChain->SetBranchAddress("t34", t34, &b_t34); */
/*    fChain->SetBranchAddress("t41", t41, &b_t41); */
/*    fChain->SetBranchAddress("t42", t42, &b_t42); */
/*    fChain->SetBranchAddress("t43", t43, &b_t43); */
/*    fChain->SetBranchAddress("t44", t44, &b_t44); */
   fChain->SetBranchAddress("ecalRecHitSumEtConeDR04", ecalRecHitSumEtConeDR04, &b_ecalRecHitSumEtConeDR04);
   fChain->SetBranchAddress("hcalTowerSumEtConeDR04", hcalTowerSumEtConeDR04, &b_hcalTowerSumEtConeDR04);
   fChain->SetBranchAddress("hcalDepth1TowerSumEtConeDR04", hcalDepth1TowerSumEtConeDR04, &b_hcalDepth1TowerSumEtConeDR04);
   fChain->SetBranchAddress("hcalDepth2TowerSumEtConeDR04", hcalDepth2TowerSumEtConeDR04, &b_hcalDepth2TowerSumEtConeDR04);
   fChain->SetBranchAddress("trkSumPtSolidConeDR04", trkSumPtSolidConeDR04, &b_trkSumPtSolidConeDR04);
   fChain->SetBranchAddress("trkSumPtHollowConeDR04", trkSumPtHollowConeDR04, &b_trkSumPtHollowConeDR04);
   fChain->SetBranchAddress("nTrkSolidConeDR04", nTrkSolidConeDR04, &b_nTrkSolidConeDR04);
   fChain->SetBranchAddress("nTrkHollowConeDR04", nTrkHollowConeDR04, &b_nTrkHollowConeDR04);
   fChain->SetBranchAddress("ecalRecHitSumEtConeDR03", ecalRecHitSumEtConeDR03, &b_ecalRecHitSumEtConeDR03);
   fChain->SetBranchAddress("hcalTowerSumEtConeDR03", hcalTowerSumEtConeDR03, &b_hcalTowerSumEtConeDR03);
   fChain->SetBranchAddress("hcalDepth1TowerSumEtConeDR03", hcalDepth1TowerSumEtConeDR03, &b_hcalDepth1TowerSumEtConeDR03);
   fChain->SetBranchAddress("hcalDepth2TowerSumEtConeDR03", hcalDepth2TowerSumEtConeDR03, &b_hcalDepth2TowerSumEtConeDR03);
   fChain->SetBranchAddress("trkSumPtSolidConeDR03", trkSumPtSolidConeDR03, &b_trkSumPtSolidConeDR03);
   fChain->SetBranchAddress("trkSumPtHollowConeDR03", trkSumPtHollowConeDR03, &b_trkSumPtHollowConeDR03);
   fChain->SetBranchAddress("nTrkSolidConeDR03", nTrkSolidConeDR03, &b_nTrkSolidConeDR03);
   fChain->SetBranchAddress("nTrkHollowConeDR03", nTrkHollowConeDR03, &b_nTrkHollowConeDR03);
   fChain->SetBranchAddress("hasConversionTracks", hasConversionTracks, &b_hasConversionTracks);
   fChain->SetBranchAddress("hasPixelSeed", hasPixelSeed, &b_hasPixelSeed);
   fChain->SetBranchAddress("isLoose", isLoose, &b_isLoose);
   fChain->SetBranchAddress("isTight", isTight, &b_isTight);
   fChain->SetBranchAddress("nTracks", nTracks, &b_nTracks);
   fChain->SetBranchAddress("isConverted", isConverted, &b_isConverted);
/*    fChain->SetBranchAddress("convPairInvariantMass", convPairInvariantMass, &b_convPairInvariantMass); */
/*    fChain->SetBranchAddress("convTrack2InnerMomentumPerp", convTrack2InnerMomentumPerp, &b_convTrack2InnerMomentumPerp); */
/*    fChain->SetBranchAddress("convTrack2InnerMomentumEta", convTrack2InnerMomentumEta, &b_convTrack2InnerMomentumEta); */
/*    fChain->SetBranchAddress("convTrack2InnerMomentumPhi", convTrack2InnerMomentumPhi, &b_convTrack2InnerMomentumPhi); */
/*    fChain->SetBranchAddress("convTrack2OuterMomentumPerp", convTrack2OuterMomentumPerp, &b_convTrack2OuterMomentumPerp); */
/*    fChain->SetBranchAddress("convTrack2OuterMomentumEta", convTrack2OuterMomentumEta, &b_convTrack2OuterMomentumEta); */
/*    fChain->SetBranchAddress("convTrack2OuterMomentumPhi", convTrack2OuterMomentumPhi, &b_convTrack2OuterMomentumPhi); */
/*    fChain->SetBranchAddress("convTrack1InnerMomentumPerp", convTrack1InnerMomentumPerp, &b_convTrack1InnerMomentumPerp); */
/*    fChain->SetBranchAddress("convTrack1InnerMomentumEta", convTrack1InnerMomentumEta, &b_convTrack1InnerMomentumEta); */
/*    fChain->SetBranchAddress("convTrack1InnerMomentumPhi", convTrack1InnerMomentumPhi, &b_convTrack1InnerMomentumPhi); */
/*    fChain->SetBranchAddress("convTrack1OuterMomentumPerp", convTrack1OuterMomentumPerp, &b_convTrack1OuterMomentumPerp); */
/*    fChain->SetBranchAddress("convTrack1OuterMomentumEta", convTrack1OuterMomentumEta, &b_convTrack1OuterMomentumEta); */
/*    fChain->SetBranchAddress("convTrack1OuterMomentumPhi", convTrack1OuterMomentumPhi, &b_convTrack1OuterMomentumPhi); */
/*    fChain->SetBranchAddress("convpairCotThetaSeparation", convpairCotThetaSeparation, &b_convpairCotThetaSeparation); */
/*    fChain->SetBranchAddress("convPairMomentumMag", convPairMomentumMag, &b_convPairMomentumMag); */
/*    fChain->SetBranchAddress("convPairMomentumPerp", convPairMomentumPerp, &b_convPairMomentumPerp); */
/*    fChain->SetBranchAddress("convPairMomentumPhi", convPairMomentumPhi, &b_convPairMomentumPhi); */
/*    fChain->SetBranchAddress("convPairMomentumEta", convPairMomentumEta, &b_convPairMomentumEta); */
/*    fChain->SetBranchAddress("convPairMomentumX", convPairMomentumX, &b_convPairMomentumX); */
/*    fChain->SetBranchAddress("convPairMomentumY", convPairMomentumY, &b_convPairMomentumY); */
/*    fChain->SetBranchAddress("convPairMomentumZ", convPairMomentumZ, &b_convPairMomentumZ); */
/*    fChain->SetBranchAddress("convPairRefittedMomentumMag", convPairRefittedMomentumMag, &b_convPairRefittedMomentumMag); */
/*    fChain->SetBranchAddress("convPairRefittedMomentumPerp", convPairRefittedMomentumPerp, &b_convPairRefittedMomentumPerp); */
/*    fChain->SetBranchAddress("convPairRefittedMomentumPhi", convPairRefittedMomentumPhi, &b_convPairRefittedMomentumPhi); */
/*    fChain->SetBranchAddress("convPairRefittedMomentumEta", convPairRefittedMomentumEta, &b_convPairRefittedMomentumEta); */
/*    fChain->SetBranchAddress("convPairRefittedMomentumX", convPairRefittedMomentumX, &b_convPairRefittedMomentumX); */
/*    fChain->SetBranchAddress("convPairRefittedMomentumY", convPairRefittedMomentumY, &b_convPairRefittedMomentumY); */
/*    fChain->SetBranchAddress("convPairRefittedMomentumZ", convPairRefittedMomentumZ, &b_convPairRefittedMomentumZ); */
/*    fChain->SetBranchAddress("convDistOfMinimumApproach", convDistOfMinimumApproach, &b_convDistOfMinimumApproach); */
/*    fChain->SetBranchAddress("convDPhiTracksAtVtx", convDPhiTracksAtVtx, &b_convDPhiTracksAtVtx); */
/*    fChain->SetBranchAddress("convDPhiTracksAtEcal", convDPhiTracksAtEcal, &b_convDPhiTracksAtEcal); */
/*    fChain->SetBranchAddress("convDEtaTracksAtEcal", convDEtaTracksAtEcal, &b_convDEtaTracksAtEcal); */
/*    fChain->SetBranchAddress("convVtxValid", convVtxValid, &b_convVtxValid); */
/*    fChain->SetBranchAddress("convVtxEta", convVtxEta, &b_convVtxEta); */
/*    fChain->SetBranchAddress("convVtxPhi", convVtxPhi, &b_convVtxPhi); */
/*    fChain->SetBranchAddress("convVtxR", convVtxR, &b_convVtxR); */
/*    fChain->SetBranchAddress("convVtxX", convVtxX, &b_convVtxX); */
/*    fChain->SetBranchAddress("convVtxY", convVtxY, &b_convVtxY); */
/*    fChain->SetBranchAddress("convVtxZ", convVtxZ, &b_convVtxZ); */
/*    fChain->SetBranchAddress("convMVALikelihood", convMVALikelihood, &b_convMVALikelihood); */
/*    fChain->SetBranchAddress("convVtxChi2Prob", convVtxChi2Prob, &b_convVtxChi2Prob); */
/*    fChain->SetBranchAddress("convEoverP", convEoverP, &b_convEoverP); */
/*    fChain->SetBranchAddress("convzOfPrimaryVertexFromTracks", convzOfPrimaryVertexFromTracks, &b_convzOfPrimaryVertexFromTracks); */
   Notify();
}

Bool_t LoopEntries::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void LoopEntries::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t LoopEntries::Cut(Long64_t entry)
{
  
   return 1;
}
#endif // #ifdef LoopEntries_cxx
