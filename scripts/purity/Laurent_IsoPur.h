//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 20 12:08:20 2009 by ROOT version 5.22/00a
// from TChain /
//////////////////////////////////////////////////////////

#ifndef Laurent_IsoPur_h
#define Laurent_IsoPur_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class Laurent_IsoPur {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nAna;
   Int_t           run;
   Int_t           event;
   Float_t         simVertexX;
   Float_t         simVertexY;
   Float_t         simVertexZ;
   Int_t           nPhotons;
   Float_t         ptHat;
   Float_t         momentumX[8];   //[nPhotons]
   Float_t         momentumY[8];   //[nPhotons]
   Float_t         momentumZ[8];   //[nPhotons]
   Float_t         vertexX[8];   //[nPhotons]
   Float_t         vertexY[8];   //[nPhotons]
   Float_t         vertexZ[8];   //[nPhotons]
   Float_t         p[8];   //[nPhotons]
   Float_t         energy[8];   //[nPhotons]
   Float_t         et[8];   //[nPhotons]
   Float_t         caloPositionX[8];   //[nPhotons]
   Float_t         caloPositionY[8];   //[nPhotons]
   Float_t         caloPositionZ[8];   //[nPhotons]
   Float_t         caloPositionEta[8];   //[nPhotons]
   Float_t         caloPositionPhi[8];   //[nPhotons]
   Int_t           isEB[8];   //[nPhotons]
   Int_t           isEE[8];   //[nPhotons]
   Int_t           isEBGap[8];   //[nPhotons]
   Int_t           isEEGap[8];   //[nPhotons]
   Int_t           isEBEEGap[8];   //[nPhotons]
   Float_t         hadronicOverEm[8];   //[nPhotons]
   Float_t         hadronicDepth1OverEm[8];   //[nPhotons]
   Float_t         hadronicDepth2OverEm[8];   //[nPhotons]
   Float_t         e1x5[8];   //[nPhotons]
   Float_t         e2x5[8];   //[nPhotons]
   Float_t         e3x3[8];   //[nPhotons]
   Float_t         maxEnergyXtal[8];   //[nPhotons]
   Float_t         sigmaEtaEta[8];   //[nPhotons]
   Float_t         sigmaIetaIeta[8];   //[nPhotons]
   Float_t         r1x5[8];   //[nPhotons]
   Float_t         r2x5[8];   //[nPhotons]
   Float_t         r9[8];   //[nPhotons]
   Float_t         ecalRecHitSumEtConeDR04[8];   //[nPhotons]
   Float_t         hcalTowerSumEtConeDR04[8];   //[nPhotons]
   Float_t         hcalDepth1TowerSumEtConeDR04[8];   //[nPhotons]
   Float_t         hcalDepth2TowerSumEtConeDR04[8];   //[nPhotons]
   Float_t         trkSumPtSolidConeDR04[8];   //[nPhotons]
   Float_t         trkSumPtHollowConeDR04[8];   //[nPhotons]
   Int_t           nTrkSolidConeDR04[8];   //[nPhotons]
   Int_t           nTrkHollowConeDR04[8];   //[nPhotons]
   Float_t         ecalRecHitSumEtConeDR03[8];   //[nPhotons]
   Float_t         hcalTowerSumEtConeDR03[8];   //[nPhotons]
   Float_t         hcalDepth1TowerSumEtConeDR03[8];   //[nPhotons]
   Float_t         hcalDepth2TowerSumEtConeDR03[8];   //[nPhotons]
   Float_t         trkSumPtSolidConeDR03[8];   //[nPhotons]
   Float_t         trkSumPtHollowConeDR03[8];   //[nPhotons]
   Int_t           nTrkSolidConeDR03[8];   //[nPhotons]
   Int_t           nTrkHollowConeDR03[8];   //[nPhotons]
   Float_t         rawEnergy[8];   //[nPhotons]
   Float_t         preshowerEnergy[8];   //[nPhotons]
   Float_t         phiWidth[8];   //[nPhotons]
   Float_t         etaWidth[8];   //[nPhotons]
   Int_t           clustersSize[8];   //[nPhotons]
   Int_t           nConversions[8];   //[nPhotons]
   Int_t           nGenMatch[8];   //[nPhotons]
   Int_t           genMatchPdgId[8][9];   //[nPhotons]
   Int_t           genMatchStatus[8][9];   //[nPhotons]
   Float_t         genMatchPt[8][9];   //[nPhotons]
   Float_t         genMatchEta[8][9];   //[nPhotons]
   Float_t         genMatchPhi[8][9];   //[nPhotons]

   // List of branches
   TBranch        *b_nAna;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_simVertexX;   //!
   TBranch        *b_simVertexY;   //!
   TBranch        *b_simVertexZ;   //!
   TBranch        *b_nPhotons;   //!
   TBranch        *b_ptHat;   //!
   TBranch        *b_momentumX;   //!
   TBranch        *b_momentumY;   //!
   TBranch        *b_momentumZ;   //!
   TBranch        *b_vertexX;   //!
   TBranch        *b_vertexY;   //!
   TBranch        *b_vertexZ;   //!
   TBranch        *b_p;   //!
   TBranch        *b_energy;   //!
   TBranch        *b_et;   //!
   TBranch        *b_caloPositionX;   //!
   TBranch        *b_caloPositionY;   //!
   TBranch        *b_caloPositionZ;   //!
   TBranch        *b_caloPositionEta;   //!
   TBranch        *b_caloPositionPhi;   //!
   TBranch        *b_isEB;   //!
   TBranch        *b_isEE;   //!
   TBranch        *b_isEBGap;   //!
   TBranch        *b_isEEGap;   //!
   TBranch        *b_isEBEEGap;   //!
   TBranch        *b_hadronicOverEm;   //!
   TBranch        *b_hadronicDepth1OverEm;   //!
   TBranch        *b_hadronicDepth2OverEm;   //!
   TBranch        *b_e1x5;   //!
   TBranch        *b_e2x5;   //!
   TBranch        *b_e3x3;   //!
   TBranch        *b_maxEnergyXtal;   //!
   TBranch        *b_sigmaEtaEta;   //!
   TBranch        *b_sigmaIetaIeta;   //!
   TBranch        *b_r1x5;   //!
   TBranch        *b_r2x5;   //!
   TBranch        *b_r9;   //!
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
   TBranch        *b_rawEnergy;   //!
   TBranch        *b_preshowerEnergy;   //!
   TBranch        *b_phiWidth;   //!
   TBranch        *b_etaWidth;   //!
   TBranch        *b_clustersSize;   //!
   TBranch        *b_nConversions;   //!
   TBranch        *b_nGenMatch;   //!
   TBranch        *b_genMatchPdgId;   //!
   TBranch        *b_genMatchStatus;   //!
   TBranch        *b_genMatchPt;   //!
   TBranch        *b_genMatchEta;   //!
   TBranch        *b_genMatchPhi;   //!

   Laurent_IsoPur(TTree *tree=0);
   virtual ~Laurent_IsoPur();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Int_t choiceEB=1);
   Bool_t IsLoosePhoton(bool isQCD, Int_t ipho);
   Bool_t CheckFiducial(Int_t choiceEB, Int_t ipho);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Double_t purCCerr(Double_t effB, Double_t effB_N, Double_t effS, Double_t effS_N, Double_t effSnB, Double_t effSnB_N);
   virtual Double_t purMCerr(Double_t fillS, Double_t scaleS, Double_t fillB, Double_t scaleB);
};

#endif

#ifdef Laurent_IsoPur_cxx
Laurent_IsoPur::Laurent_IsoPur(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TChain * chain = new TChain("","");
      
      chain->Add("/mc/Laurent/PJ_15_Data.root/photonTree");
      chain->Add("/mc/Laurent/PJ_30_Data.root/photonTree");
      chain->Add("/mc/Laurent/PJ_80_Data.root/photonTree");
      chain->Add("/mc/Laurent/PJ_170_Data.root/photonTree");
      chain->Add("/mc/Laurent/PJ_300_Data.root/photonTree");
      chain->Add("/mc/Laurent/PJ_470_Data.root/photonTree");
      chain->Add("/mc/Laurent/PJ_800_Data.root/photonTree");
      chain->Add("/mc/Laurent/QCD_15_Data.root/photonTree");
      chain->Add("/mc/Laurent/QCD_30_Data.root/photonTree");
      chain->Add("/mc/Laurent/QCD_80_Data.root/photonTree");
      chain->Add("/mc/Laurent/QCD_170_Data.root/photonTree");
      chain->Add("/mc/Laurent/QCD_300_Data.root/photonTree");
      chain->Add("/mc/Laurent/QCD_470_Data.root/photonTree");
      chain->Add("/mc/Laurent/QCD_800_Data.root/photonTree");
      
      chain->Add("/mc/Laurent/PJ_15_Template.root/photonTree");
      chain->Add("/mc/Laurent/PJ_30_Template.root/photonTree");
      chain->Add("/mc/Laurent/PJ_80_Template.root/photonTree");
      chain->Add("/mc/Laurent/PJ_170_Template.root/photonTree");
      chain->Add("/mc/Laurent/PJ_300_Template.root/photonTree");
      chain->Add("/mc/Laurent/PJ_470_Template.root/photonTree");
      chain->Add("/mc/Laurent/PJ_800_Template.root/photonTree");
      chain->Add("/mc/Laurent/QCD_15_Template.root/photonTree");
      chain->Add("/mc/Laurent/QCD_30_Template.root/photonTree");
      chain->Add("/mc/Laurent/QCD_80_Template.root/photonTree");
      chain->Add("/mc/Laurent/QCD_170_Template.root/photonTree");
      chain->Add("/mc/Laurent/QCD_300_Template.root/photonTree");
      chain->Add("/mc/Laurent/QCD_470_Template.root/photonTree");
      chain->Add("/mc/Laurent/QCD_800_Template.root/photonTree");
      
      tree = chain;
   }
   Init(tree);
}

Laurent_IsoPur::~Laurent_IsoPur()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Laurent_IsoPur::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Laurent_IsoPur::LoadTree(Long64_t entry)
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

void Laurent_IsoPur::Init(TTree *tree)
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

   fChain->SetBranchAddress("caloPositionEta", caloPositionEta, &b_caloPositionEta);
   fChain->SetBranchAddress("caloPositionPhi", caloPositionPhi, &b_caloPositionPhi);
   fChain->SetBranchAddress("ecalRecHitSumEtConeDR04", ecalRecHitSumEtConeDR04, &b_ecalRecHitSumEtConeDR04);
   fChain->SetBranchAddress("hcalTowerSumEtConeDR04", hcalTowerSumEtConeDR04, &b_hcalTowerSumEtConeDR04);
   fChain->SetBranchAddress("trkSumPtSolidConeDR04", trkSumPtSolidConeDR04, &b_trkSumPtSolidConeDR04);
   fChain->SetBranchAddress("trkSumPtHollowConeDR04", trkSumPtHollowConeDR04, &b_trkSumPtHollowConeDR04);
   fChain->SetBranchAddress("ecalRecHitSumEtConeDR03", ecalRecHitSumEtConeDR03, &b_ecalRecHitSumEtConeDR03);
   fChain->SetBranchAddress("hcalTowerSumEtConeDR03", hcalTowerSumEtConeDR03, &b_hcalTowerSumEtConeDR03);
   fChain->SetBranchAddress("r9", r9, &b_r9);
   fChain->SetBranchAddress("trkSumPtSolidConeDR03", trkSumPtSolidConeDR03, &b_trkSumPtSolidConeDR03);
   fChain->SetBranchAddress("trkSumPtHollowConeDR03", trkSumPtHollowConeDR03, &b_trkSumPtHollowConeDR03);
   fChain->SetBranchAddress("ptHat", &ptHat, &b_ptHat);
   fChain->SetBranchAddress("isEB", isEB, &b_isEB);
   fChain->SetBranchAddress("isEE", isEE, &b_isEE);
   fChain->SetBranchAddress("et", et, &b_et);

   fChain->SetBranchAddress("nPhotons", &nPhotons, &b_nPhotons);
   fChain->SetBranchAddress("isEBGap", isEBGap, &b_isEBGap);
   fChain->SetBranchAddress("isEEGap", isEEGap, &b_isEEGap);
   fChain->SetBranchAddress("isEBEEGap", isEBEEGap, &b_isEBEEGap);
   fChain->SetBranchAddress("hadronicOverEm", hadronicOverEm, &b_hadronicOverEm);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   
   /*fChain->SetBranchAddress("nAna", &nAna, &b_nAna);
   fChain->SetBranchAddress("simVertexX", &simVertexX, &b_simVertexX);
   fChain->SetBranchAddress("simVertexY", &simVertexY, &b_simVertexY);
   fChain->SetBranchAddress("simVertexZ", &simVertexZ, &b_simVertexZ);
   fChain->SetBranchAddress("momentumX", momentumX, &b_momentumX);
   fChain->SetBranchAddress("momentumY", momentumY, &b_momentumY);
   fChain->SetBranchAddress("momentumZ", momentumZ, &b_momentumZ);
   fChain->SetBranchAddress("vertexX", vertexX, &b_vertexX);
   fChain->SetBranchAddress("vertexY", vertexY, &b_vertexY);
   fChain->SetBranchAddress("vertexZ", vertexZ, &b_vertexZ);
   fChain->SetBranchAddress("p", p, &b_p);
   fChain->SetBranchAddress("energy", energy, &b_energy);
   fChain->SetBranchAddress("caloPositionX", caloPositionX, &b_caloPositionX);
   fChain->SetBranchAddress("caloPositionY", caloPositionY, &b_caloPositionY);
   fChain->SetBranchAddress("caloPositionZ", caloPositionZ, &b_caloPositionZ);
   fChain->SetBranchAddress("hadronicDepth1OverEm", hadronicDepth1OverEm, &b_hadronicDepth1OverEm);
   fChain->SetBranchAddress("hadronicDepth2OverEm", hadronicDepth2OverEm, &b_hadronicDepth2OverEm);
   fChain->SetBranchAddress("e1x5", e1x5, &b_e1x5);
   fChain->SetBranchAddress("e2x5", e2x5, &b_e2x5);
   fChain->SetBranchAddress("e3x3", e3x3, &b_e3x3);
   fChain->SetBranchAddress("maxEnergyXtal", maxEnergyXtal, &b_maxEnergyXtal);
   fChain->SetBranchAddress("sigmaEtaEta", sigmaEtaEta, &b_sigmaEtaEta);
   fChain->SetBranchAddress("sigmaIetaIeta", sigmaIetaIeta, &b_sigmaIetaIeta);
   fChain->SetBranchAddress("r1x5", r1x5, &b_r1x5);
   fChain->SetBranchAddress("r2x5", r2x5, &b_r2x5);
   fChain->SetBranchAddress("hcalDepth1TowerSumEtConeDR04", hcalDepth1TowerSumEtConeDR04, &b_hcalDepth1TowerSumEtConeDR04);
   fChain->SetBranchAddress("hcalDepth2TowerSumEtConeDR04", hcalDepth2TowerSumEtConeDR04, &b_hcalDepth2TowerSumEtConeDR04);
   fChain->SetBranchAddress("nTrkSolidConeDR04", nTrkSolidConeDR04, &b_nTrkSolidConeDR04);
   fChain->SetBranchAddress("nTrkHollowConeDR04", nTrkHollowConeDR04, &b_nTrkHollowConeDR04);
   fChain->SetBranchAddress("hcalDepth1TowerSumEtConeDR03", hcalDepth1TowerSumEtConeDR03, &b_hcalDepth1TowerSumEtConeDR03);
   fChain->SetBranchAddress("hcalDepth2TowerSumEtConeDR03", hcalDepth2TowerSumEtConeDR03, &b_hcalDepth2TowerSumEtConeDR03);
   fChain->SetBranchAddress("nTrkSolidConeDR03", nTrkSolidConeDR03, &b_nTrkSolidConeDR03);
   fChain->SetBranchAddress("nTrkHollowConeDR03", nTrkHollowConeDR03, &b_nTrkHollowConeDR03);
   fChain->SetBranchAddress("rawEnergy", rawEnergy, &b_rawEnergy);
   fChain->SetBranchAddress("preshowerEnergy", preshowerEnergy, &b_preshowerEnergy);
   fChain->SetBranchAddress("phiWidth", phiWidth, &b_phiWidth);
   fChain->SetBranchAddress("etaWidth", etaWidth, &b_etaWidth);
   fChain->SetBranchAddress("clustersSize", clustersSize, &b_clustersSize);
   fChain->SetBranchAddress("nConversions", nConversions, &b_nConversions);
   fChain->SetBranchAddress("nGenMatch", nGenMatch, &b_nGenMatch);
   fChain->SetBranchAddress("genMatchPdgId", genMatchPdgId, &b_genMatchPdgId);
   fChain->SetBranchAddress("genMatchStatus", genMatchStatus, &b_genMatchStatus);
   fChain->SetBranchAddress("genMatchPt", genMatchPt, &b_genMatchPt);
   fChain->SetBranchAddress("genMatchEta", genMatchEta, &b_genMatchEta);
   fChain->SetBranchAddress("genMatchPhi", genMatchPhi, &b_genMatchPhi);*/
   Notify();
}

Bool_t Laurent_IsoPur::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Laurent_IsoPur::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Laurent_IsoPur::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Laurent_IsoPur_cxx
