//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb 18 17:41:20 2010 by ROOT version 5.22/00d
// from TTree photonTree/photonTree
// found on file: ../PJ_15.root
// Author: Shin-Shan Eiko Yu, National Central University, Taiwan
//////////////////////////////////////////////////////////

#ifndef produceHisto_h
#define produceHisto_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string>

using namespace std;

const Int_t MAXNPHOTONS  = 20;
const Int_t MAXNGENMATCH = 9;
const Int_t NMC          = 7;


class produceHisto {
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
   Float_t         momentumX[MAXNPHOTONS];   //[nPhotons]
   Float_t         momentumY[MAXNPHOTONS];   //[nPhotons]
   Float_t         momentumZ[MAXNPHOTONS];   //[nPhotons]
   Float_t         vertexX[MAXNPHOTONS];   //[nPhotons]
   Float_t         vertexY[MAXNPHOTONS];   //[nPhotons]
   Float_t         vertexZ[MAXNPHOTONS];   //[nPhotons]
   Float_t         p[MAXNPHOTONS];   //[nPhotons]
   Float_t         energy[MAXNPHOTONS];   //[nPhotons]
   Float_t         et[MAXNPHOTONS];   //[nPhotons]
   Float_t         caloPositionX[MAXNPHOTONS];   //[nPhotons]
   Float_t         caloPositionY[MAXNPHOTONS];   //[nPhotons]
   Float_t         caloPositionZ[MAXNPHOTONS];   //[nPhotons]
   Float_t         caloPositionEta[MAXNPHOTONS];   //[nPhotons]
   Float_t         caloPositionPhi[MAXNPHOTONS];   //[nPhotons]
   Int_t           isEB[MAXNPHOTONS];   //[nPhotons]
   Int_t           isEE[MAXNPHOTONS];   //[nPhotons]
   Int_t           isEBGap[MAXNPHOTONS];   //[nPhotons]
   Int_t           isEEGap[MAXNPHOTONS];   //[nPhotons]
   Int_t           isEBEEGap[MAXNPHOTONS];   //[nPhotons]
   Float_t         hadronicOverEm[MAXNPHOTONS];   //[nPhotons]
   Float_t         hadronicDepth1OverEm[MAXNPHOTONS];   //[nPhotons]
   Float_t         hadronicDepth2OverEm[MAXNPHOTONS];   //[nPhotons]
   Float_t         e1x5[MAXNPHOTONS];   //[nPhotons]
   Float_t         e2x5[MAXNPHOTONS];   //[nPhotons]
   Float_t         e3x3[MAXNPHOTONS];   //[nPhotons]
   Float_t         maxEnergyXtal[MAXNPHOTONS];   //[nPhotons]
   Float_t         sigmaEtaEta[MAXNPHOTONS];   //[nPhotons]
   Float_t         sigmaIetaIeta[MAXNPHOTONS];   //[nPhotons]
   Float_t         r1x5[MAXNPHOTONS];   //[nPhotons]
   Float_t         r2x5[MAXNPHOTONS];   //[nPhotons]
   Float_t         r9[MAXNPHOTONS];   //[nPhotons]
   Float_t         ecalRecHitSumEtConeDR04[MAXNPHOTONS];   //[nPhotons]
   Float_t         hcalTowerSumEtConeDR04[MAXNPHOTONS];   //[nPhotons]
   Float_t         hcalDepth1TowerSumEtConeDR04[MAXNPHOTONS];   //[nPhotons]
   Float_t         hcalDepth2TowerSumEtConeDR04[MAXNPHOTONS];   //[nPhotons]
   Float_t         trkSumPtSolidConeDR04[MAXNPHOTONS];   //[nPhotons]
   Float_t         trkSumPtHollowConeDR04[MAXNPHOTONS];   //[nPhotons]
   Int_t           nTrkSolidConeDR04[MAXNPHOTONS];   //[nPhotons]
   Int_t           nTrkHollowConeDR04[MAXNPHOTONS];   //[nPhotons]
   Float_t         ecalRecHitSumEtConeDR03[MAXNPHOTONS];   //[nPhotons]
   Float_t         hcalTowerSumEtConeDR03[MAXNPHOTONS];   //[nPhotons]
   Float_t         hcalDepth1TowerSumEtConeDR03[MAXNPHOTONS];   //[nPhotons]
   Float_t         hcalDepth2TowerSumEtConeDR03[MAXNPHOTONS];   //[nPhotons]
   Float_t         trkSumPtSolidConeDR03[MAXNPHOTONS];   //[nPhotons]
   Float_t         trkSumPtHollowConeDR03[MAXNPHOTONS];   //[nPhotons]
   Int_t           nTrkSolidConeDR03[MAXNPHOTONS];   //[nPhotons]
   Int_t           nTrkHollowConeDR03[MAXNPHOTONS];   //[nPhotons]
   Float_t         rawEnergy[MAXNPHOTONS];   //[nPhotons]
   Float_t         preshowerEnergy[MAXNPHOTONS];   //[nPhotons]
   Float_t         phiWidth[MAXNPHOTONS];   //[nPhotons]
   Float_t         etaWidth[MAXNPHOTONS];   //[nPhotons]
   Int_t           clustersSize[MAXNPHOTONS];   //[nPhotons]
   Int_t           nConversions[MAXNPHOTONS];   //[nPhotons]
   Int_t           nGenMatch[MAXNPHOTONS];   //[nPhotons]
   Int_t           genMatchPdgId[MAXNPHOTONS][MAXNGENMATCH];   //[nPhotons]
   Int_t           genMatchStatus[MAXNPHOTONS][MAXNGENMATCH];   //[nPhotons]
   Float_t         genMatchPt[MAXNPHOTONS][MAXNGENMATCH];   //[nPhotons]
   Float_t         genMatchEta[MAXNPHOTONS][MAXNGENMATCH];   //[nPhotons]
   Float_t         genMatchPhi[MAXNPHOTONS][MAXNGENMATCH];   //[nPhotons]

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

   produceHisto(std::string filename);
   std::string    inputFile_;
   virtual ~produceHisto();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   void             SetUseLeadPhoton(bool value){useLeadingPhotonOnly_ = value;}
   void             SetMatching(bool value){matching_ = value;}
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);   
   bool useLeadingPhotonOnly_;
   bool matching_;
};

#endif

#ifdef produceHisto_cxx
produceHisto::produceHisto(std::string filename):
  useLeadingPhotonOnly_(true),
  matching_(true)
{
  inputFile_ = filename;
  TTree* tree;
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename.data());
  if (!f) {
    f = new TFile(filename.data());
  }
  tree = (TTree*)gDirectory->Get("skimLaurent/photonTree");
  if(!tree)
    tree = (TTree*)gDirectory->Get("photonTree");
  Init(tree);
}

produceHisto::~produceHisto()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t produceHisto::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t produceHisto::LoadTree(Long64_t entry)
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

void produceHisto::Init(TTree *tree)
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

   fChain->SetBranchAddress("nAna", &nAna, &b_nAna);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("simVertexX", &simVertexX, &b_simVertexX);
   fChain->SetBranchAddress("simVertexY", &simVertexY, &b_simVertexY);
   fChain->SetBranchAddress("simVertexZ", &simVertexZ, &b_simVertexZ);
   fChain->SetBranchAddress("nPhotons", &nPhotons, &b_nPhotons);
   fChain->SetBranchAddress("ptHat", &ptHat, &b_ptHat);
   fChain->SetBranchAddress("momentumX", momentumX, &b_momentumX);
   fChain->SetBranchAddress("momentumY", momentumY, &b_momentumY);
   fChain->SetBranchAddress("momentumZ", momentumZ, &b_momentumZ);
   fChain->SetBranchAddress("vertexX", vertexX, &b_vertexX);
   fChain->SetBranchAddress("vertexY", vertexY, &b_vertexY);
   fChain->SetBranchAddress("vertexZ", vertexZ, &b_vertexZ);
   fChain->SetBranchAddress("p", p, &b_p);
   fChain->SetBranchAddress("energy", energy, &b_energy);
   fChain->SetBranchAddress("et", et, &b_et);
   fChain->SetBranchAddress("caloPositionX", caloPositionX, &b_caloPositionX);
   fChain->SetBranchAddress("caloPositionY", caloPositionY, &b_caloPositionY);
   fChain->SetBranchAddress("caloPositionZ", caloPositionZ, &b_caloPositionZ);
   fChain->SetBranchAddress("caloPositionEta", caloPositionEta, &b_caloPositionEta);
   fChain->SetBranchAddress("caloPositionPhi", caloPositionPhi, &b_caloPositionPhi);
   fChain->SetBranchAddress("isEB", isEB, &b_isEB);
   fChain->SetBranchAddress("isEE", isEE, &b_isEE);
   fChain->SetBranchAddress("isEBGap", isEBGap, &b_isEBGap);
   fChain->SetBranchAddress("isEEGap", isEEGap, &b_isEEGap);
   fChain->SetBranchAddress("isEBEEGap", isEBEEGap, &b_isEBEEGap);
   fChain->SetBranchAddress("hadronicOverEm", hadronicOverEm, &b_hadronicOverEm);
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
   fChain->SetBranchAddress("r9", r9, &b_r9);
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
   fChain->SetBranchAddress("genMatchPhi", genMatchPhi, &b_genMatchPhi);
   Notify();
}

Bool_t produceHisto::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void produceHisto::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t produceHisto::Cut(Long64_t entry)
{
  return 1;
}
#endif // #ifdef produceHisto_cxx
