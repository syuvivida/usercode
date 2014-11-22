//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Nov 22 16:48:17 2014 by ROOT version 5.34/00
// from TTree tree/tree
// found on file: /data1/delpanjNtuples/delpanj_v4_WW_pythia_filtered.root
//////////////////////////////////////////////////////////

#ifndef skimTree_h
#define skimTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
// Header file for the classes stored in the TTree if any.
#include <vector>
using namespace std;
// Fixed size dimensions of array or collections stored in the TTree if any.

class skimTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         pu_nTrueInt;
   Int_t           pu_nPUVert;
   Int_t           info_isData;
   Int_t           info_eventId;
   Int_t           info_runId;
   Int_t           info_lumiSection;
   Int_t           info_bunchXing;
   Int_t           info_nVtx;
   vector<float>   *info_vx;
   vector<float>   *info_vy;
   vector<float>   *info_vz;
   Float_t         ptHat;
   Float_t         mcWeight;
   Int_t           nGenPar;
   vector<float>   *genParE;
   vector<float>   *genParPt;
   vector<float>   *genParEta;
   vector<float>   *genParPhi;
   vector<float>   *genParM;
   vector<int>     *genParQ;
   vector<int>     *genParId;
   vector<int>     *genParSt;
   vector<int>     *genMomParId;
   vector<int>     *genParIndex;
   vector<int>     *genNMo;
   vector<int>     *genNDa;
   vector<int>     *genMo1;
   vector<int>     *genMo2;
   vector<int>     *genDa1;
   vector<int>     *genDa2;
   Int_t           nGenJet;
   vector<float>   *genJetE;
   vector<float>   *genJetPt;
   vector<float>   *genJetEta;
   vector<float>   *genJetPhi;
   Int_t           nMu;
   vector<float>   *muPt;
   vector<float>   *muEta;
   vector<float>   *muPhi;
   vector<float>   *muM;
   vector<float>   *muTrkIso;
   vector<float>   *muCorrTrkIso;
   vector<float>   *muHcalIso;
   vector<float>   *muEcalIso;
   vector<float>   *muCharge;
   vector<int>     *muPassID;
   vector<float>   *muChHadIso;
   vector<float>   *muNeHadIso;
   vector<float>   *muGamIso;
   vector<float>   *muPUPt;
   vector<float>   *muCorrPfIso;
   vector<float>   *muTauCorrPfIso;
   vector<int>     *isGlobalMuon;
   vector<int>     *isTrackerMuon;
   vector<float>   *muPtErrx;
   vector<float>   *mudxy;
   vector<float>   *mudz;
   vector<int>     *muTrkLayers;
   vector<int>     *muPixelHits;
   vector<int>     *muHits;
   vector<int>     *muMatches;
   Int_t           nEle;
   Float_t         eleRho;
   vector<float>   *eleEt;
   vector<float>   *eleEnergy;
   vector<float>   *elePt;
   vector<float>   *eleEta;
   vector<float>   *elePhi;
   vector<float>   *eleM;
   vector<float>   *eleScEta;
   vector<float>   *eleSigIhIh;
   vector<float>   *eleDelEtaIn;
   vector<float>   *eleDelPhiIn;
   vector<float>   *eleHoE;
   vector<float>   *eleTrkIso;
   vector<float>   *eleHcalIso;
   vector<float>   *eleEcalIso;
   vector<float>   *eleUserTrkIso;
   vector<float>   *eleUserCalIso;
   vector<float>   *eleCorrPfIso;
   vector<float>   *eleTauCorrPfIso;
   vector<float>   *eleEoverP;
   vector<float>   *eleDxy;
   vector<float>   *eleDz;
   vector<float>   *eleChHadIso;
   vector<float>   *eleNeHadIso;
   vector<float>   *eleGamIso;
   vector<float>   *eleMissingHits;
   vector<float>   *eleInBarrel;
   vector<float>   *eleInEndcap;
   vector<int>     *elePassConv;
   vector<int>     *elePassID;
   Float_t         pfMetCorrPt;
   Float_t         pfMetCorrPhi;
   Float_t         pfMetCorrSumEt;
   Float_t         pfMetCorrSig;
   Float_t         pfMetRawPt;
   Float_t         pfMetRawPhi;
   Float_t         pfMetRawSumEt;
   Float_t         pfMetRawCov00;
   Float_t         pfMetRawCov01;
   Float_t         pfMetRawCov10;
   Float_t         pfMetRawCov11;
   Int_t           CA8nJet;
   vector<float>   *CA8jetPt;
   vector<float>   *CA8jetEta;
   vector<float>   *CA8jetPhi;
   vector<float>   *CA8jetMass;
   vector<float>   *CA8jetEn;
   vector<float>   *CA8jetCorrUncUp;
   vector<float>   *CA8jetCorrUncDown;
   vector<int>     *CA8jetCharge;
   vector<int>     *CA8jetPartonFlavor;
   vector<int>     *CA8jetPassID;
   vector<float>   *CA8jetSSV;
   vector<float>   *CA8jetCSV;
   vector<float>   *CA8jetTCHP;
   vector<float>   *CA8jetTCHE;
   vector<float>   *CA8jetJP;
   vector<float>   *CA8jetJBP;
   vector<float>   *CA8jetTau1;
   vector<float>   *CA8jetTau2;
   vector<float>   *CA8jetTau3;
   vector<float>   *CA8jetTau4;
   vector<float>   *CA8jetMuEF;
   vector<float>   *CA8jetPhoEF;
   vector<float>   *CA8jetCEmEF;
   vector<float>   *CA8jetCHadEF;
   vector<float>   *CA8jetNEmEF;
   vector<float>   *CA8jetNHadEF;
   vector<float>   *CA8jetCMulti;
   vector<float>   *CA8jetPrunedPt;
   vector<float>   *CA8jetPrunedEta;
   vector<float>   *CA8jetPrunedPhi;
   vector<float>   *CA8jetPrunedMass;
   vector<float>   *CA8jetPrunedEn;
   vector<float>   *CA8jetPrunedCorrUncUp;
   vector<float>   *CA8jetPrunedCorrUncDown;
   vector<int>     *CA8jetPrunedCharge;
   vector<int>     *CA8jetPrunedPartonFlavor;
   vector<float>   *CA8jetPrunedSSV;
   vector<float>   *CA8jetPrunedCSV;
   vector<float>   *CA8jetPrunedTCHP;
   vector<float>   *CA8jetPrunedTCHE;
   vector<float>   *CA8jetPrunedJP;
   vector<float>   *CA8jetPrunedJBP;
   vector<int>     *CA8nSubPrunedJet;
   vector<vector<int> > *CA8subjetMotherIndex;
   vector<vector<float> > *CA8subjetPrunedPt;
   vector<vector<float> > *CA8subjetPrunedEta;
   vector<vector<float> > *CA8subjetPrunedPhi;
   vector<vector<float> > *CA8subjetPrunedMass;
   vector<vector<float> > *CA8subjetPrunedEn;
   vector<vector<int> > *CA8subjetPrunedCharge;
   vector<vector<int> > *CA8subjetPrunedPartonFlavor;
   vector<vector<float> > *CA8subjetPrunedCSV;
   Int_t           AK5nJet;
   vector<float>   *AK5jetPt;
   vector<float>   *AK5jetEta;
   vector<float>   *AK5jetPhi;
   vector<float>   *AK5jetMass;
   vector<float>   *AK5jetEn;
   vector<float>   *AK5jetCorrUncUp;
   vector<float>   *AK5jetCorrUncDown;
   vector<int>     *AK5jetCharge;
   vector<int>     *AK5jetPartonFlavor;
   vector<int>     *AK5jetPassID;
   vector<float>   *AK5jetSSV;
   vector<float>   *AK5jetCSV;
   vector<float>   *AK5jetTCHP;
   vector<float>   *AK5jetTCHE;
   vector<float>   *AK5jetJP;
   vector<float>   *AK5jetJBP;
   vector<float>   *AK5jetTau1;
   vector<float>   *AK5jetTau2;
   vector<float>   *AK5jetTau3;
   vector<float>   *AK5jetTau4;
   vector<float>   *AK5jetMuEF;
   vector<float>   *AK5jetPhoEF;
   vector<float>   *AK5jetCEmEF;
   vector<float>   *AK5jetCHadEF;
   vector<float>   *AK5jetNEmEF;
   vector<float>   *AK5jetNHadEF;
   vector<float>   *AK5jetCMulti;
   Int_t           hlt_nTrigs;
   vector<int>     *hlt_trigResult;
   vector<string>  *hlt_trigName;

   // List of branches
   TBranch        *b_pu_nTrueInt;   //!
   TBranch        *b_pu_nPUVert;   //!
   TBranch        *b_info_isData;   //!
   TBranch        *b_info_eventId;   //!
   TBranch        *b_info_runId;   //!
   TBranch        *b_info_lumiSection;   //!
   TBranch        *b_info_bunchXing;   //!
   TBranch        *b_info_nVtx;   //!
   TBranch        *b_info_vx;   //!
   TBranch        *b_info_vy;   //!
   TBranch        *b_info_vz;   //!
   TBranch        *b_ptHat;   //!
   TBranch        *b_mcWeight;   //!
   TBranch        *b_nGenPar;   //!
   TBranch        *b_genParE;   //!
   TBranch        *b_genParPt;   //!
   TBranch        *b_genParEta;   //!
   TBranch        *b_genParPhi;   //!
   TBranch        *b_genParM;   //!
   TBranch        *b_genParQ;   //!
   TBranch        *b_genParId;   //!
   TBranch        *b_genParSt;   //!
   TBranch        *b_genMomParId;   //!
   TBranch        *b_genParIndex;   //!
   TBranch        *b_genNMo;   //!
   TBranch        *b_genNDa;   //!
   TBranch        *b_genMo1;   //!
   TBranch        *b_genMo2;   //!
   TBranch        *b_genDa1;   //!
   TBranch        *b_genDa2;   //!
   TBranch        *b_nGenJet;   //!
   TBranch        *b_genJetE;   //!
   TBranch        *b_genJetPt;   //!
   TBranch        *b_genJetEta;   //!
   TBranch        *b_genJetPhi;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_muPt;   //!
   TBranch        *b_muEta;   //!
   TBranch        *b_muPhi;   //!
   TBranch        *b_muM;   //!
   TBranch        *b_muTrkIso;   //!
   TBranch        *b_muCorrTrkIso;   //!
   TBranch        *b_muHcalIso;   //!
   TBranch        *b_muEcalIso;   //!
   TBranch        *b_muCharge;   //!
   TBranch        *b_muPassID;   //!
   TBranch        *b_muChHadIso;   //!
   TBranch        *b_muNeHadIso;   //!
   TBranch        *b_muGamIso;   //!
   TBranch        *b_muPUPt;   //!
   TBranch        *b_muCorrPfIso;   //!
   TBranch        *b_muTauCorrPfIso;   //!
   TBranch        *b_isGlobalMuon;   //!
   TBranch        *b_isTrackerMuon;   //!
   TBranch        *b_muPtErrx;   //!
   TBranch        *b_mudxy;   //!
   TBranch        *b_mudz;   //!
   TBranch        *b_muTrkLayers;   //!
   TBranch        *b_muPixelHits;   //!
   TBranch        *b_muHits;   //!
   TBranch        *b_muMatches;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_eleRho;   //!
   TBranch        *b_eleEt;   //!
   TBranch        *b_eleEnergy;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleM;   //!
   TBranch        *b_eleScEta;   //!
   TBranch        *b_eleSigIhIh;   //!
   TBranch        *b_eleDelEtaIn;   //!
   TBranch        *b_eleDelPhiIn;   //!
   TBranch        *b_eleHoE;   //!
   TBranch        *b_eleTrkIso;   //!
   TBranch        *b_eleHcalIso;   //!
   TBranch        *b_eleEcalIso;   //!
   TBranch        *b_eleUserTrkIso;   //!
   TBranch        *b_eleUserCalIso;   //!
   TBranch        *b_eleCorrPfIso;   //!
   TBranch        *b_eleTauCorrPfIso;   //!
   TBranch        *b_eleEoverP;   //!
   TBranch        *b_eleDxy;   //!
   TBranch        *b_eleDz;   //!
   TBranch        *b_eleChHadIso;   //!
   TBranch        *b_eleNeHadIso;   //!
   TBranch        *b_eleGamIso;   //!
   TBranch        *b_eleMissingHits;   //!
   TBranch        *b_eleInBarrel;   //!
   TBranch        *b_eleInEndcap;   //!
   TBranch        *b_elePassConv;   //!
   TBranch        *b_elePassID;   //!
   TBranch        *b_pfMetCorrPt;   //!
   TBranch        *b_pfMetCorrPhi;   //!
   TBranch        *b_pfMetCorrSumEt;   //!
   TBranch        *b_pfMetCorrSig;   //!
   TBranch        *b_pfMetRawPt;   //!
   TBranch        *b_pfMetRawPhi;   //!
   TBranch        *b_pfMetRawSumEt;   //!
   TBranch        *b_pfMetRawCov00;   //!
   TBranch        *b_pfMetRawCov01;   //!
   TBranch        *b_pfMetRawCov10;   //!
   TBranch        *b_pfMetRawCov11;   //!
   TBranch        *b_CA8nJet;   //!
   TBranch        *b_CA8jetPt;   //!
   TBranch        *b_CA8jetEta;   //!
   TBranch        *b_CA8jetPhi;   //!
   TBranch        *b_CA8jetMass;   //!
   TBranch        *b_CA8jetEn;   //!
   TBranch        *b_CA8jetCorrUncUp;   //!
   TBranch        *b_CA8jetCorrUncDown;   //!
   TBranch        *b_CA8jetCharge;   //!
   TBranch        *b_CA8jetPartonFlavor;   //!
   TBranch        *b_CA8jetPassID;   //!
   TBranch        *b_CA8jetSSV;   //!
   TBranch        *b_CA8jetCSV;   //!
   TBranch        *b_CA8jetTCHP;   //!
   TBranch        *b_CA8jetTCHE;   //!
   TBranch        *b_CA8jetJP;   //!
   TBranch        *b_CA8jetJBP;   //!
   TBranch        *b_CA8jetTau1;   //!
   TBranch        *b_CA8jetTau2;   //!
   TBranch        *b_CA8jetTau3;   //!
   TBranch        *b_CA8jetTau4;   //!
   TBranch        *b_CA8jetMuEF;   //!
   TBranch        *b_CA8jetPhoEF;   //!
   TBranch        *b_CA8jetCEmEF;   //!
   TBranch        *b_CA8jetCHadEF;   //!
   TBranch        *b_CA8jetNEmEF;   //!
   TBranch        *b_CA8jetNHadEF;   //!
   TBranch        *b_CA8jetCMulti;   //!
   TBranch        *b_CA8jetPrunedPt;   //!
   TBranch        *b_CA8jetPrunedEta;   //!
   TBranch        *b_CA8jetPrunedPhi;   //!
   TBranch        *b_CA8jetPrunedMass;   //!
   TBranch        *b_CA8jetPrunedEn;   //!
   TBranch        *b_CA8jetPrunedCorrUncUp;   //!
   TBranch        *b_CA8jetPrunedCorrUncDown;   //!
   TBranch        *b_CA8jetPrunedCharge;   //!
   TBranch        *b_CA8jetPrunedPartonFlavor;   //!
   TBranch        *b_CA8jetPrunedSSV;   //!
   TBranch        *b_CA8jetPrunedCSV;   //!
   TBranch        *b_CA8jetPrunedTCHP;   //!
   TBranch        *b_CA8jetPrunedTCHE;   //!
   TBranch        *b_CA8jetPrunedJP;   //!
   TBranch        *b_CA8jetPrunedJBP;   //!
   TBranch        *b_CA8nSubPrunedJet;   //!
   TBranch        *b_CA8subjetMotherIndex;   //!
   TBranch        *b_CA8subjetPrunedPt;   //!
   TBranch        *b_CA8subjetPrunedEta;   //!
   TBranch        *b_CA8subjetPrunedPhi;   //!
   TBranch        *b_CA8subjetPrunedMass;   //!
   TBranch        *b_CA8subjetPrunedEn;   //!
   TBranch        *b_CA8subjetPrunedCharge;   //!
   TBranch        *b_CA8subjetPrunedPartonFlavor;   //!
   TBranch        *b_CA8subjetPrunedCSV;   //!
   TBranch        *b_AK5nJet;   //!
   TBranch        *b_AK5jetPt;   //!
   TBranch        *b_AK5jetEta;   //!
   TBranch        *b_AK5jetPhi;   //!
   TBranch        *b_AK5jetMass;   //!
   TBranch        *b_AK5jetEn;   //!
   TBranch        *b_AK5jetCorrUncUp;   //!
   TBranch        *b_AK5jetCorrUncDown;   //!
   TBranch        *b_AK5jetCharge;   //!
   TBranch        *b_AK5jetPartonFlavor;   //!
   TBranch        *b_AK5jetPassID;   //!
   TBranch        *b_AK5jetSSV;   //!
   TBranch        *b_AK5jetCSV;   //!
   TBranch        *b_AK5jetTCHP;   //!
   TBranch        *b_AK5jetTCHE;   //!
   TBranch        *b_AK5jetJP;   //!
   TBranch        *b_AK5jetJBP;   //!
   TBranch        *b_AK5jetTau1;   //!
   TBranch        *b_AK5jetTau2;   //!
   TBranch        *b_AK5jetTau3;   //!
   TBranch        *b_AK5jetTau4;   //!
   TBranch        *b_AK5jetMuEF;   //!
   TBranch        *b_AK5jetPhoEF;   //!
   TBranch        *b_AK5jetCEmEF;   //!
   TBranch        *b_AK5jetCHadEF;   //!
   TBranch        *b_AK5jetNEmEF;   //!
   TBranch        *b_AK5jetNHadEF;   //!
   TBranch        *b_AK5jetCMulti;   //!
   TBranch        *b_hlt_nTrigs;   //!
   TBranch        *b_hlt_trigResult;   //!
   TBranch        *b_hlt_trigName;   //!

   skimTree(std::string inputFile,TTree *tree=0);
   virtual ~skimTree();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   string inputFile_;

};

#endif

#ifdef skimTree_cxx
skimTree::skimTree(std::string inputFile, TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  inputFile_= inputFile;
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(inputFile.data());
    if (!f || !f->IsOpen()) {
      f = new TFile(inputFile.data());
    }
    TDirectory * dir = (TDirectory*)f->Get(Form("%s:/tree",inputFile.data()));
    dir->GetObject("tree",tree);

  }
  Init(tree);
}


skimTree::~skimTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t skimTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t skimTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void skimTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   info_vx = 0;
   info_vy = 0;
   info_vz = 0;
   genParE = 0;
   genParPt = 0;
   genParEta = 0;
   genParPhi = 0;
   genParM = 0;
   genParQ = 0;
   genParId = 0;
   genParSt = 0;
   genMomParId = 0;
   genParIndex = 0;
   genNMo = 0;
   genNDa = 0;
   genMo1 = 0;
   genMo2 = 0;
   genDa1 = 0;
   genDa2 = 0;
   genJetE = 0;
   genJetPt = 0;
   genJetEta = 0;
   genJetPhi = 0;
   muPt = 0;
   muEta = 0;
   muPhi = 0;
   muM = 0;
   muTrkIso = 0;
   muCorrTrkIso = 0;
   muHcalIso = 0;
   muEcalIso = 0;
   muCharge = 0;
   muPassID = 0;
   muChHadIso = 0;
   muNeHadIso = 0;
   muGamIso = 0;
   muPUPt = 0;
   muCorrPfIso = 0;
   muTauCorrPfIso = 0;
   isGlobalMuon = 0;
   isTrackerMuon = 0;
   muPtErrx = 0;
   mudxy = 0;
   mudz = 0;
   muTrkLayers = 0;
   muPixelHits = 0;
   muHits = 0;
   muMatches = 0;
   eleEt = 0;
   eleEnergy = 0;
   elePt = 0;
   eleEta = 0;
   elePhi = 0;
   eleM = 0;
   eleScEta = 0;
   eleSigIhIh = 0;
   eleDelEtaIn = 0;
   eleDelPhiIn = 0;
   eleHoE = 0;
   eleTrkIso = 0;
   eleHcalIso = 0;
   eleEcalIso = 0;
   eleUserTrkIso = 0;
   eleUserCalIso = 0;
   eleCorrPfIso = 0;
   eleTauCorrPfIso = 0;
   eleEoverP = 0;
   eleDxy = 0;
   eleDz = 0;
   eleChHadIso = 0;
   eleNeHadIso = 0;
   eleGamIso = 0;
   eleMissingHits = 0;
   eleInBarrel = 0;
   eleInEndcap = 0;
   elePassConv = 0;
   elePassID = 0;
   CA8jetPt = 0;
   CA8jetEta = 0;
   CA8jetPhi = 0;
   CA8jetMass = 0;
   CA8jetEn = 0;
   CA8jetCorrUncUp = 0;
   CA8jetCorrUncDown = 0;
   CA8jetCharge = 0;
   CA8jetPartonFlavor = 0;
   CA8jetPassID = 0;
   CA8jetSSV = 0;
   CA8jetCSV = 0;
   CA8jetTCHP = 0;
   CA8jetTCHE = 0;
   CA8jetJP = 0;
   CA8jetJBP = 0;
   CA8jetTau1 = 0;
   CA8jetTau2 = 0;
   CA8jetTau3 = 0;
   CA8jetTau4 = 0;
   CA8jetMuEF = 0;
   CA8jetPhoEF = 0;
   CA8jetCEmEF = 0;
   CA8jetCHadEF = 0;
   CA8jetNEmEF = 0;
   CA8jetNHadEF = 0;
   CA8jetCMulti = 0;
   CA8jetPrunedPt = 0;
   CA8jetPrunedEta = 0;
   CA8jetPrunedPhi = 0;
   CA8jetPrunedMass = 0;
   CA8jetPrunedEn = 0;
   CA8jetPrunedCorrUncUp = 0;
   CA8jetPrunedCorrUncDown = 0;
   CA8jetPrunedCharge = 0;
   CA8jetPrunedPartonFlavor = 0;
   CA8jetPrunedSSV = 0;
   CA8jetPrunedCSV = 0;
   CA8jetPrunedTCHP = 0;
   CA8jetPrunedTCHE = 0;
   CA8jetPrunedJP = 0;
   CA8jetPrunedJBP = 0;
   CA8nSubPrunedJet = 0;
   CA8subjetMotherIndex = 0;
   CA8subjetPrunedPt = 0;
   CA8subjetPrunedEta = 0;
   CA8subjetPrunedPhi = 0;
   CA8subjetPrunedMass = 0;
   CA8subjetPrunedEn = 0;
   CA8subjetPrunedCharge = 0;
   CA8subjetPrunedPartonFlavor = 0;
   CA8subjetPrunedCSV = 0;
   AK5jetPt = 0;
   AK5jetEta = 0;
   AK5jetPhi = 0;
   AK5jetMass = 0;
   AK5jetEn = 0;
   AK5jetCorrUncUp = 0;
   AK5jetCorrUncDown = 0;
   AK5jetCharge = 0;
   AK5jetPartonFlavor = 0;
   AK5jetPassID = 0;
   AK5jetSSV = 0;
   AK5jetCSV = 0;
   AK5jetTCHP = 0;
   AK5jetTCHE = 0;
   AK5jetJP = 0;
   AK5jetJBP = 0;
   AK5jetTau1 = 0;
   AK5jetTau2 = 0;
   AK5jetTau3 = 0;
   AK5jetTau4 = 0;
   AK5jetMuEF = 0;
   AK5jetPhoEF = 0;
   AK5jetCEmEF = 0;
   AK5jetCHadEF = 0;
   AK5jetNEmEF = 0;
   AK5jetNHadEF = 0;
   AK5jetCMulti = 0;
   hlt_trigResult = 0;
   hlt_trigName = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("pu_nTrueInt", &pu_nTrueInt, &b_pu_nTrueInt);
   fChain->SetBranchAddress("pu_nPUVert", &pu_nPUVert, &b_pu_nPUVert);
   fChain->SetBranchAddress("info_isData", &info_isData, &b_info_isData);
   fChain->SetBranchAddress("info_eventId", &info_eventId, &b_info_eventId);
   fChain->SetBranchAddress("info_runId", &info_runId, &b_info_runId);
   fChain->SetBranchAddress("info_lumiSection", &info_lumiSection, &b_info_lumiSection);
   fChain->SetBranchAddress("info_bunchXing", &info_bunchXing, &b_info_bunchXing);
   fChain->SetBranchAddress("info_nVtx", &info_nVtx, &b_info_nVtx);
   fChain->SetBranchAddress("info_vx", &info_vx, &b_info_vx);
   fChain->SetBranchAddress("info_vy", &info_vy, &b_info_vy);
   fChain->SetBranchAddress("info_vz", &info_vz, &b_info_vz);
   fChain->SetBranchAddress("ptHat", &ptHat, &b_ptHat);
   fChain->SetBranchAddress("mcWeight", &mcWeight, &b_mcWeight);
   fChain->SetBranchAddress("nGenPar", &nGenPar, &b_nGenPar);
   fChain->SetBranchAddress("genParE", &genParE, &b_genParE);
   fChain->SetBranchAddress("genParPt", &genParPt, &b_genParPt);
   fChain->SetBranchAddress("genParEta", &genParEta, &b_genParEta);
   fChain->SetBranchAddress("genParPhi", &genParPhi, &b_genParPhi);
   fChain->SetBranchAddress("genParM", &genParM, &b_genParM);
   fChain->SetBranchAddress("genParQ", &genParQ, &b_genParQ);
   fChain->SetBranchAddress("genParId", &genParId, &b_genParId);
   fChain->SetBranchAddress("genParSt", &genParSt, &b_genParSt);
   fChain->SetBranchAddress("genMomParId", &genMomParId, &b_genMomParId);
   fChain->SetBranchAddress("genParIndex", &genParIndex, &b_genParIndex);
   fChain->SetBranchAddress("genNMo", &genNMo, &b_genNMo);
   fChain->SetBranchAddress("genNDa", &genNDa, &b_genNDa);
   fChain->SetBranchAddress("genMo1", &genMo1, &b_genMo1);
   fChain->SetBranchAddress("genMo2", &genMo2, &b_genMo2);
   fChain->SetBranchAddress("genDa1", &genDa1, &b_genDa1);
   fChain->SetBranchAddress("genDa2", &genDa2, &b_genDa2);
   fChain->SetBranchAddress("nGenJet", &nGenJet, &b_nGenJet);
   fChain->SetBranchAddress("genJetE", &genJetE, &b_genJetE);
   fChain->SetBranchAddress("genJetPt", &genJetPt, &b_genJetPt);
   fChain->SetBranchAddress("genJetEta", &genJetEta, &b_genJetEta);
   fChain->SetBranchAddress("genJetPhi", &genJetPhi, &b_genJetPhi);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("muPt", &muPt, &b_muPt);
   fChain->SetBranchAddress("muEta", &muEta, &b_muEta);
   fChain->SetBranchAddress("muPhi", &muPhi, &b_muPhi);
   fChain->SetBranchAddress("muM", &muM, &b_muM);
   fChain->SetBranchAddress("muTrkIso", &muTrkIso, &b_muTrkIso);
   fChain->SetBranchAddress("muCorrTrkIso", &muCorrTrkIso, &b_muCorrTrkIso);
   fChain->SetBranchAddress("muHcalIso", &muHcalIso, &b_muHcalIso);
   fChain->SetBranchAddress("muEcalIso", &muEcalIso, &b_muEcalIso);
   fChain->SetBranchAddress("muCharge", &muCharge, &b_muCharge);
   fChain->SetBranchAddress("muPassID", &muPassID, &b_muPassID);
   fChain->SetBranchAddress("muChHadIso", &muChHadIso, &b_muChHadIso);
   fChain->SetBranchAddress("muNeHadIso", &muNeHadIso, &b_muNeHadIso);
   fChain->SetBranchAddress("muGamIso", &muGamIso, &b_muGamIso);
   fChain->SetBranchAddress("muPUPt", &muPUPt, &b_muPUPt);
   fChain->SetBranchAddress("muCorrPfIso", &muCorrPfIso, &b_muCorrPfIso);
   fChain->SetBranchAddress("muTauCorrPfIso", &muTauCorrPfIso, &b_muTauCorrPfIso);
   fChain->SetBranchAddress("isGlobalMuon", &isGlobalMuon, &b_isGlobalMuon);
   fChain->SetBranchAddress("isTrackerMuon", &isTrackerMuon, &b_isTrackerMuon);
   fChain->SetBranchAddress("muPtErrx", &muPtErrx, &b_muPtErrx);
   fChain->SetBranchAddress("mudxy", &mudxy, &b_mudxy);
   fChain->SetBranchAddress("mudz", &mudz, &b_mudz);
   fChain->SetBranchAddress("muTrkLayers", &muTrkLayers, &b_muTrkLayers);
   fChain->SetBranchAddress("muPixelHits", &muPixelHits, &b_muPixelHits);
   fChain->SetBranchAddress("muHits", &muHits, &b_muHits);
   fChain->SetBranchAddress("muMatches", &muMatches, &b_muMatches);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("eleRho", &eleRho, &b_eleRho);
   fChain->SetBranchAddress("eleEt", &eleEt, &b_eleEt);
   fChain->SetBranchAddress("eleEnergy", &eleEnergy, &b_eleEnergy);
   fChain->SetBranchAddress("elePt", &elePt, &b_elePt);
   fChain->SetBranchAddress("eleEta", &eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", &elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleM", &eleM, &b_eleM);
   fChain->SetBranchAddress("eleScEta", &eleScEta, &b_eleScEta);
   fChain->SetBranchAddress("eleSigIhIh", &eleSigIhIh, &b_eleSigIhIh);
   fChain->SetBranchAddress("eleDelEtaIn", &eleDelEtaIn, &b_eleDelEtaIn);
   fChain->SetBranchAddress("eleDelPhiIn", &eleDelPhiIn, &b_eleDelPhiIn);
   fChain->SetBranchAddress("eleHoE", &eleHoE, &b_eleHoE);
   fChain->SetBranchAddress("eleTrkIso", &eleTrkIso, &b_eleTrkIso);
   fChain->SetBranchAddress("eleHcalIso", &eleHcalIso, &b_eleHcalIso);
   fChain->SetBranchAddress("eleEcalIso", &eleEcalIso, &b_eleEcalIso);
   fChain->SetBranchAddress("eleUserTrkIso", &eleUserTrkIso, &b_eleUserTrkIso);
   fChain->SetBranchAddress("eleUserCalIso", &eleUserCalIso, &b_eleUserCalIso);
   fChain->SetBranchAddress("eleCorrPfIso", &eleCorrPfIso, &b_eleCorrPfIso);
   fChain->SetBranchAddress("eleTauCorrPfIso", &eleTauCorrPfIso, &b_eleTauCorrPfIso);
   fChain->SetBranchAddress("eleEoverP", &eleEoverP, &b_eleEoverP);
   fChain->SetBranchAddress("eleDxy", &eleDxy, &b_eleDxy);
   fChain->SetBranchAddress("eleDz", &eleDz, &b_eleDz);
   fChain->SetBranchAddress("eleChHadIso", &eleChHadIso, &b_eleChHadIso);
   fChain->SetBranchAddress("eleNeHadIso", &eleNeHadIso, &b_eleNeHadIso);
   fChain->SetBranchAddress("eleGamIso", &eleGamIso, &b_eleGamIso);
   fChain->SetBranchAddress("eleMissingHits", &eleMissingHits, &b_eleMissingHits);
   fChain->SetBranchAddress("eleInBarrel", &eleInBarrel, &b_eleInBarrel);
   fChain->SetBranchAddress("eleInEndcap", &eleInEndcap, &b_eleInEndcap);
   fChain->SetBranchAddress("elePassConv", &elePassConv, &b_elePassConv);
   fChain->SetBranchAddress("elePassID", &elePassID, &b_elePassID);
   fChain->SetBranchAddress("pfMetCorrPt", &pfMetCorrPt, &b_pfMetCorrPt);
   fChain->SetBranchAddress("pfMetCorrPhi", &pfMetCorrPhi, &b_pfMetCorrPhi);
   fChain->SetBranchAddress("pfMetCorrSumEt", &pfMetCorrSumEt, &b_pfMetCorrSumEt);
   fChain->SetBranchAddress("pfMetCorrSig", &pfMetCorrSig, &b_pfMetCorrSig);
   fChain->SetBranchAddress("pfMetRawPt", &pfMetRawPt, &b_pfMetRawPt);
   fChain->SetBranchAddress("pfMetRawPhi", &pfMetRawPhi, &b_pfMetRawPhi);
   fChain->SetBranchAddress("pfMetRawSumEt", &pfMetRawSumEt, &b_pfMetRawSumEt);
   fChain->SetBranchAddress("pfMetRawCov00", &pfMetRawCov00, &b_pfMetRawCov00);
   fChain->SetBranchAddress("pfMetRawCov01", &pfMetRawCov01, &b_pfMetRawCov01);
   fChain->SetBranchAddress("pfMetRawCov10", &pfMetRawCov10, &b_pfMetRawCov10);
   fChain->SetBranchAddress("pfMetRawCov11", &pfMetRawCov11, &b_pfMetRawCov11);
   fChain->SetBranchAddress("CA8nJet", &CA8nJet, &b_CA8nJet);
   fChain->SetBranchAddress("CA8jetPt", &CA8jetPt, &b_CA8jetPt);
   fChain->SetBranchAddress("CA8jetEta", &CA8jetEta, &b_CA8jetEta);
   fChain->SetBranchAddress("CA8jetPhi", &CA8jetPhi, &b_CA8jetPhi);
   fChain->SetBranchAddress("CA8jetMass", &CA8jetMass, &b_CA8jetMass);
   fChain->SetBranchAddress("CA8jetEn", &CA8jetEn, &b_CA8jetEn);
   fChain->SetBranchAddress("CA8jetCorrUncUp", &CA8jetCorrUncUp, &b_CA8jetCorrUncUp);
   fChain->SetBranchAddress("CA8jetCorrUncDown", &CA8jetCorrUncDown, &b_CA8jetCorrUncDown);
   fChain->SetBranchAddress("CA8jetCharge", &CA8jetCharge, &b_CA8jetCharge);
   fChain->SetBranchAddress("CA8jetPartonFlavor", &CA8jetPartonFlavor, &b_CA8jetPartonFlavor);
   fChain->SetBranchAddress("CA8jetPassID", &CA8jetPassID, &b_CA8jetPassID);
   fChain->SetBranchAddress("CA8jetSSV", &CA8jetSSV, &b_CA8jetSSV);
   fChain->SetBranchAddress("CA8jetCSV", &CA8jetCSV, &b_CA8jetCSV);
   fChain->SetBranchAddress("CA8jetTCHP", &CA8jetTCHP, &b_CA8jetTCHP);
   fChain->SetBranchAddress("CA8jetTCHE", &CA8jetTCHE, &b_CA8jetTCHE);
   fChain->SetBranchAddress("CA8jetJP", &CA8jetJP, &b_CA8jetJP);
   fChain->SetBranchAddress("CA8jetJBP", &CA8jetJBP, &b_CA8jetJBP);
   fChain->SetBranchAddress("CA8jetTau1", &CA8jetTau1, &b_CA8jetTau1);
   fChain->SetBranchAddress("CA8jetTau2", &CA8jetTau2, &b_CA8jetTau2);
   fChain->SetBranchAddress("CA8jetTau3", &CA8jetTau3, &b_CA8jetTau3);
   fChain->SetBranchAddress("CA8jetTau4", &CA8jetTau4, &b_CA8jetTau4);
   fChain->SetBranchAddress("CA8jetMuEF", &CA8jetMuEF, &b_CA8jetMuEF);
   fChain->SetBranchAddress("CA8jetPhoEF", &CA8jetPhoEF, &b_CA8jetPhoEF);
   fChain->SetBranchAddress("CA8jetCEmEF", &CA8jetCEmEF, &b_CA8jetCEmEF);
   fChain->SetBranchAddress("CA8jetCHadEF", &CA8jetCHadEF, &b_CA8jetCHadEF);
   fChain->SetBranchAddress("CA8jetNEmEF", &CA8jetNEmEF, &b_CA8jetNEmEF);
   fChain->SetBranchAddress("CA8jetNHadEF", &CA8jetNHadEF, &b_CA8jetNHadEF);
   fChain->SetBranchAddress("CA8jetCMulti", &CA8jetCMulti, &b_CA8jetCMulti);
   fChain->SetBranchAddress("CA8jetPrunedPt", &CA8jetPrunedPt, &b_CA8jetPrunedPt);
   fChain->SetBranchAddress("CA8jetPrunedEta", &CA8jetPrunedEta, &b_CA8jetPrunedEta);
   fChain->SetBranchAddress("CA8jetPrunedPhi", &CA8jetPrunedPhi, &b_CA8jetPrunedPhi);
   fChain->SetBranchAddress("CA8jetPrunedMass", &CA8jetPrunedMass, &b_CA8jetPrunedMass);
   fChain->SetBranchAddress("CA8jetPrunedEn", &CA8jetPrunedEn, &b_CA8jetPrunedEn);
   fChain->SetBranchAddress("CA8jetPrunedCorrUncUp", &CA8jetPrunedCorrUncUp, &b_CA8jetPrunedCorrUncUp);
   fChain->SetBranchAddress("CA8jetPrunedCorrUncDown", &CA8jetPrunedCorrUncDown, &b_CA8jetPrunedCorrUncDown);
   fChain->SetBranchAddress("CA8jetPrunedCharge", &CA8jetPrunedCharge, &b_CA8jetPrunedCharge);
   fChain->SetBranchAddress("CA8jetPrunedPartonFlavor", &CA8jetPrunedPartonFlavor, &b_CA8jetPrunedPartonFlavor);
   fChain->SetBranchAddress("CA8jetPrunedSSV", &CA8jetPrunedSSV, &b_CA8jetPrunedSSV);
   fChain->SetBranchAddress("CA8jetPrunedCSV", &CA8jetPrunedCSV, &b_CA8jetPrunedCSV);
   fChain->SetBranchAddress("CA8jetPrunedTCHP", &CA8jetPrunedTCHP, &b_CA8jetPrunedTCHP);
   fChain->SetBranchAddress("CA8jetPrunedTCHE", &CA8jetPrunedTCHE, &b_CA8jetPrunedTCHE);
   fChain->SetBranchAddress("CA8jetPrunedJP", &CA8jetPrunedJP, &b_CA8jetPrunedJP);
   fChain->SetBranchAddress("CA8jetPrunedJBP", &CA8jetPrunedJBP, &b_CA8jetPrunedJBP);
   fChain->SetBranchAddress("CA8nSubPrunedJet", &CA8nSubPrunedJet, &b_CA8nSubPrunedJet);
   fChain->SetBranchAddress("CA8subjetMotherIndex", &CA8subjetMotherIndex, &b_CA8subjetMotherIndex);
   fChain->SetBranchAddress("CA8subjetPrunedPt", &CA8subjetPrunedPt, &b_CA8subjetPrunedPt);
   fChain->SetBranchAddress("CA8subjetPrunedEta", &CA8subjetPrunedEta, &b_CA8subjetPrunedEta);
   fChain->SetBranchAddress("CA8subjetPrunedPhi", &CA8subjetPrunedPhi, &b_CA8subjetPrunedPhi);
   fChain->SetBranchAddress("CA8subjetPrunedMass", &CA8subjetPrunedMass, &b_CA8subjetPrunedMass);
   fChain->SetBranchAddress("CA8subjetPrunedEn", &CA8subjetPrunedEn, &b_CA8subjetPrunedEn);
   fChain->SetBranchAddress("CA8subjetPrunedCharge", &CA8subjetPrunedCharge, &b_CA8subjetPrunedCharge);
   fChain->SetBranchAddress("CA8subjetPrunedPartonFlavor", &CA8subjetPrunedPartonFlavor, &b_CA8subjetPrunedPartonFlavor);
   fChain->SetBranchAddress("CA8subjetPrunedCSV", &CA8subjetPrunedCSV, &b_CA8subjetPrunedCSV);
   fChain->SetBranchAddress("AK5nJet", &AK5nJet, &b_AK5nJet);
   fChain->SetBranchAddress("AK5jetPt", &AK5jetPt, &b_AK5jetPt);
   fChain->SetBranchAddress("AK5jetEta", &AK5jetEta, &b_AK5jetEta);
   fChain->SetBranchAddress("AK5jetPhi", &AK5jetPhi, &b_AK5jetPhi);
   fChain->SetBranchAddress("AK5jetMass", &AK5jetMass, &b_AK5jetMass);
   fChain->SetBranchAddress("AK5jetEn", &AK5jetEn, &b_AK5jetEn);
   fChain->SetBranchAddress("AK5jetCorrUncUp", &AK5jetCorrUncUp, &b_AK5jetCorrUncUp);
   fChain->SetBranchAddress("AK5jetCorrUncDown", &AK5jetCorrUncDown, &b_AK5jetCorrUncDown);
   fChain->SetBranchAddress("AK5jetCharge", &AK5jetCharge, &b_AK5jetCharge);
   fChain->SetBranchAddress("AK5jetPartonFlavor", &AK5jetPartonFlavor, &b_AK5jetPartonFlavor);
   fChain->SetBranchAddress("AK5jetPassID", &AK5jetPassID, &b_AK5jetPassID);
   fChain->SetBranchAddress("AK5jetSSV", &AK5jetSSV, &b_AK5jetSSV);
   fChain->SetBranchAddress("AK5jetCSV", &AK5jetCSV, &b_AK5jetCSV);
   fChain->SetBranchAddress("AK5jetTCHP", &AK5jetTCHP, &b_AK5jetTCHP);
   fChain->SetBranchAddress("AK5jetTCHE", &AK5jetTCHE, &b_AK5jetTCHE);
   fChain->SetBranchAddress("AK5jetJP", &AK5jetJP, &b_AK5jetJP);
   fChain->SetBranchAddress("AK5jetJBP", &AK5jetJBP, &b_AK5jetJBP);
   fChain->SetBranchAddress("AK5jetTau1", &AK5jetTau1, &b_AK5jetTau1);
   fChain->SetBranchAddress("AK5jetTau2", &AK5jetTau2, &b_AK5jetTau2);
   fChain->SetBranchAddress("AK5jetTau3", &AK5jetTau3, &b_AK5jetTau3);
   fChain->SetBranchAddress("AK5jetTau4", &AK5jetTau4, &b_AK5jetTau4);
   fChain->SetBranchAddress("AK5jetMuEF", &AK5jetMuEF, &b_AK5jetMuEF);
   fChain->SetBranchAddress("AK5jetPhoEF", &AK5jetPhoEF, &b_AK5jetPhoEF);
   fChain->SetBranchAddress("AK5jetCEmEF", &AK5jetCEmEF, &b_AK5jetCEmEF);
   fChain->SetBranchAddress("AK5jetCHadEF", &AK5jetCHadEF, &b_AK5jetCHadEF);
   fChain->SetBranchAddress("AK5jetNEmEF", &AK5jetNEmEF, &b_AK5jetNEmEF);
   fChain->SetBranchAddress("AK5jetNHadEF", &AK5jetNHadEF, &b_AK5jetNHadEF);
   fChain->SetBranchAddress("AK5jetCMulti", &AK5jetCMulti, &b_AK5jetCMulti);
   fChain->SetBranchAddress("hlt_nTrigs", &hlt_nTrigs, &b_hlt_nTrigs);
   fChain->SetBranchAddress("hlt_trigResult", &hlt_trigResult, &b_hlt_trigResult);
   fChain->SetBranchAddress("hlt_trigName", &hlt_trigName, &b_hlt_trigName);
   Notify();
}

Bool_t skimTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void skimTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
#endif // #ifdef skimTree_cxx
