//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan  3 18:12:00 2012 by ROOT version 5.27/06b
// from TTree EventTree/Event data
// found on file: /home/cdxfe2/Summer2011Ntuples/YJPhotonjetanaQCD15.root
//////////////////////////////////////////////////////////

#ifndef yjqcd_h
#define yjqcd_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string>
#include <iostream>

const int MAXN=50;


using namespace std;

class yjqcd {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           orbit;
   Int_t           bx;
   Int_t           lumis;
   Bool_t          isData;
   Float_t         bspotPos[3];
   Int_t           nVtxGood;
   Int_t           nVtxNotFake;
   Float_t         vertexX[MAXN];   //[nVtxNotFake]
   Float_t         vertexY[MAXN];   //[nVtxNotFake]
   Float_t         vertexZ[MAXN];   //[nVtxNotFake]
   Float_t         vertexXError[MAXN];   //[nVtxNotFake]
   Float_t         vertexYError[MAXN];   //[nVtxNotFake]
   Float_t         vertexZError[MAXN];   //[nVtxNotFake]
   Float_t         vertexChi2[MAXN];   //[nVtxNotFake]
   Float_t         vertexNormChi2[MAXN];   //[nVtxNotFake]
   Float_t         vertexNdof[MAXN];   //[nVtxNotFake]
   Float_t         vertexNTrk[MAXN];   //[nVtxNotFake]
   Float_t         vertexNTrkWeight05[MAXN];   //[nVtxNotFake]
   Int_t           nPatJet;
   Float_t         PatJetEn[MAXN];   //[nPatJet]
   Float_t         PatJetPt[MAXN];   //[nPatJet]
   Float_t         PatJetEta[MAXN];   //[nPatJet]
   Float_t         PatJetPhi[MAXN];   //[nPatJet]
   Float_t         PatJetEt[MAXN];   //[nPatJet]
   Float_t         PatJetRawPt[MAXN];   //[nPatJet]
   Float_t         PatJetRawEn[MAXN];   //[nPatJet]
   Float_t         PatJetCHF[MAXN];   //[nPatJet]
   Float_t         PatJetNHF[MAXN];   //[nPatJet]
   Float_t         PatJetCEF[MAXN];   //[nPatJet]
   Float_t         PatJetNEF[MAXN];   //[nPatJet]
   Int_t           PatJetNCH[MAXN];   //[nPatJet]
   Int_t           PatJetPartonID[MAXN];   //[nPatJet]
   Int_t           PatJetGenJetIndex[MAXN];   //[nPatJet]
   Float_t         PatJetGenJetEn[MAXN];   //[nPatJet]
   Float_t         PatJetGenJetPt[MAXN];   //[nPatJet]
   Float_t         PatJetGenJetEta[MAXN];   //[nPatJet]
   Float_t         PatJetGenJetPhi[MAXN];   //[nPatJet]
   Int_t           PatJetGenPartonID[MAXN];   //[nPatJet]
   Float_t         PatJetGenEn[MAXN];   //[nPatJet]
   Float_t         PatJetGenPt[MAXN];   //[nPatJet]
   Float_t         PatJetGenEta[MAXN];   //[nPatJet]
   Float_t         PatJetGenPhi[MAXN];   //[nPatJet]
   Int_t           nPho;
   Float_t         PhoP[MAXN];   //[nPho]
   Float_t         PhoEt[MAXN];   //[nPho]
   Float_t         PhoEnergy[MAXN];   //[nPho]
   Float_t         PhoPx[MAXN];   //[nPho]
   Float_t         PhoPy[MAXN];   //[nPho]
   Float_t         PhoPz[MAXN];   //[nPho]
   Float_t         PhoPt[MAXN];   //[nPho]
   Float_t         PhoEta[MAXN];   //[nPho]
   Float_t         PhoPhi[MAXN];   //[nPho]
   Float_t         PhoR9[MAXN];   //[nPho]
   Float_t         PhoPhiWidth[MAXN];   //[nPho]
   Float_t         PhoEtaWidth[MAXN];   //[nPho]
   Float_t         PhoScPhi[MAXN];   //[nPho]
   Float_t         PhoScEta[MAXN];   //[nPho]
   Float_t         PhoESRatio[MAXN];   //[nPho]
   Float_t         PhoSigmaIetaIeta[MAXN];   //[nPho]
   Float_t         PhoSigmaIphiIphi[MAXN];   //[nPho]
   Float_t         PhoSeedTime[MAXN];   //[nPho]
   Int_t           PhoseedSeverity[MAXN];   //[nPho]
   Float_t         PhoE2overe9[MAXN];   //[nPho]
   Float_t         PhohadronicOverEm[MAXN];   //[nPho]
   Float_t         PhoecalRecHitSumEtConeDR04[MAXN];   //[nPho]
   Float_t         PhohcalTowerSumEtConeDR04[MAXN];   //[nPho]
   Float_t         PhotrkSumPtHollowConeDR04[MAXN];   //[nPho]
   Float_t         PhoisConverted[MAXN];   //[nPho]
   Float_t         PhohasPixelSeed[MAXN];   //[nPho]
   Float_t         PhoisGenMatched[MAXN];   //[nPho]
   Float_t         PhogenMomId[MAXN];   //[nPho]
   Float_t         PhogenGrandMomId[MAXN];   //[nPho]
   Float_t         PhogenNSiblings[MAXN];   //[nPho]
   Float_t         PhogenMatchedE[MAXN];   //[nPho]
   Float_t         PhogenMatchedPx[MAXN];   //[nPho]
   Float_t         PhogenMatchedPy[MAXN];   //[nPho]
   Float_t         PhogenMatchedPz[MAXN];   //[nPho]
   Float_t         PhogenMatchedPt[MAXN];   //[nPho]
   Float_t         PhogenMatchedEta[MAXN];   //[nPho]
   Float_t         PhogenMatchedPhi[MAXN];   //[nPho]
   Float_t         PhogenCalIsoDR04[MAXN];   //[nPho]
   Float_t         PhogenTrkIsoDR04[MAXN];   //[nPho]
   Float_t         PhogenIsoDR04[MAXN];   //[nPho]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_orbit;   //!
   TBranch        *b_bx;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_bspotPos;   //!
   TBranch        *b_nVtxGood;   //!
   TBranch        *b_nVtxNotFake;   //!
   TBranch        *b_vertexX;   //!
   TBranch        *b_vertexY;   //!
   TBranch        *b_vertexZ;   //!
   TBranch        *b_vertexXError;   //!
   TBranch        *b_vertexYError;   //!
   TBranch        *b_vertexZError;   //!
   TBranch        *b_vertexChi2;   //!
   TBranch        *b_vertexNormChi2;   //!
   TBranch        *b_vertexNdof;   //!
   TBranch        *b_vertexNTrk;   //!
   TBranch        *b_vertexNTrkWeight05;   //!
   TBranch        *b_nPatJet;   //!
   TBranch        *b_PatJetEn;   //!
   TBranch        *b_PatJetPt;   //!
   TBranch        *b_PatJetEta;   //!
   TBranch        *b_PatJetPhi;   //!
   TBranch        *b_PatJetEt;   //!
   TBranch        *b_PatJetRawPt;   //!
   TBranch        *b_PatJetRawEn;   //!
   TBranch        *b_PatJetCHF;   //!
   TBranch        *b_PatJetNHF;   //!
   TBranch        *b_PatJetCEF;   //!
   TBranch        *b_PatJetNEF;   //!
   TBranch        *b_PatJetNCH;   //!
   TBranch        *b_PatJetPartonID;   //!
   TBranch        *b_PatJetGenJetIndex;   //!
   TBranch        *b_PatJetGenJetEn;   //!
   TBranch        *b_PatJetGenJetPt;   //!
   TBranch        *b_PatJetGenJetEta;   //!
   TBranch        *b_PatJetGenJetPhi;   //!
   TBranch        *b_PatJetGenPartonID;   //!
   TBranch        *b_PatJetGenEn;   //!
   TBranch        *b_PatJetGenPt;   //!
   TBranch        *b_PatJetGenEta;   //!
   TBranch        *b_PatJetGenPhi;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_PhoP;   //!
   TBranch        *b_PhoEt;   //!
   TBranch        *b_PhoEnergy;   //!
   TBranch        *b_PhoPx;   //!
   TBranch        *b_PhoPy;   //!
   TBranch        *b_PhoPz;   //!
   TBranch        *b_PhoPt;   //!
   TBranch        *b_PhoEta;   //!
   TBranch        *b_PhoPhi;   //!
   TBranch        *b_PhoR9;   //!
   TBranch        *b_PhoPhiWidth;   //!
   TBranch        *b_PhoEtaWidth;   //!
   TBranch        *b_PhoScPhi;   //!
   TBranch        *b_PhoScEta;   //!
   TBranch        *b_PhoESRatio;   //!
   TBranch        *b_PhoSigmaIetaIeta;   //!
   TBranch        *b_PhoSigmaIphiIphi;   //!
   TBranch        *b_PhoSeedTime;   //!
   TBranch        *b_PhoseedSeverity;   //!
   TBranch        *b_PhoE2overe9;   //!
   TBranch        *b_PhohadronicOverEm;   //!
   TBranch        *b_PhoecalRecHitSumEtConeDR04;   //!
   TBranch        *b_PhohcalTowerSumEtConeDR04;   //!
   TBranch        *b_PhotrkSumPtHollowConeDR04;   //!
   TBranch        *b_PhoisConverted;   //!
   TBranch        *b_PhohasPixelSeed;   //!
   TBranch        *b_PhoisGenMatched;   //!
   TBranch        *b_PhogenMomId;   //!
   TBranch        *b_PhogenGrandMomId;   //!
   TBranch        *b_PhogenNSiblings;   //!
   TBranch        *b_PhogenMatchedE;   //!
   TBranch        *b_PhogenMatchedPx;   //!
   TBranch        *b_PhogenMatchedPy;   //!
   TBranch        *b_PhogenMatchedPz;   //!
   TBranch        *b_PhogenMatchedPt;   //!
   TBranch        *b_PhogenMatchedEta;   //!
   TBranch        *b_PhogenMatchedPhi;   //!
   TBranch        *b_PhogenCalIsoDR04;   //!
   TBranch        *b_PhogenTrkIsoDR04;   //!
   TBranch        *b_PhogenIsoDR04;   //!

   yjqcd(std::string inputFileName);
   virtual ~yjqcd();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   
   std::string      _inputFile;

};

#endif

#ifdef yjqcd_cxx
yjqcd::yjqcd(std::string inputFileName)
{
  _inputFile= inputFileName;
  
  cout << "Reading in file " << _inputFile << endl;
  
  TTree* tree=0;
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(_inputFile.data());
      if (!f) {
	f = new TFile(_inputFile.data());
	f->cd(Form("%s:/photonjetana",_inputFile.data()));
      }
      tree = (TTree*)gDirectory->Get("EventTree");

   }
   Init(tree);
}

yjqcd::~yjqcd()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t yjqcd::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t yjqcd::LoadTree(Long64_t entry)
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

void yjqcd::Init(TTree *tree)
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
   fChain->SetBranchAddress("bx", &bx, &b_bx);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("bspotPos", bspotPos, &b_bspotPos);
   fChain->SetBranchAddress("nVtxGood", &nVtxGood, &b_nVtxGood);
   fChain->SetBranchAddress("nVtxNotFake", &nVtxNotFake, &b_nVtxNotFake);
   fChain->SetBranchAddress("vertexX", vertexX, &b_vertexX);
   fChain->SetBranchAddress("vertexY", vertexY, &b_vertexY);
   fChain->SetBranchAddress("vertexZ", vertexZ, &b_vertexZ);
   fChain->SetBranchAddress("vertexXError", vertexXError, &b_vertexXError);
   fChain->SetBranchAddress("vertexYError", vertexYError, &b_vertexYError);
   fChain->SetBranchAddress("vertexZError", vertexZError, &b_vertexZError);
   fChain->SetBranchAddress("vertexChi2", vertexChi2, &b_vertexChi2);
   fChain->SetBranchAddress("vertexNormChi2", vertexNormChi2, &b_vertexNormChi2);
   fChain->SetBranchAddress("vertexNdof", vertexNdof, &b_vertexNdof);
   fChain->SetBranchAddress("vertexNTrk", vertexNTrk, &b_vertexNTrk);
   fChain->SetBranchAddress("vertexNTrkWeight05", vertexNTrkWeight05, &b_vertexNTrkWeight05);
   fChain->SetBranchAddress("nPatJet", &nPatJet, &b_nPatJet);
   fChain->SetBranchAddress("PatJetEn", PatJetEn, &b_PatJetEn);
   fChain->SetBranchAddress("PatJetPt", PatJetPt, &b_PatJetPt);
   fChain->SetBranchAddress("PatJetEta", PatJetEta, &b_PatJetEta);
   fChain->SetBranchAddress("PatJetPhi", PatJetPhi, &b_PatJetPhi);
   fChain->SetBranchAddress("PatJetEt", PatJetEt, &b_PatJetEt);
   fChain->SetBranchAddress("PatJetRawPt", PatJetRawPt, &b_PatJetRawPt);
   fChain->SetBranchAddress("PatJetRawEn", PatJetRawEn, &b_PatJetRawEn);
   fChain->SetBranchAddress("PatJetCHF", PatJetCHF, &b_PatJetCHF);
   fChain->SetBranchAddress("PatJetNHF", PatJetNHF, &b_PatJetNHF);
   fChain->SetBranchAddress("PatJetCEF", PatJetCEF, &b_PatJetCEF);
   fChain->SetBranchAddress("PatJetNEF", PatJetNEF, &b_PatJetNEF);
   fChain->SetBranchAddress("PatJetNCH", PatJetNCH, &b_PatJetNCH);
   fChain->SetBranchAddress("PatJetPartonID", PatJetPartonID, &b_PatJetPartonID);
   fChain->SetBranchAddress("PatJetGenJetIndex", PatJetGenJetIndex, &b_PatJetGenJetIndex);
   fChain->SetBranchAddress("PatJetGenJetEn", PatJetGenJetEn, &b_PatJetGenJetEn);
   fChain->SetBranchAddress("PatJetGenJetPt", PatJetGenJetPt, &b_PatJetGenJetPt);
   fChain->SetBranchAddress("PatJetGenJetEta", PatJetGenJetEta, &b_PatJetGenJetEta);
   fChain->SetBranchAddress("PatJetGenJetPhi", PatJetGenJetPhi, &b_PatJetGenJetPhi);
   fChain->SetBranchAddress("PatJetGenPartonID", PatJetGenPartonID, &b_PatJetGenPartonID);
   fChain->SetBranchAddress("PatJetGenEn", PatJetGenEn, &b_PatJetGenEn);
   fChain->SetBranchAddress("PatJetGenPt", PatJetGenPt, &b_PatJetGenPt);
   fChain->SetBranchAddress("PatJetGenEta", PatJetGenEta, &b_PatJetGenEta);
   fChain->SetBranchAddress("PatJetGenPhi", PatJetGenPhi, &b_PatJetGenPhi);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("PhoP", PhoP, &b_PhoP);
   fChain->SetBranchAddress("PhoEt", PhoEt, &b_PhoEt);
   fChain->SetBranchAddress("PhoEnergy", PhoEnergy, &b_PhoEnergy);
   fChain->SetBranchAddress("PhoPx", PhoPx, &b_PhoPx);
   fChain->SetBranchAddress("PhoPy", PhoPy, &b_PhoPy);
   fChain->SetBranchAddress("PhoPz", PhoPz, &b_PhoPz);
   fChain->SetBranchAddress("PhoPt", PhoPt, &b_PhoPt);
   fChain->SetBranchAddress("PhoEta", PhoEta, &b_PhoEta);
   fChain->SetBranchAddress("PhoPhi", PhoPhi, &b_PhoPhi);
   fChain->SetBranchAddress("PhoR9", PhoR9, &b_PhoR9);
   fChain->SetBranchAddress("PhoPhiWidth", PhoPhiWidth, &b_PhoPhiWidth);
   fChain->SetBranchAddress("PhoEtaWidth", PhoEtaWidth, &b_PhoEtaWidth);
   fChain->SetBranchAddress("PhoScPhi", PhoScPhi, &b_PhoScPhi);
   fChain->SetBranchAddress("PhoScEta", PhoScEta, &b_PhoScEta);
   fChain->SetBranchAddress("PhoESRatio", PhoESRatio, &b_PhoESRatio);
   fChain->SetBranchAddress("PhoSigmaIetaIeta", PhoSigmaIetaIeta, &b_PhoSigmaIetaIeta);
   fChain->SetBranchAddress("PhoSigmaIphiIphi", PhoSigmaIphiIphi, &b_PhoSigmaIphiIphi);
   fChain->SetBranchAddress("PhoSeedTime", PhoSeedTime, &b_PhoSeedTime);
   fChain->SetBranchAddress("PhoseedSeverity", PhoseedSeverity, &b_PhoseedSeverity);
   fChain->SetBranchAddress("PhoE2overe9", PhoE2overe9, &b_PhoE2overe9);
   fChain->SetBranchAddress("PhohadronicOverEm", PhohadronicOverEm, &b_PhohadronicOverEm);
   fChain->SetBranchAddress("PhoecalRecHitSumEtConeDR04", PhoecalRecHitSumEtConeDR04, &b_PhoecalRecHitSumEtConeDR04);
   fChain->SetBranchAddress("PhohcalTowerSumEtConeDR04", PhohcalTowerSumEtConeDR04, &b_PhohcalTowerSumEtConeDR04);
   fChain->SetBranchAddress("PhotrkSumPtHollowConeDR04", PhotrkSumPtHollowConeDR04, &b_PhotrkSumPtHollowConeDR04);
   fChain->SetBranchAddress("PhoisConverted", PhoisConverted, &b_PhoisConverted);
   fChain->SetBranchAddress("PhohasPixelSeed", PhohasPixelSeed, &b_PhohasPixelSeed);
   fChain->SetBranchAddress("PhoisGenMatched", PhoisGenMatched, &b_PhoisGenMatched);
   fChain->SetBranchAddress("PhogenMomId", PhogenMomId, &b_PhogenMomId);
   fChain->SetBranchAddress("PhogenGrandMomId", PhogenGrandMomId, &b_PhogenGrandMomId);
   fChain->SetBranchAddress("PhogenNSiblings", PhogenNSiblings, &b_PhogenNSiblings);
   fChain->SetBranchAddress("PhogenMatchedE", PhogenMatchedE, &b_PhogenMatchedE);
   fChain->SetBranchAddress("PhogenMatchedPx", PhogenMatchedPx, &b_PhogenMatchedPx);
   fChain->SetBranchAddress("PhogenMatchedPy", PhogenMatchedPy, &b_PhogenMatchedPy);
   fChain->SetBranchAddress("PhogenMatchedPz", PhogenMatchedPz, &b_PhogenMatchedPz);
   fChain->SetBranchAddress("PhogenMatchedPt", PhogenMatchedPt, &b_PhogenMatchedPt);
   fChain->SetBranchAddress("PhogenMatchedEta", PhogenMatchedEta, &b_PhogenMatchedEta);
   fChain->SetBranchAddress("PhogenMatchedPhi", PhogenMatchedPhi, &b_PhogenMatchedPhi);
   fChain->SetBranchAddress("PhogenCalIsoDR04", PhogenCalIsoDR04, &b_PhogenCalIsoDR04);
   fChain->SetBranchAddress("PhogenTrkIsoDR04", PhogenTrkIsoDR04, &b_PhogenTrkIsoDR04);
   fChain->SetBranchAddress("PhogenIsoDR04", PhogenIsoDR04, &b_PhogenIsoDR04);
   Notify();
}

Bool_t yjqcd::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void yjqcd::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t yjqcd::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef yjqcd_cxx
