//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun 28 20:35:17 2012 by ROOT version 5.27/06b
// from TTree tree/tree
// found on file: ../debug/hzz.root
//////////////////////////////////////////////////////////

#ifndef mjj_why_h
#define mjj_why_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxEvtInfo_VertexX = 1;
const Int_t kMaxEvtInfo_VertexY = 1;
const Int_t kMaxEvtInfo_VertexZ = 1;
const Int_t kMaxptHat = 1;
const Int_t kMaxmcWeight = 1;
const Int_t kMaxgenParE = 1;
const Int_t kMaxgenParPt = 1;
const Int_t kMaxgenParEta = 1;
const Int_t kMaxgenParPhi = 1;
const Int_t kMaxgenParM = 1;
const Int_t kMaxgenParQ = 1;
const Int_t kMaxgenParId = 1;
const Int_t kMaxgenParSt = 1;
const Int_t kMaxgenMomParId = 1;
const Int_t kMaxgenParIndex = 1;
const Int_t kMaxgenJetE = 1;
const Int_t kMaxgenJetPt = 1;
const Int_t kMaxgenJetEta = 1;
const Int_t kMaxgenJetPhi = 1;



class mjj_why {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
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
   vector<double>  *genParM_;
   vector<int>     *genParQ_;
   vector<int>     *genParId_;
   vector<int>     *genParSt_;
   vector<int>     *genMomParId_;
   vector<int>     *genParIndex_;
   vector<double>  *genJetE_;
   vector<double>  *genJetPt_;
   vector<double>  *genJetEta_;
   vector<double>  *genJetPhi_;
   Int_t           nAllHCand;
   Double_t        metSig;
   Double_t        metSigNoPU;
   vector<double>  *higgsPt;
   vector<double>  *higgsEta;
   vector<double>  *higgsPhi;
   vector<double>  *higgsM;
   vector<double>  *higgsMRefit;
   vector<double>  *zllPt;
   vector<double>  *zllEta;
   vector<double>  *zllPhi;
   vector<double>  *zllM;
   vector<double>  *zlldR;
   vector<double>  *zjjPt;
   vector<double>  *zjjEta;
   vector<double>  *zjjPhi;
   vector<double>  *zjjM;
   vector<double>  *zjjMRefit;
   vector<double>  *zjjdR;
   vector<int>     *jetIndex;
   vector<int>     *jetHiggsIndex;
   vector<double>  *jetE;
   vector<double>  *jetPt;
   vector<double>  *jetEta;
   vector<double>  *jetPhi;
   vector<double>  *jetGenE;
   vector<double>  *jetGenPt;
   vector<double>  *jetGenEta;
   vector<double>  *jetGenPhi;
   vector<double>  *jetPartonE;
   vector<double>  *jetPartonPt;
   vector<double>  *jetPartonEta;
   vector<double>  *jetPartonPhi;
   vector<double>  *heliLD;
   vector<double>  *heliLDRefit;
   vector<int>     *nBTags;
   vector<int>     *lepType;
   vector<int>     *passBit;

   // List of branches
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
   TBranch        *b_genParM_;   //!
   TBranch        *b_genParQ_;   //!
   TBranch        *b_genParId_;   //!
   TBranch        *b_genParSt_;   //!
   TBranch        *b_genMomParId_;   //!
   TBranch        *b_genParIndex_;   //!
   TBranch        *b_genJetE_;   //!
   TBranch        *b_genJetPt_;   //!
   TBranch        *b_genJetEta_;   //!
   TBranch        *b_genJetPhi_;   //!
   TBranch        *b_nAllHCand;   //!
   TBranch        *b_metSig;   //!
   TBranch        *b_metSigNoPU;   //!
   TBranch        *b_higgsPt;   //!
   TBranch        *b_higgsEta;   //!
   TBranch        *b_higgsPhi;   //!
   TBranch        *b_higgsM;   //!
   TBranch        *b_higgsMRefit;   //!
   TBranch        *b_zllPt;   //!
   TBranch        *b_zllEta;   //!
   TBranch        *b_zllPhi;   //!
   TBranch        *b_zllM;   //!
   TBranch        *b_zlldR;   //!
   TBranch        *b_zjjPt;   //!
   TBranch        *b_zjjEta;   //!
   TBranch        *b_zjjPhi;   //!
   TBranch        *b_zjjM;   //!
   TBranch        *b_zjjMRefit;   //!
   TBranch        *b_zjjdR;   //!
   TBranch        *b_jetIndex;   //!
   TBranch        *b_jetHiggsIndex;   //!
   TBranch        *b_jetE;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetGenE;   //!
   TBranch        *b_jetGenPt;   //!
   TBranch        *b_jetGenEta;   //!
   TBranch        *b_jetGenPhi;   //!
   TBranch        *b_jetPartonE;   //!
   TBranch        *b_jetPartonPt;   //!
   TBranch        *b_jetPartonEta;   //!
   TBranch        *b_jetPartonPhi;   //!
   TBranch        *b_heliLD;   //!
   TBranch        *b_heliLDRefit;   //!
   TBranch        *b_nBTags;   //!
   TBranch        *b_lepType;   //!
   TBranch        *b_passBit;   //!

   mjj_why(string filename, TTree *tree=0);
   virtual ~mjj_why();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int DEBUG=0);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Bool_t   matchGenToParton(Int_t igen, Int_t ijet);
   virtual Bool_t   matchGenToPFJet(Int_t igen, Int_t ijet);

   string _inputFileName;
};

#endif

#ifdef mjj_why_cxx
mjj_why::mjj_why(string filename, TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename.data());
      if (!f) {
	f = new TFile(filename.data());
	f->cd(Form("%s:/tree",filename.data()));
      }
      tree = (TTree*)gDirectory->Get("tree");

   }

   _inputFileName = filename;
   Init(tree);
}

mjj_why::~mjj_why()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t mjj_why::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t mjj_why::LoadTree(Long64_t entry)
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

void mjj_why::Init(TTree *tree)
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
   genParM_ = 0;
   genParQ_ = 0;
   genParId_ = 0;
   genParSt_ = 0;
   genMomParId_ = 0;
   genParIndex_ = 0;
   genJetE_ = 0;
   genJetPt_ = 0;
   genJetEta_ = 0;
   genJetPhi_ = 0;
   higgsPt = 0;
   higgsEta = 0;
   higgsPhi = 0;
   higgsM = 0;
   higgsMRefit = 0;
   zllPt = 0;
   zllEta = 0;
   zllPhi = 0;
   zllM = 0;
   zlldR = 0;
   zjjPt = 0;
   zjjEta = 0;
   zjjPhi = 0;
   zjjM = 0;
   zjjMRefit = 0;
   zjjdR = 0;
   jetIndex = 0;
   jetHiggsIndex = 0;
   jetE = 0;
   jetPt = 0;
   jetEta = 0;
   jetPhi = 0;
   jetGenE = 0;
   jetGenPt = 0;
   jetGenEta = 0;
   jetGenPhi = 0;
   jetPartonE = 0;
   jetPartonPt = 0;
   jetPartonEta = 0;
   jetPartonPhi = 0;
   heliLD = 0;
   heliLDRefit = 0;
   nBTags = 0;
   lepType = 0;
   passBit = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

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
   fChain->SetBranchAddress("genParM_", &genParM_, &b_genParM_);
   fChain->SetBranchAddress("genParQ_", &genParQ_, &b_genParQ_);
   fChain->SetBranchAddress("genParId_", &genParId_, &b_genParId_);
   fChain->SetBranchAddress("genParSt_", &genParSt_, &b_genParSt_);
   fChain->SetBranchAddress("genMomParId_", &genMomParId_, &b_genMomParId_);
   fChain->SetBranchAddress("genParIndex_", &genParIndex_, &b_genParIndex_);
   fChain->SetBranchAddress("genJetE_", &genJetE_, &b_genJetE_);
   fChain->SetBranchAddress("genJetPt_", &genJetPt_, &b_genJetPt_);
   fChain->SetBranchAddress("genJetEta_", &genJetEta_, &b_genJetEta_);
   fChain->SetBranchAddress("genJetPhi_", &genJetPhi_, &b_genJetPhi_);
   fChain->SetBranchAddress("nAllHCand", &nAllHCand, &b_nAllHCand);
   fChain->SetBranchAddress("metSig", &metSig, &b_metSig);
   fChain->SetBranchAddress("metSigNoPU", &metSigNoPU, &b_metSigNoPU);
   fChain->SetBranchAddress("higgsPt", &higgsPt, &b_higgsPt);
   fChain->SetBranchAddress("higgsEta", &higgsEta, &b_higgsEta);
   fChain->SetBranchAddress("higgsPhi", &higgsPhi, &b_higgsPhi);
   fChain->SetBranchAddress("higgsM", &higgsM, &b_higgsM);
   fChain->SetBranchAddress("higgsMRefit", &higgsMRefit, &b_higgsMRefit);
   fChain->SetBranchAddress("zllPt", &zllPt, &b_zllPt);
   fChain->SetBranchAddress("zllEta", &zllEta, &b_zllEta);
   fChain->SetBranchAddress("zllPhi", &zllPhi, &b_zllPhi);
   fChain->SetBranchAddress("zllM", &zllM, &b_zllM);
   fChain->SetBranchAddress("zlldR", &zlldR, &b_zlldR);
   fChain->SetBranchAddress("zjjPt", &zjjPt, &b_zjjPt);
   fChain->SetBranchAddress("zjjEta", &zjjEta, &b_zjjEta);
   fChain->SetBranchAddress("zjjPhi", &zjjPhi, &b_zjjPhi);
   fChain->SetBranchAddress("zjjM", &zjjM, &b_zjjM);
   fChain->SetBranchAddress("zjjMRefit", &zjjMRefit, &b_zjjMRefit);
   fChain->SetBranchAddress("zjjdR", &zjjdR, &b_zjjdR);
   fChain->SetBranchAddress("jetIndex", &jetIndex, &b_jetIndex);
   fChain->SetBranchAddress("jetHiggsIndex", &jetHiggsIndex, &b_jetHiggsIndex);
   fChain->SetBranchAddress("jetE", &jetE, &b_jetE);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetGenE", &jetGenE, &b_jetGenE);
   fChain->SetBranchAddress("jetGenPt", &jetGenPt, &b_jetGenPt);
   fChain->SetBranchAddress("jetGenEta", &jetGenEta, &b_jetGenEta);
   fChain->SetBranchAddress("jetGenPhi", &jetGenPhi, &b_jetGenPhi);
   fChain->SetBranchAddress("jetPartonE", &jetPartonE, &b_jetPartonE);
   fChain->SetBranchAddress("jetPartonPt", &jetPartonPt, &b_jetPartonPt);
   fChain->SetBranchAddress("jetPartonEta", &jetPartonEta, &b_jetPartonEta);
   fChain->SetBranchAddress("jetPartonPhi", &jetPartonPhi, &b_jetPartonPhi);
   fChain->SetBranchAddress("heliLD", &heliLD, &b_heliLD);
   fChain->SetBranchAddress("heliLDRefit", &heliLDRefit, &b_heliLDRefit);
   fChain->SetBranchAddress("nBTags", &nBTags, &b_nBTags);
   fChain->SetBranchAddress("lepType", &lepType, &b_lepType);
   fChain->SetBranchAddress("passBit", &passBit, &b_passBit);
   Notify();
}

Bool_t mjj_why::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void mjj_why::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t mjj_why::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef mjj_why_cxx
