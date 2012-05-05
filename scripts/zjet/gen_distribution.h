//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr  3 17:56:46 2012 by ROOT version 5.27/06b
// from TTree tree/tree
// found on file: /scratch/syu/test_sherpa/genonly/zjets_mc_1_1_WVP.root
//////////////////////////////////////////////////////////

#ifndef gen_distribution_h
#define gen_distribution_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream>
#include <string>

using namespace std;

   const Int_t kMaxEvtInfo_VertexX = 1;
   const Int_t kMaxEvtInfo_VertexY = 1;
   const Int_t kMaxEvtInfo_VertexZ = 1;
   const Int_t kMaxptHat = 1;
   const Int_t kMaxmcWeight = 1;
   const Int_t kMaxgenParPx = 1;
   const Int_t kMaxgenParPy = 1;
   const Int_t kMaxgenParPz = 1;
   const Int_t kMaxgenParE = 1;
   const Int_t kMaxgenParP = 1;
   const Int_t kMaxgenParTheta = 1;
   const Int_t kMaxgenParPt = 1;
   const Int_t kMaxgenParEta = 1;
   const Int_t kMaxgenParPhi = 1;
   const Int_t kMaxgenParEt = 1;
   const Int_t kMaxgenParQ = 1;
   const Int_t kMaxgenParId = 1;
   const Int_t kMaxgenParSt = 1;
   const Int_t kMaxgenJetPx = 1;
   const Int_t kMaxgenJetPy = 1;
   const Int_t kMaxgenJetPz = 1;
   const Int_t kMaxgenJetE = 1;
   const Int_t kMaxgenJetP = 1;
   const Int_t kMaxgenJetTheta = 1;
   const Int_t kMaxgenJetPt = 1;
   const Int_t kMaxgenJetEta = 1;
   const Int_t kMaxgenJetPhi = 1;
   const Int_t kMaxgenJetEt = 1;
   const Int_t kMaxgenJetQ = 1;

class gen_distribution {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        PU_weight;
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
   //   vector<double>  *genParPx_;
   //   vector<double>  *genParPy_;
   //   vector<double>  *genParPz_;
   vector<double>  *genParE_;
   //   vector<double>  *genParP_;
   //   vector<double>  *genParTheta_;
   vector<double>  *genParPt_;
   vector<double>  *genParEta_;
   vector<double>  *genParPhi_;
   //   vector<double>  *genParEt_;
   vector<double>  *genParQ_;
   vector<double>  *genParId_;
   vector<double>  *genParSt_;
   //   vector<double>  *genJetPx_;
   //   vector<double>  *genJetPy_;
   //   vector<double>  *genJetPz_;
   vector<double>  *genJetE_;
   //   vector<double>  *genJetP_;
   //   vector<double>  *genJetTheta_;
   vector<double>  *genJetPt_;
   vector<double>  *genJetEta_;
   vector<double>  *genJetPhi_;
   //   vector<double>  *genJetEt_;
   //   vector<double>  *genJetQ_;

   // List of branches
   TBranch        *b_PU_weight;   //!
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
   //   TBranch        *b_genParPx_;   //!
   //   TBranch        *b_genParPy_;   //!
   //   TBranch        *b_genParPz_;   //!
   TBranch        *b_genParE_;   //!
   //   TBranch        *b_genParP_;   //!
   //   TBranch        *b_genParTheta_;   //!
   TBranch        *b_genParPt_;   //!
   TBranch        *b_genParEta_;   //!
   TBranch        *b_genParPhi_;   //!
   //   TBranch        *b_genParEt_;   //!
   TBranch        *b_genParQ_;   //!
   TBranch        *b_genParId_;   //!
   TBranch        *b_genParSt_;   //!
   //   TBranch        *b_genJetPx_;   //!
   //   TBranch        *b_genJetPy_;   //!
   //   TBranch        *b_genJetPz_;   //!
   TBranch        *b_genJetE_;   //!
   //   TBranch        *b_genJetP_;   //!
   //   TBranch        *b_genJetTheta_;   //!
   TBranch        *b_genJetPt_;   //!
   TBranch        *b_genJetEta_;   //!
   TBranch        *b_genJetPhi_;   //!
   //   TBranch        *b_genJetEt_;   //!
   //   TBranch        *b_genJetQ_;   //!

   gen_distribution(std::string filename,TTree *tree=0);
   virtual ~gen_distribution();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int lepID, bool applyWeight=true, bool exclusive=false, int DEBUG=0); // 13: muon, 11: electron
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   std::string _inputFileName;


};

#endif

#ifdef gen_distribution_cxx
gen_distribution::gen_distribution(std::string filename, TTree *tree)
{
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename.data());
    if (!f) {
      f = new TFile(filename.data());
      f->cd(Form("%s:/tree",filename.data()));
    }
    tree = (TTree*)gDirectory->Get("tree");

  }
  Init(tree);
  _inputFileName = filename;


}

gen_distribution::~gen_distribution()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gen_distribution::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gen_distribution::LoadTree(Long64_t entry)
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

void gen_distribution::Init(TTree *tree)
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
   //   genParPx_ = 0;
   //   genParPy_ = 0;
   //   genParPz_ = 0;
   genParE_ = 0;
   //   genParP_ = 0;
   //   genParTheta_ = 0;
   genParPt_ = 0;
   genParEta_ = 0;
   genParPhi_ = 0;
   //   genParEt_ = 0;
   genParQ_ = 0;
   genParId_ = 0;
   genParSt_ = 0;
   //   genJetPx_ = 0;
   //   genJetPy_ = 0;
   //   genJetPz_ = 0;
   genJetE_ = 0;
   //   genJetP_ = 0;
   //   genJetTheta_ = 0;
   genJetPt_ = 0;
   genJetEta_ = 0;
   genJetPhi_ = 0;
   //   genJetEt_ = 0;
   //   genJetQ_ = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("PU_weight", &PU_weight, &b_PU_weight);
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
   //   fChain->SetBranchAddress("genParPx_", &genParPx_, &b_genParPx_);
   //   fChain->SetBranchAddress("genParPy_", &genParPy_, &b_genParPy_);
   //   fChain->SetBranchAddress("genParPz_", &genParPz_, &b_genParPz_);
   fChain->SetBranchAddress("genParE_", &genParE_, &b_genParE_);
   //   fChain->SetBranchAddress("genParP_", &genParP_, &b_genParP_);
   //   fChain->SetBranchAddress("genParTheta_", &genParTheta_, &b_genParTheta_);
   fChain->SetBranchAddress("genParPt_", &genParPt_, &b_genParPt_);
   fChain->SetBranchAddress("genParEta_", &genParEta_, &b_genParEta_);
   fChain->SetBranchAddress("genParPhi_", &genParPhi_, &b_genParPhi_);
   //   fChain->SetBranchAddress("genParEt_", &genParEt_, &b_genParEt_);
   fChain->SetBranchAddress("genParQ_", &genParQ_, &b_genParQ_);
   fChain->SetBranchAddress("genParId_", &genParId_, &b_genParId_);
   fChain->SetBranchAddress("genParSt_", &genParSt_, &b_genParSt_);
   //   fChain->SetBranchAddress("genJetPx_", &genJetPx_, &b_genJetPx_);
   //   fChain->SetBranchAddress("genJetPy_", &genJetPy_, &b_genJetPy_);
   //   fChain->SetBranchAddress("genJetPz_", &genJetPz_, &b_genJetPz_);
   fChain->SetBranchAddress("genJetE_", &genJetE_, &b_genJetE_);
   //   fChain->SetBranchAddress("genJetP_", &genJetP_, &b_genJetP_);
   //   fChain->SetBranchAddress("genJetTheta_", &genJetTheta_, &b_genJetTheta_);
   fChain->SetBranchAddress("genJetPt_", &genJetPt_, &b_genJetPt_);
   fChain->SetBranchAddress("genJetEta_", &genJetEta_, &b_genJetEta_);
   fChain->SetBranchAddress("genJetPhi_", &genJetPhi_, &b_genJetPhi_);
   //   fChain->SetBranchAddress("genJetEt_", &genJetEt_, &b_genJetEt_);
   //   fChain->SetBranchAddress("genJetQ_", &genJetQ_, &b_genJetQ_);
   Notify();
}

Bool_t gen_distribution::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gen_distribution::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gen_distribution::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gen_distribution_cxx
