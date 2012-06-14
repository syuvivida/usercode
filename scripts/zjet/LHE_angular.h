//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun 14 19:01:57 2012 by ROOT version 5.27/06b
// from TTree tree/LHE tree
// found on file: gen.root
//////////////////////////////////////////////////////////

#ifndef LHE_angular_h
#define LHE_angular_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream>
#include <string>

using namespace std;


   const Int_t kMaxgenParId = 1;
   const Int_t kMaxgenParSt = 1;
   const Int_t kMaxgenParPx = 1;
   const Int_t kMaxgenParPy = 1;
   const Int_t kMaxgenParPz = 1;
   const Int_t kMaxgenParE = 1;

class LHE_angular {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   vector<double>  *genParId_;
   vector<double>  *genParSt_;
   vector<double>  *genParPx_;
   vector<double>  *genParPy_;
   vector<double>  *genParPz_;
   vector<double>  *genParE_;

   // List of branches
   TBranch        *b_genParId_;   //!
   TBranch        *b_genParSt_;   //!
   TBranch        *b_genParPx_;   //!
   TBranch        *b_genParPy_;   //!
   TBranch        *b_genParPz_;   //!
   TBranch        *b_genParE_;   //!

   LHE_angular(std::string filename,TTree *tree=0);
   virtual ~LHE_angular();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int lepID, bool exclusive=true);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   std::string _inputFileName;
};

#endif

#ifdef LHE_angular_cxx
LHE_angular::LHE_angular(std::string filename,TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->
	FindObject(Form("%s",filename.data()));
      if (!f) {
	f = new TFile(Form("%s",filename.data()));
	f->cd(Form("%s:/dummy",filename.data()));
      }
      tree = (TTree*)gDirectory->Get("tree");

   }
   Init(tree);
   _inputFileName = filename;
}

LHE_angular::~LHE_angular()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t LHE_angular::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t LHE_angular::LoadTree(Long64_t entry)
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

void LHE_angular::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   genParId_ = 0;
   genParSt_ = 0;
   genParPx_ = 0;
   genParPy_ = 0;
   genParPz_ = 0;
   genParE_ = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("genParId_", &genParId_, &b_genParId_);
   fChain->SetBranchAddress("genParSt_", &genParSt_, &b_genParSt_);
   fChain->SetBranchAddress("genParPx_", &genParPx_, &b_genParPx_);
   fChain->SetBranchAddress("genParPy_", &genParPy_, &b_genParPy_);
   fChain->SetBranchAddress("genParPz_", &genParPz_, &b_genParPz_);
   fChain->SetBranchAddress("genParE_", &genParE_, &b_genParE_);
   Notify();
}

Bool_t LHE_angular::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void LHE_angular::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t LHE_angular::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef LHE_angular_cxx
