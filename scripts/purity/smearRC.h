//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 22 22:22:18 2010 by ROOT version 5.22/00d
// from TTree photons/reco photon info
// found on file: randomConeNtuple_data.root
//////////////////////////////////////////////////////////

#ifndef smearRC_h
#define smearRC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class smearRC {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         et;
   Float_t         eta;
   Float_t         phi;
   Float_t         tkiso;
   Float_t         eciso;
   Float_t         hciso;
   Float_t         h1iso;
   Float_t         h2iso;

   // List of branches
   TBranch        *b_et;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_tkiso;   //!
   TBranch        *b_eciso;   //!
   TBranch        *b_hciso;   //!
   TBranch        *b_h1iso;   //!
   TBranch        *b_h2iso;   //!

   smearRC(TTree *tree=0);
   virtual ~smearRC();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int ieta=0);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef smearRC_cxx
smearRC::smearRC(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("randomConeNtuple_data.root");
      if (!f) {
         f = new TFile("randomConeNtuple_data.root");
         f->cd("randomConeNtuple_data.root:/ana");
      }
      tree = (TTree*)gDirectory->Get("photons");

   }
   Init(tree);
}

smearRC::~smearRC()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t smearRC::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t smearRC::LoadTree(Long64_t entry)
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

void smearRC::Init(TTree *tree)
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

   fChain->SetBranchAddress("et", &et, &b_et);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("tkiso", &tkiso, &b_tkiso);
   fChain->SetBranchAddress("eciso", &eciso, &b_eciso);
   fChain->SetBranchAddress("hciso", &hciso, &b_hciso);
   fChain->SetBranchAddress("h1iso", &h1iso, &b_h1iso);
   fChain->SetBranchAddress("h2iso", &h2iso, &b_h2iso);
   Notify();
}

Bool_t smearRC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void smearRC::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t smearRC::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef smearRC_cxx
