//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Sep  3 11:03:19 2010 by ROOT version 5.22/00d
// from TTree photons/reco photon info
// found on file: mc/randomPhotonNtuple_mc.root
//////////////////////////////////////////////////////////

#ifndef makeRCHisto_h
#define makeRCHisto_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <string>

using namespace std;

class makeRCHisto {
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
   Float_t         nVtx;
   Float_t         nGoodVtx3;
   Float_t         nGoodVtx4;
   Float_t         nPix;
   Float_t         hfP;
   Float_t         hfM;
   Float_t         hf;
   Float_t         run;

   // List of branches
   TBranch        *b_et;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_tkiso;   //!
   TBranch        *b_eciso;   //!
   TBranch        *b_hciso;   //!
   TBranch        *b_h1iso;   //!
   TBranch        *b_h2iso;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_nGoodVtx3;   //!
   TBranch        *b_nGoodVtx4;   //!
   TBranch        *b_nPix;   //!
   TBranch        *b_hfP;   //!
   TBranch        *b_hfM;   //!
   TBranch        *b_hf;   //!
   TBranch        *b_run;   //!

/*    makeRCHisto(TTree *tree=0); */
   makeRCHisto(std::string inputDirName, std::string dirName="randomPhotonAnalyzer", TTree *tree=0);
   std::string    _outputName;
   std::string    _inputDirName;
   virtual ~makeRCHisto();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef makeRCHisto_cxx
makeRCHisto::makeRCHisto(std::string inputDirName, std::string dirName, TTree *tree)
{
  _outputName   = dirName;
  _inputDirName = inputDirName;
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    std::string inputFileName = _inputDirName + "/randomPhotonNtuple.root";
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(inputFileName.data());
    //    if (!f) {
      f = new TFile(inputFileName.data());
      std::string tempName = inputFileName+":/" + dirName;
      std::cout << tempName << std::endl;
      f->cd(tempName.data());
      //    }
    tree = (TTree*)gDirectory->Get("photons");
  }
  Init(tree);

}

makeRCHisto::~makeRCHisto()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t makeRCHisto::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t makeRCHisto::LoadTree(Long64_t entry)
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

void makeRCHisto::Init(TTree *tree)
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
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("nGoodVtx3", &nGoodVtx3, &b_nGoodVtx3);
   fChain->SetBranchAddress("nGoodVtx4", &nGoodVtx4, &b_nGoodVtx4);
   fChain->SetBranchAddress("nPix", &nPix, &b_nPix);
   fChain->SetBranchAddress("hfP", &hfP, &b_hfP);
   fChain->SetBranchAddress("hfM", &hfM, &b_hfM);
   fChain->SetBranchAddress("hf", &hf, &b_hf);
   fChain->SetBranchAddress("run", &run, &b_run);
   Notify();
}

Bool_t makeRCHisto::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void makeRCHisto::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t makeRCHisto::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef makeRCHisto_cxx
