//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Aug 25 17:51:06 2011 by ROOT version 5.30/00
// from TTree tree/tree
// found on file: patTree.root
//////////////////////////////////////////////////////////

#ifndef vectorAngular_h
#define vectorAngular_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string>
#include <iostream>

using namespace std;
   const Int_t kMaxpatJetPFlowPt = 1;
   const Int_t kMaxpatJetPFlowEta = 1;
   const Int_t kMaxpatJetPFlowPhi = 1;
   const Int_t kMaxpatJetPFlowM = 1;
   const Int_t kMaxpatJetPFlowRapidity = 1;
   const Int_t kMaxpatJetPFlowPx = 1;
   const Int_t kMaxpatJetPFlowPy = 1;
   const Int_t kMaxpatJetPFlowPz = 1;
   const Int_t kMaxpatJetPFlowEn = 1;
   const Int_t kMaxpatJetPFlowUnCorrPt = 1;
   const Int_t kMaxpatJetPFlowUnCorrPx = 1;
   const Int_t kMaxpatJetPFlowUnCorrPy = 1;
   const Int_t kMaxpatJetPFlowUnCorrPz = 1;
   const Int_t kMaxpatJetPFlowUnCorrEnt = 1;
   const Int_t kMaxpatJetPFlowPhotEn = 1;
   const Int_t kMaxpatJetPFlowElecEn = 1;
   const Int_t kMaxpatJetPFlowMuonEn = 1;
   const Int_t kMaxpatJetPFlowHfHadEn = 1;
   const Int_t kMaxpatJetPFlowHfEmEn = 1;
   const Int_t kMaxpatJetPFlowCharHadE = 1;
   const Int_t kMaxpatJetPFlowNeutHadE = 1;
   const Int_t kMaxpatJetPFlowCharEmE = 1;
   const Int_t kMaxpatJetPFlowCharMuE = 1;
   const Int_t kMaxpatJetPFlowNeutEmE = 1;
   const Int_t kMaxpatJetPFlowMuonMulti = 1;
   const Int_t kMaxpatJetPFlowNeutMulti = 1;
   const Int_t kMaxpatJetPFlowCharMulti = 1;
   const Int_t kMaxpatJetPFlowCharHadMulti = 1;
   const Int_t kMaxpatJetPFlowNeutHadMulti = 1;
   const Int_t kMaxpatJetPFlowPhotMulti = 1;
   const Int_t kMaxpatJetPFlowElecMulti = 1;
   const Int_t kMaxpatJetPFlowHfHadMulti = 1;
   const Int_t kMaxpatJetPFlowHfEmMulti = 1;
   const Int_t kMaxpatJetPFlowPhotEnFr = 1;
   const Int_t kMaxpatJetPFlowMuonEnFr = 1;
   const Int_t kMaxpatJetPFlowHfHadEnFr = 1;
   const Int_t kMaxpatJetPFlowHfEmEnFr = 1;
   const Int_t kMaxpatJetPFlowNeutEmEFr = 1;
   const Int_t kMaxpatJetPFlowCharHadEFr = 1;
   const Int_t kMaxpatJetPFlowNeutHadEFr = 1;
   const Int_t kMaxpatJetPFlowCharEmEFr = 1;
   const Int_t kMaxpatJetPFlowCharMuEFr = 1;

class vectorAngular {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           EvtInfo_EventNum;
   Int_t           EvtInfo_RunNum;
   Int_t           EvtInfo_LumiSection;
   Int_t           EvtInfo_BunchXing;
   vector<double>  *genMuPx;
   vector<double>  *genMuPy;
   vector<double>  *genMuPz;
   vector<double>  *genMuE;
   vector<double>  *genMuP;
   vector<double>  *genMuTheta;
   vector<double>  *genMuPt;
   vector<double>  *genMuEta;
   vector<double>  *genMuPhi;
   vector<double>  *genMuEt;
   vector<double>  *genMuQ;
   vector<double>  *genElePx;
   vector<double>  *genElePy;
   vector<double>  *genElePz;
   vector<double>  *genEleE;
   vector<double>  *genEleP;
   vector<double>  *genEleTheta;
   vector<double>  *genElePt;
   vector<double>  *genEleEta;
   vector<double>  *genElePhi;
   vector<double>  *genEleEt;
   vector<double>  *genEleQ;
   vector<double>  *genPhoPx;
   vector<double>  *genPhoPy;
   vector<double>  *genPhoPz;
   vector<double>  *genPhoE;
   vector<double>  *genPhoP;
   vector<double>  *genPhoTheta;
   vector<double>  *genPhoPt;
   vector<double>  *genPhoEta;
   vector<double>  *genPhoPhi;
   vector<double>  *genPhoEt;
   vector<int>     *genPhoMomPID;
   vector<double>  *genJetPx;
   vector<double>  *genJetPy;
   vector<double>  *genJetPz;
   vector<double>  *genJetE;
   vector<double>  *genJetP;
   vector<double>  *genJetTheta;
   vector<double>  *genJetPt;
   vector<double>  *genJetEta;
   vector<double>  *genJetPhi;
   vector<double>  *genJetEt;
   vector<double>  *genJetQ;
   vector<double>  *patJetPFlowPt_;
   vector<double>  *patJetPFlowEta_;
   vector<double>  *patJetPFlowPhi_;
   vector<double>  *patJetPFlowM_;
   vector<double>  *patJetPFlowRapidity_;
   vector<double>  *patJetPFlowPx_;
   vector<double>  *patJetPFlowPy_;
   vector<double>  *patJetPFlowPz_;
   vector<double>  *patJetPFlowEn_;
   vector<double>  *patJetPFlowUnCorrPt_;
   vector<double>  *patJetPFlowUnCorrPx_;
   vector<double>  *patJetPFlowUnCorrPy_;
   vector<double>  *patJetPFlowUnCorrPz_;
   vector<double>  *patJetPFlowUnCorrEnt_;
   vector<double>  *patJetPFlowPhotEn_;
   vector<double>  *patJetPFlowElecEn_;
   vector<double>  *patJetPFlowMuonEn_;
   vector<double>  *patJetPFlowHfHadEn_;
   vector<double>  *patJetPFlowHfEmEn_;
   vector<double>  *patJetPFlowCharHadE_;
   vector<double>  *patJetPFlowNeutHadE_;
   vector<double>  *patJetPFlowCharEmE_;
   vector<double>  *patJetPFlowCharMuE_;
   vector<double>  *patJetPFlowNeutEmE_;
   vector<double>  *patJetPFlowMuonMulti_;
   vector<double>  *patJetPFlowNeutMulti_;
   vector<double>  *patJetPFlowCharMulti_;
   vector<double>  *patJetPFlowCharHadMulti_;
   vector<double>  *patJetPFlowNeutHadMulti_;
   vector<double>  *patJetPFlowPhotMulti_;
   vector<double>  *patJetPFlowElecMulti_;
   vector<double>  *patJetPFlowHfHadMulti_;
   vector<double>  *patJetPFlowHfEmMulti_;
   vector<double>  *patJetPFlowPhotEnFr_;
   vector<double>  *patJetPFlowMuonEnFr_;
   vector<double>  *patJetPFlowHfHadEnFr_;
   vector<double>  *patJetPFlowHfEmEnFr_;
   vector<double>  *patJetPFlowNeutEmEFr_;
   vector<double>  *patJetPFlowCharHadEFr_;
   vector<double>  *patJetPFlowNeutHadEFr_;
   vector<double>  *patJetPFlowCharEmEFr_;
   vector<double>  *patJetPFlowCharMuEFr_;
   vector<int>     *trigResults;
   vector<string>  *trigName;

   // List of branches
   TBranch        *b_EvtInfo_EventNum;   //!
   TBranch        *b_EvtInfo_RunNum;   //!
   TBranch        *b_EvtInfo_LumiSection;   //!
   TBranch        *b_EvtInfo_BunchXing;   //!
   TBranch        *b_genMuPx;   //!
   TBranch        *b_genMuPy;   //!
   TBranch        *b_genMuPz;   //!
   TBranch        *b_genMuE;   //!
   TBranch        *b_genMuP;   //!
   TBranch        *b_genMuTheta;   //!
   TBranch        *b_genMuPt;   //!
   TBranch        *b_genMuEta;   //!
   TBranch        *b_genMuPhi;   //!
   TBranch        *b_genMuEt;   //!
   TBranch        *b_genMuQ;   //!
   TBranch        *b_genElePx;   //!
   TBranch        *b_genElePy;   //!
   TBranch        *b_genElePz;   //!
   TBranch        *b_genEleE;   //!
   TBranch        *b_genEleP;   //!
   TBranch        *b_genEleTheta;   //!
   TBranch        *b_genElePt;   //!
   TBranch        *b_genEleEta;   //!
   TBranch        *b_genElePhi;   //!
   TBranch        *b_genEleEt;   //!
   TBranch        *b_genEleQ;   //!
   TBranch        *b_genPhoPx;   //!
   TBranch        *b_genPhoPy;   //!
   TBranch        *b_genPhoPz;   //!
   TBranch        *b_genPhoE;   //!
   TBranch        *b_genPhoP;   //!
   TBranch        *b_genPhoTheta;   //!
   TBranch        *b_genPhoPt;   //!
   TBranch        *b_genPhoEta;   //!
   TBranch        *b_genPhoPhi;   //!
   TBranch        *b_genPhoEt;   //!
   TBranch        *b_genPhoMomPID;   //!
   TBranch        *b_genJetPx;   //!
   TBranch        *b_genJetPy;   //!
   TBranch        *b_genJetPz;   //!
   TBranch        *b_genJetE;   //!
   TBranch        *b_genJetP;   //!
   TBranch        *b_genJetTheta;   //!
   TBranch        *b_genJetPt;   //!
   TBranch        *b_genJetEta;   //!
   TBranch        *b_genJetPhi;   //!
   TBranch        *b_genJetEt;   //!
   TBranch        *b_genJetQ;   //!
   TBranch        *b_patJetPFlowPt_;   //!
   TBranch        *b_patJetPFlowEta_;   //!
   TBranch        *b_patJetPFlowPhi_;   //!
   TBranch        *b_patJetPFlowM_;   //!
   TBranch        *b_patJetPFlowRapidity_;   //!
   TBranch        *b_patJetPFlowPx_;   //!
   TBranch        *b_patJetPFlowPy_;   //!
   TBranch        *b_patJetPFlowPz_;   //!
   TBranch        *b_patJetPFlowEn_;   //!
   TBranch        *b_patJetPFlowUnCorrPt_;   //!
   TBranch        *b_patJetPFlowUnCorrPx_;   //!
   TBranch        *b_patJetPFlowUnCorrPy_;   //!
   TBranch        *b_patJetPFlowUnCorrPz_;   //!
   TBranch        *b_patJetPFlowUnCorrEnt_;   //!
   TBranch        *b_patJetPFlowPhotEn_;   //!
   TBranch        *b_patJetPFlowElecEn_;   //!
   TBranch        *b_patJetPFlowMuonEn_;   //!
   TBranch        *b_patJetPFlowHfHadEn_;   //!
   TBranch        *b_patJetPFlowHfEmEn_;   //!
   TBranch        *b_patJetPFlowCharHadE_;   //!
   TBranch        *b_patJetPFlowNeutHadE_;   //!
   TBranch        *b_patJetPFlowCharEmE_;   //!
   TBranch        *b_patJetPFlowCharMuE_;   //!
   TBranch        *b_patJetPFlowNeutEmE_;   //!
   TBranch        *b_patJetPFlowMuonMulti_;   //!
   TBranch        *b_patJetPFlowNeutMulti_;   //!
   TBranch        *b_patJetPFlowCharMulti_;   //!
   TBranch        *b_patJetPFlowCharHadMulti_;   //!
   TBranch        *b_patJetPFlowNeutHadMulti_;   //!
   TBranch        *b_patJetPFlowPhotMulti_;   //!
   TBranch        *b_patJetPFlowElecMulti_;   //!
   TBranch        *b_patJetPFlowHfHadMulti_;   //!
   TBranch        *b_patJetPFlowHfEmMulti_;   //!
   TBranch        *b_patJetPFlowPhotEnFr_;   //!
   TBranch        *b_patJetPFlowMuonEnFr_;   //!
   TBranch        *b_patJetPFlowHfHadEnFr_;   //!
   TBranch        *b_patJetPFlowHfEmEnFr_;   //!
   TBranch        *b_patJetPFlowNeutEmEFr_;   //!
   TBranch        *b_patJetPFlowCharHadEFr_;   //!
   TBranch        *b_patJetPFlowNeutHadEFr_;   //!
   TBranch        *b_patJetPFlowCharEmEFr_;   //!
   TBranch        *b_patJetPFlowCharMuEFr_;   //!
   TBranch        *b_trigResults;   //!
   TBranch        *b_trigName;   //!

   vectorAngular(std::string inputFile);
   virtual ~vectorAngular();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Float_t pstarmin=-9999.0, Float_t pstarmax=9999.0,
			 Float_t ybmin = 0.0, Float_t ybmax=9999.0);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   std::string  inputFile_;

};

#endif

#ifdef vectorAngular_cxx
vectorAngular::vectorAngular(std::string inputFile)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   cout << "Reading " << inputFile.data() << endl;
   TTree *tree;
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(inputFile.data());
   if (!f || !f->IsOpen()) {
     f = new TFile(inputFile.data());
   }
   f->GetObject("tree",tree);
   Init(tree);
   inputFile_ = inputFile;
   
}

vectorAngular::~vectorAngular()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t vectorAngular::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t vectorAngular::LoadTree(Long64_t entry)
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

void vectorAngular::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   genMuPx = 0;
   genMuPy = 0;
   genMuPz = 0;
   genMuE = 0;
   genMuP = 0;
   genMuTheta = 0;
   genMuPt = 0;
   genMuEta = 0;
   genMuPhi = 0;
   genMuEt = 0;
   genMuQ = 0;
   genElePx = 0;
   genElePy = 0;
   genElePz = 0;
   genEleE = 0;
   genEleP = 0;
   genEleTheta = 0;
   genElePt = 0;
   genEleEta = 0;
   genElePhi = 0;
   genEleEt = 0;
   genEleQ = 0;
   genPhoPx = 0;
   genPhoPy = 0;
   genPhoPz = 0;
   genPhoE = 0;
   genPhoP = 0;
   genPhoTheta = 0;
   genPhoPt = 0;
   genPhoEta = 0;
   genPhoPhi = 0;
   genPhoEt = 0;
   genPhoMomPID = 0;
   genJetPx = 0;
   genJetPy = 0;
   genJetPz = 0;
   genJetE = 0;
   genJetP = 0;
   genJetTheta = 0;
   genJetPt = 0;
   genJetEta = 0;
   genJetPhi = 0;
   genJetEt = 0;
   genJetQ = 0;
   patJetPFlowPt_ = 0;
   patJetPFlowEta_ = 0;
   patJetPFlowPhi_ = 0;
   patJetPFlowM_ = 0;
   patJetPFlowRapidity_ = 0;
   patJetPFlowPx_ = 0;
   patJetPFlowPy_ = 0;
   patJetPFlowPz_ = 0;
   patJetPFlowEn_ = 0;
   patJetPFlowUnCorrPt_ = 0;
   patJetPFlowUnCorrPx_ = 0;
   patJetPFlowUnCorrPy_ = 0;
   patJetPFlowUnCorrPz_ = 0;
   patJetPFlowUnCorrEnt_ = 0;
   patJetPFlowPhotEn_ = 0;
   patJetPFlowElecEn_ = 0;
   patJetPFlowMuonEn_ = 0;
   patJetPFlowHfHadEn_ = 0;
   patJetPFlowHfEmEn_ = 0;
   patJetPFlowCharHadE_ = 0;
   patJetPFlowNeutHadE_ = 0;
   patJetPFlowCharEmE_ = 0;
   patJetPFlowCharMuE_ = 0;
   patJetPFlowNeutEmE_ = 0;
   patJetPFlowMuonMulti_ = 0;
   patJetPFlowNeutMulti_ = 0;
   patJetPFlowCharMulti_ = 0;
   patJetPFlowCharHadMulti_ = 0;
   patJetPFlowNeutHadMulti_ = 0;
   patJetPFlowPhotMulti_ = 0;
   patJetPFlowElecMulti_ = 0;
   patJetPFlowHfHadMulti_ = 0;
   patJetPFlowHfEmMulti_ = 0;
   patJetPFlowPhotEnFr_ = 0;
   patJetPFlowMuonEnFr_ = 0;
   patJetPFlowHfHadEnFr_ = 0;
   patJetPFlowHfEmEnFr_ = 0;
   patJetPFlowNeutEmEFr_ = 0;
   patJetPFlowCharHadEFr_ = 0;
   patJetPFlowNeutHadEFr_ = 0;
   patJetPFlowCharEmEFr_ = 0;
   patJetPFlowCharMuEFr_ = 0;
   trigResults = 0;
   trigName = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EvtInfo_EventNum", &EvtInfo_EventNum, &b_EvtInfo_EventNum);
   fChain->SetBranchAddress("EvtInfo_RunNum", &EvtInfo_RunNum, &b_EvtInfo_RunNum);
   fChain->SetBranchAddress("EvtInfo_LumiSection", &EvtInfo_LumiSection, &b_EvtInfo_LumiSection);
   fChain->SetBranchAddress("EvtInfo_BunchXing", &EvtInfo_BunchXing, &b_EvtInfo_BunchXing);
   fChain->SetBranchAddress("genMuPx", &genMuPx, &b_genMuPx);
   fChain->SetBranchAddress("genMuPy", &genMuPy, &b_genMuPy);
   fChain->SetBranchAddress("genMuPz", &genMuPz, &b_genMuPz);
   fChain->SetBranchAddress("genMuE", &genMuE, &b_genMuE);
   fChain->SetBranchAddress("genMuP", &genMuP, &b_genMuP);
   fChain->SetBranchAddress("genMuTheta", &genMuTheta, &b_genMuTheta);
   fChain->SetBranchAddress("genMuPt", &genMuPt, &b_genMuPt);
   fChain->SetBranchAddress("genMuEta", &genMuEta, &b_genMuEta);
   fChain->SetBranchAddress("genMuPhi", &genMuPhi, &b_genMuPhi);
   fChain->SetBranchAddress("genMuEt", &genMuEt, &b_genMuEt);
   fChain->SetBranchAddress("genMuQ", &genMuQ, &b_genMuQ);
   fChain->SetBranchAddress("genElePx", &genElePx, &b_genElePx);
   fChain->SetBranchAddress("genElePy", &genElePy, &b_genElePy);
   fChain->SetBranchAddress("genElePz", &genElePz, &b_genElePz);
   fChain->SetBranchAddress("genEleE", &genEleE, &b_genEleE);
   fChain->SetBranchAddress("genEleP", &genEleP, &b_genEleP);
   fChain->SetBranchAddress("genEleTheta", &genEleTheta, &b_genEleTheta);
   fChain->SetBranchAddress("genElePt", &genElePt, &b_genElePt);
   fChain->SetBranchAddress("genEleEta", &genEleEta, &b_genEleEta);
   fChain->SetBranchAddress("genElePhi", &genElePhi, &b_genElePhi);
   fChain->SetBranchAddress("genEleEt", &genEleEt, &b_genEleEt);
   fChain->SetBranchAddress("genEleQ", &genEleQ, &b_genEleQ);
   fChain->SetBranchAddress("genPhoPx", &genPhoPx, &b_genPhoPx);
   fChain->SetBranchAddress("genPhoPy", &genPhoPy, &b_genPhoPy);
   fChain->SetBranchAddress("genPhoPz", &genPhoPz, &b_genPhoPz);
   fChain->SetBranchAddress("genPhoE", &genPhoE, &b_genPhoE);
   fChain->SetBranchAddress("genPhoP", &genPhoP, &b_genPhoP);
   fChain->SetBranchAddress("genPhoTheta", &genPhoTheta, &b_genPhoTheta);
   fChain->SetBranchAddress("genPhoPt", &genPhoPt, &b_genPhoPt);
   fChain->SetBranchAddress("genPhoEta", &genPhoEta, &b_genPhoEta);
   fChain->SetBranchAddress("genPhoPhi", &genPhoPhi, &b_genPhoPhi);
   fChain->SetBranchAddress("genPhoEt", &genPhoEt, &b_genPhoEt);
   fChain->SetBranchAddress("genPhoMomPID", &genPhoMomPID, &b_genPhoMomPID);
   fChain->SetBranchAddress("genJetPx", &genJetPx, &b_genJetPx);
   fChain->SetBranchAddress("genJetPy", &genJetPy, &b_genJetPy);
   fChain->SetBranchAddress("genJetPz", &genJetPz, &b_genJetPz);
   fChain->SetBranchAddress("genJetE", &genJetE, &b_genJetE);
   fChain->SetBranchAddress("genJetP", &genJetP, &b_genJetP);
   fChain->SetBranchAddress("genJetTheta", &genJetTheta, &b_genJetTheta);
   fChain->SetBranchAddress("genJetPt", &genJetPt, &b_genJetPt);
   fChain->SetBranchAddress("genJetEta", &genJetEta, &b_genJetEta);
   fChain->SetBranchAddress("genJetPhi", &genJetPhi, &b_genJetPhi);
   fChain->SetBranchAddress("genJetEt", &genJetEt, &b_genJetEt);
   fChain->SetBranchAddress("genJetQ", &genJetQ, &b_genJetQ);
   fChain->SetBranchAddress("patJetPFlowPt_", &patJetPFlowPt_, &b_patJetPFlowPt_);
   fChain->SetBranchAddress("patJetPFlowEta_", &patJetPFlowEta_, &b_patJetPFlowEta_);
   fChain->SetBranchAddress("patJetPFlowPhi_", &patJetPFlowPhi_, &b_patJetPFlowPhi_);
   fChain->SetBranchAddress("patJetPFlowM_", &patJetPFlowM_, &b_patJetPFlowM_);
   fChain->SetBranchAddress("patJetPFlowRapidity_", &patJetPFlowRapidity_, &b_patJetPFlowRapidity_);
   fChain->SetBranchAddress("patJetPFlowPx_", &patJetPFlowPx_, &b_patJetPFlowPx_);
   fChain->SetBranchAddress("patJetPFlowPy_", &patJetPFlowPy_, &b_patJetPFlowPy_);
   fChain->SetBranchAddress("patJetPFlowPz_", &patJetPFlowPz_, &b_patJetPFlowPz_);
   fChain->SetBranchAddress("patJetPFlowEn_", &patJetPFlowEn_, &b_patJetPFlowEn_);
   fChain->SetBranchAddress("patJetPFlowUnCorrPt_", &patJetPFlowUnCorrPt_, &b_patJetPFlowUnCorrPt_);
   fChain->SetBranchAddress("patJetPFlowUnCorrPx_", &patJetPFlowUnCorrPx_, &b_patJetPFlowUnCorrPx_);
   fChain->SetBranchAddress("patJetPFlowUnCorrPy_", &patJetPFlowUnCorrPy_, &b_patJetPFlowUnCorrPy_);
   fChain->SetBranchAddress("patJetPFlowUnCorrPz_", &patJetPFlowUnCorrPz_, &b_patJetPFlowUnCorrPz_);
   fChain->SetBranchAddress("patJetPFlowUnCorrEnt_", &patJetPFlowUnCorrEnt_, &b_patJetPFlowUnCorrEnt_);
   fChain->SetBranchAddress("patJetPFlowPhotEn_", &patJetPFlowPhotEn_, &b_patJetPFlowPhotEn_);
   fChain->SetBranchAddress("patJetPFlowElecEn_", &patJetPFlowElecEn_, &b_patJetPFlowElecEn_);
   fChain->SetBranchAddress("patJetPFlowMuonEn_", &patJetPFlowMuonEn_, &b_patJetPFlowMuonEn_);
   fChain->SetBranchAddress("patJetPFlowHfHadEn_", &patJetPFlowHfHadEn_, &b_patJetPFlowHfHadEn_);
   fChain->SetBranchAddress("patJetPFlowHfEmEn_", &patJetPFlowHfEmEn_, &b_patJetPFlowHfEmEn_);
   fChain->SetBranchAddress("patJetPFlowCharHadE_", &patJetPFlowCharHadE_, &b_patJetPFlowCharHadE_);
   fChain->SetBranchAddress("patJetPFlowNeutHadE_", &patJetPFlowNeutHadE_, &b_patJetPFlowNeutHadE_);
   fChain->SetBranchAddress("patJetPFlowCharEmE_", &patJetPFlowCharEmE_, &b_patJetPFlowCharEmE_);
   fChain->SetBranchAddress("patJetPFlowCharMuE_", &patJetPFlowCharMuE_, &b_patJetPFlowCharMuE_);
   fChain->SetBranchAddress("patJetPFlowNeutEmE_", &patJetPFlowNeutEmE_, &b_patJetPFlowNeutEmE_);
   fChain->SetBranchAddress("patJetPFlowMuonMulti_", &patJetPFlowMuonMulti_, &b_patJetPFlowMuonMulti_);
   fChain->SetBranchAddress("patJetPFlowNeutMulti_", &patJetPFlowNeutMulti_, &b_patJetPFlowNeutMulti_);
   fChain->SetBranchAddress("patJetPFlowCharMulti_", &patJetPFlowCharMulti_, &b_patJetPFlowCharMulti_);
   fChain->SetBranchAddress("patJetPFlowCharHadMulti_", &patJetPFlowCharHadMulti_, &b_patJetPFlowCharHadMulti_);
   fChain->SetBranchAddress("patJetPFlowNeutHadMulti_", &patJetPFlowNeutHadMulti_, &b_patJetPFlowNeutHadMulti_);
   fChain->SetBranchAddress("patJetPFlowPhotMulti_", &patJetPFlowPhotMulti_, &b_patJetPFlowPhotMulti_);
   fChain->SetBranchAddress("patJetPFlowElecMulti_", &patJetPFlowElecMulti_, &b_patJetPFlowElecMulti_);
   fChain->SetBranchAddress("patJetPFlowHfHadMulti_", &patJetPFlowHfHadMulti_, &b_patJetPFlowHfHadMulti_);
   fChain->SetBranchAddress("patJetPFlowHfEmMulti_", &patJetPFlowHfEmMulti_, &b_patJetPFlowHfEmMulti_);
   fChain->SetBranchAddress("patJetPFlowPhotEnFr_", &patJetPFlowPhotEnFr_, &b_patJetPFlowPhotEnFr_);
   fChain->SetBranchAddress("patJetPFlowMuonEnFr_", &patJetPFlowMuonEnFr_, &b_patJetPFlowMuonEnFr_);
   fChain->SetBranchAddress("patJetPFlowHfHadEnFr_", &patJetPFlowHfHadEnFr_, &b_patJetPFlowHfHadEnFr_);
   fChain->SetBranchAddress("patJetPFlowHfEmEnFr_", &patJetPFlowHfEmEnFr_, &b_patJetPFlowHfEmEnFr_);
   fChain->SetBranchAddress("patJetPFlowNeutEmEFr_", &patJetPFlowNeutEmEFr_, &b_patJetPFlowNeutEmEFr_);
   fChain->SetBranchAddress("patJetPFlowCharHadEFr_", &patJetPFlowCharHadEFr_, &b_patJetPFlowCharHadEFr_);
   fChain->SetBranchAddress("patJetPFlowNeutHadEFr_", &patJetPFlowNeutHadEFr_, &b_patJetPFlowNeutHadEFr_);
   fChain->SetBranchAddress("patJetPFlowCharEmEFr_", &patJetPFlowCharEmEFr_, &b_patJetPFlowCharEmEFr_);
   fChain->SetBranchAddress("patJetPFlowCharMuEFr_", &patJetPFlowCharMuEFr_, &b_patJetPFlowCharMuEFr_);
   fChain->SetBranchAddress("trigResults", &trigResults, &b_trigResults);
   fChain->SetBranchAddress("trigName", &trigName, &b_trigName);
   Notify();
}

Bool_t vectorAngular::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void vectorAngular::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t vectorAngular::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef vectorAngular_cxx
