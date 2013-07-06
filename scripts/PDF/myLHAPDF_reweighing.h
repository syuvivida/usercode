//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul  5 11:08:22 2013 by ROOT version 5.27/06b
// from TTree tree/tree
// found on file: /scratch/syu/submit_PDF/CMSSW_4_4_4/src/CRAB/zjets_mc.root
//////////////////////////////////////////////////////////

#ifndef myLHAPDF_reweighing_h
#define myLHAPDF_reweighing_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <iostream>
#include <string>
/*
#include <TBranch.h>
#include <TLeaf.h>
#include <TSystemDirectory.h>
#include <TList.h>
*/
#include "LHAPDF/LHAPDF.h"
using namespace std;


const Int_t kMaxEvtInfo_VertexX = 1;
const Int_t kMaxEvtInfo_VertexY = 1;
const Int_t kMaxEvtInfo_VertexZ = 1;
const Int_t kMaxptHat = 1;
const Int_t kMaxmcWeight = 1;
const Int_t kMaxpdfInfo = 1;
const Int_t kMaxcteq66PDFw = 1;
const Int_t kMaxcteq6l1PDFw = 1;
const Int_t kMaxmstw2008PDFw = 1;
const Int_t kMaxgenParE = 1;
const Int_t kMaxgenParPt = 1;
const Int_t kMaxgenParEta = 1;
const Int_t kMaxgenParPhi = 1;
const Int_t kMaxgenParQ = 1;
const Int_t kMaxgenParId = 1;
const Int_t kMaxgenParSt = 1;
const Int_t kMaxgenMomParId = 1;
const Int_t kMaxgenParIndex = 1;
const Int_t kMaxgenLepE = 1;
const Int_t kMaxgenLepPt = 1;
const Int_t kMaxgenLepEta = 1;
const Int_t kMaxgenLepPhi = 1;
const Int_t kMaxgenLepPhoE = 1;
const Int_t kMaxgenLepPhoPt = 1;
const Int_t kMaxgenLepPhoEta = 1;
const Int_t kMaxgenLepPhoPhi = 1;
const Int_t kMaxgenLepId = 1;
const Int_t kMaxgenJetE = 1;
const Int_t kMaxgenJetPt = 1;
const Int_t kMaxgenJetEta = 1;
const Int_t kMaxgenJetPhi = 1;


class myLHAPDF_reweighing {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
  Double_t        PU_weight;
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
  vector<double>  *pdfInfo_;
  vector<double>  *cteq66PDFw_;
  vector<double>  *cteq6l1PDFw_;
  vector<double>  *mstw2008PDFw_;
  vector<double>  *genParE_;
  vector<double>  *genParPt_;
  vector<double>  *genParEta_;
  vector<double>  *genParPhi_;
  vector<int>     *genParQ_;
  vector<int>     *genParId_;
  vector<int>     *genParSt_;
  vector<int>     *genMomParId_;
  vector<int>     *genParIndex_;
  vector<double>  *genLepE_;
  vector<double>  *genLepPt_;
  vector<double>  *genLepEta_;
  vector<double>  *genLepPhi_;
  vector<double>  *genLepPhoE_;
  vector<double>  *genLepPhoPt_;
  vector<double>  *genLepPhoEta_;
  vector<double>  *genLepPhoPhi_;
  vector<int>     *genLepId_;
  vector<double>  *genJetE_;
  vector<double>  *genJetPt_;
  vector<double>  *genJetEta_;
  vector<double>  *genJetPhi_;

  // List of branches
  TBranch        *b_PU_weight;   //!
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
  TBranch        *b_pdfInfo_;   //!
  TBranch        *b_cteq66PDFw_;   //!
  TBranch        *b_cteq6l1PDFw_;   //!
  TBranch        *b_mstw2008PDFw_;   //!
  TBranch        *b_genParE_;   //!
  TBranch        *b_genParPt_;   //!
  TBranch        *b_genParEta_;   //!
  TBranch        *b_genParPhi_;   //!
  TBranch        *b_genParQ_;   //!
  TBranch        *b_genParId_;   //!
  TBranch        *b_genParSt_;   //!
  TBranch        *b_genMomParId_;   //!
  TBranch        *b_genParIndex_;   //!
  TBranch        *b_genLepE_;   //!
  TBranch        *b_genLepPt_;   //!
  TBranch        *b_genLepEta_;   //!
  TBranch        *b_genLepPhi_;   //!
  TBranch        *b_genLepPhoE_;   //!
  TBranch        *b_genLepPhoPt_;   //!
  TBranch        *b_genLepPhoEta_;   //!
  TBranch        *b_genLepPhoPhi_;   //!
  TBranch        *b_genLepId_;   //!
  TBranch        *b_genJetE_;   //!
  TBranch        *b_genJetPt_;   //!
  TBranch        *b_genJetEta_;   //!
  TBranch        *b_genJetPhi_;   //!

  myLHAPDF_reweighing(std::string filename, TTree *tree=0);
  std::string _inputFileName;
  virtual ~myLHAPDF_reweighing();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop(int LepID);
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);

};


class MyPDF {
 public:
  int nIndex;
  MyPDF(const char* pdfname, int n) {
    nIndex = n;	
    LHAPDF::initPDFSet(nIndex,pdfname);
    LHAPDF::initPDFSet(3,"cteq6ll.LHpdf");
  }
  virtual ~MyPDF(){};

  double weight(vector<double>* Info_, unsigned int member)
  {
    int id1 = Info_->at(0); 
    int id2 = Info_->at(1); 
    double x1 = Info_->at(2); 
    double x2 = Info_->at(3); 
    double Q = Info_->at(6); 

    LHAPDF::usePDFMember(nIndex,0);
    double pdf1 = LHAPDF::xfx(nIndex,x1,Q,id1)*LHAPDF::xfx(nIndex,x2,Q,id2);

    LHAPDF::usePDFMember(3,member);
    double pdf2 = LHAPDF::xfx(3,x1,Q,id1)*LHAPDF::xfx(3,x2,Q,id2);

    if (pdf2>0) {
      return pdf1/pdf2;
    } else {
      printf ("pdf1 = %e, pdf2 = %e\n", pdf1, pdf2);
      return 1.;
    }

  }
}; // end of myPDF definitions


#endif

#ifdef myLHAPDF_reweighing_cxx
myLHAPDF_reweighing::myLHAPDF_reweighing(std::string filename, TTree *tree)
{
  _inputFileName = filename;
  cout << "Input file = " << _inputFileName << endl;
  TChain* chain = new TChain("tree/tree");
  chain->Add(_inputFileName.data());
  //  for(int i=1;i<=4;i++)      
  //    chain->Add(Form("/data2/syu/zjet_vectorNtuple/PDF/PDF_DYJets_madgraph_%d.root",i));  
  /*
  TSystemDirectory *base = new TSystemDirectory("root","root");
  string dirName = "root://eoscms//eos/cms/store/user/syu/PDF_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola";
  //    string dirName = "root://eoscms//eos/cms/store/user/syu/PDF_GJets_TuneZ2_40_HT_100_7TeV-madgraph";
  base->SetDirectory(dirName.data());
  TList *listOfFiles = base->GetListOfFiles();
  TIter fileIt(listOfFiles);
  TFile *fileH = new TFile();
  int nfile=0;
  while(fileH = (TFile*)fileIt()) {
    std::string fileN = fileH->GetName();
    std::string baseString = "zjets_mc";
    if( fileH->IsFolder())  continue;
    if(fileN.find(baseString) == std::string::npos)continue;
    cout << fileN.data() << endl;
    nfile++;
    chain->Add(Form("%s/%s",dirName.data(),fileN.data()));
  }
  
  cout << "Opening " << nfile << " files " << endl;
  */

  tree = chain;
  Init(tree);
  cout << "Number of entries are " << tree->GetEntries() << endl;
}

myLHAPDF_reweighing::~myLHAPDF_reweighing()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t myLHAPDF_reweighing::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t myLHAPDF_reweighing::LoadTree(Long64_t entry)
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

void myLHAPDF_reweighing::Init(TTree *tree)
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
  pdfInfo_ = 0;
  cteq66PDFw_ = 0;
  cteq6l1PDFw_ = 0;
  mstw2008PDFw_ = 0;
  genParE_ = 0;
  genParPt_ = 0;
  genParEta_ = 0;
  genParPhi_ = 0;
  genParQ_ = 0;
  genParId_ = 0;
  genParSt_ = 0;
  genMomParId_ = 0;
  genParIndex_ = 0;
  genLepE_ = 0;
  genLepPt_ = 0;
  genLepEta_ = 0;
  genLepPhi_ = 0;
  genLepPhoE_ = 0;
  genLepPhoPt_ = 0;
  genLepPhoEta_ = 0;
  genLepPhoPhi_ = 0;
  genLepId_ = 0;
  genJetE_ = 0;
  genJetPt_ = 0;
  genJetEta_ = 0;
  genJetPhi_ = 0;
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("PU_weight", &PU_weight, &b_PU_weight);
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
  fChain->SetBranchAddress("pdfInfo_", &pdfInfo_, &b_pdfInfo_);
  fChain->SetBranchAddress("cteq66PDFw_", &cteq66PDFw_, &b_cteq66PDFw_);
  fChain->SetBranchAddress("cteq6l1PDFw_", &cteq6l1PDFw_, &b_cteq6l1PDFw_);
  fChain->SetBranchAddress("mstw2008PDFw_", &mstw2008PDFw_, &b_mstw2008PDFw_);
  fChain->SetBranchAddress("genParE_", &genParE_, &b_genParE_);
  fChain->SetBranchAddress("genParPt_", &genParPt_, &b_genParPt_);
  fChain->SetBranchAddress("genParEta_", &genParEta_, &b_genParEta_);
  fChain->SetBranchAddress("genParPhi_", &genParPhi_, &b_genParPhi_);
  fChain->SetBranchAddress("genParQ_", &genParQ_, &b_genParQ_);
  fChain->SetBranchAddress("genParId_", &genParId_, &b_genParId_);
  fChain->SetBranchAddress("genParSt_", &genParSt_, &b_genParSt_);
  fChain->SetBranchAddress("genMomParId_", &genMomParId_, &b_genMomParId_);
  fChain->SetBranchAddress("genParIndex_", &genParIndex_, &b_genParIndex_);
  fChain->SetBranchAddress("genLepE_", &genLepE_, &b_genLepE_);
  fChain->SetBranchAddress("genLepPt_", &genLepPt_, &b_genLepPt_);
  fChain->SetBranchAddress("genLepEta_", &genLepEta_, &b_genLepEta_);
  fChain->SetBranchAddress("genLepPhi_", &genLepPhi_, &b_genLepPhi_);
  fChain->SetBranchAddress("genLepPhoE_", &genLepPhoE_, &b_genLepPhoE_);
  fChain->SetBranchAddress("genLepPhoPt_", &genLepPhoPt_, &b_genLepPhoPt_);
  fChain->SetBranchAddress("genLepPhoEta_", &genLepPhoEta_, &b_genLepPhoEta_);
  fChain->SetBranchAddress("genLepPhoPhi_", &genLepPhoPhi_, &b_genLepPhoPhi_);
  fChain->SetBranchAddress("genLepId_", &genLepId_, &b_genLepId_);
  fChain->SetBranchAddress("genJetE_", &genJetE_, &b_genJetE_);
  fChain->SetBranchAddress("genJetPt_", &genJetPt_, &b_genJetPt_);
  fChain->SetBranchAddress("genJetEta_", &genJetEta_, &b_genJetEta_);
  fChain->SetBranchAddress("genJetPhi_", &genJetPhi_, &b_genJetPhi_);
  Notify();
}

Bool_t myLHAPDF_reweighing::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void myLHAPDF_reweighing::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t myLHAPDF_reweighing::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef myLHAPDF_reweighing_cxx
