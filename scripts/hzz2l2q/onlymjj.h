//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Nov 28 15:14:38 2012 by ROOT version 5.32/00
// from TTree tree/tree
// found on file: /data4/syu/new_aod/GluGluToHToZZTo2L2Q_M-300_8TeV-powheg-pythia6.root
//////////////////////////////////////////////////////////

#ifndef onlymjj_h
#define onlymjj_h

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
const Int_t kMaxgenNMo = 1;
const Int_t kMaxgenNDa = 1;
const Int_t kMaxgenMo1 = 1;
const Int_t kMaxgenMo2 = 1;
const Int_t kMaxgenDa1 = 1;
const Int_t kMaxgenDa2 = 1;
const Int_t kMaxgenJetE = 1;
const Int_t kMaxgenJetPt = 1;
const Int_t kMaxgenJetEta = 1;
const Int_t kMaxgenJetPhi = 1;
const Int_t kMaxpatMuonPt = 1;
const Int_t kMaxpatMuonEta = 1;
const Int_t kMaxpatMuonPhi = 1;
const Int_t kMaxpatMuonM = 1;
const Int_t kMaxpatMuonTrkIso = 1;
const Int_t kMaxpatMuonHcalIso = 1;
const Int_t kMaxpatMuonEcalIso = 1;
const Int_t kMaxpatMuonCharge = 1;
const Int_t kMaxpatMuonNumChambers = 1;
const Int_t kMaxpatMuonNumMatches = 1;
const Int_t kMaxpatMuonStationMask = 1;
const Int_t kMaxpatMuonNumSegments = 1;
const Int_t kMaxpatMuonChi2Ndoff = 1;
const Int_t kMaxpatMuonNhits = 1;
const Int_t kMaxpatMuonDxy = 1;
const Int_t kMaxpatMuonDz = 1;
const Int_t kMaxpatMuonPassID = 1;
const Int_t kMaxpatMuonChHadIso = 1;
const Int_t kMaxpatMuonNeHadIso = 1;
const Int_t kMaxpatMuonGamIso = 1;
const Int_t kMaxpatMuonPUPt = 1;
const Int_t kMaxpatElecEt = 1;
const Int_t kMaxpatElecEnergy = 1;
const Int_t kMaxpatElecPt = 1;
const Int_t kMaxpatElecEta = 1;
const Int_t kMaxpatElecPhi = 1;
const Int_t kMaxpatElecM = 1;
const Int_t kMaxpatElecScEta = 1;
const Int_t kMaxpatElecSigIhIh = 1;
const Int_t kMaxpatElecDelEtaIn = 1;
const Int_t kMaxpatElecDelPhiIn = 1;
const Int_t kMaxpatElecHoE = 1;
const Int_t kMaxpatElecTrkIso = 1;
const Int_t kMaxpatElecHcalIso = 1;
const Int_t kMaxpatElecEcalIso = 1;
const Int_t kMaxpatElecEoverP = 1;
const Int_t kMaxpatElecDxy = 1;
const Int_t kMaxpatElecDz = 1;
const Int_t kMaxpatElecChHadIso = 1;
const Int_t kMaxpatElecNeHadIso = 1;
const Int_t kMaxpatElecGamIso = 1;
const Int_t kMaxpatElecMissingHits = 1;
const Int_t kMaxpatElecDist = 1;
const Int_t kMaxpatElecDeltaCotTheta = 1;
const Int_t kMaxpatElecInBarrel = 1;
const Int_t kMaxpatElecInEndcap = 1;
const Int_t kMaxpatElecHasConv = 1;
const Int_t kMaxpatElecPassID = 1;
const Int_t kMaxpatJetPfAk05Pt = 1;
const Int_t kMaxpatJetPfAk05Eta = 1;
const Int_t kMaxpatJetPfAk05Phi = 1;
const Int_t kMaxpatJetPfAk05M = 1;
const Int_t kMaxpatJetPfAk05En = 1;
const Int_t kMaxpatJetPfAk05Beta = 1;
const Int_t kMaxpatJetPfAk05PassID = 1;
const Int_t kMaxpatJetPfAk05CharMulti = 1;
const Int_t kMaxpatJetPfAk05NeutEmEFr = 1;
const Int_t kMaxpatJetPfAk05CharHadEFr = 1;
const Int_t kMaxpatJetPfAk05NeutHadEFr = 1;
const Int_t kMaxpatJetPfAk05CharEmEFr = 1;
const Int_t kMaxjetHiggsIndex = 1;

class onlymjj {
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
  vector<int>     *genNMo_;
  vector<int>     *genNDa_;
  vector<int>     *genMo1_;
  vector<int>     *genMo2_;
  vector<int>     *genDa1_;
  vector<int>     *genDa2_;
  vector<double>  *genJetE_;
  vector<double>  *genJetPt_;
  vector<double>  *genJetEta_;
  vector<double>  *genJetPhi_;
  Double_t        NumMu;
  vector<double>  *patMuonPt_;
  vector<double>  *patMuonEta_;
  vector<double>  *patMuonPhi_;
  vector<double>  *patMuonM_;
  vector<double>  *patMuonTrkIso_;
  vector<double>  *patMuonHcalIso_;
  vector<double>  *patMuonEcalIso_;
  vector<double>  *patMuonCharge_;
  vector<double>  *patMuonNumChambers_;
  vector<double>  *patMuonNumMatches_;
  vector<double>  *patMuonStationMask_;
  vector<double>  *patMuonNumSegments_;
  vector<double>  *patMuonChi2Ndoff_;
  vector<double>  *patMuonNhits_;
  vector<double>  *patMuonDxy_;
  vector<double>  *patMuonDz_;
  vector<int>     *patMuonPassID_;
  vector<double>  *patMuonChHadIso_;
  vector<double>  *patMuonNeHadIso_;
  vector<double>  *patMuonGamIso_;
  vector<double>  *patMuonPUPt_;
  vector<double>  *patElecEt_;
  vector<double>  *patElecEnergy_;
  vector<double>  *patElecPt_;
  vector<double>  *patElecEta_;
  vector<double>  *patElecPhi_;
  vector<double>  *patElecM_;
  vector<double>  *patElecScEta_;
  vector<double>  *patElecSigIhIh_;
  vector<double>  *patElecDelEtaIn_;
  vector<double>  *patElecDelPhiIn_;
  vector<double>  *patElecHoE_;
  vector<double>  *patElecTrkIso_;
  vector<double>  *patElecHcalIso_;
  vector<double>  *patElecEcalIso_;
  vector<double>  *patElecEoverP_;
  vector<double>  *patElecDxy_;
  vector<double>  *patElecDz_;
  vector<double>  *patElecChHadIso_;
  vector<double>  *patElecNeHadIso_;
  vector<double>  *patElecGamIso_;
  vector<double>  *patElecMissingHits_;
  vector<double>  *patElecDist_;
  vector<double>  *patElecDeltaCotTheta_;
  vector<double>  *patElecInBarrel_;
  vector<double>  *patElecInEndcap_;
  vector<int>     *patElecHasConv_;
  vector<int>     *patElecPassID_;
  vector<double>  *patJetPfAk05Pt_;
  vector<double>  *patJetPfAk05Eta_;
  vector<double>  *patJetPfAk05Phi_;
  vector<double>  *patJetPfAk05M_;
  vector<double>  *patJetPfAk05En_;
  vector<double>  *patJetPfAk05Beta_;
  vector<int>     *patJetPfAk05PassID_;
  vector<double>  *patJetPfAk05CharMulti_;
  vector<double>  *patJetPfAk05NeutEmEFr_;
  vector<double>  *patJetPfAk05CharHadEFr_;
  vector<double>  *patJetPfAk05NeutHadEFr_;
  vector<double>  *patJetPfAk05CharEmEFr_;
  Double_t        eleRho;
  Double_t        muoRho;
  Double_t        metSig;
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
  vector<int>     *jetHiggsIndex_;
  vector<double>  *jetE;
  vector<double>  *jetPt;
  vector<double>  *jetEta;
  vector<double>  *jetPhi;
  vector<double>  *jetEtInAnnulus0105;
  vector<double>  *jetMaxDist;
  vector<int>     *jetNConst;
  vector<double>  *jetConstPt;
  vector<double>  *jetConstEtaPhiSpread;
  vector<double>  *jetArea;
  vector<double>  *jetPileup;
  vector<double>  *jetEtaEtaMoment;
  vector<double>  *jetPhiPhiMoment;
  vector<double>  *jetEtaPhiMoment;
  vector<double>  *heliLD;
  vector<double>  *costhetaNT1;
  vector<double>  *costhetaNT2;
  vector<double>  *phiNT;
  vector<double>  *phiNT1;
  vector<double>  *costhetastarNT;
  vector<double>  *heliLDRefit;
  vector<double>  *costhetaNT1Refit;
  vector<double>  *costhetaNT2Refit;
  vector<double>  *phiNTRefit;
  vector<double>  *phiNT1Refit;
  vector<double>  *costhetastarNTRefit;
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
  TBranch        *b_genNMo_;   //!
  TBranch        *b_genNDa_;   //!
  TBranch        *b_genMo1_;   //!
  TBranch        *b_genMo2_;   //!
  TBranch        *b_genDa1_;   //!
  TBranch        *b_genDa2_;   //!
  TBranch        *b_genJetE_;   //!
  TBranch        *b_genJetPt_;   //!
  TBranch        *b_genJetEta_;   //!
  TBranch        *b_genJetPhi_;   //!
  TBranch        *b_ElectronNumMu;   //!
  TBranch        *b_patMuonPt_;   //!
  TBranch        *b_patMuonEta_;   //!
  TBranch        *b_patMuonPhi_;   //!
  TBranch        *b_patMuonM_;   //!
  TBranch        *b_patMuonTrkIso_;   //!
  TBranch        *b_patMuonHcalIso_;   //!
  TBranch        *b_patMuonEcalIso_;   //!
  TBranch        *b_patMuonCharge_;   //!
  TBranch        *b_patMuonNumChambers_;   //!
  TBranch        *b_patMuonNumMatches_;   //!
  TBranch        *b_patMuonStationMask_;   //!
  TBranch        *b_patMuonNumSegments_;   //!
  TBranch        *b_patMuonChi2Ndoff_;   //!
  TBranch        *b_patMuonNhits_;   //!
  TBranch        *b_patMuonDxy_;   //!
  TBranch        *b_patMuonDz_;   //!
  TBranch        *b_patMuonPassID_;   //!
  TBranch        *b_patMuonChHadIso_;   //!
  TBranch        *b_patMuonNeHadIso_;   //!
  TBranch        *b_patMuonGamIso_;   //!
  TBranch        *b_patMuonPUPt_;   //!
  TBranch        *b_patElecEt_;   //!
  TBranch        *b_patElecEnergy_;   //!
  TBranch        *b_patElecPt_;   //!
  TBranch        *b_patElecEta_;   //!
  TBranch        *b_patElecPhi_;   //!
  TBranch        *b_patElecM_;   //!
  TBranch        *b_patElecScEta_;   //!
  TBranch        *b_patElecSigIhIh_;   //!
  TBranch        *b_patElecDelEtaIn_;   //!
  TBranch        *b_patElecDelPhiIn_;   //!
  TBranch        *b_patElecHoE_;   //!
  TBranch        *b_patElecTrkIso_;   //!
  TBranch        *b_patElecHcalIso_;   //!
  TBranch        *b_patElecEcalIso_;   //!
  TBranch        *b_patElecEoverP_;   //!
  TBranch        *b_patElecDxy_;   //!
  TBranch        *b_patElecDz_;   //!
  TBranch        *b_patElecChHadIso_;   //!
  TBranch        *b_patElecNeHadIso_;   //!
  TBranch        *b_patElecGamIso_;   //!
  TBranch        *b_patElecMissingHits_;   //!
  TBranch        *b_patElecDist_;   //!
  TBranch        *b_patElecDeltaCotTheta_;   //!
  TBranch        *b_patElecInBarrel_;   //!
  TBranch        *b_patElecInEndcap_;   //!
  TBranch        *b_patElecHasConv_;   //!
  TBranch        *b_patElecPassID_;   //!
  TBranch        *b_patJetPfAk05Pt_;   //!
  TBranch        *b_patJetPfAk05Eta_;   //!
  TBranch        *b_patJetPfAk05Phi_;   //!
  TBranch        *b_patJetPfAk05M_;   //!
  TBranch        *b_patJetPfAk05En_;   //!
  TBranch        *b_patJetPfAk05Beta_;   //!
  TBranch        *b_patJetPfAk05PassID_;   //!
  TBranch        *b_patJetPfAk05CharMulti_;   //!
  TBranch        *b_patJetPfAk05NeutEmEFr_;   //!
  TBranch        *b_patJetPfAk05CharHadEFr_;   //!
  TBranch        *b_patJetPfAk05NeutHadEFr_;   //!
  TBranch        *b_patJetPfAk05CharEmEFr_;   //!
  TBranch        *b_eleRho;   //!
  TBranch        *b_muoRho;   //!
  TBranch        *b_metSig;   //!
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
  TBranch        *b_jetHiggsIndex_;   //!
  TBranch        *b_jetE;   //!
  TBranch        *b_jetPt;   //!
  TBranch        *b_jetEta;   //!
  TBranch        *b_jetPhi;   //!
  TBranch        *b_jetEtInAnnulus0105;   //!
  TBranch        *b_jetMaxDist;   //!
  TBranch        *b_jetNConst;   //!
  TBranch        *b_jetConstPt;   //!
  TBranch        *b_jetConstEtaPhiSpread;   //!
  TBranch        *b_jetArea;   //!
  TBranch        *b_jetPileup;   //!
  TBranch        *b_jetEtaEtaMoment;   //!
  TBranch        *b_jetPhiPhiMoment;   //!
  TBranch        *b_jetEtaPhiMoment;   //!
  TBranch        *b_heliLD;   //!
  TBranch        *b_costhetaNT1;   //!
  TBranch        *b_costhetaNT2;   //!
  TBranch        *b_phiNT;   //!
  TBranch        *b_phiNT1;   //!
  TBranch        *b_costhetastarNT;   //!
  TBranch        *b_heliLDRefit;   //!
  TBranch        *b_costhetaNT1Refit;   //!
  TBranch        *b_costhetaNT2Refit;   //!
  TBranch        *b_phiNTRefit;   //!
  TBranch        *b_phiNT1Refit;   //!
  TBranch        *b_costhetastarNTRefit;   //!
  TBranch        *b_nBTags;   //!
  TBranch        *b_lepType;   //!
  TBranch        *b_passBit;   //!

  onlymjj(int mass, TTree *tree=0);
  virtual ~onlymjj();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop(int DEBUG=0, 
			bool reweighSignal=true,
			bool reweighPU=true);
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);

  virtual Bool_t   matchGenJetToParton(Int_t igen, Int_t ijet);
  virtual Bool_t   matchPFJetToParton(Int_t igen, Int_t ijet);
  int     _higgsMass;

};

#endif

#ifdef onlymjj_cxx
onlymjj::onlymjj(int mass, TTree *tree) : fChain(0)
{
  string filename = Form("/data4/syu/new_aod/GluGluToHToZZTo2L2Q_M-%d_8TeV-powheg-pythia6.root",mass);
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

  Init(tree);
  _higgsMass = mass;
}

onlymjj::~onlymjj()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t onlymjj::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t onlymjj::LoadTree(Long64_t entry)
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

void onlymjj::Init(TTree *tree)
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
  genNMo_ = 0;
  genNDa_ = 0;
  genMo1_ = 0;
  genMo2_ = 0;
  genDa1_ = 0;
  genDa2_ = 0;
  genJetE_ = 0;
  genJetPt_ = 0;
  genJetEta_ = 0;
  genJetPhi_ = 0;
  patMuonPt_ = 0;
  patMuonEta_ = 0;
  patMuonPhi_ = 0;
  patMuonM_ = 0;
  patMuonTrkIso_ = 0;
  patMuonHcalIso_ = 0;
  patMuonEcalIso_ = 0;
  patMuonCharge_ = 0;
  patMuonNumChambers_ = 0;
  patMuonNumMatches_ = 0;
  patMuonStationMask_ = 0;
  patMuonNumSegments_ = 0;
  patMuonChi2Ndoff_ = 0;
  patMuonNhits_ = 0;
  patMuonDxy_ = 0;
  patMuonDz_ = 0;
  patMuonPassID_ = 0;
  patMuonChHadIso_ = 0;
  patMuonNeHadIso_ = 0;
  patMuonGamIso_ = 0;
  patMuonPUPt_ = 0;
  patElecEt_ = 0;
  patElecEnergy_ = 0;
  patElecPt_ = 0;
  patElecEta_ = 0;
  patElecPhi_ = 0;
  patElecM_ = 0;
  patElecScEta_ = 0;
  patElecSigIhIh_ = 0;
  patElecDelEtaIn_ = 0;
  patElecDelPhiIn_ = 0;
  patElecHoE_ = 0;
  patElecTrkIso_ = 0;
  patElecHcalIso_ = 0;
  patElecEcalIso_ = 0;
  patElecEoverP_ = 0;
  patElecDxy_ = 0;
  patElecDz_ = 0;
  patElecChHadIso_ = 0;
  patElecNeHadIso_ = 0;
  patElecGamIso_ = 0;
  patElecMissingHits_ = 0;
  patElecDist_ = 0;
  patElecDeltaCotTheta_ = 0;
  patElecInBarrel_ = 0;
  patElecInEndcap_ = 0;
  patElecHasConv_ = 0;
  patElecPassID_ = 0;
  patJetPfAk05Pt_ = 0;
  patJetPfAk05Eta_ = 0;
  patJetPfAk05Phi_ = 0;
  patJetPfAk05M_ = 0;
  patJetPfAk05En_ = 0;
  patJetPfAk05Beta_ = 0;
  patJetPfAk05PassID_ = 0;
  patJetPfAk05CharMulti_ = 0;
  patJetPfAk05NeutEmEFr_ = 0;
  patJetPfAk05CharHadEFr_ = 0;
  patJetPfAk05NeutHadEFr_ = 0;
  patJetPfAk05CharEmEFr_ = 0;
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
  jetHiggsIndex_ = 0;
  jetE = 0;
  jetPt = 0;
  jetEta = 0;
  jetPhi = 0;
  jetEtInAnnulus0105 = 0;
  jetMaxDist = 0;
  jetNConst = 0;
  jetConstPt = 0;
  jetConstEtaPhiSpread = 0;
  jetArea = 0;
  jetPileup = 0;
  jetEtaEtaMoment = 0;
  jetPhiPhiMoment = 0;
  jetEtaPhiMoment = 0;
  heliLD = 0;
  costhetaNT1 = 0;
  costhetaNT2 = 0;
  phiNT = 0;
  phiNT1 = 0;
  costhetastarNT = 0;
  heliLDRefit = 0;
  costhetaNT1Refit = 0;
  costhetaNT2Refit = 0;
  phiNTRefit = 0;
  phiNT1Refit = 0;
  costhetastarNTRefit = 0;
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
  fChain->SetBranchAddress("genNMo_", &genNMo_, &b_genNMo_);
  fChain->SetBranchAddress("genNDa_", &genNDa_, &b_genNDa_);
  fChain->SetBranchAddress("genMo1_", &genMo1_, &b_genMo1_);
  fChain->SetBranchAddress("genMo2_", &genMo2_, &b_genMo2_);
  fChain->SetBranchAddress("genDa1_", &genDa1_, &b_genDa1_);
  fChain->SetBranchAddress("genDa2_", &genDa2_, &b_genDa2_);
  fChain->SetBranchAddress("genJetE_", &genJetE_, &b_genJetE_);
  fChain->SetBranchAddress("genJetPt_", &genJetPt_, &b_genJetPt_);
  fChain->SetBranchAddress("genJetEta_", &genJetEta_, &b_genJetEta_);
  fChain->SetBranchAddress("genJetPhi_", &genJetPhi_, &b_genJetPhi_);
  /*
  fChain->SetBranchAddress("NumMu", &NumMu, &b_ElectronNumMu);
  fChain->SetBranchAddress("patMuonPt_", &patMuonPt_, &b_patMuonPt_);
  fChain->SetBranchAddress("patMuonEta_", &patMuonEta_, &b_patMuonEta_);
  fChain->SetBranchAddress("patMuonPhi_", &patMuonPhi_, &b_patMuonPhi_);
  fChain->SetBranchAddress("patMuonM_", &patMuonM_, &b_patMuonM_);
  fChain->SetBranchAddress("patMuonTrkIso_", &patMuonTrkIso_, &b_patMuonTrkIso_);
  fChain->SetBranchAddress("patMuonHcalIso_", &patMuonHcalIso_, &b_patMuonHcalIso_);
  fChain->SetBranchAddress("patMuonEcalIso_", &patMuonEcalIso_, &b_patMuonEcalIso_);
  fChain->SetBranchAddress("patMuonCharge_", &patMuonCharge_, &b_patMuonCharge_);
  fChain->SetBranchAddress("patMuonNumChambers_", &patMuonNumChambers_, &b_patMuonNumChambers_);
  fChain->SetBranchAddress("patMuonNumMatches_", &patMuonNumMatches_, &b_patMuonNumMatches_);
  fChain->SetBranchAddress("patMuonStationMask_", &patMuonStationMask_, &b_patMuonStationMask_);
  fChain->SetBranchAddress("patMuonNumSegments_", &patMuonNumSegments_, &b_patMuonNumSegments_);
  fChain->SetBranchAddress("patMuonChi2Ndoff_", &patMuonChi2Ndoff_, &b_patMuonChi2Ndoff_);
  fChain->SetBranchAddress("patMuonNhits_", &patMuonNhits_, &b_patMuonNhits_);
  fChain->SetBranchAddress("patMuonDxy_", &patMuonDxy_, &b_patMuonDxy_);
  fChain->SetBranchAddress("patMuonDz_", &patMuonDz_, &b_patMuonDz_);
  fChain->SetBranchAddress("patMuonPassID_", &patMuonPassID_, &b_patMuonPassID_);
  fChain->SetBranchAddress("patMuonChHadIso_", &patMuonChHadIso_, &b_patMuonChHadIso_);
  fChain->SetBranchAddress("patMuonNeHadIso_", &patMuonNeHadIso_, &b_patMuonNeHadIso_);
  fChain->SetBranchAddress("patMuonGamIso_", &patMuonGamIso_, &b_patMuonGamIso_);
  fChain->SetBranchAddress("patMuonPUPt_", &patMuonPUPt_, &b_patMuonPUPt_);
  fChain->SetBranchAddress("patElecEt_", &patElecEt_, &b_patElecEt_);
  fChain->SetBranchAddress("patElecEnergy_", &patElecEnergy_, &b_patElecEnergy_);
  fChain->SetBranchAddress("patElecPt_", &patElecPt_, &b_patElecPt_);
  fChain->SetBranchAddress("patElecEta_", &patElecEta_, &b_patElecEta_);
  fChain->SetBranchAddress("patElecPhi_", &patElecPhi_, &b_patElecPhi_);
  fChain->SetBranchAddress("patElecM_", &patElecM_, &b_patElecM_);
  fChain->SetBranchAddress("patElecScEta_", &patElecScEta_, &b_patElecScEta_);
  fChain->SetBranchAddress("patElecSigIhIh_", &patElecSigIhIh_, &b_patElecSigIhIh_);
  fChain->SetBranchAddress("patElecDelEtaIn_", &patElecDelEtaIn_, &b_patElecDelEtaIn_);
  fChain->SetBranchAddress("patElecDelPhiIn_", &patElecDelPhiIn_, &b_patElecDelPhiIn_);
  fChain->SetBranchAddress("patElecHoE_", &patElecHoE_, &b_patElecHoE_);
  fChain->SetBranchAddress("patElecTrkIso_", &patElecTrkIso_, &b_patElecTrkIso_);
  fChain->SetBranchAddress("patElecHcalIso_", &patElecHcalIso_, &b_patElecHcalIso_);
  fChain->SetBranchAddress("patElecEcalIso_", &patElecEcalIso_, &b_patElecEcalIso_);
  fChain->SetBranchAddress("patElecEoverP_", &patElecEoverP_, &b_patElecEoverP_);
  fChain->SetBranchAddress("patElecDxy_", &patElecDxy_, &b_patElecDxy_);
  fChain->SetBranchAddress("patElecDz_", &patElecDz_, &b_patElecDz_);
  fChain->SetBranchAddress("patElecChHadIso_", &patElecChHadIso_, &b_patElecChHadIso_);
  fChain->SetBranchAddress("patElecNeHadIso_", &patElecNeHadIso_, &b_patElecNeHadIso_);
  fChain->SetBranchAddress("patElecGamIso_", &patElecGamIso_, &b_patElecGamIso_);
  fChain->SetBranchAddress("patElecMissingHits_", &patElecMissingHits_, &b_patElecMissingHits_);
  fChain->SetBranchAddress("patElecDist_", &patElecDist_, &b_patElecDist_);
  fChain->SetBranchAddress("patElecDeltaCotTheta_", &patElecDeltaCotTheta_, &b_patElecDeltaCotTheta_);
  fChain->SetBranchAddress("patElecInBarrel_", &patElecInBarrel_, &b_patElecInBarrel_);
  fChain->SetBranchAddress("patElecInEndcap_", &patElecInEndcap_, &b_patElecInEndcap_);
  fChain->SetBranchAddress("patElecHasConv_", &patElecHasConv_, &b_patElecHasConv_);
  fChain->SetBranchAddress("patElecPassID_", &patElecPassID_, &b_patElecPassID_);
  fChain->SetBranchAddress("patJetPfAk05Pt_", &patJetPfAk05Pt_, &b_patJetPfAk05Pt_);
  fChain->SetBranchAddress("patJetPfAk05Eta_", &patJetPfAk05Eta_, &b_patJetPfAk05Eta_);
  fChain->SetBranchAddress("patJetPfAk05Phi_", &patJetPfAk05Phi_, &b_patJetPfAk05Phi_);
  fChain->SetBranchAddress("patJetPfAk05M_", &patJetPfAk05M_, &b_patJetPfAk05M_);
  fChain->SetBranchAddress("patJetPfAk05En_", &patJetPfAk05En_, &b_patJetPfAk05En_);
  fChain->SetBranchAddress("patJetPfAk05Beta_", &patJetPfAk05Beta_, &b_patJetPfAk05Beta_);
  fChain->SetBranchAddress("patJetPfAk05PassID_", &patJetPfAk05PassID_, &b_patJetPfAk05PassID_);
  fChain->SetBranchAddress("patJetPfAk05CharMulti_", &patJetPfAk05CharMulti_, &b_patJetPfAk05CharMulti_);
  fChain->SetBranchAddress("patJetPfAk05NeutEmEFr_", &patJetPfAk05NeutEmEFr_, &b_patJetPfAk05NeutEmEFr_);
  fChain->SetBranchAddress("patJetPfAk05CharHadEFr_", &patJetPfAk05CharHadEFr_, &b_patJetPfAk05CharHadEFr_);
  fChain->SetBranchAddress("patJetPfAk05NeutHadEFr_", &patJetPfAk05NeutHadEFr_, &b_patJetPfAk05NeutHadEFr_);
  fChain->SetBranchAddress("patJetPfAk05CharEmEFr_", &patJetPfAk05CharEmEFr_, &b_patJetPfAk05CharEmEFr_);
  */
  fChain->SetBranchAddress("eleRho", &eleRho, &b_eleRho);
  fChain->SetBranchAddress("muoRho", &muoRho, &b_muoRho);
  fChain->SetBranchAddress("metSig", &metSig, &b_metSig);
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
  fChain->SetBranchAddress("jetHiggsIndex_", &jetHiggsIndex_, &b_jetHiggsIndex_);
  fChain->SetBranchAddress("jetE", &jetE, &b_jetE);
  fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
  fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
  fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
  fChain->SetBranchAddress("jetEtInAnnulus0105", &jetEtInAnnulus0105, &b_jetEtInAnnulus0105);
  fChain->SetBranchAddress("jetMaxDist", &jetMaxDist, &b_jetMaxDist);
  fChain->SetBranchAddress("jetNConst", &jetNConst, &b_jetNConst);
  fChain->SetBranchAddress("jetConstPt", &jetConstPt, &b_jetConstPt);
  fChain->SetBranchAddress("jetConstEtaPhiSpread", &jetConstEtaPhiSpread, &b_jetConstEtaPhiSpread);
  fChain->SetBranchAddress("jetArea", &jetArea, &b_jetArea);
  fChain->SetBranchAddress("jetPileup", &jetPileup, &b_jetPileup);
  fChain->SetBranchAddress("jetEtaEtaMoment", &jetEtaEtaMoment, &b_jetEtaEtaMoment);
  fChain->SetBranchAddress("jetPhiPhiMoment", &jetPhiPhiMoment, &b_jetPhiPhiMoment);
  fChain->SetBranchAddress("jetEtaPhiMoment", &jetEtaPhiMoment, &b_jetEtaPhiMoment);
  fChain->SetBranchAddress("heliLD", &heliLD, &b_heliLD);
  fChain->SetBranchAddress("costhetaNT1", &costhetaNT1, &b_costhetaNT1);
  fChain->SetBranchAddress("costhetaNT2", &costhetaNT2, &b_costhetaNT2);
  fChain->SetBranchAddress("phiNT", &phiNT, &b_phiNT);
  fChain->SetBranchAddress("phiNT1", &phiNT1, &b_phiNT1);
  fChain->SetBranchAddress("costhetastarNT", &costhetastarNT, &b_costhetastarNT);
  fChain->SetBranchAddress("heliLDRefit", &heliLDRefit, &b_heliLDRefit);
  fChain->SetBranchAddress("costhetaNT1Refit", &costhetaNT1Refit, &b_costhetaNT1Refit);
  fChain->SetBranchAddress("costhetaNT2Refit", &costhetaNT2Refit, &b_costhetaNT2Refit);
  fChain->SetBranchAddress("phiNTRefit", &phiNTRefit, &b_phiNTRefit);
  fChain->SetBranchAddress("phiNT1Refit", &phiNT1Refit, &b_phiNT1Refit);
  fChain->SetBranchAddress("costhetastarNTRefit", &costhetastarNTRefit, &b_costhetastarNTRefit);
  fChain->SetBranchAddress("nBTags", &nBTags, &b_nBTags);
  fChain->SetBranchAddress("lepType", &lepType, &b_lepType);
  fChain->SetBranchAddress("passBit", &passBit, &b_passBit);
  Notify();
}

Bool_t onlymjj::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void onlymjj::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t onlymjj::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef onlymjj_cxx
