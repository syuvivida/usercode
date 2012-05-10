//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May  1 15:27:18 2012 by ROOT version 5.27/06b
// from TTree tree/tree
// found on file: /data4/syu/Gjet_vectorNtuple/GJets_TuneZ2_40_HT_100_7TeV-madgraph.root
//////////////////////////////////////////////////////////

#ifndef vector_angular_h
#define vector_angular_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <string>
#include <TH2.h>

using namespace std;

   const Int_t kMaxEvtInfo_VertexX = 1;
   const Int_t kMaxEvtInfo_VertexY = 1;
   const Int_t kMaxEvtInfo_VertexZ = 1;
   const Int_t kMaxptHat = 1;
   const Int_t kMaxmcWeight = 1;
   const Int_t kMaxgenParE = 1;
   const Int_t kMaxgenParPt = 1;
   const Int_t kMaxgenParEta = 1;
   const Int_t kMaxgenParPhi = 1;
   const Int_t kMaxgenParQ = 1;
   const Int_t kMaxgenParId = 1;
   const Int_t kMaxgenParSt = 1;
   const Int_t kMaxgenMomParId = 1;
   const Int_t kMaxgenParIndex = 1;
   const Int_t kMaxgenJetE = 1;
   const Int_t kMaxgenJetPt = 1;
   const Int_t kMaxgenJetEta = 1;
   const Int_t kMaxgenJetPhi = 1;
   const Int_t kMaxpatJetPfAk05Pt = 1;
   const Int_t kMaxpatJetPfAk05Eta = 1;
   const Int_t kMaxpatJetPfAk05Phi = 1;
   const Int_t kMaxpatJetPfAk05M = 1;
   const Int_t kMaxpatJetPfAk05Rapidity = 1;
   const Int_t kMaxpatJetPfAk05Px = 1;
   const Int_t kMaxpatJetPfAk05Py = 1;
   const Int_t kMaxpatJetPfAk05Pz = 1;
   const Int_t kMaxpatJetPfAk05En = 1;
   const Int_t kMaxpatJetPfAk05UnCorrPt = 1;
   const Int_t kMaxpatJetPfAk05UnCorrPx = 1;
   const Int_t kMaxpatJetPfAk05UnCorrPy = 1;
   const Int_t kMaxpatJetPfAk05UnCorrPz = 1;
   const Int_t kMaxpatJetPfAk05UnCorrEnt = 1;
   const Int_t kMaxpatJetPfAk05PhotEn = 1;
   const Int_t kMaxpatJetPfAk05ElecEn = 1;
   const Int_t kMaxpatJetPfAk05MuonEn = 1;
   const Int_t kMaxpatJetPfAk05HfHadEn = 1;
   const Int_t kMaxpatJetPfAk05HfEmEn = 1;
   const Int_t kMaxpatJetPfAk05CharHadE = 1;
   const Int_t kMaxpatJetPfAk05NeutHadE = 1;
   const Int_t kMaxpatJetPfAk05CharEmE = 1;
   const Int_t kMaxpatJetPfAk05CharMuE = 1;
   const Int_t kMaxpatJetPfAk05NeutEmE = 1;
   const Int_t kMaxpatJetPfAk05MuonMulti = 1;
   const Int_t kMaxpatJetPfAk05NeutMulti = 1;
   const Int_t kMaxpatJetPfAk05CharMulti = 1;
   const Int_t kMaxpatJetPfAk05CharHadMulti = 1;
   const Int_t kMaxpatJetPfAk05NeutHadMulti = 1;
   const Int_t kMaxpatJetPfAk05PhotMulti = 1;
   const Int_t kMaxpatJetPfAk05ElecMulti = 1;
   const Int_t kMaxpatJetPfAk05HfHadMulti = 1;
   const Int_t kMaxpatJetPfAk05HfEmMulti = 1;
   const Int_t kMaxpatJetPfAk05PhotEnFr = 1;
   const Int_t kMaxpatJetPfAk05MuonEnFr = 1;
   const Int_t kMaxpatJetPfAk05HfHadEnFr = 1;
   const Int_t kMaxpatJetPfAk05HfEmEnFr = 1;
   const Int_t kMaxpatJetPfAk05NeutEmEFr = 1;
   const Int_t kMaxpatJetPfAk05CharHadEFr = 1;
   const Int_t kMaxpatJetPfAk05NeutHadEFr = 1;
   const Int_t kMaxpatJetPfAk05CharEmEFr = 1;
   const Int_t kMaxpatJetPfAk05CharMuEFr = 1;
   const Int_t kMaxpatJetPfAk05N60 = 1;
   const Int_t kMaxpatJetPfAk05N90 = 1;
   const Int_t kMaxpatJetPfAk05EmFr = 1;
   const Int_t kMaxpatJetPfAk05HadFr = 1;
   const Int_t kMaxpatJetPfAk05EmEbEn = 1;
   const Int_t kMaxpatJetPfAk05EmEeEn = 1;
   const Int_t kMaxpatJetPfAk05EmHfEn = 1;
   const Int_t kMaxpatJetPfAk05HadHbEn = 1;
   const Int_t kMaxpatJetPfAk05HadHeEn = 1;
   const Int_t kMaxpatJetPfAk05HadHfEn = 1;
   const Int_t kMaxpatJetPfAk05HadHoEn = 1;
   const Int_t kMaxpatJetPfAk05TotalEm = 1;
   const Int_t kMaxpatJetPfAk05TotalHad = 1;
   const Int_t kMaxpatJetPfAk05EmEbFr = 1;
   const Int_t kMaxpatJetPfAk05EmEeFr = 1;
   const Int_t kMaxpatJetPfAk05EmHfFr = 1;
   const Int_t kMaxpatJetPfAk05HadHbFr = 1;
   const Int_t kMaxpatJetPfAk05HadHeFr = 1;
   const Int_t kMaxpatJetPfAk05HadHfFr = 1;
   const Int_t kMaxpatJetPfAk05HadHoFr = 1;
   const Int_t kMaxpatJetPfAk05JesUncert = 1;
   const Int_t kMaxpatPhotonNumPh = 1;

class vector_angular {
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
   vector<double>  *genParE_;
   vector<double>  *genParPt_;
   vector<double>  *genParEta_;
   vector<double>  *genParPhi_;
   vector<int>     *genParQ_;
   vector<int>     *genParId_;
   vector<int>     *genParSt_;
   vector<int>     *genMomParId_;
   vector<int>     *genParIndex_;
   vector<double>  *genJetE_;
   vector<double>  *genJetPt_;
   vector<double>  *genJetEta_;
   vector<double>  *genJetPhi_;
   vector<double>  *patJetPfAk05Pt_;
   vector<double>  *patJetPfAk05Eta_;
   vector<double>  *patJetPfAk05Phi_;
   vector<double>  *patJetPfAk05M_;
   vector<double>  *patJetPfAk05Rapidity_;
   vector<double>  *patJetPfAk05Px_;
   vector<double>  *patJetPfAk05Py_;
   vector<double>  *patJetPfAk05Pz_;
   vector<double>  *patJetPfAk05En_;
   vector<double>  *patJetPfAk05UnCorrPt_;
   vector<double>  *patJetPfAk05UnCorrPx_;
   vector<double>  *patJetPfAk05UnCorrPy_;
   vector<double>  *patJetPfAk05UnCorrPz_;
   vector<double>  *patJetPfAk05UnCorrEnt_;
   vector<double>  *patJetPfAk05PhotEn_;
   vector<double>  *patJetPfAk05ElecEn_;
   vector<double>  *patJetPfAk05MuonEn_;
   vector<double>  *patJetPfAk05HfHadEn_;
   vector<double>  *patJetPfAk05HfEmEn_;
   vector<double>  *patJetPfAk05CharHadE_;
   vector<double>  *patJetPfAk05NeutHadE_;
   vector<double>  *patJetPfAk05CharEmE_;
   vector<double>  *patJetPfAk05CharMuE_;
   vector<double>  *patJetPfAk05NeutEmE_;
   vector<double>  *patJetPfAk05MuonMulti_;
   vector<double>  *patJetPfAk05NeutMulti_;
   vector<double>  *patJetPfAk05CharMulti_;
   vector<double>  *patJetPfAk05CharHadMulti_;
   vector<double>  *patJetPfAk05NeutHadMulti_;
   vector<double>  *patJetPfAk05PhotMulti_;
   vector<double>  *patJetPfAk05ElecMulti_;
   vector<double>  *patJetPfAk05HfHadMulti_;
   vector<double>  *patJetPfAk05HfEmMulti_;
   vector<double>  *patJetPfAk05PhotEnFr_;
   vector<double>  *patJetPfAk05MuonEnFr_;
   vector<double>  *patJetPfAk05HfHadEnFr_;
   vector<double>  *patJetPfAk05HfEmEnFr_;
   vector<double>  *patJetPfAk05NeutEmEFr_;
   vector<double>  *patJetPfAk05CharHadEFr_;
   vector<double>  *patJetPfAk05NeutHadEFr_;
   vector<double>  *patJetPfAk05CharEmEFr_;
   vector<double>  *patJetPfAk05CharMuEFr_;
   vector<double>  *patJetPfAk05N60_;
   vector<double>  *patJetPfAk05N90_;
   vector<double>  *patJetPfAk05EmFr_;
   vector<double>  *patJetPfAk05HadFr_;
   vector<double>  *patJetPfAk05EmEbEn_;
   vector<double>  *patJetPfAk05EmEeEn_;
   vector<double>  *patJetPfAk05EmHfEn_;
   vector<double>  *patJetPfAk05HadHbEn_;
   vector<double>  *patJetPfAk05HadHeEn_;
   vector<double>  *patJetPfAk05HadHfEn_;
   vector<double>  *patJetPfAk05HadHoEn_;
   vector<double>  *patJetPfAk05TotalEm_;
   vector<double>  *patJetPfAk05TotalHad_;
   vector<double>  *patJetPfAk05EmEbFr_;
   vector<double>  *patJetPfAk05EmEeFr_;
   vector<double>  *patJetPfAk05EmHfFr_;
   vector<double>  *patJetPfAk05HadHbFr_;
   vector<double>  *patJetPfAk05HadHeFr_;
   vector<double>  *patJetPfAk05HadHfFr_;
   vector<double>  *patJetPfAk05HadHoFr_;
   vector<double>  *patJetPfAk05JesUncert_;
   Double_t        patPhotonrho25;
   Double_t        patPhotonNumPh_;
   vector<double>  *patPhotonPt;
   vector<double>  *patPhotonEta;
   vector<double>  *patPhotonPhi;
   vector<double>  *patPhotonEt;
   vector<double>  *patPhotonEnergy;
   vector<double>  *patPhotonPx;
   vector<double>  *patPhotonPy;
   vector<double>  *patPhotonPz;
   vector<double>  *patPhotonR9;
   vector<double>  *patPhotonPhiWidth;
   vector<double>  *patPhotonEtaWidth;
   vector<double>  *patPhotonScPhi;
   vector<double>  *patPhotonScEta;
   vector<double>  *patPhotonSigmaIetaIeta;
   vector<double>  *patPhotonhadronicOverEm;
   vector<double>  *patPhotonecalRecHitSumEtConeDR04;
   vector<double>  *patPhotonhcalTowerSumEtConeDR04;
   vector<double>  *patPhotontrkSumPtHollowConeDR04;
   vector<double>  *patPhotonhasPixelSeed;
   Double_t        patPhotonMCpthat;
   vector<int>     *patPhotonisGenMatched;
   vector<double>  *patPhotongenMomId;
   vector<double>  *patPhotongenGrandMomId;
   vector<double>  *patPhotongenNSiblings;
   vector<double>  *patPhotongenMatchedE;
   vector<double>  *patPhotongenMatchedPx;
   vector<double>  *patPhotongenMatchedPy;
   vector<double>  *patPhotongenMatchedPz;
   vector<double>  *patPhotongenMatchedPt;
   vector<double>  *patPhotongenMatchedEta;
   vector<double>  *patPhotongenMatchedPhi;
   vector<double>  *patPhotongenCalIsoDR04;
   vector<double>  *patPhotongenTrkIsoDR04;
   vector<double>  *patPhotongenIsoDR04;

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
   TBranch        *b_genParE_;   //!
   TBranch        *b_genParPt_;   //!
   TBranch        *b_genParEta_;   //!
   TBranch        *b_genParPhi_;   //!
   TBranch        *b_genParQ_;   //!
   TBranch        *b_genParId_;   //!
   TBranch        *b_genParSt_;   //!
   TBranch        *b_genMomParId_;   //!
   TBranch        *b_genParIndex_;   //!
   TBranch        *b_genJetE_;   //!
   TBranch        *b_genJetPt_;   //!
   TBranch        *b_genJetEta_;   //!
   TBranch        *b_genJetPhi_;   //!
   TBranch        *b_patJetPfAk05Pt_;   //!
   TBranch        *b_patJetPfAk05Eta_;   //!
   TBranch        *b_patJetPfAk05Phi_;   //!
   TBranch        *b_patJetPfAk05M_;   //!
   TBranch        *b_patJetPfAk05Rapidity_;   //!
   TBranch        *b_patJetPfAk05Px_;   //!
   TBranch        *b_patJetPfAk05Py_;   //!
   TBranch        *b_patJetPfAk05Pz_;   //!
   TBranch        *b_patJetPfAk05En_;   //!
   TBranch        *b_patJetPfAk05UnCorrPt_;   //!
   TBranch        *b_patJetPfAk05UnCorrPx_;   //!
   TBranch        *b_patJetPfAk05UnCorrPy_;   //!
   TBranch        *b_patJetPfAk05UnCorrPz_;   //!
   TBranch        *b_patJetPfAk05UnCorrEnt_;   //!
   TBranch        *b_patJetPfAk05PhotEn_;   //!
   TBranch        *b_patJetPfAk05ElecEn_;   //!
   TBranch        *b_patJetPfAk05MuonEn_;   //!
   TBranch        *b_patJetPfAk05HfHadEn_;   //!
   TBranch        *b_patJetPfAk05HfEmEn_;   //!
   TBranch        *b_patJetPfAk05CharHadE_;   //!
   TBranch        *b_patJetPfAk05NeutHadE_;   //!
   TBranch        *b_patJetPfAk05CharEmE_;   //!
   TBranch        *b_patJetPfAk05CharMuE_;   //!
   TBranch        *b_patJetPfAk05NeutEmE_;   //!
   TBranch        *b_patJetPfAk05MuonMulti_;   //!
   TBranch        *b_patJetPfAk05NeutMulti_;   //!
   TBranch        *b_patJetPfAk05CharMulti_;   //!
   TBranch        *b_patJetPfAk05CharHadMulti_;   //!
   TBranch        *b_patJetPfAk05NeutHadMulti_;   //!
   TBranch        *b_patJetPfAk05PhotMulti_;   //!
   TBranch        *b_patJetPfAk05ElecMulti_;   //!
   TBranch        *b_patJetPfAk05HfHadMulti_;   //!
   TBranch        *b_patJetPfAk05HfEmMulti_;   //!
   TBranch        *b_patJetPfAk05PhotEnFr_;   //!
   TBranch        *b_patJetPfAk05MuonEnFr_;   //!
   TBranch        *b_patJetPfAk05HfHadEnFr_;   //!
   TBranch        *b_patJetPfAk05HfEmEnFr_;   //!
   TBranch        *b_patJetPfAk05NeutEmEFr_;   //!
   TBranch        *b_patJetPfAk05CharHadEFr_;   //!
   TBranch        *b_patJetPfAk05NeutHadEFr_;   //!
   TBranch        *b_patJetPfAk05CharEmEFr_;   //!
   TBranch        *b_patJetPfAk05CharMuEFr_;   //!
   TBranch        *b_patJetPfAk05N60_;   //!
   TBranch        *b_patJetPfAk05N90_;   //!
   TBranch        *b_patJetPfAk05EmFr_;   //!
   TBranch        *b_patJetPfAk05HadFr_;   //!
   TBranch        *b_patJetPfAk05EmEbEn_;   //!
   TBranch        *b_patJetPfAk05EmEeEn_;   //!
   TBranch        *b_patJetPfAk05EmHfEn_;   //!
   TBranch        *b_patJetPfAk05HadHbEn_;   //!
   TBranch        *b_patJetPfAk05HadHeEn_;   //!
   TBranch        *b_patJetPfAk05HadHfEn_;   //!
   TBranch        *b_patJetPfAk05HadHoEn_;   //!
   TBranch        *b_patJetPfAk05TotalEm_;   //!
   TBranch        *b_patJetPfAk05TotalHad_;   //!
   TBranch        *b_patJetPfAk05EmEbFr_;   //!
   TBranch        *b_patJetPfAk05EmEeFr_;   //!
   TBranch        *b_patJetPfAk05EmHfFr_;   //!
   TBranch        *b_patJetPfAk05HadHbFr_;   //!
   TBranch        *b_patJetPfAk05HadHeFr_;   //!
   TBranch        *b_patJetPfAk05HadHfFr_;   //!
   TBranch        *b_patJetPfAk05HadHoFr_;   //!
   TBranch        *b_patJetPfAk05JesUncert_;   //!
   TBranch        *b_patPhotonrho25;   //!
   TBranch        *b_patPhotonNumPh_;   //!
   TBranch        *b_patPhotonPt;   //!
   TBranch        *b_patPhotonEta;   //!
   TBranch        *b_patPhotonPhi;   //!
   TBranch        *b_patPhotonEt;   //!
   TBranch        *b_patPhotonEnergy;   //!
   TBranch        *b_patPhotonPx;   //!
   TBranch        *b_patPhotonPy;   //!
   TBranch        *b_patPhotonPz;   //!
   TBranch        *b_patPhotonR9;   //!
   TBranch        *b_patPhotonPhiWidth;   //!
   TBranch        *b_patPhotonEtaWidth;   //!
   TBranch        *b_patPhotonScPhi;   //!
   TBranch        *b_patPhotonScEta;   //!
   TBranch        *b_patPhotonSigmaIetaIeta;   //!
   TBranch        *b_patPhotonhadronicOverEm;   //!
   TBranch        *b_patPhotonecalRecHitSumEtConeDR04;   //!
   TBranch        *b_patPhotonhcalTowerSumEtConeDR04;   //!
   TBranch        *b_patPhotontrkSumPtHollowConeDR04;   //!
   TBranch        *b_patPhotonhasPixelSeed;   //!
   TBranch        *b_patPhotonMCpthat;   //!
   TBranch        *b_patPhotonisGenMatched;   //!
   TBranch        *b_patPhotongenMomId;   //!
   TBranch        *b_patPhotongenGrandMomId;   //!
   TBranch        *b_patPhotongenNSiblings;   //!
   TBranch        *b_patPhotongenMatchedE;   //!
   TBranch        *b_patPhotongenMatchedPx;   //!
   TBranch        *b_patPhotongenMatchedPy;   //!
   TBranch        *b_patPhotongenMatchedPz;   //!
   TBranch        *b_patPhotongenMatchedPt;   //!
   TBranch        *b_patPhotongenMatchedEta;   //!
   TBranch        *b_patPhotongenMatchedPhi;   //!
   TBranch        *b_patPhotongenCalIsoDR04;   //!
   TBranch        *b_patPhotongenTrkIsoDR04;   //!
   TBranch        *b_patPhotongenIsoDR04;   //!

   vector_angular(std::string filename, TTree *tree=0);
   virtual ~vector_angular();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(bool onlyOneJet=false,bool DEBUG=false);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Bool_t   isGoodPho(Int_t ipho, bool applyPileupCorr=true);
   virtual Double_t phoEcalIso(Int_t ipho, bool applyPileupCorr=true);
   virtual Double_t phoHcalIso(Int_t ipho, bool applyPileupCorr=true);
   virtual Double_t phoTrkIso(Int_t ipho, bool applyPileupCorr=true);
   virtual Bool_t   isGoodLooseJet(Int_t ijet);
   virtual Bool_t   isGoodMediumJet(Int_t ijet);
   virtual Bool_t   isGoodTightJet(Int_t ijet);
   virtual Int_t    rec_phoDecCode(Int_t ipho);
   virtual Int_t    gen_phoDecCode(Int_t ipho);
   virtual Int_t    phoDecCode(double eta);
   virtual Bool_t   rec_isFidPho (Int_t ipho);
   virtual Bool_t   gen_isFidPho (Int_t ipho);
   virtual Bool_t   isFidPho(double pt, double eta);
   virtual Bool_t   rec_isFidJet(Int_t ijet);
   virtual Bool_t   gen_isFidJet(Int_t ijet);
   virtual Bool_t   isFidJet (double pt, double eta);
   virtual Int_t    matchRecoToGenJet(Int_t ijet); // return matched index in genjet
   virtual Int_t    matchGenToRecoJet(Int_t ijet); // return matched index in recojet

   std::string _inputFileName;

};

#endif

#ifdef vector_angular_cxx
vector_angular::vector_angular(std::string filename, TTree *tree)
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
  Init(tree);
  _inputFileName = filename;

}

vector_angular::~vector_angular()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t vector_angular::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t vector_angular::LoadTree(Long64_t entry)
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

void vector_angular::Init(TTree *tree)
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
   genParQ_ = 0;
   genParId_ = 0;
   genParSt_ = 0;
   genMomParId_ = 0;
   genParIndex_ = 0;
   genJetE_ = 0;
   genJetPt_ = 0;
   genJetEta_ = 0;
   genJetPhi_ = 0;
   patJetPfAk05Pt_ = 0;
   patJetPfAk05Eta_ = 0;
   patJetPfAk05Phi_ = 0;
   patJetPfAk05M_ = 0;
   patJetPfAk05Rapidity_ = 0;
   patJetPfAk05Px_ = 0;
   patJetPfAk05Py_ = 0;
   patJetPfAk05Pz_ = 0;
   patJetPfAk05En_ = 0;
   patJetPfAk05UnCorrPt_ = 0;
   patJetPfAk05UnCorrPx_ = 0;
   patJetPfAk05UnCorrPy_ = 0;
   patJetPfAk05UnCorrPz_ = 0;
   patJetPfAk05UnCorrEnt_ = 0;
   patJetPfAk05PhotEn_ = 0;
   patJetPfAk05ElecEn_ = 0;
   patJetPfAk05MuonEn_ = 0;
   patJetPfAk05HfHadEn_ = 0;
   patJetPfAk05HfEmEn_ = 0;
   patJetPfAk05CharHadE_ = 0;
   patJetPfAk05NeutHadE_ = 0;
   patJetPfAk05CharEmE_ = 0;
   patJetPfAk05CharMuE_ = 0;
   patJetPfAk05NeutEmE_ = 0;
   patJetPfAk05MuonMulti_ = 0;
   patJetPfAk05NeutMulti_ = 0;
   patJetPfAk05CharMulti_ = 0;
   patJetPfAk05CharHadMulti_ = 0;
   patJetPfAk05NeutHadMulti_ = 0;
   patJetPfAk05PhotMulti_ = 0;
   patJetPfAk05ElecMulti_ = 0;
   patJetPfAk05HfHadMulti_ = 0;
   patJetPfAk05HfEmMulti_ = 0;
   patJetPfAk05PhotEnFr_ = 0;
   patJetPfAk05MuonEnFr_ = 0;
   patJetPfAk05HfHadEnFr_ = 0;
   patJetPfAk05HfEmEnFr_ = 0;
   patJetPfAk05NeutEmEFr_ = 0;
   patJetPfAk05CharHadEFr_ = 0;
   patJetPfAk05NeutHadEFr_ = 0;
   patJetPfAk05CharEmEFr_ = 0;
   patJetPfAk05CharMuEFr_ = 0;
   patJetPfAk05N60_ = 0;
   patJetPfAk05N90_ = 0;
   patJetPfAk05EmFr_ = 0;
   patJetPfAk05HadFr_ = 0;
   patJetPfAk05EmEbEn_ = 0;
   patJetPfAk05EmEeEn_ = 0;
   patJetPfAk05EmHfEn_ = 0;
   patJetPfAk05HadHbEn_ = 0;
   patJetPfAk05HadHeEn_ = 0;
   patJetPfAk05HadHfEn_ = 0;
   patJetPfAk05HadHoEn_ = 0;
   patJetPfAk05TotalEm_ = 0;
   patJetPfAk05TotalHad_ = 0;
   patJetPfAk05EmEbFr_ = 0;
   patJetPfAk05EmEeFr_ = 0;
   patJetPfAk05EmHfFr_ = 0;
   patJetPfAk05HadHbFr_ = 0;
   patJetPfAk05HadHeFr_ = 0;
   patJetPfAk05HadHfFr_ = 0;
   patJetPfAk05HadHoFr_ = 0;
   patJetPfAk05JesUncert_ = 0;
   patPhotonPt = 0;
   patPhotonEta = 0;
   patPhotonPhi = 0;
   patPhotonEt = 0;
   patPhotonEnergy = 0;
   patPhotonPx = 0;
   patPhotonPy = 0;
   patPhotonPz = 0;
   patPhotonR9 = 0;
   patPhotonPhiWidth = 0;
   patPhotonEtaWidth = 0;
   patPhotonScPhi = 0;
   patPhotonScEta = 0;
   patPhotonSigmaIetaIeta = 0;
   patPhotonhadronicOverEm = 0;
   patPhotonecalRecHitSumEtConeDR04 = 0;
   patPhotonhcalTowerSumEtConeDR04 = 0;
   patPhotontrkSumPtHollowConeDR04 = 0;
   patPhotonhasPixelSeed = 0;
   patPhotonisGenMatched = 0;
   patPhotongenMomId = 0;
   patPhotongenGrandMomId = 0;
   patPhotongenNSiblings = 0;
   patPhotongenMatchedE = 0;
   patPhotongenMatchedPx = 0;
   patPhotongenMatchedPy = 0;
   patPhotongenMatchedPz = 0;
   patPhotongenMatchedPt = 0;
   patPhotongenMatchedEta = 0;
   patPhotongenMatchedPhi = 0;
   patPhotongenCalIsoDR04 = 0;
   patPhotongenTrkIsoDR04 = 0;
   patPhotongenIsoDR04 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("PU_weight", &PU_weight, &b_PU_weight);
   // fChain->SetBranchAddress("EvtInfo_EventNum", &EvtInfo_EventNum, &b_EvtInfo_EventNum);
   // fChain->SetBranchAddress("EvtInfo_RunNum", &EvtInfo_RunNum, &b_EvtInfo_RunNum);
   // fChain->SetBranchAddress("EvtInfo_LumiSection", &EvtInfo_LumiSection, &b_EvtInfo_LumiSection);
   // fChain->SetBranchAddress("EvtInfo_BunchXing", &EvtInfo_BunchXing, &b_EvtInfo_BunchXing);
   fChain->SetBranchAddress("EvtInfo_NumVtx", &EvtInfo_NumVtx, &b_EvtInfo_NumVtx);
   // fChain->SetBranchAddress("EvtInfo_VertexX_", &EvtInfo_VertexX_, &b_EvtInfo_VertexX_);
   // fChain->SetBranchAddress("EvtInfo_VertexY_", &EvtInfo_VertexY_, &b_EvtInfo_VertexY_);
   // fChain->SetBranchAddress("EvtInfo_VertexZ_", &EvtInfo_VertexZ_, &b_EvtInfo_VertexZ_);
   fChain->SetBranchAddress("ptHat_", &ptHat_, &b_ptHat_);
   fChain->SetBranchAddress("mcWeight_", &mcWeight_, &b_mcWeight_);
   fChain->SetBranchAddress("genParE_", &genParE_, &b_genParE_);
   fChain->SetBranchAddress("genParPt_", &genParPt_, &b_genParPt_);
   fChain->SetBranchAddress("genParEta_", &genParEta_, &b_genParEta_);
   fChain->SetBranchAddress("genParPhi_", &genParPhi_, &b_genParPhi_);
   fChain->SetBranchAddress("genParQ_", &genParQ_, &b_genParQ_);
   fChain->SetBranchAddress("genParId_", &genParId_, &b_genParId_);
   fChain->SetBranchAddress("genParSt_", &genParSt_, &b_genParSt_);
   fChain->SetBranchAddress("genMomParId_", &genMomParId_, &b_genMomParId_);
   fChain->SetBranchAddress("genParIndex_", &genParIndex_, &b_genParIndex_);
   fChain->SetBranchAddress("genJetE_", &genJetE_, &b_genJetE_);
   fChain->SetBranchAddress("genJetPt_", &genJetPt_, &b_genJetPt_);
   fChain->SetBranchAddress("genJetEta_", &genJetEta_, &b_genJetEta_);
   fChain->SetBranchAddress("genJetPhi_", &genJetPhi_, &b_genJetPhi_);
   fChain->SetBranchAddress("patJetPfAk05Pt_", &patJetPfAk05Pt_, &b_patJetPfAk05Pt_);
   fChain->SetBranchAddress("patJetPfAk05Eta_", &patJetPfAk05Eta_, &b_patJetPfAk05Eta_);
   fChain->SetBranchAddress("patJetPfAk05Phi_", &patJetPfAk05Phi_, &b_patJetPfAk05Phi_);
   fChain->SetBranchAddress("patJetPfAk05M_", &patJetPfAk05M_, &b_patJetPfAk05M_);
   fChain->SetBranchAddress("patJetPfAk05Rapidity_", &patJetPfAk05Rapidity_, &b_patJetPfAk05Rapidity_);
   // fChain->SetBranchAddress("patJetPfAk05Px_", &patJetPfAk05Px_, &b_patJetPfAk05Px_);
   // fChain->SetBranchAddress("patJetPfAk05Py_", &patJetPfAk05Py_, &b_patJetPfAk05Py_);
   // fChain->SetBranchAddress("patJetPfAk05Pz_", &patJetPfAk05Pz_, &b_patJetPfAk05Pz_);
   fChain->SetBranchAddress("patJetPfAk05En_", &patJetPfAk05En_, &b_patJetPfAk05En_);
   // fChain->SetBranchAddress("patJetPfAk05UnCorrPt_", &patJetPfAk05UnCorrPt_, &b_patJetPfAk05UnCorrPt_);
   // fChain->SetBranchAddress("patJetPfAk05UnCorrPx_", &patJetPfAk05UnCorrPx_, &b_patJetPfAk05UnCorrPx_);
   // fChain->SetBranchAddress("patJetPfAk05UnCorrPy_", &patJetPfAk05UnCorrPy_, &b_patJetPfAk05UnCorrPy_);
   // fChain->SetBranchAddress("patJetPfAk05UnCorrPz_", &patJetPfAk05UnCorrPz_, &b_patJetPfAk05UnCorrPz_);
   // fChain->SetBranchAddress("patJetPfAk05UnCorrEnt_", &patJetPfAk05UnCorrEnt_, &b_patJetPfAk05UnCorrEnt_);
   // fChain->SetBranchAddress("patJetPfAk05PhotEn_", &patJetPfAk05PhotEn_, &b_patJetPfAk05PhotEn_);
   // fChain->SetBranchAddress("patJetPfAk05ElecEn_", &patJetPfAk05ElecEn_, &b_patJetPfAk05ElecEn_);
   // fChain->SetBranchAddress("patJetPfAk05MuonEn_", &patJetPfAk05MuonEn_, &b_patJetPfAk05MuonEn_);
   // fChain->SetBranchAddress("patJetPfAk05HfHadEn_", &patJetPfAk05HfHadEn_, &b_patJetPfAk05HfHadEn_);
   // fChain->SetBranchAddress("patJetPfAk05HfEmEn_", &patJetPfAk05HfEmEn_, &b_patJetPfAk05HfEmEn_);
   // fChain->SetBranchAddress("patJetPfAk05CharHadE_", &patJetPfAk05CharHadE_, &b_patJetPfAk05CharHadE_);
   // fChain->SetBranchAddress("patJetPfAk05NeutHadE_", &patJetPfAk05NeutHadE_, &b_patJetPfAk05NeutHadE_);
   // fChain->SetBranchAddress("patJetPfAk05CharEmE_", &patJetPfAk05CharEmE_, &b_patJetPfAk05CharEmE_);
   // fChain->SetBranchAddress("patJetPfAk05CharMuE_", &patJetPfAk05CharMuE_, &b_patJetPfAk05CharMuE_);
   // fChain->SetBranchAddress("patJetPfAk05NeutEmE_", &patJetPfAk05NeutEmE_, &b_patJetPfAk05NeutEmE_);
   // fChain->SetBranchAddress("patJetPfAk05MuonMulti_", &patJetPfAk05MuonMulti_, &b_patJetPfAk05MuonMulti_);
   // fChain->SetBranchAddress("patJetPfAk05NeutMulti_", &patJetPfAk05NeutMulti_, &b_patJetPfAk05NeutMulti_);
   fChain->SetBranchAddress("patJetPfAk05CharMulti_", &patJetPfAk05CharMulti_, &b_patJetPfAk05CharMulti_);
   // fChain->SetBranchAddress("patJetPfAk05CharHadMulti_", &patJetPfAk05CharHadMulti_, &b_patJetPfAk05CharHadMulti_);
   // fChain->SetBranchAddress("patJetPfAk05NeutHadMulti_", &patJetPfAk05NeutHadMulti_, &b_patJetPfAk05NeutHadMulti_);
   // fChain->SetBranchAddress("patJetPfAk05PhotMulti_", &patJetPfAk05PhotMulti_, &b_patJetPfAk05PhotMulti_);
   // fChain->SetBranchAddress("patJetPfAk05ElecMulti_", &patJetPfAk05ElecMulti_, &b_patJetPfAk05ElecMulti_);
   // fChain->SetBranchAddress("patJetPfAk05HfHadMulti_", &patJetPfAk05HfHadMulti_, &b_patJetPfAk05HfHadMulti_);
   // fChain->SetBranchAddress("patJetPfAk05HfEmMulti_", &patJetPfAk05HfEmMulti_, &b_patJetPfAk05HfEmMulti_);
   // fChain->SetBranchAddress("patJetPfAk05PhotEnFr_", &patJetPfAk05PhotEnFr_, &b_patJetPfAk05PhotEnFr_);
   // fChain->SetBranchAddress("patJetPfAk05MuonEnFr_", &patJetPfAk05MuonEnFr_, &b_patJetPfAk05MuonEnFr_);
   fChain->SetBranchAddress("patJetPfAk05HfHadEnFr_", &patJetPfAk05HfHadEnFr_, &b_patJetPfAk05HfHadEnFr_);
   fChain->SetBranchAddress("patJetPfAk05HfEmEnFr_", &patJetPfAk05HfEmEnFr_, &b_patJetPfAk05HfEmEnFr_);
   fChain->SetBranchAddress("patJetPfAk05NeutEmEFr_", &patJetPfAk05NeutEmEFr_, &b_patJetPfAk05NeutEmEFr_);
   fChain->SetBranchAddress("patJetPfAk05CharHadEFr_", &patJetPfAk05CharHadEFr_, &b_patJetPfAk05CharHadEFr_);
   fChain->SetBranchAddress("patJetPfAk05NeutHadEFr_", &patJetPfAk05NeutHadEFr_, &b_patJetPfAk05NeutHadEFr_);
   fChain->SetBranchAddress("patJetPfAk05CharEmEFr_", &patJetPfAk05CharEmEFr_, &b_patJetPfAk05CharEmEFr_);
   fChain->SetBranchAddress("patJetPfAk05CharMuEFr_", &patJetPfAk05CharMuEFr_, &b_patJetPfAk05CharMuEFr_);
   fChain->SetBranchAddress("patJetPfAk05N60_", &patJetPfAk05N60_, &b_patJetPfAk05N60_);
   fChain->SetBranchAddress("patJetPfAk05N90_", &patJetPfAk05N90_, &b_patJetPfAk05N90_);
   fChain->SetBranchAddress("patJetPfAk05EmFr_", &patJetPfAk05EmFr_, &b_patJetPfAk05EmFr_);
   fChain->SetBranchAddress("patJetPfAk05HadFr_", &patJetPfAk05HadFr_, &b_patJetPfAk05HadFr_);
   // fChain->SetBranchAddress("patJetPfAk05EmEbEn_", &patJetPfAk05EmEbEn_, &b_patJetPfAk05EmEbEn_);
   // fChain->SetBranchAddress("patJetPfAk05EmEeEn_", &patJetPfAk05EmEeEn_, &b_patJetPfAk05EmEeEn_);
   // fChain->SetBranchAddress("patJetPfAk05EmHfEn_", &patJetPfAk05EmHfEn_, &b_patJetPfAk05EmHfEn_);
   // fChain->SetBranchAddress("patJetPfAk05HadHbEn_", &patJetPfAk05HadHbEn_, &b_patJetPfAk05HadHbEn_);
   // fChain->SetBranchAddress("patJetPfAk05HadHeEn_", &patJetPfAk05HadHeEn_, &b_patJetPfAk05HadHeEn_);
   // fChain->SetBranchAddress("patJetPfAk05HadHfEn_", &patJetPfAk05HadHfEn_, &b_patJetPfAk05HadHfEn_);
   // fChain->SetBranchAddress("patJetPfAk05HadHoEn_", &patJetPfAk05HadHoEn_, &b_patJetPfAk05HadHoEn_);
   // fChain->SetBranchAddress("patJetPfAk05TotalEm_", &patJetPfAk05TotalEm_, &b_patJetPfAk05TotalEm_);
   // fChain->SetBranchAddress("patJetPfAk05TotalHad_", &patJetPfAk05TotalHad_, &b_patJetPfAk05TotalHad_);
   // fChain->SetBranchAddress("patJetPfAk05EmEbFr_", &patJetPfAk05EmEbFr_, &b_patJetPfAk05EmEbFr_);
   // fChain->SetBranchAddress("patJetPfAk05EmEeFr_", &patJetPfAk05EmEeFr_, &b_patJetPfAk05EmEeFr_);
   // fChain->SetBranchAddress("patJetPfAk05EmHfFr_", &patJetPfAk05EmHfFr_, &b_patJetPfAk05EmHfFr_);
   // fChain->SetBranchAddress("patJetPfAk05HadHbFr_", &patJetPfAk05HadHbFr_, &b_patJetPfAk05HadHbFr_);
   // fChain->SetBranchAddress("patJetPfAk05HadHeFr_", &patJetPfAk05HadHeFr_, &b_patJetPfAk05HadHeFr_);
   // fChain->SetBranchAddress("patJetPfAk05HadHfFr_", &patJetPfAk05HadHfFr_, &b_patJetPfAk05HadHfFr_);
   // fChain->SetBranchAddress("patJetPfAk05HadHoFr_", &patJetPfAk05HadHoFr_, &b_patJetPfAk05HadHoFr_);
   // fChain->SetBranchAddress("patJetPfAk05JesUncert_", &patJetPfAk05JesUncert_, &b_patJetPfAk05JesUncert_);
   fChain->SetBranchAddress("patPhotonrho25", &patPhotonrho25, &b_patPhotonrho25);
   fChain->SetBranchAddress("patPhotonNumPh_", &patPhotonNumPh_, &b_patPhotonNumPh_);
   fChain->SetBranchAddress("patPhotonPt", &patPhotonPt, &b_patPhotonPt);
   fChain->SetBranchAddress("patPhotonEta", &patPhotonEta, &b_patPhotonEta);
   fChain->SetBranchAddress("patPhotonPhi", &patPhotonPhi, &b_patPhotonPhi);
   fChain->SetBranchAddress("patPhotonEt", &patPhotonEt, &b_patPhotonEt);
   fChain->SetBranchAddress("patPhotonEnergy", &patPhotonEnergy, &b_patPhotonEnergy);
   // fChain->SetBranchAddress("patPhotonPx", &patPhotonPx, &b_patPhotonPx);
   // fChain->SetBranchAddress("patPhotonPy", &patPhotonPy, &b_patPhotonPy);
   // fChain->SetBranchAddress("patPhotonPz", &patPhotonPz, &b_patPhotonPz);
   // fChain->SetBranchAddress("patPhotonR9", &patPhotonR9, &b_patPhotonR9);
   // fChain->SetBranchAddress("patPhotonPhiWidth", &patPhotonPhiWidth, &b_patPhotonPhiWidth);
   // fChain->SetBranchAddress("patPhotonEtaWidth", &patPhotonEtaWidth, &b_patPhotonEtaWidth);
   fChain->SetBranchAddress("patPhotonScPhi", &patPhotonScPhi, &b_patPhotonScPhi);
   fChain->SetBranchAddress("patPhotonScEta", &patPhotonScEta, &b_patPhotonScEta);
   fChain->SetBranchAddress("patPhotonSigmaIetaIeta", &patPhotonSigmaIetaIeta, &b_patPhotonSigmaIetaIeta);
   fChain->SetBranchAddress("patPhotonhadronicOverEm", &patPhotonhadronicOverEm, &b_patPhotonhadronicOverEm);
   fChain->SetBranchAddress("patPhotonecalRecHitSumEtConeDR04", &patPhotonecalRecHitSumEtConeDR04, &b_patPhotonecalRecHitSumEtConeDR04);
   fChain->SetBranchAddress("patPhotonhcalTowerSumEtConeDR04", &patPhotonhcalTowerSumEtConeDR04, &b_patPhotonhcalTowerSumEtConeDR04);
   fChain->SetBranchAddress("patPhotontrkSumPtHollowConeDR04", &patPhotontrkSumPtHollowConeDR04, &b_patPhotontrkSumPtHollowConeDR04);
   fChain->SetBranchAddress("patPhotonhasPixelSeed", &patPhotonhasPixelSeed, &b_patPhotonhasPixelSeed);
   fChain->SetBranchAddress("patPhotonMCpthat", &patPhotonMCpthat, &b_patPhotonMCpthat);
   fChain->SetBranchAddress("patPhotonisGenMatched", &patPhotonisGenMatched, &b_patPhotonisGenMatched);
   fChain->SetBranchAddress("patPhotongenMomId", &patPhotongenMomId, &b_patPhotongenMomId);
   fChain->SetBranchAddress("patPhotongenGrandMomId", &patPhotongenGrandMomId, &b_patPhotongenGrandMomId);
   fChain->SetBranchAddress("patPhotongenNSiblings", &patPhotongenNSiblings, &b_patPhotongenNSiblings);
   fChain->SetBranchAddress("patPhotongenMatchedE", &patPhotongenMatchedE, &b_patPhotongenMatchedE);
   fChain->SetBranchAddress("patPhotongenMatchedPx", &patPhotongenMatchedPx, &b_patPhotongenMatchedPx);
   fChain->SetBranchAddress("patPhotongenMatchedPy", &patPhotongenMatchedPy, &b_patPhotongenMatchedPy);
   fChain->SetBranchAddress("patPhotongenMatchedPz", &patPhotongenMatchedPz, &b_patPhotongenMatchedPz);
   fChain->SetBranchAddress("patPhotongenMatchedPt", &patPhotongenMatchedPt, &b_patPhotongenMatchedPt);
   fChain->SetBranchAddress("patPhotongenMatchedEta", &patPhotongenMatchedEta, &b_patPhotongenMatchedEta);
   fChain->SetBranchAddress("patPhotongenMatchedPhi", &patPhotongenMatchedPhi, &b_patPhotongenMatchedPhi);
   fChain->SetBranchAddress("patPhotongenCalIsoDR04", &patPhotongenCalIsoDR04, &b_patPhotongenCalIsoDR04);
   fChain->SetBranchAddress("patPhotongenTrkIsoDR04", &patPhotongenTrkIsoDR04, &b_patPhotongenTrkIsoDR04);
   fChain->SetBranchAddress("patPhotongenIsoDR04", &patPhotongenIsoDR04, &b_patPhotongenIsoDR04);
   Notify();
}

Bool_t vector_angular::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void vector_angular::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t vector_angular::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef vector_angular_cxx
