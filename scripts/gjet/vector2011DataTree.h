//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jan 14 16:12:15 2012 by ROOT version 5.27/06b
// from TTree tree/tree
// found on file: /data4/yunju/DData/May10rereco/TreeMay10ReReco_v9/TreeMay10ReReco_v9_11to12.root
//////////////////////////////////////////////////////////

#ifndef vector2011DataTree_h
#define vector2011DataTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TChain.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <string>
#include <iostream>

using namespace std;


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
   const Int_t kMaxPhotonNumPh = 1;

class vector2011DataTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           EvtInfo_EventNum;
   Int_t           EvtInfo_RunNum;
   Int_t           EvtInfo_LumiSection;
   Int_t           EvtInfo_BunchXing;
   Int_t           EvtInfo_nVtxGood;
   Int_t           EvtInfo_nVtxNotFake;
   vector<double>  *EvtInfo_vertexX;
   vector<double>  *EvtInfo_vertexY;
   vector<double>  *EvtInfo_vertexZ;
   vector<double>  *EvtInfo_vertexXError;
   vector<double>  *EvtInfo_vertexYError;
   vector<double>  *EvtInfo_vertexZError;
   vector<double>  *EvtInfo_vertexChi2;
   vector<double>  *EvtInfo_vertexNormChi2;
   vector<double>  *EvtInfo_vertexNdof;
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
   vector<int>     *trigResults;
   vector<int>     *trigPrescale;
   vector<string>  *trigName;
   Double_t        Photonrho25;
   Double_t        Photonrho44;
   Double_t        PhotonNumPh_;
   vector<double>  *PhotonPt;
   vector<double>  *PhotonEta;
   vector<double>  *PhotonPhi;
   vector<double>  *PhotonEt;
   vector<double>  *PhotonEnergy;
   vector<double>  *PhotonPx;
   vector<double>  *PhotonPy;
   vector<double>  *PhotonPz;
   vector<double>  *PhotonR9;
   vector<double>  *PhotonPhiWidth;
   vector<double>  *PhotonEtaWidth;
   vector<double>  *PhotonScPhi;
   vector<double>  *PhotonScEta;
   vector<double>  *PhotonSigmaIetaIeta;
   vector<double>  *PhotonSeedTime;
   vector<double>  *PhotonseedSeverity;
   vector<double>  *PhotonhadronicOverEm;
   vector<double>  *PhotonecalRecHitSumEtConeDR04;
   vector<double>  *PhotonhcalTowerSumEtConeDR04;
   vector<double>  *PhotontrkSumPtHollowConeDR04;
   vector<double>  *PhotonhasPixelSeed;

   // List of branches
   TBranch        *b_EvtInfo_EventNum;   //!
   TBranch        *b_EvtInfo_RunNum;   //!
   TBranch        *b_EvtInfo_LumiSection;   //!
   TBranch        *b_EvtInfo_BunchXing;   //!
   TBranch        *b_EvtInfo_nVtxGood;   //!
   TBranch        *b_EvtInfo_nVtxNotFake;   //!
   TBranch        *b_EvtInfo_vertexX;   //!
   TBranch        *b_EvtInfo_vertexY;   //!
   TBranch        *b_EvtInfo_vertexZ;   //!
   TBranch        *b_EvtInfo_vertexXError;   //!
   TBranch        *b_EvtInfo_vertexYError;   //!
   TBranch        *b_EvtInfo_vertexZError;   //!
   TBranch        *b_EvtInfo_vertexChi2;   //!
   TBranch        *b_EvtInfo_vertexNormChi2;   //!
   TBranch        *b_EvtInfo_vertexNdof;   //!
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
   TBranch        *b_trigResults;   //!
   TBranch        *b_trigPrescale;   //!
   TBranch        *b_trigName;   //!
   TBranch        *b_Photonrho25;   //!
   TBranch        *b_Photonrho44;   //!
   TBranch        *b_PhotonNumPh_;   //!
   TBranch        *b_PhotonPt;   //!
   TBranch        *b_PhotonEta;   //!
   TBranch        *b_PhotonPhi;   //!
   TBranch        *b_PhotonEt;   //!
   TBranch        *b_PhotonEnergy;   //!
   TBranch        *b_PhotonPx;   //!
   TBranch        *b_PhotonPy;   //!
   TBranch        *b_PhotonPz;   //!
   TBranch        *b_PhotonR9;   //!
   TBranch        *b_PhotonPhiWidth;   //!
   TBranch        *b_PhotonEtaWidth;   //!
   TBranch        *b_PhotonScPhi;   //!
   TBranch        *b_PhotonScEta;   //!
   TBranch        *b_PhotonSigmaIetaIeta;   //!
   TBranch        *b_PhotonSeedTime;   //!
   TBranch        *b_PhotonseedSeverity;   //!
   TBranch        *b_PhotonhadronicOverEm;   //!
   TBranch        *b_PhotonecalRecHitSumEtConeDR04;   //!
   TBranch        *b_PhotonhcalTowerSumEtConeDR04;   //!
   TBranch        *b_PhotontrkSumPtHollowConeDR04;   //!
   TBranch        *b_PhotonhasPixelSeed;   //!

   vector2011DataTree(std::string dirName, TTree *tree=0);
   virtual ~vector2011DataTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Bool_t   isGoodLooseJet(Long64_t entry, Int_t ijet);
   virtual Bool_t   isGoodMediumJet(Long64_t entry, Int_t ijet);
   virtual Bool_t   isGoodTightJet(Long64_t entry, Int_t ijet);
   virtual Bool_t   isFidJet (Long64_t entry, Int_t ijet);

   std::string _inputDirName;

};

#endif

#ifdef vector2011DataTree_cxx
vector2011DataTree::vector2011DataTree(std::string dirName, TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     
     TChain* pho = new TChain("tree");

     TSystemDirectory *base = new TSystemDirectory("root","root");
     base->SetDirectory(dirName.data());
     TList *listOfFiles = base->GetListOfFiles();
     TIter fileIt(listOfFiles);
     TFile *fileH = new TFile();
     int nfile=0;
     while(fileH = (TFile*)fileIt()) {
       std::string fileN = fileH->GetName();
       std::string baseString = "TreeMay";
       if( fileH->IsFolder())  continue;
       if(fileN.find(baseString) == std::string::npos)continue;
       cout << fileN.data() << endl;
       nfile++;
       pho->Add(fileN.data());
     }

     cout << "Opening " << nfile << " files " << endl;
     tree = pho;
     _inputDirName = dirName;
   }


   cout << "Number of entries are " << tree->GetEntries() << endl;
   Init(tree);
}

vector2011DataTree::~vector2011DataTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t vector2011DataTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t vector2011DataTree::LoadTree(Long64_t entry)
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

void vector2011DataTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   EvtInfo_vertexX = 0;
   EvtInfo_vertexY = 0;
   EvtInfo_vertexZ = 0;
   EvtInfo_vertexXError = 0;
   EvtInfo_vertexYError = 0;
   EvtInfo_vertexZError = 0;
   EvtInfo_vertexChi2 = 0;
   EvtInfo_vertexNormChi2 = 0;
   EvtInfo_vertexNdof = 0;
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
   trigResults = 0;
   trigPrescale = 0;
   trigName = 0;
   PhotonPt = 0;
   PhotonEta = 0;
   PhotonPhi = 0;
   PhotonEt = 0;
   PhotonEnergy = 0;
   PhotonPx = 0;
   PhotonPy = 0;
   PhotonPz = 0;
   PhotonR9 = 0;
   PhotonPhiWidth = 0;
   PhotonEtaWidth = 0;
   PhotonScPhi = 0;
   PhotonScEta = 0;
   PhotonSigmaIetaIeta = 0;
   PhotonSeedTime = 0;
   PhotonseedSeverity = 0;
   PhotonhadronicOverEm = 0;
   PhotonecalRecHitSumEtConeDR04 = 0;
   PhotonhcalTowerSumEtConeDR04 = 0;
   PhotontrkSumPtHollowConeDR04 = 0;
   PhotonhasPixelSeed = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EvtInfo_EventNum", &EvtInfo_EventNum, &b_EvtInfo_EventNum);
   fChain->SetBranchAddress("EvtInfo_RunNum", &EvtInfo_RunNum, &b_EvtInfo_RunNum);
   fChain->SetBranchAddress("EvtInfo_LumiSection", &EvtInfo_LumiSection, &b_EvtInfo_LumiSection);
   fChain->SetBranchAddress("EvtInfo_BunchXing", &EvtInfo_BunchXing, &b_EvtInfo_BunchXing);
   fChain->SetBranchAddress("EvtInfo_nVtxGood", &EvtInfo_nVtxGood, &b_EvtInfo_nVtxGood);
   fChain->SetBranchAddress("EvtInfo_nVtxNotFake", &EvtInfo_nVtxNotFake, &b_EvtInfo_nVtxNotFake);
   fChain->SetBranchAddress("EvtInfo_vertexX", &EvtInfo_vertexX, &b_EvtInfo_vertexX);
   fChain->SetBranchAddress("EvtInfo_vertexY", &EvtInfo_vertexY, &b_EvtInfo_vertexY);
   fChain->SetBranchAddress("EvtInfo_vertexZ", &EvtInfo_vertexZ, &b_EvtInfo_vertexZ);
   fChain->SetBranchAddress("EvtInfo_vertexXError", &EvtInfo_vertexXError, &b_EvtInfo_vertexXError);
   fChain->SetBranchAddress("EvtInfo_vertexYError", &EvtInfo_vertexYError, &b_EvtInfo_vertexYError);
   fChain->SetBranchAddress("EvtInfo_vertexZError", &EvtInfo_vertexZError, &b_EvtInfo_vertexZError);
   fChain->SetBranchAddress("EvtInfo_vertexChi2", &EvtInfo_vertexChi2, &b_EvtInfo_vertexChi2);
   fChain->SetBranchAddress("EvtInfo_vertexNormChi2", &EvtInfo_vertexNormChi2, &b_EvtInfo_vertexNormChi2);
   fChain->SetBranchAddress("EvtInfo_vertexNdof", &EvtInfo_vertexNdof, &b_EvtInfo_vertexNdof);
   fChain->SetBranchAddress("patJetPfAk05Pt_", &patJetPfAk05Pt_, &b_patJetPfAk05Pt_);
   fChain->SetBranchAddress("patJetPfAk05Eta_", &patJetPfAk05Eta_, &b_patJetPfAk05Eta_);
   fChain->SetBranchAddress("patJetPfAk05Phi_", &patJetPfAk05Phi_, &b_patJetPfAk05Phi_);
   fChain->SetBranchAddress("patJetPfAk05M_", &patJetPfAk05M_, &b_patJetPfAk05M_);
   fChain->SetBranchAddress("patJetPfAk05Rapidity_", &patJetPfAk05Rapidity_, &b_patJetPfAk05Rapidity_);
   fChain->SetBranchAddress("patJetPfAk05Px_", &patJetPfAk05Px_, &b_patJetPfAk05Px_);
   fChain->SetBranchAddress("patJetPfAk05Py_", &patJetPfAk05Py_, &b_patJetPfAk05Py_);
   fChain->SetBranchAddress("patJetPfAk05Pz_", &patJetPfAk05Pz_, &b_patJetPfAk05Pz_);
   fChain->SetBranchAddress("patJetPfAk05En_", &patJetPfAk05En_, &b_patJetPfAk05En_);
   fChain->SetBranchAddress("patJetPfAk05UnCorrPt_", &patJetPfAk05UnCorrPt_, &b_patJetPfAk05UnCorrPt_);
   fChain->SetBranchAddress("patJetPfAk05UnCorrPx_", &patJetPfAk05UnCorrPx_, &b_patJetPfAk05UnCorrPx_);
   fChain->SetBranchAddress("patJetPfAk05UnCorrPy_", &patJetPfAk05UnCorrPy_, &b_patJetPfAk05UnCorrPy_);
   fChain->SetBranchAddress("patJetPfAk05UnCorrPz_", &patJetPfAk05UnCorrPz_, &b_patJetPfAk05UnCorrPz_);
   fChain->SetBranchAddress("patJetPfAk05UnCorrEnt_", &patJetPfAk05UnCorrEnt_, &b_patJetPfAk05UnCorrEnt_);
   fChain->SetBranchAddress("patJetPfAk05PhotEn_", &patJetPfAk05PhotEn_, &b_patJetPfAk05PhotEn_);
   fChain->SetBranchAddress("patJetPfAk05ElecEn_", &patJetPfAk05ElecEn_, &b_patJetPfAk05ElecEn_);
   fChain->SetBranchAddress("patJetPfAk05MuonEn_", &patJetPfAk05MuonEn_, &b_patJetPfAk05MuonEn_);
   fChain->SetBranchAddress("patJetPfAk05HfHadEn_", &patJetPfAk05HfHadEn_, &b_patJetPfAk05HfHadEn_);
   fChain->SetBranchAddress("patJetPfAk05HfEmEn_", &patJetPfAk05HfEmEn_, &b_patJetPfAk05HfEmEn_);
   fChain->SetBranchAddress("patJetPfAk05CharHadE_", &patJetPfAk05CharHadE_, &b_patJetPfAk05CharHadE_);
   fChain->SetBranchAddress("patJetPfAk05NeutHadE_", &patJetPfAk05NeutHadE_, &b_patJetPfAk05NeutHadE_);
   fChain->SetBranchAddress("patJetPfAk05CharEmE_", &patJetPfAk05CharEmE_, &b_patJetPfAk05CharEmE_);
   fChain->SetBranchAddress("patJetPfAk05CharMuE_", &patJetPfAk05CharMuE_, &b_patJetPfAk05CharMuE_);
   fChain->SetBranchAddress("patJetPfAk05NeutEmE_", &patJetPfAk05NeutEmE_, &b_patJetPfAk05NeutEmE_);
   fChain->SetBranchAddress("patJetPfAk05MuonMulti_", &patJetPfAk05MuonMulti_, &b_patJetPfAk05MuonMulti_);
   fChain->SetBranchAddress("patJetPfAk05NeutMulti_", &patJetPfAk05NeutMulti_, &b_patJetPfAk05NeutMulti_);
   fChain->SetBranchAddress("patJetPfAk05CharMulti_", &patJetPfAk05CharMulti_, &b_patJetPfAk05CharMulti_);
   fChain->SetBranchAddress("patJetPfAk05CharHadMulti_", &patJetPfAk05CharHadMulti_, &b_patJetPfAk05CharHadMulti_);
   fChain->SetBranchAddress("patJetPfAk05NeutHadMulti_", &patJetPfAk05NeutHadMulti_, &b_patJetPfAk05NeutHadMulti_);
   fChain->SetBranchAddress("patJetPfAk05PhotMulti_", &patJetPfAk05PhotMulti_, &b_patJetPfAk05PhotMulti_);
   fChain->SetBranchAddress("patJetPfAk05ElecMulti_", &patJetPfAk05ElecMulti_, &b_patJetPfAk05ElecMulti_);
   fChain->SetBranchAddress("patJetPfAk05HfHadMulti_", &patJetPfAk05HfHadMulti_, &b_patJetPfAk05HfHadMulti_);
   fChain->SetBranchAddress("patJetPfAk05HfEmMulti_", &patJetPfAk05HfEmMulti_, &b_patJetPfAk05HfEmMulti_);
   fChain->SetBranchAddress("patJetPfAk05PhotEnFr_", &patJetPfAk05PhotEnFr_, &b_patJetPfAk05PhotEnFr_);
   fChain->SetBranchAddress("patJetPfAk05MuonEnFr_", &patJetPfAk05MuonEnFr_, &b_patJetPfAk05MuonEnFr_);
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
   fChain->SetBranchAddress("patJetPfAk05EmEbEn_", &patJetPfAk05EmEbEn_, &b_patJetPfAk05EmEbEn_);
   fChain->SetBranchAddress("patJetPfAk05EmEeEn_", &patJetPfAk05EmEeEn_, &b_patJetPfAk05EmEeEn_);
   fChain->SetBranchAddress("patJetPfAk05EmHfEn_", &patJetPfAk05EmHfEn_, &b_patJetPfAk05EmHfEn_);
   fChain->SetBranchAddress("patJetPfAk05HadHbEn_", &patJetPfAk05HadHbEn_, &b_patJetPfAk05HadHbEn_);
   fChain->SetBranchAddress("patJetPfAk05HadHeEn_", &patJetPfAk05HadHeEn_, &b_patJetPfAk05HadHeEn_);
   fChain->SetBranchAddress("patJetPfAk05HadHfEn_", &patJetPfAk05HadHfEn_, &b_patJetPfAk05HadHfEn_);
   fChain->SetBranchAddress("patJetPfAk05HadHoEn_", &patJetPfAk05HadHoEn_, &b_patJetPfAk05HadHoEn_);
   fChain->SetBranchAddress("patJetPfAk05TotalEm_", &patJetPfAk05TotalEm_, &b_patJetPfAk05TotalEm_);
   fChain->SetBranchAddress("patJetPfAk05TotalHad_", &patJetPfAk05TotalHad_, &b_patJetPfAk05TotalHad_);
   fChain->SetBranchAddress("patJetPfAk05EmEbFr_", &patJetPfAk05EmEbFr_, &b_patJetPfAk05EmEbFr_);
   fChain->SetBranchAddress("patJetPfAk05EmEeFr_", &patJetPfAk05EmEeFr_, &b_patJetPfAk05EmEeFr_);
   fChain->SetBranchAddress("patJetPfAk05EmHfFr_", &patJetPfAk05EmHfFr_, &b_patJetPfAk05EmHfFr_);
   fChain->SetBranchAddress("patJetPfAk05HadHbFr_", &patJetPfAk05HadHbFr_, &b_patJetPfAk05HadHbFr_);
   fChain->SetBranchAddress("patJetPfAk05HadHeFr_", &patJetPfAk05HadHeFr_, &b_patJetPfAk05HadHeFr_);
   fChain->SetBranchAddress("patJetPfAk05HadHfFr_", &patJetPfAk05HadHfFr_, &b_patJetPfAk05HadHfFr_);
   fChain->SetBranchAddress("patJetPfAk05HadHoFr_", &patJetPfAk05HadHoFr_, &b_patJetPfAk05HadHoFr_);
   fChain->SetBranchAddress("trigResults", &trigResults, &b_trigResults);
   fChain->SetBranchAddress("trigPrescale", &trigPrescale, &b_trigPrescale);
   fChain->SetBranchAddress("trigName", &trigName, &b_trigName);
   fChain->SetBranchAddress("Photonrho25", &Photonrho25, &b_Photonrho25);
   fChain->SetBranchAddress("Photonrho44", &Photonrho44, &b_Photonrho44);
   fChain->SetBranchAddress("PhotonNumPh_", &PhotonNumPh_, &b_PhotonNumPh_);
   fChain->SetBranchAddress("PhotonPt", &PhotonPt, &b_PhotonPt);
   fChain->SetBranchAddress("PhotonEta", &PhotonEta, &b_PhotonEta);
   fChain->SetBranchAddress("PhotonPhi", &PhotonPhi, &b_PhotonPhi);
   fChain->SetBranchAddress("PhotonEt", &PhotonEt, &b_PhotonEt);
   fChain->SetBranchAddress("PhotonEnergy", &PhotonEnergy, &b_PhotonEnergy);
   fChain->SetBranchAddress("PhotonPx", &PhotonPx, &b_PhotonPx);
   fChain->SetBranchAddress("PhotonPy", &PhotonPy, &b_PhotonPy);
   fChain->SetBranchAddress("PhotonPz", &PhotonPz, &b_PhotonPz);
   fChain->SetBranchAddress("PhotonR9", &PhotonR9, &b_PhotonR9);
   fChain->SetBranchAddress("PhotonPhiWidth", &PhotonPhiWidth, &b_PhotonPhiWidth);
   fChain->SetBranchAddress("PhotonEtaWidth", &PhotonEtaWidth, &b_PhotonEtaWidth);
   fChain->SetBranchAddress("PhotonScPhi", &PhotonScPhi, &b_PhotonScPhi);
   fChain->SetBranchAddress("PhotonScEta", &PhotonScEta, &b_PhotonScEta);
   fChain->SetBranchAddress("PhotonSigmaIetaIeta", &PhotonSigmaIetaIeta, &b_PhotonSigmaIetaIeta);
   fChain->SetBranchAddress("PhotonSeedTime", &PhotonSeedTime, &b_PhotonSeedTime);
   fChain->SetBranchAddress("PhotonseedSeverity", &PhotonseedSeverity, &b_PhotonseedSeverity);
   fChain->SetBranchAddress("PhotonhadronicOverEm", &PhotonhadronicOverEm, &b_PhotonhadronicOverEm);
   fChain->SetBranchAddress("PhotonecalRecHitSumEtConeDR04", &PhotonecalRecHitSumEtConeDR04, &b_PhotonecalRecHitSumEtConeDR04);
   fChain->SetBranchAddress("PhotonhcalTowerSumEtConeDR04", &PhotonhcalTowerSumEtConeDR04, &b_PhotonhcalTowerSumEtConeDR04);
   fChain->SetBranchAddress("PhotontrkSumPtHollowConeDR04", &PhotontrkSumPtHollowConeDR04, &b_PhotontrkSumPtHollowConeDR04);
   fChain->SetBranchAddress("PhotonhasPixelSeed", &PhotonhasPixelSeed, &b_PhotonhasPixelSeed);
   Notify();
}

Bool_t vector2011DataTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void vector2011DataTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t vector2011DataTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef vector2011DataTree_cxx
