/*
  Anil Singh
  Panjab University
*/

#include "DelPanj/TreeMaker/interface/genInfoTree.h"

//---------------------------------------------------------------
// Add Branches to the genTree
//---------------------------------------------------------------
void
genInfoTree::AddBranch(double* x, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),x,(brName+"/D").c_str());
}


//---------------------------------------------------
//---------------------------------------------------
void
genInfoTree::AddBranch(std::vector<double>* vec, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);
}

//---------------------------------------------------------------
//---------------------------------------------------------------
void
genInfoTree::AddBranch(std::vector<int>* vec, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);
}

genInfoTree::genInfoTree(std::string name, TTree* tree, const edm::ParameterSet& iConfig)
{
  tree_=tree; 
  genPartLabel_ = iConfig.getParameter<edm::InputTag>("genPartLabel");
  SetBranches();
}


genInfoTree::~genInfoTree()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  delete tree_;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
genInfoTree::Fill(const edm::Event& iEvent)
{
  Clear();
  using namespace edm;
  edm::Handle<reco::GenParticleCollection> genParticleHandle;
  if(not iEvent.getByLabel(genPartLabel_, genParticleHandle))
    {
      std::cout<<
	"GenAnalyzer: Generator Level Information not found\n"
	       <<std::endl;
    }

  edm::Handle<GenEventInfoProduct>    genEventScale;

  if (iEvent.getByLabel("generator", genEventScale)) {
    if (genEventScale->hasBinningValues())
      ptHat_ = genEventScale->binningValues()[0];
      
    mcWeight_ = genEventScale->weight();

  }

  unsigned int genIndex=0;
  const reco::GenParticleCollection* genColl= &(*genParticleHandle);
  reco::GenParticleCollection::const_iterator geni = genColl->begin();
  for(; geni!=genColl->end() && genIndex < 200;geni++){
    reco::GenParticle gen = *geni;
    
    genIndex++;
    //Look out for the GenMuons/Electrons
//     if(gen.status()!=1)continue; // remove this condition so to save electrons/muons before 
//                             FSR, 2012/04/25, SSY
    
    double pid = fabs(gen.pdgId());
    if(!( (gen.status()==3) || (pid==11)||(pid==13) || (pid==23) || (pid==24)))continue;
    genParE_.push_back(gen.energy());
    genParPt_.push_back(gen.pt());
    genParEta_.push_back(gen.eta());
    genParPhi_.push_back(gen.phi());
    genParM_.push_back(gen.mass());
    genParQ_.push_back(gen.charge());
    genParId_.push_back(gen.pdgId());
    genParSt_.push_back(gen.status());

    int mompid = -9999;
    if( gen.numberOfMothers() ==1 ) 
      mompid = gen.mother()->pdgId();
    else
      mompid = 10000+gen.numberOfMothers();

    genMomParId_.push_back(mompid);

    genParIndex_.push_back(genIndex);
      
  }

  edm::Handle<reco::GenJetCollection> genJetsHandle;
  if( not iEvent.getByLabel("ak5GenJets",genJetsHandle)){ 
    edm::LogInfo("GenAnalyzer") << "genJets not found, "
      "skipping event"; 
    return;
  }
  const reco::GenJetCollection* genJetColl = &(*genJetsHandle);
  reco::GenJetCollection::const_iterator gjeti = genJetColl->begin();
   
  for(; gjeti!=genJetColl->end();gjeti++){
    reco::GenParticle gjet = *gjeti;
    if(gjet.pt()<=20)continue;
    if(fabs(gjet.eta())>3.0)continue;

    genJetE_.push_back(gjet.energy()); 
    genJetPt_.push_back(gjet.pt());
    genJetEta_.push_back(gjet.eta());
    genJetPhi_.push_back(gjet.phi());
  }
   
}



void  
genInfoTree::SetBranches(){
  AddBranch(&ptHat_, "ptHat_");
  AddBranch(&mcWeight_, "mcWeight_");
  AddBranch(&genParE_, "genParE_");
  AddBranch(&genParPt_, "genParPt_");
  AddBranch(&genParEta_,"genParEta_");
  AddBranch(&genParPhi_,"genParPhi_");
  AddBranch(&genParM_,"genParM_");
  AddBranch(&genParQ_,"genParQ_");
  AddBranch(&genParId_,"genParId_");
  AddBranch(&genParSt_,"genParSt_");
  AddBranch(&genMomParId_,"genMomParId_");
  AddBranch(&genParIndex_,"genParIndex_");
  
  AddBranch(&genJetE_, "genJetE_");
  AddBranch(&genJetPt_,"genJetPt_");
  AddBranch(&genJetEta_,"genJetEta_");
  AddBranch(&genJetPhi_,"genJetPhi_");

}


void  
genInfoTree::Clear(){

  ptHat_ = -9999.0;
  mcWeight_ = -9999.0; 

  genParE_.clear();
  genParPt_.clear();
  genParEta_.clear();
  genParPhi_.clear();
  genParM_.clear();
  genParQ_.clear();
  genParId_.clear();
  genParSt_.clear();
  genMomParId_.clear();
  genParIndex_.clear();
  genJetE_.clear();
  genJetPt_.clear();
  genJetEta_.clear(); 
  genJetPhi_.clear(); 
}



