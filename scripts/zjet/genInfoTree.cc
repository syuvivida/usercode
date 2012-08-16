#include "DelPanj/TreeMaker/interface/genInfoTree.h"
#include <TLorentzVector.h>

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
    //                             FSR, 2012/04/25, SSY
    
    int pid = abs(gen.pdgId());
    if( gen.status()!=3 
	&& pid!=11 
	&& pid!=13 
	&& pid!=22
	&& pid!=23
	&& pid!=24
	)
      continue;

    // for dressing leptons
    
    if( (pid==11 || pid==13) && gen.status()==1 && gen.pt() > 10.0 
	&& fabs(gen.eta()) < 3.0
	){
      
      // check if this lepton is really from Z decays directly
      
      // for Madgraph MC
      if(gen.numberOfMothers() ==1 && 
	 gen.mother()->pdgId() != gen.pdgId())continue;

      // for Sherpa MC
      if(gen.numberOfMothers() ==2 && 
	 (abs(gen.mother(0)->pdgId()) != abs(gen.pdgId()) ||
	  abs(gen.mother(1)->pdgId()) != abs(gen.pdgId()))
	 )continue;
      
      TLorentzVector genLep(0,0,0,0);
      genLep.SetPtEtaPhiE(gen.pt(),
			  gen.eta(),
			  gen.phi(),
			  gen.energy());

      // now check photons nearby

      TLorentzVector genPho(0,0,0,0);
      const reco::GenParticleCollection* genColl2= &(*genParticleHandle);
      reco::GenParticleCollection::const_iterator geni2 = genColl2->begin();
      

      for(; geni2!=genColl2->end();geni2++){
	reco::GenParticle gen2 = *geni2;

	if(gen2.pdgId()!=22)continue;
	if(gen2.energy()<1e-6)continue;

	
	TLorentzVector thisPho(0,0,0,0);	
	thisPho.SetPtEtaPhiE(gen2.pt(),
			     gen2.eta(),
			     gen2.phi(),
			     gen2.energy());

	double dR = genLep.DeltaR(thisPho);

	if(dR > 0.1)continue;

 	if(gen2.numberOfMothers() !=1 && 
	   gen2.numberOfMothers() !=2)continue;

	// for Madgraph MC
	if(gen2.numberOfMothers() ==1 && 
	   gen2.mother()->pdgId() != gen.pdgId())continue;

	// for Sherpa MC
	if(gen2.numberOfMothers() ==2 && 
	   (abs(gen2.mother(0)->pdgId()) != abs(gen.pdgId()) ||
	    abs(gen2.mother(1)->pdgId()) != abs(gen.pdgId()))
	   )continue;
	
// 	std::cout << "lepton index: " << genIndex << "\t lepton Id: " << gen.pdgId() << "\t photon numMother: " << gen2.numberOfMothers()  << "\t";

// 	if(gen2.numberOfMothers() ==1)std::cout << "Mother PDGID = " << 
// 	  gen2.mother()->pdgId() << std::endl;
// 	if(gen2.numberOfMothers() ==2)std::cout << "Mother PDGIDs = " << 
// 	  gen2.mother(0)->pdgId() << "\t" <<
// 	  gen2.mother(1)->pdgId() << std::endl;

	
	genPho += thisPho;
	
	
      } // end of loop over gen particles

      genLepE_.push_back(gen.energy());
      genLepPt_.push_back(gen.pt());
      genLepEta_.push_back(gen.eta());
      genLepPhi_.push_back(gen.phi());
      genLepId_.push_back(gen.pdgId());

      
      if(genPho.E()>1e-6){
	genLepPhoE_.push_back  (genPho.E());
	genLepPhoPt_.push_back (genPho.Pt());
	genLepPhoEta_.push_back(genPho.Eta());
	genLepPhoPhi_.push_back(genPho.Phi());
      }

      else
	{
	  genLepPhoE_.push_back  (1e-6);
	  genLepPhoPt_.push_back (1e-6);
	  genLepPhoEta_.push_back(1e-6);
	  genLepPhoPhi_.push_back(1e-6);
	}

    }


    if( pid==22 ) continue;

    genParE_.push_back(gen.energy());
    genParPt_.push_back(gen.pt());
    genParEta_.push_back(gen.eta());
    genParPhi_.push_back(gen.phi());
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
  AddBranch(&genParQ_,"genParQ_");
  AddBranch(&genParId_,"genParId_");
  AddBranch(&genParSt_,"genParSt_");
  AddBranch(&genMomParId_,"genMomParId_");
  AddBranch(&genParIndex_,"genParIndex_");

  AddBranch(&genLepE_, "genLepE_");
  AddBranch(&genLepPt_, "genLepPt_");
  AddBranch(&genLepEta_,"genLepEta_");
  AddBranch(&genLepPhi_,"genLepPhi_");

  AddBranch(&genLepPhoE_, "genLepPhoE_");
  AddBranch(&genLepPhoPt_, "genLepPhoPt_");
  AddBranch(&genLepPhoEta_,"genLepPhoEta_");
  AddBranch(&genLepPhoPhi_,"genLepPhoPhi_");

  AddBranch(&genLepId_,"genLepId_");

  
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
  genParQ_.clear();
  genParId_.clear();
  genParSt_.clear();

  genLepE_.clear();
  genLepPt_.clear();
  genLepEta_.clear();
  genLepPhi_.clear();

  genLepPhoE_.clear();
  genLepPhoPt_.clear();
  genLepPhoEta_.clear();
  genLepPhoPhi_.clear();

  genLepId_.clear();

  genMomParId_.clear();
  genParIndex_.clear();
  genJetE_.clear();
  genJetPt_.clear();
  genJetEta_.clear(); 
  genJetPhi_.clear(); 


}



