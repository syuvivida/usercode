// -----------------------------------------------------
// Zee.cc 
// -- a very simple example analyzer for RECO data
// ----------------------------------------------------- 
// Shin-Shan Yu


#include <TROOT.h>
#include <TSystem.h>
#include <TObject.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/UtilAlgos/interface/TFileDirectory.h"
#include "CommonTools/Utils/interface/PtComparator.h"


#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

// this file contains the format of lepton, photon, and event structures
#include "zeeformat.hh" 
#include <map>                     

using namespace edm;
using namespace reco;
using namespace std;
using namespace math;
using namespace ROOT;
using namespace zee;

typedef std::map<reco::GsfElectronCollection::const_iterator,reco::PhotonCollection::const_iterator> elePhoMap;
typedef std::map<reco::GsfElectronCollection::const_iterator,reco::GenParticleCollection::const_iterator> eleGenMap;



class Zee : public edm::EDAnalyzer {
public:
  explicit Zee(const edm::ParameterSet&) ;
  ~Zee();  

    
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void dumpGenInfo(const edm::Event&); 
  virtual void analyze(const edm::Event&, const edm::EventSetup&) ;
  virtual void endJob() ;
  bool isLoosePhoton(const reco::Photon& it_ph);
  void MatchElectronToPhoton(const edm::Event&);
  void MatchElectronToGenp(const edm::Event&);
  

  TTree *root;
  EvtInfoBranches  EvtInfo;
  PhoInfoBranches  PhoInfo;
  GenInfoBranches  GenInfo;
  elePhoMap        myElePhoMap;
  eleGenMap        myEleGenMap;

  bool             _dumpHEP;
  int              _nIn;
  int              _nOut;


  // histogram
  TH1F* h_trkisodeno;
  TH1F* h_trkisonumr;
  TH1F* h_sigietadeno;
  TH1F* h_sigietanumr;
  TH2F* h_xy;
  TH2F* h_rz;
  TH1F* h_dE;



};


Zee::Zee(const edm::ParameterSet& iConfig):
  _nIn(0), _nOut(0)

{  
  _dumpHEP = iConfig.getUntrackedParameter<bool>("dumpHEP", false);
 
}


Zee::~Zee()
{
}

void Zee::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  TFileDirectory results = TFileDirectory( fs->mkdir("Zee") );


  h_trkisodeno = fs->make<TH1F>("h_trkisodeno","Track isolation denominator",
				200,0,100);
  h_trkisonumr = fs->make<TH1F>("h_trkisonumr","Track isolation numerator",
				200,0,100);
  h_sigietadeno = fs->make<TH1F>("h_sigietadeno","Sigma ieta ieta denominator",
				 200,0,1);
  h_sigietanumr = fs->make<TH1F>("h_sigietanumr","Sigma ieta ieta numerator",
				 200,0,1);
  
  h_xy          = fs->make<TH2F>("h_xy","",100,-10.0,10.0,100,-10.0,10.0);
  h_rz          = fs->make<TH2F>("h_rz","",100,0,10.0,100,-10.0,10.0);
  h_dE          = fs->make<TH1F>("h_dE","",100,0,2);


  root = new TTree("root","root");
  EvtInfo.Register(root);  
  PhoInfo.Register(root);
  GenInfo.Register(root);
  EvtInfo.Initialize();  
  PhoInfo.Initialize();
  GenInfo.Initialize();
}

void Zee::endJob() 
{
  std::cout << "Zee has " << _nIn << " input events and " << _nOut  << " events" << endl;
}

// dump generator-level information

void Zee::dumpGenInfo(const edm::Event& iEvent)
{
  edm::Handle<reco::GenParticleCollection> GenHandle;  
  bool hasGenParticle = iEvent.getByLabel("genParticles", GenHandle);
  if(!hasGenParticle)return;

  std::cout << "============================ " << 
    "Run " <<  iEvent.id().run() << " Evt " <<  iEvent.id().event() <<
    " ============================= " << std::endl << std::endl;

  int genIndex = 0;

  printf("GenIndex  ");
  printf("PDG  ");
  printf("Status ");
  printf("Mother PID");
  printf("   ");
  printf("Mass ");
  printf("Energy ");
  printf("Pt ");
  printf("Eta ");
  printf("\n");


  for( GenParticleCollection::const_iterator it_gen = GenHandle->begin(); 
       it_gen != GenHandle->end(); it_gen++ ) {

    printf("%4i", genIndex);
    printf("%7i", it_gen->pdgId());
    printf("%6i", it_gen->status());
    if(it_gen->mother())
      printf("%7i", it_gen->mother()->pdgId());
    else
      printf("%7i", -1);
    printf("   ");
    printf("%9.3f",it_gen->p4().mass());
    printf("%9.3f",it_gen->energy());
    printf("%9.3f",it_gen->pt());
    printf("%9.3f",it_gen->eta());
    printf("\n");
    genIndex ++;
  }

  std::cout << std::endl 
	    << "==========================================================="
            << "===============" 
	    << std::endl;

}


// analyzing reconstructed electrons, muons, and photons

void Zee::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  _nIn++;
  
   if(_dumpHEP)
    dumpGenInfo(iEvent);

   MatchElectronToPhoton(iEvent);
   MatchElectronToGenp(iEvent);

  
  // Start to fill the main root branches
  // initialize the variables of ntuples
  EvtInfo.Initialize();  
  EvtInfo.RunNo  = iEvent.id().run();
  EvtInfo.EvtNo  = iEvent.id().event();


  // fill generator photons information
  edm::Handle<reco::GenParticleCollection> GenHandle;  
  bool hasGenParticle = iEvent.getByLabel("genParticles", GenHandle); 
  GenInfo.Initialize();
  if(hasGenParticle){
    for( std::vector<GenParticle>::const_iterator it_gen = 
	   GenHandle->begin(); 
	 it_gen != GenHandle->end() && GenInfo.Size < 30; it_gen++ ) {
   
      if(it_gen->pdgId()==23 && GenInfo.Size < 10)EvtInfo.GenZMass = it_gen->mass();
      GenInfo.PID[GenInfo.Size] = it_gen->pdgId();
      GenInfo.MPID[GenInfo.Size] = it_gen->mother()?
	it_gen->mother()->pdgId():-999;
      GenInfo.Mass[GenInfo.Size] = it_gen->mass();
      GenInfo.Pt[GenInfo.Size] = it_gen->pt();
      GenInfo.Eta[GenInfo.Size] = it_gen->eta();
      GenInfo.Phi[GenInfo.Size] = it_gen->phi();
      GenInfo.Size ++;

    } // check generator-level
  } // end of if hasGenParticle

  
  for(eleGenMap::iterator it_e=myEleGenMap.begin(); 
      it_e!= myEleGenMap.end(); ++it_e)
    {
      // find match to photon
      if(myElePhoMap.find(it_e->first)==myElePhoMap.end())continue;

      for(eleGenMap::const_iterator jt_e = myEleGenMap.begin();
	  jt_e != myEleGenMap.end(); ++jt_e){
	if(jt_e->first == it_e->first)continue;

      // find match to photon
      if(myElePhoMap.find(jt_e->first)==myElePhoMap.end())continue;

      EvtInfo.RecZMass = (it_e->first->p4() 
			  + jt_e->first->p4()).M();


      reco::PhotonCollection::const_iterator phoMatch[2];
      phoMatch[0] = myElePhoMap[it_e->first];
      phoMatch[1] = myElePhoMap[jt_e->first];

    // initialize the variables of ntuples
      PhoInfo.Initialize();

      for (unsigned int ey=0; ey < 2; ey++){
    
	reco::PhotonCollection::const_iterator it_ph = phoMatch[ey];
    

	float ecalIso = it_ph->ecalRecHitSumEtConeDR04();
	float hcalIso = it_ph->hcalTowerSumEtConeDR04();
	float trkIso  = it_ph->trkSumPtHollowConeDR04();
	
	int location = -1; 
	bool isBarrel = it_ph->isEB();
	bool isEndCap = it_ph->isEE();
	bool inAnyGap = it_ph->isEBEEGap() || 
	  (it_ph->isEB() && it_ph->isEBGap()) || 
	  (it_ph->isEE() && it_ph->isEEGap());

	if(inAnyGap)location=0;
	else if(isBarrel)location=1;
	else if(isEndCap)location=2;

	PhoInfo.Location[ey]     = location;
	PhoInfo.E    [ey] 	 = it_ph->energy(); 
	PhoInfo.Et   [ey] 	 = it_ph->et();     
	PhoInfo.Pz   [ey] 	 = it_ph->pz();     
	PhoInfo.Eta  [ey] 	 = it_ph->eta();    
	PhoInfo.Phi  [ey] 	 = it_ph->phi();    
	PhoInfo.R9   [ey] 	 = it_ph->r9();     
	PhoInfo.TrkPtSum[ey]     = trkIso;
	PhoInfo.EcalRecHitEtSum[ey]= ecalIso;
	PhoInfo.HcalTowerEtSum[ey] = hcalIso;
	PhoInfo.HoverE[ey]       = it_ph->hadronicOverEm();

	PhoInfo.SCE[ey]          = it_ph->superCluster()->energy();
	PhoInfo.SCEta[ey]        = it_ph->superCluster()->eta();
	PhoInfo.SCPhi[ey]        = it_ph->superCluster()->phi();
	PhoInfo.SCEtaWidth[ey]   = it_ph->superCluster()->etaWidth();
	PhoInfo.SCPhiWidth[ey]   = it_ph->superCluster()->phiWidth();
	
	float scet = it_ph->superCluster()->energy()/cosh(it_ph->superCluster()->eta());
	PhoInfo.SCEt[ey]       = scet;

	PhoInfo.SCNCrystal[ey] = it_ph->superCluster()->size();

	PhoInfo.SigEta[ey]       = it_ph->sigmaEtaEta();
	PhoInfo.SigIEta[ey]      = it_ph->sigmaIetaIeta();
      
	bool passOffline = isLoosePhoton(*it_ph);
	PhoInfo.IsLoose[ey] = passOffline? 1: 0;
 

      } // if number of photons > 0
    
      root->Fill();

      } // end of loop over jt_e
  } // end of loop over it_e

   
  _nOut++;
  
}


void Zee::MatchElectronToPhoton(const edm::Event& iEvent)
{
  myElePhoMap.clear();

  // look for electron collection
  Handle<reco::GsfElectronCollection> electronColl;
  iEvent.getByLabel("gsfElectrons", electronColl);
  bool eleIsValid = electronColl.isValid();
  if(!eleIsValid){cout << "Electron is not valid" << endl; return;}

  // look for photon collection
  Handle<reco::PhotonCollection> photonColl;
  iEvent.getByLabel("photons", photonColl);
  bool phoIsValid = photonColl.isValid();
  if(!phoIsValid){cout << "Photon is not valid" << endl; return;}
  

  // first loop over electrons
  unsigned int count=0;
  for(reco::GsfElectronCollection::const_iterator it_e = electronColl->begin();
      it_e != electronColl->end(); it_e++){

    if(it_e->pt()>10)count ++ ;
    double ptMax   = 0;
    bool findMatch = false;

    // now loop over photons
    reco::PhotonCollection::const_iterator thisPho;
    for (reco::PhotonCollection::const_iterator it_ph = photonColl->begin(); 
	 it_ph!=photonColl->end(); it_ph++){
  
      // find the highest pt photon that has the same supercluster seed as 
      // electron
      double thisPhoEnergy = it_ph->superCluster()->energy();
      double distx = fabs(it_e->superCluster()->seed()->x() -
			  it_ph->superCluster()->seed()->x());
      double disty = fabs(it_e->superCluster()->seed()->y() -
			  it_ph->superCluster()->seed()->y());
      double distrphi = sqrt(distx*distx + disty*disty);

      double distz = fabs(it_e->superCluster()->seed()->z() -
			  it_ph->superCluster()->seed()->z());

      h_xy->Fill(distx,disty);
      h_rz->Fill(distrphi,distz);
   
      double relE  = it_e->superCluster()->seed()->energy()>1e-3? 
	fabs(it_ph->superCluster()->seed()->energy() - 
	     it_e->superCluster()->seed()->energy())
	/it_e->superCluster()->seed()->energy(): -999;
      
      h_dE->Fill(relE);

      if( distx < 3.0 && disty < 3.0 && distz < 3.0 && relE < 0.2	 
	 && thisPhoEnergy > ptMax)
	{
	  ptMax   = thisPhoEnergy;
	  thisPho = it_ph;
	  findMatch = true;
	} // if matched
//       else
// 	{
// 	  cout << "photon seed energy = " << 
// 	    it_ph->superCluster()->seed()->energy() << "\t" 
// 	       << "electron seed energy = " << 
// 	    it_e->superCluster()->seed()->energy() << endl;
// 	  cout << "photon seed position = (" << 
// 	    it_ph->superCluster()->seed()->x() << ", " << 
// 	    it_ph->superCluster()->seed()->y() << ", " << 
// 	    it_ph->superCluster()->seed()->z() << ") \t " 
// 	       << "electron seed position = (" << 
// 	    it_e->superCluster()->seed()->x() << ", " << 
// 	    it_e->superCluster()->seed()->y() << ", " << 
// 	    it_e->superCluster()->seed()->z() << ") \t " << endl;
// 	}
      
    } // end of loop over photon

    if(findMatch)
      myElePhoMap.insert(std::pair<reco::GsfElectronCollection::const_iterator,
			 reco::PhotonCollection::const_iterator>(it_e,thisPho));
    
  
  } // end of loop over electron

  if(myElePhoMap.size()!=count)
    cout << "There are " << count  << " electrons" << endl;

  return;

}


void Zee::MatchElectronToGenp(const edm::Event& iEvent)
{

  myEleGenMap.clear();

  // look for electron collection
  Handle<reco::GsfElectronCollection> electronColl;
  iEvent.getByLabel("gsfElectrons", electronColl);
  bool eleIsValid = electronColl.isValid();
  if(!eleIsValid){cout << "Electron is not valid" << endl; return;}

  // look for Gen particle collection
  edm::Handle<reco::GenParticleCollection> GenHandle;  
  bool hasGenParticle = iEvent.getByLabel("genParticles", GenHandle);
  if(!hasGenParticle)return;


  // loop over electron first
  for(reco::GsfElectronCollection::const_iterator it_e = electronColl->begin();
      it_e != electronColl->end(); it_e++){

    double ptMax   = 0;
    bool findMatch = false;

    // find the highest et gen particle from Z matched to electron 
    GenParticleCollection::const_iterator thisGen;
    for(GenParticleCollection::const_iterator it_gen = GenHandle->begin(); 
	it_gen != GenHandle->end(); it_gen++ ) {

      double thisGenEnergy = it_gen->pt();

      if(abs(it_gen->pdgId())!=11 || it_gen->status()!=1)continue;
      if(it_gen->mother() && abs(it_gen->mother()->pdgId())!=11)continue;
      if(it_gen->mother()->mother() 
	 && abs(it_gen->mother()->mother()->pdgId()!=23))continue;
      
      float dR = reco::deltaR(it_e->momentum(), it_gen->momentum());
      float relPt = fabs(it_e->pt()-it_gen->pt())/it_gen->pt();

      if(dR<0.5 && relPt < 1.0 && thisGenEnergy > ptMax)
	{
	  ptMax = thisGenEnergy;
	  thisGen = it_gen;
	  findMatch = true;
	} // if find a match
    } // end of loop over gen particle

    if(findMatch)
      myEleGenMap.insert(std::pair<reco::GsfElectronCollection::const_iterator,
			 reco::GenParticleCollection::const_iterator>
			 (it_e,thisGen));

  } // end of loop over electron
  
  return;

}


bool Zee::isLoosePhoton(const reco::Photon& it_ph)
{
  float et = it_ph.et();
  float eta = it_ph.eta();

  float ecalIso = it_ph.ecalRecHitSumEtConeDR04();
  float hcalIso = it_ph.hcalTowerSumEtConeDR04();
  float trkIso  = it_ph.trkSumPtHollowConeDR04();

  bool isBarrel = it_ph.isEB();
  bool isEndCap = it_ph.isEE();
  bool inAnyGap = it_ph.isEBEEGap() || 
    (it_ph.isEB() && it_ph.isEBGap()) || 
    (it_ph.isEE() && it_ph.isEEGap());

  if(inAnyGap)return false;

  // Barrel cuts
  if(isBarrel && ecalIso  > 5.0 + 0.0045*et)return false;
  if(isBarrel && hcalIso  > 5.0)return false;
  if(isBarrel && trkIso > 9.0 )return false;
  if(isBarrel && it_ph.hadronicOverEm() > 0.15)return false;

  // Endcap cuts
  if(isEndCap && ecalIso  > 5.0 + 0.02*et)return false;
  if(isEndCap && hcalIso  > 7.0)return false;
  if(isEndCap && trkIso > 9.0 )return false;
  if(isEndCap && it_ph.hadronicOverEm() > 0.15)return false;

  if(et        < 10.0)return false;
  if(fabs(eta) > 1.44 && fabs(eta)<1.56)return false;
  if(fabs(eta) > 2.5)return false;

  return true;

}

//define this as a plug-in
DEFINE_FWK_MODULE(Zee);


