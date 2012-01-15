// -----------------------------------------------------
// GJetAngular.cc 
// -- a very simple example analyzer for RECO data
// ----------------------------------------------------- 
// Shin-Shan Yu

#include "syu/MyObjectCounter/header/GJetAngular.hh" 
#include "DataFormats/Math/interface/deltaPhi.h"
#include <algorithm>


GJetAngular::GJetAngular(const edm::ParameterSet& iConfig):
  _nIn(0), _nOut(0), _alg(iConfig),
  _etaMax(iConfig.getUntrackedParameter<double>("etaMax",9999.0)),
  _ptMin(iConfig.getUntrackedParameter<double>("ptMin",0.0))

{  
}


GJetAngular::~GJetAngular()
{
}

void GJetAngular::beginJob()
{
  TFileDirectory results = TFileDirectory( fs->mkdir("GJetAngular") );

  fPtHat = fs->make<TH1F>("fPtHat","Event pthat", 200,0.0,200.0);
  fPtPho = fs->make<TH1F>("fPtPho","Event pt of stable photon", 200,0.0,100.0);
  fPtHard = fs->make<TH1F>("fPtHard","Event pt of hard scattering photon", 200,0.0,100.0);
  fPStar = fs->make<TH1F>("fPStar","Event momentum", 400,0.0,400.0);

  fCosThetaStar=fs->make<TH1F>("fCosThetaStar","cos#theta^{*}",100,0.0,1.0);

  fTanH=fs->make<TH1F>("fTanH","cos#theta^{*}",100,0.0,1.0);
  
  fPtHatVsCosThetaStar=fs->make<TH2F>("fPtHatVsCosThetaStar","pthat vs cos#theta^{*}",100,0.0,1.0,
				      200,0.0,200.0);
  fPStarVsCosThetaStar=fs->make<TH2F>("fPStarVsCosThetaStar","pstar vs cos#theta^{*}",100,0.0,1.0,
				      400,0.0,400.0);
  pfPStarVsCosThetaStar=fs->make<TProfile>("pfPStarVsCosThetaStar","pstar vs cos#theta^{*}",100,0.0,1.0,
                                      0.0,400.0);
  pfPtHatVsCosThetaStar=fs->make<TProfile>("pfPtHatVsCosThetaStar","pthat vs cos#theta^{*}",100,0.0,1.0,
                                      0.0,400.0);


  std::cout << "Focusing on photon with etaMax = " << 
    _etaMax << " and ptMin = " << _ptMin << std::endl;
}

void GJetAngular::endJob() 
{

  std::cout << "GJetAngular has " << _nIn << " input events and " << _nOut  << " output events" << endl;
}


// analyzing reconstructed electrons, muons, and photons

// void GJetAngular::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
bool GJetAngular::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // initialization
  bool accept = false;
  _nIn++;
  _alg.init(iEvent, false, false, false); 
  float ptHat = _alg.ptHat(); 
  fPtHat->Fill(ptHat);


  GenParticleCollection::const_iterator gamma_hard;
  GenParticleCollection::const_iterator parton_hard;
  bool hasOut1=false;
  bool hasOut2=false;

  // loop over GenParticles
  edm::Handle<reco::GenParticleCollection> GenHandle = _alg.getGenHandle();
  if(GenHandle.isValid()){

    int Index = -1;
    for( GenParticleCollection::const_iterator it_gen = 
	   GenHandle->begin(); 
	 it_gen != GenHandle->end(); it_gen++ ) {

      Index ++;
      
      // look for hard scattering parton and photon
      if( (Index == 6 || Index==7) && it_gen->status()==3 
	  && it_gen->pdgId()==22 && it_gen->pt() > _ptMin
	  && fabs(it_gen->eta())<_etaMax){
	gamma_hard = it_gen;
	hasOut1   = true;
      }
      else if((Index == 6 || Index==7) && it_gen->status()==3
	      && it_gen->pdgId()!=22){
	parton_hard = it_gen;
	hasOut2   = true;
      }
    } // loop over gen particle

  }
  
  if(!hasOut1 || !hasOut2){accept=false; return accept;}

  accept = true;

  reco::Candidate::LorentzVector p4_pho = gamma_hard->p4();
  reco::Candidate::LorentzVector p4_parton = parton_hard->p4();
  reco::Candidate::LorentzVector p4_all = p4_pho+p4_parton;

  reco::Candidate::Vector boostVector= p4_all.BoostToCM();

  fPtPho->Fill(p4_pho.pt());
  fPtHard->Fill(p4_parton.pt());
  //  cout << "before boost photon pt = " << p4_pho.pt() << " and jet pt = "
  // << p4_parton.pt() << endl;

  //  cout << "before boost photon ( " << p4_pho.px() << ", " << p4_pho.py() << ", " <<
  //    p4_pho.pz() << " )" << endl;
  //  cout << "before boost parton ( " << p4_parton.px() << ", " << p4_parton.py() << ", " <<
  //    p4_parton.pz() << " )" << endl;
  
  //  ROOT::Math::Boost boost(boostVector.x(), boostVector.y(), boostVector.z());
  ROOT::Math::Boost boost(boostVector);
  //  reco::Candidate::LorentzVector p4_all_afterboost(boost(p4_all));
  reco::Candidate::LorentzVector p4_pho_afterboost(boost(p4_pho));
  reco::Candidate::LorentzVector p4_parton_afterboost(boost(p4_parton));


  double ptCM  = p4_pho_afterboost.pt();
  double pstar = p4_pho_afterboost.P();
  double cosThetaStar = fabs(cos(p4_pho_afterboost.Theta()));
  double tanH = tanh(fabs(p4_pho.Rapidity()-p4_parton.Rapidity())*0.5);

  fPStar->Fill(pstar);
  fCosThetaStar->Fill(cosThetaStar);
  fTanH->Fill(tanH);
  fPtHatVsCosThetaStar->Fill(cosThetaStar, ptHat);
  fPStarVsCosThetaStar->Fill(cosThetaStar, pstar);
  pfPStarVsCosThetaStar->Fill(cosThetaStar, pstar);
  pfPtHatVsCosThetaStar->Fill(cosThetaStar, ptCM);

  //  cout << "after boost photon pt = " << p4_pho_afterboost.pt() << " and jet pt = "
  //    << p4_parton_afterboost.pt() << endl;  


  //  cout << "after boost photon ( " << p4_pho_afterboost.px() << ", " << p4_pho_afterboost.py() << ", " <<
  //    p4_pho_afterboost.pz() << " )" << endl;
  //  cout << "before boost parton ( " << p4_parton_afterboost.px() << ", " << p4_parton_afterboost.py() << ", " <<p4_parton_afterboost.pz() << " )" << endl;


    
  
  if(accept){_nOut++;}
  return accept;
  
}



//define this as a plug-in
DEFINE_FWK_MODULE(GJetAngular);


