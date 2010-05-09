// -----------------------------------------------------
// QCDGenIso.cc 
// -- a very simple example analyzer for RECO data
// ----------------------------------------------------- 
// Shin-Shan Yu

#include "syu/MyObjectCounter/header/QCDGenIso.hh" 
#include "DataFormats/Math/interface/deltaPhi.h"
#include <algorithm>


QCDGenIso::QCDGenIso(const edm::ParameterSet& iConfig):
  _nIn(0), _nOut(0), _alg(iConfig),
  _pdgCode(iConfig.getUntrackedParameter<int>("pdgCode",22)),
  _etaMax(iConfig.getUntrackedParameter<double>("etaMax",2.5)),
  _ptMin(iConfig.getUntrackedParameter<double>("ptMin",15.0)),
  _isoMaxCut(iConfig.getUntrackedParameter<double>("isoMax",2.0))
{  
}


QCDGenIso::~QCDGenIso()
{
}

void QCDGenIso::beginJob()
{
  TFileDirectory results = TFileDirectory( fs->mkdir("QCDGenIso") );

  fPtHat = fs->make<TH1F>("fPtHat","Event pthat", 200,0.0,200.0);
  fPtHard = fs->make<TH1F>("fPtHard","Event pt of hard scattering partons", 200,0.0,100.0);
  fPtStablePho = fs->make<TH1F>("fPtStablePho","Event pt of stable photon", 200,0.0,100.0);
  fPtBremPho = fs->make<TH1F>("fPtBremPho","Event pt of brem photons", 200,0.0,100.0);
  fPtFragPho = fs->make<TH1F>("fPtFragPho","Event pt of fragmentation photons", 200,0.0,100.0);
  fPtDecayPho = fs->make<TH1F>("fPtDecayPho","Event pt of decay photons", 200,0.0,100.0);

  histIsoDR04_all = fs->make<TH1F>("histIsoDR04_all","Isolation for all photons",2000,0.0,200.0) ;
  histIsoDR04_brem = fs->make<TH1F>("histIsoDR04_brem","Isolation for brem photons",2000,0.0,200.0) ;
  histIsoDR04_frag = fs->make<TH1F>("histIsoDR04_frag","Isolation for fragmentation photons",2000,0.0,200.0) ;
  histIsoDR04_decay = fs->make<TH1F>("histIsoDR04_decay","Isolation for decay photons",2000,0.0,200.0) ;


  // isolation histograms
  const int nptbin   = 50;
  const float xptmin = 0.0;
  const float xptmax = 200.0;


  fTrkIsoPt00 = fs->make<TH1F>("fTrkIsoPt00", "Photon pt before track"
			       " isolation cuts", nptbin, xptmin, xptmax);
  fTrkIsoPt01 = fs->make<TH1F>("fTrkIsoPt01", "Photon pt after track isolation"
			       " pt < 2 GeV/c cuts", nptbin, xptmin, xptmax);
  fEffTrkIso  = fs->make<TH1F>("fEffTrkIso", "Photon track isolation "
			       "efficiency", nptbin, xptmin, xptmax);
  fTrkIsoPt00->Sumw2();
  fTrkIsoPt01->Sumw2();
  fEffTrkIso->Sumw2();



  fIsoPt00 = fs->make<TH1F>("fIsoPt00", "Photon pt before isolation cuts", nptbin, xptmin, xptmax);  
  fIsoPt01 = fs->make<TH1F>("fIsoPt01", "Photon pt after isolation Et < 2 GeV cuts", nptbin, xptmin, xptmax);
  fEffIso  = fs->make<TH1F>("fEffIso", "Photon isolation efficiency", nptbin, xptmin, xptmax);
  fIsoPt00->Sumw2();
  fIsoPt01->Sumw2();
  fEffIso->Sumw2();


  fIsoCone0300 = fs->make<TH1F>("fIsoCone0300", "Photon pt before isolation cuts", nptbin, xptmin, xptmax);  
  fIsoCone0300->Sumw2();
  fIsoCone0301 = fs->make<TH1F>("fIsoCone0301", "Photon pt after isolation Et < 2 GeV (cone 0.3) cuts", nptbin, xptmin, xptmax);
  fIsoCone0301->Sumw2();
  fEffIsoCone03  = fs->make<TH1F>("fEffIsoCone03", "Photon isolation (cone 0.3) efficiency", nptbin, xptmin, xptmax);
  fEffIsoCone03->Sumw2();
 
//   TGraphAsymmErrors* fEffIsoCone03 = fs->make<TGraphAsymmErrors>(fIsoCone0300->GetNbinsX());


  // isolation Et < 2 GeV vs eta
  const int netabin = 40;
  const float xetamin = -4.0;
  const float xetamax =  4.0;

  fIsoEta00 = fs->make<TH1F>("fIsoEta00", "Photon Eta before isolation cuts", netabin, xetamin, xetamax);
  fIsoEta01 = fs->make<TH1F>("fIsoEta01", "Photon Eta after isolation Et < 2 GeV cuts", netabin, xetamin, xetamax);
  fEffIsoEta = fs->make<TH1F>("fEffIsoEta", "Photon isolation efficiency", netabin, xetamin, xetamax);
  fIsoEta00->Sumw2();
  fIsoEta01->Sumw2();
  fEffIsoEta->Sumw2();


  std::cout << "Focusing on PDG code = " << _pdgCode << " with etaMax = " << 
    _etaMax << " and ptMin = " << _ptMin << std::endl;
}

void QCDGenIso::endJob() 
{
//   fEffIsoCone03->BayesDivide(fIsoCone0301,fIsoCone0300);
  fEffIso->Divide(fIsoPt01,fIsoPt00,1.0,1.0,"B");
  fEffIsoEta->Divide(fIsoEta01, fIsoEta00,1.0,1.0,"B");

  fEffTrkIso->Divide(fTrkIsoPt01,fTrkIsoPt00,1.0,1.0,"B");
  fEffIsoCone03->Divide(fIsoCone0301,fIsoCone0300,1.0,1.0,"B");

  std::cout << "QCDGenIso has " << _nIn << " input events and " << _nOut  << " output events" << endl;
}


// analyzing reconstructed electrons, muons, and photons

// void QCDGenIso::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
bool QCDGenIso::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // initialization
  bool accept = false;
  _nIn++;
  _alg.init(iEvent, false, false, false); 
  float ptHat = _alg.ptHat(); 
  fPtHat->Fill(ptHat);

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
	  && it_gen->pdgId()<22){
	fPtHard->Fill(it_gen->pt());
      }

      // check if there's a stable photon
      float genpt   = it_gen->pt();
      float geneta  = it_gen->eta();
      int genMomPID = it_gen->mother()? it_gen->mother()->pdgId():_pdgCode;
      if(it_gen->pdgId()!=_pdgCode || it_gen->status()!=1)continue;
   
      float gen_caliso04 = _alg.getGenCalIso(it_gen,0.0,0.4);
      float gen_trkiso04 = _alg.getGenTrkIso(it_gen,0.0,0.4);

      float gen_caliso03 = _alg.getGenCalIso(it_gen,0.0,0.3);
      float gen_trkiso03 = _alg.getGenTrkIso(it_gen,0.0,0.3);

      fPtStablePho->Fill(genpt);
      histIsoDR04_all->Fill(gen_caliso04);

      // if this is a photon from BREM or fragmentation
      if(fabs(genMomPID) <= 21)
	{
	      // fragmentation
	  if(it_gen->mother()->numberOfDaughters()>2)
	    {	
	      fPtFragPho      ->Fill(genpt);
	      if(genpt>20.0 && fabs(geneta) < 2.5)
		histIsoDR04_frag->Fill(gen_caliso04);
	    }
	  // brem
	  else
	    {
	      fPtBremPho      ->Fill(genpt);
	      if(genpt>20.0 && fabs(geneta) < 2.5)
		histIsoDR04_brem->Fill(gen_caliso04);		    
	    }
	} 
      // if it's a photon not from BREM or fragmentation
      else
	{
	  
	  fPtDecayPho      ->Fill(genpt);
	  if(genpt>20.0 && fabs(geneta) < 2.5)
	    histIsoDR04_decay->Fill(gen_caliso04);		    

	}

    } // end of filling generator-level information    
  
  } // if have generator-level information


  if(accept){_nOut++;}
  return accept;
  
}



//define this as a plug-in
DEFINE_FWK_MODULE(QCDGenIso);


