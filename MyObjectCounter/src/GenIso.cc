// -----------------------------------------------------
// GenIso.cc 
// -- a very simple example analyzer for RECO data
// ----------------------------------------------------- 
// Shin-Shan Yu

#include "syu/MyObjectCounter/header/GenIso.hh" 
#include "DataFormats/Math/interface/deltaPhi.h"
#include <algorithm>


GenIso::GenIso(const edm::ParameterSet& iConfig):
  _nIn(0), _nOut(0), _alg(iConfig),
  _pdgCode(iConfig.getUntrackedParameter<int>("pdgCode",22)),
  _etaMax(iConfig.getUntrackedParameter<double>("etaMax",2.5)),
  _ptMin(iConfig.getUntrackedParameter<double>("ptMin",0.0)),
  _isoMaxCut(iConfig.getUntrackedParameter<double>("isoMax",2.0))
{  
}


GenIso::~GenIso()
{
}

void GenIso::beginJob(const edm::EventSetup&)
{
  TFileDirectory results = TFileDirectory( fs->mkdir("GenIso") );

  fPtHat = fs->make<TH1F>("fPtHat","Event pthat", 200,0.0,200.0);

  // angular separation between photons and jets

  const int   ndphibin = 50;
  const float dphimin  = 0.0;
  const float dphimax  = TMath::Pi();
  
  dphigj_match  = fs->make<TH1F>("dphigj_match","for photon pt close to pthat" 
				 " within 10%", ndphibin, dphimin, dphimax);
  dphigj_lowpt  = fs->make<TH1F>("dphigj_lowpt","for photon pt not close to pthat" 
				 " within 50%", ndphibin, dphimin, dphimax);
  dphigj_mismatch  = fs->make<TH1F>("dphigj_mismatch","for photon pt "
				    "much lower than hard-scattering photon pt", 
				    ndphibin, dphimin, dphimax);


  const int   ndetabin = 100;
  const float detamin  = 0;
  const float detamax  = 10.0;

  detagj_match  = fs->make<TH1F>("detagj_match","for photon pt close to pthat" 
				 " within 10%", ndetabin, detamin, detamax);
  detagj_lowpt  = fs->make<TH1F>("detagj_lowpt","for photon pt not close to pthat" 
				 " within 50%", ndetabin, detamin, detamax);
  detagj_mismatch  = fs->make<TH1F>("detagj_mismatch","for photon pt "
				    "much lower than hard-scattering photon pt", 
				    ndetabin, detamin, detamax);

  drgj_match  = fs->make<TH1F>("dr_match","for photon pt close to pthat"
                                 " within 10%", ndetabin, detamin, detamax);
  drgj_lowpt  = fs->make<TH1F>("dr_lowpt","for photon pt not close to pthat"
                                 " within 50%", ndetabin, detamin, detamax);
  drgj_mismatch  = fs->make<TH1F>("dr_mismatch","for photon pt "
                                    "much lower than hard-scattering photon pt",
                                    ndetabin, detamin, detamax);
  

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

  fIsoPt00_match = fs->make<TH1F>("fIsoPt00_match", "(Matched to pthat) Photon pt before isolation cuts", nptbin, xptmin, xptmax);  
  fIsoPt01_match = fs->make<TH1F>("fIsoPt01_match", "(Matched to pthat) Photon pt after isolation Et < 2 GeV cuts", nptbin, xptmin, xptmax);
  fEffIso_match  = fs->make<TH1F>("fEffIso_match", "(Matched to pthat) Photon isolation efficiency", nptbin, xptmin, xptmax);
  fIsoPt00_match->Sumw2();
  fIsoPt01_match->Sumw2();
  fEffIso_match->Sumw2();

  fIsoPt00_lowpt = fs->make<TH1F>("fIsoPt00_lowpt", "(pthat-ptmin diff>50%) Photon pt before isolation cuts", nptbin, xptmin, xptmax);  
  fIsoPt01_lowpt = fs->make<TH1F>("fIsoPt01_lowpt", "(pthat-ptmin diff>50%) Photon pt after isolation Et < 2 GeV cuts", nptbin, xptmin, xptmax);
  fEffIso_lowpt  = fs->make<TH1F>("fEffIso_lowpt", "(pthat-ptmin diff>50%) Photon isolation efficiency", nptbin, xptmin, xptmax);
  fIsoPt00_lowpt->Sumw2();
  fIsoPt01_lowpt->Sumw2();
  fEffIso_lowpt->Sumw2();

  fIsoPt00_mismatch = fs->make<TH1F>("fIsoPt00_mismatch", "(Mismatched to hardscattering) Photon pt before isolation cuts", nptbin, xptmin, xptmax);  
  fIsoPt01_mismatch = fs->make<TH1F>("fIsoPt01_mismatch", "(Mismatched to hardscattering) Photon pt after isolation Et < 2 GeV cuts", nptbin, xptmin, xptmax);
  fEffIso_mismatch  = fs->make<TH1F>("fEffIso_mismatch", "(Mismatched to hardscattering) Photon isolation efficiency", nptbin, xptmin, xptmax);
  fIsoPt00_mismatch->Sumw2();
  fIsoPt01_mismatch->Sumw2();
  fEffIso_mismatch->Sumw2();


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

  fIsoEta00_match = fs->make<TH1F>("fIsoEta00_match", "(Matched to pthat) Photon eta before isolation cuts", netabin, xetamin, xetamax);  
  fIsoEta01_match = fs->make<TH1F>("fIsoEta01_match", "(Matched to pthat) Photon eta after isolation Et < 2 GeV cuts", netabin, xetamin, xetamax);
  fEffIsoEta_match  = fs->make<TH1F>("fEffIsoEta_match", "(Matched to pthat) Photon isolation efficiency", netabin, xetamin, xetamax);
  fIsoEta00_match->Sumw2();
  fIsoEta01_match->Sumw2();
  fEffIsoEta_match->Sumw2();

  fIsoEta00_lowpt = fs->make<TH1F>("fIsoEta00_lowpt", "(pthat-ptmin diff>50%) Photon eta before isolation cuts", netabin, xetamin, xetamax);  
  fIsoEta01_lowpt = fs->make<TH1F>("fIsoEta01_lowpt", "(pthat-ptmin diff>50%) Photon eta after isolation Et < 2 GeV cuts", netabin, xetamin, xetamax);
  fEffIsoEta_lowpt  = fs->make<TH1F>("fEffIsoEta_lowpt", "(pthat-ptmin diff>50%) Photon isolation efficiency", netabin, xetamin, xetamax);
  fIsoEta00_lowpt->Sumw2();
  fIsoEta01_lowpt->Sumw2();
  fEffIsoEta_lowpt->Sumw2();

  fIsoEta00_mismatch = fs->make<TH1F>("fIsoEta00_mismatch", "(Mismatched to hardscattering) Photon eta before isolation cuts", netabin, xetamin, xetamax);  
  fIsoEta01_mismatch = fs->make<TH1F>("fIsoEta01_mismatch", "(Mismatched to hardscattering) Photon eta after isolation Et < 2 GeV cuts", netabin, xetamin, xetamax);
  fEffIsoEta_mismatch  = fs->make<TH1F>("fEffIsoEta_mismatch", "(Mismatched to hardscattering) Photon isolation efficiency", netabin, xetamin, xetamax);
  fIsoEta00_mismatch->Sumw2();
  fIsoEta01_mismatch->Sumw2();
  fEffIsoEta_mismatch->Sumw2();


  std::cout << "Focusing on PDG code = " << _pdgCode << " with etaMax = " << 
    _etaMax << " and ptMin = " << _ptMin << std::endl;
}

void GenIso::endJob() 
{
//   fEffIsoCone03->BayesDivide(fIsoCone0301,fIsoCone0300);
  fEffIso->Divide(fIsoPt01,fIsoPt00,1.0,1.0,"B");
  fEffIsoEta->Divide(fIsoEta01, fIsoEta00,1.0,1.0,"B");

  fEffIso_match->Divide(fIsoPt01_match,fIsoPt00_match,1.0,1.0,"B");
  fEffIsoEta_match->Divide(fIsoEta01_match, fIsoEta00_match,1.0,1.0,"B");

  fEffIso_mismatch->Divide(fIsoPt01_mismatch,fIsoPt00_mismatch,1.0,1.0,"B");
  fEffIsoEta_mismatch->Divide(fIsoEta01_mismatch, fIsoEta00_mismatch,1.0,1.0,"B");

  fEffIso_lowpt->Divide(fIsoPt01_lowpt,fIsoPt00_lowpt,1.0,1.0,"B");
  fEffIsoEta_lowpt->Divide(fIsoEta01_lowpt, fIsoEta00_lowpt,1.0,1.0,"B");

  fEffTrkIso->Divide(fTrkIsoPt01,fTrkIsoPt00,1.0,1.0,"B");
  fEffIsoCone03->Divide(fIsoCone0301,fIsoCone0300,1.0,1.0,"B");

  std::cout << "GenIso has " << _nIn << " input events and " << _nOut  << " output events" << endl;
}


// analyzing reconstructed electrons, muons, and photons

// void GenIso::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
bool GenIso::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // initialization
  bool accept = false;
  _nIn++;
  _alg.init(iEvent, false, false, false); 
  float ptHat = _alg.ptHat(); 
  fPtHat->Fill(ptHat);

  float stable_gammapt = -1;
  float stable_gammaeta= -999.;
  float stable_gammacaliso04 = -999.;
  float stable_gammatrkiso04 = -999.;
  float stable_gammacaliso03 = -999.;
  float stable_gammatrkiso03 = -999.;


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
	  && it_gen->pdgId()==22){
	gamma_hard = it_gen;
	hasOut1   = true;
      }
      else if((Index == 6 || Index==7) && it_gen->status()==3
	      && it_gen->pdgId()!=22){
	parton_hard = it_gen;
	hasOut2   = true;
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

      // if this is a hard scattering photon
      if(genMomPID ==_pdgCode && it_gen->pdgId() == _pdgCode && 
	 it_gen->mother()->status() == 3)
	{
	  // assign variables for later use
	  if(genpt > stable_gammapt)
	    {
	      stable_gammapt       = genpt;
	      stable_gammaeta      = geneta; 
	      stable_gammacaliso04  = gen_caliso04;
	      stable_gammatrkiso04  = gen_trkiso04;
	      stable_gammacaliso03  = gen_caliso03;
	      stable_gammatrkiso03  = gen_trkiso03;
	    }

	} // if it's a hard scattering photon

    } // end of filling generator-level information    
  
    // check only the first 8 particles
    bool ptMinMatchedToPtHat=false;
    bool ptMinDiffFromPtHat =false;
    bool ptMinDiffFromPtHard=false;
    
    if(hasOut1 && hasOut2 && stable_gammapt > 0){
      float dphi = fabs(reco::deltaPhi(gamma_hard->phi(),
				       parton_hard->phi()));
      float deta = fabs(gamma_hard->eta()-parton_hard->eta());
      
      float dr   = sqrt(dphi*dphi+deta*deta);
 
      float diffpt =  fabs( (ptHat-stable_gammapt)/ptHat);

      //      cout << "diffpt = " << diffpt << endl;
      
      if(diffpt < 0.1)
	{
	  //	  cout << "matched pthat!!" << endl;
	  accept = true;
	  dphigj_match->Fill(dphi);
	  detagj_match->Fill(deta);	  
	  drgj_match  ->Fill(dr);
	  ptMinMatchedToPtHat = true;

	}
      
      else if(diffpt > 0.5)
	{
	  //	  cout << "low pt" << endl;
	  dphigj_lowpt->Fill(dphi);
	  detagj_lowpt->Fill(deta);
	  drgj_lowpt  ->Fill(dr);
	  ptMinDiffFromPtHat = true;
	  
	}

      if( fabs(stable_gammapt - gamma_hard->pt())/gamma_hard->pt() > 0.4)
	{
	  //	  cout << "mismatched pthat" << endl;	  
	  dphigj_mismatch->Fill(dphi);
	  detagj_mismatch->Fill(deta);
	  drgj_mismatch  ->Fill(dr);
	  ptMinDiffFromPtHard = true;
	}

      // now plot efficiency figures

      if(fabs(stable_gammaeta)< _etaMax)
	{
	  fIsoPt00 ->Fill(stable_gammapt);
	  if(stable_gammacaliso04 < _isoMaxCut)
	    fIsoPt01 ->Fill(stable_gammapt);

	  if(ptMinMatchedToPtHat){
	    fIsoPt00_match ->Fill(stable_gammapt);
	    if(stable_gammacaliso04 < _isoMaxCut)
	      fIsoPt01_match ->Fill(stable_gammapt);	  
	  }
	  
	  if(ptMinDiffFromPtHat){
	    fIsoPt00_lowpt ->Fill(stable_gammapt);
	    if(stable_gammacaliso04 < _isoMaxCut)
	      fIsoPt01_lowpt ->Fill(stable_gammapt);	  
	  }
	  

	  if(ptMinDiffFromPtHard){
	    fIsoPt00_mismatch ->Fill(stable_gammapt);
	    if(stable_gammacaliso04 < _isoMaxCut)
	      fIsoPt01_mismatch ->Fill(stable_gammapt);	  
	  }
	  


	  fIsoCone0300->Fill(stable_gammapt);
	  if(stable_gammacaliso03 < _isoMaxCut)
	    fIsoCone0301->Fill(stable_gammapt);
	  
	  fTrkIsoPt00 ->Fill(stable_gammapt);
	  if(stable_gammatrkiso04 < _isoMaxCut)
	    fTrkIsoPt01 ->Fill(stable_gammapt);
	}
      

      if(stable_gammapt>_ptMin)
	{
	  fIsoEta00->Fill(stable_gammaeta);
	  if(stable_gammacaliso04 < _isoMaxCut)
		fIsoEta01->Fill(stable_gammaeta);

	  if(ptMinMatchedToPtHat){
	    fIsoEta00_match ->Fill(stable_gammaeta);
	    if(stable_gammacaliso04 < _isoMaxCut)
	      fIsoEta01_match ->Fill(stable_gammaeta);	  
	  }
	  
	  if(ptMinDiffFromPtHat){
	    fIsoEta00_lowpt ->Fill(stable_gammaeta);
	    if(stable_gammacaliso04 < _isoMaxCut)
	      fIsoEta01_lowpt ->Fill(stable_gammaeta);	  
	  }
	  

	  if(ptMinDiffFromPtHard){
	    fIsoEta00_mismatch ->Fill(stable_gammaeta);
	    if(stable_gammacaliso04 < _isoMaxCut)
	      fIsoEta01_mismatch ->Fill(stable_gammaeta);	  
	  }
	  

	}
     

    }


  } // if have generator-level information


  if(accept){_nOut++;}
  return accept;
  
}



//define this as a plug-in
DEFINE_FWK_MODULE(GenIso);


