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
  _ptMin(iConfig.getUntrackedParameter<double>("ptMin",15.0)),
  _isoMaxCut(iConfig.getUntrackedParameter<double>("isoMax",2.0))
{  
}


GenIso::~GenIso()
{
}

void GenIso::beginJob()
{
  TFileDirectory results = TFileDirectory( fs->mkdir("GenIso") );

  fPtHat = fs->make<TH1F>("fPtHat","Event pthat", 200,0.0,100.0);
  fPtPho = fs->make<TH1F>("fPtPho","Event pt of stable photon", 200,0.0,100.0);
  fPtStablePho = fs->make<TH1F>("fPtStablePho","Event pt of stable photon", 200,0.0,100.0);
  fPtHard = fs->make<TH1F>("fPtHard","Event pt of hard scattering photon", 200,0.0,100.0);

  histIsoDR04_all = fs->make<TH1F>("histIsoDR04_all","Isolation for all photons",200,0.0,20.0) ;
  histIsoDR04_goodpt = fs->make<TH1F>("histIsoDR04_goodpt","Isolation for photons with good pt",200,0.0,20.0);


  // angular separation between photons and jets

  const int   ndphibin = 50;
  const float dphimin  = 0.0;
  const float dphimax  = TMath::Pi();
  
  dphigj_match  = fs->make<TH1F>("dphigj_match","for photon pt close to pthat" 
				 " within 10%", ndphibin, dphimin, dphimax);
  dphigj_lopt  = fs->make<TH1F>("dphigj_lopt","for photon pt not close to pthat" 
				 " within -50%", ndphibin, dphimin, dphimax);
  dphigj_hipt  = fs->make<TH1F>("dphigj_hipt","for photon pt not close to pthat" 
				 " within +50%", ndphibin, dphimin, dphimax);
  dphigj_mismatch  = fs->make<TH1F>("dphigj_mismatch","for photon pt "
				    "much lower than hard-scattering photon pt", 
				    ndphibin, dphimin, dphimax);


  const int   ndetabin = 100;
  const float detamin  = 0;
  const float detamax  = 10.0;

  detagj_match  = fs->make<TH1F>("detagj_match","for photon pt close to pthat" 
				 " within 10%", ndetabin, detamin, detamax);
  detagj_lopt  = fs->make<TH1F>("detagj_lopt","for photon pt not close to pthat" 
				 " within -50%", ndetabin, detamin, detamax);
  detagj_hipt  = fs->make<TH1F>("detagj_hipt","for photon pt not close to pthat" 
				 " within +50%", ndetabin, detamin, detamax);
  detagj_mismatch  = fs->make<TH1F>("detagj_mismatch","for photon pt "
				    "much lower than hard-scattering photon pt", 
				    ndetabin, detamin, detamax);

  drgj_match  = fs->make<TH1F>("dr_match","for photon pt close to pthat"
                                 " within 10%", ndetabin, detamin, detamax);
  drgj_lopt  = fs->make<TH1F>("dr_lopt","for photon pt not close to pthat"
                                 " within -50%", ndetabin, detamin, detamax);
  drgj_hipt  = fs->make<TH1F>("dr_hipt","for photon pt not close to pthat"
                                 " within +50%", ndetabin, detamin, detamax);
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

  fIsoPt00_lopt = fs->make<TH1F>("fIsoPt00_lopt", "(pthat-ptmin diff>50%) Photon pt before isolation cuts", nptbin, xptmin, xptmax);  
  fIsoPt01_lopt = fs->make<TH1F>("fIsoPt01_lopt", "(pthat-ptmin diff>50%) Photon pt after isolation Et < 2 GeV cuts", nptbin, xptmin, xptmax);
  fEffIso_lopt  = fs->make<TH1F>("fEffIso_lopt", "(pthat-ptmin diff>50%) Photon isolation efficiency", nptbin, xptmin, xptmax);
  fIsoPt00_lopt->Sumw2();
  fIsoPt01_lopt->Sumw2();
  fEffIso_lopt->Sumw2();

  fIsoPt00_hipt = fs->make<TH1F>("fIsoPt00_hipt", "(pthat-ptmin diff>50%) Photon pt before isolation cuts", nptbin, xptmin, xptmax);  
  fIsoPt01_hipt = fs->make<TH1F>("fIsoPt01_hipt", "(pthat-ptmin diff>50%) Photon pt after isolation Et < 2 GeV cuts", nptbin, xptmin, xptmax);
  fEffIso_hipt  = fs->make<TH1F>("fEffIso_hipt", "(pthat-ptmin diff>50%) Photon isolation efficiency", nptbin, xptmin, xptmax);
  fIsoPt00_hipt->Sumw2();
  fIsoPt01_hipt->Sumw2();
  fEffIso_hipt->Sumw2();


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

  fIsoEta00_lopt = fs->make<TH1F>("fIsoEta00_lopt", "(pthat-ptmin diff>50%) Photon eta before isolation cuts", netabin, xetamin, xetamax);  
  fIsoEta01_lopt = fs->make<TH1F>("fIsoEta01_lopt", "(pthat-ptmin diff>50%) Photon eta after isolation Et < 2 GeV cuts", netabin, xetamin, xetamax);
  fEffIsoEta_lopt  = fs->make<TH1F>("fEffIsoEta_lopt", "(pthat-ptmin diff>50%) Photon isolation efficiency", netabin, xetamin, xetamax);
  fIsoEta00_lopt->Sumw2();
  fIsoEta01_lopt->Sumw2();
  fEffIsoEta_lopt->Sumw2();


  fIsoEta00_hipt = fs->make<TH1F>("fIsoEta00_hipt", "(pthat-ptmin diff>50%) Photon eta before isolation cuts", netabin, xetamin, xetamax);  
  fIsoEta01_hipt = fs->make<TH1F>("fIsoEta01_hipt", "(pthat-ptmin diff>50%) Photon eta after isolation Et < 2 GeV cuts", netabin, xetamin, xetamax);
  fEffIsoEta_hipt  = fs->make<TH1F>("fEffIsoEta_hipt", "(pthat-ptmin diff>50%) Photon isolation efficiency", netabin, xetamin, xetamax);
  fIsoEta00_hipt->Sumw2();
  fIsoEta01_hipt->Sumw2();
  fEffIsoEta_hipt->Sumw2();


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

  fEffIso_lopt->Divide(fIsoPt01_lopt,fIsoPt00_lopt,1.0,1.0,"B");
  fEffIsoEta_lopt->Divide(fIsoEta01_lopt, fIsoEta00_lopt,1.0,1.0,"B");

  fEffIso_hipt->Divide(fIsoPt01_hipt,fIsoPt00_hipt,1.0,1.0,"B");
  fEffIsoEta_hipt->Divide(fIsoEta01_hipt, fIsoEta00_hipt,1.0,1.0,"B");

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
	fPtHard->Fill(it_gen->pt());
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
	  fPtPho->Fill(it_gen->pt());
 	  if(it_gen->pt()>20.0 && fabs(it_gen->eta()) < 2.5)
	    histIsoDR04_all->Fill(gen_caliso04);
	  
	  // assign variables for later use, look for the hardest photon
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
  

    fPtStablePho->Fill(stable_gammapt);
    // check only the first 8 particles
    bool ptMinMatchedToPtHat=false;
    bool ptMinDiffFromPtHat_plus =false;
    bool ptMinDiffFromPtHat_minus =false;
    bool ptMinDiffFromPtHard=false;
    
    if(hasOut1 && hasOut2 && stable_gammapt > 0){
      float dphi = fabs(reco::deltaPhi(gamma_hard->phi(),
				       parton_hard->phi()));
      float deta = fabs(gamma_hard->eta()-parton_hard->eta());
      
      float dr   = sqrt(dphi*dphi+deta*deta);
 
      float diffpt =  (stable_gammapt-ptHat)/ptHat ;

      //      cout << "diffpt = " << diffpt << endl;
      
      if(fabs(diffpt) < 0.2)
	{
	  //	  cout << "matched pthat!!" << endl;
	  accept = true;
	  dphigj_match->Fill(dphi);
	  detagj_match->Fill(deta);	  
	  drgj_match  ->Fill(dr);
	  ptMinMatchedToPtHat = true;
 	  if(stable_gammapt>20.0 && fabs(stable_gammaeta) < 2.5)
	    histIsoDR04_goodpt->Fill(stable_gammacaliso04);
	}
      
      else if(diffpt < -0.5)
	{
	  //	  cout << "low pt" << endl;
	  dphigj_lopt->Fill(dphi);
	  detagj_lopt->Fill(deta);
	  drgj_lopt  ->Fill(dr);
	  ptMinDiffFromPtHat_minus = true;
	  
	}


      else if(diffpt > 0.5)
	{
	  //	  cout << "low pt" << endl;
	  dphigj_hipt->Fill(dphi);
	  detagj_hipt->Fill(deta);
	  drgj_hipt  ->Fill(dr);
	  ptMinDiffFromPtHat_plus = true;
	  
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
	  
	  if(ptMinDiffFromPtHat_minus){
	    fIsoPt00_lopt ->Fill(stable_gammapt);
	    if(stable_gammacaliso04 < _isoMaxCut)
	      fIsoPt01_lopt ->Fill(stable_gammapt);	  
	  }
	  
	  if(ptMinDiffFromPtHat_plus){
	    fIsoPt00_hipt ->Fill(stable_gammapt);
	    if(stable_gammacaliso04 < _isoMaxCut)
	      fIsoPt01_hipt ->Fill(stable_gammapt);	  
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
	  
	  if(ptMinDiffFromPtHat_minus){
	    fIsoEta00_lopt ->Fill(stable_gammaeta);
	    if(stable_gammacaliso04 < _isoMaxCut)
	      fIsoEta01_lopt ->Fill(stable_gammaeta);	  
	  }
	  
	  if(ptMinDiffFromPtHat_plus){
	    fIsoEta00_hipt ->Fill(stable_gammaeta);
	    if(stable_gammacaliso04 < _isoMaxCut)
	      fIsoEta01_hipt ->Fill(stable_gammaeta);	  
	  }

	  if(ptMinDiffFromPtHard){
	    fIsoEta00_mismatch ->Fill(stable_gammaeta);
	    if(stable_gammacaliso04 < _isoMaxCut)
	      fIsoEta01_mismatch ->Fill(stable_gammaeta);	  
	  }
	  

	}
     

    }


  } // if have generator-level information


  // the following code is ported from Vasu
//   using namespace HepMC;  
//   edm::Handle<HepMCProduct> hepMCProduct = _alg.getHepMCHandle();
//   const GenEvent* genEventInfo = hepMCProduct->GetEvent();
  

//   for(GenEvent::vertex_const_iterator vertexIter=genEventInfo->vertices_begin();vertexIter!=genEventInfo->vertices_end(); vertexIter++){
//     HepMC::GenVertex *genVtxTemp = (*vertexIter);

//     if(genVtxTemp->particles_in_size()!=2)
//       continue;
    
//     if(genVtxTemp->particles_out_size()!=2)
//       continue;   

//     iBarCode =genVtxTemp->barcode();
//     iVerticiesLeft++;
 
//   }

//   if(iVerticiesLeft>1)
//     std::cout<<"************************************************** More than one vertex left"<<std::endl;

//   HepMC::GenVertex *genVtxHard = genEventInfo->barcode_to_vertex(iBarCode);

//   if(genVtxHard&&iVerticiesLeft==1){
//     for(GenVertex::particles_out_const_iterator particleOut=genVtxHard->particles_out_const_begin();particleOut!=genVtxHard->particles_out_const_end(); particleOut++){
//       HepMC::GenParticle* genParticle = (*particleOut);
      
//       if(genParticle->pdg_id()!=22&&genParticle->status()!=1)
// 	continue;
//   if(accept)    
//   genEventInfo->print();
//     }
//   }
 



  if(accept){_nOut++;}
  return accept;
  
}



//define this as a plug-in
DEFINE_FWK_MODULE(GenIso);


