// -----------------------------------------------------
// GenIso.cc 
// -- a very simple example analyzer for RECO data
// ----------------------------------------------------- 
// Shin-Shan Yu

#include "syu/MyObjectCounter/header/GenIso.hh" 
#include <algorithm>


GenIso::GenIso(const edm::ParameterSet& iConfig):
  _nIn(0), _nOut(0), _alg(iConfig)
{  
  _pdgCode = iConfig.getUntrackedParameter<int>("pdgCode",22);
  _etaMax = iConfig.getUntrackedParameter<double>("etaMax",2.5);  
  _ptMin  = iConfig.getUntrackedParameter<double>("ptMin",0.0);  
}


GenIso::~GenIso()
{
}

void GenIso::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  TFileDirectory results = TFileDirectory( fs->mkdir("GenIso") );

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

void GenIso::endJob() 
{
//   fEffIso->BayesDivide(fIsoPt01,fIsoPt00, "w");
//   fEffIsoEta->BayesDivide(fIsoEta01,fIsoEta00, "w");
  fEffIso->Divide(fIsoPt01,fIsoPt00,1.0,1.0,"B");
  fEffTrkIso->Divide(fTrkIsoPt01,fTrkIsoPt00,1.0,1.0,"B");
  fEffIsoEta->Divide(fIsoEta01, fIsoEta00,1.0,1.0,"B");

  std::cout << "GenIso has " << _nIn << " input events and " << _nOut  << " output events" << endl;
}


// analyzing reconstructed electrons, muons, and photons

void GenIso::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  _nIn++;
  _alg.init(iEvent, false, false, false); 
 
  edm::Handle<reco::GenParticleCollection> GenHandle = _alg.getGenHandle();
  if(GenHandle.isValid()){
    for( GenParticleCollection::const_iterator it_gen = 
	   GenHandle->begin(); 
	 it_gen != GenHandle->end(); it_gen++ ) {

      float genpt   = it_gen->pt();
      float geneta  = it_gen->eta();
      int genMomPID = it_gen->mother()? it_gen->mother()->pdgId():_pdgCode;
      if(it_gen->pdgId()!=_pdgCode || it_gen->pt() < 2. 
	 || it_gen->status()!=1)continue;
   
      float gen_caliso = _alg.getGenCalIso(it_gen);
      float gen_trkiso = _alg.getGenTrkIso(it_gen);

      if(genMomPID ==22 && it_gen->pdgId() == 22 && 
	 it_gen->mother()->status() == 3)
	{

	  if(fabs(geneta)< _etaMax)
	    {
	      fIsoPt00 ->Fill(genpt);
	      if(gen_caliso < 2.0)
		fIsoPt01 ->Fill(genpt);
	      
	      fTrkIsoPt00 ->Fill(genpt);
	      if(gen_trkiso < 2.0)
		fTrkIsoPt01 ->Fill(genpt);
	    }


	  if(genpt>_ptMin)
	    {
	      fIsoEta00->Fill(geneta);
	      if(gen_caliso < 2.0)
		fIsoEta01->Fill(geneta);
	    }
	
	} // if it's a hard scattering photon

    } // end of filling generator-level information    
  } // if have generator-level information

  _nOut++;
  
}



//define this as a plug-in
DEFINE_FWK_MODULE(GenIso);


