/// @file
/// File containing the definition of the methods associated to the class.
///

#include "DelPanj/TreeMaker/interface/ZZTree.hh"
#include "DelPanj/TreeMaker/interface/cutvalues.h"
#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
#include "DelPanj/TreeMaker/interface/MuonEffectiveArea.h"

// system include files

#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
// #include "DataFormats/RecoCandidate/interface/IsoDeposit.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "PhysicsTools/CandUtils/interface/CenterOfMassBooster.h"
#include "PhysicsTools/CandUtils/interface/Booster.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// ROOT classes
#include <TMath.h>
#include <TVector3.h>
#include <algorithm>
#include <fstream>
#include <Math/VectorUtil.h>
#include <TLegend.h>
#include <TCanvas.h>

typedef std::vector< edm::Handle< edm::ValueMap<double> > >             
IsoDepositVals;


void
ZZTree::AddBranch(int* x, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),x,(brName+"/I").c_str());
}


void
ZZTree::AddBranch(double* x, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),x,(brName+"/D").c_str());
}

//---------------------------------------------------
//---------------------------------------------------
void
ZZTree::AddBranch(std::vector<double>* vec, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);
}

void
ZZTree::AddBranch(std::vector<int>* vec, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);
}


//---------------------------------------------------
//---------------------------------------------------
void 
ZZTree::AddBranchArray(const int arraySize, double* x, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),x,Form("%s[%d]/D",brName.data(),arraySize));
}

//---------------------------------------------------------------
//---------------------------------------------------------------

ZZTree::ZZTree(std::string name, TTree* tree, const edm::ParameterSet& iConfig):
  e2012ID_ ( iConfig.getParameter<edm::ParameterSet>("e2012IDSet")),
  mu2012ID_ ( iConfig.getParameter<edm::ParameterSet>("mu2012IDSet")),
  hzzeejj_(iConfig.getParameter<edm::InputTag>("hzzeejjTag")),
  hzzmmjj_ (iConfig.getParameter<edm::InputTag>("hzzmmjjTag")),
  eleRhoIsoInputTag_(iConfig.getParameter<edm::InputTag>("eleRhoIso")),
  muoRhoIsoInputTag_(iConfig.getParameter<edm::InputTag>("muoRhoIso"))//,
	      //   primaryVertexInputTag_(iConfig.getParameter<edm::InputTag>("primaryVertex")),

{
  tree_=tree; 
  SetBranches();

  
  // the second argument is the random seed, any reason to set it 
  // differently or the same for the 3 taggers

  btsfutiltch = new BTagSFUtil("TCHE", 13); 
  btsfutilcsv = new BTagSFUtil("CSV", 13);
  btsfutiljp = new BTagSFUtil("JP", 13);

}


ZZTree::~ZZTree()
{
  delete tree_;
  delete btsfutiltch;
  delete btsfutilcsv;
  delete btsfutiljp;
}


// ------------ method called to for each event  ------------
void ZZTree::Fill(const edm::Event& iEvent)
{
  Clear();
  
  bool isData = iEvent.isRealData();

  // vertices
  //   edm::Handle<reco::VertexCollection> vtx_h;
  //   iEvent.getByLabel(primaryVertexInputTag_, vtx_h);
  //   const reco::Vertex &pv = (*vtx_h)[0];

  // missing Et significance
  edm::Handle<pat::METCollection> met_H;
  iEvent.getByLabel("patMETsAK5", met_H);
  metSig_ = (&(met_H->front()))->significance();

  edm::Handle<pat::METCollection> metNoPU_H;
  iEvent.getByLabel("patMETsAK5NoPUSub", metNoPU_H);
  metSigNoPU_ = (&(metNoPU_H->front()))->significance();

  // rho for electron
  edm::Handle<double> ele_rho_event;
  iEvent.getByLabel(eleRhoIsoInputTag_,ele_rho_event);
  double ele_rho= *(ele_rho_event.product());
  e2012ID_.SetData(isData);
  e2012ID_.SetRho(ele_rho);

  // rho for muon
  edm::Handle<double> muo_rho_event;
  iEvent.getByLabel(muoRhoIsoInputTag_,muo_rho_event);
  double muo_rho= *(muo_rho_event.product());
  mu2012ID_.SetData(isData);
  mu2012ID_.SetRho(muo_rho);

  //initialize variables
  int hcand = -1;

  double best_mZjj = 9999999.0;
  int bestHCandIndex = -1;

  //LOOP IN 2 KINDS OF LEPTONS: ELECTRONS AND MUONS 
  for (int ilep=0; ilep<2; ilep++) {
    hcand=-1;

    //GET THE HIGGS->ZZ->LLJJ COLLECTION
    Handle<std::vector<pat::CompositeCandidate> > hzzlljj;
    if(ilep==0) iEvent.getByLabel(hzzeejj_, hzzlljj);
    else iEvent.getByLabel(hzzmmjj_, hzzlljj);

    // LOOP OVER HIGGS CANDIDATES

    for(unsigned i=0; i<hzzlljj->size(); i++){
      const pat::CompositeCandidate & h = (*hzzlljj)[i];	

      // LOOK FOR TWO GOOD CHARGED LEPTONS
      int nLepPtHi=0;
      int nLepPtLo=0;

      for(unsigned int il=0; il < 2; il++)
	{

	  const reco::Candidate*  lep = h.daughter(LEPZ)->daughter(il)->masterClone().get();
	  double pt = lep->pt();
	  
 	  if(pt >MIN_LEPPT1) nLepPtHi++;
 	  if(pt >MIN_LEPPT2) nLepPtLo++;

	}

      if(nLepPtHi < 1)continue;
      if(nLepPtLo < 2)continue;
      
      bool OppCharge = h.daughter(LEPZ)->daughter(0)->masterClone().get()->charge()*
	h.daughter(LEPZ)->daughter(1)->masterClone().get()->charge() < 0 ? true: false;

      if(!OppCharge)continue;


      // they need to pass ID cuts
      int nPassID=0;

      if(ilep==0){ // electron
	nPassID=0;

	for(unsigned int iele=0; iele < 2; iele++){
	  
	  const pat::Electron* myEle
	    = dynamic_cast<const pat::Electron*>(h.daughter(LEPZ)->daughter(iele)->masterClone().get());

	  std::map<std::string, bool> Pass = e2012ID_.CutRecord(*myEle); 
   	  if(!PassAll(Pass))continue; // 2012 loose electron ID	  
	  nPassID++;
	  
	  /*
	    eID01.push_back(myEle->pt());
	    eID02.push_back(myEle->superCluster()->eta());
	    eID03.push_back(myEle->deltaEtaSuperClusterTrackAtVtx());
	    eID04.push_back(myEle->deltaPhiSuperClusterTrackAtVtx());
	    eID05.push_back(myEle->sigmaIetaIeta());
	    eID06.push_back(myEle->hadronicOverEm());
	    eID07.push_back(myEle->userFloat("dxy"));
	    eID08.push_back(myEle->userFloat("dz"));
	    eID09.push_back(fabs(1.0/myEle->ecalEnergy() - 
	    myEle->eSuperClusterOverP()/myEle->ecalEnergy()));
	    eID10.push_back(myEle->convDcot());
	    eID11.push_back(myEle->convDist());
	    eID12.push_back(myEle->userFloat("hasMatchConv"));
	    eID13.push_back(myEle->gsfTrack().get()->trackerExpectedHitsInner().numberOfHits());
	  
	    double iso1 = myEle->chargedHadronIso();
	    double iso2 = myEle->neutralHadronIso();
	    double iso3 = myEle->photonIso();
	  
	    ElectronEffectiveArea::ElectronEffectiveAreaTarget effAreaTarget_ = 
	    isData? ElectronEffectiveArea::kEleEAData2011:
	    ElectronEffectiveArea::kEleEAFall11MC;
      
	    ElectronEffectiveArea::ElectronEffectiveAreaType effAreaType_ =
	    ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03;      
 
	    double AEff = ElectronEffectiveArea::GetElectronEffectiveArea
	    (effAreaType_, fabs(myEle->superCluster()->eta()), effAreaTarget_);

	    double iso4 = iso1 + std::max(0.0, iso2+iso3-ele_rho*AEff);
	  
	    eID14.push_back(iso1);
	    eID15.push_back(iso2);
	    eID16.push_back(iso3);
	    eID17.push_back(iso4);
	  */

	}
      } // if it's an electron type

      else if(ilep==1){ // muon

	nPassID=0;
	for(unsigned int imuo=0; imuo < 2; imuo++){
	 
	  const pat::Muon* myMuo
	    = dynamic_cast<const pat::Muon*>(h.daughter(LEPZ)->daughter(imuo)->masterClone().get());
	  std::map<std::string, bool> Pass = mu2012ID_.CutRecord(*myMuo);
	  if(!PassAll(Pass))continue; // 2012 Tight Muon ID
	  nPassID++;
	  	  
	  /*
	    double pt    =  myMuo->pt();
	    double eta   =  myMuo->eta();
	    double iso1  =  myMuo->pfIsolationR04().sumChargedHadronPt;
	    double iso2  =  myMuo->pfIsolationR04().sumNeutralHadronEt;
	    double iso3  =  myMuo->pfIsolationR04().sumPhotonEt;
	    double isoPU =  myMuo->pfIsolationR04().sumPUPt;    
	    double iso4Beta = iso1 + std::max(iso3+iso2-0.5*isoPU,0.);
	    MuonEffectiveArea::MuonEffectiveAreaTarget effAreaTarget_ = 
 	    isData? MuonEffectiveArea::kMuEAData2012:
 	    MuonEffectiveArea::kMuEAData2012;
	    // 	MuonEffectiveArea::kMuEAFall11MC;
	    MuonEffectiveArea::MuonEffectiveAreaType effAreaType_= 
 	    MuonEffectiveArea::kMuGammaAndNeutralHadronIso04;
	    double Area = MuonEffectiveArea::GetMuonEffectiveArea(
	    effAreaType_, fabs(eta), effAreaTarget_);

	    double iso4Rho =  iso1 + std::max(iso3+iso2-muo_rho*Area,0.);
	    double isoBeta = iso4Beta/pt;	  
	    double isoRho  = iso4Rho/pt;


	    muID01.push_back(isoBeta);
	    // 	  muID01.push_back(myMuo->isGlobalMuon());
	    muID02.push_back(myMuo->isPFMuon());
	    muID03.push_back(myMuo->isTrackerMuon());
	    if(myMuo->isTrackerMuon() && myMuo->isGlobalMuon()){
	    muID04.push_back(myMuo->globalTrack()->normalizedChi2());
	    muID05.push_back(myMuo->globalTrack()->hitPattern().numberOfValidMuonHits());
	    muID06.push_back(myMuo->numberOfMatchedStations());
	    muID07.push_back(myMuo->dB());
	    double dzV =myMuo->userFloat("dzVtx") ;
	    muID08.push_back(dzV);
	    muID09.push_back(myMuo->innerTrack()->hitPattern().numberOfValidPixelHits());
	    muID10.push_back(myMuo->innerTrack()->hitPattern().trackerLayersWithMeasurement());
	    }
	    muID11.push_back(myMuo->pfIsolationR04().sumChargedHadronPt);
	    muID12.push_back(myMuo->pfIsolationR04().sumNeutralHadronEt);
	    muID13.push_back(myMuo->pfIsolationR04().sumPhotonEt);
	    muID14.push_back(myMuo->pfIsolationR04().sumPUPt);
	    muID15.push_back(isoRho);
	  */

	} // end of loop over muon
      } // if is a muon type
     
      
      if(nPassID < 2)continue;


      // LOOK FOR 2 JETS PASSING BETA CUTS
      // STILL NEED TO ADD LOOSE JET ID CUTS

      int nGoodJets=0;
      int nLooseBTags=0;
      int nMediumBTags=0;

      for(unsigned int ijet=0; ijet < 2; ijet++){	 
	
	const pat::Jet * myJet =
	  dynamic_cast<const pat::Jet *>(h.daughter(HADZ)->daughter(ijet)->
					 masterClone().get());

	double pt  = myJet->pt();
 	if(pt < MIN_JETPT)continue;
	

	// to suppress jets from pileups
	double puBeta = myJet->userFloat("puBeta");
 	if(puBeta < MIN_JETBETA)continue;


	bool separatedFromLepton=false;
	for(unsigned int il=0; il < 2; il++)
	  {
	    const reco::Candidate*  lep = h.daughter(LEPZ)->daughter(il)->masterClone().get();

	    if(deltaR(lep->eta(), lep->phi(),
		      myJet->eta(), myJet->phi()) > MIN_DR_JETLEP)
	      separatedFromLepton=true;

	  }

	if(!separatedFromLepton)continue;
	  

	double eta = myJet->eta();
	bool isLoose  = false;
	bool isMedium = false;
	double jpBTag = myJet->bDiscriminator("jetProbabilityBJetTags");
	int flavor    = myJet->partonFlavour();
	
	if(jpBTag > MIN_BTAG_JP_LOOSE)isLoose=true;
	if(jpBTag > MIN_BTAG_JP_MEDIUM)isMedium=true;
	
	if(!isData)
	  btsfutiljp->modifyBTagsWithSF_fast(isLoose, isMedium, pt, eta, 
					     flavor, "mean");
	if(isLoose)  nLooseBTags++;
	if(isMedium) nMediumBTags++;

	nGoodJets++; 
      } // end of loop over jets
      
      // there must be at least two jets
      if(nGoodJets < 2)continue;

      // number of btags
      int nBTags = 0;
      if(nMediumBTags >= 1 && nLooseBTags == 2)nBTags=2;
      else if(nMediumBTags ==0 && nLooseBTags >=1)nBTags=1;
      else nBTags=0;

      hcand = i;
	
      //     }// end of loop over H->ZZ candidates
  
      //     if(hcand <0)continue; // there is no higgs candidate for this lepton type
  
  
      const pat::CompositeCandidate & goodH = (*hzzlljj)[hcand];
      const reco::Candidate * Zll = goodH.daughter(LEPZ);
      const reco::Candidate * Zjj = goodH.daughter(HADZ);
      //Final selection in categories

      double higgsPt_local  = goodH.pt();
      double higgsEta_local = goodH.eta();
      double higgsPhi_local = goodH.phi();
      double higgsM_local   = goodH.mass();
      double higgsM_refit_local = goodH.userFloat("HZZRefitMass");

      double zllPt_local  = Zll->pt();
      double zllEta_local = Zll->eta();
      double zllPhi_local = Zll->phi();
      double zllM_local   = Zll->mass();

      double zjjPt_local  = Zjj->pt();
      double zjjEta_local = Zjj->eta();
      double zjjPhi_local = Zjj->phi();
      double zjjM_local   = Zjj->mass();
      double zjjM_refit_local = goodH.userFloat("ZjjRefitMass");

      double heliLD_local = goodH.userFloat("helyLD");
      double heliLD_refit_local = goodH.userFloat("helyLDRefit");

      // check the dilepton Z mass
      if(zllM_local < MIN_MZ_LL)continue;
      if(zllM_local > MAX_MZ_LL)continue;

      // check the mass of the two jets
      if(zjjM_local < LOOSE_MIN_MZ_JJ)continue;
      if(zjjM_local > LOOSE_MAX_MZ_JJ)continue;

      bestHCandIndex++;
      nAllHCand_++;
      
      if(fabs(zjjM_local - MZ_PDG) < fabs(best_mZjj - MZ_PDG))
	{
	  best_mZjj = zjjM_local;
	  bestHCand_ = bestHCandIndex;
	}
    
      higgsPt_.push_back(higgsPt_local);
      higgsEta_.push_back(higgsEta_local);
      higgsPhi_.push_back(higgsPhi_local);
      higgsM_.push_back(higgsM_local);
      higgsMRefit_.push_back(higgsM_refit_local);

      zllPt_.push_back(zllPt_local);
      zllEta_.push_back(zllEta_local);
      zllPhi_.push_back(zllPhi_local);
      zllM_.push_back(zllM_local);

      zjjPt_.push_back(zjjPt_local);
      zjjEta_.push_back(zjjEta_local);
      zjjPhi_.push_back(zjjPhi_local);
      zjjM_.push_back(zjjM_local);
      zjjMRefit_.push_back(zjjM_refit_local);


      heliLD_.push_back(heliLD_local);
      heliLDRefit_.push_back(heliLD_refit_local);

      nBTags_.push_back(nBTags);
      lepType_.push_back(ilep);

      // CHECK WHICH CUTS THIS HIGGS CANDIDATE PASSES
      int thisBit = 0;

      if(zjjM_local > MIN_MZ_JJ && zjjM_local < MAX_MZ_JJ)
	thisBit |= MZJJ_SIGNAL;
      else
	thisBit |= MZJJ_SIDEBAND;

      if(metSig_ < MAX_MET_SIG)
	thisBit |= PFMET_SIG;

      double heliLD_cutValue = -9999;
      heliLD_cutValue = MIN_HELD_CONSTANT[nBTags] + 
	MIN_HELD_SLOPE[nBTags]* higgsM_refit_local;
      
      if(heliLD_refit_local > heliLD_cutValue)
	thisBit |= HELI_LD;

      double qgld= 99999.0; // needs to be fixed
      if(qgld > MIN_QUARK_GLUON_LD_0BTAG)
	thisBit |= QG_LD;
      
      // need to implement cuts on 4-body mass, need to be fixed
      double mH_input = higgsM_refit_local;
      if(higgsM_refit_local > MIN_MH_RATIO*mH_input &&
	 higgsM_refit_local < MAX_MH_RATIO*mH_input)
	thisBit |= MH_SIGNAL;


      // now check all cuts

      if( (thisBit & MZJJ_SIGNAL) &&
	  (thisBit & MH_SIGNAL)   &&
	  (thisBit & PFMET_SIG)   &&	      
	  (thisBit & HELI_LD)     &&
	  (thisBit & QG_LD))
	thisBit |= ALL_SIGNAL;
      
      else if( (thisBit & MZJJ_SIDEBAND) &&
	       (thisBit & MH_SIGNAL)   &&
	       (thisBit & PFMET_SIG)   &&	      
	       (thisBit & HELI_LD)     &&
	       (thisBit & QG_LD))
	thisBit |= ALL_SIDEBAND;
      
      passBit_.push_back(thisBit);

      if((thisBit & ALL_SIGNAL) ||
	 (thisBit & ALL_SIDEBAND))
	nGoodHCand_++;


      /*


      for(unsigned int ijet=0; ijet < 2; ijet++){	 
      const pat::Jet * thisJet = 
      dynamic_cast<const pat::Jet *>(goodH.daughter(HADZ)->daughter(ijet)->
      masterClone().get());

      jetPt_[ijet]  = thisJet->pt();
      jetEta_[ijet] = thisJet->eta();
      jetPhi_[ijet] = thisJet->phi();
      jetE_[ijet]   = thisJet->energy();

      jetRefitPt_[ijet]  = goodH.userFloat(Form("j%dRefitPt",ijet+1));
      jetRefitEta_[ijet] = goodH.userFloat(Form("j%dRefitEta",ijet+1));
      jetRefitPhi_[ijet] = goodH.userFloat(Form("j%dRefitPhi",ijet+1));
      jetRefitE_[ijet]   = goodH.userFloat(Form("j%dRefitE",ijet+1));


      } // end of loop over jets



      lepType_ = ilep;
    
      for(unsigned int il=0; il < 2; il++){	 

      const reco::Candidate* thisLep = 
      (goodH.daughter(LEPZ)->daughter(il)->
      masterClone().get());
      
      lepPt_[il]  = thisLep->pt();
      lepEta_[il] = thisLep->eta();
      lepPhi_[il] = thisLep->phi();
      lepE_[il]   = thisLep->energy();

      } // end of loop over jets
      */

    } // end of loop over Higgs candidates  
    
  } // end of looping over lepton types

} // end of Fill()

  //-----------------------------------------------------------------------
void  
ZZTree::SetBranches(){

  AddBranch(&nGoodHCand_, "nGoodHCand");
  AddBranch(&nAllHCand_, "nAllHCand");
  AddBranch(&bestHCand_, "bestHCand");
  AddBranch(&metSig_, "metSig");
  AddBranch(&metSigNoPU_, "metSigNoPU");

  AddBranch(&higgsPt_,"higgsPt");
  AddBranch(&higgsEta_,"higgsEta");
  AddBranch(&higgsPhi_,"higgsPhi");
  AddBranch(&higgsM_,"higgsM");
  AddBranch(&higgsMRefit_,"higgsMRefit");

  AddBranch(&zllPt_,"zllPt");
  AddBranch(&zllEta_,"zllEta");
  AddBranch(&zllPhi_,"zllPhi");
  AddBranch(&zllM_,"zllM");

  AddBranch(&zjjPt_,"zjjPt");
  AddBranch(&zjjEta_,"zjjEta");
  AddBranch(&zjjPhi_,"zjjPhi");
  AddBranch(&zjjM_,"zjjM");
  AddBranch(&zjjMRefit_,"zjjMRefit");

  AddBranch(&heliLD_,"heliLD");
  AddBranch(&heliLDRefit_,"heliLDRefit");

  AddBranch(&nBTags_,"nBTags");
  AddBranch(&lepType_,"lepType");
  AddBranch(&passBit_,"passBit");


  /*

  AddBranch(&higgsPt_,  "higgsPt");
  AddBranch(&higgsEta_, "higgsEta");
  AddBranch(&higgsPhi_, "higgsPhi");
  AddBranch(&higgsM_,   "higgsM");

  AddBranch(&zllPt_,  "zllPt");
  AddBranch(&zllEta_, "zllEta");
  AddBranch(&zllPhi_, "zllPhi");
  AddBranch(&zllM_,   "zllM");

  AddBranch(&zjjPt_,  "zjjPt");
  AddBranch(&zjjEta_, "zjjEta");
  AddBranch(&zjjPhi_, "zjjPhi");
  AddBranch(&zjjM_,   "zjjM");

  int arraySize = sizeof(jetPt_)/sizeof(jetPt_[0]);
  AddBranchArray(arraySize, jetPt_,  "jetPt");
  AddBranchArray(arraySize, jetEta_, "jetEta");
  AddBranchArray(arraySize, jetPhi_, "jetPhi");
  AddBranchArray(arraySize, jetE_,   "jetE");
  AddBranchArray(arraySize, jetRefitPt_,  "jetRefitPt");
  AddBranchArray(arraySize, jetRefitEta_, "jetRefitEta");
  AddBranchArray(arraySize, jetRefitPhi_, "jetRefitPhi");
  AddBranchArray(arraySize, jetRefitE_, "jetRefitE");

  
  AddBranch(&lepType_,   "lepType");

  arraySize = sizeof(lepPt_)/sizeof(lepPt_[0]);
  AddBranchArray(arraySize, lepPt_,  "lepPt");
  AddBranchArray(arraySize, lepEta_, "lepEta");
  AddBranchArray(arraySize, lepPhi_, "lepPhi");
  AddBranchArray(arraySize, lepE_,   "lepE");

  AddBranch(&muID01, "muID01");
  AddBranch(&muID02, "muID02");
  AddBranch(&muID03, "muID03");
  AddBranch(&muID04, "muID04");
  AddBranch(&muID05, "muID05");
  AddBranch(&muID06, "muID06");
  AddBranch(&muID07, "muID07");
  AddBranch(&muID08, "muID08");
  AddBranch(&muID09, "muID09");
  AddBranch(&muID10, "muID10");
  AddBranch(&muID11, "muID11");
  AddBranch(&muID12, "muID12");
  AddBranch(&muID13, "muID13");
  AddBranch(&muID14, "muID14");
  AddBranch(&muID15, "muID15");


  AddBranch(&eID01, "eID01");
  AddBranch(&eID02, "eID02");
  AddBranch(&eID03, "eID03");
  AddBranch(&eID04, "eID04");
  AddBranch(&eID05, "eID05");
  AddBranch(&eID06, "eID06");
  AddBranch(&eID07, "eID07");
  AddBranch(&eID08, "eID08");
  AddBranch(&eID09, "eID09");
  AddBranch(&eID10, "eID10");
  AddBranch(&eID11, "eID11");
  AddBranch(&eID12, "eID12");
  AddBranch(&eID13, "eID13");
  AddBranch(&eID14, "eID14");
  AddBranch(&eID15, "eID15");
  AddBranch(&eID16, "eID16");
  AddBranch(&eID17, "eID17");
  */

}


void  
ZZTree::Clear(){

  nGoodHCand_ = 0;
  nAllHCand_ = 0;
  bestHCand_ = -1;

  metSig_ = -99999.;
  metSigNoPU_= -99999.;

  higgsPt_.clear();
  higgsEta_.clear();
  higgsPhi_.clear();
  higgsM_.clear();
  higgsMRefit_.clear();

  zllPt_.clear();
  zllEta_.clear();
  zllPhi_.clear();
  zllM_.clear();

  zjjPt_.clear();
  zjjEta_.clear();
  zjjPhi_.clear();
  zjjM_.clear();
  zjjMRefit_.clear();

  heliLD_.clear();
  heliLDRefit_.clear();

  nBTags_.clear();
  lepType_.clear();
  passBit_.clear();

  /*

  higgsPt_  = -99999.0;
  higgsEta_ = -99999.0;
  higgsPhi_ = -99999.0;
  higgsM_   = -99999.0;

  zllPt_  = -99999.0;
  zllEta_ = -99999.0;
  zllPhi_ = -99999.0;
  zllM_   = -99999.0;

  zjjPt_  = -99999.0;
  zjjEta_ = -99999.0;
  zjjPhi_ = -99999.0;
  zjjM_   = -99999.0;


  int arraySize = sizeof(jetPt_)/sizeof(jetPt_[0]);

  for(int i=0; i<arraySize;i++)
  {
  jetPt_[i] =-99999.0;
  jetEta_[i]=-99999.0;
  jetPhi_[i]=-99999.0;
  jetE_[i]  =-99999.0;
  jetRefitPt_[i] =-99999.0;
  jetRefitEta_[i]=-99999.0;
  jetRefitPhi_[i]=-99999.0;
  jetRefitE_[i]=-99999.0;

  }
  lepType_ = -1;
  
  arraySize = sizeof(lepPt_)/sizeof(lepPt_[0]);

  for(int i=0; i<arraySize;i++)
  {
  lepPt_[i] =-99999.0;
  lepEta_[i]=-99999.0;
  lepPhi_[i]=-99999.0;
  lepE_[i]  =-99999.0;

  }

  muID01.clear();
  muID02.clear();
  muID03.clear();
  muID04.clear();
  muID05.clear();
  muID06.clear();
  muID07.clear();
  muID08.clear();
  muID09.clear();
  muID10.clear();
  muID11.clear();
  muID12.clear();
  muID13.clear();
  muID14.clear();
  muID15.clear();


  eID01.clear();
  eID02.clear();
  eID03.clear();
  eID04.clear();
  eID05.clear();
  eID06.clear();
  eID07.clear();
  eID08.clear();
  eID09.clear();
  eID10.clear();
  eID11.clear();
  eID12.clear();
  eID13.clear();
  eID14.clear();
  eID15.clear();
  eID16.clear();
  eID17.clear();
  */

}
