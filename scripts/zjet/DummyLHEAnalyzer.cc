#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <iostream>
#include <vector>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TTree.h>

using namespace std;
using namespace edm;
using namespace lhef;

const int LEPPID=13;
    
const double minZPt  = 40.0;
const double minPartonPt=30.0;
const double maxPartonEta=2.4;
const double mindR=0.5;

const double minElePt = 20.0;
const double minEleBarrelEta = 0.0;
const double maxEleBarrelEta = 1.442;

const double minEleEndcapEta = 1.566;
const double maxEleEndcapEta = 2.1;

const double minMuoPt = 20.0;
const double maxMuoEta = 2.1;

const double minM = 76.0;
const double maxM =106.0;


class DummyLHEAnalyzer : public EDAnalyzer {
private: 
  bool dumpLHE_;
  bool checkPDG_;

  std::vector<double> genParId_;
  std::vector<double> genParSt_;
  std::vector<double> genParPx_;
  std::vector<double> genParPy_;
  std::vector<double> genParPz_;
  std::vector<double> genParE_;

  TTree *tree;
  
  TH1D* h_mZ;
  TH1D* h_zpt;
  TH1D* h_jetpt;
  TH1D* h_zy;
  TH1D* h_jety; 
  TH1D* h_yB; 
  TH1D* h_ystar;

  

public:
  explicit DummyLHEAnalyzer( const ParameterSet & cfg ) : 
    src_( cfg.getParameter<InputTag>( "src" ) )
  {
  }
  
  void beginJob(){
    
    Service<TFileService> fs;
    h_mZ = fs->make<TH1D>("h_mZ", "",200,20.0,220.0);
    h_mZ->SetXTitle("M_{ll} [GeV/c^{2}");
    h_mZ->Sumw2();

    h_zpt  = fs->make<TH1D>("h_zpt","",40,0,400);
    h_zpt->SetXTitle("p_{T}(Z) [GeV]");
    h_zpt->Sumw2();

    h_jetpt = fs->make<TH1D>("h_jetpt","",40,0,400);
    h_jetpt->SetXTitle("p_{T}(jet) [GeV]");
    h_jetpt->Sumw2();

    h_zy = fs->make<TH1D>("h_zy","",15,0.0,3.0); 
    h_zy->SetXTitle("|y_{Z}|");
    h_zy->Sumw2();

    h_jety= fs->make<TH1D>("h_jety","",15,0.0,3.0); 
    h_jety->SetXTitle("|y_{jet}|");
    h_jety->Sumw2();

    h_yB  = fs->make<TH1D>("h_yB","",15,0.0,3.0); 
    h_yB->SetXTitle("0.5|y_{Z}+y_{jet}|");
    h_yB->Sumw2();

    h_ystar  = fs->make<TH1D>("h_ystar","",15,0.0,3.0); 
    h_ystar->SetXTitle("0.5|y_{Z}-y_{jet}|");
    h_ystar->Sumw2();
  
    tree = fs->make<TTree>("tree", "LHE tree");
    
    tree->Branch("genParId_",&genParId_);
    tree->Branch("genParSt_",&genParSt_);
    tree->Branch("genParPx_",&genParPx_);
    tree->Branch("genParPy_",&genParPy_);
    tree->Branch("genParPz_",&genParPz_);
    tree->Branch("genParE_",&genParE_);


  }

  void Clear(){

    genParId_.clear();
    genParSt_.clear();
    genParPx_.clear();
    genParPy_.clear();
    genParPz_.clear();
    genParE_.clear();
  }

private:
  void analyze( const Event & iEvent, const EventSetup & iSetup ) {

    Clear();

    Handle<LHEEventProduct> evt;
    iEvent.getByLabel( src_, evt );

    const lhef::HEPEUP hepeup_ = evt->hepeup();

    const int nup_ = hepeup_.NUP; 
    const std::vector<int> idup_ = hepeup_.IDUP;
    const std::vector<int> istup_ = hepeup_.ISTUP;

    const std::vector<lhef::HEPEUP::FiveVector> pup_ = hepeup_.PUP;

    std::cout << "Number of particles = " << nup_ << std::endl;

    TLorentzVector lep_p(0,0,0,0);
    TLorentzVector lep_m(0,0,0,0);

    for ( unsigned int icount = 0 ; icount < (unsigned int)nup_; icount++ ) {

//       std::cout << "# " << icount << " " << idup_[icount] << " " 
// 		<< istup_[icount] << " " 
// 		<< hepeup_.MOTHUP[icount].first << " " << hepeup_.MOTHUP[icount].second << " "
// 		<< (pup_[icount])[0] << " " << (pup_[icount])[1] << " " << (pup_[icount])[2] << " " << (pup_[icount])[3] << " " << (pup_[icount])[4] << std::endl;
      int PID    = idup_[icount];
      int status = istup_[icount];
      double px = (pup_[icount])[0];
      double py = (pup_[icount])[1];
      double pz = (pup_[icount])[2];
      double e  = (pup_[icount])[3];

      genParId_.push_back(PID);
      genParSt_.push_back(status);
      genParPx_.push_back(px);
      genParPy_.push_back(py);
      genParPz_.push_back(pz);
      genParE_.push_back(e);
      
      if(status!=1)continue;
      if(PID==LEPPID)
	lep_m.SetPxPyPzE(px,py,pz,e);
      else if(PID== -LEPPID)
	lep_p.SetPxPyPzE(px,py,pz,e);
      
    } // end of loop over particles

    if(lep_m.Pt() < minMuoPt)return;
    if(lep_p.Pt() < minMuoPt)return;

    if(fabs(lep_m.Eta()) > maxMuoEta)return;
    if(fabs(lep_p.Eta()) > maxMuoEta)return;

    TLorentzVector l4_z = lep_m+lep_p;

    double mz = l4_z.M();   
    h_mZ->Fill(mz);

    if(mz < minM || mz > maxM)return;

    double zpt = l4_z.Pt();
    double zy  = l4_z.Rapidity();

    if(zpt < minZPt)return;

    cout << "zpt=" << zpt << endl;

    TLorentzVector l4_j(0,0,0,0);
    unsigned int nPartons = 0;
    double maxPartonPt = -9999;
    
    for ( unsigned int icount = 0 ; icount < (unsigned int)nup_; icount++ ) {

      int PID    = abs(idup_[icount]);
      int status = istup_[icount];
      double px = (pup_[icount])[0];
      double py = (pup_[icount])[1];
      double pz = (pup_[icount])[2];
      double e  = (pup_[icount])[3];

      if(status!=1)continue;

      if(PID > 5 && PID!=21)continue;

      
      TLorentzVector thisParton_l4(0,0,0,0);
      thisParton_l4.SetPxPyPzE(px,py,pz,e);

      double pt = thisParton_l4.Pt();
      double eta= thisParton_l4.Eta();

      if(pt < minPartonPt)continue;
      if(fabs(eta) > maxPartonEta)continue;

      double dr_ep = lep_p.DeltaR(thisParton_l4);
      double dr_em = lep_m.DeltaR(thisParton_l4);

      
      if(dr_ep < mindR)continue;
      if(dr_em < mindR)continue;
      nPartons++; 

      if(pt > maxPartonPt)
	{
	  maxPartonPt = pt;

	  l4_j.SetPxPyPzE(px,py,pz,e);
	}
	

    } // end of loop over particles

    
    
    if(nPartons!=1)return;

    double jetpt = l4_j.Pt();
    double jety  = l4_j.Rapidity();
    double  yB   = 0.5*(zy+jety);
    double  ystar= 0.5*(zy-jety);

    h_zpt->Fill(zpt);
    h_jetpt->Fill(jetpt);
    h_zy->Fill(fabs(zy));
    h_jety->Fill(fabs(jety));
    h_yB->Fill(fabs(yB));
    h_ystar->Fill(fabs(ystar));


    tree->Fill();

  }
  InputTag src_;
}; 

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( DummyLHEAnalyzer );


