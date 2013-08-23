#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include <iostream>
#include <string>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TLorentzVector.h>


#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "TLorentzVector.h"
#include "TVector3.h"

using namespace std;
using namespace edm;
using namespace lhef;


class HTransformToHelicityFrame {

public:

  HTransformToHelicityFrame(const Double_t& energy, const Double_t& mass=0.938272){
    HbeamPartMass=mass;
    SetCmsEnergy(energy);
  }

  ~HTransformToHelicityFrame() { };

  void SetCmsEnergy(const Double_t& cms_energy)
  {
    if(HcmsEnergy!=cms_energy) {
      HcmsEnergy = cms_energy;
      //       cout << "HTransformToHelicityFrame::SetCmsEnergy: Set center of mass energy: " << HcmsEnergy << " GeV" <<  endl;
      //      cout << "HTransformToHelicityFrame::SetCmsEnergy: Mass of beam particles:    " << HbeamPartMass << endl;
      
      PF_org.SetPxPyPzE(0.,0., sqrt(0.25*HcmsEnergy*HcmsEnergy-HbeamPartMass*HbeamPartMass),0.5*HcmsEnergy);
      PW_org.SetPxPyPzE(0.,0.,-sqrt(0.25*HcmsEnergy*HcmsEnergy-HbeamPartMass*HbeamPartMass),0.5*HcmsEnergy);
    }
  }

  void TransformToHelicityFrame(TLorentzVector muplus,TLorentzVector muminus){

    Wv= muplus+muminus;// this is the Z boson 4vector
    b = Wv.BoostVector();
    muplus.Boost(-b);
    muminus.Boost(-b);
  
    PF = PF_org;
    PW = PW_org;

    PF.Boost(-b);
    PW.Boost(-b);

    Htheta_star = muminus.Vect().Angle(Wv.Vect());
    Hcos_theta_star = TMath::Cos(Htheta_star);
    // choose what to call proton and what anti-proton
    if(Wv.Angle(PF.Vect())>Wv.Angle(PW.Vect()))
      {
	Hproton_Z_plane_3v = Wv.Vect().Cross(PF.Vect()).Unit();
      }
    else
      {
	Hproton_Z_plane_3v = Wv.Vect().Cross(PW.Vect()).Unit();
      }
    HMu_Z_plane_3v = (Wv.Vect().Cross(muminus.Vect())).Unit();
    Hphi_star = HMu_Z_plane_3v.Angle(Hproton_Z_plane_3v);
    if(Hproton_Z_plane_3v.Dot(muminus.Vect()) < 0.0) {
      Hphi_star = -Hphi_star;
    }
  }


  const Double_t& GetTheta() { return Htheta_star; }
  
  const Double_t& GetCosTheta() { return Hcos_theta_star; }

  const Double_t& GetPhi() { return Hphi_star; }

private:

  Double_t HcmsEnergy;
  Double_t HbeamPartMass;

  Double_t Htheta_star;
  Double_t Hcos_theta_star;
  Double_t Hphi_star;

  TLorentzVector PF;
  TLorentzVector PW;
  TLorentzVector PF_org;
  TLorentzVector PW_org;
  TLorentzVector Wv;
  TVector3 b;
  TVector3 Hproton_Z_plane_3v;
  TVector3 HMu_Z_plane_3v;

}; // transform



class DummyLHEAnalyzer : public EDAnalyzer {
private: 

  TFile * output;
  TH1D* h_phistar;
  TH1D* h_thetastar;
  TH2D* h_phithetastar;
  TH1D* h_ms;
  TH1D* h_mZ;
  TH1D* h_yZ;
  TH1D* h_phiZ;
  TH1D* h_thetaZ;
  TH1D* h_ptZ;
  TH1D* h_jetpt;
  TH1D* h_zy;
  TH1D* h_jety; 
  TH1D* h_yB; 
  TH1D* h_ystar;

  TH1D* h_egluon;
  TH1D* h_equark;

  TH2D* h_ini;
  TH2D* h_fin;
  TH2D* h_flep;
  TH1D* h_nparton;

public:
  explicit DummyLHEAnalyzer( const edm::ParameterSet & cfg ) : 
    src_( cfg.getParameter<InputTag>( "src" )),
    fileName_(cfg.getUntrackedParameter<std::string>("histoutputFile"))
  {
  }
  void beginJob(){
    //       cout << "beginJob" <<endl;
    //     edm::ESHandle<TrackerGeometry> tkgeom;
    //    c.get<TrackerDigiGeometryRecord>().get( tkgeom );
    //     tracker=&(* tkgeom);
    output = new TFile(fileName_.data(), "RECREATE");


    h_phistar = new TH1D("h_phistar","",50,-TMath::Pi(),TMath::Pi());
    h_phistar->SetXTitle("#phi^{*}");
    h_thetastar = new TH1D("h_thetastar","",50,-1,1);
    h_thetastar->SetXTitle("#theta^{*}");
    h_phithetastar = new TH2D("h_phithetastar","",50,-TMath::Pi(),TMath::Pi(),
			      50,-1,-1);

    h_phithetastar->SetXTitle("#phi^{*}");
    h_phithetastar->SetYTitle("#theta^{*}");

    h_ini = new TH2D("h_ini","Initial state partons",
		     28,-6.5,21.5,28,-6.5,21.5);    
    h_ini->SetXTitle("Particle ID 1");
    h_ini->SetYTitle("Particle ID 2");

    h_nparton = new TH1D("h_nparton","",6,-0.5,5.5);
    h_nparton->SetXTitle("Number of partons");

    h_fin = new TH2D("h_fin","Final state partons",28,-6.5,21.5,28,-6.5,21.5);
    h_fin->SetXTitle("Particle ID 1");
    h_fin->SetYTitle("Particle ID 2");

    h_flep = new TH2D("h_flep","Final state leptons",33,-16.5,16.5,33,-16.5,16.5);
    h_flep->SetXTitle("Particle ID 1");
    h_flep->SetYTitle("Particle ID 2");

    h_equark = new TH1D("h_equark","",500,0,500);
    h_equark ->SetXTitle("Energy of quark [GeV]");
 
    h_egluon = new TH1D("h_egluon","",500,0,500);
    h_egluon ->SetXTitle("Energy of gluon [GeV]");

    h_ptZ   = new TH1D("h_ptZ","",200,0,200);
    h_ptZ   ->SetXTitle("p_{T} of Z [GeV]");
    h_ptZ   ->Sumw2();
    
    h_yZ   = new TH1D("h_yZ","",200,-2.0,2.0);
    h_yZ   ->SetXTitle("Rapidity of Z");
    h_yZ   ->Sumw2();

    h_phiZ   = new TH1D("h_phiZ","",200,-TMath::Pi(),TMath::Pi());
    h_phiZ   ->SetXTitle("phi of Z");
    h_phiZ   ->Sumw2();

    h_thetaZ = new TH1D("h_thetaZ","",200,0.0,TMath::Pi());
    h_thetaZ  ->SetXTitle("#theta of Z");
    h_thetaZ  ->Sumw2();


    h_mZ    = new TH1D("h_mZ", "",240,60.0,120.0);
    h_mZ->SetXTitle("M_{ll} [GeV/c^{2}]");
    h_mZ->SetYTitle("Candidates per 0.25 GeV");
    h_mZ->Sumw2();

    h_ms    = new TH1D("h_ms", "",200,20.0,220.0);
    h_ms->SetXTitle("Scale [GeV]");
    h_ms->Sumw2();


    h_jetpt = new TH1D("h_jetpt","",40,0,400);
    h_jetpt->SetXTitle("p_{T}(jet) [GeV]");
    h_jetpt->Sumw2();

    h_zy    = new TH1D("h_zy","",15,0.0,3.0); 
    h_zy->SetXTitle("|y_{Z}|");
    h_zy->Sumw2();

    h_jety  = new TH1D("h_jety","",15,0.0,3.0); 
    h_jety->SetXTitle("|y_{jet}|");
    h_jety->Sumw2();

    h_yB    = new TH1D("h_yB","",15,0.0,3.0); 
    h_yB->SetXTitle("0.5|y_{Z}+y_{jet}|");
    h_yB->Sumw2();

    h_ystar = new TH1D("h_ystar","",15,0.0,3.0); 
    h_ystar->SetXTitle("0.5|y_{Z}-y_{jet}|");
    h_ystar->Sumw2();
  }

  ~DummyLHEAnalyzer(){
    delete output;
  }

private:
  void analyze( const Event & iEvent, const EventSetup & iSetup ) {

    Handle<LHEEventProduct> evt;
    iEvent.getByLabel( src_, evt );

    const lhef::HEPEUP hepeup_ = evt->hepeup();
    const int nup_ = hepeup_.NUP; 
    const std::vector<int> idup_ = hepeup_.IDUP;
    const std::vector<lhef::HEPEUP::FiveVector> pup_ = hepeup_.PUP;
    const std::vector<int> istup_ = hepeup_.ISTUP;
    double scale = hepeup_.SCALUP;
    
    h_ms->Fill(scale);

    TLorentzVector l4_z_real(0,0,0,0);
    TLorentzVector l4_mup(0,0,0,0);
    TLorentzVector l4_mum(0,0,0,0);

    double ine[2]={-999,-999};
    int inPID[2]={-999,-999};
    int outPID[2]={-999,-999};
    int nParton=0;
    int lepPID[2]={-999,-999};
    
    for ( unsigned int icount = 0 ; icount < (unsigned int)nup_; icount++ ) {

      int PID    = idup_[icount];
      int status = istup_[icount];
      double px = (pup_[icount])[0];
      double py = (pup_[icount])[1];
      double pz = (pup_[icount])[2];
      double e  = (pup_[icount])[3];

      if(PID==23)
	l4_z_real.SetPxPyPzE(px,py,pz,e);

      if(status==-1 && inPID[0]<-998){
	inPID[0]=PID;
	ine[0] = e;
      }
      else if(status==-1 && inPID[1]<-998){
	inPID[1]=PID;
	ine[1] = e;
      }

      if(status==1 && (abs(PID)==11 || abs(PID)==13 || abs(PID)==15))
	{
	  if(PID>0){
	    lepPID[0]=PID;
	    l4_mum.SetPxPyPzE(px,py,pz,e);
	  }
	  else{
	    lepPID[1]=PID;
	    l4_mup.SetPxPyPzE(px,py,pz,e);
	  }
	}

      else if(status==1 && (abs(PID)<=6 || abs(PID)==21))
	{
	  nParton++;
	  if(outPID[0]<-998)outPID[0]=PID;
	  else if(outPID[1]<-998)outPID[1]=PID;
	}
      
      
    }// end of loop over particles


    h_nparton->Fill(nParton);

    if(inPID[0]==21 && abs(inPID[1])<=6)
      {
	h_equark->Fill(ine[1]);
	h_egluon->Fill(ine[0]);
      }
    else if(abs(inPID[0])<6 && inPID[1]==21)
      {
	h_equark->Fill(ine[0]);
	h_egluon->Fill(ine[1]);
      }

    if(inPID[0]>-999 && inPID[1]>-999)
      h_ini->Fill(inPID[0],inPID[1]);

    if(lepPID[0]>-999 && lepPID[1]>-999 && abs(lepPID[0])==abs(lepPID[1]) )
      {
	
	double mz = l4_z_real.M();   
	h_mZ->Fill(mz);
	h_yZ->Fill(l4_z_real.Rapidity());
	h_flep->Fill(lepPID[0],lepPID[1]);
	
	HTransformToHelicityFrame myHTransformToHelicityFrame(centerEnergy_);
	myHTransformToHelicityFrame.TransformToHelicityFrame(l4_mum,l4_mup);

	double phiStar = myHTransformToHelicityFrame.GetPhi();
	double cosThetaStar = cos(myHTransformToHelicityFrame.GetTheta());
	h_thetastar->Fill(cosThetaStar);

	if(nParton>0){
	  h_phiZ->Fill(l4_z_real.Phi());
	  h_thetaZ->Fill(l4_z_real.Theta());
	  h_ptZ->Fill(l4_z_real.Pt());	
	  h_phistar->Fill(phiStar);
	  h_phithetastar->Fill(phiStar,cosThetaStar);
	}

      }
  

    if(outPID[0]>-999 && outPID[1]>-999)
      h_fin->Fill(outPID[0],outPID[1]);



  } // end of function


    
      
  void endJob(){


    h_phistar->Write();
    h_thetastar->Write();
    h_phithetastar->Write();


    h_ini->Write();
    h_nparton->Write();
    h_fin->Write();
    h_flep->Write();

    h_egluon->Write();
    h_equark->Write();

    h_yZ->Write();
    h_phiZ->Write();
    h_thetaZ->Write();
    h_ptZ->Write();
    h_mZ->Write();
    h_ms->Write();
    h_jetpt->Write();
    h_zy->Write();
    h_jety->Write();
    h_yB->Write();
    h_ystar->Write();
    output->Write();
    output->Close();
  }   

  void beginRun(edm::Run const& iRun, edm::EventSetup const& es){

    Handle<LHERunInfoProduct> run;
    iRun.getByLabel( src_, run );
    
    const lhef::HEPRUP thisHeprup_ = run->heprup();
    centerEnergy_ = thisHeprup_.EBMUP.first + thisHeprup_.EBMUP.second;

    /*
      std::cout << "HEPRUP \n" << std::endl;
      std::cout << "IDBMUP " << std::setw(14) << std::fixed << thisHeprup_.IDBMUP.first 
      << std::setw(14) << std::fixed << thisHeprup_.IDBMUP.second << std::endl; 
      std::cout << "EBMUP  " << std::setw(14) << std::fixed << thisHeprup_.EBMUP.first 
      << std::setw(14) << std::fixed << thisHeprup_.EBMUP.second << std::endl; 
      std::cout << "PDFGUP " << std::setw(14) << std::fixed << thisHeprup_.PDFGUP.first 
      << std::setw(14) << std::fixed << thisHeprup_.PDFGUP.second << std::endl; 
      std::cout << "PDFSUP " << std::setw(14) << std::fixed << thisHeprup_.PDFSUP.first 
      << std::setw(14) << std::fixed << thisHeprup_.PDFSUP.second << std::endl; 
      std::cout << "IDWTUP " << std::setw(14) << std::fixed << thisHeprup_.IDWTUP << std::endl; 
      std::cout << "NPRUP  " << std::setw(14) << std::fixed << thisHeprup_.NPRUP << std::endl; 
      std::cout << "        XSECUP " << std::setw(14) << std::fixed 
      << "        XERRUP " << std::setw(14) << std::fixed 
      << "        XMAXUP " << std::setw(14) << std::fixed 
      << "        LPRUP  " << std::setw(14) << std::fixed << std::endl;
      for ( unsigned int iSize = 0 ; iSize < thisHeprup_.XSECUP.size() ; iSize++ ) {
      std::cout  << std::setw(14) << std::fixed << thisHeprup_.XSECUP[iSize]
      << std::setw(14) << std::fixed << thisHeprup_.XERRUP[iSize]
      << std::setw(14) << std::fixed << thisHeprup_.XMAXUP[iSize]
      << std::setw(14) << std::fixed << thisHeprup_.LPRUP[iSize] 
      << std::endl;
      }
      std::cout << " " << std::endl;
    */
  }

  InputTag src_;
  std::string fileName_;
  double centerEnergy_;
};

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( DummyLHEAnalyzer );


