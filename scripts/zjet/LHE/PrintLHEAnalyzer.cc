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
#include <TH1F.h>
#include <TFile.h>
using namespace std;
using namespace edm;
using namespace lhef;


class PrintLHEAnalyzer : public EDAnalyzer {
private: 
  TFile * output;
  TH1F* h_mb_ini;
  TH1F* h_mc_ini;
  TH1F* h_mb_out;
  TH1F* h_mc_out;

public:
  explicit PrintLHEAnalyzer( const edm::ParameterSet & cfg ) : 
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
    
    h_mb_ini    = new TH1F("h_mb_ini", "",200,-10.0,10.0);
    h_mb_ini    ->SetXTitle("M_{b} [GeV/c^{2}");
    h_mb_ini    ->Sumw2();

    h_mc_ini    = new TH1F("h_mc_ini", "",200,-10.0,10.0);
    h_mc_ini    ->SetXTitle("M_{c} [GeV/c^{2}");
    h_mc_ini    ->Sumw2();

    h_mb_out    = new TH1F("h_mb_out", "",200,-10.0,10.0);
    h_mb_out    ->SetXTitle("M_{b} [GeV/c^{2}");
    h_mb_out    ->Sumw2();

    h_mc_out    = new TH1F("h_mc_out", "",200,-10.0,10.0);
    h_mc_out    ->SetXTitle("M_{c} [GeV/c^{2}");
    h_mc_out    ->Sumw2();

  }

  ~PrintLHEAnalyzer(){
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
    
    for ( unsigned int icount = 0 ; icount < (unsigned int)nup_; icount++ ) {

      int PID    = idup_[icount];
      int status = istup_[icount];
      double px = (pup_[icount])[0];
      double py = (pup_[icount])[1];
      double pz = (pup_[icount])[2];
      double e  = (pup_[icount])[3];
      double m  = (e*e-px*px-py*py-pz*pz)<0? -99999: 
	sqrt(e*e-px*px-py*py-pz*pz);
      
      if(abs(PID)==5 || abs(PID)==4)
	cout << "ID: " << std::setw(14) << std::fixed << PID
	     << "\t Px: " << std::setw(14) << std::fixed << px 
	     << "\t Py: " << std::setw(14) << std::fixed << py 
	     << "\t Pz: " << std::setw(14) << std::fixed << pz 
	     << "\t E: " << std::setw(14) << std::fixed << e 
	     << "\t M: " << std::setw(14) << std::fixed 
	     << m << endl;

      if(status==-1 && abs(PID)==5)
	h_mb_ini->Fill(m);
      else if(status==-1 && abs(PID)==4)
	h_mc_ini->Fill(m);      
      else if(status== 1 && abs(PID)==5)
	h_mb_out->Fill(m);
      else if(status== 1 && abs(PID)==4)
	h_mc_out->Fill(m);

      
    } // end of loop over particles

  }    
  void endJob(){
    h_mb_ini->Write();
    h_mc_ini->Write();
    h_mb_out->Write();
    h_mc_out->Write();
    output->Write();
    output->Close();
  }   

  InputTag src_;
  std::string fileName_;
};

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( PrintLHEAnalyzer );


