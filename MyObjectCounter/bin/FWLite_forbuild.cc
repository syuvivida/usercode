#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include <TH1.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>
#include <TLorentzVector.h>

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

char FILENAME[100];


void parseCommandLine(int argc, char* argv[]) {
  for( int arg=1; arg<argc; arg++) {
   if( strcmp(argv[arg], "-file") == 0) {
      if (argc < arg+1)
        {
	  std::cout << "Missing data file name" << std::endl;
          exit(-1);
        }
      strcpy(FILENAME,argv[arg+1]);
      arg++;
   }
  }
 }


int main(int argc, char* argv[]) 
{
  parseCommandLine(argc,argv);
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();
  TH1F* h1 = new TH1F("h1","",5,-0.5,4.5);
  TH1F* h2 = new TH1F("h2","",5,-0.5,4.5);
  TH1F* h3 = new TH1F("h3","",5,-0.5,4.5);
  TH1F* het = new TH1F("het","",100,0,200);
  TH1F* hdR = new TH1F("hdR","",50,0,0.02);
  TH1F* hdptrel = new TH1F("hdptrel","",50,-1.0,1.0);
  TH1F* hstatus = new TH1F("hstatus","",3,0.5,3.5);
  TH1F* hpdg = new TH1F("hpdg","",61,-30.5,30.5); 

  TFile* infile = TFile::Open(FILENAME);
  std::cout << "Opening " << FILENAME << std::endl;

  fwlite::Event event(infile);

  std::cout << "Inside fwlight" << std::endl;

  unsigned int iEvent = 0;

  for(event.toBegin(); !event.atEnd();
      ++event, ++iEvent)
    {

      //    if(iEvent==1000)break;
      //      if(iEvent%10==0)
	std::cout << "Processing event " << iEvent << std::endl;


      fwlite::Handle<std::vector<pat::Jet> > jets;
      jets.getByLabel(event, "cleanLayer1Jets");

      fwlite::Handle<std::vector<pat::Electron> > electrons;
      electrons.getByLabel(event, "cleanLayer1Electrons");

      h1->Fill(jets->size());
      
      unsigned int njet_afterpat=0;
      unsigned int njet_after=0;
      // electron loop
      for(unsigned int rece=0; rece < electrons->size(); rece++)
	{
	  pat::Electron thisEle = (*electrons)[rece];
	  het->Fill(thisEle.p4().Et());
	  if(!thisEle.genLepton())continue;
	  // only check those electrons which has a matched genp electron
	  TLorentzVector genp_p4(thisEle.genLepton()->p4().Px(),
				 thisEle.genLepton()->p4().Py(),
				 thisEle.genLepton()->p4().Pz(),
				 thisEle.genLepton()->p4().E());

	  TLorentzVector reco_p4(thisEle.p4().Px(),
				 thisEle.p4().Py(),
				 thisEle.p4().Pz(),
				 thisEle.p4().E());
	  int pid = thisEle.genLepton()->pdgId();
	  std::cout << "PDG ID " << pid << "\t";
	  
	  std::cout << "Generator lepton charge " << thisEle.genLepton()->charge() << " ("
		    << thisEle.genLepton()->p4().Px() << ","
		    << thisEle.genLepton()->p4().Py() << ","
		    << thisEle.genLepton()->p4().Pz() << ","
		    << thisEle.genLepton()->p4().E() << ")\t";

	  std::cout << "Reconstructed lepton charge " << thisEle.charge() << " ("
		    << thisEle.p4().Px() << ","
		    << thisEle.p4().Py() << ","
		    << thisEle.p4().Pz() << ","
		    << thisEle.p4().E() << ")" << std::endl;

	  hpdg->Fill(pid);
	  hstatus->Fill(thisEle.genLepton()->status());

	  Double_t dR = genp_p4.DeltaR(reco_p4);
	  Double_t dptrel = genp_p4.Pt()>1e-6? 
	    (reco_p4.Pt()-genp_p4.Pt())/genp_p4.Pt() : 9999.0;
	  if(fabs(dptrel)>0.5 || dR>0.5)std::cout << "dptrel = " << dptrel << ", dR=" << dR << std::endl; 
	  hdR->Fill(dR);
	  hdptrel->Fill(dptrel);
	}
      // jet loop
      for(unsigned int i=0; i < jets->size(); i++)
	{
	  bool pat_matched = false;
	  pat::Jet thisJet = (*jets)[i];
	  pat_matched = thisJet.hasOverlaps("electrons");
	  //	  std::cout << "pat_matched = " << pat_matched << std::endl;
	  bool matched=false;
	  for(unsigned int j=0; j < electrons->size(); j++)
	    {
	      pat::Electron thisEle = (*electrons)[j];

	      double dR = reco::deltaR(thisJet.eta(),
				 thisJet.phi(),
				 thisEle.eta(),
				 thisEle.phi());
	      
	      if(dR < 0.2)
		matched=true;

	    } // end of loop over electrons
	  if(!pat_matched)njet_afterpat++;
	  if(!matched)njet_after++;

	} // end of loop over jets

      h2->Fill(njet_after);
      h3->Fill(njet_afterpat);

    } // end of loop over events

  infile->Close();

  TFile outFile("histo.root","recreate");
  h1->Write();
  h2->Write();
  h3->Write();
  het->Write();
  hdR->Write();
  hdptrel->Write();
  hstatus->Write();
  hpdg->Write();
  outFile.Close();
  delete h1;
  delete h2;
  delete h3;
  return 0;

}
