#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TChain.h>
#include <TBranch.h>
#include <iostream.h>
#include <fstream.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystemDirectory.h>
#include <TProfile.h>					
#include <TList.h>
#include <TF1.h>
#include <vector>
#include <TLorentzVector.h>
#include "NCUNtupleClasses.hh"
// #include "TrigEff_AsymmetryErrors.C"
// #include <TGraphAsymmErrors.h>

using namespace SY_NT;
using namespace std;

void TrigEff_withTriggerBit(std::string dirname)
{

  TChain* pho = new TChain("ncuAnalyzerKit/EventTree");

  TSystemDirectory *base = new TSystemDirectory("root","root");
  std::string filename = "/mc/NCUHEP/"+ dirname;

  base->SetDirectory(filename.data());
  TList *listOfFiles = base->GetListOfFiles();
  TIter fileIt(listOfFiles);
  TFile *fileH = new TFile();
  int nfile=0;
  while(fileH = (TFile*)fileIt()) {
    std::string fileN = fileH->GetName();
    std::string baseString = "ncutree";
    if( fileH->IsFolder())  continue;
    if(fileN.find(baseString) == std::string::npos)continue;
    cout << fileN.data() << endl;
    nfile++;
    pho->Add(fileN.data());
  }

  std::cout << "Opened " << nfile << " files" << std::endl;
  

  Long64_t nentries = (Long64_t)pho->GetEntries();

  TH1D* h_getdeno = new TH1D("h_getdeno","Reconstructed and matched photon Et before photon"
			     " trigger cuts", 100,0,200);
  TH1D* h_getnumr = new TH1D("h_getnumr","Reconstructed and matched photon Et after photon"
			     " trigger cuts", 100,0,200);


  TH1D* h_recgetdeno = new TH1D("h_recgetdeno","Reconstructed photon Et before photon"
			     " trigger cuts", 100,0,200);
  TH1D* h_recgetnumr = new TH1D("h_recgetnumr","Reconstructed photon Et after photon"
			     " trigger cuts", 100,0,200);


  TH1D* h_gengetdeno = new TH1D("h_gengetdeno","Generated photon Et before photon"
			     " trigger cuts", 100,0,200);
  TH1D* h_gengetnumr = new TH1D("h_gengetnumr","Generated photon Et after photon"
			     " trigger cuts", 100,0,200);




  // loop over event entries
  for( Long64_t j=0; j< nentries; j++){

    // this line is very important!!
    pho->GetEntry(j);

    std::vector<PhoObj> PhoVect;
    LoadPhos(pho, PhoVect);
    EvtObj thisEvent;
    LoadEvt(pho, thisEvent);
    std::vector<GenObj> GenVect;
    LoadGens(pho, GenVect);


    int* HLT = thisEvent.HLT();


//     bool passL1EGTrigger = true;
    bool passMuonTrigger = (HLT[25]==1);
    bool passL1EGTrigger = (HLT[38]==1);
    bool passPhotonTrigger = (HLT[56]==1);

    //    bool passPhotonTrigger = (HLT[56]==1);
//     bool passL1EGTrigger = (HLT[54]==1);
//     bool passPhotonTrigger = (HLT[58]==1);

    bool hasLoosePhoton = false;

    for(int i=0; i<PhoVect.size(); i++)
      {
	PhoObj thisPho = PhoVect[i];
	double thisEt = thisPho.Et();
	
	// has to contain a loose photon
 	if(!thisPho.IsLoosePhoton())continue;
	hasLoosePhoton = true;

	TLorentzVector rec_pho(0,0,0,0);
	rec_pho.SetPtEtaPhiE(
			     thisPho.Et(),
			     thisPho.Eta(),
			     thisPho.Phi(),
			     thisPho.E());
	

	// matched to a genp photon
	bool isMatched = thisPho.GenIndex()>=0 
	  && thisPho.GenMomPID()==22;
	
	TLorentzVector  gen_pho(0,0,0,0);
	// plotted vs genp photon Et

	// look for matched genp object
	for(int ig=0; ig<GenVect.size(); ig++)
	  {
	    GenObj thisGenObj = GenVect[ig];
	    
	    if(thisGenObj.PID()!=22 ||  abs(thisGenObj.MomPID())!=22)continue;

	    gen_pho.SetPtEtaPhiM(
				 thisGenObj.Pt(),
				 thisGenObj.Eta(),
				 thisGenObj.Phi(),
				 thisGenObj.Mass());

	    double dR = rec_pho.DeltaR(gen_pho);
		
	    if(dR<0.5)
	      break;
	    
	  } // find matched genp object
	    

	// denominator
// 	if(hasLoosePhoton && passL1EGTrigger)
	if(hasLoosePhoton && passMuonTrigger)
	  {
	    h_recgetdeno->Fill(thisEt);
	    if(isMatched)
	      h_getdeno->Fill(thisEt);
	    if(gen_pho.E()>0)
	      h_gengetdeno->Fill(gen_pho.Et());

	    // numerator
	    if(passPhotonTrigger)
	      {
		h_recgetnumr->Fill(thisEt);
		if(isMatched)
		  h_getnumr->Fill(thisEt);
		if(gen_pho.E()>0)
		  h_gengetnumr->Fill(gen_pho.Et());
	      }

	    break;
	  
	  } // if there is a loose photon that passes base trigger	    
	
      } // end of loop over photon objects

    
  } // end of loop over event entries

//   TGraphAsymmErrors *heff = MyDivide(h_getdeno,h_getnumr);  
//   heff->GetXaxis()->SetTitle("Offline photon E_{T} [GeV]");
//   heff->GetYaxis()->SetTitle("Trigger efficiency");                               
//   heff->SetTitle("");
//   heff->Draw("ap");

  std::string histoFile = "/mc/NCUHEP/scripts/"+ dirname + "_histo.root";
  TFile* outFile = new TFile(histoFile.data(),"recreate");

  h_getdeno->Write();
  h_getnumr->Write();
  h_recgetdeno->Write();
  h_recgetnumr->Write();
  h_gengetdeno->Write();
  h_gengetnumr->Write();

  outFile->Close();
  
  

}


