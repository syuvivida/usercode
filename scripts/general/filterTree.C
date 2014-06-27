#define filterTree_cxx
#include "filterTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>

// skim mode: 0 for both electron and muon, 
//            1 for electron
//            2 for muon

void filterTree::Loop(int skimMode)
{
  if (fChain == 0) return;
  cout << "Input file is " << inputFile_ << endl;
  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries = " << nentries << endl;

  if(skimMode>2 || skimMode < 0) 
    {
      cout << "Skim mode must be 0, 1, or 2. Program exited!" << endl;
      return;
    }

  std::string endfix = "_filteredtree.root";
    
  switch(skimMode){
  
  case 0:
    cout << "filtering both double electron and double muon" << endl;
    break;
  case 1: 
    endfix = "_electron.root";
    cout << "filtering double elctron" << endl;
    break;
  case 2: 
    endfix = "_muon.root";
    cout << "filtering double muon" << endl;
    break;
  }


  // rename the output file
  std::string remword=".root";
  size_t pos = inputFile_.find(remword);
  std::string forOutput = inputFile_;  

  if(pos!= std::string::npos)
    forOutput.swap(forOutput.erase(pos,remword.length()));   


  std::string outputFile = forOutput + endfix;
  // now open new root file
  TFile* newfile_data = new TFile(outputFile.data(),"recreate");
  cout << "Output file " << outputFile << endl;

  // clone tree
  TTree* newtree = fChain->CloneTree(0);
  newtree->SetMaxTreeSize(4000000000);
  cout << "Saving tree" << endl;

  Long64_t nPassEvt=0;
   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if (jentry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jentry + 1, nentries);
    
      // look for two electrons with pt > 20 GeV, |eta| <2.5
      Int_t    nGoodEle         = 0;

      for(Int_t i=0; i < nEle; i++)
	{
	  if(elePt->at(i)>20.0 && fabs(eleSCEta->at(i)) < 1.4442)
	    nGoodEle++;
	  else if(elePt->at(i)>20.0 && fabs(eleSCEta->at(i)) > 1.566 && fabs(eleSCEta->at(i))<2.5)
	  nGoodEle++;
	}

    // look for two muons with pt > 20 GeV, |eta| <2.4
      Int_t   nGoodMuo        = 0;
      for(Int_t i=0; i < nMu; i++)
	{
	  if(muPt->at(i)>20.0 && fabs(muEta->at(i))<2.4 && (muType->at(i) & 4))
	    nGoodMuo++;
	}

      bool Pass = false;
      switch(skimMode){
      case 0:
	if(nGoodEle >=2 || nGoodMuo >=2)Pass=true;
	else Pass = false;
	break;
      case 1:
	if(nGoodEle >=2)Pass=true;
	else Pass = false;
	break;
      case 2:
	if(nGoodMuo >=2)Pass=true;
      else Pass = false;
	break;

      }

    if(!Pass)continue;

    newtree->Fill();
    nPassEvt++;

   } // end of loop over event entries
 
  newtree->Print();
  newtree->AutoSave();
  delete newfile_data;


  cout << "Number of passed events = " << nPassEvt << endl;
}
