#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <string>
#include <TROOT.h>
#include <algorithm>
#include <vector>
#include "untuplizer.h"

using namespace std;

// skim mode: 0 for both electron and muon, 
//            1 for electron
//            2 for muon

void skimTree(std::string inputFile_, int skimMode=0)
{
  TreeReader data(inputFile_.data()); // v5.3.12
  cout << "Input file is " << inputFile_ << endl;
  Long64_t nentries = data.GetEntriesFast();
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
  TTree* newtree = data.GetTree()->CloneTree(0);
  newtree->SetMaxTreeSize(6000000000);
  cout << "Saving tree" << endl;



  Long64_t nPassEvt=0;

  for (Long64_t ev = 0; ev < nentries; ev++) {

    // print progress
    if (ev % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());
    
    data.GetEntry(ev);

    // look for two electrons with pt > 20 GeV, |eta| <2.5
    Int_t    nGoodEle         = 0;
    Int_t    nEle             = data.GetInt("nEle");
    Float_t* elePt            = data.GetPtrFloat("elePt");
    Float_t* eleSCEta         = data.GetPtrFloat("eleSCEta");

    for(Int_t i=0; i < nEle; i++)
      {
	if(elePt[i]>20.0 && fabs(eleSCEta[i]) < 1.4442)
	  nGoodEle++;
	else if(elePt[i]>20.0 && fabs(eleSCEta[i]) > 1.566 && fabs(eleSCEta[i])<2.5)
	  nGoodEle++;
      }

    // look for two muons with pt > 20 GeV, |eta| <2.4
    Int_t   nGoodMuo        = 0;
    Int_t    nMu             = data.GetInt("nMu");
    Float_t* muPt            = data.GetPtrFloat("muPt");
    Float_t* muEta           = data.GetPtrFloat("muEta");

    for(Int_t i=0; i < nMu; i++)
      {
	if(muPt[i]>20.0 && fabs(muEta[i])<2.4)
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

  }
 
  newtree->Print();
  newtree->AutoSave();
  delete newfile_data;


  cout << "Number of passed events = " << nPassEvt << endl;
}
