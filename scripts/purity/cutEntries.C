#include <TFile.h>
#include <TH1.h>
#include <TChain.h>
#include <TBranch.h>
#include <iostream>
#include <fstream>
#include <string>
#include <TSystemDirectory.h>
#include <TList.h>
#include "tightSelection.h"

void cutEntries(std::string filename,
		std::string treename="skimLaurent/photonTree")
{

  // chain in all the ncu ntuples in the same directory
  TChain* pho = new TChain(treename.data());

  pho->Add(filename.data());
  ofstream fout;
  fout.open("event.dat", ios::out | ios::app);
  int ntotal = pho->Draw("run",basicCut);
  cout << "For " << filename.data() << endl; 
  cout << "All photons = " << ntotal << endl;
  fout << ntotal << " "; 
  
  int npass = 0; 
//   if(filename.find("qcd") == std::string::npos){
    npass= pho->Draw("run",basicCut && realCut); 
//   }
  cout << "Real photons = " << npass << endl;
  fout << npass << endl;
  fout.close();

}


