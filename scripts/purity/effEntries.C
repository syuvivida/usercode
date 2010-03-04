#include <TFile.h>
#include <TH1.h>
#include <TChain.h>
#include <TBranch.h>
#include <iostream>
#include <fstream>
#include <string>
#include <TSystemDirectory.h>
#include <TList.h>
#include "Selection.h"

void effEntries(std::string filename,
		std::string treename="ncuAnalyzerKit/EventTree")
{

  // chain in all the ncu ntuples in the same directory
  TChain* pho = new TChain(treename.data());

  pho->Add(filename.data());
  ofstream fout;
  fout.open("event.dat", ios::out | ios::app);
  int ntotal = pho->Draw("run",basicCut);
  fout << ntotal << " "; 
  
  int npass = 0; 
  npass= pho->Draw("run",basicCut && IDCut); 
  fout << npass << " ";

  int ntotal_real = pho->Draw("run",basicCut && realCut);
  fout << ntotal_real << " "; 
  
  int npass_real = 0; 
  npass_real= pho->Draw("run",basicCut && IDCut && realCut); 
  fout << npass_real << endl;

  fout.close();

}


