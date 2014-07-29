#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TChain.h>

void count(char *infname, char* treeName="tree/tree")
{
  TChain* t = new TChain(treeName);
  t->Add(infname);
  cout <<infname<<" # "<<t->GetEntries()<<endl;
   
}
