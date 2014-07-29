#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

void countHisto(char *infname)
{
   TFile *inf = new TFile(infname);
   TH1F *h = (TH1F*)inf->FindObjectAny("hEvents");
//   TH1F *h = (TH1F*)inf->Get("compleSuperCluster/NoE");
//   cout <<infname<<" # "<<h->GetEntries()<<endl;
   cout <<infname<<" # "<<h->GetBinContent(1)<<endl;
   h->GetBinContent(1);
}
