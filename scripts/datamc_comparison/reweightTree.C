#define reweightTree_cxx
#include "reweightTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

const int NBIN=42;
float nhitRatio[NBIN]={
//   0.462101,
//   0.684404,
//   0.765442,
//   0.888395,
//   1.53893,
//   3.27762,
//   11.9284,
//   58.1366,
//   65.1835
  0,
  0,
  0.375337,
  0.464066,
  0.472002,
  0.50328,
  0.796283,
  0.552819,
  0.736838,
  0.64276,
  0.645026,
  0.740528,
  0.892373,
  0.747445,
  0.750035,
  0.678839,
  0.671351,
  0.862722,
  0.713806,
  0.908073,
  0.871133,
  0.879219,
  0.932858,
  1.13504,
  1.02422,
  1.48307,
  1.5291,
  1.90471,
  1.90781,
  1.9915,
  2.98713,
  2.34516,
  2.59062,
  4.2468,
  4.34556,
  5.31387,
  8.14793,
  13.3939,
  11.031,
  11.4071,
  20.2293,
  21.1846
};


void reweightTree::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   TFile* newfile = new TFile("small.root","recreate");
   TTree* newtree = fChain->CloneTree(0);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      bool pass = false;
      int Bin = nHits/10;
      if(Bin < 0 || Bin > NBIN-1)pass=true;
      else
	{
	  float ref = nhitRatio[Bin]/nhitRatio[NBIN-1];
	  if( ref > gRandom->Rndm(jentry))pass=true;
	  else pass=false;
	}
      if(pass)newtree->Fill();

   }

   newtree->Print();
   newtree->AutoSave();
}
