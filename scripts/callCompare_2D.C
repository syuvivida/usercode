#include <TH1.h>
#include <TFile.h>
#include <iostream>
#include <TString.h>
#include <TPaveText.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TChain.h>
#include <TROOT.h>
#include <TMath.h>

void callCompare_2D()
{

  std::string treeName = "NTuples/Analysis";

  TChain* dataF = new TChain(treeName.data());
  dataF->Add("jan23data_all.root");
  TChain* mcF = new TChain(treeName.data());
  mcF->Add("jan21mc_900GeV.root");

  gROOT->ProcessLine(".L /home/syu/testYenJie/CMSSW_3_3_6_patch3/src/CRAB/scripts/compare2D.C");

  int decCode = 0;

  compare2D(dataF,mcF,"PHOLEAD_scPhi", "PHOLEAD_scEta",decCode,"scetavsscphi",24,-TMath::Pi(),TMath::Pi(),50,-3.0,3.0,"SC #phi","SC #eta");

  compare2D(dataF,mcF,"PHOLEAD_phi", "PHOLEAD_eta",decCode,"etavsscphi",24,-TMath::Pi(),TMath::Pi(),50,-3.0,3.0,"#phi(#gamma)","#eta(#gamma)");

  compare2D(dataF,mcF,"PHOLEAD_pt","PHOLEAD_eta",decCode,"etavspt",80,0,20.0,120,-3.0,3.0,"p_{T}(#gamma)","#eta(#gamma)");


}
