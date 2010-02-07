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
#include "compare2D.C"

void callCompare_2D(std::string dataName, std::string mcName, bool applyWeight=false)
{

  std::string treeName = "NTuples/Analysis";

  TChain* dataF = new TChain(treeName.data());
  dataF->Add(dataName.data());

  std::string treeName2 = applyWeight? "Analysis": treeName;
  TChain* mcF = new TChain(treeName2.data());
  mcF->Add(mcName.data());

  //  gROOT->ProcessLine(".L $PWD/compare2D.C");

  int decCode = 0;

  std::string tempName = "_900GeV";
  if(dataName.find("2360") != std::string::npos)
    tempName = "_2360GeV";

  compare2D(dataF,mcF,"PHOLEAD_scPhi", "PHOLEAD_scEta",decCode,"scetavsscphi" + tempName,120,-TMath::Pi(),TMath::Pi(),115,-3.0,3.0,"SC #phi","SC #eta");

  compare2D(dataF,mcF,"PHOLEAD_phi", "PHOLEAD_eta",decCode,"etavsphi" + tempName,120,-TMath::Pi(),TMath::Pi(),115,-3.0,3.0,"#phi(#gamma)","#eta(#gamma)");

  compare2D(dataF,mcF,"PHOLEAD_pt","PHOLEAD_eta",decCode,"etavspt" + tempName,80,0,10.0,115,-3.0,3.0,"p_{T}(#gamma)","#eta(#gamma)");


}
