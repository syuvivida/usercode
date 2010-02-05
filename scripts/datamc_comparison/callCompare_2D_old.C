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

void callCompare_2D_old(std::string dataName, std::string mcName)
{

  std::string treeName = "NTuples/Analysis";

  TChain* dataF = new TChain(treeName.data());

 int nfile = 0;
  TSystemDirectory *base = new TSystemDirectory("root","root");
  base->SetDirectory("$CMSSW_BASE/src/CRAB/allOldData");
  TList *listOfFiles = base->GetListOfFiles();
  TIter fileIt(listOfFiles);
  TFile *fileH = new TFile();
  while(fileH = (TFile*)fileIt()) {
    std::string fileN = fileH->GetName();
    std::string baseString = "minbias_data";
    if( fileH->IsFolder())  continue;
    if(fileN.find(baseString) == std::string::npos)continue;
    cout << fileN.data() << endl;
    nfile++;
    dataF->Add(fileN.data());
  }
  std::cout << "Opened " << nfile << " files" << std::endl;


  TChain* mcF = new TChain(treeName.data());
  nfile = 0;

  TList *listOfFiles2 = base->GetListOfFiles();
  TIter fileIt2(listOfFiles2);
  TFile *fileH2 = new TFile();
  while(fileH2 = (TFile*)fileIt2()) {
    std::string fileN = fileH2->GetName();
    std::string baseString = "minbias_mc";
    if( fileH2->IsFolder())  continue;
    if(fileN.find(baseString) == std::string::npos)continue;
    cout << fileN.data() << endl;
    nfile++;
    mcF->Add(fileN.data());
  }

  std::cout << "Opened " << nfile << " files" << std::endl;


  gROOT->ProcessLine(".L $PWD/compare2D.C");

  int decCode = 0;

  std::string tempName = "_900GeV";
  if(dataName.find("2360") != std::string::npos)
    tempName = "_2360GeV";

  compare2D(dataF,mcF,"PHOLEAD_scPhi", "PHOLEAD_scEta",decCode,"scetavsscphi" + tempName,120,-TMath::Pi(),TMath::Pi(),115,-3.0,3.0,"SC #phi","SC #eta");

  compare2D(dataF,mcF,"PHOLEAD_phi", "PHOLEAD_eta",decCode,"etavsphi" + tempName,120,-TMath::Pi(),TMath::Pi(),115,-3.0,3.0,"#phi(#gamma)","#eta(#gamma)");

  compare2D(dataF,mcF,"PHOLEAD_pt","PHOLEAD_eta",decCode,"etavspt" + tempName,80,0,10.0,115,-3.0,3.0,"p_{T}(#gamma)","#eta(#gamma)");


}
