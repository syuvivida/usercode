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
#include "overlayData.C"

void callOnlyData_endcap(std::string dataName)
{

  std::string treeName = "NTuples/Analysis";

  TChain* dataF = new TChain(treeName.data());
  dataF->Add(dataName.data());

  std::string tempName = "_900GeV";
  if(dataName.find("2360") != std::string::npos)
    tempName = "_2360GeV";  
  

  //  gROOT->ProcessLine(".L $PWD/overlayData.C");

  int decCode = 2;

  overlayData(dataF,"PHOLEAD_preshowerEnergy",decCode,"esenergy"+tempName,10,0,10.0,"Preshower E(#gamma) [GeV]","","",true);

  overlayData(dataF,"PHOLEAD_energy",decCode,"energy"+tempName, 25,0., 100,"E(#gamma) [GeV]","","",true);


  overlayData(dataF,"PHOLEAD_rawEnergy*sin( (2*atan(exp(-PHOLEAD_scEta))))",decCode,"rawet"+tempName,10,0.,20.,"raw E_{T}(#gamma) [GeV]","","",true);


  overlayData(dataF,"PHOLEAD_eta",decCode,"eta"+tempName,20,-3.,7.,"#eta(#gamma) ","","");


  overlayData(dataF,"PHOLEAD_phi",decCode,"phi"+tempName,10,-1.5*TMath::Pi(),3.5*TMath::Pi(),"#phi(#gamma) ","","");


  overlayData(dataF,"PHOLEAD_scEta",decCode,"sceta"+tempName,20,-3.,7.,"SC #eta","","");


  overlayData(dataF,"PHOLEAD_scPhi",decCode,"scphi"+tempName,10,-1.5*TMath::Pi(),3.5*TMath::Pi(),"SC #phi","","");

  TCut etaCut = "PHOLEAD_scEta > 0";
  overlayData(dataF,"PHOLEAD_scPhi",decCode,"scpphi"+tempName,10,-1.5*TMath::Pi(),3.5*TMath::Pi(),"SC #phi (#eta >0)","",etaCut);

  etaCut = "PHOLEAD_scEta < 0";
  overlayData(dataF,"PHOLEAD_scPhi",decCode,"scnphi"+tempName,10,-1.5*TMath::Pi(),3.5*TMath::Pi(),"SC #phi (#eta <0)","",etaCut);


  overlayData(dataF,"PHOLEAD_rawEnergy",decCode,"rawenergy"+tempName, 25,0., 100,"SC raw energy [GeV]","","",true);

  
  overlayData(dataF,"PHOLEAD_r9",decCode,"r9"+tempName,10,0,2.0,"R9(#gamma) ","","");

  overlayData(dataF,"PHOLEAD_clustersSize",decCode,"nbcinsc"+tempName,11,-0.5,10.5,"Number of basic clusters in SC","","");


  overlayData(dataF,"PHOLEAD_phiWidth",decCode,"scphiwidth"+tempName,10,0,0.25,"SC #phi width","","");

  overlayData(dataF,"PHOLEAD_etaWidth",decCode,"scetawidth"+tempName,20,0,0.25,"SC #eta width","","");

  overlayData(dataF,"PHOLEAD_eMax",decCode,"sceMax"+tempName,15,0,60.0,"Energy of most energetic crystal [GeV]","","",true);

  overlayData(dataF,"PHOLEAD_e2nd",decCode,"sce2nd"+tempName, 5,0,20.0,"Energy of second energetic crystal [GeV]","","",true);


  overlayData(dataF,"PHOLEAD_eMax/PHOLEAD_e3x3", decCode,"sceSpikeRatio", 20,0.0,2.0,"eMax/e3x3","","");

  overlayData(dataF,"PHOLEAD_e2x2",decCode,"sce2x2"+tempName, 25,0, 100.,"e2x2 [GeV]","","",true);

  overlayData(dataF,"PHOLEAD_e3x3",decCode,"sce3x3"+tempName, 25,0, 100.,"e3x3 [GeV]","","",true);

  overlayData(dataF,"PHOLEAD_e4x4",decCode,"sce4x4"+tempName, 25,0, 100.,"e4x4 [GeV]","","",true);

  overlayData(dataF,"PHOLEAD_e5x5",decCode,"sce5x5"+tempName, 25,0, 100.,"e5x5 [GeV]","","",true);

  overlayData(dataF,"PHOLEAD_covPhiPhi",decCode,"covphiphi"+tempName,10,0,0.02,"covPhiPhi","","");

  overlayData(dataF,"PHOLEAD_covEtaPhi",decCode,"covetaphi"+tempName, 5, -0.005,0.02,"covEtaPhi","","");

  overlayData(dataF,"PHOLEAD_covEtaEta",decCode,"covetaeta"+tempName, 10,0,0.02,"covEtaEta","","");

  overlayData(dataF,"PHOLEAD_sigmaEtaEta",decCode,"sigetaeta"+tempName, 10,0,0.25,"#sigma_{#eta#eta}","","");

  overlayData(dataF,"PHOLEAD_sigmaIetaIeta",decCode,"sigietaieta"+tempName, 10,0,0.25,"#sigma_{i#etai#eta}","","");


  overlayData(dataF,"PHOLEAD_hadronicOverEm",decCode,"hadem"+tempName,16,0,0.8,"E_{had}/E_{EM}","","",true);


  overlayData(dataF,"PHOLEAD_ecalRecHitSumEtConeDR04",decCode,"ecaliso04"+tempName,10,-2.0,6.0,"ecalRecHitSumEtConeDR04 [GeV]","","");

  overlayData(dataF,"PHOLEAD_hasConversionTracks",decCode,"hasconv"+tempName,2,-0.5,1.5,"hasConversionTracks","","");

  overlayData(dataF,"PHOLEAD_hasPixelSeed",decCode,"haspixel"+tempName,2,-0.5,1.5,"hasPixelSeed","","");

  overlayData(dataF,"PHOLEAD_trkSumPtSolidConeDR04",decCode,"close_solidtrkiso04"+tempName,12,0,24.0,"trkSumPtSolidConeDR04 [GeV/c]","","",true);

  overlayData(dataF,"PHOLEAD_trkSumPtHollowConeDR04",decCode,"close_hollowtrkiso04"+tempName,12,0,24.0,"trkSumPtHollowConeDR04 [GeV/c]","","",true);


  overlayData(dataF,"PHOLEAD_nTrkSolidConeDR04",decCode,"solidntrk04"+tempName,11,-0.5,10.5,"nTrkSolidConeDR04","","");

  overlayData(dataF,"PHOLEAD_nTrkHollowConeDR04",decCode,"hollowntrk04"+tempName,11,-0.5,10.5,"nTrkHollowConeDR04","","");


  overlayData(dataF,"PHOLEAD_hcalTowerSumEtConeDR04",decCode,"hcaliso04"+tempName, 5,-1.0, 9.0,"hcalTowerSumEtConeDR04 [GeV]","","",true);

}

