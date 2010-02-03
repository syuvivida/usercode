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


void callOnlyData_barrel(std::string dataName)
{

  std::string treeName = "NTuples/Analysis";

  TChain* dataF = new TChain(treeName.data());
  dataF->Add(dataName.data());

 
  std::string tempName = "_900GeV";
  if(dataName.find("2360") != std::string::npos)
    tempName = "_2360GeV";  

  gROOT->ProcessLine(".L /home/syu/testYenJie/CMSSW_3_3_6_patch3/src/CRAB/scripts/overlayData.C");


  int decCode = 1;


  decCode = 0;

  overlayData(dataF,"nHits",decCode,"npixelhits"+tempName,30,-0.5,599.5,"Number of pixel hits in each event ","","");

  overlayData(dataF,"PHOLEAD_eta",decCode,"eta"+tempName,55,-3.,8.,"#eta(#gamma) ","","");

  overlayData(dataF,"PHOLEAD_scEta",decCode,"sceta"+tempName,55,-3.,8.,"SC #eta","","");
  
  decCode = 1;

  overlayData(dataF,"vtxX",decCode,"vx"+tempName,50,-0.1,1.5,"Highest #Sigma p_{T} vertex x position [cm]","","");

  overlayData(dataF,"vtxY",decCode,"vy"+tempName,50,-0.3,2.0,"Highest #Sigma p_{T} vertex y position [cm]","","");

  overlayData(dataF,"vtxZ",decCode,"vz"+tempName,50,-20.,30.,"Highest #Sigma p_{T} vertex z position [cm]","","");
  
  overlayData(dataF,"PHOLEAD_energy",decCode,"energy"+tempName,25,0.,25.,"E(#gamma) [GeV]","","");

  overlayData(dataF,"PHOLEAD_pt",decCode,"pt"+tempName,20,0.,20.,"p_{T}(#gamma) [GeV/c]","","");


  overlayData(dataF,"PHOLEAD_eta",decCode,"eta"+tempName,40,-3.,5.,"#eta(#gamma) ","","");

  overlayData(dataF,"PHOLEAD_scEta",decCode,"sceta"+tempName,40,-3.,5.,"SC #eta","","");


  overlayData(dataF,"PHOLEAD_phi",decCode,"phi"+tempName,60,-1.5*TMath::Pi(),3.5*TMath::Pi(),"#phi(#gamma) ","","");



  overlayData(dataF,"PHOLEAD_scPhi",decCode,"scphi"+tempName,60,-1.5*TMath::Pi(),3.5*TMath::Pi(),"SC #phi","","");


  overlayData(dataF,"PHOLEAD_rawEnergy",decCode,"rawenergy"+tempName,25,0.,25.,"SC raw energy [GeV]","","");

  
  overlayData(dataF,"PHOLEAD_r9",decCode,"r9"+tempName,40,0,2.0,"R9(#gamma) ","","");

  overlayData(dataF,"PHOLEAD_clustersSize",decCode,"nbcinsc"+tempName,11,-0.5,10.5,"Number of basic clusters in SC","","");


  overlayData(dataF,"PHOLEAD_phiWidth",decCode,"scphiwidth"+tempName,32,0,0.16,"SC #phi width","","");

  overlayData(dataF,"PHOLEAD_etaWidth",decCode,"scetawidth"+tempName,25,0,0.04,"SC #eta width","","");

  overlayData(dataF,"PHOLEAD_eMax",decCode,"sceMax"+tempName,25,0,25.0,"Energy of most energetic crystal [GeV]","","");

  overlayData(dataF,"PHOLEAD_e2nd",decCode,"sce2nd"+tempName,10,0,5.0,"Energy of second energetic crystal [GeV]","","");


  overlayData(dataF,"PHOLEAD_eMax/PHOLEAD_e3x3", decCode,"sceSpikeRatio",200,0.0,2.0,"eMax/e3x3","","");

  overlayData(dataF,"PHOLEAD_e2x2",decCode,"sce2x2"+tempName,30,0,30.0,"e2x2 [GeV]","","");

  overlayData(dataF,"PHOLEAD_e3x3",decCode,"sce3x3"+tempName,30,0,30.0,"e3x3 [GeV]","","");

  overlayData(dataF,"PHOLEAD_e4x4",decCode,"sce4x4"+tempName,30,0,30.0,"e4x4 [GeV]","","");

  overlayData(dataF,"PHOLEAD_e5x5",decCode,"sce5x5"+tempName,30,0,30.0,"e5x5 [GeV]","","");

  overlayData(dataF,"PHOLEAD_covPhiPhi",decCode,"covphiphi"+tempName,20,0,0.001,"covPhiPhi","","");

  overlayData(dataF,"PHOLEAD_covEtaPhi",decCode,"covetaphi"+tempName,30, -0.001,0.002,"covEtaPhi","","");

  overlayData(dataF,"PHOLEAD_covEtaEta",decCode,"covetaeta"+tempName,20,0,0.001,"covEtaEta","","");

  overlayData(dataF,"PHOLEAD_sigmaEtaEta",decCode,"sigetaeta"+tempName,20,0,0.05,"#sigma_{#eta#eta}","","");
  overlayData(dataF,"PHOLEAD_sigmaIetaIeta",decCode,"sigietaieta"+tempName,20,0,0.05,"#sigma_{i#etai#eta}","","");

  overlayData(dataF,"PHOLEAD_hadronicOverEm",decCode,"hadem"+tempName,60,0,0.6,"E_{had}/E_{EM}","","");


  overlayData(dataF,"PHOLEAD_caloIso",decCode,"patcaloiso"+tempName,55,-1.0,10.0,"PAT calIso [GeV]","","");

  overlayData(dataF,"PHOLEAD_ecalIso",decCode,"patecaliso"+tempName,55,-1.0,10.0,"PAT ecalIso [GeV]","","");


  overlayData(dataF,"PHOLEAD_ecalRecHitSumEtConeDR04",decCode,"ecaliso04"+tempName,55,-1.0,10.0,"ecalRecHitSumEtConeDR04 [GeV]","","");


  overlayData(dataF,"PHOLEAD_hasConversionTracks",decCode,"hasconv"+tempName,2,-0.5,1.5,"hasConversionTracks","","");

  overlayData(dataF,"PHOLEAD_hasPixelSeed",decCode,"haspixel"+tempName,2,-0.5,1.5,"hasPixelSeed","","");


  overlayData(dataF,"PHOLEAD_trackIso",decCode,"pattrkiso"+tempName,60,0.0,30.0,"PAT trackIso [GeV/c]","","");
 
  overlayData(dataF,"PHOLEAD_trkSumPtSolidConeDR04",decCode,"solidtrkiso04"+tempName,60,0,30.0,"trkSumPtSolidConeDR04 [GeV/c]","","");

  overlayData(dataF,"PHOLEAD_trkSumPtHollowConeDR04",decCode,"hollowtrkiso04"+tempName,60,0,30.0,"trkSumPtHollowConeDR04 [GeV/c]","","");
  

  overlayData(dataF,"PHOLEAD_trackIso",decCode,"close_pattrkiso"+tempName,20,0.0,10.0,"PAT trackIso [GeV/c]","","");
 
  overlayData(dataF,"PHOLEAD_trkSumPtSolidConeDR04",decCode,"close_solidtrkiso04"+tempName,20,0,10.0,"trkSumPtSolidConeDR04 [GeV/c]","","");

  overlayData(dataF,"PHOLEAD_trkSumPtHollowConeDR04",decCode,"close_hollowtrkiso04"+tempName,20,0,10.0,"trkSumPtHollowConeDR04 [GeV/c]","","");
  

  overlayData(dataF,"PHOLEAD_nTrkSolidConeDR04",decCode,"solidntrk04"+tempName,11,-0.5,10.5,"nTrkSolidConeDR04","","");

  overlayData(dataF,"PHOLEAD_nTrkHollowConeDR04",decCode,"hollowntrk04"+tempName,11,-0.5,10.5,"nTrkHollowConeDR04","","");

  overlayData(dataF,"PHOLEAD_hcalIso",decCode,"pathcaliso"+tempName,55,-1.0,10.0,"PAT hcalIso [GeV]","","");

  overlayData(dataF,"PHOLEAD_hcalTowerSumEtConeDR04",decCode,"hcaliso04"+tempName,55,-1.0,10.0,"hcalTowerSumEtConeDR04 [GeV]","","");

}
