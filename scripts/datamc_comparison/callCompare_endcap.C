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
#include "compareTwoTrees.C"

void callCompare_endcap(std::string dataName, std::string mcName, bool applyWeight=false)
{

  std::string treeName = "NTuples/Analysis";

  TChain* dataF = new TChain(treeName.data());
  dataF->Add(dataName.data());

  std::string treeName2 = applyWeight? "Analysis": treeName;
  TChain* mcF = new TChain(treeName2.data());
  mcF->Add(mcName.data());

  std::string tempName = "_900GeV";
  if(dataName.find("2360") != std::string::npos)
    tempName = "_2360GeV";  
  

  //  gROOT->ProcessLine(".L $PWD/compareTwoTrees.C");

  int decCode = 2;

  compareTwoTrees(dataF,mcF,"PHOLEAD_preshowerEnergy",decCode,"esenergy"+tempName,100,0,10.0,"Preshower E(#gamma) [GeV]","","",8.0);

  compareTwoTrees(dataF,mcF,"vtxX",decCode,"vx"+tempName,50,-0.1,1.5,"Highest #Sigma p_{T} vertex x position [cm]","","",-1000);

  compareTwoTrees(dataF,mcF,"vtxY",decCode,"vy"+tempName,50,-0.3,2.0,"Highest #Sigma p_{T} vertex y position [cm]","","",-1000);

  compareTwoTrees(dataF,mcF,"vtxZ",decCode,"vz"+tempName,50,-20.,30.,"Highest #Sigma p_{T} vertex z position [cm]","","",-1000);

  compareTwoTrees(dataF,mcF,"PHOLEAD_ncry",decCode,"scncry"+tempName,25,-0.5,49.5,"Number of hits used in SC","","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_energy",decCode,"energy"+tempName,100,0.,100.,"E(#gamma) [GeV]","","",6.5);


  compareTwoTrees(dataF,mcF,"PHOLEAD_pt",decCode,"pt"+tempName,20,0.,20.,"p_{T}(#gamma) [GeV/c]","","");

  // fake seed pt                                                                             
  compareTwoTrees(dataF,mcF,"PHOLEAD_eMax*sin( (2*atan(exp(-PHOLEAD_scEta))))",decCode,"pseudoseedpt"+tempName,20
		  ,0.,20.,"eMax*sin(SC #theta) [GeV]","","",10.0);


  compareTwoTrees(dataF,mcF,"PHOLEAD_eta",decCode,"eta"+tempName,55,-3.,8.,"#eta(#gamma) ","","");


  compareTwoTrees(dataF,mcF,"PHOLEAD_phi",decCode,"phi"+tempName,60,-1.5*TMath::Pi(),3.5*TMath::Pi(),"#phi(#gamma) ","","");


  compareTwoTrees(dataF,mcF,"PHOLEAD_scEta",decCode,"sceta"+tempName,55,-3.,8.,"SC #eta","","");


  compareTwoTrees(dataF,mcF,"PHOLEAD_scPhi",decCode,"scphi"+tempName,60,-1.5*TMath::Pi(),3.5*TMath::Pi(),"SC #phi","","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_rawEnergy",decCode,"rawenergy"+tempName,100,0.,100.,"SC raw energy [GeV]","","",6.5);

  compareTwoTrees(dataF,mcF,"PHOLEAD_rawEnergy*sin( (2*atan(exp(-PHOLEAD_scEta))))",decCode,"rawet"+tempName,25,0.,25.,"SC raw E_{T} [GeV]","","",6.5);
  
  compareTwoTrees(dataF,mcF,"PHOLEAD_r9",decCode,"r9"+tempName,40,0,2.0,"R9(#gamma) ","","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_clustersSize",decCode,"nbcinsc"+tempName,11,-0.5,10.5,"Number of basic clusters in SC","","");


  compareTwoTrees(dataF,mcF,"PHOLEAD_phiWidth",decCode,"scphiwidth"+tempName,50,0,0.25,"SC #phi width","","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_etaWidth",decCode,"scetawidth"+tempName,50,0,0.2,"SC #eta width","","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_eMax",decCode,"sceMax"+tempName,60,0,60.0,"Energy of most energetic crystal [GeV]","","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_e2nd",decCode,"sce2nd"+tempName,35,0,35.0,"Energy of second energetic crystal [GeV]","","",6.0);


  compareTwoTrees(dataF,mcF,"PHOLEAD_eMax/PHOLEAD_e3x3", decCode,"sceSpikeRatio"+tempName,100,0.0,2.0,"eMax/e3x3","","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_e2x2",decCode,"sce2x2"+tempName,100,0,100.0,"e2x2 [GeV]","","",6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_e3x3",decCode,"sce3x3"+tempName,100,0,100.0,"e3x3 [GeV]","","",6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_e4x4",decCode,"sce4x4"+tempName,100,0,100.0,"e4x4 [GeV]","","",6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_e5x5",decCode,"sce5x5"+tempName,100,0,100.0,"e5x5 [GeV]","","",6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_covPhiPhi",decCode,"covphiphi"+tempName,20,0,0.01,"covPhiPhi","","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_covEtaPhi",decCode,"covetaphi"+tempName,60, -0.005,0.01,"covEtaPhi","","",6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_covEtaEta",decCode,"covetaeta"+tempName,20,0,0.01,"covEtaEta","","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_sigmaEtaEta",decCode,"sigetaeta"+tempName,25,0,0.25,"#sigma_{#eta#eta}","","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_sigmaIetaIeta",decCode,"sigietaieta"+tempName,25,0,0.25,"#sigma_{i#etai#eta}","","");


  compareTwoTrees(dataF,mcF,"PHOLEAD_hadronicOverEm",decCode,"hadem"+tempName,60,0,0.6,"E_{had}/E_{EM}","","");


  compareTwoTrees(dataF,mcF,"PHOLEAD_caloIso",decCode,"patcaloiso"+tempName,55,-1.0,10.0,"PAT calIso [GeV]","","",6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_ecalIso",decCode,"patecaliso"+tempName,45,-2.0,7.0,"PAT ecalIso [GeV]","","",6.0);


  compareTwoTrees(dataF,mcF,"PHOLEAD_ecalRecHitSumEtConeDR04",decCode,"ecaliso04"+tempName,45,-2.0,7.0,"ecalRecHitSumEtConeDR04 [GeV]","","",6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_hasConversionTracks",decCode,"hasconv"+tempName,2,-0.5,1.5,"hasConversionTracks","","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_hasPixelSeed",decCode,"haspixel"+tempName,2,-0.5,1.5,"hasPixelSeed","","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_trackIso",decCode,"pattrkiso"+tempName,60,0.0,30.0,"PAT trackIso [GeV/c]","","",30.0);
 

  compareTwoTrees(dataF,mcF,"PHOLEAD_trkSumPtSolidConeDR04",decCode,"solidtrkiso04"+tempName,60,0,30.0,"trkSumPtSolidConeDR04 [GeV/c]","","",30.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_trkSumPtHollowConeDR04",decCode,"hollowtrkiso04"+tempName,60,0,30.0,"trkSumPtHollowConeDR04 [GeV/c]","","",30.0);


  compareTwoTrees(dataF,mcF,"PHOLEAD_trackIso",decCode,"close_pattrkiso"+tempName,20,0.0,10.0,"PAT trackIso [GeV/c]","","",30.0);
 
  compareTwoTrees(dataF,mcF,"PHOLEAD_trkSumPtSolidConeDR04",decCode,"close_solidtrkiso04"+tempName,20,0,10.0,"trkSumPtSolidConeDR04 [GeV/c]","","",30.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_trkSumPtHollowConeDR04",decCode,"close_hollowtrkiso04"+tempName,20,0,10.0,"trkSumPtHollowConeDR04 [GeV/c]","","",30.0);


  compareTwoTrees(dataF,mcF,"PHOLEAD_nTrkSolidConeDR04",decCode,"solidntrk04"+tempName,11,-0.5,10.5,"nTrkSolidConeDR04","","",30.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_nTrkHollowConeDR04",decCode,"hollowntrk04"+tempName,11,-0.5,10.5,"nTrkHollowConeDR04","","",30.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_hcalIso",decCode,"pathcaliso"+tempName,55,-1.0,10.0,"PAT hcalIso [GeV]","","",6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_hcalTowerSumEtConeDR04",decCode,"hcaliso04"+tempName,55,-1.0,10.0,"hcalTowerSumEtConeDR04 [GeV]","","",6.0);

  // ========================================================================================
  // EE+
  //  
  // ========================================================================================

  decCode = 12;

  TCut etaCut = "PHOLEAD_scEta > 0";
  

  compareTwoTrees(dataF,mcF,"PHOLEAD_preshowerEnergy",decCode,"esenergy"+tempName,100,0,10.0,"Preshower E(#gamma) [GeV]","",etaCut,8.0);

  compareTwoTrees(dataF,mcF,"vtxX",decCode,"vx"+tempName,50,-0.1,1.5,"Highest #Sigma p_{T} vertex x position [cm]","",etaCut,-1000);

  compareTwoTrees(dataF,mcF,"vtxY",decCode,"vy"+tempName,50,-0.3,2.0,"Highest #Sigma p_{T} vertex y position [cm]","",etaCut,-1000);

  compareTwoTrees(dataF,mcF,"vtxZ",decCode,"vz"+tempName,50,-20.,30.,"Highest #Sigma p_{T} vertex z position [cm]","",etaCut,-1000);


  compareTwoTrees(dataF,mcF,"PHOLEAD_energy",decCode,"energy"+tempName,100,0.,100.,"E(#gamma) [GeV]","",etaCut,6.5);


  compareTwoTrees(dataF,mcF,"PHOLEAD_pt",decCode,"pt"+tempName,20,0.,20.,"p_{T}(#gamma) [GeV/c]","",etaCut);

  // fake seed pt                                                                             
  compareTwoTrees(dataF,mcF,"PHOLEAD_eMax*sin( (2*atan(exp(-PHOLEAD_scEta))))",decCode,"pseudoseedpt"+tempName,20
		  ,0.,20.,"eMax*sin(SC #theta) [GeV]","",etaCut,10.0);


  compareTwoTrees(dataF,mcF,"PHOLEAD_eta",decCode,"eta"+tempName,55,-3.,8.,"#eta(#gamma) ","",etaCut);


  compareTwoTrees(dataF,mcF,"PHOLEAD_phi",decCode,"phi"+tempName,60,-1.5*TMath::Pi(),3.5*TMath::Pi(),"#phi(#gamma) ","",etaCut);


  compareTwoTrees(dataF,mcF,"PHOLEAD_scEta",decCode,"sceta"+tempName,55,-3.,8.,"SC #eta","",etaCut);


  compareTwoTrees(dataF,mcF,"PHOLEAD_scPhi",decCode,"scphi"+tempName,60,-1.5*TMath::Pi(),3.5*TMath::Pi(),"SC #phi","",etaCut);

  compareTwoTrees(dataF,mcF,"PHOLEAD_rawEnergy",decCode,"rawenergy"+tempName,100,0.,100.,"SC raw energy [GeV]","",etaCut,6.5);

  compareTwoTrees(dataF,mcF,"PHOLEAD_rawEnergy*sin( (2*atan(exp(-PHOLEAD_scEta))))",decCode,"rawet"+tempName,25,0.,25.,"SC raw E_{T} [GeV]","",etaCut,6.5);
  
  compareTwoTrees(dataF,mcF,"PHOLEAD_r9",decCode,"r9"+tempName,40,0,2.0,"R9(#gamma) ","",etaCut);

  compareTwoTrees(dataF,mcF,"PHOLEAD_clustersSize",decCode,"nbcinsc"+tempName,11,-0.5,10.5,"Number of basic clusters in SC","",etaCut);


  compareTwoTrees(dataF,mcF,"PHOLEAD_phiWidth",decCode,"scphiwidth"+tempName,50,0,0.25,"SC #phi width","",etaCut);

  compareTwoTrees(dataF,mcF,"PHOLEAD_etaWidth",decCode,"scetawidth"+tempName,50,0,0.2,"SC #eta width","",etaCut);

  compareTwoTrees(dataF,mcF,"PHOLEAD_eMax",decCode,"sceMax"+tempName,60,0,60.0,"Energy of most energetic crystal [GeV]","",etaCut);

  compareTwoTrees(dataF,mcF,"PHOLEAD_e2nd",decCode,"sce2nd"+tempName,35,0,35.0,"Energy of second energetic crystal [GeV]","",etaCut,6.0);


  compareTwoTrees(dataF,mcF,"PHOLEAD_eMax/PHOLEAD_e3x3", decCode,"sceSpikeRatio"+tempName,100,0.0,2.0,"eMax/e3x3","",etaCut);

  compareTwoTrees(dataF,mcF,"PHOLEAD_e2x2",decCode,"sce2x2"+tempName,100,0,100.0,"e2x2 [GeV]","",etaCut,6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_e3x3",decCode,"sce3x3"+tempName,100,0,100.0,"e3x3 [GeV]","",etaCut,6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_e4x4",decCode,"sce4x4"+tempName,100,0,100.0,"e4x4 [GeV]","",etaCut,6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_e5x5",decCode,"sce5x5"+tempName,100,0,100.0,"e5x5 [GeV]","",etaCut,6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_covPhiPhi",decCode,"covphiphi"+tempName,20,0,0.01,"covPhiPhi","",etaCut);

  compareTwoTrees(dataF,mcF,"PHOLEAD_covEtaPhi",decCode,"covetaphi"+tempName,60, -0.005,0.01,"covEtaPhi","",etaCut,6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_covEtaEta",decCode,"covetaeta"+tempName,20,0,0.01,"covEtaEta","",etaCut);

  compareTwoTrees(dataF,mcF,"PHOLEAD_sigmaEtaEta",decCode,"sigetaeta"+tempName,25,0,0.25,"#sigma_{#eta#eta}","",etaCut);

  compareTwoTrees(dataF,mcF,"PHOLEAD_sigmaIetaIeta",decCode,"sigietaieta"+tempName,25,0,0.25,"#sigma_{i#etai#eta}","",etaCut);


  compareTwoTrees(dataF,mcF,"PHOLEAD_hadronicOverEm",decCode,"hadem"+tempName,60,0,0.6,"E_{had}/E_{EM}","",etaCut);


  compareTwoTrees(dataF,mcF,"PHOLEAD_caloIso",decCode,"patcaloiso"+tempName,55,-1.0,10.0,"PAT calIso [GeV]","",etaCut,6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_ecalIso",decCode,"patecaliso"+tempName,45,-2.0,7.0,"PAT ecalIso [GeV]","",etaCut,6.0);


  compareTwoTrees(dataF,mcF,"PHOLEAD_ecalRecHitSumEtConeDR04",decCode,"ecaliso04"+tempName,45,-2.0,7.0,"ecalRecHitSumEtConeDR04 [GeV]","",etaCut,6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_hasConversionTracks",decCode,"hasconv"+tempName,2,-0.5,1.5,"hasConversionTracks","",etaCut);

  compareTwoTrees(dataF,mcF,"PHOLEAD_hasPixelSeed",decCode,"haspixel"+tempName,2,-0.5,1.5,"hasPixelSeed","",etaCut);

  compareTwoTrees(dataF,mcF,"PHOLEAD_trackIso",decCode,"pattrkiso"+tempName,60,0.0,30.0,"PAT trackIso [GeV/c]","",etaCut,30.0);
 

  compareTwoTrees(dataF,mcF,"PHOLEAD_trkSumPtSolidConeDR04",decCode,"solidtrkiso04"+tempName,60,0,30.0,"trkSumPtSolidConeDR04 [GeV/c]","",etaCut,30.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_trkSumPtHollowConeDR04",decCode,"hollowtrkiso04"+tempName,60,0,30.0,"trkSumPtHollowConeDR04 [GeV/c]","",etaCut,30.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_trackIso",decCode,"close_pattrkiso"+tempName,20,0.0,10.0,"PAT trackIso [GeV/c]","",etaCut,30.0);
 
  compareTwoTrees(dataF,mcF,"PHOLEAD_trkSumPtSolidConeDR04",decCode,"close_solidtrkiso04"+tempName,20,0,10.0,"trkSumPtSolidConeDR04 [GeV/c]","",etaCut,30.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_trkSumPtHollowConeDR04",decCode,"close_hollowtrkiso04"+tempName,20,0,10.0,"trkSumPtHollowConeDR04 [GeV/c]","",etaCut,30.0);


  compareTwoTrees(dataF,mcF,"PHOLEAD_nTrkSolidConeDR04",decCode,"solidntrk04"+tempName,11,-0.5,10.5,"nTrkSolidConeDR04","",etaCut,30.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_nTrkHollowConeDR04",decCode,"hollowntrk04"+tempName,11,-0.5,10.5,"nTrkHollowConeDR04","",etaCut,30.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_hcalIso",decCode,"pathcaliso"+tempName,55,-1.0,10.0,"PAT hcalIso [GeV]","",etaCut,6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_hcalTowerSumEtConeDR04",decCode,"hcaliso04"+tempName,55,-1.0,10.0,"hcalTowerSumEtConeDR04 [GeV]","",etaCut,6.0);

  // ========================================================================================
  // EE-
  //  
  // ========================================================================================

  decCode = 22;

  etaCut = "PHOLEAD_scEta < 0";
  

  compareTwoTrees(dataF,mcF,"PHOLEAD_preshowerEnergy",decCode,"esenergy"+tempName,100,0,10.0,"Preshower E(#gamma) [GeV]","",etaCut,8.0);

  compareTwoTrees(dataF,mcF,"vtxX",decCode,"vx"+tempName,50,-0.1,1.5,"Highest #Sigma p_{T} vertex x position [cm]","",etaCut,-1000);

  compareTwoTrees(dataF,mcF,"vtxY",decCode,"vy"+tempName,50,-0.3,2.0,"Highest #Sigma p_{T} vertex y position [cm]","",etaCut,-1000);

  compareTwoTrees(dataF,mcF,"vtxZ",decCode,"vz"+tempName,50,-20.,30.,"Highest #Sigma p_{T} vertex z position [cm]","",etaCut,-1000);


  compareTwoTrees(dataF,mcF,"PHOLEAD_energy",decCode,"energy"+tempName,100,0.,100.,"E(#gamma) [GeV]","",etaCut,6.5);


  compareTwoTrees(dataF,mcF,"PHOLEAD_pt",decCode,"pt"+tempName,20,0.,20.,"p_{T}(#gamma) [GeV/c]","",etaCut);

  // fake seed pt                                                                             
  compareTwoTrees(dataF,mcF,"PHOLEAD_eMax*sin( (2*atan(exp(-PHOLEAD_scEta))))",decCode,"pseudoseedpt"+tempName,20
		  ,0.,20.,"eMax*sin(SC #theta) [GeV]","",etaCut,10.0);


  compareTwoTrees(dataF,mcF,"PHOLEAD_eta",decCode,"eta"+tempName,55,-3.,8.,"#eta(#gamma) ","",etaCut);


  compareTwoTrees(dataF,mcF,"PHOLEAD_phi",decCode,"phi"+tempName,60,-1.5*TMath::Pi(),3.5*TMath::Pi(),"#phi(#gamma) ","",etaCut);


  compareTwoTrees(dataF,mcF,"PHOLEAD_scEta",decCode,"sceta"+tempName,55,-3.,8.,"SC #eta","",etaCut);


  compareTwoTrees(dataF,mcF,"PHOLEAD_scPhi",decCode,"scphi"+tempName,60,-1.5*TMath::Pi(),3.5*TMath::Pi(),"SC #phi","",etaCut);

  compareTwoTrees(dataF,mcF,"PHOLEAD_rawEnergy",decCode,"rawenergy"+tempName,100,0.,100.,"SC raw energy [GeV]","",etaCut,6.5);

  compareTwoTrees(dataF,mcF,"PHOLEAD_rawEnergy*sin( (2*atan(exp(-PHOLEAD_scEta))))",decCode,"rawet"+tempName,25,0.,25.,"SC raw E_{T} [GeV]","",etaCut,6.5);
  
  compareTwoTrees(dataF,mcF,"PHOLEAD_r9",decCode,"r9"+tempName,40,0,2.0,"R9(#gamma) ","",etaCut);

  compareTwoTrees(dataF,mcF,"PHOLEAD_clustersSize",decCode,"nbcinsc"+tempName,11,-0.5,10.5,"Number of basic clusters in SC","",etaCut);


  compareTwoTrees(dataF,mcF,"PHOLEAD_phiWidth",decCode,"scphiwidth"+tempName,50,0,0.25,"SC #phi width","",etaCut);

  compareTwoTrees(dataF,mcF,"PHOLEAD_etaWidth",decCode,"scetawidth"+tempName,50,0,0.2,"SC #eta width","",etaCut);

  compareTwoTrees(dataF,mcF,"PHOLEAD_eMax",decCode,"sceMax"+tempName,60,0,60.0,"Energy of most energetic crystal [GeV]","",etaCut);

  compareTwoTrees(dataF,mcF,"PHOLEAD_e2nd",decCode,"sce2nd"+tempName,35,0,35.0,"Energy of second energetic crystal [GeV]","",etaCut,6.0);


  compareTwoTrees(dataF,mcF,"PHOLEAD_eMax/PHOLEAD_e3x3", decCode,"sceSpikeRatio"+tempName,100,0.0,2.0,"eMax/e3x3","",etaCut);

  compareTwoTrees(dataF,mcF,"PHOLEAD_e2x2",decCode,"sce2x2"+tempName,100,0,100.0,"e2x2 [GeV]","",etaCut,6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_e3x3",decCode,"sce3x3"+tempName,100,0,100.0,"e3x3 [GeV]","",etaCut,6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_e4x4",decCode,"sce4x4"+tempName,100,0,100.0,"e4x4 [GeV]","",etaCut,6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_e5x5",decCode,"sce5x5"+tempName,100,0,100.0,"e5x5 [GeV]","",etaCut,6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_covPhiPhi",decCode,"covphiphi"+tempName,20,0,0.01,"covPhiPhi","",etaCut);

  compareTwoTrees(dataF,mcF,"PHOLEAD_covEtaPhi",decCode,"covetaphi"+tempName,60, -0.005,0.01,"covEtaPhi","",etaCut,6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_covEtaEta",decCode,"covetaeta"+tempName,20,0,0.01,"covEtaEta","",etaCut);

  compareTwoTrees(dataF,mcF,"PHOLEAD_sigmaEtaEta",decCode,"sigetaeta"+tempName,25,0,0.25,"#sigma_{#eta#eta}","",etaCut);

  compareTwoTrees(dataF,mcF,"PHOLEAD_sigmaIetaIeta",decCode,"sigietaieta"+tempName,25,0,0.25,"#sigma_{i#etai#eta}","",etaCut);


  compareTwoTrees(dataF,mcF,"PHOLEAD_hadronicOverEm",decCode,"hadem"+tempName,60,0,0.6,"E_{had}/E_{EM}","",etaCut);


  compareTwoTrees(dataF,mcF,"PHOLEAD_caloIso",decCode,"patcaloiso"+tempName,55,-1.0,10.0,"PAT calIso [GeV]","",etaCut,6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_ecalIso",decCode,"patecaliso"+tempName,45,-2.0,7.0,"PAT ecalIso [GeV]","",etaCut,6.0);


  compareTwoTrees(dataF,mcF,"PHOLEAD_ecalRecHitSumEtConeDR04",decCode,"ecaliso04"+tempName,45,-2.0,7.0,"ecalRecHitSumEtConeDR04 [GeV]","",etaCut,6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_hasConversionTracks",decCode,"hasconv"+tempName,2,-0.5,1.5,"hasConversionTracks","",etaCut);

  compareTwoTrees(dataF,mcF,"PHOLEAD_hasPixelSeed",decCode,"haspixel"+tempName,2,-0.5,1.5,"hasPixelSeed","",etaCut);

  compareTwoTrees(dataF,mcF,"PHOLEAD_trackIso",decCode,"pattrkiso"+tempName,60,0.0,30.0,"PAT trackIso [GeV/c]","",etaCut,30.0);
 

  compareTwoTrees(dataF,mcF,"PHOLEAD_trkSumPtSolidConeDR04",decCode,"solidtrkiso04"+tempName,60,0,30.0,"trkSumPtSolidConeDR04 [GeV/c]","",etaCut,30.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_trkSumPtHollowConeDR04",decCode,"hollowtrkiso04"+tempName,60,0,30.0,"trkSumPtHollowConeDR04 [GeV/c]","",etaCut,30.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_trackIso",decCode,"close_pattrkiso"+tempName,20,0.0,10.0,"PAT trackIso [GeV/c]","",etaCut,30.0);
 
  compareTwoTrees(dataF,mcF,"PHOLEAD_trkSumPtSolidConeDR04",decCode,"close_solidtrkiso04"+tempName,20,0,10.0,"trkSumPtSolidConeDR04 [GeV/c]","",etaCut,30.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_trkSumPtHollowConeDR04",decCode,"close_hollowtrkiso04"+tempName,20,0,10.0,"trkSumPtHollowConeDR04 [GeV/c]","",etaCut,30.0);


  compareTwoTrees(dataF,mcF,"PHOLEAD_nTrkSolidConeDR04",decCode,"solidntrk04"+tempName,11,-0.5,10.5,"nTrkSolidConeDR04","",etaCut,30.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_nTrkHollowConeDR04",decCode,"hollowntrk04"+tempName,11,-0.5,10.5,"nTrkHollowConeDR04","",etaCut,30.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_hcalIso",decCode,"pathcaliso"+tempName,55,-1.0,10.0,"PAT hcalIso [GeV]","",etaCut,6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_hcalTowerSumEtConeDR04",decCode,"hcaliso04"+tempName,55,-1.0,10.0,"hcalTowerSumEtConeDR04 [GeV]","",etaCut,6.0);


}

