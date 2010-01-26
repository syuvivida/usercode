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

void callCompare_endcap()
{

  std::string treeName = "NTuples/Analysis";

  TChain* dataF = new TChain(treeName.data());
  dataF->Add("jan23data_all.root");
  TChain* mcF = new TChain(treeName.data());
  mcF->Add("jan21mc_900GeV.root");

  gROOT->ProcessLine(".L /home/syu/testYenJie/CMSSW_3_3_6_patch3/src/CRAB/scripts/compareTwoTrees.C");

  int decCode = 2;

  compareTwoTrees(dataF,mcF,"PHOLEAD_preshowerEnergy",decCode,"esenergy_900GeV",100,0,10.0,"Preshower E(#gamma) [GeV]","");

  compareTwoTrees(dataF,mcF,"vtxX",decCode,"vx_900GeV",50,-0.1,1.5,"Highest #Sigma p_{T} vertex x position [cm]","");

  compareTwoTrees(dataF,mcF,"vtxY",decCode,"vy_900GeV",50,-0.3,2.0,"Highest #Sigma p_{T} vertex y position [cm]","");

  compareTwoTrees(dataF,mcF,"vtxZ",decCode,"vz_900GeV",50,-20.,30.,"Highest #Sigma p_{T} vertex z position [cm]","");


  compareTwoTrees(dataF,mcF,"PHOLEAD_energy",decCode,"energy_900GeV",100,0.,100.,"E(#gamma) [GeV]","");


  compareTwoTrees(dataF,mcF,"PHOLEAD_pt",decCode,"pt_900GeV",20,0.,20.,"p_{T}(#gamma) [GeV/c]","");


  compareTwoTrees(dataF,mcF,"PHOLEAD_eta",decCode,"eta_900GeV",55,-3.,8.,"#eta(#gamma) ","");


  compareTwoTrees(dataF,mcF,"PHOLEAD_phi",decCode,"phi_900GeV",60,-1.5*TMath::Pi(),3.5*TMath::Pi(),"#phi(#gamma) ","");

  
  compareTwoTrees(dataF,mcF,"PHOLEAD_r9",decCode,"r9_900GeV",40,0,2.0,"R9(#gamma) ","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_clustersSize",decCode,"nbcinsc_900GeV",11,-0.5,10.5,"Number of basic clusters in SC","");


  compareTwoTrees(dataF,mcF,"PHOLEAD_phiWidth",decCode,"scphiwidth_900GeV",50,0,0.25,"SC #phi width","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_etaWidth",decCode,"scetawidth_900GeV",50,0,0.25,"SC #eta width","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_eMax",decCode,"sceMax_900GeV",60,0,60.0,"Energy of most energetic crystal [GeV]","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_e2nd",decCode,"sce2nd_900GeV",35,0,35.0,"Energy of second energetic crystal [GeV]","");


  compareTwoTrees(dataF,mcF,"PHOLEAD_eMax/PHOLEAD_e3x3", decCode,"sceSpikeRatio",50,0.0,2.0,"eMax/e3x3","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_e2x2",decCode,"sce2x2_900GeV",100,0,100.0,"e2x2 [GeV]","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_e3x3",decCode,"sce3x3_900GeV",100,0,100.0,"e3x3 [GeV]","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_e4x4",decCode,"sce4x4_900GeV",100,0,100.0,"e4x4 [GeV]","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_e5x5",decCode,"sce5x5_900GeV",100,0,100.0,"e5x5 [GeV]","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_covPhiPhi",decCode,"covphiphi_900GeV",20,0,0.03,"covPhiPhi","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_covEtaPhi",decCode,"covetaphi_900GeV",60, -0.02,0.04,"covEtaPhi","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_covEtaEta",decCode,"covetaeta_900GeV",20,0,0.03,"covEtaEta","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_sigmaEtaEta",decCode,"sigetaeta_900GeV",25,0,0.25,"#sigma_{#eta#eta}","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_sigmaIetaIeta",decCode,"sigietaieta_900GeV",25,0,0.25,"#sigma_{i#etai#eta}","");


  compareTwoTrees(dataF,mcF,"PHOLEAD_hadronicOverEm",decCode,"hadem_900GeV",60,0,0.6,"E_{had}/E_{EM}","");


  compareTwoTrees(dataF,mcF,"PHOLEAD_caloIso",decCode,"patcaloiso_900GeV",55,-1.0,10.0,"PAT calIso [GeV]","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_ecalIso",decCode,"patecaliso_900GeV",45,-2.0,7.0,"PAT ecalIso [GeV]","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_hcalIso",decCode,"pathcaliso_900GeV",55,-1.0,10.0,"PAT hcalIso [GeV]","");


  compareTwoTrees(dataF,mcF,"PHOLEAD_ecalRecHitSumEtConeDR04",decCode,"ecaliso04_900GeV",45,-2.0,7.0,"ecalRecHitSumEtConeDR04 [GeV]","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_hcalTowerSumEtConeDR04",decCode,"hcaliso04_900GeV",55,-1.0,10.0,"hcalTowerSumEtConeDR04 [GeV]","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_hasConversionTracks",decCode,"hasconv_900GeV",2,-0.5,1.5,"hasConversionTracks","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_hasPixelSeed",decCode,"haspixel_900GeV",2,-0.5,1.5,"hasPixelSeed","");


  compareTwoTrees(dataF,mcF,"PHOLEAD_trackIso",decCode,"pattrkiso_900GeV",60,0.0,30.0,"PAT trackIso [GeV/c]","");
 

  compareTwoTrees(dataF,mcF,"PHOLEAD_trkSumPtSolidConeDR04",decCode,"solidtrkiso04_900GeV",60,0,30.0,"trkSumPtSolidConeDR04 [GeV/c]","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_trkSumPtHollowConeDR04",decCode,"hollowtrkiso04_900GeV",60,0,30.0,"trkSumPtHollowConeDR04 [GeV/c]","");
  

  compareTwoTrees(dataF,mcF,"PHOLEAD_nTrkSolidConeDR04",decCode,"solidntrk04_900GeV",11,-0.5,10.5,"nTrkSolidConeDR04","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_nTrkHollowConeDR04",decCode,"hollowntrk04_900GeV",11,-0.5,10.5,"nTrkHollowConeDR04","");



  // asking only track isolation  > 0

  TCut newcut = "PHOLEAD_trackIso > 0.0";
  compareTwoTrees(dataF,mcF,"PHOLEAD_trackIso",decCode,"pattrkisoNo0_900GeV",60,0,30.,"PAT trackIso [GeV/c]","",newcut); 
  newcut = "PHOLEAD_trkSumPtSolidConeDR04 > 0.0";
  compareTwoTrees(dataF,mcF,"PHOLEAD_trkSumPtSolidConeDR04",decCode,"solidtrkiso04No0_900GeV",60,0,30.,"trkSumPtSolidConeDR04 [GeV/c]","",newcut);
  compareTwoTrees(dataF,mcF,"PHOLEAD_nTrkSolidConeDR04",decCode,"solidntrk04No0_900GeV",11,-0.5,10.5,"nTrkSolidConeDR04","",newcut);

  newcut = "PHOLEAD_trkSumPtHollowConeDR04 > 0.0";
  compareTwoTrees(dataF,mcF,"PHOLEAD_trkSumPtHollowConeDR04",decCode,"hollowtrkiso04No0_900GeV",60,0,30.,"trkSumPtHollowConeDR04 [GeV/c]","",newcut);
  compareTwoTrees(dataF,mcF,"PHOLEAD_nTrkHollowConeDR04",decCode,"hollowntrk04No0_900GeV",11,-0.5,10.5,"nTrkHollowConeDR04","",newcut);





}
  // 1463,15233
  // 1370,14976
