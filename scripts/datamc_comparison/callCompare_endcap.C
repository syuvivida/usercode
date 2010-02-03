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
  

  gROOT->ProcessLine(".L /home/syu/testYenJie/CMSSW_3_3_6_patch3/src/CRAB/scripts/compareTwoTrees.C");

  int decCode = 2;

  compareTwoTrees(dataF,mcF,"PHOLEAD_preshowerEnergy",decCode,"esenergy"+tempName,100,0,10.0,"Preshower E(#gamma) [GeV]","","",8.0);

  compareTwoTrees(dataF,mcF,"vtxX",decCode,"vx"+tempName,50,-0.1,1.5,"Highest #Sigma p_{T} vertex x position [cm]","","",-1000);

  compareTwoTrees(dataF,mcF,"vtxY",decCode,"vy"+tempName,50,-0.3,2.0,"Highest #Sigma p_{T} vertex y position [cm]","","",-1000);

  compareTwoTrees(dataF,mcF,"vtxZ",decCode,"vz"+tempName,50,-20.,30.,"Highest #Sigma p_{T} vertex z position [cm]","","",-1000);


  compareTwoTrees(dataF,mcF,"PHOLEAD_energy",decCode,"energy"+tempName,100,0.,100.,"E(#gamma) [GeV]","","",6.5);


  compareTwoTrees(dataF,mcF,"PHOLEAD_pt",decCode,"pt"+tempName,20,0.,20.,"p_{T}(#gamma) [GeV/c]","","");


  compareTwoTrees(dataF,mcF,"PHOLEAD_eta",decCode,"eta"+tempName,55,-3.,8.,"#eta(#gamma) ","","");


  compareTwoTrees(dataF,mcF,"PHOLEAD_phi",decCode,"phi"+tempName,60,-1.5*TMath::Pi(),3.5*TMath::Pi(),"#phi(#gamma) ","","");


  compareTwoTrees(dataF,mcF,"PHOLEAD_scEta",decCode,"sceta"+tempName,55,-3.,8.,"SC #eta","","");


  compareTwoTrees(dataF,mcF,"PHOLEAD_scPhi",decCode,"scphi"+tempName,60,-1.5*TMath::Pi(),3.5*TMath::Pi(),"SC #phi","","");

  TCut etaCut = "PHOLEAD_scEta > 0";
  compareTwoTrees(dataF,mcF,"PHOLEAD_scPhi",decCode,"scpphi"+tempName,60,-1.5*TMath::Pi(),3.5*TMath::Pi(),"SC #phi (#eta >0)","",etaCut);

  etaCut = "PHOLEAD_scEta < 0";
  compareTwoTrees(dataF,mcF,"PHOLEAD_scPhi",decCode,"scnphi"+tempName,60,-1.5*TMath::Pi(),3.5*TMath::Pi(),"SC #phi (#eta <0)","",etaCut);


  compareTwoTrees(dataF,mcF,"PHOLEAD_rawEnergy",decCode,"rawenergy"+tempName,100,0.,100.,"SC raw energy [GeV]","","",6.5);

  
  compareTwoTrees(dataF,mcF,"PHOLEAD_r9",decCode,"r9"+tempName,40,0,2.0,"R9(#gamma) ","","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_clustersSize",decCode,"nbcinsc"+tempName,11,-0.5,10.5,"Number of basic clusters in SC","","");


  compareTwoTrees(dataF,mcF,"PHOLEAD_phiWidth",decCode,"scphiwidth"+tempName,50,0,0.25,"SC #phi width","","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_etaWidth",decCode,"scetawidth"+tempName,50,0,0.25,"SC #eta width","","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_eMax",decCode,"sceMax"+tempName,60,0,60.0,"Energy of most energetic crystal [GeV]","","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_e2nd",decCode,"sce2nd"+tempName,35,0,35.0,"Energy of second energetic crystal [GeV]","","",6.0);


  compareTwoTrees(dataF,mcF,"PHOLEAD_eMax/PHOLEAD_e3x3", decCode,"sceSpikeRatio",200,0.0,2.0,"eMax/e3x3","","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_e2x2",decCode,"sce2x2"+tempName,100,0,100.0,"e2x2 [GeV]","","",6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_e3x3",decCode,"sce3x3"+tempName,100,0,100.0,"e3x3 [GeV]","","",6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_e4x4",decCode,"sce4x4"+tempName,100,0,100.0,"e4x4 [GeV]","","",6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_e5x5",decCode,"sce5x5"+tempName,100,0,100.0,"e5x5 [GeV]","","",6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_covPhiPhi",decCode,"covphiphi"+tempName,30,0,0.03,"covPhiPhi","","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_covEtaPhi",decCode,"covetaphi"+tempName,60, -0.02,0.04,"covEtaPhi","","");

  compareTwoTrees(dataF,mcF,"PHOLEAD_covEtaEta",decCode,"covetaeta"+tempName,30,0,0.03,"covEtaEta","","");

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


  /*
  // asking only track isolation  > 0
  TCut newcut = "PHOLEAD_trackIso > 0.0";
  compareTwoTrees(dataF,mcF,"PHOLEAD_trackIso",decCode,"pattrkisoNo0"+tempName,20,0,10.,"PAT trackIso [GeV/c]","",newcut,30.0); 
  newcut = "PHOLEAD_trkSumPtSolidConeDR04 > 0.0";
  compareTwoTrees(dataF,mcF,"PHOLEAD_trkSumPtSolidConeDR04",decCode,"solidtrkiso04No0"+tempName,20,0,10.,"trkSumPtSolidConeDR04 [GeV/c]","",newcut,30.0);
  compareTwoTrees(dataF,mcF,"PHOLEAD_nTrkSolidConeDR04",decCode,"solidntrk04No0"+tempName,11,-0.5,10.5,"nTrkSolidConeDR04","",newcut,30.0);

  newcut = "PHOLEAD_trkSumPtHollowConeDR04 > 0.0";
  compareTwoTrees(dataF,mcF,"PHOLEAD_trkSumPtHollowConeDR04",decCode,"hollowtrkiso04No0"+tempName,20,0,10.,"trkSumPtHollowConeDR04 [GeV/c]","",newcut,30.0);
  compareTwoTrees(dataF,mcF,"PHOLEAD_nTrkHollowConeDR04",decCode,"hollowntrk04No0"+tempName,11,-0.5,10.5,"nTrkHollowConeDR04","",newcut,30.0);
  */


  compareTwoTrees(dataF,mcF,"PHOLEAD_hcalIso",decCode,"pathcaliso"+tempName,55,-1.0,10.0,"PAT hcalIso [GeV]","","",6.0);

  compareTwoTrees(dataF,mcF,"PHOLEAD_hcalTowerSumEtConeDR04",decCode,"hcaliso04"+tempName,55,-1.0,10.0,"hcalTowerSumEtConeDR04 [GeV]","","",6.0);

  /*
  TCut noZCut = "PHOLEAD_hcalIso>0";

  compareTwoTrees(dataF,mcF,"PHOLEAD_hcalIso",decCode,"pathcalisoNo0"+tempName,55,-1.0,10.0,"PAT hcalIso [GeV]","",noZCut,6.0);

  noZCut = "PHOLEAD_hcalTowerSumEtConeDR04 > 0";

  compareTwoTrees(dataF,mcF,"PHOLEAD_hcalTowerSumEtConeDR04",decCode,"hcaliso04No0"+tempName,55,-1.0,10.0,"hcalTowerSumEtConeDR04 [GeV]","",noZCut,6.0);
  */

}

