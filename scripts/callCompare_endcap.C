{

  std::string dataF = "eiko123596.root";
  std::string mcF   = "meridian_mc.root";
  std::string treeName = "NTuples/Analysis";
  std::string cut   = "nVtx>0 && PHOLEAD_isEBEEGap==0 && PHOLEAD_isEBGap==0 && PHOLEAD_isEEGap==0 && PHOLEAD_isTransGap==0 && abs(PHOLEAD_eta)>1.56";
  gROOT->ProcessLine(".L compareTwoTrees.C");

  int decCode = 2;

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_energy",decCode,"energy",45,0.,45.,"E(#gamma) [GeV]","",cut);
  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_pt",decCode,"pt",20,0.,20.,"p_{T}(#gamma) [GeV/c]","",cut);
  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_eta",decCode,"eta",24,-3.,3.,"#eta(#gamma) ","",cut);
  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_phi",decCode,"phi",24,-TMath::Pi(),TMath::Pi(),"#phi(#gamma) ","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_r9",decCode,"r9",24,0,1.2,"R9(#gamma) ","",cut);
  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_preshowerEnergy",decCode,"esenergy",25,0,0.04,"Preshower E(#gamma) [GeV]","",cut);
  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_clustersSize",decCode,"nbcinsc",7,-0.5,6.5,"Number of basic clusters in SC","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_phiWidth",decCode,"scphiwidth",25,0,0.16,"SC #phi width","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_etaWidth",decCode,"scetawidth",25,0,0.09,"SC #eta width","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_eMax",decCode,"sceMax",30,0,30.0,"Energy of most energetic crystal [GeV]","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_e2nd",decCode,"sce2nd",10,0,10.0,"Energy of second energetic crystal [GeV]","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_e2x2",decCode,"sce2x2",40,0,40.0,"e2x2 [GeV]","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_e3x3",decCode,"sce3x3",40,0,40.0,"e3x3 [GeV]","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_e4x4",decCode,"sce4x4",40,0,40.0,"e4x4 [GeV]","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_e5x5",decCode,"sce5x5",40,0,40.0,"e5x5 [GeV]","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_covPhiPhi",decCode,"covphiphi",20,0,0.012,"covPhiPhi","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_covEtaPhi",decCode,"covetaphi",20, -0.006,0.006,"covEtaPhi","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_covEtaEta",decCode,"covetaeta",20,0,0.014,"covEtaEta","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_sigmaEtaEta",decCode,"sigetaeta",20,0,0.1,"#sigma_{#eta#eta}","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_sigmaIetaIeta",decCode,"sigietaieta",20,0,0.1,"#sigma_{i#etai#eta}","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_hadronicOverEm",decCode,"hadem",50,0,0.5,"E_{had}/E_{EM}","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_caloIso",decCode,"patcaloiso",20,-1.0,6.0,"PAT calIso [GeV]","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_ecalIso",decCode,"patecaliso",20,-1.0,3.0,"PAT ecalIso [GeV]","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_hcalIso",decCode,"pathcaliso",20,-1.0,5.0,"PAT hcalIso [GeV]","",cut);


  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_trackIso",decCode,"pattrkiso",20,0,2.0,"PAT trackIso [GeV/c]","",cut);

 
  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_trkSumPtSolidConeDR04",decCode,"solidtrkiso04",20,0,2.0,"trkSumPtSolidConeDR04 [GeV/c]","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_trkSumPtHollowConeDR04",decCode,"hollowtrkiso04",20,0,2.0,"trkSumPtHollowConeDR04 [GeV/c]","",cut);


  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_nTrkSolidConeDR04",decCode,"solidntrk04",5,-0.5,4.5,"nTrkSolidConeDR04","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_nTrkHollowConeDR04",decCode,"hollowntrk04",5,-0.5,4.5,"nTrkHollowConeDR04","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_ecalRecHitSumEtConeDR04",decCode,"ecaliso04",20,-1.0,3.0,"ecalRecHitSumEtConeDR04 [GeV]","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_hcalTowerSumEtConeDR04",decCode,"hcaliso04",20,-1.0,5.0,"hcalTowerSumEtConeDR04 [GeV]","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_hasConversionTracks",decCode,"hasconv",2,-0.5,1.5,"hasConversionTracks","",cut);

  compareTwoTrees(dataF,mcF,treeName,"PHOLEAD_hasPixelSeed",decCode,"hasconv",2,-0.5,1.5,"hasPixelSeed","",cut);

}
