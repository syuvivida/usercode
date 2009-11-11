{
// void displayTwoGraphs(std::string file1, std::string file2, std::string histo1,
// 		      std::string xtitle="", std::string ytitle="",
// 		      std::string histo2="", std::string title="",std::string leg1="",
// 		      std::string leg2="", std::string output="test", float xmin=-1, float xmax=-1)


  gROOT->ProcessLine(".L ~/scripts/displayTwoGraphs.C");
  
  displayTwoGraphs("../triggerdata/PhotonJet_15/PhotonJet_15_bit_default_trigeff.root","../triggerdata/PhotonJet_15/PhotonJet_15_bit_R9_trigeff.root","heff01", "","","","HLT_Photon15_L1R","Real #gamma:No R9 cut","Real #gamma:R9 > 0.93","PhotonJet15_realgamma_R9");

  displayTwoGraphs("../triggerdata/PhotonJet_15/PhotonJet_15_bit_default_trigeff.root","../triggerdata/PhotonJet_15/PhotonJet_15_bit_default_trigeff.root","heff00", "","","heff06","HLT_Photon15_L1R","All #gamma: /HLT_L1EG5","All #gamma: /HLT_MU5","PhotonJet15_realgamma_mutrig",13,70);


  displayTwoGraphs("../triggerdata/PhotonJet_15/PhotonJet_15_bit_default_trigeff.root","../triggerdata/PhotonJet_15/PhotonJet_15_bitEG8_default_trigeff.root","heff01", "","","","HLT_Photon15_L1R","/HLT_L1EG5","/HLT_L1EG8","PhotonJet15_realgamma_EG58",0,50);


  displayTwoGraphs("../triggerdata/PhotonJet_15/PhotonJet_15_bit_default_trigeff.root","../triggerdata/PhotonJet_15/PhotonJet_15_bit_R9_trigeff.root","heff01", "","","","HLT_Photon15_L1R","Real #gamma:No R9 cut","Real #gamma:R9 > 0.93","PhotonJet15_realgamma_R9_close",0,50);


  displayTwoGraphs("../triggerdata/PhotonJet_15/PhotonJet_15_bit20Iso_default_trigeff.root","../triggerdata/PhotonJet_15/PhotonJet_15_bit20Iso_R9_trigeff.root","heff01", "","","","HLT_Photon20_Iso","Real #gamma:No R9 cut","Real #gamma:R9 > 0.93","PhotonJet15_20Iso_realgamma_R9_close",10,70);


  displayTwoGraphs("../triggerdata/PhotonJet_15/PhotonJet_15_bit25R_default_trigeff.root","../triggerdata/PhotonJet_15/PhotonJet_15_bit25R_R9_trigeff.root","heff01", "","","","HLT_Photon25_L1R","Real #gamma:No R9 cut","Real #gamma:R9 > 0.93","PhotonJet15_25R_realgamma_R9_close",10,70);



  displayTwoGraphs("../triggerdata/PhotonJet_15/PhotonJet_15_bit_default_trigeff.root","../triggerdata/QCDEMEnriched_20_30_reco/QCDEMEnriched_20_30_reco_bit_default_trigeff.root","heff01", "","","heff02","HLT_Photon15_L1R","Real #gamma","Jets","PhotonJet15_QCD_NoR9close",13,70);

  displayTwoGraphs("../triggerdata/PhotonJet_15/PhotonJet_15_bit_R9_trigeff.root","../triggerdata/QCDEMEnriched_20_30_reco/QCDEMEnriched_20_30_reco_bit_R9_trigeff.root","heff01", "","","heff02","HLT_Photon15_L1R (R9>0.93)","Real #gamma","Jets","PhotonJet15_QCD_R9close",13,70);


  displayTwoGraphs("../triggerdata/PhotonJet_15/PhotonJet_15_bit25R_default_trigeff.root","../triggerdata/QCDEMEnriched_20_30_reco/QCDEMEnriched_20_30_reco_bit25R_default_trigeff.root","heff01", "","","heff02","HLT_Photon25_L1R","Real #gamma","Jets","PhotonJet15_QCD_25RNoR9close",20,70);

  displayTwoGraphs("../triggerdata/PhotonJet_15/PhotonJet_15_bit25R_R9_trigeff.root","../triggerdata/QCDEMEnriched_20_30_reco/QCDEMEnriched_20_30_reco_bit25R_R9_trigeff.root","heff01", "","","heff02","HLT_Photon25_L1R(R9>0.93)","Real #gamma","Jets","PhotonJet15_QCD_25RR9close",20,70);


  displayTwoGraphs("../triggerdata/PhotonJet_15/PhotonJet_15_bit20Iso_default_trigeff.root","../triggerdata/QCDEMEnriched_20_30_reco/QCDEMEnriched_20_30_reco_bit20Iso_default_trigeff.root","heff01", "","","heff02","HLT_Photon20_Iso","Real #gamma","Jets","PhotonJet15_QCD_20IsoNoR9close",15,70);

  displayTwoGraphs("../triggerdata/PhotonJet_15/PhotonJet_15_bit20Iso_R9_trigeff.root","../triggerdata/QCDEMEnriched_20_30_reco/QCDEMEnriched_20_30_reco_bit20Iso_R9_trigeff.root","heff01", "","","heff00","HLT_Photon20_Iso(R9>0.93)","Real #gamma","Jets","PhotonJet15_QCD_20IsoR9close",15,70);


  displayTwoGraphs("../triggerdata/PhotonJet_15/PhotonJet_15_bit25R_default_trigeff.root","../triggerdata/singleGamma/singleGamma_match_default_trigeff.root","heff01", "","","heff07","HLT_Photon25_L1R","RECO","HLTDEBUG: single #gamma","PhotonJet15_SingleGamma_bit",10,70);


  displayTwoGraphs("../triggerdata/PhotonJet_15/PhotonJet_15_bit25R_default_trigeff.root","../triggerdata/singleElectron/singleElectron_match_default_trigeff.root","heff01", "","","heff07","HLT_Photon25_L1R","RECO","HLTDEBUG: single e","PhotonJet15_SingleElectron_bit",10,70);


  displayTwoGraphs("../triggerdata/PhotonJet_15/PhotonJet_15_bit20Iso_default_trigeff.root","../triggerdata/QCDFlat/QCDFlat_bit_default_trigeff.root","heff01", "","","heff08","HLT_Photon20_Iso","RECO","HLTDEBUG: single #gamma","PhotonJet15_QCDFlat_bit20Iso",15,70);



  displayTwoGraphs("../triggerdata/PhotonJet_15/PhotonJet_15_bit_default_trigeff.root","../triggerdata/QCDEMEnriched_20_30_reco/QCDEMEnriched_20_30_reco_bit_default_trigeff.root","heff01", "","","heff03","HLT_Photon15_L1R","Real #gamma","Quark fragmentation","PhotonJet15_QCD_QuarkNoR9close",13,70);


  displayTwoGraphs("../triggerdata/PhotonJet_15/PhotonJet_15_bit_default_trigeff.root","../triggerdata/QCDEMEnriched_20_30_reco/QCDEMEnriched_20_30_reco_bit_default_trigeff.root","heff01", "","","heff04","HLT_Photon15_L1R","Real #gamma","Gluon fragmentation","PhotonJet15_QCD_GluonNoR9close",13,70);

  displayTwoGraphs("../triggerdata/PhotonJet_15/PhotonJet_15_bit20Iso_default_trigeff.root","../triggerdata/QCDEMEnriched_20_30_reco/QCDEMEnriched_20_30_reco_bit20Iso_default_trigeff.root","heff01", "","","heff03","HLT_Photon20_Iso","Real #gamma","Quark fragmentation","PhotonJet15_QCD_20IsoQuarkNoR9close",15,70);


  displayTwoGraphs("../triggerdata/PhotonJet_15/PhotonJet_15_bit20Iso_default_trigeff.root","../triggerdata/QCDEMEnriched_20_30_reco/QCDEMEnriched_20_30_reco_bit20Iso_default_trigeff.root","heff01", "","","heff04","HLT_Photon20_Iso","Real #gamma","Gluon fragmentation","PhotonJet15_QCD_20IsoGluonNoR9close",15,70);



  displayTwoGraphs("../triggerdata/PhotonJet_15/PhotonJet_15_new_default_trigeff.root","../triggerdata/PhotonJet_15/PhotonJet_15_new_R9_trigeff.root","heff10", "","","heff10","HLT_L1SingleEG5","All #gamma: No R9","All #gamma: R9>0.93","PhotonJet15_L1EG5R9close",0,50);


  displayTwoGraphs("../triggerdata/QCDEMEnriched_20_30_reco/QCDEMEnriched_20_30_reco_new_default_trigeff.root","../triggerdata/QCDEMEnriched_20_30_reco/QCDEMEnriched_20_30_reco_new_R9_trigeff.root","heff10", "","","heff10","HLT_L1SingleEG5","All #gamma: No R9","All #gamma: R9>0.93","QCD_L1EG5R9close",0,50);


  displayTwoGraphs("../triggerdata/QCDEMEnriched_20_30_reco/QCDEMEnriched_20_30_reco_new_default_trigeff.root","../triggerdata/QCDEMEnriched_20_30_reco/QCDEMEnriched_20_30_reco_new_R9_trigeff.root","heff00", "","","heff06","HLT_Photon15_L1R","All #gamma: /HLT_L1EG5","All #gamma: /HLT_MU5","QCD_mutrig",13,70);



  displayTwoGraphs("../triggerdata/PhotonJet_15/PhotonJet_15_barrel_default_trigeff.root","../triggerdata/PhotonJet_15/PhotonJet_15_barrel_R9_trigeff.root","heff10", "","","heff10","HLT_L1SingleEG5","Barrel #gamma: No R9","Barrel #gamam: R9>0.93","PhotonJet15_BarrelL1EG8R9close",0,50);



}
