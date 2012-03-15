{
  TStopwatch* myWatch = new TStopwatch();
  myWatch->Start();
  gROOT->ProcessLine(".L combineMC_eff.C++");

  // jetpt/photon pt
  gROOT->ProcessLine(
"combineMC_eff(\"h_ptratio_EB\",\"p_{T}(jet)/p_{T}(#gamma) in barrel\",5)"
);

  gROOT->ProcessLine(
"combineMC_eff(\"h_ptratio_EE\",\"p_{T}(jet)/p_{T}(#gamma) in endcap\",5)"
);

  //nvertex

 
  gROOT->ProcessLine(
"combineMC_eff(\"h_nvtx_EB\",\"number of vertex in barrel\",1,0.5,20.5)"
);

  gROOT->ProcessLine(
"combineMC_eff(\"h_nvtx_EE\",\"number of vertex in endcap\",1,0.5,20.5)"
);

 
   
  // jet pt and eta
  gROOT->ProcessLine(
"combineMC_eff(\"h_jetpt_eff_EB\",\"p_{T}(jet)\",4,30,200)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_jetpt_eff_EE\",\"p_{T}(jet)\",4,30,200)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_jeteta_eff_EB\",\"#eta(jet)\",4,-2.2,2.2)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_jeteta_eff_EE\",\"#eta(jet)\",4,-2.2,2.2)");


  gROOT->ProcessLine(
"combineMC_eff(\"h_phopt_eff_EB\",\"p_{T}(#gamma) in barrel\",4,40,200)"
);

  gROOT->ProcessLine(
"combineMC_eff(\"h_phopt_eff_EE\",\"p_{T}(#gamma) in endcap\",4,40,200)"
);

  gROOT->ProcessLine(
"combineMC_eff(\"h_phoeta_eff_40_45\",\"#eta(#gamma)\",4,-2.5,2.5)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_phoeta_eff_45_50\",\"#eta(#gamma)\",4,-2.5,2.5)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_phoeta_eff_50_55\",\"#eta(#gamma)\",4,-2.5,2.5)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_phoeta_eff_55_60\",\"#eta(#gamma)\",4,-2.5,2.5)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_phoeta_eff_60_65\",\"#eta(#gamma)\",4,-2.5,2.5)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_phoeta_eff_65_70\",\"#eta(#gamma)\",4,-2.5,2.5)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_phoeta_eff_70_75\",\"#eta(#gamma)\",4,-2.5,2.5)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_phoeta_eff_75_85\",\"#eta(#gamma)\",4,-2.5,2.5)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_phoeta_eff_85_95\",\"#eta(#gamma)\",4,-2.5,2.5)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_phoeta_eff_95_110\",\"#eta(#gamma)\",4,-2.5,2.5)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_phoeta_eff_110_130\",\"#eta(#gamma)\",4,-2.5,2.5)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_phoeta_eff_130_160\",\"#eta(#gamma)\",4,-2.5,2.5)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_phoeta_eff_160_200\",\"#eta(#gamma)\",4,-2.5,2.5)");


  gROOT->ProcessLine(
"combineMC_eff(\"h_phoeta_eff_allpt\",\"#eta(#gamma)\",4,-2.2,2.2)");


  // pstar in barrel
  gROOT->ProcessLine(
		     "combineMC_eff(\"h_pstar_EB_allpt\",\"p^{*} in barrel\",4,40,300)");

  gROOT->ProcessLine(
		     "combineMC_eff(\"h_pstar_COM3D_EB_allpt\",\"p^{*} in barrel\",4,40,300)");

  gROOT->ProcessLine(
		     "combineMC_eff(\"h_pstar_COMZ_EB_allpt\",\"p^{*} in barrel\",4,40,300)");

  // pstar in endcap
  gROOT->ProcessLine(
		     "combineMC_eff(\"h_pstar_EE_allpt\",\"p^{*} in endcap\",4,40,300)");

  gROOT->ProcessLine(
		     "combineMC_eff(\"h_pstar_COM3D_EE_allpt\",\"p^{*} in endcap\",4,40,300)");

  gROOT->ProcessLine(
		     "combineMC_eff(\"h_pstar_COMZ_EE_allpt\",\"p^{*} in endcap\",4,40,300)");

  // yBoost in barrel

  gROOT->ProcessLine(
"combineMC_eff(\"h_yB_EB_allpt\",\"y^{boost} in barrel\",4,-2.4,2.4)");

  // yboost in endcap

  gROOT->ProcessLine(
"combineMC_eff(\"h_yB_EE_allpt\",\"y^{boost} in endcap\",4,-2.4,2.4)");


  // costheta in barrel
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EB_40_45\",\"cos#theta^{*} in barrel\",4,0.0,0.95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EB_45_50\",\"cos#theta^{*} in barrel\",4,0.0,0.95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EB_50_55\",\"cos#theta^{*} in barrel\",4,0.0,0.95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EB_55_60\",\"cos#theta^{*} in barrel\",4,0.0,0.95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EB_60_65\",\"cos#theta^{*} in barrel\",4,0.0,0.95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EB_65_70\",\"cos#theta^{*} in barrel\",4,0.0,0.95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EB_70_75\",\"cos#theta^{*} in barrel\",4,0.0,0.95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EB_75_85\",\"cos#theta^{*} in barrel\",4,0.0,0.95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EB_85_95\",\"cos#theta^{*} in barrel\",4,0.0,0.95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EB_95_110\",\"cos#theta^{*} in barrel\",4,0.0,0.95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EB_110_130\",\"cos#theta^{*} in barrel\",4,0.0,0.95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EB_130_160\",\"cos#theta^{*} in barrel\",4,0.0,95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EB_160_200\",\"cos#theta^{*} in barrel\",4,0.0,0.95)");
 
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_EB_allpt\",\"cos#theta^{*} in barrel\",4,0.0,0.95)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COM3D_EB_allpt\",\"cos#theta^{*} in barrel\",4,0.0,0.95)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EB_allpt\",\"cos#theta^{*} in barrel\",4,0.0,0.95)");

  // costheta in endcap
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EE_40_45\",\"cos#theta^{*} in endcap\",4,0.0,0.95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EE_45_50\",\"cos#theta^{*} in endcap\",4,0.0,0.95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EE_50_55\",\"cos#theta^{*} in endcap\",4,0.0,0.95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EE_55_60\",\"cos#theta^{*} in endcap\",4,0.0,0.95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EE_60_65\",\"cos#theta^{*} in endcap\",4,0.0,0.95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EE_65_70\",\"cos#theta^{*} in endcap\",4,0.0,0.95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EE_70_75\",\"cos#theta^{*} in endcap\",4,0.0,0.95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EE_75_85\",\"cos#theta^{*} in endcap\",4,0.0,0.95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EE_85_95\",\"cos#theta^{*} in endcap\",4,0.0,0.95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EE_95_110\",\"cos#theta^{*} in endcap\",4,0.0,0.95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EE_110_130\",\"cos#theta^{*} in endcap\",4,0.0,0.95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EE_130_160\",\"cos#theta^{*} in endcap\",4,0.0,0.95)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EE_160_200\",\"cos#theta^{*} in endcap\",4,0.0,0.95)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_EE_allpt\",\"cos#theta^{*} in endcap\",4,0.0,0.95)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COM3D_EE_allpt\",\"cos#theta^{*} in endcap\",4,0.0,0.95)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EE_allpt\",\"cos#theta^{*} in endcap\",4,0.0,0.95)");


   // deltaphi in barrel
  gROOT->ProcessLine(
"combineMC_eff(\"h_dphi_EB_allpt\",\"#Delta #phi in barrel\",4)");

  // deltaphi in endcap

  gROOT->ProcessLine(
"combineMC_eff(\"h_dphi_EE_allpt\",\"#Delta #phi in endcap\",4)");


 
  // y at the center of mass frame in barrel
  gROOT->ProcessLine(
"combineMC_eff(\"h_ystar_EB_allpt\",\"y^{*} in barrel\",4,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_ystar_COM3D_EB_allpt\",\"y^{*} in barrel\",4,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_ystar_COMZ_EB_allpt\",\"y^{*} in barrel\",4,-2.4,2.4)");

  // y at the center of mass frame in endcap

  gROOT->ProcessLine(
"combineMC_eff(\"h_ystar_EE_allpt\",\"y^{*} in endcap\",4,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_ystar_COM3D_EE_allpt\",\"y^{*} in endcap\",4,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_ystar_COMZ_EE_allpt\",\"y^{*} in endcap\",4,-2.4,2.4)");


  // new variables
  gROOT->ProcessLine(
"combineMC_eff(\"h_njet_EB\",\"Jet multiplicity in barrel\",1, 0.5,5.5)");
  gROOT->ProcessLine(
"combineMC_eff(\"h_njet_EE\",\"Jet multiplicity in endcap\",1, 0.5,5.5)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_sumJet_EB_allpt\",\"cos#theta^{*} in barrel (using all jets)\",4,0.0,0.95)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_sumJet_EE_allpt\",\"cos#theta^{*} in barrel (using all jets)\",4,0.0,0.95)");


  gROOT->ProcessLine(
"combineMC_eff(\"h_pstar_sumJet_EB_allpt\",\"p^{*} in barrel (using all jets)\",4,40,300)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_pstar_sumJet_EE_allpt\",\"p^{*} in endcap (using all jets)\",4,40,300)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_ystar_sumJet_EB_allpt\",\"y^{*} in barrel (using all jets)\",4,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_ystar_sumJet_EE_allpt\",\"y^{*} in endcap (using all jets)\",4,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_sumjetpt_eff_EB\",\"p_{T}(jet^{all}) with barrel photon\",4,30,200)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_sumjetpt_eff_EB\",\"p_{T}(jet^{all}) with barrel photon\",4,30,200)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_sumjeteta_eff_EE\",\"#eta(jet^{all}) with endcap photon\",4,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_sumjeteta_eff_EE\",\"#eta(jet^{all}) with endcap photon\",4,-2.4,2.4)");


  myWatch->Stop();
  cout << myWatch->RealTime() << " seconds has passed! " << endl; 

}
