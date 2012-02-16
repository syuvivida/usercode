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
"combineMC_eff(\"h_jetpt_eff\",\"p_{T}(jet)\",2,30,300)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_jeteta_eff\",\"#eta(jet)\",2,-2.4,2.4)");


  // photon pt and eta
  gROOT->ProcessLine(
"combineMC_eff(\"h_pt_eff_EB\",\"p_{T}(#gamma) in barrel\",2,85,200)"
);

  gROOT->ProcessLine(
"combineMC_eff(\"h_pt_eff_EE\",\"p_{T}(#gamma) in endcap\",2,85,200)"
);

 
  gROOT->ProcessLine(
"combineMC_eff(\"h_eta_eff_85_95\",\"#eta(#gamma)\",2,-2.5,2.5)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_eta_eff_95_110\",\"#eta(#gamma)\",2,-2.5,2.5)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_eta_eff_110_130\",\"#eta(#gamma)\",2,-2.5,2.5)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_eta_eff_130_160\",\"#eta(#gamma)\",2,-2.5,2.5)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_eta_eff_160_200\",\"#eta(#gamma)\",2,-2.5,2.5)");


  gROOT->ProcessLine(
"combineMC_eff(\"h_eta_eff_allpt\",\"#eta(#gamma)\",2,-2.5,2.5)");


  // pstar in barrel
  gROOT->ProcessLine(
 		     "combineMC_eff(\"h_pstar_EB_85_95\",\"p^{*} in barrel\",2,85,300)");

  gROOT->ProcessLine(
		     "combineMC_eff(\"h_pstar_EB_95_110\",\"p^{*} in barrel\",2,85,300)");

  gROOT->ProcessLine(
		     "combineMC_eff(\"h_pstar_EB_110_130\",\"p^{*} in barrel\",2,85,300)");

  gROOT->ProcessLine(
		     "combineMC_eff(\"h_pstar_EB_130_160\",\"p^{*} in barrel\",2,85,300)");

  gROOT->ProcessLine(
		     "combineMC_eff(\"h_pstar_EB_160_200\",\"p^{*} in barrel\",2,85,300)");

  gROOT->ProcessLine(
		     "combineMC_eff(\"h_pstar_EB_allpt\",\"p^{*} in barrel\",4,85,300)");

  gROOT->ProcessLine(
		     "combineMC_eff(\"h_pstar_COM3D_EB_allpt\",\"p^{*} in barrel\",4,85,300)");

  gROOT->ProcessLine(
		     "combineMC_eff(\"h_pstar_COMZ_EB_allpt\",\"p^{*} in barrel\",4,85,300)");

  // pstar in endcap
  gROOT->ProcessLine(
		     "combineMC_eff(\"h_pstar_EE_85_95\",\"p^{*} in endcap\",2,85,300)");

  gROOT->ProcessLine(
		     "combineMC_eff(\"h_pstar_EE_95_110\",\"p^{*} in endcap\",2,85,300)");

  gROOT->ProcessLine(
		     "combineMC_eff(\"h_pstar_EE_110_130\",\"p^{*} in endcap\",2,85,300)");

  gROOT->ProcessLine(
		     "combineMC_eff(\"h_pstar_EE_130_160\",\"p^{*} in endcap\",2,85,300)");

  gROOT->ProcessLine(
		     "combineMC_eff(\"h_pstar_EE_160_200\",\"p^{*} in endcap\",2,85,300)");

  gROOT->ProcessLine(
		     "combineMC_eff(\"h_pstar_EE_allpt\",\"p^{*} in endcap\",4,85,300)");


  gROOT->ProcessLine(
		     "combineMC_eff(\"h_pstar_COM3D_EE_allpt\",\"p^{*} in endcap\",4,85,300)");

  gROOT->ProcessLine(
		     "combineMC_eff(\"h_pstar_COMZ_EE_allpt\",\"p^{*} in endcap\",4,85,300)");

  // yBoost in barrel
  gROOT->ProcessLine(
"combineMC_eff(\"h_yB_EB_85_95\",\"y^{boost} in barrel\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_yB_EB_95_110\",\"y^{boost} in barrel\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_yB_EB_110_130\",\"y^{boost} in barrel\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_yB_EB_130_160\",\"y^{boost} in barrel\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_yB_EB_160_200\",\"y^{boost} in barrel\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_yB_EB_allpt\",\"y^{boost} in barrel\",2,-2.4,2.4)");

  // yboost in endcap
  gROOT->ProcessLine(
"combineMC_eff(\"h_yB_EE_85_95\",\"y^{boost} in endcap\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_yB_EE_95_110\",\"y^{boost} in endcap\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_yB_EE_110_130\",\"y^{boost} in endcap\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_yB_EE_130_160\",\"y^{boost} in endcap\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_yB_EE_160_200\",\"y^{boost} in endcap\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_yB_EE_allpt\",\"y^{boost} in endcap\",2,-2.4,2.4)");


  // costheta in barrel
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_EB_85_95\",\"cos#theta^{*} in barrel\",2,0.0,1.0)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_EB_95_110\",\"cos#theta^{*} in barrel\",2,0.0,1.0)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_EB_110_130\",\"cos#theta^{*} in barrel\",2,0.0,1.0)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_EB_130_160\",\"cos#theta^{*} in barrel\",2,0.0,1.0)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_EB_160_200\",\"cos#theta^{*} in barrel\",2,0.0,1.0)");
 
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_EB_allpt\",\"cos#theta^{*} in barrel\",4,0.0,1.0)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COM3D_EB_allpt\",\"cos#theta^{*} in barrel\",4,0.0,1.0)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EB_allpt\",\"cos#theta^{*} in barrel\",4,0.0,1.0)");

  // costheta in endcap
  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_EE_85_95\",\"cos#theta^{*} in endcap\",2,0.0,1.0)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_EE_95_110\",\"cos#theta^{*} in endcap\",2,0.0,1.0)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_EE_110_130\",\"cos#theta^{*} in endcap\",2,0.0,1.0)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_EE_130_160\",\"cos#theta^{*} in endcap\",2,0.0,1.0)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_EE_160_200\",\"cos#theta^{*} in endcap\",2,0.0,1.0)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_EE_allpt\",\"cos#theta^{*} in endcap\",4,0.0,1.0)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COM3D_EE_allpt\",\"cos#theta^{*} in endcap\",4,0.0,1.0)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_cost_COMZ_EE_allpt\",\"cos#theta^{*} in endcap\",4,0.0,1.0)");


   // deltaphi in barrel
  gROOT->ProcessLine(
"combineMC_eff(\"h_dphi_EB_85_95\",\"#Delta #phi in barrel\",2)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_dphi_EB_95_110\",\"#Delta #phi in barrel\",2)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_dphi_EB_110_130\",\"#Delta #phi in barrel\",2)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_dphi_EB_130_160\",\"#Delta #phi in barrel\",2)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_dphi_EB_160_200\",\"#Delta #phi in barrel\",2)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_dphi_EB_allpt\",\"#Delta #phi in barrel\",4)");

  // deltaphi in endcap
  gROOT->ProcessLine(
"combineMC_eff(\"h_dphi_EE_85_95\",\"#Delta #phi in endcap\",2)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_dphi_EE_95_110\",\"#Delta #phi in endcap\",2)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_dphi_EE_110_130\",\"#Delta #phi in endcap\",2)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_dphi_EE_130_160\",\"#Delta #phi in endcap\",2)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_dphi_EE_160_200\",\"#Delta #phi in endcap\",2)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_dphi_EE_allpt\",\"#Delta #phi in endcap\",4)");


 
  // y at the center of mass frame in barrel
  gROOT->ProcessLine(
"combineMC_eff(\"h_ystar_EB_85_95\",\"y^{*} in barrel\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_ystar_EB_95_110\",\"y^{*} in barrel\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_ystar_EB_110_130\",\"y^{*} in barrel\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_ystar_EB_130_160\",\"y^{*} in barrel\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_ystar_EB_160_200\",\"y^{*} in barrel\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_ystar_EB_allpt\",\"y^{*} in barrel\",4,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_ystar_COM3D_EB_allpt\",\"y^{*} in barrel\",4,-2.4,2.4)");

  // y at the center of mass frame in endcap
  gROOT->ProcessLine(
"combineMC_eff(\"h_ystar_EE_85_95\",\"y^{*} in endcap\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_ystar_EE_95_110\",\"y^{*} in endcap\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_ystar_EE_110_130\",\"y^{*} in endcap\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_ystar_EE_130_160\",\"y^{*} in endcap\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_ystar_EE_160_200\",\"y^{*} in endcap\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_ystar_EE_allpt\",\"y^{*} in endcap\",4,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_ystar_COM3D_EE_allpt\",\"y^{*} in endcap\",4,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff(\"h_ystar_COMZ_EE_allpt\",\"y^{*} in endcap\",4,-2.4,2.4)");

  myWatch->Stop();
  cout << myWatch->RealTime() << " seconds has passed! " << endl; 

}
