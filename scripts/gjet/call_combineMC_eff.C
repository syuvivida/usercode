{
  TStopwatch* myWatch = new TStopwatch();
  myWatch->Start();
  gROOT->ProcessLine(".L combineMC_eff.C++");

  // jetpt/photon pt
  gROOT->ProcessLine(
"combineMC_eff\(\"h_ptratio_EB\",\"p_{T}(jet)/p_{T}(#gamma) in barrel\",2)"
);

  gROOT->ProcessLine(
"combineMC_eff\(\"h_ptratio_EE\",\"p_{T}(jet)/p_{T}(#gamma) in endcap\",2)"
);

  //nvertex

 
  gROOT->ProcessLine(
"combineMC_eff\(\"h_nvtx_eff_EB\",\"number of vertex in barrel\")"
);

  gROOT->ProcessLine(
"combineMC_eff\(\"h_nvtx_eff_EE\",\"number of vertex in endcap\")"
);

 
   
  // jet pt and eta
  gROOT->ProcessLine(
"combineMC_eff\(\"h_jetpt_eff\",\"p_{T}(jet)\",2,30,300)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_jeteta_eff\",\"#eta(jet)\",2,-2.4,2.4)");


  // photon pt and eta
  gROOT->ProcessLine(
"combineMC_eff\(\"h_pt_eff_EB\",\"p_{T}(#gamma) in barrel\",2,85,200)"
);

  gROOT->ProcessLine(
"combineMC_eff\(\"h_pt_eff_EE\",\"p_{T}(#gamma) in endcap\",2,85,200)"
);

 
  gROOT->ProcessLine(
"combineMC_eff\(\"h_eta_eff_85_95\",\"#eta(#gamma)\",2,-2.5,2.5)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_eta_eff_95_110\",\"#eta(#gamma)\",2,-2.5,2.5)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_eta_eff_110_130\",\"#eta(#gamma)\",2,-2.5,2.5)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_eta_eff_130_160\",\"#eta(#gamma)\",2,-2.5,2.5)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_eta_eff_160_200\",\"#eta(#gamma)\",2,-2.5,2.5)");

  // pstar in barrel
  gROOT->ProcessLine(
 		     "combineMC_eff\(\"h_pstar_EB_85_95\",\"p^{*}\",2,85,300)");

  gROOT->ProcessLine(
		     "combineMC_eff\(\"h_pstar_EB_95_110\",\"p^{*}\",2,85,300)");

  gROOT->ProcessLine(
		     "combineMC_eff\(\"h_pstar_EB_110_130\",\"p^{*}\",2,85,300)");

  gROOT->ProcessLine(
		     "combineMC_eff\(\"h_pstar_EB_130_160\",\"p^{*}\",2,85,300)");

  gROOT->ProcessLine(
		     "combineMC_eff\(\"h_pstar_EB_160_200\",\"p^{*}\",2,85,300)");

  // pstar in endcap
  gROOT->ProcessLine(
		     "combineMC_eff\(\"h_pstar_EE_85_95\",\"p^{*}\",2,85,300)");

  gROOT->ProcessLine(
		     "combineMC_eff\(\"h_pstar_EE_95_110\",\"p^{*}\",2,85,300)");

  gROOT->ProcessLine(
		     "combineMC_eff\(\"h_pstar_EE_110_130\",\"p^{*}\",2,85,300)");

  gROOT->ProcessLine(
		     "combineMC_eff\(\"h_pstar_EE_130_160\",\"p^{*}\",2,85,300)");

  gROOT->ProcessLine(
		     "combineMC_eff\(\"h_pstar_EE_160_200\",\"p^{*}\",2,85,300)");


  // yBoost in barrel
  gROOT->ProcessLine(
"combineMC_eff\(\"h_yB_EB_85_95\",\"y^{boost}\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_yB_EB_95_110\",\"y^{boost}\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_yB_EB_110_130\",\"y^{boost}\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_yB_EB_130_160\",\"y^{boost}\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_yB_EB_160_200\",\"y^{boost}\",2,-2.4,2.4)");

  // yboost in endcap
  gROOT->ProcessLine(
"combineMC_eff\(\"h_yB_EE_85_95\",\"y^{boost}\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_yB_EE_95_110\",\"y^{boost}\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_yB_EE_110_130\",\"y^{boost}\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_yB_EE_130_160\",\"y^{boost}\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_yB_EE_160_200\",\"y^{boost}\",2,-2.4,2.4)");

  // costheta in barrel
  gROOT->ProcessLine(
"combineMC_eff\(\"h_cost_EB_85_95\",\"cos#theta^{*}\",2,0.0,1.0)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_cost_EB_95_110\",\"cos#theta^{*}\",2,0.0,1.0)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_cost_EB_110_130\",\"cos#theta^{*}\",2,0.0,1.0)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_cost_EB_130_160\",\"cos#theta^{*}\",2,0.0,1.0)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_cost_EB_160_200\",\"cos#theta^{*}\",2,0.0,1.0)");

  // costheta in endcap
  gROOT->ProcessLine(
"combineMC_eff\(\"h_cost_EE_85_95\",\"cos#theta^{*}\",2,0.0,1.0)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_cost_EE_95_110\",\"cos#theta^{*}\",2,0.0,1.0)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_cost_EE_110_130\",\"cos#theta^{*}\",2,0.0,1.0)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_cost_EE_130_160\",\"cos#theta^{*}\",2,0.0,1.0)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_cost_EE_160_200\",\"cos#theta^{*}\",2,0.0,1.0)");

   // deltaphi in barrel
  gROOT->ProcessLine(
"combineMC_eff\(\"h_dphi_EB_85_95\",\"#Delta #phi\",2)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_dphi_EB_95_110\",\"#Delta #phi\",2)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_dphi_EB_110_130\",\"#Delta #phi\",2)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_dphi_EB_130_160\",\"#Delta #phi\",2)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_dphi_EB_160_200\",\"#Delta #phi\",2)");

  // deltaphi in endcap
  gROOT->ProcessLine(
"combineMC_eff\(\"h_dphi_EE_85_95\",\"#Delta #phi\",2)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_dphi_EE_95_110\",\"#Delta #phi\",2)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_dphi_EE_110_130\",\"#Delta #phi\",2)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_dphi_EE_130_160\",\"#Delta #phi\",2)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_dphi_EE_160_200\",\"#Delta #phi\",2)");


 
  // y at the center of mass frame in barrel
  gROOT->ProcessLine(
"combineMC_eff\(\"h_yCOM_EB_85_95\",\"y^{boost}\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_yCOM_EB_95_110\",\"y^{boost}\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_yCOM_EB_110_130\",\"y^{boost}\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_yCOM_EB_130_160\",\"y^{boost}\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_yCOM_EB_160_200\",\"y^{boost}\",2,-2.4,2.4)");

  // y at the center of mass frame in endcap
  gROOT->ProcessLine(
"combineMC_eff\(\"h_yCOM_EE_85_95\",\"y^{boost}\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_yCOM_EE_95_110\",\"y^{boost}\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_yCOM_EE_110_130\",\"y^{boost}\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_yCOM_EE_130_160\",\"y^{boost}\",2,-2.4,2.4)");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_yCOM_EE_160_200\",\"y^{boost}\",2,-2.4,2.4)");


  myWatch->Stop();
  cout << myWatch->RealTime() << " seconds has passed! " << endl; 

}
