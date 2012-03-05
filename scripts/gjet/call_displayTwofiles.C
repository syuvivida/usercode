
void call_displayTwofiles()
{
  TStopwatch* myWatch = new TStopwatch();
  myWatch->Start();
  // first set of root files
  // first histogram requires recreating root files
  gROOT->ProcessLine(".L /afs/cern.ch/user/s/syu/scripts/displayTwofiles.C");
 
  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_genciso_all\",\"\",0,50,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_gentiso_all\",\"\",0,50,true,true)");


  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_ptpho_all\",\"\",85,200,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_ptjet_all\",\"\",30,400,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_jetpt_eff_0_all\",\"\",30,400,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_etapho_all\",\"\",-2.5,2.5)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_etapho_all\",\"\",-2.3,2.3,false,false,\"h_etapho_restricted\")");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_etajet_all\",\"\",-2.4,2.4)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_jeteta_eff_0_all\",\"\",-2.4,2.4)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_sieie_leadingEB_all\",\"\",0,0.02,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_eciso_leadingEB_all\",\"\",-1.5,50,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_hciso_leadingEB_all\",\"\",-1.5,50,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_tkiso_leadingEB_all\",\"\",-1.5,50,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_hovere_leadingEB_all\",\"\",0,0.5,true,true)");
  
  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_pixel_leadingEB_all\")");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_sieie_leadingEE_all\",\"\",0,0.08,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_eciso_leadingEE_all\",\"\",-1.5,50,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_hciso_leadingEE_all\",\"\",-1.5,50,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_tkiso_leadingEE_all\",\"\",-1.5,50,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_hovere_leadingEE_all\",\"\",0,0.5,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_pixel_leadingEE_all\")");
 
  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_ptratio_EB_0_all\",\"\",-9999,-9999,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_ptratio_EE_0_all\",\"\",-9999,-9999,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_nvtx_EB_0_all\")");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_nvtx_EE_0_all\")");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_dphi_EB_allpt_0_all\",\"\",-9999,-9999,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_dphi_EB_85_95_0_all\",\"\",-9999,-9999,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_dphi_EE_allpt_0_all\",\"\",-9999,-9999,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_dphi_EE_85_95_0_all\",\"\",-9999,-9999,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_cost_COMZ_EB_allpt_0_all\",\"\",-9999,-9999,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_cost_COMZ_EE_allpt_0_all\",\"\",-9999,-9999,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_pstar_COMZ_EB_allpt_0_all\",\"\",85,500,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_pstar_COMZ_EE_allpt_0_all\",\"\",85,500,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_yB_EB_allpt_0_all\",\"\",-2,2)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_yB_EE_allpt_0_all\",\"\",-2,2)");


  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_ystar_COMZ_EB_allpt_0_all\",\"\",-3.0,3.0)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_ystar_COMZ_EE_allpt_0_all\",\"\",-3.0,3.0)");


  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_njetraw_EB_all\",\"\",-0.5,10.5,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_njetraw_EE_all\",\"\",-0.5,10.5,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_ptclosestjet_EB_all\",\"\",30,400,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_ptclosestjet_EE_all\",\"\",30,400,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_etaclosestjet_EB_all\",\"\",-2.4,2.4,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_etaclosestjet_EE_all\",\"\",-2.4,2.4,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_dRLeadingPhoJet_EB_all\",\"\",-9999,-9999,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_dRLeadingPhoJet_EE_all\",\"\",-9999,-9999,false,true)");


  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_dRclosestjet_EB_all\",\"\",-9999,-9999,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_dRclosestjet_EE_all\",\"\",-9999,-9999,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_alljetpt_eff_EB_0_all\",\"\", 30,400,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_alljetpt_eff_EE_0_all\",\"\", 30,400,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_alljeteta_eff_EB_0_all\",\"\",-2.4,2.4,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_alljeteta_eff_EE_0_all\",\"\",-2.4,2.4,false,true)");


  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_cost_COMZ_alljet_EB_0_all\",\"\",-9999,-9999,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_cost_COMZ_alljet_EE_0_all\",\"\",-9999,-9999,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_pstar_COMZ_alljet_EB_0_all\",\"\",85,500,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_pstar_COMZ_alljet_EE_0_all\",\"\",85,500,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_ystar_COMZ_alljet_EB_0_all\",\"\",-3.0,3.0)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_ystar_COMZ_alljet_EE_0_all\",\"\",-3.0,3.0)");


  myWatch->Stop();
  cout << myWatch->RealTime() << " seconds has passed! " << endl; 

}
