
void call_displayTwofiles()
{
  TStopwatch* myWatch = new TStopwatch();
  myWatch->Start();
  // first set of root files
  // first histogram requires recreating root files
  gROOT->ProcessLine(".L /afs/cern.ch/user/s/syu/scripts/displayTwofiles.C");
 
  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_genciso_all\",\"\",0,50,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_gentiso_all\",\"\",0,50,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_jeteta_eff_EB_0_all\",\"\",-2.4,2.4)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_jeteta_eff_EE_0_all\",\"\",-2.4,2.4)");
  
  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_phoeta_eff_allpt_0_all\",\"\",-2.3,2.3)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_jetpt_eff_EB_0_all\",\"\", 30,300)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_jetpt_eff_EE_0_all\",\"\", 30,300)");
  
  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_phopt_eff_EB_0_all\",\"\",40,200)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_phopt_eff_EE_0_all\",\"\",40,200)");


  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_sieie_EB_all\",\"\",0,0.02,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_eciso_EB_all\",\"\",-1.5,50,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_hciso_EB_all\",\"\",-1.5,50,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_tkiso_EB_all\",\"\",-1.5,50,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_hovere_EB_all\",\"\",0,0.5,true,true)");
  
  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_pixel_EB_all\")");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_sieie_EE_all\",\"\",0,0.08,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_eciso_EE_all\",\"\",-1.5,50,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_hciso_EE_all\",\"\",-1.5,50,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_tkiso_EE_all\",\"\",-1.5,50,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_hovere_EE_all\",\"\",0,0.5,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_pixel_EE_all\")");
 
  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_ptratio_EB_0_all\",\"\",-9999,-9999,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_ptratio_EE_0_all\",\"\",-9999,-9999,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_nvtx_EB_0_all\")");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_nvtx_EE_0_all\")");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_gendphi_EB_0_all\",\"\",-9999,-9999,true,true)");
  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_gendphi_EB_1_all\",\"\",-9999,-9999,true,true)");
  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_dphi_EB_allpt_0_all\",\"\",-9999,-9999,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_dphi_EB_85_95_0_all\",\"\",-9999,-9999,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_gendphi_EE_0_all\",\"\",-9999,-9999,true,true)");
  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_gendphi_EE_1_all\",\"\",-9999,-9999,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_dphi_EE_allpt_0_all\",\"\",-9999,-9999,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_dphi_EE_85_95_0_all\",\"\",-9999,-9999,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_cost_COMZ_EB_allpt_0_all\",\"\",-9999,-9999,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_cost_COMZ_EE_allpt_0_all\",\"\",-9999,-9999,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_pstar_COMZ_EB_allpt_0_all\",\"\",85,500,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_pstar_COMZ_EE_allpt_0_all\",\"\",85,500,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_yB_EB_allpt_0_all\",\"\",-2,2)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_yB_EE_allpt_0_all\",\"\",-2,2)");


  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_ystar_COMZ_EB_allpt_0_all\",\"\",-3.0,3.0)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_ystar_COMZ_EE_allpt_0_all\",\"\",-3.0,3.0)");


  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_njet_EB_0_all\",\"\",-0.5,10.5,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_njet_EE_0_all\",\"\",-0.5,10.5,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_ngenjet_all\",\"\",-0.5,10.5,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_nrecjet_all\",\"\",-0.5,10.5,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_ptclosestjet_EB_all\",\"\",30,400,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_ptclosestjet_EE_all\",\"\",30,400,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_etaclosestjet_EB_all\",\"\",-2.4,2.4,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_etaclosestjet_EE_all\",\"\",-2.4,2.4,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_dRLeadingPhoJet_EB_all\",\"\",-9999,-9999,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_dRLeadingPhoJet_EE_all\",\"\",-9999,-9999,false,true)");


  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_dRclosestjet_EB_all\",\"\",-9999,-9999,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_dRclosestjet_EE_all\",\"\",-9999,-9999,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_sumjetpt_eff_EB_0_all\",\"\", 30,300,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_sumjetpt_eff_EE_0_all\",\"\", 30,300,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_sumjeteta_eff_EB_0_all\",\"\",-2.4,2.4,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_sumjeteta_eff_EE_0_all\",\"\",-2.4,2.4,false,true)");


  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_cost_sumJet_EB_allpt_0_all\",\"\",-9999,-9999,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_cost_sumJet_EE_allpt_0_all\",\"\",-9999,-9999,false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_pstar_sumJet_EB_allpt_0_all\",\"\",85,500,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_pstar_sumJet_EE_allpt_0_all\",\"\",85,500,true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_ystar_sumJet_EB_allpt_0_all\",\"\",-3.0,3.0)");

  gROOT->ProcessLine("displayTwofiles(\"combined_pythia_fragmentation.root\",\"combined_pythia_direct.root\",\"h_ystar_sumJet_EE_allpt_0_all\",\"\",-3.0,3.0)");


  myWatch->Stop();
  cout << myWatch->RealTime() << " seconds has passed! " << endl; 

}
