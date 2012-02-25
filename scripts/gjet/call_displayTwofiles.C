
void call_displayTwofiles()
{
  TStopwatch* myWatch = new TStopwatch();
  myWatch->Start();
  // first set of root files
  // first histogram requires recreating root files
  gROOT->ProcessLine(".L /afs/cern.ch/user/s/syu/scripts/displayTwofiles.C");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_ptpho_all\",\"\",85,200)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_ptjet_all\",\"\",30,400,\"\",false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_etapho_all\",\"\",-2.5,2.5)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_etapho_all\",\"\",-2.3,2.3,\"h_etapho_restricted\")");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_etajet_all\",\"\",-2.4,2.4)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_sieie_leadingEB_all\",\"\",0,0.02,\"\",true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_eciso_leadingEB_all\",\"\",-1.5,50,\"\",true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_hciso_leadingEB_all\",\"\",-1.5,50,\"\",true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_tkiso_leadingEB_all\",\"\",-1.5,50,\"\",true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_hovere_leadingEB_all\",\"\",0,0.5,\"\",true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_pixel_leadingEB_all\"");
 
  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_sieie_leadingEE_all\",\"\",0,0.08,\"\",true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_eciso_leadingEE_all\",\"\",-1.5,50,\"\",true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_hciso_leadingEE_all\",\"\",-1.5,50,\"\",true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_tkiso_leadingEE_all\",\"\",-1.5,50,\"\",true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_hovere_leadingEE_all\",\"\",0,0.5,\"\",true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_pixel_leadingEE_all\"");
 
  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_ptratio_EB_0_all\",\"\",-9999,-9999,\"\",false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_ptratio_EE_0_all\",\"\",-9999,-9999,\"\",false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_nvtx_EB_0_all\"");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_nvtx_EE_0_all\"");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_dphi_EB_allpt_0_all\",\"\",-9999,-9999,\"\",true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_cost_COMZ_EB_allpt_0_all\",\"\",-9999,-9999,\"\",false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_pstar_COMZ_EB_allpt_0_all\",\"\",85,500,\"\",true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_yB_EB_allpt_0_all\",\"\",-2,2)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_ystar_COMZ_EB_allpt_0_all\",\"\",-3.0,3.0)");



  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_dphi_EE_allpt_0_all\",\"\",-9999,-9999,\"\",true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_cost_COMZ_EE_allpt_0_all\",\"\",-9999,-9999,\"\",false,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_pstar_COMZ_EE_allpt_0_all\",\"\",85,500,\"\",true,true)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_yB_EE_allpt_0_all\",\"\",-2,2)");

  gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_ystar_COMZ_EE_allpt_0_all\",\"\",-3.0,3.0)");



//   gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_njet_EB_0_all\",\"\",-9999,-9999,true)");

//   gROOT->ProcessLine("displayTwofiles(\"combined_madgraph.root\",\"combined_pythia.root\",\"h_njet_EE_0_all\",\"\",-9999,-9999,true)");


myWatch->Stop();
  cout << myWatch->RealTime() << " seconds has passed! " << endl; 

}
