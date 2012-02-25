void call_combineMC_spectrum(std::string fileName)
{
  TStopwatch* myWatch = new TStopwatch();
  myWatch->Start();
  gROOT->ProcessLine(".L combineMC_spectrum.C++");

  // first set of root files
  // first histogram requires recreating root files
  gROOT->ProcessLine(Form(
			  "combineMC_spectrum(\"h_pthat\",\"#hat{p_{T}} [GeV]\",\"%s\")", fileName.data()));
  gROOT->ProcessLine(Form(
			  "combineMC_spectrum(\"h_ptpho\",\"p_{T}(#gamma) [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_ptjet\",\"p_{T}(jet) [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_etapho\",\"#eta(#gamma)\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_etajet\",\"#eta(jet)\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_sieie_leadingEB\",\"Barrel #sigma_{i#eta i#eta}\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_eciso_leadingEB\",\"Barrel ISO_{ECAL} [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_hciso_leadingEB\",\"Barrel ISO_{HCAL} [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_tkiso_leadingEB\",\"Barrel ISO_{TRK} [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_hovere_leadingEB\",\"Barrel H/E\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_pixel_leadingEB\",\"Barrel hasPixelSeed\",\"%s\",true)", fileName.data()));
 

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_sieie_leadingEE\",\"Endcap #sigma_{i#eta i#eta}\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_eciso_leadingEE\",\"Endcap ISO_{ECAL} [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_hciso_leadingEE\",\"Endcap ISO_{HCAL} [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_tkiso_leadingEE\",\"Endcap ISO_{TRK} [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_hovere_leadingEE\",\"Endcap H/E\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_pixel_leadingEE\",\"Endcap hasPixelSeed\",\"%s\",true)", fileName.data()));
 
  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_ptratio_EB_0\",\"Barrel p_{T}(jet)/p_{T}(#gamma)\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_ptratio_EE_0\",\"Endcap p_{T}(jet)/p_{T}(#gamma)\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_nvtx_EB_0\",\"Barrel number of good vertices\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_nvtx_EE_0\",\"Endcap number of good vertices\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_dphi_EB_allpt_0\",\"#Delta #phi(barrel #gamma,jet) [radian]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_cost_COMZ_EB_allpt_0\",\"Barrel cos#theta^{*}\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_pstar_COMZ_EB_allpt_0\",\"Barrel p^{*} [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_yB_EB_allpt_0\",\"Barrel 0.5#times(y_{1}+y_{2})\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_ystar_COMZ_EB_allpt_0\",\"Barrel y^{*}\",\"%s\",true)", fileName.data()));


  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_dphi_EE_allpt_0\",\"#Delta #phi(endcap #gamma,jet) [radian]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_cost_COMZ_EE_allpt_0\",\"Endcap cos#theta^{*}\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_pstar_COMZ_EE_allpt_0\",\"Endcap p^{*} [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_yB_EE_allpt_0\",\"Endcap 0.5#times(y_{1}+y_{2})\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_ystar_COMZ_EE_allpt_0\",\"Endcap y^{*}\",\"%s\",true)", fileName.data()));

//   gROOT->ProcessLine(Form(
// "combineMC_spectrum(\"h_njet_EB_0\",\"Jet multiplicity with barrel photon\",\"%s\",true)", fileName.data()));

//   gROOT->ProcessLine(Form(
// "combineMC_spectrum(\"h_njet_EE_0\",\"Jet multiplicity with endcap photon\",\"%s\",true)", fileName.data()));


  myWatch->Stop();
  cout << myWatch->RealTime() << " seconds has passed! " << endl; 

}
