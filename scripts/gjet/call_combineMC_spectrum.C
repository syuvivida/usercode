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
			  "combineMC_spectrum(\"h_genciso\",\"Gen-level calorimeter isolation [GeV]\",\"%s\",true)", fileName.data()));
  gROOT->ProcessLine(Form(
			  "combineMC_spectrum(\"h_gentiso\",\"Gen-level track isolation [GeV]\",\"%s\",true)", fileName.data()));


  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_jetpt_eff_EB_0\", \"p_{T}(jet) with barrel photon [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_jetpt_eff_EE_0\", \"p_{T}(jet) with endcap photon [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_jeteta_eff_EB_0\",\"#eta(jet) with barrel photon\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_jeteta_eff_EE_0\",\"#eta(jet) with endcap photon\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_phopt_eff_EB_0\", \"p_{T}(barrel #gamma) [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_phopt_eff_EE_0\", \"p_{T}(endcap #gamma) [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_phoeta_eff_allpt_0\", \"#eta(#gamma)\",\"%s\",true)", fileName.data()));


  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_sieie_EB\",\"Barrel #sigma_{i#eta i#eta}\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_eciso_EB\",\"Barrel ISO_{ECAL} [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_hciso_EB\",\"Barrel ISO_{HCAL} [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_tkiso_EB\",\"Barrel ISO_{TRK} [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_hovere_EB\",\"Barrel H/E\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_pixel_EB\",\"Barrel hasPixelSeed\",\"%s\",true)", fileName.data()));
 

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_sieie_EE\",\"Endcap #sigma_{i#eta i#eta}\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_eciso_EE\",\"Endcap ISO_{ECAL} [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_hciso_EE\",\"Endcap ISO_{HCAL} [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_tkiso_EE\",\"Endcap ISO_{TRK} [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_hovere_EE\",\"Endcap H/E\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_pixel_EE\",\"Endcap hasPixelSeed\",\"%s\",true)", fileName.data()));
 
  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_ptratio_EB_0\",\"Barrel p_{T}(jet)/p_{T}(#gamma)\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_ptratio_EE_0\",\"Endcap p_{T}(jet)/p_{T}(#gamma)\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_nvtx_EB_0\",\"Barrel number of good vertices\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_nvtx_EE_0\",\"Endcap number of good vertices\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_gendphi_EB_0\",\"Generator-level #Delta #phi(barrel #gamma,jet) [radian]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_gendphi_EB_1\",\"Generator-level #Delta #phi(barrel #gamma,jet) [radian]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_dphi_EB_allpt_0\",\"#Delta #phi(barrel #gamma,jet) [radian]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_dphi_EB_85_95_0\",\"#Delta #phi(barrel #gamma,jet) [radian]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_cost_COMZ_EB_allpt_0\",\"Barrel cos#theta^{*}\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_pstar_COMZ_EB_allpt_0\",\"Barrel p^{*} [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_yB_EB_allpt_0\",\"Barrel 0.5#times(y_{1}+y_{2})\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_ystar_COMZ_EB_allpt_0\",\"Barrel y^{*}\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_gendphi_EE_0\",\"Generator-level #Delta #phi(endcap #gamma,jet) [radian]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_gendphi_EE_1\",\"Generator-level #Delta #phi(endcap #gamma,jet) [radian]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_dphi_EE_allpt_0\",\"#Delta #phi(endcap #gamma,jet) [radian]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_dphi_EE_85_95_0\",\"#Delta #phi(endcap #gamma,jet) [radian]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_cost_COMZ_EE_allpt_0\",\"Endcap cos#theta^{*}\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_pstar_COMZ_EE_allpt_0\",\"Endcap p^{*} [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_yB_EE_allpt_0\",\"Endcap 0.5#times(y_{1}+y_{2})\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_ystar_COMZ_EE_allpt_0\",\"Endcap y^{*}\",\"%s\",true)", fileName.data()));


  gROOT->ProcessLine(Form(
			  "combineMC_spectrum(\"h_ngenjet\",\"Generator-level jet multiplicity\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
			  "combineMC_spectrum(\"h_nrecjet\",\"Reconstruction-level jet multiplicity\",\"%s\",true)", fileName.data()));


  gROOT->ProcessLine(Form(
			  "combineMC_spectrum(\"h_njet_EB_0\",\"Fiducial jet multiplicity with barrel photon\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
			  "combineMC_spectrum(\"h_njet_EE_0\",\"Fiducial jet multiplicity with endcap photon\",\"%s\",true)", fileName.data()));


  gROOT->ProcessLine(Form(
			  "combineMC_spectrum(\"h_ptclosestjet_EB\",\"p_{T} of the closest jet to barrel photon\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
			  "combineMC_spectrum(\"h_ptclosestjet_EE\",\"p_{T} of the closest jet to endcap photon\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
			  "combineMC_spectrum(\"h_etaclosestjet_EB\",\"#eta of the closest jet to barrel photon\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
			  "combineMC_spectrum(\"h_etaclosestjet_EE\",\"#eta of the closest jet to endcap photon\",\"%s\",true)", fileName.data()));


  gROOT->ProcessLine(Form(
			  "combineMC_spectrum(\"h_dRclosestjet_EB\",\"#Delta R of the closest jet to barrel photon\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
			  "combineMC_spectrum(\"h_dRclosestjet_EE\",\"#Delta R of the closest jet to endcap photon\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
			  "combineMC_spectrum(\"h_dRLeadingPhoJet_EB\",\"#Delta R of the leading photon (barrel) and jet\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
			  "combineMC_spectrum(\"h_dRLeadingPhoJet_EE\",\"#Delta R of the leading photon (endcap) and jet\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
			  "combineMC_spectrum(\"h_sumjetpt_eff_EB_0\",\"p_{T}(jet^{all}) with barrel photon\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
			  "combineMC_spectrum(\"h_sumjetpt_eff_EE_0\",\"p_{T}(jet^{all}) with endcap photon\",\"%s\",true)", fileName.data()));


  gROOT->ProcessLine(Form(
			  "combineMC_spectrum(\"h_sumjeteta_eff_EB_0\",\"#eta(jet^{all}) with barrel photon\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
			  "combineMC_spectrum(\"h_sumjeteta_eff_EE_0\",\"#eta(jet^{all}) with endcap photon\",\"%s\",true)", fileName.data()));


  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_cost_sumJet_EB_allpt_0\",\"Barrel cos#theta^{*} using all jets\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_pstar_sumJet_EB_allpt_0\",\"Barrel p^{*} using all jets [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_ystar_sumJet_EB_allpt_0\",\"Barrel y^{*} using all jets\",\"%s\",true)", fileName.data()));


  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_cost_sumJet_EE_allpt_0\",\"Endcap cos#theta^{*} using all jets\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_pstar_sumJet_EE_allpt_0\",\"Endcap p^{*} using all jets [GeV]\",\"%s\",true)", fileName.data()));

  gROOT->ProcessLine(Form(
"combineMC_spectrum(\"h_ystar_sumJet_EE_allpt_0\",\"Endcap y^{*} using all jets\",\"%s\",true)", fileName.data()));

  myWatch->Stop();
  cout << myWatch->RealTime() << " seconds has passed! " << endl; 

}
