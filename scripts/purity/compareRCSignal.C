#include "/afs/cern.ch/user/s/syu/scripts/displayTwofiles.C"

void compareRCSignal()
{
  // barrel
  displayTwofiles("randomConeHists_Barrel_v5.root","template_comb3Iso.root",
		  "hsumisomc","h_EB_comb3Iso_sig_pt15",
		  "Iso_{ECAL}+Iso_{HCAL}+Iso_{TRK} (GeV)","","RCMC_signalMCpt15_EB_comb3Iso");


  displayTwofiles("randomConeHists_Barrel_v5.root","template_ecalIso.root",
		  "hecisomc","h_EB_ecalIso_sig_pt15",
		  "Iso_{ECAL} (GeV)","","RCMC_signalMCpt15_EB_ecalIso");
  

  displayTwofiles("randomConeHists_Barrel_v5.root","template_hcalIso.root",
		  "hhcisomc","h_EB_hcalIso_sig_pt15",
		  "Iso_{HCAL} (GeV)","","RCMC_signalMCpt15_EB_hcalIso");
  

  displayTwofiles("randomConeHists_Barrel_v5.root","template_hcalIso.root",
		  "hhcisomc","h_EB_hcalIso_sig_pt15",
		  "Iso_{HCAL} (GeV)","","RCMC_signalMCpt15_EB_hcalIso_log",
		  true);
 

  displayTwofiles("randomConeHists_Barrel_v5.root","template_trkIso.root",
		  "htkisomc","h_EB_trkIso_sig_pt15",
		  "Iso_{TRK} (GeV)","","RCMC_signalMCpt15_EB_trkIso");
  

  displayTwofiles("randomConeHists_Barrel_v5.root","template_trkIso.root",
		  "htkisomc","h_EB_trkIso_sig_pt15",
		  "Iso_{TRK} (GeV)","","RCMC_signalMCpt15_EB_trkIso_log",true);
  
  // endcap


  displayTwofiles("randomConeHists_Endcap_v5.root","template_comb3Iso.root",
		  "hsumisomc","h_EE_comb3Iso_sig_pt15",
		  "Iso_{ECAL}+Iso_{HCAL}+Iso_{TRK} (GeV)","","RCMC_signalMCpt15_EE_comb3Iso");


  displayTwofiles("randomConeHists_Endcap_v5.root","template_ecalIso.root",
		  "hecisomc","h_EE_ecalIso_sig_pt15",
		  "Iso_{ECAL} (GeV)","","RCMC_signalMCpt15_EE_ecalIso");
  

  displayTwofiles("randomConeHists_Endcap_v5.root","template_hcalIso.root",
		  "hhcisomc","h_EE_hcalIso_sig_pt15",
		  "Iso_{HCAL} (GeV)","","RCMC_signalMCpt15_EE_hcalIso");
  

  displayTwofiles("randomConeHists_Endcap_v5.root","template_hcalIso.root",
		  "hhcisomc","h_EE_hcalIso_sig_pt15",
		  "Iso_{HCAL} (GeV)","","RCMC_signalMCpt15_EE_hcalIso_log",
		  true);
 

  displayTwofiles("randomConeHists_Endcap_v5.root","template_trkIso.root",
		  "htkisomc","h_EE_trkIso_sig_pt15",
		  "Iso_{TRK} (GeV)","","RCMC_signalMCpt15_EE_trkIso");
  

  displayTwofiles("randomConeHists_Endcap_v5.root","template_trkIso.root",
		  "htkisomc","h_EE_trkIso_sig_pt15",
		  "Iso_{TRK} (GeV)","","RCMC_signalMCpt15_EE_trkIso_log",true);
  


  
}
