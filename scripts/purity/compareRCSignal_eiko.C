#include "/afs/cern.ch/user/s/syu/scripts/displayTwofiles.C"

void compareRCSignal_eiko()
{
  // barrel
  displayTwofiles("randomConeHists_Barrel_v5.root","TightEESB_proj_comb.root",
		  "hsumisomc","h_EB_comb3Iso_PhoJet30_SIG",
		  "Iso_{ECAL}+Iso_{HCAL}+Iso_{TRK} (GeV)","","RCMC_PhoJetPt30_EB_comb3Iso");


  displayTwofiles("randomConeHists_Barrel_v5.root","TightEESB_proj_comb.root",
		  "hecisomc","h_EB_ecalIso_PhoJet30_SIG",
		  "Iso_{ECAL} (GeV)","","RCMC_PhoJetPt30_EB_ecalIso");
  

  displayTwofiles("randomConeHists_Barrel_v5.root","TightEESB_proj_comb.root",
		  "hhcisomc","h_EB_hcalIso_PhoJet30_SIG",
		  "Iso_{HCAL} (GeV)","","RCMC_PhoJetPt30_EB_hcalIso");
  

  displayTwofiles("randomConeHists_Barrel_v5.root","TightEESB_proj_comb.root",
		  "hhcisomc","h_EB_hcalIso_PhoJet30_SIG",
		  "Iso_{HCAL} (GeV)","","RCMC_PhoJetPt30_EB_hcalIso_log",
		  true);
 

  displayTwofiles("randomConeHists_Barrel_v5.root","TightEESB_proj_comb.root",
		  "htkisomc","h_EB_trkIso_PhoJet30_SIG",
		  "Iso_{TRK} (GeV)","","RCMC_PhoJetPt30_EB_trkIso");
  

  displayTwofiles("randomConeHists_Barrel_v5.root","TightEESB_proj_comb.root",
		  "htkisomc","h_EB_trkIso_PhoJet30_SIG",
		  "Iso_{TRK} (GeV)","","RCMC_PhoJetPt30_EB_trkIso_log",true);
  
  // endcap


  displayTwofiles("randomConeHists_Endcap_v5.root","TightEESB_proj_comb.root",
		  "hsumisomc","h_EE_comb3Iso_PhoJet30_SIG",
		  "Iso_{ECAL}+Iso_{HCAL}+Iso_{TRK} (GeV)","","RCMC_PhoJetPt30_EE_comb3Iso");


  displayTwofiles("randomConeHists_Endcap_v5.root","TightEESB_proj_comb.root",
		  "hecisomc","h_EE_ecalIso_PhoJet30_SIG",
		  "Iso_{ECAL} (GeV)","","RCMC_PhoJetPt30_EE_ecalIso");
  

  displayTwofiles("randomConeHists_Endcap_v5.root","TightEESB_proj_comb.root",
		  "hhcisomc","h_EE_hcalIso_PhoJet30_SIG",
		  "Iso_{HCAL} (GeV)","","RCMC_PhoJetPt30_EE_hcalIso");
  

  displayTwofiles("randomConeHists_Endcap_v5.root","TightEESB_proj_comb.root",
		  "hhcisomc","h_EE_hcalIso_PhoJet30_SIG",
		  "Iso_{HCAL} (GeV)","","RCMC_PhoJetPt30_EE_hcalIso_log",
		  true);
 

  displayTwofiles("randomConeHists_Endcap_v5.root","TightEESB_proj_comb.root",
		  "htkisomc","h_EE_trkIso_PhoJet30_SIG",
		  "Iso_{TRK} (GeV)","","RCMC_PhoJetPt30_EE_trkIso");
  

  displayTwofiles("randomConeHists_Endcap_v5.root","TightEESB_proj_comb.root",
		  "htkisomc","h_EE_trkIso_PhoJet30_SIG",
		  "Iso_{TRK} (GeV)","","RCMC_PhoJetPt30_EE_trkIso_log",true);
  


  
}
