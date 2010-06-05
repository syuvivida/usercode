#include "/afs/cern.ch/user/s/syu/scripts/displayTwofiles.C"

void compareRCDataMC()
{
  // Barrel first
  displayTwofiles("randomConeHists_Barrel.root","randomConeHists_Barrel.root",
		  "hsumiso","hsumisomc",
		  "Iso_{ECAL}+Iso_{HCAL}+Iso_{TRK} (GeV)","","RC_EB_comb3Iso");


  displayTwofiles("randomConeHists_Barrel.root","randomConeHists_Barrel.root",
		  "heciso","hecisomc",
		  "Iso_{ECAL} (GeV)","","RC_EB_ecalIso");
  

  displayTwofiles("randomConeHists_Barrel.root","randomConeHists_Barrel.root",
		  "hhciso","hhcisomc",
		  "Iso_{HCAL} (GeV)","","RC_EB_hcalIso");
  

  displayTwofiles("randomConeHists_Barrel.root","randomConeHists_Barrel.root",
		  "hhciso","hhcisomc",
		  "Iso_{HCAL} (GeV)","","RC_EB_hcalIso_log",true);
 

  displayTwofiles("randomConeHists_Barrel.root","randomConeHists_Barrel.root",
		  "htkiso","htkisomc",
		  "Iso_{TRK} (GeV)","","RC_EB_trkIso");
  

  displayTwofiles("randomConeHists_Barrel.root","randomConeHists_Barrel.root",
		  "htkiso","htkisomc",
		  "Iso_{TRK} (GeV)","","RC_EB_trkIso_log",true);
 
  c1->SetLogy(1);



  // endcap second

  displayTwofiles("randomConeHists_Endcap.root","randomConeHists_Endcap.root",
		  "hsumiso","hsumisomc",
		  "Iso_{ECAL}+Iso_{HCAL}+Iso_{TRK} (GeV)","","RC_EE_comb3Iso");


  displayTwofiles("randomConeHists_Endcap.root","randomConeHists_Endcap.root",
		  "heciso","hecisomc",
		  "Iso_{ECAL} (GeV)","","RC_EE_ecalIso");
  

  displayTwofiles("randomConeHists_Endcap.root","randomConeHists_Endcap.root",
		  "hhciso","hhcisomc",
		  "Iso_{HCAL} (GeV)","","RC_EE_hcalIso");
  

  displayTwofiles("randomConeHists_Endcap.root","randomConeHists_Endcap.root",
		  "hhciso","hhcisomc",
		  "Iso_{HCAL} (GeV)","","RC_EE_hcalIso_log",true);
 

  displayTwofiles("randomConeHists_Endcap.root","randomConeHists_Endcap.root",
		  "htkiso","htkisomc",
		  "Iso_{TRK} (GeV)","","RC_EE_trkIso");
  

  displayTwofiles("randomConeHists_Endcap.root","randomConeHists_Endcap.root",
		  "htkiso","htkisomc",
		  "Iso_{TRK} (GeV)","","RC_EE_trkIso_log",true);
 
  c1->SetLogy(1);

  
}
