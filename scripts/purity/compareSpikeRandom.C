#include "/afs/cern.ch/user/s/syu/scripts/displayTwofiles.C"

void compareSpikeRandom()
{
  displayTwofiles("spike_withID.root","randomConeHists_Barrel.root",
		  "h_EB_comb3Iso_EGdata_SIG","hsumiso",
		  "Iso_{ECAL}+Iso_{HCAL}+Iso_{TRK} (GeV)","","spike_random_EB_comb3Iso");


  displayTwofiles("spike_withID.root","randomConeHists_Barrel.root",
		  "h_EB_ecalIso_EGdata_SIG","heciso",
		  "Iso_{ECAL} (GeV)","","spike_random_EB_ecalIso");
  

  displayTwofiles("spike_withID.root","randomConeHists_Barrel.root",
		  "h_EB_hcalIso_EGdata_SIG","hhciso",
		  "Iso_{HCAL} (GeV)","","spike_random_EB_hcalIso");
  

  displayTwofiles("spike_withID.root","randomConeHists_Barrel.root",
		  "h_EB_hcalIso_EGdata_SIG","hhciso",
		  "Iso_{HCAL} (GeV)","","spike_random_EB_hcalIso_log",true);
 

  displayTwofiles("spike_withID.root","randomConeHists_Barrel.root",
		  "h_EB_trkIso_EGdata_SIG","htkiso",
		  "Iso_{TRK} (GeV)","","spike_random_EB_trkIso");
  

  displayTwofiles("spike_withID.root","randomConeHists_Barrel.root",
		  "h_EB_trkIso_EGdata_SIG","htkiso",
		  "Iso_{TRK} (GeV)","","spike_random_EB_trkIso_log",true);
 
  c1->SetLogy(1);

  
}
