#include "/afs/cern.ch/user/s/syu/scripts/displayTwofiles.C"

void compareSpikeSignal()
{
  displayTwofiles("spike_withID.root","template_comb3Iso.root",
		  "h_EB_comb3Iso_EGdata_SIG","h_EB_comb3Iso_sig_pt20",
		  "Iso_{ECAL}+Iso_{HCAL}+Iso_{TRK} (GeV)","","spike_signalMCpt20_EB_comb3Iso");


  displayTwofiles("spike_withID.root","template_ecalIso.root",
		  "h_EB_ecalIso_EGdata_SIG","h_EB_ecalIso_sig_pt20",
		  "Iso_{ECAL} (GeV)","","spike_signalMCpt20_EB_ecalIso");
  

  displayTwofiles("spike_withID.root","template_hcalIso.root",
		  "h_EB_hcalIso_EGdata_SIG","h_EB_hcalIso_sig_pt20",
		  "Iso_{HCAL} (GeV)","","spike_signalMCpt20_EB_hcalIso");
  

  displayTwofiles("spike_withID.root","template_hcalIso.root",
		  "h_EB_hcalIso_EGdata_SIG","h_EB_hcalIso_sig_pt20",
		  "Iso_{HCAL} (GeV)","","spike_signalMCpt20_EB_hcalIso_log",
		  true);
 

  displayTwofiles("spike_withID.root","template_trkIso.root",
		  "h_EB_trkIso_EGdata_SIG","h_EB_trkIso_sig_pt20",
		  "Iso_{TRK} (GeV)","","spike_signalMCpt20_EB_trkIso");
  

  displayTwofiles("spike_withID.root","template_trkIso.root",
		  "h_EB_trkIso_EGdata_SIG","h_EB_trkIso_sig_pt20",
		  "Iso_{TRK} (GeV)","","spike_signalMCpt20_EB_trkIso_log",true);
  


  
}
