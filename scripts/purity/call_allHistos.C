
void call_allHistos(std::string prefix)
{
  gROOT->ProcessLine(".L allHistos.C++");
//   allHistos(prefix+"_ecalIso","ecalRecHitSumEtConeDR04",20,-1,11.0);
//   allHistos(prefix+"_hcalIso","hcalTowerSumEtConeDR04",20,-1,11.0);
//   allHistos(prefix+"_trackIso","trkSumPtHollowConeDR04",20,-1,11.0);
  allHistos(prefix+"_sumIso","ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04",20,-1,11.0);

}
