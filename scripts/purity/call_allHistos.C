
void call_allHistos(std::string prefix)
{
  gROOT->ProcessLine(".L allHistos.C++");
//   allHistos(prefix+"_ecalIso","ecalRecHitSumEtConeDR04",40,-1,19.0);
//   allHistos(prefix+"_hcalIso","hcalTowerSumEtConeDR04",40,-1,19.0);
//   allHistos(prefix+"_trackIso","trkSumPtHollowConeDR04",40,-1,19.0);
  allHistos(prefix+"_sumIso","ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04",20,-1,9.0);

}
