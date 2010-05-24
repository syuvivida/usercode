
void call_fragDirect(std::string prefix)
{
  gROOT->ProcessLine(".L fragDirect.C++");
  fragDirect(prefix+"_sumIso","ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04",20,-1,11.0,5);

}
