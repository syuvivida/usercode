
void call_scanGenIso(std::string prefix)
{
  gROOT->ProcessLine(".L scanGenIso.C++");
  for(int i=1;i<=10;i++)
    scanGenIso(prefix+"_sumIso","ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04",24,-1,11.0,i);

}
