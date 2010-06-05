/*==============================================================================================
  Produce histograms for combined isolation
    
  string prefix = "FileName";// File name for output
  
  root -q -b call_allHistos.C\(\"test\"\)
  Within call_allHistos.C, the main file is allHistos.C
 
  allHistos(std::string outputName="",                                                                               
  std::string var="(ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04)",                        
  int nbin=20,                                                                                                       
  double xmin=-1.0, double xmax=11.0, double split = 1.0,bool normalize=false)                                                             

  ==============================================================================================*/


void call_allHistos(std::string prefix)
{
  gROOT->ProcessLine(".L allHistos.C++");
  // if you are using split MC for producing templates
//   allHistos(prefix+"_sumIso","ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04",20,-1,11.0,2.0);
  // if you are using total MC for producing templates
  allHistos(prefix+"_sumIso","ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04",20,-1,11.0,1.0);

}
