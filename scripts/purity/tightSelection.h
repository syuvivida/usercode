#include <TCut.h>


TCut basicCut    = "phoEt > 20 && phoEt < 80"
  " && phoPos==0"
  " && phoTrkIsoHollowDR04 < 3.4"
  " && phoEcalIsoDR04 < (3+0.0013*phoEt) "
  " && phoHcalIsoDR04 < 3.0"
  " && phoHoverE < 0.054"
  //  " && phoSigmaIEtaIEta < 0.01"
  ;


/*
TCut basicCut    = "phoEt > 20 && phoEt < 80"                                 
  " && phoPos==1" 
  " && abs(phoEta)<2.5"
  " && phoTrkIsoHollowDR03 < 1.5+0.02*phoEt"
  " && phoEcalIsoDR03 < 1.0+0.05*phoEt "
  " && phoHcalIsoDR03 < 0.5+0.04*phoEt"
  " && phoHoverE < 0.15"                                                       
  ;                                                                             
*/

/*
TCut basicCut    = "phoEt > 20 && phoEt < 30"                                 
  " && phoPos==0"                                                             
  " && phoR9 > 0.93"
  " && phoTrkIsoHollowDR03 < 9.0"
  " && phoEcalIsoDR03 < (5.0+0.04*phoEt) "                                  
  " && phoHcalIsoDR03 < 5.0"                                       
  " && phoHoverE < 0.15"                                        
;                                                                             
*/

TCut realCut = "phoGenIndex>=0 && phoGenMomPID==22";

/*

x-sections

photon+jets 

20-30 5.718e+4 110000
30-50 1.652e+4 104060
50-80 2.723e+3 104060

QCD

20-30 0.2355e+9 0.0073 33305352
30-80 0.0593e+9 0.059  40899325


*/
