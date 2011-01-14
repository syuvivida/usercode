#include <TCut.h>


TCut eventCut     = "(nVtxNotFake>0 && vtxNdof[0] > 4 && abs(vtxZ[0]) <= 24 && sqrt(vtxX[0]*vtxX[0]+vtxY[0]*vtxY[0]) <= 2.0 )"; 


/*
TCut trigCut = 
  "(pt >=25 && pt < 35 && HLT_Photon20_Cleaned_L1R) || "
  "(pt >=35 && pt < 55 && (HLT_Photon20_Cleaned_L1R || HLT_Photon30_Cleaned_L1R)) ||"
  "(pt >=55 && pt < 80 && (HLT_Photon20_Cleaned_L1R || HLT_Photon30_Cleaned_L1R || HLT_Photon50_Cleaned_L1R_v1)) ||"
  "(pt >=80 && (HLT_Photon20_Cleaned_L1R || HLT_Photon30_Cleaned_L1R || HLT_Photon50_Cleaned_L1R_v1 || HLT_Photon70_Cleaned_L1R_v1))";
*/


TCut trig20 = "(HLT_Photon20_Cleaned_L1R  && run>=138564 && run<=143962)";
TCut trig30 = "(HLT_Photon30_Cleaned_L1R && run>=144010 && run<=147116)";
TCut trig50 = "(HLT_Photon50_Cleaned_L1R_v1 && run>=147196 && run<=148058)";
TCut trig70 = "(HLT_Photon70_Cleaned_L1R_v1 && run>=148822 && run<=149294)";

TCut pt2535 = "(pt >=25 && pt < 35)";
TCut pt3555 = "(pt >=35 && pt < 55)";
TCut pt5580 = "(pt >=55 && pt < 80)";
TCut pt80   = "(pt >=80)";

TCut trigCut = 
  ( pt2535 &&  trig20)                      || 
  ( pt3555 && (trig20 || trig30))           ||
  ( pt5580 && (trig20 || trig30 || trig50)) ||
  ( pt80   && (trig20 || trig30 || trig50 || trig70));


TCut removeSpikeCut = "((abs(scEta)<1.4442 && seedSeverity!=3 && sigmaIetaIeta > 0.001 && sigmaIphiIphi >0.001 ) || (abs(scEta) > 1.566 && abs(scEta) < 2.5))";

TCut rsCutEB = "pt > 25.0 && abs(scEta)<1.4442 && hadronicOverEm <  0.05 && !hasPixelSeed && sigmaIetaIeta < 0.01 && (ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04) < 20.0 ";

TCut rsCutEE = "pt > 25.0 && abs(scEta) > 1.566 && abs(scEta) < 2.5 && hadronicOverEm < 0.05 && !hasPixelSeed && sigmaIetaIeta < 0.028 && (ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04) < 20.0";

/* TCut rsCutEB = "pt > 25.0 && abs(scEta)<1.4442 && hadronicOverEm <  0.05 && !hasPixelSeed";  */

/* TCut rsCutEE = "pt > 25.0 && abs(scEta) > 1.566 && abs(scEta) < 2.5 && hadronicOverEm < 0.05 && !hasPixelSeed";  */


TCut isolationCut2 = "(ecalRecHitSumEtConeDR04 < 4.2 + 0.003 * pt) && (hcalTowerSumEtConeDR04 <  2.2 + 0.001*pt) && (trkSumPtHollowConeDR04 < 2.0 + 0.001 * pt)";

TCut spikeCut = !removeSpikeCut + "isEB";
 
TCut rsCutSidebandEB = "pt > 25.0 && abs(scEta)<1.4442 && hadronicOverEm <  0.05 && !hasPixelSeed && sigmaIetaIeta > 0.011 && sigmaIetaIeta < 0.0115 && (ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04) < 20.0 ";

TCut rsCutSidebandEE = "pt > 25.0 && abs(scEta) > 1.566 && abs(scEta) < 2.5 && hadronicOverEm < 0.05 && !hasPixelSeed && sigmaIetaIeta > 0.038 && (ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04) < 20.0";

TCut hardScatterCut = "isGenMatched && genMomId==22 && genIsoDR04 < 5.0"; 
TCut fragCut = "isGenMatched && abs(genMomId)<22 && genNSiblings>2 && genCalIsoDR04 < 5.0";
TCut decayCut = "(isGenMatched && abs(genMomId)>50) || !isGenMatched"; 

TCut sigCut = "isGenMatched && abs(genMomId) <= 22 && genIsoDR04 < 5.0"; 
TCut bkgCut = !sigCut;





