#include <TCut.h>

TCut removeSpikeCut = "((isEB && (seedSeverity!=3 && seedSeverity!=4 ) && (seedRecoFlag != 2) && sigmaIetaIeta > 0.002 ) || isEE)";

TCut eventCut     = "( !TTBit[36] && !TTBit[37] && !TTBit[38] && !TTBit[39] && !vtxIsFake && vtxNdof > 4 && abs(vtxZ) <= 15) && HLT_Photon15_L1R";

TCut rsCutEB = "pt > 15.0 && abs(eta)<1.45 && hadronicOverEm <  0.05 && sigmaIetaIeta < 0.01 && (ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04) < 11.0"; 
TCut rsCutEE = "pt > 15.0 && abs(eta) > 1.7 && abs(eta) < 2.5 && hadronicOverEm <  0.05 && abs(ESRatio) > 0.1 && (ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04) < 11.0"; 

/*
TCut rsCutEB = "pt > 15.0 && abs(eta)<1.45 && hadronicOverEm <  0.05 && sigmaIetaIeta < 0.01 && (ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04) < 5.0";
TCut rsCutEE = "pt > 15.0 && abs(eta) > 1.7 && abs(eta) < 2.5 && hadronicOverEm <  0.05 && abs(ESRatio) > 0.1 && (ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04) < 5.0";
*/

TCut spikeCut = !removeSpikeCut + "isEB";


TCut rsCutSidebandEB = "pt > 15.0 && abs(eta)<1.45 && hadronicOverEm <  0.05 && sigmaIetaIeta > 0.011 &&sigmaIetaIeta<0.012 && (ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04) < 11.0";

TCut rsCutSidebandEE = "pt > 15.0 && abs(eta) > 1.7 && abs(eta) < 2.5 && hadronicOverEm <  0.05 && sigmaIetaIeta > 0.05 &&sigmaIetaIeta<0.06 && (ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04) < 11.0";

TCut hardScatterCut = "isGenMatched && genMomId==22 && genCalIsoDR04 < 5.0"; 
TCut fragCut = "isGenMatched && abs(genMomId)<22 && genNSiblings>2 && genCalIsoDR04 < 5.0";
TCut decayCut = "(isGenMatched && abs(genMomId)>50) || !isGenMatched"; 

TCut sigCut = "isGenMatched && abs(genMomId) <= 22 && genCalIsoDR04 < 5.0"; 
TCut bkgCut = !sigCut;





