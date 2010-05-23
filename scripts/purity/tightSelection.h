#include <TCut.h>


TCut loosePhoton356 = "pt > 5 && "
  "( (isEB==1 && isEE==0 && "
  " ecalRecHitSumEtConeDR04 < (2.5+0.004*pt) && "
  " hcalTowerSumEtConeDR04  < 3.0 && "
  " trkSumPtHollowConeDR04  < 9.0 && "
  " hadronicOverEm < 0.15) || "
  "  (isEB==0 && isEE==1 && "
  " ecalRecHitSumEtConeDR04 < (3.0+0.0021*pt) && "
  " hcalTowerSumEtConeDR04  < 3.5 && "
  " trkSumPtHollowConeDR04  < 9.0 && "
  " hadronicOverEm < 0.15) ) ";


TCut removeSpikeCut = "((isEB && (seedSeverity!=3 && seedSeverity!=4 ) && (seedRecoFlag != 2) ) || isEE)";

TCut eventCut     = "( !TTBit[36] && !TTBit[37] && !TTBit[38] && !TTBit[39] && !vtxIsFake && vtxNdof > 4 && abs(vtxZ) <= 15) && HLT_Photon15_L1R";

TCut rsCut = "( pt > 15.0 && abs(eta)<1.45 && hadronicOverEm <  0.05 && sigmaIetaIeta < 0.01 ) || "
  "(pt > 15.0 && abs(eta) > 1.7 && abs(eta) < 2.5 && hadronicOverEm <  0.05 && abs(ESRatio) > 0.1)";


TCut RSCut_EB = "pt > 15.0 && abs(eta)<1.45 &&"
  " hadronicOverEm <  0.05 && "
  " sigmaIetaIeta < 0.01";

TCut RSCut_EE = "pt > 15.0 && abs(eta) > 1.7 && abs(eta) < 2.5 && "
  //  "(ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04) < 1.45 &&"
  " hadronicOverEm <  0.05 && "
  " abs(ESRatio) > 0.1";


TCut hardScatterCut = "isGenMatched && genMomId==22 && genCalIsoDR04 < 5.0"; 
TCut fragCut = "isGenMatched && abs(genMomId)<22 && genNSiblings>2 && genCalIsoDR04 < 5.0";
TCut decayCut = "(isGenMatched && abs(genMomId)>50) || !isGenMatched"; 

TCut IDCut = rsCut;

TCut sigCut = "isGenMatched && abs(genMomId) <= 22 && genCalIsoDR04 < 10.0"; 
TCut bkgCut = !sigCut;

TCut sigCut1 = "isGenMatched && abs(genMomId)==22 && genCalIsoDR04 < 10.0";
TCut sigCut2 = "isGenMatched && abs(genMomId)<= 22 && genCalIsoDR04 < 10.0";


/*
TCut sigCut1 = hardScatterCut;
TCut sigCut2 = "isGenMatched && abs(genMomId)<= 22";
TCut sigCut3 = "isGenMatched && genMomId==22 && genCalIsoDR04 < 5.0";
TCut sigCut4 = "isGenMatched && abs(genMomId)<=22 && genCalIsoDR04 < 5.0";
TCut sigCut5 = "isGenMatched && genCalIsoDR04 < 5.0";
*/

TCut simpleCut = "hadronicOverEm < 0.15&& abs(eta)<2.5 && pt>17"; 

TCut etaCut1 = "isEB && !isEE";

TCut etaCut2 = "isEE && !isEB";




