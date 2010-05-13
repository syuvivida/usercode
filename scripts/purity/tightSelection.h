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

TCut allRSCut = "(isEB && !isEE && ecalRecHitSumEtConeDR04/et < 0.095 && "
  "hcalTowerSumEtConeDR04/et < 0.045 && "
  "trkSumPtHollowConeDR04/et < 0.045 && "
  "hadronicOverEm < 0.045 && sigmaIetaIeta < 0.01) || "
  "(isEE && !isEB && ecalRecHitSumEtConeDR04/et < 0.07  && "
  "hcalTowerSumEtConeDR04/et < 0.025 && "
  "trkSumPtHollowConeDR04/et < 0.025 && "
  "hadronicOverEm < 0.015 && (ESRatio > 0.36 || ESRatio==-1))";


TCut rsCut = "(isEB && !isEE && isLoose && "
  "hadronicOverEm < 0.045 && sigmaIetaIeta < 0.01) || "
  "(isEE && !isEB && isLoose && "
  "hadronicOverEm < 0.015 && (ESRatio > 0.36 || ESRatio==-1))";

TCut hardScatterCut = "isGenMatched && genMomId==22"; 
TCut fragCut = "isGenMatched==1 && abs(genMomId)<22 && genNSiblings>2";
TCut decayCut = "(isGenMatched && abs(genMomId)>50) || !isGenMatched"; 

TCut basicCut = "isEB &&pt>15 && pt < 30 && !isEBGap && !isEEGap && !isEBEEGap" + rsCut;

//TCut sigCut = "isGenMatched && genCalIsoDR04 < 5.0";  
TCut sigCut = "isGenMatched && abs(genMomId) <= 22 && genCalIsoDR04 < 5.0"; 
TCut bkgCut = !sigCut;

TCut simpleCut = "hadronicOverEm < 0.15&& abs(eta)<2.5 && pt>17"; 

TCut etaCut1 = "isEB && !isEE";

TCut etaCut2 = "isEE && !isEB";




