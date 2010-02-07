#include <TCut.h>

TCut basicCut    = "vtxIsFake==0 && vtxNdof > 4"
  " && ((PHOLEAD_isEB==1 && PHOLEAD_isEE==0 && (PHOLEAD_eMax/PHOLEAD_e3x3<0.9)) || (PHOLEAD_isEB==0 && PHOLEAD_isEE==1))"
  " && nHfTowersP>=1 && nHfTowersN>=1 && clusVtxCut==1"
  " && PHOLEAD_pt > 2.0";

TCut fiducialCut = "PHOLEAD_isEBEEGap==0 && PHOLEAD_isEBGap==0"
  " && PHOLEAD_isEEGap==0 && PHOLEAD_isTransGap==0"
  " && abs(PHOLEAD_scEta) < 2.4"
  ;

TCut myCut       = basicCut + fiducialCut;

TCut looseCut    = basicCut;

TCut Cut_900GeV  = "run<124120";
TCut Cut_2360GeV = "run==124120"; 

TCut barrelCut   = "PHOLEAD_isEB==1 && PHOLEAD_isEE==0";
TCut endcapCut   = "PHOLEAD_isEB==0 && PHOLEAD_isEE==1";

TCut runCut = "run == 123596 ||"
  "run == 123615 ||"
  "run == 123732 ||"
  "run == 123815 ||"
  "run == 123818 ||"
  "run == 123906 ||"
  "run == 123908 ||"
  "run == 123909 ||"
  "run == 124009 ||"
  "run == 124020 ||"
  "run == 124022 ||"
  "run == 124023 ||"
  "run == 124024 ||"
  "run == 124025 ||"
  "run == 124027 ||"
  "run == 124030 ||"
  "run == 124120 ||"
  "run == 124230";

