#include <TCut.h>

TCut myCut = "vtxIsFake==0 && PHOLEAD_isEBEEGap==0 && PHOLEAD_isEBGap==0"
  " && PHOLEAD_isEEGap==0 && PHOLEAD_isTransGap==0"
  " && !(abs(PHOLEAD_eMax/PHOLEAD_e3x3-1)<0.01)"
  " && nHfTowersP>=1 && nHfTowersN>=1 && clusVtxCut==1";

TCut Cut_900GeV = "run!=124120";
TCut Cut_2360GeV = "run==124120";

TCut barrelCut = "PHOLEAD_isEB==1 && PHOLEAD_isEE==0";
TCut endcapCut = "PHOLEAD_isEB==0 && PHOLEAD_isEE==1";

