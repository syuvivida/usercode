#ifndef TRIGFORMAT_HH
#define TRIGFORMAT_HH

namespace TRIGGER{
  enum {    
    HLT_L1SingleEG5                                = 0x1 << 1,                                               
    HLT_L1SingleEG8                                = 0x1 << 2,                                               
    HLT_Photon10_L1R                               = 0x1 << 3,                                              
    HLT_Photon10_LooseEcalIso_TrackIso_L1R         = 0x1 << 4,                                               
    HLT_Photon15_L1R                               = 0x1 << 5,            
    
    HLT_Photon15_TrackIso_L1R                      = 0x1 << 6, 

    HLT_Photon15_LooseEcalIso_L1R                  = 0x1 << 7,

    HLT_Photon20_LooseEcalIso_TrackIso_L1R         = 0x1 << 8, 

    HLT_Photon25_L1R                               = 0x1 << 9, 

    HLT_Mu5                                        = 0x1 << 10                   
  };

}

#endif
