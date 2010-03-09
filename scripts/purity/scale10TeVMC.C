#include <iostream>
#include <string>
#include "/mc/Laurent/scripts/scaleHistos.C"

using namespace std;

const Int_t NMC          = 7;
const Float_t scaleInfo[2][NMC][3]=
  {
    // xsec in pb, filter efficiency, number of generated events
    // photon + jet signal
    {
      {288700, 1.0, 1073270},
      {32220,  1.0, 1006118},
      {1010,   1.0, 1003673},
      {51.43,  1.0, 1496900},
      {4.193,  1.0, 1048016},
      {0.4515, 1.0, 1014413},
      {0.02003, 1.0, 1241880}
    },

    // QCD background
    {
      {1458000000, 1.0, 5172530},
      {109000000, 1.0, 13643440},
      {1936000, 1.0, 2147800},
      {62510, 1.0, 2189700},
      {3669, 1.0, 2087000},
      {315.3, 1.0, 2131360},
      {11.94, 1.0, 2092640}
    }
  };


void scale10TeVMC(std::string inputFile_)
{

  // check if this is a QCD or a photon+jet samples
  // should check the particle list in GenpBlock when they are available
  Int_t isQCD = 0;
  if(inputFile_.find("QCD") != std::string::npos)isQCD=1;

  // check the pthat sample
  const Int_t pthat[NMC+1]={15, 30, 80, 170, 300, 470, 800, 10000000};
  Char_t name[300];
  Float_t PTHAT_MIN = -999.;
  Float_t PTHAT_MAX = -999.;
  Int_t thisMC = -1;
  for(Int_t i=0; i < NMC; i++)
    {
      sprintf(name,"%d",pthat[i]);
      if(inputFile_.find(name) != std::string::npos){
	PTHAT_MIN = pthat[i];
	PTHAT_MAX = pthat[i+1];
	thisMC = i;
      }
    }
 

  cout << "Restricted to pthat = " << PTHAT_MIN << " to " << PTHAT_MAX 
       << " GeV" << endl;

  Float_t scaleFactor = 1.0;
  scaleFactor = scaleInfo[isQCD][thisMC][0]*scaleInfo[isQCD][thisMC][1]/scaleInfo[isQCD][thisMC][2];

  cout << "For this MC, the x-section = " << 
    scaleInfo[isQCD][thisMC][0] << " pb, the filter efficiency = " << 
    scaleInfo[isQCD][thisMC][1] << ", the number of generated events = " << 
    scaleInfo[isQCD][thisMC][2] << ", and the scale factor = " << scaleFactor << endl;

  scaleHistos(inputFile_,scaleFactor);
  
}
