#ifndef signalShapeReWeighting_h
#define signalShapeReWeighting_h

/*
   \author Matthias Mozer, modified by Shin-Shan Yu
*/

#include "TH1.h"
#include "TFile.h"
#include <string>
#include <vector>
#include <TROOT.h>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

class signalShapeReWeighting {
 public:
    
  signalShapeReWeighting(std::string inputfile); 
  virtual ~signalShapeReWeighting();
  double weight( double genHiggsMass, int mode); 
  // 0: central, -1: minus, +1: plus

 protected:
  std::vector<double> bincenters_;
  std::vector<double> weight_;
  std::vector<double> weightP_;
  std::vector<double> weightM_;
  
  bool FILE_NOT_EXIST;

};


#endif

