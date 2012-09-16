/**
   \author Matthias Mozer, modified by Shin-Shan Yu
  
*/
#ifndef signalShapeReWeighting_cxx
#define signalShapeReWeighting_cxx
#include "TH1.h"
#include "TFile.h"
#include <string>
#include "signalShapeReWeighting.h"
#include <algorithm>


signalShapeReWeighting::signalShapeReWeighting(std::string inputFile):
  FILE_NOT_EXIST(false)
{
  std::ifstream ifs(inputFile.data());
  double bincenter, initial, pow, powp, powm, out, outp, outm;

  if(!ifs){
    cout << "Error opening file " << inputFile.data() << endl;
    FILE_NOT_EXIST=true;
  }
  else
    cout << "Reading input file " << inputFile.data() << endl;

  while( ifs.good() ) {
    ifs >> bincenter >> initial >> pow >> powp >>  powm >> out>> outp >> outm;
    
    bincenters_.push_back(bincenter);
    if(initial > 0){
      weight_.push_back(   TMath::Max(0.,out/initial) );
      weightP_.push_back(  TMath::Max(0.,outp/initial) );
      weightM_.push_back(  TMath::Max(0.,outm/initial) );
    }else{
      //weights are not defined if initial distribution is 0 => set weight to 0
      weight_.push_back( 0. );
      weightP_.push_back( 0. );
      weightM_.push_back( 0. );
    }

  }


}

signalShapeReWeighting::~signalShapeReWeighting()
{
}



double signalShapeReWeighting::weight(double m, int mode ) {

  if(FILE_NOT_EXIST)return 1.0;
  
  if( m < bincenters_.front() || m >  bincenters_.back() ){ 
    // set weights to 0 if out of range
    return 0.0;
  }

  std::vector<double>::iterator low;
  low=lower_bound( bincenters_.begin(), bincenters_.end(),m ); 
  int lowindex=(low-  bincenters_.begin());

  double weight_lo   = 0.0;
  double weight_hi   = 0.0;
  double weight_this = 0.0;

  switch (mode)
    {
    case 0: 
      weight_lo = weight_[lowindex-1];
      weight_hi = weight_[lowindex];
      break;
    case 1: 
      weight_lo = weightP_[lowindex-1];
      weight_hi = weightP_[lowindex];
      break;
    case -1: 
      weight_lo = weightM_[lowindex-1];
      weight_hi = weightM_[lowindex];
      break;
    default:
      weight_lo = 0.0;
      weight_hi = 0.0;
      break;
    } // end of switch                    

  if(m == *low )weight_this = weight_hi;
  else //linear interpolation
    weight_this = weight_lo + (m -bincenters_[lowindex-1] ) * 
      (weight_hi - weight_lo)/(bincenters_[lowindex]-bincenters_[lowindex-1]);
  
  return weight_this;
}


#endif
