#ifndef CHI2NBINS_H
#define CHI2NBINS_H

#include <iostream>
#include <TF1.h>
#include <TH1.h>

using namespace std;

void chi2Nbins( TF1 *func, const TH1F *hist, double& chi2, int& nbins ){

  printf( "Chi2 Calculation: \n");

  chi2 = 0;
  nbins = 0;

  int nBins = 0;
  int evtThresh = -1;

  double theoryErr, theory, lowEdge, highEdge, binChi2;

  double  binWidth    = hist -> GetXaxis() -> GetBinWidth(1);

  for ( int i = 1; i <= hist -> GetNbinsX(); i++){
 
    double nEvts = hist -> GetBinContent(i);

    if ( nEvts < evtThresh) continue;

    lowEdge   = hist -> GetXaxis() -> GetBinLowEdge( i );
    highEdge  = hist -> GetXaxis() -> GetBinUpEdge( i );
    theory    = func -> Integral( lowEdge, highEdge ) / binWidth;
    theoryErr    = sqrt(theory);

    binChi2 = (nEvts - theory) / theoryErr;
    binChi2 *= binChi2;
    
    chi2 += binChi2;

    printf( "%d) [%d, %d] [%f, %f] data: %f theo: %f chi2: %f total: %f \n ", 
	    nBins, i, i, lowEdge, highEdge,nEvts, theory, binChi2, chi2);
    
    nBins++;

  }
  

  int npar = func -> GetNpar();

  double parmin, parmax;
  int fixed = 0;

  for ( int i = 0; i < npar; i++){
    func -> GetParLimits(i , parmin, parmax );
    if (parmin > parmax)
      fixed++;
  }

  printf("Function has %d parameters, %d are fixed, really: %d pars.\n", npar, fixed, npar - fixed);

  int NDF = nBins;

  printf("Fit chi2/NDF = %f/%d, prob: %f\n", chi2, NDF, TMath::Prob(chi2,NDF)*100);
  nbins = NDF;
}

void chi2NbinsHisto( const TH1F* htheory, const TH1F *hist, double& chi2, int& nbins ){

  printf( "Chi2 Calculation: \n");

  chi2 = 0;
  nbins = 0;

  int nBins = 0;
  int evtThresh = -1;

  double theoryErr, theory, lowEdge, highEdge, binChi2;

  double  binWidth    = hist -> GetXaxis() -> GetBinWidth(1);

  for ( int i = 1; i <= hist -> GetNbinsX(); i++){
 
    double nEvts = hist -> GetBinContent(i);

    if ( nEvts < evtThresh) continue;
    lowEdge   = hist -> GetXaxis() -> GetBinLowEdge( i );
    highEdge  = hist -> GetXaxis() -> GetBinUpEdge( i );

    theory       = htheory->GetBinContent(i);
    theoryErr    = sqrt(theory);
    
    binChi2 = (nEvts - theory) / theoryErr;
    binChi2 *= binChi2;
    
    chi2 += binChi2;

    printf( "%d) [%d, %d] [%f, %f] data: %f theo: %f chi2: %f total: %f \n ", 
	    nBins, i, i, lowEdge, highEdge,nEvts, theory, binChi2, chi2);
    
    nBins++;

  }
  

  int NDF = nBins;

  printf("Fit chi2/NDF = %f/%d, prob: %f\n", chi2, NDF, TMath::Prob(chi2,NDF)*100);
  nbins = NDF;
}




#endif
