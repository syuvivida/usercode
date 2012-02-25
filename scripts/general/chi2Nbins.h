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
  int evtThresh = 25;

  double theoryErr, theory, lowEdge, highEdge, binChi2;

  double  binWidth    = hist -> GetXaxis() -> GetBinWidth(1);

  for ( int i = 1; i <= hist -> GetNbinsX(); i++){
 
    double nEvts = hist -> GetBinContent(i);
    cout << "bin " << i << ": " << nEvts << endl;
    if ( nEvts < evtThresh) continue;

    lowEdge   = hist -> GetXaxis() -> GetBinLowEdge( i );
    highEdge  = hist -> GetXaxis() -> GetBinUpEdge( i );
    theory    = func -> Integral( lowEdge, highEdge ) / binWidth;
    theoryErr    = sqrt(theory);

    binChi2 = (nEvts - theory) / theoryErr;
    binChi2 *= binChi2;
    
    chi2 += binChi2;

/*     printf( "%d) [%d, %d] [%f, %f] data: %f theo: %f chi2: %f total: %f \n ",  */
/* 	    nBins, i, i, lowEdge, highEdge,nEvts, theory, binChi2, chi2); */
    
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

    if(theoryErr < 1e-6)continue;
    
    binChi2 = (nEvts - theory) / theoryErr;
    binChi2 *= binChi2;
    
    chi2 += binChi2;

/*     printf( "%d) [%d, %d] [%f, %f] data: %f theo: %f theoErr: %f chi2: %f total: %f \n ",  */
/* 	    nBins, i, i, lowEdge, highEdge,nEvts, theory, theoryErr, binChi2, chi2); */
    
    nBins++;

  }
  

  int NDF = nBins;

  printf("Fit chi2/NDF = %f/%d, prob: %f\n", chi2, NDF, TMath::Prob(chi2,NDF)*100);
  nbins = NDF;
}

void chi2NbinsCompare( const TH1F* h1, const TH1F *h2, double& chi2, int& nbins , int binLo=-1, int binHi=-1){

  printf( "Chi2 Calculation: \n");

  chi2 = 0;
  nbins = 0;

  int nBins = 0;
  double evtThresh = 1e-6;

  double h1Err, h1Value, h2Err, h2Value, lowEdge, highEdge, binChi2;

  double  binWidth    = h1 -> GetXaxis() -> GetBinWidth(1);

  if(binLo < 0 || binHi < 0)
    {
      binLo = 1;
      binHi = h1->GetNbinsX();
    }

  for ( int i = binLo ; i <= binHi; i++){
 
    h1Value = h1 -> GetBinContent(i);
    h1Err   = h1 -> GetBinError(i);

    h2Value = h2 -> GetBinContent(i);
    h2Err   = h2 -> GetBinError(i);

    if ( h1Value < evtThresh && h2Value  < evtThresh) continue;
    if ( h1Err < 1e-6 && h2Err  < 1e-6) continue;
    lowEdge   = h1 -> GetXaxis() -> GetBinLowEdge( i );
    highEdge  = h1 -> GetXaxis() -> GetBinUpEdge( i );

    binChi2 = (h1Value - h2Value) / sqrt(h1Err*h1Err+h2Err*h2Err);
    binChi2 *= binChi2;
    
    chi2 += binChi2;

    printf( "%d) [%d, %d] [%f, %f] h1: %f h2: %f h1Err: %f h2Err: %f chi2: %f total: %f \n ",  
   	    nBins, i, i, lowEdge, highEdge,h1Value, h2Value, h1Err, h2Err, binChi2, chi2); 
    
    nBins++;

  }
  

  int NDF = nBins;

  printf("Fit chi2/NDF = %f/%d, prob: %f\n", chi2, NDF, TMath::Prob(chi2,NDF)*100);
  nbins = NDF;
}



#endif
