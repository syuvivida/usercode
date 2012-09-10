/**
   \class    standalone_LumiReWeighting standalone_LumiReWeighting.h "PhysicsTools/Utilities/interface/standalone_LumiReWeighting.h"
   \brief    Class to provide lumi weighting for analyzers to weight "flat-to-N" MC samples to data

   This class will trivially take two histograms:
   1. The generated "flat-to-N" distributions from a given processing (or any other generated input)
   2. A histogram generated from the "estimatePileup" macro here:

   https://twiki.cern.ch/twiki/bin/view/CMS/LumiCalc#How_to_use_script_estimatePileup

   and produce weights to convert the input distribution (1) to the latter (2).

   \author Salvatore Rappoccio, modified by Mike Hildreth
  
*/
#ifndef standalone_LumiReWeighting_cxx
#define standalone_LumiReWeighting_cxx
#include "TH1.h"
#include "TFile.h"
#include <string>
#include "standalone_LumiReWeighting.h"

//=======================================================
// For 2012 Data and MC
//=======================================================

Double_t Summer2012[60] = {
  2.344E-05,
  2.344E-05,
  2.344E-05,
  2.344E-05,
  4.687E-04,
  4.687E-04,
  7.032E-04,
  9.414E-04,
  1.234E-03,
  1.603E-03,
  2.464E-03,
  3.250E-03,
  5.021E-03,
  6.644E-03,
  8.502E-03,
  1.121E-02,
  1.518E-02,
  2.033E-02,
  2.608E-02,
  3.171E-02,
  3.667E-02,
  4.060E-02,
  4.338E-02,
  4.520E-02,
  4.641E-02,
  4.735E-02,
  4.816E-02,
  4.881E-02,
  4.917E-02,
  4.909E-02,
  4.842E-02,
  4.707E-02,
  4.501E-02,
  4.228E-02,
  3.896E-02,
  3.521E-02,
  3.118E-02,
  2.702E-02,
  2.287E-02,
  1.885E-02,
  1.508E-02,
  1.166E-02,
  8.673E-03,
  6.190E-03,
  4.222E-03,
  2.746E-03,
  1.698E-03,
  9.971E-04,
  5.549E-04,
  2.924E-04,
  1.457E-04,
  6.864E-05,
  3.054E-05,
  1.282E-05,
  5.081E-06,
  1.898E-06,
  6.688E-07,
  2.221E-07,
  6.947E-08,
  2.047E-08
};  

// this was the central value produced with buggy pileup JSON file
double Data2012Buggy[60]={
  8994.65,
  14398.2,
  97362.2,
  4.13664e+06,
  5.76045e+06,
  2.27729e+06,
  1.12598e+07,
  2.92359e+07,
  5.95919e+07,
  1.0363e+08,
  1.64679e+08,
  2.37665e+08,
  3.05316e+08,
  3.63753e+08,
  4.1495e+08,
  4.41045e+08,
  4.15773e+08,
  3.57062e+08,
  3.02194e+08,
  2.63202e+08,
  2.35093e+08,
  2.11978e+08,
  1.90889e+08,
  1.70897e+08,
  1.51584e+08,
  1.32755e+08,
  1.14521e+08,
  9.71723e+07,
  8.10282e+07,
  6.63574e+07,
  5.33436e+07,
  4.20779e+07,
  3.25602e+07,
  2.47114e+07,
  1.83915e+07,
  1.34213e+07,
  9.6024e+06,
  6.73474e+06,
  4.62982e+06,
  3.11927e+06,
  2.05931e+06,
  1.33199e+06,
  843963,
  523736,
  318272,
  189370,
  110302,
  62886.9,
  35089.9,
  19160.3,
  10237.2,
  5351.54,
  2736.95,
  1369.37,
  670.216,
  320.875,
  150.269,
  68.8335,
  30.8404,
  13.5152
};

double Data2012[60]={
  0.223516,
  25.8075,
  2351.55,
  108612,
  268471,
  1.90027e+06,
  1.04981e+07,
  2.71424e+07,
  5.54473e+07,
  9.6718e+07,
  1.51618e+08,
  2.13609e+08,
  2.65407e+08,
  3.00605e+08,
  3.26343e+08,
  3.45545e+08,
  3.4915e+08,
  3.3622e+08,
  3.14297e+08,
  2.90056e+08,
  2.68352e+08,
  2.52931e+08,
  2.44914e+08,
  2.40746e+08,
  2.32957e+08,
  2.14909e+08,
  1.85442e+08,
  1.47536e+08,
  1.06784e+08,
  6.9313e+07,
  3.98601e+07,
  2.01579e+07,
  8.93546e+06,
  3.46912e+06,
  1.18145e+06,
  354734,
  94898.1,
  23026.4,
  5197.28,
  1122.76,
  237.463,
  49.5914,
  10.1747,
  2.02972,
  0.390504,
  0.0721848,
  0.0128062,
  0.00217838,
  0.000354464,
  5.49583e-05,
  8.0806e-06,
  1.1215e-06,
  1.4629e-07,
  1.78541e-08,
  2.05299e-09,
  2.02801e-10,
  0,
  0,
  0,
  0
};


double Data2012Up[60]={
  0.148557,
  18.2664,
  1226.84,
  60640.9,
  255564,
  919867,
  6.60171e+06,
  1.91216e+07,
  4.0348e+07,
  7.27575e+07,
  1.16497e+08,
  1.71426e+08,
  2.252e+08,
  2.66543e+08,
  2.94373e+08,
  3.16521e+08,
  3.314e+08,
  3.31462e+08,
  3.18047e+08,
  2.97704e+08,
  2.75697e+08,
  2.55982e+08,
  2.41664e+08,
  2.33841e+08,
  2.29966e+08,
  2.24179e+08,
  2.10302e+08,
  1.86398e+08,
  1.54038e+08,
  1.17349e+08,
  8.14035e+07,
  5.08171e+07,
  2.83063e+07,
  1.40061e+07,
  6.14548e+06,
  2.39118e+06,
  826872,
  255536,
  71323.4,
  18284.3,
  4404.03,
  1021.55,
  232.599,
  52.3634,
  11.6066,
  2.51038,
  0.525855,
  0.106269,
  0.0206927,
  0.0038799,
  0.000699421,
  0.000120863,
  1.99424e-05,
  3.12883e-06,
  4.64965e-07,
  6.52845e-08,
  8.56987e-09,
  1.08471e-09,
  8.82189e-11,
  0
};

double Data2012Down[60]={
  0.344462,
  39.6011,
  4753.77,
  174365,
  343330,
  3.75385e+06,
  1.60453e+07,
  3.8747e+07,
  7.61792e+07,
  1.28897e+08,
  1.96074e+08,
  2.59746e+08,
  3.05442e+08,
  3.3612e+08,
  3.59999e+08,
  3.68222e+08,
  3.56389e+08,
  3.32807e+08,
  3.05993e+08,
  2.81983e+08,
  2.65346e+08,
  2.57148e+08,
  2.52424e+08,
  2.4178e+08,
  2.18367e+08,
  1.82255e+08,
  1.38358e+08,
  9.40039e+07,
  5.62515e+07,
  2.93008e+07,
  1.32075e+07,
  5.14151e+06,
  1.72973e+06,
  505109,
  129355,
  29601.8,
  6225.4,
  1243.6,
  242.433,
  46.5942,
  8.77441,
  1.59975,
  0.279923,
  0.0468395,
  0.00748728,
  0.00114162,
  0.00016547,
  2.26834e-05,
  2.92449e-06,
  3.52748e-07,
  3.96928e-08,
  4.11229e-09,
  3.77351e-10,
  2.04314e-11,
  0,
  0,
  0,
  0,
  0,
  0
};


//=======================================================
// For 2011 Data and MC
//=======================================================

Double_t Fall2011[50] = {
  0.003388501,
  0.010357558,
  0.024724258,
  0.042348605,
  0.058279812,
  0.068851751,
  0.072914824,
  0.071579609,
  0.066811668,
  0.060672356,
  0.054528356,
  0.04919354,
  0.044886042,
  0.041341896,
  0.0384679,
  0.035871463,
  0.03341952,
  0.030915649,
  0.028395374,
  0.025798107,
  0.023237445,
  0.020602754,
  0.0180688,
  0.015559693,
  0.013211063,
  0.010964293,
  0.008920993,
  0.007080504,
  0.005499239,
  0.004187022,
  0.003096474,
  0.002237361,
  0.001566428,
  0.001074149,
  0.000721755,
  0.000470838,
  0.00030268,
  0.000184665,
  0.000112883,
  6.74043E-05,
  3.82178E-05,
  2.22847E-05,
  1.20933E-05,
  6.96173E-06,
  3.4689E-06,
  1.96172E-06,
  8.49283E-07,
  5.02393E-07,
  2.15311E-07,
  9.56938E-08
};


double Data2011[50]={
  5.07024e+06,
  1.25241e+06,
  9.03288e+06,
  1.15106e+08,
  3.56113e+08,
  5.05445e+08,
  5.15838e+08,
  4.69266e+08,
  4.24162e+08,
  3.98722e+08,
  3.71911e+08,
  3.46841e+08,
  3.27692e+08,
  3.01809e+08,
  2.58952e+08,
  1.98323e+08,
  1.31701e+08,
  7.46387e+07,
  3.57587e+07,
  1.44475e+07,
  4.97014e+06,
  1.4923e+06,
  405908,
  104272,
  26235.1,
  6600.02,
  1659.27,
  415.404,
  109.906,
  41.2309,
  33.2132,
  43.8802,
  63.9808,
  91.6263,
  126.102,
  166.165,
  209.506,
  252.713,
  291.616,
  321.941,
  340.153,
  343.94,
  332.511,
  307.736,
  272.51,
  230.858,
  187.096,
  145.067,
  107.618,
  76.3918
};

double Data2011Up[50]={
  5.01089e+06,
  1.15407e+06,
  5.97158e+06,
  7.88451e+07,
  2.90835e+08,
  4.60041e+08,
  4.98014e+08,
  4.64675e+08,
  4.19087e+08,
  3.89226e+08,
  3.67961e+08,
  3.42229e+08,
  3.22879e+08,
  3.0469e+08,
  2.7688e+08,
  2.32927e+08,
  1.75234e+08,
  1.15356e+08,
  6.56347e+07,
  3.20356e+07,
  1.33956e+07,
  4.84198e+06,
  1.54594e+06,
  450181,
  123968,
  33354.9,
  8955.57,
  2407.16,
  643.666,
  174.949,
  56.3123,
  31.9003,
  34.8057,
  48.5878,
  69.5207,
  96.7445,
  129.626,
  166.92,
  206.5,
  245.404,
  280.142,
  307.214,
  323.747,
  327.857,
  318.729,
  297.838,
  267.405,
  230.625,
  191.063,
  152.057
};

double Data2011Down[50]={
  5.13071e+06,
  1.41095e+06,
  1.4471e+07,
  1.64152e+08,
  4.268e+08,
  5.46143e+08,
  5.288e+08,
  4.73097e+08,
  4.32847e+08,
  4.06601e+08,
  3.75696e+08,
  3.52948e+08,
  3.28549e+08,
  2.87695e+08,
  2.2503e+08,
  1.51263e+08,
  8.55545e+07,
  4.02287e+07,
  1.5665e+07,
  5.10339e+06,
  1.43125e+06,
  360929,
  85884.8,
  20077,
  4699.5,
  1096.48,
  256.996,
  69.7713,
  35.4979,
  39.7273,
  58.0342,
  85.6779,
  121.581,
  164.519,
  212.022,
  260.171,
  303.964,
  338.138,
  358.297,
  361.584,
  347.27,
  317.767,
  276.884,
  229.704,
  181.44,
  136.467,
  97.7427,
  66.6703,
  43.3105,
  26.7966
};


standalone_LumiReWeighting::standalone_LumiReWeighting(int year,int mode) {

  std::cout << "=======================================================================" << std::endl;
  
  std::vector<double> MC_distr;
  std::vector<double> Lumi_distr;

  MC_distr.clear();
  Lumi_distr.clear();
  std::cout << "Year " << year << std::endl;
  if(year!=2011 && year!=2012)
    {
      std::cout << "The year is not correct!! Reset to year 2012." << 
	std::endl;
      year=2012;
      std::cout << "Year " << year << std::endl;
    }
  switch (mode)
    {
    case 0:
      std::cout << "Using central value " << std::endl;
      break;
    case 1:
      std::cout << "Using +1 sigma 5% value " << std::endl;
      break;
    case -1:
      std::cout << "Using -1 sigma 5% value " << std::endl;
      break;
    case 2:
      std::cout << "Using old pileup JSON value" << std::endl;
      break;
    default:
      std::cout << "Using central value " << std::endl;
      break;
    } // end of switch

  Int_t NBins = 60;
  if(year==2011) NBins = 50;
  
  for( int i=0; i< NBins; ++i) {
    if(year==2011)
      {
	switch (mode){
	case 0:
	  Lumi_distr.push_back(Data2011[i]);
	  break;
	case 1:
	  Lumi_distr.push_back(Data2011Up[i]);
	  break;
	case -1:
	  Lumi_distr.push_back(Data2011Down[i]);
	  break;
	default:
	  Lumi_distr.push_back(Data2011[i]);
	  break;
	} // end of switch
	MC_distr.push_back(Fall2011[i]);
      }

    else if(year==2012)
      {
	switch (mode){
	case 0:
	  Lumi_distr.push_back(Data2012[i]);
	  break;
	case 1:
	  Lumi_distr.push_back(Data2012Up[i]);
	  break;
	case -1:
	  Lumi_distr.push_back(Data2012Down[i]);
	  break;
	case 2:
	  Lumi_distr.push_back(Data2012Buggy[i]);
	  break;
	default:
	  Lumi_distr.push_back(Data2012[i]);
	  break;
	} // end of switch
	MC_distr.push_back(Summer2012[i]);
      }

  } // end of loop over bins

  // no histograms for input: use vectors
  
  // now, make histograms out of them:

  // first, check they are the same size...

  if( MC_distr.size() != Lumi_distr.size() ){   
    std::cout << "MC_distr.size() = " << MC_distr.size() << std::endl;
    std::cout << "Lumi_distr.size() = " << Lumi_distr.size() << std::endl;
    std::cerr <<"ERROR: standalone_LumiReWeighting: input vectors have different sizes. Quitting... \n";

  }


  weights_ = new TH1D(Form("luminumer_%d",mode),
 		      Form("luminumer_%d",mode),
 		      NBins,0.0, double(NBins));

  TH1D* den = new TH1D(Form("lumidenom_%d",mode),
 		       Form("lumidenom_%d",mode),
 		       NBins,0.0, double(NBins));


  
  for(int ibin = 1; ibin<NBins+1; ++ibin ) {
    weights_->SetBinContent(ibin, Lumi_distr[ibin-1]);
    den->SetBinContent(ibin,MC_distr[ibin-1]);
  }

  std::cout << "Data Input " << std::endl;
  for(int ibin = 1; ibin<NBins+1; ++ibin){
    std::cout << "   " << ibin-1 << " " << weights_->GetBinContent(ibin) << std::endl;
  }
  std::cout << "MC Input " << std::endl;
  for(int ibin = 1; ibin<NBins+1; ++ibin){
    std::cout << "   " << ibin-1 << " " << den->GetBinContent(ibin) << std::endl;
  }

  // check integrals, make sure things are normalized

  double deltaH = weights_->Integral();
  if(fabs(1.0 - deltaH) > 0.02 ) { //*OOPS*...
    weights_->Scale( 1.0/ weights_->Integral() );
  }
  double deltaMC = den->Integral();
  if(fabs(1.0 - deltaMC) > 0.02 ) {
    den->Scale(1.0/ den->Integral());
  }

  weights_->Divide( den );  // so now the average weight should be 1.0    

  std::cout << "Reweighting: Computed Weights per In-Time Nint " << std::endl;


  for(int ibin = 1; ibin<NBins+1; ++ibin){
    std::cout << "   " << ibin-1 << " " << weights_->GetBinContent(ibin) << std::endl;
  }

  std::cout << "=======================================================================" << std::endl;

}

standalone_LumiReWeighting::~standalone_LumiReWeighting()
{
}



double standalone_LumiReWeighting::weight( double npv ) {
  int bin = weights_->GetXaxis()->FindBin( npv );
  return weights_->GetBinContent( bin );
}


#endif
