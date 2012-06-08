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

/* Lovedeep's original array
double Data[50]={2.44414e+06,1.25579e+06,3.42957e+06,5.18618e+07,2.54054e+08,4.36586e+08,4.86031e+08,4.63551e+08,4.18993e+08,3.84891e+08,3.65724e+08,3.41505e+08,3.22711e+08,3.0886e+08,2.87693e+08,2.5129e+08,1.99438e+08,1.40551e+08,8.66577e+07,4.6234e+07,2.12058e+07,8.37396e+06,2.88178e+06,882886,246537,63900.8,15571.6,3628.24,840.61,211.248,66.7507,32.1445,26.92,33.3738,47.1181,66.9628,92.522,123.311,158.278,195.606,232.736,266.603,294.026,312.195,319.142,314.095,297.617,271.501,238.455,1.93299e+07};
*/

double Data[50]={
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

double DataUp[50]={
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

double DataDown[50]={
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


standalone_LumiReWeighting::standalone_LumiReWeighting(int mode) {

  std::cout << "=======================================================================" << std::endl;
  
  std::vector<float> MC_distr;
  std::vector<float> Lumi_distr;

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
    default:
      std::cout << "Using central value " << std::endl;
      break;
    } // end of switch

  for( int i=0; i<50; ++i) {
    switch (mode){
    case 0:
      Lumi_distr.push_back(Data[i]);
      break;
    case 1:
      Lumi_distr.push_back(DataUp[i]);
      break;
    case -1:
      Lumi_distr.push_back(DataDown[i]);
      break;
    default:
      Lumi_distr.push_back(Data[i]);
      break;
    } // end of switch
    MC_distr.push_back(Fall2011[i]);
  }


  // no histograms for input: use vectors
  
  // now, make histograms out of them:

  // first, check they are the same size...

  if( MC_distr.size() != Lumi_distr.size() ){   
    std::cout << "MC_distr.size() = " << MC_distr.size() << std::endl;
    std::cout << "Lumi_distr.size() = " << Lumi_distr.size() << std::endl;
    std::cerr <<"ERROR: standalone_LumiReWeighting: input vectors have different sizes. Quitting... \n";

  }

  Int_t NBins = MC_distr.size();

  weights_ = new TH1F(Form("luminumer_%d",mode),
		      Form("luminumer_%d",mode),
		      NBins,-0.5, float(NBins)-0.5);

  TH1F* den = new TH1F(Form("lumidenom_%d",mode),
		       Form("lumidenom_%d",mode),
		       NBins,-0.5, float(NBins)-0.5);

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

  float deltaH = weights_->Integral();
  if(fabs(1.0 - deltaH) > 0.02 ) { //*OOPS*...
    weights_->Scale( 1.0/ weights_->Integral() );
  }
  float deltaMC = den->Integral();
  if(fabs(1.0 - deltaMC) > 0.02 ) {
    den->Scale(1.0/ den->Integral());
  }

  weights_->Divide( den );  // so now the average weight should be 1.0    

  std::cout << "Reweighting: Computed Weights per In-Time Nint " << std::endl;


  for(int ibin = 1; ibin<NBins+1; ++ibin){
    std::cout << "   " << ibin-1 << " " << weights_->GetBinContent(ibin) << std::endl;
  }

  //   weightOOT_init();
  std::cout << "=======================================================================" << std::endl;

}

standalone_LumiReWeighting::~standalone_LumiReWeighting()
{
}



double standalone_LumiReWeighting::weight( int npv ) {
  int bin = weights_->GetXaxis()->FindBin( npv );
  return weights_->GetBinContent( bin );
}

void standalone_LumiReWeighting::weightOOT_init() {

  // The following are poisson distributions with different means, where the maximum
  // of the function has been normalized to weight 1.0
  // These are used to reweight the out-of-time pileup to match the in-time distribution.
  // The total event weight is the product of the in-time weight, the out-of-time weight,
  // and a residual correction to fix the distortions caused by the fact that the out-of-time
  // distribution is not flat.

  static double weight_24[25] = {
    0,
    0,
    0,
    0,
    2.46277e-06,
    2.95532e-05,
    0.000104668,
    0.000401431,
    0.00130034,
    0.00342202,
    0.00818132,
    0.0175534,
    0.035784,
    0.0650836,
    0.112232,
    0.178699,
    0.268934,
    0.380868,
    0.507505,
    0.640922,
    0.768551,
    0.877829,
    0.958624,
    0.99939,
    1
  };

  static double weight_23[25] = {
    0,
    1.20628e-06,
    1.20628e-06,
    2.41255e-06,
    1.20628e-05,
    6.39326e-05,
    0.000252112,
    0.000862487,
    0.00244995,
    0.00616527,
    0.0140821,
    0.0293342,
    0.0564501,
    0.100602,
    0.164479,
    0.252659,
    0.36268,
    0.491427,
    0.627979,
    0.75918,
    0.873185,
    0.957934,
    0.999381,
    1,
    0.957738
  };

  static double weight_22[25] = {
    0,
    0,
    0,
    5.88636e-06,
    3.0609e-05,
    0.000143627,
    0.000561558,
    0.00173059,
    0.00460078,
    0.0110616,
    0.0238974,
    0.0475406,
    0.0875077,
    0.148682,
    0.235752,
    0.343591,
    0.473146,
    0.611897,
    0.748345,
    0.865978,
    0.953199,
    0.997848,
    1,
    0.954245,
    0.873688
  };

  static double weight_21[25] = {
    0,
    0,
    1.15381e-06,
    8.07665e-06,
    7.1536e-05,
    0.000280375,
    0.00107189,
    0.00327104,
    0.00809396,
    0.0190978,
    0.0401894,
    0.0761028,
    0.13472,
    0.216315,
    0.324649,
    0.455125,
    0.598241,
    0.739215,
    0.861866,
    0.953911,
    0.998918,
    1,
    0.956683,
    0.872272,
    0.76399
  };
 
 
  static double weight_20[25] = {
    0,
    0,
    1.12532e-06,
    2.58822e-05,
    0.000145166,
    0.000633552,
    0.00215048,
    0.00592816,
    0.0145605,
    0.0328367,
    0.0652649,
    0.11893,
    0.19803,
    0.305525,
    0.436588,
    0.581566,
    0.727048,
    0.8534,
    0.949419,
    0.999785,
    1,
    0.953008,
    0.865689,
    0.753288,
    0.62765
  }; 
  static double weight_19[25] = {
    0,
    0,
    1.20714e-05,
    5.92596e-05,
    0.000364337,
    0.00124994,
    0.00403953,
    0.0108149,
    0.025824,
    0.0544969,
    0.103567,
    0.17936,
    0.283532,
    0.416091,
    0.562078,
    0.714714,
    0.846523,
    0.947875,
    1,
    0.999448,
    0.951404,
    0.859717,
    0.742319,
    0.613601,
    0.48552
  };

  static double weight_18[25] = {
    0,
    3.20101e-06,
    2.88091e-05,
    0.000164319,
    0.000719161,
    0.00250106,
    0.00773685,
    0.0197513,
    0.0443693,
    0.0885998,
    0.159891,
    0.262607,
    0.392327,
    0.543125,
    0.69924,
    0.837474,
    0.943486,
    0.998029,
    1,
    0.945937,
    0.851807,
    0.729309,
    0.596332,
    0.467818,
    0.350434
  };

 
  static double weight_17[25] = {
    1.03634e-06,
    7.25437e-06,
    4.97443e-05,
    0.000340956,
    0.00148715,
    0.00501485,
    0.0143067,
    0.034679,
    0.0742009,
    0.140287,
    0.238288,
    0.369416,
    0.521637,
    0.682368,
    0.828634,
    0.939655,
    1,
    0.996829,
    0.94062,
    0.841575,
    0.716664,
    0.582053,
    0.449595,
    0.331336,
    0.234332
  };

 
  static double weight_16[25] = {
    4.03159e-06,
    2.41895e-05,
    0.000141106,
    0.00081942,
    0.00314565,
    0.00990662,
    0.026293,
    0.0603881,
    0.120973,
    0.214532,
    0.343708,
    0.501141,
    0.665978,
    0.820107,
    0.938149,
    1,
    0.99941,
    0.940768,
    0.837813,
    0.703086,
    0.564023,
    0.42928,
    0.312515,
    0.216251,
    0.14561
  };
 
 
  static double weight_15[25] = {
    9.76084e-07,
    5.07564e-05,
    0.000303562,
    0.00174036,
    0.00617959,
    0.0188579,
    0.047465,
    0.101656,
    0.189492,
    0.315673,
    0.474383,
    0.646828,
    0.809462,
    0.934107,
    0.998874,
    1,
    0.936163,
    0.827473,
    0.689675,
    0.544384,
    0.40907,
    0.290648,
    0.198861,
    0.12951,
    0.0808051
  };
 
 
  static double weight_14[25] = {
    1.13288e-05,
    0.000124617,
    0.000753365,
    0.00345056,
    0.0123909,
    0.0352712,
    0.0825463,
    0.16413,
    0.287213,
    0.44615,
    0.625826,
    0.796365,
    0.930624,
    0.999958,
    1,
    0.934414,
    0.816456,
    0.672939,
    0.523033,
    0.386068,
    0.269824,
    0.180342,
    0.114669,
    0.0698288,
    0.0406496
  };

 
  static double weight_13[25] = {
    2.54296e-05,
    0.000261561,
    0.00167018,
    0.00748083,
    0.0241308,
    0.0636801,
    0.138222,
    0.255814,
    0.414275,
    0.600244,
    0.779958,
    0.92256,
    0.999155,
    1,
    0.927126,
    0.804504,
    0.651803,
    0.497534,
    0.35976,
    0.245834,
    0.160904,
    0.0991589,
    0.0585434,
    0.0332437,
    0.0180159
  };

  static double weight_12[25] = {
    5.85742e-05,
    0.000627706,
    0.00386677,
    0.0154068,
    0.0465892,
    0.111683,
    0.222487,
    0.381677,
    0.5719,
    0.765001,
    0.915916,
    1,
    0.999717,
    0.921443,
    0.791958,
    0.632344,
    0.475195,
    0.334982,
    0.223666,
    0.141781,
    0.0851538,
    0.048433,
    0.0263287,
    0.0133969,
    0.00696683
  };

 
  static double weight_11[25] = {
    0.00015238,
    0.00156064,
    0.00846044,
    0.0310939,
    0.0856225,
    0.187589,
    0.343579,
    0.541892,
    0.74224,
    0.909269,
    0.998711,
    1,
    0.916889,
    0.77485,
    0.608819,
    0.447016,
    0.307375,
    0.198444,
    0.121208,
    0.070222,
    0.0386492,
    0.0201108,
    0.0100922,
    0.00484937,
    0.00222458
  };

  static double weight_10[25] = {
    0.000393044,
    0.00367001,
    0.0179474,
    0.060389,
    0.151477,
    0.302077,
    0.503113,
    0.720373,
    0.899568,
    1,
    0.997739,
    0.909409,
    0.75728,
    0.582031,
    0.415322,
    0.277663,
    0.174147,
    0.102154,
    0.0566719,
    0.0298642,
    0.0147751,
    0.00710995,
    0.00319628,
    0.00140601,
    0.000568796
  };

 
  static double weight_9[25] = {
    0.00093396,
    0.00854448,
    0.0380306,
    0.113181,
    0.256614,
    0.460894,
    0.690242,
    0.888781,
    1,
    0.998756,
    0.899872,
    0.735642,
    0.552532,
    0.382726,
    0.246114,
    0.147497,
    0.0825541,
    0.0441199,
    0.0218157,
    0.0103578,
    0.00462959,
    0.0019142,
    0.000771598,
    0.000295893,
    0.000111529
  };

 
  static double weight_8[25] = {
    0.00240233,
    0.0192688,
    0.0768653,
    0.205008,
    0.410958,
    0.65758,
    0.875657,
    0.999886,
    1,
    0.889476,
    0.711446,
    0.517781,
    0.345774,
    0.212028,
    0.121208,
    0.0644629,
    0.0324928,
    0.0152492,
    0.00673527,
    0.0028547,
    0.00117213,
    0.000440177,
    0.000168471,
    5.80689e-05,
    1.93563e-05
  };

  static double weight_7[25] = {
    0.00617233,
    0.0428714,
    0.150018,
    0.350317,
    0.612535,
    0.856525,
    0.999923,
    1,
    0.87544,
    0.679383,
    0.478345,
    0.303378,
    0.176923,
    0.0950103,
    0.0476253,
    0.0222211,
    0.00972738,
    0.00392962,
    0.0015258,
    0.000559168,
    0.000183928,
    6.77983e-05,
    1.67818e-05,
    7.38398e-06,
    6.71271e-07
  };
 
  static double weight_6[25] = {
    0.0154465,
    0.0923472,
    0.277322,
    0.55552,
    0.833099,
    0.999035,
    1,
    0.855183,
    0.641976,
    0.428277,
    0.256804,
    0.139798,
    0.0700072,
    0.0321586,
    0.0137971,
    0.00544756,
    0.00202316,
    0.000766228,
    0.000259348,
    8.45836e-05,
    1.80362e-05,
    8.70713e-06,
    3.73163e-06,
    6.21938e-07,
    0
  };
 
 
  static double weight_5[25] = {
    0.0382845,
    0.191122,
    0.478782,
    0.797314,
    1,
    0.997148,
    0.831144,
    0.59461,
    0.371293,
    0.205903,
    0.103102,
    0.0471424,
    0.0194997,
    0.00749415,
    0.00273709,
    0.000879189,
    0.000286049,
    0.000102364,
    1.70606e-05,
    3.98081e-06,
    2.27475e-06,
    0,
    0,
    0,
    0
  };
 
 
  static double weight_4[25] = {
    0.0941305,
    0.373824,
    0.750094,
    1,
    0.997698,
    0.800956,
    0.532306,
    0.304597,
    0.152207,
    0.0676275,
    0.0270646,
    0.00975365,
    0.00326077,
    0.00101071,
    0.000301781,
    7.41664e-05,
    1.58563e-05,
    3.58045e-06,
    1.02299e-06,
    0,
    5.11493e-07,
    0,
    0,
    0,
    0
  };
 
 
  static double weight_3[25] = {
    0.222714,
    0.667015,
    1,
    0.999208,
    0.750609,
    0.449854,
    0.224968,
    0.0965185,
    0.0361225,
    0.012084,
    0.00359618,
    0.000977166,
    0.000239269,
    6.29422e-05,
    1.16064e-05,
    1.78559e-06,
    0,
    4.46398e-07,
    0,
    0,
    0,
    0,
    0,
    0,
    0
  };
 
  static double weight_2[25] = {
    0.499541,
    0.999607,
    1,
    0.666607,
    0.333301,
    0.13279,
    0.0441871,
    0.0127455,
    0.00318434,
    0.00071752,
    0.000132204,
    2.69578e-05,
    5.16999e-06,
    2.21571e-06,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0
  };
 
  static double weight_1[25] = {
    0.999165,
    1,
    0.499996,
    0.166868,
    0.0414266,
    0.00831053,
    0.00137472,
    0.000198911,
    2.66302e-05,
    2.44563e-06,
    2.71737e-07,
    2.71737e-07,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0
  };
 
  static double weight_0[25] = {
    1,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0
  };


  double* WeightPtr = 0;

  for(int iint = 0; iint<25; ++iint){
    if(iint ==0) WeightPtr = weight_0;
    if(iint ==1) WeightPtr = weight_1;
    if(iint ==2) WeightPtr = weight_2;
    if(iint ==3) WeightPtr = weight_3;
    if(iint ==4) WeightPtr = weight_4;
    if(iint ==5) WeightPtr = weight_5;
    if(iint ==6) WeightPtr = weight_6;
    if(iint ==7) WeightPtr = weight_7;
    if(iint ==8) WeightPtr = weight_8;
    if(iint ==9) WeightPtr = weight_9;
    if(iint ==10) WeightPtr = weight_10;
    if(iint ==11) WeightPtr = weight_11;
    if(iint ==12) WeightPtr = weight_12;
    if(iint ==13) WeightPtr = weight_13;
    if(iint ==14) WeightPtr = weight_14;
    if(iint ==15) WeightPtr = weight_15;
    if(iint ==16) WeightPtr = weight_16;
    if(iint ==17) WeightPtr = weight_17;
    if(iint ==18) WeightPtr = weight_18;
    if(iint ==19) WeightPtr = weight_19;
    if(iint ==20) WeightPtr = weight_20;
    if(iint ==21) WeightPtr = weight_21;
    if(iint ==22) WeightPtr = weight_22;
    if(iint ==23) WeightPtr = weight_23;
    if(iint ==24) WeightPtr = weight_24;

  }

}

#endif
