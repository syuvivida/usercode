#ifndef CONVERSIONRATE_HH
#define CONVERSIONRATE_HH

#define NErrCode -9999;
#define PErrCode  9999;

namespace syu{


  double EtaToTheta(const double eta)
  {
    
    return (2*atan(exp(-eta)));
  }


  double EtaToSintheta(const double eta)
  {
    return sin(EtaToTheta(eta));


  }

//____________________________________________________________________________
  double scaleEDep(const double epho)
  {

    double truepp, ratpp, exppp = 1.4564199;
    truepp = pairCrossSection( 13., epho );
    ratpp = truepp/exppp;

    return ratpp;

  }

//____________________________________________________________________________
  double conversionRate(const double npho, 
			const double epho,
			const double sth,
			const double mat)
  {
    double epsilon = 1e-6;
    double mysth = fabs(sth);
    if(mysth < epsilon)return NErrCode;

    double ratpp = scaleEDep(epho);
    double rate = 1. - exp((-7.*mat*npho*ratpp)/(9.*mysth));

    return rate;
  }


  //____________________________________________________________________________
  double nPho(const double cRate, const double epho, const double sth, const double mat, const double cErr, double& err)
  {
    double epsilon = 1e-6;
    double mysth = fabs(sth);
    if(mysth < epsilon)return NErrCode;
    if(cRate-1.0 > epsilon)return PErrCode;
    if(mat < epsilon)return NErrCode;

    double ratpp = scaleEDep(epho);
    if(ratpp < epsilon)return NErrCode;

    double a = (7.*mat*ratpp)/(9.*mysth);
    double n = -log(1-cRate)/a;

    err = cErr/a/(1-cRate);
    
    return n;
  }



  //____________________________________________________________________________
  double matInX0(const double cRate, const double epho, const double sth, 
		 const double cErr, double & err, const double npho=1.0)
  {
    double epsilon = 1e-6;
    double mysth = fabs(sth);
    if(mysth < epsilon)return NErrCode;
    if(cRate-1.0 > epsilon)return PErrCode;
    if(npho < epsilon)return NErrCode;

    double ratpp = scaleEDep(epho);
    if(ratpp < epsilon)return NErrCode;

    double a = (7.*npho*ratpp)/(9.*mysth);
    double mat = -log(1-cRate)/a; 

    err = cErr/a/(1-cRate);

    return mat;
  }



  //_____________________________________________________________________________
  double pairCrossSection( const double z, const double energy  ) {
    // by Bob Blair
    // ported to c++ 2/2001
    // The tables below are basically a transcription of table III.5
    //    in Rev. Mod. Phys., Vol. 46, No. 4 p815 (Oct. 1974)
    //    which summarize the pair production cross section
    //    The article was written by Yung-Su Tsai.
    const double tsaidel[14][9] = {
      0.011,0.028,0.039,0.079,0.126,0.174,0.222,0.323,0.441
      ,0.012,0.023,0.030,0.058,0.091,0.128,0.166,0.253,0.367
      ,0.004,0.024,0.034,0.073,0.113,0.154,0.195,0.283,0.391
      ,0.003,0.020,0.029,0.064,0.101,0.139,0.178,0.263,0.370
      ,0.002,0.016,0.023,0.053,0.087,0.122,0.158,0.238,0.343
      ,0.002,0.015,0.022,0.050,0.082,0.116,0.151,0.230,0.334
      ,0.002,0.012,0.019,0.044,0.073,0.105,0.137,0.213,0.315
      ,0.002,0.011,0.017,0.040,0.068,0.098,0.129,0.202,0.302
      ,0.001,0.009,0.014,0.033,0.057,0.084,0.112,0.180,0.275
      ,0.001,0.009,0.013,0.032,0.056,0.082,0.110,0.177,0.272
      ,0.001,0.008,0.012,0.029,0.051,0.075,0.102,0.166,0.259
      ,0.001,0.007,0.011,0.028,0.049,0.073,0.099,0.162,0.256
      ,0.001,0.007,0.011,0.028,0.049,0.072,0.098,0.162,0.257
      ,0.001,0.007,0.011,0.028,0.048,0.072,0.098,0.162,0.258};
    const double tsaiz[14]= { 1, 2, 3, 4, 6, 7,10,13,26,29,50,74,82,92};
    const double tsaisig[14] ={20.73,55.06,108.8,179.4,361.5,473.8,896.1,
			       1443.,5182.,6343.,17276,34869,41720,50870};
    const double tsaie[9] ={ 100.,10.,6.,2.,1.,0.6,0.4,0.2,0.1 };
    const double emass = 0.000511;
    const double rad = 1.0093;

    if(energy<=2.0*emass) return 0.0;
    double sig = 0.0;

    // find energy bin
    int je = 0;
    while(energy<tsaie[je] && je<8) je++;

    // find nearest z
    double zdif = 1.0e10;
    int jz;
    for(int j=0; j<14; j++) {
      if(fabs(tsaiz[j]-z) < fabs(zdif) ) {
	zdif = z-tsaiz[j];
	jz = j;
      }
    }

    // Z<1 HUH?
    if(jz<1 && z<tsaiz[0]) jz=1;
    // Z>92 well maybe.
    if(jz>12 && z>tsaiz[13]) jz=12;

    // interpolate z if needed
    double siginf;
    if( zdif > 0.99 ) {
      siginf = tsaisig[jz] + (z*z - tsaiz[jz]*tsaiz[jz])*(tsaisig[jz+1] -
							  tsaisig[jz])/(tsaiz[jz+1]*tsaiz[jz+1] - tsaiz[jz]*tsaiz[jz]);
    } else if( zdif < -.99 ) {
      siginf = tsaisig[jz] + (z*z-tsaiz[jz]*tsaiz[jz])*(tsaisig[jz-1] -
							tsaisig[jz])/(tsaiz[jz-1]*tsaiz[jz-1] - tsaiz[jz]*tsaiz[jz]);
    } else {
      siginf = tsaisig[jz];
    }

    // interpolate in 1/log(e)  (at low energy vs. log(e))
    double sigpl,sigmn,xpl,xmn,xx;
    if(je == 9 ) {
      // at the low end interpolate between threshold and the first point
      sigpl=tsaidel[jz][8];
      sigmn=1.0;
      xpl=log(tsaie[8]/emass);
      xmn=log(2.);
      xx=log(energy/emass);
    } else if( je == 0 ) {
      // in the high energy end interpolate to infinity via 1/log(e)
      sigpl=0.0;
      sigmn=tsaidel[jz][0];
      xmn=1.0/log(tsaie[0]/emass);
      xpl=0.0;
      xx=1.0/log(energy/emass);
    } else {
      sigpl=tsaidel[jz][je-1];
      sigmn=tsaidel[jz][je];
      xmn=1.0/log(tsaie[je]/emass);
      xpl=1.0/log(tsaie[je-1]/emass);
      xx=1.0/log(energy/emass);
    }

    // calculate adjustment to value at infinity
    sig = sigmn + (xx-xmn)*(sigpl-sigmn)/(xpl-xmn);
    // include cross section at infinity and radiative correction (rad)
    sig=(1.0-sig)*siginf*rad;
    return sig/1000.0;
    
  } // pair cross section





};

#endif
