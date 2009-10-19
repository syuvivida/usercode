#include<math.h>
#include<float.h>
#include<assert.h>
#include "TF1.h"
#include "TMath.h"
#include <iostream.h>

double incompletebeta(double a, double b, double x2);

double inverseBeta(double p, double alpha, double beta, double A=0, double B=1);

void MyEff(double n, double N, double& eff, double& errp, double& errn);


static double ib(double a,double b,double x) {
  double n=a+b, d=a+1, prod=1, sum=1;
  while (prod>sum*DBL_EPSILON)
    sum += (prod *= n++*x/d++);
  return exp(a*log(x)+b*log(1-x)+lgamma(a+b)-lgamma(a+1)-lgamma(b))*sum;
}

double incompletebeta(double a,double b,double x) {
  assert(x>=0);assert(x<=1);assert(a>0);assert(b>0);
  return (x==0||x==1) ? x : (x*(a+b+2)<a+1) ? ib(a,b,x) : 1-ib(b,a,1-x);
}


double inverseBeta(double p, double alpha, double beta, double A, double B)
{

    double x = 0;
    double a = 0;
    double b = 1;

    A=0;
    B=1;


    double precision = DBL_EPSILON; // converge until there is 6 decimal places precision

    while ((b - a) > precision)
    {

        x = (a + b) / 2;

        if (incompletebeta(alpha, beta, x) > p)
        {

            b = x;

        }
        else
        {

            a = x;

        }

    }

    if ((B > 0) && (A > 0))
    {

        x = x * (B - A) + A;

    }
    return x;

}


void MyEff(double n, double N, double& eff, double& errp, double& errn)
{

  TF1* gaus = new TF1("gaus","TMath::Gaus(x,0,1,1)");
  double lowerB = gaus->Integral(-1,0);
  double upperB = 1- lowerB;

  eff = n/N;
  if(n > N)
    {
      errp = errn = -999;
      return;
    }

  double a = n+1;
  double b = N-n+1;
  
  double integralPeak = incompletebeta(a, b, eff);

  if(integralPeak >= lowerB && integralPeak <= upperB)
    {
      double xmin = inverseBeta(integralPeak - lowerB, a, b);
      double xmax = inverseBeta(integralPeak + lowerB, a, b);
      errp = xmax-eff;
      errn = eff-xmin;
      return;

    }
  else if(integralPeak < lowerB)
    {
      double xmax = inverseBeta(2*lowerB, a, b);      
      errp = xmax-eff;
      errn = eff;
      return;

    }

  else if(integralPeak > upperB)
    {
      double xmin = inverseBeta(1-2*lowerB, a, b);
      errp = 1-eff;
      errn = eff-xmin;
      return;      
    }

  else{

    cout << "There is a bug !!!" << endl;
    errp = errn = -999;
    return;

  }


}
