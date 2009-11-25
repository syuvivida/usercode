#include </home/syu/Photon_Code/CMSSW_3_1_2/src/syu/MyObjectCounter/test/triggerscripts/beta_asymmetry.C>
#include <TGraphAsymmErrors.h>
#include <TH1.h>

namespace SY_NT{

TGraphAsymmErrors* MyDivide(TH1F* hdeno, TH1F* hnum);

TGraphAsymmErrors* MyDivide(TH1F* hdeno, TH1F* hnum)
{

  int nbin = hdeno->GetNbinsX();
  double x[1000]={0};
  double y[1000]={0};
  double ex_l[1000]={0};
  double ex_h[1000]={0};
  double ey_l[1000]={0};
  double ey_h[1000]={0};

  for(int i=1; i<= nbin; i++)
    {
      double ntotal= hdeno->GetBinContent(i);
      double nsig  = hnum ->GetBinContent(i);
      if(ntotal==0)continue;
      double eff = nsig/ntotal;
      double errp=-999, errn=-999;
      MyEff(nsig, ntotal, eff, errp, errn);
      x[i-1] = hdeno->GetBinCenter(i);
      ex_l[i-1]= fabs(x[i-1]-hdeno->GetBinLowEdge(i));
      ex_h[i-1]= fabs(x[i-1]-hdeno->GetBinLowEdge(i+1));
      y[i-1] = eff;
      ey_l[i-1] = errn;
      ey_h[i-1] = errp;

      printf("%2.0lf to %2.0lf GeV : deno: %f, numr: %f, %1.6lf + %1.6lf - %1.6lf \n",
             hdeno->GetBinLowEdge(i),
             hdeno->GetBinLowEdge(i+1),
	     hdeno->GetBinContent(i),
	     hnum->GetBinContent(i),
             y[i-1], ey_h[i-1], ey_l[i-1]);
    }

  TGraphAsymmErrors* hanswer= new TGraphAsymmErrors(nbin,x,y,ex_l,ex_h,ey_l,ey_h);
  return hanswer;
}

}
