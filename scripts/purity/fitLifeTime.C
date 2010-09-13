Double_t exp_conv (Double_t *v, Double_t *par)
{
  Double_t ctau = par[1];
  Double_t sigma = par[3]; //using narrow prompt width
  Double_t x = v[0]-par[2];

  Double_t arg1 = TMath::Exp( 0.5*sigma*sigma/ctau/ctau - x/ctau );
  Double_t arg2 = 1.0 - TMath::Freq( (sigma/ctau - x/sigma) );
  Double_t func = par[0]/ctau * arg1 * arg2;

  if (func<=0) func=1e-10;
  return func;
}


void fitLifeTime(const char* histoName)
{
  const float fit_lo_edge =  0.1;
  const float fit_hi_edge =  3.0;

  TH1F* hfit = (TH1F*)gROOT->FindObject(histoName);
  TF1* f1 = new TF1("f1",exp_conv, fit_lo_edge, fit_hi_edge,4);
  
  f1->SetParameters(hfit->Integral(), 1., 0.5, 0.3);
  f1->SetNpx(2500);
  f1->SetLineColor(2);

  hfit->Fit("f1","","",fit_lo_edge,fit_hi_edge);
  hfit->Fit("f1","","",fit_lo_edge,fit_hi_edge);

  cout << " ==================== Fit result ==================" << endl;
  for(int i=0; i< f1-> GetNpar(); i++)
    cout << " Parameter " << i << " = " << f1->GetParameter(i) << endl;
  
  cout << " ==================================================" << endl;

}

