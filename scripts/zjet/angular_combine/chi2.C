void computeChi2(TH1F* hdata, TH1F* hmc, TH1F* hscale, int nbin, char* psname)
{
  // scale the area of the cloned histogram to 1
  char name[100];
  float fmc = 1/(float)hmc->Integral();
  float fdata = 1/(float)hdata->Integral();


  hmc->Sumw2();
  hdata->Sumw2();

  TH1F *hmcc = (TH1F*)hmc->Clone();
  hmcc->SetName("hmcc");
  hmcc->Sumw2();
  hmcc->Scale(fmc);

  TH1F *hdatac = (TH1F*)hdata->Clone();
  hdatac->SetName("hdatac");
  hdatac->Sumw2();
  hdatac->Scale(fdata);
  
  int n= hmc->GetNbinsX();
  double low = hmc->GetBinLowEdge(1);
  double high= hmc->GetBinLowEdge(1+n);
  double chi2=0;
  int realbin =0;
 

  
  // loop over n bins and compute chi2, also fill hscale
 for(int i=0;i<n;i++){

    double nmc=hmcc->GetBinContent(i+1);
    double ndata=hdatac->GetBinContent(i+1); 
    double nmcerr=hmcc->GetBinError(i+1);
    double ndataerr=hdatac->GetBinError(i+1); 
    
    if(hmc->GetBinContent(i+1)<=0 || hdata->GetBinContent(i+1)<0){
      hmcc->SetBinContent(i+1,0);
      hdatac->SetBinContent(i+1,0);
      hmcc->SetBinError(i+1,0);
      hdatac->SetBinError(i+1,0);
      continue;
    }
    
    hscale->SetBinContent(i+1,ndata/nmc);
    double err = 
	(ndata/nmc)*sqrt(pow(nmcerr/nmc,2)+pow(ndataerr/ndata,2));
    hscale->SetBinError(i+1,err);

    double chi2ndef = (nmc-ndata)*(nmc-ndata)/
      ( nmcerr*nmcerr+ ndataerr*ndataerr);
    chi2 += chi2ndef;

    cout << hmcc->GetBinCenter(i+1) << " : " << ndata;
    cout << " " << chi2ndef << endl;
    realbin++;

  }

 //  cout << "n = " << n << endl;
 cout << "real bin = " << realbin << endl;
  
 cout << "chi2=" << chi2 << endl;
 cout << "chi2/N=" << (float)chi2/(float)(realbin-1) << endl;
 double prob = hdata->KolmogorovTest(hmcc);
 cout << "prob KS = " << prob << endl;
 prob = TMath::Prob(chi2,realbin-1);
 cout << "prob = " << prob << endl;
 
  

}




