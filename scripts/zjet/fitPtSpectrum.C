

void fitPtSpectrum(std::string filename, std::string histoName)
{  
  const double fBinsPt[]={30,40,55,75,105,150,210,315,500};
  const int nPtBins = sizeof(fBinsPt)/sizeof(fBinsPt[0])-1;
  TF1* f1 = new TF1("f1","[0]*pow(x,-[1])");

  TH1F* h1;
  TFile *fmcfm = TFile::Open(filename.data());
  h1 =  (TH1F*)(fmcfm->Get(histoName.data()));

  for(int i=0; i < nPtBins; i++){

    f1->SetParameters(10000,2.0);   
    h1->Fit("f1","q","",fBinsPt[i],fBinsPt[i+1]);
    
    double p1 = f1->GetParameter(1);
    
    cout << "Bin pt = " << fBinsPt[i] << " ~ " << fBinsPt[i+1] << " GeV \t slope=" 
	 << p1 << endl;


  }




}
