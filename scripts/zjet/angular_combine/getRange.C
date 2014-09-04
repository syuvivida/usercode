void getRange(std::string inFile, std::string histo, double xmin=-10000,
	      double xmax=-10000)
{


  double max = -999999.0;
  double min =  999999.0;

  int binmax = -1;
  int binmin = -1;
  
  TH1D* h1;


  TFile f(inFile.data());
  if(f.IsZombie()){
    cout << endl << "Error opening file" << f.GetName() << endl << endl;
    return;
  }
  else
    cout << endl << "Opened " << f.GetName() << endl << endl;

  h1 = (TH1D*)(f.Get(histo.data()));
  h1 -> SetName("h1");


  int binLo = -1;
  int binHi = -1;
  int nbins = h1->GetNbinsX();
  if(xmin>-9999.0 && xmax>-9999.0)
    {

      binLo = h1->FindBin(xmin);
      binHi = h1->FindBin(xmax)-1;

    }

  else
    {
      binLo = 1;
      binHi = nbins;
      xmin = h1->GetBinLowEdge(1);
      xmax = h1->GetBinLowEdge(nbins+1);
    }

  cout << "BinLo = " << binLo << "\t BinHi = " << binHi << endl;

  for(int i=binLo; i<= binHi; i++){
    
    if(h1->GetBinContent(i)<1e-6)continue;
    double relErr = h1->GetBinError(i) / h1->GetBinContent(i);

    if(relErr > max){
      
      binmax = i;
      max = relErr;

    }

    if(relErr < min){

      binmin = i;
      min = relErr;

    }

  }

  cout << "Bin " << binmin << " has minimum error " << min*100 << "%" << endl;
  cout << "Bin " << binmax << " has maximum error " << max*100 << "%" << endl;
  cout << endl;
  cout << min*100 << "% ~ " << max*100 << "%" << endl;
}
