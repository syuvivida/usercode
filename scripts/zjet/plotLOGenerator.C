
void plotLOGenerator(std::string file="weighted_genHisto_electron_genOnly_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola.root")
{

  const double fBinsPt[]={30,40,55,75,105,150,210,315,500};
  const int nPtBins = sizeof(fBinsPt)/sizeof(fBinsPt[0])-1;

  cout << "There are " << nPtBins << " bins." << endl;

  const double fBinsY[]={0.0,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4};
  const int nYBins = sizeof(fBinsY)/sizeof(fBinsY[0])-1;

  cout << "There are " << nYBins << " bins." << endl;
    

  const int nJets=4;
  TH1D* h_diff_mc_jetpt[nJets];
  TH1D* h_diff_mc_jety[nJets]; 
  
  TH1D* h_mc_jetpt[nJets];
  TH1D* h_mc_jety[nJets];


  TH1D* h_data_jetpt[nJets];
  TH1D* h_data_jety[nJets];

  TFile *fmc = TFile::Open(file.data());

  double xsec = 3048.0*1000;
  double ngen = 2.29809910000000000e+07; // madgraph, 3.38254373785480205e+06 sherpa
  
  
  for(int ij=0; ij< nJets; ij++){

    h_mc_jetpt[ij] =  (TH1D*)(fmc->Get(Form("h_mc_jetpt%02i",ij+1)));
    h_mc_jetpt[ij] -> SetName(Form("h_mc_jetpt%02i",ij+1));
    h_mc_jety[ij]  =  (TH1D*)(fmc->Get(Form("h_mc_jety%02i",ij+1)));
    h_mc_jety[ij]  -> SetName(Form("h_mc_jety%02i",ij+1));

    h_mc_jetpt[ij]->Sumw2();
    h_mc_jetpt[ij]->Scale(xsec/ngen);
    h_mc_jety[ij]->Sumw2();
    h_mc_jety[ij]->Scale(xsec/ngen);


    h_diff_mc_jetpt[ij] =  new TH1D(Form("h_diff_mc_jetpt%02i",ij+1),"",nPtBins,fBinsPt);
    h_diff_mc_jety[ij]  =  new TH1D(Form("h_diff_mc_jety%02i",ij+1),"",nYBins,fBinsY);

  }


  for(int ij=0; ij < nJets; ij++){
    for(int k=1; k<= h_mc_jetpt[ij]->GetNbinsX(); k++)
    {
      double binWidth = h_mc_jetpt[ij]->GetBinWidth(k);
      double xsec = h_mc_jetpt[ij]->GetBinContent(k);
      double diff_xsec = xsec/binWidth;

      double xsec_err = h_mc_jetpt[ij]->GetBinError(k);
      double diff_xsec_err = xsec_err/binWidth;
//       cout << "jet " << ij+1 << " pt bin " << k << ": diff_xsec = " << diff_xsec << endl;
      h_diff_mc_jetpt[ij]->SetBinContent(k,diff_xsec);
      h_diff_mc_jetpt[ij]->SetBinError(k,diff_xsec_err);
    }
  }


  for(int ij=0; ij < nJets; ij++){
    for(int k=1; k<= h_mc_jety[ij]->GetNbinsX(); k++)
    {
      double binWidth = h_mc_jety[ij]->GetBinWidth(k);
      double xsec = h_mc_jety[ij]->GetBinContent(k);
      double diff_xsec = xsec/binWidth;

      double xsec_err = h_mc_jety[ij]->GetBinError(k);
      double diff_xsec_err = xsec_err/binWidth;
      
//       cout << "jet " << ij+1 << " y bin " << k << ": diff_xsec = " << diff_xsec << endl;

      h_diff_mc_jety[ij]->SetBinContent(k,diff_xsec);
      h_diff_mc_jety[ij]->SetBinError(k,diff_xsec_err);
    }
  }
  

  cout << "Final cross check" << endl;
  cout << "as a function of pt" << endl;

  for(int ij=0; ij < 4; ij++)
    {
      cout << "=============================================" << endl;
      cout << "jet " << (ij+1) << endl;
      double total_xsec = 0;

      for(int k=1; k<= h_diff_mc_jetpt[ij]->GetNbinsX();k++){
	cout << "Bin " << k  << ": " 
	     << h_diff_mc_jetpt[ij] -> GetBinContent(k)
	     << " +- " << h_diff_mc_jetpt[ij] -> GetBinError(k)
	     << " fb/GeV" << endl;

	total_xsec += h_diff_mc_jetpt[ij] -> GetBinContent(k)*
	  h_diff_mc_jetpt[ij] ->GetBinWidth(k);
	
      }
      cout << "total xsec = " << total_xsec << " fb" << endl;
      cout << "=============================================" << endl;
    }


  cout << "as a function of Y" << endl;
  for(int ij=0; ij < 4; ij++)
    {
      cout << "=============================================" << endl;
      cout << "jet " << ij+1 << endl;
      double total_xsec = 0;
      for(int k=1; k<= h_diff_mc_jety[ij]->GetNbinsX(); k++){
	cout << "Bin " << k  << ": " 
	     << h_diff_mc_jety[ij] -> GetBinContent(k)
	     << " +- " << h_diff_mc_jety[ij] -> GetBinError(k)
	     << " fb/GeV" << endl;

	total_xsec += h_diff_mc_jety[ij] -> GetBinContent(k)*
	  h_diff_mc_jety[ij] ->GetBinWidth(k);
	
      }
      cout << "total xsec = " << total_xsec << " fb" << endl;
      cout << "=============================================" << endl;
    }

  TFile* outFile = new TFile(Form("20120423_%s",file.data()),"recreate");       
  
  for(int ij=0;ij<nJets;ij++)
    {
       h_diff_mc_jetpt[ij]->Write();
       h_diff_mc_jety[ij]->Write();

    }

  
  outFile->Close();




}
