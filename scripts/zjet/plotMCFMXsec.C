
void plotMCFMXsec(std::string file="ZEE_atleast1jet_tota_mstw200_eiko.root", 
		  std::string histoName="id1")
{

  const double fBinsPt[]={30,40,55,75,105,150,210,315,500};
  const int nPtBins = sizeof(fBinsPt)/sizeof(fBinsPt[0])-1;

  cout << "There are " << nPtBins << " bins." << endl;

  const double fBinsY[]={0.0,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4};
  const int nYBins = sizeof(fBinsY)/sizeof(fBinsY[0])-1;

  cout << "There are " << nYBins << " bins." << endl;
    

  TH1D* h_mcfm_jetpt = new TH1D("h_mcfm_jetpt", "",nPtBins, fBinsPt);
  
  TH1F* h_mcfm_original;

  TH1D* h_data_jetpt;

  TFile *fmcfm = TFile::Open(file.data());
  TFile *fdata = TFile::Open("DiffCrossSection_1515.root");
  
  h_mcfm_original =  (TH1F*)(fmcfm->Get(histoName.data()));
  h_mcfm_original -> SetName("h_mcfm_original");

  h_data_jetpt = (TH1D*)(fdata->Get("r15JetPt_Z1jet"));
//   h_data_jetpt = (TH1D*)(fdata->Get("FirstJetPt_Z1jet"));
  h_data_jetpt -> SetName("h_data_jetpt");

  double default_width = h_mcfm_original->GetBinWidth(1);

  for(int i=1; i<= nPtBins-1; i++)
    {
      cout << "i = " << i << endl;
      int BinStart = h_mcfm_original->FindBin(fBinsPt[i-1]);
      int BinEnd   = h_mcfm_original->FindBin(fBinsPt[i])-1;

      double xStart = h_mcfm_original->GetBinLowEdge(BinStart);
      double xEnd   = h_mcfm_original->GetBinLowEdge(BinEnd+1);

      double Width = xEnd-xStart;

      cout << "Bin started at index " << BinStart << " at " 
	   << xStart << endl;

      cout << "Bin ended at index " << BinEnd << " at " 
	   << xEnd << endl;

      cout << "Width = " << Width << " GeV" << endl;
      
      double mcfm_diff = h_mcfm_original->Integral(BinStart,BinEnd);

      mcfm_diff *= default_width;

      cout << "Integrated cross section = " << mcfm_diff << " fb" << endl;

      if(Width<1e-6){
	cout << "Incorrect width" << endl;
	continue;
      }
      mcfm_diff /= Width;

      cout << "Differential cross section = " << mcfm_diff << " fb/GeV" << endl;

      h_mcfm_jetpt->SetBinContent(i,mcfm_diff);
      h_mcfm_jetpt->SetBinError(i,1e-4);


    }

  h_mcfm_jetpt->SetLineColor(2);
  h_data_jetpt->SetMarkerSize(1);
  h_data_jetpt->SetMarkerStyle(8);
  h_data_jetpt->Draw("e");
  h_mcfm_jetpt->Draw("same");
  h_data_jetpt->Draw("esame");

  TFile* outFile = new TFile(Form("1515data_mcfm_%s",file.data()),"recreate");       
  
  h_mcfm_original->Write();
  h_mcfm_jetpt->Write();
  h_data_jetpt->Write();
 
  outFile->Close();
}
