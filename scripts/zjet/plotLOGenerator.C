
void plotLOGenerator(std::string file="weighted_genHisto_electron_genOnly_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola.root")
{

  const double fBinsPt[]={30,40,55,75,105,150,210,315,500};
  const int nPtBins = sizeof(fBinsPt)/sizeof(fBinsPt[0])-1;

  cout << "There are " << nPtBins << " bins." << endl;

  const double fBinsY[]={0.0,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7};
  const int nYBins = sizeof(fBinsY)/sizeof(fBinsY[0])-1;

  cout << "There are " << nYBins << " bins." << endl;
    

  TH1D* h_mc_jetpt;
  TH1D* h_mc_jety; 
  TH1D* h_diff_mc_jetpt;
  TH1D* h_diff_mc_jety; 
  

  TH1D* h_data_jetpt;
  TH1D* h_data_jety;

  TFile *fmc = TFile::Open(file.data());
  TFile *fdata = TFile::Open("DiffCrossSection.root");
  
  h_mc_jetpt =  (TH1D*)(fmc->Get("h_mc_jetpt"));
  h_mc_jetpt -> SetName("h_mc_jetpt");
  h_mc_jety =  (TH1D*)(fmc->Get("h_mc_jety"));
  h_mc_jety  -> SetName("h_mc_jety");

  h_data_jetpt = (TH1D*)(fdata->Get("JetPt_Z1jet"));
  h_data_jetpt -> SetName("h_data_jetpt");

  h_data_jety  = (TH1D*)(fdata->Get("JetY_Z1jet"));
  h_data_jety -> SetName("h_data_jety");

  TH1D* h_diff_mc_jetpt = (TH1D*)h_mc_jetpt->Clone("h_diff_mc_jetpt");
  TH1D* h_diff_mc_jety  = (TH1D*)h_mc_jety->Clone("h_diff_mc_jety"); 


  double xsec = 3048.0*1000;
  double ngen = 2.29809910000000000e+07; // madgraph, 3.38254373785480205e+06 sherpa
  h_mc_jetpt->Sumw2();
  h_mc_jety->Sumw2();
  
  h_mc_jetpt->Scale(xsec/ngen);
  h_mc_jety->Scale(xsec/ngen);

  for(int i=1; i<= h_mc_jetpt->GetNbinsX(); i++)
    {

      double xsec = h_mc_jetpt->GetBinContent(i);
      double diff_xsec = xsec/h_mc_jetpt->GetBinWidth(i);

      h_diff_mc_jetpt->SetBinContent(i,diff_xsec);
      h_diff_mc_jetpt->SetBinError(i,1e-4);
    }

  for(int i=1; i<= h_mc_jety->GetNbinsX(); i++)
    {

      double xsec = h_mc_jety->GetBinContent(i);
      double diff_xsec = xsec/h_mc_jety->GetBinWidth(i);

      h_diff_mc_jety->SetBinContent(i,diff_xsec);
      h_diff_mc_jety->SetBinError(i,1e-4);
    }

  TCanvas* c1 = new TCanvas("c1","",0,0,500,500);

  h_diff_mc_jetpt->SetLineColor(2);
  h_data_jetpt->SetMarkerSize(1);
  h_data_jetpt->SetMarkerStyle(8);
  h_data_jetpt->Draw("e");
  h_diff_mc_jetpt->Draw("same");
  h_data_jetpt->Draw("esame");


  TCanvas* c2 = new TCanvas("c2","",500,0,500,500);

  h_diff_mc_jety ->SetLineColor(2);
  h_data_jety ->SetMarkerSize(1);
  h_data_jety ->SetMarkerStyle(8);
  h_data_jety ->Draw("e");
  h_diff_mc_jety ->Draw("same");
  h_data_jety ->Draw("esame");

  TFile* outFile = new TFile(Form("20120420data_mc_%s",file.data()),"recreate");       
  
  h_diff_mc_jety->Write();
  h_diff_mc_jetpt->Write();
  h_data_jetpt->Write();
  h_data_jety->Write();
 
  outFile->Close();
}
