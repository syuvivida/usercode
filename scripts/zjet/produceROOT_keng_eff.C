
void produceROOT_keng_eff(std::string var, bool update=false)
{
  TH1F* hcentral;
  TH1F* hraw;
  TH1F* heff_mu;

  std::string xtitle;
  std::string kengName;
  std::string effMuonName;

  if(var=="h_ystar")
    {
      xtitle = "0.5|Y_{Z}-Y_{jet}|";
      kengName = "DEta_per2_Z1jets_BE";
      effMuonName = "EffCorr_Ydiff";
    }

  else if(var=="h_yB")
    {
      xtitle = "0.5|Y_{Z}+Y_{jet}|";
      kengName = "SumEta_per2_Z1jets_BE";
      effMuonName = "EffCorr_Ysum";
    }
  else if(var=="h_jety")
    {
      xtitle = "|Y(jet)|";
      kengName = "Z1jets_1jeta_BE";
      effMuonName = "EffCorr_Yjet";
    }
  else if(var=="h_zy")
    {
      xtitle = "|Y(Z)|";
      kengName ="dimuoneta1jet_BE";
      effMuonName = "EffCorr_Yz";
    }

 
  // first get the histogram files
  TFile *fmc  = TFile::Open("/afs/cern.ch/work/s/syu/angular_combine/mainCore/DYJetsToLL_EffCorr_091712.root");  
  cout << "Opening " << fmc->GetName() << endl;

  std::string var1=var;
  std::string remword1 ="h_";
  size_t pos1  = var1.find(remword1);
  if(pos1!= std::string::npos)
    var1.replace(pos1,remword1.length(),"");

  hcentral    = (TH1F*)(fmc->Get(kengName.data()));
  hcentral    -> SetName(Form("hcentral_%s",var1.data()));
  hcentral    -> SetXTitle(xtitle.data());
  hcentral    -> Sumw2();


  TFile *fraw = TFile::Open("/afs/cern.ch/work/s/syu/angular_combine/mainCore/"
			    "DYJetsToLL_WithOutID_CorrErr_091612.root");
  cout << "Opening " << fraw->GetName() << endl;
  hraw  = (TH1F*)(fraw->Get(kengName.data()));
  hraw  -> SetName(Form("hraw_%s",var1.data()));
  hraw  -> SetXTitle(xtitle.data());
  hraw  -> Sumw2();

  TFile *f_eff_muon = TFile::Open("/afs/cern.ch/work/s/syu/angular_combine/mainCore/EffCorr_root_091912.root");

  heff_mu = (TH1F*)(f_eff_muon->Get(effMuonName.data()));
  heff_mu -> SetName("heff_mu");


  int binLo = -1;
  int binHi = -1;
  int nbins = hcentral->GetNbinsX();
  binLo = 1;
  binHi = nbins;

  cout << "Before scaling = " << endl;
  cout << "Integral of hcentral = " << hcentral->Integral(binLo, binHi) << endl;
  cout << "Integral of hraw = " << hraw->Integral(binLo, binHi) << endl;

  double scale_central=1.0/hcentral->Integral(binLo, binHi);  
  double scale_raw    =1.0/hraw->Integral(binLo, binHi);


  TH1F* r_corr = (TH1F*) heff_mu->Clone(Form("eff_%s_corr",var1.data()));
  r_corr->Reset();
  r_corr->SetTitle("");
  r_corr->SetYTitle("1/#epsilon Correction");
  r_corr->SetXTitle(xtitle.data());
  r_corr->Divide(hcentral,hraw,scale_central,scale_raw);


  for(int i=binLo; i <= binHi; i++){
    

    if(hcentral->GetBinContent(i)<1e-6)continue;
    if(heff_mu->GetBinContent(i)<1e-6)continue;

    double relative_err = heff_mu->GetBinError(i)/heff_mu->GetBinContent(i);
    double error = relative_err * r_corr->GetBinContent(i);
    r_corr->SetBinError(i,error);

  }

 
  std::string command = "recreate";
  if(update)command = "update";
  TFile* outFile = new TFile(Form("muon_forSteve.root"),
			     command.data());       
  r_corr->Write();

  outFile->Close();

}
		     
