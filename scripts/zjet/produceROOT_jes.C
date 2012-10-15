void produceROOT_jes(std::string var, bool update=false)
{
  TH1D* h1;
  TH1D* hup;
  TH1D* hdown;
  TH1D* hcentral;
  TH1D* hraw;

  std::string xtitle;
  std::string kengName;
  std::string rawName;

  if(var=="h_ystar")
    {
      kengName = "DEta_per2_Z1jets_BE";
      rawName  = "YDiff";
      xtitle = "0.5|Y_{Z}-Y_{jet}|";	    
    }
  else if(var=="h_yB")
    {
      kengName = "SumEta_per2_Z1jets_BE";
      rawName  = "YSum";
      xtitle = "0.5|Y_{Z}+Y_{jet}|";
    }
  else if(var=="h_jety")
    {
      kengName = "Z1jets_1jeta_BE";
      rawName  = "YJet";
      xtitle = "|Y_{jet}|";
    }
  else if(var=="h_zy")
    {
      kengName ="dimuoneta1jet_BE";
      rawName  ="YZ";
      xtitle = "|Y_{Z}|";
    }
  // first get the histogram files
  TFile *fkeng  = TFile::Open("/afs/cern.ch/work/s/syu/angular_combine/mainCore/"
			      "DoubleMu2011_JESuncertainty_JetY_061712.root");  
  cout << "Opening " << fkeng->GetName() << endl;

  std::string var1=var;
  std::string remword1 ="h_";
  size_t pos1  = var1.find(remword1);
  if(pos1!= std::string::npos)
    var1.replace(pos1,remword1.length(),"");

  hcentral    = (TH1D*)(fkeng->Get(kengName.data()));
  hcentral    -> SetName(Form("hcentral_%s",var1.data()));
  hcentral    -> SetXTitle(xtitle.data());
  hcentral    -> Sumw2();


  TFile *fraw = TFile::Open("/afs/cern.ch/work/s/syu/angular_combine/mainCore/"
			    "DoubleMu2011_RawJet_083112.root");
  cout << "Opening " << fraw->GetName() << endl;
  hraw  = (TH1D*)(fraw->Get(rawName.data()));
  hraw  -> SetName(Form("hraw_%s",var1.data()));
  hraw  -> SetXTitle(xtitle.data());
  hraw  -> Sumw2();


  TFile *fmc1 = TFile::Open("/afs/cern.ch/work/s/syu/angular_combine/mainCore/cts_CorrectedPlotsZCut.root");
  TFile *fmcup = TFile::Open("/afs/cern.ch/work/s/syu/angular_combine/mainCore/cts_CorrectedPlotsZCut_JesUp_NoBkgSub.root");
  TFile *fmcdown = TFile::Open("/afs/cern.ch/work/s/syu/angular_combine/mainCore/cts_CorrectedPlotsZCut_JesDn_NoBkgSub.root");

  h1    = (TH1D*)(fmc1->Get(var.data()));
  hup   = (TH1D*)(fmcup->Get(var.data()));
  hdown = (TH1D*)(fmcdown->Get(var.data()));

  int binLo = -1;
  int binHi = -1;
  int nbins = h1->GetNbinsX();
  binLo = 1;
  binHi = nbins;



  double scale_central=1.0/hcentral->Integral(binLo, binHi);
  hcentral->Scale(scale_central);
  
  double scale_raw    =1.0/hraw->Integral(binLo, binHi);
  hraw->Scale(scale_raw);

  TH1D* r_corr = (TH1D*) hcentral->Clone(Form("jes_%s_corr_stat",var1.data()));
  r_corr->Reset();
  r_corr->SetXTitle(xtitle.data());
  r_corr->SetYTitle("JES correction with stat. err");
  r_corr->Divide(hcentral,hraw,1.0,1.0,"B");

  TH1D* r_sys = (TH1D*) hcentral->Clone(Form("jes_%s_corr_sys",var1.data()));
  r_sys->Reset();
  r_sys->SetXTitle(xtitle.data());
  r_sys->SetYTitle("JES correction with syst. err");

  double scale=h1->Integral(binLo, binHi)/h1->Integral(binLo, binHi);
  h1->Sumw2();
  h1->Scale(scale);
  
  scale=h1->Integral(binLo, binHi)/hup->Integral(binLo, binHi);
  hup->Sumw2();
  hup->Scale(scale);

  scale=h1->Integral(binLo, binHi)/hdown->Integral(binLo, binHi);
  hdown->Sumw2();
  hdown->Scale(scale);

  cout << "Integral of h1 = " << h1->Integral(binLo, binHi) << endl;
  cout << "Integral of hup = " << hup->Integral(binLo, binHi) << endl;
  cout << "Integral of hdown = " << hdown->Integral(binLo, binHi) << endl;

  TH1D* r_up =(TH1D*) h1->Clone("r_up");
  r_up->Reset();
  r_up->Divide(hup,h1,1.0,1.0,"B");

  TH1D* r_down =(TH1D*) h1->Clone("r_down");
  r_down->Reset();
  r_down->Divide(hdown,h1,1.0,1.0,"B");

  for(int i=binLo; i <= binHi; i++){
    
    if(h1->GetBinContent(i)<1e-6)continue;

    double tempcorr= r_corr->GetBinContent(i);
    double tempUp = fabs(r_up->GetBinContent(i)-1.0);
    double tempDn = fabs(r_down->GetBinContent(i)-1.0);
    double max_relative_err = TMath::Max(tempUp,tempDn);
    double error = tempcorr * max_relative_err;

    r_sys->SetBinContent(i,tempcorr);
    r_sys->SetBinError(i, error);

  }

  std::string command = "recreate";
  if(update)command = "update";
  TFile* outFile = new TFile(Form("electron_forSteve.root"),
			     command.data());       
  
  r_corr->Write();
  r_sys->Write();
  outFile->Close();

}
		     
