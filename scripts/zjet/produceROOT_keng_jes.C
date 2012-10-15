#include "/afs/cern.ch/user/s/syu/scripts/histoLib.h"

void produceROOT_keng_jes(std::string var, bool update=false)
{
  TH1D* hcentral;
  TH1D* hup;
  TH1D* hdown;
  TH1D* hraw;

  std::string kengName;
  std::string rawName;
  std::string xtitle;

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
  TFile *fmc  = TFile::Open("/afs/cern.ch/work/s/syu/angular_combine/mainCore/"
			    "DoubleMu2011_JESuncertainty_JetY_061712.root");  
  cout << "Opening " << fmc->GetName() << endl;

  std::string var1=var;
  std::string remword1 ="h_";
  size_t pos1  = var1.find(remword1);
  if(pos1!= std::string::npos)
    var1.replace(pos1,remword1.length(),"");

  hcentral    = (TH1D*)(fmc->Get(kengName.data()));
  hcentral    -> SetName(Form("hcentral_%s",var1.data()));
  hcentral    -> SetXTitle(xtitle.data());
  hcentral    -> Sumw2();

  hup         = (TH1D*)(fmc->Get(Form("%sUp",kengName.data())));
  hup         -> SetName(Form("hup_%s",var1.data()));
  hup         -> SetXTitle(xtitle.data());
  hup         -> Sumw2();

  hdown       = (TH1D*)(fmc->Get(Form("%sDn",kengName.data())));
  hdown       -> SetName(Form("hdown_%s",var1.data()));
  hdown       -> SetXTitle(xtitle.data());
  hdown       -> Sumw2();

  TFile *fraw = TFile::Open("/afs/cern.ch/work/s/syu/angular_combine/mainCore/"
			    "DoubleMu2011_RawJet_083112.root");
  cout << "Opening " << fraw->GetName() << endl;
  hraw  = (TH1D*)(fraw->Get(rawName.data()));
  hraw  -> SetName(Form("hraw_%s",var1.data()));
  hraw  -> SetXTitle(xtitle.data());
  hraw  -> Sumw2();

  int binLo = -1;
  int binHi = -1;
  int nbins = hcentral->GetNbinsX();
  binLo = 1;
  binHi = nbins;


  cout << "Before scaling = " << endl;
  cout << "Integral of hcentral = " << hcentral->Integral(binLo, binHi) << endl;  cout << "Integral of hup = " << hup->Integral(binLo, binHi) << endl;
  cout << "Integral of hdown = " << hdown->Integral(binLo, binHi) << endl;
  cout << "Integral of hraw = " << hraw->Integral(binLo, binHi) << endl;

  double scale_central=1.0/hcentral->Integral(binLo, binHi);
  hcentral->Scale(scale_central);
  
  double scale_up     =1.0/hup->Integral(binLo, binHi);
  hup->Scale(scale_up);

  double scale_down   =1.0/hdown->Integral(binLo, binHi);
  hdown->Scale(scale_down);

  double scale_raw    =1.0/hraw->Integral(binLo, binHi);
  hraw->Scale(scale_raw);

  cout << "After scaling = " << endl;
  cout << "Integral of hcentral = " << hcentral->Integral(binLo, binHi) << endl;  cout << "Integral of hup = " << hup->Integral(binLo, binHi) << endl;
  cout << "Integral of hdown = " << hdown->Integral(binLo, binHi) << endl;
  cout << "Integral of hraw = " << hraw->Integral(binLo, binHi) << endl;

  TH1D* r_up =(TH1D*) hcentral->Clone("r_up");
  r_up->Reset();
  r_up->SetXTitle(xtitle.data());
  r_up->Divide(hup,hcentral,1.0,1.0,"B");

  TH1D* r_down =(TH1D*) hcentral->Clone("r_down");
  r_down->Reset();
  r_down->SetXTitle(xtitle.data());
  r_down->Divide(hdown,hcentral,1.0,1.0,"B");

  TH1D* r_corr = (TH1D*) hcentral->Clone(Form("jes_%s_corr_stat",var1.data()));
  r_corr->Reset();
  r_corr->SetXTitle(xtitle.data());
  r_corr->SetYTitle("JES correction with stat. err");
  r_corr->Divide(hcentral,hraw,1.0,1.0,"B");

  TH1D* r_sys = (TH1D*) hcentral->Clone(Form("jes_%s_corr_sys",var1.data()));
  r_sys->Reset();
  r_sys->SetXTitle(xtitle.data());
  r_sys->SetYTitle("JES correction with syst. err");

  for(int i=binLo; i <= binHi; i++){
    
    if(hcentral->GetBinContent(i)<1e-6)continue;

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
  TFile* outFile = new TFile(Form("muon_forSteve.root"),
			     command.data());       
  
  r_corr->Write();
  r_sys->Write();
  outFile->Close();


}
		     
