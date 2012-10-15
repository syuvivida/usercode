#include "/afs/cern.ch/user/s/syu/scripts/histoLib.h"

void produceROOT_pu(std::string var, bool update=false)
{
  TH1D* hcentral;
  TH1D* hup;
  TH1D* hdown;
  TH1D* hraw;
  std::string xtitle;

  if(var=="h_ystar")
    xtitle = "0.5|Y_{Z}-Y_{jet}|";	    
  else if(var=="h_yB")
    xtitle = "0.5|Y_{Z}+Y_{jet}|";
  else if(var=="h_jety")
    xtitle = "|Y_{jet}|";
  else if(var=="h_zy")
    xtitle = "|Y_{Z}|";

  // first get the histogram files
  std::string mcfile = "pileup/puweight_sys_ZPt40_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola.root";
  TFile *fmc  = TFile::Open(mcfile.data());  
  cout << "Opening " << fmc->GetName() << endl;

  std::string var1=var;
  std::string remword1 ="h_";
  size_t pos1  = var1.find(remword1);
  if(pos1!= std::string::npos)
    var1.replace(pos1,remword1.length(),"");

  hcentral    = (TH1D*)(fmc->Get(Form("%s_central_1",var.data())));
  hcentral    -> SetName(Form("hcentral_%s",var1.data()));
  hcentral    -> Sumw2();

  hup         = (TH1D*)(fmc->Get(Form("%s_up_1",var.data())));
  hup         -> SetName(Form("hup_%s",var1.data()));
  hup         -> Sumw2();

  hdown       = (TH1D*)(fmc->Get(Form("%s_down_1",var.data())));
  hdown       -> SetName(Form("hdown_%s",var1.data()));
  hdown       -> Sumw2();

  std::string rawfilename = mcfile;
  std::string remword  ="puweight";
  size_t pos  = rawfilename.find(remword);
  if(pos!= std::string::npos)
    rawfilename.replace(pos,remword.length(),"raw");

  TFile *fraw = TFile::Open(rawfilename.data());
  cout << "Opening " << fraw->GetName() << endl;
  hraw  = (TH1D*)(fraw->Get(Form("%s_central_1",var.data())));
  hraw  -> SetName(Form("hraw_%s",var1.data()));
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
  
  double scale_up     =1.0/hup->Integral(binLo, binHi);

  double scale_down   =1.0/hdown->Integral(binLo, binHi);

  double scale_raw    =1.0/hraw->Integral(binLo, binHi);

  cout << "After scaling = " << endl;
  cout << "Integral of hcentral = " << hcentral->Integral(binLo, binHi) << endl;  cout << "Integral of hup = " << hup->Integral(binLo, binHi) << endl;
  cout << "Integral of hdown = " << hdown->Integral(binLo, binHi) << endl;
  cout << "Integral of hraw = " << hraw->Integral(binLo, binHi) << endl;

  TH1D* r_up =(TH1D*) hcentral->Clone(Form("r_%s_up",var1.data()));
  r_up->Reset();
  r_up->Divide(hup,hcentral,scale_up,scale_central);

  TH1D* r_down =(TH1D*) hcentral->Clone(Form("r_%s_down",var1.data()));
  r_down->Reset();
  r_down->Divide(hdown,hcentral,scale_down,scale_central);

  TH1D* r_corr = (TH1D*) hcentral->Clone(Form("pu_%s_corr_stat",var1.data()));
  r_corr->Reset();
  r_corr->SetTitle("");
  r_corr->Divide(hcentral,hraw,scale_central,scale_raw);
  r_corr->SetXTitle(xtitle.data());
  r_corr->SetYTitle("PU correction with stat. error");

  TH1D* r_sys = (TH1D*) hcentral->Clone(Form("pu_%s_corr_sys",var1.data()));
  r_sys->Reset();
  r_sys->SetTitle("");
  r_sys->SetXTitle(xtitle.data());
  r_sys->SetYTitle("PU correction with stat. error");

  // now reset the error to get the spread of weights for r_corr

  getWeightedHistoErrors(r_corr, hraw, hcentral);

  
  for(int i=binLo; i <= binHi; i++){
    
    if(r_corr->GetBinContent(i)<1e-6)continue;

    double tempcorr = r_corr->GetBinContent(i);

    double tempUp = fabs(r_up->GetBinContent(i)-1.0);
    double tempDn = fabs(r_down->GetBinContent(i)-1.0);

    double max_relative_err = TMath::Max(tempUp,tempDn);

    double error = max_relative_err * tempcorr;

    r_sys->SetBinContent(i,tempcorr);
    r_sys->SetBinError(i,error);

  }



  
  std::string command = "recreate";
  if(update)command = "update";
  TFile* outFile = new TFile(Form("electron_forSteve.root"),
			     command.data());       
  r_corr->Write();
  r_sys->Write();

  outFile->Close();

}
		     
