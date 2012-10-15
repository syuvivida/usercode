#include "/afs/cern.ch/user/s/syu/scripts/histoLib.h"

void produceROOT_keng_stat(std::string var, bool update=false)
{
  TH1D* hraw;

 
  // first get the histogram files
  std::string var1=var;
  std::string remword1 ="h_";
  size_t pos1  = var1.find(remword1);
  if(pos1!= std::string::npos)
    var1.replace(pos1,remword1.length(),"");

  std::string xtitle;
  std::string muName;

  if(var=="h_ystar")
    {
      xtitle = "0.5|Y_{Z}-Y_{jet}|";
      muName = "YDiff";
    }

  else if(var=="h_yB")
    {
      xtitle = "0.5|Y_{Z}+Y_{jet}|";
      muName = "YSum";
    }
  else if(var=="h_jety")
    {
      xtitle = "|Y(jet)|";
      muName = "Yjet";
    }
  else if(var=="h_zy")
    {
      xtitle = "|Y(Z)|";
      muName = "YZ";
    }


  TFile *fraw   = TFile::Open("/afs/cern.ch/work/s/syu/angular_combine/mainCore/DoubleMu2011_EffCorr_091812.root");
  cout << "Opening " << fraw->GetName() << endl;

  hraw    = (TH1D*)(fraw->Get(muName.data()));
  hraw    -> SetName(Form("hraw_%s",var1.data()));
  hraw    -> SetXTitle(xtitle.data());
  hraw    -> Sumw2();


  int binLo = -1;
  int binHi = -1;
  int nbins = hraw->GetNbinsX();
  binLo = 1;
  binHi = nbins;



  TH1D* r_corr =(TH1D*) hraw->Clone(Form("stat_%s_corr",var1.data()));
  r_corr->Reset();
  r_corr->SetYTitle("statistical");
  r_corr->SetXTitle(xtitle.data());
  r_corr->SetTitle("");

  for(int i=binLo;i<= binHi; i++){
    
    if(hraw->GetBinContent(i)<1e-6)continue;
    r_corr->SetBinContent(i,1.0);

    double relative_error = hraw->GetBinError(i)/hraw->GetBinContent(i);

    r_corr->SetBinError(i,relative_error);


  }

  
  std::string command = "recreate";
  if(update)command = "update";
  TFile* outFile = new TFile(Form("muon_forSteve.root"),
			     command.data());       
  r_corr->Write();

  outFile->Close();

}
		     
