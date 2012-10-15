#include "/afs/cern.ch/user/s/syu/scripts/histoLib.h"

void produceROOT_eff(std::string var, bool update=false)
{
  TH1D* heff;
  TH1D* hup;
  TH1D* hdown;

  TH1D* hraw;
  TH1D* hcentral;

  std::string rawName;
  std::string anilName;
  std::string xtitle;

  if(var=="h_ystar")
    {
      rawName  = "Ydiff_Z1Jet";
      anilName = "EfficienyVsYdif";
      xtitle = "0.5|Y_{Z}-Y_{jet}|";	    
    }
  else if(var=="h_yB")
    {
      rawName  = "Ysum_Z1Jet";
      anilName = "EfficienyVsYsum";
      xtitle = "0.5|Y_{Z}+Y_{jet}|";
    }
  else if(var=="h_jety")
    {
      rawName  = "FirstJetY_1jet";
      anilName = "EfficienyVsYjet";
      xtitle = "|Y_{jet}|";
    }
  else if(var=="h_zy")
    {
      rawName  ="ZRap_1jet";
      anilName ="EfficienyVsYz";
      xtitle = "|Y_{Z}|";
    }
 
  // first get the histogram files
  std::string var1=var;
  std::string remword1 ="h_";
  size_t pos1  = var1.find(remword1);
  if(pos1!= std::string::npos)
    var1.replace(pos1,remword1.length(),"");

  TFile *fc   = TFile::Open("/afs/cern.ch/work/s/syu/anil_rootfile/"
			    "cts_CorrectedPlotsZCut.root");
  cout << "Opening " << fc->GetName() << endl;

  hcentral    = (TH1D*)(fc->Get(var.data()));
  hcentral    -> SetName(Form("hcentral_%s",var1.data()));
  hcentral    -> SetXTitle(xtitle.data());
  hcentral    -> Sumw2();


  TFile *fraw   = TFile::Open("/afs/cern.ch/work/s/syu/anil_rootfile/"
			      "VJets_Analysed_ElecData_ZCut_Jes0.root");
  cout << "Opening " << fraw->GetName() << endl;

  hraw    = (TH1D*)(fraw->Get(rawName.data()));
  hraw    -> SetName(Form("hraw_%s",var1.data()));
  hraw    -> SetXTitle(xtitle.data());
  hraw    -> Sumw2();


  TFile *fmc  = TFile::Open("/afs/cern.ch/work/s/syu/anil_rootfile/"
			    "EffiZCut_Jer0.root");  
  cout << "Opening " << fmc->GetName() << endl;

  heff    = (TH1D*)(fmc->Get(anilName.data()));
  heff    -> SetName(Form("heff_%s",var1.data()));
  heff    -> SetXTitle(xtitle.data());
  heff    -> Sumw2();

  TFile *fup  = TFile::Open("/afs/cern.ch/work/s/syu/anil_rootfile/"
			    "EffiZCut_Jer0_sfUp.root");  
  cout << "Opening " << fup->GetName() << endl;

  hup         = (TH1D*)(fup->Get(anilName.data()));
  hup         -> SetName(Form("hup_%s",var1.data()));
  hup         -> SetXTitle(xtitle.data());
  hup         -> Sumw2();

  TFile *fdown  = TFile::Open("/afs/cern.ch/work/s/syu/anil_rootfile/"
			    "EffiZCut_Jer0_sfDn.root");  
  cout << "Opening " << fdown->GetName() << endl;

  hdown         = (TH1D*)(fdown->Get(anilName.data()));
  hdown         -> SetName(Form("hdown_%s",var1.data()));
  hdown         -> SetXTitle(xtitle.data());
  hdown         -> Sumw2();


  int binLo = -1;
  int binHi = -1;
  int nbins = heff->GetNbinsX();
  binLo = 1;
  binHi = nbins;



  double scale_eff    =1.0/heff->Integral(binLo, binHi);  
  double scale_up     =1.0/hup->Integral(binLo, binHi);
  double scale_down   =1.0/hdown->Integral(binLo, binHi);

  double scale_central=1.0/hcentral->Integral(binLo,binHi);
  double scale_raw    =1.0/hraw->Integral(binLo,binHi);
  

  TH1D* r_corr =(TH1D*) heff->Clone(Form("eff_%s_corr",var1.data()));
  r_corr->Reset();
  r_corr->SetYTitle("1/#epsilon Correction");
  r_corr->SetXTitle(xtitle.data());
  r_corr->SetTitle("");
  r_corr->Divide(hcentral,hraw,scale_central,scale_raw);

  for(int i=binLo;i<= binHi; i++){
    
    if(heff->GetBinContent(i) < 1e-6)continue;
    double ratio = r_corr->GetBinContent(i);

    double relative_err = heff->GetBinError(i)/heff->GetBinContent(i);

    double error = ratio*relative_err;

    r_corr->SetBinError(i,error);
  }

  
  std::string command = "recreate";
  if(update)command = "update";
  TFile* outFile = new TFile(Form("electron_forSteve.root"),
			     command.data());       
  r_corr->Write();

  outFile->Close();

}
		     
