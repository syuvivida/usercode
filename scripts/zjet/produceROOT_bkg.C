#include "/afs/cern.ch/user/s/syu/scripts/histoLib.h"

void produceROOT_bkg(std::string var, bool update=false)
{
  TH1D* hcentral;
  TH1D* hup;

  std::string rawName;
  std::string xtitle;

  if(var=="h_ystar")
    {
      rawName  = "Ydiff_Z1Jet";
      xtitle = "0.5|Y_{Z}-Y_{jet}|";	    
    }
  else if(var=="h_yB")
    {
      rawName  = "Ysum_Z1Jet";
      xtitle = "0.5|Y_{Z}+Y_{jet}|";
    }
  else if(var=="h_jety")
    {
      rawName  = "FirstJetY_1jet";
      xtitle = "|Y_{jet}|";
    }
  else if(var=="h_zy")
    {
      rawName  ="ZRap_1jet";
      xtitle = "|Y_{Z}|";
    }
 
  // first get the histogram files
  TFile *fmc  = TFile::Open("/afs/cern.ch/work/s/syu/anil_rootfile/"
			    "VJets_Analysed_ElecData_ZCut_Jes0.root");  
  cout << "Opening " << fmc->GetName() << endl;

  std::string var1=var;
  std::string remword1 ="h_";
  size_t pos1  = var1.find(remword1);
  if(pos1!= std::string::npos)
    var1.replace(pos1,remword1.length(),"");

  hcentral    = (TH1D*)(fmc->Get(rawName.data()));
  hcentral    -> SetName(Form("hcentral_%s",var1.data()));
  hcentral    -> SetXTitle(xtitle.data());
  hcentral    -> Sumw2();

  TFile *fup = TFile::Open("/afs/cern.ch/work/s/syu/anil_rootfile/"
			   "BkgSubSignalZCut_BkgUp.root");
  cout << "Opening " << fup->GetName() << endl;
  hup  = (TH1D*)(fup->Get(rawName.data()));
  hup  -> SetName(Form("hup_%s",var1.data()));
  hup  -> SetXTitle(xtitle.data());
  hup  -> Sumw2();

  int binLo = -1;
  int binHi = -1;
  int nbins = hcentral->GetNbinsX();
  binLo = 1;
  binHi = nbins;


  cout << "Before scaling = " << endl;
  cout << "Integral of hcentral = " << hcentral->Integral(binLo, binHi) << endl;  cout << "Integral of hup = " << hup->Integral(binLo, binHi) << endl;

  double scale_central=1.0/hcentral->Integral(binLo, binHi);  
  double scale_up     =1.0/hup->Integral(binLo, binHi);

  cout << "After scaling = " << endl;
  cout << "Integral of hcentral = " << hcentral->Integral(binLo, binHi) << endl;  cout << "Integral of hup = " << hup->Integral(binLo, binHi) << endl;

  TH1D* r_corr =(TH1D*) hcentral->Clone(Form("bkg_%s_corr",var1.data()));
  r_corr->Reset();
  r_corr->SetXTitle(xtitle.data());
  r_corr->SetYTitle("Background correction");

  // correct for the scaling
  for(int ib=binLo; ib<= binHi; ib++)
    {
      if(hcentral->GetBinContent(ib)<1e-6)continue;
      r_corr->SetBinContent(ib,1.0);
      
      double original_value = (hup->GetBinContent(ib)*scale_up)/
	(hcentral->GetBinContent(ib)*scale_central)-1.0;
      double error = fabs(original_value);

      r_corr->SetBinError(ib, error);
    }

  
  std::string command = "recreate";
  if(update)command = "update";
  TFile* outFile = new TFile(Form("electron_forSteve.root"),
			     command.data());       
  r_corr->Write();
  outFile->Close();

}
		     
