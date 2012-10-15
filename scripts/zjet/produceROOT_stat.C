#include "/afs/cern.ch/user/s/syu/scripts/histoLib.h"

void produceROOT_stat(std::string var, bool update=false)
{
  TH1D* hraw;
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



  TFile *fraw   = TFile::Open("/afs/cern.ch/work/s/syu/anil_rootfile/"
			      "VJets_Analysed_ElecData_ZCut_Jes0.root");
  cout << "Opening " << fraw->GetName() << endl;

  hraw    = (TH1D*)(fraw->Get(rawName.data()));
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
  TFile* outFile = new TFile(Form("electron_forSteve.root"),
			     command.data());       
  r_corr->Write();

  outFile->Close();

}
		     
