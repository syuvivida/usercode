#include "/afs/cern.ch/user/s/syu/scripts/histoLib.h"

void getAnilEffSysRange(std::string var, bool update=false, double xmin=-9999.0, double xmax=-9999.0)
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
  if(xmin>-9999.0 && xmax>-9999.0)
    {

      binLo = heff->FindBin(xmin);
      binHi = heff->FindBin(xmax)-1;

    }

  else
    {
      binLo = 1;
      binHi = nbins;
      xmin = heff->GetBinLowEdge(1);
      xmax = heff->GetBinLowEdge(nbins+1);
    }


  cout << "Before scaling = " << endl;
  cout << "Integral of heff = " << heff->Integral(binLo, binHi) << endl;  cout << "Integral of hup = " << hup->Integral(binLo, binHi) << endl;
  cout << "Integral of hdown = " << hdown->Integral(binLo, binHi) << endl;

  double scale_eff    =1.0/heff->Integral(binLo, binHi);  
  double scale_up     =1.0/hup->Integral(binLo, binHi);
  double scale_down   =1.0/hdown->Integral(binLo, binHi);

  double scale_central=1.0/hcentral->Integral(binLo,binHi);
  double scale_raw    =1.0/hraw->Integral(binLo,binHi);
  

  cout << "After scaling = " << endl;
  cout << "Integral of heff = " << heff->Integral(binLo, binHi) << endl;  cout << "Integral of hup = " << hup->Integral(binLo, binHi) << endl;
  cout << "Integral of hdown = " << hdown->Integral(binLo, binHi) << endl;

  TH1D* r_corr =(TH1D*) heff->Clone(Form("r_%s_corr",var1.data()));
  r_corr->Reset();
  r_corr->SetXTitle(xtitle.data());
  r_corr->SetTitle("");
  r_corr->Divide(hcentral,hraw,scale_central,scale_raw);

  for(int i=1;i<= r_corr->GetNbinsX(); i++){
    
    if(heff->GetBinContent(i) < 1e-6)continue;
    double ratio = r_corr->GetBinContent(i);

    double relative_err = heff->GetBinError(i)/heff->GetBinContent(i);

    double error = ratio*relative_err;

    r_corr->SetBinError(i,error);
  }

  // note, unlike the other uncertainty, here we are showing the difference, 
  // not the ratio
  TH1D* r_up   =(TH1D*) heff->Clone(Form("r_%s_up",var1.data()));
  r_up->Reset();
  r_up->SetXTitle(xtitle.data());

  TH1D* r_down   =(TH1D*) heff->Clone(Form("r_%s_down",var1.data()));
  r_down->Reset();
  r_down->SetXTitle(xtitle.data());

  for(int ib=1; ib<= r_corr->GetNbinsX(); ib++){

    if(heff->GetBinContent(ib)<1e-6)continue;

    double upDiff = hup->GetBinContent(ib)*scale_up - 
      heff->GetBinContent(ib)*scale_eff;

    double dnDiff = hdown->GetBinContent(ib)*scale_down - 
      heff->GetBinContent(ib)*scale_eff;

    double statErr = heff->GetBinError(ib);

    double ErrUp = 
      sqrt(upDiff*upDiff+ statErr*statErr)/heff->GetBinContent(ib);

    double ErrDn = 
      sqrt(dnDiff*dnDiff+ statErr*statErr)/heff->GetBinContent(ib);

    r_up  ->SetBinContent(ib, ErrUp);

    r_down->SetBinContent(ib, ErrDn);

  }


  double minCorr  =  99999.0;
  double maxCorr  = -99999.0;

  double minSysUp =  99999.0;
  double maxSysUp = -99999.0;
  double minSysDn =  99999.0;
  double maxSysDn = -99999.0;

  cout << "binLo = " << binLo << "\t binHi = " << binHi << endl;
  
  for(int i=binLo; i <= binHi; i++){
    
    if(r_corr->GetBinContent(i)<1e-6)continue;

    double tempcorr = r_corr->GetBinContent(i)-1.0;
    if(tempcorr < minCorr)minCorr=tempcorr;
    if(tempcorr > maxCorr)maxCorr=tempcorr;

    double tempUp = r_up->GetBinContent(i);
    double tempDn = r_down->GetBinContent(i);
    
    if(tempUp < minSysUp)minSysUp=tempUp;
    if(tempUp > maxSysUp)maxSysUp=tempUp;

    if(tempDn < minSysDn)minSysDn=tempDn;
    if(tempDn > maxSysDn)maxSysDn=tempDn;

  }

  cout << "maxCorr = " << maxCorr << "\t minCorr = " << minCorr << endl;
  cout << "The range of correction for " << var << " is " << 
    fabs(maxCorr-minCorr)*100 << " % " << endl;

  cout << "The range of systematic uncertainty for up " << var << " is " << 
    minSysUp*100 << " %" << " -- " << maxSysUp*100 << " %" << endl;
  cout << "The total range of systematic uncertainty for up " << var << " is " << 
    fabs(maxSysUp-minSysUp)*100 << " %" << endl;

  cout << "The range of systematic uncertainty for down " << var << " is " << 
    minSysDn*100 << " %" << " -- " << maxSysDn*100 << " %" << endl;
  cout << "The total range of systematic uncertainty for down " << var << " is " << 
    fabs(maxSysDn-minSysDn)*100 << " %" << endl;

  cout << endl;

  
  std::string command = "recreate";
  if(update)command = "update";
  TFile* outFile = new TFile(Form("anilEffSys.root"),
			     command.data());       
  r_down->Write();
  r_up->Write();
  r_corr->Write();

  heff->Write();
  hup->Write();
  hdown->Write();

  hcentral->Write();
  hraw->Write();

  outFile->Close();

}
		     
