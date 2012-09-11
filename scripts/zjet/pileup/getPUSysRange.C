#include "/afs/cern.ch/user/s/syu/scripts/histoLib.h"

void getPUSysRange(std::string mcfile,std::string var, bool update=false, double xmin=-9999.0, double xmax=-9999.0)
{
  TH1D* hcentral;
  TH1D* hup;
  TH1D* hdown;
  TH1D* hraw;

  // first get the histogram files
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
  if(xmin>-9999.0 && xmax>-9999.0)
    {

      binLo = hcentral->FindBin(xmin);
      binHi = hcentral->FindBin(xmax)-1;

    }

  else
    {
      binLo = 1;
      binHi = nbins;
      xmin = hcentral->GetBinLowEdge(1);
      xmax = hcentral->GetBinLowEdge(nbins+1);
    }


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

  TH1D* r_corr = (TH1D*) hcentral->Clone(Form("r_%s_corr",var1.data()));
  r_corr->Reset();
  r_corr->Divide(hcentral,hraw,scale_central,scale_raw);

  // now reset the error to get the spread of weights for r_corr

  getWeightedHistoErrors(r_corr, hraw, hcentral);
  getWeightedHistoErrors(r_up, hraw, hup);
  getWeightedHistoErrors(r_down, hraw, hdown);

  double minCorr  =  99999.0;
  double maxCorr  = -99999.0;

  double minSysUp =  99999.0;
  double maxSysUp = -99999.0;
  double minSysDn =  99999.0;
  double maxSysDn = -99999.0;

  double firstSysUp = -99999.0;
  double lastSysUp = -99999.0;

  double firstSysDn = -99999.0;
  double lastSysDn = -99999.0;

  cout << "binLo = " << binLo << "\t binHi = " << binHi << endl;
  
  for(int i=binLo; i <= binHi; i++){
    
    if(r_up->GetBinContent(i)<1e-6)continue;
    if(r_down->GetBinContent(i)<1e-6)continue;

    double tempcorr = r_corr->GetBinContent(i);
    if(tempcorr < minCorr)minCorr=tempcorr;
    if(tempcorr > maxCorr)maxCorr=tempcorr;

    double tempUp = r_up->GetBinContent(i)-1.0;
    double tempDn = r_down->GetBinContent(i)-1.0;
    
    if(tempUp < minSysUp)minSysUp=tempUp;
    if(tempUp > maxSysUp)maxSysUp=tempUp;

    if(tempDn < minSysDn)minSysDn=tempDn;
    if(tempDn > maxSysDn)maxSysDn=tempDn;

    if(i==binLo)firstSysUp=tempUp;
    if(i==binLo)firstSysDn=tempDn;

    if(i==binHi)lastSysUp=tempUp;
    if(i==binHi)lastSysDn=tempDn;


  }

  cout << "The range of correction for " << var << " is " << 
    fabs(maxCorr-minCorr)*100 << " % " << endl;

  cout << "The x-axis range of correction for " << var << " is " << 
    r_corr->GetBinContent(binLo) << " --- " << 
    r_corr->GetBinContent(binHi) << endl;
  cout << endl;

  cout << "The range of systematic uncertainty for up " << var << " is " << 
    minSysUp*100 << " %" << " -- " << maxSysUp*100 << " %" << endl;
  cout << "The total range of systematic uncertainty for up " << var << " is " << 
    fabs(maxSysUp-minSysUp)*100 << " %" << endl;

  cout << "The range of systematic uncertainty for down " << var << " is " << 
    minSysDn*100 << " %" << " -- " << maxSysDn*100 << " %" << endl;
  cout << "The total range of systematic uncertainty for down " << var << " is " << 
    fabs(maxSysDn-minSysDn)*100 << " %" << endl;

  cout << endl;

  cout << "The x-range of systematic uncertainty for up: " << var << " is " << 
    firstSysUp*100 << " %" << " -- " << lastSysUp*100 << " %" << endl;

  cout << "The total range of systematic uncertainty for up " << var << " is " << 
    fabs(firstSysUp-lastSysUp)*100 << " %" << endl;

  cout << "The x-range of systematic uncertainty for down: " << var << " is " << 
    firstSysDn*100 << " %" << " -- " << lastSysDn*100 << " %" << endl;

  cout << "The total range of systematic uncertainty for down: " << var << " is " << 
    fabs(firstSysDn-lastSysDn)*100 << " %" << endl;



  
  std::string command = "recreate";
  if(update)command = "update";
  TFile* outFile = new TFile(Form("myPUSys.root"),
			     command.data());       
  r_down->Write();
  r_up->Write();
  r_corr->Write();

  hcentral->Write();
  hup->Write();
  hdown->Write();
  hraw->Write();
  outFile->Close();

}
		     
