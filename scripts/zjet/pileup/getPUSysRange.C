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


  cout << "binLo = " << binLo << "\t binHi = " << binHi << endl;
  
  cout << "Before scaling = " << endl;
  cout << "Integral of hcentral = " << hcentral->Integral(binLo, binHi) << endl;
  cout << "Integral of hup = " << hup->Integral(binLo, binHi) << endl;
  cout << "Integral of hdown = " << hdown->Integral(binLo, binHi) << endl;
  cout << "Integral of hraw = " << hraw->Integral(binLo, binHi) << endl;

  double scale=1.0/hcentral->Integral(binLo, binHi);
  hcentral->Scale(scale);
  
  scale=1.0/hup->Integral(binLo, binHi);
  hup->Scale(scale);

  scale=1.0/hdown->Integral(binLo, binHi);
  hdown->Scale(scale);

  scale=1.0/hraw->Integral(binLo, binHi);
  hraw->Scale(scale);

  cout << "Integral of hcentral = " << hcentral->Integral(binLo, binHi) << endl;
  cout << "Integral of hup = " << hup->Integral(binLo, binHi) << endl;
  cout << "Integral of hdown = " << hdown->Integral(binLo, binHi) << endl;
  cout << "Integral of hraw = " << hraw->Integral(binLo, binHi) << endl;

  TH1D* r_up =(TH1D*) hcentral->Clone(Form("r_%s_up",var1.data()));
  r_up->Reset();
  r_up->Divide(hup,hcentral,1.0,1.0,"B");

  TH1D* r_down =(TH1D*) hcentral->Clone(Form("r_%s_down",var1.data()));
  r_down->Reset();
  r_down->Divide(hdown,hcentral,1.0,1.0,"B");

  TH1D* r_corr = (TH1D*) hcentral->Clone(Form("r_%s_corr",var1.data()));
  r_corr->Reset();
  r_corr->Divide(hcentral,hraw,1.0,1.0,"B");

  double minSys =  99999.0;
  double maxSys = -99999.0;

  double firstSys = -1.0;
  double lastSys = -1.0;
  for(int i=binLo; i <= binHi; i++){
    
    if(r_up->GetBinContent(i)<1e-6)continue;
    if(r_down->GetBinContent(i)<1e-6)continue;

    double tempUp = fabs(r_up->GetBinContent(i)-1.0);
    double tempDn = fabs(r_down->GetBinContent(i)-1.0);
    
    if(tempUp < minSys)minSys=tempUp;
    if(tempUp > maxSys)maxSys=tempUp;

    if(tempDn < minSys)minSys=tempDn;
    if(tempDn > maxSys)maxSys=tempDn;

    if(i==binLo && tempUp > firstSys)firstSys=tempUp;
    if(i==binLo && tempDn > firstSys)firstSys=tempDn;

    if(i==binHi && tempUp > lastSys)lastSys=tempUp;
    if(i==binHi && tempDn > lastSys)lastSys=tempDn;


  }

  cout << "The range of systematic uncertainty for " << var << " is " << 
    minSys*100 << " %" << " -- " << maxSys*100 << " %" << endl;

  cout << "The x-range of systematic uncertainty for " << var << " is " << 
    firstSys*100 << " %" << " -- " << lastSys*100 << " %" << endl;
  
  
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
		     
