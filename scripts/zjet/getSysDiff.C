void getSysDiff(std::string default_file, 
		std::string var, bool update=false, 
		double xmin=-9999.0, double xmax=-9999.0)
{

  std::string var1=var;
  std::string remword1 ="h_";
  size_t pos1  = var1.find(remword1);
  if(pos1!= std::string::npos)
    var1.replace(pos1,remword1.length(),"");

  TH1D* hcentral;
  TH1D* hup;
  TH1D* hdown;

  TFile *fd = TFile::Open(default_file.data());
  cout << "Opening " << fd->GetName() << endl;
  hcentral  = (TH1D*)(fd->Get(var.data()));
  hcentral  -> SetName(Form("hcentral_%s",var1.data()));
  hcentral  -> Sumw2();

  std::string sys_file_up = default_file ;
  std::string remword  =".root";
  size_t pos  = sys_file_up.find(remword);
  if(pos!= std::string::npos)
    sys_file_up.replace(pos,remword.length(),"_sfUp.root");  
//   if(pos!= std::string::npos)
//     sys_file_up.replace(pos,remword.length(),"_BkgUp.root");  
  // first get the histogram files
  TFile *fmc_up  = TFile::Open(sys_file_up.data());  
  cout << "Opening " << fmc_up->GetName() << endl;
  hup         = (TH1D*)(fmc_up->Get(var.data()));
  hup         -> SetName(Form("hup_%s",var1.data()));
  hup         -> Sumw2();


  std::string sys_file_down = default_file ;
  size_t pos2  = sys_file_down.find(remword);
  if(pos2!= std::string::npos)
    sys_file_down.replace(pos2,remword.length(),"_sfDn.root");
  TFile *fmc_down  = TFile::Open(sys_file_down.data());  
  cout << "Opening " << fmc_down->GetName() << endl;
  hdown       = (TH1D*)(fmc_down->Get(var.data()));
  hdown       -> SetName(Form("hdown_%s",var1.data()));
  hdown       -> Sumw2();

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

  double scale=1.0/hcentral->Integral(binLo, binHi);
  hcentral->Scale(scale);
  
  scale=1.0/hup->Integral(binLo, binHi);
  hup->Scale(scale);

  scale=1.0/hdown->Integral(binLo, binHi);
  hdown->Scale(scale);

  cout << "Integral of hcentral = " << hcentral->Integral(binLo, binHi) << endl;
  cout << "Integral of hup = " << hup->Integral(binLo, binHi) << endl;
  cout << "Integral of hdown = " << hdown->Integral(binLo, binHi) << endl;

  TH1D* r_up =(TH1D*) hcentral->Clone(Form("r_%s_up",var1.data()));
  r_up->Reset();
  r_up->Divide(hup,hcentral,1.0,1.0,"B");

  TH1D* r_down =(TH1D*) hcentral->Clone(Form("r_%s_down",var1.data()));
  r_down->Reset();
  r_down->Divide(hdown,hcentral,1.0,1.0,"B");



  double minSysUp =  99999.0;
  double maxSysUp = -99999.0;
  double minSysDn =  99999.0;
  double maxSysDn = -99999.0;

  double firstSysUp = -99999.0;
  double lastSysUp = -99999.0;

  double firstSysDn = -99999.0;
  double lastSysDn = -99999.0;

  for(int i=binLo; i <= binHi; i++){
    
    if(r_up->GetBinContent(i)<1e-6)continue;
    if(r_down->GetBinContent(i)<1e-6)continue;

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
  TFile* outFile = new TFile(Form("effSys.root"),
			     command.data());       
  r_down->Write();
  r_up->Write();

  hcentral->Write();
  hup->Write();
  hdown->Write();
  outFile->Close();

}
		     
