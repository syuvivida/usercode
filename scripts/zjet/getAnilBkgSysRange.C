#include "/afs/cern.ch/user/s/syu/scripts/histoLib.h"

void getAnilBkgSysRange(std::string var, bool update=false, 
			double xmin=-9999.0, double xmax=-9999.0)
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

  double scale_central=1.0/hcentral->Integral(binLo, binHi);  
  double scale_up     =1.0/hup->Integral(binLo, binHi);

  cout << "After scaling = " << endl;
  cout << "Integral of hcentral = " << hcentral->Integral(binLo, binHi) << endl;  cout << "Integral of hup = " << hup->Integral(binLo, binHi) << endl;

  TH1D* r_up =(TH1D*) hcentral->Clone(Form("r_%s_up",var1.data()));
  r_up->Reset();
  r_up->SetXTitle(xtitle.data());
  r_up->Divide(hup,hcentral,1.0,1.0,"B");

  // correct for the scaling
  for(int ib=1; ib<= r_up->GetNbinsX(); ib++)
    {
      double original_value = r_up->GetBinContent(ib);
      double original_error = r_up->GetBinError(ib);

      cout << "bin " << ib << " has: " << 
	original_value << " +- " << original_error << endl;
      
      double new_value = (scale_up/scale_central)*original_value;
      double new_error = (scale_up/scale_central)*original_error;

      r_up->SetBinContent(ib, new_value);
      r_up->SetBinError(ib, new_error);
    }

  cout << endl << " After scaling " << endl;

  for(int ib=1; ib<= r_up->GetNbinsX(); ib++)
    cout << r_up->GetBinContent(ib) << "\t" << 
      r_up->GetBinError(ib) << endl;


  // now reset the error to get the spread of weights for r_corr

  double minSysUp =  99999.0;
  double maxSysUp = -99999.0;


  cout << "binLo = " << binLo << "\t binHi = " << binHi << endl;
  
  for(int i=binLo; i <= binHi; i++){
    
    if(r_up->GetBinContent(i)<1e-6)continue;


    double tempUp = r_up->GetBinContent(i)-1.0;
    
    if(tempUp < minSysUp)minSysUp=tempUp;
    if(tempUp > maxSysUp)maxSysUp=tempUp;

  }


  cout << "The range of systematic uncertainty for up " << var << " is " << 
    minSysUp*100 << " %" << " -- " << maxSysUp*100 << " %" << endl;
  cout << "The total range of systematic uncertainty for up " << var << " is " << 
    fabs(maxSysUp-minSysUp)*100 << " %" << endl;

  
  std::string command = "recreate";
  if(update)command = "update";
  TFile* outFile = new TFile(Form("anilBkgSys.root"),
			     command.data());       
  r_up->Write();

  hcentral->Write();
  hup->Write();
  outFile->Close();

}
		     
