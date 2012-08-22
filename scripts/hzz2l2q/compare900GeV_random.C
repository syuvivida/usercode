#include "setTDRStyle.C"
#include "cutvalues.h"

void compare900GeV_random(
			  std::string var="h_mh_stable_random", 
			  bool logScale=false,
			  float xmin=-9999.0, float xmax=-9999.0
	      )
{
  
  setTDRStyle();
  gStyle->SetOptStat(0);

  const int NHISTOS=4;
  TH1D* h[NHISTOS];

  const int NTEMPS=6;
  TH1D* htemp[NTEMPS];

  char tempName[300];
  TFile *fmc;
  
  std::string mcfile = "dijetmass_study/studymjj_AOD_GluGluToHToZZTo2L2Q_M-900_8TeV-powheg-pythia6.root";

  std::string mcName[NHISTOS]={
    "separated",
    "merged",
    "1-quark",
    "no match"
  };


  fmc = TFile::Open(mcfile.data());

  for(int i=0; i<NTEMPS; i++){

    htemp[i]   = (TH1D*)(fmc->Get(Form("%s%d",var.data(),i)));

  }

  h[0] = (TH1D*)htemp[0]->Clone("h0");
  h[1] = (TH1D*)htemp[1]->Clone("h1");
  h[1] -> Reset();
  h[1] -> Sumw2();
  h[1] -> Add(htemp[1],htemp[2],1.0,1.0);

  h[2] = (TH1D*)htemp[3]->Clone("h2");
  h[2] -> Reset();
  h[2] -> Sumw2();
  h[2] -> Add(htemp[3],htemp[4],1.0,1.0);

  h[3] = (TH1D*)htemp[5]->Clone("h3");

  int COLOR[NHISTOS]={4,1,2, kOrange-1};
  int MARKERSTYLE[NHISTOS]={24,21,29,8};

  for(int i=0; i < NHISTOS; i++){


    h[i]->SetTitle("mH=900 GeV");
    h[i]->GetXaxis()->SetNdivisions(5);
    h[i]->GetYaxis()->SetNdivisions(5);
//     h[i]->GetYaxis()->SetDecimals();
    h[i]->SetLineColor(COLOR[i]);
    h[i]->SetMarkerColor(COLOR[i]);
    h[i]->SetMarkerSize(1);
    h[i]->SetMarkerStyle(MARKERSTYLE[i]);
    h[i]->SetTitleOffset(1.0,"Y");

  }

  // if normalizing to the same area, set the scale 

  int binLo = -1;
  int binHi = -1;
  int nbins = h[0]->GetNbinsX();
  if(xmin>-9999.0 && xmax>-9999.0)
    {

      binLo = h[0]->FindBin(xmin);
      binHi = h[0]->FindBin(xmax)-1;

    }

  else
    {
      binLo = 1;
      binHi = nbins;
      xmin = h[0]->GetBinLowEdge(1);
      xmax = h[0]->GetBinLowEdge(nbins+1);
    }

  cout << "Before scaling" << endl;
  for(int ih=0; ih < NHISTOS; ih++)
    cout << "histogram " << ih << " entries  = " << h[ih]->GetEntries() << endl;
  for(int ih=0; ih < NHISTOS; ih++)
    cout << "histogram " << ih << " integral = " << h[ih]->Integral() << endl;

  float scaleFactor[NHISTOS]={1};

  for(int ih=0; ih < NHISTOS; ih++){
    
    float integral = h[ih]->Integral(binLo,binHi);
    scaleFactor[ih] = integral > 0? 1.0/integral: 1;
    h[ih]->Sumw2();
    h[ih]->Scale(scaleFactor[ih]);

  }

  for(int ih=0; ih < NHISTOS; ih++)
    cout << "histogram " << ih << " integral = " << h[ih]->Integral() << endl;

  for(int ih=0; ih < NHISTOS; ih++)
    cout << "histogram " << ih << " integral = " << h[ih]->Integral(binLo, 
								    binHi)
	 << endl;


  vector<float> maxArray;
  maxArray.clear();

  for(int ih=0; ih < NHISTOS; ih++){

    h[ih]->GetXaxis()->SetRangeUser(xmin,xmax);
    float max_this  = h[ih]->GetBinError(h[ih]->GetMaximumBin()) + h[ih]->GetMaximum();
    maxArray.push_back(max_this);

  }


  float max = *max_element(maxArray.begin(),maxArray.end());
  cout << "Max = " << max << endl;

  for(int ih=0; ih < NHISTOS; ih++){

    h[ih]->SetMaximum(1.1*max);

  }

  cout << "here" << endl;
  TCanvas* c1 = new TCanvas("c1","",600,500);  
  if(logScale)
    gPad->SetLogy(1);
  gStyle->SetTitleX(0.1);
  h[0]->Draw("e");
  for(int ih=0; ih < NHISTOS-1; ih++)
    h[ih]->Draw("esame");

  h[NHISTOS-1]->Draw("hesame");

  TLine mLocation(MZ_PDG,0.0,MZ_PDG,h[0]->GetMaximum());
  mLocation.SetLineStyle(2);
  mLocation.SetLineColor(kMagenta-6);
  mLocation.SetLineWidth(3);

  if(var.find("mjj")!=std::string::npos || var.find("mll")!=std::string::npos)
    mLocation.Draw("same");
  
  float x1NDC = 0.744966;
  float y1NDC = 0.764831;
  float x2NDC = 0.964765;
  float y2NDC = 0.930085;

  TLegend* leg = new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  for(int i=0; i < NHISTOS; i++)
    leg->AddEntry(h[i], mcName[i].data());
  leg->Draw("same");


  string dirName = "20120822_3histograms";
  gSystem->mkdir(dirName.data());

  std::string remword  ="h_";
  std::string remword2 ="h";
  size_t pos  = var.find(remword);
  if(pos!= std::string::npos)
    var.replace(pos,remword2.length(),"");

  std::string filename;
  std::string psname = dirName + "/mh900" + var;
  
  if(logScale)
    psname += "_log";

  filename = psname + ".eps";
  c1->Print(filename.data());
  filename = psname + ".gif";
  c1->Print(filename.data());
  filename = psname + ".pdf";
  c1->Print(filename.data());
  //   c1->Close();
}
		     
