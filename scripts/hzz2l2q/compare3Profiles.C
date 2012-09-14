#include "setTDRStyle.C"
#include "cutvalues.h"

void compare3Profiles(
		      std::string var="pf_dR_Rm_gen", 
		      float xmin=0.5, float xmax=3.0
		      )
{
  
  setTDRStyle();
  gStyle->SetOptStat(0);

  const int NHISTOS=3;
  TProfile* h[NHISTOS];

  char tempName[300];
  TFile *fmc[NHISTOS];
  
  std::string mcfile[NHISTOS]={
    "dijetmass_study_20120913/studymjj_ak5_GluGluToHToZZTo2L2Q_M-300_8TeV-powheg-pythia6.root",
    "dijetmass_study_20120913/studymjj_ak5_GluGluToHToZZTo2L2Q_M-600_8TeV-powheg-pythia6.root",
    "dijetmass_study_20120913/studymjj_ak5_GluGluToHToZZTo2L2Q_M-900_8TeV-powheg-pythia6.root"
  };

  std::string mcName[NHISTOS]={
    "M_{H}=300 GeV",
    "M_{H}=600 GeV",
    "M_{H}=900 GeV"

  };


  for(int i=0; i<NHISTOS; i++){

    fmc[i] = TFile::Open(mcfile[i].data());
    h[i]   = (TProfile*)(fmc[i]->Get(var.data()));

  }


  int COLOR[NHISTOS]={4,1, 2};
  int MARKERSTYLE[NHISTOS]={24,21,29};

  for(int i=0; i < NHISTOS; i++){


//     h[i]->SetTitle("");
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



  vector<float> maxArray;
  maxArray.clear();

  vector<float> minArray;
  minArray.clear();

  for(int ih=0; ih < NHISTOS; ih++){

    h[ih]->GetXaxis()->SetRangeUser(xmin,xmax);
    float max_this  =  h[ih]->GetBinError(h[ih]->GetMaximumBin()) + h[ih]->GetMaximum();
    maxArray.push_back(max_this);

    float min_this  = -h[ih]->GetBinError(h[ih]->GetMinimumBin()) + h[ih]->GetMinimum();
    minArray.push_back(min_this);

  }


  float max = *max_element(maxArray.begin(),maxArray.end());
  cout << "Max = " << max << endl;

  float min = *min_element(minArray.begin(),minArray.end());
  cout << "Min = " << min << endl;

  bool isRatio=false;
  if(var.find("Rpt")!=std::string::npos || 
     var.find("Rm")!=std::string::npos)
    isRatio=true;

  for(int ih=0; ih < NHISTOS; ih++){

//     h[ih]->SetMaximum(1.2*max);
//     if(min>0)
//       h[ih]->SetMinimum(0.8*min);
//     else
//       h[ih]->SetMinimum(1.2*min);
    if(!isRatio)
      {
	h[ih]->SetMaximum(20);
	h[ih]->SetMinimum(-15);
      }
    else
      {
// 	h[ih]->SetMaximum(1.5);
// 	h[ih]->SetMinimum(0.5);
	h[ih]->SetMaximum(1.2);
	h[ih]->SetMinimum(0.8);
      }
  }

  cout << "here" << endl;
  TCanvas* c1 = new TCanvas("c1","",600,500);  
  gStyle->SetTitleX(0.1);
  h[0]->Draw("e");
  for(int ih=0; ih < NHISTOS-1; ih++)
    h[ih]->Draw("esame");

  h[NHISTOS-1]->Draw("hesame");


  float yLine = isRatio? 1.0: 0.0;
  TLine level(xmin,yLine,xmax,yLine);
  level.SetLineWidth(3);
  level.SetLineStyle(2);
  level.SetLineColor(kGray+2);
  level.Draw("same");
 
  float x1NDC = 0.192953;
  float y1NDC = 0.747881;
  float x2NDC = 0.412752;
  float y2NDC = 0.913136;

  if(var.find("dpt")!=std::string::npos)
    {
      x1NDC = 0.669355;
      y1NDC = 0.758475;
      x2NDC = 0.889113;
      y2NDC = 0.925847;
    }
  std::string headertitle = "H^{0}#rightarrow ZZ#rightarrow 2l2q";
  TLegend* leg = new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
//   leg->SetHeader(headertitle.data());
  for(int i=0; i < NHISTOS; i++)
    leg->AddEntry(h[i], mcName[i].data());
  leg->Draw("same");


  string dirName = "20120913_3profiles";
  gSystem->mkdir(dirName.data());

  std::string filename;
  std::string psname = dirName + "/ak5_" + var;
  
  filename = psname + ".eps";
  c1->Print(filename.data());
  filename = psname + ".gif";
  c1->Print(filename.data());
  filename = psname + ".pdf";
  c1->Print(filename.data());
  //   c1->Close();
}
		     
