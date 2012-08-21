#include "setTDRStyle.C"

void compare3(
		    std::string var="h_mh_parton", 
		    bool logScale=false,
		    float xmin=-9999.0, float xmax=-9999.0,
		    std::string output="test"
		    )
{
  
  setTDRStyle();
  gStyle->SetOptStat(0);

  const int NHISTOS=3;
  TH1F* h[NHISTOS];

  char tempName[300];
  TFile *fmc[NHISTOS];
  
  std::string mcfile[NHISTOS]={
    "biggerBin_best/mjj_hzz2l2q_M300.root",
    "biggerBin_best/mjj_hzz2l2q_M600.root",
    "biggerBin_best/mjj_hzz2l2q_M900.root"
  };

  std::string mcName[NHISTOS]={
    "M_{H}=300 GeV/c^{2}",
    "M_{H}=600 GeV/c^{2}",
    "M_{H}=900 GeV/c^{2}"

  };


  for(int i=0; i<NHISTOS; i++){

    fmc[i] = TFile::Open(mcfile[i].data());
    h[i]   = (TH1F*)(fmc[i]->Get(var.data()));

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
  TCanvas* c1 = new TCanvas("c1","",500,500);  
  if(logScale)
    gPad->SetLogy(1);

  h[0]->Draw("e");
  for(int ih=0; ih < NHISTOS-1; ih++)
    h[ih]->Draw("esame");

  h[NHISTOS-1]->Draw("hesame");

  float x1NDC = 0.620806;
  float y1NDC = 0.75136;
  float x2NDC = 0.840565;
  float y2NDC = 0.917373;

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


  string dirName = "20120629_3histograms";
  gSystem->mkdir(dirName.data());

  std::string filename;
  std::string psname = dirName + "/" + var;
  if(output !="test")
    psname = dirName+ "/" + output;
  else
    psname = dirName+ "/" + var;
  filename = psname + ".eps";
  c1->Print(filename.data());
  filename = psname + ".gif";
  c1->Print(filename.data());
  filename = psname + ".pdf";
  c1->Print(filename.data());
  //   c1->Close();
}
		     