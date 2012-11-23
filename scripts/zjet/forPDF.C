#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"

void forPDF(std::string var="h_ystar", 
	    std::string mcfile1="PDF_study_electron.root", 
	    bool logScale=false,
	    std::string headertitle="Z#rightarrow ee + 1 jet",
	    std::string mcName1="MSTW2008 NLO",
	    std::string mcName2="CTEQ 66m",
	    std::string mcName3="CTEQ 6L1"
	    )
{
  
  setTDRStyle();
  gStyle->SetOptStat(0);

  const int LINEWIDTH=3;
  const int NHISTOS=3;
  TH1F* h[NHISTOS];

  char tempName[300];

  std::string var1 = var + "_mstw2008";
  std::string var2 = var + "_cteq66m";
  std::string var3 = var + "_cteq6l1";

  std::string xtitle;
  std::string output;
  std::string theoryName;

  std::string remword3  ="h_";
  std::string corrName = var;
  size_t pos3  = corrName.find(remword3);
  if(pos3!= std::string::npos)
    corrName.replace(pos3,remword3.length(),"");


  double ymax=1.22;
  double ymin=0.78;

  double xmax=3.0;
  double xmin=0.0;

  if(var=="h_ystar")
    {
//       xtitle = "0.5|Y_{Z}-Y_{jet}|";
      xtitle = "Y_{dif}";
      output = "DifYEle";
      theoryName = "Ydif";
      xmax = 1.7999;
    }
  else if(var=="h_yB")
    {
//       xtitle = "0.5|Y_{Z}+Y_{jet}|";
      xtitle = "Y_{sum}";
      output = "SumYEle";
      theoryName = "Ysum";
      xmax = 2.1999;
    }
  else if(var=="h_jety")
    {
      xtitle = "|Y_{jet}|";
      output = "YJetEle";
      theoryName = "Yjet";
      xmax = 2.3999;

    }
  else if(var=="h_zy")
    {
      xtitle = "|Y_{Z}|";
      output = "YZedEle";
      theoryName = "Yzed";
      xmax = 2.1999;
    }

  

  // first get the histogram files
  TFile *fmc1   = TFile::Open(mcfile1.data());

  cout << "Opening " << var1 << endl;
  cout << "Opening " << var2 << endl;
  cout << "Opening " << var3 << endl;

  h[0]    = (TH1F*)(fmc1->Get(var1.data()));
  h[1]    = (TH1F*)(fmc1->Get(var2.data()));
  h[2]    = (TH1F*)(fmc1->Get(var3.data()));

  TH1D* hscale[NHISTOS];

  int COLOR[NHISTOS]={1,4,2};
  int MARKERSTYLE[NHISTOS]={8,24,21};
  int LINESTYLE[NHISTOS]={2,4,1};

  for(int i=0; i < NHISTOS; i++){

    hscale[i]   =(TH1D*) h[0]->Clone(Form("hscale%02i",i));
    hscale[i]   ->SetYTitle(Form("Ratio to %s",mcName3.data()));
    hscale[i]   ->SetXTitle(xtitle.data());
//     hscale[i]   ->GetXaxis()->SetNdivisions(510);
//     hscale[i]->GetYaxis()->SetNdivisions(5);
    hscale[i]   ->GetXaxis()->SetDecimals();
    hscale[i]   ->GetYaxis()->SetDecimals();

    hscale[i]->SetLineColor(COLOR[i]);
    hscale[i]->SetLineWidth(LINEWIDTH);
    hscale[i]->SetLineStyle(LINESTYLE[i]);
    hscale[i]->SetMarkerColor(COLOR[i]);
    hscale[i]->SetMarkerStyle(MARKERSTYLE[i]);
    hscale[i]->SetMarkerSize(1);

    hscale[i]->SetTitle("");
    hscale[i]->SetMaximum(ymax);
    hscale[i]->SetMinimum(ymin);

    hscale[i]->SetTitleOffset(1.2,"X");
    hscale[i]->SetTitleOffset(1.2,"Y");

    h[i]->SetTitle("");
//     h[i]->GetXaxis()->SetNdivisions(510);
//     h[i]->GetYaxis()->SetNdivisions(5);
    h[i]->SetLineStyle(LINESTYLE[i]);
    h[i]->GetXaxis()->SetDecimals();
    h[i]->GetYaxis()->SetDecimals();
    h[i]->SetMarkerSize(1);
    h[i]->SetLineColor(COLOR[i]);
    h[i]->SetLineWidth(LINEWIDTH);
    h[i]->SetMarkerColor(COLOR[i]);
    h[i]->SetMarkerStyle(MARKERSTYLE[i]);
    h[i]->SetTitleOffset(1.2,"Y");

  }
  h[0]->SetLineWidth(1);
  hscale[0]->SetLineWidth(1);

  // if normalizing to the same area, set the scale 

  int binLo = -1;
  int binHi = -1;
  int nbins = h[0]->GetNbinsX();
//   if(xmin>-9999.0 && xmax>-9999.0)
//     {

//       binLo = h[0]->FindBin(xmin);
//       binHi = h[0]->FindBin(xmax)-1;

//     }

//   else
//     {
      binLo = 1;
      binHi = nbins;
//       xmin = h[0]->GetBinLowEdge(1);
//       xmax = h[0]->GetBinLowEdge(nbins+1);
//     }


  float scaleFactor[NHISTOS]={1};

  for(int ih=0; ih < NHISTOS; ih++){
    
    float integral = h[ih]->Integral(binLo,binHi);
    scaleFactor[ih] = integral > 0? 1.0/integral: 1;
    h[ih]->Sumw2();
    h[ih]->Scale(scaleFactor[ih]);

  }

  for(int ih=0; ih < NHISTOS; ih++)
    cout << "histogram " << ih << " integral = " << h[ih]->Integral() << endl;


  // get the ratio

  for(int ih=0; ih < NHISTOS-1; ih++){
    cout << "===================================================" << endl;
    cout << "For histogram " << ih << endl;
    hscale[ih]->Divide(h[ih], h[NHISTOS-1]);
    hscale[ih]->SetMaximum(ymax);
    hscale[ih]->SetMinimum(ymin);

  } // end of loop over histograms

  vector<float> maxArray;
  maxArray.clear();

  for(int ih=0; ih < NHISTOS; ih++){

    h[ih]->GetXaxis()->SetRangeUser(xmin,xmax);
    hscale[ih]->GetXaxis()->SetRangeUser(xmin,xmax);
    float max_this  = h[ih]->GetBinError(h[ih]->GetMaximumBin()) + h[ih]->GetMaximum();
    maxArray.push_back(max_this);

  }


  float max = *max_element(maxArray.begin(),maxArray.end());
  cout << "Max = " << max << endl;

  for(int ih=0; ih < NHISTOS; ih++){

    h[ih]->SetMaximum(1.1*max);
    if(!logScale)
      h[ih]->SetMinimum(-0.015);
    else
      h[ih]->SetMinimum(5e-6);
  }

  cout << "here" << endl;
  TCanvas* c1 = new TCanvas("c1","",700,1000);  
  c1->Divide(1,2,0.01,0);
  c1->cd(1);
  if(logScale)
    gPad->SetLogy(1);
  gPad->SetTopMargin(0.01);
  gPad->SetBottomMargin(0);
  gPad->SetRightMargin(0.04);

  
  cout << "here1" << endl;
//   h[0]->SetYTitle("Arbitrary Unit");
  h[0]->SetYTitle("1/#sigma d#sigma/dY");
  h[0]->Draw("chist");
  for(int ih=1; ih < NHISTOS; ih++)
    h[ih]->Draw("chistsame");
  h[0]->Draw("chistsame");
  cout << "here2" << endl;


  cout << "here3" << endl;

  float x1NDC = 0.609697;
  float y1NDC = 0.636937;
  float x2NDC = 0.819518;
  float y2NDC = 0.948603;

  TLegend* leg = new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  leg->SetBorderSize(0);
  leg->SetHeader(headertitle.data());
  leg->AddEntry(h[0], mcName1.data(),"l");
  leg->AddEntry(h[1], mcName2.data(),"l");
  leg->AddEntry(h[2], mcName3.data(),"l");
  leg->Draw("same");


  c1->cd(2);
  gStyle->SetStatW       (0.3);
  gStyle->SetStatH       (0.3);
  gStyle->SetStatX       (0.879447);
  gStyle->SetStatY       (0.939033);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatBorderSize(0);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0);
  gPad->SetBottomMargin(0.2);
  gPad->SetTickx();
  gStyle->SetOptFit(1);

  cout << "here4" << endl;

  hscale[0]->Draw("chist");
  for(int ih=0; ih < NHISTOS-1; ih++){
    hscale[ih]->Draw("chistsame");
  }
  hscale[0]->Draw("chistsame");

  cout << "here5" << endl;
  TLine* l2 = new TLine(xmin,1.,xmax,1.);
  l2->SetLineColor(2);
  l2->SetLineStyle(1);
  l2->Draw("same");

  
  string dirName = "forPDF";
  gSystem->mkdir(dirName.data());

  std::string filename;
  std::string psname = dirName + "/" + output;
  filename = psname + ".eps";
  c1->Print(filename.data());
  filename = psname + ".gif";
  c1->Print(filename.data());
  filename = psname + ".pdf";
  c1->Print(filename.data());
  //   c1->Close();
}
		     
