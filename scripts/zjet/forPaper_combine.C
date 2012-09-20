#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"

void forPaper_combine(std::string var1="h_ystar", 
		      float xmin=-9999.0, float xmax=-9999.0,
		      bool logScale=false,
		      std::string datafile="darko_root/goldenWithMCerror_bigmatrix.root",
		      std::string mcfile1="darko_root/bare_exclusive1Jet_zPt40_both_dressed_DYToLL_M-50_1jEnh2_2jEnh35_3jEnh40_4jEnh50_7TeV-sherpa.root", 
		      std::string mcfile2="darko_root/bare_exclusive1Jet_zPt40_both_dressed_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola.root", 
		      std::string mcfile3="unified_angular_distributions/rebinnings/zpt40/Z_1jet_tota_cteq66_1___1___ex_m34.root",
		 std::string var2="", 
		 std::string headertitle="Z(#rightarrow ee, #mu#mu)+1 jet",
		 std::string dataName="Data",
		 std::string mcName1="Sherpa",
		 std::string mcName2="Madgraph",
		 std::string mcName3="MCFM"
		 )
{
  
  setTDRStyle();
  gStyle->SetOptStat(0);

  const int NHISTOS=4;
  TH1F* h[NHISTOS];

  char tempName[300];
  if(var2 ==  "" )var2=var1;

  std::string xtitle;
  std::string output;
  std::string datavar;
  std::string var3;
  std::string theoryName;

  std::string remword3  ="h_";
  std::string corrName = var1;
  size_t pos3  = corrName.find(remword3);
  if(pos3!= std::string::npos)
    corrName.replace(pos3,remword3.length(),"");

  datavar = "h_combine_" + corrName;

  if(var1=="h_ystar")
    {
      var3  = "id4";
      xtitle = "0.5|Y_{Z}-Y_{jet}|";
      output = "DifYAll";
      theoryName = "Ydif";
    }
  else if(var1=="h_yB")
    {
      var3 = "id3";
      xtitle = "0.5|Y_{Z}+Y_{jet}|";
      output = "SumYAll";
      theoryName = "Ysum";
    }
  else if(var1=="h_jety")
    {
      var3 = "id1";
      xtitle = "|Y_{jet}|";
      output = "YJetAll";
      datafile = "darko_root/goldenWithMCerror_jety_bigmatrix.root";
      theoryName = "Yjet";
    }
  else if(var1=="h_zy")
    {
      var3 = "id2";
      xtitle = "|Y_{Z}|";
      output = "YZedAll";
      theoryName = "Yzed";
    }

  

  // first get the histogram files
  TFile *fdata  = TFile::Open(datafile.data());
  TFile *fmc1   = TFile::Open(mcfile1.data());
  TFile *fmc2   = TFile::Open(mcfile2.data());
  TFile *fmc3   = TFile::Open(mcfile3.data());

  cout << "Opening " << fdata->GetName() << endl;
  cout << "Opening " << fmc1->GetName() << endl;
  cout << "Opening " << fmc2->GetName() << endl;
  cout << "Opening " << fmc3->GetName() << endl;

  h[0] = (TH1F*)(fdata->Get(datavar.data()));
  h[1]    = (TH1F*)(fmc1->Get(var1.data()));
  h[2]    = (TH1F*)(fmc2->Get(var2.data()));
  h[3]    = (TH1F*)(fmc3->Get(var3.data()));

  TH1D* hscale[NHISTOS];

  int COLOR[NHISTOS]={1,4,2,kOrange-1};
  int MARKERSTYLE[NHISTOS]={8,24,21,29};

  for(int i=0; i < NHISTOS; i++){

    hscale[i]   =(TH1D*) h[0]->Clone(Form("hscale%02i",i));
    hscale[i]   ->SetYTitle(Form("Ratio to %s",mcName3.data()));
    hscale[i]   ->SetXTitle(xtitle.data());
//     hscale[i]   ->GetXaxis()->SetNdivisions(5);
    hscale[i]   ->GetXaxis()->SetDecimals();
    hscale[i]   ->GetYaxis()->SetDecimals();

    hscale[i]->SetLineColor(COLOR[i]);
    hscale[i]->SetMarkerColor(COLOR[i]);
    hscale[i]->SetMarkerStyle(MARKERSTYLE[i]);

    hscale[i]->SetTitle("");
    hscale[i]->SetMaximum(2.05);
    hscale[i]->SetMinimum(0.0);

//     hscale[i]->SetMaximum(1.45);
//     hscale[i]->SetMinimum(0.4);
    hscale[i]->SetTitleOffset(1.2,"X");
    hscale[i]->SetTitleOffset(1.2,"Y");

    h[i]->SetTitle("");
//     h[i]->GetXaxis()->SetNdivisions(5);
    h[i]->GetXaxis()->SetDecimals();
    h[i]->GetYaxis()->SetDecimals();

    h[i]->SetLineColor(COLOR[i]);
    h[i]->SetMarkerColor(COLOR[i]);
    h[i]->SetMarkerSize(1);
    h[i]->SetMarkerStyle(MARKERSTYLE[i]);
    h[i]->SetTitleOffset(1.2,"Y");

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
//     scaleFactor[ih] = integral > 0? (float)h[0]->Integral(binLo,binHi)/integral: 1;
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
    hscale[ih]->SetMaximum(2.05);
    hscale[ih]->SetMinimum(0.0);

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
  h[0]->SetYTitle("Arbitrary Unit");
  h[0]->Draw("e");
  for(int ih=0; ih < NHISTOS-1; ih++)
    h[ih]->Draw("esame");

  cout << "here2" << endl;

  h[NHISTOS-1]->Draw("hesame");

  cout << "here3" << endl;

  float x1NDC = 0.682679;
  float y1NDC = 0.560219;
  float x2NDC = 0.892501;
  float y2NDC = 0.871885;

  TLegend* leg = new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  leg->SetBorderSize(0);
  leg->SetHeader(headertitle.data());
  leg->AddEntry(h[0], dataName.data());
  leg->AddEntry(h[1], mcName1.data());
  leg->AddEntry(h[2], mcName2.data());
  leg->AddEntry(h[3], mcName3.data());
  leg->Draw("same");

  TLatex *lar = new TLatex(0.48, 0.91, "CMS #sqrt{s} = 7 TeV, L_{int} = 4.9~5.1 fb^{-1}");
  lar->SetNDC(kTRUE);
  lar->SetTextSize(0.05);
  lar->Draw();

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

  hscale[0]->Draw("e1");
  for(int ih=0; ih < NHISTOS-1; ih++){
    hscale[ih]->Draw("e1same");
  }

  cout << "here5" << endl;
  TLine* l2 = new TLine(xmin,1.,xmax,1.);
  l2->SetLineColor(4);
  l2->SetLineStyle(3);
  l2->Draw("same");

  gROOT->ProcessLine(".L ~/scripts/theoryError.c");
  theoryError(theoryName.data());
  
  string dirName = "forPaper";
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
		     
