#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"
#include <algorithm>

void forPaper_approval(std::string var1="h_ystar", 
			 bool logScale=false,
			 std::string datafile="darko_root/goldenWithMCerror_withJESCorr_bigmatrix_corr.root",
			 std::string mcfile1="darko_root/bare_exclusive1Jet_zPt40_both_dressed_DYToLL_M-50_1jEnh2_2jEnh35_3jEnh40_4jEnh50_7TeV-sherpa.root", 
			 std::string mcfile2="darko_root/bare_exclusive1Jet_zPt40_both_dressed_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola.root", 
			 std::string mcfile3="unified_angular_distributions/rebinnings/zpt40/Z_1jet_tota_cteq66_1___1___ex_m34.root",
			 std::string var2="", 
			 std::string headertitle="Z + 1 jet",
			 std::string dataName="CMS Data",
			 std::string mcName1="Sherpa",
			 std::string mcName2="Madgraph",
			 std::string mcName3="MCFM"
			 )
{
  
  setTDRStyle();
  gStyle->SetOptStat(0);

  const double LABELSIZE = 25.0;
  const int LINEWIDTH=3;
  const int NHISTOS=4;
  TH1F* h[NHISTOS];

  char tempName[300];
  if(var2 ==  "" )var2=var1;

  std::string xtitle;
  std::string output;
  std::string datavar;
  std::string var3;
  std::string theoryName;
  std::string indexName;

  std::string remword3  ="h_";
  std::string corrName = var1;
  size_t pos3  = corrName.find(remword3);
  if(pos3!= std::string::npos)
    corrName.replace(pos3,remword3.length(),"");

  datavar = "h_combine_" + corrName;

  double ymax=1.22;
  double ymin=0.78;

  double xmax=3.0;
  double xmin=0.0;

  if(var1=="h_ystar")
    {
      var3  = "id4";
      xtitle = "Y_{dif}";
      output = "DifYAll";
      theoryName = "Ydif";
      ymax = 1.45;
      ymin = 0.55; 
      xmax = 1.7999;
      indexName = "(d)";
    }
  else if(var1=="h_yB")
    {
      var3 = "id3";
      xtitle = "Y_{sum}";
      output = "SumYAll";
      theoryName = "Ysum";
      ymax = 1.45;
      ymin = 0.55;
      xmax = 2.1999;
      indexName = "(c)";
    }
  else if(var1=="h_jety")
    {
      var3 = "id1";
      xtitle = "|Y_{jet}|";
      output = "YJetAll";
      datafile = "darko_root/goldenWithMCerror_withJESCorr_jety_bigmatrix.root";
      theoryName = "Yjet";
      xmax = 2.3999;
      indexName = "(b)";
    }
  else if(var1=="h_zy")
    {
      var3 = "id2";
      xtitle = "|Y_{Z}|";
      output = "YZedAll";
      theoryName = "Yzed";
      xmax = 2.1999;
      indexName = "(a)";
    }

  

  // first get the histogram files
  TFile *fdata  = TFile::Open(datafile.data());
  cout << "Opening " << fdata->GetName() << endl;

  TFile *fmc1   = TFile::Open(mcfile1.data());
  cout << "Opening " << fmc1->GetName() << endl;

  TFile *fmc2   = TFile::Open(mcfile2.data());
  cout << "Opening " << fmc2->GetName() << endl;

  TFile *fmc3   = TFile::Open(mcfile3.data());
  cout << "Opening " << fmc3->GetName() << endl;

  h[0] = (TH1F*)(fdata->Get(datavar.data()));
  h[1] = (TH1F*)(fmc1->Get(var1.data()));
  h[2] = (TH1F*)(fmc2->Get(var2.data()));
  h[3] = (TH1F*)(fmc3->Get(var3.data()));

  if(var1=="h_jety")
    {
      double value = h[3]->GetBinContent(12);
      value *= 0.9;
      h[3]->SetBinContent(12,value);
    }


  TH1D* hscale[NHISTOS];

  int COLOR[NHISTOS]={1,4,2,kOrange-1};
  int MARKERSTYLE[NHISTOS]={8,24,21,29};
  int MARKERSIZE[NHISTOS]={1,0,0,0};
  int LINESTYLE[NHISTOS]={1,1,2,6};
  //  int FILLSTYLE[NHISTOS]={1,3345,3436,1};
  int FILLSTYLE[NHISTOS]={1,3345,3354,1};

  for(int i=0; i < NHISTOS; i++){

    hscale[i]   =(TH1D*) h[0]->Clone(Form("hscale%02i",i));
    hscale[i]   ->SetYTitle(Form("Ratio to %s",mcName3.data()));
    hscale[i]   ->SetXTitle(xtitle.data());
    hscale[i]   ->GetXaxis()->SetDecimals();
    hscale[i]   ->GetYaxis()->SetDecimals();

    hscale[i]->SetLineColor(COLOR[i]);
    hscale[i]->SetLineWidth(LINEWIDTH);
    hscale[i]->SetLineStyle(1);
    hscale[i]->SetMarkerColor(COLOR[i]);
    hscale[i]->SetMarkerStyle(MARKERSTYLE[i]);
    hscale[i]->SetMarkerSize(MARKERSIZE[i]);
    hscale[i]->SetFillColor(COLOR[i]);
    hscale[i]->SetFillStyle(FILLSTYLE[i]);

    hscale[i]->SetTitle("");
    hscale[i]->SetMaximum(ymax);
    hscale[i]->SetMinimum(ymin);

    hscale[i]->SetTitleOffset(1.2,"X");
    hscale[i]->SetTitleOffset(1.2,"Y");

    h[i]->SetTitle("");
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
  binLo = 1;
  binHi = nbins;

  double scaleFactor[NHISTOS]={1};

  for(int ih=0; ih < NHISTOS; ih++){
    
    double integral = h[ih]->Integral(binLo,binHi);
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

  vector<double> maxArray;
  maxArray.clear();
  
  for(int ih=0; ih < NHISTOS; ih++){

    h[ih]->GetXaxis()->SetRangeUser(xmin,xmax);
    hscale[ih]->GetXaxis()->SetRangeUser(xmin,xmax);
    double max_this  = h[ih]->GetBinError(h[ih]->GetMaximumBin()) + h[ih]->GetMaximum();
    maxArray.push_back(max_this);

  }


  double max = *(std::max_element(maxArray.begin(),maxArray.end()));
  cout << "Max = " << max << endl;

  for(int ih=0; ih < NHISTOS; ih++){

    h[ih]->SetMaximum(1.5*max);
    if(!logScale)
      h[ih]->SetMinimum(-0.015);
    else
      h[ih]->SetMinimum(5e-6);
  }


  //////////////////////////////////////////////////////////////////////////
  /// 
  ///   Making final figures and save the canvas in files
  ///
  //////////////////////////////////////////////////////////////////////////

  TCanvas* c1 = new TCanvas("c1","",700,1000);  
  c1->Divide(1,2,0.01,0);
  c1->cd(1);
  if(logScale)
    gPad->SetLogy(1);
  gPad->SetTopMargin(0.01);
  gPad->SetBottomMargin(0);
  gPad->SetRightMargin(0.04);
  gStyle->SetStatFontSize(0.05);

  // pad 1 and pad 2 have different sizes, need to rescale to 
  // make the label size the same
  double temp1_pad = gPad->GetWh()*gPad->GetAbsHNDC();
  cout << "pad 1 size in pixels = " << temp1_pad << endl;
  double label1_size = LABELSIZE/temp1_pad;
  for(int ih=0; ih < NHISTOS; ih++)
    {
      h[ih]->GetYaxis()->SetLabelSize(label1_size);
      h[ih]->GetXaxis()->SetLabelSize(label1_size);
    }

  h[0]->SetYTitle("1/#sigma d#sigma/dY");
  h[0]->Draw("9e1");
  for(int ih=1; ih < NHISTOS; ih++)
    h[ih]->Draw("9histsame");
  h[0]->Draw("9e1same");


  double x1NDC = 0.647708;
  double y1NDC = 0.512271;
  double x2NDC = 0.857529;
  double y2NDC = 0.823937;

  TLegend* leg = new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.055);
  leg->SetBorderSize(0);
  leg->SetHeader(headertitle.data());
  leg->AddEntry(h[0], dataName.data());
  leg->AddEntry(h[1], mcName1.data(),"l");
  leg->AddEntry(h[2], mcName2.data(),"l");
  leg->AddEntry(h[3], Form("%s (NLO)",mcName3.data()),"l");
  leg->Draw("same");

  TLatex *lar = new TLatex(0.35, 0.89, "CMS,   #sqrt{s} = 7 TeV, L_{int} = 5 fb^{-1}");
  lar->SetNDC(kTRUE); 
  lar->SetTextSize(0.06);
  lar->Draw();


  TLatex *larIndex = new TLatex(0.20, 0.89, indexName.data());
  larIndex->SetNDC(kTRUE);
  larIndex->SetTextSize(0.06);
  larIndex->Draw();

  TPaveText *pavetex = new TPaveText(0.17029, 
				     0.0495665, 
				     0.478939, 
				     0.5,"NDCBR");
  pavetex->SetBorderSize(0);
  pavetex->SetFillColor(0);
  pavetex->SetFillStyle(0);
  pavetex->SetLineWidth(3);
  pavetex->SetTextAlign(12);
  pavetex->SetTextSize(0.05);
  pavetex->AddText("76 < M_{ll} < 106 GeV");
  pavetex->AddText("p_{T}^{ll} > 40 GeV"); 
  pavetex->AddText("p_{T}^{l} > 20 GeV, |#eta^{l}| < 2.1");
  pavetex->AddText("N_{jet}=1, p_{T}^{jet} > 30 GeV, |#eta^{jet}| < 2.4");
  pavetex->AddText("#Delta R(l,jet)>0.5");
  if(var1=="h_zy"){
    pavetex->Draw();
  }
  cout << "pad 1 label size = " << h[0]->GetYaxis()->GetLabelSize() << endl;

  //////////////////////////////////////////////////////////////////////  

  c1->cd(2);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0);
  gPad->SetBottomMargin(0.2);
  gPad->SetTickx();

  // pad 1 and pad 2 have different sizes, need to rescale to 
  // make the label size the same
  double temp2_pad = gPad->GetWh()*gPad->GetAbsHNDC();
  cout << "pad 2 size in pixels = " << temp2_pad << endl;
  double label2_size = LABELSIZE/temp2_pad;
  for(int ih=0; ih < NHISTOS; ih++)
    {
      hscale[ih]->GetYaxis()->SetLabelSize(label2_size);
      hscale[ih]->GetXaxis()->SetLabelSize(label2_size);
    }


  hscale[0]->Draw("9e1");
  for(int ih=1; ih < NHISTOS-1; ih++){
    hscale[ih]->Draw("9e3same");
  }
  hscale[0]->Draw("9e1same");

  TLine* l2 = new TLine(xmin,1.,xmax,1.);
  l2->SetLineColor(kOrange-1);
  l2->SetLineStyle(1);
  l2->Draw("same");


  TLegend* leg2 = new TLegend(0.17333,0.244648,0.392274,0.476857);
  
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.05);
  leg2->SetBorderSize(0);
  leg2->AddEntry(hscale[1], "Stat. uncertainty on Sherpa","f");
  leg2->AddEntry(hscale[2], "Stat. uncertainty on Madgraph","f");
  leg2->Draw("same");


  gROOT->ProcessLine(".L /afs/cern.ch/user/s/syu/scripts/theoryErrorZed.c");
  theoryErrorZed(theoryName.data());

  hscale[0]->Draw("9e1same");
  for(int ih=1; ih < NHISTOS-1; ih++){
    hscale[ih]->Draw("9e3same");
  }
  hscale[0]->Draw("9e1same");
  l2->Draw("same");

  cout << "pad 2 label size = " << hscale[0]->GetYaxis()->GetLabelSize() << endl;
  
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
		     
