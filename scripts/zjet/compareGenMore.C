#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"

void compareGenMore(std::string mcfile1, std::string mcfile2, std::string mcfile3,
		    std::string var1, std::string var2="", std::string var3="",
		    std::string output="test",
		    std::string headertitle="Z(#rightarrow ee)+1 jet",
		    std::string mcName1="Sherpa",
		    std::string mcName2="Madgraph",
		    std::string mcName3="MCFM",
		    float xmin=-9999.0, float xmax=-9999.0,
		    bool logScale=false)
{
  
  setTDRStyle();
  gStyle->SetOptStat(0);

  TH1F* h1;
  TH1F* h2;
  TH1F* h3;

  char tempName[300];
  if(var2 ==  "" )var2=var1;
  if(var3 ==  "" )var3=var1;
  
  // first get the histogram files
  TFile *fmc1 = TFile::Open(mcfile1.data());
  TFile *fmc2   = TFile::Open(mcfile2.data());
  TFile *fmc3   = TFile::Open(mcfile3.data());

  h1  = (TH1F*)(fmc1->Get(var1.data()));
  h2    = (TH1F*)(fmc2->Get(var2.data()));
  h3    = (TH1F*)(fmc3->Get(var3.data()));

  TH1D* hscale1_3 =(TH1D*) h1->Clone("hscale1_3");
  hscale1_3->SetYTitle(Form("Ratio to %s",mcName3.data()));

  TH1D* hscale2_3 =(TH1D*) h1->Clone("hscale2_3");
  hscale2_3->SetYTitle(Form("Ratio to %s",mcName3.data()));

  h1->GetXaxis()->SetNdivisions(5);
  h1->GetYaxis()->SetDecimals();

  h2->GetXaxis()->SetNdivisions(5);
  h2->GetYaxis()->SetDecimals();

  h3->GetXaxis()->SetNdivisions(5);
  h3->GetYaxis()->SetDecimals();

  hscale1_3->GetXaxis()->SetNdivisions(5);
  hscale1_3->GetYaxis()->SetDecimals();

  hscale2_3->GetXaxis()->SetNdivisions(5);
  hscale2_3->GetYaxis()->SetDecimals();


  int COLOR[3]={4,2,kOrange-1};
  int MARKERSTYLE[3]={24,21,29};

  h1->SetLineColor(COLOR[0]);
  h1->SetMarkerColor(COLOR[0]);
  h1->SetMarkerSize(1);
  h1->SetMarkerStyle(MARKERSTYLE[0]);


  h2->SetLineColor(COLOR[1]);
  h2->SetMarkerColor(COLOR[1]);
  h2->SetMarkerSize(1);
  h2->SetMarkerStyle(MARKERSTYLE[1]);

  h3->SetLineColor(COLOR[2]);
  h3->SetMarkerColor(COLOR[2]);
  h3->SetMarkerSize(1);
  h3->SetMarkerStyle(MARKERSTYLE[2]);

  // if normalizing to the same area, set the scale 

  int binLo = -1;
  int binHi = -1;
  int nbins = h1->GetNbinsX();
  if(xmin>-9999.0 && xmax>-9999.0)
    {

      binLo = h1->FindBin(xmin);
      binHi = h1->FindBin(xmax)-1;

    }

  else
    {
      binLo = 1;
      binHi = nbins;
      xmin = h1->GetBinLowEdge(1);
      xmax = h1->GetBinLowEdge(nbins+1);
    }


//   cout << "binLo = " << binLo << ", binHi = " << binHi << endl;
//    cout << "xmin = " << xmin << "xmax = " << xmax << endl;

  float scale_mc2= (float)h1->Integral(binLo,binHi)/(float)h2->Integral(binLo,binHi);

  h2->Sumw2();
  h2->Scale(scale_mc2);

  float scale_mc3= (float)h1->Integral(binLo,binHi)/(float)h3->Integral(binLo,binHi);

  h3->Sumw2();
  h3->Scale(scale_mc3);

  cout << "h1 integral = "   << h1->Integral() << endl;;
  cout << "h2 integral = "   << h2->Integral() << endl;
  cout << "h3 integral = "   << h3->Integral() << endl;

  // get the ratio
  for(int i=1;i<= nbins;i++){

    double nref=h3->GetBinContent(i);
    double nreferr=h3->GetBinError(i);

    double nmc1=h1->GetBinContent(i); 
    double nmc1err=h1->GetBinError(i); 

    double nmc2=h2->GetBinContent(i); 
    double nmc2err=h2->GetBinError(i); 

    
    if(nref<=0 || nmc1<=0 || nmc2<=0)continue;    
    
    if(nreferr<=0 || nmc1err<=0 || nmc2err<=0)
      continue;

    // now calculate the ratio 1
    cout << "Bin " << i << " ratio1 = " << nmc1/nref << "\t";
    double err1_3 = 0;
    err1_3 =
      (nmc1/nref)*sqrt(pow(nreferr/nref,2)+pow(nmc1err/nmc1,2));
    hscale1_3->SetBinContent(i,nmc1/nref);
    hscale1_3->SetBinError(i,err1_3);

    // now calculate the ratio 2
    cout << nmc2/nref << endl;
    double err2_3 = 0;
    err2_3 =
      (nmc2/nref)*sqrt(pow(nreferr/nref,2)+pow(nmc2err/nmc2,2));
    hscale2_3->SetBinContent(i,nmc2/nref);
    hscale2_3->SetBinError(i,err2_3);

  }

  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h2->GetXaxis()->SetRangeUser(xmin,xmax);
  h3->GetXaxis()->SetRangeUser(xmin,xmax);
  hscale1_3->GetXaxis()->SetRangeUser(xmin,xmax);
  hscale2_3->GetXaxis()->SetRangeUser(xmin,xmax);


  TCanvas* c1 = new TCanvas("c1","",700,1000);  
  c1->Divide(1,2,0.01,0);
  c1->cd(1);
  if(logScale)
    gPad->SetLogy(1);
  gPad->SetTopMargin(0.01);
  gPad->SetBottomMargin(0);
  gPad->SetRightMargin(0.04);

  
  float max_mc1   = h1->GetBinError(h1->GetMaximumBin()) + h1->GetMaximum();
  float max_mc2   = h2->GetBinError(h2->GetMaximumBin()) + h2->GetMaximum();
  float max_mc3   = h3->GetBinError(h3->GetMaximumBin()) + h3->GetMaximum();


  vector<float> maxArray;
  maxArray.push_back(max_mc1);
  maxArray.push_back(max_mc2);
  maxArray.push_back(max_mc3);

  float max = *max_element(maxArray.begin(),maxArray.end());
  cout << "Max = " << max << endl;

  h1->SetMaximum(1.1*max);
  h2->SetMaximum(1.1*max);
  h3->SetMaximum(1.1*max);

  h1->Draw("e");
  h2->Draw("esame");
  h3->Draw("hesame");

  float x1NDC = 0.735894;
  float y1NDC = 0.620155;
  float x2NDC = 0.956358;
  float y2NDC = 0.934218;

  TLegend* leg = new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->SetHeader(headertitle.data());
  leg->AddEntry(h1, mcName1.data());
  leg->AddEntry(h2, mcName2.data());
  leg->AddEntry(h3, mcName3.data());
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
  hscale1_3->SetTitle("");
  hscale1_3->SetMaximum(3.0);
  hscale1_3->SetMinimum(0.1);
  hscale1_3->SetTitleOffset(1.2,"Y");
  hscale1_3->SetLineColor(COLOR[0]);
  hscale1_3->SetMarkerColor(COLOR[0]);
  hscale1_3->SetMarkerStyle(MARKERSTYLE[0]);
  hscale1_3->Draw("e1");

  hscale2_3->SetTitle("");
  hscale2_3->SetMaximum(3.0);
  hscale2_3->SetMinimum(0.1);
  hscale2_3->SetTitleOffset(1.2,"Y");
  hscale2_3->SetLineColor(COLOR[1]);
  hscale2_3->SetMarkerColor(COLOR[1]);
  hscale2_3->SetMarkerStyle(MARKERSTYLE[1]);
  hscale2_3->Draw("e1same");


  TLine* l2 = new TLine(xmin,1.,xmax,1.);
  l2->SetLineColor(4);
  l2->SetLineStyle(3);
  l2->Draw("same");


  string dirName = "compareGen";
  gSystem->mkdir(dirName.data());

  std::string filename;
  std::string psname = dirName + "/" + var1;
  if(output !="test")
    psname = dirName+ "/" + output;
  else
    psname = dirName+ "/" + var1;
  filename = psname + ".eps";
  c1->Print(filename.data());
  filename = psname + ".gif";
  c1->Print(filename.data());
  filename = psname + ".pdf";
  c1->Print(filename.data());
  //   c1->Close();
}
		     
