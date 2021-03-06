#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"

void compareGen(std::string mcfile1, std::string mcfile2, 
		std::string var1, std::string var2="",
		std::string xtitle, 
		std::string output="test",
		std::string header="Z(#rightarrow ee)+#geq 1 jet",
		std::string mcName1="Data",
		std::string mcName2="Madgraph",
		float xmin=-9999.0, float xmax=-9999.0,
		bool logScale=true
		)
{
  
  setTDRStyle();
  gStyle->SetOptStat(0);

  TH1F* h1;
  TH1F* h2;

  char tempName[300];
  if(var2 ==  "" )var2=var1;
  
  // first get the histogram files
  TFile *fmc1 = TFile::Open(mcfile1.data());
  TFile *fmc2   = TFile::Open(mcfile2.data());

  h1  = (TH1F*)(fmc1->Get(var1.data()));
  h2    = (TH1F*)(fmc2->Get(var2.data()));

  TH1D* hscale =(TH1D*) h1->Clone("hscale");
  hscale->SetYTitle(Form("%s/%s",mcName1.data(),mcName2.data()));

  h1->GetXaxis()->SetNdivisions(5);
  h1->GetYaxis()->SetDecimals();

  h2->GetXaxis()->SetNdivisions(5);
  h2->GetYaxis()->SetDecimals();

  hscale->GetXaxis()->SetNdivisions(5);
  hscale->GetYaxis()->SetDecimals();


  h1->SetLineColor(2);
  h1->SetMarkerColor(2);
  h1->SetMarkerSize(1);
  h1->SetMarkerStyle(24);


  h2->SetLineColor(4);
  h2->SetMarkerColor(4);
  h2->SetMarkerSize(1);
  h2->SetMarkerStyle(21);

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


  float scale_mc = (float)h1->Integral(binLo,binHi)/(float)h2->Integral(binLo,binHi);
//   cout << "binLo = " << binLo << ", binHi = " << binHi << endl;
//    cout << "xmin = " << xmin << "xmax = " << xmax << endl;

//   h2->Sumw2();
// //   scale_mc = 1000.0*4.890*3048.0/2.29809910000000000e+07;
//   h2->Scale(scale_mc);

  cout << "h2 integral = " << h2->Integral() << endl;
  cout << "h1 integral = "   << h1->Integral() << endl;;

  // get the ratio
  double chi2 = 0;
  int realbin = 0;
  for(int i=1;i<= nbins;i++){

    double nmc=h2->GetBinContent(i);
    double ndata=h1->GetBinContent(i); 
    double nmcerr=h2->GetBinError(i);
    double ndataerr=h1->GetBinError(i); 
    
    if(nmc<0 || ndata<0)continue;    
    
    if(nmcerr==0 && ndataerr==0)continue;

    if(nmc==0 && ndata==0)continue;

    double chi2ndef = (nmc-ndata)*(nmc-ndata)/
      ( nmcerr*nmcerr+ ndataerr*ndataerr);
    chi2 += chi2ndef;
    realbin++;

    cout << "Bin " << i << " : " << ndata << ", " << nmc;
    cout << " " << chi2ndef << endl;


    // now calculate the ratio
    if(nmc==0 || nmcerr==0 || ndata==0 || ndataerr==0)
      {
	hscale->SetBinContent(i,-9999);
	hscale->SetBinError(i,1e-4);
	continue;
      }

    cout << "Bin " << i << " ratio = " << ndata/nmc << endl;
    hscale->SetBinContent(i,ndata/nmc);
    double err = 0;
    err=
      (ndata/nmc)*sqrt(pow(nmcerr/nmc,2)+pow(ndataerr/ndata,2));
    hscale->SetBinError(i,err);

  }

  for(int i=1;i<=hscale->GetNbinsX();i++)
    cout << i << ": " << hscale->GetBinContent(i) << endl;

  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h2->GetXaxis()->SetRangeUser(xmin,xmax);
  hscale->GetXaxis()->SetRangeUser(xmin,xmax);


  TCanvas* c1 = new TCanvas("c1","",700,1000);  
  c1->Divide(1,2,0.01,0);
  c1->cd(1);
  if(logScale)
    gPad->SetLogy(1);
  gPad->SetTopMargin(0.01);
  gPad->SetBottomMargin(0);
  gPad->SetRightMargin(0.04);

  
  float max_data  = h1->GetBinError(h1->GetMaximumBin()) + h1->GetMaximum();
  float max_mc    = h2->GetBinError(h2->GetMaximumBin()) + h2->GetMaximum();

  if(max_data > max_mc)
    {
      h1->Draw("e");
      h2->Draw("hesame");
    }
  else
    { h2->Draw("he");
      h1->Draw("esame");
    }


  float x1NDC = 0.725;
  float y1NDC = 0.615;
  float x2NDC = 0.928;
  float y2NDC = 0.951;

  TLegend* leg = new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  
  leg->SetHeader(header.data());
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->AddEntry(h1, mcName1.data());
  leg->AddEntry(h2, mcName2.data());
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
  hscale->SetTitle("");
  hscale->GetXaxis()->SetTitle(xtitle.data());
//   hscale->SetMaximum(3.0);
//   hscale->SetMinimum(0.1);
  hscale->SetMaximum(2.0);
  hscale->SetMinimum(0.5);
  hscale->SetTitleOffset(1.2,"Y");
  hscale->Draw("e1");
  TF1* fline = new TF1("fline","pol1");
  TLine* l2 = new TLine(xmin,1.,xmax,1.);
  l2->SetLineColor(4);
  l2->SetLineStyle(3);
  fline->SetLineWidth(3);
  fline->SetLineColor(6);
  fline->SetNpx(2500);
//   hscale->Fit("fline","","");
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
		     
