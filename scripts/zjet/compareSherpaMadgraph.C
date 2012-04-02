#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"

void compareSherpaMadgraph(std::string sherpafile, std::string madgraphfile, 
			   std::string var,
			   float xmin=-9999.0, float xmax=-9999.0,
			   bool logScale=false, 
			   std::string output="test")
{
  
  setTDRStyle();
  gStyle->SetOptStat(0);

  TH1F* hsherpa;
  TH1F* hmadgraph;

  char tempName[300];


  // first get the histogram files
  TFile *fsherpa = TFile::Open(sherpafile.data());
  TFile *fmadgraph   = TFile::Open(madgraphfile.data());

  hsherpa  = (TH1F*)(fsherpa->Get(var.data()));
  hmadgraph    = (TH1F*)(fmadgraph->Get(var.data()));

  TH1D* hscale =(TH1D*) hsherpa->Clone("hscale");
  hscale->SetYTitle("SHERPA/MADGRAPH");

  int nREBIN=2;
  if(var.find("pt")!= std::string::npos)
    nREBIN=4;

  hsherpa->GetXaxis()->SetNdivisions(5);
  hsherpa->GetYaxis()->SetDecimals();
  hsherpa->Rebin(nREBIN);

  hmadgraph->GetXaxis()->SetNdivisions(5);
  hmadgraph->GetYaxis()->SetDecimals();
  hmadgraph->Rebin(nREBIN);

  hscale->GetXaxis()->SetNdivisions(5);
  hscale->GetYaxis()->SetDecimals();
  hscale->Rebin(nREBIN);

  hsherpa->SetLineColor(2);
  hsherpa->SetMarkerColor(2);
  hsherpa->SetMarkerSize(1);
  hsherpa->SetMarkerStyle(24);


  hmadgraph->SetLineColor(4);
  hmadgraph->SetMarkerColor(4);
  hmadgraph->SetMarkerSize(1);
  hmadgraph->SetMarkerStyle(21);

  // if normalizing to the same area, set the scale 

  int binLo = -1;
  int binHi = -1;
  int nbins = hsherpa->GetNbinsX();
  if(xmin>-9999.0 && xmax>-9999.0)
    {

      binLo = hsherpa->FindBin(xmin);
      binHi = hsherpa->FindBin(xmax)-1;

    }

  else
    {
      binLo = 1;
      binHi = nbins;
      xmin = hsherpa->GetBinLowEdge(1);
      xmax = hsherpa->GetBinLowEdge(nbins+1);
    }


   float scale_mc = (float)hsherpa->Integral(binLo,binHi)/(float)hmadgraph->Integral(binLo,binHi);
   cout << "binLo = " << binLo << ", binHi = " << binHi << endl;
   cout << "xmin = " << xmin << "xmax = " << xmax << endl;

//    hmadgraph->Sumw2();
//    hmadgraph->Scale(scale_mc);

   float scale_sherpa   = 1000.0*4.890*3048.0/30607750;
   float scale_madgraph = 1000.0*4.890*3048.0/36209629;

   hsherpa->Sumw2();
   hsherpa->Scale(scale_sherpa);

   hmadgraph->Sumw2();
   hmadgraph->Scale(scale_madgraph);

   cout << "hmadgraph integral = " << hmadgraph->Integral() << endl;
   cout << "hsherpa integral = "   << hsherpa->Integral() << endl;;

  // get the ratio
  double chi2 = 0;
  int realbin = 0;
  for(int i=1;i<= nbins;i++){

    double nmc=hmadgraph->GetBinContent(i);
    double ndata=hsherpa->GetBinContent(i); 
    double nmcerr=hmadgraph->GetBinError(i);
    double ndataerr=hsherpa->GetBinError(i); 
    
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
    if(nmc==0 || nmcerr==0 || ndata==0 || ndataerr==0)continue;
    cout << "Bin " << i << " ratio = " << ndata/nmc << endl;
    hscale->SetBinContent(i,ndata/nmc);
    double err = 0;
    err=
      (ndata/nmc)*sqrt(pow(nmcerr/nmc,2)+pow(ndataerr/ndata,2));
    hscale->SetBinError(i,err);

  }

  
  hsherpa->GetXaxis()->SetRangeUser(xmin,xmax);
  hmadgraph->GetXaxis()->SetRangeUser(xmin,xmax);
  hscale->GetXaxis()->SetRangeUser(xmin,xmax);


  TCanvas* c1 = new TCanvas("c1","",700,1000);  
  c1->Divide(1,2,0.01,0);
  c1->cd(1);
  if(logScale)
    gPad->SetLogy(1);
  gPad->SetTopMargin(0.01);
  gPad->SetBottomMargin(0);
  gPad->SetRightMargin(0.04);

  
  float max_data  = hsherpa->GetBinError(hsherpa->GetMaximumBin()) + hsherpa->GetMaximum();
  float max_mc    = hmadgraph->GetBinError(hmadgraph->GetMaximumBin()) + hmadgraph->GetMaximum();

  if(max_data > max_mc)
    {
      hsherpa->Draw("e");
      hmadgraph->Draw("hesame");
    }
  else
    { hmadgraph->Draw("he");
      hsherpa->Draw("esame");
    }


  float x1NDC = 0.691;
  float y1NDC = 0.757;
  float x2NDC = 0.894;
  float y2NDC = 0.973;

  TLegend* leg = new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->AddEntry(hsherpa, "SHERPA");
  leg->AddEntry(hmadgraph, "MADGRAPH");
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
  hscale->SetMaximum(3.0);
  hscale->SetMinimum(-0.5);
  hscale->SetTitleOffset(1.2,"Y");
  hscale->Draw("e1");
  TF1* fline = new TF1("fline","pol1");
  TLine* l2 = new TLine(xmin,1.,xmax,1.);
  l2->SetLineColor(4);
  l2->SetLineStyle(3);
  fline->SetLineWidth(3);
  fline->SetLineColor(6);
  fline->SetNpx(2500);
  hscale->Fit("fline","","");
  l2->Draw("same");


  string dirName = "compareSherpaMadgraph";
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
		     
