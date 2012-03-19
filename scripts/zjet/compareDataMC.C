#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"

void compareDataMC(std::string datafile, std::string mcfile, 
		   std::string var,
		   float xmin=-9999.0, float xmax=-9999.0,
		   bool logScale=false, 
		   std::string output="test")
{
  
  setTDRStyle();
  gStyle->SetOptStat(0);

  TH1F* hdata;
  TH1F* hmc;

  char tempName[300];


  // first get the histogram files
  TFile *fdata = TFile::Open(datafile.data());
  TFile *fmc   = TFile::Open(mcfile.data());

  hdata  = (TH1F*)(fdata->Get(var.data()));
  hmc    = (TH1F*)(fmc->Get(var.data()));

  TH1D* hscale =(TH1D*) hdata->Clone("hscale");
  hscale->SetYTitle("Data/MC");

  int nREBIN=2;
  if(var.find("pt")!= std::string::npos)
    nREBIN=4;

  hdata->GetXaxis()->SetNdivisions(5);
  hdata->GetYaxis()->SetDecimals();
  hdata->Rebin(nREBIN);

  hmc->GetXaxis()->SetNdivisions(5);
  hmc->GetYaxis()->SetDecimals();
  hmc->Rebin(nREBIN);

  hscale->GetXaxis()->SetNdivisions(5);
  hscale->GetYaxis()->SetDecimals();
  hscale->Rebin(nREBIN);

  hdata->SetLineColor(2);
  hdata->SetMarkerColor(2);
  hdata->SetMarkerSize(1);
  hdata->SetMarkerStyle(24);


  hmc->SetLineColor(4);
  hmc->SetMarkerColor(4);
  hmc->SetMarkerSize(1);
  hmc->SetMarkerStyle(21);

  // if normalizing to the same area, set the scale 

  int binLo = -1;
  int binHi = -1;
  int nbins = hdata->GetNbinsX();
  if(xmin>-9999.0 && xmax>-9999.0)
    {

      binLo = hdata->FindBin(xmin);
      binHi = hdata->FindBin(xmax)-1;

    }

  else
    {
      binLo = 1;
      binHi = nbins;
      xmin = hdata->GetBinLowEdge(1);
      xmax = hdata->GetBinLowEdge(nbins+1);
    }


  float scale_mc = (float)hdata->Integral(binLo,binHi)/(float)hmc->Integral(binLo,binHi);
  cout << "binLo = " << binLo << ", binHi = " << binHi << endl;
  cout << "xmin = " << xmin << "xmax = " << xmax << endl;

  hmc->Sumw2();
  hmc->Scale(scale_mc);

  // get the ratio
  double chi2 = 0;
  int realbin = 0;
  for(int i=1;i<= nbins;i++){

    double nmc=hmc->GetBinContent(i);
    double ndata=hdata->GetBinContent(i); 
    double nmcerr=hmc->GetBinError(i);
    double ndataerr=hdata->GetBinError(i); 
    
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

  
  hdata->GetXaxis()->SetRangeUser(xmin,xmax);
  hmc->GetXaxis()->SetRangeUser(xmin,xmax);
  hscale->GetXaxis()->SetRangeUser(xmin,xmax);


  TCanvas* c1 = new TCanvas("c1","",700,1000);  
  c1->Divide(1,2,0.01,0);
  c1->cd(1);
  if(logScale)
    gPad->SetLogy(1);
  gPad->SetTopMargin(0.01);
  gPad->SetBottomMargin(0);
  gPad->SetRightMargin(0.04);

  
  float max_data  = hdata->GetBinError(hdata->GetMaximumBin()) + hdata->GetMaximum();
  float max_mc    = hmc->GetBinError(hmc->GetMaximumBin()) + hmc->GetMaximum();

  if(max_data > max_mc)
    {
      hdata->Draw("e");
      hmc->Draw("hesame");
    }
  else
    { hmc->Draw("he");
      hdata->Draw("esame");
    }


  float x1NDC = 0.616;
  float y1NDC = 0.687;
  float x2NDC = 0.818;
  float y2NDC = 0.901;

  TLegend* leg = new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->AddEntry(hdata, "2011 data");
  leg->AddEntry(hmc, "DY+Jets MC");
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


  string dirName = "compareDataMC";
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
		     
