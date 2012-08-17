#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"
#include <string>
#include <iostream>

using namespace std;

void compareDressed_darko(std::string mcfilePostfix, 
			  std::string var1="h_ystar", 
			  float xmin=-9999, float xmax=-9999,
			  float ymin=0.8, float ymax=1.1,
			  bool logScale=false
			  )
{
  
  setTDRStyle();
  gStyle->SetOptStat(0);

  TH1D* h_numr;
  TH1D* h_deno;

  TH1D* h_both;

  char tempName[300];
 
  std::string mcfile_numr = "bare_exclusive1Jet_zPt40_" + mcfilePostfix;
  std::string mcfile_deno = "dressed_exclusive1Jet_zPt40_" + mcfilePostfix;
  std::string mcfile_both = "both_exclusive1Jet_zPt40_" + mcfilePostfix;

  std::string mcName_numr="Bare";
  std::string mcName_deno="Dressed";

  std::string header;
  std::string output;
  if(mcfile_numr.find("madgraph")!=std::string::npos)
    {
      header="Madgraph";
      output="madgraph";
    }
  else if(mcfile_numr.find("sherpa")!=std::string::npos)
    {
      header="Sherpa";
      output="sherpa";
    }
      
  if(mcfile_numr.find("electron")!=std::string::npos)
    {
      header+= " e";
      output+= "E";
    }
  else if(mcfile_numr.find("muon")!=std::string::npos)
    {
      header+= " #mu";
      output+= "Mu";
    }

 

  // first get the histogram files
  TFile *fmc1 = TFile::Open(mcfile_numr.data());
  cout << "Reading file 1: " << fmc1->GetName() << endl;

  TFile *fmc2 = TFile::Open(mcfile_deno.data());
  cout << "Reading file 2: " << fmc2->GetName() << endl;

  TFile *fmc3 = TFile::Open(mcfile_both.data());
  cout << "Reading file 3: " << fmc3->GetName() << endl;

  h_numr  = (TH1D*)(fmc1->FindObjectAny(var1.data()));
  h_deno  = (TH1D*)(fmc2->FindObjectAny(var1.data()));
  h_both  = (TH1D*)(fmc3->FindObjectAny(var1.data()));

  TH1D* hratio =(TH1D*) h_numr->Clone("hratio");
  hratio->SetYTitle(Form("%s/%s",mcName_numr.data(),mcName_deno.data()));
  hratio->SetLineColor(1);
  hratio->SetMarkerColor(1);

  h_numr->GetXaxis()->SetNdivisions(5);
  h_numr->GetYaxis()->SetDecimals();
  h_numr->SetTitleOffset(1.2,"Y");

  h_deno->GetXaxis()->SetNdivisions(5);
  h_deno->GetYaxis()->SetDecimals();
  h_deno->SetTitleOffset(1.2,"Y");

  hratio->GetXaxis()->SetNdivisions(5);
  hratio->GetYaxis()->SetDecimals();


  h_numr->SetLineColor(2);
  h_numr->SetMarkerColor(2);
  h_numr->SetMarkerSize(1);
  h_numr->SetMarkerStyle(24);


  h_deno->SetLineColor(4);
  h_deno->SetMarkerColor(4);
  h_deno->SetMarkerSize(1);
  h_deno->SetMarkerStyle(21);

  // if normalizing to the same area, set the scale 

  int binLo = -1;
  int binHi = -1;
  int nbins = h_numr->GetNbinsX();
  if(xmin>-9999.0 && xmax>-9999.0)
    {

      binLo = h_numr->FindBin(xmin);
      binHi = h_numr->FindBin(xmax)-1;

    }

  else
    {
      binLo = 1;
      binHi = nbins;
      xmin = h_numr->GetBinLowEdge(1);
      xmax = h_numr->GetBinLowEdge(nbins+1);
    }




  cout << "h_numr integral = "   << h_numr->Integral() << endl;;
  cout << "h_deno integral = " << h_deno->Integral() << endl;

  float area_h_deno = h_deno->Integral(binLo, binHi);

  h_numr->Sumw2();
  h_numr->Scale(1.0/area_h_deno);

  h_deno->Sumw2();
  h_deno->Scale(1.0/area_h_deno);

  h_both->Sumw2();
  h_both->Scale(1.0/area_h_deno);

  // now use formulas similar to Darko's for error calculation
  hratio->Reset();
  for(int i=1; i<=hratio->GetNbinsX();i++){
    
    double n_n = h_numr->GetBinContent(i);
    double n_d = h_deno->GetBinContent(i);
    
    if(n_d<1e-6)continue;

    double diff = n_n-n_d;

    double err_n = h_numr->GetBinError(i);
    double err_d = h_deno->GetBinError(i);

    double err_b = h_both->GetBinError(i);
    
    double err_y = sqrt(err_n*err_n - err_b*err_b);
    double err_z = sqrt(err_d*err_d - err_b*err_b);

    double value = n_n/n_d;
    
    double variance = 
      err_b*err_b*diff*diff/n_d/n_d + 
      err_y*err_y + 
      err_z*err_z*value*value;

    variance /= (n_d*n_d);

    if(variance < 0){cout << "bin " << i << " has an error in variance" << endl; continue;}
    
    variance = sqrt(variance);
    
    hratio->SetBinContent(i, value);
    hratio->SetBinError(i, variance);

  }


  for(int i=1;i<=hratio->GetNbinsX();i++)
    cout << "Bin " << i << " ( " << hratio->GetBinLowEdge(i) << "~" << hratio->GetBinLowEdge(i+1) << " ): " 
	 << h_numr->GetBinContent(i) << "/" << h_deno->GetBinContent(i) << " = " 
	 << hratio->GetBinContent(i) << " +- " << hratio->GetBinError(i) << endl;

  h_numr->GetXaxis()->SetRangeUser(xmin,xmax);
  h_deno->GetXaxis()->SetRangeUser(xmin,xmax);
  hratio->GetXaxis()->SetRangeUser(xmin,xmax);


  TCanvas* c1 = new TCanvas("c1","",700,1000);  
  c1->Divide(1,2,0.01,0);
  c1->cd(1);
  if(logScale)
    gPad->SetLogy(1);
  gPad->SetTopMargin(0.01);
  gPad->SetBottomMargin(0);
  gPad->SetRightMargin(0.04);

  
  float max_data  = h_numr->GetBinError(h_numr->GetMaximumBin()) + h_numr->GetMaximum();
  float max_mc    = h_deno->GetBinError(h_deno->GetMaximumBin()) + h_deno->GetMaximum();

  if(max_data > max_mc)
    {
      h_numr->Draw("e");
      h_deno->Draw("hesame");
    }
  else
    { h_deno->Draw("he");
      h_numr->Draw("esame");
    }


  float x1NDC = 0.7;
  float y1NDC = 0.620;
  float x2NDC = 0.9;
  float y2NDC = 0.956;

  TLegend* leg = new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  
  leg->SetHeader(header.data());
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->AddEntry(h_numr, mcName_numr.data());
  leg->AddEntry(h_deno, mcName_deno.data());
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
  hratio->SetTitle("");
  hratio->SetMaximum(ymax);
  hratio->SetMinimum(ymin);
  hratio->SetTitleOffset(1.2,"Y");
  hratio->Draw("e1");
  TF1* fline = new TF1("fline","pol0");
  TLine* l2 = new TLine(xmin,1.,xmax,1.);
  l2->SetLineColor(4);
  l2->SetLineStyle(3);
  fline->SetLineWidth(3);
  fline->SetLineColor(kMagenta);
  fline->SetNpx(2500);
  if(var1.find("mZ")== std::string::npos && 
     var1.find("zpt")== std::string::npos && 
     var1.find("jetpt")== std::string::npos)
    hratio->Fit("fline","","",0,2.0);
  l2->Draw("same");


  string dirName = "compareDressed";
  gSystem->mkdir(dirName.data());

  std::string filename;
  std::string remword  ="h_";
  std::string remword2 ="h";
  size_t pos  = var1.find(remword);
  if(pos!= std::string::npos)
    var1.replace(pos,remword2.length(),"");

  std::string psname = dirName + "/" + var1;
  if(output !="test")
    psname = dirName+ "/" + output + var1;
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
		     
