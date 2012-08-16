#include "/afs/cern.ch/user/s/syu/scripts/chi2Nbins.h"
#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"
#include <string>
#include <iostream>
#include <TSystem.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TH2.h>
#include <TFile.h>
#include <TLegend.h>
#include <TMath.h>

using namespace std;

void displayTwofiles(std::string file1, std::string file2, 
 		     std::string var1, 
		     std::string var2="",
  		     std::string output="test",
		     float xmin=-9999.0, float xmax=-9999.0,
		     std::string xtitle="Number of interactions", 
		     std::string ytitle="Arbitrary Unit",		     
		     bool logScale=true ,bool reBin=false
		     )
{
  if(var2 ==  "" )var2=var1;
  
  setTDRStyle();
  gStyle->SetOptStat(0);

  TH1F* h1;
  TH1F* h2;

  char tempName[300];

  TCanvas* c1 = new TCanvas("c1","",0, 0, 500,500);
  if(logScale)c1->SetLogy(1);
  else c1->SetLogy(0);

  // first get the histogram files
  TFile *f1 = TFile::Open(file1.data());
  TFile *f2 = TFile::Open(file2.data());

  h1  = (TH1F*)(f1->Get(var1.data()));
  h1->SetXTitle(xtitle.data());  
  if(reBin)
    h1->Rebin(4);
  h1->SetYTitle(ytitle.data());  
  h1->SetTitleOffset(1.1,"Y");

  h2  = (TH1F*)(f2->Get(var2.data()));
  h2->SetXTitle(xtitle.data());
  if(reBin)
    h2->Rebin(4);
  h2->SetYTitle(ytitle.data());  
  h2->SetTitleOffset(1.1,"Y");

  h1->GetYaxis()->SetDecimals();
  h2->GetYaxis()->SetDecimals();
  h1->GetXaxis()->SetNdivisions(5);
  h2->GetXaxis()->SetNdivisions(5);
  if(xmin>-9999.0 && xmax>-9999.0)
    {
      h1->GetXaxis()->SetRangeUser(xmin,xmax);
      h2->GetXaxis()->SetRangeUser(xmin,xmax);
    }


  h2->SetLineColor(4);
  h2->SetMarkerColor(4);
  h2->SetMarkerSize(1);
  h2->SetMarkerStyle(21);
  h1->Draw();
  cout << "h1 mean = " << h1->GetMean() << " and width = " << h1->GetRMS() << 
    " and entries = " << h1->GetEntries() << " and integral = " << 
    h1->Integral() << endl;

  h1->SetLineColor(2);
  h1->SetMarkerColor(2);
  h1->SetMarkerSize(1);
  h1->SetMarkerStyle(24);
  h2->Draw();
  cout << "h2 mean = " << h2->GetMean() << " and width = " << h2->GetRMS() << 
    " and entries = " << h2->GetEntries() << " and integral = " << 
    h2->Integral() << endl;

  // if normalizing to the same area, set the scale 
  float scale1 = 1.0/(float)h1->Integral();
  float scale2 = 1.0/(float)h2->Integral();


  int binLo = -1;
  int binHi = -1;
  if(xmin>-9999.0 && xmax>-9999.0)
    {

      binLo = h1->FindBin(xmin);
      binHi = h1->FindBin(xmax)-1;

      scale1 = 1.0/(float)h1->Integral(binLo,binHi);
      scale2 = 1.0/(float)h2->Integral(binLo,binHi);
      cout << "binLo = " << binLo << ", binHi = " << binHi << endl;

    }


  h1->SetTitle("");
  h2->SetTitle("");
  h1->Sumw2();
  h1->Scale(scale1);
  h2->Sumw2();
  h2->Scale(scale2);

  float max1   = h1->GetBinError(h1->GetMaximumBin()) + h1->GetMaximum();
  float max2   = h2->GetBinError(h2->GetMaximumBin()) + h2->GetMaximum();

  if(max1 > max2)
    {
      h1->Draw("h");
      h2->Draw("hsame");
    }
  else
    { h2->Draw("h");
      h1->Draw("hsame");
    }
  cout << "Now h1 integral = " << h1->Integral() << " and h2 integral = " << h2->Integral() << endl;
  cout << "difference is " << (float)h1->Integral()/(float)h2->Integral() << endl;

//   double chi2=0.;
//   int nbins=0;
//   chi2NbinsCompare(h1,h2,chi2,nbins,binLo,binHi);
//   double kstestProb1 = h1->KolmogorovTest(h2,"X");
//   double kstestProb2 = h2->KolmogorovTest(h1,"X");

//   cout << "KS test prob 1 = " << kstestProb1 << endl;
//   cout << "KS test prob 2 = " << kstestProb2 << endl;

  float x1NDC = 0.232;
  float y1NDC = 0.318;
  float x2NDC = 0.434;
  float y2NDC = 0.531;

  if(!logScale)
    {
      x1NDC = 0.457661;
      y1NDC = 0.680085;
      x2NDC = 0.659274;
      y2NDC = 0.891949;

    }

  TLegend* leg = new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  
  //   leg->SetHeader(Form("#chi^{2}/NDF=%.1f/%d"
  // 		      ,chi2,nbins));
  //   leg->SetHeader(Form("#chi^{2} Prob=%.2f"
  // 		      ,TMath::Prob(chi2,nbins)));
  //   //    leg->SetHeader(var1.data());
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  if(var2==var1)
    {
      leg->AddEntry(h1, file1.data());
      leg->AddEntry(h2, file2.data());
    }
  else
    {
      leg->AddEntry(h1,h1->GetName());
      leg->AddEntry(h2,h2->GetName());
    }
  leg->Draw("same");

  gSystem->mkdir("compareMC");

  std::string filename;
  std::string psname = "compareMC/" + var1;
  if(output !="test")
    psname = "compareMC/"+ output;
  if(logScale)
    psname += "_logScale";
  else
    psname += "_linScale";
  filename = psname + ".eps";
  c1->Print(filename.data());
  filename = psname + ".gif";
  c1->Print(filename.data());
  filename = psname + ".pdf";
  c1->Print(filename.data());
  //   c1->Close();

  //   TCanvas* c2 = new TCanvas("c2","",500,0,500,500);
  //   TH1F* hratio = (TH1F*)h1->Clone("hratio");
  //   hratio->Reset();

  //   for(int i=1; i<= h1->GetNbinsX(); i++){
    
  //     double r = h2->GetBinContent(i)>0? 
  //       h1->GetBinContent(i)/h2->GetBinContent(i):-1;
  //     hratio->SetBinContent(i,r);
  //     cout << "bin " << i << " ratio = " << r << endl;

  //   }


  //   std::string remword ="_flat.root";
  //   std::string epsname = file2;
  //   size_t pos  = epsname.find(remword);

  //   if(pos!= std::string::npos)
  //     epsname.swap(epsname.erase(pos,remword.length()));

  //   hratio->SetMinimum(0.8);
  //   hratio->SetMaximum(1.2);
  //   hratio->Draw();
  //   hratio->SetXTitle("p_{T}(#gamma)");
  //   leg->Clear();
  //   leg->SetHeader(epsname.data());
  //   leg->AddEntry(hratio,"Double ratio: x^{-2}/x^{0} correction");
  //   leg->Draw("same");
  //   TLine* a = new TLine(h1->GetBinLowEdge(1),1.0,
  // 		       h1->GetBinLowEdge(h1->GetNbinsX()+1),1.0);
  //   a->Draw();


  //   epsname += ".eps";

  //   c2->Print(epsname.data());

}
		     
