#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TF1.h>
#include <iostream.h>
#include <fstream.h>
#include <TFile.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMinuit.h>
#include <vector>
#include <TPaveStats.h>
#include <TPad.h>
#include "setTDRStyle.C"

void computeChi2New(TH1F* hdata, TH1F* hmc, TH1F* hscale, std::string psname, int decCode, bool logy, float ymax)
{
  // check if one should use ndivision5 or not
  bool useNDivison5=true;
  if(psname.find("sceta") != std::string::npos || psname.find("scphi") != std::string::npos)
    useNDivison5=false;
  if(psname.find("scetawidth") != std::string::npos || psname.find("scphiwidth") != std::string::npos)
    useNDivison5=true;


  // scale the area of the cloned histogram to 1
  char name[100];
  float fmc = (float)hdata->GetEntries()/(float)hmc->GetEntries();
  float fdata = 1.0;

  cout << "Data total entries = " << hdata->GetEntries() << endl;
  cout << "MC total entries = " << hmc->GetEntries() << endl;
  cout << "fmc = " << fmc << endl;
  cout << "overall scale factor = " << hscale->GetBinContent(1) << endl;
  
  // if we want to use an overall scale factor, uncommented the following line
  fmc = hscale->GetBinContent(1) ;
  cout << "Now fmc = " << fmc << endl;

  TH1F *hmcc = (TH1F*)hmc->Clone();
  hmcc->SetName("hmcc");
  hmcc->Sumw2();
  hmcc->Scale(fmc);

  TH1F *hdatac = (TH1F*)hdata->Clone();
  hdatac->SetName("hdatac");
  hdatac->Sumw2();
  hdatac->Scale(fdata);
  
  int n= hmc->GetNbinsX();
  double low = hmc->GetBinLowEdge(1);
  double high= hmc->GetBinLowEdge(1+n);
  double chi2=0;
  int realbin =0;
 

  hscale->Reset();
  // loop over n bins and compute chi2, also fill hscale
 for(int i=1;i<=n;i++){

    double nmc=hmcc->GetBinContent(i);
    double ndata=hdatac->GetBinContent(i); 
    double nmcerr=hmcc->GetBinError(i);
    double ndataerr=hdatac->GetBinError(i); 
    
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

 //  cout << "n = " << n << endl;
 cout << "real bin = " << realbin << endl;
  
 cout << "chi2=" << chi2 << endl;
 
 double prob = hdatac->KolmogorovTest(hmcc);
 cout << "prob KS = " << prob << endl;
 if(realbin > 1)
   {
     cout << "chi2/N=" << (float)chi2/(float)(realbin-1) << endl;
     prob = TMath::Prob(chi2,realbin-1);
     cout << "prob = " << prob << endl;
   }
 
  int dataColor=1;
//  int mcColor = 64;
  int mcColor = kYellow;
  
  TCanvas* c1 = new TCanvas("c1","",700,1000);  
  c1->Divide(1,2,0.01,0);
  c1->cd(1);
  if(logy)
    gPad->SetLogy(1);
  gPad->SetTopMargin(0.01);
  gPad->SetBottomMargin(0);
  gPad->SetRightMargin(0.04);
  hdatac->SetMarkerColor(dataColor);
  hdatac->SetLineColor(dataColor);
  hmcc->SetMarkerSize(0);
  hmcc->SetMarkerColor(mcColor);
  hmcc->SetLineColor(mcColor);
  hmcc->SetFillColor(mcColor);
  hmcc->SetFillStyle(1001);
  hmcc->SetXTitle("");
  hmcc->GetYaxis()->SetDecimals();
  hdatac->GetYaxis()->SetDecimals();
  if(useNDivison5)
    {
      hmcc->GetXaxis()->SetNdivisions(5);
      hdatac->GetXaxis()->SetNdivisions(5);
    }
  hdatac->SetXTitle("");
  setTDRStyle();
 
  if(!logy)
    {
      float mini = -hdatac->GetMaximum()*0.05;
      hmcc->SetMinimum(mini);
      hdatac->SetMinimum(mini);
    }

  float datamaxi = hdatac->GetMaximum()+hdatac->GetBinError(hdatac->GetMaximumBin());
  float mcmaxi   = hmcc->GetMaximum()+hmcc->GetBinError(hmcc->GetMaximumBin());

  if(mcmaxi > datamaxi)
    {
      hmcc->Draw("hist");
      hdatac->Draw("e1,same");
    }
  else{
      hdatac->Draw("e1");
      hmcc->Draw("hist,same");
  }

  hdatac->Draw("e,same");
  TLatex latexLabel;
  latexLabel.SetTextSize(0.055);
  latexLabel.SetNDC();
  latexLabel.DrawLatex(0.56, 0.87, "CMS Preliminary 2009");
  if(psname.find("2360") != std::string::npos)
    latexLabel.DrawLatex(0.56, 0.80, "#sqrt{s} = 2.36 TeV");
  else
    latexLabel.DrawLatex(0.56, 0.80, "#sqrt{s} = 900 GeV");

  TLine* l1 = new TLine(low,0,high,0);
  l1->SetLineColor(2);
  l1->SetLineStyle(3);
  l1->Draw("same");

  sprintf(name,"#chi^{2}/N = %d/%d",chi2,realbin-1 );

  TLegend* leg4 = new TLegend(0.60302,0.455274,0.80987,0.714113);
  leg4->SetFillColor(0);
  leg4->SetFillStyle(0);
  leg4->SetTextSize(0.055);
  leg4->SetBorderSize(0);
  std::string dataname = "Data";
  if(decCode%10==1) dataname += " (barrel)";
  else if(decCode%10==2) dataname += " (endcap)";
  leg4->AddEntry(hdatac, dataname.data());
  leg4->AddEntry(hmcc,"MC");
  leg4->Draw("same");

  TPaveText *pt = new TPaveText(0.60302,0.3,0.80987,0.4,"NDCBR");
  pt->SetBorderSize(0);
  pt->SetFillColor(0);
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetLineWidth(3);
  pt->SetTextAlign(12);
  //  pt->SetTextFont(42);
  pt->SetTextSize(0.055);

  if(decCode >= 10 && decCode < 20){
    pt->AddText("SC #eta > 0");
    pt->Draw();
  }
  else if(decCode >= 20 && decCode < 30){
    pt->AddText("SC #eta < 0");
    pt->Draw();
  }


  c1->cd(2);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0);
  gPad->SetBottomMargin(0.2);
  gPad->SetTickx();
  gStyle->SetOptFit(1);
  hscale->SetTitle("");
  hscale->GetYaxis()->SetDecimals();
  if(useNDivison5)
    hscale->GetXaxis()->SetNdivisions(5);
  hscale->SetMinimum(-0.5);
  float maximum = hscale->GetMaximum()+2.0;
  if(ymax>-999.)hscale->SetMaximum(ymax);
  hscale->Draw("e1");
  TF1* fline = new TF1("fline","pol1");
  TLine* l2 = new TLine(low,1.,high,1.);
  l2->SetLineColor(4);
  l2->SetLineStyle(3);
  fline->SetLineWidth(3);
  fline->SetLineColor(6);
  fline->SetNpx(2500);
  //  hscale->Fit("fline","","");
  l2->Draw("same");
  std::string dirname = logy ? 
    "$CMSSW_BASE/src/CRAB/preprod_figures/logy":
    "$CMSSW_BASE/src/CRAB/preprod_figures/";

  if(decCode%10==1) psname += "_barrel";
  else if(decCode%10==2) psname += "_endcap";
  else if(decCode%10==0) psname += "_all";

  if(decCode>=10 && decCode<20) psname += "_possceta";
  else if(decCode>=20 && decCode<30) psname += "_negsceta";
  

  std::string filename = dirname + psname + ".eps";
  c1->Print(filename.data());
//   filename = dirname + psname + ".gif";  
//   c1->Print(filename.data());

  delete c1;
  delete l1;
  delete l2;
  delete leg4;
  delete fline;
}




