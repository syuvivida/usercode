#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"
#include <string>
#include <iostream>

using namespace std;

void compareDressed(std::string mcfilePostfix, 
		    std::string var1="h_ystar", 
		    float xmin=-9999, float xmax=-9999,
		    float ymin=0.9, float ymax=1.2,
		    std::string var2="",
		    std::string mcName1="Dressed",
		    std::string mcName2="Bare",
		    bool logScale=false
		    )
{
  
  setTDRStyle();
  gStyle->SetOptStat(0);

  TH1F* h1;
  TH1F* h2;

  char tempName[300];
  if(var2 ==  "" )var2=var1;
 
  std::string mcfile1 = "dressed_exclusive1Jet_zPt40_" + mcfilePostfix;
  std::string mcfile2 = "weighted_exclusive1Jet_zPt40_" + mcfilePostfix;

  std::string header;
  std::string output;
  if(mcfile1.find("madgraph")!=std::string::npos)
    {
      header="Madgraph";
      output="madgraph";
    }
  else if(mcfile1.find("sherpa")!=std::string::npos)
    {
      header="Sherpa";
      output="sherpa";
    }
      
  if(mcfile1.find("electron")!=std::string::npos)
    {
      header+= " e";
      output+= "E";
    }
  else if(mcfile1.find("muon")!=std::string::npos)
    {
      header+= " #mu";
      output+= "Mu";
    }

 

  // first get the histogram files
  TFile *fmc1 = TFile::Open(mcfile1.data());
  cout << "Reading file 1: " << fmc1->GetName() << endl;

  TFile *fmc2   = TFile::Open(mcfile2.data());
  cout << "Reading file 2: " << fmc2->GetName() << endl;


  h1  = (TH1F*)(fmc1->FindObjectAny(var1.data()));
  h2    = (TH1F*)(fmc2->FindObjectAny(var2.data()));

  TH1D* hratio =(TH1D*) h1->Clone("hratio");
  hratio->SetYTitle(Form("%s/%s",mcName1.data(),mcName2.data()));
  hratio->SetLineColor(2);
  hratio->SetMarkerColor(2);

  h1->GetXaxis()->SetNdivisions(5);
  h1->GetYaxis()->SetDecimals();
  h1->SetTitleOffset(1.2,"Y");

  h2->GetXaxis()->SetNdivisions(5);
  h2->GetYaxis()->SetDecimals();
  h2->SetTitleOffset(1.2,"Y");

  hratio->GetXaxis()->SetNdivisions(5);
  hratio->GetYaxis()->SetDecimals();


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




  cout << "h1 integral = "   << h1->Integral() << endl;;
  cout << "h2 integral = " << h2->Integral() << endl;

  hratio->Divide(h1, h2, 1,1,"B");

  for(int i=1;i<=hratio->GetNbinsX();i++)
    cout << "Bin " << i << " ( " << hratio->GetBinLowEdge(i) << "~" << hratio->GetBinLowEdge(i+1) << " ): " 
	 << h1->GetBinContent(i) << "/" << h2->GetBinContent(i) << " = " 
	 << hratio->GetBinContent(i) << " +- " << hratio->GetBinError(i) << endl;

  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h2->GetXaxis()->SetRangeUser(xmin,xmax);
  hratio->GetXaxis()->SetRangeUser(xmin,xmax);


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
  hratio->SetTitle("");
  hratio->SetMaximum(ymax);
  hratio->SetMinimum(ymin);
  hratio->SetTitleOffset(1.2,"Y");
  hratio->Draw("e1");
  TF1* fline = new TF1("fline","pol1");
  TLine* l2 = new TLine(xmin,1.,xmax,1.);
  l2->SetLineColor(4);
  l2->SetLineStyle(3);
  fline->SetLineWidth(3);
  fline->SetLineColor(6);
  fline->SetNpx(2500);
  //   hratio->Fit("fline","","");
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
		     
