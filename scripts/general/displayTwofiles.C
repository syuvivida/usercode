#include "/afs/cern.ch/user/s/syu/scripts/chi2Nbins.h"
#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"

void displayTwofiles(std::string file1, std::string file2, 
 		     std::string var1, std::string var2="",
		     float xmin=-9999.0, float xmax=-9999.0,
		     //  		     std::string xtitle="", std::string ytitle="",
		     std::string output="test", bool logScale=false, bool reBin=false)
{
  if(var2 ==  "" )var2=var1;
  
  setTDRStyle();
  gStyle->SetOptStat(0);

  TH1F* h1;
  TH1F* h2;

  char tempName[300];

  TCanvas* c1 = new TCanvas("c1","",500,500);
  if(logScale)c1->SetLogy(1);
  else c1->SetLogy(0);

  // first get the histogram files
  TFile *f1 = TFile::Open(file1.data());
  TFile *f2 = TFile::Open(file2.data());

  h1  = (TH1F*)(f1->Get(var1.data()));
  //   h1->SetXTitle(xtitle.data());  
  if(reBin)
    h1->Rebin(4);
  //   h1->SetYTitle(ytitle.data());  

  h2  = (TH1F*)(f2->Get(var2.data()));
  //   h2->SetXTitle(xtitle.data());
  if(reBin)
    h2->Rebin(4);
  //   h2->SetYTitle(ytitle.data());  

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

  double chi2=0.;
  int nbins=0;
  chi2NbinsCompare(h1,h2,chi2,nbins,binLo,binHi);
  double kstestProb1 = h1->KolmogorovTest(h2,"X");
  double kstestProb2 = h2->KolmogorovTest(h1,"X");

  cout << "KS test prob 1 = " << kstestProb1 << endl;
  cout << "KS test prob 2 = " << kstestProb2 << endl;

  float x1NDC = 0.415;
  float y1NDC = 0.701;
  float x2NDC = 0.617;
  float y2NDC = 0.915;

  if(var1.find("eta")!= std::string::npos || 
     var1.find("h_sieie_leadingEB")!= std::string::npos )
     
    {
      x1NDC = 0.333;
      y1NDC = 0.182;
      x2NDC = 0.536;
      y2NDC = 0.396;
    }

  TLegend* leg = new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  
//   leg->SetHeader(Form("#chi^{2}/NDF=%.1f/%d"
// 		      ,chi2,nbins));
  leg->SetHeader(Form("#chi^{2} Prob=%.2f"
		      ,TMath::Prob(chi2,nbins)));
  //    leg->SetHeader(var1.data());
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  if(var2==var1)
    {
      //        leg->AddEntry(h1,"Data");
      //        leg->AddEntry(h2,"MC");
      //        leg->AddEntry(h1, "Default");
      //        leg->AddEntry(h2, "No IP cuts");
      leg->AddEntry(h1, file1.data());
      leg->AddEntry(h2, file2.data());
      //       leg->AddEntry(h1, "q#bar{q} annilation");
      //       leg->AddEntry(h2, "q-gluon scattering");
      //        leg->AddEntry(h1,"Spikes");
      //        leg->AddEntry(h2,"Signal MC 20 < p_{T}(#gamma) < 30");
    }
  else
    {
      //        leg->AddEntry(h1,"Smeared Random Cone Data");
      //        leg->AddEntry(h2,"Signal MC 20 < p_{T}(#gamma) < 30");
      //        leg->AddEntry(h1,"Spikes");
      //        leg->AddEntry(h2,"Smeared Random Cone Data");
      //        leg->AddEntry(h2,"PhotonJet_Pt30");
      //        leg->AddEntry(h2,"Random Cone");
      //        leg->AddEntry(h1,"Data");
      //        leg->AddEntry(h2,"Background MC");   
      leg->AddEntry(h1,h1->GetName());
      leg->AddEntry(h2,h2->GetName());
      //        leg->AddEntry(h1,"Data");
      //        leg->AddEntry(h2,"MC");
    }
  leg->Draw("same");

  gSystem->mkdir("compareMC");

  std::string filename;
  std::string psname = "compareMC/" + var1;
  if(output.data()!="")
    psname = "compareMC/"+ output;
  filename = psname + ".eps";
  c1->Print(filename.data());
  filename = psname + ".gif";
  c1->Print(filename.data());
  filename = psname + ".pdf";
  c1->Print(filename.data());
}
		     
