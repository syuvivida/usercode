void overLay(std::string output, std::string file1, std::string name1, std::string file2, std::string name2, std::string file3, std::string name3, std::string xtitle="", float xmin=-1, float xmax=-1, bool log=true, float YMAX=0.5)
{
  TH1F* h1;
  TH1F* h2;
  TH1F* h3;
  gStyle->SetOptStat(0);
  TFile *f_unweighted = TFile::Open(file1.data());
  h1 =(TH1F*)(f_unweighted->Get(name1.data()));

  TFile *f_weighted = TFile::Open(file2.data());
  h2 =(TH1F*)(f_weighted->Get(name2.data()));

  TFile *f_weighted = TFile::Open(file3.data());
  h3 =(TH1F*)(f_weighted->Get(name3.data()));
 
  h1->SetTitle("");
  h2->SetTitle("");
  h3->SetTitle("");

  
  h1->SetLineColor(4);
  h2->SetLineColor(2);
  h3->SetLineColor(13);

  h1->SetLineWidth(4);
  h2->SetLineWidth(4);
  h3->SetLineWidth(4);

  h1->SetMarkerSize(0);
  h2->SetMarkerSize(0);
  h3->SetMarkerSize(0);

  h1->Sumw2();
  h2->Sumw2();
  h3->Sumw2();

  float scale1=1.0,scale2=1.0,scale3=1.0;
  if(h1->Integral()>0)scale1= 1.0/(float)h1->Integral();
  if(h2->Integral()>0)scale2= 1.0/(float)h2->Integral();
  if(h3->Integral()>0)scale3= 1.0/(float)h3->Integral();

  h1->Scale(scale1);
  h2->Scale(scale2);
  h3->Scale(scale3);

  if(xmin < xmax){
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h2->GetXaxis()->SetRangeUser(xmin,xmax);
  h3->GetXaxis()->SetRangeUser(xmin,xmax);
  }

  h1->SetMaximum(YMAX);
  h2->SetMaximum(YMAX);
  h3->SetMaximum(YMAX);

  h1->SetXTitle(xtitle.data());
  h1->SetYTitle("Normalized entries");
  h1->SetTitleOffset(1.5,"X");
  h1->SetTitleOffset(1.5,"Y");
  h1->SetTitleFont(42,"X");
  h1->SetTitleSize(0.06  ,"X");
  h1->SetTitleFont(42,"Y");
  h1->SetTitleSize(0.06  ,"Y");

  TCanvas* c1 = new TCanvas("c1","",500,500);

  h1->Draw("hist");
  h2->Draw("histsame");
  h3->Draw("histsame");
  
  if(log==true)
    c1->SetLogy(1);

  //TLegend* leg = new TLegend(0.48,0.57,0.77,0.83);
  TLegend* leg = new TLegend(0.39,0.65,0.68,0.91);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.043);
  leg->SetTextFont(42);
  leg->AddEntry(h1,"PYTHIA 6.4");
  leg->AddEntry(h2,"WGAMMA_NLO born-level");
  leg->AddEntry(h3,"WGAMMA");
  leg->Draw("same");

  std::string outfile = output + ".eps";
  c1->Print(outfile.data());
  outfile = output + ".gif";
  c1->Print(outfile.data());
  outfile = output + ".C";
  c1->Print(outfile.data());

}
