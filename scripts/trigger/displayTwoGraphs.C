void displayTwoGraphs(std::string file1, std::string file2, std::string histo1,
		      std::string xtitle="", std::string ytitle="",
		      std::string histo2="", std::string title="",std::string leg1="",
		      std::string leg2="", std::string output="test", float xmin=-1, float xmax=-1, float ymin=-1, float ymax=-1)
{
  if(histo2 == ""  )histo2 = histo1;
  if(file2  == "")file2=file1;
  TGraph* h1;
  TGraph* h2;
  TCanvas* c1 = new TCanvas("c1");
  TH1F *hr = c1->DrawFrame(0.0,0.0,60.0,2.0);
  TFile *f1 = new TFile(file1.data());
  cout << "Opening " << file1.data() << endl;
  h1 = dynamic_cast<TGraph*>(gDirectory->Get(histo1.data()));
  h1->Draw();
  hr->SetFillColor(0);
  hr->SetXTitle(xtitle.data());
  hr->SetYTitle(ytitle.data());
  h1->SetTitle("");
  if(xmin<xmax)
    h1->GetXaxis()->SetRangeUser(xmin,xmax);

  h1->SetMarkerColor(2);
  h1->SetLineColor(2);
  h1->GetYaxis()->SetDecimals();
  h1->SetLineWidth(1);

  TFile *f2 = new TFile(file2.data());
  cout << "Opening " << file2.data() << endl;
  h2 = dynamic_cast<TGraph*>(gDirectory->Get(histo2.data()));
  h2->Draw();
  h2->SetMarkerColor(4);
  h2->SetMarkerStyle(30);
  h2->SetLineColor(4);
  h2->SetLineWidth(1);
  h2->GetYaxis()->SetDecimals();
  h2->GetXaxis()->SetTitleSize(0.5);
  h2->GetYaxis()->SetTitleSize(0.5);
  h2->GetXaxis()->SetTitleOffset(1.0);
  h2->GetYaxis()->SetTitleOffset(1.0);
  h2->GetXaxis()->SetRangeUser(0.0,60.0);
  h2->SetTitle("");

  if(ymin<ymax)
    {
      h1->SetMaximum(ymax);
      h1->SetMinimum(ymin);
    }
  h1->Draw("ap");
  h2->Draw("psame");

   TLegend* leg = new TLegend(0.4612,0.315678,0.613506,0.53178) ;
//   TLegend* leg = new TLegend(0.2612,0.315678,0.613506,0.53178) ;
 
  leg->SetFillColor(0);
//   leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetHeader(title.data());
 
  if(leg1=="")
    leg->AddEntry(h1,file1.data());
  else 
    leg->AddEntry(h1,leg1.data());

  if(leg2=="")
    leg->AddEntry(h2,file2.data());
  else
    leg->AddEntry(h2,leg2.data());
  leg->Draw("same");

  std::string temp= output + ".eps";
  c1->Print(temp.data());
  temp= output + ".gif";
  c1->Print(temp.data());
  temp= output + ".C";
  c1->Print(temp.data());

}
