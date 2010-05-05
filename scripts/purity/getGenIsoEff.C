#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>

void getGenIsoEff(std::string inputfile, const float thre = 0.90)
{
  TFile *inf = new TFile(inputfile.data());
  TH1F *hiso=(TH1F*)inf->FindObjectAny("histIsoDR04_all");
  hiso->SetTitle("");
  hiso->SetName("hiso");
  TH1F *heff= hiso->Clone();
  heff->Reset();
  heff->SetTitle("");
  heff->SetName("heff");
  //  heff->SetTitle("gen pho pt > 20 GeV, |eta| < 2.5 ");
  heff->SetXTitle("Generator isolation cut value [GeV]");
  heff->SetYTitle("Generator isolation cut efficiency");
  heff->SetTitleOffset( 1.400,"Y");
  heff->GetYaxis()->SetDecimals();


  TH1F *hdiff= hiso->Clone();
  hdiff->Reset();
  hdiff->SetName("hdiff");
  //  hdiff->SetTitle("gen pho pt > 20 GeV, |eta| < 2.5 ");
  hdiff->SetXTitle("Generator isolation cut value [GeV]");
  hdiff->SetYTitle(Form("Generator isolation cut efficiency - %.2f",thre));  
  hdiff->SetTitleOffset( 1.600,"Y");
  hdiff->GetYaxis()->SetDecimals();

  float total = hiso->GetEntries();      
  for(int i=1; i<=hiso->GetNbinsX(); i++)
    {
      float temp_eff = hiso->Integral(0,i)/total;
      float temp_diff = fabs(temp_eff - thre);
      
      heff->SetBinContent(i,temp_eff);
      hdiff->SetBinContent(i,temp_diff);
    } 

  TCanvas* c1=new TCanvas("c1","",1000,500);
  c1->Divide(2,1);
  c1->cd(1);
  heff->Draw();
  c1->cd(2);
  hdiff->Draw();
  int minBin = hdiff->GetMinimumBin();
  float minEff = hdiff->GetBinLowEdge(minBin+1);
  TPaveText *pt = new TPaveText(0.25302,0.7,0.80987,0.9,"NDCBR");
  pt->SetBorderSize(0);
  pt->SetFillColor(0);
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetLineWidth(3);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.04);

  cout << thre*100 << "% efficiency point is " << minEff << endl;
  pt->AddText("Gen #gamma pt > 20 GeV, |eta| < 2.5"); 
  pt->AddText(Form("%.2f efficiency point: iso < %.1f",thre,minEff));
  pt->Draw();

  c1->Print(Form("~/public/%s_%.2f.gif",inputfile.data(),thre));  



}
