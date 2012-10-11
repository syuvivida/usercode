
#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"
void displayTwoTreesRatio(std::string file, 
			  std::string var,
			  int nbin=100,
			  double xmin=0,
			  double xmax=100,
			  bool logScale=false,
			  std::string output="test",
			  std::string cut="")
{
  setTDRStyle();
  gStyle->SetOptStat(0);
 
  TChain* t1 = new TChain("tree");
  std::string file1 = "/data4/syu/52X_533_validation/533_"+file;
  t1->Add(file1.data());

  TChain* t2 = new TChain("tree");
  std::string file2 = "/data4/syu/52X_533_validation/52X_"+file;
  t2->Add(file2.data());

  std::string endfix = "MC";
  std::string legname= "MC";

  if(file.find("DoubleElectron")!= std::string::npos)
    {
      legname = "e Data";
      endfix  = "DiEle";
    }

  else if(file.find("DoubleMu")!= std::string::npos)
    {
      legname = "#mu Data";
      endfix  = "DiMuo";
    }


  
  
  TH1F* h1 = new TH1F("h1","",nbin,xmin,xmax);
  TH1F* h2 = new TH1F("h2","",nbin,xmin,xmax);
  std::string temp = var  + Form(">>h1",nbin,xmin,xmax);

  t1->Draw(temp.data(), cut.data());
  h1->SetLineColor(4);
  h1->SetName("h1");

  cout << "h1 mean = " << h1->GetMean() << " and width = " << h1->GetRMS() << 
    " and entries = " << h1->GetEntries() << endl;

  temp = var  + Form(">>h2",nbin,xmin,xmax);

  t2->Draw(temp.data(), cut.data());
  h2->SetLineColor(2);
  h2->SetName("h2");

  cout << "h2 mean = " << h2->GetMean() << " and width = " << h2->GetRMS() << 
    " and entries = " << h2->GetEntries() << endl;


  
  TCanvas* c1 = new TCanvas("c1","",600,1000);  
  c1->Divide(1,2,0.01,0);
  c1->cd(1);
  if(logScale)
    gPad->SetLogy(1);
 
  gPad->SetTopMargin(0.01);
  gPad->SetBottomMargin(0);
  gPad->SetRightMargin(0.04);

  h1->GetXaxis()->SetNdivisions(5);
  h2->GetXaxis()->SetNdivisions(5);
  h1->SetMinimum(0.1);
  h2->SetMinimum(0.1);
 
  float max_data  = h1->GetBinError(h1->GetMaximumBin()) + h1->GetMaximum();
  float max_mc    = h2->GetBinError(h2->GetMaximumBin()) + h2->GetMaximum();

  if(max_data > max_mc)
    {
      h1->Draw("hist");
      h2->Draw("histsame");
    }
  else
    { h2->Draw("hist");
      h1->Draw("histsame");
    }


  float x1NDC = 0.704977;
  float y1NDC = 0.775988;
  float x2NDC = 0.976553;
  float y2NDC = 0.991757;

  TLegend* leg = new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->AddEntry(h1, Form("53X %s", legname.data()),"l");
  leg->AddEntry(h2, Form("52X %s", legname.data()),"l");
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

  TH1F* hscale = new TH1F("hscale","",nbin,xmin,xmax);
  hscale->Divide(h1,h2,1,1,"B");

  hscale->SetTitle("");
  hscale->SetMaximum(2.0);
  hscale->SetMinimum(0.0);
  hscale->SetXTitle(var.data());
  hscale->SetYTitle("53X/52X");
  hscale->GetXaxis()->SetNdivisions(5);
  hscale->GetYaxis()->SetDecimals();
  hscale->SetTitleOffset(1.2,"Y");
  hscale->Draw("hist");


  string dirName = "validation";
  gSystem->mkdir(dirName.data());

  std::string filename;
  std::string psname = dirName + "/" + endfix + "_" + var ;
  if(output!="test")
    psname = dirName + "/" + endfix + "_" + output;
  filename = psname + ".eps";
  c1->Print(filename.data());
  filename = psname + ".gif";
  c1->Print(filename.data());
  filename = psname + ".pdf";
  c1->Print(filename.data());


  gSystem->Sleep(1000);


}
		     
