#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"

void acceptCorr(std::string filepostfix,
		std::string var, 
		bool update=false,
		float xmin=-9999.0, float xmax=-9999.0,
		bool logScale=false
		)
{
  
  setTDRStyle();
  gStyle->SetOptStat(0);

  const int NFILES=2;

  TH1F* h[NFILES];

  std::string mcfile[NFILES]={
    "darko_root/excludedECALCrack_exclusive1Jet_zPt40_electron_"+filepostfix,
    "darko_root/bare_exclusive1Jet_zPt40_electron_"+filepostfix
  };
  
  std::string mcName[NFILES]={
    "Crack excluded",
    "Crack included"
  };

  TFile* fmc[NFILES];

  // first get the histogram files
  for(int ifile=0; ifile < NFILES; ifile++){

    fmc[ifile] = TFile::Open(mcfile[ifile] .data());
    cout << "opened " << fmc[ifile]->GetName() << endl;
    h[ifile]  = (TH1F*)(fmc[ifile]->Get(var.data()));
    h[ifile]->SetName(Form("h%d",ifile));

  }

  std::string remword  ="h_";
  size_t pos  = var.find(remword);
  if(pos!= std::string::npos)
    var.replace(pos,remword.length(),"");

  TH1D* hscale =(TH1D*) h[0]->Clone(Form("hscale_%s",var.data()));
  hscale->SetYTitle(Form("%s/%s",mcName[0].data(),mcName[1].data()));

  int COLORCODE[]={2,4};
  int MARKERCODE[]={24,21};

  for(int i=0; i < NFILES; i++){
    h[i]->GetXaxis()->SetNdivisions(5);
    h[i]->GetYaxis()->SetDecimals();
    h[i]->SetLineColor(COLORCODE[i]);
    h[i]->SetMarkerColor(COLORCODE[i]);
    h[i]->SetMarkerSize(1);
    h[i]->SetMarkerStyle(MARKERCODE[i]);
  };

  hscale->GetXaxis()->SetNdivisions(5);
  hscale->GetYaxis()->SetDecimals();

  // if normalizing to the same area, set the scale 

  int binLo = -1;
  int binHi = -1;
  int nbins = hscale->GetNbinsX();
  if(xmin>-9999.0 && xmax>-9999.0)
    {

      binLo = hscale->FindBin(xmin);
      binHi = hscale->FindBin(xmax)-1;

    }

  else
    {
      binLo = 1;
      binHi = nbins;
      xmin = hscale->GetBinLowEdge(1);
      xmax = hscale->GetBinLowEdge(nbins+1);
    }

  float area_h = h[1]->Integral(binLo, binHi);

  for(int ifile=0; ifile < NFILES; ifile++){
    h[ifile]->Sumw2();
    h[ifile]->Scale(1.0/area_h);
    cout << "h[" << ifile << "] integral = " << h[ifile]->Integral() << endl;
  }
    

  
  hscale->Divide(h[0], h[1], 1,1,"B");

  for(int i=1;i<=hscale->GetNbinsX();i++)
    cout << "Bin " << i << " ( " << hscale->GetBinLowEdge(i) << "~" 
	 << hscale->GetBinLowEdge(i+1) << " ): " 
	 << h[0]->GetBinContent(i) << "/" << h[1]->GetBinContent(i) << " = " 
	 << hscale->GetBinContent(i) << " +- " << hscale->GetBinError(i) 
	 << endl;


  h[0]->GetXaxis()->SetRangeUser(xmin,xmax);
  h[1]->GetXaxis()->SetRangeUser(xmin,xmax);
  hscale->GetXaxis()->SetRangeUser(xmin,xmax);

  ///////////////////////////////////////////////////////////////////////////
  //
  //      Saving histograms
  //
  ///////////////////////////////////////////////////////////////////////////

  string dirName = "acceptCorr";
  gSystem->mkdir(dirName.data());

  std::string command = update? "update": "recreate";
  TFile* outFile = new TFile(Form("%s/acceptCorr_%s", dirName.data(),
				  filepostfix.data()),
			     command.data());       

  //   h[0]->Write();
  //   h[1]->Write();
  hscale->Write();
  outFile->Close();

  ///////////////////////////////////////////////////////////////////////////
  //
  //      Making histograms
  //
  ///////////////////////////////////////////////////////////////////////////

  TCanvas* c1 = new TCanvas("c1","",700,1000);  
  c1->Divide(1,2,0.01,0);
  c1->cd(1);
  if(logScale)
    gPad->SetLogy(1);
  gPad->SetTopMargin(0.01);
  gPad->SetBottomMargin(0);
  gPad->SetRightMargin(0.04);

  
  h[1]->Draw("he");
  h[0]->Draw("hesame");

  float x1NDC = 0.631;
  float y1NDC = 0.776;
  float x2NDC = 0.770;
  float y2NDC = 0.941;

  TLegend* leg = new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  for(int ifile=0; ifile<NFILES; ifile++)
    leg->AddEntry(h[ifile], mcName[ifile].data());
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
  hscale->SetMaximum(1.04);
  hscale->SetMinimum(0.65);
  hscale->SetTitleOffset(1.2,"Y");
  hscale->Draw("e1");
  TLine* l2 = new TLine(xmin,1.,xmax,1.);
  l2->SetLineColor(4);
  l2->SetLineStyle(3);
  l2->Draw("same");

  std::string filename;
  std::string psname = dirName + "/" + var;

  filename = psname + ".eps";
  c1->Print(filename.data());
  filename = psname + ".gif";
  c1->Print(filename.data());
  filename = psname + ".pdf";
  c1->Print(filename.data());
}
		     
