#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"

void compareGenMore(
		    std::string var="h_mH_parton", 
		    float xmin=-9999.0, float xmax=-9999.0,
		    bool logScale=true,
		    std::string mcfile1="histo_ggH_mc_cut_900.root",
		    std::string mcfile2="histo_ggH_mc_cut_525.root",
		    std::string mcfile3="histo_ggH_mc_cut_300.root",
		    std::string mcfile4="histo_ggH_mc_cut_200.root",
		    std::string output="test",
		    std::string headertitle="H^{0}#rightarrow ZZ#rightarrow 2l2q",
		    std::string mcName1="M_{H}=900 GeV/c^{2}",
		    std::string mcName2="M_{H}=525 GeV/c^{2}",
		    std::string mcName3="M_{H}=300 GeV/c^{2}",
		    std::string mcName4="M_{H}=200 GeV/c^{2}"
		    )
{
  
  setTDRStyle();
  gStyle->SetOptStat(0);

  const int NHISTOS=4;
  TH1F* h[NHISTOS];

  char tempName[300];
  
  // first get the histogram files
  TFile *fmc1   = TFile::Open(mcfile1.data());
  TFile *fmc2   = TFile::Open(mcfile2.data());
  TFile *fmc3   = TFile::Open(mcfile3.data());
  TFile *fmc4   = TFile::Open(mcfile4.data());

  h[0]    = (TH1F*)(fmc1->Get(var.data()));
  h[1]    = (TH1F*)(fmc2->Get(var.data()));
  h[2]    = (TH1F*)(fmc3->Get(var.data()));
  h[3]    = (TH1F*)(fmc4->Get(var.data()));

  TH1D* hscale[NHISTOS];

  int COLOR[NHISTOS]={4,2,kOrange-1,1};
  int MARKERSTYLE[NHISTOS]={24,21,29,8};

  for(int i=0; i < NHISTOS; i++){

    hscale[i]   =(TH1D*) h[0]->Clone(Form("hscale%02i",i));
    hscale[i]   ->SetYTitle(Form("Ratio to %s",mcName3.data()));
    
    hscale[i]   ->GetXaxis()->SetNdivisions(5);
    hscale[i]   ->GetYaxis()->SetDecimals();

    hscale[i]->SetLineColor(COLOR[i]);
    hscale[i]->SetMarkerColor(COLOR[i]);
    hscale[i]->SetMarkerStyle(MARKERSTYLE[i]);

    hscale[i]->SetTitle("");
    hscale[i]->SetMaximum(2.3);
    hscale[i]->SetMinimum(0.1);
    hscale[i]->SetTitleOffset(1.2,"Y");
    hscale[i]->SetTitleOffset(1.2,"X");

    h[i]->SetTitle("");
    h[i]->GetXaxis()->SetNdivisions(5);
    h[i]->GetYaxis()->SetDecimals();

    h[i]->SetLineColor(COLOR[i]);
    h[i]->SetMarkerColor(COLOR[i]);
    h[i]->SetMarkerSize(1);
    h[i]->SetMarkerStyle(MARKERSTYLE[i]);
    h[i]->SetTitleOffset(1.2,"Y");

  }

  // if normalizing to the same area, set the scale 

  int binLo = -1;
  int binHi = -1;
  int nbins = h[0]->GetNbinsX();
  if(xmin>-9999.0 && xmax>-9999.0)
    {

      binLo = h[0]->FindBin(xmin);
      binHi = h[0]->FindBin(xmax)-1;

    }

  else
    {
      binLo = 1;
      binHi = nbins;
      xmin = h[0]->GetBinLowEdge(1);
      xmax = h[0]->GetBinLowEdge(nbins+1);
    }


  float scaleFactor[NHISTOS]={1};

  for(int ih=0; ih < NHISTOS; ih++){
    
    float integral = h[ih]->Integral(binLo,binHi);
    scaleFactor[ih] = integral > 0? 1.0/integral: 1;
    h[ih]->Sumw2();
    h[ih]->Scale(scaleFactor[ih]);

  }

  for(int ih=0; ih < NHISTOS; ih++)
    cout << "histogram " << ih << " integral = " << h[ih]->Integral() << endl;

  for(int ih=0; ih < NHISTOS; ih++)
    cout << "histogram " << ih << " integral = " << h[ih]->Integral(binLo, 
								    binHi)
	 << endl;


  // get the ratio

  for(int ih=0; ih < NHISTOS-1; ih++){
    cout << "===================================================" << endl;
    cout << "For histogram " << ih << endl;
    for(int ib=1;ib<= nbins;ib++){

      double nref    =h[NHISTOS-1]->GetBinContent(ib);
      double nreferr =h[NHISTOS-1]->GetBinError(ib);

      double ncomp   =h[ih]->GetBinContent(ib); 
      double ncomperr=h[ih]->GetBinError(ib); 

      double ratio = -9999;
      double err   = 1e-4;

      hscale[ih]->SetBinContent(ib,ratio);
      hscale[ih]->SetBinError(ib,err);
    
      if(ncomp<=0 || nref<=0)continue;    
    
      if(ncomperr<=0 || nreferr<=0)
	continue;
      
      // now calculate the ratio 1
      ratio = ncomp/nref;
      err =
	(ratio)*sqrt(pow(nreferr/nref,2)+pow(ncomperr/ncomp,2));

      hscale[ih]->SetBinContent(ib,ratio);
      hscale[ih]->SetBinError(ib,err);

      cout << "Bin " << ib << " ratio = " << ratio << " +- " << err << endl;

      cout << "===================================================" << endl;
    } // end of loop over bins
  } // end of loop over histograms

  vector<float> maxArray;
  maxArray.clear();

  for(int ih=0; ih < NHISTOS; ih++){

    h[ih]->GetXaxis()->SetRangeUser(xmin,xmax);
    hscale[ih]->GetXaxis()->SetRangeUser(xmin,xmax);
    float max_this  = h[ih]->GetBinError(h[ih]->GetMaximumBin()) + h[ih]->GetMaximum();
    maxArray.push_back(max_this);

  }


  float max = *max_element(maxArray.begin(),maxArray.end());
  cout << "Max = " << max << endl;

  for(int ih=0; ih < NHISTOS; ih++){

    h[ih]->SetMaximum(1.1*max);

  }

  cout << "here" << endl;
  TCanvas* c1 = new TCanvas("c1","",700,1000);  
  c1->Divide(1,2,0.01,0);
  c1->cd(1);
  if(logScale)
    gPad->SetLogy(1);
  gPad->SetTopMargin(0.01);
  gPad->SetBottomMargin(0);
  gPad->SetRightMargin(0.04);

  
  cout << "here1" << endl;
  h[0]->Draw("e");
  for(int ih=0; ih < NHISTOS-1; ih++)
    h[ih]->Draw("esame");

  cout << "here2" << endl;

  h[NHISTOS-1]->Draw("hesame");

  cout << "here3" << endl;
  float x1NDC = 0.697883;
  float y1NDC = 0.656116;
  float x2NDC = 0.918347;
  float y2NDC = 0.972577;

  TLegend* leg = new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->SetHeader(headertitle.data());
  leg->AddEntry(h[0], mcName1.data());
  leg->AddEntry(h[1], mcName2.data());
  leg->AddEntry(h[2], mcName3.data());
  leg->AddEntry(h[3], mcName4.data());
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

  cout << "here4" << endl;

  hscale[0]->Draw("e1");
  for(int ih=0; ih < NHISTOS-1; ih++){
    hscale[ih]->Draw("e1same");
  }

  cout << "here5" << endl;
  TLine* l2 = new TLine(xmin,1.,xmax,1.);
  l2->SetLineColor(4);
  l2->SetLineStyle(3);
  l2->Draw("same");


  string dirName = "compareGenMore";
  gSystem->mkdir(dirName.data());

  std::string filename;
  std::string psname = dirName + "/" + var;
  if(output !="test")
    psname = dirName+ "/" + output;
  else
    psname = dirName+ "/" + var;
  filename = psname + ".eps";
  c1->Print(filename.data());
  filename = psname + ".gif";
  c1->Print(filename.data());
  filename = psname + ".pdf";
  c1->Print(filename.data());
  //   c1->Close();
}
		     
