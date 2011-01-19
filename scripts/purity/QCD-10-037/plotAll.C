void plotAll(int fileStart=0, int fileEnd=23,
	     std::string inputRoot="projectedHistos/projected_bkgMC_pt_histo.root",
	     std::string histoPrefix="TemplateB_Eta_0_25000_",
	     std::string histoEndfix="0")
{

  int REBINNINGS=3;
  
  TFile* inf = new TFile(inputRoot.data());

  TCanvas* c3 = new TCanvas("c3");

  TH1D* hsum;

  hsum = (TH1D*)inf->FindObjectAny(Form("%s%s",
					histoPrefix,histoEndfix));
  hsum->Rebin(REBINNINGS);
  hsum->Draw();

  
  TCanvas* c2 = new TCanvas("c2");

  const int NFILE=fileEnd-fileStart+1;

  TH1D* h[100];

  TH1D* hTotal;

  cout << inputRoot << endl;
  cout << histoPrefix << endl;
  
  for(int i=fileStart; i<= fileEnd; i++){

    
    h[i] = (TH1D*)inf->FindObjectAny(Form("%s%02i_%s",
					       histoPrefix,i,histoEndfix));
    h[i] -> SetLineColor((i%6)+1);

    h[i] -> Rebin(REBINNINGS);
    
    if(i==fileStart)
      {
	hTotal = (TH1D*)h[i]->Clone();
	hTotal->Reset();
	hTotal->SetLineStyle(4);
	hTotal->SetLineColor(4);
	h[i] -> Draw();
// 	h[i] -> GetXaxis()->SetRangeUser(0,300   );
      }
    else 
      h[i] -> Draw("same");

    hTotal->Add(h[i]);

  }

  hTotal->Draw("hist");
  for(int i=fileStart; i<= fileEnd; i++)
    h[i]->Draw("same");


  TCanvas* c4 = new TCanvas("c4");
  double scale = 1.0/hsum->Integral();
  hsum->Scale(scale);
  scale = 1.0/hTotal->Integral();
  hTotal->Scale(scale);
  hsum->Draw();
  hTotal->Draw("hist,same");
}
