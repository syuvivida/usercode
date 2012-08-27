#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"

void compareAcceptance(std::string file1,
		       std::string file2,
		       std::string var, 
		       float xmin=-9999.0, float xmax=-9999.0,
		       bool update=false
		       )
{
  
  setTDRStyle();
  gStyle->SetOptStat(0);

  const int NFILES=2;

  TH1F* h[NFILES];

  std::string mcfile[NFILES]={
    file1,
    file2
  };
  
  std::string mcName[NFILES]={
    "SHERPA",
    "MADGRAPH"
  };

  TFile* fmc[NFILES];

  // first get the histogram files
  for(int ifile=0; ifile < NFILES; ifile++){

    fmc[ifile] = TFile::Open(mcfile[ifile] .data());
    cout << "opened " << fmc[ifile]->GetName() << endl;
    h[ifile]  = (TH1F*)(fmc[ifile]->Get(var.data()));
    h[ifile]->SetName(Form("h%d",ifile));

  }

  std::string remword  ="hscale_";
  size_t pos  = var.find(remword);
  if(pos!= std::string::npos)
    var.replace(pos,remword.length(),"");

  TH1D* have =(TH1D*) h[0]->Clone(Form("have_%s",var.data()));

  int COLORCODE[]={2,4};
  int MARKERCODE[]={24,21};

  for(int i=0; i < NFILES; i++){
    h[i]->GetXaxis()->SetNdivisions(5);
    h[i]->GetYaxis()->SetDecimals();
    h[i]->SetTitleOffset(1.1,"Y");
    h[i]->SetLineColor(COLORCODE[i]);
    h[i]->SetMarkerColor(COLORCODE[i]);
    h[i]->SetMarkerSize(1);
    h[i]->SetMarkerStyle(MARKERCODE[i]);
  };

  have->GetXaxis()->SetNdivisions(5);
  have->GetYaxis()->SetDecimals();

  // if normalizing to the same area, set the scale 

  int binLo = -1;
  int binHi = -1;
  int nbins = have->GetNbinsX();
  if(xmin>-9999.0 && xmax>-9999.0)
    {

      binLo = have->FindBin(xmin);
      binHi = have->FindBin(xmax)-1;

    }

  else
    {
      binLo = 1;
      binHi = nbins;
      xmin = have->GetBinLowEdge(1);
      xmax = have->GetBinLowEdge(nbins+1);
    }

  for(int i=1; i <= have->GetNbinsX(); i++)
    {
      double value[2];
      double sigma[2];

      for(int j=0; j<2; j++)
	{
	  value[j] = h[j]->GetBinContent(i);
	  sigma[j] = h[j]->GetBinError(i);
	}

      if(value[0]==0 || value[1]==0)continue;
      if(sigma[0]==0 || sigma[1]==0)continue;

      double sigmaAll = 0.0;
      double valueAll = 0.0;

      for(int j=0; j<2; j++)
	{
	  sigmaAll += 1.0/pow(sigma[j],2);
	  valueAll += value[j]/pow(sigma[j],2);
	}

      valueAll /= sigmaAll;
      sigmaAll = sqrt(1.0/sigmaAll);
      
      have->SetBinContent(i,valueAll);
      have->SetBinError(i,sigmaAll);
      
    }
  
  for(int i=1; i <= have->GetNbinsX(); i++)
    {
      cout << "Bin " << i << ": (" << 
	h[0]->GetBinContent(i) << "+-" << h[0]->GetBinError(i) << ") + (" << 
	h[1]->GetBinContent(i) << "+-" << h[1]->GetBinError(i) << ") = (" << 
	have->GetBinContent(i) << "+-" << have->GetBinError(i) << ")" << endl;
    }
	

  h[0]->GetXaxis()->SetRangeUser(xmin,xmax);
  h[1]->GetXaxis()->SetRangeUser(xmin,xmax);
  have->GetXaxis()->SetRangeUser(xmin,xmax);



  ///////////////////////////////////////////////////////////////////////////
  //
  //      Saving histograms
  //
  ///////////////////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////////////////
  //
  //      Making histograms
  //
  ///////////////////////////////////////////////////////////////////////////

  TCanvas* c1 = new TCanvas("c1","",600,500);  
  
  h[1]->Draw("e");
  h[0]->Draw("esame");
  have->Draw("e1same");

  float x1NDC = 0.198;
  float y1NDC = 0.256;
  float x2NDC = 0.3356;
  float y2NDC = 0.4216;

  TLegend* leg = new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  for(int ifile=0; ifile<NFILES; ifile++)
    leg->AddEntry(h[ifile], mcName[ifile].data());
  leg->AddEntry(have,"Average");
  leg->Draw("same");


  string dirName = "acceptCorr";
  gSystem->mkdir(dirName.data());


  std::string filename;
  std::string psname = dirName + "/ave_" + var;

  filename = psname + ".eps";
  c1->Print(filename.data());
  filename = psname + ".gif";
  c1->Print(filename.data());
  filename = psname + ".pdf";
  c1->Print(filename.data());

  std::string command = update? "update": "recreate";
  TFile* outFile = new TFile(Form("%s/ave_sherpamadgraph.root", dirName.data()),
			     command.data());       

  have->Write();
  outFile->Close();




}
		     
