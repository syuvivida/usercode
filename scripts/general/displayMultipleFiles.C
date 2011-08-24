#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"
#include <vector>
#include <string>

void displayMultipleFiles(std::string inputTextFile, 
			  std::string var, 
			  std::string xtitle="", std::string ytitle="",
			  std::string output="test", bool logScale=false)
{
  
  setTDRStyle();
  gStyle->SetOptStat(0);

  
  FILE *fTable = fopen(inputTextFile.data(),"r");
   
  int flag=1;   

  std::vector<std::string> dataFile;
  std::vector<std::string> legendName;

  char filename[300];

  while (flag!=-1){
    // first reading input file
    flag=fscanf(fTable,"%s",filename);
    std::string tempFile = filename;

    flag=fscanf(fTable,"%s",filename);
    std::string tempLegend = filename;

    if (flag!=-1) {
      dataFile.push_back(tempFile);
      cout << "Reading " << tempFile.data() << endl;
      legendName.push_back(tempLegend);
      cout << "Legend is " << tempLegend.data() << endl;
    }
  }	 

  int nfile=dataFile.size();
	
  cout << "There are " << nfile << " files to compare" << endl;

  TH1F* h[20];
  TFile* f[20];

  const int NMAXFILE=10;

  if(nfile > NMAXFILE){
    cout << "Too many files!!" << endl;
    return;
  }

  TCanvas* c1 = new TCanvas("c1","",500,500);
  if(logScale)c1->SetLogy(1);
  else c1->SetLogy(0);

  float max=-9999.0;
  int maxHisto=-1;
  
  for(int i=0; i < nfile; i++){
    f[i] = TFile::Open(dataFile[i].data());
    h[i] = (TH1F*)(f[i]->Get(var.data()));
    h[i]->SetXTitle(xtitle.data());  
    h[i]->SetYTitle(ytitle.data());  
    h[i]->GetYaxis()->SetDecimals();
    h[i]->GetYaxis()->SetNdivisions(5);
    h[i]->GetXaxis()->SetNdivisions(5);
    h[i]->SetLineColor(1+i);
    h[i]->SetMarkerColor(1+i);
    h[i]->SetMarkerSize(1);
    h[i]->SetMarkerStyle(21);
    h[i]->Rebin(2);

    cout << "h[" << i << "] mean = " << h[i]->GetMean() << " and width = " << h[i]->GetRMS() << 
      " and entries = " << h[i]->GetEntries() << " and integral = " << 
      h[i]->Integral() << endl;
    

    float scale = 1.0/(float)h[i]->Integral();
    h[i]->Sumw2();
    h[i]->Scale(scale);

    float max1   = h[i]->GetBinError(h[i]->GetMaximumBin()) + h[i]->GetMaximum();

    if(max1>max){max = max1; maxHisto=i;}
  }

//   h[maxHisto]->SetMaximum(0.15);
  h[maxHisto]->SetMaximum(0.1);
  h[maxHisto]->Draw("histe");
  for(int i=0; i<nfile; i++)h[i]->Draw("histesame");
    

  TLegend* leg = new TLegend(0.31,0.425,0.51,0.85);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->SetHeader("|y_{B}| < 0.5");
  for(int i=0; i< nfile; i++)
    leg->AddEntry(h[i], legendName[i].data());
  leg->Draw("same");

  std::string outputFile;
  std::string psname = output;
  outputFile = psname + ".eps";
  c1->Print(outputFile.data());
  outputFile = psname + ".gif";
  c1->Print(outputFile.data());
}
		     
