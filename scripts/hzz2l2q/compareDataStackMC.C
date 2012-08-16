#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <TStyle.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1.h>
#include <THStack.h>

using namespace std;

struct MCFile{
  
  std::string filename;
  double      scaleFactor;

};


void compareDataStackMC(std::string histoName,                     
			std::string inputFile="inputFile.txt", bool update=false,
			double xmin=-9999.0, double xmax=-9999.0,int rebin=1)
{
  
  //---------------------------------------------------------------------------
  // read in the MC root files
  //---------------------------------------------------------------------------
  std::string mcName[]={
    "Dimuon Data",
    "DY+jets",
    "t#bar{t} +jets",
    "WW",
    "WZ"
  };

  std::string leptonName="test";
  double lumi =5000;

  if(inputFile.find("electron")!= std::string::npos){
    leptonName="electron";
    lumi = 5200.7;
    mcName[0]="DiElectron Data";
  }
  else if(inputFile.find("muon")!= std::string::npos){
    leptonName="muon";
    lumi = 5204.7;
  }

  cout << "Data luminosity = " << lumi << endl;

  std::vector<MCFile> myMCFiles;
  std::string dataFile;

  FILE *fTable = fopen(inputFile.data(),"r");
   
  int flag=1;   
  while (flag!=-1){
    // first reading input file
    char filename[500];
    flag=fscanf(fTable,"%s",filename);

    char tmp[1000];
    // read in x-section
    flag=fscanf(fTable,"%s",tmp);
    double cross=atof(tmp);

    // read in number of events
    flag=fscanf(fTable,"%s",tmp);
    double nevt=atof(tmp);

    double scale =lumi*cross/nevt;

    MCFile tempMCfile;
    tempMCfile.filename = filename;
    tempMCfile.scaleFactor = scale;

    if(flag!=-1 && cross>1e-6)
      {
	myMCFiles.push_back(tempMCfile);
        cout << tempMCfile.filename << "\t" << tempMCfile.scaleFactor << endl;
      }
    else if(flag!=-1) // is data
      {
	dataFile = filename;
      }
  }	 
	
  const unsigned int nfiles = myMCFiles.size();
  cout << "Data file = " << dataFile << endl;
  cout << "Reading " << nfiles << " files" << endl;

  if(nfiles==0)
    {
      cout << "There is no file for me to process!" << endl;
      return;
    }

  //-------------------------------------------------------------------------------
  // Declaring histograms to be used
  //-------------------------------------------------------------------------------

  // local histograms
  THStack *hs = new THStack(Form("%s_stack",histoName.data()),"");

  // for local debugging
  TH1D* h_deno[nfiles];
  int COLORCODE[]={
    kRed-7,
    kRed-10,
    kRed-6,
    kMagenta-2,
    kMagenta-6,
    kBlue-7,
    kBlue-9
  };

  cout << "opening " << myMCFiles[0].filename << endl;
  TFile *f1 = TFile::Open(myMCFiles[0].filename.data());

  TH1D* h_template = (TH1D*)(f1->Get(Form("%s",histoName.data())));
  h_template->Reset();
  
  TH1D* h_all = (TH1D*)h_template->Clone(Form("%s_all",histoName.data()));
  h_all -> Reset();

  TFile *fdata = TFile::Open(dataFile.data());

  TH1D* h_data= (TH1D*)(fdata->Get(Form("%s",histoName.data())));
  h_data-> SetName(Form("%s_data",histoName.data()));

  h_data->Draw();

  //-------------------------------------------------------------------------------
  // combine
  //-------------------------------------------------------------------------------


  for(int ifile=nfiles-1; ifile>=0; ifile--){

    cout << "File " << ifile << endl;
    TFile *f_temp = TFile::Open(myMCFiles[ifile].filename.data());
    cout << "opening " << myMCFiles[ifile].filename << endl;


    h_deno[ifile] = (TH1D*)(f_temp->Get(Form("%s",histoName.data())));
    h_deno[ifile] -> SetName(Form("h_deno_%d",ifile));
    h_deno[ifile] -> SetTitle(f_temp->GetName());
    h_deno[ifile] -> Rebin(rebin);
    h_deno[ifile] -> SetLineColor(COLORCODE[ifile]);
    h_deno[ifile] -> SetFillColor(COLORCODE[ifile]);
    h_deno[ifile] -> SetFillStyle(1001);
    h_deno[ifile] -> SetMarkerColor(COLORCODE[ifile]);
    h_deno[ifile] -> Sumw2();
    double weight = myMCFiles[ifile].scaleFactor;
    h_deno[ifile] -> Scale(weight);
    if(ifile==nfiles-1)
      {
	h_all    -> Rebin(rebin);
	h_all    -> Sumw2();
      }

    h_all    -> Add(h_deno[ifile]);
    hs       -> Add(h_deno[ifile]);

    // to be used with TEfficiency methods
    cout << h_deno[ifile]->GetEntries() << endl;
  } // end of loop over files

  float x1NDC = 0.620968;
  float y1NDC = 0.684322;
  float x2NDC = 0.762097;
  float y2NDC = 0.898305;

  TLegend* leg = new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->AddEntry(h_data,mcName[0].data());
  for(int i=0; i < nfiles; i++)
    leg->AddEntry(h_deno[i], mcName[i+1].data());

  int maxBin = h_data->GetMaximumBin();
  double max = h_data->GetMaximum()+h_data->GetBinError(maxBin);
  double maxmc = hs->GetMaximum();

  if(maxmc>max)max=maxmc;

  h_data->SetMaximum(1.1*max);

  TCanvas* c1 = new TCanvas("c1",inputFile.data(),0,0,500,500);
  h_data->SetMarkerSize(1);
  h_data->SetMarkerStyle(24);
  h_data->SetTitle("");
  h_data->Draw("e");
  if(xmin>-9999 && xmax>-9999){
    h_data->GetXaxis()->SetRangeUser(xmin,xmax);
  }
  cout << "Data integral = " << h_data->Integral() << endl;
  hs->Draw("histsame");
  h_data->Draw("esame");
  leg->Draw("same");
  

  std::string dirName = "compareDataMC_" + leptonName;
  gSystem->mkdir(dirName.data());

  std::string filename;
  std::string psname ;
  psname = dirName+ "/overlay_" + histoName;
  filename = psname + ".eps";
  c1->Print(filename.data());
  filename = psname + ".gif";
  c1->Print(filename.data());
  filename = psname + ".pdf";
  c1->Print(filename.data());
  

  // study the ratios

  TCanvas* c2 = new TCanvas("c2","",700,0,700,1000);  
  c2->Divide(1,2,0.01,0);
  c2->cd(1);
  gPad->SetTopMargin(0.01);
  gPad->SetBottomMargin(0);
  gPad->SetRightMargin(0.04);

  h_data->Draw("e");
  if(xmin>-9999 && xmax>-9999){
    h_data->GetXaxis()->SetRangeUser(xmin,xmax);
  }
  cout << h_data->GetName() << " integral = " << h_data->Integral() << endl;
  hs->Draw("histsame");
  h_data->Draw("esame");
  leg->Draw("same");


  c2->cd(2);
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

  TH1D* hratio = (TH1D*)h_data->Clone("hratio");
  hratio->Reset();
  hratio->Divide(h_data,h_all,1.0,1.0);
  hratio->SetTitle("");
  hratio->SetMaximum(1.5);
  hratio->SetMinimum(0.5);
  hratio->SetTitleOffset(1.2,"Y");
  hratio->GetXaxis()->SetTitle(h_data->GetXaxis()->GetTitle());
  hratio->Draw("e1");

  cout << "( " << h_data->GetBinContent(maxBin) << "+-" << h_data->GetBinError(maxBin) << " )/(" 
       << h_all->GetBinContent(maxBin) << "+-" << h_all->GetBinError(maxBin) <<  ")= " 
       << hratio->GetBinContent(maxBin) << "+-" << hratio->GetBinError(maxBin) << endl;

  psname = dirName+ "/ratio_" + histoName;
  filename = psname + ".eps";
  c2->Print(filename.data());
  filename = psname + ".gif";
  c2->Print(filename.data());
  filename = psname + ".pdf";
  c2->Print(filename.data());
  

  std::string remword  =".txt";
  size_t pos  = inputFile.find(remword);

  if(pos!= std::string::npos)
    inputFile.replace(pos,remword.length(),"");

  std::string command = "recreate";

  if(update)command ="update";

  TFile* outFile = new TFile(Form("combined_%s.root",inputFile.data()),command.data());               
  h_all->Write();
  hs->Write();
  h_data->Write();
  outFile->Close();     


}
