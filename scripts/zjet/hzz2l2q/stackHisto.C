#include <vector>
#include <string>
#include <iostream>
#include <fstream>

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


void stackHisto(std::string histoName,                     
		std::string inputFile="inputFile.txt", bool update=false,
		int rebin=1, double xmin=-9999.0, double xmax=-9999.0)
{
  
  //---------------------------------------------------------------------------
  // read in the MC root files
  //---------------------------------------------------------------------------

  std::vector<MCFile> myMCFiles;

  double lumi =  5204.7; // for muon, for electron: 5200.7
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

    if(flag!=-1)
      {
	myMCFiles.push_back(tempMCfile);
        cout << tempMCfile.filename << "\t" << tempMCfile.scaleFactor << endl;
      }

      
  }	 
	
  const unsigned int nfiles = myMCFiles.size();
   
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
  TH1D* h_all;
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

  h_all = (TH1D*)h_template->Clone(Form("%s_all",histoName.data()));
  h_all -> Reset();



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

  
  TCanvas* c1 = new TCanvas("c1",inputFile.data(),0,0,500,500);
  h_all->SetMarkerSize(1);
  h_all->SetMarkerStyle(24);
  h_all->Draw("e");
  //   h_all->SetXTitle(xtitle.data());
  h_all->GetXaxis()->SetRangeUser(xmin,xmax);
  cout << histoName << " integral = " << h_all->Integral() << endl;
  for(int ifile=0; ifile < nfiles; ifile++){
    h_deno[ifile]->Draw("hist,same");
  }
  
  std::string remword  =".txt";
  size_t pos  = inputFile.find(remword);

  if(pos!= std::string::npos)
    inputFile.replace(pos,remword.length(),"");

  std::string command = "recreate";

  if(update)command ="update";

  TCanvas* c2 = new TCanvas("c2",inputFile.data(),500,0,500,500);
  hs->Draw("hist");
  hs->GetXaxis()->SetRangeUser(xmin,xmax);
  hs->GetXaxis()->SetTitle(h_all->GetXaxis()->GetTitle());
  //   hs->GetXaxis()->SetTitle(xtitle.data());
  hs->Draw("hist");

  

  TFile* outFile = new TFile(Form("combined_%s.root",inputFile.data()),command.data());               
  h_all->Write();
  hs->Write();
  outFile->Close();     


}
