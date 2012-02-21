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
#include <TEfficiency.h>
#include <TList.h>
#include <TGraphAsymmErrors.h>


using namespace std;

struct MCFile{
  
  std::string filename;
  double      scaleFactor;

};


void combineMC_spectrum(std::string histoName, std::string xtitle, int rebin=1, double xmin=-9999.0, double xmax=-9999.0)
{
  
  //-------------------------------------------------------------------------------
  // read in the MC root files
  //-------------------------------------------------------------------------------
  

  std::vector<MCFile> myMCFiles;

  double lumi = 1.0;
  FILE *fTable = fopen("inputFile.txt","r");
   
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
	
  const int nfiles = myMCFiles.size();
   
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

  // for local debugging
  TH1D* h_deno[nfiles];

  cout << "opening " << myMCFiles[0].filename << endl;
  TFile *f1 = TFile::Open(myMCFiles[0].filename.data());

  TH1D* h_template = (TH1D*)(f1->Get(Form("%s",histoName.data())));
  h_template->Reset();

  h_all = (TH1D*)h_template->Clone("h_all");

  h_all -> Reset();

  //-------------------------------------------------------------------------------
  // combine
  //-------------------------------------------------------------------------------

  

  for(unsigned int ifile=0; ifile< nfiles; ifile++){

    cout << "File " << ifile << endl;
    TFile *f_temp = TFile::Open(myMCFiles[ifile].filename.data());
    cout << "opening " << myMCFiles[ifile].filename << endl;


    h_deno[ifile] = (TH1D*)(f_temp->Get(Form("%s",histoName.data())));
    h_deno[ifile] -> SetName(Form("h_deno_%d",ifile));
    h_deno[ifile] -> Rebin(rebin);
    h_deno[ifile] -> SetLineColor(2*ifile+2);
    h_deno[ifile] -> Sumw2();
    double weight = myMCFiles[ifile].scaleFactor;
    h_deno[ifile] -> Scale(weight);
    if(ifile==0)
      {
	h_all    -> Rebin(rebin);
	h_all    -> Sumw2();
      }

    h_all    -> Add(h_deno[ifile]);


    // to be used with TEfficiency methods
    cout << h_deno[ifile]->GetEntries() << endl;
  } // end of loop over files

  
  TCanvas* c1 = new TCanvas("c1","After combining all MC samples",500,500);
  h_all->SetMarkerSize(1);
  h_all->SetMarkerStyle(24);
  h_all->Draw("e");
  h_all->SetXTitle(xtitle.data());
  h_all->GetXaxis()->SetRangeUser(xmin,xmax);
  for(unsigned int ifile=0; ifile < nfiles; ifile++){
    h_deno[ifile]->Draw("hist,same");
  }
  


}
