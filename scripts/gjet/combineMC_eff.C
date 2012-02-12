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


void combineMC_eff(std::string histoName, std::string xtitle, int rebin=1, double xmin=-9999.0, double xmax=-9999.0)
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
  TH1D* h_all_numr;
  TH1D* h_all_deno;
  TH1D* h_eff;

  // for local debugging
  TH1D* h_numr[nfiles];
  TH1D* h_deno[nfiles];

  cout << "opening " << myMCFiles[0].filename << endl;
  TFile *f1 = TFile::Open(myMCFiles[0].filename.data());

  TH1D* h_template = (TH1D*)(f1->Get(Form("%s_0",histoName.data())));
  h_template->Reset();

  h_all_numr = (TH1D*)h_template->Clone("h_all_numr");
  h_all_deno = (TH1D*)h_template->Clone("h_all_deno");
  h_eff      = (TH1D*)h_template->Clone("h_eff");
  h_eff      -> SetXTitle(xtitle.data());
  h_eff      -> SetYTitle("Efficiency");
  h_eff      -> SetTitle("Darko method");
  h_eff      -> SetLineWidth(2);

  if(h_all_numr->GetEntries()!=0)
    {
      cout << "h_all_numr->GetEntries()!=0" << endl;
      return;
    }

  if(h_all_deno->GetEntries()!=0)
    {
      cout << "h_all_deno->GetEntries()!=0" << endl;
      return;
    }


  //-------------------------------------------------------------------------------
  // Prepare TEfficiency
  //-------------------------------------------------------------------------------

  TEfficiency* eff[nfiles];
  
  double weight_temp[nfiles];

  for(unsigned int ifile=0; ifile< nfiles; ifile++){

    cout << "File " << ifile << endl;
    TFile *f_temp = TFile::Open(myMCFiles[ifile].filename.data());
    double weight = myMCFiles[ifile].scaleFactor;
    weight_temp[ifile] = weight;

    // numerator
    h_numr[ifile] = (TH1D*)(f_temp->Get(Form("%s_1",histoName.data())));
    h_numr[ifile] -> SetName(Form("h_numr_%d",ifile));
    h_all_numr    -> Add(h_numr[ifile], weight);

    h_numr[ifile] -> Rebin(rebin);
    if(ifile==0)
      h_all_numr    -> Rebin(rebin);


    // denominator
    h_deno[ifile] = (TH1D*)(f_temp->Get(Form("%s_0",histoName.data())));
    h_deno[ifile] -> SetName(Form("h_deno_%d",ifile));
    h_all_deno    -> Add(h_deno[ifile], weight);
    h_deno[ifile] -> Rebin(rebin);
    if(ifile==0)
      h_all_deno    -> Rebin(rebin);

    // to be used with TEfficiency methods
    cout << h_numr[ifile]->GetEntries() << "\t" << h_deno[ifile]->GetEntries() << endl;
    eff[ifile] = new TEfficiency(*(h_numr[ifile]),*(h_deno[ifile]));
    eff[ifile]->SetTitle(xtitle.data());
  } // end of loop over files


  TEfficiency* eff_combine = new TEfficiency();

  TList* effList = new TList();
  for(int i=0; i<nfiles; i++)
    effList->Add(eff[i]);
  TGraphAsymmErrors* eff_final = eff_combine->Combine(effList, "v", nfiles, weight_temp);
  eff_final->GetXaxis()->SetTitle(xtitle.data());
  eff_final->GetYaxis()->SetTitle("Efficiency");
  eff_final->SetTitle("TEfficiency");

  h_eff->Rebin(rebin);
  int nbins = h_eff->GetNbinsX();
  cout << "There are " << nbins << " bins" << endl;

  double* y_eff;
  double* y_eff_h;
  double* y_eff_l;

  y_eff    = eff_final->GetY();
  y_eff_h  = eff_final->GetEYhigh();
  y_eff_l  = eff_final->GetEYlow();

  //-------------------------------------------------------------------------------
  // Darko's efficiency calculation
  //-------------------------------------------------------------------------------

  for(int ibin=1; ibin<= nbins; ibin++)
    {
      if(h_all_deno->GetBinContent(ibin) < 1e-6)continue;
      if(h_all_numr->GetBinContent(ibin) < 1e-6)continue;
      double central_value =h_all_numr->GetBinContent(ibin)/h_all_deno->GetBinContent(ibin);
      h_eff->SetBinContent(ibin, central_value);
    }

  // Darko's efficiency error calculation
  for(int ibin=1; ibin<= nbins; ibin++)
    {
      if(h_all_deno->GetBinContent(ibin) < 1e-6)continue;
      if(h_all_numr->GetBinContent(ibin) < 1e-6)continue;
      double central_value = h_eff->GetBinContent(ibin);      
      double sigma2_numerator = 0;

      for(int ifile=0; ifile < nfiles; ifile++)
	{
	  double ntotal =  h_deno[ifile]->GetBinContent(ibin);
	  if(ntotal < 10)continue;
	  double npass  =  h_numr[ifile]->GetBinContent(ibin);
	  double nfail  =  ntotal -npass;
	  double weight = myMCFiles[ifile].scaleFactor;
	  sigma2_numerator += weight*weight*((1.-central_value)*(1.-central_value)*npass
					     + central_value*central_value*nfail);
					     
	}
      double sigma = sqrt(sigma2_numerator)/h_all_deno->GetBinContent(ibin);
      h_eff->SetBinError(ibin, sigma);
    }


  // compare the two methods

  //-------------------------------------------------------------------------------
  // Final efficiency comparison and display
  //-------------------------------------------------------------------------------
  

  if(xmin > -9999.0 && xmax > -9999.0)
    {
      h_eff->GetXaxis()->SetRangeUser(xmin,xmax);
      eff_final->GetXaxis()->SetRangeUser(xmin,xmax);
    }

  h_eff->GetYaxis()->SetRangeUser(0,1.0);
  eff_final->GetYaxis()->SetRangeUser(0,1.0);
  h_eff->GetXaxis()->SetNdivisions(5);

//    for(int ibin=1; ibin<= nbins; ibin++)
//      {
//        cout << "Efficiency bin " << ibin << ": " << eff_temp[0].GetEfficiency(ibin) 
// 	    << " -" << eff_temp[0].GetEfficiencyErrorLow(ibin)
//  	   << " +" << eff_temp[0].GetEfficiencyErrorUp(ibin) << "\t" 
//  	   << h_eff->GetBinContent(ibin) << " +- " << h_eff->GetBinError(ibin) 
//  	   << endl;
 
//      }

   std::string remword="h_";

   size_t pos = histoName.find(remword);
  
   if(pos!= std::string::npos)
     histoName.swap(histoName.erase(pos,remword.length()));


   gSystem->mkdir("effDataFiles");
   ofstream fout;
   fout.open(Form("effDataFiles/%s.dat",histoName.data()));
   for(int ibin=1; ibin<= nbins; ibin++)
     {
       cout << "Efficiency bin " << ibin << ": " << y_eff[ibin-1]
	    << " -" << y_eff_l[ibin-1]
	    << " +" << y_eff_h[ibin-1] << "\t"
	    << h_eff->GetBinContent(ibin) << " +- " << h_eff->GetBinError(ibin) 
	    << endl;

       fout << ibin << " " << y_eff[ibin-1]
	    << " " << y_eff_l[ibin-1]
	    << " " << y_eff_h[ibin-1] << endl;
 
     }

   fout.close();

   TCanvas* c1 = new TCanvas("c1","",1000,500);
   c1->Divide(2,1);
   c1->cd(1);
   eff_final->Draw("ap");
   c1->cd(2);
   h_eff->Draw("e1");

   gSystem->mkdir("effHistos");
   c1->Print(Form("effHistos/eff_%s.eps",histoName.data()));
   c1->Print(Form("effHistos/eff_%s.pdf",histoName.data()));
   c1->Print(Form("effHistos/eff_%s.gif",histoName.data()));

//    if(nfiles>1){
//      TCanvas* c2 = new TCanvas("c2","",1000,1000);
//      c2->Divide(2,2);
//      for(int i=0;i<3;i++){
//        c2->cd(i+1);
//        eff[i]->Draw();
//      }
//      c2->cd(4);
//      eff_final->Draw("ap");

//      c2->Print(Form("effHistos/eff2_%s.eps",histoName.data()));
//      c2->Print(Form("effHistos/eff2_%s.pdf",histoName.data()));
//      c2->Print(Form("effHistos/eff2_%s.gif",histoName.data()));
     

//    }





}
