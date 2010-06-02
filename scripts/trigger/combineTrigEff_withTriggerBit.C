#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TChain.h>
#include <TBranch.h>
#include <iostream.h>
#include <fstream.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystemDirectory.h>
#include <TProfile.h>					
#include <TList.h>
#include <TF1.h>
#include <vector>
#include <TLorentzVector.h>
#include "TrigEff_AsymmetryErrors.C"
#include <TGraphAsymmErrors.h>

using namespace SY_NT;
using namespace std;

void combineTrigEff_withTriggerBit()
{
  std::string dirname[100];

  int NFILES = 0;

  ifstream fin;
  fin.open("InputFile");
  while(!fin.eof())
    {
      fin >> dirname[NFILES];
      if(!fin.eof())NFILES++;
    }

  cout << "NFILES = " << NFILES << endl;
  for(int i=0; i < NFILES; i++)
    cout << dirname[i] << endl;

  TH1D* horiginal[NFILES][2];
  TH1D* hscale[NFILES][2];
  TH1D* hdeno;
  TH1D* hnumr;

  TFile* histoFile[NFILES];
  
  // assign files
  for(int i=0; i < NFILES; i++)
    {     
      std::string fileName = dirname[i] + "_HLT15L1RtoL1EG5.root";
      histoFile[i] = TFile::Open(fileName.data());
      cout << "Opened " << fileName << endl;
    }
  
  char name[300];
  // assign histograms
  for(int i=0; i< NFILES; i++)
    {

       horiginal[i][0] = (TH1D*)(histoFile[i]->Get("h_recgetdeno"));
       horiginal[i][1] = (TH1D*)(histoFile[i]->Get("h_recgetnumr"));
//        horiginal[i][0] = (TH1D*)(histoFile[i]->Get("h_getdeno"));
//        horiginal[i][1] = (TH1D*)(histoFile[i]->Get("h_getnumr"));

      ifstream fin;
      double xsec, ngen, filtereff;
      std::string location="/mc/NCUHEP/" + dirname[i] + "/geninfo";
      fin.open(location.data());
      fin >> xsec >> ngen >> filtereff;
      fin.close();
     
      double scale = 50.0*xsec*filtereff/ngen;

      cout << dirname[i] << " : " << xsec << "\t" << ngen << 
	"\t" << filtereff << endl;
      
      for(int k=0; k < 2; k++)	
	{
	  sprintf(name,"horiginal%d%d",i,k);
	  horiginal[i][k]->SetName(name);

	  hscale[i][k] = (TH1D*)horiginal[i][k]->Clone();
	  hscale[i][k] -> Sumw2();

	  sprintf(name,"hscale%d%d",i,k);
	  hscale[i][k]->SetName(name);
//       	  hscale[i][k]->Scale(scale);

	}

      cout << "File " << i << " has " << scale*horiginal[i][0]->Integral() << 
	" denominator events in 50/pb of 10 TeV data" << endl;

      cout << "File " << i << " has " << scale*horiginal[i][1]->Integral() << 
	" numerator events in 50/pb of 10 TeV data" << endl;

    }


  

   hdeno = (TH1D*)hscale[0][0]->Clone();
   hdeno->SetName("hdeno");
   hdeno->Reset();
   hnumr = (TH1D*)hscale[0][1]->Clone();
   hnumr->SetName("hnumr");
   hnumr->Reset();

   for(int i=0; i< NFILES; i++)
     {
       hdeno->Add(hscale[i][0]);
       hnumr->Add(hscale[i][1]);
     }
  
   hnumr->Draw();
   
   TGraphAsymmErrors *heff = MyDivide(hdeno,hnumr);  
   heff->GetXaxis()->SetTitle("Offline photon E_{T} [GeV]");
   heff->GetYaxis()->SetTitle("Trigger efficiency");                               
   heff->SetTitle("");
   heff->Draw("ap");


}


