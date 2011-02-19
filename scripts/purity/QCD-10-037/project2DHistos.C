#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <string>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include <TSystem.h>

#include <string>
#include <iostream>

using namespace std;

// const double fBinsPt[]={20,30,50,80,120,180,240,500};
// const double fBinsPt[]={21.,23.,26.,30.,35.,40.,45.,50.,60.,85.,120.,300};

const double fBinsPt[]={25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 100, 120, 200,300,400,500,600,700,800,900,1000};
const int nPtBin = sizeof(fBinsPt)/sizeof(fBinsPt[0])-1;

void project2DHistos(std::string inputFile_)
{

  //Create directory, output file
  // dump histogram to a root file
  std::string outputFile_ = "projected_" + inputFile_;
  gSystem->mkdir("projectedHistos");
  Char_t strFullPath[500];
  sprintf(strFullPath,"projectedHistos/%s",outputFile_.data());

  TFile* fOutFile=new TFile(strFullPath,"RECREATE");

  cout << "Target path: " << fOutFile->GetPath() << endl;
  TString path( (char*)strstr( fOutFile->GetPath(), ":" ) );
  path.Remove( 0, 2 );

  TFile *_file = TFile::Open(inputFile_.data());
  _file->cd(path);

  TDirectory *current_sourcedir = gDirectory;

  // loop over all keys in this directory
  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key;
  while (key = (TKey*)nextkey() ) {

    TObject *obj = key->ReadObj();
    if ( obj->IsA()->InheritsFrom( "TH2" ) ) {

      // obtain bins for pt projections
  
      int binIndex[nPtBin]={0};
      TH2 *h2 = (TH2*)obj;

      for(int ibin=0; ibin<= nPtBin; ibin++)
	{
	  binIndex[ibin] = 
	    h2->GetYaxis()->FindBin(fBinsPt[ibin]);
	}

      for(int ibin=0; ibin< nPtBin; ibin++)
	{
	  cout << "Projecting 2D histogram " << h2->GetName() << " onto pt=" 
	       << fBinsPt[ibin] << " --- " << fBinsPt[ibin+1] << endl;

	  TH1 *h1 = (TH1*)h2->ProjectionX(Form("%s_%d",key->GetName(),
						(int)fBinsPt[ibin]),
					   binIndex[ibin], binIndex[ibin+1]-1);     
	  h1->SetName(Form("%s_%d",key->GetName(),(int)fBinsPt[ibin]));
	  cout << "projected histogram name: "<< h1->GetName() << " with"
	       << " pt bins = " << binIndex[ibin] << "--" << 
	    binIndex[ibin+1]-1 << endl;


	  fOutFile->cd();
	  h1->Write(h1->GetName() );

	} // end of loop over pt bins
      
    } // if the object is a histogram
  } // loop over keys


  //  fOutFile->Write();
  fOutFile->Close();









}
