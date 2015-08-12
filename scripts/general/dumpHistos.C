#include <TROOT.h>
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include <string.h>
#include <iostream.h>
#include <TSystem.h>

void dumpHistos(std::string inputFile_)
{

  TString endfix=gSystem->GetFromPipe(Form("file=%s; test=${file##*Diboson/}; echo \"${test%%_narrow.root*}\"",inputFile_.data()));

  cout << endfix << endl;
  std::string dirName = "/afs/cern.ch/user/s/syu/www/public/xtozh/priority0_histos/" + endfix;
  gSystem->mkdir(dirName.data());

  TCanvas* c1 = new TCanvas("c1","",500,500);

  TFile *_file = TFile::Open(inputFile_.data());
  _file->cd();
  TDirectory *current_sourcedir = gDirectory;

  // loop over all keys in this directory
  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key;
  while (key = (TKey*)nextkey() ) {

    TObject *obj = key->ReadObj();
    if ( obj->IsA()->InheritsFrom( "TH1" ) ) {
      // descendant of TH1 -> scale it

      cout << "outputing histogram " << obj->GetName() << endl;
      TH1 *h1 = (TH1*)obj;
      h1->SetMarkerStyle(8);
      h1->SetMarkerSize(1);
      h1->SetLineWidth(3);
      h1->SetLineColor(4);
      h1->Draw();
      c1->Print(Form("%s/%s.gif",dirName.data(),obj->GetName()));
    } // if the object is a histogram
  } // loop over keys



}
