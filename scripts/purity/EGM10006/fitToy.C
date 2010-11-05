#include <TFile.h>
#include <TChain.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <iostream>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TKey.h>
#include <TH1.h>

using namespace std;


void fitToy(std::string dirname, std::string baseString="fixsigtail")
{


  TSystemDirectory *base = new TSystemDirectory("root","root");

  base->SetDirectory(dirname.data());
  TList *listOfFiles = base->GetListOfFiles();
  TIter fileIt(listOfFiles);
  TFile *fileH = new TFile();
  int nfile=0;
    cout << "Hello " << endl;
  while(fileH = (TFile*)fileIt()) {
    cout << "Hello " << endl;
    std::string fileN = fileH->GetName();

    cout << "Hello " << endl;
    if( fileH->IsFolder())  continue;
    if(fileN.find(baseString) == std::string::npos)continue;

    cout << fileN << endl;
  // rename the output file
   std::string remword=".root";
   size_t pos = fileN.find(remword);
   std::string forOutput = fileN;  
   if(pos!= std::string::npos)
     forOutput.swap(forOutput.erase(pos,remword.length()));   
   TFile* _file = TFile::Open(fileN.data());  
    
  TDirectory *current_sourcedir = gDirectory;

  // loop over all keys in this directory
  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key;
  cout << "Ma" << endl;
  while (key = (TKey*)nextkey() ) {
    cout << "hey" << endl;
    TObject *obj = key->ReadObj();
    if ( obj->IsA()->InheritsFrom( "TH1" ) ) {
      // descendant of TH1 -> scale it

      cout << "fitting histogram " << obj->GetName() << endl;

      TCanvas* c1 = new TCanvas("c1");      
      TH1 *h1 = (TH1*)obj;

      h1->Fit("gaus");


      c1->Print(Form("%s_%s.png",forOutput.data(),key->GetName()));
      delete c1;
      
    } // if the object is a histogram
  } // loop over keys




  }






}
