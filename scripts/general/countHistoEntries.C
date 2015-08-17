#include <TROOT.h>
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include <string>
#include <iostream>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TList.h>
using namespace std;

void countHistoEntries(std::string dirname, std::string rootname="NCU")
{

  TSystemDirectory *base = new TSystemDirectory("root","root");

  base->SetDirectory(dirname.data());
  TList *listOfFiles = base->GetListOfFiles();
  TIter fileIt(listOfFiles);
  TFile *fileH = new TFile();
  Long64_t jEntry=0;
  while(fileH = (TFile*)fileIt()) {
    std::string fileN = fileH->GetName();
    std::string baseString = rootname.data();
    if( fileH->IsFolder())  continue;
    if(fileN.find(baseString) == std::string::npos)continue;
    TFile *inf = new TFile(Form("%s/%s",dirname.data(),fileN.data()));
    TH1F *h = (TH1F*)inf->FindObjectAny("totalEvents");
    jEntry += h->GetEntries();

  }

  cout << "Total number of entries = " << jEntry << endl;


}
