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
    // cout << fileN << endl;
    if(fileN.find("fail") != std::string::npos)continue;
    if(fileH->IsFolder()){
     std::string newDir=dirname+fileN;
      base->SetDirectory(newDir.data());
      TList *listOfFiles2 = base->GetListOfFiles();
      TIter fileIt2(listOfFiles2);
      TFile *fileH2 = new TFile();  
      while(fileH2 = (TFile*)fileIt2()) {
	std::string fileN2 = fileH2->GetName();
	// cout << fileN2 << endl;
	if(fileH2->IsFolder())continue;
 	if(fileN2.find("fail") != std::string::npos)continue;
 	if(fileN2.find(baseString) == std::string::npos)continue;
	// cout << fileN2.data() << endl;
	TFile *inf = TFile::Open(Form("%s/%s",newDir.data(),fileN2.data()));
	TH1F *h = (TH1F*)inf->FindObjectAny("totalEvents");
	jEntry += h->GetEntries();
	inf->Close();
      }
    } 

  }

  cout << "Total number of entries = " << jEntry << endl;


}

