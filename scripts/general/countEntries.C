#include <TFile.h>
#include <TChain.h>
#include <TBranch.h>
#include <iostream>
#include <TSystemDirectory.h>
#include <TList.h>

void countEntries(std::string dirname, std::string rootname="MultiPhotonAnalyzer",
		  std::string treename="NTuples/Analysis")
{

  // chain in all the ncu ntuples in the same directory
  TChain* pho = new TChain(treename.data());

  TSystemDirectory *base = new TSystemDirectory("root","root");
  std::string filename =  dirname;

  base->SetDirectory(filename.data());
  TList *listOfFiles = base->GetListOfFiles();
  TIter fileIt(listOfFiles);
  TFile *fileH = new TFile();
  int nfile=0;
  while(fileH = (TFile*)fileIt()) {
    std::string fileN = fileH->GetName();
    std::string baseString = rootname.data();
    if( fileH->IsFolder())  continue;
    if(fileN.find(baseString) == std::string::npos)continue;
    cout << fileN.data() << endl;
    nfile++;
    pho->Add(Form("%s/%s",dirname.data(),fileN.data()));
  }

  std::cout << "Opened " << nfile << " files" << std::endl;
  
  
  // histogram definition
  Long64_t nentries = (Long64_t)pho->GetEntries();
  cout << "The total number of entries is " << nentries << endl;
  

}


