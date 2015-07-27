#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TKey.h>
#include <string>
#include <iostream>
#include <TChain.h>

using namespace std;

void loopOverBranches(std::string inputFile_,std::string treeName="tree/treeMaker")
{

  TString endfix=gSystem->GetFromPipe(Form("file=%s; echo \"${file%%.root*}\"",inputFile_.data()));
  endfix += ".txt";
  cout << endfix << endl;
  TChain* fTree = new TChain(treeName.data());
  fTree->Add(inputFile_.data());
  ofstream fout; 
  fout.open(endfix.Data());
  TObjArray* all=fTree->GetListOfBranches();

  for(unsigned int n=0; n < all->GetSize(); n++){
    
    TBranch* thisBranch = (TBranch*)all->At(n);
    cout << "branch " << n << " = " << thisBranch->GetName() << endl;
    fTree->Draw(Form("%s>>h1",thisBranch->GetName()));
    cout << h1->GetMean() << " " << h1->GetRMS() << endl;
    fout << thisBranch->GetName() << " " << h1->GetMean() << " " << h1->GetRMS() << endl;
  }
  
  fout.close();

}
