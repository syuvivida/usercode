#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <string>

using namespace std;

void separateTree(std::string inputFile_,int mode=0)
{

  // get the correct tree from input file
   TFile *inf = new TFile(inputFile_.data());
   TTree *fChain = (TTree*)inf->FindObjectAny("Analysis");
   cout << "Input file is " << inputFile_ << endl;

   // rename the output file
   std::string remword=".root";
   size_t pos = inputFile_.find(remword);
   std::string forOutput = inputFile_;  
   if(pos!= std::string::npos)
     forOutput.swap(forOutput.erase(pos,remword.length()));   
   std::string endfix = mode == 0 ? "_data.root" : "_template.root";
   std::string outputFile = forOutput + endfix;

   // now open new root file
   TFile* newfile_data = new TFile(outputFile.data(),"recreate");
   cout << "Output file " << outputFile << endl;

   // clone tree
   TTree* newtree = fChain->CloneTree(0);
   newtree->SetMaxTreeSize(4000000000);
   cout << "Saving "  << endfix << " tree" << endl;

   Long64_t nentries = fChain->GetEntries();
   cout << "nentries = " << nentries << endl;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
//       cout << jentry << endl;
     fChain->GetEntry(jentry);
      // if (Cut(ientry) < 0) continue;
      if(jentry%2==mode)
	newtree->Fill();
   }

 
   newtree->Print();
   newtree->AutoSave();
   delete newfile_data;


}
