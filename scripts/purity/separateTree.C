#define separateTree_cxx
#include "separateTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void separateTree::Loop(int mode)
{
   if (fChain == 0) return;

   std::string remword=".root";
   size_t pos = inputFile_.find(remword);
  
   if(pos!= std::string::npos)
     inputFile_.swap(inputFile_.erase(pos,remword.length()));   
   std::string endfix = mode == 0 ? "_data.root" : "_template.root";
   std::string outputFile_data = inputFile_ + endfix;
   cout << "Saving "  << endfix << " tree" << endl;
   TFile* newfile_data = new TFile(outputFile_data.data(),"recreate");
   cout << "Saving "  << endfix << " tree" << endl;

   Long64_t nentries = fChain->GetEntries();
   TTree* newtree = fChain->CloneTree(0);

   Long64_t nbytes = 0, nb = 0;
   cout << "nentries = " << nentries << endl;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if(jentry%2==mode)
	newtree->Fill();
   }


   newtree->Print();
   newtree->AutoSave();
   delete newfile_data;


}
