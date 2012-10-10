#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <string>
#include <TROOT.h>
#include <algorithm>
#include <vector>


using namespace std;

class runInfo
{
public:
  int run;
  int evt;
  bool operator==( const runInfo& other ) const 
  {
    return ((this->run== other.run) && 
	    (this->evt== other.evt));
  };

};


void skimTree(std::string inputFile_)
{

  
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(inputFile_.data());
  if (!f) {
    f = new TFile(inputFile_.data());
    f->cd(Form("%s:/tree",inputFile_.data()));
  }
  TTree* fChain = (TTree*)gDirectory->Get("tree");
  cout << "Input file is " << inputFile_ << endl;

  // rename the output file
  std::string remword=".root";
  size_t pos = inputFile_.find(remword);
  std::string forOutput = inputFile_;  
  if(pos!= std::string::npos)
    forOutput.swap(forOutput.erase(pos,remword.length()));   
  std::string endfix = "_filteredtree.root";
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

  Int_t           eventNo=-1;
  Int_t           runNo=-1;

  vector<int>*     trigResults;
  vector<string>*  trigName;


  bool isData=false;
  if(forOutput.find("DoubleElectron")!= std::string::npos || forOutput.find("DoubleMu")!= std::string::npos)
    isData=true;


  fChain->SetBranchAddress("EvtInfo_EventNum",&eventNo);
  fChain->SetBranchAddress("EvtInfo_RunNum",&runNo);
  if(isData){
    fChain->SetBranchAddress("trigResults", &trigResults);
    fChain->SetBranchAddress("trigName", &trigName);
  }

  Long64_t nlines=0;

  vector<runInfo> myList;

  ifstream fin;

  if(forOutput.find("DoubleElectron")!= std::string::npos)
    fin.open("/data4/syu/52X_533_validation/DoubleElectron_common.txt");
  else if(forOutput.find("DoubleMu")!= std::string::npos)
    fin.open("/data4/syu/52X_533_validation/DoubleMu_common.txt");
  else if(forOutput.find("GluGlu")!= std::string::npos)
    fin.open("/data4/syu/52X_533_validation/MC_common.txt");
  runInfo tempInfo;
  fin >> tempInfo.run >> tempInfo.evt;
  while(!fin.eof())
    {
      nlines++;
      myList.push_back(tempInfo);
      fin >> tempInfo.run >> tempInfo.evt;
    }
  fin.close();

  cout << "There are " << nlines << " lines" << endl;

  ofstream fout;
  fout.open(Form("%s_updated",forOutput.data()));

  Long64_t nPassEvt=0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    
    eventNo = -1;
    runNo = -1;
    trigName = 0;
    trigResults = 0;
    
    fChain->GetEntry(jentry);


    bool passTrigger=false;
    if(isData){

      for(int it=0; it< trigName->size(); it++)
	{
	  std::string thisTrig= trigName->at(it);
	  int results = trigResults->at(it);
//   	  cout << thisTrig << "\t" << results << endl;

	  if(forOutput.find("DoubleElectron")!= std::string::npos && 
	     thisTrig.find("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL")
	     != std::string::npos && results==1)
	    {
	      passTrigger=true;
// 	      cout << "Find HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL" << endl;
	      break;
	    }	  

	  if(forOutput.find("DoubleMu")!= std::string::npos && thisTrig.find("HLT_Mu17_Mu8")!= std::string::npos && results==1)
	    {
	      passTrigger=true;
// 	      cout << "Find HLT_Mu17_Mu8" << endl;
	      break;
	    }

	  if(forOutput.find("DoubleMu")!= std::string::npos && thisTrig.find("HLT_Mu17_TkMu8")!= std::string::npos && results==1)
	    {
	      passTrigger=true;
// 	      cout << "Find HLT_Mu17_TkMu8" << endl;
	      break;
	    }

	}

    }
    if(!passTrigger && isData)continue;

    bool pass=false;

    runInfo thisInfo;
    thisInfo.run = runNo;
    thisInfo.evt = eventNo;

    vector<runInfo>::const_iterator location = std::find(myList.begin(), myList.end(), thisInfo);
    if(location != myList.end()) pass=true;

    if(!pass)continue;
    newtree->Fill();
    nPassEvt++;
    fout << runNo << "\t" << eventNo << endl;

    if (jentry%100==0)
      printf("%4.1f%% done.\r",(float)jentry/(float)nentries*100.);		

  }
 
  newtree->Print();
  newtree->AutoSave();
  delete newfile_data;

  fout.close();

  cout << "Number of passed events = " << nPassEvt << endl;
}
