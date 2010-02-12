#include "tightSelection.h"
#include <map>
#include <string>
#include <TChain.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <iostream>
#include <algorithm>
#include <fstream>

using namespace std;

typedef map<int,int> RunMap;

void dumpRunNumber(std::string dataname)
{

  TChain* t1 = new TChain("NTuples/Analysis");
  t1->Add(dataname.data());

  RunMap runs;

  // loop over event entries
  Long64_t nentries = (Long64_t)t1->GetEntries();
  for( Long64_t j=0; j< nentries; j++){

    // this line is very important!!
    t1->GetEntry(j);
    int thisrun = (int)t1->GetBranch("run")->GetLeaf("run")->GetValue();
    
    int vtxIsFake = (int)t1->GetBranch("vtxIsFake")->GetLeaf("vtxIsFake")->GetValue();

    int nHfTowersP = (int)t1->GetBranch("nHfTowersP")->GetLeaf("nHfTowersP")->GetValue();
    int nHfTowersN = (int)t1->GetBranch("nHfTowersN")->GetLeaf("nHfTowersN")->GetValue();
    if(vtxIsFake==1)continue;
    if(nHfTowersP < 1)continue;
    if(nHfTowersN < 1)continue;
    
    if(runs.find(thisrun)!=runs.end())
      runs[thisrun]++;
    else
      runs.insert(std::pair<int,int>(thisrun,1));

  }
  
  ofstream fout;
  fout.open("runs.txt");

  int ngood = 0;
  for(RunMap::iterator i=runs.begin(); i!=runs.end(); ++i)
    {
      if(i->second < 400)continue;
      cout << i->first << " " << i->second << endl;
      fout << i->first << endl;
      ngood ++;
    }
  cout << "There are " << ngood << " good runs" << endl;
  fout.close();

}
		     
