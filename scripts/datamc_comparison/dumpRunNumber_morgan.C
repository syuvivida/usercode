#include "tightSelection.h"
#include <map>
#include <string>
#include <TChain.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <iostream>
#include <fstream>

using namespace std;


typedef map<int,int> RunMap;

void dumpRunNumber(std::string dataname)
{

  TChain* t1 = new TChain("NTuples/Analysis");
  t1->Add(dataname.data());

  RunMap runs;
  RunMap minruns;
  RunMap maxruns;

  // loop over event entries
  Long64_t nentries = (Long64_t)t1->GetEntries();
  for( Long64_t j=0; j< nentries; j++){

    // this line is very important!!
    t1->GetEntry(j);
    int thisrun = (int)t1->GetBranch("run")->GetLeaf("run")->GetValue();
    int lumi    = (int)t1->GetBranch("luminosityBlock")->GetLeaf("luminosityBlock")->GetValue();
    
    int vtxIsFake = (int)t1->GetBranch("vtxIsFake")->GetLeaf("vtxIsFake")->GetValue();

    int nHfTowersP = (int)t1->GetBranch("nHfTowersP")->GetLeaf("nHfTowersP")->GetValue();
    int nHfTowersN = (int)t1->GetBranch("nHfTowersN")->GetLeaf("nHfTowersN")->GetValue();
//     if(vtxIsFake==1)continue;
//     if(nHfTowersP < 1)continue;
//     if(nHfTowersN < 1)continue;
    
//     if ( ! (thisrun==123591 || thisrun==123596 || thisrun==123615 || thisrun==123732
// 	    || thisrun==123818 || thisrun==123906 || thisrun==123908 || thisrun==123970
// 	    || thisrun==123976 || thisrun==123977 || thisrun==123978 || thisrun==123987
// 	    || thisrun==124006 || thisrun==124008 || thisrun==124009 || thisrun==124020
// 	    || thisrun==124022 || thisrun==124023 || thisrun==124024 || thisrun==124025
// 	    || thisrun==124027 || thisrun==124030
// 	    ) ) continue;

    if(!(thisrun == 123596 ||
	 thisrun == 123615 ||
	 thisrun == 123732 ||
	 thisrun == 123818 ||
	 thisrun == 123906 ||
	 thisrun == 123908 ||
	 thisrun == 123909 ||
	 thisrun == 124009 ||
	 thisrun == 124020 ||
	 thisrun == 124022 ||
	 thisrun == 124023 ||
	 thisrun == 124024 ||
	 thisrun == 124025 ||
	 thisrun == 124027 ||
	 thisrun == 124030))continue;


   if(runs.find(thisrun)!=runs.end())
      runs[thisrun]++;
    else
      runs.insert(std::pair<int,int>(thisrun,1));

    if(minruns.find(thisrun)!= 
       minruns.end())
      {
	if(lumi < minruns[thisrun])minruns[thisrun]=lumi;
      }
    else
      minruns.insert(std::pair<int,int>(thisrun,lumi));


    if(maxruns.find(thisrun)!= 
       maxruns.end())
      {
	if(lumi > maxruns[thisrun])maxruns[thisrun]=lumi;
      }
    else
      maxruns.insert(std::pair<int,int>(thisrun,lumi));

  }
  
//   ofstream fout;
//   fout.open("runs.txt");

//   int ngood = 0;
//   for(RunMap::iterator i=runs.begin(); i!=runs.end(); ++i)
//     {
//       if(i->second < 0)continue;
//       cout << i->first << " " << i->second << endl;
//       fout << i->first << endl;
//       ngood ++;
//     }
//   cout << "There are " << ngood << " good runs" << endl;
//   fout.close();

  for(RunMap::iterator i=minruns.begin(); i!=minruns.end(); ++i)
    {
      cout << i->first << " " << i->second << " -- " << 
	maxruns[i->first] << endl;
    }

}
		     
