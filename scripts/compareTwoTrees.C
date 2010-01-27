#include <TCut.h>
#include "selection.h"

void compareTwoTrees(TChain* t1, TChain* t2, 
		     std::string var1,
		     int decCode, 
		     std::string psname,
		     int nbin, float xmin, float xmax, 
		     std::string xtitle="", std::string ytitle="",
		     TCut cut1="")
{
  // for cuts the same for barrel and endcap
  TCut allCut = cut1 + myCut;

  if(decCode==1)
    allCut += barrelCut;
  else if(decCode==2)
    allCut += endcapCut;

//   // only for data
  TCut dataCut;
  if(psname.find("2360") != std::string::npos)
    dataCut= Cut_2360GeV;
  else
    dataCut= Cut_900GeV;

  
  gStyle->SetOptStat(0);


  float binwidth = (xmax-xmin)/(float)nbin;

  char tempName[300];
  sprintf(tempName, "Candidates per %1.3lf", binwidth);
  if(binwidth < 0.001)
  sprintf(tempName, "Candidates per %1.4lf", binwidth);
  if(binwidth < 0.0001)
  sprintf(tempName, "Candidates per %1.5lf", binwidth);
  if(ytitle == "")ytitle = tempName;

  
  TH1F* h1 = new TH1F("h1","",nbin,xmin,xmax);
  h1->SetXTitle(xtitle.data());  
  h1->SetYTitle(ytitle.data());  

  TH1F* h2 = new TH1F("h2","",nbin,xmin,xmax);
  h2->SetXTitle(xtitle.data());
  h2->SetYTitle(ytitle.data());  

  TH1F* hscale = new TH1F("hscale","",nbin,xmin,xmax);
  hscale->SetXTitle(xtitle.data());
  hscale->SetYTitle("Data/MC");

  std::string temp = var1 + ">>h1";

  t1->Draw(temp.data(), allCut && dataCut);
  h1->SetLineColor(4);
  h1->Draw();
  cout << "h1 mean = " << h1->GetMean() << " and width = " << h1->GetRMS() << 
    " and entries = " << h1->GetEntries() << endl;

  temp = var1 + ">>h2";
  t2->Draw(temp.data(), allCut);
  h2->SetLineColor(2);
  h2->Draw();
  cout << "h2 mean = " << h2->GetMean() << " and width = " << h2->GetRMS() << 
    " and entries = " << h2->GetEntries() << endl;

  
  gROOT->ProcessLine(".L /home/syu/testYenJie/CMSSW_3_3_6_patch3/src/CRAB/scripts/computeChi2New.C");
  computeChi2New(h1,h2,hscale,psname,decCode,false);
  computeChi2New(h1,h2,hscale,psname,decCode,true);

  delete h1;
  delete h2;
  delete hscale;


}
		     
