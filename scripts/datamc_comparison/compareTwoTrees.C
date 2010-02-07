#include <TCut.h>
#include "tightSelection.h"
#include "computeChi2New.C"

void compareTwoTrees(TChain* t1, TChain* t2, 
		     std::string var1,
		     int decCode, 
		     std::string psname,
		     int nbin, float xmin, float xmax, 
		     std::string xtitle="", std::string ytitle="",
		     TCut cut1="",float yMAX=3.8)
{
  // for cuts the same for barrel and endcap
  TCut allCut = cut1 + myCut;

  // cuts for getting overall scale factor
  TCut scaleCut = allCut;

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

  

  float binwidth = (xmax-xmin)/(float)nbin;

  char tempName[300];
  sprintf(tempName, "Candidates per %1.3lf", binwidth);
  if(binwidth < 0.001)
  sprintf(tempName, "Candidates per %1.4lf", binwidth);
  if(binwidth < 0.0001)
  sprintf(tempName, "Candidates per %1.5lf", binwidth);
  if(ytitle == "")ytitle = tempName;
  
  float runmin = 123586;
  float runmax = 124240;
  TProfile* hf1 = new TProfile("hf1","",16,runmin,runmax);
  std::string temp_variable = var1 + ": run >>hf1";
  t1->Draw(temp_variable.data(), allCut && dataCut,"prof");

  TCanvas* c2 = new TCanvas("c2");
  hf1->SetTitle("");
  hf1->SetXTitle("Run");
  hf1->SetMarkerColor(2);
  hf1->SetMarkerSize(1);
  hf1->SetYTitle(xtitle.data());
  hf1->Draw("");

  std::string decName;
  if(decCode==0)decName = "_all";
  else if(decCode==1)decName = "_barrel";
  else if(decCode==2)decName = "_endcap";
  
  std::string runName = "$CMSSW_BASE/src/CRAB/preprod_figures/runDep_" + psname + decName + ".eps";
  c2->Print(runName.data());
  runName = "$CMSSW_BASE/src/CRAB/preprod_figures/runDep_" + psname + decName + ".gif";
  c2->Print(runName.data());
  
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


  // below is to get the overall scale factor

  TH1F* scaleh1 = new TH1F("scaleh1","",nbin,xmin,xmax);
  TH1F* scaleh2 = new TH1F("scaleh2","",nbin,xmin,xmax);
  temp = var1 + ">>scaleh1";
  t1->Draw(temp.data(), scaleCut && dataCut);
  temp = var1 + ">>scaleh2";
  t2->Draw(temp.data(), scaleCut);
  float scaleFactor = (float)scaleh1->GetEntries()/(float)scaleh2->GetEntries();
  hscale->SetBinContent(1, scaleFactor);
  hscale->SetBinError(1,0.001);
  //  gROOT->ProcessLine(".L $PWD/computeChi2New.C");
  // in linear scale
  computeChi2New(h1,h2,hscale,psname,decCode,false,yMAX);
  // in log scale
  hscale->Reset();
  hscale->SetBinContent(1, scaleFactor);
  hscale->SetBinError(1,0.001);
  computeChi2New(h1,h2,hscale,psname,decCode,true,yMAX);



  delete h1;
  delete h2;
  delete hscale;
  delete scaleh1;
  delete scaleh2;
  delete hf1;
  delete c2;

}
		     
