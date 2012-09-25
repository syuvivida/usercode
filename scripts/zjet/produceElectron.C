#include <TH1.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TFile.h>
#include <iostream>
#include <string>

using namespace std;

void produceElectron(std::string eikoName="h_jety", 
		     bool update=false,
		     bool acceptanceCorr=true,
		     double correlation=1.0, 
		     bool logScale=false
		     )
{
  std::string xtitle;
  if(eikoName=="h_ystar")
    xtitle = "0.5|Y_{Z}-Y_{jet}|";
  else if(eikoName=="h_yB")
    xtitle = "0.5|Y_{Z}+Y_{jet}|";
  else if(eikoName=="h_jety")
    xtitle = "|Y(jet)|";
  else if(eikoName=="h_zy")
    xtitle = "|Y(Z)|";
  
  // declare histograms
  TH1D* h_e;

  TH1D* h_raw;

  std::string remword3  ="h_";
  std::string corrName = eikoName;
  size_t pos3  = corrName.find(remword3);
  if(pos3!= std::string::npos)
    corrName.replace(pos3,remword3.length(),"");

  // acceptance correction for electrons
  TFile f_crack("mainCore/ave_sherpamadgraph.root");
  if (f_crack.IsZombie()) {
    cout << endl << "Error opening file" << f_crack.GetName() << endl;
    return;
  }
  else
    cout << endl << "Opened " << f_crack.GetName() << endl << endl;

  h_corr = (TH1D*)(f_crack.Get(Form("have_%s",corrName.data())));

  // for debugging
  if( !acceptanceCorr )
    h_corr->Reset();
  for(int i=1; i<= h_corr->GetNbinsX(); i++)
    {
      if( !acceptanceCorr )
	{
	  h_corr->SetBinContent(i,1.0);
	  h_corr->SetBinError(i,1e-6);
	}
      cout << "Correction for bin " << i << " = " 
	   << h_corr->GetBinContent(i) << " +- " << h_corr->GetBinError(i) 
	   << endl;
    }
  
  // electron channel

  TFile f_e("mainCore/cts_CorrectedPlotsZCut.root");
  if (f_e.IsZombie()) {
    cout << endl << "Error opening file" << f_e.GetName() << endl << endl;
    return;
  }
  else
    cout << endl << "Opened " << f_e.GetName() << endl << endl;
  h_raw = (TH1D*)(f_e.Get(eikoName.data()));
  h_raw  -> SetName(Form("h_raw_%s",corrName.data()));


  h_e  = (TH1D*)(f_e.Get(eikoName.data()));
  h_e  -> SetName(eikoName.data());
  h_e  -> SetTitle(""); 


  double darko[12][2]={
    {13218.9 , 155.4},
    {13038.2 , 159.7},
    {12338.0 , 155.2},
    {12004.6 , 153.5},
    {11173.1 , 150.4},
    { 9959.5 , 142.6},
    { 9266.5 , 137.1},
    { 7725.4 , 126.6},
    { 6897.7 , 120.1},
    { 5876.9 , 110.5},
    { 4814.4 , 101.1},
    { 3277.8 , 80.8}
  };

  if(eikoName=="h_jety")
    {
      h_raw  -> Reset();
      h_e    -> Reset();

      for(int i=1;i<=12;i++){ 

	h_raw->SetBinContent(i, darko[i-1][0]);
	h_raw->SetBinError(i, darko[i-1][1]);

	h_e->SetBinContent(i, darko[i-1][0]);
	h_e->SetBinError(i, darko[i-1][1]);
	
      } 
    }

  //===================================================
  // 2012/09/10, New!! crack acceptance correction
  //===================================================
  h_e  -> Sumw2();
  h_e  -> Divide(h_corr);
  h_e  -> SetYTitle("Arbitrary Unit");
  h_e  -> SetTitleOffset(2.0,"Y");
  h_e  -> GetYaxis()->SetDecimals();
  h_e  -> GetXaxis()->SetDecimals();
  h_e  -> SetLineColor(kBlue-7);
  h_e  -> SetMarkerColor(kBlue-7);
  h_e  -> SetMarkerSize(1);
  h_e  -> SetMarkerStyle(24);

  cout << "h_e integral = " << h_e->Integral() << endl;

  // to get the JES of electron channel

  h_e      ->SetXTitle(xtitle.data());


  // save the original electron and muon root files and the combined 
  // result in a ROOT file

  std::string command = "recreate";
  if(update)command="update";
  TFile* outFile = new TFile("corrected_electron.root", command.data());
  h_e      ->Write();
  h_raw    ->Write();

}
