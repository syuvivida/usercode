#define signalEff_pileupSys_cxx
#include "signalEff_pileupSys.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <standalone_LumiReWeighting.cc>
#include <cutvalues.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <fstream>
#include <iomanip>

void signalEff_pileupSys::Loop(int lepCode)
{
  if (fChain == 0) return;
  std::string leptonName = lepCode==0? "electron":"muon";

  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries = " << nentries << endl;
  standalone_LumiReWeighting LumiWeights_central(2012,0);
  standalone_LumiReWeighting LumiWeights_up(2012,+1);
  standalone_LumiReWeighting LumiWeights_down(2012,-1);

  Long64_t nbytes = 0, nb = 0;
  
  double nPassEvts[3]={0};
  double nPassEvts_central[3]={0};
  double nPassEvts_up[3]={0};
  double nPassEvts_down[3]={0};

  Long64_t nPassTotal = 0;

  const int nPUBin = 60;
  TH1D* h_nint_template = new TH1D("h_nint_template","",nPUBin,
 				   -0.5,nPUBin-0.5);
  h_nint_template->SetXTitle("Number of true interactions");
  h_nint_template->Sumw2();


  TH1D* h_input_nint_data   = (TH1D*)h_nint_template->Clone("h_input_nint_data");

  TH1D* h_input_nint_mc   = (TH1D*)h_nint_template->Clone("h_input_nint_mc");


  for(int i=0; i < nPUBin; i++){
    h_input_nint_data->SetBinContent(i+1,Data2012[i]);
    h_input_nint_mc  ->SetBinContent(i+1,Summer2012[i]);
  }

  const int nPUs= 3;

  TH1D* h_output_nint_mc[nPUs];

  std::string titleName[nPUs]={"central", "+5%", "-5%"};

  for(int i=0; i < nPUs; i++)
    {
      h_output_nint_mc[i]
	= (TH1D*)h_nint_template->Clone(Form("h_output_nint_mc%d",i));
      h_output_nint_mc[i]
	-> SetTitle(titleName[i].data());
    }

  ofstream fout;
  fout.open(Form("%s.dat",leptonName.data()));

  ofstream fout0;
  fout0.open(Form("%s_0btag.dat",leptonName.data()));

  ofstream fout1;
  fout1.open(Form("%s_1btag.dat",leptonName.data()));

  ofstream fout2;
  fout2.open(Form("%s_2btag.dat",leptonName.data()));

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    double PU_weight_central =  LumiWeights_central.weight(PU_nTrueInt);
    double PU_weight_up      =  LumiWeights_up.weight(PU_nTrueInt);
    double PU_weight_down    =  LumiWeights_down.weight(PU_nTrueInt);

    int myBest = -1;
    double best_mZjj = 9999999.0;
    double best_mZll = 9999999.0;

    h_output_nint_mc[0]->Fill(PU_nTrueInt, PU_weight_central); 
    h_output_nint_mc[1]->Fill(PU_nTrueInt, PU_weight_up); 
    h_output_nint_mc[2]->Fill(PU_nTrueInt, PU_weight_down); 


    int NBTAGMAX = -1;
    for(int ih=0; ih < lepType->size(); ih++){

      int lepton = lepType->at(ih);
      if(lepton!= lepCode)continue;
    
      int bitmap = passBit->at(ih);

      bool Pass=false;
    
      if((bitmap & MZJJ_SIGNAL) &&
	 (bitmap & PFMET_SIG) &&
	 (bitmap & HELI_LD)
	 )
	Pass=true;
      if(!Pass)continue;
	 
      double zjjMass = zjjM->at(ih);
      double zllMass = zllM->at(ih);

      if(zllMass < MIN_MZ_LL || zllMass > MAX_MZ_LL)continue; 

      int nbtag = nBTags->at(ih);

      if(nbtag > NBTAGMAX)
	{
	  myBest    = ih;
	  best_mZjj = zjjMass;
	  best_mZll = zllMass;
	  NBTAGMAX  = nbtag;
	}
      else if(nbtag == NBTAGMAX)
	{
	  if( ( fabs(zjjMass - MZ_PDG) + fabs(zllMass-MZ_PDG) ) < 
	      ( fabs(best_mZjj - MZ_PDG)+fabs(best_mZll-MZ_PDG))
	      )
 	    {
	      myBest    = ih;
	      best_mZjj = zjjMass;
	      best_mZll = zllMass;
	      NBTAGMAX  = nbtag;	      
	    }
	}      
    } // loop over candidates

    if(myBest<0)continue;
    if(NBTAGMAX<0)continue;

    nPassTotal ++;

    fout << EvtInfo_EventNum << endl;

    if(NBTAGMAX==0)
      fout0 << EvtInfo_EventNum << endl;
    else if(NBTAGMAX==1)
      fout1 << EvtInfo_EventNum << endl;
    else if(NBTAGMAX==2)
      fout2 << EvtInfo_EventNum << endl;

      
    nPassEvts[NBTAGMAX] += 1.0;
    nPassEvts_central[NBTAGMAX]+= PU_weight_central;
    nPassEvts_up[NBTAGMAX]+= PU_weight_up;
    nPassEvts_down[NBTAGMAX]+= PU_weight_down;
    
  } // loop over entries
    
  fout.close();

  fout0.close();
  fout1.close();
  fout2.close();


  double sysP[3]={0};
  double sysM[3]={0};
  for(int ib=0; ib<3; ib++)
    {

      double changeUp = (nPassEvts_up[ib] - nPassEvts_central[ib])/
	nPassEvts_central[ib];

      double changeDown = (nPassEvts_down[ib] - nPassEvts_central[ib])/
	nPassEvts_central[ib];

      cout << "Relative systematic = " << changeUp << "\t" << changeDown << endl;
      if(changeUp>0)
	{
	  sysP[ib] = changeUp*100;
	  sysM[ib] = changeDown*100;
	}
      else
	{
	  sysP[ib] = changeDown*100;
	  sysM[ib] = changeUp*100;
	}
    }

//   std::string remword  ="/home/syu/HZZ/CMSSW_5_2_3_patch2/src/runJob/2012BTag_Sychronized/";

//   size_t pos  = _inputFile.find(remword);

//   if(pos!= std::string::npos)
//     _inputFile.swap(_inputFile.erase(pos,remword.length()));
//   else
//     _inputFile = "test.root";


//   TFile* outFile = new TFile(Form("pusys_%s_%s",leptonName.data(),
// 				  _inputFile.data()),"recreate");   

 
//   h_input_nint_data->Write();
//   h_input_nint_mc  ->Write();

//   for(int i=0; i< nPUs; i++)
//     h_output_nint_mc[i]->Write();

//   outFile->Close();     

  cout << "Total number of events = " << nPassTotal << endl;

  cout << "Separated into different number of btags" << endl;
  for(int ib=0; ib<3; ib++)
    {
      cout << ib << " btag: " << nPassEvts[ib] << endl;
    }
  cout << endl;
  cout << endl;

  for(int ib=0; ib<3; ib++)
    {

      cout << " & ";
      cout << nPassEvts[ib];
    }
  cout << " \\\\" << endl;
  cout << endl;


  for(int ib=0; ib<3; ib++)
    {
      cout << " & ";
      cout << fixed;
      cout << "${+ " << setprecision(2) << sysP[ib] << " \\atop " << 
	sysM[ib] << "}$";     
    }

  cout << " \\\\" << endl;
  cout << endl;
  cout.unsetf(ios_base::fixed);
  cout.precision(6);
    
}

  
