#define pileup_cxx
#include "pileup.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <standalone_LumiReWeighting.cc>
#include <cutvalues.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <fstream>

void pileup::Loop(int lepCode)
{
  if (fChain == 0) return;
  std::string leptonName = lepCode==0? "electron":"muon";
  
  bool isData=false;
  if( _inputFile.find("Run2012")!= std::string::npos) isData=true;

  Long64_t nentries = fChain->GetEntriesFast();
  cout << "isData: " << isData << endl;
  cout << "nentries = " << nentries << endl;
  standalone_LumiReWeighting LumiWeights_central(2012,0);
  standalone_LumiReWeighting LumiWeights_up(2012,+1);
  standalone_LumiReWeighting LumiWeights_down(2012,-1);
  
  // for comparison with data
  standalone_LumiReWeighting LumiWeights_inflat(2012,2);

  Long64_t nbytes = 0, nb = 0;
  
  double nPassEvts[3]={0};
  double nPassEvts_central[3]={0};
  double nPassEvts_up[3]={0};
  double nPassEvts_down[3]={0};

  Long64_t nPassTotal = 0;
  Long64_t nCount[20]={0};

  TH1D* h_weight = new TH1D("h_weight","",1000,-10,10);

  TH1D* h_nvtx_template = new TH1D("h_nvtx_template","",50,0.5,50.5);
  h_nvtx_template->SetXTitle("Number of reconstructed good vertices");
  h_nvtx_template->Sumw2();

  TH1D* h_eleRho_template = new TH1D("h_eleRho_template","",50,0,60);
  h_eleRho_template->SetXTitle("#rho from kt6PFJetsForIso");
  h_eleRho_template->Sumw2();

  TH1D* h_muoRho_template = new TH1D("h_muoRho_template","",50,0,10);
  h_muoRho_template->SetXTitle("#rho from kt6PFJetsCentralNeutral");
  h_muoRho_template->Sumw2();

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


  const int nPUs= 4;

  TH1D* h_output_nint_mc[nPUs];
  TH1D* h_nvtx[nPUs];
  TH1D* h_eleRho[nPUs];
  TH1D* h_muoRho[nPUs];

  std::string titleName[nPUs]={"central", "+5%", "-5%", "73.5 mb"};

  for(int i=0; i < nPUs; i++)
    {
      h_output_nint_mc[i]
	= (TH1D*)h_nint_template->Clone(Form("h_output_nint_mc%d",i));
      h_output_nint_mc[i]
	-> SetTitle(titleName[i].data());

      h_nvtx[i]
	= (TH1D*)h_nvtx_template->Clone(Form("h_nvtx%d",i));
      h_nvtx[i]
	-> SetTitle(titleName[i].data());

      h_eleRho[i]
	= (TH1D*)h_eleRho_template->Clone(Form("h_eleRho%d",i));
      h_eleRho[i]
	-> SetTitle(titleName[i].data());

      h_muoRho[i]
	= (TH1D*)h_muoRho_template->Clone(Form("h_muoRho%d",i));
      h_muoRho[i]
	-> SetTitle(titleName[i].data());
    }

  ofstream fout;
  fout.open(Form("%s.dat",leptonName.data()));

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    nCount[0]++;
    // if it is data, need to do trigger validation
    bool passTrigger = false;

    #ifdef ISDATA
    if(isData){
      
      for(int it=0; it< trigName->size(); it++)
	{
	  std::string thisTrig= trigName->at(it);
	  int results = trigResults->at(it);
//  	  cout << thisTrig << "\t" << results << endl;

	  if(lepCode==0 && thisTrig.find("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL")
	     != std::string::npos && results==1)
	    {
	      passTrigger=true;
	      break;
	    }	  

	  if(lepCode==1 && thisTrig.find("HLT_Mu17_Mu8")!= std::string::npos && results==1)
	    {
	      passTrigger=true;
	      break;
	    }

	  if(lepCode==1 && thisTrig.find("HLT_Mu17_TkMu8")!= std::string::npos && results==1)
	    {
	      passTrigger=true;
	      break;
	    }

	}

    }
    #endif
    if(!passTrigger && isData)continue;
    nCount[1]++;

    double PU_weight_central =  isData? 1.0: LumiWeights_central.weight(PU_nTrueInt);
    double PU_weight_up      =  isData? 1.0: LumiWeights_up.weight(PU_nTrueInt);
    double PU_weight_down    =  isData? 1.0: LumiWeights_down.weight(PU_nTrueInt);
    double PU_weight_inflat  =  isData? 1.0: LumiWeights_inflat.weight(PU_nTrueInt);

    h_weight->Fill(PU_weight_central);
    double weight[4] = {PU_weight_central,
			PU_weight_up,
			PU_weight_down,
			PU_weight_inflat};

    for(int ip=0; ip < nPUs; ip++)    
      h_output_nint_mc[ip]->Fill(PU_nTrueInt, weight[ip]); 


    int NBTAGMAX = -1;
    int myBest = -1;
    double best_mZjj = 9999999.0;

    for(int ih=0; ih < lepType->size(); ih++){

      int lepton = lepType->at(ih);
      if(lepton!= lepCode)continue;
    
      int bitmap = passBit->at(ih);

      bool Pass=false;
    
//       if(
// 	 (bitmap & MZJJ_SIGNAL) 
//       //  	 && (bitmap & PFMET_SIG) 
//       //  	 && (bitmap & HELI_LD)
//  	 )
	Pass=true;
      if(!Pass)continue;
	 
      double zjjMass = zjjM->at(ih);

      int nbtag = nBTags->at(ih);

      if(nbtag > NBTAGMAX)
	{
	  myBest    = ih;
	  best_mZjj = zjjMass;
	  NBTAGMAX  = nbtag;
	}
      else if(nbtag == NBTAGMAX)
	{
	  if(fabs(zjjMass - MZ_PDG) < fabs(best_mZjj - MZ_PDG))
 	    {
	      myBest    = ih;
	      best_mZjj = zjjMass;
	      NBTAGMAX  = nbtag;	      
	    }
	}      
    } // loop over candidates

    nCount[2]++;
    if(myBest<0)continue;
    if(NBTAGMAX<0)continue;

    nCount[3]++;

    nPassTotal ++;

    fout << EvtInfo_EventNum << endl;
    
    for(int ip=0; ip < nPUs; ip++){

      h_nvtx[ip]->Fill( EvtInfo_NumVtx, weight[ip]);
      h_eleRho[ip]->Fill( eleRho, weight[ip]);
      h_muoRho[ip]->Fill( muoRho, weight[ip]);

    }     

    nPassEvts[NBTAGMAX] += 1.0;
    nPassEvts_central[NBTAGMAX]+= weight[0];
    nPassEvts_up[NBTAGMAX]+= weight[1];
    nPassEvts_down[NBTAGMAX]+= weight[2];
    
  } // loop over entries
    
  fout.close();


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

  std::string remword  ="/data4/syu/hzzTrees/";

  size_t pos  = _inputFile.find(remword);

  if(pos!= std::string::npos)
    _inputFile.swap(_inputFile.erase(pos,remword.length()));
  else
    _inputFile = "test.root";


  TFile* outFile = new TFile(Form("loosemjj_pileup_%s_%s",leptonName.data(),
				  _inputFile.data()),"recreate");   

  h_input_nint_data->Write();
  h_input_nint_mc  ->Write();

  h_weight->Write();

  for(int i=0; i< nPUs; i++)
    {
      h_output_nint_mc[i]->Write();
      h_nvtx[i]->Write();
      h_eleRho[i]->Write();
      h_muoRho[i]->Write();
    }


  outFile->Close();     

  cout << "Total number of events = " << nPassTotal << endl;

  cout << "Separated into different number of btags" << endl;
  for(int ib=0; ib<3; ib++)
    {
      cout << ib << " btag: " << nPassEvts[ib] << endl;
    }
  cout << endl;
  cout << endl;
  for(int icount=0; icount<20;icount++)
    if(nCount[icount]>0)
      cout << "nCount[" << icount << "]=" << nCount[icount] << endl;

}

  
