#define signalEff_pileupSys_cxx
#include "signalEff_pileupSys.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <standalone_LumiReWeighting.cc>
#include <cutvalues.h>
#include <TLorentzVector.h>
#include <TMath.h>

void signalEff_pileupSys::Loop(int lepCode)
{
  if (fChain == 0) return;

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

  TH1D* hsys[3];
  TH1D* hpass[3];
  for(int i=0; i<3; i++)
    {
      hsys[i] = new TH1D(Form("hsys%d",i),Form("nbtag=%d",i),4,0.5,4.5);
      hpass[i] = new TH1D(Form("hpass%d",i),Form("nbtag=%d",i),3,-1.5,1.5);
    }

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

    h_output_nint_mc[0]->Fill(PU_nTrueInt, PU_weight_central); 
    h_output_nint_mc[1]->Fill(PU_nTrueInt, PU_weight_up); 
    h_output_nint_mc[2]->Fill(PU_nTrueInt, PU_weight_down); 

//     cout << "Looking for best candidate" << endl;

    for(int ih=0; ih < lepType->size(); ih++){

      int lepton = lepType->at(ih);
      if(lepton!= lepCode)continue;
    
      int bitmap = passBit->at(ih);

      bool Pass=false;
      if(bitmap & ALL_SIGNAL)Pass=true;
      if(!Pass)continue;

      double zjjMass = zjjM->at(ih);
      if(fabs(zjjMass - MZ_PDG) < fabs(best_mZjj - MZ_PDG))
	{
	  best_mZjj = zjjMass;
	  myBest    = ih;
// 	  cout << "myBest = " << myBest << "\t zjjM = " << best_mZjj << 
// 	    "\t bestHCand = " << bestHCand << endl; 
	}
      
    }

    if(myBest<0)continue;

    int nbtag=nBTags->at(myBest);
      
    if(nbtag>=0){
      nPassEvts[nbtag] += 1.0;
      nPassEvts_central[nbtag]+= PU_weight_central;
      nPassEvts_up[nbtag]+= PU_weight_up;
      nPassEvts_down[nbtag]+= PU_weight_down;

    }

  } // loop over entries
    
  for(int ib=0; ib<3; ib++)
    {
      cout << "npass raw with " << ib << " btag = " << nPassEvts[ib] << endl;
      cout << "npass with " << ib << " btag = " << nPassEvts_central[ib] << 
	"\t" << nPassEvts_up[ib] << "\t" << nPassEvts_down[ib] << endl;

      hpass[ib]->SetBinContent(1, nPassEvts_down[ib]);
      hpass[ib]->SetBinContent(2, nPassEvts_central[ib]);
      hpass[ib]->SetBinContent(3, nPassEvts_up[ib]);

      double changeUp = (nPassEvts_up[ib] - nPassEvts_central[ib])/
	nPassEvts_central[ib];

      hsys[ib]->SetBinContent(1, changeUp);

      double changeDown = (nPassEvts_down[ib] - nPassEvts_central[ib])/
	nPassEvts_central[ib];

      hsys[ib]->SetBinContent(2, changeDown);

      cout << "Relative systematic = " << changeUp << "\t" << changeDown << endl;
//       cout << "7 TeV style Relative systematic = " << changeUpInt << "\t" << 
// 	changeDownInt << endl;

	
    }

  std::string remword  ="/home/syu/HZZ/CMSSW_5_2_3_patch2/src/runJob/match_dRrelPt/";

  size_t pos  = _inputFile.find(remword);

  if(pos!= std::string::npos)
    _inputFile.swap(_inputFile.erase(pos,remword.length()));
  else
    _inputFile = "test.root";

  std::string leptonName = lepCode==0? "electron":"muon";
  TFile* outFile = new TFile(Form("pusys_%s_%s",leptonName.data(),
				  _inputFile.data()),"recreate");   

  for(int ib=0; ib<3; ib++)
    {
      hsys[ib]->Write();
      hpass[ib]->Write();
    }
 
  h_input_nint_data->Write();
  h_input_nint_mc  ->Write();

  for(int i=0; i< nPUs; i++)
    h_output_nint_mc[i]->Write();

  outFile->Close();     



}

  
