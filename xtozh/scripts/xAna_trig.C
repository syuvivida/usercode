//root -q -b -l juwu.C++\(\"inputFile\"\,\"outputFile\"\)


#include <vector>
#include <iostream>
#include <algorithm>
#include <map>
#include <TH1D.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TFile.h>
#include "untuplizer.h"
#include <TGraphAsymmErrors.h>
// #include "TrigEff_AsymmetryErrors.C"

using namespace std;
void xAna_trig(std::string inputFile){

  //get TTree from file ...
  TreeReader data(inputFile.data());

  //histogram anoucement
  float xmin[12]={0, 0.05,0.1,
		  0.15,0.2,
		  0.25,0.3,
		  0.35,0.4,
		  0.45,0.5,
		  1.5};
  TH1F* h_mZ   = new TH1F("h_mZ","",100,50,150);
  TH1F* h_dR   = new TH1F("h_dR","",11,xmin);
  h_dR->SetXTitle("#Delta R between generator-level muons");
  TH1F* h_dR_deno = (TH1F*)h_dR->Clone("h_dR_deno");
  TH1F* h_dR_numr_muonReco[3];
  TH1F* h_dR_numr_muonTrig[3];

  std::string title[3]={"Two Global muon", "1 Global + 1 Trk", "2 Trks"};
  std::string muTrigNames[3]={"HLT_Mu17_Mu8",
			      "HLT_Mu17_TkMu8",
			      "HLT_Mu30_TkMu11"};

  for(int i=0; i < 3; i++)
    {
      h_dR_numr_muonReco[i] = (TH1F*)h_dR->Clone(Form("h_dR_numr_muonReco%d",i));
      h_dR_numr_muonReco[i] -> SetTitle(title[i].data());
      h_dR_numr_muonTrig[i] = (TH1F*)h_dR->Clone(Form("h_dR_numr_muonTrig%d",i));
      h_dR_numr_muonTrig[i] -> SetTitle(muTrigNames[i].data());
    }
  
  //Event loop
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){

    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());

    data.GetEntry(jEntry);


    std::string* trigName = data.GetPtrString("trigName");
    Int_t* trigResult = data.GetPtrInt("hlt_trigResult");

    Int_t nGenPar        = data.GetInt("nGenPar");
    Int_t* genParId      = data.GetPtrInt("genParId");
    Int_t* genParSt      = data.GetPtrInt("genParSt");
    Int_t* genMomParId   = data.GetPtrInt("genMomParId");
    Float_t* genParPt    = data.GetPtrFloat("genParPt");
    Float_t* genParEta    = data.GetPtrFloat("genParEta");
    Float_t* genParPhi    = data.GetPtrFloat("genParPhi");
    Float_t* genParM    = data.GetPtrFloat("genParM");

    int muIndex[2]={-1,-1};
    for(int ig=0; ig < nGenPar; ig++){

      if(abs(genParId[ig])!=13)continue;
      if(genParSt[ig]!=1)continue;
      if(genMomParId[ig] != 23 && genMomParId[ig] != genParId[ig])continue;

      if(genParPt[ig]<20)continue;
      if(abs(genParEta[ig])>2.4)continue;

      if(muIndex[0]<0)muIndex[0]=ig;
      else if(muIndex[1]<0)muIndex[1]=ig;
      if(muIndex[0]>=0 && muIndex[1]>=0)break;
    }


    if(muIndex[0]<0 || muIndex[1]<0)continue;

    if(genParPt[muIndex[0]] < 40 && genParPt[muIndex[1]] < 40)continue;

    TLorentzVector l4_mu[2];
    for(int i=0;i<2;i++)
      {
	l4_mu[i].SetPtEtaPhiM(
			      genParPt[muIndex[i]],
			      genParEta[muIndex[i]],
			      genParPhi[muIndex[i]],
			      genParM[muIndex[i]]
			      );

      }
    float gendR = l4_mu[0].DeltaR(l4_mu[1]);
    h_mZ->Fill((l4_mu[0]+l4_mu[1]).M());
    h_dR_deno->Fill(gendR);



    Float_t* muPt        = data.GetPtrFloat("muPt");
    Float_t* muEta       = data.GetPtrFloat("muEta");
    Float_t* muPhi       = data.GetPtrFloat("muPhi");
    Float_t* muM         = data.GetPtrFloat("muM");
    Int_t*   isGlobalMuon = data.GetPtrInt("isGlobalMuon");
    Int_t*   isTrackerMuon = data.GetPtrInt("isTrackerMuon");
    Int_t    nMu   = data.GetInt("nMu"); 

    bool Matched[2][2];

    for(int imu=0; imu <2; imu++)
      {
	float dRMin[2];
	for(int icut=0;icut<2; icut++){
	  Matched[imu][icut]=false;
	  dRMin[icut]=0.4;
	}
	
	for(int i=0; i < nMu; i++){
	  
	  bool cut[2]={
	    (isGlobalMuon[i]==1),
	    (isTrackerMuon[i]==1)};
	    
	  
	  TLorentzVector thisMu(0,0,0,0);
	  thisMu.SetPtEtaPhiM
	    (
	     muPt[i],
	     muEta[i],
	     muPhi[i],
	     muM[i]
	     );
	  float thisDR= thisMu.DeltaR(l4_mu[imu]);
	  for(int icut=0; icut<2;icut++){
	    if(thisDR<dRMin[icut] && cut[icut])
	    {
	      dRMin[icut] = thisDR;
	      Matched[imu][icut]=true;

	    }
	  }// end of loop over cuts

	}

      }



    if(Matched[0][0] && Matched[1][0])
      h_dR_numr_muonReco[0]->Fill(gendR);
    if(Matched[0][1] && Matched[1][1])
      h_dR_numr_muonReco[2]->Fill(gendR); 
    if( (Matched[0][0] && Matched[1][0]) ||
	(Matched[0][0] && Matched[1][1]) ||
	(Matched[0][1] && Matched[1][0]) 
	)
      h_dR_numr_muonReco[1]->Fill(gendR);
     
    if( ( Matched[0][0] || Matched[0][1]) &&
	( Matched[1][0] || Matched[1][1])) {

	  const Int_t nsize = data.GetPtrStringSize();
	  bool passTrigger[3]={false,false,false};
	  for(int it=0; it< nsize; it++)
	    {
	      std::string thisTrig= trigName[it];
	      int results = trigResult[it];
	  
	      for(int i=0;i<3;i++){
		if(thisTrig.find(muTrigNames[i].data())!= std::string::npos 
		   && results==1)
		  passTrigger[i]=true;
	      }
	    }

	  for(int i=0;i<3;i++)
	    {
	      if(passTrigger[i])
		h_dR_numr_muonTrig[i]->Fill(gendR);
	
	    }
	}  

  } // end of loop over entries


  TGraphAsymmErrors* htrigeff[3];
  TGraphAsymmErrors* hrecoeff[3];

  //save output
  TString endfix=gSystem->GetFromPipe(Form("file=%s; echo \"${file##*/}\"",inputFile.data()));
  TString outputFile = "trigHisto_" + endfix;
  TFile* outFile = new TFile(outputFile.Data(),"recreate");
  h_dR_deno ->Write();
  h_mZ->Write();
  for(int i=0; i<3; i++){
  h_dR_numr_muonReco[i]->Write();
  h_dR_numr_muonTrig[i]->Write();
  htrigeff[i]=new TGraphAsymmErrors(h_dR_numr_muonTrig[i], h_dR_deno);
				    // "cl=0.683 b(1,1) mode");
  htrigeff[i]->SetName(Form("htrigeff%02i",i));
  htrigeff[i]->GetXaxis()->SetTitle("#Delta R between generator-level muons");
  htrigeff[i]->GetYaxis()->SetTitle("Efficiency");

  hrecoeff[i]=new TGraphAsymmErrors(h_dR_numr_muonReco[i], h_dR_deno);
  hrecoeff[i]->SetName(Form("hrecoeff%02i",i));
  hrecoeff[i]->GetXaxis()->SetTitle("#Delta R between generator-level muons");
  hrecoeff[i]->GetYaxis()->SetTitle("Reconstruction Efficiency");


  htrigeff[i]->Write();
  hrecoeff[i]->Write();
  }
  outFile->Close();
}
