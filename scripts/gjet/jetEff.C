#define jetEff_cxx
#include "jetEff.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TProfile.h>

const Float_t ENDCAP_MINETA=2.6;
const Float_t ENDCAP_MAXETA=3.0;
const Int_t   NETA=2;
const Float_t fEtaBin[]={0,ENDCAP_MINETA,ENDCAP_MAXETA};


void jetEff::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   // declare and define histograms
   TH1F* h_dec = new TH1F("h_dec","",2,fEtaBin);

   TH1F* h_nvtx_template = new TH1F("h_nvtx_template","",30,-0.5,29.5);
   TH1F* h_ngood  = (TH1F*)h_nvtx_template->Clone("h_ngood");


   TH1F* h_eta_template = new TH1F("h_eta_template","",60,-3.0,3.0);
   h_eta_template->SetXTitle("#eta(jet)");
   TH1F* h_recetaall = (TH1F*)h_eta_template->Clone("h_recetaall");
   TH1F* h_recetapass = (TH1F*)h_eta_template->Clone("h_recetapass");
   TH1F* h_recetapass_medium = (TH1F*)h_eta_template->Clone("h_recetapass_medium");
   TH1F* h_recetapass_tight = (TH1F*)h_eta_template->Clone("h_recetapass_tight");
   TGraphAsymmErrors* h_ideff_eta = new TGraphAsymmErrors();
   h_ideff_eta->SetNameTitle("ideff_eta", 
			       "Jet loose ID efficiency as a function of eta");
   h_ideff_eta->GetXaxis()->SetTitle("#eta(jet)");
   h_ideff_eta->GetYaxis()->SetTitle("Loose jet ID efficiency");

   TGraphAsymmErrors* h_ideff_eta_medium = new TGraphAsymmErrors();
   h_ideff_eta_medium->SetNameTitle("ideff_eta_medium", 
			       "Jet medium ID efficiency as a function of eta");
   h_ideff_eta_medium->GetXaxis()->SetTitle("#eta(jet)");
   h_ideff_eta_medium->GetYaxis()->SetTitle("Medium jet ID efficiency");

   TGraphAsymmErrors* h_ideff_eta_tight = new TGraphAsymmErrors();
   h_ideff_eta_tight->SetNameTitle("ideff_eta_tight", 
			       "Jet tight ID efficiency as a function of eta");
   h_ideff_eta_tight->GetXaxis()->SetTitle("#eta(jet)");
   h_ideff_eta_tight->GetYaxis()->SetTitle("Tight jet ID efficiency");



   TH1F* h_pt_template = new TH1F("h_pt_template","",50,0,200);
   h_pt_template->SetXTitle("p_{T}(jet) [GeV]");

   TH1F* h_genpt = (TH1F*)h_pt_template->Clone("h_genpt");

   //0: barrel, 1: endcap, for ID efficiency
   TH1F* h_recptall[NETA];
   TH1F* h_recptpass[NETA];
   TH1F* h_recptpass_medium[NETA];
   TH1F* h_recptpass_tight[NETA];
   TGraphAsymmErrors* h_ideff_pt[NETA];
   TGraphAsymmErrors* h_ideff_pt_medium[NETA];
   TGraphAsymmErrors* h_ideff_pt_tight[NETA];

   TH1F* h_ngood_all[NETA];
   TH1F* h_ngood_pass[NETA];
   TGraphAsymmErrors* h_ideff_ngood[NETA];


   // profile template
   TProfile* pf_pt_template = new TProfile("pf_pt_template","",50,0,200,-10.0,10.0);
   pf_pt_template->SetXTitle("p_{T}(jet) [GeV]");
   pf_pt_template->SetYTitle("p_{T}^{RECO}(jet)/p_{T}^{GEN}(jet)"); 

   TProfile* pf_alljet_ratio = (TProfile*)pf_pt_template->Clone("pf_alljet_ratio");
   TProfile* pf_1stjet_ratio = (TProfile*)pf_pt_template->Clone("pf_1stjet_ratio");
   TProfile* pf_2ndjet_ratio = (TProfile*)pf_pt_template->Clone("pf_2ndjet_ratio");

   TProfile* pf_eta_template = new TProfile("pf_eta_template","",60,-3.0,3.0,-10.0,10.0);
   pf_eta_template->SetXTitle("#eta(jet)");
   pf_eta_template->SetYTitle("p_{T}^{GEN}(jet)/p_{T}^{RAW}(jet)"); 

   TProfile* pf_alljet_corr_eta = (TProfile*)pf_eta_template->Clone("pf_alljet_corr_eta");
   TProfile* pf_1stjet_corr_eta = (TProfile*)pf_eta_template->Clone("pf_1stjet_corr_eta");
   TProfile* pf_2ndjet_corr_eta = (TProfile*)pf_eta_template->Clone("pf_2ndjet_corr_eta");
					   

   for(int ieta=0; ieta< NETA; ieta++){
   
     h_recptall[ieta]  = (TH1F*)h_pt_template->Clone(Form("h_recptall%d",ieta));
     h_recptpass[ieta] = (TH1F*)h_pt_template->Clone(Form("h_recptpass%d",ieta));
     h_recptpass_medium[ieta] = (TH1F*)h_pt_template->Clone(Form("h_recptpass_medium%d",ieta));
     h_recptpass_tight[ieta] = (TH1F*)h_pt_template->Clone(Form("h_recptpass_tight%d",ieta));
   
     h_ideff_pt[ieta]  = new TGraphAsymmErrors();
     h_ideff_pt[ieta]->SetNameTitle(Form("ideff_pt%d",ieta),
			       "Jet loose ID efficiency as a function of pT");
     h_ideff_pt[ieta]->GetXaxis()->SetTitle("p_{T}(jet) [GeV]");
     h_ideff_pt[ieta]->GetYaxis()->SetTitle("Loose jet ID efficiency");


     h_ideff_pt_medium[ieta]  = new TGraphAsymmErrors();
     h_ideff_pt_medium[ieta]->SetNameTitle(Form("ideff_pt_medium%d",ieta),
			       "Jet medium ID efficiency as a function of pT");
     h_ideff_pt_medium[ieta]->GetXaxis()->SetTitle("p_{T}(jet) [GeV]");
     h_ideff_pt_medium[ieta]->GetYaxis()->SetTitle("Medium jet ID efficiency");


     h_ideff_pt_tight[ieta]  = new TGraphAsymmErrors();
     h_ideff_pt_tight[ieta]->SetNameTitle(Form("ideff_pt_tight%d",ieta),
			       "Jet tight ID efficiency as a function of pT");
     h_ideff_pt_tight[ieta]->GetXaxis()->SetTitle("p_{T}(jet) [GeV]");
     h_ideff_pt_tight[ieta]->GetYaxis()->SetTitle("Tight jet ID efficiency");

     // as a function of the number of good vertices
     
     h_ngood_all[ieta] = (TH1F*)h_nvtx_template->Clone(
						       Form("h_ngood_all%d",ieta));
     h_ngood_pass[ieta] = (TH1F*)h_nvtx_template->Clone(
						       Form("h_ngood_pass%d",ieta));
     h_ideff_ngood[ieta] = new TGraphAsymmErrors();
     h_ideff_ngood[ieta]->SetNameTitle(Form("ideff_ngood%d",ieta),
				      "ID efficiency as a function of nvertex");
     h_ideff_ngood[ieta]->GetXaxis()->SetTitle("Number of good vertices");
     h_ideff_ngood[ieta]->GetYaxis()->SetTitle("Efficiency");


   } // end of loop over eta bins


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      // check the number of good vertices per event
      Int_t ngood_vtx=nGoodVtx(ientry);
      h_ngood->Fill(ngood_vtx);

      // first check which reco jet is the one from the highest and 
      // second et gen jet
      Int_t leadingJetIndex=-1;
      Float_t genJetMaxPt=-9999.;

      Int_t secondLeadingJetIndex=-1;
      Float_t secondJetMaxPt=-9999.;

      for(int ijet=0; ijet < nJet; ijet++){

	if(jetGenJetIndex[ijet]<1)continue;
	if(jetGenJetPt[ijet] > genJetMaxPt)
	  {
	    secondJetMaxPt = genJetMaxPt;
	    secondLeadingJetIndex = leadingJetIndex;
	    genJetMaxPt = jetGenJetPt[ijet];
	    leadingJetIndex= ijet;
	  }
	else if(jetGenJetPt[ijet] > secondJetMaxPt)
	  {
	    secondJetMaxPt = jetGenJetPt[ijet];
	    secondLeadingJetIndex = ijet;
	  }


      } // end of loop over jets

      
      // now check the reconstruction level
      for(int ijet=0; ijet< nJet; ijet++){

	  // matched or not to generator-level jets
 	  if(jetGenJetIndex[ijet]<1)continue;
  	  if(!isFidJet(ientry, ijet))continue;

	  h_genpt->Fill(jetGenJetPt[ijet]);

 	  Float_t jetPtRatio = jetPt[ijet]/jetGenJetPt[ijet];

	  Float_t jetRawPtCorr = jetGenJetPt[ijet]/jetRawPt[ijet];

	  pf_alljet_ratio->Fill(jetGenJetPt[ijet],jetPtRatio);
	  pf_alljet_corr_eta->Fill(jetEta[ijet],jetRawPtCorr);

	  if(ijet==leadingJetIndex)
	    {
	      pf_1stjet_ratio->Fill(jetGenJetPt[ijet],jetPtRatio);
	      pf_1stjet_corr_eta->Fill(jetEta[ijet],jetRawPtCorr);
	    }
	  else if(ijet==secondLeadingJetIndex)
	    {
	      pf_2ndjet_ratio->Fill(jetGenJetPt[ijet],jetPtRatio);
	      pf_2ndjet_corr_eta->Fill(jetEta[ijet],jetRawPtCorr);
	    }

	  
	  int etaIndex=-1;
	  etaIndex = h_dec->GetXaxis()->FindBin(fabs(jetEta[ijet]))-1;
	  if(etaIndex < 0 || etaIndex > 1)continue;

	  h_recetaall->Fill(jetEta[ijet]);
	  h_recptall[etaIndex]->Fill(jetPt[ijet]);

	  h_ngood_all[etaIndex]->Fill(ngood_vtx);

	  if(!isGoodLooseJet(ientry, ijet))continue;

	  h_recetapass->Fill(jetEta[ijet]);
	  h_recptpass[etaIndex]->Fill(jetPt[ijet]);

	  if(!isGoodMediumJet(ientry, ijet))continue;
	  h_recetapass_medium->Fill(jetEta[ijet]);
	  h_recptpass_medium[etaIndex]->Fill(jetPt[ijet]);

	  if(!isGoodTightJet(ientry, ijet))continue;
	  h_recetapass_tight->Fill(jetEta[ijet]);
	  h_recptpass_tight[etaIndex]->Fill(jetPt[ijet]);
	  h_ngood_pass[etaIndex]->Fill(ngood_vtx);	
      }      
      
   } // end of loop over entries

   // now save the histogram in a root file
 
   
   h_ideff_eta->BayesDivide(h_recetapass,h_recetaall,"v");
   h_ideff_eta_medium->BayesDivide(h_recetapass_medium,h_recetaall,"v");
   h_ideff_eta_tight->BayesDivide(h_recetapass_tight,h_recetaall,"v");

   for(int ieta=0; ieta < NETA; ieta++)
     {
       h_ideff_pt[ieta]->BayesDivide(h_recptpass[ieta],h_recptall[ieta],"v");
       h_ideff_pt_medium[ieta]->BayesDivide(h_recptpass_medium[ieta],h_recptall[ieta],"v");
       h_ideff_pt_tight[ieta]->BayesDivide(h_recptpass_tight[ieta],h_recptall[ieta],"v");

       h_ideff_ngood[ieta]->BayesDivide(h_ngood_pass[ieta],
					h_ngood_all[ieta], "v");
     }
  
   TFile* outFile = new TFile(Form("jeteff_%s",inputFile_.data()),"recreate");               
   h_ngood->Write();
   h_recetaall ->Write();
   h_recetapass ->Write();
   h_recetapass_medium ->Write();
   h_recetapass_tight ->Write();
   h_ideff_eta->Write();
   h_ideff_eta_medium->Write();
   h_ideff_eta_tight->Write();
   pf_alljet_ratio->Write();
   pf_1stjet_ratio->Write();
   pf_2ndjet_ratio->Write();
   pf_alljet_corr_eta->Write();
   pf_1stjet_corr_eta->Write();
   pf_2ndjet_corr_eta->Write();
   h_genpt->Write();

   for(int ieta=0;ieta < NETA; ieta++){
     h_recptall[ieta] ->Write();
     h_recptpass[ieta] ->Write();
     h_recptpass_medium[ieta] ->Write();
     h_recptpass_tight[ieta] ->Write();
     h_ideff_pt[ieta]->Write();
     h_ideff_pt_medium[ieta]->Write();
     h_ideff_pt_tight[ieta]->Write();
     h_ngood_all[ieta]->Write();
     h_ngood_pass[ieta]->Write();
     h_ideff_ngood[ieta]->Write();

   }
   outFile->Close();      
}


Bool_t jetEff::isFidJet (Long64_t entry, Int_t ijet)
{
//   if(jetPt[ijet] < 10.0)return false;
  if(jetRawPt[ijet] < 30.0)return false;
  if(fabs(jetEta[ijet]) > 3.0)return false;
  return true;

}


// check if this reco-jet is a good loose jet
Bool_t jetEff::isGoodLooseJet(Long64_t entry, Int_t ijet)
{
  if(jetPt[ijet] < 10.0)return false;
  if(fabs(jetEta[ijet]) > 3.0)return false;
  if(jetNConstituents[ijet] <= 1)return false;
  if(jetNHF[ijet] >= 0.99)return false;
  if(jetNEF[ijet] >= 0.99)return false;

  // for the tracker region
  if(fabs(jetEta[ijet])<2.4 && jetCHF[ijet] <= 0.0)return false;
  if(fabs(jetEta[ijet])<2.4 && jetCEF[ijet] >= 0.99)return false;
  if(fabs(jetEta[ijet])<2.4 && jetNCH[ijet] <= 0)return false;

  return true;

}


// check if this reco-jet is a good medium jet
Bool_t jetEff::isGoodMediumJet(Long64_t entry, Int_t ijet)
{
  if(jetPt[ijet] < 10.0)return false;
  if(fabs(jetEta[ijet]) > 3.0)return false;
  if(jetNConstituents[ijet] <= 1)return false;
  if(jetNHF[ijet] >= 0.95)return false;
  if(jetNEF[ijet] >= 0.95)return false;

  // for the tracker region
  if(fabs(jetEta[ijet])<2.4 && jetCHF[ijet] <= 0.0)return false;
  if(fabs(jetEta[ijet])<2.4 && jetCEF[ijet] >= 0.99)return false;
  if(fabs(jetEta[ijet])<2.4 && jetNCH[ijet] <= 0)return false;

  return true;

}

// check if this reco-jet is a good medium jet
Bool_t jetEff::isGoodTightJet(Long64_t entry, Int_t ijet)
{
  if(jetPt[ijet] < 10.0)return false;
  if(fabs(jetEta[ijet]) > 3.0)return false;
  if(jetNConstituents[ijet] <= 1)return false;
  if(jetNHF[ijet] >= 0.90)return false;
  if(jetNEF[ijet] >= 0.90)return false;

  // for the tracker region
  if(fabs(jetEta[ijet])<2.4 && jetCHF[ijet] <= 0.0)return false;
  if(fabs(jetEta[ijet])<2.4 && jetCEF[ijet] >= 0.99)return false;
  if(fabs(jetEta[ijet])<2.4 && jetNCH[ijet] <= 0)return false;

  return true;

}



Int_t jetEff::nGoodVtx (Long64_t entry)
{
  Int_t ngood=0;
  for(int iv=0; iv < nVtx; iv++)
    {
      if(vtxNDF[iv] > 4 && fabs(vtxD0[iv]) < 2.0 && fabs(vtx[iv][2])<24.0)
	ngood++;
    }
  
  return ngood;


}
