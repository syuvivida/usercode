#define vectorJetEff_cxx
#include "vectorJetEff.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TProfile.h>
#include <TMath.h>

const double ENDCAP_MINETA=2.6;
const double ENDCAP_MAXETA=3.0;
const Int_t   NETA=2;
const double fEtaBin[]={0,ENDCAP_MINETA,ENDCAP_MAXETA};

void vectorJetEff::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   // declare and define histograms
   TH1D* h_dec = new TH1D("h_dec","",2,fEtaBin);


   // adding efficiency 
   TH1D* h_eta_template = new TH1D("h_eta_template","",100,-5.0,5.0);
   h_eta_template->SetXTitle("#eta(jet)");
   TH1D* h_recetaall = (TH1D*)h_eta_template->Clone("h_recetaall");
   TH1D* h_recetapass = (TH1D*)h_eta_template->Clone("h_recetapass");
   TH1D* h_recetapass_medium = (TH1D*)h_eta_template->Clone("h_recetapass_medium");
   TH1D* h_recetapass_tight = (TH1D*)h_eta_template->Clone("h_recetapass_tight");
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


   TH1D* h_pt_template = new TH1D("h_pt_template","",50,0,200);
   h_pt_template->SetXTitle("p_{T}(jet) [GeV]");

   TH1D* h_genpt = (TH1D*)h_pt_template->Clone("h_genpt");

   //0: barrel, 1: endcap, for ID efficiency
   TH1D* h_recptall[NETA];
   TH1D* h_recptpass[NETA];
   TH1D* h_recptpass_medium[NETA];
   TH1D* h_recptpass_tight[NETA];
   TGraphAsymmErrors* h_ideff_pt[NETA];
   TGraphAsymmErrors* h_ideff_pt_medium[NETA];
   TGraphAsymmErrors* h_ideff_pt_tight[NETA];


   // profile template
   TProfile* pf_pt_template = new TProfile("pf_pt_template","",400,0,200,-10.0,10.0);
   pf_pt_template->SetXTitle("p_{T}(jet) [GeV]");
   pf_pt_template->SetYTitle("p_{T}^{RECO}(jet)/p_{T}^{GEN}(jet)"); 

   TProfile* pf_alljet_ratio = (TProfile*)pf_pt_template->Clone("pf_alljet_ratio");
   TProfile* pf_1stjet_ratio = (TProfile*)pf_pt_template->Clone("pf_1stjet_ratio");
   TProfile* pf_2ndjet_ratio = (TProfile*)pf_pt_template->Clone("pf_2ndjet_ratio");

   TProfile* pf_eta_template = new TProfile("pf_eta_template","",100,-5.0,5.0,-10.0,10.0);
   pf_eta_template->SetXTitle("#eta(jet)");
   pf_eta_template->SetYTitle("p_{T}^{GEN}(jet)/p_{T}^{RAW}(jet)"); 

   TProfile* pf_alljet_corr_eta = (TProfile*)pf_eta_template->Clone("pf_alljet_corr_eta");
   TProfile* pf_1stjet_corr_eta = (TProfile*)pf_eta_template->Clone("pf_1stjet_corr_eta");
   TProfile* pf_2ndjet_corr_eta = (TProfile*)pf_eta_template->Clone("pf_2ndjet_corr_eta");

   TProfile* pf_alljet_ratio_eta = (TProfile*)pf_eta_template->Clone("pf_alljet_ratio_eta");
   pf_alljet_ratio_eta->SetYTitle("p_{T}^{RECO}(jet)/p_{T}^{GEN}(jet)"); 
   TProfile* pf_1stjet_ratio_eta = (TProfile*)pf_eta_template->Clone("pf_1stjet_ratio_eta");
   pf_1stjet_ratio_eta->SetYTitle("p_{T}^{RECO}(jet)/p_{T}^{GEN}(jet)"); 
   TProfile* pf_2ndjet_ratio_eta = (TProfile*)pf_eta_template->Clone("pf_2ndjet_ratio_eta");
   pf_2ndjet_ratio_eta->SetYTitle("p_{T}^{RECO}(jet)/p_{T}^{GEN}(jet)"); 					   

   for(int ieta=0; ieta< NETA; ieta++){
   
     h_recptall[ieta]  = (TH1D*)h_pt_template->Clone(Form("h_recptall%d",ieta));
     h_recptpass[ieta] = (TH1D*)h_pt_template->Clone(Form("h_recptpass%d",ieta));
     h_recptpass_medium[ieta] = (TH1D*)h_pt_template->Clone(Form("h_recptpass_medium%d",ieta));
     h_recptpass_tight[ieta] = (TH1D*)h_pt_template->Clone(Form("h_recptpass_tight%d",ieta));
   
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


   } // end of loop over eta bins


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;


      // first check which reco jet is the one from the highest and 
      // second et gen jet
      Int_t leadingJetIndex=-1;
      double genJetMaxPt=-9999.;

      Int_t secondLeadingJetIndex=-1;
      double secondJetMaxPt=-9999.;

      for(unsigned int k=0; k < genJetPt->size(); k++){

	if(genJetPt->at(k) > genJetMaxPt)
	  {
	    secondJetMaxPt = genJetMaxPt;
	    secondLeadingJetIndex = leadingJetIndex;
	    genJetMaxPt = genJetPt->at(k);
	    leadingJetIndex= k;
	  }
	else if(genJetPt->at(k) > secondJetMaxPt)
	  {
	    secondJetMaxPt = genJetPt->at(k);
	    secondLeadingJetIndex = k;
	  }


      } // end of loop over jets


     
      
      // now check the reconstruction level
      for(int ijet=0; ijet< patJetPFlowPt_->size(); ijet++){

	  // matched or not to generator-level jets
	  if(!isFidJet(ientry, ijet))continue;
	  
	  int genIndex =  matchedGenJet(ientry, ijet);
	  if(genIndex <0)continue;

	  double etaRec = patJetPFlowEta_->at(ijet);
	  double ptRec  = patJetPFlowPt_->at(ijet);
	  double ptTruth = genJetPt->at(genIndex);
	  double etaTruth = genJetEta->at(genIndex);
	  
	  if(ptTruth<20.0)continue;
	  if(fabs(etaTruth)>3.0)continue;
	  

	  h_genpt->Fill(ptTruth);


 	  double jetPtRatio = patJetPFlowPt_->at(ijet)/ptTruth;

	  double jetRawPtCorr = ptTruth/patJetPFlowUnCorrPt_->at(ijet);


	  pf_alljet_ratio->Fill(ptTruth,jetPtRatio);
	  pf_alljet_ratio_eta->Fill(etaTruth,jetPtRatio);

	  pf_alljet_corr_eta->Fill(etaRec,jetRawPtCorr);

	  if(genIndex==leadingJetIndex)
	    {
	      pf_1stjet_ratio->Fill(ptTruth,jetPtRatio);
	      pf_1stjet_ratio_eta->Fill(etaTruth,jetPtRatio);
	      pf_1stjet_corr_eta->Fill(etaRec,jetRawPtCorr);
	    }
	  else if(genIndex==secondLeadingJetIndex)
	    {
	      pf_2ndjet_ratio->Fill(ptTruth,jetPtRatio);
	      pf_2ndjet_ratio_eta->Fill(etaTruth,jetPtRatio);
	      pf_2ndjet_corr_eta->Fill(etaRec,jetRawPtCorr);
	    }

	  
	  int etaIndex=-1;
	  etaIndex = h_dec->GetXaxis()->FindBin(fabs(etaRec))-1;
	  if(etaIndex < 0 || etaIndex > 1)continue;

	  h_recetaall->Fill(etaRec);
	  h_recptall[etaIndex]->Fill(ptRec);

	  if(!isGoodLooseJet(ientry, ijet))continue;

	  h_recetapass->Fill(etaRec);
	  h_recptpass[etaIndex]->Fill(ptRec);

	  if(!isGoodMediumJet(ientry, ijet))continue;
	  h_recetapass_medium->Fill(etaRec);
	  h_recptpass_medium[etaIndex]->Fill(ptRec);

	  if(!isGoodTightJet(ientry, ijet))continue;
	  h_recetapass_tight->Fill(etaRec);
	  h_recptpass_tight[etaIndex]->Fill(ptRec);
      }      // end of loop over jets
      

      
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

     }
  
   std::string remword="/data2/syu/anilNtuple/";
   size_t pos = inputFile_.find(remword);
  
   if(pos!= std::string::npos)
     inputFile_.swap(inputFile_.erase(pos,remword.length()));
   
   TFile* outFile = new TFile(Form("vectorJetEff_%s", 
				   inputFile_.data()),"recreate");        

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
   pf_alljet_ratio_eta->Write();
   pf_1stjet_ratio_eta->Write();
   pf_2ndjet_ratio_eta->Write();
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

   }
   outFile->Close();  

}


Bool_t vectorJetEff::isFidJet (Long64_t entry, Int_t ijet)
{
  if(patJetPFlowPt_->at(ijet) < 10.0)return false;
  if(fabs(patJetPFlowEta_->at(ijet)) > 3.0)return false;
  return true;

}


// check if this reco-jet is a good loose jet
Bool_t vectorJetEff::isGoodLooseJet(Long64_t entry, Int_t ijet)
{
  if(patJetPFlowPt_->at(ijet) < 10.0)return false;
  if(fabs(patJetPFlowEta_->at(ijet)) > 3.0)return false;
  if(patJetPFlowNConstituents_->at(ijet) <= 1)return false;
  if(patJetPFlowNeutHadEFr_->at(ijet) >= 0.99)return false;
  if(patJetPFlowNeutEmEFr_->at(ijet) >= 0.99)return false;

//   // for the tracker region
  if(fabs(patJetPFlowEta_->at(ijet))<2.4 && patJetPFlowCharHadEFr_->at(ijet) <= 0.0)return false;
  if(fabs(patJetPFlowEta_->at(ijet))<2.4 && patJetPFlowCharEmEFr_->at(ijet) >= 0.99)return false;
  if(fabs(patJetPFlowEta_->at(ijet))<2.4 && patJetPFlowCharMulti_->at(ijet) <= 0)return false;

  return true;

}


// check if this reco-jet is a good medium jet
Bool_t vectorJetEff::isGoodMediumJet(Long64_t entry, Int_t ijet)
{
  if(patJetPFlowPt_->at(ijet) < 10.0)return false;
  if(fabs(patJetPFlowEta_->at(ijet)) > 3.0)return false;
  if(patJetPFlowNConstituents_->at(ijet) <= 1)return false;
  if(patJetPFlowNeutHadEFr_->at(ijet) >= 0.95)return false;
  if(patJetPFlowNeutEmEFr_->at(ijet) >= 0.95)return false;

//   // for the tracker region
  if(fabs(patJetPFlowEta_->at(ijet))<2.4 && patJetPFlowCharHadEFr_->at(ijet) <= 0.0)return false;
  if(fabs(patJetPFlowEta_->at(ijet))<2.4 && patJetPFlowCharEmEFr_->at(ijet) >= 0.99)return false;
  if(fabs(patJetPFlowEta_->at(ijet))<2.4 && patJetPFlowCharMulti_->at(ijet) <= 0)return false;


  return true;

}

// check if this reco-jet is a good medium jet
Bool_t vectorJetEff::isGoodTightJet(Long64_t entry, Int_t ijet)
{
  if(patJetPFlowPt_->at(ijet) < 10.0)return false;
  if(fabs(patJetPFlowEta_->at(ijet)) > 3.0)return false;
  if(patJetPFlowNConstituents_->at(ijet) <= 1)return false;
  if(patJetPFlowNeutHadEFr_->at(ijet) >= 0.90)return false;
  if(patJetPFlowNeutEmEFr_->at(ijet) >= 0.90)return false;


//   // for the tracker region
  if(fabs(patJetPFlowEta_->at(ijet))<2.4 && patJetPFlowCharHadEFr_->at(ijet) <= 0.0)return false;
  if(fabs(patJetPFlowEta_->at(ijet))<2.4 && patJetPFlowCharEmEFr_->at(ijet) >= 0.99)return false;
  if(fabs(patJetPFlowEta_->at(ijet))<2.4 && patJetPFlowCharMulti_->at(ijet) <= 0)return false;

  return true;

}


Int_t  vectorJetEff::matchedGenJet(Long64_t entry, Int_t ijet)
{
  
  int matchedGenIndex = -1;
  for(int k=0; k< genJetPt->size(); k++)
    {

      double deta = genJetEta->at(k)-patJetPFlowEta_->at(ijet);
      double dphi = genJetPhi->at(k)-patJetPFlowPhi_->at(ijet);      
      while (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
      while (dphi <= -TMath::Pi()) dphi += 2*TMath::Pi();

      double dR = sqrt(deta*deta+dphi*dphi);

      double relPt = genJetPt->at(k)>1e-6? 
	fabs(genJetPt->at(k)-patJetPFlowPt_->at(ijet))/genJetPt->at(k): -9999.0;

      if(dR<0.4 && relPt < 3.0)
	{
	  matchedGenIndex = k;
	  break;
	}
      
    } // end of loop over generator-level jets
      
  return matchedGenIndex;
}

