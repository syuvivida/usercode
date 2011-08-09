#define phoEff_cxx
#include "phoEff.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TProfile.h>

const Float_t BARREL_MAXETA=1.4442;
const Float_t ENDCAP_MINETA=1.566;
const Float_t ENDCAP_MAXETA=2.5;
const Int_t   NETA=2;
const Float_t fEtaBin[]={0,ENDCAP_MINETA,ENDCAP_MAXETA};
const Float_t rhoCorr_higgs[]={0.175,0.215};
const Float_t rhoCorr[]={0.267,0.267};

void phoEff::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TH1F* h_dec = new TH1F("h_dec","",2,fEtaBin);

   TH1F* h_nvtx_template = new TH1F("h_nvtx_template","",30,-0.5,29.5);
   TH1F* h_ngood  = (TH1F*)h_nvtx_template->Clone("h_ngood");

   TH1F* h_eta_template = new TH1F("h_eta_template","",60,-3.0,3.0);
   h_eta_template->SetXTitle("#eta(#gamma)");
   TH1F* h_genetaall = (TH1F*)h_eta_template->Clone("h_genetaall");
   TH1F* h_genetapass = (TH1F*)h_eta_template->Clone("h_genetapass");
   TGraphAsymmErrors* h_recoeff_eta = new TGraphAsymmErrors();
   h_recoeff_eta->SetNameTitle("recoeff_eta", 
			       "Reconstruction efficiency as a function of eta");
   h_recoeff_eta->GetXaxis()->SetTitle("#eta(#gamma)");
   h_recoeff_eta->GetYaxis()->SetTitle("Efficiency");

   TH1F* h_et_template = new TH1F("h_et_template","",100,0,500);
   h_et_template->SetXTitle("E_{T}(#gamma) [GeV]");

   // 0: barrel, 1: endcap, for reconstruction efficiency
   TH1F* h_genetall[NETA];
   TH1F* h_genetpass[NETA];
   TGraphAsymmErrors* h_recoeff_et[NETA];

   //0: barrel, 1: endcap, for ID efficiency
   TH1F* h_recetall[NETA];
   TH1F* h_recetpass[NETA];
   TGraphAsymmErrors* h_ideff_et[NETA];


   // checking efficiency vs # of good vertices
   TH1F* h_ngood_all[NETA];
   TH1F* h_ngood_pass[NETA];
   TGraphAsymmErrors* h_ideff_ngood[NETA];

   // checking mean isolation vs # of good vertices
   TProfile* pf_nvtx_template = 
     new TProfile("pf_nvtx_template","",30,-0.5,29.5,-20,50);
//     new TProfile("pf_nvtx_template","",30,-0.5,29.5,-10,20);
   pf_nvtx_template->SetXTitle("Number of good vertices");
   TProfile* pf_eciso04[NETA];
   TProfile* pf_eciso03[NETA];
   TProfile* pf_hciso04[NETA];
   TProfile* pf_hciso03[NETA];
   TProfile* pf_tkiso04[NETA];
   TProfile* pf_tkiso03[NETA];
   TProfile* pf_siso[NETA];
   TProfile* pf_sisocorr[NETA];
   TProfile* pf_siso_higgs[NETA];
   TProfile* pf_sisocorr_higgs[NETA];

   TH1F* h_iso_template = new TH1F("h_iso_template","",100,-2.0,23.0);
   h_iso_template->SetXTitle("Combined isolation [GeV]");
   TH1F* h_iso[NETA];

   for(int ieta=0; ieta< NETA; ieta++){

     // for the isolation template
     h_iso[ieta] = (TH1F*)h_iso_template->Clone(Form("h_iso%d",ieta));
   
     // as a function of Et
     h_genetall[ieta]   = 
       (TH1F*)h_et_template->Clone(Form("h_genetall%d",ieta));
     h_genetpass[ieta] = 
       (TH1F*)h_et_template->Clone(Form("h_genetpass%d",ieta));
   
     h_recoeff_et[ieta]  = new TGraphAsymmErrors();
     h_recoeff_et[ieta]->SetNameTitle(Form("recoeff_et%d",ieta),
				      "Reconstruction efficiency as a function of ET");
     h_recoeff_et[ieta]->GetXaxis()->SetTitle("E_{T}(#gamma) [GeV]");
     h_recoeff_et[ieta]->GetYaxis()->SetTitle("Efficiency");

     h_recetall[ieta]  = (TH1F*)h_et_template->Clone(Form("h_recetall%d",ieta));
     h_recetpass[ieta] = (TH1F*)h_et_template->Clone(Form("h_recetpass%d",ieta));
   
     h_ideff_et[ieta]  = new TGraphAsymmErrors();
     h_ideff_et[ieta]->SetNameTitle(Form("ideff_et%d",ieta),
				      "ID efficiency as a function of ET");
     h_ideff_et[ieta]->GetXaxis()->SetTitle("E_{T}(#gamma) [GeV]");
     h_ideff_et[ieta]->GetYaxis()->SetTitle("Efficiency");

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
     

     pf_eciso04[ieta] = (TProfile*)pf_nvtx_template->Clone(
			     Form("pf_eciso04%d",ieta));
     pf_eciso04[ieta] -> SetYTitle("Average ECAL #DeltaR=0.4 isolation [GeV]");

     pf_eciso03[ieta] = (TProfile*)pf_nvtx_template->Clone(
			     Form("pf_eciso03%d",ieta));
     pf_eciso03[ieta] -> SetYTitle("Average ECAL #DeltaR=0.3 isolation [GeV]");

     pf_hciso04[ieta] = (TProfile*)pf_nvtx_template->Clone(
			     Form("pf_hciso04%d",ieta));
     pf_hciso04[ieta] -> SetYTitle("Average HCAL #DeltaR=0.4 isolation [GeV]");

     pf_hciso03[ieta] = (TProfile*)pf_nvtx_template->Clone(
			     Form("pf_hciso03%d",ieta));
     pf_hciso03[ieta] -> SetYTitle("Average HCAL #DeltaR=0.3 isolation [GeV]");


     pf_tkiso04[ieta] = (TProfile*)pf_nvtx_template->Clone(
			     Form("pf_tkiso04%d",ieta));
     pf_tkiso04[ieta] -> SetYTitle("Average tracker #DeltaR=0.4 isolation [GeV]");

     pf_tkiso03[ieta] = (TProfile*)pf_nvtx_template->Clone(
			     Form("pf_tkiso03%d",ieta));
     pf_tkiso03[ieta] -> SetYTitle("Average tracker #DeltaR=0.3 isolation [GeV]");

     pf_siso[ieta] = (TProfile*)pf_nvtx_template->Clone(
	    	     Form("pf_siso%d",ieta));
     pf_siso[ieta]-> SetYTitle("Sum isolation [GeV]");

     pf_sisocorr[ieta] = (TProfile*)pf_nvtx_template->Clone(
				  Form("pf_sisocorr%d",ieta));
     pf_sisocorr[ieta]-> SetYTitle("Corrected sum isolation [GeV]");

     pf_siso_higgs[ieta] = (TProfile*)pf_nvtx_template->Clone(
			     Form("pf_siso_higgs%d",ieta));
     pf_siso_higgs[ieta]-> SetYTitle("Sum isolation [GeV]");

     pf_sisocorr_higgs[ieta] = (TProfile*)pf_nvtx_template->Clone(
								  Form("pf_sisocorr_higgs%d",ieta));
     pf_sisocorr_higgs[ieta]-> SetYTitle("Corrected sum isolation [GeV]");

   } // end of loop over eta bins

   float eff_deno[2]={0};
   float eff_numr[2]={0};

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      int ngood_vtx=nGoodVtx(ientry);
      if(ngood_vtx<1)continue;
      h_ngood->Fill(ngood_vtx);

      
      // calculate reconstruction efficiency first
      for(int imc=0; imc < nMC; imc++){
	
	// check if this is a real photon and that this photon is from 
	// Higgs decays
	if(mcPID[imc]!=22 || mcMomPID[imc]!=22)continue;

	// the photon should be away from the detector boundary
	if(fabs(mcEta[imc]) > 1.44
	   && fabs(mcEta[imc]) < 1.57)continue;
	if(fabs(mcEta[imc]) > ENDCAP_MAXETA)continue;

	// at least 15 GeV
	if(mcPt[imc]<15.0)continue;

	int etaIndex=-1;
	etaIndex = h_dec->GetXaxis()->FindBin(fabs(mcEta[imc]))-1;
	if(etaIndex < 0 || etaIndex > 1)continue;

	// check the denominator
	h_genetaall->Fill(mcEta[imc]);
 	h_genetall[etaIndex]->Fill(mcEt[imc]);

	// now check the reconstruction level
	bool find=false;
	for(int ipho=0; ipho< nPho; ipho++){

	  // matched or not
 	  if(phoGenIndex[ipho]!=mcIndex[imc])continue;
 	  if(!isFidPho(ientry, ipho))continue;
	  find=true;

 	  h_genetapass->Fill(mcEta[imc]);	  
 	  h_genetpass[etaIndex]->Fill(mcEt[imc]);	  

	  h_recetall[etaIndex]->Fill(phoEt[ipho]);
	  h_ngood_all[etaIndex]->Fill(ngood_vtx);

	  pf_eciso04[etaIndex]->Fill(ngood_vtx,phoEcalIsoDR04[ipho]);
	  pf_eciso03[etaIndex]->Fill(ngood_vtx,phoEcalIsoDR03[ipho]);
	  pf_hciso04[etaIndex]->Fill(ngood_vtx,phoHcalIsoDR04[ipho]);
	  pf_hciso03[etaIndex]->Fill(ngood_vtx,phoHcalIsoDR03[ipho]);
	  pf_tkiso04[etaIndex]->Fill(ngood_vtx,phoTrkIsoHollowDR04[ipho]);
	  pf_tkiso03[etaIndex]->Fill(ngood_vtx,phoTrkIsoHollowDR03[ipho]);

	  Float_t sumIso_higgs = phoEcalIsoDR03[ipho]+phoHcalIsoDR04[ipho]
	    + phoCiCTrkIsoDR04[ipho][0];

	  Float_t sumIso = phoEcalIsoDR04[ipho]+phoHcalIsoDR04[ipho]+
	    phoTrkIsoHollowDR04[ipho];

	  pf_siso_higgs[etaIndex]->Fill(ngood_vtx,sumIso_higgs);	      
// 	  pf_sisocorr_higgs[etaIndex]->Fill(ngood_vtx,sumIso_higgs-rhoCorr_higgs[etaIndex]*rho25);
	  pf_sisocorr_higgs[etaIndex]->Fill(ngood_vtx,rho25);
	  
	  pf_siso[etaIndex]->Fill(ngood_vtx,sumIso); 
	  pf_sisocorr[etaIndex]->Fill(ngood_vtx,sumIso-rhoCorr[etaIndex]*rho25);


	  if(isGoodPho(ientry, ipho))
	    {
	      h_recetpass[etaIndex]->Fill(phoEt[ipho]);
	      h_ngood_pass[etaIndex]->Fill(ngood_vtx);	  
	      h_iso[etaIndex]->Fill(sumIso);
	      eff_deno[etaIndex] += 1.0;

	      if(sumIso<5.0)
		eff_numr[etaIndex] += 1.0;
	      
	    }

	  if(find)break;
	  
	} // end of loop over reco photons
      
      } // end of loop over MC particles


   } // end of loop over entries

   // now save the histogram in a root file
 
   
   h_recoeff_eta->BayesDivide(h_genetapass,h_genetaall,"v");
   for(int ieta=0; ieta < NETA; ieta++)
     {
       h_recoeff_et[ieta]->BayesDivide(h_genetpass[ieta],h_genetall[ieta],"v");      
       h_ideff_et[ieta]->BayesDivide(h_recetpass[ieta],h_recetall[ieta],"v");

       h_ideff_ngood[ieta]->BayesDivide(h_ngood_pass[ieta],
					h_ngood_all[ieta], "v");
     }
  
   TFile* outFile = new TFile(Form("phoeff_hovere_%s",inputFile_.data()),"recreate");               
   h_ngood->Write();
   h_genetaall ->Write();
   h_genetapass ->Write();
   h_recoeff_eta->Write();

   for(int ieta=0;ieta < NETA; ieta++){
     
     h_iso[ieta]->Write();

     h_genetall[ieta] ->Write();
     h_genetpass[ieta] ->Write();
     h_recoeff_et[ieta]->Write();

     h_recetall[ieta] ->Write();
     h_recetpass[ieta] ->Write();
     h_ideff_et[ieta]->Write();

     h_ngood_all[ieta]->Write();
     h_ngood_pass[ieta]->Write();
     h_ideff_ngood[ieta]->Write();


     pf_eciso04[ieta]->Write();
     pf_eciso03[ieta]->Write();
     pf_hciso04[ieta]->Write();
     pf_hciso03[ieta]->Write();
     pf_tkiso04[ieta]->Write();
     pf_tkiso03[ieta]->Write();

     pf_siso[ieta]->Write();
     pf_sisocorr[ieta]->Write();

     pf_siso_higgs[ieta]->Write();
     pf_sisocorr_higgs[ieta]->Write();
     cout << "efficiency " << ieta << " = " << eff_numr[ieta]/eff_deno[ieta]<< endl;

   }
   outFile->Close();      


}

// check if this reco-photon is a good photon
Bool_t phoEff::isGoodPho(Long64_t entry, Int_t ipho)
{


  Bool_t isEB=false;
  Bool_t isEE=false;
  if(fabs(phoSCEta[ipho]) < BARREL_MAXETA)isEB=true;
  if(fabs(phoSCEta[ipho]) > ENDCAP_MINETA && 
     fabs(phoSCEta[ipho]) < ENDCAP_MAXETA) isEE=true;

  if(!isEB && !isEE)return false;
  if(phoEt[ipho] < 15.0)return false;
  if(phoHoverE[ipho] > 0.05)return false;
  if(phohasPixelSeed[ipho]==1)return false;
//   if(phoEcalIsoDR04[ipho] > 4.2 +0.006 * phoEt[ipho])return false;
//   if(phoHcalIsoDR04[ipho] > 2.2 +0.0025* phoEt[ipho])return false;
//   if(phoTrkIsoHollowDR04[ipho] > 2.0 +0.001* phoEt[ipho])return false;
  if(phoEcalIsoDR04[ipho] > 4.2 +0.003 * phoEt[ipho])return false;
  if(phoHcalIsoDR04[ipho] > 2.2 +0.001* phoEt[ipho])return false;
  if(phoTrkIsoHollowDR04[ipho] > 2.0 +0.001* phoEt[ipho])return false;
  
//   Float_t sumIso = phoEcalIsoDR04[ipho]+phoHcalIsoDR04[ipho]+
//     phoTrkIsoHollowDR04[ipho];

//   sumIso -= isEB? (rhoCorr[0]*rho25): (rhoCorr[1]*rho25);

//   if(sumIso > 5.0)return false;
//   if(isEB && phoSigmaIEtaIEta[ipho] > 0.013)return false;
//   if(isEE && phoSigmaIEtaIEta[ipho] > 0.030)return false;

  return true;

}

// check if this reco-photon is a good photon
Bool_t phoEff::isFidPho(Long64_t entry, Int_t ipho)
{


  bool isEB=false;
  bool isEE=false;
  if(fabs(phoSCEta[ipho]) < BARREL_MAXETA)isEB=true;
  if(fabs(phoSCEta[ipho]) > ENDCAP_MINETA && 
     fabs(phoSCEta[ipho]) < ENDCAP_MAXETA) isEE=true;

  if(!isEB && !isEE)return false;
  if(phoEt[ipho] < 15.0)return false;

  return true;

}


Int_t phoEff::nGoodVtx (Long64_t entry)
{
  Int_t ngood=0;
  for(int iv=0; iv < nVtx; iv++)
    {
      if(vtxNDF[iv] > 4 && fabs(vtxD0[iv]) < 2.0 && fabs(vtx[iv][2])<24.0)
	ngood++;
    }
  
  return ngood;


}
