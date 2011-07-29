#define phoEff_cxx
#include "phoEff.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>

const Float_t BARREL_MAXETA=1.4442;
const Float_t ENDCAP_MINETA=1.566;
const Float_t ENDCAP_MAXETA=2.5;
const Int_t   NETA=2;
const Float_t fEtaBin[]={0,ENDCAP_MINETA,ENDCAP_MAXETA};

void phoEff::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TH1F* h_dec = new TH1F("h_dec","",2,fEtaBin);
   TH1F* h_ngood = new TH1F("h_ngood","",25,-0.5,24.5);

   TH1F* h_eta_template = new TH1F("h_eta_template","",60,-3.0,3.0);
   h_eta_template->SetXTitle("#eta(#gamma)");

   TH1F* h_genetaall = (TH1F*)h_eta_template->Clone("h_genetaall");
   TH1F* h_genetapass = (TH1F*)h_eta_template->Clone("h_genetapass");
   TGraphAsymmErrors* h_recoeff_eta = new TGraphAsymmErrors();
   h_recoeff_eta->SetNameTitle("recoeff_eta", 
			       "Reconstruction efficiency as a function of eta");
   h_recoeff_eta->GetXaxis()->SetTitle("#eta(#gamma)");
   h_recoeff_eta->GetYaxis()->SetTitle("Efficiency");

   TH1F* h_et_template = new TH1F("h_et_template","",50,0,200);
   h_et_template->SetXTitle("E_{T}(#gamma) [GeV]");

   // 0: barrel, 1: endcap, for reconstruction efficiency
   TH1F* h_genetall[NETA];
   TH1F* h_genetpass[NETA];
   TGraphAsymmErrors* h_recoeff_et[NETA];

   //0: barrel, 1: endcap, for ID efficiency
   TH1F* h_recetall[NETA];
   TH1F* h_recetpass[NETA];
   TGraphAsymmErrors* h_ideff_et[NETA];

   for(int ieta=0; ieta< NETA; ieta++){
   
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

   } // end of loop over eta bins


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      int ngood_vtx=nGoodVtx(ientry);
      h_ngood->Fill(ngood_vtx);
//       if(ngood_vtx<=2)continue;

      
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
	  if(isGoodPho(ientry, ipho))
	    h_recetpass[etaIndex]->Fill(phoEt[ipho]);

	  if(find)break;
	  
	} // end of loop over reco photons
      
      } // end of loop over MC particles


   } // end of loop over entries

   // now save the histogram in a root file
 
   
   h_recoeff_eta->BayesDivide(h_genetapass,h_genetaall,"v");
   for(int ieta=0; ieta < NETA; ieta++)
     h_recoeff_et[ieta]->BayesDivide(h_genetpass[ieta],h_genetall[ieta],"v");

   for(int ieta=0; ieta < NETA; ieta++)
     h_ideff_et[ieta]->BayesDivide(h_recetpass[ieta],h_recetall[ieta],"v");
  
   TFile* outFile = new TFile("output_eff.root","recreate");               
   h_ngood->Write();
   h_genetaall ->Write();
   h_genetapass ->Write();
   h_recoeff_eta->Write();

   for(int ieta=0;ieta < NETA; ieta++){
     h_genetall[ieta] ->Write();
     h_genetpass[ieta] ->Write();
     h_recoeff_et[ieta]->Write();

     h_recetall[ieta] ->Write();
     h_recetpass[ieta] ->Write();
     h_ideff_et[ieta]->Write();

   }
   outFile->Close();      

}

// check if this reco-photon is a good photon
Bool_t phoEff::isGoodPho(Long64_t entry, Int_t ipho)
{


  bool isEB=false;
  bool isEE=false;
  if(fabs(phoSCEta[ipho]) < BARREL_MAXETA)isEB=true;
  if(fabs(phoSCEta[ipho]) > ENDCAP_MINETA && 
     fabs(phoSCEta[ipho]) < ENDCAP_MAXETA) isEE=true;

  if(!isEB && !isEE)return false;
  if(phoEt[ipho] < 15.0)return false;
  if(phoHoverE[ipho] > 0.05)return false;
  if(phohasPixelSeed[ipho]==1)return false;
  if(phoEcalIsoDR04[ipho] > 4.2 +0.006 * phoEt[ipho])return false;
  if(phoHcalIsoDR04[ipho] > 2.2 +0.0025* phoEt[ipho])return false;
  if(phoTrkIsoHollowDR04[ipho] > 2.0 +0.001* phoEt[ipho])return false;
  if(isEB && phoSigmaIEtaIEta[ipho] > 0.013)return false;
  if(isEE && phoSigmaIEtaIEta[ipho] > 0.030)return false;

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
