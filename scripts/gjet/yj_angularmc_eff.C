#define yj_angularmc_eff_cxx
#include "yj_angularmc_eff.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include "myLib.h"

// const double BARREL_MAXETA=1.4442;
// const double ENDCAP_MINETA=1.566;
const double BARREL_MAXETA=1.44;
const double ENDCAP_MINETA=1.57;
const double ENDCAP_MAXETA=2.5;

void yj_angularmc_eff::Loop()
{
  cout << "This is version 0" << endl;
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  cout << "There are " << nentries << " entries" << endl;

  double ptbound[]={85., 95., 110., 130., 160., 200.};
  const int nPtBins = sizeof(ptbound)/sizeof(ptbound[0])-1;
  cout << "There are " << nPtBins << " pt bins" << endl;

  TH1D* h_pt_template = new TH1D("h_pt_template","",500,0,500);
  TH1D* h_pthat       = (TH1D*)h_pt_template->Clone("h_pthat");
  TH1D* h_ptpho       = (TH1D*)h_pt_template->Clone("h_ptpho");
  TH1D* h_ptjet       = (TH1D*)h_pt_template->Clone("h_ptjet");

  TH1D* h_eta_template = new TH1D("h_eta_template","",100,-5.0,5.0);
  TH1D* h_etapho      = (TH1D*)h_eta_template->Clone("h_etapho");
  TH1D* h_etajet      = (TH1D*)h_eta_template->Clone("h_etajet");

  TH1D* h_ptbin_template= new TH1D("h_ptbin_template","", nPtBins, ptbound);

  TH1D* h_zgamma[2][nPtBins][2];
  TH1D* h_dphi[2][nPtBins][2];
  TH1D* h_cost[2][nPtBins][2];
  TH1D* h_chi[2][nPtBins][2];
  TH1D* h_pstar[2][nPtBins][2];
  TH1D* h_yB[2][nPtBins][2];

  TH1D* h_zgamma_template = new TH1D("h_zgamma_template","",
				     200,-4.0,4.0);
  
  TH1D* h_dphi_template = new TH1D("h_dphi_template","",
				   100,0.0,TMath::Pi());

  TH1D* h_cost_template = new TH1D("h_cost_template","",
				   100,0.0,1.0);

  TH1D* h_chi_template = new TH1D("h_chi_template","",
				  100,0.0,20.0);      
  
  TH1D* h_pstar_template = new TH1D("h_pstar_template","",500,0,500);

  TH1D* h_yB_template = new TH1D("h_yB_template","",100,-5.0,5.0);

  // debugging histograms
  TH1D* h_sieie_template = new TH1D("h_sieie_template","", 100, 0,0.1);
  
  TH1D* h_eciso_template = new TH1D("h_eciso_template","", 250, 0,50.0);

  TH1D* h_hciso_template = new TH1D("h_hciso_template","", 250, 0,50.0);

  TH1D* h_tkiso_template = new TH1D("h_tkiso_template","", 250, 0,50.0);

  TH1D* h_hovere_template = new TH1D("h_hovere_template","", 250, 0,0.5);

  TH1D* h_pixel_template = new TH1D("h_pixel_template","", 3, -0.5,2.5);


  // creating histograms
  char* decName[4]={"EB","EE","leadingEB","leadingEE"};

  TH1D* h_sieie[4];
  TH1D* h_eciso[4];
  TH1D* h_hciso[4];
  TH1D* h_tkiso[4];
  TH1D* h_hovere[4];
  TH1D* h_pixel[4];

  for(int idec=0; idec<4; idec++)
    {
      // EB:idec=0, EE:idec=1, leading phtoon only idec=2
      h_sieie[idec] = (TH1D*)h_sieie_template->Clone(Form("h_sieie_%s", decName[idec]));
      h_eciso[idec] = (TH1D*)h_eciso_template->Clone(Form("h_eciso_%s", decName[idec]));
      h_hciso[idec] = (TH1D*)h_hciso_template->Clone(Form("h_hciso_%s", decName[idec]));
      h_tkiso[idec] = (TH1D*)h_tkiso_template->Clone(Form("h_tkiso_%s", decName[idec]));
      h_hovere[idec] = (TH1D*)h_hovere_template->Clone(Form("h_hovere_%s", decName[idec]));
      h_pixel[idec] = (TH1D*)h_pixel_template->Clone(Form("h_pixel_%s", decName[idec]));    
    }

  for(int idec=0; idec < 2; idec ++){

    for(int ipt=0; ipt < nPtBins; ipt++){
     
      for(int ip=0; ip < 2; ip++){ // 0: before cut, 1: after cut
 
	h_zgamma[idec][ipt][ip] = (TH1D*)h_zgamma_template->Clone(Form("h_zgamma_%s_%d_%d_%d", decName[idec],
								       (int)ptbound[ipt], (int)ptbound[ipt+1], ip));
      
	h_dphi[idec][ipt][ip] = (TH1D*)h_dphi_template->Clone(Form("h_dphi_%s_%d_%d_%d", decName[idec],
								   (int)ptbound[ipt], (int)ptbound[ipt+1], ip));

	h_cost[idec][ipt][ip] = (TH1D*)h_cost_template->Clone(Form("h_cost_%s_%d_%d_%d", decName[idec],
								   (int)ptbound[ipt], (int)ptbound[ipt+1], ip));

	h_chi[idec][ipt][ip] = (TH1D*)h_chi_template->Clone(Form("h_chi_%s_%d_%d_%d", decName[idec],
								 (int)ptbound[ipt], (int)ptbound[ipt+1], ip));

	h_pstar[idec][ipt][ip] = (TH1D*)h_pstar_template->Clone(Form("h_pstar_%s_%d_%d_%d", decName[idec],
								     (int)ptbound[ipt], (int)ptbound[ipt+1], ip));

	h_yB[idec][ipt][ip] = (TH1D*)h_yB_template->Clone(Form("h_yB_%s_%d_%d_%d", decName[idec],
							       (int)ptbound[ipt], (int)ptbound[ipt+1], ip));
	
      } // end of loop over processes
   
    } // end of loop over pt bins

  } // loop over barrel and endcaps


  

  Long64_t nPass[30]={0};

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
//     if(jentry > 50000) break;
    nPass[0]++;
    h_pthat->Fill(PhotonMCpthat); // make sure we are checking the right MC samples


    //vertex selection
    if(EvtInfo_nVtxGood<1) continue;
    nPass[1]++;

    int leadingPhotonIndex = -1;
    double phoMaxPt = -9999.;

     
    // now find a good leading photon that is matched to hard scattering photon
    for(unsigned int ipho=0; ipho < PhotonPt->size(); ipho++)
      {	

	int decIndex = phoDecCode(ientry, ipho);
	if(!isFidPho(ientry,ipho))continue;

	h_sieie[decIndex] -> Fill(PhotonSigmaIetaIeta->at(ipho));
	h_eciso[decIndex] -> Fill(PhotonecalRecHitSumEtConeDR04->at(ipho));
	h_hciso[decIndex] -> Fill(PhotonhcalTowerSumEtConeDR04->at(ipho));
	h_tkiso[decIndex] -> Fill(PhotontrkSumPtHollowConeDR04->at(ipho));
	h_hovere[decIndex]-> Fill(PhotonhadronicOverEm->at(ipho));
	h_pixel[decIndex] -> Fill(PhotonhasPixelSeed->at(ipho));

	// 	if(!PhotonisGenMatched->at(ipho))continue;
	
	// 	cout << PhotongenMomId->size() << "\t" << PhotonPt->size() << endl;
	// 	if(fabs(PhotongenMomId->at(ipho)-22)>1e-6)continue; // not prompt photon

	double thisPhoPt= PhotonEt->at(ipho);      
	if(thisPhoPt > phoMaxPt)
	  {
	    phoMaxPt = thisPhoPt;
	    leadingPhotonIndex= ipho;
	  }

      }
    // end of leading photon search

    if(leadingPhotonIndex<0)continue;
    nPass[2]++;

    double leadingPhotonEt = PhotonEt->at(leadingPhotonIndex);
    int phoPtBinIndex =  h_ptbin_template->GetXaxis()->FindBin(leadingPhotonEt)-1; 

    double leadingPhotonScEta = PhotonScEta->at(leadingPhotonIndex);
    int phoDecBinIndex = phoDecCode(ientry, leadingPhotonIndex);

    h_sieie[phoDecBinIndex+2] -> Fill(PhotonSigmaIetaIeta->at(leadingPhotonIndex));
    h_eciso[phoDecBinIndex+2] -> Fill(PhotonecalRecHitSumEtConeDR04->at(leadingPhotonIndex));
    h_hciso[phoDecBinIndex+2] -> Fill(PhotonhcalTowerSumEtConeDR04->at(leadingPhotonIndex));
    h_tkiso[phoDecBinIndex+2] -> Fill(PhotontrkSumPtHollowConeDR04->at(leadingPhotonIndex));
    h_hovere[phoDecBinIndex+2]-> Fill(PhotonhadronicOverEm->at(leadingPhotonIndex));
    h_pixel[phoDecBinIndex+2] -> Fill(PhotonhasPixelSeed->at(leadingPhotonIndex));

    h_ptpho->Fill(leadingPhotonEt);
    h_etapho->Fill(PhotonEta->at(leadingPhotonIndex));

    if(leadingPhotonIndex>0)nPass[10]++;

    // assign 4-vector of leading photon
    TLorentzVector l4_pho(0,0,0,0);
    l4_pho.SetPtEtaPhiE(
			PhotonEt->at(leadingPhotonIndex),
			PhotonEta->at(leadingPhotonIndex),
			PhotonPhi->at(leadingPhotonIndex),
			PhotonEnergy->at(leadingPhotonIndex)
			);

    // find the reco jet that has the highest matched gen jet pt

    Int_t leadingJetIndex=-1;
    double genJetMaxPt=-9999.;
    
    for(int ijet=0; ijet < patJetPfAk05Pt_->size(); ijet++){

      int genJetIndex = matchedGenJet(ientry, ijet);
      if(genJetIndex < 0)continue;


      double thisGenJetPt = genJetPt_->at(genJetIndex);

      if(thisGenJetPt > genJetMaxPt)
	{
	  genJetMaxPt = thisGenJetPt;
	  leadingJetIndex= ijet;
	}


    } // end of loop over jets

    
    if(leadingPhotonIndex<0 || leadingJetIndex<0)continue;

    nPass[3]++;

    h_ptjet->Fill(patJetPfAk05Pt_->at(leadingJetIndex));
    h_etajet->Fill(patJetPfAk05Eta_->at(leadingJetIndex));


    // counting entries
    if(phoPtBinIndex < 0 || phoPtBinIndex > (nPtBins-1))continue;
    if(phoDecBinIndex < 0)continue;

    for(int ipt=0; ipt < nPtBins; ipt++)nPass[phoPtBinIndex+4]++;

    
    TLorentzVector l4_1stjet(0,0,0,0);
    l4_1stjet.SetPtEtaPhiE(
			   patJetPfAk05Pt_->at(leadingJetIndex),
			   patJetPfAk05Eta_->at(leadingJetIndex),
			   patJetPfAk05Phi_->at(leadingJetIndex),
			   patJetPfAk05En_->at(leadingJetIndex)
			   );
    

    double gj_zgamma = eiko::zgamma(l4_pho, l4_1stjet);
    double gj_dphi   = eiko::deltaPhi(l4_pho, l4_1stjet);
    double gj_cost   = eiko::cosThetaStar(l4_pho, l4_1stjet);
    double gj_chi    = eiko::chiPair(l4_pho, l4_1stjet);
    double gj_pstar  = eiko::pstar(l4_pho, l4_1stjet);
    double gj_yB     = eiko::yB(l4_pho, l4_1stjet);    

    // before applying cuts
    h_zgamma[phoDecBinIndex][phoPtBinIndex][0]->Fill(gj_zgamma);
    h_dphi[phoDecBinIndex][phoPtBinIndex][0]->Fill(gj_dphi);
    h_cost[phoDecBinIndex][phoPtBinIndex][0]->Fill(gj_cost);
    h_chi[phoDecBinIndex][phoPtBinIndex][0]->Fill(gj_chi);
    h_pstar[phoDecBinIndex][phoPtBinIndex][0]->Fill(gj_pstar);
    h_yB[phoDecBinIndex][phoPtBinIndex][0]->Fill(gj_yB);		   

    // after applying cuts
    if(!eiko::separated(l4_pho, l4_1stjet))continue;

    if(!isGoodPho(ientry, leadingPhotonIndex))continue;

    if(!isGoodLooseJet(ientry, leadingJetIndex))continue;

    // after applying cuts
    h_zgamma[phoDecBinIndex][phoPtBinIndex][1]->Fill(gj_zgamma);
    h_dphi[phoDecBinIndex][phoPtBinIndex][1]->Fill(gj_dphi);
    h_cost[phoDecBinIndex][phoPtBinIndex][1]->Fill(gj_cost);
    h_chi[phoDecBinIndex][phoPtBinIndex][1]->Fill(gj_chi);
    h_pstar[phoDecBinIndex][phoPtBinIndex][1]->Fill(gj_pstar);
    h_yB[phoDecBinIndex][phoPtBinIndex][1]->Fill(gj_yB);		   
   
  } // end of loop over entries

  for(int i=0;i<30;i++)if(nPass[i]>0)cout << "nPass[" << i << "] = " << nPass[i] << endl;

  // write out output histogram files to calculate efficiency
  
  std::string remword="/data4/yunju/NewtreeT/MC/";

  size_t pos = _inputDirName.find(remword);
  
  if(pos!= std::string::npos)
    _inputDirName.swap(_inputDirName.erase(pos,remword.length()));
     

  TFile* outFile = new TFile(Form("/home/syu/CVSCode/eff_%s.root",_inputDirName.data()),"recreate");               

  h_pthat->Write();
  h_ptpho->Write();
  h_ptjet->Write();
  h_etapho->Write();
  h_etajet->Write();

  for(int idec=0; idec<4; idec++){

    h_sieie[idec]->Write();
    h_eciso[idec]->Write();
    h_hciso[idec]->Write();
    h_tkiso[idec]->Write();
    h_hovere[idec]->Write();
    h_pixel[idec]->Write();

  }

  for(int idec=0; idec<2; idec++){
    for(int ipt=0; ipt< nPtBins; ipt++){
      for(int ip=0; ip<2; ip++){

	h_zgamma[idec][ipt][ip]->Write();
	h_dphi[idec][ipt][ip]->Write();
	h_cost[idec][ipt][ip]->Write();
	h_chi[idec][ipt][ip]->Write();
	h_pstar[idec][ipt][ip]->Write();
	h_yB[idec][ipt][ip]->Write();				   

      }
    }
  }


  outFile->Close();     

}

Int_t  yj_angularmc_eff::phoDecCode(Long64_t entry, Int_t ipho)
{
  double eta = PhotonScEta->at(ipho);
  if(fabs(eta) < BARREL_MAXETA)return 0;
  if(fabs(eta) > ENDCAP_MINETA && 
     fabs(eta) < ENDCAP_MAXETA)return 1;
  return -1;

}

Bool_t yj_angularmc_eff::isFidPho (Long64_t entry, Int_t ipho)
{
  if(phoDecCode(entry,ipho)<0)return false;
  return true;
}



Bool_t yj_angularmc_eff::isGoodPho(Long64_t entry, Int_t ipho)
{
  double et  = PhotonEt   ->at(ipho);

  if(!isFidPho(entry,ipho))return false;
  if(PhotonhadronicOverEm->at(ipho) > 0.05)return false;
  if(PhotonhasPixelSeed->at(ipho)   > 1e-6)return false; // this should be saved as bool
  if(PhotonecalRecHitSumEtConeDR04->at(ipho) > 4.2 +0.003*et)return false;
  if(PhotonhcalTowerSumEtConeDR04->at(ipho)  > 2.2 +0.001*et)return false;
  if(PhotontrkSumPtHollowConeDR04->at(ipho)  > 2.0 +0.001*et)return false;
 

  return true;
}


Bool_t yj_angularmc_eff::isFidJet (Long64_t entry, Int_t ijet)
{
  if(fabs(patJetPfAk05Eta_->at(ijet)) > 2.4)return false;
  return true;
}


// check if this reco-jet is a good loose jet
Bool_t yj_angularmc_eff::isGoodLooseJet(Long64_t entry, Int_t ijet)
{
  if(!isFidJet(entry,ijet))return false;
  if(patJetPfAk05Pt_->at(ijet) < 30.0)return false;
  if(patJetPfAk05JetNConstituents_->at(ijet) <= 1)return false;
  if(patJetPfAk05NeutHadEFr_->at(ijet) >= 0.99)return false;
  if(patJetPfAk05NeutEmEFr_->at(ijet) >= 0.99)return false;

  //   // for the tracker region
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharHadEFr_->at(ijet) <= 0.0)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharEmEFr_->at(ijet) >= 0.99)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharMulti_->at(ijet) <= 0)return false;

  return true;

}


// check if this reco-jet is a good medium jet
Bool_t yj_angularmc_eff::isGoodMediumJet(Long64_t entry, Int_t ijet)
{
  if(patJetPfAk05Pt_->at(ijet) < 10.0)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet)) > 3.0)return false;
  if(patJetPfAk05JetNConstituents_->at(ijet) <= 1)return false;
  if(patJetPfAk05NeutHadEFr_->at(ijet) >= 0.95)return false;
  if(patJetPfAk05NeutEmEFr_->at(ijet) >= 0.95)return false;

  //   // for the tracker region
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharHadEFr_->at(ijet) <= 0.0)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharEmEFr_->at(ijet) >= 0.99)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharMulti_->at(ijet) <= 0)return false;


  return true;

}

// check if this reco-jet is a good medium jet
Bool_t yj_angularmc_eff::isGoodTightJet(Long64_t entry, Int_t ijet)
{
  if(patJetPfAk05Pt_->at(ijet) < 10.0)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet)) > 3.0)return false;
  if(patJetPfAk05JetNConstituents_->at(ijet) <= 1)return false;
  if(patJetPfAk05NeutHadEFr_->at(ijet) >= 0.90)return false;
  if(patJetPfAk05NeutEmEFr_->at(ijet) >= 0.90)return false;


  //   // for the tracker region
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharHadEFr_->at(ijet) <= 0.0)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharEmEFr_->at(ijet) >= 0.99)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharMulti_->at(ijet) <= 0)return false;

  return true;

}


Int_t  yj_angularmc_eff::matchedGenJet(Long64_t entry, Int_t ijet)
{
  
  int matchedGenIndex = -1;
  for(int k=0; k< genJetPt_->size(); k++)
    {

      double deta = genJetEta_->at(k)-patJetPfAk05Eta_->at(ijet);
      double dphi = genJetPhi_->at(k)-patJetPfAk05Phi_->at(ijet);      
      while (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
      while (dphi <= -TMath::Pi()) dphi += 2*TMath::Pi();

      double dR = sqrt(deta*deta+dphi*dphi);

      double relPt = genJetPt_->at(k)>1e-6? 
	fabs(genJetPt_->at(k)-patJetPfAk05Pt_->at(ijet))/genJetPt_->at(k): -9999.0;

      if(dR<0.4 && relPt < 3.0)
	{
	  matchedGenIndex = k;
	  break;
	}
      
    } // end of loop over generator-level jets
      
  return matchedGenIndex;
}

