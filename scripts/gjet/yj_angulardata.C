#define yj_angulardata_cxx
#include "yj_angulardata.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

const double BARREL_MAXETA=1.4442;
const double ENDCAP_MINETA=1.566;
const double ENDCAP_MAXETA=2.5;

void yj_angulardata::Loop()
{
  cout << "This is version 0" << endl;
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  cout << "There are " << nentries << " entries" << endl;

  Long64_t nPass[30]={0};

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    

    nPass[0]++;
    

    // first count only the entries

    if(EvtInfo_nVtxGood<1) continue;
    //vertex selection
    nPass[1]++;

    bool fireHLT_Photon75_CaloIdVL = false;
    for(int itrig=0; itrig < trigName->size(); itrig++)
      {
	if(trigName->at(itrig).find("HLT_Photon75_CaloIdVL_v")!= string::npos)
	  {
	    fireHLT_Photon75_CaloIdVL = true;
	    break;
	  }
      }

    //
    if(!fireHLT_Photon75_CaloIdVL)continue;
      
    nPass[2]++;


    int leadingPhotonIndex = -1;
    double phoMaxPt = -9999.;
     
    // now find a good leading photon
    for(unsigned int ipho=0; ipho < PhotonPt->size(); ipho++)
      {	
	if(!isGoodPho(ientry,ipho))continue;

	double thisPhoPt= PhotonEt->at(ipho);      
	if(thisPhoPt > phoMaxPt)
	  {
	    phoMaxPt = thisPhoPt;
	    leadingPhotonIndex= ipho;
	  }

	
      } // end of leading photon search

    if(leadingPhotonIndex<0)continue;

    nPass[3]++;

    TLorentzVector l4_pho(0,0,0,0);
    l4_pho.SetPtEtaPhiE(
			PhotonEt->at(leadingPhotonIndex),
			PhotonEta->at(leadingPhotonIndex),
			PhotonPhi->at(leadingPhotonIndex),
			PhotonEnergy->at(leadingPhotonIndex)
			);



    // first check which reco jet is the one from the highest and 
    // second et gen jet
    Int_t leadingJetIndex=-1;
    double jetMaxPt=-9999.;

    Int_t secondLeadingJetIndex=-1;
    double secondJetMaxPt=-9999.;

    

    for(int ijet=0; ijet < patJetPfAk05Pt_->size(); ijet++){

      TLorentzVector l4_tempjet(0,0,0,0);
      l4_tempjet.SetPtEtaPhiE(
			      patJetPfAk05Pt_->at(ijet),
			      patJetPfAk05Eta_->at(ijet),
			      patJetPfAk05Phi_->at(ijet),
			      patJetPfAk05En_->at(ijet)
			      );

      double dR = l4_tempjet.DeltaR(l4_pho);

      if(dR<0.5)continue;
      
      if(!isGoodLooseJet(ientry,ijet))continue;

      double thisJetPt = patJetPfAk05Pt_->at(ijet);

      if(thisJetPt > jetMaxPt)
	{
	  secondJetMaxPt = jetMaxPt;
	  secondLeadingJetIndex = leadingJetIndex;
	  jetMaxPt = thisJetPt;
	  leadingJetIndex= ijet;
	}
      else if(thisJetPt > secondJetMaxPt)
	{
	  secondJetMaxPt = thisJetPt;
	  secondLeadingJetIndex = ijet;
	}


    } // end of loop over jets

    
    
    if(leadingPhotonIndex<0 || leadingJetIndex<0)continue;
    nPass[4]++;
    
    double leadingPhotonEt = PhotonEt->at(leadingPhotonIndex);

    
    if(leadingPhotonEt >= 85 && leadingPhotonEt < 95)nPass[5]++;
    
    if(leadingPhotonEt >= 95 && leadingPhotonEt < 110)nPass[6]++;
    
    if(leadingPhotonEt >= 110 && leadingPhotonEt < 130)nPass[7]++;
    
    if(leadingPhotonEt >= 130 && leadingPhotonEt < 160)nPass[8]++;
    
    if(leadingPhotonEt >= 160 && leadingPhotonEt < 200)nPass[9]++;

  } // end of loop over entries


  for(int i=0;i<30;i++)if(nPass[i]>0)cout << "nPass[" << i << "] = " << nPass[i] << endl;


}


Bool_t yj_angulardata::isGoodPho(Long64_t entry, Int_t ipho)
{
  Bool_t isEB=false;
  Bool_t isEE=false;
  double eta = PhotonScEta->at(ipho);
  double et  = PhotonEt   ->at(ipho);

  if(fabs(eta) < BARREL_MAXETA)isEB=true;
  if(fabs(eta) > ENDCAP_MINETA && 
     fabs(eta) < ENDCAP_MAXETA) isEE=true;

  //   if(!isEB && !isEE)return false;
  if(fabs(eta) > 1.44)return false;

  if(PhotonhadronicOverEm->at(ipho) > 0.05)return false;
  if(PhotonhasPixelSeed->at(ipho)   > 1e-6)return false; // this should be saved as bool
  if(PhotonecalRecHitSumEtConeDR04->at(ipho) > 4.2 +0.003*et)return false;
  if(PhotonhcalTowerSumEtConeDR04->at(ipho)  > 2.2 +0.001*et)return false;
  if(PhotontrkSumPtHollowConeDR04->at(ipho)  > 2.0 +0.001*et)return false;
 

  return true;
}


Bool_t yj_angulardata::isFidJet (Long64_t entry, Int_t ijet)
{
  if(patJetPfAk05Pt_->at(ijet) < 10.0)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet)) > 3.0)return false;
  return true;
}


// check if this reco-jet is a good loose jet
Bool_t yj_angulardata::isGoodLooseJet(Long64_t entry, Int_t ijet)
{
  if(patJetPfAk05Pt_->at(ijet) < 30.0)return false;
  //  if(fabs(patJetPfAk05Eta_->at(ijet)) > 3.0)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet)) > 2.4)return false;
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
Bool_t yj_angulardata::isGoodMediumJet(Long64_t entry, Int_t ijet)
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
Bool_t yj_angulardata::isGoodTightJet(Long64_t entry, Int_t ijet)
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

