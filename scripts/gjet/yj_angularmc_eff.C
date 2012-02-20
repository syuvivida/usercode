#define yj_angularmc_eff_cxx
#include "yj_angularmc_eff.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include "myLib.h"

// const double BARREL_MAXETA=1.4442;
// const double ENDCAP_MINETA=1.566;
const double BARREL_MAXETA=1.44;
const double ENDCAP_MINETA=1.57;
const double ENDCAP_MAXETA=2.5;
const int nDECs = 2;

double ptbound[]={85., 95., 110., 130., 160., 200.};
const int nPtBins = sizeof(ptbound)/sizeof(ptbound[0])-1;

// effective area for correcting pileup, ECAL/HCAL/Tracker from EWK-11-251
const double aeff[2][3]={
  {0.183, 0.062, 0.0167},
  {0.090, 0.180, 0.032}
};




void yj_angularmc_eff::Loop(bool applyCOMCut, bool applyPileUpCorr)
{
  cout << "This is version 3" << endl;
  cout << "There are " << nPtBins << " pt bins" << endl;
  cout << "applyCOMCut: " << applyCOMCut << "\t applyPileUpCorr: " << applyPileUpCorr << endl;

//const double ystarMax = 1.4;  // if photon pt min~25 GeV
  const double ystarMax = 1.0;  // if photon pt min~ 85 GeV
  const double pstarMin = ptbound[0]       *TMath::CosH(ystarMax);
  const double pstarMax = ptbound[nPtBins];
  
  cout << "ystar range: |y*| < " << ystarMax << endl;
  cout << "pstar range: " << pstarMin << " < p* < " << pstarMax << " GeV" << endl;

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  cout << "There are " << nentries << " entries" << endl;


  TH2D* h2_ystar_pstar_template = new TH2D("h2_ystar_pstar_template","",100.0,-5.0,5.0,500,0,500);


  TH1D* h_pt_template = new TH1D("h_pt_template","",500,0,500);
  TH1D* h_pthat       = (TH1D*)h_pt_template->Clone("h_pthat");
  TH1D* h_ptpho       = (TH1D*)h_pt_template->Clone("h_ptpho");
  TH1D* h_ptjet       = (TH1D*)h_pt_template->Clone("h_ptjet");

  TH1D* h_eta_template = new TH1D("h_eta_template","",100,-5.0,5.0);
  TH1D* h_etapho      = (TH1D*)h_eta_template->Clone("h_etapho");
  TH1D* h_etajet      = (TH1D*)h_eta_template->Clone("h_etajet");

  TH1D* h_ptbin_template= new TH1D("h_ptbin_template","", nPtBins, ptbound);
  
  TH1D* h_ptratio_template = new TH1D("h_ptratio_template","",500,0.0,2.5);

  TH1D* h_nvtx_template = new TH1D("h_nvtx_template","",40,0.5,40.5);

  TProfile* pf_nvtx_eciso_template = new TProfile("pf_nvtx_eciso_template","",
						  40,0.5,40.5,-10.0,50.0);

  TProfile* pf_nvtx_hciso_template = new TProfile("pf_nvtx_hciso_template","",
						  40,0.5,40.5,-10.0,50.0);

  TProfile* pf_nvtx_tkiso_template = new TProfile("pf_nvtx_tkiso_template","",
						  40,0.5,40.5,-10.0,50.0);

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

  TH1D* h_ystar_template = new TH1D("h_ystar_template","",100,-5.0,5.0);

  TH1D* h_sieie_template = new TH1D("h_sieie_template","", 100, 0,0.1);
  
  TH1D* h_eciso_template = new TH1D("h_eciso_template","", 240, -10.0, 50.0);

  TH1D* h_hciso_template = new TH1D("h_hciso_template","", 240, -10.0, 50.0);

  TH1D* h_tkiso_template = new TH1D("h_tkiso_template","", 240, -10.0, 50.0);

  TH1D* h_hovere_template = new TH1D("h_hovere_template","", 250, 0,0.5);

  TH1D* h_pixel_template = new TH1D("h_pixel_template","", 3, -0.5,2.5);


  TH2D* h2_ystar_pstar[2]; // in EB and EE

  TH1D* h_jetpt_eff[2];
  TH1D* h_jeteta_eff[2];

  TProfile* pf_nvtx_eciso[2][2];       // in EB and EE, before and after pileupcorr
  TProfile* pf_nvtx_hciso[2][2];       // in EB and EE, before and after pileupcorr
  TProfile* pf_nvtx_tkiso[2][2];       // in EB and EE, before and after pileupcorr

  TH1D* h_pt_eff[nDECs][2];        // in EB and EE
  TH1D* h_nvtx_eff[nDECs][2];        // in EB and EE
  TH1D* h_eta_eff[nPtBins+1][2]; // in various pt bins
  TH1D* h_ptratio[nDECs][2];
  TH1D* h_zgamma[nDECs][nPtBins+1][2];
  TH1D* h_dphi[nDECs][nPtBins+1][2];
  TH1D* h_cost[nDECs][nPtBins+1][2]; // tanH(0.5(y1-y2))
  TH1D* h_cost_COM3D[nDECs][nPtBins+1][2]; // boost to the COM using total momentum
  TH1D* h_cost_COMZ[nDECs][nPtBins+1][2]; // boost to the COM using z-momentum
  TH1D* h_chi[nDECs][nPtBins+1][2];
  TH1D* h_pstar[nDECs][nPtBins+1][2];
  TH1D* h_pstar_COM3D[nDECs][nPtBins+1][2]; // boost to the COM using total momentum
  TH1D* h_pstar_COMZ[nDECs][nPtBins+1][2]; // boost to the COM using z-momentum
  TH1D* h_yB[nDECs][nPtBins+1][2];
  TH1D* h_ystar[nDECs][nPtBins+1][2];
  TH1D* h_ystar_COM3D[nDECs][nPtBins+1][2]; // boost to the COM using total momentum
  TH1D* h_ystar_COMZ[nDECs][nPtBins+1][2]; // boost to the COM using the z-momentum


  // creating histograms
  char* decName[6]={"EB","EE","leadingEB","leadingEE",
		    "leadingEBPileupCorr","leadingEEPileupCorr"};
  char* puName[2]={"raw","pucorr"};

  // debugging histograms
  TH1D* h_sieie[6];
  TH1D* h_eciso[6];
  TH1D* h_hciso[6];
  TH1D* h_tkiso[6];
  TH1D* h_hovere[6];
  TH1D* h_pixel[6];

  for(int idec=0; idec<6; idec++)
    {
      // EB:idec=0, EE:idec=1, leading phtoon only idec=2
      h_sieie[idec] = (TH1D*)h_sieie_template->Clone(Form("h_sieie_%s", decName[idec]));
      h_eciso[idec] = (TH1D*)h_eciso_template->Clone(Form("h_eciso_%s", decName[idec]));
      h_hciso[idec] = (TH1D*)h_hciso_template->Clone(Form("h_hciso_%s", decName[idec]));
      h_tkiso[idec] = (TH1D*)h_tkiso_template->Clone(Form("h_tkiso_%s", decName[idec]));
      h_hovere[idec] = (TH1D*)h_hovere_template->Clone(Form("h_hovere_%s", decName[idec]));
      h_pixel[idec] = (TH1D*)h_pixel_template->Clone(Form("h_pixel_%s", decName[idec]));    
    }


  // histograms for efficiency

  for(int ip=0; ip<2; ip++){

    h_jetpt_eff[ip] = (TH1D*)h_pt_template->Clone(Form("h_jetpt_eff_%d",ip));
    h_jeteta_eff[ip] = (TH1D*)h_eta_template->Clone(Form("h_jeteta_eff_%d",ip));


    for(int idec=0; idec< nDECs; idec++)
      {
	h_pt_eff[idec][ip] = (TH1D*)h_pt_template->Clone(Form("h_pt_eff_%s_%d",
							      decName[idec], ip));
	h_ptratio[idec][ip]= (TH1D*)h_ptratio_template->Clone(Form("h_ptratio_%s_%d",
								   decName[idec],ip));
	h_nvtx_eff[idec][ip]=(TH1D*)h_nvtx_template->Clone(Form("h_nvtx_%s_%d",
								decName[idec],ip));
      }
    for(int ipt=0; ipt < nPtBins+1; ipt++)
      
      if(ipt == (nPtBins))
	h_eta_eff[ipt][ip] = (TH1D*)h_eta_template->Clone(Form("h_eta_eff_allpt_%d",ip));
      else
	h_eta_eff[ipt][ip] = (TH1D*)h_eta_template->Clone(Form("h_eta_eff_%d_%d_%d",
							       (int)ptbound[ipt], (int)ptbound[ipt+1], ip));
  }// end of loop over process, ip=0 before cut, 1 after cut
  
  
  for(int idec=0; idec < nDECs; idec ++){


    h2_ystar_pstar[idec] = (TH2D*)h2_ystar_pstar_template->Clone(Form("h2_ystar_pstar_%s",decName[idec]));

    for(int ip=0; ip<2; ip++){
      pf_nvtx_eciso[idec][ip] = (TProfile*)pf_nvtx_eciso_template->Clone
	(Form("pf_nvtx_eciso_%s_%s",decName[idec],puName[ip]));
      
      pf_nvtx_hciso[idec][ip] = (TProfile*)pf_nvtx_hciso_template->Clone
	(Form("pf_nvtx_hciso_%s_%s",decName[idec],puName[ip]));
      
      pf_nvtx_tkiso[idec][ip] = (TProfile*)pf_nvtx_tkiso_template->Clone
	(Form("pf_nvtx_tkiso_%s_%s",decName[idec],puName[ip]));
    }

    for(int ipt=0; ipt < nPtBins+1; ipt++){
     
      for(int ip=0; ip < 2; ip++){ // 0: before cut, 1: after cut
 
	if(ipt == (nPtBins))
	  {

	    h_zgamma[idec][ipt][ip] = (TH1D*)h_zgamma_template->Clone(Form("h_zgamma_%s_allpt_%d", decName[idec],ip));
      
	    h_dphi[idec][ipt][ip] = (TH1D*)h_dphi_template->Clone(Form("h_dphi_%s_allpt_%d", decName[idec],ip));

	    h_cost[idec][ipt][ip] = (TH1D*)h_cost_template->Clone(Form("h_cost_%s_allpt_%d", decName[idec], ip));

	    h_cost_COM3D[idec][ipt][ip] = (TH1D*)h_cost_template->Clone(Form("h_cost_COM3D_%s_allpt_%d", decName[idec], ip));

	    h_cost_COMZ[idec][ipt][ip] = (TH1D*)h_cost_template->Clone(Form("h_cost_COMZ_%s_allpt_%d", decName[idec], ip));

	    h_chi[idec][ipt][ip] = (TH1D*)h_chi_template->Clone(Form("h_chi_%s_allpt_%d", decName[idec], ip));

	    h_pstar[idec][ipt][ip] = (TH1D*)h_pstar_template->Clone(Form("h_pstar_%s_allpt_%d", decName[idec], ip));

	    h_pstar_COM3D[idec][ipt][ip] = (TH1D*)h_pstar_template->Clone(Form("h_pstar_COM3D_%s_allpt_%d", decName[idec], ip));

	    h_pstar_COMZ[idec][ipt][ip] = (TH1D*)h_pstar_template->Clone(Form("h_pstar_COMZ_%s_allpt_%d", decName[idec], ip));

	    h_yB[idec][ipt][ip] = (TH1D*)h_yB_template->Clone(Form("h_yB_%s_allpt_%d", decName[idec], ip));

	    h_ystar[idec][ipt][ip] = (TH1D*)h_ystar_template->Clone(Form("h_ystar_%s_allpt_%d", decName[idec], ip));

	    h_ystar_COM3D[idec][ipt][ip] = (TH1D*)h_ystar_template->Clone(Form("h_ystar_COM3D_%s_allpt_%d", decName[idec], ip));

	    h_ystar_COMZ[idec][ipt][ip] = (TH1D*)h_ystar_template->Clone(Form("h_ystar_COMZ_%s_allpt_%d", decName[idec], ip));

	  }

	else
	  {

	    h_zgamma[idec][ipt][ip] = (TH1D*)h_zgamma_template->Clone(Form("h_zgamma_%s_%d_%d_%d", decName[idec],
									   (int)ptbound[ipt], (int)ptbound[ipt+1], ip));
      
	    h_dphi[idec][ipt][ip] = (TH1D*)h_dphi_template->Clone(Form("h_dphi_%s_%d_%d_%d", decName[idec],
								       (int)ptbound[ipt], (int)ptbound[ipt+1], ip));

	    h_cost[idec][ipt][ip] = (TH1D*)h_cost_template->Clone(Form("h_cost_%s_%d_%d_%d", decName[idec],
								       (int)ptbound[ipt], (int)ptbound[ipt+1], ip));

	    h_cost_COM3D[idec][ipt][ip] = (TH1D*)h_cost_template->Clone(Form("h_cost_COM3D_%s_%d_%d_%d", decName[idec],
									     (int)ptbound[ipt], (int)ptbound[ipt+1], ip));

	    h_cost_COMZ[idec][ipt][ip] = (TH1D*)h_cost_template->Clone(Form("h_cost_COMZ_%s_%d_%d_%d", decName[idec],
									     (int)ptbound[ipt], (int)ptbound[ipt+1], ip));

	    h_chi[idec][ipt][ip] = (TH1D*)h_chi_template->Clone(Form("h_chi_%s_%d_%d_%d", decName[idec],
								     (int)ptbound[ipt], (int)ptbound[ipt+1], ip));

	    h_pstar[idec][ipt][ip] = (TH1D*)h_pstar_template->Clone(Form("h_pstar_%s_%d_%d_%d", decName[idec],
									 (int)ptbound[ipt], (int)ptbound[ipt+1], ip));

	    h_pstar_COM3D[idec][ipt][ip] = (TH1D*)h_pstar_template->Clone(Form("h_pstar_COM3D_%s_%d_%d_%d", decName[idec],
									     (int)ptbound[ipt], (int)ptbound[ipt+1], ip));

	    h_pstar_COMZ[idec][ipt][ip] = (TH1D*)h_pstar_template->Clone(Form("h_pstar_COMZ_%s_%d_%d_%d", decName[idec],
									     (int)ptbound[ipt], (int)ptbound[ipt+1], ip));

	    h_yB[idec][ipt][ip] = (TH1D*)h_yB_template->Clone(Form("h_yB_%s_%d_%d_%d", decName[idec],
								   (int)ptbound[ipt], (int)ptbound[ipt+1], ip));

	    h_ystar[idec][ipt][ip] = (TH1D*)h_ystar_template->Clone(Form("h_ystar_%s_%d_%d_%d", decName[idec],
								       (int)ptbound[ipt], (int)ptbound[ipt+1], ip));

	    h_ystar_COM3D[idec][ipt][ip] = (TH1D*)h_ystar_template->Clone(Form("h_ystar_COM3D_%s_%d_%d_%d", decName[idec],
									     (int)ptbound[ipt], (int)ptbound[ipt+1], ip));

	    h_ystar_COMZ[idec][ipt][ip] = (TH1D*)h_ystar_template->Clone(Form("h_ystar_COMZ_%s_%d_%d_%d", decName[idec],
									     (int)ptbound[ipt], (int)ptbound[ipt+1], ip));

	  }
	
      } // end of loop over processes
   
    } // end of loop over pt bins

  } // loop over barrel and endcaps

  //---------------------------------------------------------------------------------------------------------------------
  //
  //   Now we start loop over events
  // 
  //---------------------------------------------------------------------------------------------------------------------

  Long64_t nPass[30]={0};

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
//     if(jentry > 50000 ) break;
    nPass[0]++;
    h_pthat->Fill(PhotonMCpthat); // make sure we are checking the right MC samples


    //vertex selection
    if(EvtInfo_nVtxGood<1) continue;
    nPass[1]++;

    int leadingPhotonIndex = -1;
    double phoMaxPt = -9999.;

    //-------------------------------------------------------------------------
    // now find a good leading photon that is matched to hard scattering photon
    //-------------------------------------------------------------------------
    for(unsigned int ipho=0; ipho < PhotonPt->size(); ipho++)
      {	

	int decIndex = phoDecCode(ientry, ipho);
	if(!isFidPho(ientry,ipho))continue;

	h_sieie[decIndex] -> Fill(PhotonSigmaIetaIeta->at(ipho));
	h_eciso[decIndex] -> Fill(phoEcalIso(ientry,ipho,false));
	h_hciso[decIndex] -> Fill(phoHcalIso(ientry,ipho,false));
	h_tkiso[decIndex] -> Fill(phoTrkIso(ientry,ipho,false));
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

    double leadingPhotonEta = PhotonEta->at(leadingPhotonIndex);
    int phoDecBinIndex = phoDecCode(ientry, leadingPhotonIndex);

    h_sieie[phoDecBinIndex+2] -> Fill(PhotonSigmaIetaIeta->at(leadingPhotonIndex));
    h_hovere[phoDecBinIndex+2]-> Fill(PhotonhadronicOverEm->at(leadingPhotonIndex));
    h_pixel[phoDecBinIndex+2] -> Fill(PhotonhasPixelSeed->at(leadingPhotonIndex));

    h_eciso[phoDecBinIndex+2] -> Fill(phoEcalIso(ientry,leadingPhotonIndex,false));
    h_hciso[phoDecBinIndex+2] -> Fill(phoHcalIso(ientry,leadingPhotonIndex,false));
    h_tkiso[phoDecBinIndex+2] -> Fill(phoTrkIso(ientry,leadingPhotonIndex,false));


   // after pileup correction
    h_eciso[phoDecBinIndex+4] -> Fill(phoEcalIso(ientry,leadingPhotonIndex,true));
    h_hciso[phoDecBinIndex+4] -> Fill(phoHcalIso(ientry,leadingPhotonIndex,true));
    h_tkiso[phoDecBinIndex+4] -> Fill(phoTrkIso(ientry,leadingPhotonIndex,true));


    h_ptpho->Fill(leadingPhotonEt);
    h_etapho->Fill(leadingPhotonEta);

    if(leadingPhotonIndex>0)nPass[10]++;

    // assign 4-vector of leading photon
    TLorentzVector l4_pho(0,0,0,0);
    l4_pho.SetPtEtaPhiE(
			PhotonEt->at(leadingPhotonIndex),
			PhotonEta->at(leadingPhotonIndex),
			PhotonPhi->at(leadingPhotonIndex),
			PhotonEnergy->at(leadingPhotonIndex)
			);

    //-------------------------------------------------------------------------
    // find the reco jet that has the highest matched gen jet pt
    //-------------------------------------------------------------------------

    Int_t leadingJetIndex=-1;
    double genJetMaxPt=-9999.;
    
    for(int ijet=0; ijet < patJetPfAk05Pt_->size(); ijet++){

      if(!isFidJet(ientry,ijet))continue;
      int genJetIndex = matchedGenJet(ientry, ijet);
      if(genJetIndex < 0)continue;


      double thisGenJetPt = genJetPt_->at(genJetIndex);

      if(thisGenJetPt > genJetMaxPt)
	{
	  genJetMaxPt = thisGenJetPt;
	  leadingJetIndex= ijet;
	}


    } // end of loop over jets


    // need to have a leading photon and a leading jet
    if(leadingPhotonIndex<0 || leadingJetIndex<0)continue;

    double leadingJetPt  = patJetPfAk05Pt_->at(leadingJetIndex);
    double leadingJetEta = patJetPfAk05Eta_->at(leadingJetIndex);

    h_ptjet->Fill(leadingJetPt);
    h_etajet->Fill(leadingJetEta);

    nPass[3]++;


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



    //-------------------------------------------------------------------------
    // Preselections                                    
    //-------------------------------------------------------------------------
    
    if(!eiko::separated(l4_pho, l4_1stjet))continue;

    double gj_pstar_comZ = eiko::pstar_ZBoostToCM(l4_pho, l4_1stjet);
    double gj_ystar_comZ = eiko::ystar_ZBoostToCM(l4_pho, l4_1stjet);

    bool passCOMCut = gj_pstar_comZ > pstarMin && gj_pstar_comZ < pstarMax
      && fabs(gj_ystar_comZ) < ystarMax;

    if(applyCOMCut && !passCOMCut)continue;


    h2_ystar_pstar[phoDecBinIndex]->Fill(gj_ystar_comZ, gj_pstar_comZ);


    //-------------------------------------------------------------------------
    // Now we plot distributions before and after cuts
    //-------------------------------------------------------------------------

    double gj_zgamma = eiko::zgamma(l4_pho, l4_1stjet);
    double gj_dphi   = eiko::deltaPhi(l4_pho, l4_1stjet);

    double gj_cost   = eiko::cosThetaStar(l4_pho, l4_1stjet);
    double gj_cost_com3D = eiko::cosThetaStar_BoostToCM(l4_pho, l4_1stjet);
    double gj_cost_comZ = eiko::cosThetaStar_ZBoostToCM(l4_pho, l4_1stjet);

    double gj_chi    = eiko::chiPair(l4_pho, l4_1stjet);

    double gj_pstar  = eiko::pstar(l4_pho, l4_1stjet);
    double gj_pstar_com3D = eiko::pstar_BoostToCM(l4_pho, l4_1stjet);

    double gj_yB     = eiko::yB(l4_pho, l4_1stjet);    

    double gj_ystar   = eiko::ystar(l4_pho, l4_1stjet);
    double gj_ystar_com3D = eiko::ystar_BoostToCM(l4_pho, l4_1stjet);

    double ptRatio   = leadingJetPt/leadingPhotonEt;
    // before applying cuts

    h_nvtx_eff[phoDecBinIndex][0]->Fill(EvtInfo_nVtxGood);

    h_pt_eff[phoDecBinIndex][0]->Fill(leadingPhotonEt);
    h_eta_eff[phoPtBinIndex][0]->Fill(leadingPhotonEta);

    h_jetpt_eff[0]->Fill(leadingJetPt);
    h_jeteta_eff[0]->Fill(leadingJetEta);

    h_ptratio[phoDecBinIndex][0]->Fill(ptRatio);

    h_zgamma[phoDecBinIndex][phoPtBinIndex][0]->Fill(gj_zgamma);
    h_dphi[phoDecBinIndex][phoPtBinIndex][0]->Fill(gj_dphi);

    h_cost[phoDecBinIndex][phoPtBinIndex][0]->Fill(gj_cost);
    h_cost_COM3D[phoDecBinIndex][phoPtBinIndex][0]->Fill(gj_cost_com3D);
    h_cost_COMZ[phoDecBinIndex][phoPtBinIndex][0]->Fill(gj_cost_comZ);

    h_chi[phoDecBinIndex][phoPtBinIndex][0]->Fill(gj_chi);

    h_pstar[phoDecBinIndex][phoPtBinIndex][0]->Fill(gj_pstar);
    h_pstar_COM3D[phoDecBinIndex][phoPtBinIndex][0]->Fill(gj_pstar_com3D);
    h_pstar_COMZ[phoDecBinIndex][phoPtBinIndex][0]->Fill(gj_pstar_comZ);

    h_yB[phoDecBinIndex][phoPtBinIndex][0]->Fill(gj_yB);		   
    h_ystar[phoDecBinIndex][phoPtBinIndex][0]->Fill(gj_ystar);		   
    h_ystar_COM3D[phoDecBinIndex][phoPtBinIndex][0]->Fill(gj_ystar_com3D);
    h_ystar_COMZ[phoDecBinIndex][phoPtBinIndex][0]->Fill(gj_ystar_comZ);


    // lump all pt together
    h_eta_eff[nPtBins][0]->Fill(leadingPhotonEta);
    h_zgamma[phoDecBinIndex][nPtBins][0]->Fill(gj_zgamma);
    h_dphi[phoDecBinIndex][nPtBins][0]->Fill(gj_dphi);

    h_cost[phoDecBinIndex][nPtBins][0]->Fill(gj_cost);
    h_cost_COM3D[phoDecBinIndex][nPtBins][0]->Fill(gj_cost_com3D);
    h_cost_COMZ[phoDecBinIndex][nPtBins][0]->Fill(gj_cost_comZ);

    h_chi[phoDecBinIndex][nPtBins][0]->Fill(gj_chi);

    h_pstar[phoDecBinIndex][nPtBins][0]->Fill(gj_pstar);
    h_pstar_COM3D[phoDecBinIndex][nPtBins][0]->Fill(gj_pstar_com3D);
    h_pstar_COMZ[phoDecBinIndex][nPtBins][0]->Fill(gj_pstar_comZ);

    h_yB[phoDecBinIndex][nPtBins][0]->Fill(gj_yB);		   
    h_ystar[phoDecBinIndex][nPtBins][0]->Fill(gj_ystar);		   
    h_ystar_COM3D[phoDecBinIndex][nPtBins][0]->Fill(gj_ystar_com3D);
    h_ystar_COMZ[phoDecBinIndex][nPtBins][0]->Fill(gj_ystar_comZ);


    // checking the isolation behaviour as a function of nvertex
    pf_nvtx_eciso[phoDecBinIndex][0]->Fill(EvtInfo_nVtxGood,
					   phoEcalIso(ientry,leadingPhotonIndex,false));
    pf_nvtx_hciso[phoDecBinIndex][0]->Fill(EvtInfo_nVtxGood,
					   phoHcalIso(ientry,leadingPhotonIndex,false));
    pf_nvtx_tkiso[phoDecBinIndex][0]->Fill(EvtInfo_nVtxGood,
					   phoTrkIso(ientry,leadingPhotonIndex,false));

    pf_nvtx_eciso[phoDecBinIndex][1]->Fill(EvtInfo_nVtxGood,
					   phoEcalIso(ientry,leadingPhotonIndex,true));
    pf_nvtx_hciso[phoDecBinIndex][1]->Fill(EvtInfo_nVtxGood,
					   phoHcalIso(ientry,leadingPhotonIndex,true));
    pf_nvtx_tkiso[phoDecBinIndex][1]->Fill(EvtInfo_nVtxGood,
					   phoTrkIso(ientry,leadingPhotonIndex,true));

    // after applying cuts

    //=============================================================
    if(!eiko::separated(l4_pho, l4_1stjet))continue;
    if(!isGoodPho(ientry, leadingPhotonIndex, applyPileUpCorr))continue;
    if(!isGoodLooseJet(ientry, leadingJetIndex))continue;
    //=============================================================

    h_nvtx_eff[phoDecBinIndex][1]->Fill(EvtInfo_nVtxGood);

    h_pt_eff[phoDecBinIndex][1]->Fill(leadingPhotonEt);
    h_eta_eff[phoPtBinIndex][1]->Fill(leadingPhotonEta);

    h_jetpt_eff[1]->Fill(leadingJetPt);
    h_jeteta_eff[1]->Fill(leadingJetEta);

    h_ptratio[phoDecBinIndex][1]->Fill(ptRatio);

    h_zgamma[phoDecBinIndex][phoPtBinIndex][1]->Fill(gj_zgamma);
    h_dphi[phoDecBinIndex][phoPtBinIndex][1]->Fill(gj_dphi);

    h_cost[phoDecBinIndex][phoPtBinIndex][1]->Fill(gj_cost);
    h_cost_COM3D[phoDecBinIndex][phoPtBinIndex][1]->Fill(gj_cost_com3D);
    h_cost_COMZ[phoDecBinIndex][phoPtBinIndex][1]->Fill(gj_cost_comZ);

    h_chi[phoDecBinIndex][phoPtBinIndex][1]->Fill(gj_chi);

    h_pstar[phoDecBinIndex][phoPtBinIndex][1]->Fill(gj_pstar);
    h_pstar_COM3D[phoDecBinIndex][phoPtBinIndex][1]->Fill(gj_pstar_com3D);
    h_pstar_COMZ[phoDecBinIndex][phoPtBinIndex][1]->Fill(gj_pstar_comZ);

    h_yB[phoDecBinIndex][phoPtBinIndex][1]->Fill(gj_yB);		   
    h_ystar[phoDecBinIndex][phoPtBinIndex][1]->Fill(gj_ystar);		   
    h_ystar_COM3D[phoDecBinIndex][phoPtBinIndex][1]->Fill(gj_ystar_com3D);
    h_ystar_COMZ[phoDecBinIndex][phoPtBinIndex][1]->Fill(gj_ystar_comZ);


    // lump all pt together
    h_eta_eff[nPtBins][1]->Fill(leadingPhotonEta);
    h_zgamma[phoDecBinIndex][nPtBins][1]->Fill(gj_zgamma);
    h_dphi[phoDecBinIndex][nPtBins][1]->Fill(gj_dphi);

    h_cost[phoDecBinIndex][nPtBins][1]->Fill(gj_cost);
    h_cost_COM3D[phoDecBinIndex][nPtBins][1]->Fill(gj_cost_com3D);
    h_cost_COMZ[phoDecBinIndex][nPtBins][1]->Fill(gj_cost_comZ);

    h_chi[phoDecBinIndex][nPtBins][1]->Fill(gj_chi);

    h_pstar[phoDecBinIndex][nPtBins][1]->Fill(gj_pstar);
    h_pstar_COM3D[phoDecBinIndex][nPtBins][1]->Fill(gj_pstar_com3D);
    h_pstar_COMZ[phoDecBinIndex][nPtBins][1]->Fill(gj_pstar_comZ);

    h_yB[phoDecBinIndex][nPtBins][1]->Fill(gj_yB);		   
    h_ystar[phoDecBinIndex][nPtBins][1]->Fill(gj_ystar);		   
    h_ystar_COM3D[phoDecBinIndex][nPtBins][1]->Fill(gj_ystar_com3D);
    h_ystar_COMZ[phoDecBinIndex][nPtBins][1]->Fill(gj_ystar_comZ);



  } // end of loop over entries

  for(int i=0;i<30;i++)if(nPass[i]>0)cout << "nPass[" << i << "] = " << nPass[i] << endl;


  //---------------------------------------------------------------------------------------
  // Write these histograms to an output root file, the output will be used for efficiency
  //---------------------------------------------------------------------------------------
     
  
  std::string remword="/data4/yunju/NewtreeT/MC/";

  size_t pos = _inputDirName.find(remword);

  if(pos!= std::string::npos)
    _inputDirName.swap(_inputDirName.erase(pos,remword.length()));

  if(applyCOMCut)_inputDirName += "_withCOMCut";
  if(applyPileUpCorr)_inputDirName += "_pileCorr";
  TFile* outFile = new TFile(Form("/home/syu/CVSCode/eff_%s.root",_inputDirName.data()),"recreate");               
  h_pthat->Write();
  h_ptpho->Write();
  h_ptjet->Write();
  h_etapho->Write();
  h_etajet->Write();

  for(int idec=0; idec<6; idec++){

    h_sieie[idec]->Write();
    h_eciso[idec]->Write();
    h_hciso[idec]->Write();
    h_tkiso[idec]->Write();
    h_hovere[idec]->Write();
    h_pixel[idec]->Write();

  }


  for(int ip=0; ip<2; ip++){

    h_jetpt_eff[ip]->Write();
    h_jeteta_eff[ip]->Write();

    for(int idec=0; idec<nDECs; idec++) h_pt_eff[idec][ip]->Write();
    for(int idec=0; idec<nDECs; idec++) h_ptratio[idec][ip]->Write();
    for(int idec=0; idec<nDECs; idec++) h_nvtx_eff[idec][ip]->Write();
    for(int ipt=0; ipt < nPtBins+1; ipt++) h_eta_eff[ipt][ip]->Write();
  }// end of loop over process, ip=0 before cut, 1 after cut

  for(int idec=0; idec< nDECs; idec++){

    h2_ystar_pstar[idec]->Write();

    for(int ip=0; ip<2; ip++){
      pf_nvtx_eciso[idec][ip]->Write();
      pf_nvtx_hciso[idec][ip]->Write();
      pf_nvtx_tkiso[idec][ip]->Write();
    }


    for(int ipt=0; ipt< nPtBins+1; ipt++){
      for(int ip=0; ip<2; ip++){
	
	h_zgamma[idec][ipt][ip]->Write();
	h_dphi[idec][ipt][ip]->Write();
	h_cost[idec][ipt][ip]->Write();
	h_cost_COM3D[idec][ipt][ip]->Write();
	h_cost_COMZ[idec][ipt][ip]->Write();
	h_chi[idec][ipt][ip]->Write();
	h_pstar[idec][ipt][ip]->Write();
	h_pstar_COM3D[idec][ipt][ip]->Write();
	h_pstar_COMZ[idec][ipt][ip]->Write();
	h_yB[idec][ipt][ip]->Write();				   
	h_ystar[idec][ipt][ip]->Write();
	h_ystar_COM3D[idec][ipt][ip]->Write();
	h_ystar_COMZ[idec][ipt][ip]->Write();

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
  double et  = PhotonEt   ->at(ipho);  
  if(et < ptbound[0])return false;
  if(et > ptbound[nPtBins])return false;
  return true;
}



Bool_t yj_angularmc_eff::isGoodPho(Long64_t entry, Int_t ipho, bool applyPileUpCorr)
{
  double et  = PhotonEt   ->at(ipho);

  if(!isFidPho(entry,ipho))return false;
  if(PhotonhadronicOverEm->at(ipho) > 0.05)return false;
  if(PhotonhasPixelSeed->at(ipho)   > 1e-6)return false; // this should be saved as bool

  double eciso = phoEcalIso(entry, ipho, applyPileUpCorr); 
  double hciso = phoHcalIso(entry, ipho, applyPileUpCorr);
  double tkiso = phoTrkIso(entry, ipho, applyPileUpCorr);

  if(eciso > 4.2 +0.006*et)return false;
  if(hciso > 2.2 +0.0025*et)return false;
  if(tkiso > 2.0 +0.001*et)return false;

  return true;
}

Double_t  yj_angularmc_eff::phoEcalIso(Long64_t entry, Int_t ipho, bool applyPileUpCorr)
{
  if(!isFidPho(entry,ipho))return -9999.0;

  double eciso = PhotonecalRecHitSumEtConeDR04->at(ipho);
  int decBinIndex = phoDecCode(entry,ipho);

  if(applyPileUpCorr)
    eciso -= aeff[decBinIndex][0]*Photonrho25; 
  return eciso;
}
 
Double_t yj_angularmc_eff::phoHcalIso(Long64_t entry, Int_t ipho, bool applyPileUpCorr)
{
  if(!isFidPho(entry,ipho))return -9999.0;

  double hciso = PhotonhcalTowerSumEtConeDR04->at(ipho);
  int decBinIndex = phoDecCode(entry,ipho);

  if(applyPileUpCorr)
    hciso -= aeff[decBinIndex][1]*Photonrho25; 
  return hciso;

}

Double_t yj_angularmc_eff::phoTrkIso(Long64_t entry, Int_t ipho, bool applyPileUpCorr)
{
  if(!isFidPho(entry,ipho))return -9999.0;

  double tkiso = PhotontrkSumPtHollowConeDR04->at(ipho);
  int decBinIndex = phoDecCode(entry,ipho);

  if(applyPileUpCorr)
    tkiso -= aeff[decBinIndex][2]*Photonrho25; 
  return tkiso;

}


Bool_t yj_angularmc_eff::isFidJet (Long64_t entry, Int_t ijet)
{
  if(fabs(patJetPfAk05Eta_->at(ijet)) > 2.4)return false;
  if(patJetPfAk05Pt_->at(ijet) < 30.0)return false;
  return true;
}


// check if this reco-jet is a good loose jet
Bool_t yj_angularmc_eff::isGoodLooseJet(Long64_t entry, Int_t ijet)
{
  if(!isFidJet(entry,ijet))return false;
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

