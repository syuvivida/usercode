#define gg_angularmc_eff_cxx
#include "gg_angularmc_eff.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include "myLib.h"

const double BARREL_MAXETA=1.44;
const double ENDCAP_MINETA=1.57;
const double ENDCAP_MAXETA=2.5;
const int nDECs = 2;

double ptbound[]={85., 95., 110., 130., 160., 200.};
// double ptbound[]={40., 45., 50., 55., 60., 65., 70., 75., 85., 95., 110., 130., 160., 200.};
const int nPtBins = sizeof(ptbound)/sizeof(ptbound[0])-1;

// effective area for correcting pileup, ECAL/HCAL/Tracker from EWK-11-251
const double aeff[2][3]={
  {0.183, 0.062, 0.0167},
  {0.090, 0.180, 0.032}
};


void gg_angularmc_eff::Loop(bool applyCOMCut, bool applyPileUpCorr)
{
  cout << "This is version 4" << endl;
  cout << "There are " << nPtBins << " pt bins" << endl;
  cout << "applyCOMCut: " << applyCOMCut << "\t applyPileUpCorr: " << applyPileUpCorr << endl;

//const double ystarMax = 1.4;  // if photon pt min~25 GeV
  const double ystarMax = 1.0;  // if photon pt min~ 85 GeV
  const double pstarMin = ptbound[0]       *TMath::CosH(ystarMax);
  const double pstarMax = ptbound[nPtBins];
  
  if(applyCOMCut){
    cout << "ystar range: |y*| < " << ystarMax << endl;
    cout << "pstar range: " << pstarMin << " < p* < " << pstarMax << " GeV" << endl;
  }
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  cout << "There are " << nentries << " entries" << endl;

  // defining templates of various histograms

  TH2D* h2_ystar_pstar_template = new TH2D("h2_ystar_pstar_template","",100.0,-5.0,5.0,500,0,500);


  TH1D* h_pt_template = new TH1D("h_pt_template","",500,0,500);

  TH1I* h_njet_template = new TH1I("h_njet_template","",50,-0.5,49.5);

  TH1D* h_eta_template = new TH1D("h_eta_template","",100,-5.0,5.0);

  TH1D* h_ptbin_template= new TH1D("h_ptbin_template","", nPtBins, ptbound);

  TH1D* h_ptratio_template = new TH1D("h_ptratio_template","",600,0.0,3.0);

  TH1D* h_dR_template = new TH1D("h_dR_template","",200,0.0,2.0);

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

  // now declare the necessary histograms

  TH1D* h_sieie[6];
  TH1D* h_eciso[6];
  TH1D* h_hciso[6];
  TH1D* h_tkiso[6];
  TH1D* h_hovere[6];
  TH1D* h_pixel[6];

  TH1D* h_recgenpt_jet;
  TH1D* h_dR_jet;
  TH1D* h_recgenpt_pho;
  TH1D* h_dR_pho;
  TH1D* h_dR_gj;


  TH2D* h2_ystar_pstar[2]; // in EB and EE

  TH1D* h_jetpt_eff[2];
  TH1D* h_jeteta_eff[2];

  TH1I* h_njet[2][2];                  // in EB and EE, before and after ID selections
  TH1I* h_njetraw[2];                  // jet multiplicity before perselections in barrel and endcap
  TH1D* h_ptclosestjet[2];             // pt of the closest jet in deltaR to the leading photon in barrel and endcap
  TH1D* h_etaclosestjet[2];            // eta of the closest jet in deltaR to the leading photon in barrel and endcap
  TH1D* h_dRclosestjet[2];             // deltaR of the closest jet in deltaR to the leading photon in barrel and endcap


  TProfile* pf_nvtx_eciso[2][2];       // in EB and EE, before and after pileupcorr
  TProfile* pf_nvtx_hciso[2][2];       // in EB and EE, before and after pileupcorr
  TProfile* pf_nvtx_tkiso[2][2];       // in EB and EE, before and after pileupcorr

  TH1D* h_pt_eff[nDECs][2];        // in EB and EE
  TH1D* h_eta_eff[nPtBins+1][2]; // in various pt bins
  TH1D* h_nvtx_eff[nDECs][2];        // in EB and EE
  TH1D* h_ptratio[nDECs][2];

  // angular and COM distributions
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


  // histograms for photon and the vector sum of all good jets
  TH1D* h_alljetpt_eff[nDECs][2];   // in EB and EE before and after the selections
  TH1D* h_alljeteta_eff[nDECs][2];
  TH1D* h_cost_COMZ_alljet[nDECs][2]; 
  TH1D* h_pstar_COMZ_alljet[nDECs][2];
  TH1D* h_ystar_COMZ_alljet[nDECs][2];


  // creating histograms
  char* decName[6]={"EB","EE","leadingEB","leadingEE",
		    "leadingEBPileupCorr","leadingEEPileupCorr"};
  char* puName[2]={"raw","pucorr"};

  // debugging histograms
  TH1D* h_pthat       = (TH1D*)h_pt_template->Clone("h_pthat");

  TH1D* h_ptpho       = (TH1D*)h_pt_template->Clone("h_ptpho");
  TH1D* h_ptjet       = (TH1D*)h_pt_template->Clone("h_ptjet");
  TH1D* h_etapho      = (TH1D*)h_eta_template->Clone("h_etapho");
  TH1D* h_etajet      = (TH1D*)h_eta_template->Clone("h_etajet");

  TH1I* h_jetparton   = new TH1I("h_jetparton","",150,-99.5,50.5);
  TH1I* h_leadingjetparton   = new TH1I("h_leadingjetparton","",150,-99.5,50.5);

  h_recgenpt_jet = (TH1D*)h_ptratio_template->Clone("h_recgenpt_jet");
  h_recgenpt_pho = (TH1D*)h_ptratio_template->Clone("h_recgenpt_pho");
  h_dR_jet       = (TH1D*)h_dR_template->Clone("h_dR_jet");
  h_dR_pho       = (TH1D*)h_dR_template->Clone("h_dR_pho");
  h_dR_gj        = (TH1D*)h_dR_template->Clone("h_dR_gj");
  

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

    h_ptclosestjet[idec]  = (TH1D*)h_pt_template->Clone(Form("h_ptclosestjet_%s",decName[idec]));
    h_etaclosestjet[idec] = (TH1D*)h_eta_template->Clone(Form("h_etaclosestjet_%s",decName[idec])); 
    h_dRclosestjet[idec]  = (TH1D*)h_dR_template->Clone(Form("h_dRclosestjet_%s",decName[idec])); 

    h_njetraw[idec] = (TH1I*)h_njet_template->Clone(Form("h_njetraw_%s",decName[idec]));

    for(int ip=0; ip<2; ip++){

      h_alljetpt_eff[idec][ip] = (TH1D*)h_pt_template->Clone(Form("h_alljetpt_eff_%s_%d",decName[idec],ip));
      h_alljeteta_eff[idec][ip] = (TH1D*)h_eta_template->Clone(Form("h_alljeteta_eff_%s_%d",decName[idec],ip));
      h_cost_COMZ_alljet[idec][ip] = (TH1D*)h_cost_template->Clone(Form("h_cost_COMZ_alljet_%s_%d",decName[idec],ip));
      h_pstar_COMZ_alljet[idec][ip] = (TH1D*)h_pstar_template->Clone(Form("h_pstar_COMZ_alljet_%s_%d",decName[idec],ip));
      h_ystar_COMZ_alljet[idec][ip] = (TH1D*)h_ystar_template->Clone(Form("h_ystar_COMZ_alljet_%s_%d",decName[idec],ip));


      h_njet[idec][ip] = (TH1I*)h_njet_template->Clone(Form("h_njet_%s_%d",decName[idec],ip));

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
    h_pthat->Fill(pthat); // make sure we are checking the right MC samples


    //vertex selection
    if(IsVtxGood==0) continue;
    nPass[1]++;

    int nVtxGood = 0;
    for(int iv=0; iv< nVtx; iv++){

      Float_t d0 = vtx[iv][0]*vtx[iv][0] + vtx[iv][1]*vtx[iv][1];
      d0 = sqrt(d0);
      if (fabs(vtx[iv][2]) <= 24 
	  && d0) nVtxGood++;

    }

    //-------------------------------------------------------------------------
    // now find a good leading photon that is matched to hard scattering photon
    //-------------------------------------------------------------------------

    int leadingPhotonIndex = -1;
    double phoMaxPt = -9999.;

    for(unsigned int ipho=0; ipho < nPho; ipho++)
      {	

	int decIndex = phoDecCode(ipho);
	if(!isFidPho(ipho))continue;

	h_sieie[decIndex] -> Fill(phoSigmaIEtaIEta[ipho]);
	h_eciso[decIndex] -> Fill(phoEcalIso(ipho,false));
	h_hciso[decIndex] -> Fill(phoHcalIso(ipho,false));
	h_tkiso[decIndex] -> Fill(phoTrkIso(ipho,false));
	h_hovere[decIndex]-> Fill(phoHoverE[ipho]);
	h_pixel[decIndex] -> Fill(phohasPixelSeed[ipho]);

	// make sure this is the signal photon we care about
 	if(phoGenIndex[ipho]<0)continue;
 	if(phoGenMomPID[ipho]!=22)continue;

	double thisPhoPt= phoEt[ipho]; 

	if(thisPhoPt > phoMaxPt)
	  {
	    phoMaxPt = thisPhoPt;
	    leadingPhotonIndex= ipho;
	  }

      }
    // end of leading photon search

    if(leadingPhotonIndex<0)continue;
    nPass[2]++;

    double leadingPhotonEt = phoEt[leadingPhotonIndex];
    int phoPtBinIndex =  h_ptbin_template->GetXaxis()->FindBin(leadingPhotonEt)-1; 

    double leadingPhotonEta = phoEta[leadingPhotonIndex];
    int phoDecBinIndex = phoDecCode(leadingPhotonIndex);

    h_sieie[phoDecBinIndex+2] -> Fill(phoSigmaIEtaIEta[leadingPhotonIndex]);
    h_hovere[phoDecBinIndex+2]-> Fill(phoHoverE[leadingPhotonIndex]);
    h_pixel[phoDecBinIndex+2] -> Fill(phohasPixelSeed[leadingPhotonIndex]);

    h_eciso[phoDecBinIndex+2] -> Fill(phoEcalIso(leadingPhotonIndex,false));
    h_hciso[phoDecBinIndex+2] -> Fill(phoHcalIso(leadingPhotonIndex,false));
    h_tkiso[phoDecBinIndex+2] -> Fill(phoTrkIso(leadingPhotonIndex,false));


   // after pileup correction
    h_eciso[phoDecBinIndex+4] -> Fill(phoEcalIso(leadingPhotonIndex,true));
    h_hciso[phoDecBinIndex+4] -> Fill(phoHcalIso(leadingPhotonIndex,true));
    h_tkiso[phoDecBinIndex+4] -> Fill(phoTrkIso(leadingPhotonIndex,true));


    h_ptpho->Fill(leadingPhotonEt);
    h_etapho->Fill(leadingPhotonEta);


    // assign 4-vector of leading photon
    TLorentzVector l4_pho(0,0,0,0);
    l4_pho.SetPtEtaPhiE(
			phoEt[leadingPhotonIndex],
			phoEta[leadingPhotonIndex],
			phoPhi[leadingPhotonIndex],
		        phoE[leadingPhotonIndex]
			);


    //-------------------------------------------------------------------------
    // find the reco jet that has the highest matched gen jet pt
    //-------------------------------------------------------------------------

    Int_t leadingJetIndex=-1;
    double genJetMaxPt=-9999.;
    int nGoodJet=0;
    double closestPt = -9999.;
    double closestEta = -9999.;
    double closestdR = 9999.;
     
    TLorentzVector l4_alljet(0,0,0,0);

    for(int ijet=0; ijet < nJet; ijet++){

      if(!isFidJet(ijet))continue;
      if(jetGenJetIndex[ijet] < 1)continue;

      TLorentzVector l4_thisjet(0,0,0,0);

      l4_thisjet.SetPtEtaPhiE(
			      jetPt[ijet],
			      jetEta[ijet],
			      jetPhi[ijet],
			      jetEn[ijet]
			      );

      double thisdR = l4_thisjet.DeltaR(l4_pho);
      h_dR_gj->Fill(thisdR);
      if(thisdR < 0.5)continue; // remove the overlap between photon and jet


      l4_alljet += l4_thisjet;

      
      nGoodJet++;
      h_jetparton->Fill(jetGenPartonID[ijet]);
      
      double thisGenJetPt = jetGenJetPt[ijet];
      // study the properties of the closest jet

      if(thisdR < closestdR)
	{
	  closestPt  = jetPt[ijet];
	  closestEta = jetEta[ijet];
	  closestdR  = thisdR;
	}



      // check which matched recojet has the highest genjet pt
      if(thisGenJetPt > genJetMaxPt)
 	{
 	  genJetMaxPt = thisGenJetPt;
 	  leadingJetIndex= ijet;
 	}


    } // end of loop over jets

    h_njetraw[phoDecBinIndex]->Fill(nGoodJet);

    // need to have a leading photon and a leading jet
    if(leadingPhotonIndex<0 || leadingJetIndex<0)continue;
    nPass[3]++;


    double leadingJetPt  = jetPt[leadingJetIndex];
    double leadingJetEta = jetEta[leadingJetIndex];

    h_ptjet->Fill(leadingJetPt);
    h_etajet->Fill(leadingJetEta);
    h_leadingjetparton->Fill(jetGenPartonID[leadingJetIndex]);

    h_ptclosestjet[phoDecBinIndex]->Fill(closestPt);
    h_etaclosestjet[phoDecBinIndex]->Fill(closestEta);
    h_dRclosestjet[phoDecBinIndex]->Fill(closestdR);

    
    // assign 4-vector

    TLorentzVector l4_1stjet(0,0,0,0);
    l4_1stjet.SetPtEtaPhiE(
			   jetPt[leadingJetIndex],
			   jetEta[leadingJetIndex],
			   jetPhi[leadingJetIndex],
			   jetEn[leadingJetIndex]
			   );


    //-------------------------------------------------------------------------
    // Preselections                                    
    //-------------------------------------------------------------------------

    if(phoPtBinIndex < 0 || phoPtBinIndex > (nPtBins-1))continue;
    if(phoDecBinIndex < 0)continue;
    nPass[phoPtBinIndex+4]++;

    
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

    double allgj_cost_comZ = eiko::cosThetaStar_ZBoostToCM(l4_pho, l4_alljet);
    double allgj_pstar_comZ = eiko::pstar_ZBoostToCM(l4_pho, l4_alljet);
    double allgj_ystar_comZ = eiko::ystar_ZBoostToCM(l4_pho, l4_alljet);

    // before applying cuts

    h_nvtx_eff[phoDecBinIndex][0]->Fill(nVtxGood);

    h_pt_eff[phoDecBinIndex][0]->Fill(leadingPhotonEt);
    h_eta_eff[phoPtBinIndex][0]->Fill(leadingPhotonEta);

    h_jetpt_eff[0]->Fill(leadingJetPt);
    h_jeteta_eff[0]->Fill(leadingJetEta);

    h_njet[phoDecBinIndex][0]->Fill(nGoodJet);

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

    // lumping all jets together
    h_alljetpt_eff[phoDecBinIndex][0]->Fill(l4_alljet.Pt());
    h_alljeteta_eff[phoDecBinIndex][0]->Fill(l4_alljet.Eta());
    h_cost_COMZ_alljet[phoDecBinIndex][0]->Fill(allgj_cost_comZ);
    h_pstar_COMZ_alljet[phoDecBinIndex][0]->Fill(allgj_pstar_comZ);
    h_ystar_COMZ_alljet[phoDecBinIndex][0]->Fill(allgj_ystar_comZ);


    // checking the isolation behaviour as a function of nvertex
    pf_nvtx_eciso[phoDecBinIndex][0]->Fill(nVtxGood,
					   phoEcalIso(leadingPhotonIndex,false));
    pf_nvtx_hciso[phoDecBinIndex][0]->Fill(nVtxGood,
					   phoHcalIso(leadingPhotonIndex,false));
    pf_nvtx_tkiso[phoDecBinIndex][0]->Fill(nVtxGood,
					   phoTrkIso(leadingPhotonIndex,false));

    pf_nvtx_eciso[phoDecBinIndex][1]->Fill(nVtxGood,
					   phoEcalIso(leadingPhotonIndex,true));
    pf_nvtx_hciso[phoDecBinIndex][1]->Fill(nVtxGood,
					   phoHcalIso(leadingPhotonIndex,true));
    pf_nvtx_tkiso[phoDecBinIndex][1]->Fill(nVtxGood,
					   phoTrkIso(leadingPhotonIndex,true));




    // after applying cuts

    //=============================================================
    if(!eiko::separated(l4_pho, l4_1stjet))continue;
    if(!isGoodPho(leadingPhotonIndex, applyPileUpCorr))continue;
    if(!isGoodLooseJet(leadingJetIndex))continue;
    //=============================================================

    h_nvtx_eff[phoDecBinIndex][1]->Fill(nVtxGood);

    h_pt_eff[phoDecBinIndex][1]->Fill(leadingPhotonEt);
    h_eta_eff[phoPtBinIndex][1]->Fill(leadingPhotonEta);

    h_jetpt_eff[1]->Fill(leadingJetPt);
    h_jeteta_eff[1]->Fill(leadingJetEta);

    h_njet[phoDecBinIndex][1]->Fill(nGoodJet);

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


    // lumping all jets together
    h_alljetpt_eff[phoDecBinIndex][1]->Fill(l4_alljet.Pt());
    h_alljeteta_eff[phoDecBinIndex][1]->Fill(l4_alljet.Eta());
    h_cost_COMZ_alljet[phoDecBinIndex][1]->Fill(allgj_cost_comZ);
    h_pstar_COMZ_alljet[phoDecBinIndex][1]->Fill(allgj_pstar_comZ);
    h_ystar_COMZ_alljet[phoDecBinIndex][1]->Fill(allgj_ystar_comZ);



  } // end of loop over entries

  for(int i=0;i<30;i++)if(nPass[i]>0)cout << "nPass[" << i << "] = " << nPass[i] << endl;


  //---------------------------------------------------------------------------------------
  // Write these histograms to an output root file, the output will be used for efficiency
  //---------------------------------------------------------------------------------------
     
  
  std::string remword  ="/data4/syu/7TeV_pythiaMC/";
  std::string remword2 ="/data2/syu/7TeV_madgraphMC/";

  size_t pos  = _inputFileName.find(remword);
  size_t pos2 = _inputFileName.find(remword2);

  if(pos!= std::string::npos)
    _inputFileName.swap(_inputFileName.erase(pos,remword.length()));

  else if(pos2!= std::string::npos)
    _inputFileName.swap(_inputFileName.erase(pos2,remword2.length()));

  else
    _inputFileName = "test.root";
      


  if(applyPileUpCorr)_inputFileName = "pileCorr_" + _inputFileName;
  if(applyCOMCut)_inputFileName = "withCOMCut_" + _inputFileName ;
  TFile* outFile = new TFile(Form("/home/syu/CVSCode/ggeff_%s",_inputFileName.data()),"recreate");               
  h_jetparton->Write();
  h_leadingjetparton->Write();
  h_pthat->Write();
  h_ptpho->Write();
  h_ptjet->Write();
  h_etapho->Write();
  h_etajet->Write();

  h_dR_gj->Write();
//   h_recgenpt_jet -> Write();
//   h_recgenpt_pho -> Write(); 
//   h_dR_jet       -> Write(); 
//   h_dR_pho       -> Write();


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

    h_ptclosestjet[idec]->Write();             
    h_etaclosestjet[idec]->Write();            
    h_dRclosestjet[idec]->Write();             
    h_njetraw[idec]->Write();

    for(int ip=0; ip<2; ip++){
      h_njet[idec][ip]->Write();
      pf_nvtx_eciso[idec][ip]->Write();
      pf_nvtx_hciso[idec][ip]->Write();
      pf_nvtx_tkiso[idec][ip]->Write();

      h_alljetpt_eff[idec][ip]->Write();   
      h_alljeteta_eff[idec][ip]->Write();
      h_cost_COMZ_alljet[idec][ip]->Write();
      h_pstar_COMZ_alljet[idec][ip]->Write();
      h_ystar_COMZ_alljet[idec][ip]->Write();

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

Int_t  gg_angularmc_eff::phoDecCode(Int_t ipho)
{
  double eta = phoSCEta[ipho];
  if(fabs(eta) < BARREL_MAXETA)return 0;
  if(fabs(eta) > ENDCAP_MINETA && 
     fabs(eta) < ENDCAP_MAXETA)return 1;
  return -1;

}

Bool_t gg_angularmc_eff::isFidPho (Int_t ipho)
{
  if(phoDecCode(ipho)<0)return false;
  double et  = phoEt[ipho];  
  if(et < ptbound[0])return false;
  if(et > ptbound[nPtBins])return false;
  return true;
}



Bool_t gg_angularmc_eff::isGoodPho(Int_t ipho, bool applyPileUpCorr)
{
  double et  = phoEt[ipho];

  if(!isFidPho(ipho))return false;
  if(phoHoverE[ipho] > 0.05)return false;
  if(phohasPixelSeed[ipho]==1)return false; // this should be saved as bool

  double eciso = phoEcalIso(ipho, applyPileUpCorr); 
  double hciso = phoHcalIso(ipho, applyPileUpCorr);
  double tkiso = phoTrkIso(ipho, applyPileUpCorr);

  if(eciso > 4.2 +0.006*et)return false;
  if(hciso > 2.2 +0.0025*et)return false;
  if(tkiso > 2.0 +0.001*et)return false;

  return true;
}

Double_t  gg_angularmc_eff::phoEcalIso(Int_t ipho, bool applyPileUpCorr)
{
  if(!isFidPho(ipho))return -9999.0;

  double eciso = phoEcalIsoDR04[ipho];
  int decBinIndex = phoDecCode(ipho);

  if(applyPileUpCorr)
    eciso -= aeff[decBinIndex][0]*rho25; 
  return eciso;
}
 
Double_t gg_angularmc_eff::phoHcalIso(Int_t ipho, bool applyPileUpCorr)
{
  if(!isFidPho(ipho))return -9999.0;

  double hciso = phoHcalIsoDR04[ipho];
  int decBinIndex = phoDecCode(ipho);

  if(applyPileUpCorr)
    hciso -= aeff[decBinIndex][1]*rho25; 
  return hciso;

}

Double_t gg_angularmc_eff::phoTrkIso(Int_t ipho, bool applyPileUpCorr)
{
  if(!isFidPho(ipho))return -9999.0;

  double tkiso = phoTrkIsoHollowDR04[ipho];
  int decBinIndex = phoDecCode(ipho);

  if(applyPileUpCorr)
    tkiso -= aeff[decBinIndex][2]*rho25; 
  return tkiso;

}


Bool_t gg_angularmc_eff::isFidJet (Int_t ijet)
{
  if(fabs(jetEta[ijet]) > 2.4)return false;
  if(jetPt[ijet] < 30.0)return false;
  return true;
}

// check if this reco-jet is a good loose jet
Bool_t gg_angularmc_eff::isGoodLooseJet(Int_t ijet)
{
  if(fabs(jetEta[ijet]) > 2.4)return false;
  if(jetPt[ijet] < 30.0)return false;
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
Bool_t gg_angularmc_eff::isGoodMediumJet(Int_t ijet)
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
Bool_t gg_angularmc_eff::isGoodTightJet(Int_t ijet)
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

