#define yj_angularmc_eff_cxx
#include "yj_angularmc_eff.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include "myLib.h"

// const double BARREL_MAXETA=1.4442;
// const double ENDCAP_MINETA=1.566;
// const double ENDCAP_MAXETA=2.5;

const double BARREL_MAXETA=1.44;
const double ENDCAP_MINETA=1.57;
const double ENDCAP_MAXETA=2.3;
const int nDECs = 2;

// double ptbound[]={85., 95., 110., 130., 160., 200.};
double ptbound[]={40., 45., 50., 55., 60., 65., 70., 75., 85., 95., 110., 130., 160., 200.};
const int nPtBins = sizeof(ptbound)/sizeof(ptbound[0])-1;

// effective area for correcting pileup, ECAL/HCAL/Tracker from EWK-11-251
const double aeff[2][3]={
  {0.183, 0.062, 0.0167},
  {0.090, 0.180, 0.032}
};




void yj_angularmc_eff::Loop(bool applyCOMCut, bool applyPileUpCorr)
{
  if (fChain == 0) return;
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

  bool isDirPho = false;
  bool isFraPho = false;
  if(_inputFileName.find("G_Pt")!= std::string::npos || _inputFileName.find("GJets")!= std::string::npos)
    isDirPho = true;
  else if(_inputFileName.find("QCD")!= std::string::npos)
    isFraPho = true;

  if(isDirPho)
    cout << "This is a gamma+jet direct photon MC sample." << endl;
  if(isFraPho)
    cout << "This is a dijet fragmentation photon MC sample." << endl;

  Long64_t nentries = fChain->GetEntriesFast();
  cout << "There are " << nentries << " entries" << endl;

  //---------------------------------------------------------------------------------------------------------------------
  //
  //   Defining templates for each type of physics variable
  // 
  //---------------------------------------------------------------------------------------------------------------------

  TH2D* h2_ystar_pstar_template = new TH2D("h2_ystar_pstar_template","",100.0,-5.0,5.0,500,0,500);

  TH1D* h_njet_template = new TH1D("h_njet_template","",21,-0.5,20.5);

  TH1D* h_pt_template = new TH1D("h_pt_template","",500,0,500);
  TH1D* h_eta_template = new TH1D("h_eta_template","",100,-5.0,5.0);
  TH1D* h_ptbin_template= new TH1D("h_ptbin_template","", nPtBins, ptbound);
  TH1D* h_ptratio_template = new TH1D("h_ptratio_template","",500,0.0,2.5);
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

  //---------------------------------------------------------------------------------------------------------------------
  //
  //   Declaring and defining the histograms to be plotted
  // 
  //---------------------------------------------------------------------------------------------------------------------

  string decName[4];
  decName[0] = "EB";
  decName[1] = "EE";
  decName[2] = "EBPileupCorr";
  decName[3] = "EEPileupCorr";

  // debugging histograms
  TH1D* h_pthat       = (TH1D*)h_pt_template->Clone("h_pthat");
  h_pthat   -> SetXTitle("#hat{p_{T}} of the MC sample [GeV]");
  h_pthat   -> SetTitle("Before applying any selections");

  TH1D* h_ngenjet     = (TH1D*)h_njet_template->Clone("h_ngenjet");
  h_ngenjet -> SetXTitle("Number of generator-level jets");
  h_ngenjet -> SetTitle("Before applying any selections");

  TH1D* h_nrecjet     = (TH1D*)h_njet_template->Clone("h_nrecjet");
  h_nrecjet -> SetXTitle("Number of reconstruction-level jets");
  h_nrecjet -> SetTitle("Before applying any selections");

  TH1D* h_ptpho       = (TH1D*)h_pt_template->Clone("h_ptpho");
  h_ptpho   -> SetXTitle("Reconstructed leading photon p_{T} [GeV]");
  h_ptpho   -> SetTitle("Before applying ID selections");

  TH1D* h_etapho      = (TH1D*)h_eta_template->Clone("h_etapho");
  h_etapho   -> SetXTitle("Reconstructed leading photon #eta");
  h_etapho   -> SetTitle("Before applying ID selections");

  TH1D* h_genciso     = (TH1D*)h_eciso_template->Clone("h_genciso");
  h_genciso -> SetXTitle("Generator-level calorimetric isolation [GeV]");
  h_genciso -> SetTitle("Already reconstructed generator-level leading photons");

  TH1D* h_gentiso     = (TH1D*)h_tkiso_template->Clone("h_gentiso");
  h_gentiso -> SetXTitle("Generator-level track isolation [GeV]");
  h_gentiso -> SetTitle("Already reconstructed generator-level leading photons");


  TH1D* h_sieie[4];
  TH1D* h_eciso[4];
  TH1D* h_hciso[4];
  TH1D* h_tkiso[4];
  TH1D* h_hovere[4];
  TH1D* h_pixel[4];

  for(int ip=0; ip<4; ip++)
    {
      // EB:ip=0, EE:ip=1, pileup corrected ip=2 and 3
      h_sieie[ip] = (TH1D*)h_sieie_template->Clone(Form("h_sieie_%s", decName[ip].data()));
      h_sieie[ip] -> SetXTitle("#sigma_{i#eta i#eta}");

      h_eciso[ip] = (TH1D*)h_eciso_template->Clone(Form("h_eciso_%s", decName[ip].data()));
      h_eciso[ip] -> SetXTitle("ISO_{ECAL} [GeV]");

      h_hciso[ip] = (TH1D*)h_hciso_template->Clone(Form("h_hciso_%s", decName[ip].data()));
      h_hciso[ip] -> SetXTitle("ISO_{HCAL} [GeV]");

      h_tkiso[ip] = (TH1D*)h_tkiso_template->Clone(Form("h_tkiso_%s", decName[ip].data()));
      h_tkiso[ip] -> SetXTitle("ISO_{TRK} [GeV]");

      h_hovere[ip] = (TH1D*)h_hovere_template->Clone(Form("h_hovere_%s", decName[ip].data()));
      h_hovere[ip]-> SetXTitle("E_{HAD}/E_{EM}");

      h_pixel[ip] = (TH1D*)h_pixel_template->Clone(Form("h_pixel_%s", decName[ip].data()));    
      h_pixel[ip] -> SetXTitle("Matched to pixel hits");

      
      std::string tempTitle;
      if(ip<2)
	tempTitle = "Leading photon, before applying #rho-correction";
      else
	tempTitle = "Leading photon, after applying #rho-correction";

      h_sieie[ip]  -> SetTitle(tempTitle.data());
      h_eciso[ip]  -> SetTitle(tempTitle.data());
      h_hciso[ip]  -> SetTitle(tempTitle.data());
      h_tkiso[ip]  -> SetTitle(tempTitle.data());
      h_hovere[ip] -> SetTitle(tempTitle.data());
      h_pixel[ip]  -> SetTitle(tempTitle.data());
    

    }


  TProfile* pf_nvtx_eciso[2][2];       // in EB and EE, before and after pileupcorr
  TProfile* pf_nvtx_hciso[2][2];       // in EB and EE, before and after pileupcorr
  TProfile* pf_nvtx_tkiso[2][2];       // in EB and EE, before and after pileupcorr

  string puName[2];
  puName[0] = "raw";
  puName[1] = "pucorr";

  for(int idec=0; idec < nDECs; idec++){
    for(int ip=0; ip<2; ip++){      

      pf_nvtx_eciso[idec][ip] = (TProfile*)pf_nvtx_eciso_template->Clone
	(Form("pf_nvtx_eciso_%s_%s",decName[idec].data(),puName[ip].data()));
      pf_nvtx_eciso[idec][ip]->SetYTitle("Average ISO_{ECAL} [GeV]");
      
      pf_nvtx_hciso[idec][ip] = (TProfile*)pf_nvtx_hciso_template->Clone
	(Form("pf_nvtx_hciso_%s_%s",decName[idec].data(),puName[ip].data()));
      pf_nvtx_hciso[idec][ip]->SetYTitle("Average ISO_{HCAL} [GeV]");
      
      pf_nvtx_tkiso[idec][ip] = (TProfile*)pf_nvtx_tkiso_template->Clone
	(Form("pf_nvtx_tkiso_%s_%s",decName[idec].data(),puName[ip].data()));
      pf_nvtx_tkiso[idec][ip]->SetYTitle("Average ISO_{TRK} [GeV]");

      std::string tempXTitle = "Number of good vertices";
      pf_nvtx_eciso[idec][ip]->SetXTitle(tempXTitle.data());
      pf_nvtx_hciso[idec][ip]->SetXTitle(tempXTitle.data());
      pf_nvtx_tkiso[idec][ip]->SetXTitle(tempXTitle.data());

      std::string tempTitle;
      if(ip ==0) tempTitle = "Before #rho correction";
      else tempTitle = "After #rho correction";

      pf_nvtx_eciso[idec][ip]->SetTitle(tempTitle.data());
      pf_nvtx_hciso[idec][ip]->SetTitle(tempTitle.data());
      pf_nvtx_tkiso[idec][ip]->SetTitle(tempTitle.data());
      
    } 

  }



  TH1D* h_ptjet       = (TH1D*)h_pt_template->Clone("h_ptjet");
  h_ptjet   -> SetXTitle("Reconstructed jet p_{T} [GeV]");
  h_ptjet   -> SetTitle("Before applying ID selections");

  TH1D* h_etajet      = (TH1D*)h_eta_template->Clone("h_etajet");
  h_etajet   -> SetXTitle("Reconstructed jet #eta");
  h_etajet   -> SetTitle("Before applying ID selections");



  TH2D* h2_ystar_pstar[2]; // in EB and EE
  TH1D* h_ptclosestjet[2];             // pt of the closest jet in deltaR to the leading photon in barrel and endcap
  TH1D* h_etaclosestjet[2];            // eta of the closest jet in deltaR to the leading photon in barrel and endcap
  TH1D* h_dRclosestjet[2];             // deltaR of the closest jet in deltaR to the leading photon in barrel and endcap
  TH1D* h_dRLeadingPhoJet[2];          // deltaR of the leading photon and leading jet

  for(int idec=0; idec < nDECs; idec ++){
    

    h2_ystar_pstar[idec] = (TH2D*)h2_ystar_pstar_template->Clone(Form("h2_ystar_pstar_%s",decName[idec].data()));
    h2_ystar_pstar[idec] -> SetXTitle("y^{*}");
    h2_ystar_pstar[idec] -> SetYTitle("p^{*} [GeV]");

    h_ptclosestjet[idec]  = (TH1D*)h_pt_template->Clone(Form("h_ptclosestjet_%s",decName[idec].data()));
    h_ptclosestjet[idec] -> SetXTitle("p_{T} of the jet closest to leading photon [GeV]");

    h_etaclosestjet[idec] = (TH1D*)h_eta_template->Clone(Form("h_etaclosestjet_%s",decName[idec].data())); 
    h_etaclosestjet[idec]-> SetXTitle("#eta of the jet closest to leading photon"); 

    h_dRclosestjet[idec]  = (TH1D*)h_dR_template->Clone(Form("h_dRclosestjet_%s",decName[idec].data())); 
    h_dRclosestjet[idec] -> SetXTitle("#Delta R between the jet closest to and leading photon"); 

    h_dRLeadingPhoJet[idec] = (TH1D*)h_dR_template->Clone(Form("h_dRLeadingPhoJet_%s",decName[idec].data()));
    h_dRLeadingPhoJet[idec]-> SetXTitle("#Delta R between leading jet and leading photon");
    

    std::string tempTitle = "Before applying ID selections";
    h2_ystar_pstar[idec] -> SetTitle(tempTitle.data());
    h_ptclosestjet[idec] -> SetTitle(tempTitle.data());
    h_etaclosestjet[idec]-> SetTitle(tempTitle.data());
    h_dRclosestjet[idec] -> SetTitle(tempTitle.data());
    h_dRLeadingPhoJet[idec]-> SetTitle(tempTitle.data());

  }

  // histograms for efficiency

  TH1D* h_njet[2][2];                  // in EB and EE, before and after ID selections
  TH1D* h_jetpt_eff[2][2];            
  TH1D* h_jeteta_eff[2][2];
  TH1D* h_sumjetpt_eff[2][2];            
  TH1D* h_sumjeteta_eff[2][2];

  TH1D* h_nvtx_eff[nDECs][2];        // in EB and EE

  // for photon
  TH1D* h_phopt_eff[nDECs][2];        // in EB and EE
  TH1D* h_phoeta_eff[nPtBins+1][2]; // in various pt bins

  TH1D* h_ptratio[nDECs][2];
  TH1D* h_ptratio_sumJet[nDECs][2];

  TH1D* h_zgamma[nDECs][nPtBins+1][2];
  TH1D* h_chi[nDECs][nPtBins+1][2];
  TH1D* h_yB[nDECs][nPtBins+1][2];

  TH1D* h_dphi[nDECs][nPtBins+1][2];
  TH1D* h_dphi_sumJet[nDECs][nPtBins+1][2];

  TH1D* h_cost[nDECs][nPtBins+1][2]; // tanH(0.5(y1-y2))
  TH1D* h_cost_COM3D[nDECs][nPtBins+1][2]; // boost to the COM using total momentum
  TH1D* h_cost_COMZ[nDECs][nPtBins+1][2]; // boost to the COM using z-momentum
  TH1D* h_cost_sumJet[nDECs][nPtBins+1][2]; // boost to the COM using z-momentum, use all fiducial jets

  TH1D* h_pstar[nDECs][nPtBins+1][2];
  TH1D* h_pstar_COM3D[nDECs][nPtBins+1][2]; // boost to the COM using total momentum
  TH1D* h_pstar_COMZ[nDECs][nPtBins+1][2]; // boost to the COM using z-momentum
  TH1D* h_pstar_sumJet[nDECs][nPtBins+1][2]; // boost to the COM using z-momentum, use all fiducial jets

  TH1D* h_ystar[nDECs][nPtBins+1][2];
  TH1D* h_ystar_COM3D[nDECs][nPtBins+1][2]; // boost to the COM using total momentum
  TH1D* h_ystar_COMZ[nDECs][nPtBins+1][2]; // boost to the COM using the z-momentum
  TH1D* h_ystar_sumJet[nDECs][nPtBins+1][2]; // boost to the COM using z-momentum, use all fiducial jets



  for(int ip=0; ip<2; ip++){

    for(int idec=0; idec< nDECs; idec++)
      {

	h_njet[idec][ip] = (TH1D*)h_njet_template->Clone(Form("h_njet_%s_%d", decName[idec].data(), ip));
	h_njet[idec][ip] -> SetXTitle("Number of reconstructed jets");

	h_jetpt_eff[idec][ip] = (TH1D*)h_pt_template->Clone(Form("h_jetpt_eff_%s_%d", decName[idec].data(), ip));
	h_jetpt_eff[idec][ip] -> SetXTitle("Reconstructed jet p_{T} [GeV]");

	h_jeteta_eff[idec][ip] = (TH1D*)h_eta_template->Clone(Form("h_jeteta_eff_%s_%d", decName[idec].data(), ip));
	h_jeteta_eff[idec][ip] -> SetXTitle("Reconstructed jet #eta");

	h_sumjetpt_eff[idec][ip] = (TH1D*)h_pt_template->Clone(Form("h_sumjetpt_eff_%s_%d", decName[idec].data(), ip));
	h_sumjetpt_eff[idec][ip] -> SetXTitle("Reconstructed jet^{sum} p_{T} [GeV]");

	h_sumjeteta_eff[idec][ip] = (TH1D*)h_eta_template->Clone(Form("h_sumjeteta_eff_%s_%d", decName[idec].data(), ip));
	h_sumjeteta_eff[idec][ip] -> SetXTitle("Reconstructed jet^{sum} #eta");

	h_nvtx_eff[idec][ip]=(TH1D*)h_nvtx_template->Clone(Form("h_nvtx_%s_%d",
								decName[idec].data(),ip));
	h_nvtx_eff[idec][ip]   -> SetXTitle("Number of good vertices");

	h_phopt_eff[idec][ip] = (TH1D*)h_pt_template->Clone(Form("h_phopt_eff_%s_%d",
								 decName[idec].data(), ip));
	h_phopt_eff[idec][ip]  -> SetXTitle("Leading photon p_{T} [GeV]");


	h_ptratio[idec][ip]= (TH1D*)h_ptratio_template->Clone(Form("h_ptratio_%s_%d",
								   decName[idec].data(),ip));
	h_ptratio[idec][ip]    -> SetXTitle("p_{T}(jet^{1st})/p_{T}(#gamma)");

	h_ptratio_sumJet[idec][ip]= (TH1D*)h_ptratio_template->Clone(Form("h_ptratio_sumJet%s_%d",
								   decName[idec].data(),ip));
	h_ptratio_sumJet[idec][ip]    -> SetXTitle("p_{T}(jet^{sum})/p_{T}(#gamma)");

      } // end of loop over barrel and endcap

    for(int ipt=0; ipt < nPtBins+1; ipt++){      
      if(ipt == (nPtBins))
	h_phoeta_eff[ipt][ip] = (TH1D*)h_eta_template->Clone(Form("h_phoeta_eff_allpt_%d",ip));
      else
	h_phoeta_eff[ipt][ip] = (TH1D*)h_eta_template->Clone(Form("h_phoeta_eff_%d_%d_%d",
								  (int)ptbound[ipt], (int)ptbound[ipt+1], ip));
      h_phoeta_eff[ipt][ip]  -> SetXTitle("Leading photon #eta");
    } // end of loop over pt bins
  
  }// end of loop over process, ip=0 before cut, 1 after cut
  
  
  for(int idec=0; idec < nDECs; idec ++){
    for(int ipt=0; ipt < nPtBins+1; ipt++){
      for(int ip=0; ip < 2; ip ++){
	double ptmin = -1;
	double ptmax = -1;
	
 	if(ipt == (nPtBins))
 	  {
	  
 	    h_zgamma[idec][ipt][ip] = createHisto(h_zgamma_template, Form("h_zgamma_%s", decName[idec].data()), "z^{#gamma,jet^{1st}}",ip);
	    
 	    h_chi[idec][ipt][ip] = createHisto(h_chi_template, Form("h_chi_%s", decName[idec].data()),"#chi(#gamma,jet^{1st}",ip);
	    
 	    h_yB[idec][ipt][ip] = createHisto(h_yB_template, Form("h_yB_%s", decName[idec].data()), "0.5 (y^{gamma}+y^{1stjet})",ip);
	    
 	    h_dphi[idec][ipt][ip] = createHisto(h_dphi_template, Form("h_dphi_%s", decName[idec].data()), "#Delta#phi(#gamma,jet^{1st})",ip);
	    
	    h_dphi_sumJet[idec][ipt][ip] = createHisto(h_dphi_template, Form("h_dphi_sumJet_%s", decName[idec].data()), 
						       "#Delta#phi(#gamma,jet^{sum})",ip);

 	    h_cost[idec][ipt][ip] = createHisto(h_cost_template, Form("h_cost_%s", decName[idec].data()),
						"tanh[0.5(y^{gamma}-y^{1stjet})]",ip);

 	    h_cost_COM3D[idec][ipt][ip] = createHisto(h_cost_template, Form("h_cost_COM3D_%s", decName[idec].data()),
						      "cos#theta^{*} (in the frame where #gamma and jet^{1st} back-to-back in 3D)",ip);

 	    h_cost_COMZ[idec][ipt][ip] = createHisto(h_cost_template, Form("h_cost_COMZ_%s", decName[idec].data()),
						     "cos#theta^{*} (in the frame where #gamma and jet^{1st} back-to-back in z)",ip);

	    h_cost_sumJet[idec][ipt][ip] = createHisto(h_cost_template, Form("h_cost_sumJet_%s", decName[idec].data()),
						       "cos#theta^{*} (in the frame where #gamma and jet^{sum} back-to-back in z)",ip);

	    h_pstar[idec][ipt][ip] = createHisto(h_pstar_template, Form("h_pstar_%s", decName[idec].data()),
						 "p_{T}(#gamma)cosh(y^{*})", ip);

	    h_pstar_COM3D[idec][ipt][ip] = createHisto(h_pstar_template, Form("h_pstar_COM3D_%s", decName[idec].data()), 
						       "p^{*} (in the frame where #gamma and jet^{1st} back-to-back in 3D)",ip);
	    
	    h_pstar_COMZ[idec][ipt][ip] = createHisto(h_pstar_template, Form("h_pstar_COMZ_%s", decName[idec].data()),
						      "p^{*} (in the frame where #gamma and jet^{1st} back-to-back in z)",ip);

 	    h_pstar_sumJet[idec][ipt][ip] = createHisto(h_pstar_template, Form("h_pstar_sumJet_%s", decName[idec].data()),
							"p^{*} (in the frame where #gamma and jet^{sum} back-to-back in z)",ip);
	    
	    h_ystar[idec][ipt][ip] = createHisto(h_ystar_template, Form("h_ystar_%s", decName[idec].data()),
						 "0.5(y^{#gamma} - y^{1stjet})", ip);

	    h_ystar_COM3D[idec][ipt][ip] = createHisto(h_ystar_template, Form("h_ystar_COM3D_%s", decName[idec].data()), 
						       "y^{*} (in the frame where #gamma and jet^{1st} back-to-back in 3D)",ip); 
			
	    h_ystar_COMZ[idec][ipt][ip] = createHisto(h_ystar_template, Form("h_ystar_COMZ_%s", decName[idec].data()),
						      "y^{*} (in the frame where #gamma and jet^{1st} back-to-back in z)", ip);

	    h_ystar_sumJet[idec][ipt][ip] = createHisto(h_ystar_template, Form("h_ystar_sumJet_%s", decName[idec].data()),
							"y^{*} (in the frame where #gamma and jet^{sum} back-to-back in z)", ip);
	    
 	  }

 	else
 	  {
	    ptmin = ptbound[ipt];
	    ptmax = ptbound[ipt+1];

	    h_zgamma[idec][ipt][ip] = createHisto(h_zgamma_template, Form("h_zgamma_%s", decName[idec].data()), "z^{#gamma,jet^{1st}}",ip, ptmin, ptmax);
	    
	    h_chi[idec][ipt][ip] = createHisto(h_chi_template, Form("h_chi_%s", decName[idec].data()),"#chi(#gamma,jet^{1st}",ip, ptmin, ptmax);
	    
	    h_yB[idec][ipt][ip] = createHisto(h_yB_template, Form("h_yB_%s", decName[idec].data()), "0.5 (y^{gamma}+y^{1stjet})",ip, ptmin, ptmax);
	    
 	    h_dphi[idec][ipt][ip] = createHisto(h_dphi_template, Form("h_dphi_%s", decName[idec].data()), "#Delta#phi(#gamma,jet^{1st})",ip, ptmin, ptmax);	    
 	    h_dphi_sumJet[idec][ipt][ip] = createHisto(h_dphi_template, Form("h_dphi_sumJet_%s", decName[idec].data()), 
						       "#Delta#phi(#gamma,jet^{sum})",ip, ptmin, ptmax);

 	    h_cost[idec][ipt][ip] = createHisto(h_cost_template, Form("h_cost_%s", decName[idec].data()),
						"tanh[0.5(y^{gamma}-y^{1stjet})]",ip, ptmin, ptmax);

 	    h_cost_COM3D[idec][ipt][ip] = createHisto(h_cost_template, Form("h_cost_COM3D_%s", decName[idec].data()),
						      "cos#theta^{*} (in the frame where #gamma and jet^{1st} back-to-back in 3D)",
						      ip, ptmin, ptmax);

 	    h_cost_COMZ[idec][ipt][ip] = createHisto(h_cost_template, Form("h_cost_COMZ_%s", decName[idec].data()),
						     "cos#theta^{*} (in the frame where #gamma and jet^{1st} back-to-back in z)",
						     ip, ptmin, ptmax);
	    
 	    h_cost_sumJet[idec][ipt][ip] = createHisto(h_cost_template, Form("h_cost_sumJet_%s", decName[idec].data()),
						       "cos#theta^{*} (in the frame where #gamma and jet^{sum} back-to-back in z)",
						       ip, ptmin, ptmax);

 	    h_pstar[idec][ipt][ip] = createHisto(h_pstar_template, Form("h_pstar_%s", decName[idec].data()),
						 "p_{T}(#gamma)cosh(y^{*})", ip, ptmin, ptmax);

 	    h_pstar_COM3D[idec][ipt][ip] = createHisto(h_pstar_template, Form("h_pstar_COM3D_%s", decName[idec].data()), 
						       "p^{*} (in the frame where #gamma and jet^{1st} back-to-back in 3D)",ip, ptmin, ptmax);
	    
 	    h_pstar_COMZ[idec][ipt][ip] = createHisto(h_pstar_template, Form("h_pstar_COMZ_%s", decName[idec].data()),
						      "p^{*} (in the frame where #gamma and jet^{1st} back-to-back in z)",ip, ptmin, ptmax);

	    h_pstar_sumJet[idec][ipt][ip] = createHisto(h_pstar_template, Form("h_pstar_sumJet_%s", decName[idec].data()),
							"p^{*} (in the frame where #gamma and jet^{sum} back-to-back in z)",ip, ptmin, ptmax);
	    
	    h_ystar[idec][ipt][ip] = createHisto(h_ystar_template, Form("h_ystar_%s", decName[idec].data()),
						 "0.5(y^{#gamma} - y^{1stjet})", ip, ptmin, ptmax);

	    h_ystar_COM3D[idec][ipt][ip] = createHisto(h_ystar_template, Form("h_ystar_COM3D_%s", decName[idec].data()), 
						       "y^{*} (in the frame where #gamma and jet^{1st} back-to-back in 3D)",ip, ptmin, ptmax); 
	    
	    h_ystar_COMZ[idec][ipt][ip] = createHisto(h_ystar_template, Form("h_ystar_COMZ_%s", decName[idec].data()),
						      "y^{*} (in the frame where #gamma and jet^{1st} back-to-back in z)", ip, ptmin, ptmax);

	    h_ystar_sumJet[idec][ipt][ip] = createHisto(h_ystar_template, Form("h_ystar_sumJet_%s", decName[idec].data()),
							"y^{*} (in the frame where #gamma and jet^{sum} back-to-back in z)", ip, ptmin, ptmax);

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
    h_ngenjet->Fill(genJetPt_->size());
    h_nrecjet->Fill(patJetPfAk05Pt_->size());

    // Find a good vertex first
    if(EvtInfo_nVtxGood<1) continue;
    nPass[1]++;


    int leadingPhotonIndex = -1;
    double phoMaxPt = -9999.;

    //-------------------------------------------------------------------------
    // now find a good leading photon that is matched to hard scattering photon
    //-------------------------------------------------------------------------
    for(unsigned int ipho=0; ipho < PhotonPt->size(); ipho++)
      {	
	// first require this reconstruction-level photon to be fiducial
        // and has pt within the boundary we are studying
	if(!isFidPho(ipho))continue;

	// then check if this photon is matched to a gen photon
 	if(!PhotonisGenMatched->at(ipho))continue;
	
// 	cout << "gen mom Id = " << PhotongenMomId->at(ipho) << endl;
	// require it to match to a prompt photon
 	if(fabs(PhotongenMomId->at(ipho))>22.1)continue; // not prompt photon 
// 	cout << "is a prompt photon" << endl;

	// find the leading photon, usually there's only one in the MC
	double thisPhoPt= PhotonEt->at(ipho);      
	if(thisPhoPt > phoMaxPt)
	  {
	    phoMaxPt = thisPhoPt;
	    leadingPhotonIndex= ipho;
	  }

      }
    // end of leading photon search
//     cout << "leadingPhotonIndex = " << leadingPhotonIndex << endl;
    if(leadingPhotonIndex<0)continue;
    nPass[2]++;

    if(leadingPhotonIndex>0)nPass[3]++;

    double leadingPhotonEt = PhotonEt->at(leadingPhotonIndex);
    int phoPtBinIndex =  h_ptbin_template->GetXaxis()->FindBin(leadingPhotonEt)-1; 
    h_ptpho->Fill(leadingPhotonEt);

    double leadingPhotonEta = PhotonEta->at(leadingPhotonIndex);
    int phoDecBinIndex = phoDecCode(leadingPhotonIndex);
    h_etapho->Fill(leadingPhotonEta);

    h_genciso->Fill(PhotongenCalIsoDR04->at(leadingPhotonIndex));
    h_gentiso->Fill(PhotongenTrkIsoDR04->at(leadingPhotonIndex));

    // dummy proof, but already make the same requirement in isFidPho
    if(phoPtBinIndex < 0 || phoPtBinIndex > (nPtBins-1))continue;
    if(phoDecBinIndex < 0)continue;

    // check the photon ID variable before making cuts on them
    h_sieie[phoDecBinIndex] -> Fill(PhotonSigmaIetaIeta->at(leadingPhotonIndex));
    h_hovere[phoDecBinIndex]-> Fill(PhotonhadronicOverEm->at(leadingPhotonIndex));
    h_pixel[phoDecBinIndex] -> Fill(PhotonhasPixelSeed->at(leadingPhotonIndex));


    double eciso = phoEcalIso(leadingPhotonIndex,false);
    double hciso = phoHcalIso(leadingPhotonIndex,false);
    double tkiso = phoTrkIso(leadingPhotonIndex,false);

    double eciso_corr = phoEcalIso(leadingPhotonIndex,true);
    double hciso_corr = phoHcalIso(leadingPhotonIndex,true);
    double tkiso_corr = phoTrkIso(leadingPhotonIndex,true);


    h_eciso[phoDecBinIndex] -> Fill(eciso);
    h_hciso[phoDecBinIndex] -> Fill(hciso);
    h_tkiso[phoDecBinIndex] -> Fill(tkiso);

    // checking the isolation behaviour as a function of nvertex
    pf_nvtx_eciso[phoDecBinIndex][0]->Fill(EvtInfo_nVtxGood,eciso);
    pf_nvtx_hciso[phoDecBinIndex][0]->Fill(EvtInfo_nVtxGood,hciso);
    pf_nvtx_tkiso[phoDecBinIndex][0]->Fill(EvtInfo_nVtxGood,tkiso);

    // after pileup correction
    h_eciso[phoDecBinIndex+2] -> Fill(eciso_corr);
    h_hciso[phoDecBinIndex+2] -> Fill(hciso_corr);
    h_tkiso[phoDecBinIndex+2] -> Fill(tkiso_corr);


    // checking the isolation behaviour as a function of nvertex
    pf_nvtx_eciso[phoDecBinIndex][1]->Fill(EvtInfo_nVtxGood,eciso_corr);
    pf_nvtx_hciso[phoDecBinIndex][1]->Fill(EvtInfo_nVtxGood,hciso_corr);
    pf_nvtx_tkiso[phoDecBinIndex][1]->Fill(EvtInfo_nVtxGood,tkiso_corr);


    // assign 4-vector of leading photon
    TLorentzVector l4_pho(0,0,0,0);
    l4_pho.SetPtEtaPhiE(
			PhotonEt->at(leadingPhotonIndex),
			PhotonEta->at(leadingPhotonIndex),
			PhotonPhi->at(leadingPhotonIndex),
			PhotonEnergy->at(leadingPhotonIndex)
			);

    //-------------------------------------------------------------------------
    // find the gen jet index that has the highest pt and also find the matched
    // reco jet, this gen jet needs to have pt > 10 GeV and |eta| < 3.0
    //-------------------------------------------------------------------------

    double genJetMaxPt     = -9999.0;
    int leadingGenJetIndex = -1;

    // loop over generator-level jets
    for(unsigned int ij=0; ij < genJetPt_->size(); ij++)
      {
	double thisJetPt = genJetPt_->at(ij);
	double thisJetEta= genJetEta_->at(ij);

	if(thisJetPt < 10.0)continue;
	if(fabs(thisJetEta) > 3.0)continue;

	// try to find the genJet that the highest gen pt
	if(thisJetPt > genJetMaxPt)
	  {
	    genJetMaxPt                   = thisJetPt;
	    leadingGenJetIndex            = ij;
	  }


      } // end of loop over genJets

    if(leadingGenJetIndex <0)continue; // couldn't find a gen jet with pt > 10 GeV and |eta| < 3.0
    nPass[4]++;


    // find the reconstruction-level jet that is matched to this highest pt genJet
    int leadingJetIndex = matchGenToRecoJet(leadingGenJetIndex);
    // need to have a leading photon and a leading jet
    if(leadingPhotonIndex<0 || leadingJetIndex<0)continue;
    nPass[5]++;

    if( !isFidJet(leadingJetIndex) )continue;
    nPass[6]++;

    double leadingJetPt  = patJetPfAk05Pt_->at(leadingJetIndex);
    h_ptjet->Fill(leadingJetPt);
    double leadingJetEta = patJetPfAk05Eta_->at(leadingJetIndex);
    h_etajet->Fill(leadingJetEta);
    
    // assign 4-vector of leading jet
    TLorentzVector l4_1stjet(0,0,0,0);
    l4_1stjet.SetPtEtaPhiE(
			   patJetPfAk05Pt_->at(leadingJetIndex),
			   patJetPfAk05Eta_->at(leadingJetIndex),
			   patJetPfAk05Phi_->at(leadingJetIndex),
			   patJetPfAk05En_->at(leadingJetIndex)
			   );



    //-------------------------------------------------------------------------
    // now loop over all reco jets and find the sum momentum of the jets that 
    // satisfy isFidJet 
    // and is matched to genjet?
    //-------------------------------------------------------------------------

    TLorentzVector l4_sumjet(0,0,0,0);
    int nFiducialJets = 0; // number of jets passing pt > 30 GeV, |eta| < 2.4, and matched to gen jet
    double closestPt = -9999.;
    double closestEta = -9999.;
    double closestdR = 9999.;
    bool findJetComponent = false;
    
    for(unsigned int ijet=0; ijet < patJetPfAk05Pt_->size(); ijet++){

      if(!isFidJet(ijet))continue;
      int genJetIndex = matchRecoToGenJet(ijet);
      if(genJetIndex < 0)continue; // is this needed

      nFiducialJets++;

      TLorentzVector l4_thisJet(0,0,0,0);
      l4_thisJet.SetPtEtaPhiE(
			      patJetPfAk05Pt_->at(ijet),
			      patJetPfAk05Eta_->at(ijet),
			      patJetPfAk05Phi_->at(ijet),
			      patJetPfAk05En_->at(ijet)
			      );
      
      double thisdR = l4_thisJet.DeltaR(l4_pho);
      if(thisdR < 0.5)continue; // remove the overlap between photon and jet
      l4_sumjet += l4_thisJet;
      findJetComponent = true;

      // study the properties of the closest jet
      if(thisdR < closestdR)
	{
	  closestPt  = patJetPfAk05Pt_->at(ijet);
	  closestEta = patJetPfAk05Eta_->at(ijet);
	  closestdR  = thisdR;
	}
      
    } // end of loop over reconstructed jets

    double sumJetPt  = findJetComponent? l4_sumjet.Pt(): -9999.0;
    double sumJetEta = findJetComponent? l4_sumjet.Eta(): -9999.0;

    bool findSumJet = false;
    if(findJetComponent && sumJetPt > 30.0 && fabs(sumJetEta) < 2.4)
      findSumJet = true;

    
    //-------------------------------------------------------------------------
    // Preselections                                    
    //-------------------------------------------------------------------------
    
    if(!eiko::separated(l4_pho, l4_1stjet))continue;

    double gj_pstar_comZ = eiko::pstar_ZBoostToCM(l4_pho, l4_1stjet);
    double gj_ystar_comZ = eiko::ystar_ZBoostToCM(l4_pho, l4_1stjet);

    bool passCOMCut = gj_pstar_comZ > pstarMin && gj_pstar_comZ < pstarMax
      && fabs(gj_ystar_comZ) < ystarMax;


    if(applyCOMCut && !passCOMCut)continue;


    //-------------------------------------------------------------------------
    // Now we plot distributions before and after cuts
    //-------------------------------------------------------------------------


    double gj_zgamma = eiko::zgamma(l4_pho, l4_1stjet);
    double gj_dphi   = eiko::deltaPhi(l4_pho, l4_1stjet);
    double gj_dR     = l4_pho.DeltaR(l4_1stjet);

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

    double ptRatio_sumJet   = sumJetPt/leadingPhotonEt;
    double gsumj_dphi       = eiko::deltaPhi(l4_pho, l4_sumjet);
    double gsumj_cost_comZ  = eiko::cosThetaStar_ZBoostToCM(l4_pho, l4_sumjet);
    double gsumj_pstar_comZ = eiko::pstar_ZBoostToCM(l4_pho, l4_sumjet);
    double gsumj_ystar_comZ = eiko::ystar_ZBoostToCM(l4_pho, l4_sumjet);
	

    // before applying cuts
    h2_ystar_pstar[phoDecBinIndex]->Fill(gj_ystar_comZ, gj_pstar_comZ);
    h_ptclosestjet[phoDecBinIndex]->Fill(closestPt);
    h_etaclosestjet[phoDecBinIndex]->Fill(closestEta);
    h_dRclosestjet[phoDecBinIndex]->Fill(closestdR);
    h_dRLeadingPhoJet[phoDecBinIndex]->Fill(gj_dR);


    h_njet[phoDecBinIndex][0]->Fill(nFiducialJets);
    h_jetpt_eff[phoDecBinIndex][0]->Fill(leadingJetPt);
    h_jeteta_eff[phoDecBinIndex][0]->Fill(leadingJetEta);

    h_nvtx_eff[phoDecBinIndex][0]->Fill(EvtInfo_nVtxGood);

    h_phopt_eff[phoDecBinIndex][0]->Fill(leadingPhotonEt);
    h_phoeta_eff[phoPtBinIndex][0]->Fill(leadingPhotonEta);

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

    // lumping all fiducial jets together
    if(findSumJet){
      h_sumjetpt_eff[phoDecBinIndex][0]->Fill(sumJetPt);
      h_sumjeteta_eff[phoDecBinIndex][0]->Fill(sumJetEta);
      h_ptratio_sumJet[phoDecBinIndex][0]->Fill(ptRatio_sumJet);
      h_dphi_sumJet[phoDecBinIndex][phoPtBinIndex][0]->Fill(gsumj_dphi);
      h_cost_sumJet[phoDecBinIndex][phoPtBinIndex][0]->Fill(gsumj_cost_comZ);
      h_pstar_sumJet[phoDecBinIndex][phoPtBinIndex][0]->Fill(gsumj_pstar_comZ);
      h_ystar_sumJet[phoDecBinIndex][phoPtBinIndex][0]->Fill(gsumj_ystar_comZ);
    }

    // lump all pt together
    h_phoeta_eff[nPtBins][0]->Fill(leadingPhotonEta);
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

    // lumping all fiducial jets and all pt together
    if(findSumJet){
      h_dphi_sumJet[phoDecBinIndex][nPtBins][0]->Fill(gsumj_dphi);
      h_cost_sumJet[phoDecBinIndex][nPtBins][0]->Fill(gsumj_cost_comZ);
      h_pstar_sumJet[phoDecBinIndex][nPtBins][0]->Fill(gsumj_pstar_comZ);
      h_ystar_sumJet[phoDecBinIndex][nPtBins][0]->Fill(gsumj_ystar_comZ);
    }
    
    // after applying cuts

    //=============================================================
    if(!eiko::separated(l4_pho, l4_1stjet))continue;
    if(!isGoodPho(leadingPhotonIndex, applyPileUpCorr))continue;
    if(!isGoodLooseJet(leadingJetIndex))continue;
    //=============================================================

    h_njet[phoDecBinIndex][1]->Fill(nFiducialJets);
    h_jetpt_eff[phoDecBinIndex][1]->Fill(leadingJetPt);
    h_jeteta_eff[phoDecBinIndex][1]->Fill(leadingJetEta);

    h_nvtx_eff[phoDecBinIndex][1]->Fill(EvtInfo_nVtxGood);

    h_phopt_eff[phoDecBinIndex][1]->Fill(leadingPhotonEt);
    h_phoeta_eff[phoPtBinIndex][1]->Fill(leadingPhotonEta);


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


    // lumping all fiducial jets together
    if(findSumJet){
      h_sumjetpt_eff[phoDecBinIndex][1]->Fill(sumJetPt);
      h_sumjeteta_eff[phoDecBinIndex][1]->Fill(sumJetEta);
      h_ptratio_sumJet[phoDecBinIndex][1]->Fill(ptRatio_sumJet);
      h_dphi_sumJet[phoDecBinIndex][phoPtBinIndex][1]->Fill(gsumj_dphi);
      h_cost_sumJet[phoDecBinIndex][phoPtBinIndex][1]->Fill(gsumj_cost_comZ);
      h_pstar_sumJet[phoDecBinIndex][phoPtBinIndex][1]->Fill(gsumj_pstar_comZ);
      h_ystar_sumJet[phoDecBinIndex][phoPtBinIndex][1]->Fill(gsumj_ystar_comZ);
    }

    // lump all pt together
    h_phoeta_eff[nPtBins][1]->Fill(leadingPhotonEta);
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

    // lumping all fiducial jets and all pt together
    if(findSumJet){
      h_dphi_sumJet[phoDecBinIndex][nPtBins][1]->Fill(gsumj_dphi);
      h_cost_sumJet[phoDecBinIndex][nPtBins][1]->Fill(gsumj_cost_comZ);
      h_pstar_sumJet[phoDecBinIndex][nPtBins][1]->Fill(gsumj_pstar_comZ);
      h_ystar_sumJet[phoDecBinIndex][nPtBins][1]->Fill(gsumj_ystar_comZ);
    }

  } // end of loop over entries

  for(int i=0;i<30;i++)if(nPass[i]>0)cout << "nPass[" << i << "] = " << nPass[i] << endl;


  //---------------------------------------------------------------------------------------
  // Write these histograms to an output root file, the output will be used for efficiency
  //---------------------------------------------------------------------------------------
     
  std::string remword  ="/data4/syu/7TeV_vectorNtuple/pythia/";
  std::string remword2 ="/data4/syu/7TeV_vectorNtuple/madgraph/";

  size_t pos  = _inputFileName.find(remword);
  size_t pos2 = _inputFileName.find(remword2);

  if(pos!= std::string::npos)
    _inputFileName.swap(_inputFileName.erase(pos,remword.length()));

  else if(pos2!= std::string::npos)
    _inputFileName.swap(_inputFileName.erase(pos2,remword2.length()));

  else
    _inputFileName = "test.root";

  std::string prefix = "/home/syu/CVSCode/eff_";
  if(applyCOMCut)prefix += "withCOMCut_";
  if(applyPileUpCorr)prefix += "pileCorr_";
  TFile* outFile = new TFile(Form("%s%s",prefix.data(),_inputFileName.data()),"recreate");               

  h_pthat->Write();
  h_ngenjet->Write();
  h_nrecjet->Write();

  h_ptpho->Write();
  h_etapho->Write();

  h_genciso->Write();
  h_gentiso->Write();


  for(int ip=0; ip<4; ip++){

    h_sieie[ip]->Write();
    h_eciso[ip]->Write();
    h_hciso[ip]->Write();
    h_tkiso[ip]->Write();
    h_hovere[ip]->Write();
    h_pixel[ip]->Write();

  }

  h_ptjet->Write();
  h_etajet->Write();

  for(int idec=0; idec< nDECs; idec++){

    h2_ystar_pstar[idec]->Write();
    h_ptclosestjet[idec]->Write();
    h_etaclosestjet[idec]->Write();
    h_dRclosestjet[idec]->Write();
    h_dRLeadingPhoJet[idec]->Write();

    for(int ip=0; ip<2; ip++){
      pf_nvtx_eciso[idec][ip]->Write();
      pf_nvtx_hciso[idec][ip]->Write();
      pf_nvtx_tkiso[idec][ip]->Write();
      h_njet[idec][ip]->Write();

      h_jetpt_eff[idec][ip]->Write();
      h_jeteta_eff[idec][ip]->Write();

      h_sumjetpt_eff[idec][ip]->Write();
      h_sumjeteta_eff[idec][ip]->Write();

      h_nvtx_eff[idec][ip]->Write();
      h_phopt_eff[idec][ip]->Write();
      h_ptratio[idec][ip]->Write();
      h_ptratio_sumJet[idec][ip]->Write();
    } // end of loop over process, ip=0 before cut, 1 after cut

  } // end of loop over barrel and endcap


  for(int ipt=0; ipt< nPtBins+1; ipt++){
    for(int ip=0; ip<2; ip++){
      h_phoeta_eff[ipt][ip]->Write();
      for(int idec=0; idec < nDECs; idec++){
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

	h_dphi_sumJet[idec][ipt][ip]->Write();
	h_cost_sumJet[idec][ipt][ip]->Write();
	h_pstar_sumJet[idec][ipt][ip]->Write();
	h_ystar_sumJet[idec][ipt][ip]->Write();

      }
    }
  }


  outFile->Close();     

}

Int_t  yj_angularmc_eff::phoDecCode(Int_t ipho)
{
  double eta = PhotonScEta->at(ipho);
  if(fabs(eta) < BARREL_MAXETA)return 0;
  if(fabs(eta) > ENDCAP_MINETA && 
     fabs(eta) < ENDCAP_MAXETA)return 1;
  return -1;

}

Bool_t yj_angularmc_eff::isFidPho (Int_t ipho)
{
  if(phoDecCode(ipho)<0)return false;
  double et  = PhotonEt   ->at(ipho);  
  if(et < ptbound[0])return false;
  if(et > ptbound[nPtBins])return false;
  return true;
}



Bool_t yj_angularmc_eff::isGoodPho(Int_t ipho, bool applyPileUpCorr)
{
  double et  = PhotonEt   ->at(ipho);

  if(!isFidPho(ipho))return false;
  if(PhotonhadronicOverEm->at(ipho) > 0.05)return false;
  if(PhotonhasPixelSeed->at(ipho)   > 1e-6)return false; // this should be saved as bool

  double eciso = phoEcalIso(ipho, applyPileUpCorr); 
  double hciso = phoHcalIso(ipho, applyPileUpCorr);
  double tkiso = phoTrkIso(ipho, applyPileUpCorr);

  if(eciso > 4.2 +0.006*et)return false;
  if(hciso > 2.2 +0.0025*et)return false;
  if(tkiso > 2.0 +0.001*et)return false;

  return true;
}

Double_t  yj_angularmc_eff::phoEcalIso(Int_t ipho, bool applyPileUpCorr)
{
  if(!isFidPho(ipho))return -9999.0;

  double eciso = PhotonecalRecHitSumEtConeDR04->at(ipho);
  int decBinIndex = phoDecCode(ipho);

  if(applyPileUpCorr)
    eciso -= aeff[decBinIndex][0]*Photonrho25; 
  return eciso;
}
 
Double_t yj_angularmc_eff::phoHcalIso(Int_t ipho, bool applyPileUpCorr)
{
  if(!isFidPho(ipho))return -9999.0;

  double hciso = PhotonhcalTowerSumEtConeDR04->at(ipho);
  int decBinIndex = phoDecCode(ipho);

  if(applyPileUpCorr)
    hciso -= aeff[decBinIndex][1]*Photonrho25; 
  return hciso;

}

Double_t yj_angularmc_eff::phoTrkIso(Int_t ipho, bool applyPileUpCorr)
{
  if(!isFidPho(ipho))return -9999.0;

  double tkiso = PhotontrkSumPtHollowConeDR04->at(ipho);
  int decBinIndex = phoDecCode(ipho);

  if(applyPileUpCorr)
    tkiso -= aeff[decBinIndex][2]*Photonrho25; 
  return tkiso;

}


Bool_t yj_angularmc_eff::isFidJet (Int_t ijet)
{
  if(fabs(patJetPfAk05Eta_->at(ijet)) > 2.4)return false;
  if(patJetPfAk05Pt_->at(ijet) < 30.0)return false;
  return true;
}


// check if this reco-jet is a good loose jet
Bool_t yj_angularmc_eff::isGoodLooseJet(Int_t ijet)
{
  if(!isFidJet(ijet))return false;
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
Bool_t yj_angularmc_eff::isGoodMediumJet(Int_t ijet)
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
Bool_t yj_angularmc_eff::isGoodTightJet(Int_t ijet)
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


Int_t yj_angularmc_eff::matchRecoToGenJet(Int_t ijet)
{  
  int matchedGenIndex = -1;
  for(unsigned int k=0; k< genJetPt_->size(); k++)
    {

      double dR = eiko::deltaR(genJetEta_->at(k),genJetPhi_->at(k),
			       patJetPfAk05Eta_->at(ijet),patJetPfAk05Phi_->at(ijet)); 

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


Int_t yj_angularmc_eff::matchGenToRecoJet(Int_t ijet)
{
  int matchedRecoIndex = -1;
  for(unsigned int k=0; k< patJetPfAk05Pt_->size(); k++)
    {

      double dR = eiko::deltaR(genJetEta_->at(ijet),genJetPhi_->at(ijet),
			       patJetPfAk05Eta_->at(k),patJetPfAk05Phi_->at(k)); 

      double relPt = genJetPt_->at(ijet)>1e-6? 
	fabs(genJetPt_->at(ijet)-patJetPfAk05Pt_->at(k))/genJetPt_->at(ijet): -9999.0;

      if(dR<0.4 && relPt < 3.0)
	{
	  matchedRecoIndex = k;
	  break;
	}
      
    } // end of loop over generator-level jets
      
  return matchedRecoIndex;

}

TH1D* yj_angularmc_eff::createHisto(TH1D* htemplate, string histoName, string xtitle, int ip, double ptmin, double ptmax)
{
  TH1D* output;
  if(ptmin < 0 || ptmax < 0){
    std::string tempTitle = Form("%d < p_{T}(#gamma) < %d GeV",(int)ptbound[0],(int)ptbound[nPtBins]);	    
    output = (TH1D*)htemplate->Clone(Form("%s_allpt_%d", histoName.data(),ip));     
    output -> SetTitle(tempTitle.data());
    output -> SetXTitle(xtitle.data());
  } // if this is lumping all pt bins together

  else{
    std::string tempTitle = Form("%d < p_{T}(#gamma) < %d GeV",(int)ptmin,(int)ptmax);	    
    output = (TH1D*)htemplate->Clone(Form("%s_%d_%d_%d", histoName.data(),(int)ptmin, (int)ptmax, ip)); 
    output -> SetTitle(tempTitle.data());
    output -> SetXTitle(xtitle.data());    
  }

  return output;

}
