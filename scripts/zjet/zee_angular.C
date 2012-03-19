#define zee_angular_cxx
#include "zee_angular.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include "myLib.h"


void zee_angular::Loop()
{
  if (fChain == 0) return;

  //---------------------------------------------------------------------------------------------------------------------
  //
  //   Specify the corresponding trigger path for each run period
  // 
  //---------------------------------------------------------------------------------------------------------------------

  const int nTriggers = 11;
  string zeeTriggerNames[nTriggers];

  zeeTriggerNames[0]="HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1";
  zeeTriggerNames[1]="HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2";
  zeeTriggerNames[2]="HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3";
  zeeTriggerNames[3]="HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4";
  zeeTriggerNames[4]="HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5";
  zeeTriggerNames[5]="HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6";
  zeeTriggerNames[6]="HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6";
  zeeTriggerNames[7]="HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7";
  zeeTriggerNames[8]="HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8";
  zeeTriggerNames[9]="HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9";
  zeeTriggerNames[10]="HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10";

  int runLo[nTriggers]=
    {
      160404,
      161217,
      163270,
      165088,
      165970,
      167039,
      170722,
      170826,
      173236,
      178420,
      179959
    };
  int runHi[nTriggers]=
    {
      161176,
      163261,
      163869,
      165633,
      166967,
      167913,
      170759,
      173198,
      175875,
      179889, 
      180252 
    };


  //---------------------------------------------------------------------------------------------------------------------
  //
  //   Defining templates for each type of physics variable
  // 
  //---------------------------------------------------------------------------------------------------------------------

  TH2D* h2d_mass_template = new TH2D("h2d_mass_template","",500,0,500,500,0,500);
  h2d_mass_template->SetXTitle("M_{jj} [GeV/c^{2}]");

  TH1D* h_mass_template = new TH1D("h_mass_template","",200,20,220);
  h_mass_template->SetXTitle("M_{ee} [GeV/c^{2}]");

  TH1D* h_cost_template = new TH1D("h_cost_template","", 100,0.0,1.0);
  h_cost_template->SetXTitle("cos#theta^{*}");

  TH1D* h_ystar_template = new TH1D("h_ystar_template","",100,-5.0,5.0);
  h_ystar_template->SetXTitle("y^{*}");

  TH1D* h_pt_template    = new TH1D("h_pt_template","",500,0,500);
  TH1D* h_y_template     = new TH1D("h_y_template","",100,-5.0,5.0);
  TH1D* h_njet_template = new TH1D("h_njet_template","",21,-0.5,20.5);

  //---------------------------------------------------------------------------------------------------------------------
  //
  //   Defining histograms
  // 
  //---------------------------------------------------------------------------------------------------------------------


  TH2D* h2d_mjj_mzjj = (TH2D*)h2d_mass_template->Clone("h2d_mjj_mzjj");
  h2d_mjj_mzjj->SetYTitle("M_{Zjj} [GeV/c^{2}]");

  TH2D* h2d_mjj_mdiff = (TH2D*)h2d_mass_template->Clone("h2d_mjj_mdiff");
  h2d_mjj_mdiff->SetYTitle("M_{Zjj}-M_{jj} [GeV/c^{2}]");

  TH1D* h_zmass_raw     = (TH1D*)h_mass_template->Clone("h_zmass_raw");
  h_zmass_raw->SetTitle("Before any ID selections, after zee filter");

  TH1D* h_zmass_ID      = (TH1D*)h_mass_template->Clone("h_zmass_ID");
  h_zmass_ID->SetTitle("After WP80VBTF11 electron ID selections");

   
  // 0: inclusive, at least one jet, n>0: exclusive n jet (up to 6)
  const int NMAXJETS = 7; 

  // correlation histograms
  TH1D* h_cost_COM3D[NMAXJETS];
  TH1D* h_cost_COMZ[NMAXJETS];
  TH1D* h_cost_sumJet[NMAXJETS];
  TH1D* h_ystar[NMAXJETS];
  TH1D* h_ystar_COM3D[NMAXJETS];
  TH1D* h_ystar_COMZ[NMAXJETS];
  TH1D* h_ystar_sumJet[NMAXJETS];
  
  TH1D* h_zpt[NMAXJETS];
  TH1D* h_zy[NMAXJETS];
  
  TH1D* h_leadingjetpt[NMAXJETS];
  TH1D* h_leadingjety[NMAXJETS];

  TH1D* h_sumjetpt[NMAXJETS];
  TH1D* h_sumjety[NMAXJETS];

  string jetprefix;
  string threedprefix;
  string onedprefix;
  string sumjetprefix;

  for(int ij=0; ij < NMAXJETS; ij++)
    {
      if(ij==0)jetprefix = "#geq 1 jets";
      else jetprefix = Form("= %d jets",ij);
     
      threedprefix = jetprefix + ", boosted to the COM frame in 3D";
      onedprefix   = jetprefix + ", boosted to the COM frame in z-direction";
      sumjetprefix = jetprefix + ", sum all fiducial jets, boosted to the COM frame in z-direction";

      h_cost_COM3D[ij]  = (TH1D*)h_cost_template->Clone(Form("h_cost_COM3D_%d",ij));
      h_cost_COMZ [ij]  = (TH1D*)h_cost_template->Clone(Form("h_cost_COMZ_%d",ij));
      h_cost_sumJet[ij] = (TH1D*)h_cost_template->Clone(Form("h_cost_sumJet_%d",ij));

      h_ystar      [ij]  = (TH1D*)h_ystar_template->Clone(Form("h_ystar_%d",ij));
      h_ystar      [ij]  -> SetXTitle("0.5(y_{Z}-y_{jet^{1st}})");

      h_ystar_COM3D[ij]  = (TH1D*)h_ystar_template->Clone(Form("h_ystar_COM3D_%d",ij));
      h_ystar_COMZ [ij]  = (TH1D*)h_ystar_template->Clone(Form("h_ystar_COMZ_%d",ij));
      h_ystar_sumJet[ij] = (TH1D*)h_ystar_template->Clone(Form("h_ystar_sumJet_%d",ij));

      h_zpt[ij]          = (TH1D*)h_pt_template->Clone(Form("h_zpt_%d",ij));
      h_zpt[ij]          -> SetXTitle("p_{T}(Z) [GeV]");

      h_zy[ij]           = (TH1D*)h_y_template->Clone(Form("h_zy_%d",ij));
      h_zy[ij]           -> SetXTitle("y(Z)");

      h_leadingjetpt[ij]          = (TH1D*)h_pt_template->Clone(Form("h_leadingjetpt_%d",ij));
      h_leadingjetpt[ij]          -> SetXTitle("p_{T}(jet^{1st}) [GeV]");

      h_leadingjety[ij]           = (TH1D*)h_y_template->Clone(Form("h_leadingjety_%d",ij));
      h_leadingjety[ij]           -> SetXTitle("y(jet^{1st})");

      h_sumjetpt[ij]              = (TH1D*)h_pt_template->Clone(Form("h_sumjetpt_%d",ij));;
      h_sumjetpt[ij]              -> SetXTitle("p_{T}(jet^{sum}) [GeV]");

      h_sumjety[ij]               = (TH1D*)h_y_template->Clone(Form("h_sumjety_%d",ij));
      h_sumjety[ij]               -> SetXTitle("y(jet^{sum})");

      h_cost_COM3D[ij]   -> SetTitle(threedprefix.data());
      h_ystar_COM3D[ij]  -> SetTitle(threedprefix.data());

      h_cost_COMZ [ij]   -> SetTitle(onedprefix.data());
      h_ystar_COMZ [ij]  -> SetTitle(onedprefix.data());

      h_cost_sumJet[ij]  -> SetTitle(sumjetprefix.data());
      h_ystar_sumJet[ij] -> SetTitle(sumjetprefix.data());

      h_ystar      [ij]  -> SetTitle(jetprefix.data());
      h_zpt[ij]          -> SetTitle(jetprefix.data());
      h_zy[ij]           -> SetTitle(jetprefix.data());

      h_leadingjetpt[ij] -> SetTitle(jetprefix.data());
      h_leadingjety[ij]  -> SetTitle(jetprefix.data());
      h_sumjetpt[ij]     -> SetTitle(jetprefix.data());
      h_sumjety[ij]      -> SetTitle(jetprefix.data());


      
    }


  TH1D* h_njet = (TH1D*)h_njet_template->Clone("h_njet");
  h_njet->SetXTitle("Number of good loose PF jets");


  Long64_t nPass[50]={0};
  
  Long64_t nentries = fChain->GetEntriesFast();

  cout << "There are " << nentries << " entries." << endl;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
//     if(jentry >1000 ) break;

    nPass[0]++;

    //---------------------------------------------------------------------------------------------------------------------
    //
    //   Apply trigger selections
    // 
    //---------------------------------------------------------------------------------------------------------------------
    
    bool fireTrigger = false;
    Int_t run = EvtInfo_RunNum;

    for(int iPath=0; iPath < nTriggers; iPath++)
      {
	if(run < runLo[iPath])continue;
	if(run > runHi[iPath])continue;	  

	for(unsigned int itrig=0; itrig < trigName->size(); itrig++)
	  {
	    if(trigName->at(itrig).find(zeeTriggerNames[iPath].data())!= string::npos 
	       && trigResults->at(itrig) ==1)
	      {
		fireTrigger = true;
		break;
	      }

	  } // end of loop over names
      } // end of loop over trigger paths

    // if this is a MC sample
    if(run == 1)fireTrigger = true;

    if(!fireTrigger)continue;
    nPass[1]++;


    
    //---------------------------------------------------------------------------------------------------------------------
    //
    //   Look for a Z with two good electron legs
    // 
    //---------------------------------------------------------------------------------------------------------------------
    

    bool findAZ = false;
    TLorentzVector l4_ele1(0,0,0,0);
    TLorentzVector l4_ele2(0,0,0,0);

    for(unsigned int iele=0; iele < patElecEnergy_->size(); iele++){

      if(!isWP80VBTF11Ele(iele))continue;

      for(unsigned int jele=0; jele < iele; jele++){

 	if(!isWP80VBTF11Ele(jele))continue;


	l4_ele1.SetPtEtaPhiE(patElecPt_->at(iele),
			     patElecEta_->at(iele),
			     patElecPhi_->at(iele),
			     patElecEnergy_->at(iele));

	l4_ele2.SetPtEtaPhiE(patElecPt_->at(jele),
			     patElecEta_->at(jele),
			     patElecPhi_->at(jele),
			     patElecEnergy_->at(jele));


 	findAZ = true;

	break;

  	// require oppositely charged?
	
//   	double product_charge = patElecCharge_->at(iele)
//   	  *patElecCharge_->at(jele);

//   	// if one of the electron has zero charge or the pair is same signed
//   	if(fabs(product_charge)<1e-6 || product_charge > 1e-3)
//   	  continue;

//   	h_zmass_ID_oppQ->Fill(mZ);

      
      } // end of loop over second electron
    } // end of loop over first  electron

    
    if(!findAZ)continue;
    nPass[2]++;
    
    
    //---------------------------------------------------------------------------------------------------------------------
    //
    //   Look for a loose PF jet and find the leading jet
    // 
    //---------------------------------------------------------------------------------------------------------------------
    TLorentzVector l4_sumjet(0,0,0,0);
    int nGoodLooseJets = 0; // number of jets passing loose jet PF ID
    bool findJetComponent = false;
    int leadingJetIndex = -1;
    double maxJetPt = -9999.0;
    int secondLeadingJetIndex = -1;
    double secondJetMaxPt = -9999.0;

    for(unsigned int ijet=0; ijet < patJetPfAk05Pt_->size(); ijet++){
      
      if(!isGoodLooseJet(ijet))continue;
    
      TLorentzVector l4_thisJet(0,0,0,0);
      l4_thisJet.SetPtEtaPhiE(
  			      patJetPfAk05Pt_->at(ijet),
  			      patJetPfAk05Eta_->at(ijet),
  			      patJetPfAk05Phi_->at(ijet),
  			      patJetPfAk05En_->at(ijet)
  			      );
      // avoid overlap with electrons
      if(!eiko::separated(l4_thisJet,l4_ele1,0.3))continue;
      if(!eiko::separated(l4_thisJet,l4_ele2,0.3))continue;
    
      nGoodLooseJets++;
     
      l4_sumjet += l4_thisJet;
      findJetComponent = true;

      double thisJetPt = patJetPfAk05Pt_->at(ijet);

      // find the highest pt jet
	if(thisJetPt > maxJetPt)
	  {
	    secondJetMaxPt = maxJetPt;
	    secondLeadingJetIndex = leadingJetIndex;

	    maxJetPt = thisJetPt;
	    leadingJetIndex= ijet;
	  }
	else if(thisJetPt > secondJetMaxPt)
	  {
	    secondJetMaxPt = thisJetPt;
	    secondLeadingJetIndex = ijet;
	  }


    } // end of loop over reconstructed jets

    h_njet->Fill(nGoodLooseJets);

    if(leadingJetIndex < 0)continue;
    nPass[3]++;

    if(nGoodLooseJets < 1)continue; 
    nPass[4]++;

    TLorentzVector l4_1stjet(0,0,0,0);
    l4_1stjet.SetPtEtaPhiE(
  			   patJetPfAk05Pt_->at(leadingJetIndex),
  			   patJetPfAk05Eta_->at(leadingJetIndex),
  			   patJetPfAk05Phi_->at(leadingJetIndex),
  			   patJetPfAk05En_->at(leadingJetIndex)
  			   );


    
    
    //---------------------------------------------------------------------------------------------------------------------
    //
    //   start making angular histogram for at least one jet case
    // 
    //--------------------------------------------------------------------------------------------------------------------- 

    double mZ = (l4_ele1+l4_ele2).M();


    h_zmass_ID->Fill(mZ);

    if(mZ < 76 || mZ > 106)continue;

    TLorentzVector l4_z = l4_ele1+l4_ele2;

    double zpt = l4_z.Pt();
    h_zpt[0]->Fill(zpt);

    double zy  = l4_z.Rapidity();
    h_zy[0]->Fill(zy);

    double jetpt = l4_1stjet.Pt();
    h_leadingjetpt[0]->Fill(jetpt);

    double jety  = l4_1stjet.Rapidity();
    h_leadingjety[0]->Fill(jety);


    double sumjetpt = l4_sumjet.Pt();
    h_sumjetpt[0]->Fill(sumjetpt);

    double sumjety  = l4_sumjet.Rapidity();
    h_sumjety[0]->Fill(sumjety);

    double zj_cost_com3D = eiko::cosThetaStar_BoostToCM(l4_z, l4_1stjet);
    h_cost_COM3D[0]      ->Fill(zj_cost_com3D);

    double zj_cost_comZ  = eiko::cosThetaStar_ZBoostToCM(l4_z, l4_1stjet);
    h_cost_COMZ[0]       ->Fill(zj_cost_comZ);

    double zsumj_cost_comZ  = eiko::cosThetaStar_ZBoostToCM(l4_z, l4_sumjet);
    h_cost_sumJet[0]     ->Fill(zsumj_cost_comZ);

    double zj_ystar       = eiko::ystar(l4_z, l4_1stjet);
    h_ystar[0]           ->Fill(zj_ystar);

    double zj_ystar_com3D = eiko::ystar_BoostToCM(l4_z, l4_1stjet);
    h_ystar_COM3D[0]     ->Fill(zj_ystar_com3D);

    double zj_ystar_comZ  = eiko::ystar_ZBoostToCM(l4_z, l4_1stjet);
    h_ystar_COMZ[0]      ->Fill(zj_ystar_comZ);

    double zsumj_ystar_comZ = eiko::ystar_ZBoostToCM(l4_z, l4_sumjet);
    h_ystar_sumJet[0]    ->Fill(zsumj_ystar_comZ);

	

    //---------------------------------------------------------------------------------------------------------------------
    //
    //   start making angular histogram for exclusive n jet case, up to NMAXJETS
    // 
    //--------------------------------------------------------------------------------------------------------------------- 
    if(nGoodLooseJets > NMAXJETS-1) continue; // don't want array to go out of bound


    h_zpt[nGoodLooseJets]->Fill(zpt);
    h_zy[nGoodLooseJets]->Fill(zy);
    h_leadingjetpt[nGoodLooseJets]->Fill(jetpt);
    h_leadingjety[nGoodLooseJets]->Fill(jety);
    h_sumjetpt[nGoodLooseJets]->Fill(sumjetpt);
    h_sumjety[nGoodLooseJets]->Fill(sumjety);

    h_cost_COM3D[nGoodLooseJets]      ->Fill(zj_cost_com3D);
    h_cost_COMZ[nGoodLooseJets]       ->Fill(zj_cost_comZ);
    h_cost_sumJet[nGoodLooseJets]     ->Fill(zsumj_cost_comZ);
    h_ystar[nGoodLooseJets]           ->Fill(zj_ystar);
    h_ystar_COM3D[nGoodLooseJets]     ->Fill(zj_ystar_com3D);
    h_ystar_COMZ[nGoodLooseJets]      ->Fill(zj_ystar_comZ);
    h_ystar_sumJet[nGoodLooseJets]    ->Fill(zsumj_ystar_comZ);

    if(secondLeadingJetIndex < 0)continue;


    TLorentzVector l4_2ndjet(0,0,0,0);
    l4_2ndjet.SetPtEtaPhiE(
  			   patJetPfAk05Pt_->at(secondLeadingJetIndex),
  			   patJetPfAk05Eta_->at(secondLeadingJetIndex),
  			   patJetPfAk05Phi_->at(secondLeadingJetIndex),
  			   patJetPfAk05En_->at(secondLeadingJetIndex)
  			   );


    
    double mjj = (l4_1stjet + l4_2ndjet).M();
    double mzjj = (l4_z + l4_1stjet + l4_2ndjet).M();

    h2d_mjj_mzjj->Fill(mjj,mzjj);
    h2d_mjj_mdiff->Fill(mjj,mzjj-mjj);


  } // end of loop over entries


  std::string remword  ="/data4/syu/Zjet_DelPanjNtuple/CTRLPTrigON/";

  size_t pos  = _inputFile.find(remword);

  if(pos!= std::string::npos)
    _inputFile.swap(_inputFile.erase(pos,remword.length()));


  std::string dirName = "zjet_histos";
  gSystem->mkdir(dirName.data());
  
  std::string prefix = dirName + "/histo_";

  TFile* outFile = new TFile(Form("%s%s",prefix.data(),_inputFile.data()),"recreate");               

  h2d_mjj_mzjj->Write();
  h2d_mjj_mdiff->Write();
  h_zmass_raw->Write();
  h_zmass_ID->Write();

  for(int ij=0; ij< NMAXJETS; ij++){

    h_cost_COM3D[ij]->Write();
    h_cost_COMZ[ij]->Write();
    h_cost_sumJet[ij]->Write();
    h_ystar[ij]->Write();
    h_ystar_COM3D[ij]->Write();
    h_ystar_COMZ[ij]->Write();
    h_ystar_sumJet[ij]->Write();
  
    h_zpt[ij]->Write();
    h_zy[ij]->Write();
  
    h_leadingjetpt[ij]->Write();
    h_leadingjety[ij]->Write();
    h_sumjetpt[ij]->Write();
    h_sumjety[ij]->Write();

  }

  h_njet->Write();
  outFile->Close();     
 
  for(int i=0;i<50;i++)if(nPass[i]>0)cout << "nPass[" << i << "] = " << nPass[i] << endl;
  

}


Bool_t zee_angular::isFidEle (Int_t iele)
{
  double pt = patElecPt_->at(iele);
  if(pt < 20.0)return false;
  
  double eta = patElecEta_->at(iele); // should use supercluster eta, but this variable is not filled
  bool ingap= fabs(eta)>1.446 && fabs(eta)< 1.566;
  if(ingap)return false;

  if(fabs(eta) > 2.4)return false;

  return true;
}

Bool_t zee_angular::isWP80VBTF11Ele(Int_t iele)
{
  if(!isFidEle(iele))return false;

  // // should have applied electron ID, but not the ID variables are not saved in Lovedeep's ntuple
  //   if(patElecMissingHits_->at(iele)>0)return false;

  //   if(fabs(patElecDist_->at(iele)) < 0.02)return false;

  //   if(fabs(patElecDeltaCotTheta_->at(iele)) < 0.02)return false;

  //   bool isEB = fabs(patElecInBarrel_->at(iele)-1.0)<1e-6;
  //   bool isEE = fabs(patElecInEndcap_->at(iele)-1.0)<1e-6;

  //   double sieie = patElecSigIhIh_->at(iele);
  //   double dphi_in  = fabs(patElecDelPhiIn_->at(iele));
  //   double deta_in  = fabs(patElecDelEtaIn_->at(iele));
  //   double hadem = patElecHoE_->at(iele);

  //   // shower shape
  //   // barrel
  //   if( isEB && sieie > 0.01)return false;
  //   if( isEB && dphi_in > 0.06)return false;
  //   if( isEB && deta_in > 0.004)return false;
  //   if( isEB && hadem > 0.04)return false;
  
  //   // endcap
  //   if( isEE && sieie > 0.03)return false;
  //   if( isEE && dphi_in > 0.03)return false;
  //   if( isEE && deta_in > 0.007)return false;
  //   if( isEE && hadem > 0.15)return false;
  
  // pf isolation
  //   double pfiso = patElecChHadIso_->at(iele) +
  //     patElecNeHadIso_->at(iele) + patElecGamIso_->at(iele);
  //   double pt = patElecPt_->at(iele);
  
  //   double relative_pfiso = pfiso/pt;

  //   if( relative_pfiso > 0.2)return false;
   
  
  return true;
}


Bool_t zee_angular::isFidJet (Int_t ijet)
{
  if(fabs(patJetPfAk05Eta_->at(ijet)) > 2.4)return false;
  if(patJetPfAk05Pt_->at(ijet) < 30.0)return false;
  return true;
}


// check if this reco-jet is a good loose jet
Bool_t zee_angular::isGoodLooseJet(Int_t ijet)
{
  if(!isFidJet(ijet))return false;
  //   if(patJetPfAk05JetNConstituents_->at(ijet) <= 1)return false;
  if(patJetPfAk05NeutHadEFr_->at(ijet) >= 0.99)return false;
  if(patJetPfAk05NeutEmEFr_->at(ijet) >= 0.99)return false;

  //   // for the tracker region
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharHadEFr_->at(ijet) <= 0.0)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharEmEFr_->at(ijet) >= 0.99)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharMulti_->at(ijet) <= 0)return false;

  return true;

}

