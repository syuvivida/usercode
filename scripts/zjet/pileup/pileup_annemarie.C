#define pileup_annemarie_cxx
#include "pileup_annemarie.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include "myLib.h"


const double minZPt  =40.0;

const double minJetPt=30.0;
const double maxJetEta=2.4;
const double mindR=0.5;

const double minElePt = 20.0;
const double minEleBarrelEta = 0.0;
const double maxEleBarrelEta = 1.442;

const double minEleEndcapEta = 1.566;
const double maxEleEndcapEta = 2.1;
const double minMee = 76.0;
const double maxMee =106.0;

void pileup_annemarie::Loop(bool match)
{

  if (fChain == 0) return;

  //---------------------------------------------------------------------------------------------------------------------
  //
  //   Defining templates for each type of physics variable
  // 
  //---------------------------------------------------------------------------------------------------------------------

  TH1D* h_nvtx_template = new TH1D("h_nvtx_template","",31,-0.5,30.5);
  h_nvtx_template->Sumw2();

  TH1D* h_mass_template = new TH1D("h_mass_template","",500, 0, 500);
  h_mass_template->SetXTitle("M_{ee} [GeV/c^{2}]");
  h_mass_template->Sumw2();

  TH1D* h_pt_template    = new TH1D("h_pt_template","", 40,0,400);
  h_pt_template->Sumw2();

  TH1D* h_y_template     = new TH1D("h_y_template","",15,0.0,3.0);
  h_y_template->Sumw2();

  TH1D* h_njet_template = new TH1D("h_njet_template","",21,-0.5,20.5);
  h_njet_template->Sumw2();

  //---------------------------------------------------------------------------
  //
  //   Defining histograms
  // 
  //---------------------------------------------------------------------------

  TH1D* h_nvtx_before = (TH1D*)h_nvtx_template->Clone("h_nvtx_before");
  h_nvtx_before->SetXTitle("Number of vertices");
  h_nvtx_before->SetTitle("Before pile-up reweighting");

  TH1D* h_nvtx_after = (TH1D*)h_nvtx_template->Clone("h_nvtx_after");
  h_nvtx_after->SetXTitle("Number of vertices");
  h_nvtx_after->SetTitle("After pile-up reweighting");

  TH1D* h_zmass_ID      = (TH1D*)h_mass_template->Clone("h_zmass_ID");
  h_zmass_ID->SetTitle("After WP80VBTF11 electron ID selections");
   
  // 0: inclusive, at least one jet, n>0: exclusive n jet (up to 6)
  const int NMAXJETS = 7; 

  // correlation histograms

  TH1D* h_yB[NMAXJETS];
  TH1D* h_ystar[NMAXJETS];  
  TH1D* h_zpt[NMAXJETS];
  TH1D* h_zy[NMAXJETS];
  TH1D* h_jetpt[NMAXJETS];
  TH1D* h_jety[NMAXJETS];

  string jetprefix;

  for(int ij=0; ij < NMAXJETS; ij++)
    {
      if(ij==0)jetprefix = "#geq 1 jets";
      else if(ij==1)jetprefix = Form("= %d jet",ij);
      else jetprefix = Form("= %d jets",ij);
     
      h_ystar      [ij]  = (TH1D*)h_y_template->Clone(Form("h_ystar_%d",ij));
      h_ystar      [ij]  -> SetXTitle("0.5|Y_{Z}-Y_{jet}|");

      h_yB         [ij]  = (TH1D*)h_y_template->Clone(Form("h_yB_%d",ij));
      h_yB         [ij]  -> SetXTitle("0.5|Y_{Z}+Y_{jet}|");

      h_zpt[ij]          = (TH1D*)h_pt_template->Clone(Form("h_zpt_%d",ij));
      h_zpt[ij]          -> SetXTitle("p_{T}(Z) [GeV/c]");

      h_zy[ij]           = (TH1D*)h_y_template->Clone(Form("h_zy_%d",ij));
      h_zy[ij]           -> SetXTitle("y(Z)");

      h_jetpt[ij]          = (TH1D*)h_pt_template->Clone(Form("h_jetpt_%d",ij));
      h_jetpt[ij]          -> SetXTitle("p_{T}(jet) [GeV/c]");

      h_jety[ij]           = (TH1D*)h_y_template->Clone(Form("h_jety_%d",ij));
      h_jety[ij]           -> SetXTitle("y(jet^{1st})");

      h_yB[ij]           -> SetTitle(jetprefix.data());
      h_ystar[ij]        -> SetTitle(jetprefix.data());
      h_zpt[ij]          -> SetTitle(jetprefix.data());
      h_zy[ij]           -> SetTitle(jetprefix.data());

      h_jetpt[ij]        -> SetTitle(jetprefix.data());
      h_jety[ij]         -> SetTitle(jetprefix.data());
      
    }


  TH1D* h_njet = (TH1D*)h_njet_template->Clone("h_njet");
  h_njet->SetXTitle("Number of good loose PF jets");


  Long64_t nPass[50]={0};
  
  Long64_t nentries = fChain->GetEntriesFast();

  cout << " has " << nentries << " entries." << endl;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
//     if(jentry >1000 ) break;

    nPass[0]++;

    double eventWeight = 1;
    if(PU_weight >= 0.0)eventWeight *= PU_weight;
    if(mcWeight_>0)eventWeight *= mcWeight_;

//     cout << "eventWeight = " << eventWeight << endl;

    
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

	double mZ = (l4_ele1+l4_ele2).M();
	h_zmass_ID->Fill(mZ, eventWeight);

	if(mZ < minMee || mZ > maxMee)continue;

	double zPt = (l4_ele1+l4_ele2).Pt();
	if(zPt < minZPt)continue;

 	findAZ = true;

	break;

      
      } // end of loop over second electron
    } // end of loop over first  electron

    
    if(!findAZ)continue;
    nPass[1]++;
    
    
    //---------------------------------------------------------------------------------------------------------------------
    //
    //   Look for a loose PF jet and find the leading jet
    // 
    //---------------------------------------------------------------------------------------------------------------------
    int nGoodLooseJets = 0; // number of jets passing loose jet PF ID
    int leadingJetIndex = -1;
    double maxJetPt = -9999.0;
    
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
      if(!eiko::separated(l4_thisJet,l4_ele1,mindR))continue;
      if(!eiko::separated(l4_thisJet,l4_ele2,mindR))continue;
           
      double thisJetPt = patJetPfAk05Pt_->at(ijet);

//       int PID = abs(patJetPfAk05GenPartonID_->at(ijet));
      
//       if(PID > 5 && PID!=21)continue;

      nGoodLooseJets++;


      // find the highest pt jet
      if(thisJetPt > maxJetPt)
	{
	  maxJetPt = thisJetPt;
	  leadingJetIndex= ijet;
	}


    } // end of loop over reconstructed jets

    h_njet->Fill(nGoodLooseJets, eventWeight);

    if(leadingJetIndex < 0)continue;

    if(nGoodLooseJets < 1)continue; 
    nPass[2]++;

    TLorentzVector l4_1stjet(0,0,0,0);
    l4_1stjet.SetPtEtaPhiE(
  			   patJetPfAk05Pt_->at(leadingJetIndex),
  			   patJetPfAk05Eta_->at(leadingJetIndex),
  			   patJetPfAk05Phi_->at(leadingJetIndex),
  			   patJetPfAk05En_->at(leadingJetIndex)
  			   );

//     int partonMyIndex = matchRecoToParton(leadingJetIndex);
//     if(match && partonMyIndex < 0)continue;
//     cout << "After partonPID = " << genParId_->at(partonMyIndex) << endl;
    
    int genJetMyIndex = matchRecoToGenJet(leadingJetIndex);
    if(match && genJetMyIndex < 0)continue;

    int partonPID = abs(patJetPfAk05GenPartonID_->at(leadingJetIndex));
//     cout << "leadingJetIndex =" << leadingJetIndex << endl;
//     cout << "Before partonPID = " << partonPID << endl;

//     if(match && partonPID>5 && partonPID!=21)continue;
//     cout << "After partonPID = " << partonPID << endl;
    nPass[3]++;
    
    
    //-------------------------------------------------------------------------
    //
    //   start making angular histogram for at least one jet case
    // 
    //-------------------------------------------------------------------------



    h_nvtx_before->Fill(EvtInfo_NumVtx,1.0);
    h_nvtx_after ->Fill(EvtInfo_NumVtx,eventWeight);

    TLorentzVector l4_z = l4_ele1+l4_ele2;

    double zpt = l4_z.Pt();
    h_zpt[0]->Fill(zpt, eventWeight);

    double zy  = fabs(l4_z.Rapidity());
    h_zy[0]->Fill(zy, eventWeight);

    double jetpt = l4_1stjet.Pt();
    h_jetpt[0]->Fill(jetpt, eventWeight);

    double jety  = fabs(l4_1stjet.Rapidity());
    h_jety[0]->Fill(jety, eventWeight);

    double zj_yB          = fabs(eiko::yB(l4_z, l4_1stjet));
    h_yB[0]              ->Fill(zj_yB, eventWeight);

    double zj_ystar       = fabs(eiko::ystar(l4_z, l4_1stjet));
    h_ystar[0]           ->Fill(zj_ystar, eventWeight);


    //---------------------------------------------------------------------------------------------------------------------
    //
    //   start making angular histogram for exclusive n jet case, up to NMAXJETS
    // 
    //--------------------------------------------------------------------------------------------------------------------- 
    if(nGoodLooseJets > NMAXJETS-1) continue; // don't want array to go out of bound

    h_zpt[nGoodLooseJets]->Fill(zpt, eventWeight);
    h_zy[nGoodLooseJets]->Fill(zy, eventWeight);
    h_jetpt[nGoodLooseJets]->Fill(jetpt, eventWeight);
    h_jety[nGoodLooseJets]->Fill(jety, eventWeight);

    h_yB[nGoodLooseJets]              ->Fill(zj_yB, eventWeight);
    h_ystar[nGoodLooseJets]           ->Fill(zj_ystar, eventWeight);


  } // end of loop over entries


  std::string suffix = "nomatch";
  if(match) suffix = "yesmatch";

  TFile* outFile = new TFile(Form("anne-marie_myZPt%02i_%s.root",(int)minZPt,suffix.data()),"recreate");               

  h_nvtx_before->Write();
  h_nvtx_after->Write();

  h_zmass_ID->Write();

  for(int ij=0; ij< NMAXJETS; ij++){

    h_yB[ij]->Write();
    h_ystar[ij]->Write();
  
    h_zpt[ij]->Write();
    h_zy[ij]->Write();
  
    h_jetpt[ij]->Write();
    h_jety[ij]->Write();

  }

  h_njet->Write();
  outFile->Close();     
 
  for(int i=0;i<50;i++)if(nPass[i]>0)cout << "nPass[" << i << "] = " << nPass[i] << endl;
  

}


Bool_t pileup_annemarie::isFidEle (Int_t iele)
{
  double pt = patElecPt_->at(iele);
  if(pt < minElePt )return false;
  
  double eta = patElecEta_->at(iele); // should use supercluster eta, but this variable is not filled
  bool ingap= fabs(eta)>maxEleBarrelEta && fabs(eta)< minEleEndcapEta;
  if(ingap)return false;

  if(fabs(eta) > maxEleEndcapEta)return false;

  return true;
}

Bool_t pileup_annemarie::isWP80VBTF11Ele(Int_t iele)
{
  if(!isFidEle(iele))return false;

  // // should have applied electron ID, but not the ID variables are not saved in ntuples
  if(patElecMissingHits_->at(iele)>0)return false;

  if(fabs(patElecDist_->at(iele)) < 0.02 
     && fabs(patElecDeltaCotTheta_->at(iele)) < 0.02)
    return false;

  bool isEB = fabs(patElecInBarrel_->at(iele)-1.0)<1e-6;
  bool isEE = fabs(patElecInEndcap_->at(iele)-1.0)<1e-6;

  double sieie = patElecSigIhIh_->at(iele);
  double dphi_in  = fabs(patElecDelPhiIn_->at(iele));
  double deta_in  = fabs(patElecDelEtaIn_->at(iele));
  double hadem = patElecHoE_->at(iele);

  // shower shape
  // barrel
  if( isEB && sieie > 0.01)return false;
  if( isEB && dphi_in > 0.06)return false;
  if( isEB && deta_in > 0.004)return false;
  if( isEB && hadem > 0.04)return false;
  
     // endcap
  if( isEE && sieie > 0.03)return false;
  if( isEE && dphi_in > 0.03)return false;
  if( isEE && deta_in > 0.007)return false;
  if( isEE && hadem > 0.15)return false;
  
  // pf isolation
  double pfiso = patElecChHadIso_->at(iele) +
    patElecNeHadIso_->at(iele) + patElecGamIso_->at(iele);
  double pt = patElecPt_->at(iele);
  
  double relative_pfiso = pfiso/pt;
  
  if( relative_pfiso > 0.2)return false;
   
  
  return true;
}


Bool_t pileup_annemarie::isFidJet (Int_t ijet)
{
  if(fabs(patJetPfAk05Eta_->at(ijet)) > maxJetEta)return false;
  if(patJetPfAk05Pt_->at(ijet) < minJetPt)return false;
  return true;
}


// check if this reco-jet is a good loose jet
Bool_t pileup_annemarie::isGoodLooseJet(Int_t ijet)
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

Int_t pileup_annemarie::matchRecoToParton(Int_t ijet)
{  
  int matchedGenIndex = -1;
  for(unsigned int k=0; k< genParPt_->size(); k++)
    {
      int status    = genParSt_->at(k);
      int PID       = abs(genParId_->at(k));
      int motherPID = genMomParId_->at(k);

      if(status!=3)continue;
      if(PID>5 && PID!=21)continue;
      if(motherPID!=10002)continue;

      double dR = eiko::deltaR(genParEta_->at(k),genParPhi_->at(k),
			       patJetPfAk05Eta_->at(ijet),
			       patJetPfAk05Phi_->at(ijet)); 

      double relPt = genParPt_->at(k)>1e-6? 
	fabs(genParPt_->at(k)-patJetPfAk05Pt_->at(ijet))/genParPt_->at(k)
	: -9999.0;

      if(dR<0.4 && relPt < 3.0)
	{
	  matchedGenIndex = k;
	  break;
	}
      
    } // end of loop over generator-level information
      
  return matchedGenIndex;
}

Int_t pileup_annemarie::matchRecoToGenJet(Int_t ijet)
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

