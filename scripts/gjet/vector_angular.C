#define vector_angular_cxx
#include "vector_angular.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include "myLib.h"

const double BARREL_MAXETA=1.4442;
const double ENDCAP_MINETA=1.566;
const double ENDCAP_MAXETA=2.5;

// const double BARREL_MAXETA=1.44;
// const double ENDCAP_MINETA=1.57;
// const double ENDCAP_MAXETA=2.3;
const int nDECs = 2;

double ptbound[]={40., 45., 50., 55., 60., 65., 70., 75., 85., 95., 110., 130., 160., 200.};
const int nPtBins = sizeof(ptbound)/sizeof(ptbound[0])-1;

// effective area for correcting pileup, ECAL/HCAL/Tracker from EWK-11-251
const double aeff[2][3]={
  {0.183, 0.062, 0.0167},
  {0.090, 0.180, 0.032}
};




void vector_angular::Loop(bool onlyOneJet, bool DEBUG)
{
  if (fChain == 0) return;
  cout << "This is vector_angular version 0" << endl;
  cout << "There are " << nPtBins << " pt bins" << endl;
  cout << "require there is no jet within deltaR=0.5 of photon" << endl;
  cout << "require exclusive 1 jet: " << onlyOneJet << endl;


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
  TH1D* h_pt_template = new TH1D("h_pt_template","",500,0,500); h_pt_template->Sumw2();
  TH1D* h_njet_template = new TH1D("h_njet_template","",21,-0.5,20.5); h_njet_template->Sumw2();
  TH1D* h_ptbin_template= new TH1D("h_ptbin_template","", nPtBins, ptbound);  
  TH1D* h_yB_template_oneside = new TH1D("h_yB_template_oneside","",18,0,1.8); h_yB_template_oneside->Sumw2();
  TH1D* h_ystar_template_oneside = new TH1D("h_ystar_template_oneside","",18,0,1.8); h_ystar_template_oneside->Sumw2();

  //---------------------------------------------------------------------------------------------------------------------
  //
  //   Declaring and defining the histograms to be plotted
  // 
  //---------------------------------------------------------------------------------------------------------------------

  string decName[4];
  decName[0] = "EB";
  decName[1] = "EE";

  // debugging histograms
  TH1D* h_pthat       = (TH1D*)h_pt_template->Clone("h_pthat");
  h_pthat   -> SetXTitle("#hat{p_{T}} of the MC sample [GeV]");
  h_pthat   -> SetTitle("Before applying any selections");
  h_pthat   -> Sumw2();

  TH1D* h_ngenjet     = (TH1D*)h_njet_template->Clone("h_ngenjet");
  h_ngenjet -> SetXTitle("Number of generator-level jets");
  h_ngenjet -> SetTitle("Before applying any selections");
  h_ngenjet -> Sumw2();

  TH1D* h_nrecjet     = (TH1D*)h_njet_template->Clone("h_nrecjet");
  h_nrecjet -> SetXTitle("Number of reconstruction-level jets");
  h_nrecjet -> SetTitle("Before applying any selections");
  h_nrecjet -> Sumw2();

// 0: generator-level, 
// 1: generator-level, after requiring a gen jet, a gen photon, and a reco-photon, 
// 2: generator-level, after requiring a gen jet, a gen photon, a reco-photon, and a reco-jet
// 3: generator-level, after requiring also the ID cuts
// 4: reconstruction-level distribution after only requiring reco-objects
// 5: reconstruction-level distribution after requiring reco-objects and ID

  const int NPROCS = 6;
  TH1D* h_genystar_oneside[nDECs][NPROCS]; 
  TH1D* h_genyB_oneside[nDECs][NPROCS];


  for(int ip=0; ip< NPROCS; ip++){

    for(int idec=0; idec< nDECs; idec++)
      {
	h_genystar_oneside[idec][ip] = (TH1D*)h_ystar_template_oneside->Clone(Form("h_genystar_%s_%d", decName[idec].data(),ip));
	h_genystar_oneside[idec][ip] -> SetXTitle("Generator-level 0.5 | y^{#gamma} - y^{1stjet}| ");
        h_genystar_oneside[idec][ip] -> Sumw2();

	h_genyB_oneside[idec][ip] = (TH1D*)h_yB_template_oneside->Clone(Form("h_genyB_%s_%d", decName[idec].data(),ip));
	h_genyB_oneside[idec][ip] -> SetXTitle("Generator-level 0.5|y^{#gamma} + y^{1stjet}|");
	h_genyB_oneside[idec][ip] -> Sumw2();
      } // end of loop over barrel and endcap

  }
  
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
    //    if(jentry > 5000 ) break;
    // weight
    double event_weight = 1;
    if(PU_weight >0) event_weight *= PU_weight;
    if(mcWeight_ >0) event_weight *= mcWeight_;

    nPass[0]++;
    h_pthat->Fill(patPhotonMCpthat, event_weight); // make sure we are checking the right MC samples

    //-------------------------------------------------------------------------
    // generator-level information
    //-------------------------------------------------------------------------


    //-------------------------------------------------------------------------
    // now find generator-level photon
    //-------------------------------------------------------------------------
    int genPhotonIndex = -1;
    for(unsigned int ipart=0; ipart < genParE_->size(); ipart++)
      {
        if(genParSt_->at(ipart)!=1)continue;
        if(genParId_->at(ipart)!=22)continue;
        int momId = genMomParId_->at(ipart);
        if(isDirPho && momId!=22)continue;
        if(isFraPho && momId > 21)continue;
        if(!gen_isFidPho(ipart))continue;
        if(genPhotonIndex< 0)
	  genPhotonIndex = ipart;
        else break;
      }

    if(genPhotonIndex < 0)continue; 
    nPass[2]++;
    TLorentzVector l4_genpho(0,0,0,0);
   
    l4_genpho.SetPtEtaPhiE(genParPt_->at(genPhotonIndex),
                           genParEta_->at(genPhotonIndex),
                           genParPhi_->at(genPhotonIndex),
                           genParE_->at(genPhotonIndex)
			   );


    //-------------------------------------------------------------------------
    // find the gen jet index that has the highest pt and also find the matched
    // reco jet, this gen jet needs to have pt > 30 GeV and |eta| < 3.0, delta
    // R < 0.5 with respect to gen photon
    //-------------------------------------------------------------------------


    double genJetMaxPt     = -9999.0;
    int leadingGenJetIndex = -1;
    int nGenJet = 0;

    // loop over generator-level jets
    for(unsigned int ijet=0; ijet < genJetPt_->size(); ijet++)
      {
        if(! gen_isFidJet(ijet))continue; 
	TLorentzVector l4_thisGenJet(0,0,0,0);
	l4_thisGenJet.SetPtEtaPhiE(
				   genJetPt_->at(ijet),
				   genJetEta_->at(ijet),
				   genJetPhi_->at(ijet),
				   genJetE_->at(ijet)
				   );
 	if(! eiko::separated(l4_genpho, l4_thisGenJet))continue;
        
        nGenJet++;

    	// try to find the genJet that the highest gen pt
	double thisJetPt = genJetPt_->at(ijet);
	if(thisJetPt > genJetMaxPt)
	  {
	    genJetMaxPt                   = thisJetPt;
	    leadingGenJetIndex            = ijet;
	  }


      } // end of loop over genJets
    h_ngenjet->Fill(nGenJet, event_weight);
    if(nGenJet < 1)continue;
    nPass[3]++;
    
    if(leadingGenJetIndex <0)continue; // couldn't find a gen jet with pt > 10 GeV and |eta| < 3.0
    nPass[4]++;

    if(onlyOneJet && isDirPho && nGenJet !=1)continue;
    if(onlyOneJet && isFraPho && nGenJet !=2)continue;
    nPass[5]++;
    

    // now plot generator-level ystar and yb

    TLorentzVector l4_genjet(0,0,0,0);
    l4_genjet.SetPtEtaPhiE(
			   genJetPt_->at(leadingGenJetIndex),
			   genJetEta_->at(leadingGenJetIndex),
			   genJetPhi_->at(leadingGenJetIndex),
			   genJetE_->at(leadingGenJetIndex)
			   );
    
    double genystar = eiko::ystar(l4_genpho,l4_genjet);
    double genyB    = eiko::yB(l4_genpho,l4_genjet);
  
    int gen_phoDecBinIndex = gen_phoDecCode(genPhotonIndex);
    if(DEBUG && gen_phoDecBinIndex <0) cout << "gen_phoDecBinIndex = " <<  gen_phoDecBinIndex << endl;

    h_genystar_oneside[gen_phoDecBinIndex][0]->Fill(fabs(genystar), event_weight);
    h_genyB_oneside[gen_phoDecBinIndex][0]->Fill(fabs(genyB), event_weight);
    

    //-------------------------------------------------------------------------
    // now find reconstructed information
    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    // now find a good leading photon that is matched to hard scattering photon
    //-------------------------------------------------------------------------

    int leadingPhotonIndex = -1;
    double phoMaxPt = -9999.;

    for(unsigned int ipho=0; ipho < patPhotonPt->size(); ipho++)
      {	
	// first require this reconstruction-level photon to be fiducial
        // and has pt within the boundary we are studying
	if(!rec_isFidPho(ipho))continue;

	// then check if this photon is matched to a gen photon
 	if(patPhotonisGenMatched->at(ipho)!=1)continue;
	
	// require it to match to a prompt photon
 	if(fabs(patPhotongenMomId->at(ipho)-22)>1e-6 && isDirPho)continue; // not prompt photon 
 	if(fabs(patPhotongenMomId->at(ipho))>21.1 && isFraPho)continue; // not prompt photon 
//  	cout << "gen mom Id = " << patPhotongenMomId->at(ipho) << endl;
// 	cout << "is a prompt photon" << endl;

	// find the leading photon, usually there's only one in the MC
	double thisPhoPt= patPhotonEt->at(ipho);      
	if(thisPhoPt > phoMaxPt)
	  {
	    phoMaxPt = thisPhoPt;
	    leadingPhotonIndex= ipho;
	  }

      }
    // end of leading photon search
    if(leadingPhotonIndex<0)continue;
    nPass[6]++;

    if(leadingPhotonIndex>0)nPass[7]++;

    double leadingPhotonEt = patPhotonEt->at(leadingPhotonIndex);
    int phoPtBinIndex =  h_ptbin_template->GetXaxis()->FindBin(leadingPhotonEt)-1; 

    int phoDecBinIndex = rec_phoDecCode(leadingPhotonIndex);

    // dummy proof, but already make the same requirement in isFidPho
    if(phoPtBinIndex < 0 || phoPtBinIndex > (nPtBins-1))
      {
        if(DEBUG)
          cout << "phoPtBinIndex = " << phoPtBinIndex << endl;
        continue;
      }
    nPass[1]++;
    if(phoDecBinIndex < 0)
      {
        if(DEBUG)
          cout << "phoDecBinIndex = " << phoDecBinIndex << endl;
        continue;
      }
    nPass[8]++;

    // assign 4-vector of leading photon
    TLorentzVector l4_pho(0,0,0,0);
    l4_pho.SetPtEtaPhiE(
			patPhotonEt->at(leadingPhotonIndex),
			patPhotonEta->at(leadingPhotonIndex),
			patPhotonPhi->at(leadingPhotonIndex),
			patPhotonEnergy->at(leadingPhotonIndex)
			);

    h_genystar_oneside[gen_phoDecBinIndex][1]->Fill(fabs(genystar), event_weight);
    h_genyB_oneside[gen_phoDecBinIndex][1]->Fill(fabs(genyB), event_weight);    

    //-------------------------------------------------------------------------
    // now loop over all reco jets and find the leading jet without matching 
    //-------------------------------------------------------------------------

    int nFiducialRecoJets = 0; // number of jets passing pt > 30 GeV, |eta| < 2.4, and matched to gen jet
    for(unsigned int ijet=0; ijet < patJetPfAk05Pt_->size(); ijet++){

      if(!rec_isFidJet(ijet))continue;
      //      int genJetIndex = matchRecoToGenJet(ijet);
      //      if(genJetIndex < 0)continue; // is this needed

      TLorentzVector l4_thisJet(0,0,0,0);
      l4_thisJet.SetPtEtaPhiE(
			      patJetPfAk05Pt_->at(ijet),
			      patJetPfAk05Eta_->at(ijet),
			      patJetPfAk05Phi_->at(ijet),
			      patJetPfAk05En_->at(ijet)
			      );
      
      if(! eiko::separated(l4_pho,l4_thisJet))continue;
      nFiducialRecoJets++;

      
    } // end of loop over reconstructed jets

    h_nrecjet->Fill(nFiducialRecoJets, event_weight);

    // find the reconstruction-level jet that is matched to this highest pt genJet
    int leadingJetIndex = matchGenToRecoJet(leadingGenJetIndex);
    if(leadingPhotonIndex<0 || leadingJetIndex<0)continue;
    nPass[9]++;
    // need to have a leading photon and a leading jet
    if( ! rec_isFidJet(leadingJetIndex) )continue;

    
    // assign 4-vector of leading jet
    TLorentzVector l4_1stjet(0,0,0,0);
    l4_1stjet.SetPtEtaPhiE(
			   patJetPfAk05Pt_->at(leadingJetIndex),
			   patJetPfAk05Eta_->at(leadingJetIndex),
			   patJetPfAk05Phi_->at(leadingJetIndex),
			   patJetPfAk05En_->at(leadingJetIndex)
			   );

    if(!eiko::separated(l4_pho, l4_1stjet))continue;

    nPass[10]++;

    h_genystar_oneside[gen_phoDecBinIndex][2]->Fill(fabs(genystar), event_weight);
    h_genyB_oneside[gen_phoDecBinIndex][2]->Fill(fabs(genyB), event_weight);
    

    double gj_yB     = eiko::yB(l4_pho, l4_1stjet);    
    double gj_ystar   = eiko::ystar(l4_pho, l4_1stjet);
    h_genystar_oneside[phoDecBinIndex][4]->Fill(fabs(gj_ystar), event_weight);
    h_genyB_oneside[phoDecBinIndex][4]->Fill(fabs(gj_yB), event_weight);

    // apply ID cuts

    if(!isGoodPho(leadingPhotonIndex))continue;  
    if(!isGoodLooseJet(leadingJetIndex))continue;
    //=============================================================
    nPass[11]++;

    h_genystar_oneside[gen_phoDecBinIndex][3]->Fill(fabs(genystar), event_weight);
    h_genyB_oneside[gen_phoDecBinIndex][3]->Fill(fabs(genyB), event_weight);
    
    h_genystar_oneside[phoDecBinIndex][5]->Fill(fabs(gj_ystar), event_weight);
    h_genyB_oneside[phoDecBinIndex][5]->Fill(fabs(gj_yB), event_weight);


  } // end of loop over entries

  for(int i=0;i<30;i++)if(nPass[i]>0)cout << "nPass[" << i << "] = " << nPass[i] << endl;


  //---------------------------------------------------------------------------------------
  // Write these histograms to an output root file, the output will be used for efficiency
  //---------------------------------------------------------------------------------------
     
  std::string remword  ="/data4/syu/Gjet_vectorNtuple/";
  size_t pos  = _inputFileName.find(remword);

  if(pos!= std::string::npos)
    _inputFileName.swap(_inputFileName.erase(pos,remword.length()));
  else
    _inputFileName = "test.root";

  std::string prefix = "/home/syu/CVSCode/eff_";
  if(onlyOneJet)     prefix += "exclusiveOneJet_";
  TFile* outFile = new TFile(Form("%s%s",prefix.data(),_inputFileName.data()),"recreate");               

  h_pthat->Write();
  h_ngenjet->Write();
  h_nrecjet->Write();

  for(int idec=0; idec< nDECs; idec++){
    for(int ip=0; ip< NPROCS; ip++){
      h_genystar_oneside[idec][ip]->Write();
      h_genyB_oneside[idec][ip]->Write();
    } // end of loop over process, ip=0 before cut, 1 after cut

  } // end of loop over barrel and endcap

  outFile->Close();     
}


Int_t  vector_angular::rec_phoDecCode(Int_t ipho)
{
  return phoDecCode(patPhotonScEta->at(ipho));
}

Int_t  vector_angular::gen_phoDecCode(Int_t ipho)
{
  return phoDecCode(genParEta_->at(ipho));
}

Int_t  vector_angular::phoDecCode(double eta)
{
  if(fabs(eta) < BARREL_MAXETA)return 0;
  if(fabs(eta) > ENDCAP_MINETA && 
     fabs(eta) < ENDCAP_MAXETA)return 1;
  return -1;
}

Bool_t vector_angular::rec_isFidPho (Int_t ipho)
{
  return isFidPho(patPhotonPt->at(ipho), patPhotonEta->at(ipho));
}

Bool_t vector_angular::gen_isFidPho (Int_t ipart)
{
  return isFidPho(genParPt_->at(ipart), genParEta_->at(ipart));
}

Bool_t vector_angular::isFidPho (double pt, double eta)
{
  if(phoDecCode(eta)<0)return false;
  if(pt < ptbound[0])return false;
  if(pt > ptbound[nPtBins])return false;
  return true;
}



Bool_t vector_angular::isGoodPho(Int_t ipho, bool applyPileUpCorr)
{
  double et  = patPhotonEt   ->at(ipho);

  if(!rec_isFidPho(ipho))return false;
  if(patPhotonhadronicOverEm->at(ipho) > 0.05)return false;
  if(patPhotonhasPixelSeed->at(ipho)   > 1e-6)return false; // this should be saved as bool

  double eciso = phoEcalIso(ipho, applyPileUpCorr); 
  double hciso = phoHcalIso(ipho, applyPileUpCorr);
  double tkiso = phoTrkIso(ipho, applyPileUpCorr);

  if(eciso > 4.2 +0.006*et)return false;
  if(hciso > 2.2 +0.0025*et)return false;
  if(tkiso > 2.0 +0.001*et)return false;

  return true;
}

Double_t  vector_angular::phoEcalIso(Int_t ipho, bool applyPileUpCorr)
{
  if(!rec_isFidPho(ipho))return -9999.0;

  double eciso = patPhotonecalRecHitSumEtConeDR04->at(ipho);
  int decBinIndex = phoDecCode(ipho);

  if(applyPileUpCorr)
    eciso -= aeff[decBinIndex][0]*patPhotonrho25; 
  return eciso;
}
 
Double_t vector_angular::phoHcalIso(Int_t ipho, bool applyPileUpCorr)
{
  if(!rec_isFidPho(ipho))return -9999.0;

  double hciso = patPhotonhcalTowerSumEtConeDR04->at(ipho);
  int decBinIndex = phoDecCode(ipho);

  if(applyPileUpCorr)
    hciso -= aeff[decBinIndex][1]*patPhotonrho25; 
  return hciso;

}

Double_t vector_angular::phoTrkIso(Int_t ipho, bool applyPileUpCorr)
{
  if(!rec_isFidPho(ipho))return -9999.0;

  double tkiso = patPhotontrkSumPtHollowConeDR04->at(ipho);
  int decBinIndex = phoDecCode(ipho);

  if(applyPileUpCorr)
    tkiso -= aeff[decBinIndex][2]*patPhotonrho25; 
  return tkiso;

}


Bool_t vector_angular::rec_isFidJet (Int_t ijet)
{
  return isFidJet(patJetPfAk05Pt_->at(ijet), patJetPfAk05Eta_->at(ijet));
}


Bool_t vector_angular::gen_isFidJet (Int_t ijet)
{
  return isFidJet(genJetPt_->at(ijet), genJetEta_->at(ijet));
}


Bool_t vector_angular::isFidJet (double pt, double eta)
{
  if(pt < 30.0)return false;
  if(fabs(eta) > 2.4)return false;
  return true;
}


// check if this reco-jet is a good loose jet
Bool_t vector_angular::isGoodLooseJet(Int_t ijet)
{
  if(!rec_isFidJet(ijet))return false;
  //  if(patJetPfAk05JetNConstituents_->at(ijet) <= 1)return false;
  if(patJetPfAk05NeutHadEFr_->at(ijet) >= 0.99)return false;
  if(patJetPfAk05NeutEmEFr_->at(ijet) >= 0.99)return false;

  //   // for the tracker region
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharHadEFr_->at(ijet) <= 0.0)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharEmEFr_->at(ijet) >= 0.99)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharMulti_->at(ijet) <= 0)return false;

  return true;

}


// check if this reco-jet is a good medium jet
Bool_t vector_angular::isGoodMediumJet(Int_t ijet)
{
  if(!rec_isFidJet(ijet))return false;
  if(patJetPfAk05Pt_->at(ijet) < 10.0)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet)) > 3.0)return false;
  //  if(patJetPfAk05JetNConstituents_->at(ijet) <= 1)return false;
  if(patJetPfAk05NeutHadEFr_->at(ijet) >= 0.95)return false;
  if(patJetPfAk05NeutEmEFr_->at(ijet) >= 0.95)return false;

  //   // for the tracker region
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharHadEFr_->at(ijet) <= 0.0)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharEmEFr_->at(ijet) >= 0.99)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharMulti_->at(ijet) <= 0)return false;


  return true;

}

// check if this reco-jet is a good medium jet
Bool_t vector_angular::isGoodTightJet(Int_t ijet)
{
  if(!rec_isFidJet(ijet))return false;
  if(patJetPfAk05Pt_->at(ijet) < 10.0)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet)) > 3.0)return false;
  //  if(patJetPfAk05JetNConstituents_->at(ijet) <= 1)return false;
  if(patJetPfAk05NeutHadEFr_->at(ijet) >= 0.90)return false;
  if(patJetPfAk05NeutEmEFr_->at(ijet) >= 0.90)return false;


  //   // for the tracker region
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharHadEFr_->at(ijet) <= 0.0)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharEmEFr_->at(ijet) >= 0.99)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharMulti_->at(ijet) <= 0)return false;

  return true;

}


Int_t vector_angular::matchRecoToGenJet(Int_t ijet)
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


Int_t vector_angular::matchGenToRecoJet(Int_t ijet)
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

