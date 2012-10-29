#define real_VBF_cxx
#include "real_VBF.h"
#include <TH2.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <map>


double deltaPhi(double phi1, double phi2)
{
  double dphi = phi1-phi2;
  while (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
  while (dphi <= -TMath::Pi()) dphi += 2*TMath::Pi();
  return dphi;
}

double deltaEta(double eta1, double eta2)
{
  double deta = eta1-eta2;
  return deta;
}


double deltaR(double eta1, double phi1, double eta2, double phi2)
{    
  double deta = deltaEta(eta1,eta2);
  double dphi = deltaPhi(phi1,phi2);
  double dR = sqrt(deta*deta+dphi*dphi);
  return dR;
}

TLorentzVector real_VBF::setParticleL4(Int_t igen)
{
  TLorentzVector lpart(0,0,0,0);
  lpart.SetPtEtaPhiM(
		     genParPt_->at(igen),
		     genParEta_->at(igen),
		     genParPhi_->at(igen),
		     genParM_->at(igen)
		     );
  return lpart;
}

TLorentzVector real_VBF::setGenJetL4(Int_t ijet)
{
  TLorentzVector ljet(0,0,0,0);
  
  ljet.SetPtEtaPhiE(
		    genJetPt_->at(ijet),
		    genJetEta_->at(ijet),
		    genJetPhi_->at(ijet),
		    genJetE_->at(ijet));
  return ljet;

}


void real_VBF::Loop(int DEBUG)
{
  if (fChain == 0) return;

  //================================================================
  //                    Template histograms
  //================================================================


  TH1D* h_mh_template = new TH1D("h_mh_template","",300,0.0,1500);
  h_mh_template->Sumw2();
  h_mh_template->SetXTitle("M_{lljj} [GeV/c^{2}]");
  h_mh_template->SetYTitle(Form("Events per %d GeV/c^{2}",
				(int)h_mh_template->GetBinWidth(1)));


  TH1D* h_mll_template = new TH1D("h_mll_template","",100,0.0,200.0);
  h_mll_template->Sumw2();
  h_mll_template->SetXTitle("M_{ll} [GeV/c^{2}]");
  h_mll_template->SetYTitle(Form("Events per %d GeV/c^{2}",
				 (int)h_mll_template->GetBinWidth(1)));

  TH1D* h_mjj_template = new TH1D("h_mjj_template","",200,0.0,1000.0);
  h_mjj_template->Sumw2();
  h_mjj_template->SetXTitle("M_{jj} [GeV/c^{2}]");
  h_mjj_template->SetYTitle(Form("Events per %d GeV/c^{2}",
				 (int)h_mjj_template->GetBinWidth(1)));

  TProfile* pf_dR_template = new TProfile("pf_dR_template","",50,0.0,5.0,-1000, 1000);
  pf_dR_template->Sumw2();
  pf_dR_template->SetXTitle("#Delta R(q_{1},q_{2})");


  TH1D* h_dR_template = new TH1D("h_dR_template","",200,0,10.0);
  h_dR_template->Sumw2();
  h_dR_template->SetXTitle("#DeltaR");
  h_dR_template->SetYTitle(Form("Events per %.2f",
				h_dR_template->GetBinWidth(1)));

  TH1D* h_dphi_template = new TH1D("h_dphi_template","",100,0,TMath::Pi());
  h_dphi_template->Sumw2();
  h_dphi_template->SetXTitle("#Delta#phi");
  h_dphi_template->SetYTitle(Form("Events per %.2f",
				  h_dphi_template->GetBinWidth(1)));

  TH1D* h_deta_template = new TH1D("h_deta_template","",200,0,10.0);
  h_deta_template->Sumw2();
  h_deta_template->SetXTitle("#Delta#eta");
  h_deta_template->SetYTitle(Form("Events per %.2f",
				  h_deta_template->GetBinWidth(1)));

  TH1D* h_pt_template = new TH1D("h_pt_template","",100,0.0,200.0);
  h_pt_template->Sumw2();
  h_pt_template->SetXTitle("p_{T} [GeV/c]");
  h_pt_template->SetYTitle(Form("Events per %d GeV/c",
				(int)h_pt_template->GetBinWidth(1)));

  TH1D* h_eta_template = new TH1D("h_eta_template","",100,-5,5);
  h_eta_template->Sumw2();
  h_eta_template->SetXTitle("#eta");
  h_eta_template->SetYTitle(Form("Events per %.2f",
				 h_eta_template->GetBinWidth(1)));

  TH1D* h_phi_template = new TH1D("h_phi_template","",50,-TMath::Pi(),
				  TMath::Pi());
  h_phi_template->Sumw2();
  h_phi_template->SetXTitle("#phi [Radians]");
  h_phi_template->SetYTitle(Form("Events per %.2f radians",
				 h_phi_template->GetBinWidth(1)));


  //================================================================
  //                    status=3 histograms
  //================================================================
  
  TH1D* h_mh_parton_mother = 
    (TH1D*)h_mh_template->Clone("h_mh_parton_mother");
  h_mh_parton_mother->SetTitle("status==3, PID==25");
  h_mh_parton_mother->SetXTitle("M_{H} [GeV/c^{2}]");

  TH1D* h_mh_parton_mother_weighted = 
    (TH1D*)h_mh_template->Clone("h_mh_parton_mother_weighted");
  h_mh_parton_mother_weighted->SetTitle("Reweighed, status==3, PID==25");
  h_mh_parton_mother_weighted->SetXTitle("M_{H} [GeV/c^{2}]");

  TH1D* h_mh_parton_daughter = 
    (TH1D*)h_mh_template->Clone("h_mh_parton_daughter");
  h_mh_parton_daughter->SetTitle("status==3");
  h_mh_parton_daughter->SetXTitle("M_{llqq} [GeV/c^{2}]");

  TH1D* h_mll_parton_daughter = 
    (TH1D*)h_mll_template->Clone("h_mll_parton_daughter");
  h_mll_parton_daughter->SetTitle("status==3");

  TH1D* h_mjj_parton_daughter = 
    (TH1D*)h_mjj_template->Clone("h_mjj_parton_daughter");
  h_mjj_parton_daughter->SetTitle("status==3");
  h_mjj_parton_daughter->SetXTitle("M_{qq} [GeV/c^{2}]");

  TH1D* h_dR_qq = 
    (TH1D*)h_dR_template->Clone("h_dR_qq");
  h_dR_qq->SetXTitle("#Delta R(q_{1},q_{2})");

  //================================================================
  //                    status=1 histograms
  //================================================================


  TH1D* h_vbf_mjj = (TH1D*)h_mjj_template->Clone("h_vbf_mjj");
  TH1D* h_vbf_dRjj= (TH1D*)h_dR_template->Clone("h_vbf_dRjj");

  TH1D* h_vbf_dphijj= (TH1D*)h_dphi_template->Clone("h_vbf_dphijj");
  TH1D* h_vbf_detajj= (TH1D*)h_deta_template->Clone("h_vbf_detajj");

  TH1D* h_vbf_ptjj= (TH1D*)h_pt_template->Clone("h_vbf_ptjj");
  h_vbf_ptjj->SetXTitle("p_{T}(jj) [GeV/c]");


  TH1D* h_jetpt[2];
  TH1D* h_jeteta[2];
  TH1D* h_jetphi[2];
  std::string vbftitle[2]={"Leading jet", "Subleading jet"};

  for(int ij=0; ij<2; ij++)
    {
      h_jetpt[ij] = (TH1D*)h_pt_template->Clone(Form("h_jetpt%d",ij));
      h_jeteta[ij] = (TH1D*)h_eta_template->Clone(Form("h_jeteta%d",ij));
      h_jetphi[ij] = (TH1D*)h_phi_template->Clone(Form("h_jetphi%d",ij));

      h_jetpt[ij]->SetTitle(vbftitle[ij].data());
      h_jeteta[ij]->SetTitle(vbftitle[ij].data());
      h_jetphi[ij]->SetTitle(vbftitle[ij].data());

    }

 
  int nPass[50]={0};
  typedef map<double, int,  std::greater<double> > myMap;
  myMap sorted_partonPtMap;
  typedef myMap::iterator mapIter;
  
  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    sorted_partonPtMap.clear();

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    //     if(jentry > 50)break;
    nPass[0]++;
    //================================================================
    // STATUS=3 LEVEL
    //================================================================

    TLorentzVector lep1(0,0,0,0);
    TLorentzVector lep2(0,0,0,0);
    int lep1Index=-1;
    int lep2Index=-1;

    TLorentzVector lep1_post(0,0,0,0);
    TLorentzVector lep2_post(0,0,0,0);
    int lep1PostIndex=-1;
    int lep2PostIndex=-1;

    TLorentzVector q1(0,0,0,0);
    TLorentzVector q2(0,0,0,0);
    int q1Index=-1;
    int q2Index=-1;
    const int ZPID=23;
    int LEPTYPE = -1;
    int lep1PID = -1;
    int lep1PostPID = -1;
    int q1PID = -999;

    double genHiggsMass = 0.0;
    double signal_weight= 1.0;
    
    
    for(unsigned int igen=0; igen < genParId_->size(); igen++)
      {
	
	int status = genParSt_->at(igen);
	int PID    = genParId_->at(igen);
	int momPID = genMomParId_->at(igen);

	// first status=3 lepton from Z
	if(status==3 && lep1.E() < 1e-6 && (abs(PID) == 11 ||
					    abs(PID) == 13)
	   && momPID == ZPID
	   )
	  {
	    lep1 = setParticleL4(igen);
	    LEPTYPE = abs(PID);
	    lep1Index = igen;
	    lep1PID   = PID;
	  }

	// second status=3 lepton from Z
	else if (status==3 && lep2.E() < 1e-6 && (abs(PID) == 11 ||
						  abs(PID) == 13)
		 && momPID == ZPID		 
		 && PID == -(lep1PID)
		 )
	  {
	    lep2 = setParticleL4(igen);
	    lep2Index = igen;
	  }

	// first status =1 lepton from Z
	else if(status==1 && lep1_post.E() < 1e-6 && (abs(PID) == 11 ||
						      abs(PID) == 13)
		&& abs(momPID) == LEPTYPE
		&& abs(PID) == LEPTYPE
		)
	  {
	    lep1_post = setParticleL4(igen);
	    lep1PostIndex = igen;
	    lep1PostPID = PID;
	  }

	// second status=1 lepton from Z
	else if (status==1 && lep2_post.E() < 1e-6 && (abs(PID) == 11 ||
						       abs(PID) == 13)
		 && abs(momPID) == LEPTYPE
		 && PID == -(lep1PostPID)
		 )
	  {
	    lep2_post = setParticleL4(igen);
	    lep2PostIndex = igen;
	  }

	// first parton from Z
	else if(status==3 && q1.E() < 1e-6 && abs(PID) < 6 
		&& momPID == ZPID)
	  {
	    q1 = setParticleL4(igen);
	    q1Index = igen;
	    q1PID   = PID;
	  }

	// second parton from Z
	else if(status==3 && q2.E() < 1e-6 && abs(PID) < 6 
		&& momPID == ZPID 
		&& PID == -(q1PID)
		)
	  {
	    q2 = setParticleL4(igen);
	    q2Index = igen;
	  }

			    
      } // end of loop over gen particles
    
  
    if(DEBUG==1){
      cout << "===========================================" << endl;
      cout << "Event = " << EvtInfo_EventNum << endl;
      cout << "LEPTYPE = " << LEPTYPE << endl;
      cout << "lep1 = " << lep1Index << "\t lep2 = " << lep2Index << endl;
      cout << "lep1Post =" << lep1PostIndex << "\t lep2Post = " << lep2PostIndex << endl;
      cout << "q1 =" << q1Index << "\t q2Index = " << q2Index << endl;
    }
    if(LEPTYPE==-1)continue;
  
    /////////////////////////////////////////////////////////
    ///
    /// Check only Higgs
    ///
    /////////////////////////////////////////////////////////

    
    int higgs_mom1 = -1;
    int higgs_mom2 = -1;
    for(unsigned int igen=0; igen < genParId_->size(); igen++)
      {	
	int status = genParSt_->at(igen);
	int PID    = genParId_->at(igen);
	// higgs
	if(status==3 && PID==25)
	  {
	    genHiggsMass = genParM_->at(igen);
	    h_mh_parton_mother->Fill(genHiggsMass);
	    h_mh_parton_mother_weighted->Fill(genHiggsMass,signal_weight);

	    higgs_mom1 = genMo1_->at(igen);
	    higgs_mom2 = genMo2_->at(igen);
	    break;
	  }
      }
    //     if(DEBUG==1)
    //       cout << "higgs mother 1 = " << higgs_mom1 << "\t" << higgs_mom2 << endl;
    for(unsigned int igen=0; igen < genParId_->size(); igen++)
      {	
	int status = genParSt_->at(igen);
	int PID    = genParId_->at(igen);
	int mom1   = genMo1_->at(igen);
	int mom2   = genMo2_->at(igen);
	double pt =  genParPt_->at(igen);

	if(status==3 && (abs(PID)<=6 || abs(PID)==21)  
	   && mom1 > 0 
	   && (mom1 == higgs_mom1) 
	   && (mom2 == higgs_mom2)) // 
	  {
	    sorted_partonPtMap.insert(std::pair<double, int>(pt,igen));  
	  }
      }

    ///////////////////////////////////////////////////////////

    if(lep1Index<0 || lep2Index<0)continue;
    if(lep1PostIndex<0 || lep2PostIndex<0)continue;
    if(q1Index<0 || q2Index<0)continue;

    nPass[1]++;

    double mH_parton = (lep1+lep2+q1+q2).M();
    double mZll_parton = (lep1+lep2).M();
    double mZjj_parton = (q1+q2).M();
    double dR_parton = q1.DeltaR(q2);

    h_mh_parton_daughter->Fill(mH_parton, signal_weight);
    h_mll_parton_daughter->Fill(mZll_parton, signal_weight);
    h_mjj_parton_daughter->Fill(mZjj_parton, signal_weight);
    h_dR_qq->Fill(dR_parton,signal_weight);
    nPass[2]++;



    ///////////////////////////////////////////////////////////
    //
    //      Look for real_VBF partons
    //
    //////////////////////////////////////////////////////////

    if(sorted_partonPtMap.size()<2)continue;
    nPass[3]++;

    // do the pt ordering of quarks
    if(DEBUG==1)
      cout << "q1 pt = " << q1.Pt() << "\t q2 pt = " << q2.Pt() << endl;

    TLorentzVector qTemp(0,0,0,0);
    int qTempIndex = -1;
    if(q2.Pt() > q1.Pt()){
      qTemp = q1;       qTempIndex = q1Index;
      q1    = q2;       q1Index    = q2Index;
      q2    = qTemp;    q2Index    = qTempIndex;
    }

    if(DEBUG==1)
      cout << "After sorting q1 pt = " << q1.Pt() << "\t q2 pt = " << q2.Pt() << endl << 
	" q1 index = " << q1Index << "\t q2 index = " << q2Index << endl;


    
  
    //================================================================
    // STATUS =1       LEVEL
    //================================================================


    double lep1posteta = lep1_post.Eta();
    double lep1postphi = lep1_post.Phi();

    double lep2posteta = lep2_post.Eta();
    double lep2postphi = lep2_post.Phi();

    int vbfjet1Index=-1, vbfjet2Index=-1;
    int vbfparton1Index=-1, vbfparton2Index=-1;
    bool findAPair=false;

    //================================================================
    // Plot only those that are matched to VBF partons
    //================================================================

    for(unsigned int ijet =0; ijet < genJetPt_->size(); ijet++)
      {
	double etai = genJetEta_->at(ijet);
	double phii = genJetPhi_->at(ijet);
	double pti  = genJetPt_->at(ijet);

  	if(deltaR(etai,phii, lep1posteta, lep1postphi)<0.5)continue;
  	if(deltaR(etai,phii, lep2posteta, lep2postphi)<0.5)continue;
	if(pti < 10.0)continue;
	if(fabs(etai) > 5.0)continue;
	
   	for(unsigned int kjet =0; kjet < ijet; kjet++)
	  {
	    
	    double etak = genJetEta_->at(kjet);
	    double phik = genJetPhi_->at(kjet);
	    double ptk  = genJetPt_->at(kjet);

  	    if(deltaR(etak,phik, lep1posteta, lep1postphi)<0.5)continue;
  	    if(deltaR(etak,phik, lep2posteta, lep2postphi)<0.5)continue;
	    if(ptk < 10.0)continue;
	    if(fabs(etak) > 5.0)continue;

	    for (mapIter it_part= sorted_partonPtMap.begin();
		 it_part != sorted_partonPtMap.end();  ++it_part)
	      {
		int parton_index = it_part->second;
		if(matchPartonToGenJet(parton_index,kjet))
		  vbfparton1Index=parton_index;
		if(matchPartonToGenJet(parton_index,ijet))
		  vbfparton2Index=parton_index;
	      }


	    if(vbfparton1Index < 0 || vbfparton2Index <0)continue;
	    if(vbfparton1Index == vbfparton2Index)continue;

	    findAPair = true;
	    
	    vbfjet1Index = kjet;
	    vbfjet2Index = ijet;

	    if(findAPair)break;

	  } // end of loop over kjet

	if(findAPair)break;

      }// end of loop over ijet
    

    if(vbfjet1Index>=0 && vbfjet2Index>=0)
      {
	nPass[4]++; 
	// there is no match if one of the VBF jet is too forward |eta| > 6.0 
	// 1%
	// or if the VBF jet is withint deltaR 0.5 of leptons
	// 1%
	// matched to the same parton
	// 2%
	// pt > 10 GeV, |eta| < 5.0
	// 10%

	TLorentzVector vbf_jet1(0,0,0,0);
	TLorentzVector vbf_jet2(0,0,0,0);
	vbf_jet1 = setGenJetL4(vbfjet1Index);
	vbf_jet2 = setGenJetL4(vbfjet2Index);

	double dr_jj = vbf_jet1.DeltaR(vbf_jet2);
	if(//DEBUG==1 && 
	   dr_jj < 0.5)
	  {
	    cout << "Event " << EvtInfo_EventNum << endl;
	    cout << "VBF parton1 index = " << vbfparton1Index << "\t"    
		 << "VBF parton2 index = " << vbfparton2Index << endl;
	    cout << "VBF jet1 index = " << vbfjet1Index << "\t"    
		 << "VBF jet2 index = " << vbfjet2Index << endl;
	  
	  }

	h_vbf_mjj->Fill((vbf_jet1+ vbf_jet2).M(), signal_weight);
	h_vbf_dRjj->Fill(dr_jj, signal_weight);

	h_vbf_dphijj->Fill(fabs(deltaPhi(vbf_jet1.Phi(),vbf_jet2.Phi())));
	h_vbf_detajj->Fill(fabs(deltaEta(vbf_jet1.Eta(),vbf_jet2.Eta())));

	h_vbf_ptjj->Fill((vbf_jet1+ vbf_jet2).Pt(), signal_weight);
            
	h_jetpt[0]->Fill(vbf_jet1.Pt(), signal_weight);
	h_jeteta[0]->Fill(vbf_jet1.Eta(), signal_weight);
	h_jetphi[0]->Fill(vbf_jet1.Phi(), signal_weight);
      
	h_jetpt[1]->Fill(vbf_jet2.Pt(), signal_weight);
	h_jeteta[1]->Fill(vbf_jet2.Eta(), signal_weight);
	h_jetphi[1]->Fill(vbf_jet2.Phi(), signal_weight);
      }
    //     else
    //       cout << "Event: " << EvtInfo_EventNum << endl;

    nPass[5]++;


  } // end of loop over entries

  for(int i=0; i<50;i++)
    if(nPass[i]>0)
      cout << "nPass[" << i << "] = " << nPass[i] << endl;

  TFile* outFile = new TFile(Form("studyreal_VBF_%s",
				  _inputFileName.data()),
			     "recreate");            
  

  h_vbf_mjj ->Write();
  h_vbf_dRjj->Write();
  h_vbf_dphijj->Write();
  h_vbf_detajj->Write();
  h_vbf_ptjj->Write();

  for(int ij=0; ij<2; ij++)
    {
      h_jetpt[ij]  ->Write();
      h_jeteta[ij] ->Write();
      h_jetphi[ij] ->Write();
    }

  h_dR_qq->Write();

  h_mh_parton_mother->Write();
  h_mh_parton_mother_weighted->Write();
  h_mh_parton_daughter->Write();
  h_mll_parton_daughter->Write();
  h_mjj_parton_daughter->Write();


  outFile->Close();  

}


Bool_t real_VBF::matchPartonToGenJet(Int_t igen, Int_t ijet){

  Bool_t matched = false;
  double dR = deltaR(genParEta_->at(igen), genParPhi_->at(igen),
		     genJetEta_->at(ijet), genJetPhi_->at(ijet)
		     ); 

  Double_t relPt = genParPt_->at(igen)>1e-6?
    fabs(genParPt_->at(igen)-genJetPt_->at(ijet))/genParPt_->at(igen)
    : -9999.0;

  if(dR<0.5)
    matched = true;
   
  return matched;  


}

