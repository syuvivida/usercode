#define mjj_cxx
#include "mjj.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <standalone_LumiReWeighting.cc>
#include <TLorentzVector.h>
#include <TMath.h>
Double_t deltaR(double eta1, double phi1, double eta2, double phi2)
{
    
  Double_t deta = eta1-eta2;
  double dphi = phi1-phi2;
  while (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
  while (dphi <= -TMath::Pi()) dphi += 2*TMath::Pi();
  double dR = sqrt(deta*deta+dphi*dphi);
  return dR;

}

void mjj::Loop(int DEBUG)
{
  if (fChain == 0) return;

  // book histograms
  TH1D* h_mh_template = new TH1D("h_mh_template","",280,100,1500);
  h_mh_template->Sumw2();
  h_mh_template->SetXTitle("M_{lljj} [GeV/c^{2}]");
  h_mh_template->SetYTitle(Form("Candidates per %d GeV/c^{2}",
				(int)h_mh_template->GetBinWidth(1)));


  TH1D* h_mll_template = new TH1D("h_mll_template","",35,60,130);
  h_mll_template->Sumw2();
  h_mll_template->SetXTitle("M_{ll} [GeV/c^{2}]");
  h_mll_template->SetYTitle(Form("Candidates per %d GeV/c^{2}",
				 (int)h_mll_template->GetBinWidth(1)));

  TH1D* h_mjj_template = new TH1D("h_mjj_template","",35,60,130);
  h_mjj_template->Sumw2();
  h_mjj_template->SetXTitle("M_{jj} [GeV/c^{2}]");
  h_mjj_template->SetYTitle(Form("Candidates per %d GeV/c^{2}",
				 (int)h_mjj_template->GetBinWidth(1)));


  TH1D* h_jec_template = new TH1D("h_jec_template","",50,0,2.5);
  h_jec_template->Sumw2();
  h_jec_template->SetXTitle("p_{T}^{REC}(jet)/p_{T}^{GEN}(jet)");
  h_jec_template->SetYTitle(Form("Candidates per %.1f",
				 h_jec_template->GetBinWidth(1)));
   

  TH1D* h_dR_template = new TH1D("h_dR_template","",50,0,5.0);
  h_dR_template->Sumw2();
  h_dR_template->SetXTitle("#Delta R");
  h_dR_template->SetYTitle(Form("Candidates per %.1f",
				 h_dR_template->GetBinWidth(1)));

  // status=3 level
  TH1D* h_mh_parton = (TH1D*)h_mh_template->Clone("h_mh_parton");
  h_mh_parton->SetTitle("status=3");

  TH1D* h_mll_parton = (TH1D*)h_mll_template->Clone("h_mll_parton");
  h_mll_parton->SetTitle("status=3");

  TH1D* h_mjj_parton = (TH1D*)h_mjj_template->Clone("h_mjj_parton");
  h_mjj_parton->SetTitle("status=3");

  // status=1 level
  TH1D* h_mh_stable = (TH1D*)h_mh_template->Clone("h_mh_stable");
  h_mh_stable->SetTitle("status=1");

  TH1D* h_mll_stable = (TH1D*)h_mll_template->Clone("h_mll_stable");
  h_mll_stable->SetTitle("status=1");

  TH1D* h_mjj_stable = (TH1D*)h_mjj_template->Clone("h_mjj_stable");
  h_mjj_stable->SetTitle("status=1");

  // reconstruction level
  TH1D* h_mh_rec = (TH1D*)h_mh_template->Clone("h_mh_rec");
  h_mh_rec->SetTitle("reconstructed, all");

  TH1D* h_mll_rec = (TH1D*)h_mll_template->Clone("h_mll_rec");
  h_mll_rec->SetTitle("reconstructed, all");

  TH1D* h_mjj_rec = (TH1D*)h_mjj_template->Clone("h_mjj_rec");
  h_mjj_rec->SetTitle("reconstructed, all");

  TH1D* h_jec[2];
   
  std::string jetName[2]={"leading","subleading"};
  for(int i=0; i<2; i++)
    {
      h_jec[i] = (TH1D*)h_jec_template->Clone(Form("h_jec_%s",jetName[i].data()));
    }


  // separated in dR
  const double dRArray[]={0.5,1.0,1.5,2.0,2.5,3.0,3.5};
  const int NBINS = sizeof(dRArray)/sizeof(dRArray[0])-1;

  TH1D* h_dR = new TH1D("h_dR","",NBINS,dRArray);
  TH1D* h_dR_ll = (TH1D*)h_dR_template->Clone("h_dR_ll");
  h_dR_ll->SetXTitle("#DeltaR_{ll}");

  TH1D* h_dR_jj = (TH1D*)h_dR_template->Clone("h_dR_jj");
  h_dR_jj->SetXTitle("#DeltaR_{jj}");

  TH1D* h_mll_recdR[NBINS];
  TH1D* h_mjj_recdR[NBINS];
  TH1D* h_jec_dR[NBINS][2];

  for(int i=0; i<NBINS; i++){
    h_mll_recdR[i] = (TH1D*)h_mll_template->Clone(
						  Form("h_mll_recdR%d",i));
    h_mll_recdR[i]->SetTitle(Form("reconstructed, %.1f < #DeltaR_{ll} < %.1f",
				  dRArray[i],dRArray[i+1]));
    h_mjj_recdR[i] = (TH1D*)h_mjj_template->Clone(
						  Form("h_mjj_recdR%d",i));
    h_mjj_recdR[i]->SetTitle(Form("reconstructed, %.1f < #DeltaR_{jj} < %.1f",
				  dRArray[i],dRArray[i+1]));
     
    for(int j=0; j<2; j++){
      h_jec_dR[i][j] = (TH1D*)h_jec_template->Clone(
						    Form("h_jec_dR%d_%s",
							 i,jetName[j].data()));
      h_jec_dR[i][j]->SetTitle(Form("reconstructed, %.1f < #DeltaR_{jj} < %.1f",
				    dRArray[i],dRArray[i+1]));

    }

  }
  Long64_t nentries = fChain->GetEntriesFast();
  standalone_LumiReWeighting LumiWeights_central(0);

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    double PU_weight =  LumiWeights_central.weight(PU_nTrueInt);


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
	    lep1.SetPtEtaPhiM(
			      genParPt_->at(igen),
			      genParEta_->at(igen),
			      genParPhi_->at(igen),
			      genParM_->at(igen)
			      );

	    LEPTYPE = abs(PID);
	    lep1Index = igen;
	  }

	// second status=3 lepton from Z
	else if (status==3 && lep2.E() < 1e-6 && (abs(PID) == 11 ||
						  abs(PID) == 13)
		 && momPID== ZPID
		 )
	  {
	    lep2.SetPtEtaPhiM(
			      genParPt_->at(igen),
			      genParEta_->at(igen),
			      genParPhi_->at(igen),
			      genParM_->at(igen)
			      );
	    lep2Index = igen;
	  }

	// first status =1 lepton from Z
	else if(status==1 && lep1_post.E() < 1e-6 && (abs(PID) == 11 ||
						      abs(PID) == 13)	
		&& abs(momPID) == LEPTYPE
		)
	  {
	    lep1_post.SetPtEtaPhiM(
				   genParPt_->at(igen),
				   genParEta_->at(igen),
				   genParPhi_->at(igen),
				   genParM_->at(igen)
				   );
	    lep1PostIndex = igen;
	  }

	// second status=1 lepton from Z
	else if (status==1 && lep2_post.E() < 1e-6 && (abs(PID) == 11 ||
						       abs(PID) == 13)
		 && abs(momPID) == LEPTYPE

		 )
	  {
	    lep2_post.SetPtEtaPhiM(
				   genParPt_->at(igen),
				   genParEta_->at(igen),
				   genParPhi_->at(igen),
				   genParM_->at(igen)
				   );
	    lep2PostIndex = igen;
	  }

	// first parton from Z
	else if(status==3 && q1.E() < 1e-6 && abs(PID) < 6 
		&& momPID == ZPID)
	  {
	    q1.SetPtEtaPhiM(
			    genParPt_->at(igen),
			    genParEta_->at(igen),
			    genParPhi_->at(igen),
			    genParM_->at(igen)
			    );	 
	    q1Index = igen;
	  }

	// second parton from Z
	else if(status==3 && q2.E() < 1e-6 && abs(PID) < 6 
		&& momPID == ZPID)
	  {
	    q2.SetPtEtaPhiM(
			    genParPt_->at(igen),
			    genParEta_->at(igen),
			    genParPhi_->at(igen),
			    genParM_->at(igen)
			    );	    
	    q2Index = igen;
	  }


      } // end of loop over gen particles

    if(DEBUG==1){
      cout << "LEPTYPE = " << LEPTYPE << endl;
      cout << "lep1 = " << lep1Index << "\t lep2 = " << lep2Index << endl;
      cout << "lep1Post =" << lep1PostIndex << "\t lep2Post = " << lep2PostIndex << endl;
      cout << "q1 =" << q1Index << "\t q2Index = " << q2Index << endl;
      //       cout << "jet1 = " << jet1Index << "\t jet2Index = " << jet2Index << endl;
    }
    
    if(LEPTYPE==-1)continue;
    // now look for gen jets

    //================================================================
    // STATUS =1       LEVEL
    //================================================================


    TLorentzVector jet1(0,0,0,0);
    TLorentzVector jet2(0,0,0,0);
    int jet1Index=-1, jet2Index=-1;
    //     int leadingJetIndex = -1;
    //     double maxJetPt = -9999.0;
    //     int secondLeadingJetIndex = -1;
    //     double secondJetMaxPt = -9999.0;

    for(unsigned int ijet =0; ijet < genJetPt_->size(); ijet++)
      {
	

 	if(matchGenToParton(q1Index,ijet) && jet1.E()<1e-6)
 	  {

 	    jet1.SetPtEtaPhiE(
 			      genJetPt_->at(ijet),
 			      genJetEta_->at(ijet),
 			      genJetPhi_->at(ijet),
 			      genJetE_->at(ijet));
 	    jet1Index = ijet;
 	  }
 	else if(matchGenToParton(q2Index,ijet) && jet2.E()<1e-6)
 	  {

 	    jet2.SetPtEtaPhiE(
 			      genJetPt_->at(ijet),
 			      genJetEta_->at(ijet),
 			      genJetPhi_->at(ijet),
 			      genJetE_->at(ijet));
 	    jet2Index = ijet;
 	  }

	// 	double thisJetPt = genJetPt_->at(ijet);

	// 	// find the highest pt jet
	// 	if(thisJetPt > maxJetPt)
	// 	  {
	// 	    secondJetMaxPt = maxJetPt;
	// 	    secondLeadingJetIndex = leadingJetIndex;

	// 	    maxJetPt = thisJetPt;
	// 	    leadingJetIndex= ijet;
	// 	  }
	// 	else if(thisJetPt > secondJetMaxPt)
	// 	  {
	// 	    secondJetMaxPt = thisJetPt;
	// 	    secondLeadingJetIndex = ijet;
	// 	  }
	

      }

    //      jet1Index = leadingJetIndex;
    //      jet2Index = secondLeadingJetIndex;


    double mH_parton = (lep1+lep2+q1+q2).M();

    double mZll_parton = (lep1+lep2).M();
    double mZjj_parton = (q1+q2).M();

    h_mh_parton->Fill(mH_parton);
    h_mll_parton->Fill(mZll_parton);
    h_mjj_parton->Fill(mZjj_parton);


    if(jet1Index>=0 && jet2Index>=0){

      //     jet1.SetPtEtaPhiE(
      //      		      genJetPt_->at(jet1Index),
      //      		      genJetEta_->at(jet1Index),
      //      		      genJetPhi_->at(jet1Index),
      //      		      genJetE_->at(jet1Index)
      //      		      );


      //     jet2.SetPtEtaPhiE(
      //      		      genJetPt_->at(jet2Index),
      //      		      genJetEta_->at(jet2Index),
      //      		      genJetPhi_->at(jet2Index),
      //      		      genJetE_->at(jet2Index)
      //      		      );


      double mH_particle = (jet1+jet2+lep1_post+lep2_post).M();
      double mZll_particle = (lep1_post + lep2_post).M();
      double mZjj_particle = (jet1+jet2).M();

      h_mh_stable->Fill(mH_particle, PU_weight);
      h_mll_stable->Fill(mZll_particle,PU_weight);
      h_mjj_stable->Fill(mZjj_particle,PU_weight);

    }

    //================================================================
    // RECONSTRUCTION LEVEL
    //================================================================
    
    for(int ih=0; ih < higgsM->size(); ih++){
      
      int best = bestHCand;

      double mh_rec = higgsM->at(best);
      double mll_rec = zllM->at(best);
      double mjj_rec = zjjM->at(best);

      h_mh_rec->Fill(mh_rec, PU_weight);
      h_mll_rec->Fill(mll_rec, PU_weight);
      h_mjj_rec->Fill(mjj_rec, PU_weight);


      double dr_ll = zlldR->at(best);
      h_dR_ll->Fill(dr_ll, PU_weight);
      int dRLL_index = h_dR->FindBin(dr_ll)-1;

      if(dRLL_index>=0 && dRLL_index < NBINS)
	h_mll_recdR[dRLL_index]->Fill(mll_rec, PU_weight);

      double dr_jj = zjjdR->at(best);
      h_dR_jj->Fill(dr_jj, PU_weight);
      int dRJJ_index = h_dR->FindBin(dr_jj)-1;

      if(dRJJ_index>=0 && dRJJ_index < NBINS)
	h_mjj_recdR[dRJJ_index]->Fill(mjj_rec, PU_weight);

      
      for(int ijet=0; ijet < jetIndex->size(); ijet++){

	int jet_index = jetIndex->at(ijet);
	
	if(jet_index < 0 ) continue;
	if(jet_index > 1 ) continue;
	if(jetHiggsIndex->at(ijet)!=best)continue;
	if(jetGenPt->at(ijet)<1e-6) continue;

	double ratio = jetPt->at(ijet)/jetGenPt->at(ijet);
	h_jec[jet_index]->Fill(ratio, PU_weight);

	if(dRJJ_index>=0 && dRJJ_index <NBINS)
	  {
	    h_jec_dR[dRJJ_index][jet_index]->Fill(ratio, PU_weight);
	  }

	
      } // end of loop over jets


    }
  }


  TFile* outFile = new TFile(Form("mjj_%s",_inputFileName.data()),
			     "recreate");            

  h_mh_parton->Write();
  h_mll_parton->Write();
  h_mjj_parton->Write();
  h_mh_stable->Write();
  h_mll_stable->Write();
  h_mjj_stable->Write();

  h_mh_rec->Write();
  h_mll_rec->Write();
  h_mjj_rec->Write();

  h_dR_ll->Write();
  h_dR_jj->Write();

  for(int j=0; j<2; j++)
    {
      h_jec[j]->Write();
      
      for(int i=0; i<NBINS; i++){
	h_mll_recdR[i]->Write();
	h_mjj_recdR[i]->Write();
	h_jec_dR[i][j]->Write();
      }

    }

  outFile->Close();  

}


Bool_t mjj::matchGenToParton(Int_t igen, Int_t ijet){

  Bool_t matched = false;
  Double_t dR = deltaR(genParEta_->at(igen), genParPhi_->at(igen),
		       genJetEta_->at(ijet), genJetPhi_->at(ijet)
		       ); 

  Double_t relPt = genParPt_->at(igen)>1e-6? 
    fabs(genParPt_->at(igen)-genJetPt_->at(ijet))/genParPt_->at(igen)
    : -9999.0;

  if(dR<0.4 && relPt < 3.0)
    matched = true;
  
  return matched;  


}
