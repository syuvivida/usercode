#define onlymjj_cxx
#include "onlymjj.h"
#include <TH2.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include "cutvalues.h"
#include <standalone_LumiReWeighting.cc>
#include <signalShapeReWeighting.cc>

double deltaR(double eta1, double phi1, double eta2, double phi2)
{
    
  double deta = eta1-eta2;
  double dphi = phi1-phi2;
  while (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
  while (dphi <= -TMath::Pi()) dphi += 2*TMath::Pi();
  double dR = sqrt(deta*deta + dphi*dphi);
  return dR;

}

void onlymjj::Loop(int DEBUG, bool reweighSignal, bool reweighPU)
{
  if (fChain == 0) return;

  //================================================================
  //                    Template histograms
  //================================================================

  TH1D* h_pt_template = new TH1D("h_pt_template","",100,0,500);
  TH1D* h_eta_template = new TH1D("h_eta_template","",100,-5.0,5.0);
  
  TH1D* h_prob_template = new TH1D("h_prob_template","",2,-0.5,1.5);

  TH1D* h_ratio_template = new TH1D("h_ratio_template","",50,0.5,1.5);
  h_ratio_template->Sumw2();
  h_ratio_template->SetYTitle(Form("Events per %.2lf",
				   h_ratio_template->GetBinWidth(1)));

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

  TH1D* h_mjj_template = new TH1D("h_mjj_template","",100,0.0,200.0);
  h_mjj_template->Sumw2();
  h_mjj_template->SetXTitle("M_{jj} [GeV/c^{2}]");
  h_mjj_template->SetYTitle(Form("Events per %d GeV/c^{2}",
				 (int)h_mjj_template->GetBinWidth(1)));

  TProfile* pf_dR_template = new TProfile("pf_dR_template","",50,0.0,5.0,-1000, 1000);
  pf_dR_template->Sumw2();
  pf_dR_template->SetXTitle("#Delta R(q_{1},q_{2})");


  TH1D* h_dR_template = new TH1D("h_dR_template","",60,0,3.0);
  h_dR_template->Sumw2();
  h_dR_template->SetXTitle("#DeltaR");
  h_dR_template->SetYTitle(Form("Events per %.2f",
				h_dR_template->GetBinWidth(1)));

  // separated in dR
  const double dRArray[]={0.5,1.0,1.5,2.0,2.5,3.0,3.5};
  const int NBINS = sizeof(dRArray)/sizeof(dRArray[0])-1;
  TH1D* h_dR = new TH1D("h_dR","",NBINS,dRArray);


  //================================================================
  //                    status=3 histograms
  //================================================================
  
  TH1D* h_rm_gen[NBINS];
  TH1D* h_rm_rec[NBINS];
  TH1D* h_rpt_gen[NBINS][2];
  TH1D* h_rpt_rec[NBINS][2];

  TH1D* h_qpt_parton[NBINS][2];
  TH1D* h_qeta_parton[NBINS][2];

  TH1D* h_matchprob[NBINS][2];

  for(int ibin=0; ibin < NBINS; ibin++){
    
    std::string dRstring = Form("when %.1f < #DeltaR_{qq} < %.1f",
				dRArray[ibin],dRArray[ibin+1]);
    h_rm_gen[ibin] = (TH1D*)h_ratio_template->Clone(Form("h_rm_gen%02i",ibin));
    h_rm_gen[ibin]->SetXTitle(Form("M_{jj}^{GEN}/M_{qq} %s",dRstring.data()));

    h_rm_rec[ibin] = (TH1D*)h_ratio_template->Clone(Form("h_rm_rec%02i",ibin));
    h_rm_rec[ibin]->SetXTitle(Form("M_{jj}^{REC}/M_{qq} %s",dRstring.data()));

    for(int ip=0; ip<2; ip++){
      
      h_rpt_gen[ibin][ip] = (TH1D*)h_ratio_template->Clone(Form("h_rpt_gen%02i_%d",ibin,ip));
      h_rpt_gen[ibin][ip]->SetXTitle(Form("p_{T}(genJet_{%d})/p_{T}(q_{%d}) %s",ip+1,ip+1,dRstring.data()));

      h_rpt_rec[ibin][ip] = (TH1D*)h_ratio_template->Clone(Form("h_rpt_rec%02i_%d",ibin,ip));
      h_rpt_rec[ibin][ip]->SetXTitle(Form("p_{T}(PFJet_{%d})/p_{T}(q_{%d}) %s",ip+1,ip+1,dRstring.data()));

      h_qpt_parton[ibin][ip] = 
	(TH1D*)h_pt_template->Clone(Form("h_qpt_parton%02i_%d",ibin,ip));
      h_qpt_parton[ibin][ip]->SetXTitle(Form("p_{T}(q_{%d}) [GeV] %s",
					     ip+1, dRstring.data()));

      h_qeta_parton[ibin][ip] = 
	(TH1D*)h_eta_template->Clone(Form("h_qeta_parton%02i_%d",ibin,ip));
      h_qeta_parton[ibin][ip]->SetXTitle(Form("#eta(q_{%d}) %s",
					     ip+1, dRstring.data()));

      h_matchprob[ibin][ip] = 
	(TH1D*)h_prob_template->Clone(Form("h_matchprob%02i_%d",ibin,ip));
      h_matchprob[ibin][ip]->SetXTitle(Form("Matched probability of PFJet%d %s",
					    ip+1, dRstring.data()));

    }

  }




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
  
  TH1D* h_mh_parton_event;
  TH1D* h_mll_parton_event;
  TH1D* h_mjj_parton_event; 

  TH1D* h_mh_stable_event;
  TH1D* h_mll_stable_event;
  TH1D* h_mjj_stable_event;
 
  TH1D* h_mh_rec_event;
  TH1D* h_mll_rec_event;
  TH1D* h_mjj_rec_event;

  
  h_mh_parton_event = (TH1D*)h_mh_template->Clone("h_mh_parton_event");
  h_mh_parton_event->SetXTitle("M_{llqq} [GeV/c^{2}]");
  h_mh_stable_event = (TH1D*)h_mh_template->Clone("h_mh_stable_event");
  h_mh_rec_event = (TH1D*)h_mh_template->Clone("h_mh_rec_event");

  h_mll_parton_event = (TH1D*)h_mll_template->Clone("h_mll_parton_event");
  h_mll_stable_event = (TH1D*)h_mll_template->Clone("h_mll_stable_event");
  h_mll_rec_event = (TH1D*)h_mll_template->Clone("h_mll_rec_event");

  h_mjj_parton_event = (TH1D*)h_mjj_template->Clone("h_mjj_parton_event");
  h_mjj_parton_event->SetXTitle("M_{qq} [GeV/c^{2}]");
  h_mjj_stable_event = (TH1D*)h_mjj_template->Clone("h_mjj_stable_event");
  h_mjj_rec_event = (TH1D*)h_mjj_template->Clone("h_mjj_rec_event");


  //================================================================
  //                    status=1, 3 profiles 
  //================================================================
 
  TProfile* pf_dR_Rm_gen = (TProfile*)pf_dR_template->Clone("pf_dR_Rm_gen");
  pf_dR_Rm_gen->SetYTitle("M_{jj}/M_{qq} [GeV]");

  TProfile* pf_dR_Rm_rec = (TProfile*)pf_dR_template->Clone("pf_dR_Rm_rec");
  pf_dR_Rm_rec->SetYTitle("M_{jj}/M_{qq} [GeV]");

  TProfile* pf_dR_dm_gen = (TProfile*)pf_dR_template->Clone("pf_dR_dm_gen");
  pf_dR_dm_gen->SetYTitle("M_{jj}-M_{qq} [GeV]");

  TProfile* pf_dR_dm_rec = (TProfile*)pf_dR_template->Clone("pf_dR_dm_rec");
  pf_dR_dm_rec->SetYTitle("M_{jj}-M_{qq} [GeV]");


  TProfile* pf_dR_Rpt_gen[2];
  TProfile* pf_dR_Rpt_rec[2];
  TProfile* pf_dR_jetArea_rec[2];

  TProfile* pf_dR_matchProb_gen[2];
  TProfile* pf_dR_matchProb_rec[2];

  for(int ip=0; ip<2; ip++){
   
    // ratio
    pf_dR_Rpt_gen[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_Rpt_gen%d",ip));
    pf_dR_Rpt_gen[ip]->SetYTitle(Form("p_{T}(genJet_{%d})/p_{T}(q_{%d}) [GeV]",ip+1,ip+1));
  
    pf_dR_Rpt_rec[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_Rpt_rec%d",ip));
    pf_dR_Rpt_rec[ip]->SetYTitle(Form("p_{T}(PFJet_{%d})/p_{T}(q_{%d}) [GeV]",ip+1,ip+1));

    pf_dR_jetArea_rec[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_jetArea_rec%d",ip));
    pf_dR_jetArea_rec[ip] -> SetYTitle(Form("jet area for PFJet %d",ip+1));

    pf_dR_matchProb_gen[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_matchProb_gen%d",ip));
    pf_dR_matchProb_gen[ip] -> SetYTitle(Form("Matched probability for genJet %d",ip+1));

    pf_dR_matchProb_rec[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_matchProb_rec%d",ip));
    pf_dR_matchProb_rec[ip] -> SetYTitle(Form("Matched probability for PFJet %d",ip+1));

  }
 
  int nPass[50]={0};
  
  Long64_t nentries = fChain->GetEntriesFast();
  standalone_LumiReWeighting LumiWeights_central(2012,0);
  std::string reweight_file = Form("/home/syu/cvs_scripts/scripts/hzz2l2q/data/mZZ_Higgs%d_8TeV_Lineshape+Interference.txt",_higgsMass);
  cout << "reweighing signal with the input file: " << reweight_file << endl;
  signalShapeReWeighting signalWeightMethod(reweight_file.data());

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    //    if(jentry > 1000)break;
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

	// higgs
	if(status==3 && PID==25)
	  {
	    genHiggsMass = genParM_->at(igen);
	    h_mh_parton_mother->Fill(genHiggsMass);
	    if(reweighSignal)
	      signal_weight = signalWeightMethod.weight(genHiggsMass,0);
	    h_mh_parton_mother_weighted->Fill(genHiggsMass,signal_weight);
	  }

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
	    lep1PID   = PID;
	  }

	// second status=3 lepton from Z
	else if (status==3 && lep2.E() < 1e-6 && (abs(PID) == 11 ||
						  abs(PID) == 13)
		 && momPID == ZPID		 
		 && PID == -(lep1PID)
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
		&& abs(PID) == LEPTYPE
		)
	  {
	    lep1_post.SetPtEtaPhiM(
				   genParPt_->at(igen),
				   genParEta_->at(igen),
				   genParPhi_->at(igen),
				   genParM_->at(igen)
				   );
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
	    q1PID   = PID;
	  }

	// second parton from Z
	else if(status==3 && q2.E() < 1e-6 && abs(PID) < 6 
		&& momPID == ZPID 
		&& PID == -(q1PID)
		)
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
      cout << "===========================================" << endl;
      cout << "Event = " << EvtInfo_EventNum << endl;
      cout << "LEPTYPE = " << LEPTYPE << endl;
      cout << "lep1 = " << lep1Index << "\t lep2 = " << lep2Index << endl;
      cout << "lep1Post =" << lep1PostIndex << "\t lep2Post = " << lep2PostIndex << endl;
      cout << "q1 =" << q1Index << "\t q2Index = " << q2Index << endl;
      cout << "===========================================" << endl;
    }
    if(LEPTYPE==-1)continue;

    if(lep1Index<0 || lep2Index<0)continue;
    if(lep1PostIndex<0 || lep2PostIndex<0)continue;
    nPass[1]++;

    if(q1Index<0 || q2Index<0)continue;
    nPass[2]++;


    double PU_weight =  LumiWeights_central.weight(PU_nTrueInt);
    double rec_weight = reweighPU? PU_weight*signal_weight: signal_weight;
    double mH_parton = (lep1+lep2+q1+q2).M();
    double mZll_parton = (lep1+lep2).M();
    double mZjj_parton = (q1+q2).M();
    double dR_parton = q1.DeltaR(q2);

    int dRJJ_index = h_dR->FindBin(dR_parton)-1;
    if(dRJJ_index < 0)dRJJ_index = 0;
    if(dRJJ_index >= NBINS)dRJJ_index = NBINS-1; 
    if(DEBUG==1)
      cout << "dR_parton = " << dR_parton << "\t dRJJ_index = " << 
	dRJJ_index << endl;

    h_mh_parton_daughter->Fill(mH_parton, signal_weight);
    h_mll_parton_daughter->Fill(mZll_parton, signal_weight);
    h_mjj_parton_daughter->Fill(mZjj_parton, signal_weight);
    h_dR_qq->Fill(dR_parton,signal_weight);


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

    double QuarkPt[2]={q1.Pt(), q2.Pt()};
    
    h_qpt_parton[dRJJ_index][0]->Fill(q1.Pt(),signal_weight);
    h_qpt_parton[dRJJ_index][1]->Fill(q2.Pt(),signal_weight);

    h_qeta_parton[dRJJ_index][0]->Fill(q1.Eta(),signal_weight);
    h_qeta_parton[dRJJ_index][1]->Fill(q2.Eta(),signal_weight);
  
    //================================================================
    // STATUS =1       LEVEL
    //================================================================


    double lep1posteta = lep1_post.Eta();
    double lep1postphi = lep1_post.Phi();

    double lep2posteta = lep2_post.Eta();
    double lep2postphi = lep2_post.Phi();

    double mZll_particle = (lep1_post + lep2_post).M();
    
    //================================================================
    //    Starting from quarks
    //================================================================

    // Only pick up the jets that matched and construct them in the dijetmass

    TLorentzVector jet1(0,0,0,0);
    TLorentzVector jet2(0,0,0,0);
    int jet1Index=-1, jet2Index=-1;

    int nGoodJets=0;

    for(unsigned int ijet =0; ijet < genJetPt_->size(); ijet++)
      {

	double eta = genJetEta_->at(ijet);
	double phi = genJetPhi_->at(ijet);
	double pt  = genJetPt_->at(ijet);
	
	if(deltaR(eta,phi, lep1posteta, lep1postphi)<0.5)continue;
	if(deltaR(eta,phi, lep2posteta, lep2postphi)<0.5)continue;

	if(fabs(eta) > 2.4)continue;
 	if(pt < 30.0)continue;

	double matchedProb = 
	  (matchGenJetToParton(q1Index,ijet) || 
	   matchGenJetToParton(q2Index,ijet))?
	  1.0: 0.0;
		
	if(nGoodJets<2)
	  pf_dR_matchProb_gen[nGoodJets]->
	    Fill(dR_parton, matchedProb, signal_weight);

	nGoodJets++;


	if(matchGenJetToParton(q1Index,ijet) && jet1.E()<1e-6)
	  {
	    jet1.SetPtEtaPhiE(
			      genJetPt_->at(ijet),
			      genJetEta_->at(ijet),
			      genJetPhi_->at(ijet),
			      genJetE_->at(ijet));
	    jet1Index = ijet;
	  }

	if(matchGenJetToParton(q2Index,ijet) && jet2.E()<1e-6)
	  {

	    jet2.SetPtEtaPhiE(
			      genJetPt_->at(ijet),
			      genJetEta_->at(ijet),
			      genJetPhi_->at(ijet),
			      genJetE_->at(ijet));
	    jet2Index = ijet;
	  }

      }
    

    if(DEBUG==1)
      {
	cout << "jet1Index = " << jet1Index << "\t jet2Index = " << jet2Index
	     << endl;
	cout << "jet1 pt = " << jet1.Pt() << "\t jet2 pt = " << jet2.Pt() << endl;
	cout << "q1 pt = " << q1.Pt() << "\t q2 pt = " << q2.Pt() << endl;
      }
      
    
    double mZjj_particle = -999;
    if(jet1Index>=0 && jet2Index>=0 && jet1Index!=jet2Index){
      nPass[3]++;
      h_mh_parton_event->Fill(mH_parton, signal_weight);
      h_mll_parton_event->Fill(mZll_parton, signal_weight);
      h_mjj_parton_event->Fill(mZjj_parton, signal_weight);

      double mH_particle = (jet1+jet2+lep1_post+lep2_post).M();
      mZjj_particle = (jet1+jet2).M();

      h_mh_stable_event->Fill(mH_particle, signal_weight);
      h_mll_stable_event->Fill(mZll_particle, signal_weight);
      h_mjj_stable_event->Fill(mZjj_particle, signal_weight);

      double rm = mZjj_particle/mZjj_parton;
      double dm = mZjj_particle-mZjj_parton; 
      double rpt[2]={jet1.Pt()/q1.Pt(),
		     jet2.Pt()/q2.Pt()};

      if(mZll_particle > MIN_MZ_LL       && mZll_particle < MAX_MZ_LL 
//  	 && mZjj_particle > LOOSE_MIN_MZ_JJ && mZjj_particle < LOOSE_MAX_MZ_JJ
	 )
       	{
	  h_rm_gen[dRJJ_index]->Fill(rm, signal_weight);
	  for(int ip=0;ip<2;ip++)
	    h_rpt_gen[dRJJ_index][ip]->Fill(rpt[ip], signal_weight);

	  pf_dR_Rm_gen->Fill(dR_parton, rm, signal_weight);
	  pf_dR_dm_gen->Fill(dR_parton, dm, signal_weight);
		
	  for(int ip=0;ip<2;ip++)
	    pf_dR_Rpt_gen[ip]->Fill(dR_parton, rpt[ip], signal_weight);


	}

    } 
      
    // ------------------------------------------------
    // Reconstruction level
    //-------------------------------------------------    

    for(unsigned int ih=0; ih<higgsM->size(); ih++){

      bool matchedToQ1=false;
      bool matchedToQ2=false;
      
      double jetRecPt[2]={0};
      
      int pfJetIndex[2]={-1,-1};

      for(unsigned int ijet=0; ijet<jetPt->size(); ijet++){

	if(jetHiggsIndex_->at(ijet)!=ih)continue;
	if(jetPt->at(ijet)<30.0)continue;
	if(fabs(jetEta->at(ijet))>2.4)continue;

	int jet_index = jetIndex->at(ijet);
	if(jet_index<0)continue;
	if(jet_index>1)continue;
	
	matchedToQ1=matchPFJetToParton(q1Index, ijet);
	matchedToQ2=matchPFJetToParton(q2Index, ijet);


	double matchedProb =(matchedToQ1 || matchedToQ2)?
	  1.0: 0.0;

	pf_dR_matchProb_rec[jet_index]->
	  Fill(dR_parton,matchedProb, rec_weight);

	h_matchprob[dRJJ_index][jet_index]->Fill(matchedProb, rec_weight);

	if(matchedToQ1 && !matchedToQ2)
	  {
	    jetRecPt[0] = jetPt->at(ijet);
	    pfJetIndex[0] = ijet;
	  }
	else if(matchedToQ2 && !matchedToQ1)
	  {
	    jetRecPt[1] = jetPt->at(ijet);
	    pfJetIndex[1] = ijet;
	  }
      } // end of loop over jets
      
      double mH_rec  = higgsM->at(ih);
      double mll_rec = zllM->at(ih);
      double mjj_rec = zjjM->at(ih);

      if(pfJetIndex[0]>=0 && pfJetIndex[1]>=0){

	nPass[4]++;
	if(DEBUG==1)cout << 
	  "pfJetIndex1 = " << pfJetIndex[0] << "\t" <<
	  "pfJetIndex2 = " << pfJetIndex[1] << endl;	  

	h_mh_rec_event->Fill(mH_rec, rec_weight);
	h_mll_rec_event->Fill(mll_rec, rec_weight);
	h_mjj_rec_event->Fill(mjj_rec, rec_weight);
	
	if(
	   mll_rec > MIN_MZ_LL  && mll_rec < MAX_MZ_LL  
//  	   && mjj_rec > LOOSE_MIN_MZ_JJ && mjj_rec < LOOSE_MAX_MZ_JJ  
	   && mZll_particle > MIN_MZ_LL && mZll_particle < MAX_MZ_LL  
// 	   &&  mZjj_particle > LOOSE_MIN_MZ_JJ && mZjj_particle < LOOSE_MAX_MZ_JJ
 	   && jet1Index>=0 && jet2Index>=0 && jet1Index!=jet2Index
	   )
	  {

	    double rm = mjj_rec/mZjj_parton;
	    h_rm_rec[dRJJ_index]->Fill(rm, rec_weight);
	    double dm = mjj_rec-mZjj_parton;
	
	    double rpt[2];
	    for(int ip=0; ip<2; ip++)
	      {
		rpt[ip]= jetRecPt[ip]/QuarkPt[ip];
		h_rpt_rec[dRJJ_index][ip]->Fill(rpt[ip], rec_weight);
	      }

	    pf_dR_Rm_rec->Fill(dR_parton, rm, rec_weight);
	    pf_dR_dm_rec->Fill(dR_parton, dm, rec_weight);
	    for(int ieiko=0; ieiko<2; ieiko++)
	      {
		if(DEBUG==1)cout << "jetRecPt[" << ieiko << "] = " << 
		  jetRecPt[ieiko] << "\t QuarkPt[" << ieiko << "] = " << 
		  QuarkPt[ieiko] << endl;
		pf_dR_Rpt_rec[ieiko]->Fill(dR_parton, rpt[ieiko],rec_weight);
		pf_dR_jetArea_rec[ieiko]->Fill(dR_parton,
					       jetArea->at(pfJetIndex[ieiko]),
					       rec_weight);

	      }
	  } // after pass dilepton mass cuts
	    
      }

    } // end of loop over higgs candidates
  




  } // end of loop over entries

  for(int i=0; i<50;i++)
    if(nPass[i]>0)
      cout << "nPass[" << i << "] = " << nPass[i] << endl;

  TFile* outFile = new TFile(Form("debugmjj_M%d.root",
				  _higgsMass),
			     "recreate");            
  
  h_dR_qq->Write();
  h_mh_parton_mother->Write();
  h_mh_parton_mother_weighted->Write();
  h_mh_parton_daughter->Write();
  h_mll_parton_daughter->Write();
  h_mjj_parton_daughter->Write();

  pf_dR_Rm_gen->Write();
  pf_dR_Rm_rec->Write();

  pf_dR_dm_gen->Write();
  pf_dR_dm_rec->Write();

  for(int ip=0; ip<2; ip++){

    pf_dR_Rpt_gen[ip]->Write();
    pf_dR_Rpt_rec[ip]->Write();    
    pf_dR_jetArea_rec[ip]->Write();

    pf_dR_matchProb_gen[ip]->Write();
    pf_dR_matchProb_rec[ip]->Write();

  }


  for(int ib=0; ib<NBINS; ib++){      

    h_rm_gen[ib]->Write();
    h_rm_rec[ib]->Write();

    for(int ip=0; ip<2; ip++){
      h_rpt_gen[ib][ip]->Write();
      h_rpt_rec[ib][ip]->Write();
      h_qpt_parton[ib][ip]->Write();
      h_qeta_parton[ib][ip]->Write();
      h_matchprob[ib][ip]->Write();
    }
  }


    
  h_mh_parton_event->Write();
  h_mll_parton_event->Write();
  h_mjj_parton_event->Write();

  h_mh_stable_event->Write();
  h_mll_stable_event->Write();
  h_mjj_stable_event->Write();

  h_mh_rec_event->Write();
  h_mll_rec_event->Write();
  h_mjj_rec_event->Write();


  outFile->Close();  

}


Bool_t onlymjj::matchGenJetToParton(Int_t igen, Int_t ijet){

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

Bool_t onlymjj::matchPFJetToParton(Int_t igen, Int_t ijet){

  Bool_t matched = false;
  double dR = deltaR(genParEta_->at(igen), genParPhi_->at(igen),
		     jetEta->at(ijet), jetPhi->at(ijet)
		     ); 

  Double_t relPt = genParPt_->at(igen)>1e-6?
    fabs(genParPt_->at(igen)-jetPt->at(ijet))/genParPt_->at(igen)
    : -9999.0;

  if(dR<0.5)
    matched = true;
  
  return matched;  


}
