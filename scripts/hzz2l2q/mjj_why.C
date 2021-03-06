#define mjj_why_cxx
#include "mjj_why.h"
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
  double dR = sqrt(deta*deta+dphi*dphi);
  return dR;

}

void mjj_why::Loop(int DEBUG, bool reweightSignal)
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

  TH1D* h_mjj_template = new TH1D("h_mjj_template","",100,0.0,200.0);
  h_mjj_template->Sumw2();
  h_mjj_template->SetXTitle("M_{jj} [GeV/c^{2}]");
  h_mjj_template->SetYTitle(Form("Events per %d GeV/c^{2}",
				 (int)h_mjj_template->GetBinWidth(1)));

  TProfile* pf_dR_template = new TProfile("pf_dR_template","",50,0.0,5.0,-1000, 1000);
  pf_dR_template->Sumw2();
  pf_dR_template->SetXTitle("#Delta R(q_{1},q_{2})");

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


  //================================================================
  //                    status=1 histograms
  //================================================================
  const int NCASES=6;
  std::string title[NCASES]={
    "matched and separated (#geq 2 jets)",  //0 
    "matched and merged (==1 jet)",         //1
    "matched and merged (#geq 2 jets)",     //2
    "only one quark matched (==1 jet)",     //3
    "only one quark matched (#geq 2 jets)", //4
    "no match (#geq 2 jets)"                //5
  };
 
  std::string rectitle[NCASES]={
    "matched and separated (#geq 2 jets)",  //0 
    "matched and separated (#geq 2 jets)",  //1
    "matched and merged (#geq 2 jets)",     //2
    "matched and merged (#geq 2 jets)",     //3
    "only one quark matched (#geq 2 jets)", //4
    "no match (#geq 2 jets)"                //5
  };

  
  TH1D* h_mh_parton_event[NCASES];
  TH1D* h_mll_parton_event[NCASES];
  TH1D* h_mjj_parton_event[NCASES]; 

  TH1D* h_mh_stable_event[NCASES];
  TH1D* h_mll_stable_event[NCASES];
  TH1D* h_mjj_stable_event[NCASES];
 
  TH1D* h_mh_parton_random[NCASES];
  TH1D* h_mll_parton_random[NCASES];
  TH1D* h_mjj_parton_random[NCASES]; 

  TH1D* h_mh_stable_random[NCASES];
  TH1D* h_mll_stable_random[NCASES];
  TH1D* h_mjj_stable_random[NCASES];

  TH1D* h_mh_rec_random[NCASES];
  TH1D* h_mll_rec_random[NCASES];
  TH1D* h_mjj_rec_random[NCASES];

  TH1D* h_mh_rec_truth[NCASES];
  TH1D* h_mll_rec_truth[NCASES];
  TH1D* h_mjj_rec_truth[NCASES];

  for(int i=0; i< NCASES; i++){

    //================================================================
    //                    event-level histograms
    //================================================================

    h_mh_parton_event[i] = (TH1D*)h_mh_template->Clone(Form("h_mh_parton_event%d",i));
    h_mh_parton_event[i]->SetTitle(title[i].data());
    h_mh_parton_event[i]->SetXTitle("M_{llqq} [GeV/c^{2}]");

    h_mh_stable_event[i] = (TH1D*)h_mh_template->Clone(Form("h_mh_stable_event%d",i));
    h_mh_stable_event[i]->SetTitle(title[i].data());

    h_mll_parton_event[i] = (TH1D*)h_mll_template->Clone(Form("h_mll_parton_event%d",i));
    h_mll_parton_event[i]->SetTitle(title[i].data());

    h_mll_stable_event[i] = (TH1D*)h_mll_template->Clone(Form("h_mll_stable_event%d",i));
    h_mll_stable_event[i]->SetTitle(title[i].data());

    h_mjj_parton_event[i] = (TH1D*)h_mjj_template->Clone(Form("h_mjj_parton_event%d",i));
    h_mjj_parton_event[i]->SetTitle(title[i].data());
    h_mjj_parton_event[i]->SetXTitle("M_{qq} [GeV/c^{2}]");

    h_mjj_stable_event[i] = (TH1D*)h_mjj_template->Clone(Form("h_mjj_stable_event%d",i));
    h_mjj_stable_event[i]->SetTitle(title[i].data());

    //================================================================
    //                    candidate-level histograms
    //================================================================

    std::string ytitle = Form("Candidates per %d GeV/c^{2}",
			      (int)h_mh_template->GetBinWidth(1));

    h_mh_parton_random[i] = (TH1D*)h_mh_template->Clone(Form("h_mh_parton_random%d",i));
    h_mh_parton_random[i]->SetTitle(title[i].data());
    h_mh_parton_random[i]->SetXTitle("M_{llqq} [GeV/c^{2}]");
    h_mh_parton_random[i]->SetYTitle(ytitle.data());

    h_mh_stable_random[i] = (TH1D*)h_mh_template->Clone(Form("h_mh_stable_random%d",i));
    h_mh_stable_random[i]->SetTitle(title[i].data());
    h_mh_stable_random[i]->SetYTitle(ytitle.data());

    h_mh_rec_random[i] = (TH1D*)h_mh_template->Clone(Form("h_mh_rec_random%d",i));
    h_mh_rec_random[i]->SetTitle(rectitle[i].data());
    h_mh_rec_random[i]->SetYTitle(ytitle.data());

    h_mh_rec_truth[i] = (TH1D*)h_mh_template->Clone(Form("h_mh_rec_truth%d",i));
    h_mh_rec_truth[i]->SetTitle(rectitle[i].data());
    h_mh_rec_truth[i]->SetYTitle(ytitle.data());
    h_mh_rec_truth[i]->SetXTitle("M_{llqq} [GeV/c^{2}]");


    ytitle = Form("Candidates per %d GeV/c^{2}",
		  (int)h_mll_template->GetBinWidth(1));
    h_mll_parton_random[i] = (TH1D*)h_mll_template->Clone(Form("h_mll_parton_random%d",i));
    h_mll_parton_random[i]->SetTitle(title[i].data());
    h_mll_parton_random[i]->SetYTitle(ytitle.data());   

    h_mll_stable_random[i] = (TH1D*)h_mll_template->Clone(Form("h_mll_stable_random%d",i));
    h_mll_stable_random[i]->SetTitle(title[i].data());
    h_mll_stable_random[i]->SetYTitle(ytitle.data());   

    h_mll_rec_random[i] = (TH1D*)h_mll_template->Clone(Form("h_mll_rec_random%d",i));
    h_mll_rec_random[i]->SetTitle(rectitle[i].data());
    h_mll_rec_random[i]->SetYTitle(ytitle.data());   

    h_mll_rec_truth[i] = (TH1D*)h_mll_template->Clone(Form("h_mll_rec_truth%d",i));
    h_mll_rec_truth[i]->SetTitle(rectitle[i].data());
    h_mll_rec_truth[i]->SetYTitle(ytitle.data());

    ytitle = Form("Candidates per %d GeV/c^{2}",
		  (int)h_mjj_template->GetBinWidth(1));
    h_mjj_parton_random[i] = (TH1D*)h_mjj_template->Clone(Form("h_mjj_parton_random%d",i));
    h_mjj_parton_random[i]->SetTitle(title[i].data());
    h_mjj_parton_random[i]->SetXTitle("M_{qq} [GeV/c^{2}]");
    h_mjj_parton_random[i]->SetYTitle(ytitle.data()); 
 
    h_mjj_stable_random[i] = (TH1D*)h_mjj_template->Clone(Form("h_mjj_stable_random%d",i));
    h_mjj_stable_random[i]->SetTitle(title[i].data());
    h_mjj_stable_random[i]->SetYTitle(ytitle.data()); 

    h_mjj_rec_random[i] = (TH1D*)h_mjj_template->Clone(Form("h_mjj_rec_random%d",i));
    h_mjj_rec_random[i]->SetTitle(rectitle[i].data());
    h_mjj_rec_random[i]->SetYTitle(ytitle.data()); 

    if(i==3)h_mjj_rec_random[i]->SetXTitle("M_{j} [GeV/c^{2}]");

    h_mjj_rec_truth[i] = (TH1D*)h_mjj_template->Clone(Form("h_mjj_rec_truth%d",i));
    h_mjj_rec_truth[i]->SetTitle(rectitle[i].data());
    h_mjj_rec_truth[i]->SetYTitle(ytitle.data());
    h_mjj_rec_truth[i]->SetXTitle("M_{qq} [GeV/c^{2}]");

  }

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

  for(int ip=0; ip<2; ip++){
   
    // ratio
    pf_dR_Rpt_gen[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_Rpt_gen%d",ip));
    pf_dR_Rpt_gen[ip]->SetYTitle(Form("p_{T}(genJet)/p_{T}(q_{%d}) [GeV]",ip+1));
  
    pf_dR_Rpt_rec[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_Rpt_rec%d",ip));
    pf_dR_Rpt_rec[ip]->SetYTitle(Form("p_{T}(PFJet)/p_{T}(q_{%d}) [GeV]",ip+1));

  }
 
  int nPass[50]={0};
  
  Long64_t nentries = fChain->GetEntriesFast();
  standalone_LumiReWeighting LumiWeights_central(2012,0);
  signalShapeReWeighting 
    signalWeightMethod("/home/syu/old/scripts/hzz2l2q/data/mZZ_Higgs900_8TeV_Lineshape+Interference.txt");

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
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

    for(unsigned int igen=0; igen < genParId_->size(); igen++)
      {	
	int status = genParSt_->at(igen);
	int PID    = genParId_->at(igen);
	// higgs
	if(status==3 && PID==25)
	  {
	    genHiggsMass = genParM_->at(igen);
	    h_mh_parton_mother->Fill(genHiggsMass);
	    if(reweightSignal)
	      signal_weight = signalWeightMethod.weight(genHiggsMass,0);
	    h_mh_parton_mother_weighted->Fill(genHiggsMass,signal_weight);
	    break;
	  }
      }

    if(lep1Index<0 || lep2Index<0)continue;
    if(lep1PostIndex<0 || lep2PostIndex<0)continue;
    if(q1Index<0 || q2Index<0)continue;

    nPass[1]++;

    double PU_weight =  LumiWeights_central.weight(PU_nTrueInt);
    double rec_weight = PU_weight*signal_weight;
    double mH_parton = (lep1+lep2+q1+q2).M();
    double mZll_parton = (lep1+lep2).M();
    double mZjj_parton = (q1+q2).M();
    double dR_parton = q1.DeltaR(q2);

    h_mh_parton_daughter->Fill(mH_parton, signal_weight);
    h_mll_parton_daughter->Fill(mZll_parton, signal_weight);
    h_mjj_parton_daughter->Fill(mZjj_parton, signal_weight);

    nPass[2]++;


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
    
  
    //================================================================
    // STATUS =1       LEVEL
    //================================================================


    double lep1posteta = lep1_post.Eta();
    double lep1postphi = lep1_post.Phi();

    double lep2posteta = lep2_post.Eta();
    double lep2postphi = lep2_post.Phi();

    //================================================================
    // Starting from genJets, more like real data
    //================================================================

    for(unsigned int ijet =0; ijet < genJetPt_->size(); ijet++)
      {
	double etai = genJetEta_->at(ijet);
	double phii = genJetPhi_->at(ijet);
	double pti  = genJetPt_->at(ijet);

	bool   iMatchedToQ1=false;
	bool   iMatchedToQ2=false;

	if(deltaR(etai,phii, lep1posteta, lep1postphi)<0.5)continue;
	if(deltaR(etai,phii, lep2posteta, lep2postphi)<0.5)continue;
	if(fabs(etai) > 2.4)continue;
 	if(pti < 30.0)continue;

	if(matchGenToParton(q1Index,ijet)) iMatchedToQ1=true;
	if(matchGenToParton(q2Index,ijet)) iMatchedToQ2=true;

	for(unsigned int kjet =0; kjet < ijet; kjet++)
	  {
	    
	    double etak = genJetEta_->at(kjet);
	    double phik = genJetPhi_->at(kjet);
	    double ptk  = genJetPt_->at(kjet);
	    bool   kMatchedToQ1=false;
	    bool   kMatchedToQ2=false;

	    if(deltaR(etak,phik, lep1posteta, lep1postphi)<0.5)continue;
	    if(deltaR(etak,phik, lep2posteta, lep2postphi)<0.5)continue;
	    if(fabs(etak) > 2.4)continue;
	    if(ptk < 30.0)continue;
	  
	    if(matchGenToParton(q1Index,kjet)) kMatchedToQ1=true;
	    if(matchGenToParton(q2Index,kjet)) kMatchedToQ2=true;

	    TLorentzVector l4_ijet(0,0,0,0);
	    
	    l4_ijet.SetPtEtaPhiE(
				 genJetPt_->at(ijet),
				 genJetEta_->at(ijet),
				 genJetPhi_->at(ijet),
				 genJetE_->at(ijet));

	    TLorentzVector l4_kjet(0,0,0,0);
	    
	    l4_kjet.SetPtEtaPhiE(
				 genJetPt_->at(kjet),
				 genJetEta_->at(kjet),
				 genJetPhi_->at(kjet),
				 genJetE_->at(kjet));

	    double mH_random   = (l4_ijet+l4_kjet+lep1_post+lep2_post).M();
	    double mZll_random = (lep1_post + lep2_post).M();
	    double mZjj_random = (l4_ijet+l4_kjet).M();

	    // both matched to different quarks
	    if(  (  iMatchedToQ1 &&  kMatchedToQ2 &&
		    !iMatchedToQ2 && !kMatchedToQ1) || 
		 (  iMatchedToQ2 &&  kMatchedToQ1 &&
		    !iMatchedToQ1 && !kMatchedToQ2)
		 )
	      {
		if(DEBUG==1)
		  {
		    cout << "random index = " << ijet << ", \t" << kjet << endl;
		    cout << "ijet pt = " << l4_ijet.Pt() << ", \t kjet pt = " << l4_kjet.Pt() << endl;
		  }
		h_mh_parton_random[0]->Fill(mH_parton, signal_weight);
		h_mll_parton_random[0]->Fill(mZll_parton, signal_weight);
		h_mjj_parton_random[0]->Fill(mZjj_parton, signal_weight);

		h_mh_stable_random[0]->Fill(mH_random, signal_weight);
		h_mll_stable_random[0]->Fill(mZll_random, signal_weight);
		h_mjj_stable_random[0]->Fill(mZjj_random, signal_weight);

		pf_dR_Rm_gen->Fill(dR_parton, mZjj_random/mZjj_parton
				   , signal_weight);
		pf_dR_dm_gen->Fill(dR_parton, mZjj_random-mZjj_parton
				   , signal_weight);
		

		if(iMatchedToQ1)
		  {
		    pf_dR_Rpt_gen[0]->Fill(dR_parton, l4_ijet.Pt()/q1.Pt()
					   , signal_weight);	 
		    pf_dR_Rpt_gen[1]->Fill(dR_parton, l4_kjet.Pt()/q2.Pt()
					   , signal_weight);	 
		  }
		else		
		  {
		    pf_dR_Rpt_gen[0]->Fill(dR_parton, l4_kjet.Pt()/q1.Pt()
					   , signal_weight);	 
		    pf_dR_Rpt_gen[1]->Fill(dR_parton, l4_ijet.Pt()/q2.Pt()
					   , signal_weight);	 
		  }


	      }

	    else if(  ( iMatchedToQ1 &&  iMatchedToQ2 &&
			!kMatchedToQ1 && !kMatchedToQ2) || 
		      ( kMatchedToQ1 &&  kMatchedToQ2 &&
			!iMatchedToQ1 && !iMatchedToQ2)
		      )
	      {
		h_mh_parton_random[2]->Fill(mH_parton, signal_weight);
		h_mll_parton_random[2]->Fill(mZll_parton, signal_weight);
		h_mjj_parton_random[2]->Fill(mZjj_parton, signal_weight);

		h_mh_stable_random[2]->Fill(mH_random, signal_weight);
		h_mll_stable_random[2]->Fill(mZll_random, signal_weight);
		h_mjj_stable_random[2]->Fill(mZjj_random, signal_weight);
	      }

	    else if(  ( iMatchedToQ1 && !iMatchedToQ2 &&
			!kMatchedToQ1 && !kMatchedToQ2) || 

		      (!iMatchedToQ1 &&  iMatchedToQ2 &&
		       !kMatchedToQ1 && !kMatchedToQ2) || 

		      (!iMatchedToQ1 && !iMatchedToQ2 &&
		       kMatchedToQ1 && !kMatchedToQ2) || 

		      (!iMatchedToQ1 && !iMatchedToQ2 &&
		       !kMatchedToQ1 &&  kMatchedToQ2) 
		      )
	      {

		h_mh_parton_random[4]->Fill(mH_parton, signal_weight);
		h_mll_parton_random[4]->Fill(mZll_parton, signal_weight);
		h_mjj_parton_random[4]->Fill(mZjj_parton, signal_weight);

		h_mh_stable_random[4]->Fill(mH_random, signal_weight);
		h_mll_stable_random[4]->Fill(mZll_random, signal_weight);
		h_mjj_stable_random[4]->Fill(mZjj_random, signal_weight);
	      }

	    else if(  !iMatchedToQ1 && !iMatchedToQ2 &&
		      !kMatchedToQ1 && !kMatchedToQ2
		      )
	      {
		h_mh_parton_random[5]->Fill(mH_parton, signal_weight);
		h_mll_parton_random[5]->Fill(mZll_parton, signal_weight);
		h_mjj_parton_random[5]->Fill(mZjj_parton, signal_weight);

		h_mh_stable_random[5]->Fill(mH_random, signal_weight);
		h_mll_stable_random[5]->Fill(mZll_random, signal_weight);
		h_mjj_stable_random[5]->Fill(mZjj_random, signal_weight);
	      }



	  } // end of loop over kjet

      }// end of loop over ijet

    
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
	
	nGoodJets++;

	if(matchGenToParton(q1Index,ijet) && jet1.E()<1e-6)
	  {
	    jet1.SetPtEtaPhiE(
			      genJetPt_->at(ijet),
			      genJetEta_->at(ijet),
			      genJetPhi_->at(ijet),
			      genJetE_->at(ijet));
	    jet1Index = ijet;
	  }

	if(matchGenToParton(q2Index,ijet) && jet2.E()<1e-6)
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
      

    // ------------------------------------------------
    // Reconstruction level
    //-------------------------------------------------    

    for(int ih=0; ih<higgsM->size(); ih++){

      bool matchedToQ1=false;
      bool matchedToQ2=false;
      
      int nMatch=0;
      TLorentzVector mergedJet_l4(0,0,0,0);
      double jetRecPt[2]={0};

      for(int ijet=0; ijet<jetPt->size(); ijet++){

	if(jetHiggsIndex->at(ijet)!=ih)continue;
	if(jetPt->at(ijet)<30.0)continue;
	if(fabs(jetEta->at(ijet))>2.4)continue;

	int jet_index = jetIndex->at(ijet);
	if(jet_index<0)continue;
	if(jet_index>1)continue;
	
	matchedToQ1=matchGenToPFJet(q1Index, ijet);
	matchedToQ2=matchGenToPFJet(q2Index, ijet);


	if(matchedToQ1 && matchedToQ2)
	  {
	    nMatch=100;
	    mergedJet_l4.SetPtEtaPhiE(
				      jetPt->at(ijet),
				      jetEta->at(ijet),
				      jetPhi->at(ijet),
				      jetE->at(ijet)				      
				      );	    
	  }
	else if(matchedToQ1 && !matchedToQ2)
	  {
	    nMatch++;
	    jetRecPt[0] = jetPt->at(ijet);
	  }
	else if(matchedToQ2 && !matchedToQ1)
	  {
	    nMatch++;
	    jetRecPt[1] = jetPt->at(ijet);
	  }
      } // end of loop over jets
      
      double mH_rec  = higgsM->at(ih);
      double mll_rec = zllM->at(ih);
      double mjj_rec = zjjM->at(ih);

      if(nMatch==2){
	h_mh_rec_random[0]->Fill(mH_rec, rec_weight);
	h_mll_rec_random[0]->Fill(mll_rec, rec_weight);
	h_mjj_rec_random[0]->Fill(mjj_rec, rec_weight);

	h_mh_rec_truth[0]->Fill(mH_parton, rec_weight);
	h_mll_rec_truth[0]->Fill(mZll_parton, rec_weight);
	h_mjj_rec_truth[0]->Fill(mZjj_parton, rec_weight);
	
	if(jet1Index>=0 && jet2Index>=0 && jet1Index!=jet2Index){       

	  h_mh_rec_random[1]->Fill(mH_rec, rec_weight);
	  h_mll_rec_random[1]->Fill(mll_rec, rec_weight);
	  h_mjj_rec_random[1]->Fill(mjj_rec, rec_weight);

	  h_mh_rec_truth[1]->Fill(mH_parton, rec_weight);
	  h_mll_rec_truth[1]->Fill(mZll_parton, rec_weight);
	  h_mjj_rec_truth[1]->Fill(mZjj_parton, rec_weight);

	// Roberto-style profiles
	  if(mjj_rec > LOOSE_MIN_MZ_JJ && mjj_rec < LOOSE_MAX_MZ_JJ){
	    pf_dR_Rm_rec->Fill(dR_parton, mjj_rec/mZjj_parton, rec_weight);
	    pf_dR_dm_rec->Fill(dR_parton, mjj_rec-mZjj_parton, rec_weight);
	  for(int ieiko=0; ieiko<2; ieiko++)
	    {
	      pf_dR_Rpt_rec[ieiko]->Fill(dR_parton, jetRecPt[ieiko]/QuarkPt[ieiko],rec_weight);
	    }
	  }
	  
	} // if for the same event, gen jets are also matched to quarks             

	  
      } // if both jets are matched

      else if(nMatch==100){
	h_mh_rec_random[2]->Fill(mH_rec, rec_weight);
	h_mll_rec_random[2]->Fill(mll_rec, rec_weight);
	h_mjj_rec_random[2]->Fill(mjj_rec, rec_weight);

	h_mh_rec_truth[2]->Fill(mH_parton, rec_weight);
	h_mll_rec_truth[2]->Fill(mZll_parton, rec_weight);
	h_mjj_rec_truth[2]->Fill(mZjj_parton, rec_weight);

	if(mergedJet_l4.E()>1e-6)
	  h_mjj_rec_random[3]->Fill(mergedJet_l4.M(), rec_weight);

	h_mh_rec_truth[3]->Fill(mH_parton, rec_weight);
	h_mll_rec_truth[3]->Fill(mZll_parton, rec_weight);
	h_mjj_rec_truth[3]->Fill(mZjj_parton, rec_weight);


      }

      else if(nMatch==1){
	h_mh_rec_random[4]->Fill(mH_rec, rec_weight);
	h_mll_rec_random[4]->Fill(mll_rec, rec_weight);
	h_mjj_rec_random[4]->Fill(mjj_rec, rec_weight);

	h_mh_rec_truth[4]->Fill(mH_parton, rec_weight);
	h_mll_rec_truth[4]->Fill(mZll_parton, rec_weight);
	h_mjj_rec_truth[4]->Fill(mZjj_parton, rec_weight);

      }
      else if(nMatch==0){
	h_mh_rec_random[5]->Fill(mH_rec, rec_weight);
	h_mll_rec_random[5]->Fill(mll_rec, rec_weight);
	h_mjj_rec_random[5]->Fill(mjj_rec, rec_weight);

	h_mh_rec_truth[5]->Fill(mH_parton, rec_weight);
	h_mll_rec_truth[5]->Fill(mZll_parton, rec_weight);
	h_mjj_rec_truth[5]->Fill(mZjj_parton, rec_weight);

      }
      
    } // end of loop over higgs candidate



    //================================================================
    //    Event-level histograms
    //================================================================


    if(jet1Index>=0 && jet2Index>=0 && jet1Index!=jet2Index){

      h_mh_parton_event[0]->Fill(mH_parton, signal_weight);
      h_mll_parton_event[0]->Fill(mZll_parton, signal_weight);
      h_mjj_parton_event[0]->Fill(mZjj_parton, signal_weight);

      double mH_particle = (jet1+jet2+lep1_post+lep2_post).M();
      double mZll_particle = (lep1_post + lep2_post).M();
      double mZjj_particle = (jet1+jet2).M();

      h_mh_stable_event[0]->Fill(mH_particle, signal_weight);
      h_mll_stable_event[0]->Fill(mZll_particle, signal_weight);
      h_mjj_stable_event[0]->Fill(mZjj_particle, signal_weight);

      if(higgsM->size()>0)nPass[4]++;

    }
      
    else if(jet1Index>=0 && jet2Index>=0 && jet1Index==jet2Index && nGoodJets==1){

      h_mh_parton_event[1]->Fill(mH_parton, signal_weight);
      h_mll_parton_event[1]->Fill(mZll_parton, signal_weight);
      h_mjj_parton_event[1]->Fill(mZjj_parton, signal_weight);


      double mH_particle = (jet1+lep1_post+lep2_post).M();
      double mZll_particle = (lep1_post + lep2_post).M();
      double mZjj_particle = jet1.M();

      h_mh_stable_event[1]->Fill(mH_particle, signal_weight);
      h_mll_stable_event[1]->Fill(mZll_particle, signal_weight);
      h_mjj_stable_event[1]->Fill(mZjj_particle, signal_weight);

    }

    else if(jet1Index>=0 && jet2Index>=0 && jet1Index==jet2Index && nGoodJets>1){

      h_mh_parton_event[2]->Fill(mH_parton, signal_weight);
      h_mll_parton_event[2]->Fill(mZll_parton, signal_weight);
      h_mjj_parton_event[2]->Fill(mZjj_parton, signal_weight);


      double mH_particle = (jet1+lep1_post+lep2_post).M();
      double mZll_particle = (lep1_post + lep2_post).M();
      double mZjj_particle = jet1.M();

      h_mh_stable_event[2]->Fill(mH_particle, signal_weight);
      h_mll_stable_event[2]->Fill(mZll_particle, signal_weight);
      h_mjj_stable_event[2]->Fill(mZjj_particle, signal_weight);

    }
      
    else if(jet1Index>=0 && jet2Index <0 && nGoodJets==1){

      h_mh_parton_event[3]->Fill(mH_parton, signal_weight);
      h_mll_parton_event[3]->Fill(mZll_parton, signal_weight);
      h_mjj_parton_event[3]->Fill(mZjj_parton, signal_weight);

      double mH_particle = (jet1+lep1_post+lep2_post).M();
      double mZll_particle = (lep1_post + lep2_post).M();
      double mZjj_particle = jet1.M();

      h_mh_stable_event[3]->Fill(mH_particle, signal_weight);
      h_mll_stable_event[3]->Fill(mZll_particle, signal_weight);
      h_mjj_stable_event[3]->Fill(mZjj_particle, signal_weight);

    }
      
    else if(jet2Index>=0 && jet1Index <0 && nGoodJets==1){

      h_mh_parton_event[3]->Fill(mH_parton, signal_weight);
      h_mll_parton_event[3]->Fill(mZll_parton, signal_weight);
      h_mjj_parton_event[3]->Fill(mZjj_parton, signal_weight);

      double mH_particle = (jet2+lep1_post+lep2_post).M();
      double mZll_particle = (lep1_post + lep2_post).M();
      double mZjj_particle = jet2.M();

      h_mh_stable_event[3]->Fill(mH_particle, signal_weight);
      h_mll_stable_event[3]->Fill(mZll_particle, signal_weight);
      h_mjj_stable_event[3]->Fill(mZjj_particle, signal_weight);

    }
    else if(jet1Index>=0 && jet2Index <0 && nGoodJets >1){

      h_mh_parton_event[4]->Fill(mH_parton, signal_weight);
      h_mll_parton_event[4]->Fill(mZll_parton, signal_weight);
      h_mjj_parton_event[4]->Fill(mZjj_parton, signal_weight);

      double mH_particle = (jet1+lep1_post+lep2_post).M();
      double mZll_particle = (lep1_post + lep2_post).M();
      double mZjj_particle = jet1.M();

      h_mh_stable_event[4]->Fill(mH_particle, signal_weight);
      h_mll_stable_event[4]->Fill(mZll_particle, signal_weight);
      h_mjj_stable_event[4]->Fill(mZjj_particle, signal_weight);

    }
      
    else if(jet2Index>=0 && jet1Index <0 && nGoodJets >1){

      h_mh_parton_event[4]->Fill(mH_parton, signal_weight);
      h_mll_parton_event[4]->Fill(mZll_parton, signal_weight);
      h_mjj_parton_event[4]->Fill(mZjj_parton, signal_weight);

      double mH_particle = (jet2+lep1_post+lep2_post).M();
      double mZll_particle = (lep1_post + lep2_post).M();
      double mZjj_particle = jet2.M();

      h_mh_stable_event[4]->Fill(mH_particle, signal_weight);
      h_mll_stable_event[4]->Fill(mZll_particle, signal_weight);
      h_mjj_stable_event[4]->Fill(mZjj_particle, signal_weight);

    }
         
      
    else continue;

    nPass[3]++;


  } // end of loop over entries

  for(int i=0; i<50;i++)
    if(nPass[i]>0)
      cout << "nPass[" << i << "] = " << nPass[i] << endl;

  TFile* outFile = new TFile(Form("studymjj_%s",
				  _inputFileName.data()),
			     "recreate");            

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
    
  }

  for(int i=0; i<NCASES; i++){
    
    h_mh_parton_event[i]->Write();
    h_mll_parton_event[i]->Write();
    h_mjj_parton_event[i]->Write();

    h_mh_stable_event[i]->Write();
    h_mll_stable_event[i]->Write();
    h_mjj_stable_event[i]->Write();

    h_mh_parton_random[i]->Write();
    h_mll_parton_random[i]->Write();
    h_mjj_parton_random[i]->Write();

    h_mh_stable_random[i]->Write();
    h_mll_stable_random[i]->Write();
    h_mjj_stable_random[i]->Write();

    h_mh_rec_random[i]->Write();
    h_mll_rec_random[i]->Write();
    h_mjj_rec_random[i]->Write();

    h_mh_rec_truth[i]->Write();
    h_mll_rec_truth[i]->Write();
    h_mjj_rec_truth[i]->Write();


  }

  outFile->Close();  

}


Bool_t mjj_why::matchGenToParton(Int_t igen, Int_t ijet){

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

Bool_t mjj_why::matchGenToPFJet(Int_t igen, Int_t ijet){

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
