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
  const int nCases=4;
  std::string matchName[nCases]={"none","quarks from Z","associated parton","both"};
  TH1D* h_matchCase_template = new TH1D("h_matchCase_template","",nCases,-0.5,nCases-0.5);
  h_matchCase_template->Sumw2();
  for(int ib=1; ib<=nCases;ib++)
    h_matchCase_template->GetXaxis()->SetBinLabel(ib,matchName[ib-1].data());
  TH1D* h_pt_template = new TH1D("h_pt_template","",100,0,500);
  h_pt_template->Sumw2();
  TH1D* h_eta_template = new TH1D("h_eta_template","",100,-5.0,5.0);
  h_eta_template->Sumw2();

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

  //   TProfile* pf_dR_template = new TProfile("pf_dR_template","",50,0.0,5.0,-1000, 1000);
  TProfile* pf_dR_template = new TProfile("pf_dR_template","",17,0.0,3.4,-1000, 1000);
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

  TH1D* h_match[NBINS][2];

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

      h_match[ibin][ip] = 
	(TH1D*)h_matchCase_template->Clone(Form("h_match%02i_%d",ibin,ip));
      h_match[ibin][ip]->SetXTitle(Form("jet_{%d} %s", ip+1, dRstring.data()));

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
  
  TH1D* h_mh_stable;
  TH1D* h_mll_stable;
  TH1D* h_mjj_stable;
 
  TH1D* h_mh_rec;
  TH1D* h_mll_rec;
  TH1D* h_mjj_rec;

  
  h_mh_stable = (TH1D*)h_mh_template->Clone("h_mh_stable");
  h_mh_rec = (TH1D*)h_mh_template->Clone("h_mh_rec");

  h_mll_stable = (TH1D*)h_mll_template->Clone("h_mll_stable");
  h_mll_rec = (TH1D*)h_mll_template->Clone("h_mll_rec");

  h_mjj_stable = (TH1D*)h_mjj_template->Clone("h_mjj_stable");
  h_mjj_rec = (TH1D*)h_mjj_template->Clone("h_mjj_rec");


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

  // generator-level
  TProfile* pf_dR_matchProb_gen[2]; // matched to the q and qbar from Z
  TProfile* pf_dR_matchProb_gen_ass[2]; // matched to the associated parton produced with Higgs
  TProfile* pf_dR_matchProb_gen_either[2]; // matched to either q or qbar from Z or the associated parton produced with Higgs

  // after requiring mjj or mj cut
  TProfile* pf_dR_matchProb_gen_afterMjj[2]; // matched to the q and qbar from Z
  TProfile* pf_dR_matchProb_gen_afterMjj_ass[2]; // matched to the associated parton produced with Higgs
  TProfile* pf_dR_matchProb_gen_afterMjj_either[2]; // matched to either q or qbar from Z or the associated parton produced with Higgs

  // reconstruction level from Higgs candidate
  TProfile* pf_dR_matchProb_h[2]; // matched to the q and qbar from Z
  TProfile* pf_dR_matchProb_h_ass[2]; // matched to the associated parton produced with Higgs
  TProfile* pf_dR_matchProb_h_either[2]; // matched to either q or qbar from Z or the associated parton produced with Higgs
  
  // reconstruction level for leading jets
  TProfile* pf_dR_matchProb_rec[2]; // matched to the q and qbar from Z
  TProfile* pf_dR_matchProb_rec_ass[2]; // matched to the associated parton produced with Higgs
  TProfile* pf_dR_matchProb_rec_either[2]; // matched to either q or qbar from Z or the associated parton produced with Higgs

  // after requiring mjj or mj cut
  TProfile* pf_dR_matchProb_rec_afterMjj[2]; // matched to the q and qbar from Z
  TProfile* pf_dR_matchProb_rec_afterMjj_ass[2]; // matched to the associated parton produced with Higgs
  TProfile* pf_dR_matchProb_rec_afterMjj_either[2]; // matched to either q or qb
  for(int ip=0; ip<2; ip++){
   
    // ratio
    pf_dR_Rpt_gen[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_Rpt_gen%d",ip));
    pf_dR_Rpt_gen[ip]->SetYTitle(Form("p_{T}(genJet_{%d})/p_{T}(q_{%d}) [GeV]",ip+1,ip+1));
  
    pf_dR_Rpt_rec[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_Rpt_rec%d",ip));
    pf_dR_Rpt_rec[ip]->SetYTitle(Form("p_{T}(PFJet_{%d})/p_{T}(q_{%d}) [GeV]",ip+1,ip+1));

    pf_dR_jetArea_rec[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_jetArea_rec%d",ip));
    pf_dR_jetArea_rec[ip] -> SetYTitle(Form("jet area for PFJet %d",ip+1));

    pf_dR_matchProb_gen[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_matchProb_gen%d",ip));
    pf_dR_matchProb_gen[ip] -> SetYTitle(Form("Probability of matching to q from Z for genJet %d",ip+1));

    pf_dR_matchProb_gen_ass[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_matchProb_gen_ass%d",ip));
    pf_dR_matchProb_gen_ass[ip] -> SetYTitle(Form("Probability of matching to associated parton for genJet %d",ip+1));
    pf_dR_matchProb_gen_either[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_matchProb_gen_either%d",ip));
    pf_dR_matchProb_gen_either[ip] -> SetYTitle(Form("Probability of matching to either for genJet %d",ip+1));

    pf_dR_matchProb_gen_afterMjj[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_matchProb_gen_afterMjj%d",ip));
    pf_dR_matchProb_gen_afterMjj[ip] -> SetTitle("M_{j} or M_{jj} within Z mass window");
    pf_dR_matchProb_gen_afterMjj[ip] -> SetYTitle(Form("Probability of matching to q from Z for genJet %d",ip+1));

    pf_dR_matchProb_gen_afterMjj_ass[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_matchProb_gen_afterMjj_ass%d",ip));
    pf_dR_matchProb_gen_afterMjj_ass[ip] -> SetYTitle(Form("Probability of matching to associated parton for genJet %d",ip+1));
    pf_dR_matchProb_gen_afterMjj_either[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_matchProb_gen_afterMjj_either%d",ip));
    pf_dR_matchProb_gen_afterMjj_either[ip] -> SetYTitle(Form("Probability of matching to either for genJet %d",ip+1));

    /////////////////////////////////////////////////////////

    pf_dR_matchProb_h[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_matchProb_h%d",ip));
    pf_dR_matchProb_h[ip] -> SetYTitle(Form("Probability of matching to q from Z for PFJet %d",ip+1));

    pf_dR_matchProb_h_ass[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_matchProb_h_ass%d",ip));
    pf_dR_matchProb_h_ass[ip] -> SetYTitle(Form("Probability of matching to associated parton for PFJet %d",ip+1));
    pf_dR_matchProb_h_either[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_matchProb_h_either%d",ip));
    pf_dR_matchProb_h_either[ip] -> SetYTitle(Form("Probability of matching to either for PFJet %d",ip+1));

    /////////////////////////////////////////////////////////

    pf_dR_matchProb_rec[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_matchProb_rec%d",ip));
    pf_dR_matchProb_rec[ip] -> SetYTitle(Form("Probability of matching to q from Z for PFJet %d",ip+1));

    pf_dR_matchProb_rec_ass[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_matchProb_rec_ass%d",ip));
    pf_dR_matchProb_rec_ass[ip] -> SetYTitle(Form("Probability of matching to associated parton for PFJet %d",ip+1));
    pf_dR_matchProb_rec_either[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_matchProb_rec_either%d",ip));
    pf_dR_matchProb_rec_either[ip] -> SetYTitle(Form("Probability of matching to either for PFJet %d",ip+1));

    pf_dR_matchProb_rec_afterMjj[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_matchProb_rec_afterMjj%d",ip));
    pf_dR_matchProb_rec_afterMjj[ip] -> SetTitle("M_{j} or M_{jj} within Z mass window");
    pf_dR_matchProb_rec_afterMjj[ip] -> SetYTitle(Form("Probability of matching to q from Z for PFJet %d",ip+1));

    pf_dR_matchProb_rec_afterMjj_ass[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_matchProb_rec_afterMjj_ass%d",ip));
    pf_dR_matchProb_rec_afterMjj_ass[ip] -> SetYTitle(Form("Probability of matching to associated parton for PFJet %d",ip+1));
    pf_dR_matchProb_rec_afterMjj_either[ip] = (TProfile*)pf_dR_template->Clone(Form("pf_dR_matchProb_rec_afterMjj_either%d",ip));
    pf_dR_matchProb_rec_afterMjj_either[ip] -> SetYTitle(Form("Probability of matching to either for PFJet %d",ip+1));

  }
 
  int nPass[50]={0};
  
  Long64_t nentries = fChain->GetEntriesFast();
  standalone_LumiReWeighting LumiWeights_central(2012,0);
  std::string reweight_file = Form("/home/syu/cvs_scripts/hzz2l2q/data/mZZ_Higgs%d_8TeV_Lineshape+Interference.txt",_higgsMass);
  if(reweighSignal)
    cout << "reweighing signal with the input file: " << reweight_file << endl;
  signalShapeReWeighting signalWeightMethod(reweight_file.data());

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    //     if(jentry > 1000)break;
    nPass[0]++;
    //===================================================================
    // Find the quark index of particles that come from Z hadronic decays
    //===================================================================


    int lep1Index=-1;
    int lep2Index=-1;

    int lep1PostIndex=-1;
    int lep2PostIndex=-1;

    int q1Index=-1;
    int q2Index=-1;

    const int ZPID=23;
    int LEPTYPE = -1;
    int lep1PID = -1;
    int lep1PostPID = -1;
    int q1PID = -999;

    double genHiggsMass = 0.0;
    double signal_weight= 1.0;

    int higgsMo1=-1;
    int higgsMo2=-1;
    
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
	    if(reweighSignal && _higgsMass >= 400 && _higgsMass <= 1000)
	      signal_weight = signalWeightMethod.weight(genHiggsMass,0);
	    h_mh_parton_mother_weighted->Fill(genHiggsMass,signal_weight);

	    higgsMo1 = genMo1_->at(igen);
	    higgsMo2 = genMo2_->at(igen);
	  }

	// first status=3 lepton from Z
	if(status==3 && lep1Index < 0 && (abs(PID) == 11 ||
					  abs(PID) == 13)
	   && momPID == ZPID
	   )
	  {
	    LEPTYPE = abs(PID);
	    lep1Index = igen;
	    lep1PID   = PID;
	  }

	// second status=3 lepton from Z
	else if (status==3 && lep2Index < 0 && (abs(PID) == 11 ||
						abs(PID) == 13)
		 && momPID == ZPID		 
		 && PID == -(lep1PID)
		 )
	  {
	    lep2Index = igen;
	  }

	// first status =1 lepton from Z
	else if(status==1 && lep1PostIndex < 0 && (abs(PID) == 11 ||
						   abs(PID) == 13)
		&& abs(momPID) == LEPTYPE
		&& abs(PID) == LEPTYPE
		)
	  {
	    lep1PostIndex = igen;
	    lep1PostPID = PID;
	  }

	// second status=1 lepton from Z
	else if (status==1 && lep2PostIndex < 0 && (abs(PID) == 11 ||
						    abs(PID) == 13)
		 && abs(momPID) == LEPTYPE
		 && PID == -(lep1PostPID)
		 )
	  {
	    lep2PostIndex = igen;
	  }

	// first parton from Z
	else if(status==3 && q1Index < 0 && abs(PID) < 6 
		&& momPID == ZPID)
	  {
	    q1Index = igen;
	    q1PID   = PID;
	  }

	// second parton from Z
	else if(status==3 && q2Index < 0 && abs(PID) < 6 
		&& momPID == ZPID 
		&& PID == -(q1PID)
		)
	  {
	    q2Index = igen;
	  }

			    
      } // end of loop over gen particles
    
    //===================================================================    
    // look for parton that is produced with Higgs
    //===================================================================

    int assPartonIndex = -1;
    for(unsigned int igen=0; igen < genParId_->size(); igen++)
      {
	
	int status = genParSt_->at(igen);
	int PID    = genParId_->at(igen);
	int mom1   = genMo1_->at(igen);
	int mom2   = genMo2_->at(igen);

	if(status!=3)continue;
	if(PID!=21 && (abs(PID)>=6 || abs(PID)<1))continue;
	if(mom1 < 0 && mom2 < 0)continue;
	if(mom1 != higgsMo1)continue;
	if(mom2 != higgsMo2)continue;

	assPartonIndex = igen;
	if(assPartonIndex >=0)break;
      }
 
    if(DEBUG==1){
      cout << "===========================================" << endl;
      cout << "Event = " << EvtInfo_EventNum << endl;
      cout << "LEPTYPE = " << LEPTYPE << endl;
      cout << "lep1 = " << lep1Index << "\t lep2 = " << lep2Index << endl;
      cout << "lep1Post =" << lep1PostIndex << "\t lep2Post = " << lep2PostIndex << endl;
      cout << "q1 =" << q1Index << "\t q2Index = " << q2Index << endl;
      cout << "assParton = " << assPartonIndex << endl;
      cout << "===========================================" << endl;
    }
    if(LEPTYPE==-1)continue;

    if(lep1Index<0 || lep2Index<0)continue;
    if(lep1PostIndex<0 || lep2PostIndex<0)continue;
    nPass[1]++;

    if(q1Index<0 || q2Index<0)continue;
    nPass[2]++;


    //===================================================================    
    // Setting the LorentzVectors of the particles from Higgs decay
    //===================================================================

    // setting TLorentzVector of generator-level particles
    // leptons from Z decays
    TLorentzVector lep1(0,0,0,0);
    lep1.SetPtEtaPhiM(
		      genParPt_->at(lep1Index),
		      genParEta_->at(lep1Index),
		      genParPhi_->at(lep1Index),
		      genParM_->at(lep1Index)
		      );

    TLorentzVector lep2(0,0,0,0);
    lep2.SetPtEtaPhiM(
		      genParPt_->at(lep2Index),
		      genParEta_->at(lep2Index),
		      genParPhi_->at(lep2Index),
		      genParM_->at(lep2Index)
		      );

    if(DEBUG==1)
      {
	cout << "Printing information of leptonic Z before FSR" << endl;
	lep1.Print();
	lep2.Print();
      }


    // leptons after FSR from Z decays
    TLorentzVector lep1_post(0,0,0,0);
    lep1_post.SetPtEtaPhiM(
			   genParPt_->at(lep1PostIndex),
			   genParEta_->at(lep1PostIndex),
			   genParPhi_->at(lep1PostIndex),
			   genParM_->at(lep1PostIndex)
			   );
    TLorentzVector lep2_post(0,0,0,0);
    lep2_post.SetPtEtaPhiM(
			   genParPt_->at(lep2PostIndex),
			   genParEta_->at(lep2PostIndex),
			   genParPhi_->at(lep2PostIndex),
			   genParM_->at(lep2PostIndex)
			   );

    if(DEBUG==1)
      {
	cout << "Printing information of leptonic Z after FSR" << endl;
	lep1_post.Print();
	lep2_post.Print();
      }

    
    // quark and anti-quark from Z decays
    TLorentzVector q1(0,0,0,0);
    q1.SetPtEtaPhiM(
		    genParPt_->at(q1Index),
		    genParEta_->at(q1Index),
		    genParPhi_->at(q1Index),
		    genParM_->at(q1Index)
		    );	 
    TLorentzVector q2(0,0,0,0);
    q2.SetPtEtaPhiM(
		    genParPt_->at(q2Index),
		    genParEta_->at(q2Index),
		    genParPhi_->at(q2Index),
		    genParM_->at(q2Index)
		    );	 

    if(DEBUG==1)
      {
	cout << "Printing information of hadronic Z" << endl;
	q1.Print();
	q2.Print();
      }


    //===================================================================    
    // Invariant mass of particles with status = 3
    //===================================================================

    double mH_parton = (lep1+lep2+q1+q2).M();
    h_mh_parton_daughter->Fill(mH_parton, signal_weight);

    double mZll_parton = (lep1+lep2).M();
    h_mll_parton_daughter->Fill(mZll_parton, signal_weight);

    double mZjj_parton = (q1+q2).M();
    h_mjj_parton_daughter->Fill(mZjj_parton, signal_weight);

    double dR_parton = q1.DeltaR(q2);
    int dRJJ_index = h_dR->FindBin(dR_parton)-1;
    if(dRJJ_index < 0)dRJJ_index = 0;
    if(dRJJ_index >= NBINS)dRJJ_index = NBINS-1; 
    if(DEBUG==1)
      cout << "dR_parton = " << dR_parton << "\t dRJJ_index = " << 
	dRJJ_index << endl;

    h_dR_qq->Fill(dR_parton,signal_weight);


    //===================================================================    
    // Order the quarks from hadronic decay of Z in pt
    //===================================================================

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
    // STATUS =1       LEVEL invariant mass
    //================================================================


    double lep1posteta = lep1_post.Eta();
    double lep1postphi = lep1_post.Phi();

    double lep2posteta = lep2_post.Eta();
    double lep2postphi = lep2_post.Phi();

    double mZll_particle = (lep1_post + lep2_post).M();
    
    //================================================================
    // Starting from quarks, pick up the generator-level 
    // jets that matched and construct them in the dijetmass
    //================================================================

    TLorentzVector l4_genjet1(0,0,0,0);
    TLorentzVector l4_genjet2(0,0,0,0);
    int jet1Index=-1, jet2Index=-1;

    // store generator-level jets that satisfy kinematic selections
    vector<int> myGoodJets;
    myGoodJets.clear();
    nPass[10]++;
    for(unsigned int ijet =0; ijet < genJetPt_->size(); ijet++)
      {

	double eta = genJetEta_->at(ijet);
	double phi = genJetPhi_->at(ijet);
	double pt  = genJetPt_->at(ijet);
	
	if(deltaR(eta,phi, lep1posteta, lep1postphi)<0.5)continue;
	if(deltaR(eta,phi, lep2posteta, lep2postphi)<0.5)continue;

	if(fabs(eta) > 2.4)continue;
	if(pt < 30.0)continue;

	myGoodJets.push_back(ijet);

      }
    // end of loop over generator-level jets

    if(myGoodJets.size()>0)nPass[11]++;


    // check how often the good leading and second leading generator-level jets
    // are matched to hadronic Z decays and associated produced partons
    int strictJet     = -1; // store the number of jets that are outside of Z mass 
    // window in dijet mass or single-jet mass 
    for(unsigned int aj = 0; aj < myGoodJets.size(); aj++)
      {
	int ijet = myGoodJets[aj];
	int matchType = 0;
	double matchProb_q12   = 0.0;
	double matchProb_ass   = 0.0;
	double matchProb_either= 0.0;

	TLorentzVector l4_templeading(0,0,0,0);
	l4_templeading.SetPtEtaPhiE(
				    genJetPt_->at(ijet),
				    genJetEta_->at(ijet),
				    genJetPhi_->at(ijet),
				    genJetE_->at(ijet));
	
	double singleJetMass = l4_templeading.M();
	bool isMergedZ = singleJetMass > LOOSE_MIN_MZ_JJ 
	  && singleJetMass < MORELOOSE_MAX_MZ_JJ;

	if( matchGenJetToParton(q1Index,ijet) || 
	    matchGenJetToParton(q2Index,ijet)
	    )
	  {
	    matchType += 1;
	    matchProb_q12 = 1.0;
	    matchProb_either = 1.0;
	    if(aj == 0)nPass[12]++;
	  }
	if( matchGenJetToParton(assPartonIndex,ijet))
	  {
	    matchType += 2;
	    matchProb_ass = 1.0;
	    matchProb_either = 1.0;
	    if(aj == 0)nPass[13]++;
	  }	
	
	if(aj<2)
	  {
	    pf_dR_matchProb_gen[aj]->
	      Fill(dR_parton, matchProb_q12, signal_weight);

	    pf_dR_matchProb_gen_ass[aj]->
	      Fill(dR_parton, matchProb_ass, signal_weight);

	    pf_dR_matchProb_gen_either[aj]->
	      Fill(dR_parton, matchProb_either, signal_weight);

	    h_match[dRJJ_index][aj]->
	      Fill(matchType,signal_weight);
	  }


	if(matchGenJetToParton(q1Index,ijet) && l4_genjet1.E()<1e-6)
	  {
	    l4_genjet1   = l4_templeading;
	    jet1Index = ijet;
	  }

	if(matchGenJetToParton(q2Index,ijet) && l4_genjet2.E()<1e-6)
	  {

	    l4_genjet2   = l4_templeading;
	    jet2Index = ijet;
	  }

	// now make a Z mass window cut
	double doubleJetMass = -9999.0;
	bool findAZPair = false;
	for(unsigned int kj=0; kj < myGoodJets.size(); kj++)
	  {
	    if(kj==aj)continue;
	    TLorentzVector l4_tempsub(0,0,0,0);
	    l4_tempsub.SetPtEtaPhiE(
				    genJetPt_->at(myGoodJets[kj]),
				    genJetEta_->at(myGoodJets[kj]),
				    genJetPhi_->at(myGoodJets[kj]),
				    genJetE_->at(myGoodJets[kj]));

	    doubleJetMass = (l4_templeading + l4_tempsub).M();
	    if(doubleJetMass > LOOSE_MIN_MZ_JJ 
	       && doubleJetMass <  MORELOOSE_MAX_MZ_JJ)
	      {
		findAZPair = true;
		break;
	      }
	  } // end of second loop over good jets

	if(!findAZPair && !isMergedZ)continue;

	strictJet ++;

	if(strictJet<2)
	  {
	    pf_dR_matchProb_gen_afterMjj[strictJet]->
	      Fill(dR_parton, matchProb_q12, signal_weight);

	    pf_dR_matchProb_gen_afterMjj_ass[strictJet]->
	      Fill(dR_parton, matchProb_ass, signal_weight);

	    pf_dR_matchProb_gen_afterMjj_either[strictJet]->
	      Fill(dR_parton, matchProb_either, signal_weight);
	  }

      } // end of loop over generator-level jets passing kinematic selections

    
    if(DEBUG==1)
      {
	cout << "jet1Index = " << jet1Index << "\t jet2Index = " << jet2Index
	     << endl;
	cout << "jet1 pt = " << l4_genjet1.Pt() << "\t jet2 pt = " << l4_genjet2.Pt() << endl;
	cout << "q1 pt = " << q1.Pt() << "\t q2 pt = " << q2.Pt() << endl;
	
      }
      
  
    // if two different generator-level leading jets are matched to quarks from Z
    double mZjj_particle = -999;
    
    if(jet1Index>=0 && jet2Index>=0 && jet1Index!=jet2Index){
      nPass[3]++;

      double mH_particle = (l4_genjet1+l4_genjet2+lep1_post+lep2_post).M();
      mZjj_particle = (l4_genjet1+l4_genjet2).M();

      h_mh_stable->Fill(mH_particle, signal_weight);
      h_mll_stable->Fill(mZll_particle, signal_weight);
      h_mjj_stable->Fill(mZjj_particle, signal_weight);

      double rm = mZjj_particle/mZjj_parton;
      double dm = mZjj_particle-mZjj_parton; 
      double rpt[2]={l4_genjet1.Pt()/q1.Pt(),
		     l4_genjet2.Pt()/q2.Pt()};

      // study ratio of mass
      if(mZll_particle > MIN_MZ_LL       && mZll_particle < MAX_MZ_LL 
	 && mZjj_particle > LOOSE_MIN_MZ_JJ && mZjj_particle < MORELOOSE_MAX_MZ_JJ
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
    // Reconstruction level distributions
    //-------------------------------------------------    


    double PU_weight =  LumiWeights_central.weight(PU_nTrueInt);
    double rec_weight = reweighPU? PU_weight*signal_weight: signal_weight;

    // ------------------------------------------------------
    // start from patJets, requiring good ID and beta > 0.2
    //-------------------------------------------------------    

    vector<int> myGoodRecJets;
    myGoodRecJets.clear();

    for(unsigned int ipatjet=0; ipatjet < patJetPfAk05Pt_->size(); ipatjet++){

      if(patJetPfAk05PassID_->at(ipatjet)!=1)continue;

      if(patJetPfAk05Beta_->at(ipatjet)<0.2)continue;

      if(fabs(patJetPfAk05Eta_->at(ipatjet))>2.4)continue;

      if(patJetPfAk05Pt_->at(ipatjet) < 30.0)continue;

      myGoodRecJets.push_back(ipatjet);
    } // end of loop over pat jets


    int strictRecJet= -1;
    for(unsigned int ab=0; ab < myGoodRecJets.size(); ab++){

      int ipatjet = myGoodRecJets[ab];
   
      bool matchedToQ1=false;
      bool matchedToQ2=false;
      bool matchedToAss=false;
  
      matchedToQ1=matchPatJetToParton(q1Index, ipatjet);
      matchedToQ2=matchPatJetToParton(q2Index, ipatjet);
      matchedToAss=matchPatJetToParton(assPartonIndex, ipatjet);


      double matchedProb_zqq =(matchedToQ1 || matchedToQ2)?
	1.0: 0.0;

      double matchedProb_ass = matchedToAss?
	1.0: 0.0;
  
      double matchedProb_either = (matchedToQ1 || matchedToQ2 ||
				   matchedToAss)?
	1.0: 0.0;

      if(ab < 2){
	pf_dR_matchProb_rec[ab]->
	  Fill(dR_parton,matchedProb_zqq, rec_weight);
	 
	pf_dR_matchProb_rec_ass[ab]->
	  Fill(dR_parton,matchedProb_ass, rec_weight);
  
	pf_dR_matchProb_rec_either[ab]->
	  Fill(dR_parton,matchedProb_either, rec_weight);	
      }

      TLorentzVector l4_templeading(0,0,0,0);
      l4_templeading.SetPtEtaPhiE(
  				  patJetPfAk05Pt_->at(ipatjet),
  				  patJetPfAk05Eta_->at(ipatjet),
  				  patJetPfAk05Phi_->at(ipatjet),
  				  patJetPfAk05En_->at(ipatjet)
  				  );
	
      double singleJetMass = l4_templeading.M();
      bool isMergedZ = singleJetMass > LOOSE_MIN_MZ_JJ 
  	&& singleJetMass < MORELOOSE_MAX_MZ_JJ;


      // now make a Z mass window cut
      double doubleJetMass = -9999.0;
      bool findAZPair = false;
      for(unsigned int ey=0; ey < myGoodRecJets.size(); ey++)
 	{
 	  if(ey==ab)continue;
 	  TLorentzVector l4_tempsub(0,0,0,0);
 	  l4_tempsub.SetPtEtaPhiE(
 				  patJetPfAk05Pt_->at(myGoodRecJets[ey]),
 				  patJetPfAk05Eta_->at(myGoodRecJets[ey]),
 				  patJetPfAk05Phi_->at(myGoodRecJets[ey]),
 				  patJetPfAk05En_->at(myGoodRecJets[ey])
				  );

 	  doubleJetMass = (l4_templeading + l4_tempsub).M();
 	  if(doubleJetMass > LOOSE_MIN_MZ_JJ 
 	     && doubleJetMass <  MORELOOSE_MAX_MZ_JJ)
 	    {
 	      findAZPair = true;
 	      break;
 	    }
 	} // end of second loop over good reconstructed jets

      if(!findAZPair && !isMergedZ)continue;

      strictRecJet ++;
      // only fill those passing dijet or single jet Z mass selection
      if(strictRecJet<2)
	{
	  pf_dR_matchProb_rec_afterMjj[strictRecJet]->
	    Fill(dR_parton, matchedProb_zqq, rec_weight);

	  pf_dR_matchProb_rec_afterMjj_ass[strictRecJet]->
	    Fill(dR_parton, matchedProb_ass, rec_weight);

	  pf_dR_matchProb_rec_afterMjj_either[strictRecJet]->
	    Fill(dR_parton, matchedProb_either, rec_weight);
	}
    
	
    } // end of loop over good reconstructed jets
    

    //==========================================================================
    // Instead of picking the leading jets, pick the ones that give best 
    // combination after MET, dijetmass, b-tag, and helicity selections
    //==========================================================================
    int NBTAGMAX = -1;
    int myBest = -1;
    double best_mZll = 9999999.0;
    double best_mZjj = 9999999.0;

    for(unsigned int ih=0; ih < lepType->size(); ih++){
    
      int bitmap = passBit->at(ih);

      bool Pass=false;
    
      if((bitmap & PFMET_SIG) &&
	 (bitmap & HELI_LD)
	 )
	Pass=true;
      if(!Pass)continue;
	 
      double zllMass = zllM->at(ih);
      double zjjMass = zjjM->at(ih);

      if(zllMass < MIN_MZ_LL || zllMass > MAX_MZ_LL)continue; 
      if(zjjMass < LOOSE_MIN_MZ_JJ || zjjMass > LOOSE_MAX_MZ_JJ)continue;

      int nbtag = nBTags->at(ih);

      if(nbtag > NBTAGMAX)
	{
	  myBest    = ih;
	  best_mZjj = zjjMass;
	  best_mZll = zllMass;
	  NBTAGMAX  = nbtag;
	}
      else if(nbtag == NBTAGMAX)
	{
	  if( ( fabs(zjjMass - MZ_PDG) + fabs(zllMass-MZ_PDG) ) < 
	      ( fabs(best_mZjj - MZ_PDG)+fabs(best_mZll-MZ_PDG))
	      )
	    {
	      myBest    = ih;
	      best_mZjj = zjjMass;
	      best_mZll = zllMass;
	      NBTAGMAX  = nbtag;	      
	    }
	}      
    } // loop over candidates
    
    nPass[5]++;
    if(myBest<0)continue;
    if(NBTAGMAX<0)continue;
    nPass[6]++;
    
    if(DEBUG==1)cout << "myBest higgs index = " << myBest << " and NBTAG = " << NBTAGMAX << 
      "\t dijet mass = " << zjjM->at(myBest) << endl;
    
    // once we find the best candidate, starting checking thet jet pt/area/matching probability

    double jetRecPt[2]={0};
      
    int pfJetIndex[2]={-1,-1};

    for(unsigned int ijet=0; ijet<jetPt->size(); ijet++){

      if(jetHiggsIndex_->at(ijet)!= myBest)continue;
      if(jetPt->at(ijet)<30.0)continue;
      if(fabs(jetEta->at(ijet))>2.4)continue;

      int jet_index = jetIndex->at(ijet);
      if(jet_index<0)continue;
      if(jet_index>1)continue;

      bool matchedToQ1=false;
      bool matchedToQ2=false;
      bool matchedToAss=false;
    
      matchedToQ1=matchPFJetToParton(q1Index, ijet);
      matchedToQ2=matchPFJetToParton(q2Index, ijet);
      matchedToAss=matchPFJetToParton(assPartonIndex, ijet);

      double matchedProb_zqq =(matchedToQ1 || matchedToQ2)?
	1.0: 0.0;

      if(jet_index==0 && matchedProb_zqq > 0.5)nPass[7]++;

      double matchedProb_ass = matchedToAss?
	1.0: 0.0;

      if(jet_index==0 && matchedProb_ass > 0.5)nPass[8]++;
      
      double matchedProb_either = (matchedToQ1 || matchedToQ2 ||
				   matchedToAss)?
	1.0: 0.0;

      pf_dR_matchProb_h[jet_index]->
	Fill(dR_parton,matchedProb_zqq, rec_weight);

      pf_dR_matchProb_h_ass[jet_index]->
	Fill(dR_parton,matchedProb_ass, rec_weight);
	
      pf_dR_matchProb_h_either[jet_index]->
	Fill(dR_parton,matchedProb_either, rec_weight);	


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
    } // end of loop over PF jets that are from the best higgs candidate
    
    double mH_rec  = higgsM->at(myBest);
    double mll_rec = zllM->at(myBest);
    double mjj_rec = zjjM->at(myBest);

    // if both jets are matched
    if(pfJetIndex[0]>=0 && pfJetIndex[1]>=0){

      nPass[4]++;
      if(DEBUG==1)cout << 
	"pfJetIndex1 = " << pfJetIndex[0] << "\t" <<
	"pfJetIndex2 = " << pfJetIndex[1] << endl;	  

      h_mh_rec->Fill(mH_rec, rec_weight);
      h_mll_rec->Fill(mll_rec, rec_weight);
      h_mjj_rec->Fill(mjj_rec, rec_weight);
	
      if(
	 mll_rec > MIN_MZ_LL  && mll_rec < MAX_MZ_LL  
	 && mjj_rec > LOOSE_MIN_MZ_JJ && mjj_rec < LOOSE_MAX_MZ_JJ  
	 // 	 && mZll_particle > MIN_MZ_LL && mZll_particle < MAX_MZ_LL  
	 // 	 &&  mZjj_particle > LOOSE_MIN_MZ_JJ && mZjj_particle < LOOSE_MAX_MZ_JJ
	 // 	 && jet1Index>=0 && jet2Index>=0 && jet1Index!=jet2Index
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
	    
    } // if both jets are matched
    


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
    pf_dR_matchProb_gen_ass[ip]->Write();
    pf_dR_matchProb_gen_either[ip]->Write();

    pf_dR_matchProb_gen_afterMjj[ip]->Write();
    pf_dR_matchProb_gen_afterMjj_ass[ip]->Write();
    pf_dR_matchProb_gen_afterMjj_either[ip]->Write();

    pf_dR_matchProb_h[ip]->Write();
    pf_dR_matchProb_h_ass[ip]->Write();
    pf_dR_matchProb_h_either[ip]->Write();

    pf_dR_matchProb_rec[ip]->Write();
    pf_dR_matchProb_rec_ass[ip]->Write();
    pf_dR_matchProb_rec_either[ip]->Write();

    pf_dR_matchProb_rec_afterMjj[ip]->Write();
    pf_dR_matchProb_rec_afterMjj_ass[ip]->Write();
    pf_dR_matchProb_rec_afterMjj_either[ip]->Write();




  }


  for(int ib=0; ib<NBINS; ib++){      

    h_rm_gen[ib]->Write();
    h_rm_rec[ib]->Write();

    for(int ip=0; ip<2; ip++){
      h_rpt_gen[ib][ip]->Write();
      h_rpt_rec[ib][ip]->Write();
      h_qpt_parton[ib][ip]->Write();
      h_qeta_parton[ib][ip]->Write();
      h_match[ib][ip]->Write();
    }
  }


  h_mh_stable->Write();
  h_mll_stable->Write();
  h_mjj_stable->Write();

  h_mh_rec->Write();
  h_mll_rec->Write();
  h_mjj_rec->Write();


  outFile->Close();  

}


Bool_t onlymjj::matchGenJetToParton(Int_t igen, Int_t ijet){

  Bool_t matched = false;
  if(igen < 0 || ijet < 0)return matched;
  double dR = deltaR(genParEta_->at(igen), genParPhi_->at(igen),
		     genJetEta_->at(ijet), genJetPhi_->at(ijet)
		     ); 

  if(dR<0.5)
    matched = true;
  
  return matched;  


}

Bool_t onlymjj::matchPFJetToParton(Int_t igen, Int_t ijet){

  Bool_t matched = false;
  if(igen < 0 || ijet < 0)return matched;
  double dR = deltaR(genParEta_->at(igen), genParPhi_->at(igen),
		     jetEta->at(ijet), jetPhi->at(ijet)
		     ); 

  if(dR<0.5)
    matched = true;
  
  return matched;  


}

Bool_t onlymjj::matchPatJetToParton(Int_t igen, Int_t ijet){

  Bool_t matched = false;
  if(igen < 0 || ijet < 0)return matched;
  double dR = deltaR(genParEta_->at(igen), genParPhi_->at(igen),
		     patJetPfAk05Eta_->at(ijet), patJetPfAk05Phi_->at(ijet)
		     ); 

  if(dR<0.5)
    matched = true;
  
  return matched;  


}
