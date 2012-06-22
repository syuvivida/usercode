#define genHZZ2l2q_cxx
#include "genHZZ2l2q.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
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

void genHZZ2l2q::Loop(int DEBUG)
{
  if (fChain == 0) return;

  TH1D* h_mH_template = new TH1D("h_mH_template","",280,100,1500);

  TH1D* h_mH_parton   = (TH1D*)h_mH_template->Clone("h_mH_parton");
  h_mH_parton->SetXTitle("Generated M_{llqq} [GeV/c^{2}]");

  TH1D* h_mH_particle = (TH1D*)h_mH_template->Clone("h_mH_particle");
  h_mH_particle->SetXTitle("Generated M_{lljj} [GeV/c^{2}]");

  TH1D* h_mZll_template = new TH1D("h_mZll_template","",80,60,140);
  TH1D* h_mZll_parton   = (TH1D*)h_mZll_template->Clone("h_mZll_parton");
  h_mZll_parton->SetXTitle("Generated status=3 M_{ll} [GeV/c^{2}]");

  TH1D* h_mZll_particle = (TH1D*)h_mZll_template->Clone("h_mZll_particle");
  h_mZll_particle->SetXTitle("Generated status=1 M_{ll} [GeV/c^{2}]");

  TH1D* h_mZjj_template = new TH1D("h_mZjj_template","",40,60,140);

  TH1D* h_mZjj_parton   = (TH1D*)h_mZjj_template->Clone("h_mZjj_parton");
  h_mZjj_parton->SetXTitle("Generated M_{qq} [GeV/c^{2}]");

  TH1D* h_mZjj_particle = (TH1D*)h_mZjj_template->Clone("h_mZjj_particle");
  h_mZjj_particle->SetXTitle("Generated M_{jj} [GeV/c^{2}]");
   

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    //     if(jentry >10) break;

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
	double pt  = genParPt_->at(igen);
	double eta  = genParEta_->at(igen);


	//  	if(pt < 25)continue;
	//  	if(fabs(eta)>2.4)continue;
	
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

    h_mH_parton->Fill(mH_parton);
    h_mZll_parton->Fill(mZll_parton);
    h_mZjj_parton->Fill(mZjj_parton);


    if(jet1Index < 0|| jet2Index < 0)continue;

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

    h_mH_particle->Fill(mH_particle);
    h_mZll_particle->Fill(mZll_particle);
    h_mZjj_particle->Fill(mZjj_particle);



  }
   

  std::string remword  ="../runJob/cutROOT/";

  size_t pos  = _inputFile.find(remword);

  if(pos!= std::string::npos)
    _inputFile.swap(_inputFile.erase(pos,remword.length()));


  TFile* outFile = new TFile(Form("histo_%s",
				  _inputFile.data()),"recreate");       

  h_mH_parton   ->Write();
  h_mH_particle ->Write();

  h_mZll_parton   ->Write();
  h_mZll_particle ->Write();

  h_mZjj_parton   ->Write();
  h_mZjj_particle ->Write();
   


  outFile->Close();

   
}


Bool_t genHZZ2l2q::matchGenToParton(Int_t igen, Int_t ijet){

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
