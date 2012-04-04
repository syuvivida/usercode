#define gen_distribution_cxx
#include "gen_distribution.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>


void gen_distribution::Loop(int lepID, int DEBUG)
{
   if (fChain == 0) return;
   TH1D* h_y_template = new TH1D("h_y_template","",25,-2.5,2.5);

   TH1D* h_yB    = (TH1D*)h_y_template->Clone("h_yB");
   h_yB->SetXTitle("0.5(y_{Z}+y_{jet})");

   TH1D* h_ystar = (TH1D*)h_y_template->Clone("h_ystar");
   h_ystar->SetXTitle("0.5(y_{Z}-y_{jet})");



   Long64_t nentries = fChain->GetEntriesFast();

   cout << _inputFileName << " has " << nentries << " entries." << endl;
   std::string leptonName;
   if(abs(lepID)==13)leptonName = "muon";
   else if(abs(lepID)==11)leptonName = "electron";

   cout << "studying " << leptonName << endl;

   // dummy proof, in case someone put a negative number
   int leptonPID = abs(lepID); 
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
//       if(jentry > 50000)break;

      int lepPlusIndex = -1;
      int lepMinusIndex = -1;


      for(unsigned int igen=0; igen < genParId_->size(); igen++){

	if(lepPlusIndex < 0 && fabs(genParId_->at(igen)-(-leptonPID))<1e-6 && 
	   fabs(genParSt_->at(igen)-1) < 1e-6)
	  lepPlusIndex = igen;

	if(lepMinusIndex < 0 && fabs(genParId_->at(igen)-(leptonPID))<1e-6 && 
	   fabs(genParSt_->at(igen)-1) < 1e-6)
	  lepMinusIndex = igen;


	// need to add a line and see if this is a daughter of Z-> tau tau
	if(lepPlusIndex >= 0 && lepMinusIndex >=0)
	  break;
	
      }

      // do not find m+ or lep-
      if(lepPlusIndex < 0 || lepMinusIndex < 0)continue;

      if(DEBUG==1){
	cout << "lepPlusIndex = " << lepPlusIndex << "\t pt=" << 
	  genParPt_->at(lepPlusIndex) << "\t eta=" << 
	  genParEta_->at(lepPlusIndex) << endl;
	cout << "lepMinusIndex = " << lepMinusIndex << "\t pt=" << 
	  genParPt_->at(lepMinusIndex) << "\t eta=" << 
	  genParEta_->at(lepMinusIndex) << endl;
      }

      int nPt20=0;
      int nPt10=0;

      int indexNumbers[2] = {lepPlusIndex, lepMinusIndex};

      for(int ip=0; ip < 2; ip++){

	double pt  = genParPt_->at(indexNumbers[ip]);
	double eta = genParEta_->at(indexNumbers[ip]);

	if(leptonPID==13 && pt > 20.0 && fabs(eta) < 2.1)nPt20++;
	if(leptonPID==13 && pt > 10.0 && fabs(eta) < 2.1)nPt10++;

	if(leptonPID==11 && pt > 20.0 && fabs(eta) > 1.566 && 
	   fabs(eta) < 2.5)nPt20++;

	if(leptonPID==11 && pt > 20.0 && fabs(eta) > 0.0 && 
	   fabs(eta) < 1.446)nPt20++;

      }

      if(DEBUG==1)
	cout << "nPt20 = " << nPt20 << "\t nPt10 = " << nPt10 << endl;

      if(leptonPID==13 && (nPt20 < 1 || nPt10 < 2))continue;
      if(leptonPID==11 && (nPt20 < 2))continue;


      TLorentzVector l4_lepp(0,0,0,0);

      l4_lepp.SetPtEtaPhiE(genParPt_->at(lepPlusIndex),
			   genParEta_->at(lepPlusIndex),
			   genParPhi_->at(lepPlusIndex),
			   genParE_->at(lepPlusIndex));

      TLorentzVector l4_lepm(0,0,0,0);

      l4_lepm.SetPtEtaPhiE(genParPt_->at(lepMinusIndex),
			   genParEta_->at(lepMinusIndex),
			   genParPhi_->at(lepMinusIndex),
			   genParE_->at(lepMinusIndex));


      TLorentzVector l4_z = l4_lepp + l4_lepm;


      double mll = l4_z.M();

      if(leptonPID==13 && (mll < 76 || mll > 106))continue;
      if(leptonPID==11 && (mll < 70 || mll > 110))continue;

      if(DEBUG==1)
	cout << "dilepton mass = " << mll << endl;
      

      // now look for jets
      double maxGenJetPt = -9999;
      int maxGenJetIndex = -1;


      for(unsigned int ij = 0; ij < genJetPt_->size(); ij ++){

	double thisGenJetPt  = genJetPt_->at(ij);
	double thisGenJetEta = genJetEta_->at(ij);

	if(thisGenJetPt < 30.0)continue;
	if(fabs(thisGenJetEta) > 2.4)continue;

	TLorentzVector thisGenJet_l4(0,0,0,0);
	
	thisGenJet_l4.SetPtEtaPhiE(genJetPt_->at(ij),
				   genJetEta_->at(ij),
				   genJetPhi_->at(ij),
				   genJetE_->at(ij));
	
	double dr_ep = l4_lepp.DeltaR(thisGenJet_l4);
	double dr_em = l4_lepm.DeltaR(thisGenJet_l4);

	if(leptonPID==11 && dr_ep < 0.3)continue;
	if(leptonPID==11 && dr_em < 0.3)continue;


	if(thisGenJetPt > maxGenJetPt)
	  {
	    maxGenJetPt = thisGenJetPt;
	    maxGenJetIndex = ij;
	  }


      }
      
      if(DEBUG==1)
	cout << "max gen jet index = " << maxGenJetIndex << endl;

      if(maxGenJetIndex < 0)continue;

      
      TLorentzVector l4_j(0,0,0,0);
      l4_j.SetPtEtaPhiE(genJetPt_->at(maxGenJetIndex),
			genJetEta_->at(maxGenJetIndex),
			genJetPhi_->at(maxGenJetIndex),
			genJetE_->at(maxGenJetIndex));


      double yz = l4_z.Rapidity();
      double yj = l4_j.Rapidity();

      double yB = 0.5*(yz + yj);
      double ystar = 0.5*(yz-yj);

      h_yB->Fill(yB);
      h_ystar->Fill(ystar);


   } // end of loop over entries

   std::string prefix = "genHisto_";
   //   std::string remword  ="/scratch/syu/test_sherpa/genonly/";
   std::string remword  ="/scratch/syu/";

   size_t pos  = _inputFileName.find(remword);

   if(pos!= std::string::npos)
     _inputFileName.swap(_inputFileName.erase(pos,remword.length()));


   TFile* outFile = new TFile(Form("genHisto_%s_%s",leptonName.data(),
				   _inputFileName.data()),"recreate");       
        
 
   h_yB->Write();
   h_ystar->Write();

   outFile->Close();

}
