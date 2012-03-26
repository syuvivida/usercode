#define madgraph_gen_cxx
#include "madgraph_gen.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

void madgraph_gen::Loop()
{
   if (fChain == 0) return;
   TH1D* h_y_template = new TH1D("h_y_template","",25,-2.5,2.5);

   TH1D* h_yB    = (TH1D*)h_y_template->Clone("h_yB");
   h_yB->SetXTitle("0.5(y_{Z}+y_{jet})");

   TH1D* h_ystar = (TH1D*)h_y_template->Clone("h_ystar");
   h_ystar->SetXTitle("0.5(y_{Z}-y_{jet})");



   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
//       if(jentry > 50000)break;

      int muPlusIndex = -1;
      int muMinusIndex = -1;


      for(int imc=0; imc < nMC; imc++)
	{
// 	  if(muPlusIndex < 0 && mcPID[imc]==-13)
	  if(muPlusIndex < 0 && mcPID[imc]==-13 && mcStatus[imc]==1)
	    muPlusIndex = imc;

// 	  if(muMinusIndex < 0 && mcPID[imc]== 13)
	  if(muMinusIndex < 0 && mcPID[imc]== 13 && mcStatus[imc]==1)
	    muMinusIndex = imc;

	  if(muPlusIndex >= 0 && muMinusIndex >=0)
	    break;

	}
      
      // do not find m+ or mu-
      if(muPlusIndex < 0 || muMinusIndex < 0)continue;


      int nPt20=0;
      int nPt10=0;

      int indexNumbers[2] = {muPlusIndex, muMinusIndex};

//       cout << "mu + pt and eta: " << mcPt[indexNumbers[0]] << "\t" << mcEta[indexNumbers[0]] << endl;
//       cout << "mu - pt and eta: " << mcPt[indexNumbers[1]] << "\t" << mcEta[indexNumbers[1]] << endl;

      for(int ip=0; ip < 2; ip++){

	if(mcPt[indexNumbers[ip]] > 20.0 && fabs(mcEta[indexNumbers[ip]]) < 2.1)nPt20++;
	if(mcPt[indexNumbers[ip]] > 10.0 && fabs(mcEta[indexNumbers[ip]]) < 2.1)nPt10++;

      }

//       cout << "nPt20 = " << nPt20 << "\t nPt10 = " << nPt10 << endl;

      if(nPt20 < 1 || nPt10 < 2)continue;


      TLorentzVector l4_lepp(0,0,0,0);

      l4_lepp.SetPtEtaPhiE(mcPt[muPlusIndex],
			   mcEta[muPlusIndex],
			   mcPhi[muPlusIndex],
			   mcE[muPlusIndex]);

      TLorentzVector l4_lepm(0,0,0,0);

      l4_lepm.SetPtEtaPhiE(mcPt[muMinusIndex],
			   mcEta[muMinusIndex],
			   mcPhi[muMinusIndex],
			   mcE[muMinusIndex]);


      TLorentzVector l4_z = l4_lepp + l4_lepm;


      // now look for jets
      double maxGenJetPt = -9999;
      int maxGenJetIndex = -1;

      double maxGenPartonPt = -9999;
      int maxGenPartonIndex = -1;


      for(int ij = 0; ij < nJet; ij ++){
	if(jetGenJetIndex[ij] < 1)continue;

	double thisGenJetPt = jetGenJetPt[ij];
	double thisGenPartonPt = jetGenPt[ij];

	double thisGenPartonEta = jetGenEta[ij];

	if(thisGenPartonPt < 30.0)continue;
	if(fabs(thisGenPartonEta) > 2.4)continue;

	if(thisGenJetPt > maxGenJetPt)
	  {
	    maxGenJetPt = thisGenJetPt;
	    maxGenJetIndex = ij;
	  }

	if(thisGenPartonPt > maxGenPartonPt)
	  {
	    maxGenPartonPt = thisGenPartonPt;
	    maxGenPartonIndex = ij;
	  }


      }

//       cout << "max gen jet index = " << maxGenJetIndex << "\t";
//       cout << "max gen parton index = " << maxGenPartonIndex << endl;

      if(maxGenJetIndex < 0)continue;

      
      TLorentzVector l4_j(0,0,0,0);
      l4_j.SetPtEtaPhiE(jetGenPt[maxGenJetIndex],
			jetGenEta[maxGenJetIndex],
			jetGenPhi[maxGenJetIndex],
			jetGenEn[maxGenJetIndex]);


      double yz = l4_z.Rapidity();
      double yj = l4_j.Rapidity();

      double yB = 0.5*(yz + yj);
      double ystar = 0.5*(yz-yj);

      h_yB->Fill(yB);
      h_ystar->Fill(ystar);


   } // end of loop over entries

   std::string prefix = "genHisto_";
   std::string remword  ="/data4/syu/7TeV_ggNtuple/madgraph/";

   size_t pos  = _inputFileName.find(remword);

   if(pos!= std::string::npos)
     _inputFileName.swap(_inputFileName.erase(pos,remword.length()));


   TFile* outFile = new TFile(Form("genHisto_%s",_inputFileName.data()),"recreate");               
 
   h_yB->Write();
   h_ystar->Write();

   outFile->Close();

}
