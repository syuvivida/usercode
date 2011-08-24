#define cosTheta_cxx
#include "cosTheta.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include "myLib.h"

const Float_t BARREL_MAXETA=1.4442;
const Float_t ENDCAP_MINETA=1.566;
const Float_t ENDCAP_MAXETA=2.5;

const Int_t   NETA=2;
const Float_t fEtaBin[]={0,ENDCAP_MINETA,ENDCAP_MAXETA};

// 0: no filter, 1: quark-quark, 2: quark-gluon, 3: gluon-gluon
void cosTheta::Loop(Int_t mode, Float_t pstarmin, Float_t pstarmax,
		    Float_t ybmin, Float_t ybmax) 
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t n_qq=0; // quark-quark
   Long64_t n_qg=0; // quark-gluon
   Long64_t n_gg=0; // gluon-gluon

   if(mode==0)cout << "All processes are accepted" << endl;
   else if(mode==1)cout << "Filtering quark-quark interaction" << endl;
   else if(mode==2)cout << "Filtering quark-gluon interaction" << endl;
   else if(mode==3)cout << "Filtering gluon-gluon interaction" << endl;

   char* MODE[4]={"all", "qq", "qg","gg"};
   // declare and define histograms
   TH1F* h_dec = new TH1F("h_dec","",2,fEtaBin);


   // pt distribution
   TH1F* h_pt_template = new TH1F("h_pt_template","",500,0,500);
   TH1F* h_pt_dirgamma = (TH1F*) h_pt_template->Clone("h_pt_dirgamma");
   TH1F* h_pt_fraggamma = (TH1F*) h_pt_template->Clone("h_pt_fraggamma");
   TH1F* h_pt_1stjet = (TH1F*) h_pt_template->Clone("h_pt_1stjet");
   TH1F* h_pt_2ndjet = (TH1F*) h_pt_template->Clone("h_pt_2ndjet");
   
   // rapidity distribution
   TH1F* h_y_template = new TH1F("h_y_template","",100,-5,5);
   TH1F* h_y_dirgamma = (TH1F*) h_y_template->Clone("h_y_dirgamma");
   TH1F* h_y_fraggamma = (TH1F*) h_y_template->Clone("h_y_fraggamma");
   TH1F* h_y_1stjet = (TH1F*) h_y_template->Clone("h_y_1stjet");
   TH1F* h_y_2ndjet = (TH1F*) h_y_template->Clone("h_y_2ndjet");

   // some CM variables 
   TH1F* h_pstar_template = new TH1F("h_pstar_template","",500,0,500);

   TH1F* h_pstar_dirgamma1stjet = (TH1F*)h_pstar_template->Clone("h_pstar_dirgamma1stjet");
   h_pstar_dirgamma1stjet->SetXTitle("p^{*}(#gamma^{direct},"
				      "jet^{1st})");
   
   TH1F* h_pstar_dirgamma2ndjet = (TH1F*)h_pstar_template->Clone("h_pstar_dirgamma2ndjet");
   h_pstar_dirgamma2ndjet->SetXTitle("p^{*}(#gamma^{direct},"
				      "jet^{2nd})");
   
   TH1F* h_pstar_dijet = (TH1F*)h_pstar_template->Clone("h_pstar_dijet");
   h_pstar_dijet->SetXTitle("p^{*}(jet^{1st},"
					"jet^{2nd})");

   TH1F* h_pstar_fraggamma1stjet = (TH1F*)h_pstar_template->Clone("h_pstar_fraggamma1stjet");
   h_pstar_fraggamma1stjet->SetXTitle("p^{*}(#gamma^{frag},"
					"jet^{1st})");

   TH1F* h_pstar_fraggamma2ndjet = (TH1F*)h_pstar_template->Clone("h_pstar_fraggamma2ndjet");
   h_pstar_fraggamma2ndjet->SetXTitle("p^{*}(#gamma^{frag},"
				       "jet^{2nd})");

   TH1F* h_pstar_debug = (TH1F*)h_pstar_template->Clone("h_pstar_debug");
   h_pstar_debug->SetXTitle("p^{*}(#gamma^{direct},"
				      "jet^{1st})");

   ////////////////////////////////////////////////////////////////////////


   TH1F* h_yB_template = new TH1F("h_yB_template","",100,-5.0,5.0);

   TH1F* h_yB_dirgamma1stjet = (TH1F*)h_yB_template->Clone("h_yB_dirgamma1stjet");
   h_yB_dirgamma1stjet->SetXTitle("y_{B}(#gamma^{direct},"
				      "jet^{1st})");
   
   TH1F* h_yB_dirgamma2ndjet = (TH1F*)h_yB_template->Clone("h_yB_dirgamma2ndjet");
   h_yB_dirgamma2ndjet->SetXTitle("y_{B}(#gamma^{direct},"
				      "jet^{2nd})");
   
   TH1F* h_yB_dijet = (TH1F*)h_yB_template->Clone("h_yB_dijet");
   h_yB_dijet->SetXTitle("y_{B}(jet^{1st},"
					"jet^{2nd})");

   TH1F* h_yB_fraggamma1stjet = (TH1F*)h_yB_template->Clone("h_yB_fraggamma1stjet");
   h_yB_fraggamma1stjet->SetXTitle("y_{B}(#gamma^{frag},"
					"jet^{1st})");

   TH1F* h_yB_fraggamma2ndjet = (TH1F*)h_yB_template->Clone("h_yB_fraggamma2ndjet");
   h_yB_fraggamma2ndjet->SetXTitle("y_{B}(#gamma^{frag},"
				       "jet^{2nd})");

   TH1F* h_yB_debug = (TH1F*)h_yB_template->Clone("h_yB_debug");
   h_yB_debug->SetXTitle("y_{B}(#gamma^{direct},"
				      "jet^{1st})");


   // z_gammma proposed by JETPHOX
   TH1F* h_zgamma_template = new TH1F("h_zgamma_template","",
				    200,-4.0,4.0);

   TH1F* h_zgamma_dirgamma1stjet = (TH1F*)h_zgamma_template->Clone("h_zgamma_dirgamma1stjet");
   h_zgamma_dirgamma1stjet->SetXTitle("z_{#gamma}(#gamma^{direct},"
				      "jet^{1st})");
   
   TH1F* h_zgamma_dirgamma2ndjet = (TH1F*)h_zgamma_template->Clone("h_zgamma_dirgamma2ndjet");
   h_zgamma_dirgamma2ndjet->SetXTitle("z_{#gamma}(#gamma^{direct},"
				      "jet^{2nd})");
   
   TH1F* h_zgamma_dijet = (TH1F*)h_zgamma_template->Clone("h_zgamma_dijet");
   h_zgamma_dijet->SetXTitle("z_{#gamma}(jet^{1st},"
					"jet^{2nd})");

   TH1F* h_zgamma_fraggamma1stjet = (TH1F*)h_zgamma_template->Clone("h_zgamma_fraggamma1stjet");
   h_zgamma_fraggamma1stjet->SetXTitle("z_{#gamma}(#gamma^{frag},"
					"jet^{1st})");

   TH1F* h_zgamma_fraggamma2ndjet = (TH1F*)h_zgamma_template->Clone("h_zgamma_fraggamma2ndjet");
   h_zgamma_fraggamma2ndjet->SetXTitle("z_{#gamma}(#gamma^{frag},"
				       "jet^{2nd})");



   // delta phi
   TH1F* h_dPhi_template = new TH1F("h_dPhi_template","",
				    100,0.0,TMath::Pi());
   TH1F* h_dPhi_dirgamma1stjet = (TH1F*)h_dPhi_template->Clone("h_dPhi_dirgamma1stjet");
   h_dPhi_dirgamma1stjet->SetXTitle("#Delta#phi(#gamma^{direct},"
					"jet^{1st})");
   
   TH1F* h_dPhi_dirgamma2ndjet = (TH1F*)h_dPhi_template->Clone("h_dPhi_dirgamma2ndjet");
   h_dPhi_dirgamma2ndjet->SetXTitle("#Delta#phi(#gamma^{direct},"
					"jet^{2nd})");

   TH1F* h_dPhi_dijet = (TH1F*)h_dPhi_template->Clone("h_dPhi_dijet");
   h_dPhi_dijet->SetXTitle("#Delta#phi(jet^{1st},"
					"jet^{2nd})");

   TH1F* h_dPhi_fraggamma1stjet = (TH1F*)h_dPhi_template->Clone("h_dPhi_fraggamma1stjet");
   h_dPhi_fraggamma1stjet->SetXTitle("#Delta#phi(#gamma^{frag},"
					"jet^{1st})");

   TH1F* h_dPhi_fraggamma2ndjet = (TH1F*)h_dPhi_template->Clone("h_dPhi_fraggamma2ndjet");
   h_dPhi_fraggamma2ndjet->SetXTitle("#Delta#phi(#gamma^{frag},"
					"jet^{2nd})");

   
   // costheta
   TH1F* h_cosTheta_template = new TH1F("h_cosTheta_template","",
					100,0.0,1.0);
      
   TH1F* h_cosTheta_dirgamma1stjet = (TH1F*)h_cosTheta_template->Clone("h_cosTheta_dirgamma1stjet");
   h_cosTheta_dirgamma1stjet->SetXTitle("cos#theta^{*}(#gamma^{direct},"
					"jet^{1st})");
   
   TH1F* h_cosTheta_dirgamma2ndjet = (TH1F*)h_cosTheta_template->Clone("h_cosTheta_dirgamma2ndjet");
   h_cosTheta_dirgamma2ndjet->SetXTitle("cos#theta^{*}(#gamma^{direct},"
					"jet^{2nd})");

   TH1F* h_cosTheta_dijet = (TH1F*)h_cosTheta_template->Clone("h_cosTheta_dijet");
   h_cosTheta_dijet->SetXTitle("cos#theta^{*}(jet^{1st},"
					"jet^{2nd})");

   TH1F* h_cosTheta_fraggamma1stjet = (TH1F*)h_cosTheta_template->Clone("h_cosTheta_fraggamma1stjet");
   h_cosTheta_fraggamma1stjet->SetXTitle("cos#theta^{*}(#gamma^{frag},"
					"jet^{1st})");

   TH1F* h_cosTheta_fraggamma2ndjet = (TH1F*)h_cosTheta_template->Clone("h_cosTheta_fraggamma2ndjet");
   h_cosTheta_fraggamma2ndjet->SetXTitle("cos#theta^{*}(#gamma^{frag},"
					"jet^{2nd})");


   // chi
   TH1F* h_chi_template = new TH1F("h_chi_template","",
					100,0.0,20.0);      
      
   TH1F* h_chi_dirgamma1stjet = (TH1F*)h_chi_template->Clone("h_chi_dirgamma1stjet");
   h_chi_dirgamma1stjet->SetXTitle("#chi(#gamma^{direct},"
					"jet^{1st})");
   

   TH1F* h_chi_dirgamma2ndjet = (TH1F*)h_chi_template->Clone("h_chi_dirgamma2ndjet");
   h_chi_dirgamma2ndjet->SetXTitle("#chi(#gamma^{direct},"
					"jet^{2nd})");
   
   TH1F* h_chi_dijet = (TH1F*)h_chi_template->Clone("h_chi_dijet");
   h_chi_dijet->SetXTitle("#chi(jet^{1st},"
			  "jet^{2nd})");

   TH1F* h_chi_fraggamma1stjet = (TH1F*)h_chi_template->Clone("h_chi_fraggamma1stjet");
   h_chi_fraggamma1stjet->SetXTitle("#chi(#gamma^{frag},"
				    "jet^{1st})");

   TH1F* h_chi_fraggamma2ndjet = (TH1F*)h_chi_template->Clone("h_chi_fraggamma2ndjet");
   h_chi_fraggamma2ndjet->SetXTitle("#chi(#gamma^{frag},"
				    "jet^{2nd})");



   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      // check if which physics process this is

      Int_t eventMode = 0;

      if(mcPID[4]==21 && mcPID[5]==21){ n_gg +=1; eventMode=3;}

      else if(abs(mcPID[4])<=6 && abs(mcPID[5])<=6){ n_qq +=1; eventMode=1;}

      else if((abs(mcPID[4])<=6 && mcPID[5]==21) ||
	      (mcPID[4]==21 && abs(mcPID[5])<=6)){ n_qg +=1; eventMode=2;}

      // apply filter
      if(mode!=0 && eventMode!=mode)continue;


      Int_t ngood_vtx=nGoodVtx(ientry);
//       if(ngood_vtx<1)continue;
      
      // find the leading photon at the generator-level
      bool findADirPhoton=false;
      bool findAFraPhoton=false;

      TLorentzVector gen_direct_photon(0,0,0,0);
      TLorentzVector gen_fragmentation_photon(0,0,0,0);


      for(int imc=0; imc < nMC; imc++){
	
	// check if this is a real photon

	if(mcPID[imc]!=22)continue;
  	if(mcPt[imc]< 20.0)continue;
   	if(fabs(mcEta[imc]) > 2.5)continue;


	if(mcMomPID[imc]==22 && imc > 7 && !findADirPhoton)
	  {
	    gen_direct_photon.SetPtEtaPhiE(
					   mcPt[imc],
					   mcEta[imc],
					   mcPhi[imc],
					   mcE[imc]
					   );
	    findADirPhoton = true;
	  }

	if(fabs(mcMomPID[imc]) < 22 && imc > 7 && !findAFraPhoton)
	  {
	    gen_fragmentation_photon.SetPtEtaPhiE(
						  mcPt[imc],
						  mcEta[imc],
						  mcPhi[imc],
						  mcE[imc]
						  );
	    findAFraPhoton = true;
	  }


      } // loop over MC particles



      if(findADirPhoton)
	{
	  h_y_dirgamma->Fill(gen_direct_photon.Rapidity());
	  h_pt_dirgamma->Fill(gen_direct_photon.Pt());
	}

      if(findAFraPhoton)
	{
	  h_y_fraggamma->Fill(gen_fragmentation_photon.Rapidity());
	  h_pt_fraggamma->Fill(gen_fragmentation_photon.Pt());
	}

      // first check which reco jet is the one from the highest and 
      // second et gen jet
      bool findLeadingJet = false;
      Int_t leadingJetIndex=-1;
      Float_t genJetMaxPt=-9999.;

      bool findSecondLeadingJet = false;
      Int_t secondLeadingJetIndex=-1;
      Float_t secondJetMaxPt=-9999.;

      for(int ijet=0; ijet < nJet; ijet++){

	if(jetGenJetIndex[ijet]<1)continue;
	if(jetGenJetPt[ijet] > genJetMaxPt)
	  {
	    secondJetMaxPt = genJetMaxPt;
	    secondLeadingJetIndex = leadingJetIndex;
	    genJetMaxPt = jetGenJetPt[ijet];
	    leadingJetIndex= ijet;
	  }
	else if(jetGenJetPt[ijet] > secondJetMaxPt)
	  {
	    secondJetMaxPt = jetGenJetPt[ijet];
	    secondLeadingJetIndex = ijet;
	  }


      } // end of loop over jets


      // only study
      if(leadingJetIndex >=0 &&
	 (
 	  jetGenJetPt[leadingJetIndex] > 20.0 &&
	  fabs(jetGenJetEta[leadingJetIndex]) < 3.0
	  )
	 )
	findLeadingJet = true;
	 
      if(secondLeadingJetIndex >=0 &&
	 (
	  jetGenJetPt[secondLeadingJetIndex] > 20.0 &&
	  fabs(jetGenJetEta[secondLeadingJetIndex]) < 3.0
	  )
	 )
	findSecondLeadingJet = true;


      TLorentzVector gen_1stjet(0,0,0,0);
      if(findLeadingJet)
	{
	  gen_1stjet.SetPtEtaPhiE(
				  jetGenJetPt[leadingJetIndex],
				  jetGenJetEta[leadingJetIndex],
				  jetGenJetPhi[leadingJetIndex],
				  jetGenJetEn[leadingJetIndex]
				  );
	  h_y_1stjet->Fill(gen_1stjet.Rapidity());
 	  h_pt_1stjet->Fill(gen_1stjet.Pt());
	  
	}
      

      

      TLorentzVector gen_2ndjet(0,0,0,0);           
      if(findSecondLeadingJet)
	{
	  gen_2ndjet.SetPtEtaPhiE(
				  jetGenJetPt[secondLeadingJetIndex],
				  jetGenJetEta[secondLeadingJetIndex],
				  jetGenJetPhi[secondLeadingJetIndex],
				  jetGenJetEn[secondLeadingJetIndex]
				  );

	   h_y_2ndjet->Fill(gen_2ndjet.Rapidity());
	   h_pt_2ndjet->Fill(gen_2ndjet.Pt());
	}

      
      // now calculating cosTheta and chi
      if(findADirPhoton && findLeadingJet 
	 && separated(gen_direct_photon, gen_1stjet)
	 )
	{

	  Double_t pstar_temp  = pstar(gen_direct_photon,gen_1stjet);
	  Double_t yB_temp     = yB(gen_direct_photon,gen_1stjet);
	  h_pstar_dirgamma1stjet->Fill(pstar_temp);
	  h_yB_dirgamma1stjet->Fill(yB_temp);
	  h_zgamma_dirgamma1stjet->Fill(zgamma(gen_direct_photon,gen_1stjet));
	  h_dPhi_dirgamma1stjet->Fill(deltaPhi(gen_direct_photon,gen_1stjet));

 	  if(pstar_temp > pstarmin && pstar_temp < pstarmax 
	     && fabs(yB_temp) > ybmin && fabs(yB_temp) < ybmax)
	    {
	      h_cosTheta_dirgamma1stjet->Fill(cosThetaStar(gen_direct_photon, gen_1stjet));
	      h_pstar_debug->Fill(pstar_temp);
	      h_yB_debug   ->Fill(yB_temp);
	    }

	    h_chi_dirgamma1stjet->Fill(chiPair(gen_direct_photon, gen_1stjet));
	}
      
      if(findADirPhoton && findSecondLeadingJet
	 && separated(gen_direct_photon, gen_2ndjet)
	 )
	{
	  h_pstar_dirgamma2ndjet->Fill(pstar(gen_direct_photon,gen_2ndjet));
	  h_yB_dirgamma2ndjet->Fill(yB(gen_direct_photon,gen_2ndjet));
	  h_zgamma_dirgamma2ndjet->Fill(zgamma(gen_direct_photon,gen_2ndjet));
	  h_dPhi_dirgamma2ndjet->Fill(deltaPhi(gen_direct_photon,gen_2ndjet));
	  h_cosTheta_dirgamma2ndjet->Fill(cosThetaStar(gen_direct_photon, gen_2ndjet));	  
	  h_chi_dirgamma2ndjet->Fill(chiPair(gen_direct_photon,gen_2ndjet));

	}

      if(findLeadingJet && findSecondLeadingJet 
	 && separated(gen_1stjet, gen_2ndjet)
	 )
	{
	  h_pstar_dijet->Fill(pstar(gen_1stjet, gen_2ndjet));
	  h_yB_dijet->Fill(yB(gen_1stjet, gen_2ndjet));
	  h_zgamma_dijet->Fill(zgamma(gen_1stjet, gen_2ndjet));
	  h_dPhi_dijet->Fill(deltaPhi(gen_1stjet, gen_2ndjet));
	  h_cosTheta_dijet->Fill(cosThetaStar(gen_1stjet, gen_2ndjet));
	  h_chi_dijet->Fill(chiPair(gen_1stjet, gen_2ndjet));
	}


      if(findAFraPhoton && findLeadingJet 
	 && separated(gen_fragmentation_photon,gen_1stjet)
	 )
	{
	  h_pstar_fraggamma1stjet->Fill(pstar(gen_fragmentation_photon,gen_1stjet));
	  h_yB_fraggamma1stjet->Fill(yB(gen_fragmentation_photon,gen_1stjet));
	  h_zgamma_fraggamma1stjet->Fill(zgamma(gen_fragmentation_photon,gen_1stjet));
	  h_dPhi_fraggamma1stjet->Fill(deltaPhi(gen_fragmentation_photon,gen_1stjet));
	  h_cosTheta_fraggamma1stjet->Fill(cosThetaStar(gen_fragmentation_photon,gen_1stjet));
	  h_chi_fraggamma1stjet->Fill(chiPair(gen_fragmentation_photon,gen_1stjet));
	}


      if(findAFraPhoton && findSecondLeadingJet
	 && separated(gen_fragmentation_photon,gen_2ndjet)
	 )
	{

	  h_pstar_fraggamma2ndjet->Fill(pstar(gen_fragmentation_photon,gen_2ndjet));
	  h_yB_fraggamma2ndjet->Fill(yB(gen_fragmentation_photon,gen_2ndjet));
	  h_zgamma_fraggamma2ndjet->Fill(zgamma(gen_fragmentation_photon,gen_2ndjet));
	  h_dPhi_fraggamma2ndjet->Fill(deltaPhi(gen_fragmentation_photon,gen_2ndjet));
	  h_cosTheta_fraggamma2ndjet->Fill(cosThetaStar(gen_fragmentation_photon,gen_2ndjet));
	  h_chi_fraggamma2ndjet->Fill(chiPair(gen_fragmentation_photon,gen_2ndjet));
	}

  
   } // end of loop over entries


   cout << "=================== Summary ================= " << endl;
   cout << "Number of quark-quark interaction is " << n_qq << endl;
   cout << "Number of quark-gluon interaction is " << n_qg << endl;
   cout << "Number of gluon-gluon interaction is " << n_gg << endl;
   cout << "============================================= " << endl;
   
   std::string remword="/data2/syu/ggNTuple/";
   size_t pos = inputFile_.find(remword);
  
   if(pos!= std::string::npos)
     inputFile_.swap(inputFile_.erase(pos,remword.length()));
   
   TFile* outFile = new TFile(Form("jetphoxPaper_cosTheta_pstar%dto%d_"
				   "yb%.1lf""to""%.1lf_%s_%s",
				   (Int_t)pstarmin, (Int_t)pstarmax,
				   ybmin, ybmax,
				   MODE[mode],inputFile_.data()),"recreate");               
   h_y_dirgamma -> Write();
   h_y_fraggamma -> Write();
   h_y_1stjet -> Write();
   h_y_2ndjet -> Write();

   h_pt_dirgamma -> Write();
   h_pt_fraggamma -> Write();
   h_pt_1stjet -> Write();
   h_pt_2ndjet -> Write();

   h_pstar_dirgamma1stjet -> Write();
   h_pstar_dirgamma2ndjet -> Write();
   h_pstar_dijet-> Write();
   h_pstar_fraggamma1stjet-> Write();
   h_pstar_fraggamma2ndjet-> Write();
   h_pstar_debug->Write();

   h_yB_dirgamma1stjet -> Write();
   h_yB_dirgamma2ndjet -> Write();
   h_yB_dijet-> Write();
   h_yB_fraggamma1stjet-> Write();
   h_yB_fraggamma2ndjet-> Write();
   h_yB_debug->Write();

   h_zgamma_dirgamma1stjet -> Write();
   h_zgamma_dirgamma2ndjet -> Write();
   h_zgamma_dijet-> Write();
   h_zgamma_fraggamma1stjet-> Write();
   h_zgamma_fraggamma2ndjet-> Write();

   h_dPhi_dirgamma1stjet -> Write();
   h_dPhi_dirgamma2ndjet -> Write();
   h_dPhi_dijet-> Write();
   h_dPhi_fraggamma1stjet-> Write();
   h_dPhi_fraggamma2ndjet-> Write();

   h_cosTheta_dirgamma1stjet -> Write();
   h_cosTheta_dirgamma2ndjet -> Write();
   h_cosTheta_dijet-> Write();
   h_cosTheta_fraggamma1stjet-> Write();
   h_cosTheta_fraggamma2ndjet-> Write();

   h_chi_dirgamma1stjet-> Write();
   h_chi_dirgamma2ndjet-> Write();
   h_chi_dijet-> Write();
   h_chi_fraggamma1stjet-> Write();
   h_chi_fraggamma2ndjet-> Write();


   outFile->Close();      
}


Bool_t cosTheta::isFidJet (Long64_t entry, Int_t ijet)
{
  if(jetPt[ijet] < 10.0)return false;
  if(fabs(jetEta[ijet]) > 3.0)return false;
  return true;

}


// check if this reco-jet is a good loose jet
Bool_t cosTheta::isGoodLooseJet(Long64_t entry, Int_t ijet)
{
  if(jetPt[ijet] < 10.0)return false;
  if(fabs(jetEta[ijet]) > 3.0)return false;
  if(jetNConstituents[ijet] <= 1)return false;
  if(jetNHF[ijet] >= 0.99)return false;
  if(jetNEF[ijet] >= 0.99)return false;

  // for the tracker region
  if(fabs(jetEta[ijet])<2.4 && jetCHF[ijet] <= 0.0)return false;
  if(fabs(jetEta[ijet])<2.4 && jetCEF[ijet] >= 0.99)return false;
  if(fabs(jetEta[ijet])<2.4 && jetNCH[ijet] <= 0)return false;

  return true;

}


// check if this reco-jet is a good medium jet
Bool_t cosTheta::isGoodMediumJet(Long64_t entry, Int_t ijet)
{
  if(jetPt[ijet] < 10.0)return false;
  if(fabs(jetEta[ijet]) > 3.0)return false;
  if(jetNConstituents[ijet] <= 1)return false;
  if(jetNHF[ijet] >= 0.95)return false;
  if(jetNEF[ijet] >= 0.95)return false;

  // for the tracker region
  if(fabs(jetEta[ijet])<2.4 && jetCHF[ijet] <= 0.0)return false;
  if(fabs(jetEta[ijet])<2.4 && jetCEF[ijet] >= 0.99)return false;
  if(fabs(jetEta[ijet])<2.4 && jetNCH[ijet] <= 0)return false;

  return true;

}

// check if this reco-jet is a good medium jet
Bool_t cosTheta::isGoodTightJet(Long64_t entry, Int_t ijet)
{
  if(jetPt[ijet] < 10.0)return false;
  if(fabs(jetEta[ijet]) > 3.0)return false;
  if(jetNConstituents[ijet] <= 1)return false;
  if(jetNHF[ijet] >= 0.90)return false;
  if(jetNEF[ijet] >= 0.90)return false;

  // for the tracker region
  if(fabs(jetEta[ijet])<2.4 && jetCHF[ijet] <= 0.0)return false;
  if(fabs(jetEta[ijet])<2.4 && jetCEF[ijet] >= 0.99)return false;
  if(fabs(jetEta[ijet])<2.4 && jetNCH[ijet] <= 0)return false;

  return true;

}

Bool_t cosTheta::isGoodPho(Long64_t entry, Int_t ipho)
{


  Bool_t isEB=false;
  Bool_t isEE=false;
  if(fabs(phoSCEta[ipho]) < BARREL_MAXETA)isEB=true;
  if(fabs(phoSCEta[ipho]) > ENDCAP_MINETA && 
     fabs(phoSCEta[ipho]) < ENDCAP_MAXETA) isEE=true;

  if(!isEB && !isEE)return false;
  if(phoEt[ipho] < 15.0)return false;
  if(phoHoverE[ipho] > 0.05)return false;
  if(phohasPixelSeed[ipho]==1)return false;
//   if(phoEcalIsoDR04[ipho] > 4.2 +0.006 * phoEt[ipho])return false;
//   if(phoHcalIsoDR04[ipho] > 2.2 +0.0025* phoEt[ipho])return false;
//   if(phoTrkIsoHollowDR04[ipho] > 2.0 +0.001* phoEt[ipho])return false;
  if(phoEcalIsoDR04[ipho] > 4.2 +0.003 * phoEt[ipho])return false;
  if(phoHcalIsoDR04[ipho] > 2.2 +0.001* phoEt[ipho])return false;
  if(phoTrkIsoHollowDR04[ipho] > 2.0 +0.001* phoEt[ipho])return false;
  
  Float_t sumIso = phoEcalIsoDR04[ipho]+phoHcalIsoDR04[ipho]+
    phoTrkIsoHollowDR04[ipho];

//   sumIso -= isEB? (rhoCorr[0]*rho25): (rhoCorr[1]*rho25);

//   if(sumIso > 5.0)return false;
  if(isEB && phoSigmaIEtaIEta[ipho] > 0.013)return false;
  if(isEE && phoSigmaIEtaIEta[ipho] > 0.030)return false;

  return true;

}

Int_t cosTheta::nGoodVtx (Long64_t entry)
{
  Int_t ngood=0;
  for(int iv=0; iv < nVtx; iv++)
    {
      if(vtxNDF[iv] > 4 && fabs(vtxD0[iv]) < 2.0 && fabs(vtx[iv][2])<24.0)
	ngood++;
    }
  
  return ngood;


}
