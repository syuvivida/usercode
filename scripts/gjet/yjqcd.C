#define yjqcd_cxx
#include "yjqcd.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TSystem.h>
#include <myLib.h>
#include <TLorentzVector.h>

void yjqcd::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   Long64_t nPass[50]={0};

   const int NDEC=2;
   const int NREGION=5;

   TH1F* hphopt  = new TH1F("hphopt","",100,0,500);
   TH1F* hphoeta = new TH1F("hphoeta","",120,-3.0,3.0);
   TH2F* hphoetavspt = new TH2F("hphoetavspt","",100,0,500,120,-3.0,3.0);

   TH1F* hngoodpho = new TH1F("hngoodpho","",6,-0.5,5.5);
   TH1F* hngoodjet = new TH1F("hngoodjet","",30,-0.5,29.5);

   // photon ID variables
   TH1F* hsieie_template = new TH1F("hsieie_template","",
				    500,0,0.10);
   TH1F* hsieie_dec[NDEC];
   
   for(int i=0; i<NDEC; i++){

     hsieie_dec[i] = (TH1F*)hsieie_template->Clone(
						   Form("hsieie_dec%d",i));
     hsieie_dec[i]->SetXTitle("#sigma_{i#eta i#eta}");
   }



   TH1F* hdR_template = new TH1F("hdR_template","",
				 100,0,10.0);
   TH1F* hdR_dec[NDEC];
   
   for(int i=0; i<NDEC; i++){

     hdR_dec[i] = (TH1F*)hdR_template->Clone(
					     Form("hdR_dec%d",i));
     hdR_dec[i]->SetXTitle("#Delta R");
   }


   // some CM variables 
   TH1F* h_pstar_template = new TH1F("h_pstar_template","",500,0,500);

   ////////////////////////////////////////////////////////////////////////


   TH1F* h_yB_template = new TH1F("h_yB_template","",100,-5.0,5.0);

   // z_gammma proposed by JETPHOX
   TH1F* h_zgamma_template = new TH1F("h_zgamma_template","",
				    200,-4.0,4.0);

   // delta phi
   TH1F* h_dPhi_template = new TH1F("h_dPhi_template","",
				    100,0.0,TMath::Pi());
   

//    Float_t boundsieie[NREGION]={-1,0.0088,0.0092,0.01,0.0116};

//    for(int i=0; i<NDEC; i++){

//      for(int j=0; j<NREGION; j++){




       
//      } // end of loop over sieie regions

//    } // end of detector region loop


   //
   TH2F* h_dPhi_sieie_template = new TH2F("h_dPhi_sieie_template","",
					  100,0,0.05,100,0.0,TMath::Pi());
   
   TH2F* h_dPhi_sieie_dec[NDEC];
   for(int i=0; i<NDEC; i++){

     h_dPhi_sieie_dec[i] = (TH2F*)h_dPhi_sieie_template->Clone(Form("h_dPhi_sieie_dec%d",i));
     h_dPhi_sieie_dec[i] -> SetXTitle("#sigma_{i#eta i#eta}");
     h_dPhi_sieie_dec[i] -> SetYTitle("#Delta#phi(#gamma,jet)");

   }

   TH2F* h_dPhiAllJets_sieie_dec[NDEC];
   for(int i=0; i<NDEC; i++){

     h_dPhiAllJets_sieie_dec[i] = (TH2F*)h_dPhi_sieie_template->Clone(Form("h_dPhiAllJets_sieie_dec%d",i));
     h_dPhiAllJets_sieie_dec[i] -> SetXTitle("#sigma_{i#eta i#eta}");
     h_dPhiAllJets_sieie_dec[i] -> SetYTitle("#Delta#phi(#gamma,all jets)");

   }


   /////////////////////////////////////////////////////////////////////////

   // profile of delta phi vs sigmaIetaIeta

   /////////////////////////////////////////////////////////////////////////


   TProfile* pf_dPhi_sieie_template = new TProfile("pf_dPhi_sieie_template","",
					  50,0,0.10,0.0,TMath::Pi());
   
   TProfile* pf_dPhi_sieie_dec[NDEC];
   for(int i=0; i<NDEC; i++){

     pf_dPhi_sieie_dec[i] = (TProfile*)pf_dPhi_sieie_template->Clone(Form("pf_dPhi_sieie_dec%d",i));
     pf_dPhi_sieie_dec[i] -> SetXTitle("#sigma_{i#eta i#eta}");
     pf_dPhi_sieie_dec[i] -> SetYTitle("#Delta#phi(#gamma,jet)");

   }

   TProfile* pf_dPhiAllJets_sieie_dec[NDEC];
   for(int i=0; i<NDEC; i++){

     pf_dPhiAllJets_sieie_dec[i] = (TProfile*)pf_dPhi_sieie_template->Clone(Form("pf_dPhiAllJets_sieie_dec%d",i));
     pf_dPhiAllJets_sieie_dec[i] -> SetXTitle("#sigma_{i#eta i#eta}");
     pf_dPhiAllJets_sieie_dec[i] -> SetYTitle("#Delta#phi(#gamma,all jets)");

   }


   TProfile* pf_sieie_dPhi_template = new TProfile("pf_sieie_dPhi_template","",						   
						   50,0.0,TMath::Pi(),
						   0,0.10);
   
   TProfile* pf_sieie_dPhi_dec[NDEC];
   for(int i=0; i<NDEC; i++){

     pf_sieie_dPhi_dec[i] = (TProfile*)pf_sieie_dPhi_template->Clone(Form("pf_sieie_dPhi_dec%d",i));
     pf_sieie_dPhi_dec[i] -> SetXTitle("#Delta#phi(#gamma,jet)");
     pf_sieie_dPhi_dec[i] -> SetYTitle("#sigma_{i#eta i#eta}");
   }

   TProfile* pf_sieie_dPhiAllJets_dec[NDEC];
   for(int i=0; i<NDEC; i++){

     pf_sieie_dPhiAllJets_dec[i] = (TProfile*)pf_sieie_dPhi_template->Clone(Form("pf_sieie_dPhiAllJets_dec%d",i));
     pf_sieie_dPhiAllJets_dec[i] -> SetXTitle("#Delta#phi(#gamma,all jets)");
     pf_sieie_dPhiAllJets_dec[i] -> SetYTitle("#sigma_{i#eta i#eta}");

   }



   /////////////////////////////////////////////////////////////////////////

   
   // costheta
   TH1F* h_cosTheta_template = new TH1F("h_cosTheta_template","",
					100,0.0,1.0);
      

   /////////////////////////////////////////////////////////////////////////

   // profile of cosTheta vs sigmaIetaIeta

   /////////////////////////////////////////////////////////////////////////


   TProfile* pf_cosTheta_sieie_template = new TProfile("pf_cosTheta_sieie_template","",
						       50,0,0.10,0.0,1.0);


   TProfile* pf_cosTheta_sieie_dec[NDEC];
   for(int i=0; i<NDEC; i++){

     pf_cosTheta_sieie_dec[i] = (TProfile*)pf_cosTheta_sieie_template->Clone(Form("pf_cosTheta_sieie_dec%d",i));
     pf_cosTheta_sieie_dec[i] -> SetXTitle("#sigma_{i#eta i#eta}");
     pf_cosTheta_sieie_dec[i] -> SetYTitle("cos#theta^{*}(#gamma,jet)");

   }

   TProfile* pf_cosThetaAllJets_sieie_dec[NDEC];
   for(int i=0; i<NDEC; i++){

     pf_cosThetaAllJets_sieie_dec[i] = (TProfile*)pf_cosTheta_sieie_template->Clone(Form("pf_cosThetaAllJets_sieie_dec%d",i));
     pf_cosThetaAllJets_sieie_dec[i] -> SetXTitle("#sigma_{i#eta i#eta}");
     pf_cosThetaAllJets_sieie_dec[i] -> SetYTitle("cos#theta^{*}(#gamma,all jets)");

   }



   TProfile* pf_sieie_cosTheta_template = new TProfile("pf_sieie_cosTheta_template","",
						       50,0.0,1.0,0,0.10);

   TProfile* pf_sieie_cosTheta_dec[NDEC];
   for(int i=0; i<NDEC; i++){

     pf_sieie_cosTheta_dec[i] = (TProfile*)pf_sieie_cosTheta_template->Clone(Form("pf_sieie_cosTheta_dec%d",i));
     pf_sieie_cosTheta_dec[i] -> SetXTitle("cos#theta^{*}(#gamma,jet)");
     pf_sieie_cosTheta_dec[i] -> SetYTitle("#sigma_{i#eta i#eta}");
   }

   TProfile* pf_sieie_cosThetaAllJets_dec[NDEC];
   for(int i=0; i<NDEC; i++){

     pf_sieie_cosThetaAllJets_dec[i] = (TProfile*)pf_sieie_cosTheta_template->Clone(Form("pf_sieie_cosThetaAllJets_dec%d",i));
     pf_sieie_cosThetaAllJets_dec[i] -> SetXTitle("cos#theta^{*}(#gamma,all jets)");
     pf_sieie_cosThetaAllJets_dec[i] -> SetYTitle("#sigma_{i#eta i#eta}");

   }



   // chi
   TH1F* h_chi_template = new TH1F("h_chi_template","",
					100,0.0,20.0);      
      




   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      nPass[0]++;


      // event vertex selections
      if(nVtxNotFake<0)continue;
      nPass[1]++;

      Bool_t hasGoodVtx = false;
      for(int iv=0; iv< nVtxNotFake; iv++)
	{

	  if(vertexNdof[iv] <=4)continue;

	  if(fabs(vertexZ[iv]) > 24.0)continue;
	  
	  if(sqrt(vertexX[iv]*vertexX[iv]+vertexY[iv]*vertexY[iv]) > 2.0)
	    continue;

	  hasGoodVtx = true;

	}
      
      if(!hasGoodVtx)continue;

      nPass[2]++;


      // first, look for a good leading photon
      
      Int_t   nGoodPho = 0;
      Int_t   leadingPhotonIndex = -1;
      Float_t ptMax = -9999.0;
      for(int ipho=0; ipho<nPho; ipho++){

	if(PhoisGenMatched[ipho]<0.5)continue; 

	//	if(fabs(PhogenMomId[ipho]-22.0)>1e-3)continue;

	// is not matched to any generator-level photon

	// asking for the photon to satisfy the ID selections

	Float_t pt    = PhoPt[ipho];
	Float_t sceta = PhoScEta[ipho];

	
	if(pt < 30.0)continue; 
	if( (fabs(sceta) > 1.4442 && fabs(sceta) < 1.566) 
	    || fabs(sceta) > 2.5)continue;

	if(PhohadronicOverEm[ipho] > 0.05)continue;

	if(PhoecalRecHitSumEtConeDR04[ipho]  > 4.2 + 0.003 * pt)continue;

	if(PhohcalTowerSumEtConeDR04[ipho] > 2.2 + 0.001*pt)continue;

	if(PhotrkSumPtHollowConeDR04[ipho] > 2.0 + 0.001*pt)continue;

	nGoodPho ++;
	
	if(pt > ptMax)
	  {
	    ptMax                = pt;
	    leadingPhotonIndex   = ipho;
	  }

	// ask if this photon passes selections
	
      } // end of loop over photons
      if(leadingPhotonIndex < 0 || ptMax < 0) continue;
      nPass[3]++;
      
      
      Float_t scetapho = PhoScEta[leadingPhotonIndex];
      Int_t   phoDecCode = -1;
      if(fabs(scetapho) < 1.4442) phoDecCode = 0;
      else if(fabs(scetapho) > 1.566 && fabs(scetapho)<2.5) phoDecCode = 1;
      if(phoDecCode <0)continue;
      
      nPass[4]++;


      TLorentzVector l_pho(0,0,0,0);
      l_pho.SetPtEtaPhiE(
			 PhoPt[leadingPhotonIndex],
			 PhoEta[leadingPhotonIndex],
			 PhoPhi[leadingPhotonIndex],
			 PhoEnergy[leadingPhotonIndex]
			 );



      // now, look for a leading jet
      Int_t   nGoodJet = 0;
      Int_t   leadingJetIndex = -1;
      Float_t ptMaxJet = -9999.0;


      // 4-momentum for summing up all the rest of the jets
      TLorentzVector l_allOtherJets(0,0,0,0);
      

      for(int ijet=0; ijet < nPatJet; ijet++){

	Float_t jetpt  = PatJetPt[ijet];
	Float_t jeteta = PatJetEta[ijet];

	if(jetpt < 10.0)continue;;
	if(fabs(jeteta) > 3.0)continue;
	if(PatJetNHF[ijet] >= 0.99)continue;
	if(PatJetNEF[ijet] >= 0.99)continue;

	// for the tracker region
	if(fabs(jeteta)<2.4 && PatJetCHF[ijet] <= 0.0)continue;
	if(fabs(jeteta)<2.4 && PatJetCEF[ijet] >= 0.99)continue;
	if(fabs(jeteta)<2.4 && PatJetNCH[ijet] <= 0)continue;

	TLorentzVector tempJet(0,0,0,0);
	tempJet.SetPtEtaPhiE(
			     PatJetPt[ijet],
			     PatJetEta[ijet],
			     PatJetPhi[ijet],
			     PatJetEn[ijet]
			     );


	Float_t dR = tempJet.DeltaR(l_pho);
	hdR_dec[phoDecCode]->Fill(dR);

	if(dR < 0.3)continue; // photon and jets are the same
	
        l_allOtherJets += tempJet;

	nGoodJet++;

	if(jetpt < 20.0)continue;;


	if(jetpt > ptMaxJet)
	  {
	    ptMaxJet             = jetpt;
	    leadingJetIndex      = ijet;
	  }


      } // end of loop over jets
      if(leadingJetIndex < 0 || ptMaxJet < 0)continue;
      nPass[5]++;

      Float_t ptpho  = PhoPt[leadingPhotonIndex];
      Float_t etapho = PhoEta[leadingPhotonIndex];

      hphopt->Fill(ptpho);
      hphoeta->Fill(etapho);	
      hphoetavspt->Fill(ptpho,etapho);

      hngoodpho->Fill(nGoodPho);
      hngoodjet->Fill(nGoodJet);


      TLorentzVector l_jet(0,0,0,0);
      l_jet.SetPtEtaPhiE(
			 PatJetPt[leadingJetIndex],
			 PatJetEta[leadingJetIndex],
			 PatJetPhi[leadingJetIndex],
			 PatJetEn[leadingJetIndex]
			 );



      Float_t sieie = PhoSigmaIetaIeta[leadingPhotonIndex];
      Float_t dPhi  = deltaPhi(l_pho,l_jet);
      Float_t dPhi_allotherjets  = deltaPhi(l_pho,l_allOtherJets);

      Float_t cosT  = cosThetaStar(l_pho,l_jet);
      Float_t cosT_allotherjets  = cosThetaStar(l_pho,l_allOtherJets);

      h_dPhi_sieie_dec[phoDecCode]->Fill(sieie,dPhi);
      h_dPhiAllJets_sieie_dec[phoDecCode]->Fill(sieie,dPhi_allotherjets);      
      hsieie_dec[phoDecCode]->Fill(sieie);

      
      pf_dPhi_sieie_dec[phoDecCode]->Fill(sieie,dPhi);
      pf_dPhiAllJets_sieie_dec[phoDecCode]->Fill(sieie,dPhi_allotherjets);      
      pf_sieie_dPhi_dec[phoDecCode]->Fill(dPhi, sieie);
      pf_sieie_dPhiAllJets_dec[phoDecCode]->Fill(dPhi_allotherjets, sieie);      


      pf_cosTheta_sieie_dec[phoDecCode]->Fill(sieie,cosT);
      pf_cosThetaAllJets_sieie_dec[phoDecCode]->Fill(sieie,cosT_allotherjets);      
      pf_sieie_cosTheta_dec[phoDecCode]->Fill(cosT, sieie);
      pf_sieie_cosThetaAllJets_dec[phoDecCode]->Fill(cosT_allotherjets, sieie);      

      
   } // end of loop over entries


   for(int i=0; i<50; i++)
     if(nPass[i]>0)cout << "nPass[" << i << "]=" << nPass[i] << endl;

   std::string _outputFile = _inputFile;
   std::string remword="/home/cdxfe2/Summer2011Ntuples/";
   size_t pos = _outputFile.find(remword);
   if(pos!= std::string::npos)
     _outputFile.swap(_outputFile.erase(pos,remword.length()));   

   gSystem->mkdir("yj276TeVhistos");

   TFile* _file = TFile::Open(Form("yj276TeVhistos/histo_%s",_outputFile.data()),"recreate");  
   
   hphopt->Write();
   hphoeta->Write();
   hphoetavspt->Write();
   hngoodpho->Write();
   hngoodjet->Write();


   for(int i=0; i<2; i++){

     hdR_dec[i]->Write();

     h_dPhi_sieie_dec[i]->Write();     
     h_dPhiAllJets_sieie_dec[i]->Write();     
     hsieie_dec[i]->Write();

     pf_dPhi_sieie_dec[i]->Write();     
     pf_dPhiAllJets_sieie_dec[i]->Write();     

     pf_sieie_dPhi_dec[i]->Write();     
     pf_sieie_dPhiAllJets_dec[i]->Write();     


     pf_cosTheta_sieie_dec[i]->Write();
     pf_cosThetaAllJets_sieie_dec[i]->Write();
     pf_sieie_cosTheta_dec[i]->Write();
     pf_sieie_cosThetaAllJets_dec[i]->Write();

   }

   _file->Close();

} // end of Loop
