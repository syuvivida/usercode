#define vectorJetEff_cxx
#include "vectorJetEff.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include "myLib.h"

void vectorJetEff::Loop()
{
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();

   cout << _inputFileName << " has " << nentries << " entries" << endl;

   // template
   TH1D* h_pt_template = new TH1D("h_pt_template","",500,0,500);
   TH1D* h_eta_template = new TH1D("h_eta_template","",100,-5.0,5.0);

   // book histograms
   TH1D* h_genjetpt_alljets[2]; // before and after the reconstruction
   TH1D* h_genjeteta_alljets[2];

   TH1D* h_recjetpt_alljets[2]; // before and after the ID cuts
   TH1D* h_recjeteta_alljets[2];


   // here the leading jet refers to the jet with highest gen jet pt
   TH1D* h_genjetpt_1stjet[2]; // before and after the reconstruction
   TH1D* h_genjeteta_1stjet[2];

   TH1D* h_recjetpt_1stjet[2]; // before and after the ID cuts
   TH1D* h_recjeteta_1stjet[2];

   for(int ip=0; ip<2; ip++){
     
     h_genjetpt_alljets[ip] = (TH1D*)h_pt_template->Clone(
			       Form("h_genjetpt_alljets_%d",ip));
     h_genjetpt_alljets[ip] -> SetXTitle("p_{T}^{GEN}(all jets)");

     h_recjetpt_alljets[ip] = (TH1D*)h_pt_template->Clone(
                               Form("h_recjetpt_alljets_%d",ip));
     h_recjetpt_alljets[ip] -> SetXTitle("p_{T}^{REC}(all jets)");

     h_genjeteta_alljets[ip] = (TH1D*)h_eta_template->Clone(
                               Form("h_genjeteta_alljets_%d",ip));
     h_genjeteta_alljets[ip] -> SetXTitle("#eta^{GEN}(all jets)");

     h_recjeteta_alljets[ip] = (TH1D*)h_eta_template->Clone(
			       Form("h_recjeteta_alljets_%d",ip));
     h_recjeteta_alljets[ip] -> SetXTitle("#eta^{REC}(all jets)");


     h_genjetpt_1stjet[ip] = (TH1D*)h_pt_template->Clone(
			       Form("h_genjetpt_1stjet_%d",ip));
     h_genjetpt_1stjet[ip] -> SetXTitle("p_{T}^{GEN}(leading gen jet)");


     h_recjetpt_1stjet[ip] = (TH1D*)h_pt_template->Clone(
                               Form("h_recjetpt_1stjet_%d",ip));
     h_recjetpt_1stjet[ip] -> SetXTitle("p_{T}^{REC}(leading gen jet)");

     h_genjeteta_1stjet[ip] = (TH1D*)h_eta_template->Clone(
                               Form("h_genjeteta_1stjet_%d",ip));
     h_genjeteta_1stjet[ip] -> SetXTitle("#eta^{GEN}(leading gen jet)");


     h_recjeteta_1stjet[ip] = (TH1D*)h_eta_template->Clone(
			       Form("h_recjeteta_1stjet_%d",ip));
     h_recjeteta_1stjet[ip] -> SetXTitle("#eta^{REC}(leading gen jet)");


						  
   }

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
//       if(jentry > 500 ) break;
      // find the genjet that has the highest jet pt
      // and also fill the histogram for reco and id efficiency of all jets that have gen pt > 30 GeV, |eta| < 3.0
      double genJetMaxPt = -9999.0;
      int leadingGenJetIndex = -1;
      int matchedRecoLeadingGenJetIndex = -1;

      // loop over generator-level jets
      for(unsigned int ij=0; ij < genJetPt_->size(); ij++)
	{
	  double thisJetPt = genJetPt_->at(ij);
	  double thisJetEta= genJetEta_->at(ij);

	  if(thisJetPt < 30.0)continue;
	  if(fabs(thisJetEta) > 3.0)continue;

	  // default gen jet
	  h_genjetpt_alljets[0] ->Fill(thisJetPt);
	  h_genjeteta_alljets[0]->Fill(thisJetEta);

	  // for those jets that are matched to recoJet
	  int recoJetIndex = matchGenToRecoJet(ij);
// 	  cout << "genjet index = " << ij << " and matched reco index = " << recoJetIndex << endl;

 	  if(recoJetIndex>=0)
 	    {
 	      // for jet reconstruction efficiency
 	      h_genjetpt_alljets[1] ->Fill(thisJetPt);
 	      h_genjeteta_alljets[1]->Fill(thisJetEta);

 	      // then further go ahead to calculate loose jet id efficiency
	      
 	      double thisRecoJetPt = patJetPfAk05Pt_->at(recoJetIndex);
 	      double thisRecoJetEta= patJetPfAk05Eta_->at(recoJetIndex);

 	      // should I apply jet fiducial requirement to calculate ID efficiency?
	      if(isFidJet(recoJetIndex)){

		h_recjetpt_alljets[0] ->Fill(thisRecoJetPt);
		h_recjeteta_alljets[0]->Fill(thisRecoJetEta);

		if(isGoodLooseJet(recoJetIndex))
		  {
		  
		    h_recjetpt_alljets[1] ->Fill(thisRecoJetPt);
		    h_recjeteta_alljets[1]->Fill(thisRecoJetEta);

		  }
	      } // if the reconstructed jet satisfies some fiducial cuts

 	    } // if we find a matched recoJet
	  
 	  // try to find the genJet that has highest pt
 	  if(thisJetPt > genJetMaxPt)
 	    {
 	      genJetMaxPt                   = thisJetPt;
 	      leadingGenJetIndex            = ij;
 	      // index of the reconstruction-level jet that is matched to this genjet
 	      matchedRecoLeadingGenJetIndex = recoJetIndex; 
 	    }


	} // end of loop over genJets


      // now check reconstruction and id efficiency of leading genJet

      if(leadingGenJetIndex < 0 )continue;
      
      double leadingGenJetPt = genJetPt_->at(leadingGenJetIndex);
      double leadingGenJetEta= genJetEta_->at(leadingGenJetIndex);
      h_genjetpt_1stjet[0] ->Fill(leadingGenJetPt);
      h_genjeteta_1stjet[0]->Fill(leadingGenJetEta);

      if(matchedRecoLeadingGenJetIndex >= 0)
 	{
 	  h_genjetpt_1stjet[1] ->Fill(leadingGenJetPt);
 	  h_genjeteta_1stjet[1]->Fill(leadingGenJetEta);

 	  double leadingRecoJetPt = patJetPfAk05Pt_->at(matchedRecoLeadingGenJetIndex);
 	  double leadingRecoJetEta= patJetPfAk05Eta_->at(matchedRecoLeadingGenJetIndex);

 	  // should I apply jet fiducial requirement to calculate ID efficiency?
	  if(isFidJet(matchedRecoLeadingGenJetIndex)){

	    h_recjetpt_1stjet[0] ->Fill(leadingRecoJetPt);
	    h_recjeteta_1stjet[0]->Fill(leadingRecoJetEta);

	    if(isGoodLooseJet(matchedRecoLeadingGenJetIndex))
	      {
		  
		h_recjetpt_1stjet[1] ->Fill(leadingRecoJetPt);
		h_recjeteta_1stjet[1]->Fill(leadingRecoJetEta);
	      
	      }
	  
	  } // if the reconstructed jet satisfies some fiducial cuts

 	} // if the leading genjet finds a matched reco jet

      


   } // end of loop over entries

   
   // write output root file
   string remword="/data4/syu/7TeV_vectorNtuple/pythia/";

   size_t pos = _inputFileName.find(remword);

   if(pos!= string::npos)
     _inputFileName.swap(_inputFileName.erase(pos,remword.length()));

   TFile* outFile = new TFile(Form("/home/syu/CVSCode/jeteff_%s",
 				   _inputFileName.data()),"recreate");       

   
   for(int ip=0; ip<2; ip++){
     
     h_genjetpt_alljets[ip]->Write(); 
     h_genjeteta_alljets[ip]->Write();
      
     h_recjetpt_alljets[ip]->Write(); 
     h_recjeteta_alljets[ip]->Write();

     h_genjetpt_1stjet[ip]->Write(); 
     h_genjeteta_1stjet[ip]->Write();

     h_recjetpt_1stjet[ip]->Write(); 
     h_recjeteta_1stjet[ip]->Write();

   }
   outFile->Close();     
   

}


// check if a jet satisfy basic kinematic requirement
Bool_t vectorJetEff::isFidJet (Int_t ijet)
{
  if(fabs(patJetPfAk05Eta_->at(ijet)) > 3.0)return false;
  if(patJetPfAk05Pt_->at(ijet) < 20.0)return false;
  return true;
}


// check if this reco-jet is a good loose jet
Bool_t vectorJetEff::isGoodLooseJet(Int_t ijet)
{
  if(!isFidJet(ijet))return false;
  if(patJetPfAk05JetNConstituents_->at(ijet) <= 1)return false;
  if(patJetPfAk05NeutHadEFr_->at(ijet) >= 0.99)return false;
  if(patJetPfAk05NeutEmEFr_->at(ijet) >= 0.99)return false;

  //   // for the tracker region
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharHadEFr_->at(ijet) <= 0.0)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharEmEFr_->at(ijet) >= 0.99)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharMulti_->at(ijet) <= 0)return false;

  return true;

}


// check if this reco-jet is a good medium jet
Bool_t vectorJetEff::isGoodMediumJet(Int_t ijet)
{
  if(patJetPfAk05Pt_->at(ijet) < 10.0)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet)) > 3.0)return false;
  if(patJetPfAk05JetNConstituents_->at(ijet) <= 1)return false;
  if(patJetPfAk05NeutHadEFr_->at(ijet) >= 0.95)return false;
  if(patJetPfAk05NeutEmEFr_->at(ijet) >= 0.95)return false;

  //   // for the tracker region
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharHadEFr_->at(ijet) <= 0.0)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharEmEFr_->at(ijet) >= 0.99)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharMulti_->at(ijet) <= 0)return false;


  return true;

}

// check if this reco-jet is a good medium jet
Bool_t vectorJetEff::isGoodTightJet(Int_t ijet)
{
  if(patJetPfAk05Pt_->at(ijet) < 10.0)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet)) > 3.0)return false;
  if(patJetPfAk05JetNConstituents_->at(ijet) <= 1)return false;
  if(patJetPfAk05NeutHadEFr_->at(ijet) >= 0.90)return false;
  if(patJetPfAk05NeutEmEFr_->at(ijet) >= 0.90)return false;


  //   // for the tracker region
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharHadEFr_->at(ijet) <= 0.0)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharEmEFr_->at(ijet) >= 0.99)return false;
  if(fabs(patJetPfAk05Eta_->at(ijet))<2.4 && patJetPfAk05CharMulti_->at(ijet) <= 0)return false;

  return true;

}


Int_t vectorJetEff::matchRecoToGenJet(Int_t ijet)
{  
  int matchedGenIndex = -1;
  for(unsigned int k=0; k< genJetPt_->size(); k++)
    {

      double dR = eiko::deltaR(genJetEta_->at(k),genJetPhi_->at(k),
			       patJetPfAk05Eta_->at(ijet),patJetPfAk05Phi_->at(ijet)); 

      double relPt = genJetPt_->at(k)>1e-6? 
	fabs(genJetPt_->at(k)-patJetPfAk05Pt_->at(ijet))/genJetPt_->at(k): -9999.0;

      if(dR<0.4 && relPt < 3.0)
	{
	  matchedGenIndex = k;
	  break;
	}
      
    } // end of loop over generator-level jets
      
  return matchedGenIndex;
}


Int_t vectorJetEff::matchGenToRecoJet(Int_t ijet)
{
  int matchedRecoIndex = -1;
  for(unsigned int k=0; k< patJetPfAk05Pt_->size(); k++)
    {

      double dR = eiko::deltaR(genJetEta_->at(ijet),genJetPhi_->at(ijet),
			       patJetPfAk05Eta_->at(k),patJetPfAk05Phi_->at(k)); 

      double relPt = genJetPt_->at(ijet)>1e-6? 
	fabs(genJetPt_->at(ijet)-patJetPfAk05Pt_->at(k))/genJetPt_->at(ijet): -9999.0;

      if(dR<0.4 && relPt < 3.0)
	{
	  matchedRecoIndex = k;
	  break;
	}
      
    } // end of loop over generator-level jets
      
  return matchedRecoIndex;

}
