#define LHE_angular_cxx
#include "LHE_angular.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

const double minZPt  = 40.0;

const double minPartonPt=30.0;
const double maxPartonEta=2.4;
const double mindR=0.5;

const double minElePt = 20.0;
const double minEleBarrelEta = 0.0;
const double maxEleBarrelEta = 1.442;

const double minEleEndcapEta = 1.566;
const double maxEleEndcapEta = 2.1;
const double minMee = 76.0;
const double maxMee =106.0;

const double minMuoPt = 20.0;
const double maxMuoEta = 2.1;
const double minMmm = 76.0;
const double maxMmm =106.0;

void LHE_angular::Loop(int lepID, bool exclusive)
{
  if (fChain == 0) return;

  //=============================================================================
  //   Print out information
  //=============================================================================

  cout << "=========================================================" << endl;

  Long64_t nentries = fChain->GetEntriesFast();

  cout << _inputFileName << " has " << nentries << " entries." << endl;
  std::string leptonName;
  if(abs(lepID)==13)leptonName = "muon";
  else if(abs(lepID)==11)leptonName = "electron";

  cout << "Running mode: " << endl;
  cout << "exclusive=" << exclusive << endl;
  cout << "The cuts applied: " << endl;
  
  cout << " minZPt= " << minZPt << endl;
  cout << " minPartonPt= " << minPartonPt << endl;
  cout << " maxPartonEta= " << maxPartonEta << endl;
  cout << " mindR= " << mindR << endl;

  cout << "studying " << leptonName << endl;
  if(abs(lepID)==11){
    cout << " minElePt= " << minElePt << endl;
    cout << " minEleBarrelEta = " << minEleBarrelEta << endl;
    cout << " maxEleBarrelEta = " << maxEleBarrelEta << endl;
    cout << " minEleEndcapEta = " << minEleEndcapEta << endl;
    cout << " maxEleEndcapEta = " << maxEleEndcapEta << endl;
    cout << " minMee = " << minMee << endl;
    cout << " maxMee = " << maxMee << endl;
  }
  else if(abs(lepID)==13){
    cout << " minMuoPt = " << minMuoPt << endl;
    cout << " maxMuoEta = " << maxMuoEta << endl;
    cout << " minMmm = " << minMmm << endl;
    cout << " maxMmm = " << maxMmm << endl;
  }
  cout << "=========================================================" << endl;

  // dummy proof, in case someone put a negative number
  int leptonPID = abs(lepID); 

  //=============================================================================
  //   Book histograms
  //=============================================================================


  TH1D* h_mZ   = new TH1D("h_mZ","",200,20.0,220.0);
  h_mZ->SetXTitle("M_{ll} [GeV/c^{2}");
  h_mZ->Sumw2();


  TH1D* h_njet = new TH1D("h_njet","",6,-0.5,5.5);
  h_njet->SetXTitle("#geq n jet");
  h_njet->Sumw2();


  TH1D* h_zpt_template = new TH1D("h_zpt_template","",40,0,400);
  h_zpt_template->Sumw2();

  TH1D* h_zpt   = (TH1D*)h_zpt_template->Clone("h_zpt");
  h_zpt->SetXTitle("p_{T}(Z) [GeV]");

  TH1D* h_jetpt_template = new TH1D("h_jetpt_template","",40,0,400);
  h_jetpt_template->Sumw2();

  TH1D* h_jetpt = (TH1D*)h_jetpt_template->Clone("h_jetpt");
  h_jetpt->SetXTitle("p_{T}(jet) [GeV]");

  TH1D* h_ypart_template = new TH1D("h_ypart_template","",15,0.0,3.0); 
  h_ypart_template->Sumw2();

  TH1D* h_zy    = (TH1D*)h_ypart_template->Clone("h_zy");
  h_zy->SetXTitle("|y_{Z}|");

  TH1D* h_jety    = (TH1D*)h_ypart_template->Clone("h_jety");
  h_jety->SetXTitle("|y_{jet}|");

  TH1D* h_y_template = new TH1D("h_y_template","",15,0.0,3.0);
  h_y_template->Sumw2();

  TH1D* h_yB    = (TH1D*)h_y_template->Clone("h_yB");
  h_yB->SetXTitle("0.5|y_{Z}+y_{jet}|");

  TH1D* h_ystar = (TH1D*)h_y_template->Clone("h_ystar");
  h_ystar->SetXTitle("0.5|y_{Z}-y_{jet}|");


  //=============================================================================
  //   Loop over entries
  //=============================================================================


  int nPass[30];
   
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
//     if(jentry > 1000) break;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    nPass[0]++;

    double eventWeight = 1.0;
    // now loop over particles to find leptons
    int lepPlusIndex = -1;
    int lepMinusIndex = -1;

    for(unsigned int igen=0; igen < genParId_->size(); igen++){

      int PID       = genParId_->at(igen);
      int status    = genParSt_->at(igen);

      // first check leptons

      bool isLepPlus = (PID == (-leptonPID) && status==1);	
      bool isLepMinus= (PID == ( leptonPID) && status==1);


      if(!isLepPlus && !isLepMinus)continue;

      TLorentzVector thisLepton(0,0,0,0);
      thisLepton.SetPxPyPzE(genParPx_->at(igen),
			    genParPy_->at(igen),
			    genParPz_->at(igen),
			    genParE_->at(igen));

      double pt  = thisLepton.Pt();
      double eta = thisLepton.Eta();

      // muon
      if(leptonPID==13 && (pt < minMuoPt || fabs(eta) > maxMuoEta))continue;
      if(leptonPID==11 && (pt < minElePt || !( 
					      (fabs(eta) > minEleEndcapEta && 
					       fabs(eta) < maxEleEndcapEta) ||
					      (fabs(eta) > minEleBarrelEta && 
					       fabs(eta) < maxEleBarrelEta
					       )))
	 )continue;

      if(lepPlusIndex < 0  && isLepPlus)lepPlusIndex = igen;

      if(lepMinusIndex < 0 && isLepMinus)lepMinusIndex = igen;
      
	
    }

    // do not find lep+ or lep- that satisfies cuts
    if(lepPlusIndex < 0 || lepMinusIndex < 0)continue;

    nPass[1]++;


    TLorentzVector l4_lepp(0,0,0,0);
    l4_lepp.SetPxPyPzE(
		       genParPx_->at(lepPlusIndex),
		       genParPy_->at(lepPlusIndex),
		       genParPz_->at(lepPlusIndex),
		       genParE_->at(lepPlusIndex)
		       );


    TLorentzVector l4_lepm(0,0,0,0);
    l4_lepm.SetPxPyPzE(
		       genParPx_->at(lepMinusIndex),
		       genParPy_->at(lepMinusIndex),
		       genParPz_->at(lepMinusIndex),
		       genParE_->at(lepMinusIndex)
		       );



    TLorentzVector l4_z = l4_lepp + l4_lepm;


    double mll = l4_z.M();

    h_mZ->Fill(mll,eventWeight);

    if(leptonPID==13 && (mll < minMmm || mll > maxMmm))continue;
    if(leptonPID==11 && (mll < minMee || mll > maxMee))continue;

    nPass[2]++;


    // now look for partons

    int leadingPartonIndex = -1;
    unsigned int nPartons = 0;
    double maxPartonPt = -9999;

    for(unsigned int igen=0; igen < genParId_->size(); igen++){

      int status    = genParSt_->at(igen);
      int PID       = genParId_->at(igen);

      if(status!=1)continue;
      if(abs(PID)>6 && PID!=21)continue;

      cout << "PID = " << PID << endl;
      


      TLorentzVector thisParton_l4(0,0,0,0);
      thisParton_l4.SetPxPyPzE(
			       genParPx_->at(igen),
			       genParPy_->at(igen),
			       genParPz_->at(igen),
			       genParE_->at(igen)
			       );

      double thisPartonPt = thisParton_l4.Pt();
      double thisPartonEta= thisParton_l4.Eta();

      if(thisPartonPt < minPartonPt)continue;
      if(fabs(thisPartonEta) > maxPartonEta)continue;

      double dr_ep = l4_lepp.DeltaR(thisParton_l4);
      double dr_em = l4_lepm.DeltaR(thisParton_l4);

      if(dr_ep < mindR)continue;
      if(dr_em < mindR)continue;
      nPartons++; 

      if(thisPartonPt > maxPartonPt)
	{
	  maxPartonPt = thisPartonPt;
	  leadingPartonIndex = igen;
	}
	
    } // end of loop over jets
      

    for(int ik=0; ik<h_njet->GetNbinsX(); ik++)
      if(nPartons>= ik)
	h_njet->Fill(ik,eventWeight);

    if(leadingPartonIndex < 0)continue;
    
    cout << "Leading parton index = " << leadingPartonIndex << endl;
    nPass[4]++;

    if(exclusive && nPartons!=1)continue;
     
    nPass[5]++;


    double ptz = l4_z.Pt();

    if(ptz < minZPt)continue;

 
    TLorentzVector l4_j(0,0,0,0);
    l4_j.SetPxPyPzE(
		    genParPx_->at(leadingPartonIndex),
		    genParPy_->at(leadingPartonIndex),
		    genParPz_->at(leadingPartonIndex),
		    genParE_->at(leadingPartonIndex)
		    );


    double ptjet = l4_j.Pt();

    double yz = l4_z.Rapidity();
    double yj = l4_j.Rapidity();

    double yB = 0.5*(yz + yj);
    double ystar = 0.5*(yz-yj);

    h_zpt->Fill(ptz,eventWeight);
    h_jetpt->Fill(ptjet,eventWeight); 
    h_zy->Fill(fabs(yz),eventWeight);
    h_jety->Fill(fabs(yj),eventWeight);
    h_yB->Fill(fabs(yB),eventWeight);
    h_ystar->Fill(fabs(ystar),eventWeight);

  } // end of loop over entries

  std::string prefix = "LHE";
  if(exclusive)prefix += "_exclusive1Parton";
  if(minZPt>1e-6)prefix += Form("_zPt%d",(int)minZPt);

  TFile* outFile = new TFile(Form("%s_%s_%s",prefix.data(),
				  leptonName.data(),
				  _inputFileName.data()),"recreate");       
  h_mZ->Write();
  h_njet->Write();
        
  h_yB->Write();
  h_ystar->Write();
  
  h_zy->Write();
  h_zpt->Write();

  h_jety->Write();   
  h_jetpt->Write();

  outFile->Close();


}
