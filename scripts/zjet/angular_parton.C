#define angular_parton_cxx
#include "angular_parton.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <map>

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


void angular_parton::Loop(int lepID, bool exclusive, bool applyWeight, int DEBUG)
{
  if (fChain == 0) return;

  //=============================================================================
  //   Print out information
  //=============================================================================

  cout << "=========================================================" << endl;

  bool isSherpa=false;
  size_t pos_sherpa  = _inputFileName.find("sherpa");
  if(pos_sherpa!= std::string::npos)
    isSherpa=true;
  if(isSherpa)cout << "This is a sherpa MC sample" << endl;

  bool isMadgraph=false;
  size_t pos_madgraph  = _inputFileName.find("madgraph");
  if(pos_madgraph!= std::string::npos)
    isMadgraph=true;
  if(isMadgraph)cout << "This is a madgraph MC sample" << endl;

  bool isPythia=false;
  size_t pos_pythia  = _inputFileName.find("pythia");
  if(pos_pythia!= std::string::npos)
    isPythia=true;
  if(isPythia)cout << "This is a pythia MC sample" << endl;

  const double fBinsPt01[]= {30,40,55,75,105,150,210,315,500};
  const double fBinsPt02[]= {30,40,55,75,105,150,210,315,450};
  const double fBinsPt03[]= {30,40,55,75,105,150,300};
  const double fBinsPt04[]= {30,40,55,75,105,200};
  
  int nPtBins01 = sizeof(fBinsPt01)/sizeof(fBinsPt01[0])-1;
  cout << "There are " << nPtBins01 << " bins in 1st leading jet" << endl;

  int nPtBins02 = sizeof(fBinsPt02)/sizeof(fBinsPt02[0])-1;
  cout << "There are " << nPtBins02 << " bins in 2nd leading jet" << endl;

  int nPtBins03 = sizeof(fBinsPt03)/sizeof(fBinsPt03[0])-1;
  cout << "There are " << nPtBins03 << " bins in 3rd leading jet" << endl;

  int nPtBins04 = sizeof(fBinsPt04)/sizeof(fBinsPt04[0])-1;
  cout << "There are " << nPtBins04 << " bins in 4th leading jet" << endl;

  const double fBinsY[]={0.0,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4};
  const int nYBins = sizeof(fBinsY)/sizeof(fBinsY[0])-1;

  cout << "There are " << nYBins << " bins." << endl;

  Long64_t nentries = fChain->GetEntriesFast();

  cout << _inputFileName << " has " << nentries << " entries." << endl;
  std::string leptonName;
  if(abs(lepID)==13)leptonName = "muon";
  else if(abs(lepID)==11)leptonName = "electron";

  cout << "Running mode: " << endl;
  cout << "exclusive=" << exclusive 
       << "\t applyWeight=" <<  applyWeight << "\t DEBUG=" << DEBUG << endl;
  cout << endl;
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

  TH1D* h_nvtx = new TH1D("h_nvtx","",41.5,-0.5,40.5);
  h_nvtx->SetXTitle("Number of good vertices");
  h_nvtx->Sumw2();

  TH1D* h_njet = new TH1D("h_njet","",6,-0.5,5.5);
  h_njet->SetXTitle("#geq n jet");
  h_njet->Sumw2();

  const int nMAXJETS=4;

  TH1D* h_jetpt_power_template = new TH1D("h_jetpt_power_template","",500,30,530);
  h_jetpt_power_template->SetXTitle("p_{T}(jet)^{gen} [GeV]");
  h_jetpt_power_template->Sumw2();
  TH1D* h_jetpt_power[nMAXJETS+1];
  TH1D* h_jetp_power[nMAXJETS+1];
   
  TH1D* h_mc_jetpt[nMAXJETS];
  TH1D* h_mc_jety[nMAXJETS];

  TH1D* h_predict_jetpt_template = new TH1D("h_predict_jetpt_template","",nPtBins01,fBinsPt01);
  h_predict_jetpt_template->Sumw2();
  TH1D* h_leadingjet_pt = (TH1D*)h_predict_jetpt_template->Clone("h_leadingjet_pt");
   

  h_mc_jetpt[0] = new TH1D(Form("h_mc_jetpt%02i",1),"1st leading jet",nPtBins01,fBinsPt01);
  h_mc_jetpt[1] = new TH1D(Form("h_mc_jetpt%02i",2),"2nd leading jet",nPtBins02,fBinsPt02);
  h_mc_jetpt[2] = new TH1D(Form("h_mc_jetpt%02i",3),"3rd leading jet",nPtBins03,fBinsPt03);
  h_mc_jetpt[3] = new TH1D(Form("h_mc_jetpt%02i",4),"4th leading jet",nPtBins04,fBinsPt04);

  TH1D* h_predict_jety_template = new TH1D("h_predict_jety_template","",nYBins,fBinsY);
  h_predict_jety_template->Sumw2();

  TH1D* h_leadingjet_y = (TH1D*)h_predict_jety_template->Clone("h_leadingjet_y");

  for(int ij=0; ij<nMAXJETS; ij++){

    h_jetpt_power[ij] = (TH1D*)h_jetpt_power_template->Clone(Form("h_jetpt_power%02i",ij+1));
    h_jetp_power[ij]  = (TH1D*)h_jetpt_power_template->Clone(Form("h_jetp_power%02i",ij+1));
    h_jetp_power[ij]  -> SetXTitle("p(jet) [GeV]");

    h_mc_jetpt[ij]->SetXTitle("p_{T}(jet) [GeV]");
    h_mc_jetpt[ij]->Sumw2();

    h_mc_jety[ij] = (TH1D*)h_predict_jety_template->Clone(Form("h_mc_jety%02i",ij+1));
    h_mc_jety[ij]->SetXTitle("Rapidity(jet)");

  }

  h_jetpt_power[4] = (TH1D*)h_jetpt_power_template->Clone("h_jetpt_power_inclusive");

  h_jetp_power[4] = (TH1D*)h_jetpt_power_template->Clone("h_jetp_power_inclusive");
  h_jetp_power[4]->SetXTitle("p(jet) [GeV]");

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
   
  typedef map<double, int,  std::greater<double> > myMap;
  myMap sorted_partonEtMap;
  typedef myMap::iterator mapIter;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    sorted_partonEtMap.clear();
//     if(jentry > 1000) break;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    nPass[0]++;

    double eventWeight = 1;

    if(applyWeight && PU_weight >= 0.0)eventWeight *= PU_weight;
    if(applyWeight && mcWeight_>0)eventWeight *= mcWeight_;
      
    if(DEBUG==1){
      cout << "PU_weight = " << PU_weight << "\t nvertex = " << EvtInfo_NumVtx << endl;
      cout << "MCweight = " << mcWeight_ << endl;
      cout << "eventWeight = " << eventWeight << endl;
    }

    h_nvtx->Fill(EvtInfo_NumVtx, eventWeight);
 
    // now loop over particles to find leptons
    int lepPlusIndex = -1;
    int lepMinusIndex = -1;

    for(unsigned int igen=0; igen < genParId_->size(); igen++){

      int PID       = genParId_->at(igen);
      int motherPID = genMomParId_->at(igen);
      int status    = genParSt_->at(igen);

      // first check leptons
      bool isFromZBeforeFSR = ((isMadgraph && motherPID == 23) ||
			       (isPythia && motherPID == 23) ||
 			       (isSherpa && motherPID == 10002)) 
//			       (isSherpa && motherPID < -9998)) 
	&& (status == 3);   

      bool isLepPlus = (PID == (-leptonPID) && isFromZBeforeFSR);	
      bool isLepMinus= (PID == ( leptonPID) && isFromZBeforeFSR);


      if(!isLepPlus && !isLepMinus)continue;

      double pt  = genParPt_->at(igen);
      double eta = genParEta_->at(igen);

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

    if(DEBUG==1){
      cout << "lepPlusIndex = " << lepPlusIndex << "\t pt=" << 
	genParPt_->at(lepPlusIndex) << "\t eta=" << 
	genParEta_->at(lepPlusIndex) << endl;
      cout << "lepMinusIndex = " << lepMinusIndex << "\t pt=" << 
	genParPt_->at(lepMinusIndex) << "\t eta=" << 
	genParEta_->at(lepMinusIndex) << endl;
    }


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

    h_mZ->Fill(mll);

    if(leptonPID==13 && (mll < minMmm || mll > maxMmm))continue;
    if(leptonPID==11 && (mll < minMee || mll > maxMee))continue;

    if(DEBUG==1)
      cout << "dilepton mass = " << mll << endl;
      
    nPass[2]++;


    // now look for partons

    int leadingPartonIndex = -1;
    unsigned int nPartons = 0;
    double maxPartonPt = -9999;

    for(unsigned int igen=0; igen < genParId_->size(); igen++){

      int status    = genParSt_->at(igen);
      int PID       = genParId_->at(igen);
      int motherPID = genMomParId_->at(igen);

      if(status!=3)continue;
      if(abs(PID)>6 && PID!=21)continue;
      if(motherPID!=10002)continue;

//       cout << "PID = " << PID << endl;
      
      double thisPartonPt  = genParPt_->at(igen);
      double thisPartonEta = genParEta_->at(igen);

      if(thisPartonPt < minPartonPt)continue;
      if(fabs(thisPartonEta) > maxPartonEta)continue;

      TLorentzVector thisParton_l4(0,0,0,0);
	

      thisParton_l4.SetPtEtaPhiE(genParPt_->at(igen),
				 genParEta_->at(igen),
				 genParPhi_->at(igen),
				 genParE_->at(igen));


      double dr_ep = l4_lepp.DeltaR(thisParton_l4);
      double dr_em = l4_lepm.DeltaR(thisParton_l4);

      if(dr_ep < mindR)continue;
      if(dr_em < mindR)continue;
      nPartons++; 
      sorted_partonEtMap.insert(std::pair<double, int>(thisPartonPt,igen));  

      if(thisPartonPt > maxPartonPt)
	{
	  maxPartonPt = thisPartonPt;
	  leadingPartonIndex = igen;
	}

      h_jetpt_power[4]->Fill(thisPartonPt,eventWeight);
      h_jetp_power[4]->Fill(thisParton_l4.P(),eventWeight);
	
    } // end of loop over jets
      
    if(nPartons!= sorted_partonEtMap.size())
      { 
	cout << "errors in map size" << endl;
	continue;
      }

    if(DEBUG==1)
      cout << "max gen jet index = " << leadingPartonIndex << endl;

    for(unsigned int ik=0; ik<h_njet->GetNbinsX(); ik++)
      if(nPartons>= ik)
	h_njet->Fill(ik,eventWeight);

    if(leadingPartonIndex < 0)continue;

    nPass[4]++;

    if(exclusive && nPartons!=1)continue;
     
    nPass[5]++;


    double ptz = l4_z.Pt();

    if(ptz < minZPt)continue;

 
    TLorentzVector l4_j(0,0,0,0);
    l4_j.SetPtEtaPhiE(genParPt_->at(leadingPartonIndex),
		      genParEta_->at(leadingPartonIndex),
		      genParPhi_->at(leadingPartonIndex),
		      genParE_->at(leadingPartonIndex));

    double ptjet = l4_j.Pt();

    double yz = l4_z.Rapidity();
    double yj = l4_j.Rapidity();

    double yB = 0.5*(yz + yj);
    double ystar = 0.5*(yz-yj);

    h_leadingjet_pt->Fill(ptjet,eventWeight);
    h_leadingjet_y->Fill(fabs(yj),eventWeight);
      
      
    if(DEBUG==1)
      cout << "Now ordering jets" << endl;
    int countParton=0;
    for (mapIter it_part= sorted_partonEtMap.begin();
	 it_part != sorted_partonEtMap.end() && countParton< nMAXJETS; ++it_part)
      {
        double pt_mapthis = genParPt_->at(it_part->second);

	int jet_index = it_part->second;

	TLorentzVector l4_jthis(0,0,0,0);
	l4_jthis.SetPtEtaPhiE(genParPt_->at(jet_index),
			      genParEta_->at(jet_index),
			      genParPhi_->at(jet_index),
			      genParE_->at(jet_index));

	double y_mapthis  = l4_jthis.Rapidity();
	
	if(DEBUG==1)
	  cout << "pt = " << pt_mapthis << "\t y = " << y_mapthis << endl;

        h_mc_jetpt[countParton]->Fill(pt_mapthis,eventWeight);
        h_mc_jety[countParton]->Fill(fabs(y_mapthis),eventWeight);

	h_jetpt_power[countParton]->Fill(pt_mapthis,eventWeight);	
	h_jetp_power[countParton] ->Fill(l4_jthis.P(),eventWeight);

        countParton++;
      }

    if(DEBUG==1)cout << "nPartons = " << sorted_partonEtMap.size() 
		     << "\t countParton = " << countParton << endl;

    h_zpt->Fill(ptz,eventWeight);
    h_jetpt->Fill(ptjet,eventWeight); 
    h_zy->Fill(fabs(yz),eventWeight);
    h_jety->Fill(fabs(yj),eventWeight);
    h_yB->Fill(fabs(yB),eventWeight);
    h_ystar->Fill(fabs(ystar),eventWeight);

    if(DEBUG==1){
      double dR1 = l4_lepp.DeltaR(l4_j);
      double dR2 = l4_lepm.DeltaR(l4_j);
      cout << "dR1 = " << dR1 << "\t dR2=" << dR2 << endl;
    }
  } // end of loop over entries

  std::string prefix = "parton_weighted";
  if(!applyWeight)prefix = "parton_raw";
  if(exclusive)prefix += "_exclusive1Parton";
  if(minZPt>1e-6)prefix += Form("_zPt%d",(int)minZPt);

  std::string remword  ="/data2/syu/zjet_vectorNtuple/";

  size_t pos  = _inputFileName.find(remword);

  if(pos!= std::string::npos)
    _inputFileName.swap(_inputFileName.erase(pos,remword.length()));


  TFile* outFile = new TFile(Form("%s_%s_%s",prefix.data(),
				  leptonName.data(),
				  _inputFileName.data()),"recreate");       
  h_mZ->Write();
  h_nvtx->Write();
  h_njet->Write();
        
  h_leadingjet_pt->Write();
  h_leadingjet_y->Write();
  h_yB->Write();
  h_ystar->Write();
  
  h_zy->Write();
  h_zpt->Write();

  h_jety->Write();   
  h_jetpt->Write();

  for(int ij=0; ij < nMAXJETS; ij++){
    h_mc_jetpt[ij]->Write();
    h_mc_jety[ij]->Write();
  }

  for(int ij=0; ij < nMAXJETS+1; ij++){
    h_jetpt_power[ij]->Write();
    h_jetp_power[ij]->Write();
  }
  outFile->Close();

}
