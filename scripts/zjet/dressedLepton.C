#define dressedLepton_cxx
#include "dressedLepton.h"
#include "SMP-12-017.h" // for Z+jet cross section
//#include "SMP-12-004.h" // for angular analysis
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <map>



void dressedLepton::Loop(int lepID, int mode, bool exclusive, 
			 int DEBUG)
{

  //============================================================================
  //   Print out information
  //============================================================================

  if (fChain == 0) return;
  cout << "=========================================================" << endl;

  if(mode==0)
    cout << "bare leptons" << endl;
  else if(mode==1)
    cout << "dressed leptons" << endl;
  else if(mode==2)
    cout << "both bare and dressed leptons" << endl;

  std::string leptonName;
  if(abs(lepID)==13)leptonName = "muon";
  else if(abs(lepID)==11)leptonName = "electron";
  else if(abs(lepID)==0)leptonName = "both";

  std::string mcName;  

  bool isSherpa=false;
  size_t pos_sherpa  = _inputFileName.find("sherpa");
  if(pos_sherpa!= std::string::npos)
    isSherpa=true;
  if(isSherpa){mcName = "sherpa"; cout << "This is a sherpa MC sample" << endl;}

  bool isMadgraph=false;
  size_t pos_madgraph  = _inputFileName.find("madgraph");
  if(pos_madgraph!= std::string::npos)
    isMadgraph=true;
  if(isMadgraph){mcName = "madgraph"; cout << "This is a madgraph MC sample" << endl;}

  bool isPythia=false;
  size_t pos_pythia  = _inputFileName.find("pythia");
  if(pos_pythia!= std::string::npos)
    isPythia=true;
  if(isPythia)cout << "This is a pythia MC sample" << endl;

  cout << "There are " << nPtBins01 << " bins in 1st leading jet" << endl;
  cout << "There are " << nPtBins02 << " bins in 2nd leading jet" << endl;
  cout << "There are " << nPtBins03 << " bins in 3rd leading jet" << endl;
  cout << "There are " << nPtBins04 << " bins in 4th leading jet" << endl;
  cout << "There are " << nYBins << " bins." << endl;

  Long64_t nentries = fChain->GetEntriesFast();

  cout << _inputFileName << " has " << nentries << " entries." << endl;

  cout << "Running mode: " << endl;
  cout << "exclusive=" << exclusive 
       << "\t DEBUG=" << DEBUG << endl;
  cout << endl;
  cout << "The cuts applied: " << endl;
  
  cout << " minZPt= " << minZPt << endl;
  cout << " minJetPt= " << minJetPt << endl;
  cout << " maxJetEta= " << maxJetEta << endl;
  cout << " mindR= " << mindR << endl;
  cout << " minMll = " << minMll << endl;
  cout << " maxMll = " << maxMll << endl;

  cout << "studying " << leptonName << endl;
  if(abs(lepID)==11){
    cout << " minLepPt= " << minLepPt << endl;
    cout << " minEleBarrelEta = " << minEleBarrelEta << endl;
    cout << " maxEleBarrelEta = " << maxEleBarrelEta << endl;
    cout << " minEleEndcapEta = " << minEleEndcapEta << endl;
    cout << " maxEleEndcapEta = " << maxEleEndcapEta << endl;
  }
  else if(abs(lepID)==13 || abs(lepID)==0){
    cout << " minLepPt = " << minLepPt << endl;
    cout << " maxMuoEta = " << maxMuoEta << endl;
  }
  cout << "=========================================================" << endl;

  // dummy proof, in case someone put a negative number
  int leptonPID = abs(lepID); 

  //=============================================================================
  //   Book histograms
  //=============================================================================

  TH1D* h_mZ   = new TH1D("h_mZ","",200,20.0,220.0);
  h_mZ->SetXTitle("M_{ll} [GeV/c^{2}]");
  h_mZ->Sumw2();

  TH1D* h_nvtx = new TH1D("h_nvtx","",41.5,-0.5,40.5);
  h_nvtx->SetXTitle("Number of good vertices");
  h_nvtx->Sumw2();

  TH1D* h_njet = new TH1D("h_njet","",6,-0.5,5.5);
  h_njet->SetXTitle("#geq n jet");
  h_njet->Sumw2();

  TH2D* h2_ystarpstar = new TH2D("h2_ystarpstar","",60,-3.0,3.0,125,0,250);
  h2_ystarpstar->SetXTitle("0.5(Y_{Z}-Y_{jet})");
  h2_ystarpstar->SetYTitle("p_{T}cosh[0.5(Y_{Z}-Y_{jet})]");

  const int nMAXJETS=4;
   
  TH1D* h_mc_jetpt[nMAXJETS];
  TH1D* h_mc_jety[nMAXJETS];

  TH1D* h_predict_jetpt_template = new TH1D("h_predict_jetpt_template","",nPtBins01,fBinsPt01);
  h_predict_jetpt_template->Sumw2();
   

  h_mc_jetpt[0] = new TH1D(Form("h_mc_jetpt%02i",1),"1st leading jet",nPtBins01,fBinsPt01);
  h_mc_jetpt[1] = new TH1D(Form("h_mc_jetpt%02i",2),"2nd leading jet",nPtBins02,fBinsPt02);
  h_mc_jetpt[2] = new TH1D(Form("h_mc_jetpt%02i",3),"3rd leading jet",nPtBins03,fBinsPt03);
  h_mc_jetpt[3] = new TH1D(Form("h_mc_jetpt%02i",4),"4th leading jet",nPtBins04,fBinsPt04);

  TH1D* h_predict_jety_template = new TH1D("h_predict_jety_template","",nYBins,fBinsY);
  h_predict_jety_template->Sumw2();


  for(int ij=0; ij<nMAXJETS; ij++){

    h_mc_jetpt[ij]->SetXTitle("p_{T}(jet) [GeV]");
    h_mc_jetpt[ij]->Sumw2();

    h_mc_jety[ij] = (TH1D*)h_predict_jety_template->Clone(Form("h_mc_jety%02i",ij+1));
    h_mc_jety[ij]->SetXTitle("Rapidity(jet)");

  }


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


  int nPass[30];
   
  typedef map<double, int,  std::greater<double> > myMap;
  myMap sorted_genJetEtMap;
  typedef myMap::iterator mapIter;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    sorted_genJetEtMap.clear();
    //     if(jentry > 1000) break;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    nPass[0]++;

    double eventWeight = 1;

    if(mcWeight_>0)eventWeight *= mcWeight_;
      
    if(DEBUG==1){
      cout << "MCweight = " << mcWeight_ << endl;
      cout << "eventWeight = " << eventWeight << endl;
    }

    h_nvtx->Fill(EvtInfo_NumVtx, eventWeight);

    if(genLepId_->size()<2)continue;

    
    int lepPlusIndex = -1;
    int lepMinusIndex = -1;
    for(unsigned int igen=0; igen < genLepId_->size(); igen++){

      int PID       = genLepId_->at(igen);
      bool isLepPlus = (PID == (-leptonPID)) || 
	(leptonPID==0 && (PID== -11 || PID== -13));
      bool isLepMinus= (PID == ( leptonPID)) || 
	(leptonPID==0 && (PID==  11 || PID== 13));

      if(lepPlusIndex < 0 && isLepPlus)
	lepPlusIndex = igen;

      if(lepMinusIndex < 0 && isLepMinus)
	lepMinusIndex = igen;
	

      // need to add a line and see if this is a daughter of Z-> tau tau
      if(lepPlusIndex >= 0 && lepMinusIndex >=0)
	break;
	
    }
    
    // do not find m+ or lep-
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

    int nPt20=0;
    int nPt20Bare=0;
    int indexNumbers[2] = {lepPlusIndex, lepMinusIndex};

    TLorentzVector l4_lep[2]; 
    TLorentzVector l4_lepBare[2]; 
    

    for(int ip=0; ip < 2; ip++){

      l4_lep[ip].SetPtEtaPhiE(genLepPt_->at(indexNumbers[ip]),
			      genLepEta_->at(indexNumbers[ip]),
			      genLepPhi_->at(indexNumbers[ip]),
			      genLepE_->at(indexNumbers[ip])
			      );


      l4_lepBare[ip] = l4_lep[ip];
      
      double ptBare = l4_lepBare[ip].Pt();
      double etaBare = l4_lepBare[ip].Eta();
      if((leptonPID==13 || leptonPID==0) && ptBare > minLepPt && fabs(etaBare) < maxMuoEta)nPt20Bare++;
      if(leptonPID==11 && ptBare > minLepPt && ( 
					    (fabs(etaBare) > minEleEndcapEta && 
					     fabs(etaBare) < maxEleEndcapEta) ||

					    (fabs(etaBare) > minEleBarrelEta && 
					     fabs(etaBare) < maxEleBarrelEta
					     )
					    ))nPt20Bare++;


      TLorentzVector l4_pho(0,0,0,0);
      
      if( genLepPhoE_->at(indexNumbers[ip]) > 1e-5)
	l4_pho.SetPtEtaPhiE(genLepPhoPt_->at(indexNumbers[ip]),
			    genLepPhoEta_->at(indexNumbers[ip]),
			    genLepPhoPhi_->at(indexNumbers[ip]),
			    genLepPhoE_->at(indexNumbers[ip])
			    );


      if(mode>0)
	l4_lep[ip] += l4_pho;
      double pt = l4_lep[ip].Pt();
      double eta = l4_lep[ip].Eta();

      // muon
      if( (leptonPID==13 || leptonPID==0) && pt > minLepPt 
	  && fabs(eta) < maxMuoEta)nPt20++;
      if(leptonPID==11 && pt > minLepPt && ( 
					    (fabs(eta) > minEleEndcapEta && 
					     fabs(eta) < maxEleEndcapEta) ||

					    (fabs(eta) > minEleBarrelEta && 
					     fabs(eta) < maxEleBarrelEta
					     )
					    ))nPt20++;


    }

    if(DEBUG==1)
      cout << "nPt20 = " << nPt20 << endl;

    if(nPt20 < 2)continue;
    nPass[2]++;
    

    if(mode==2 && nPt20Bare < 2)continue;
    nPass[3]++;

    TLorentzVector l4_z = l4_lep[0]+l4_lep[1];
    double mll = l4_z.M();

    TLorentzVector l4_zBare = l4_lepBare[0]+l4_lepBare[1];
    double mllBare = l4_zBare.M();
    
    if(mode<2)
      h_mZ->Fill(mll,eventWeight);
    else
      {
	int binIndex0 = h_mZ->FindBin(mll);
	int binIndex1 = h_mZ->FindBin(mllBare);
	if(binIndex0 == binIndex1)
	  h_mZ->Fill(mll,eventWeight);

      }

    if(mll < minMll || mll > maxMll)continue;
    nPass[4]++;

    if(mode==2 && (mllBare < minMll || mllBare > maxMll))continue;
    nPass[5]++;

    if(DEBUG==1)
      cout << "dilepton mass = " << mll << endl;
      
    double ptz = l4_z.Pt();
    if(ptz < minZPt)continue;
    nPass[6]++;

    double ptzBare = l4_zBare.Pt();
    if(mode==2 && ptzBare < minZPt)continue;
    nPass[6]++;

    // now look for jets
    double maxGenJetPt = -9999;
    int maxGenJetIndex = -1;
    unsigned int nGenJets = 0;


    for(unsigned int ij = 0; ij < genJetPt_->size(); ij ++){

      double thisGenJetPt  = genJetPt_->at(ij);
      double thisGenJetEta = genJetEta_->at(ij);

      if(thisGenJetPt < minJetPt)continue;
      if(fabs(thisGenJetEta) > maxJetEta)continue;

      TLorentzVector thisGenJet_l4(0,0,0,0);
	
      thisGenJet_l4.SetPtEtaPhiE(genJetPt_->at(ij),
				 genJetEta_->at(ij),
				 genJetPhi_->at(ij),
				 genJetE_->at(ij));
	
      double dr_ep = l4_lep[0].DeltaR(thisGenJet_l4);
      double dr_em = l4_lep[1].DeltaR(thisGenJet_l4);

      if(dr_ep < mindR)continue;
      if(dr_em < mindR)continue;

      double dr_epBare = l4_lepBare[0].DeltaR(thisGenJet_l4);
      double dr_emBare = l4_lepBare[1].DeltaR(thisGenJet_l4);

      if(mode==2 && dr_epBare < mindR)continue;
      if(mode==2 && dr_emBare < mindR)continue;

      nGenJets++; 
      sorted_genJetEtMap.insert(std::pair<double, int>(thisGenJetPt,ij));  

      if(thisGenJetPt > maxGenJetPt)
	{
	  maxGenJetPt = thisGenJetPt;
	  maxGenJetIndex = ij;
	}

	
    } // end of loop over jets
      
    if(nGenJets!= sorted_genJetEtMap.size())
      { 
	cout << "errors in map size" << endl;
	continue;
      }

    if(DEBUG==1)
      cout << "max gen jet index = " << maxGenJetIndex << endl;

    for(unsigned int ik=0; ik<h_njet->GetNbinsX(); ik++)
      if(nGenJets>= ik)
	h_njet->Fill(ik,eventWeight);

    if(maxGenJetIndex < 0)continue;

    nPass[7]++;

    if(exclusive && nGenJets!=1)continue;
     
    nPass[8]++;

 
    TLorentzVector l4_j(0,0,0,0);
    l4_j.SetPtEtaPhiE(genJetPt_->at(maxGenJetIndex),
		      genJetEta_->at(maxGenJetIndex),
		      genJetPhi_->at(maxGenJetIndex),
		      genJetE_->at(maxGenJetIndex));

    double ptjet = l4_j.Pt();

    double yz = l4_z.Rapidity();
    double yj = l4_j.Rapidity();
    double yB = 0.5*(yz + yj);
    double ystar = 0.5*(yz-yj);

    double yzBare = l4_zBare.Rapidity();
    double yBBare = 0.5*(yzBare + yj);
    double ystarBare = 0.5*(yzBare-yj);

     
      
    if(DEBUG==1)
      cout << "Now ordering jets" << endl;
    int countGenJet=0;
    for (mapIter it_part= sorted_genJetEtMap.begin();
	 it_part != sorted_genJetEtMap.end() && countGenJet< nMAXJETS; ++it_part)
      {
        double pt_mapthis = genJetPt_->at(it_part->second);

	int jet_index = it_part->second;

	TLorentzVector l4_jthis(0,0,0,0);
	l4_jthis.SetPtEtaPhiE(genJetPt_->at(jet_index),
			      genJetEta_->at(jet_index),
			      genJetPhi_->at(jet_index),
			      genJetE_->at(jet_index));

	double y_mapthis  = l4_jthis.Rapidity();
	
	if(DEBUG==1)
	  cout << "pt = " << pt_mapthis << "\t y = " << y_mapthis << endl;

        h_mc_jetpt[countGenJet]->Fill(pt_mapthis,eventWeight);
        h_mc_jety[countGenJet]->Fill(fabs(y_mapthis),eventWeight);

        countGenJet++;
      }

    if(DEBUG==1)cout << "nGenJets = " << sorted_genJetEtMap.size() 
		     << "\t countGenJet = " << countGenJet << endl;
    if(mode<2){
      h_zpt->Fill(ptz,eventWeight);
      h_zy->Fill(fabs(yz),eventWeight);
      h_yB->Fill(fabs(yB),eventWeight);
      h_ystar->Fill(fabs(ystar),eventWeight);
    }

    else
      {
	int binIndex0 = h_zpt->FindBin(ptz);
	int binIndex1 = h_zpt->FindBin(ptzBare);
	if(binIndex0 == binIndex1)
	  h_zpt->Fill(ptz,eventWeight);

	binIndex0 = h_zy->FindBin(fabs(yz));
	binIndex1 = h_zy->FindBin(fabs(yzBare));
	if(binIndex0 == binIndex1)
	  h_zy->Fill(fabs(yz),eventWeight);

	binIndex0 = h_yB->FindBin(fabs(yB));
	binIndex1 = h_yB->FindBin(fabs(yBBare));
	if(binIndex0 == binIndex1)
	  h_yB->Fill(fabs(yB),eventWeight);


	binIndex0 = h_ystar->FindBin(fabs(ystar));
	binIndex1 = h_ystar->FindBin(fabs(ystarBare));
	if(binIndex0 == binIndex1)
	  h_ystar->Fill(fabs(ystar),eventWeight);

      }


    
    h_jetpt->Fill(ptjet,eventWeight); 
    h_jety->Fill(fabs(yj),eventWeight);

    h2_ystarpstar->Fill(ystar,ptjet*TMath::CosH(ystar),eventWeight);
 

    if(DEBUG==1){
      double dR1 = l4_lep[0].DeltaR(l4_j);
      double dR2 = l4_lep[1].DeltaR(l4_j);
      cout << "dR1 = " << dR1 << "\t dR2=" << dR2 << endl;
    }
  } // end of loop over entries


  std::string prefix;
  switch (mode)
    {
    case 0:
      prefix = "bare";
      break;
    case 1:
      prefix = "dressed";
      break;
    case 2:
      prefix = "both";
      break;
    default:
      prefix = "test";
      break;
    } // end of switch

  if(exclusive)prefix += "_exclusive1Jet";
  else prefix += "_inclusive";
  if(minZPt>1e-6)prefix += Form("_zPt%d",(int)minZPt);

  std::string remword  ="/data2/syu/zjet_vectorNtuple/";

  size_t pos  = _inputFileName.find(remword);

  if(pos!= std::string::npos)
    _inputFileName.swap(_inputFileName.erase(pos,remword.length()));


  TFile* outFile = new TFile(Form("%s_%s_%s",prefix.data(),
				  leptonName.data(),
				  _inputFileName.data()),"recreate");       

  h2_ystarpstar->Write();
  h_mZ->Write();
  h_nvtx->Write();
  h_njet->Write();
        
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

  outFile->Close();

}
