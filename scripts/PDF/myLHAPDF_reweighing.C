#define myLHAPDF_reweighing_cxx
#include "myLHAPDF_reweighing.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

const double minVPt  =40.0;
const double maxVEta = 1.45; // for photon

const double minJetPt=30.0;
const double maxJetEta=2.4;
const double mindR=0.5;

const double minLepPt = 20.0;
const double maxLepEta = 2.1;
const double minM = 76.0;
const double maxM =106.0;



void myLHAPDF_reweighing::Loop(int lepID)
{
  //===========================================================================
  //   Print out information
  //===========================================================================

  if (fChain == 0) return;
  cout << "=========================================================" << endl;

  Long64_t nentries = fChain->GetEntries();

  cout << "The tree has " << nentries << " entries." << endl;
  std::string leptonName;
  if(abs(lepID)==13)leptonName = "muon";
  else if(abs(lepID)==11)leptonName = "electron";
  else if(abs(lepID)==22)leptonName = "photon";
  else if(abs(lepID)== 0)leptonName = "all";

  cout << "Running mode: " << endl;
  cout << "The cuts applied: " << endl;
  
  cout << " minVPt= " << minVPt << endl;
  cout << " minJetPt= " << minJetPt << endl;
  cout << " maxJetEta= " << maxJetEta << endl;
  cout << " mindR= " << mindR << endl;

  cout << "studying " << leptonName << endl;
  if(abs(lepID)==0 || abs(lepID)==11 || abs(lepID)==13){
    cout << " minLepPt= " << minLepPt << endl;
    cout << " maxLepEta = " << maxLepEta << endl;
    cout << " minM = " << minM << endl;
    cout << " maxM = " << maxM << endl;
  }
  else if(abs(lepID)==22){
    cout << " maxVEta = " << maxVEta << endl;
  }
  cout << "=========================================================" << endl;

  // dummy proof, in case someone put a negative number
  int leptonPID = abs(lepID); 

  //===========================================================================
  //   Book histograms
  //===========================================================================

  TH1D* h_mZ   = new TH1D("h_mZ","",200,20.0,220.0);
  h_mZ->SetXTitle("M_{ll} [GeV/c^{2}");
  h_mZ->Sumw2();

  TH1D* h_mZ_parton   = new TH1D("h_mZ_parton","",200,20.0,220.0);
  h_mZ_parton->SetXTitle("M_{ll} [GeV/c^{2}");
  h_mZ_parton->Sumw2();

  TH1D* h_nvtx = new TH1D("h_nvtx","",41.5,-0.5,40.5);
  h_nvtx->SetXTitle("Number of good vertices");
  h_nvtx->Sumw2();

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

  const int NPDFs=3;

  std::string pdfname[NPDFs]={"cteq6l1", "cteq66m", "mstw2008"};

  TH1D* h_zpt_pdf[NPDFs];
  TH1D* h_jetpt_pdf[NPDFs];

  TH1D* h_zy_pdf[NPDFs];
  TH1D* h_jety_pdf[NPDFs];
  TH1D* h_ystar_pdf[NPDFs];
  TH1D* h_yB_pdf[NPDFs];


  for(int i=0; i < NPDFs; i++){

    h_zpt_pdf[i] = (TH1D*)
      h_zpt_template->Clone(Form("h_zpt_%s",pdfname[i].data()));
    h_jetpt_pdf[i] = (TH1D*)
      h_jetpt_template->Clone(Form("h_jetpt_%s",pdfname[i].data()));
    h_zy_pdf[i] = (TH1D*)
      h_y_template->Clone(Form("h_zy_%s",pdfname[i].data()));
    h_jety_pdf[i] = (TH1D*)
      h_y_template->Clone(Form("h_jety_%s",pdfname[i].data()));
    h_ystar_pdf[i] = (TH1D*)
      h_y_template->Clone(Form("h_ystar_%s",pdfname[i].data()));
    h_yB_pdf[i] = (TH1D*)
      h_y_template->Clone(Form("h_yB_%s",pdfname[i].data()));
  }

  int nPass[30]={0};

  MyPDF* cteq66R = new MyPDF("cteq66.LHgrid",1);
  MyPDF* mstw2008R = new MyPDF("MSTW2008nlo68cl.LHgrid",2);

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //    if(jentry > 10000) break;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    nPass[0]++;

    int lepPlusIndex = -1;
    int lepMinusIndex = -1;
    double eventWeight = 1;

    h_nvtx->Fill(EvtInfo_NumVtx, eventWeight);


    for(unsigned int igen=0; igen < genParId_->size(); igen++){

      int PID       = genParId_->at(igen);
      int motherPID = genMomParId_->at(igen);
      int status    = genParSt_->at(igen);
      bool isFromZAfterFSR  = (motherPID == PID) && (status == 1);

      
      bool isLepPlus = ((PID == (-leptonPID)) || 
			(leptonPID==0 && (PID== -11 || PID== -13)))  
	&& isFromZAfterFSR;

      bool isLepMinus= ((PID == ( leptonPID)) || 
			(leptonPID==0 && (PID==  11 || PID== 13)))
	&& isFromZAfterFSR;


      TLorentzVector temp_l4(0,0,0,0);
      temp_l4.SetPtEtaPhiE(
			   genParPt_->at(igen),
			   genParEta_->at(igen),
			   genParPhi_->at(igen),
			   genParE_->at(igen)
			   );

      if(PID == 23 && status == 3)
	h_mZ_parton->Fill(temp_l4.M());

      if(lepPlusIndex < 0 && isLepPlus)
	lepPlusIndex = igen;

      if(lepMinusIndex < 0 && isLepMinus)
	lepMinusIndex = igen;
	
      if(leptonPID< 22 && lepPlusIndex >= 0 && lepMinusIndex >=0)
	break;
      if(leptonPID== 22 && lepMinusIndex >=0)
	break;
	
    }

    // do not find m+ or lep-
    if(leptonPID < 22 && (lepPlusIndex < 0 || lepMinusIndex < 0))continue;
    if(leptonPID == 22 && lepMinusIndex < 0)continue;

    nPass[1]++;

    int nPt20=0;

    int indexNumbers[2] = {lepPlusIndex, lepMinusIndex};

    if(leptonPID < 22){
      for(int ip=0; ip < 2; ip++){
	
	double pt  = genParPt_->at(indexNumbers[ip]);
	double eta = genParEta_->at(indexNumbers[ip]);
	
	if(pt > minLepPt && fabs(eta) < maxLepEta)nPt20++;
      
      }

    if(nPt20 < 2)continue;
    }

    nPass[2]++;

    TLorentzVector l4_lepp(0,0,0,0);

    if(leptonPID<22){
      l4_lepp.SetPtEtaPhiE(genParPt_->at(lepPlusIndex),
			   genParEta_->at(lepPlusIndex),
			   genParPhi_->at(lepPlusIndex),
			   genParE_->at(lepPlusIndex));
    }

    TLorentzVector l4_lepm(0,0,0,0);

    l4_lepm.SetPtEtaPhiE(genParPt_->at(lepMinusIndex),
			 genParEta_->at(lepMinusIndex),
			 genParPhi_->at(lepMinusIndex),
			 genParE_->at(lepMinusIndex));



    TLorentzVector l4_z = leptonPID< 22? l4_lepp + l4_lepm : l4_lepm;
//     cout << EvtInfo_EventNum << "\t vector boson " << l4_z.Pt() 
// 	 << "\t" << l4_z.Eta() << "\t" << 
//       l4_z.Phi() << "\t" << l4_z.E() << endl;

    double mll = l4_z.M();

    h_mZ->Fill(mll,eventWeight);

    if(leptonPID < 22 && (mll < minM || mll > maxM))continue;

    nPass[3]++;

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
	
      double dr_ep = leptonPID < 22? l4_lepp.DeltaR(thisGenJet_l4) : 2.0;
      double dr_em = l4_lepm.DeltaR(thisGenJet_l4);

      if(dr_ep < mindR)continue;
      if(dr_em < mindR)continue;
      
      nGenJets++; 

      if(thisGenJetPt > maxGenJetPt)
	{
	  maxGenJetPt = thisGenJetPt;
	  maxGenJetIndex = ij;
	}
	
    } // end of loop over jets
      
    if(maxGenJetIndex < 0)continue;

    nPass[4]++;

    if(nGenJets!=1)continue;
     
    nPass[5]++;


    double ptz = l4_z.Pt();
    double etaz = l4_z.Eta();
    if(ptz < minVPt)continue;
    nPass[6]++;

    if(leptonPID == 22 && fabs(etaz) > maxVEta)continue;
    nPass[7]++;
 
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


    double weight_pdf[3]={
      1.0,
      cteq66R->weight(pdfInfo_,0),
      mstw2008R->weight(pdfInfo_,0)
      //      1.0,
      //      1.0
    };

    //    cout << "weight = " << weight_pdf[1] << "\t" << 
    //       cteq66PDFw_->at(0) << endl;
    //    cout << "weight = " << weight_pdf[2] << "\t" << 
    //      mstw2008PDFw_->at(0) << endl;
    h_zpt->Fill(ptz,eventWeight);
    h_jetpt->Fill(ptjet,eventWeight); 
    h_zy->Fill(fabs(yz),eventWeight);
    h_jety->Fill(fabs(yj),eventWeight);
    h_yB->Fill(fabs(yB),eventWeight);
    h_ystar->Fill(fabs(ystar),eventWeight);

    for(int i=0; i< NPDFs; i++){
      
      h_zpt_pdf[i]  ->Fill(ptz,weight_pdf[i]);
      h_jetpt_pdf[i]->Fill(ptjet,weight_pdf[i]); 
      h_zy_pdf[i]   ->Fill(fabs(yz),weight_pdf[i]);
      h_jety_pdf[i] ->Fill(fabs(yj),weight_pdf[i]);
      h_yB_pdf[i]   ->Fill(fabs(yB),weight_pdf[i]);
      h_ystar_pdf[i]->Fill(fabs(ystar),weight_pdf[i]);

    }


  } // end of loop over entries
  delete cteq66R;
  delete mstw2008R;
  for(int i=0;i<30;i++)if(nPass[i]>0)cout << "nPass[" << i << "] = " << nPass[i] << endl;

  std::string remword  ="root://eoscms//eos/cms/store/user/syu/PDF_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/zjets_mc_";

  size_t pos  = _inputFileName.find(remword);

  if(pos!= std::string::npos)
    _inputFileName.swap(_inputFileName.erase(pos,remword.length()));

  TFile* outFile = new TFile(Form("DY_PDFstudy_%s",_inputFileName.data()),"recreate");       

  h_mZ->Write();
  h_mZ_parton->Write();
  h_nvtx->Write();
        
  h_yB->Write();
  h_ystar->Write();
  
  h_zy->Write();
  h_zpt->Write();

  h_jety->Write();   
  h_jetpt->Write();

  for(int i=0;i<NPDFs;i++)
    {
      h_zpt_pdf[i]->Write();
      h_jetpt_pdf[i]->Write();

      h_zy_pdf[i]->Write();
      h_jety_pdf[i]->Write();   
      h_ystar_pdf[i]->Write();
      h_yB_pdf[i]->Write();

    }

  outFile->Close();
}
