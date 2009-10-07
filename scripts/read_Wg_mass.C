#include <iostream>
#include <fstream>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>


double ClusterMass(TLorentzVector& gamma, TLorentzVector& electron,
		   TLorentzVector& neutrino)
{
  double PTMISS = neutrino.Pt();
  double PT2AB  = pow((gamma+electron).Pt(),2);
  double MASSAB = (gamma+electron).M();
  double CTM2   = pow(sqrt(PT2AB + MASSAB*MASSAB) + PTMISS,2)
    - pow((gamma+electron+neutrino).Pt(),2);
  double CTM    = sqrt(CTM2);
  
  return CTM;
    
}


using namespace std;
// here, only leading photon is required
void read_Wg_mass(std::string inputfile, double totalxsec=1.890E-08, int baur=0)
//void read_Wg_mass(std::string inputfile, double totalxsec=1.069E-07, int baur=0)
{
  // Booking histograms
  int ntotal=0;
  int nsubtotal=0;
  int npass=0;
  int netbin=150;
  float etmin=0.0;
  float etmax=300.0;
  int netabin=50;
  float etamin=-5.0;
  float etamax=5.0;
  
  int ndRbin=100;
  float dRmin=0;
  float dRmax=10.0;

  char name[300];
  sprintf(name,"%s%s",inputfile.data(),"_masshisto.root");  
  TFile* outFile = new TFile(name,"recreate");

  TH2F* hFSR_dRgpt = new TH2F("hFSR_dRgpt","#Delta R (lepton,#gamma) vs. "
			      " E_{T}(#gamma)",netbin,etmin,etmax,ndRbin,dRmin,dRmax);
  TH2F* hISR_dRgpt = new TH2F("hISR_dRgpt","#Delta R (lepton,#gamma) vs. "
			      " E_{T}(#gamma)",netbin,etmin,etmax,ndRbin,dRmin,dRmax);
  TH2F* hFSR_lptgpt = new TH2F("hFSR_lptgpt","E_{T}(e) vs. "
			       " E_{T}(#gamma)",netbin,etmin,etmax,netbin,etmin,etmax);
  TH2F* hISR_lptgpt = new TH2F("hISR_lptgpt","E_{T}(e) vs. "
			       " E_{T}(#gamma)",netbin,etmin,etmax,netbin,etmin,etmax);


  TH1F* h_dR = new TH1F("h_dR","#Delta R (lepton,#gamma)",ndRbin,dRmin,dRmax);
  TH1F* h_mW = new TH1F("h_mW","Mass of W boson",netbin,etmin,etmax);  
  TH1F* h_mWg = new TH1F("h_mWg","Mass of W boson + photon",netbin,etmin,etmax);
  TH2F* h_mW3 = new TH2F("h_mW3","M(W#gamma) vs M(W)",netbin,etmin,etmax,netbin,etmin,etmax);
  TH1F* h_dR2 = new TH1F("h_dR2","#Delta R (lepton,#gamma)",ndRbin,dRmin,dRmax);
  TH1F* h_mW2 = new TH1F("h_mW2","Mass of W boson",netbin,etmin,etmax);  
  TH1F* h_mWg2 = new TH1F("h_mWg2","Mass of W boson + photon",netbin,etmin,etmax);
  TH2F* h_mW32 = new TH2F("h_mW32","M(W#gamma) vs M(W)",netbin,etmin,etmax,netbin,etmin,etmax);

  TH1F* hFSR_dR = new TH1F("hFSR_dR","#Delta R (lepton,#gamma)",ndRbin,dRmin,dRmax);
  TH1F* hFSR_mW = new TH1F("hFSR_mW","Mass of W boson",netbin,etmin,etmax);  
  TH1F* hFSR_mWg = new TH1F("hFSR_mWg","Mass of W boson + photon",netbin,etmin,etmax);
  TH2F* hFSR_mW3 = new TH2F("hFSR_mW3","M(W#gamma) vs M(W)",netbin,etmin,etmax,netbin,etmin,etmax);
  TH1F* hFSR_dR2 = new TH1F("hFSR_dR2","#Delta R (lepton,#gamma)",ndRbin,dRmin,dRmax);
  TH1F* hFSR_mW2 = new TH1F("hFSR_mW2","Mass of W boson",netbin,etmin,etmax);  
  TH1F* hFSR_mWg2 = new TH1F("hFSR_mWg2","Mass of W boson + photon",netbin,etmin,etmax);
  TH2F* hFSR_mW32 = new TH2F("hFSR_mW32","M(W#gamma) vs M(W)",netbin,etmin,etmax,netbin,etmin,etmax);

  TH1F* hISR_dR = new TH1F("hISR_dR","#Delta R (lepton,#gamma)",ndRbin,dRmin,dRmax);
  TH1F* hISR_mW = new TH1F("hISR_mW","Mass of W boson",netbin,etmin,etmax);  
  TH1F* hISR_mWg = new TH1F("hISR_mWg","Mass of W boson + photon",netbin,etmin,etmax);
  TH2F* hISR_mW3 = new TH2F("hISR_mW3","M(W#gamma) vs M(W)",netbin,etmin,etmax,netbin,etmin,etmax);
  TH1F* hISR_dR2 = new TH1F("hISR_dR2","#Delta R (lepton,#gamma)",ndRbin,dRmin,dRmax);
  TH1F* hISR_mW2 = new TH1F("hISR_mW2","Mass of W boson",netbin,etmin,etmax);  
  TH1F* hISR_mWg2 = new TH1F("hISR_mWg2","Mass of W boson + photon",netbin,etmin,etmax);
  TH2F* hISR_mW32 = new TH2F("hISR_mW32","M(W#gamma) vs M(W)",netbin,etmin,etmax,netbin,etmin,etmax);


  TH1F* h_gptb = new TH1F("h_gptb","p_{T}(#gamma) before cuts",netbin,etmin,etmax);
  TH1F* h_gpta = new TH1F("h_gpta","p_{T}(#gamma) after cuts",netbin,etmin,etmax);
  TH1F* h_lptb = new TH1F("h_lptb","p_{T}(e) before cuts",netbin,etmin,etmax);
  TH1F* h_lpta = new TH1F("h_lpta","p_{T}(e) after cuts",netbin,etmin,etmax);
  TH1F* h_nptb = new TH1F("h_nptb","p_{T}(#nu) before cuts",netbin,etmin,etmax);
  TH1F* h_npta = new TH1F("h_npta","p_{T}(#nu) after cuts",netbin,etmin,etmax);

  TH1F* h_getab = new TH1F("h_getab","#eta(#gamma) before cuts",netabin,etamin,etamax);
  TH1F* h_getaa = new TH1F("h_getaa","#eta(#gamma) after cuts",netabin,etamin,etamax);
  TH1F* h_letab = new TH1F("h_letab","#eta(e) before cuts",netabin,etamin,etamax);
  TH1F* h_letaa = new TH1F("h_letaa","#eta(e) after cuts",netabin,etamin,etamax);
  TH1F* h_netab = new TH1F("h_netab","#eta(#nu) before cuts",netabin,etamin,etamax);
  TH1F* h_netaa = new TH1F("h_netaa","#eta(#nu) after cuts",netabin,etamin,etamax);

  TH1F* hISR_gptb = new TH1F("hISR_gptb","ISR p_{T}(#gamma) before cuts",
			     netbin,etmin,etmax);
  TH1F* hFSR_gptb = new TH1F("hFSR_gptb","FSR p_{T}(#gamma) before cuts",
			     netbin,etmin,etmax);
  TH1F* hISR_gpta = new TH1F("hISR_gpta","ISR p_{T}(#gamma) after cuts",
			     netbin,etmin,etmax);
  TH1F* hFSR_gpta = new TH1F("hFSR_gpta","FSR p_{T}(#gamma) after cuts",
			     netbin,etmin,etmax);
  

  TH1F* hmT_FSRb = new TH1F("hmT_FSRb","M_{T}(W) before cuts", netbin, etmin, etmax);
  TH1F* hmT_FSRa = new TH1F("hmT_FSRa","M_{T}(W) after cuts", netbin, etmin, etmax);

  TH1F* hmT_ISRb = new TH1F("hmT_ISRb","M_{T}(W) before cuts", netbin, etmin, etmax);
  TH1F* hmT_ISRa = new TH1F("hmT_ISRa","M_{T}(W) after cuts", netbin, etmin, etmax);


  TH1F* hmT3 = new TH1F("hmT3","M_{T}(W#gamma)", netbin, etmin, etmax);
 

  TH1F* hmT3_FSRb = new TH1F("hmT3_FSRb","M_{T}(W#gamma) before cuts", netbin, etmin, etmax);
  TH1F* hmT3_FSRa = new TH1F("hmT3_FSRa","M_{T}(W#gamma) after cuts", netbin, etmin, etmax);

  TH1F* hmT3_ISRb = new TH1F("hmT3_ISRb","M_{T}(W#gamma) before cuts", netbin, etmin, etmax);
  TH1F* hmT3_ISRa = new TH1F("hmT3_ISRa","M_{T}(W#gamma) after cuts", netbin, etmin, etmax);



  TH2F* hmT3mT_FSRb = new TH2F("hmT3mT_FSRb","M_{T}(W#gamma) vs M_{T}(W) before cuts", netbin, etmin, etmax, netbin, etmin, etmax);
  TH2F* hmT3mT_FSRa = new TH2F("hmT3mT_FSRa","M_{T}(W#gamma) vs M_{T}(W) after cuts", netbin, etmin, etmax, netbin, etmin, etmax);

  TH2F* hmT3mT_ISRb = new TH2F("hmT3mT_ISRb","M_{T}(W#gamma) vs M_{T}(W) before cuts", netbin, etmin, etmax, netbin, etmin, etmax);
  TH2F* hmT3mT_ISRa = new TH2F("hmT3mT_ISRa","M_{T}(W#gamma) vs M_{T}(W) after cuts", netbin, etmin, etmax, netbin, etmin, etmax);




  // Looping over entries

  ifstream fin;
  sprintf(name,"%s%s",inputfile.data(),".dat");
  fin.open(name);
  cout << "Opening " << name << endl;
  
  const int NCOL=baur == 1? 17: 16; // photon, charged lepton, neutral lepton
  double cl[NCOL];

  for(int j=0;j<NCOL;j++)
    fin >> cl[j];

  while(!fin.eof()){
    ntotal ++;
    if(ntotal%100==0)cout << "ntotal = " << ntotal << endl;
    TLorentzVector gamma1(cl[0],cl[1],cl[2],cl[3]);
    TLorentzVector clepton(cl[4],cl[5],cl[6],cl[7]);
    TLorentzVector nlepton(cl[8],cl[9],cl[10],cl[11]);
    TLorentzVector gamma2(cl[12],cl[13],cl[14],cl[15]);

    double weight = 1.0;


    TLorentzVector gamma_lead = gamma1.Et() > gamma2.Et()? 
      gamma1: gamma2;

    double mT3 = ClusterMass(gamma_lead,clepton,nlepton);

    bool hasFSR= mT3 < 90.0;
    bool hasISR= mT3 > 90.0;

    double mW = (clepton+nlepton).M();
    double mWg= (clepton+nlepton+gamma_lead).M();
    double dR = gamma_lead.DeltaR(clepton);
    double gpt  =gamma_lead.Pt();
    double geta =gamma_lead.Eta();
    

    TLorentzVector met_lv(nlepton.Px(),                    
			  nlepton.Py(),
			  0,
			  nlepton.Pt());


    TLorentzVector lep_lv(clepton.Px(),
			  clepton.Py(),
			  0,
			  clepton.Pt());

    TLorentzVector g_lv(gamma_lead.Px(),
			gamma_lead.Py(),
			0,
			gamma_lead.Pt());


    double mT = (met_lv+lep_lv).M();

    double lpt  =clepton.Pt();
    double leta =clepton.Eta();
    double npt  =nlepton.Pt();
    double neta =nlepton.Eta();

    double gy = gamma_lead.Rapidity();
    double ly = clepton.Rapidity();
    double ny = nlepton.Rapidity();

    for(int j=0;j<NCOL;j++)
      fin >> cl[j];

    
    if((!hasFSR && !hasISR) || clepton.E() <1e-6 || nlepton.E()<1e-6)
      continue;

    
    nsubtotal++;
    
    if(hasISR){
      h_mW->Fill(mW, weight);
      h_mWg->Fill(mWg, weight);
      h_mW3->Fill(mW,mWg, weight);
      h_dR->Fill(dR, weight);

      hISR_dRgpt->Fill(gpt,dR);
      hISR_lptgpt->Fill(gpt,lpt);

      hISR_mW->Fill(mW, weight);
      hISR_mWg->Fill(mWg, weight);
      hISR_mW3->Fill(mW,mWg, weight);
      hISR_dR->Fill(dR, weight);

      h_gptb->Fill(gpt, weight);
      h_lptb->Fill(lpt, weight);
      h_nptb->Fill(npt, weight);

      h_getab->Fill(gy, weight);
      h_letab->Fill(ly, weight);
      h_netab->Fill(neta, weight);

      hISR_gptb->Fill(gpt, weight);

      hmT_ISRb->Fill(mT, weight);
      hmT3_ISRb->Fill(mT3, weight);
      hmT3mT_ISRb->Fill(mT,mT3, weight);
      
    }

    if(hasFSR)
      {

	h_mW->Fill(mW, weight);
	h_mWg->Fill(mWg, weight);
	h_mW3->Fill(mW,mWg, weight);
	h_dR->Fill(dR, weight);

	hFSR_mW->Fill(mW, weight);
	hFSR_mWg->Fill(mWg, weight);
	hFSR_mW3->Fill(mW,mWg, weight);
	hFSR_dR->Fill(dR, weight);

	hFSR_dRgpt->Fill(gpt,dR,weight);
	hFSR_lptgpt->Fill(gpt,lpt,weight);


	h_gptb->Fill(gpt, weight);
	h_lptb->Fill(lpt, weight);
	h_nptb->Fill(npt, weight);

	h_getab->Fill(gy, weight);
	h_letab->Fill(ly, weight);
	h_netab->Fill(neta, weight);

	hFSR_gptb->Fill(gpt, weight);

	hmT_FSRb->Fill(mT, weight);
	hmT3_FSRb->Fill(mT3, weight);
	hmT3mT_FSRb->Fill(mT,mT3, weight);

      }
          

//     // muon channel LSP's cuts
//     bool photonCut_ISR = gpt>10.0 && fabs(geta)<2.7;
//     bool photonCut_FSR = gpt_FSR > 10.0 && fabs(geta_FSR)<2.7;
// //     if(!photonCut_ISR && !photonCut_FSR)continue;
//     if(!photonCut_ISR)continue;
//     if(lpt<5.0)continue;
//     if(npt<20.0)continue;
//     if(fabs(leta)>2.7)continue;

    // Poter's cuts, supposedly
    
//     bool photonCut_ISR = gpt>20.0 && fabs(geta)<2.7 && dR>1.1;
//     bool photonCut_FSR = gpt_FSR > 20.0 && fabs(geta_FSR)<2.7 && dR_FSR>1.1;
// //      if(!photonCut_ISR && !photonCut_FSR)continue;
//     if(!photonCut_ISR)continue;
//     if(lpt<15.0)continue;
//     if(fabs(leta)>2.7)continue;
//     if(mT < 60)continue;
//     if(npt<20.0)continue;


    // ATLAS Zhijun's cuts

//      bool photonCut_ISR = gpt>10.0 && fabs(geta)<2.5 && dR>0.7;
//      if(!photonCut_ISR)continue;
//      if(lpt<10.0)continue;
//      if(fabs(leta)>2.5)continue;


    // CDF cuts
    bool photonCut_ISR = gpt> 7.0 && dR>0.7;
    if(!photonCut_ISR)continue;

//      bool photonCut_ISR = gpt> 5.0 && dR>0.1 && mT3 > 90.;
//      if(!photonCut_ISR)continue;


    npass++;

    hmT3->Fill(mT3);
    
    if(hasISR){
      h_gpta->Fill(gpt, weight);
      h_lpta->Fill(lpt, weight);
      h_npta->Fill(npt, weight);

      h_getaa->Fill(gy, weight);
      h_letaa->Fill(ly, weight);
      h_netaa->Fill(neta, weight);
    

      h_mW2->Fill(mW, weight);
      h_mWg2->Fill(mWg, weight);
      h_mW32->Fill(mW,mWg, weight);
      h_dR2->Fill(dR, weight);

      
      hISR_mW2->Fill(mW, weight);
      hISR_mWg2->Fill(mWg, weight);
      hISR_mW32->Fill(mW,mWg, weight);
      hISR_dR2->Fill(dR, weight);
      hISR_gpta->Fill(gpt, weight);
      hmT_ISRa->Fill(mT, weight);
      hmT3_ISRa->Fill(mT3, weight);
      hmT3mT_ISRa->Fill(mT,mT3, weight);
    

    }

   if(hasFSR)
      {
	h_mW2->Fill(mW, weight);
	h_mWg2->Fill(mWg, weight);
	h_mW32->Fill(mW,mWg, weight);
	h_dR2->Fill(dR, weight);

	h_gpta->Fill(gpt, weight);
	h_lpta->Fill(lpt, weight);
	h_npta->Fill(npt, weight);

	h_getaa->Fill(gy, weight);
	h_letaa->Fill(ly, weight);
	h_netaa->Fill(neta, weight);

	hFSR_mW2->Fill(mW, weight);
	hFSR_mWg2->Fill(mWg, weight);
	hFSR_mW32->Fill(mW,mWg, weight);
	hFSR_dR2->Fill(dR, weight);
	hFSR_gpta->Fill(gpt, weight);
	hmT_FSRa->Fill(mT, weight);
	hmT3_FSRa->Fill(mT3, weight);
	hmT3mT_FSRa->Fill(mT,mT3, weight);
      

      }


  }
  outFile->Write();
  outFile->Close();

  fin.close();
  fin.clear();
  cout << "Ntotal = " << ntotal << "\t Npass = " << npass << std::endl;
  cout << "Nsubtotal = " << nsubtotal << endl;
  
  double eff = (double)npass/(double)ntotal;
  double err = sqrt( (1-eff)*eff/ntotal);
  cout << "filter eff = " << eff << " +- " << err << endl;

  cout << "totalxsec = " << totalxsec << endl;
  double xsec = totalxsec*eff*1e-3/1e-12;
  double xsec_err = totalxsec*err*1e-3/1e-12;
  
  cout << "Reduced x-section from " << totalxsec*1e-3/1e-12 << " pb to " << 
    xsec << "+-" << xsec_err << " pb" << endl;


  
}
