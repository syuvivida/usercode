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
void read_Wg_baur(std::string inputfile)
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
  sprintf(name,"%s%s",inputfile.data(),"_baur.root");  
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
  
  const int NCOL_INT=30;
  const int NCOL_SPIN=6;
  const int NCOL_FLOAT=33;
  int code[NCOL_INT];
  double spin[NCOL_SPIN];
  double cl_mom[NCOL_FLOAT];

  for(int j=0;j<NCOL_INT;j++)
    fin >> code[j];
  for(int j=0;j<NCOL_SPIN;j++)
    fin >> spin[j];
  for(int j=0;j<NCOL_FLOAT;j++)
    fin >> cl_mom[j];


  while(!fin.eof()){
    ntotal ++;
    if(ntotal%100==0)cout << "ntotal = " << ntotal << endl;


    TLorentzVector gamma_lead(cl_mom[25],cl_mom[26],cl_mom[27],cl_mom[28]);
    TLorentzVector clepton(cl_mom[20],cl_mom[21],cl_mom[22],cl_mom[23]);
    TLorentzVector nlepton(cl_mom[15],cl_mom[16],cl_mom[17],cl_mom[18]);

    double mT3 = ClusterMass(gamma_lead,clepton,nlepton);

    double mW = (clepton+nlepton).M();
    double mWg= (clepton+nlepton+gamma_lead).M();
    double dR = gamma_lead.DeltaR(clepton);
    double gpt  =gamma_lead.Pt();
    double geta =gamma_lead.Eta();
    
    bool hasFSR= mT3 < 90.0;
    bool hasISR= mT3 > 90.0;
//     bool hasFSR= mWg < 90.0;
//     bool hasISR= mWg > 90.0;



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


    for(int j=0;j<NCOL_INT;j++)
      fin >> code[j];
    for(int j=0;j<NCOL_SPIN;j++)
      fin >> spin[j];
    for(int j=0;j<NCOL_FLOAT;j++)
      fin >> cl_mom[j];

    
    if((!hasFSR && !hasISR) || clepton.E() <1e-6 || nlepton.E()<1e-6)
      continue;

    
    nsubtotal++;
    
    if(hasISR){
      h_mW->Fill(mW);
      h_mWg->Fill(mWg);
      h_mW3->Fill(mW,mWg);
      h_dR->Fill(dR);

      hISR_dRgpt->Fill(gpt,dR);
      hISR_lptgpt->Fill(gpt,lpt);

      hISR_mW->Fill(mW);
      hISR_mWg->Fill(mWg);
      hISR_mW3->Fill(mW,mWg);
      hISR_dR->Fill(dR);

      h_gptb->Fill(gpt);
      h_lptb->Fill(lpt);
      h_nptb->Fill(npt);

      h_getab->Fill(gy);
      h_letab->Fill(ly);
      h_netab->Fill(neta);

      hISR_gptb->Fill(gpt);

      hmT_ISRb->Fill(mT);
      hmT3_ISRb->Fill(mT3);
      hmT3mT_ISRb->Fill(mT,mT3);
      
    }

    if(hasFSR)
      {

	h_mW->Fill(mW);
	h_mWg->Fill(mWg);
	h_mW3->Fill(mW,mWg);
	h_dR->Fill(dR);

	hFSR_dRgpt->Fill(gpt,dR);
	hFSR_lptgpt->Fill(gpt,lpt);

	hFSR_mW->Fill(mW);
	hFSR_mWg->Fill(mWg);
	hFSR_mW3->Fill(mW,mWg);
	hFSR_dR->Fill(dR);

	h_gptb->Fill(gpt);
	h_lptb->Fill(lpt);
	h_nptb->Fill(npt);

	h_getab->Fill(gy);
	h_letab->Fill(ly);
	h_netab->Fill(neta);

	hFSR_gptb->Fill(gpt);

	hmT_FSRb->Fill(mT);
	hmT3_FSRb->Fill(mT3);
	hmT3mT_FSRb->Fill(mT,mT3);

      }
          

//     // muon channel LSP's cuts
//     bool photonCut_ISR = gpt>10.0 && fabs(geta)<2.7;
//     if(!photonCut_ISR)continue;
//     if(lpt<5.0)continue;
//     if(npt<20.0)continue;
//     if(fabs(leta)>2.7)continue;

    // Poter's cuts, supposedly
    
//      bool photonCut_ISR = gpt>20.0 && fabs(geta)<2.7 && dR>1.1;
//      if(!photonCut_ISR)continue;
//      if(lpt<15.0)continue;
//      if(fabs(leta)>2.7)continue;
//      if(mT < 60)continue;
//      if(npt<20.0)continue;

//      bool photonCut_ISR = gpt>20.0 && fabs(geta)<2.7;
//      if(!photonCut_ISR)continue;
//      if(lpt<15.0)continue;
//      if(fabs(leta)>2.7)continue;
//      if(mT < 30 || mT>120)continue;
//      if(npt<20.0)continue;


    // ATLAS Zhijun's cuts

     bool photonCut_ISR = gpt>10.0 && fabs(geta)<2.5 && dR>0.7;

        if(!photonCut_ISR)continue;
     if(lpt<10.0)continue;
     if(fabs(leta)>2.5)continue;


    // CDF cuts
//     bool photonCut_ISR = gpt> 7.0 && dR>0.7;
//     if(!photonCut_ISR)continue;
//       if(mT3 < 90)continue;


    npass++;
    
    hmT3->Fill(mT3);

    
    if(hasISR){
      h_gpta->Fill(gpt);
      h_lpta->Fill(lpt);
      h_npta->Fill(npt);

      h_getaa->Fill(gy);
      h_letaa->Fill(ly);
      h_netaa->Fill(neta);
    

      h_mW2->Fill(mW);
      h_mWg2->Fill(mWg);
      h_mW32->Fill(mW,mWg);
      h_dR2->Fill(dR);

      
      hISR_mW2->Fill(mW);
      hISR_mWg2->Fill(mWg);
      hISR_mW32->Fill(mW,mWg);
      hISR_dR2->Fill(dR);
      hISR_gpta->Fill(gpt);
      hmT_ISRa->Fill(mT);
      hmT3_ISRa->Fill(mT3);
      hmT3mT_ISRa->Fill(mT,mT3);
    

    }

   if(hasFSR)
      {
	h_mW2->Fill(mW);
	h_mWg2->Fill(mWg);
	h_mW32->Fill(mW,mWg);
	h_dR2->Fill(dR);

	h_gpta->Fill(gpt);
	h_lpta->Fill(lpt);
	h_npta->Fill(npt);

	h_getaa->Fill(gy);
	h_letaa->Fill(ly);
	h_netaa->Fill(neta);

	hFSR_mW2->Fill(mW);
	hFSR_mWg2->Fill(mWg);
	hFSR_mW32->Fill(mW,mWg);
	hFSR_dR2->Fill(dR);
	hFSR_gpta->Fill(gpt);
	hmT_FSRa->Fill(mT);
	hmT3_FSRa->Fill(mT3);
	hmT3mT_FSRa->Fill(mT,mT3);
      

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


  
}
