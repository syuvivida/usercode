#include <iostream>
#include <fstream>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>

using namespace std;
void read_Zg(std::string inputfile, double totalxsec)
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
  
  int ndRbin=200;
  float dRmin=0;
  float dRmax=10.0;

  char name[300];
  sprintf(name,"%s%s",inputfile.data(),"_histo.root");  
  TFile* outFile = new TFile(name,"recreate");

  TH1F* h_diff = new TH1F("h_diff","P_{T}(#gamma^{ISR})-P_{T}(#gamma^{FSR})",
			  100,-100,100);
  TH1F* h_diff2 = new TH1F("h_diff2","P_{T}(#gamma^{ISR})-P_{T}(#gamma^{FSR})",
			  100,-100,100);

  TH1F* h_dR = new TH1F("h_dR","#Delta R (lepton,#gamma)",ndRbin,dRmin,dRmax);
  TH1F* h_mZ = new TH1F("h_mZ","Mass of Z boson",netbin,etmin,etmax);  
  TH1F* h_mZg = new TH1F("h_mZg","Mass of Z boson + photon",netbin,etmin,etmax);
  TH2F* h_mZ3 = new TH2F("h_mZ3","M(Z#gamma) vs M(Z)",netbin,etmin,etmax,netbin,etmin,etmax);
  TH1F* h_dR2 = new TH1F("h_dR2","#Delta R (lepton,#gamma)",ndRbin,dRmin,dRmax);
  TH1F* h_mZ2 = new TH1F("h_mZ2","Mass of Z boson",netbin,etmin,etmax);  
  TH1F* h_mZg2 = new TH1F("h_mZg2","Mass of Z boson + photon",netbin,etmin,etmax);
  TH2F* h_mZ32 = new TH2F("h_mZ32","M(Z#gamma) vs M(Z)",netbin,etmin,etmax,netbin,etmin,etmax);

  TH1F* hFSR_dR = new TH1F("hFSR_dR","#Delta R (lepton,#gamma)",ndRbin,dRmin,dRmax);
  TH1F* hFSR_mZ = new TH1F("hFSR_mZ","Mass of Z boson",netbin,etmin,etmax);  
  TH1F* hFSR_mZg = new TH1F("hFSR_mZg","Mass of Z boson + photon",netbin,etmin,etmax);
  TH2F* hFSR_mZ3 = new TH2F("hFSR_mZ3","M(Z#gamma) vs M(Z)",netbin,etmin,etmax,netbin,etmin,etmax);
  TH1F* hFSR_dR2 = new TH1F("hFSR_dR2","#Delta R (lepton,#gamma)",ndRbin,dRmin,dRmax);
  TH1F* hFSR_mZ2 = new TH1F("hFSR_mZ2","Mass of Z boson",netbin,etmin,etmax);  
  TH1F* hFSR_mZg2 = new TH1F("hFSR_mZg2","Mass of Z boson + photon",netbin,etmin,etmax);
  TH2F* hFSR_mZ32 = new TH2F("hFSR_mZ32","M(Z#gamma) vs M(Z)",netbin,etmin,etmax,netbin,etmin,etmax);

  TH1F* hISR_dR = new TH1F("hISR_dR","#Delta R (lepton,#gamma)",ndRbin,dRmin,dRmax);
  TH1F* hISR_mZ = new TH1F("hISR_mZ","Mass of Z boson",netbin,etmin,etmax);  
  TH1F* hISR_mZg = new TH1F("hISR_mZg","Mass of Z boson + photon",netbin,etmin,etmax);
  TH2F* hISR_mZ3 = new TH2F("hISR_mZ3","M(Z#gamma) vs M(Z)",netbin,etmin,etmax,netbin,etmin,etmax);
  TH1F* hISR_dR2 = new TH1F("hISR_dR2","#Delta R (lepton,#gamma)",ndRbin,dRmin,dRmax);
  TH1F* hISR_mZ2 = new TH1F("hISR_mZ2","Mass of Z boson",netbin,etmin,etmax);  
  TH1F* hISR_mZg2 = new TH1F("hISR_mZg2","Mass of Z boson + photon",netbin,etmin,etmax);
  TH2F* hISR_mZ32 = new TH2F("hISR_mZ32","M(Z#gamma) vs M(Z)",netbin,etmin,etmax,netbin,etmin,etmax);


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

  // after requiring only one leading photon

  

  TH1F* hlead_dR = new TH1F("hlead_dR","#Delta R (lepton,#gamma)",ndRbin,dRmin,dRmax);
  TH1F* hlead_mZ = new TH1F("hlead_mZ","Mass of Z boson",netbin,etmin,etmax);  
  TH1F* hlead_mZg = new TH1F("hlead_mZg","Mass of Z boson + photon",netbin,etmin,etmax);
  TH2F* hlead_mZ3 = new TH2F("hlead_mZ3","M(Z#gamma) vs M(Z)",netbin,etmin,etmax,netbin,etmin,etmax);
  TH1F* hlead_dR2 = new TH1F("hlead_dR2","#Delta R (lepton,#gamma)",ndRbin,dRmin,dRmax);
  TH1F* hlead_mZ2 = new TH1F("hlead_mZ2","Mass of Z boson",netbin,etmin,etmax);  
  TH1F* hlead_mZg2 = new TH1F("hlead_mZg2","Mass of Z boson + photon",netbin,etmin,etmax);
  TH2F* hlead_mZ32 = new TH2F("hlead_mZ32","M(Z#gamma) vs M(Z)",netbin,etmin,etmax,netbin,etmin,etmax);

  TH1F* hleadFSR_dR = new TH1F("hleadFSR_dR","#Delta R (lepton,#gamma)",ndRbin,dRmin,dRmax);
  TH1F* hleadFSR_mZ = new TH1F("hleadFSR_mZ","Mass of Z boson",netbin,etmin,etmax);  
  TH1F* hleadFSR_mZg = new TH1F("hleadFSR_mZg","Mass of Z boson + photon",netbin,etmin,etmax);
  TH2F* hleadFSR_mZ3 = new TH2F("hleadFSR_mZ3","M(Z#gamma) vs M(Z)",netbin,etmin,etmax,netbin,etmin,etmax);
  TH1F* hleadFSR_dR2 = new TH1F("hleadFSR_dR2","#Delta R (lepton,#gamma)",ndRbin,dRmin,dRmax);
  TH1F* hleadFSR_mZ2 = new TH1F("hleadFSR_mZ2","Mass of Z boson",netbin,etmin,etmax);  
  TH1F* hleadFSR_mZg2 = new TH1F("hleadFSR_mZg2","Mass of Z boson + photon",netbin,etmin,etmax);
  TH2F* hleadFSR_mZ32 = new TH2F("hleadhFSR_mZ32","M(Z#gamma) vs M(Z)",netbin,etmin,etmax,netbin,etmin,etmax);

  TH1F* hleadISR_dR = new TH1F("hleadISR_dR","#Delta R (lepton,#gamma)",ndRbin,dRmin,dRmax);
  TH1F* hleadISR_mZ = new TH1F("hleadISR_mZ","Mass of Z boson",netbin,etmin,etmax);  
  TH1F* hleadISR_mZg = new TH1F("hleadISR_mZg","Mass of Z boson + photon",netbin,etmin,etmax);
  TH2F* hleadISR_mZ3 = new TH2F("hleadISR_mZ3","M(Z#gamma) vs M(Z)",netbin,etmin,etmax,netbin,etmin,etmax);
  TH1F* hleadISR_dR2 = new TH1F("hleadISR_dR2","#Delta R (lepton,#gamma)",ndRbin,dRmin,dRmax);
  TH1F* hleadISR_mZ2 = new TH1F("hleadISR_mZ2","Mass of Z boson",netbin,etmin,etmax);  
  TH1F* hleadISR_mZg2 = new TH1F("hleadISR_mZg2","Mass of Z boson + photon",netbin,etmin,etmax);
  TH2F* hleadISR_mZ32 = new TH2F("hleadISR_mZ32","M(Z#gamma) vs M(Z)",netbin,etmin,etmax,netbin,etmin,etmax);


  TH1F* hlead_gptb = new TH1F("hlead_gptb","p_{T}(#gamma) before cuts",netbin,etmin,etmax);
  TH1F* hlead_gpta = new TH1F("hlead_gpta","p_{T}(#gamma) after cuts",netbin,etmin,etmax);
  TH1F* hlead_lptb = new TH1F("hlead_lptb","p_{T}(e) before cuts",netbin,etmin,etmax);
  TH1F* hlead_lpta = new TH1F("hlead_lpta","p_{T}(e) after cuts",netbin,etmin,etmax);
  TH1F* hlead_nptb = new TH1F("hlead_nptb","p_{T}(#nu) before cuts",netbin,etmin,etmax);
  TH1F* hlead_npta = new TH1F("hlead_npta","p_{T}(#nu) after cuts",netbin,etmin,etmax);

  TH1F* hlead_getab = new TH1F("hlead_getab","#eta(#gamma) before cuts",netabin,etamin,etamax);
  TH1F* hlead_getaa = new TH1F("hlead_getaa","#eta(#gamma) after cuts",netabin,etamin,etamax);
  TH1F* hlead_letab = new TH1F("hlead_letab","#eta(e) before cuts",netabin,etamin,etamax);
  TH1F* hlead_letaa = new TH1F("hlead_letaa","#eta(e) after cuts",netabin,etamin,etamax);
  TH1F* hlead_netab = new TH1F("hlead_netab","#eta(#nu) before cuts",netabin,etamin,etamax);
  TH1F* hlead_netaa = new TH1F("hlead_netaa","#eta(#nu) after cuts",netabin,etamin,etamax);

  



  // Looping over entries

  ifstream fin;
  sprintf(name,"%s%s",inputfile.data(),".dat");
  fin.open(name);
  cout << "Opening " << name << endl;
  
  const int NCOL=16; // photon, charged lepton, neutral lepton
  double cl[NCOL];

  for(int j=0;j<NCOL;j++)
    fin >> cl[j];

  while(!fin.eof()){
    ntotal ++;
    if(ntotal%100==0)cout << "ntotal = " << ntotal << endl;
    TLorentzVector gamma(cl[0],cl[1],cl[2],cl[3]);
    TLorentzVector clepton(cl[4],cl[5],cl[6],cl[7]);
    TLorentzVector nlepton(cl[8],cl[9],cl[10],cl[11]);
    TLorentzVector gamma_FSR(cl[12],cl[13],cl[14],cl[15]);

    double mZ = (clepton+nlepton).M();
    double mZg= (clepton+nlepton+gamma).M();
    double dR1= gamma.DeltaR(clepton);
    double dR2= gamma.DeltaR(nlepton);
    double gpt  =gamma.Pt();
    double geta =gamma.Eta();


    double mZg_FSR  = (clepton+nlepton+gamma_FSR).M();
    double dR_FSR1  = gamma_FSR.DeltaR(clepton);
    double dR_FSR2  = gamma_FSR.DeltaR(nlepton);
    double gpt_FSR  = gamma_FSR.Pt();
    double geta_FSR = gamma_FSR.Eta();


    TLorentzVector met_lv(nlepton.Px(),                    
			  nlepton.Py(),
			  0,
			  nlepton.Pt());




    double lpt  =clepton.Pt();
    double leta =clepton.Eta();
    double npt  =nlepton.Pt();
    double neta =nlepton.Eta();

    bool leadIsISR = gpt > gpt_FSR;
    bool leadIsFSR = gpt_FSR >gpt;

    for(int j=0;j<NCOL;j++)
      fin >> cl[j];

    bool hasFSR= gamma_FSR.E()>1e-6;
    bool hasISR= gamma.E()>1e-6;
    
    if((!hasFSR && !hasISR) || clepton.E() <1e-6 || nlepton.E()<1e-6)
      continue;

    cout << gpt << "\t" << gpt_FSR << endl;
    
    nsubtotal++;
    
    if(hasISR){
      h_mZ->Fill(mZ);
      h_mZg->Fill(mZg);
      h_mZ3->Fill(mZ,mZg);
      h_dR->Fill(dR1);
      h_dR->Fill(dR2);

      hISR_mZ->Fill(mZ);
      hISR_mZg->Fill(mZg);
      hISR_mZ3->Fill(mZ,mZg);
      hISR_dR->Fill(dR1);
      hISR_dR->Fill(dR2);

      h_gptb->Fill(gpt);
      h_lptb->Fill(lpt);
      h_nptb->Fill(npt);

      h_getab->Fill(geta);
      h_letab->Fill(leta);
      h_netab->Fill(neta);

      hISR_gptb->Fill(gpt);

      
      if(leadIsISR){
	hlead_mZ->Fill(mZ);
	hlead_mZg->Fill(mZg);
	hlead_mZ3->Fill(mZ,mZg);
	hlead_dR->Fill(dR1);
	hlead_dR->Fill(dR2);

	hleadISR_mZ->Fill(mZ);
	hleadISR_mZg->Fill(mZg);
	hleadISR_mZ3->Fill(mZ,mZg);
	hleadISR_dR->Fill(dR1);
	hleadISR_dR->Fill(dR2);

	hlead_gptb->Fill(gpt);
	hlead_lptb->Fill(lpt);
	hlead_nptb->Fill(npt);

	hlead_getab->Fill(geta);
	hlead_letab->Fill(leta);
	hlead_netab->Fill(neta);
      }

    }

    if(hasFSR)
      {

	h_mZ->Fill(mZ);
	h_mZg->Fill(mZg_FSR);
	h_mZ3->Fill(mZ,mZg_FSR);
	h_dR->Fill(dR_FSR1);
	h_dR->Fill(dR_FSR2);

	hFSR_mZ->Fill(mZ);
	hFSR_mZg->Fill(mZg_FSR);
	hFSR_mZ3->Fill(mZ,mZg_FSR);
	hFSR_dR->Fill(dR_FSR1);
	hFSR_dR->Fill(dR_FSR2);

	h_gptb->Fill(gpt_FSR);
	h_lptb->Fill(lpt);
	h_nptb->Fill(npt);

	h_getab->Fill(geta_FSR);
	h_letab->Fill(leta);
	h_netab->Fill(neta);

	hFSR_gptb->Fill(gpt_FSR);


	if(leadIsFSR){
	  hlead_mZ->Fill(mZ);
	  hlead_mZg->Fill(mZg_FSR);
	  hlead_mZ3->Fill(mZ,mZg_FSR);
	  hlead_dR->Fill(dR_FSR1);
	  hlead_dR->Fill(dR_FSR2);

	  hleadFSR_mZ->Fill(mZ);
	  hleadFSR_mZg->Fill(mZg_FSR);
	  hleadFSR_mZ3->Fill(mZ,mZg_FSR);
	  hleadFSR_dR->Fill(dR_FSR1);
	  hleadFSR_dR->Fill(dR_FSR2);

	  hlead_gptb->Fill(gpt_FSR);
	  hlead_lptb->Fill(lpt);
	  hlead_nptb->Fill(npt);

	  hlead_getab->Fill(geta_FSR);
	  hlead_letab->Fill(leta);
	  hlead_netab->Fill(neta);
	}

      }
          
    if(hasFSR && hasISR)h_diff->Fill(gpt-gpt_FSR);


    // Poter's cuts, supposedly
    
//     bool photonCut_ISR = gpt>10.0&& fabs(geta)<2.7 && dR1>0.7 && dR2>0.7;
//     bool photonCut_FSR = gpt_FSR > 10.0 && fabs(geta_FSR)<2.7 && dR_FSR1>0.7 && dR_FSR2>0.7;
// //     if(!photonCut_ISR && !photonCut_FSR)continue;
//     if(!photonCut_ISR)continue;
//     if(lpt< 3.0)continue;
//     if(npt< 3.0)continue;
//     if(fabs(leta)>2.7)continue;
//     if(fabs(neta)>2.7)continue;

    // CDF cuts
    bool photonCut_ISR = gpt> 7.0 && dR1>0.7 && dR2>0.7;
    bool photonCut_FSR = gpt_FSR > 7.0 && dR_FSR1>0.7 && dR_FSR2>0.7;
    if(!photonCut_ISR && !photonCut_FSR)continue;
//     if(!photonCut_ISR)continue;
    if(mZ<40)continue;


    npass++;
    
    if(hasISR){
      h_gpta->Fill(gpt);
      h_lpta->Fill(lpt);
      h_npta->Fill(npt);

      h_getaa->Fill(geta);
      h_letaa->Fill(leta);
      h_netaa->Fill(neta);
    

      h_mZ2->Fill(mZ);
      h_mZg2->Fill(mZg);
      h_mZ32->Fill(mZ,mZg);
      h_dR2->Fill(dR1);
      h_dR2->Fill(dR2);

      
      if(photonCut_ISR){
	hISR_mZ2->Fill(mZ);
	hISR_mZg2->Fill(mZg);
	hISR_mZ32->Fill(mZ,mZg);
	hISR_dR2->Fill(dR1);
	hISR_dR2->Fill(dR2);
	hISR_gpta->Fill(gpt);
      }


      if(leadIsISR && photonCut_ISR){
	hlead_gpta->Fill(gpt);
	hlead_lpta->Fill(lpt);
	hlead_npta->Fill(npt);

	hlead_getaa->Fill(geta);
	hlead_letaa->Fill(leta);
	hlead_netaa->Fill(neta);
    

	hlead_mZ2->Fill(mZ);
	hlead_mZg2->Fill(mZg);
	hlead_mZ32->Fill(mZ,mZg);
	hlead_dR2->Fill(dR1);
	hlead_dR2->Fill(dR2);

	hleadISR_mZ2->Fill(mZ);
	hleadISR_mZg2->Fill(mZg);
	hleadISR_mZ32->Fill(mZ,mZg);
	hleadISR_dR2->Fill(dR1);
	hleadISR_dR2->Fill(dR2);
      }

    }

   if(hasFSR)
      {
	h_mZ2->Fill(mZ);
	h_mZg2->Fill(mZg_FSR);
	h_mZ32->Fill(mZ,mZg_FSR);
	h_dR2->Fill(dR_FSR1);
	h_dR2->Fill(dR_FSR2);

	h_gpta->Fill(gpt_FSR);
	h_lpta->Fill(lpt);
	h_npta->Fill(npt);

	h_getaa->Fill(geta_FSR);
	h_letaa->Fill(leta);
	h_netaa->Fill(neta);

	if(photonCut_FSR){
	  hFSR_mZ2->Fill(mZ);
	  hFSR_mZg2->Fill(mZg_FSR);
	  hFSR_mZ32->Fill(mZ,mZg_FSR);
	  hFSR_dR2->Fill(dR_FSR1);
	  hFSR_dR2->Fill(dR_FSR2);
	  hFSR_gpta->Fill(gpt_FSR);
	}

	if(leadIsFSR && photonCut_FSR){
	  hlead_mZ2->Fill(mZ);
	  hlead_mZg2->Fill(mZg_FSR);
	  hlead_mZ32->Fill(mZ,mZg_FSR);
	  hlead_dR2->Fill(dR_FSR1);
	  hlead_dR2->Fill(dR_FSR2);

	  hleadFSR_mZ2->Fill(mZ);
	  hleadFSR_mZg2->Fill(mZg_FSR);
	  hleadFSR_mZ32->Fill(mZ,mZg_FSR);
	  hleadFSR_dR2->Fill(dR_FSR1);
	  hleadFSR_dR2->Fill(dR_FSR2);

	  hlead_gpta->Fill(gpt_FSR);
	  hlead_lpta->Fill(lpt);
	  hlead_npta->Fill(npt);

	  hlead_getaa->Fill(geta_FSR);
	  hlead_letaa->Fill(leta);
	  hlead_netaa->Fill(neta);

	}

      }

    if(hasFSR && hasISR)h_diff2->Fill(gpt-gpt_FSR);


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

//   sprintf(name,inputfile.data(),"_xsec");
//   fin.open(name);
//   fin.close();
  cout << "totalxsec = " << totalxsec << endl;
  double xsec = totalxsec*eff*1e-3/1e-12;
  double xsec_err = totalxsec*err*1e-3/1e-12;
  
  cout << "Reduced x-section from " << totalxsec*1e-3/1e-12 << " pb to " << 
    xsec << "+-" << xsec_err << " pb" << endl;


  
}
