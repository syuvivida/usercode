#include <iostream>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TH1.h>
#include <TRandom2.h>
#include <TSystem.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h> 
#include <TCanvas.h>
#include <TLegend.h>
#include "ext-mL_fit_unbinned.C"
#include "setTDRStyle.C"

using namespace std;

// the pt and eta binning
const Double_t fBinsEta[]={0,1.45,1.7,2.5};
const Double_t fBinsPt[]={15,20,30,50,80,120};
const int nPtBin = 5;//sizeof(fBinsPt)/sizeof(fBinsPt[0])-1;
const int nEtaBin = 2;//(sizeof(fBinsEta)/sizeof(fBinsEta[0]))/2;

const int REBINNINGS=1;


void compareCorrRCSpike(){
	

  char tmp[1000];
  setTDRStyle();
  // settings for purity TGraphAsymmetryErrors
  Double_t fBinsPtMidPoint[nPtBin]={0};
  Double_t fBinsPtError[nPtBin]={0};
  
  for(int ipt=0; ipt < nPtBin; ipt++)
    {
      fBinsPtMidPoint[ipt] = 0.5*(fBinsPt[ipt+1]+fBinsPt[ipt]);
      fBinsPtError[ipt] = 0.5*(fBinsPt[ipt+1]-fBinsPt[ipt]);
    }


  TH1F* hTemplate_S[nEtaBin][nPtBin];
  TH1F* hdata_S[nEtaBin][nPtBin];

  
  std::string dataFile     = "SBDataTemplate_131511_139239.root";
  std::string templateFile = "spike_131511_139239.root";

  TFile* inf_data = new TFile(dataFile.data());
  TFile* inf_template = new TFile(templateFile.data());

  
  char* dec[2] = {"EB","EE"};
  for(int ieta=0; ieta<1; ieta++){
    for(int ipt=0; ipt < nPtBin; ipt++){

      // getting histograms from data root file

      sprintf(tmp,"h_%s_comb3Iso_sig_pt%d",dec[ieta],(int)fBinsPt[ipt]);
      cout << "looking for histogram " << tmp << " in file " << 
	inf_data->GetName() << endl;
      hdata_S[ieta][ipt] = (TH1F*)inf_data->FindObjectAny(tmp);
      hdata_S[ieta][ipt]->SetName(Form("hdata_S_%d_%d",ieta,ipt));
      hdata_S[ieta][ipt]->Rebin(REBINNINGS);

      // getting histogram for template root file
      sprintf(tmp,"h_EB_comb3Iso_EGdata_SIG");
      cout << "looking for histogram " << tmp << " in file " << 
	inf_template->GetName() << endl;
      hTemplate_S[ieta][ipt] = (TH1F*)inf_template->FindObjectAny(tmp);
      hTemplate_S[ieta][ipt]->SetName(Form("hTemplate_S_%d_%d",ieta,ipt));
      hTemplate_S[ieta][ipt]->Rebin(REBINNINGS);

    }
  }


  TH1F* hfit_sig;
  TH1F* hfit_spike;

  TH1F* hsig = new TH1F("hsig","",120,-1,11);
  hsig->SetXTitle("Iso (GeV)");
  hsig->SetYTitle("A.U.");
  hsig->GetYaxis()->SetDecimals();
  hsig->SetYTitle("Iso (GeV)");

  hsig->SetLineColor(2);
  TH1F* hspike = new TH1F("hspike","",120,-1,11);
  hspike->SetLineColor(4);

  TCanvas* c1 = new TCanvas("c1","",500,500);

  for(int ieta=0; ieta< 1; ieta++){
    for(int ipt=0; ipt < nPtBin; ipt++){


      hsig->Reset();
      hspike->Reset();

      // first obtain fit for the random cone corrected values
      hfit_sig  = (TH1F*)hdata_S[ieta][ipt]->Clone();
      hfit_sig  -> SetName("hfit_sig");
      hfit_sig  -> Rebin(4);

      double sigPar[20] = {hfit_sig->GetMaximum(), 1., 0.5, 0.3,
			   1.,-.3,-1., 0.01, 0.5, 0.01, 1., 1.};

      TF1 *fsig = new TF1("fsig", exp_conv, -1., 11., 12);
      fsig->SetParameters(sigPar);
      fsig->SetNpx(2500);
      fsig->SetLineColor(2);
      hfit_sig->Fit(fsig,"","",-1,5.0);

      c1->Print(Form("hfit_sig_pt%d.gif",(int)fBinsPt[ipt]));
      
      fsig->SetParameter(0,2);
      fsig->SetParameter(1,fsig->GetParameter(1)*8.39614e-01/6.83080e-01);//correction from RC
      fsig->SetParameter(2,4.83182e-01);
      fsig->SetParameter(3,fsig->GetParameter(3)*2.33769e-01/2.26323e-01);
    
      cout << "fsig Integral = " << fsig->Integral(-1,11) << endl;
      for(int i=0;i<4;i++)
	cout << "fsig par " << i << "= " << fsig->GetParameter(i) << endl;

      // second obtain fit for the spikes
      hfit_spike  = (TH1F*)hTemplate_S[ieta][ipt]->Clone();
      hfit_spike  -> SetName("hfit_spike");
      hfit_spike  -> Rebin(4);

      TF1 *fspike = new TF1("fspike", exp_conv, -1., 11., 12);
      fspike->SetParameters(sigPar);
      fspike->SetParameter(0,hfit_spike->GetMaximum());
      fspike->SetLineColor(4);
      fspike->SetNpx(2500);
      hfit_spike->Fit(fspike,"","",-1,5.0);

      c1->Print(Form("hfit_spike_pt%d.gif",(int)fBinsPt[ipt]));
      
      fspike->SetParameter(0,2.0);
      cout << "fspike Integral = " << fspike->Integral(-1,11) << endl;
      float scale = (float)fsig->Integral(-1,11)/(float)fspike->Integral(-1,11);
      fspike->SetParameter(0,2.0*scale);
      cout << "now fspike Integral = " << fspike->Integral(-1,11) << endl;

      for(int i=0;i<4;i++)
	cout << "fspike par " << i << "= " << fspike->GetParameter(i) << endl;

    
      hsig->FillRandom("fsig",10000);
      hspike->FillRandom("fspike",10000);

      cout << "KS test for sideband data and sideband MC = " << 
	hsig->KolmogorovTest(hspike,"X") << endl;
      cout << "KS test 2 for sideband data and sideband MC = " << 
	hsig->KolmogorovTest(hspike) << endl;
     
      hsig->Reset();
      hsig->SetMaximum(1.2);
      hsig->Draw();
//       hspike->Draw("same");
      fsig->Draw("same");
      fspike->Draw("same");

      TLegend* leg = new TLegend(0.35,0.70,0.55,0.85);
      leg->SetHeader(Form("p_{T}[%d,%d]",(int)fBinsPt[ipt],(int)fBinsPt[ipt+1]));
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.04);
      leg->SetBorderSize(0);
      leg->AddEntry(fsig,"Random-cone-corrected MC","f");
      leg->AddEntry(fspike,"Spike","f");   
      leg->Draw("same");
      
      c1->Print(Form("CorrRCSpike_pt%d.eps",(int)fBinsPt[ipt]));
      c1->Print(Form("CorrRCSpike_pt%d.gif",(int)fBinsPt[ipt]));
      c1->Print(Form("CorrRCSpike_pt%d.C",(int)fBinsPt[ipt]));
      
    } 
  }

    
   

}
