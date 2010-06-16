#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"
// the pt and eta binning
const double fBinsEta[]={0,1.45,1.7,2.5};
const double fBinsPt[]={15,20,30,50,80,120};
const int nPtBin = 5;
const int nEtaBin = 2;


void fitPull()
{
  setTDRStyle();
  gStyle->SetOptStat(0);
  char tmp[1000];
  TH1F hTemp;
  char* dec[2] = {"EB","EE"};
  TCanvas* c1 = new TCanvas("c1","",500,500);
  for(int ieta=0; ieta< nEtaBin; ieta++)
    {
      TFile* inf_data = new TFile(Form("pull_%s.root",dec[ieta]));
      for(int ipt=0; ipt < nPtBin; ipt++){
      sprintf(tmp,"hpull_Eta_%.2f_%.2f_Et_%d_%d",		
		    fBinsEta[ieta*2],fBinsEta[ieta*2+1],
	      (int)fBinsPt[ipt], (int)fBinsPt[ipt+1]);
      
      hTemp = (TH1F*)inf_data->FindObjectAny(tmp);
      sprintf(tmp,"%.2f < |#eta(#gamma)| < %.2f, %d < p_{T}(#gamma) < %d GeV",
	      fBinsEta[ieta*2],fBinsEta[ieta*2+1],
	      (int)fBinsPt[ipt], (int)fBinsPt[ipt+1]);
      hTemp->SetTitle(tmp);
      hTemp->Fit("gaus");
      hTemp->SetXTitle("(n_{sig}^{fit}-n_{sig}^{input})/err^{fit}");

      sprintf(tmp,"hpull_Eta_%.2f_%.2f_Et_%d_%d",		
		    fBinsEta[ieta*2],fBinsEta[ieta*2+1],
	      (int)fBinsPt[ipt], (int)fBinsPt[ipt+1]);
      
      c1->Print(Form("%s.eps",tmp));
      c1->Print(Form("%s.gif",tmp));
      }

    }







}
