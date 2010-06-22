#define smearRC_cxx
#include "smearRC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom2.h>

void smearRC::Loop(int ieta)
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TH1F* hTemplate= new TH1F("hTemplate","",120,-1.0,11.0);
   TH1F* hsumiso = (TH1F*)hTemplate->Clone();
   hsumiso->SetName("hsumiso");
   TH1F* heciso = (TH1F*)hTemplate->Clone();
   heciso->SetName("heciso");
   TH1F* hhciso = (TH1F*)hTemplate->Clone();
   hhciso->SetName("hhciso");
   TH1F* htkiso = (TH1F*)hTemplate->Clone();
   htkiso->SetName("htkiso");
   

   double RCMean     = 5.00462e-01;
   double RCWidth    = 2.37595e-01;

   double spikeMean  = 0.451003;
   double spikeWidth = 0.325818;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (Cut(ientry) < 0) continue;
      
      if (fabs(eta)>1.45 && ieta==0)continue;
      if (fabs(eta)<1.70 && ieta==1)continue;
      if (fabs(eta)>2.50 && ieta==1)continue;

      double smear_eciso =  (ieta==0)?
   	(eciso-RCMean)*spikeWidth/RCWidth + spikeMean:
	gRandom->Gaus(eciso, 
		      sqrt(spikeWidth*spikeWidth - RCWidth*RCWidth));


      Float_t comb3Iso = tkiso+hciso+smear_eciso;
      if(comb3Iso  > 11.0)continue;

      htkiso->Fill(tkiso);
      hhciso->Fill(hciso);
      heciso->Fill(smear_eciso);
      hsumiso->Fill(comb3Iso);
   }
   char* dec[2] = {"EB","EE"};
   TFile* outFile = new TFile(Form("smearRC_%s.root",dec[ieta]),"recreate"); 
   htkiso->Write();
   heciso->Write();
   hhciso->Write();
   hsumiso->Write();
   outFile->Close();


}
