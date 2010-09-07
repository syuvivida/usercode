#define makeRCHisto_cxx
#include "makeRCHisto.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>

void makeRCHisto::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   // simple 1-D histograms
   TH1F* hTemplate= new TH1F("hTemplate","",120,-1.0,11.0);
   TH1F* hsumiso = (TH1F*)hTemplate->Clone();
   hsumiso->SetName("hsumiso");
   TH1F* heciso = (TH1F*)hTemplate->Clone();
   heciso->SetName("heciso");
   TH1F* hhciso = (TH1F*)hTemplate->Clone();
   hhciso->SetName("hhciso");
   TH1F* htkiso = (TH1F*)hTemplate->Clone();
   htkiso->SetName("htkiso");
  

   // isolation vs. number of vertices
   TProfile* nvtx_eciso = new TProfile("nvtx_eciso","",10,0.5,10.5); 
   nvtx_eciso->SetXTitle("Number of good vertices");
   nvtx_eciso->SetYTitle("<Iso_{ECAL}> (GeV)");

   TProfile* nvtx_hciso = new TProfile("nvtx_hciso","",10,0.5,10.5); 
   nvtx_hciso->SetXTitle("Number of good vertices");
   nvtx_hciso->SetYTitle("<Iso_{HCAL}> (GeV)");

   TProfile* nvtx_tkiso = new TProfile("nvtx_tkiso","",10,0.5,10.5); 
   nvtx_tkiso->SetXTitle("Number of good vertices");
   nvtx_tkiso->SetYTitle("<Iso_{TRK}> (GeV/c)");

   TProfile* nvtx_sumiso = new TProfile("nvtx_sumiso","",10,0.5,10.5); 
   nvtx_sumiso->SetXTitle("Number of good vertices");
   nvtx_sumiso->SetYTitle("<Iso> (GeV)");

   // isolation vs. number of pixel hits

   TProfile* npix_eciso = new TProfile("npix_eciso","",100,0.5,1000.5); 
   npix_eciso->SetXTitle("Number of pixel hits");
   npix_eciso->SetYTitle("<Iso_{ECAL}> (GeV)");

   TProfile* npix_hciso = new TProfile("npix_hciso","",100,0.5,1000.5); 
   npix_hciso->SetXTitle("Number of pixel hits");
   npix_hciso->SetYTitle("<Iso_{HCAL}> (GeV)");

   TProfile* npix_tkiso = new TProfile("npix_tkiso","",100,0.5,1000.5); 
   npix_tkiso->SetXTitle("Number of pixel hits");
   npix_tkiso->SetYTitle("<Iso_{TRK}> (GeV/c)");

   TProfile* npix_sumiso = new TProfile("npix_sumiso","",100,0.5,1000.5); 
   npix_sumiso->SetXTitle("Number of pixel hits");
   npix_sumiso->SetYTitle("<Iso> (GeV)");


   int ieta=0;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     
     
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if (fabs(eta)>1.4442 && ieta==0)continue;
      if (fabs(eta)<1.566 && ieta==1)continue;
      if (fabs(eta)>2.50 && ieta==1)continue;


      Float_t comb3Iso = tkiso+hciso+eciso;

      htkiso->Fill(tkiso);
      hhciso->Fill(hciso);
      heciso->Fill(eciso);

      if(comb3Iso  > 11.0)continue;

      hsumiso->Fill(comb3Iso);

      nvtx_eciso->Fill(nGoodVtx4,eciso);
      nvtx_hciso->Fill(nGoodVtx4,hciso);
      nvtx_tkiso->Fill(nGoodVtx4,tkiso);
      nvtx_sumiso->Fill(nGoodVtx4,comb3Iso);

      npix_eciso->Fill(nPix,eciso);
      npix_hciso->Fill(nPix,hciso);
      npix_tkiso->Fill(nPix,tkiso);
      npix_sumiso->Fill(nPix,comb3Iso);

      if (jentry%100==0)
	printf("%4.1f%% done.\r",(float)jentry/(float)nentries*100.); 
      

   }


   char* dec[2] = {"EB","EE"};
   TFile* outFile = new TFile(Form("%s_nGoodVtx4_%s.root",_inputDirName.data(),
				   _outputName.data(),dec[ieta]),"recreate"); 
   htkiso->Write();
   heciso->Write();
   hhciso->Write();
   hsumiso->Write();

   nvtx_eciso->Write();
   nvtx_hciso->Write();
   nvtx_tkiso->Write();
   nvtx_sumiso->Write();

   npix_eciso->Write();
   npix_hciso->Write();
   npix_tkiso->Write();
   npix_sumiso->Write();

   outFile->Close();

}
