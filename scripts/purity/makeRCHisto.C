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

   const int nvtxs = 5;
   TH1F* hsumison[nvtxs];
   TH1F* hecison[nvtxs];
   TH1F* hhcison[nvtxs];
   TH1F* htkison[nvtxs];

   for(int iv=0; iv<nvtxs; iv++){

     hsumison[iv] = (TH1F*)hTemplate->Clone();
     hsumison[iv]->SetName(Form("hsumison%d",iv+1));

     hecison[iv] = (TH1F*)hTemplate->Clone();
     hecison[iv]->SetName(Form("hecison%d",iv+1));

     hhcison[iv] = (TH1F*)hTemplate->Clone();
     hhcison[iv]->SetName(Form("hhcison%d",iv+1));

     htkison[iv] = (TH1F*)hTemplate->Clone();
     htkison[iv]->SetName(Form("htkison%d",iv+1));
			  
   }

   // isolation vs. number of vertices (is not fake)
   
   TProfile* pTemplate= new TProfile("pTemplate","",10,0.5,10.5);
   TProfile* nvtx_eciso = (TProfile*)pTemplate->Clone();
   nvtx_eciso->SetName("nvtx_eciso");
   nvtx_eciso->SetXTitle("Number of vertices (!isFake)");
   nvtx_eciso->SetYTitle("<Iso_{ECAL}> (GeV)");

   TProfile* nvtx_hciso = (TProfile*)pTemplate->Clone();
   nvtx_hciso->SetName("nvtx_hciso");
   nvtx_hciso->SetXTitle("Number of vertices (!isFake)");
   nvtx_hciso->SetYTitle("<Iso_{HCAL}> (GeV)");

   TProfile* nvtx_tkiso = (TProfile*)pTemplate->Clone();
   nvtx_tkiso->SetName("nvtx_tkiso");
   nvtx_tkiso->SetXTitle("Number of vertices (!isFake)");
   nvtx_tkiso->SetYTitle("<Iso_{TRK}> (GeV/c)");

   TProfile* nvtx_sumiso = (TProfile*)pTemplate->Clone();
   nvtx_sumiso->SetName("nvtx_sumiso");
   nvtx_sumiso->SetXTitle("Number of vertices (!isFake)");
   nvtx_sumiso->SetYTitle("<Iso> (GeV)");

   // isolation vs. number of vertices (is not fake && >=3 tracks)
   TProfile* nvtx3_eciso = (TProfile*)pTemplate->Clone();
   nvtx3_eciso->SetName("nvtx3_eciso");
   nvtx3_eciso->SetXTitle("Number of vertices (!isFake && nTrk>=3)");
   nvtx3_eciso->SetYTitle("<Iso_{ECAL}> (GeV)");

   TProfile* nvtx3_hciso = (TProfile*)pTemplate->Clone();
   nvtx3_hciso->SetName("nvtx3_hciso");
   nvtx3_hciso->SetXTitle("Number of vertices (!isFake && nTrk>=3)");
   nvtx3_hciso->SetYTitle("<Iso_{HCAL}> (GeV)");

   TProfile* nvtx3_tkiso = (TProfile*)pTemplate->Clone();
   nvtx3_tkiso->SetName("nvtx3_tkiso");
   nvtx3_tkiso->SetXTitle("Number of vertices (!isFake && nTrk>=3)");
   nvtx3_tkiso->SetYTitle("<Iso_{TRK}> (GeV/c)");

   TProfile* nvtx3_sumiso = (TProfile*)pTemplate->Clone();
   nvtx3_sumiso->SetName("nvtx3_sumiso");
   nvtx3_sumiso->SetXTitle("Number of vertices (!isFake && nTrk>=3)");
   nvtx3_sumiso->SetYTitle("<Iso> (GeV)");

   // isolation vs. number of vertices (is not fake && >=4 tracks)
   TProfile* nvtx4_eciso = (TProfile*)pTemplate->Clone();
   nvtx4_eciso->SetName("nvtx4_eciso");
   nvtx4_eciso->SetXTitle("Number of vertices (!isFake && nTrk>=4)");
   nvtx4_eciso->SetYTitle("<Iso_{ECAL}> (GeV)");

   TProfile* nvtx4_hciso = (TProfile*)pTemplate->Clone();
   nvtx4_hciso->SetName("nvtx4_hciso");
   nvtx4_hciso->SetXTitle("Number of vertices (!isFake && nTrk>=4)");
   nvtx4_hciso->SetYTitle("<Iso_{HCAL}> (GeV)");

   TProfile* nvtx4_tkiso = (TProfile*)pTemplate->Clone();
   nvtx4_tkiso->SetName("nvtx4_tkiso");
   nvtx4_tkiso->SetXTitle("Number of vertices (!isFake && nTrk>=4)");
   nvtx4_tkiso->SetYTitle("<Iso_{TRK}> (GeV/c)");

   TProfile* nvtx4_sumiso = (TProfile*)pTemplate->Clone();
   nvtx4_sumiso->SetName("nvtx4_sumiso");
   nvtx4_sumiso->SetXTitle("Number of vertices (!isFake && nTrk>=4)");
   nvtx4_sumiso->SetYTitle("<Iso> (GeV)");

   // isolation vs. number of pixel hits

   TProfile* ppTemplate= new TProfile("ppTemplate","",100,0.5,1000.5); 
   TProfile* npix_eciso = (TProfile*)ppTemplate->Clone();
   npix_eciso->SetName("npix_eciso");
   npix_eciso->SetXTitle("Number of pixel hits");
   npix_eciso->SetYTitle("<Iso_{ECAL}> (GeV)");

   TProfile* npix_hciso = (TProfile*)ppTemplate->Clone();
   npix_hciso->SetName("npix_hciso");
   npix_hciso->SetXTitle("Number of pixel hits");
   npix_hciso->SetYTitle("<Iso_{HCAL}> (GeV)");

   TProfile* npix_tkiso = (TProfile*)ppTemplate->Clone();
   npix_tkiso->SetName("npix_tkiso");
   npix_tkiso->SetXTitle("Number of pixel hits");
   npix_tkiso->SetYTitle("<Iso_{TRK}> (GeV/c)");

   TProfile* npix_sumiso = (TProfile*)ppTemplate->Clone();
   npix_sumiso->SetName("npix_sumiso");
   npix_sumiso->SetXTitle("Number of pixel hits");
   npix_sumiso->SetYTitle("<Iso> (GeV)");


   int ieta=0;
   cout << "Total entries = " << nentries << endl;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     

      if (jentry%100==0)
	printf("%4.1f%% done.\r",(float)jentry/(float)nentries*100.); 

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if (fabs(eta)>1.4442 && ieta==0)continue;
      if (fabs(eta)<1.566 && ieta==1)continue;
      if (fabs(eta)>2.50 && ieta==1)continue;

            
      Float_t comb3Iso = tkiso+hciso+eciso;

      const Float_t upper = 10.0;


      htkiso->Fill(tkiso);
      hhciso->Fill(hciso);
      heciso->Fill(eciso);
      hsumiso->Fill(comb3Iso);

      // separate based on the number of vertices
      
      int tempIndex = nGoodVtx4-1;
      if(nGoodVtx4<6){
	hsumison[tempIndex]->Fill(comb3Iso);
	hecison[tempIndex]->Fill(eciso);
	hhcison[tempIndex]->Fill(hciso);	
	htkison[tempIndex]->Fill(tkiso);
      }

      // filling profiles

      if(eciso < upper){
	nvtx_eciso->Fill(nVtx,eciso);
	nvtx3_eciso->Fill(nGoodVtx3,eciso);
	nvtx4_eciso->Fill(nGoodVtx4,eciso);
	npix_eciso->Fill(nPix,eciso);
      }

      if(hciso < upper){
	nvtx_hciso->Fill(nVtx,hciso);
	nvtx3_hciso->Fill(nGoodVtx3,hciso);
	nvtx4_hciso->Fill(nGoodVtx4,hciso);
	npix_hciso->Fill(nPix,hciso);
      }

      if(tkiso < upper){
	nvtx_tkiso->Fill(nVtx,tkiso);
	nvtx3_tkiso->Fill(nGoodVtx3,tkiso);
	nvtx4_tkiso->Fill(nGoodVtx4,tkiso);
	npix_tkiso->Fill(nPix,tkiso);
      }

      
      if(comb3Iso < upper){
	nvtx_sumiso->Fill(nVtx,comb3Iso);
	nvtx3_sumiso->Fill(nGoodVtx3,comb3Iso);
	nvtx4_sumiso->Fill(nGoodVtx4,comb3Iso);
	npix_sumiso->Fill(nPix,comb3Iso);
      }

      

   }


   char* dec[2] = {"EB","EE"};
   TFile* outFile = new TFile(Form("%s_%s.root",_inputDirName.data(),
				   _outputName.data(),dec[ieta]),"recreate"); 
   htkiso->Write();
   heciso->Write();
   hhciso->Write();
   hsumiso->Write();

   
   for(int iv=0; iv<nvtxs; iv++){

     hsumison[iv]->Write();
     hecison[iv]->Write();
     hhcison[iv]->Write();
     htkison[iv]->Write();

   }

   nvtx_eciso->Write();
   nvtx_hciso->Write();
   nvtx_tkiso->Write();
   nvtx_sumiso->Write();

   nvtx3_eciso->Write();
   nvtx3_hciso->Write();
   nvtx3_tkiso->Write();
   nvtx3_sumiso->Write();

   nvtx4_eciso->Write();
   nvtx4_hciso->Write();
   nvtx4_tkiso->Write();
   nvtx4_sumiso->Write();

   npix_eciso->Write();
   npix_hciso->Write();
   npix_tkiso->Write();
   npix_sumiso->Write();

   outFile->Close();

}
