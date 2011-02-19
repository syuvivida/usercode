#define LoopEntries_cxx
#include "LoopEntries.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>

void LoopEntries::Loop(std::string outputFile)
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
//    const double fBinsPt[]={25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 100, 120, 200,300};
 const double fBinsPt[]={25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 100, 120, 200,300,400,500,600,700,800,900,1000};
   const int nPtBin = sizeof(fBinsPt)/sizeof(fBinsPt[0])-1;

   const double fBinsEta[]={0.0,0.9,1.5,2.1,2.5};
   const int nEtaBin = sizeof(fBinsEta)/sizeof(fBinsEta[0])-1;

   TH2D* hetapt = new TH2D("hetapt","Number of entries in pt and eta bins",
			   nEtaBin,fBinsEta,nPtBin,fBinsPt);

   ofstream fout;
   fout.open(outputFile.data());
   cout << "writing out data file " << outputFile.data() << endl;

   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries ;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (Cut(ientry) < 0) continue;
      //     if (jentry%1000==0)
      //	printf("%4.1f%% started.\r",(float)jentry/(float)nentries*100.);		

      // good vertex selections
      
      if(nVtxNotFake<1)continue;

      bool hasGoodVtx = false;
 
      for(int ivtx=0; ivtx<1; ivtx++)
	{
	  if(vtxNdof[ivtx] > 4 && 
	     fabs(vtxZ[ivtx])<=24.0 && 
	     sqrt(vtxX[ivtx]*vtxX[ivtx]+vtxY[ivtx]*vtxY[ivtx])<=2.0)
	    {
	      hasGoodVtx=true;
	      break;
	    }
	}


      if(!hasGoodVtx)continue;

//      if(nVtxGood <1)continue;
//      if(!notScraping)continue;
      
//       if (jentry%100==0)
// 	printf("GoodVtx %4.1f%% done.\r",(float)jentry/(float)nentries*100.);		

      for(int ipho=0; ipho < nPhotons; ipho++){
	
	bool isBarrel = fabs(scEta[ipho]) < 1.4442;
	bool isEndcap = fabs(scEta[ipho]) < 2.50 && fabs(scEta[ipho]) > 1.566;
	

	if(!isBarrel && !isEndcap)continue;


// 	if (jentry%100==0)
// 	  printf("EBEE %4.1f%% done.\r",(float)jentry/(float)nentries*100.);		

	// remove spikes
	bool removeSpike = (isBarrel 
			    && seedSeverity[ipho]!=3
			    && sigmaIetaIeta[ipho] > 0.001 && 
			    sqrt(covPhiPhi[ipho]) >0.001);
// 	if (jentry%100==0)
// 	  printf("IDCut1 %4.1f%% done.\r",(float)jentry/(float)nentries*100.);		

	// trigger requirements
	//	if(pt[ipho] >= 25 && pt[ipho] < 35 
	//	   && !HLT_Photon20_Cleaned_L1R)continue;
	//	if(pt[ipho] >= 35 && pt[ipho] < 55
	//	   && !(HLT_Photon20_Cleaned_L1R || HLT_Photon30_Cleaned_L1R))continue;
	//	if(pt[ipho] >=55 && pt[ipho] < 80 
	//	   && !(HLT_Photon20_Cleaned_L1R || HLT_Photon30_Cleaned_L1R || HLT_Photon50_Cleaned_L1R_v1))continue;
	//	if(pt[ipho] >=80 
	//	   && !(HLT_Photon20_Cleaned_L1R || HLT_Photon30_Cleaned_L1R || HLT_Photon50_Cleaned_L1R_v1 || HLT_Photon70_Cleaned_L1R_v1))continue;

	float ptScaled = pt[ipho];
	if(isBarrel)ptScaled *= 1.005*0.9957;
	else ptScaled *= 1.014*0.996;


 	bool trig20 = HLT_Photon20_Cleaned_L1R  && run>=138564 && run<=143962;
 	bool trig30 = HLT_Photon30_Cleaned_L1R  && run>=144010 && run<=147116;
 	bool trig50 = HLT_Photon50_Cleaned_L1R_v1 && run>=147196 && run<=148058;
 	bool trig70 = HLT_Photon70_Cleaned_L1R_v1 && run>=148822 && run<=149294;

 	if(ptScaled >= 25 && ptScaled < 35 
 	   && !trig20)continue;
 	if(ptScaled >= 35 && ptScaled < 55
 	   && !(trig20 || trig30))continue;
 	if(ptScaled >=55 && ptScaled < 80 
 	   && !(trig20 || trig30 || trig50))continue;
 	if(ptScaled >=80 
 	   && !(trig20 || trig30 || trig50 || trig70))continue;
	
	// basic ID cuts
	bool IDCut = (isBarrel && removeSpike &&
		      hadronicOverEm[ipho] <  0.05 && 
		      sigmaIetaIeta[ipho] < 0.01 && 
		      !hasPixelSeed[ipho]) ||
	  (isEndcap && 
	   hadronicOverEm[ipho] < 0.05 && 
	   sigmaIetaIeta[ipho] < 0.028 && 
	   !hasPixelSeed[ipho]);


	if(!IDCut)continue;
// 	if (jentry%100==0)
// 	  printf("IDCut %4.1f%% done.\r",(float)jentry/(float)nentries*100.); 

	Float_t totalIso = 
	  ecalRecHitSumEtConeDR04[ipho]+
	  hcalTowerSumEtConeDR04[ipho]+
	  trkSumPtHollowConeDR04[ipho];

// 	if (jentry%100==0)
// 	  printf("IDCut3 %4.1f%% done.\r",(float)jentry/(float)nentries*100.);		



	if(totalIso >= 20.0)continue;

	hetapt->Fill(fabs(scEta[ipho]),ptScaled);

	if(ptScaled < 25.0 || ptScaled > 1000.0)continue;
  	fout << totalIso << " " << ptScaled << " " << scEta[ipho] << endl;

      } // end of looping over photons

      //      if (jentry%100==0)
      //	printf("Final %4.1f%% done.\r",(float)jentry/(float)nentries*100.);		
   } // end of loop over entries

   fout.close();


   // dump number of entries
   for(int i=0; i < nEtaBin; i++){
     for(int j=0; j < nPtBin; j++){
       cout << "eta " << i << " pt bin " << (int)(fBinsPt[j])
	    << "= " << hetapt->GetBinContent(i+1,j+1) << endl;
     }
   }

}
