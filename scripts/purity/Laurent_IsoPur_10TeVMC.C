#define Laurent_IsoPur_cxx
#include "Laurent_IsoPur.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TH2F.h>
#include <TLine.h>
#include <iostream> 
#include <fstream> 
#include <TFile.h> 
#include <TGraphAsymmErrors.h>
#include <TFractionFitter.h>
#include <TMath.h> 

void Laurent_IsoPur::Loop(Int_t choiceEB) // 0 all  1 Barrel  2 Transition  3 Endcap
{
  
  // constants, these may vary depending on the samples we use for weighting or fitting
  const Int_t nSignalSamples     = 7;
  const Int_t nBackgroundSamples = 7;
  const Int_t nOriginalSamples   = nSignalSamples + nBackgroundSamples;
  // 6 additional files for adding signal data, background data, total data, signal template, background template
  const Int_t nSamples = nOriginalSamples*2 + 6; 

  const Int_t SignalDataSample_start = 0;
  const Int_t SignalDataSample_end = SignalDataSample_start + nSignalSamples-1;
  const Int_t BackgroundDataSample_start = SignalDataSample_end + 1;
  const Int_t BackgroundDataSample_end = BackgroundDataSample_start + nBackgroundSamples-1;

  const Int_t SDataIndex = BackgroundDataSample_end+1;
  const Int_t BDataIndex = BackgroundDataSample_end+2;
  const Int_t SandBDataIndex = BackgroundDataSample_end+3;

  const Int_t SignalTemplateSample_start = SandBDataIndex + 1;
  const Int_t SignalTemplateSample_end = SignalTemplateSample_start + nSignalSamples-1;
  const Int_t BackgroundTemplateSample_start = SignalTemplateSample_end + 1;
  const Int_t BackgroundTemplateSample_end = BackgroundTemplateSample_start + nBackgroundSamples-1;

  const Int_t STemplateIndex = BackgroundTemplateSample_end + 1;
  const Int_t BTemplateIndex = BackgroundTemplateSample_end + 2;
  const Int_t SandBTemplateIndex = BackgroundTemplateSample_end + 3;

  cout << "Printing all constant parameters" << endl;
  cout <<   " nSignalSamples                 = " << nSignalSamples << endl;
  cout <<   " nBackgroundSamples             = " << nBackgroundSamples << endl;
  cout <<   " nOriginalSamples               = " << nOriginalSamples << endl;
  cout <<   " nSamples                       = " << nSamples << endl;
  cout <<   " SignalDataSample_start         = " << SignalDataSample_start << endl;
  cout <<   " SignalDataSample_end           = " << SignalDataSample_end << endl;
  cout <<   " BackgroundDataSample_start     = " << BackgroundDataSample_start << endl;
  cout <<   " BackgroundDataSample_end       = " << BackgroundDataSample_end << endl;
  cout <<   " SDataIndex                     = " << SDataIndex << endl;
  cout <<   " BDataIndex                     = " << BDataIndex << endl;
  cout <<   " SandBDataIndex                 = " << SandBDataIndex << endl;
  cout <<   " SignalTemplateSample_start     = " << SignalTemplateSample_start << endl; 
  cout <<   " SignalTemplateSample_end       = " << SignalTemplateSample_end << endl;                             
  cout <<   " BackgroundTemplateSample_start = " << BackgroundTemplateSample_start << endl;
  cout <<   " BackgroundTemplateSample_end   = " << BackgroundTemplateSample_end << endl;
  cout <<   " STemplateIndex                 = " << STemplateIndex << endl;   
  cout <<   " BTemplateIndex                 = " << BTemplateIndex << endl;   
  cout <<   " SandBTemplateIndex             = " << SandBTemplateIndex << endl;


  const Int_t nPtBin = 14;

  Int_t maxLoop = 5e8;
  Double_t reachS = 0.98; // How much signal do we want to keep?
  Double_t intLum = 100.; // pb-1
  Int_t nBins = 60;
  Int_t reachBins = 300;
  Int_t smooD = 0;
  Int_t smooT = 6;
  Double_t eMin = -0.1;     Double_t eMax = 0.5;
  Double_t hMin = eMin;     Double_t hMax = 4;
  Double_t tMin = eMin;     Double_t tMax = 15;

  // scaled isolation histogram
  Double_t sMin = eMin;     Double_t sMaxLow = 3.5 ; Double_t sMaxHigh = 18.5;
  Double_t sMax[nPtBin]; sMax[0] = 18;
  for(Int_t ptb=0; ptb<nPtBin; ptb++) {
    sMax[ptb] = sMaxLow + ptb*ptb/(13.*13.)*sMaxHigh;
  }
  
  
  // Defining histograms
  
  char nameEB[4][100] = {"allEta","centr","trans","endcp"};
  char printEB[4][100] = {"All |#eta|<2.5","Central |#eta|<.9","Transition .9<|#eta|<1.45","Endcap 1.55<|#eta|<2.5"};
  
  Double_t nFillA[nSamples]; // for the weights
  Double_t nFill_ForIsoSum[nSamples]; // for counting in Laurent's isolation sum
  Double_t nFill[nSamples][nPtBin]; // for the weighted isolation sum
  Double_t realFill[nSamples][nPtBin];
  
  // distributions of ecal, hcal, trackIso
  TH1F* e[nSamples]; 
  TH1F* h[nSamples];
  TH1F* t[nSamples]; 
  
  // Integral(1,x)-effS minimum indicates where effS % of the signal are
  TH1F* eEff; 
  TH1F* hEff; 
  TH1F* tEff;


  //  eff  % of the signal and background
  TH1F* eEff_sig;
  TH1F* hEff_sig;
  TH1F* tEff_sig;
  TH1F* eEff_bkg;
  TH1F* hEff_bkg;
  TH1F* tEff_bkg;

  TH1F* eEff_diff;
  TH1F* hEff_diff;
  TH1F* tEff_diff;

   
  // to check the smoothness of  distributions
  TH1F* ptHatDis[nSamples]; 
    
  // the weighted sum 
  TH1F* iso[nSamples][nPtBin]; 
  TH1F* sumEffS;
  TH1F* sumEffB;
  TH1F* sumEffSB_diff;
  TH1* outMCpredS[nPtBin];
  TH1* outMCpredB[nPtBin];

  // other stuff
  gStyle->SetPadLeftMargin(0.06);
  TString sochn1;
  TString sochnB; 
  char sochn2[100];
  Double_t yield[3][nPtBin];
  
  TLatex histLabel;
  histLabel.SetNDC();
  
  Double_t CCpurity[nPtBin];
  Double_t MCpurity[nPtBin];
  Double_t FFpurity[nPtBin];
  Double_t cutPurity[nPtBin];
  
  Long64_t entries[nSamples];
  Long64_t fileStart[nSamples];
  TFile* indepFile[nSamples];
  TTree* indepTree[nSamples];
  Double_t sc[nSamples];  // is the coefficient applied to the REAL nFill    xSec*skimEff*intLum/nAna 
  
  TGraphAsymmErrors* purMC = new TGraphAsymmErrors(nPtBin);
  TGraphAsymmErrors* purCC = new TGraphAsymmErrors(nPtBin);
  TGraphAsymmErrors* purFF = new TGraphAsymmErrors(nPtBin);
  
  TGraphAsymmErrors* errMC = new TGraphAsymmErrors(nPtBin);
  TGraphAsymmErrors* errCC = new TGraphAsymmErrors(nPtBin);
  TGraphAsymmErrors* errFF = new TGraphAsymmErrors(nPtBin);
  
  TGraphAsymmErrors* nbMC = new TGraphAsymmErrors(nPtBin);
  TGraphAsymmErrors* nbCC = new TGraphAsymmErrors(nPtBin);
  TGraphAsymmErrors* nbFF = new TGraphAsymmErrors(nPtBin);
  
  TGraphAsymmErrors* frame = new TGraphAsymmErrors(1);
  
  TFile* graphs = new TFile("/mc/QCD_31x/scripts/isoTempOutput.root","UPDATE");
    
  // Cross sections
  
  Double_t prexSec[nSamples] = { // in pb
    2.887e+05,  // 0 S15
    3.222e+04,  // 1 S30
    1.010e+03,  // 2 S80
    5.143e+01,  // 3 S170
    4.193e+00,  // 4 S300
    4.515e-01,  // 5 S470
    2.003e-02,  // 6 S800
    
    1.458e+09,  // 7  B15
    1.090e+08,  // 8  B30
    1.936e+06,  // 9  B80
    6.251e+04,  // 10 B170
    3.669e+03,  // 11 B300
    3.153e+02,  // 12 B470
    1.194e+01,  // 13 B800

    // weighted sum
    0,0,0       // 30,31,32
  };    
  
  Double_t xSec[nSamples];
  
  for(Int_t i=SignalDataSample_start;i<SignalDataSample_end;i++)   
    {xSec[i] = prexSec[i]-prexSec[i+1];}
  xSec[SignalDataSample_end] = prexSec[SignalDataSample_end];
  for(Int_t i=BackgroundDataSample_start;i<BackgroundDataSample_end;i++)  
    {xSec[i] = prexSec[i]-prexSec[i+1];}
  xSec[BackgroundDataSample_end] = prexSec[BackgroundDataSample_end];
  xSec[14]=xSec[15]=xSec[16]=0.0;
  for(Int_t i=SignalTemplateSample_start;
      i<=BackgroundTemplateSample_end;i++) {xSec[i] = xSec[i-17];}


  for(Int_t i=0;i<nSamples; i++)
    cout << "xSec[" << i << "]=" << xSec[i] << endl;

  // Skim and Crab Efficiencies
  
  Double_t skimEff[nSamples]={ // skimming efficiency
    .58,  // 0 S15
    .71,  // 1 S30
    .89,  // 2 S80
    .97,  // 3 S170
    .99,  // 4 S300
    .99,  // 5 S470
    .99,  // 6 S800
    
    .019,  // 7  B15
    .32,  // 8  B30
    .66,  // 9  B80
    .81,  // 10 B170
    .91,  // 11 B300
    .96,  // 12 B470
    .97,  // 13 B800
    0,0,0 // 14,15,16
  };    
  
  for(Int_t i=SignalTemplateSample_start;i<nSamples;i++)skimEff[i] = skimEff[i-SignalTemplateSample_start];
  for(Int_t i=0;i<nSamples; i++)
    cout << "skimEff[" << i << "]=" << skimEff[i] << endl;
  
  Double_t nGen[nSamples] = { // generated events
    600533,  // 0 S15
    648508,
    434327,
    1329973,
    996616,
    1009019,
    1235245,
    
    195358,
    2530689,
    1175392,
    1756969,
    1667122,
    1997663,
    1851280,
    
    0,0,0 // 14,15,16

  };    
  

  for(Int_t i=SignalTemplateSample_start;i<nSamples;i++)nGen[i] = nGen[i-SignalTemplateSample_start];
  for(Int_t i=0;i<nSamples; i++)
    cout << "nGen[" << i << "]=" << nGen[i] << endl;
  
  char fileName[nSamples][1000] = {
    "/mc/Laurent/PJ_15_Data.root",
    "/mc/Laurent/PJ_30_Data.root",
    "/mc/Laurent/PJ_80_Data.root",
    "/mc/Laurent/PJ_170_Data.root",
    "/mc/Laurent/PJ_300_Data.root",
    "/mc/Laurent/PJ_470_Data.root",
    "/mc/Laurent/PJ_800_Data.root",
    
    "/mc/Laurent/QCD_15_Data.root",
    "/mc/Laurent/QCD_30_Data.root",
    "/mc/Laurent/QCD_80_Data.root",
    "/mc/Laurent/QCD_170_Data.root",
    "/mc/Laurent/QCD_300_Data.root",
    "/mc/Laurent/QCD_470_Data.root",
    "/mc/Laurent/QCD_800_Data.root",
    
    "","","",
    
    "/mc/Laurent/PJ_15_Template.root",
    "/mc/Laurent/PJ_30_Template.root",
    "/mc/Laurent/PJ_80_Template.root",
    "/mc/Laurent/PJ_170_Template.root",
    "/mc/Laurent/PJ_300_Template.root",
    "/mc/Laurent/PJ_470_Template.root",
    "/mc/Laurent/PJ_800_Template.root",
    
    "/mc/Laurent/QCD_15_Template.root",
    "/mc/Laurent/QCD_30_Template.root",
    "/mc/Laurent/QCD_80_Template.root",
    "/mc/Laurent/QCD_170_Template.root",
    "/mc/Laurent/QCD_300_Template.root",
    "/mc/Laurent/QCD_470_Template.root",
    "/mc/Laurent/QCD_800_Template.root",
    
    "","",""
  };
  
  char fileLbl[nSamples][200] = {
    "S1530_Data", "S3080_Data", "S80170_Data", "S170300_Data", "S300470_Data", "S470800_Data","S800INF_Data",
    "B1530_Data", "B3080_Data", "B80170_Data", "B170300_Data", "B300470_Data", "B470800_Data","B800INF_Data",
    "S_Data", "B_Data","SandB_Data",
    "S1530_Template", "S3080_Template", "S80170_Template", "S170800_Template", "S300470_Template", "S470800_Template", "S800INF_Template",
    "B1530_Template", "B3080_Template", "B80170_Template", "B170300_Template", "B300470_Template", "B470800_Template", "B800INF_Template",
    "S_Template", "B_Template"
  };
  
  Double_t minPtHat[40] = {
    15,30,80,170,300,470,800,
    15,30,80,170,300,470,800,
    0,0,0
  };
  for(Int_t i=17;i<31;i++) {minPtHat[i] = minPtHat[i-17];}
  
  Double_t maxPtHat[40] ;
  
  for(Int_t i=0; i<6;  i++) {maxPtHat[i]=minPtHat[i+1];}
  for(Int_t i=7; i<13; i++) {maxPtHat[i]=minPtHat[i+1];}
  maxPtHat[6] = maxPtHat[13] = 1e5;
  for(Int_t i=17;i<31;i++)  {maxPtHat[i] = maxPtHat[i-17];}

 Double_t GamEt[nPtBin+1]= { 15,20,27,35,45,57,72,90,120,150,200,300,400,550,
			   1000}; 
 Double_t minGamEt[nPtBin]; 
 Double_t maxGamEt[nPtBin];
 Double_t gamEtBinCenter[nPtBin];
 for(Int_t i=0; i< nPtBin; i++){ 
   minGamEt[i]=GamEt[i]; 
   maxGamEt[i]=GamEt[i+1];
   gamEtBinCenter[i]=(minGamEt[i]+maxGamEt[i])/2.;
 }
 
  
  
  cout << endl;
  cout << "Weights for samples " << endl;
  
  cout << endl;
  for(Int_t samp=0; samp<nSamples; samp++) {
    entries[samp] = 0;
    cout << "samp = " << samp << endl;
    if (TString(fileName[samp])!=""){
      indepFile[samp] = new TFile(fileName[samp]);
      indepTree[samp] = (TTree*) indepFile[samp]->Get("photonTree");
      entries[samp] = indepTree[samp]->GetEntries();
      sprintf(sochn2,"entries[%2i] %lld    %s",samp,entries[samp],fileName[samp]);
      cout << sochn2 << endl;
      fileStart[samp] = 0;
    }
  }
    
  for(Int_t samp=0; samp<nSamples; samp++) {
    for(Int_t samp2=samp+1; samp2<nSamples; samp2++) {
      fileStart[samp2] += entries[samp];
    }
  }
  
  cout << endl;
  for(Int_t samp=0; samp<nSamples; samp++) {
    sprintf(sochn2,"fileStart[%2i] %lld   %s",samp,fileStart[samp],fileName[samp]);
    if (entries[samp]!=0) cout << sochn2 << endl;
  }
  
  // HISTOGRAMS
  
  eEff = new TH1F("eEff","eEff",reachBins,eMin,eMax);
  hEff = new TH1F("hEff","hEff",reachBins,hMin,hMax);
  tEff = new TH1F("tEff","tEff",reachBins,tMin,tMax);

  eEff_sig = new TH1F("eEff_sig","eEff_sig",reachBins,eMin,eMax);
  hEff_sig = new TH1F("hEff_sig","hEff_sig",reachBins,hMin,hMax);
  tEff_sig = new TH1F("tEff_sig","tEff_sig",reachBins,tMin,tMax);

  eEff_bkg = new TH1F("eEff_bkg","eEff_bkg",reachBins,eMin,eMax);
  hEff_bkg = new TH1F("hEff_bkg","hEff_bkg",reachBins,hMin,hMax);
  tEff_bkg = new TH1F("tEff_bkg","tEff_bkg",reachBins,tMin,tMax);

  eEff_diff = new TH1F("eEff_diff","eEff_diff",reachBins,eMin,eMax);
  hEff_diff = new TH1F("hEff_diff","hEff_diff",reachBins,hMin,hMax);
  tEff_diff = new TH1F("tEff_diff","tEff_diff",reachBins,tMin,tMax);
  
  sumEffS = new TH1F("sumEffS","sumEffS",nBins,sMin,sMaxLow);
  sumEffB = new TH1F("sumEffB","sumEffB",nBins,sMin,sMaxLow);
  sumEffSB_diff = new TH1F("sumEffSB_diff","sumEffSB_diff",nBins,sMin,sMaxLow);
  
  for (Int_t i=0; i<nSamples; i++){
    for (Int_t pt=0; pt<nPtBin; pt++) {
      nFill[i][pt] = 0;
      realFill[i][pt] = 0;
      sprintf(sochn2,"iso_%s_%i",fileLbl[i],pt);
      iso[i][pt] = new TH1F(sochn2,sochn2,nBins,sMin,sMax[pt]);
      iso[i][pt]->Sumw2();
    }
    
    sprintf(sochn2,"ptHatDis_%s",fileLbl[i]);
    ptHatDis[i] = new TH1F(sochn2,sochn2,1000,0,1000);

    sprintf(sochn2,"e_%s",fileLbl[i]);
    e[i] = new TH1F(sochn2,sochn2,reachBins,eMin,eMax);
    e[i] -> Sumw2();

    sprintf(sochn2,"h_%s",fileLbl[i]);
    h[i] = new TH1F(sochn2,sochn2,reachBins,hMin,hMax);
    h[i] -> Sumw2();
    
    sprintf(sochn2,"t_%s",fileLbl[i]);
    t[i] = new TH1F(sochn2,sochn2,reachBins,tMin,tMax);
    t[i] -> Sumw2();
  
    nFill_ForIsoSum[i] = 0;
    nFillA[i] = 0;
  }

  // FILLING WEIGHTS
  cout << endl;
  cout << "Filling weight histograms." << endl;


  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
   
  Long64_t nbytes = 0, nb = 0;
  
  // first loop over samples, and check the point of 98% isolation scaling
  for (Int_t samp=SignalDataSample_start; samp<=BackgroundDataSample_end; samp++){
    for (Long64_t jentry=fileStart[samp]; jentry<fileStart[samp]+maxLoop && jentry<fileStart[samp]+entries[samp] ;jentry++) {
      Long64_t ientry = LoadTree(jentry); if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      sochn1 = TString(fChain->GetCurrentFile()->GetName()); // file name
      if (sochn1 != sochnB) { cout << " >> looping on " << sochn1 << " ... " << endl; }
      sochnB = TString(fChain->GetCurrentFile()->GetName());


      if(ptHat < minPtHat[samp]) continue;
      if(ptHat > maxPtHat[samp]) continue;

      ptHatDis[samp]->Fill(ptHat);
      
      // check if this is a QCD sample
      Bool_t isQCD = false;
      std::string inputFile_(fileName[samp]);
      if(inputFile_.find("DiJet") != std::string::npos)isQCD=true;
      
      Bool_t findLeadingPhoton=false;
      // first loose photon ID cuts
      for (int ipho=0; ipho < nPhotons; ipho++){

	if(!IsLoosePhoton(isQCD,ipho))continue;
	if(findLeadingPhoton)continue; // only use leading photon
	findLeadingPhoton = true;
	
	if(!CheckFiducial(choiceEB,ipho))continue;
	Double_t thiset = et[ipho];	
 	if(thiset        < GamEt[0])continue;
 	if(thiset        >= GamEt[nPtBin])continue;
        e[samp]->Fill(ecalRecHitSumEtConeDR04[ipho]/et[ipho]);
        h[samp]->Fill(hcalTowerSumEtConeDR04[ipho]);
        t[samp]->Fill(trkSumPtHollowConeDR04[ipho]);
	nFill_ForIsoSum[samp] += 1.0;

      } // loop over photons
    } // loop over entries
  } // loop over samples

  
  cout << "Weight histograms filled." << endl;
  cout << endl;
  
  // Define the scaling used througcout the sochn
  cout << endl;
  for(Int_t i=0; i<nSamples; i++) {
    if (nGen[i]>0) 
      sc[i]    = intLum * xSec[i]    * skimEff[i] / nGen[i];
    else 
      sc[i] = 0;
    
  }

  for(Int_t i=0; i<nSamples; i++){
    if(sc[i]) cout << "sc[" << i << "] : " << sc[i] << endl;
  }
  
  cout << endl; // first print nFillA
  // then print nAna
  for (Int_t i=0;i<nSamples;i++) {
    if(nFill_ForIsoSum[i]) cout << "nFill_ForIsoSum[" << i << "] : " << nFill_ForIsoSum[i] << endl;
    ptHatDis[i]->Scale(sc[i]);
  }
  
  Double_t weightNormS = 0;
  Double_t weightNormB = 0;
  // now get the sum of weighted histogram
  // signal data
  e[SDataIndex]->Reset();
  h[SDataIndex]->Reset();
  t[SDataIndex]->Reset();

  e[BDataIndex]->Reset();
  h[BDataIndex]->Reset();
  t[BDataIndex]->Reset();

  for(Int_t i=SignalDataSample_start; i<=SignalDataSample_end; i++){
    Double_t scaleFactor = sc[i];
    e[SDataIndex]->Add(e[i],scaleFactor);
    h[SDataIndex]->Add(h[i],scaleFactor);
    t[SDataIndex]->Add(t[i],scaleFactor);
    weightNormS  += scaleFactor*nFill_ForIsoSum[i];
  }

  cout << "weightNormS = " << weightNormS << endl;

  // background data
  for(Int_t i=BackgroundDataSample_start; i<=BackgroundDataSample_end; i++){
    Double_t scaleFactor = sc[i];
    e[BDataIndex]->Add(e[i],scaleFactor);
    h[BDataIndex]->Add(h[i],scaleFactor);
    t[BDataIndex]->Add(t[i],scaleFactor);
    weightNormB += scaleFactor*nFill_ForIsoSum[i];
  }


  cout << "weightNormB = " << weightNormB << endl;

  e[SDataIndex]->Scale(1./weightNormS);
  h[SDataIndex]->Scale(1./weightNormS);
  t[SDataIndex]->Scale(1./weightNormS);
 
  e[BDataIndex]->Scale(1./weightNormB);
  h[BDataIndex]->Scale(1./weightNormB);
  t[BDataIndex]->Scale(1./weightNormB);
  
  cout << "e[SDataIndex] integral = " << e[SDataIndex]->Integral() << endl;
  cout << "h[SDataIndex] integral = " << h[SDataIndex]->Integral() << endl;
  cout << "t[SDataIndex] integral = " << t[SDataIndex]->Integral() << endl;


  cout << "e[BDataIndex] integral = " << e[BDataIndex]->Integral() << endl;
  cout << "h[BDataIndex] integral = " << h[BDataIndex]->Integral() << endl;
  cout << "t[BDataIndex] integral = " << t[BDataIndex]->Integral() << endl;

 
  // EFFICIENCY HISTOGRAMS
   
  for (Int_t bin=1; bin<reachBins+1; bin++) {
    eEff->SetBinContent(bin,fabs((Double_t)e[SDataIndex]->Integral(1,bin)-reachS));
    hEff->SetBinContent(bin,fabs((Double_t)h[SDataIndex]->Integral(1,bin)-reachS));
    tEff->SetBinContent(bin,fabs((Double_t)t[SDataIndex]->Integral(1,bin)-reachS));

    eEff_sig->SetBinContent(bin,(Double_t)e[SDataIndex]->Integral(1,bin));
    hEff_sig->SetBinContent(bin,(Double_t)h[SDataIndex]->Integral(1,bin));
    tEff_sig->SetBinContent(bin,(Double_t)t[SDataIndex]->Integral(1,bin));

    eEff_bkg->SetBinContent(bin,(Double_t)e[BDataIndex]->Integral(1,bin));
    hEff_bkg->SetBinContent(bin,(Double_t)h[BDataIndex]->Integral(1,bin));
    tEff_bkg->SetBinContent(bin,(Double_t)t[BDataIndex]->Integral(1,bin));


    eEff_diff->SetBinContent(bin,(Double_t)(e[SDataIndex]->Integral(1,bin)-e[BDataIndex]->Integral(1,bin)));
    hEff_diff->SetBinContent(bin,(Double_t)(h[SDataIndex]->Integral(1,bin)-h[BDataIndex]->Integral(1,bin)));
    tEff_diff->SetBinContent(bin,(Double_t)(t[SDataIndex]->Integral(1,bin)-t[BDataIndex]->Integral(1,bin)));

  } 
  
  // Find out the maximum difference in signal and background efficiency
  cout << endl;
  cout << "--- Maximum difference between signal and background --- " << endl;
  
  cout << "eCalIso/pt : " << eEff_diff->GetBinCenter(
						     eEff_diff->GetMaximumBin()) << endl;
  cout << "hcalIso : " << hEff_diff->GetBinCenter(
						  hEff_diff->GetMaximumBin()) << endl;
  cout << "trkIso : " << tEff_diff->GetBinCenter(
						 tEff_diff->GetMaximumBin()) << endl;


  // DETERMINING THE WEIGHTS
  
  cout << endl;
  cout << "-- Scaling & Efficiencies " << endl;
  Double_t eReach = (eMax-eMin)*eEff->GetMinimumBin()/reachBins + eMin;
  cout << "Cut on ecal/et signal with " << reachS*100 << "% efficiency : " << eReach << endl;
  cout << "Efficiency for ecal background : " << 100*e[BDataIndex]->Integral(1,eEff->GetMinimumBin())<< endl;
  Double_t hReach = (hMax-hMin)*hEff->GetMinimumBin()/reachBins + hMin;
  cout << "Cut on hcal signal with " << reachS*100 << "% efficiency : " << hReach << endl;
  cout << "Efficiency for hcal background : " << 100*h[BDataIndex]->Integral(1,hEff->GetMinimumBin())<< endl;
  Double_t tReach = (tMax-tMin)*tEff->GetMinimumBin()/reachBins + tMin;
  cout << "Cut on trck signal with " << reachS*100 << "% efficiency : " << tReach << endl;
  cout << "Efficiency for trck background : " << 100*t[BDataIndex]->Integral(1,tEff->GetMinimumBin())<< endl;
  cout << endl;
   
  // FILLING THE SUMS
  
  //  Double_t noSpikeCut[nSamples];
  //  for(Int_t i =0; i < nSamples; i++) noSpikeCut[i]=maxPtHat[i];
  Double_t noSpikeCut[nSamples] = {
    30.,80.,170.,300.,470.,800.,1e5,
    20.,45., 90.,300.,400.,550.,1e5,
    0,0,0,
    30.,80.,170.,300.,470.,800.,1e5,
    20.,45., 90.,300.,400.,550.,1e5,
    0,0
  };


  Double_t antiSpike[nSamples][nPtBin];
  for(Int_t i=0;i<nSamples;i++) {for(Int_t j=0; j<nPtBin; j++) {antiSpike[i][j]=0;}}
  
  for(Int_t i=0;i<nSamples;i++) cout << "noSpikeCut[" << i << "]  et[ipho] < " << noSpikeCut[i] << endl;
  
  
  cout << endl;
  cout << "Filling distribution histograms." << endl;


  if (fChain == 0) return;
  nbytes = 0; nb = 0;


  // NOW the main jobs to fill scaled isolation sums
  
  for (Int_t samp=0; samp< nSamples; samp++){
    for (Long64_t jentry=fileStart[samp]; jentry<fileStart[samp]+maxLoop && jentry<fileStart[samp]+entries[samp] ;jentry++) {
      Long64_t ientry = LoadTree(jentry); if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      sochn1 = TString(fChain->GetCurrentFile()->GetName()); // file name
      if (sochn1 != sochnB) { cout << " >> looping on " << sochn1 << " ... " << endl; }
      sochnB = TString(fChain->GetCurrentFile()->GetName());
       // check pthat

      // check if this is a QCD sample
      Bool_t isQCD = false;
      std::string inputFile_(fileName[samp]);
      if(inputFile_.find("DiJet") != std::string::npos)isQCD=true;

      if(ptHat < minPtHat[samp]) continue;
      if(ptHat > maxPtHat[samp]) continue;

      
      Bool_t findLeadingPhoton = false;      
      for (int ipho=0; ipho < nPhotons; ipho++){

	
 	if(!IsLoosePhoton(isQCD,ipho))continue;	
 	if(findLeadingPhoton)continue;
 	findLeadingPhoton = true;
	
	if(!CheckFiducial(choiceEB,ipho))continue;
	  
	Double_t thiset = et[ipho];	
	if(thiset        < GamEt[0])continue;
	if(thiset        >= GamEt[nPtBin])continue;
	int EtBin = -1;
	EtBin = TMath::BinarySearch(nPtBin+1, GamEt,thiset);

	if(thiset        >noSpikeCut[samp]
	   && thiset     <maxPtHat[samp] )
	  antiSpike[samp][EtBin] +=1.0;
        nFillA[samp]             +=1.0;
	// This is to work around peaks and yet compute the right purity
	if(EtBin<0 || EtBin >= nPtBin){cout << "Err" << endl; continue;}
	realFill[samp][EtBin]    +=1.0;

	if(thiset        > noSpikeCut[samp])continue;

	Double_t scaledIsoSum = 
	  (ecalRecHitSumEtConeDR04[ipho]/et[ipho])/eReach + 
	  hcalTowerSumEtConeDR04[ipho]/hReach + 
	  trkSumPtHollowConeDR04[ipho]/tReach;

	iso[samp][EtBin]->Fill(scaledIsoSum);
	nFill[samp][EtBin]+=1.0;
      } // end of loop over photons
    } // end of loop over entries
  } // end of loop over samples


  cout << "Distribution histograms filled." << endl;
  cout << endl;

  for (Int_t i=0;i<nSamples;i++) {
    //    if(nFillA[i]) cout << "nFillA[" << i << "] : " << nFillA[i] << endl;
    for(Int_t j=0; j<nPtBin; j++)
      if(nFill[i][j]) cout << "nFill[" << i << "]["<< j << "] : " 
			   << nFill[i][j] << endl;
  }
  cout << endl; 


  for (Int_t i=0;i<nSamples;i++) {
    bool Print = false;
    for (Int_t j=0;j<nPtBin;j++) {
      if(antiSpike[i][j] && !Print) {
	sprintf(sochn2,"aSpike[%2i] \n  ",i);
	Print = true;
      }
      sprintf(sochn2,"%.3e [%i] ",antiSpike[i][j],j);
      cout << sochn2;
    }
  }



  // Data and Template IsoSum

  
  for (Int_t samp=0; samp<nSamples; samp++) {
    if(nFillA[samp]<=0)continue;
    for (Int_t ptb=0; ptb<nPtBin; ptb++) {
      cout << "sc[" << samp << "]["<< ptb << "] = " << sc[samp] << endl;
      iso[samp][ptb]->Scale(sc[samp]);
      antiSpike[samp][ptb] *= sc[samp];
    }
  }      


  for (Int_t ptb=0;ptb<nPtBin;ptb++) {
    
    weightNormS = weightNormB = 0;
    nFill[SDataIndex][ptb] =0;
    nFill[BDataIndex][ptb] =0;
    nFill[STemplateIndex][ptb] =0;
    nFill[BTemplateIndex][ptb] =0;

    // 1st adding signal data histograms
    for(int i=SignalDataSample_start; i <= SignalDataSample_end; i++)
      {
	iso[SDataIndex][ptb]->Add(iso[i][ptb]);
   	nFill[SDataIndex][ptb] += iso[i][ptb]->Integral() + iso[i][ptb]->GetBinContent(0) + iso[i][ptb]->GetBinContent(nBins+1) + antiSpike[i][ptb];
	if(nFill[i][ptb])    weightNormS += sc[i];
      }

    // 2nd adding background data histograms
    for(int i=BackgroundDataSample_start; i <= BackgroundDataSample_end; i++)
      {
	iso[BDataIndex][ptb]->Add(iso[i][ptb]);
   	nFill[BDataIndex][ptb] += iso[i][ptb]->Integral() + iso[i][ptb]->GetBinContent(0) + iso[i][ptb]->GetBinContent(nBins+1) + antiSpike[i][ptb];
	if(nFill[i][ptb])    weightNormB += sc[i];
      }

    // 3rd adding signal template histograms
    for(int i=SignalTemplateSample_start; i <= SignalTemplateSample_end; i++)
      {
	iso[STemplateIndex][ptb]->Add(iso[i][ptb]);
   	nFill[STemplateIndex][ptb] += iso[i][ptb]->Integral() + iso[i][ptb]->GetBinContent(0) + iso[i][ptb]->GetBinContent(nBins+1) + antiSpike[i][ptb];
      }
    // 4th adding background data histograms
    for(int i=BackgroundTemplateSample_start; i<= BackgroundTemplateSample_end; i++)
      {
	iso[BTemplateIndex][ptb]->Add(iso[i][ptb]);
  	nFill[BTemplateIndex][ptb] += iso[i][ptb]->Integral() + iso[i][ptb]->GetBinContent(0) + iso[i][ptb]->GetBinContent(nBins+1) + antiSpike[i][ptb];
      }

  
    iso[BDataIndex][ptb]->Smooth(smooD);
    iso[SandBDataIndex][ptb]->Add(iso[SDataIndex][ptb],iso[BDataIndex][ptb]);
    nFill[SandBDataIndex][ptb] = nFill[SDataIndex][ptb] + nFill[BDataIndex][ptb];
    realFill[SDataIndex][ptb] = nFill[SDataIndex][ptb]/weightNormS;
    realFill[BDataIndex][ptb] = nFill[BDataIndex][ptb]/weightNormB;
    realFill[SandBDataIndex][ptb] = realFill[SDataIndex][ptb]+realFill[BDataIndex][ptb];

    realFill[STemplateIndex][ptb] = nFill[STemplateIndex][ptb]/weightNormB;
    realFill[BTemplateIndex][ptb] = nFill[BTemplateIndex][ptb]/weightNormB;
      
    iso[BTemplateIndex][ptb]->Smooth(smooT);
  }
  
  for (Int_t i=0;i<nSamples;i++) {
    for(Int_t j=0; j<nPtBin; j++)
      {
	if(nFill[i][j]) cout << "nFill[" << i << "]["<< j << "] : " 
			     << nFill[i][j] << endl;

	if(realFill[i][j]) cout << "realFill[" << i << "]["<< j << "] : " 
				<< realFill[i][j] << endl;

	if(antiSpike[i][j]) cout << "antiSpike[" << i << "]["<< j << "] : " 
				<< antiSpike[i][j] << endl;

	if(iso[i][j]->GetBinContent(0)) cout << "iso[" << i << "]["<< j << "] : " 
			     << iso[i][j]->GetBinContent(0) << endl;

      }
    
  }




  // Drawing the sumIso distribution in etBin 0 = all ets.

  Double_t cut = 1.0;  
  if(1) { 
    sprintf(sochn2,"SUMiso_%s",nameEB[choiceEB]);
    TCanvas* canIso = new TCanvas(sochn2,sochn2,800,800);
    canIso->cd(1)->SetLogy();
    Double_t intSDataIndex = iso[SDataIndex][1]->Integral()+iso[SDataIndex][1]->GetBinContent(nBins+1);
    Double_t intBDataIndex = iso[BDataIndex][1]->Integral()+iso[BDataIndex][1]->GetBinContent(nBins+1);
    
    iso[SDataIndex][1]->Scale(1./intSDataIndex);
    iso[BDataIndex][1]->Scale(1./intBDataIndex);
    
    for (Int_t i=1; i<nBins+1; i++) {
      sumEffS->SetBinContent(i,iso[SDataIndex][1]->GetMaximum()*iso[SDataIndex][1]->Integral(1,i));
      sumEffB->SetBinContent(i,iso[SDataIndex][1]->GetMaximum()*iso[BDataIndex][1]->Integral(1,i));
      
      float diffEff = iso[SDataIndex][1]->Integral(1,i) - 
	iso[BDataIndex][1]->Integral(1,i);
      sumEffSB_diff->SetBinContent(i, diffEff);
    }
    
    iso[SDataIndex][1]->GetXaxis()->SetTitle("SUM ISO");
    iso[SDataIndex][1]->SetLineColor(2);
    iso[BDataIndex][1]->SetLineColor(4);
    sumEffS->SetLineColor(2);
    sumEffS->SetLineStyle(2);
    sumEffB->SetLineColor(4);
    sumEffB->SetLineStyle(2);
    sumEffSB_diff->Scale(iso[SDataIndex][1]->GetMaximum());
    sumEffSB_diff->SetLineColor(3);
    sumEffSB_diff->SetLineStyle(2);
  
    cut = sumEffSB_diff->GetBinCenter(sumEffSB_diff->GetMaximumBin());
    cout << "The maximum difference is at " << cut << endl;
  
    iso[SDataIndex][1]->DrawCopy("hist");
    iso[BDataIndex][1]->DrawCopy("hist same");
    sumEffS->DrawCopy("hist same");
    sumEffB->DrawCopy("hist same");
    sumEffSB_diff->DrawCopy("hist same");
   
    canIso->SaveAs(".gif");
    canIso->SaveAs(".eps");
    
    iso[SDataIndex][1]->Scale(intSDataIndex);
    iso[BDataIndex][1]->Scale(intBDataIndex);
    delete canIso;
  }


    
  // DETERMINING EFFICIENCIES
  
  Int_t cutBin[nPtBin]; 
  
  for (Int_t i=0; i<3; i++) {
    for (Int_t j=0; j<nPtBin; j++) { yield[i][j]=0; }
  }
  
  for (Int_t ptb=0;ptb<nPtBin;ptb++) {

    cutBin[ptb] = TMath::Nint((cut-sMin)/(sMax[ptb]-sMin)*nBins); 
    
    // Efficiencies
    if(nFill[SandBDataIndex]) yield[0][ptb]  = iso[SandBDataIndex][ptb]->Integral(1,cutBin[ptb])/(nFill[SandBDataIndex][ptb]); // eff SnB
    if(nFill[STemplateIndex]) yield[1][ptb]  = iso[STemplateIndex][ptb]->Integral(1,cutBin[ptb])/(nFill[STemplateIndex][ptb]); // eff S
    if(nFill[BTemplateIndex]) yield[2][ptb]  = iso[BTemplateIndex][ptb]->Integral(1,cutBin[ptb])/(nFill[BTemplateIndex][ptb]); // eff B

    cutPurity[ptb] = iso[SDataIndex][ptb]->Integral(1,cutBin[ptb])/iso[SandBDataIndex][ptb]->Integral(1,cutBin[ptb]);
    
    // Purities
    CCpurity[ptb] = 0;    MCpurity[ptb] = 0; 
    
    if (yield[1][ptb]-yield[2][ptb]>1e-3)       CCpurity[ptb] = (yield[0][ptb] - yield[2][ptb])/(yield[1][ptb] - yield[2][ptb]);
    if (nFill[SDataIndex][ptb]+nFill[nPtBin][ptb]>1e-3)     MCpurity[ptb] = nFill[SDataIndex][ptb]/(nFill[SandBDataIndex][ptb]);
  }
  

  for (Int_t j=0; j<nPtBin; j++){
    for (Int_t i=0; i<3; i++) {
      cout << "yield[" << i << "][" << j << "] = " << yield[i][j] << endl; 
    }
    cout << "CCpurity[" << j << "]= " << CCpurity[j] << endl;
    cout << "MCpurity[" << j << "]= " << MCpurity[j] << endl;
  }

  Double_t errorMC[nPtBin];
  Double_t errorCC[nPtBin];
  Double_t errorFF[nPtBin];
  for(Int_t i=0;i<nPtBin;i++) {errorMC[i] = errorCC[i] = errorFF[i] = 0.;}
  
  for(Int_t i=0;i<nPtBin;i++) {
    purMC->SetPoint(i,gamEtBinCenter[i],100*MCpurity[i]);
    purMC->SetPointError(i,gamEtBinCenter[i]-minGamEt[i], TMath::Min(maxGamEt[i],1000.)-gamEtBinCenter[i],0.,0.);
    
    purCC->SetPoint(i,gamEtBinCenter[i],100*CCpurity[i]);
    errorCC[i] = 50*purCCerr(yield[2][i], realFill[BTemplateIndex][i],
                             yield[1][i], realFill[STemplateIndex][i],
                             yield[0][i], realFill[SandBDataIndex][i]);
    purCC->SetPointError(i,gamEtBinCenter[i]-minGamEt[i], TMath::Min(maxGamEt[i],1000.)-gamEtBinCenter[i],
                         TMath::Min(99.*CCpurity[i],errorCC[i]), TMath::Min(100.-100.*CCpurity[i],errorCC[i]));
    
    nbMC->SetPoint(i,gamEtBinCenter[i],nFill[SDataIndex][i]);
    nbMC->SetPointError(i,gamEtBinCenter[i]-minGamEt[i], TMath::Min(maxGamEt[i],1000.)-gamEtBinCenter[i],0.,0.);
    
    nbCC->SetPoint(i,gamEtBinCenter[i],nFill[SandBDataIndex][i]*CCpurity[i]);
    nbCC->SetPointError(i,gamEtBinCenter[i]-minGamEt[i], TMath::Min(maxGamEt[i],1000.)-gamEtBinCenter[i],
                        nFill[SandBDataIndex][i]*errorCC[i]/100.,nFill[SandBDataIndex][i]*errorCC[i]/100.);
  }


  
  frame->SetPoint(0,500.,10.);
  frame->SetPointError(0,500.-12.5,600.,0.,0.);


  // FRACTION FITTING 
  
  TFractionFitter* fitter[nPtBin];
  TObjArray *fitmc[nPtBin];
  Int_t status[nPtBin];
  TH1F* fitresult[nPtBin];
  Double_t preFFpurity[nPtBin];
  Double_t probability[nPtBin];
  
  Int_t doFit = 1;
  
  if(doFit) {
  
  for (Int_t ptb=0; ptb<nPtBin; ptb++) {
    fitmc[ptb] = new TObjArray(2);
    
    cout << "ratio = " << realFill[STemplateIndex][ptb]/nFill[STemplateIndex][ptb] << endl;
    if (nFill[STemplateIndex][ptb]) iso[STemplateIndex][ptb]->Scale(realFill[STemplateIndex][ptb]/nFill[STemplateIndex][ptb]);
    if (nFill[BTemplateIndex][ptb]) iso[BTemplateIndex][ptb]->Scale(realFill[BTemplateIndex][ptb]/nFill[BTemplateIndex][ptb]);
    if(iso[SandBDataIndex][ptb]->Integral()<=0)
      {cout << "iso["<< SandBDataIndex << "][" << ptb << "] integral <0 " << endl;continue;}
    if(iso[STemplateIndex][ptb]->Integral()<=0)
      {cout << "iso["<< STemplateIndex << "][" << ptb << "] integral <0 " << endl;continue;}
    if(iso[BTemplateIndex][ptb]->Integral()<=0)
      {cout << "iso["<< BTemplateIndex << "][" << ptb << "] integral <0 " << endl;continue;}

    fitmc[ptb]->Add(iso[STemplateIndex][ptb]);
    fitmc[ptb]->Add(iso[BTemplateIndex][ptb]);
    
    fitter[ptb] = new TFractionFitter(iso[SandBDataIndex][ptb],fitmc[ptb]);
    fitter[ptb]->Constrain(0,.002,.998);
    fitter[ptb]->Constrain(1,.002,.998);
    status[ptb] = fitter[ptb]->Fit();
    cout << "Pt bin " << ptb << " fit status: " << status[ptb] << endl;
    cout << endl; cout << endl;
    if (status[ptb] == 0) {
      fitresult[ptb] = (TH1F*) fitter[ptb]->GetPlot();
      probability[ptb] = fitter[ptb]->GetProb();
    }
    fitter[ptb]->GetResult(0,preFFpurity[ptb],errorFF[ptb]);
    cout << "preFFpurity[" << ptb << "] = " << preFFpurity[ptb] << " +- "
	 << errorFF[ptb] << endl;

    FFpurity[ptb] = preFFpurity[ptb]*(iso[SDataIndex][ptb]->Integral()+iso[BDataIndex][ptb]->Integral())/(iso[SDataIndex][ptb]->Integral()+nFill[BDataIndex][ptb]);
    errorFF[ptb] *= 50;
    purFF->SetPoint(ptb-1,gamEtBinCenter[ptb],100*FFpurity[ptb]);
    purFF->SetPointError(ptb-1,gamEtBinCenter[ptb]-minGamEt[ptb], TMath::Min(maxGamEt[ptb],1000.)-gamEtBinCenter[ptb],
                         TMath::Min(100.*FFpurity[ptb],errorFF[ptb]), TMath::Min(100.-100.*FFpurity[ptb],errorFF[ptb]));
    if (nFill[STemplateIndex][ptb]) iso[STemplateIndex][ptb]->Scale(nFill[STemplateIndex][ptb]/realFill[STemplateIndex][ptb]);
    if (nFill[BTemplateIndex][ptb]) iso[BTemplateIndex][ptb]->Scale(nFill[BTemplateIndex][ptb]/realFill[BTemplateIndex][ptb]);
    cout << endl;
    
    nbFF->SetPoint(ptb-1,gamEtBinCenter[ptb],nFill[SandBDataIndex][ptb]*FFpurity[ptb]);
    nbFF->SetPointError(ptb-1,gamEtBinCenter[ptb]-minGamEt[ptb], TMath::Min(maxGamEt[ptb],1000.)-gamEtBinCenter[ptb],
                        nFill[SandBDataIndex][ptb]*errorFF[ptb]/100.,nFill[SandBDataIndex][ptb]*errorFF[ptb]/100.);
  }
  cout << endl;
  }
  
  
  cout << endl;
  cout << "-- Isolation Summary" << endl;
  for(Int_t i=0;i<nPtBin;i++) {
    sprintf(sochn2,"Et_pho bin %i [%.1e,%.1e]         cut (%.1f) purity %.4f",i,minGamEt[i],maxGamEt[i],cut,cutPurity[i]);
    cout << sochn2 << endl;
    sprintf(sochn2,"    eff_S     %.4e",yield[1][i]);
    cout << sochn2 << endl;
    sprintf(sochn2,"    eff_B     %.4e",yield[2][i]);
    cout << sochn2 << endl;
    sprintf(sochn2,"    eff_Data  %.4e",yield[0][i]);
    cout << sochn2 << endl;
    sprintf(sochn2,"    purity MC %.4e     pm     %.4e",MCpurity[i],errorMC[i]);
    cout << sochn2 << endl;
    sprintf(sochn2,"    purity CC %.4e     pm     %.4e",CCpurity[i],errorCC[i]);
    cout << sochn2 << endl;
    sprintf(sochn2,"    purity FF %.4e     pm     %.4e",FFpurity[i],errorFF[i]);
    cout << sochn2 << endl;
    
    errMC->SetPoint(i,gamEtBinCenter[i],TMath::Min(errorMC[i],100.));
    errMC->SetPointError(i,gamEtBinCenter[i]-minGamEt[i], TMath::Min(maxGamEt[i],1000.)-gamEtBinCenter[i],0.,0.);
    errCC->SetPoint(i,gamEtBinCenter[i],TMath::Min(errorCC[i],100.));
    errCC->SetPointError(i,gamEtBinCenter[i]-minGamEt[i], TMath::Min(maxGamEt[i],1000.)-gamEtBinCenter[i],0.,0.);
    errFF->SetPoint(i,gamEtBinCenter[i],TMath::Min(errorFF[i],100.));
    errFF->SetPointError(i,gamEtBinCenter[i]-minGamEt[i], TMath::Min(maxGamEt[i],1000.)-gamEtBinCenter[i],0.,0.);
  }
  cout << endl;
  
  for(Int_t i=0; i<nPtBin; i++){
    cout << "prob[" << i << "] : " << probability[i] << endl;
  }
  
  cout << endl;
  for (Int_t i=0; i<nPtBin; i++) {
    sprintf(sochn2,"RealFill [SDataIndex] %7f [BDataIndex] %7f [SandBDataIndex] %7f",realFill[SDataIndex][i],realFill[BDataIndex][i],realFill[SandBDataIndex][i]);
    cout << sochn2 << endl;
  }
  cout << endl;
  for (Int_t i=0; i<nPtBin; i++) {
    sprintf(sochn2,"sMax[ptb] %7f",sMax[i]);
    cout << sochn2 << endl;
  }
  cout << endl;
  
  // DRAWING
    
  // Checking the smoothness of pthat distributions
  
  if (1) { 
    TCanvas* canPtHatDis = new TCanvas("PtHatDisSig","PtHatDisSig",800,800);
    canPtHatDis->SetLogy();
    ptHatDis[SignalDataSample_start]->SetMinimum(1e-10*ptHatDis[SignalDataSample_start]->GetMaximum());
    ptHatDis[SignalDataSample_start]->DrawCopy();
    for(int i=SignalDataSample_start; i <= SignalDataSample_end; i++)
      {
	int COLOR = (i%2==1)? 2: 1;
	ptHatDis[i]->SetLineColor(COLOR);
	ptHatDis[i]->DrawCopy("same");
      }
    canPtHatDis->SaveAs(".gif");
    canPtHatDis->SaveAs(".eps");
    delete canPtHatDis;

    canPtHatDis = new TCanvas("PtHatDisBkg","PtHatDisBkg",800,800);
    canPtHatDis->SetLogy();
    ptHatDis[BackgroundDataSample_start]->SetMinimum(1e-10*ptHatDis[BackgroundDataSample_start]->GetMaximum());
    ptHatDis[BackgroundDataSample_start]->DrawCopy();
    for(int i=BackgroundDataSample_start; i <= BackgroundDataSample_end; i++)
      {
	int COLOR = (i%2==1)? 2: 1;
	ptHatDis[i]->SetLineColor(COLOR);
	ptHatDis[i]->DrawCopy("same");
      }
    canPtHatDis->SaveAs(".gif");
    canPtHatDis->SaveAs(".eps");
    delete canPtHatDis;

  }
  
  // Drawing the initial distributions of e,h,t and efficiencies
  
  e[SDataIndex]->SetLineColor(2);
  e[SDataIndex]->GetXaxis()->SetTitle("ECALISO/ET");
  e[SDataIndex]->GetXaxis()->SetTitleSize(0.07);
  e[BDataIndex]->SetLineColor(4);
  h[SDataIndex]->SetLineColor(2);
  h[SDataIndex]->GetXaxis()->SetTitle("HCALISO");
  h[SDataIndex]->GetXaxis()->SetTitleSize(0.07);
  h[BDataIndex]->SetLineColor(4);
  t[SDataIndex]->SetLineColor(2);
  t[SDataIndex]->GetXaxis()->SetTitle("TRACKISO");
  t[SDataIndex]->GetXaxis()->SetTitleSize(0.07);
  t[BDataIndex]->SetLineColor(4);
  
  TLine* eLine = new TLine(eReach,0,eReach,1.1*e[SDataIndex]->GetMaximum());
  eLine->SetLineColor(2);
  TLine* hLine = new TLine(hReach,0,hReach,1.1*h[SDataIndex]->GetMaximum());
  hLine->SetLineColor(2);
  TLine* tLine = new TLine(tReach,0,tReach,1.1*t[SDataIndex]->GetMaximum());
  tLine->SetLineColor(2);
  
  if(1) { // Drawing the 3 S&B distributions
    sprintf(sochn2,"EHTiso_%s",nameEB[choiceEB]);
    TCanvas* canDist = new TCanvas(sochn2,sochn2,800,800);
    canDist->Divide(1,3); 
    canDist->cd(1)->SetLogy();
    e[SDataIndex]->DrawCopy();
    e[BDataIndex]->SetLineColor(4);
    e[BDataIndex]->DrawCopy("same");
    eLine->Draw("same");
    canDist->cd(2)->SetLogy();
    h[SDataIndex]->DrawCopy();
    h[BDataIndex]->SetLineColor(4);
    h[BDataIndex]->DrawCopy("same");
    hLine->Draw("same");
    canDist->cd(3)->SetLogy();
    t[SDataIndex]->DrawCopy();
    t[BDataIndex]->SetLineColor(4);
    t[BDataIndex]->DrawCopy("same");
    tLine->Draw("same");
       
    canDist->SaveAs(".gif");
    canDist->SaveAs(".eps");
    delete canDist;
  }
  
  if(1) { // Drawing the 3 signal efficiencies
    sprintf(sochn2,"EHTeff_%s",nameEB[choiceEB]);
    TCanvas* canEff = new TCanvas(sochn2,sochn2,1);
    canEff->Divide(1,3);
    canEff->cd(1)->SetLogy();
    eEff->DrawCopy();
    canEff->cd(2)->SetLogy();
    hEff->DrawCopy();
    canEff->cd(3)->SetLogy();
    tEff->DrawCopy();
      
    canEff->SaveAs(".gif");
    canEff->SaveAs(".eps");
    delete canEff;
  }
  

  if(1) { // Drawing the 3 signal efficiencies
    sprintf(sochn2,"EHTSBeff_%s",nameEB[choiceEB]);
    TCanvas* canEff = new TCanvas(sochn2,sochn2,1);
    canEff->Divide(1,3);

    canEff->cd(1)->SetLogy();
    eEff_sig->DrawCopy();
    eEff_bkg->SetLineColor(2);
    eEff_bkg->DrawCopy("same");

    canEff->cd(2)->SetLogy();
    hEff_sig->DrawCopy();
    hEff_bkg->SetLineColor(2);
    hEff_bkg->DrawCopy("same");

    canEff->cd(3)->SetLogy();
    tEff_sig->DrawCopy();
    tEff_bkg->SetLineColor(2);
    tEff_bkg->DrawCopy("same");      

    canEff->SaveAs(".gif");
    canEff->SaveAs(".eps");
    delete canEff;
  }


  if(1) { // Drawing the difference between signal and background efficiency
    sprintf(sochn2,"EHTSBDiffEff_%s",nameEB[choiceEB]);
    TCanvas* canEff = new TCanvas(sochn2,sochn2,1);
    canEff->Divide(1,3);

    canEff->cd(1)->SetLogy();
    eEff_diff->DrawCopy();

    canEff->cd(2)->SetLogy();
    hEff_diff->DrawCopy();

    canEff->cd(3)->SetLogy();
    tEff_diff->DrawCopy();

    canEff->SaveAs(".gif");
    canEff->SaveAs(".eps");
    delete canEff;
  }
  

  
  // PURITY PLOT 
  
  frame->SetMarkerColor(0);
  frame->SetLineColor(0);
  
  purMC->SetMarkerSize(1);
  purMC->SetMarkerColor(2);
  
  purMC->SetLineColor(2);
  
  purCC->SetMarkerSize(1);
  purCC->SetMarkerColor(4);
  purCC->SetLineColor(4);
  
  purFF->SetMarkerSize(1);
  purFF->SetMarkerColor(3);
  purFF->SetLineColor(3);
  
  Double_t legX = 0.2;
  Double_t legY = 0.65;  
  TLegend* leg1 = new TLegend(legX,legY,legX+0.35,legY+0.15);
  leg1->SetFillColor(0);      
  leg1->SetTextSize(0.03);
  leg1->AddEntry(purMC,"MC purity","l");
  leg1->AddEntry(purCC,"Two-bin purity","l");
  leg1->AddEntry(purFF,"Fitted purity","l");
  
   
  gStyle->SetPadLeftMargin(0.15);
  sprintf(sochn2,"PurVsPt_%s",nameEB[choiceEB]);
  TCanvas* canPurPt = new TCanvas(sochn2,sochn2,800,800);
  if(1) {
    canPurPt->cd()->SetLogx();
    purMC->GetXaxis()->SetTitle("Photon E_{t} (GeV)");
    sprintf(sochn2,"Photon Purity (%)",intLum);
    purMC->GetYaxis()->SetTitle(sochn2);
    purMC->SetMinimum(-1);
    purMC->SetMaximum(100);
    purMC->SetMarkerSize(.001);
    purMC->Draw("APZ");
    purCC->Draw("* same");
    purFF->Draw("* same");
    leg1->Draw("same");
    
    canPurPt->cd();
    sprintf(sochn2,"%s",printEB[choiceEB]);
    histLabel.SetTextSize(0.05);
    histLabel.DrawLatex(0.2,0.86,sochn2);
    
    canPurPt->SaveAs(".gif");
    canPurPt->SaveAs(".eps");
    canPurPt->SaveAs(".pdf");
  }
  
  // NB ISOLATED PHOTONS EVENTS Plot
  
  nbMC->SetMarkerSize(1);
  nbMC->SetMarkerColor(2);
  nbMC->GetXaxis()->SetTitle("Photon E_{t} (GeV)");
  sprintf(sochn2,"Isolated photon evts for %.0f pb^{-1}",intLum);
  nbMC->GetYaxis()->SetTitle(sochn2);
  nbMC->SetLineColor(2);
  
  nbCC->SetMarkerSize(1);
  nbCC->SetMarkerColor(4);
  nbCC->SetLineColor(4);
  
  nbFF->SetMarkerSize(1);
  nbFF->SetMarkerColor(3);
  nbFF->SetLineColor(3);
  
  legX = 0.65;
  legY = 0.65;
  TLegend* leg1b = new TLegend(legX,legY,legX+0.3,legY+0.15);
  leg1b->SetFillColor(0);
  leg1b->SetTextSize(0.03);
  leg1b->AddEntry(nbMC,"From MC","l");
  leg1b->AddEntry(nbCC,"From two bins","l");
  leg1b->AddEntry(nbFF,"From fit","l");
  
  sprintf(sochn2,"nbVsPt_%s",nameEB[choiceEB]);
  TCanvas* canNbPt = new TCanvas(sochn2,sochn2,800,800);
  if(1) { 
    canNbPt->cd()->SetLogx();
    canNbPt->cd()->SetLogy();
    nbMC->SetMaximum(1e7);
    nbMC->SetMinimum(.01);
    nbMC->SetMarkerSize(.001);
    nbMC->Draw("APZ");
    nbCC->Draw("* same");
    nbFF->Draw("* same");
    leg1b->Draw("same");
    
    canNbPt->cd();
    sprintf(sochn2,"%s",printEB[choiceEB]);
    histLabel.SetTextSize(0.05);
    histLabel.DrawLatex(0.45,0.85,sochn2);
    
    canNbPt->SaveAs(".gif");
    canNbPt->SaveAs(".eps");
    canNbPt->SaveAs(".pdf");
    gStyle->SetPadLeftMargin(0.08);
  }
  
  // PURITY ERRORS for Isolation
  
  errMC->SetMarkerSize(1);
  errMC->SetMarkerColor(2);
  errMC->GetXaxis()->SetTitle("Photon E_{t} (GeV)");
  errMC->GetYaxis()->SetTitle("Error on photon purity (%)");
  errMC->SetLineColor(2);
  errMC->SetMinimum(-1.);
  errMC->SetMaximum(101.);
  
  errCC->SetMarkerSize(1);
  errCC->SetMarkerColor(4);
  errCC->SetLineColor(4);
  
  errFF->SetMarkerSize(1);
  errFF->SetMarkerColor(3);
  errFF->SetLineColor(3);
  
  legX = 0.15;
  legY = 0.75;  
  TLegend* legErrPur = new TLegend(legX,legY,legX+0.35,legY+0.15);
  legErrPur->SetFillColor(0);      
  legErrPur->SetTextSize(0.03);
  legErrPur->SetHeader("Error on :");
  legErrPur->AddEntry(purMC,"MC purity","l");
  legErrPur->AddEntry(purCC,"Two-bin purity","l");
  legErrPur->AddEntry(purFF,"Fitted purity","l");
  
  if(1) { 
    sprintf(sochn2,"ErrPur_%s",nameEB[choiceEB]);
    TCanvas* canErrPur = new TCanvas(sochn2,sochn2,800,800);
    errMC->Draw("A*");
    errCC->Draw("* same");
    errFF->Draw("* same");
    legErrPur->Draw("same");
    
    canErrPur->cd();
    sprintf(sochn2,"%s",printEB[choiceEB]);
    histLabel.SetTextSize(0.05);
    histLabel.DrawLatex(0.5,0.8,sochn2);
    
    canErrPur->SaveAs(".gif");
    canErrPur->SaveAs(".eps");
    delete canErrPur;
  }
  
  // Distributions for the DATA
  
  for (Int_t i=0;i<nPtBin;i++) {
    iso[SDataIndex][i]->SetLineColor(2);
    iso[SDataIndex][i]->SetLineStyle(2);
    iso[SDataIndex][i]->SetMarkerSize(0.001);
    sprintf(sochn2,"SUM ISO DATA for gamEt in [%.0f,%.0f]",minGamEt[i],maxGamEt[i]);
    iso[SandBDataIndex][i]->GetXaxis()->SetTitle(sochn2);
    iso[BDataIndex][i]->SetLineColor(4);
    iso[BDataIndex][i]->SetMarkerSize(0.001);
    iso[BDataIndex][i]->SetLineStyle(2);
    iso[SandBDataIndex][i]->SetMarkerColor(3);
    iso[SandBDataIndex][i]->SetLineColor(3);
    iso[SandBDataIndex][i]->SetMarkerSize(0.5);
    iso[SandBDataIndex][i]->SetMinimum(0);
  }
  
  if(1) {
    sprintf(sochn2,"DataPtBins_%s",nameEB[choiceEB]);
    TCanvas* canDataPtBins = new TCanvas(sochn2,sochn2,800,800);
    canDataPtBins->Divide(5,3);
    
    for (Int_t i=0;i<nPtBin;i++) {
      canDataPtBins->cd(i);
      iso[SandBDataIndex][i]->DrawCopy("E1");
      iso[SDataIndex][i]->DrawCopy("hist same");
      iso[BDataIndex][i]->DrawCopy("hist same");
    }
    canDataPtBins->SaveAs(".gif");
    canDataPtBins->SaveAs(".eps");
    delete canDataPtBins;
  }
  
  // Distributions for the TEMPLATE
  
  for (Int_t i=0;i<nPtBin;i++) {
    iso[STemplateIndex][i]->SetLineColor(2);
    iso[STemplateIndex][i]->SetMarkerColor(2);
    iso[STemplateIndex][i]->SetMarkerSize(.5);
    sprintf(sochn2,"SUM ISO TEMPLATE for gamEt in [%.0f,%.0f]",minGamEt[i],maxGamEt[i]);
    iso[BTemplateIndex][i]->GetXaxis()->SetTitle(sochn2);
    iso[BTemplateIndex][i]->SetLineColor(4);
    iso[BTemplateIndex][i]->SetMarkerColor(4);
    iso[BTemplateIndex][i]->SetMarkerSize(.5);
    iso[BTemplateIndex][i]->SetMaximum(iso[SandBDataIndex][i]->GetMaximum());
    iso[BTemplateIndex][i]->SetMinimum(0);
  }
  
  if(1) {
    sprintf(sochn2,"TempPtBins_%s",nameEB[choiceEB]);
    TCanvas* canTempPtBins = new TCanvas(sochn2,sochn2,800,800);
    canTempPtBins->Divide(5,3);
    
    for (Int_t i=0;i<nPtBin;i++) {
      if(iso[BTemplateIndex][i]->Integral()) canTempPtBins->cd(i);
      else canTempPtBins->cd(i);
      iso[BTemplateIndex][i]->DrawCopy("E1");
      iso[STemplateIndex][i]->DrawCopy("E1 same");
      cout << "intSTemplateIndex " << iso[STemplateIndex][i]->Integral() << " intBTemplateIndex " << iso[BTemplateIndex][i]->Integral() << endl;
    }
    
    canTempPtBins->SaveAs(".gif");
    canTempPtBins->SaveAs(".eps");
    delete canTempPtBins;
  }
  
  char name[300];
  sprintf(sochn2,"TFracFits_%s",nameEB[choiceEB]);
//   TCanvas* canTFracFits = new TCanvas(sochn2,sochn2,800,800);
  if(doFit) { // Fitted sochn
//     canTFracFits->Divide(5,3);
    
    for (Int_t ptb=0; ptb<nPtBin; ptb++) {
      int ptlo = minGamEt[ptb];
      int pthi = maxGamEt[ptb];
      sprintf(name,"TFracFits_pt%dpt%d", ptlo, pthi);
      cout << "name = " << name << endl;
      TCanvas* canTFracFits = new TCanvas(name,name,500,500);
      if (status[ptb] == 0) {
        outMCpredS[ptb] = fitter[ptb]->GetMCPrediction(0);
        outMCpredB[ptb] = fitter[ptb]->GetMCPrediction(1);

// 	canTFracFits->cd(ptb+1);
        iso[SandBDataIndex][ptb]->SetMinimum(0);
        if(fitresult[ptb]->GetMaximum() > iso[SandBDataIndex][ptb]->GetMaximum()) iso[SandBDataIndex][ptb]->SetMaximum(1.5*fitresult[ptb]->GetMaximum());
        iso[SandBDataIndex][ptb]->DrawCopy("E1");
        fitresult[ptb]->DrawCopy("E1 same");
        iso[STemplateIndex][ptb]->SetLineStyle(2);
        iso[BTemplateIndex][ptb]->SetLineStyle(2);    

        iso[STemplateIndex][ptb]->Scale(iso[SandBDataIndex][ptb]->Integral()*preFFpurity[ptb]/nFill[STemplateIndex][ptb]);
        iso[STemplateIndex][ptb]->DrawCopy("hist same");
        iso[STemplateIndex][ptb]->Scale(nFill[STemplateIndex][ptb]/(iso[SandBDataIndex][ptb]->Integral()*preFFpurity[ptb]));

        iso[BTemplateIndex][ptb]->Scale(iso[SandBDataIndex][ptb]->Integral()*(1-preFFpurity[ptb])/iso[BTemplateIndex][ptb]->Integral());
        iso[BTemplateIndex][ptb]->DrawCopy("hist same");
        iso[BTemplateIndex][ptb]->Scale(iso[BTemplateIndex][ptb]->Integral()/(iso[SandBDataIndex][ptb]->Integral()*(1-preFFpurity[ptb])));
      }
      else // if the fit fails
	{
        iso[SandBDataIndex][ptb]->SetMinimum(0);
        iso[SandBDataIndex][ptb]->DrawCopy("E1");
        iso[STemplateIndex][ptb]->SetLineStyle(2);
        iso[BTemplateIndex][ptb]->SetLineStyle(2);    
        iso[STemplateIndex][ptb]->Scale(iso[SandBDataIndex][ptb]->Integral()*preFFpurity[ptb]/nFill[STemplateIndex][ptb]);
        iso[STemplateIndex][ptb]->DrawCopy("hist same");
        iso[STemplateIndex][ptb]->Scale(nFill[STemplateIndex][ptb]/(iso[SandBDataIndex][ptb]->Integral()*preFFpurity[ptb]));
        iso[BTemplateIndex][ptb]->Scale(iso[SandBDataIndex][ptb]->Integral()*(1-preFFpurity[ptb])/iso[BTemplateIndex][ptb]->Integral());
        iso[BTemplateIndex][ptb]->DrawCopy("hist same");
        iso[BTemplateIndex][ptb]->Scale(iso[BTemplateIndex][ptb]->Integral()/(iso[SandBDataIndex][ptb]->Integral()*(1-preFFpurity[ptb])));
	}

      canTFracFits->SaveAs(".gif");
      canTFracFits->SaveAs(".eps");
      canTFracFits->SaveAs(".pdf");
      delete canTFracFits;
    }
    
//      canTFracFits->SaveAs(".gif");
//      canTFracFits->SaveAs(".eps");
//      canTFracFits->SaveAs(".pdf");
  }

  
  
  
    // BURNING THE REMAINS
  
  graphs->cd();
  sprintf(sochn2,"Purity_MC_%s",nameEB[choiceEB]);
  purMC->Write(sochn2);
  sprintf(sochn2,"Purity_2bin_%s",nameEB[choiceEB]);
  purCC->Write(sochn2);
  sprintf(sochn2,"Purity_Tfrac_%s",nameEB[choiceEB]);
  purFF->Write(sochn2);
  sprintf(sochn2,"NbPho_MC_%s",nameEB[choiceEB]);
  nbMC->Write(sochn2);
  sprintf(sochn2,"NbPho_2bin_%s",nameEB[choiceEB]);
  nbCC->Write(sochn2);
  sprintf(sochn2,"NnPho_Tfrac_%s",nameEB[choiceEB]);
  nbFF->Write(sochn2);
  sprintf(sochn2,"Purity_plot_%s",nameEB[choiceEB]);
  canPurPt->Write(sochn2);
  sprintf(sochn2,"NbPho_plot_%s",nameEB[choiceEB]);
  canNbPt->Write(sochn2);
  sprintf(sochn2,"Fitting_plot_%s",nameEB[choiceEB]);
//   canTFracFits->Write(sochn2);
  graphs->Close();


//   delete canPurPt;
//   delete canNbPt;

  delete eLine;
  delete hLine;
  delete tLine;

  delete purMC ;
  delete purCC ;
  delete purFF ;

  delete errMC ;
  delete errCC ;
  delete errFF ;

  delete nbMC  ;
  delete nbCC  ;
  delete nbFF  ;
  delete frame ;


  delete eEff;
  delete hEff;
  delete tEff;

  delete eEff_sig;
  delete hEff_sig;
  delete tEff_sig;

  delete eEff_bkg;
  delete hEff_bkg;
  delete tEff_bkg;

  delete eEff_diff;
  delete hEff_diff;
  delete tEff_diff;

  delete sumEffS;
  delete sumEffB;
  delete sumEffSB_diff;

  delete leg1;
  delete leg1b;
  delete legErrPur;

  // the following can not be deleted
//   for(Int_t pt=0; pt<nPtBin; pt++)
//      {
//        delete fitter[pt];
//        delete fitmc[pt];
//      }
//     delete indepFile[i];


  
  for (Int_t i=0; i<nSamples; i++){

    delete ptHatDis[i];
    delete e[i];
    delete h[i];
    delete t[i];

    
    for (Int_t pt=0; pt<nPtBin; pt++) {
      delete iso[i][pt];
      
    }
  }
    
   
}


Double_t Laurent_IsoPur::purCCerr(Double_t effB, Double_t effB_N, Double_t effS, Double_t effS_N, Double_t effSnB, Double_t effSnB_N) 
{
  Double_t error, errS, errB, errSnB;
  error = errS = errB = errSnB = 0.;
  errS = effS*(1.-effS)/TMath::Max(1.,effS_N); 
  errB = effB*(1.-effB)/TMath::Max(1.,effB_N); 
  errSnB = effSnB*(1.-effSnB)/TMath::Max(1.,effSnB_N);
  error = TMath::Sqrt(errS + 2.*errB + errSnB);
  //cout << "details purCCerr | effS " << effS << " | effS_N " << effS_N << " | errS " << errS << " | effB " << effB << " | errB " << errB << " | errSnB " << errSnB << endl; 
  return error;
}

Double_t Laurent_IsoPur::purMCerr(Double_t fillS, Double_t scaleS, Double_t fillB, Double_t scaleB) 
{
  Double_t error, errS, errB, errNum, errDen;
  error = errS = errB = errNum = errDen = 0.;
  errS = TMath::Sqrt(fillS)*scaleS; errB = TMath::Sqrt(fillB)*scaleB;
  errNum = errS/TMath::Max(1.,fillS*scaleS + fillB*scaleB); errDen = fillS*scaleS*(errS + errB)/TMath::Max(1.,(fillS*scaleS + fillB*scaleB)*(fillS*scaleS + fillB*scaleB));
  error = TMath::Sqrt(errNum*errNum + errDen*errDen);
  //   cout << "fillS " <<  fillS << " scaleS " << scaleS << " fillB " << fillB << " scaleB " << scaleB << endl;
  //   cout << "errS " << errS << " errB " << errB << " errNum " << errNum << " errDen " << errDen << endl;
  //   cout << "error " << error << endl;
  return error;
}


Bool_t Laurent_IsoPur::CheckFiducial(Int_t choiceEB, Int_t ipho)
{

  if(choiceEB == 1 && !(fabs(caloPositionEta[ipho])>0.0  && fabs(caloPositionEta[ipho])<=0.9))return false;
  if(choiceEB == 2 && !(fabs(caloPositionEta[ipho])>0.9  && fabs(caloPositionEta[ipho])<=1.45))return false;
  if(choiceEB == 3 && !(fabs(caloPositionEta[ipho])>1.55 && fabs(caloPositionEta[ipho])<=2.5))return false;
  return true;

}

Bool_t Laurent_IsoPur::IsLoosePhoton(Bool_t isQCD, Int_t ipho)
{
//   if(!isGenMatched[ipho])return false;
//   // if it's a Photon+Jet events, match photon to hard scattering
//   if(!isQCD && genMomId[ipho]!=22)return false;

//   // if it's a DiJet events, match photon from hadrons
//   if(isQCD && abs(genMomId[ipho])<50)return false;
  	

  bool isBarrel = (isEB[ipho]==1);
  bool isEndCap = (isEE[ipho]==1);
  bool inAnyGap = (isEBEEGap[ipho]==1) || 
    (isEB[ipho]==1 && isEBGap[ipho]==1) || 
    (isEE[ipho]==1 && isEEGap[ipho]==1);
  Float_t thiset      = et[ipho];
  Float_t ecalIso     = ecalRecHitSumEtConeDR04[ipho];
  Float_t hcalIso     = hcalTowerSumEtConeDR04[ipho];
  Float_t trkIso      = trkSumPtHollowConeDR04[ipho];
  Float_t hadem       = hadronicOverEm[ipho];

  if(inAnyGap)return false;

  // Barrel cuts
  if(isBarrel && ecalIso  > 5.0 + 0.004*thiset)return false;
  if(isBarrel && hcalIso  > 5.0)return false;
  if(isBarrel && trkIso > 9.0 )return false;
  if(isBarrel && hadem > 0.15)return false;

  // Endcap cuts
  if(isEndCap && ecalIso  > 5.0 + 0.0021*thiset)return false;
  if(isEndCap && hcalIso  > 5.0)return false;
  if(isEndCap && trkIso > 9.0 )return false;
  if(isEndCap && hadem > 0.15)return false;

  // R9 cut is turned off now
  Double_t R9cut = 0.0;
  if(r9[ipho] < R9cut)return false;

  return true;

}
