/*==============================================================================================     
                  
  Produce histograms for combined isolation                                                                          
  
                                                                                                                     
  std::string outputName = "FileName";// prefix for output file name                                                 
  std::string var;                    // variables for making histograms    
  int nbin;                           // total number of bins
  double xmin, double xmax;           // minimum and maximum of histograms 
  double split;                      // if 1.0 for unsplit and original MC samples, 2.0 for split MC samples
  bool normalize;                   // store the histogram by normalizing the integrated area to 1.0 
      
  allHistos(std::string outputName="";
  std::string var="(ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04)",
  int nbin=20,
  double xmin=-1.0, double xmax=11.0, double split = 1.0,bool normalize=false)

  ==============================================================================================*/


#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <vector>
#include <string>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include <string>
#include "tightSelection.h"

using namespace std;
const int MAXBIN_ALLOWED=1000;
const double lumi = 0.1; // 100/nb
// the pt and eta binning
const double fBinsEta[]={0,1.45,1.7,2.5};
const double fBinsPt[]={15,20,30,50,80,120,180,250,500};
const int nPtBin = sizeof(fBinsPt)/sizeof(fBinsPt[0])-1;
const int nEtaBin = (sizeof(fBinsEta)/sizeof(fBinsEta[0]))/2;

 

// calculate errors from weighted histograms
void histoError(TH1F* h, int maxbin, double& total_value, double& total_err)
{

  total_value = total_err = 0;
  for(int i=1; i<= maxbin; i++)
    {
      total_value += h->GetBinContent(i);
      total_err   += pow(h->GetBinError(i),2);
      
    }

  total_err = sqrt(total_err);
  return;

}

// calculate errors on the ratios
void ratioErr(Double_t n1, Double_t n1err, Double_t n2, Double_t n2err,
	      Double_t& ratio, Double_t& err)
{

  ratio = (n1)/(n1+n2);
  
  err= pow(1/(n1+n2) - n1/(n1+n2)/(n1+n2),2)*n1err*n1err+
    pow(n1/(n1+n2)/(n1+n2),2)*n2err*n2err; 

  err = sqrt(err);

}



// fill trees
void fillTrees(vector<string> &fileName,vector<TTree*> &trees,string treeName)
{
   for (unsigned int i=0; i<fileName.size(); i++){
      TFile *f = new TFile(fileName[i].data());
      TTree *t = (TTree*)f->FindObjectAny(treeName.data());
      if (t==0) cout <<"Error!"<<endl;
      trees.push_back(t);
   }

}

// core function, make histograms, with pt hat cuts, and scale each histogram according to weights
void makePlot(vector<TTree*> sigTree,vector<double> sigWeight,
	      vector<int> ptHatLo, vector<int> ptHatHi,
              std::string var,TCut cut,TH1F* h,bool norm)
{
   TH1F *hRes = (TH1F*)h->Clone();
   hRes->SetName("hRes");
   hRes->Reset();
   hRes->Sumw2();

   char tmp[300];
  
   for (unsigned int i=0; i<sigTree.size(); i++)
     {
       // first determine the pthat cut
       
       if(ptHatLo[i]>0 && ptHatHi[i]>0)
	 sprintf(tmp, "ptHat >= %d && ptHat <= %d",ptHatLo[i], ptHatHi[i]);
       else
	 sprintf(tmp, "pt>0");	 
       TCut ptHatCut = tmp;
       TCut allCut = cut + ptHatCut;
       cout << "Current cut = " << allCut.GetTitle() << endl;
       TH1F *htmp = (TH1F*)h->Clone();
       htmp->SetName("htmp");
       htmp->Reset();
       sigTree[i]->Draw(Form("%s>>htmp",var.data()),allCut);
       htmp->Sumw2();
       htmp->Scale(sigWeight[i]);
       cout << "scale = " << sigWeight[i] << endl;
       cout << "After scaling htmp -> entries() " << htmp->GetEntries() << endl;
       cout << "After scaling htmp -> Integral() " << htmp->Integral(0,MAXBIN_ALLOWED) << endl;
       cout << "After scaling htmp -> GetMean()  " << htmp->GetMean() << endl;
       cout << "After scaling htmp -> GetRMS()  " << htmp->GetRMS() << endl;
       hRes->Add(htmp);
       delete htmp;
   }
   h->Sumw2();
   h->Add(hRes);
   if(norm)h->Scale(1.0/(double)h->Integral(0,MAXBIN_ALLOWED));
   cout << "After scaling h-> entries() " << h->GetEntries() << endl;
   cout << "After scaling h-> Integral() " << h->Integral(0,MAXBIN_ALLOWED) << endl;
   cout << "After scaling h -> GetMean()  " << h->GetMean() << endl;
   cout << "After scaling h -> GetRMS()  " << h->GetRMS() << endl;
   
   delete hRes;
}


// calling functions
void allHistos(std::string outputName="", 
	       std::string var="(ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04)",
	       int nbin=20,
	       double xmin=-1.0, double xmax=11.0, double split = 1.0,
	       bool normalize=false)
{
  TH1F *hTemplate = new TH1F("hTemplate","",nbin,xmin,xmax);
  char tmp[1000];

  if(outputName=="")outputName = var;
  cout << "output file prefix = " << outputName << endl;
  vector<string> mixFile;
   vector<double> mixWeight;
   vector<TTree*> mixTree;
   vector<int> mixPtHatLo;
   vector<int> mixPtHatHi;

   vector<string> phoFile;
   vector<double> phoWeight;
   vector<TTree*> phoTree;
   vector<int> phoPtHatLo;
   vector<int> phoPtHatHi;

   vector<string> jetFile;
   vector<double> jetWeight;
   vector<TTree*> jetTree;
   vector<int> jetPtHatLo;
   vector<int> jetPtHatHi;

   vector<string> dataFile;
   vector<double> dataWeight;
   vector<TTree*> dataTree;
   vector<int> dataPtHatLo;
   vector<int> dataPtHatHi;

   FILE *fTable = fopen("inputFile.txt","r");
   
   int flag=1;   
   int nfile=0;
   char filename[100];
   while (flag!=-1){
     // first reading input file
     flag=fscanf(fTable,"%s",filename);
     std::string tempFile = filename;

     bool isData      = false;
     bool isPhotonJet = false;

     if(tempFile.find("EGData")    != std::string::npos)isData     =true;
     else if(tempFile.find("PhotonJet") != std::string::npos)isPhotonJet=true;

     // read in x-section
     flag=fscanf(fTable,"%s",tmp);
     double cross=atof(tmp);
     // read in number of events
     flag=fscanf(fTable,"%s",tmp);
     double nevt=atof(tmp);
     double scale = lumi*split*cross/nevt;

     flag=fscanf(fTable,"%s",tmp);
     int ptHatLo=atof(tmp);
     
     flag=fscanf(fTable,"%s",tmp);
     int ptHatHi=atof(tmp);

     if (flag!=-1) {
       if(isData)
	 {
	   cout << "filling data" << endl;
	   dataFile.push_back(tempFile);
	   dataWeight.push_back(1.0);
	   dataPtHatLo.push_back(-1);
	   dataPtHatHi.push_back(-1);      
	 }
       else if(isPhotonJet)
	 {
	   cout << "filling photon jet" << endl;
	   phoFile.push_back(tempFile);
	   phoWeight.push_back(scale);
	   phoPtHatLo.push_back(ptHatLo);
	   phoPtHatHi.push_back(ptHatHi);      
	 }
       else
	 {
	   cout << "filling dijet" << endl;
	   jetFile.push_back(tempFile);
	   jetWeight.push_back(scale);
	   jetPtHatLo.push_back(ptHatLo);
	   jetPtHatHi.push_back(ptHatHi);      
	 }
      
       if(!isData){
	 cout << "filling mixture" << endl;
	 mixFile.push_back(tempFile);
	 mixWeight.push_back(scale);
	 mixPtHatLo.push_back(ptHatLo);
	 mixPtHatHi.push_back(ptHatHi);      
       }
       cout <<filename<<" "<<cross<<" "<<nevt<< " " << ptHatLo << " " << ptHatHi << endl;

       nfile++; 
     }	 
	
   } 

   std::string treeName = "Analysis";

   // fill data trees
   fillTrees(dataFile,dataTree,treeName);
   const unsigned int nSize_data = dataFile.size();
   if(dataTree.size()!= nSize_data){cout << "error 1"<< endl; return;}
   if(dataWeight.size()!= nSize_data){cout << "error 2"<< endl; return;}
   if(dataPtHatLo.size()!= nSize_data){cout << "error 3"<< endl; return;}
   if(dataPtHatHi.size()!= nSize_data){cout << "error 4"<< endl; return;}

   // first photon trees
   fillTrees(phoFile,phoTree,treeName);
   const unsigned int nSize_pho = phoFile.size();
   if(phoTree.size()!= nSize_pho){cout << "error 1"<< endl; return;}
   if(phoWeight.size()!= nSize_pho){cout << "error 2"<< endl; return;}
   if(phoPtHatLo.size()!= nSize_pho){cout << "error 3"<< endl; return;}
   if(phoPtHatHi.size()!= nSize_pho){cout << "error 4"<< endl; return;}

   // second fill dijet MC trees
   fillTrees(jetFile,jetTree,treeName);
   const unsigned int nSize_jet = jetFile.size();
   if(jetTree.size()!= nSize_jet){cout << "error 1"<< endl; return;}
   if(jetWeight.size()!= nSize_jet){cout << "error 2"<< endl; return;}
   if(jetPtHatLo.size()!= nSize_jet){cout << "error 3"<< endl; return;}
   if(jetPtHatHi.size()!= nSize_jet){cout << "error 4"<< endl; return;}


   // last fill mix trees
   fillTrees(mixFile,mixTree,treeName);
   const unsigned int nSize_mix = mixFile.size();
   if(mixTree.size()!= nSize_mix){cout << "error 1"<< endl; return;}
   if(mixWeight.size()!= nSize_mix){cout << "error 2"<< endl; return;}
   if(mixPtHatLo.size()!= nSize_mix){cout << "error 3"<< endl; return;}
   if(mixPtHatHi.size()!= nSize_mix){cout << "error 4"<< endl; return;}

   
   TH1F *hIsoMixSig[nEtaBin][nPtBin];
   TH1F *hIsoMixBkg[nEtaBin][nPtBin];
   TH1F *hIsoMixData[nEtaBin][nPtBin];
   TH1F *hIsoMixSBMC[nEtaBin][nPtBin];
   TH1F *hIsoMixSBData[nEtaBin][nPtBin];
   TH1F *hIsoAllSBData[nEtaBin];

   // making templates for spikes
   TH1F *hIsoSpikeData;
   hIsoSpikeData = (TH1F*)hTemplate->Clone();
   hIsoSpikeData->SetName("hIsoSpikeData");
   hIsoSpikeData->Reset();
   hIsoSpikeData->SetXTitle("ecalIso+hcalIso+trackIso [GeV]");
   cout << "making histograms for spikes " << endl;
   TCut basicSpikeCut = eventCut + rsCutEB + spikeCut;
   makePlot(dataTree,dataWeight,dataPtHatLo,dataPtHatHi,Form("%s",var.data()),
	    basicSpikeCut,hIsoSpikeData,normalize);

  
   // first looping over eta bins
   for(int ieta = 0; ieta < nEtaBin; ieta++){

     sprintf(tmp,"abs(eta)>%.2lf && abs(eta)<%.2lf",
	     fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
     TCut etaCut = tmp;
     
     hIsoAllSBData[ieta] = (TH1F*)hTemplate->Clone();
     sprintf(tmp,"hIsoAllSBData_Eta_%.2lf_%.2lf",
	     fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
     hIsoAllSBData[ieta]->SetName(tmp);
     hIsoAllSBData[ieta]->Reset();
     hIsoAllSBData[ieta]->Sumw2();
     hIsoAllSBData[ieta]->SetXTitle("ecalIso+hcalIso+trackIso [GeV]");
     
     // second, looping over pt bins
     for(int ipt=0; ipt < nPtBin; ipt++){

       sprintf(tmp, "pt >= %d && pt <  %d",(int)fBinsPt[ipt], (int)fBinsPt[ipt+1]);
       //====================== Defining cuts ====================== //

       TCut ptCut    = tmp;     
       TCut IDCut    = ieta==0? rsCutEB: rsCutEE;
       TCut IDSidebandCut = ieta==0? rsCutSidebandEB : rsCutSidebandEE;
       TCut basicCut   = eventCut+ ptCut + etaCut + IDCut + removeSpikeCut;
       TCut basicSBCut = eventCut+ ptCut + etaCut + IDSidebandCut + removeSpikeCut;

       //====================== booking histograms ====================== //

       hIsoMixSig[ieta][ipt] = (TH1F*)hTemplate->Clone();
       sprintf(tmp,"TemplateS_Eta_%.2f_%.2f_Et_%d_%d",
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1],
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1]); 

       hIsoMixSig[ieta][ipt]->SetName(tmp);
       hIsoMixSig[ieta][ipt]->Reset();
       sprintf(tmp,"%.2f < |#eta(#gamma)| < %.2f, %d < p_{T}(#gamma) < %d GeV",
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1],
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1]);
       hIsoMixSig[ieta][ipt]->SetTitle(tmp);
       hIsoMixSig[ieta][ipt]->SetXTitle("ecalIso+hcalIso+trackIso [GeV]");

       hIsoMixBkg[ieta][ipt] = (TH1F*)hTemplate->Clone(); 
       sprintf(tmp,"TemplateB_Eta_%.2f_%.2f_Et_%d_%d",
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1],
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1]);
       hIsoMixBkg[ieta][ipt]->SetName(tmp);
       hIsoMixBkg[ieta][ipt]->Reset();
       sprintf(tmp,"%.2f < |#eta(#gamma)| < %.2f, %d < p_{T}(#gamma) < %d GeV",
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1],
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1]);
       hIsoMixBkg[ieta][ipt]->SetTitle(tmp);
       hIsoMixBkg[ieta][ipt]->SetXTitle("ecalIso+hcalIso+trackIso [GeV]");

       hIsoMixSBData[ieta][ipt] = (TH1F*)hTemplate->Clone(); 
       sprintf(tmp,"TemplateSBData_Eta_%.2f_%.2f_Et_%d_%d",
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1],
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1]);
       hIsoMixSBData[ieta][ipt]->SetName(tmp);
       hIsoMixSBData[ieta][ipt]->Reset();
       sprintf(tmp,"%.2f < |#eta(#gamma)| < %.2f, %d < p_{T}(#gamma) < %d GeV",
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1],
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1]);
       hIsoMixSBData[ieta][ipt]->SetTitle(tmp);
       hIsoMixSBData[ieta][ipt]->SetXTitle("ecalIso+hcalIso+trackIso [GeV]");

       hIsoMixSBMC[ieta][ipt] = (TH1F*)hTemplate->Clone(); 
       sprintf(tmp,"TemplateSBMC_Eta_%.2f_%.2f_Et_%d_%d",
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1],
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1]);
       hIsoMixSBMC[ieta][ipt]->SetName(tmp);
       hIsoMixSBMC[ieta][ipt]->Reset();
       sprintf(tmp,"%.2f < |#eta(#gamma)| < %.2f, %d < p_{T}(#gamma) < %d GeV",
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1],
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1]);
       hIsoMixSBMC[ieta][ipt]->SetTitle(tmp);
       hIsoMixSBMC[ieta][ipt]->SetXTitle("ecalIso+hcalIso+trackIso [GeV]");


       hIsoMixData[ieta][ipt] = (TH1F*)hTemplate->Clone();
       sprintf(tmp,"Template_Eta_%.2f_%.2f_Et_%d_%d",
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1],
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1]);
       hIsoMixData[ieta][ipt]->SetName(tmp);
       hIsoMixData[ieta][ipt]->Reset();
       sprintf(tmp,"%.2f < |#eta(#gamma)| < %.2f,%d < p_{T}(#gamma) < %d GeV",
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1],
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1]);
       hIsoMixData[ieta][ipt]->SetTitle(tmp);
       hIsoMixData[ieta][ipt]->SetXTitle("ecalIso+hcalIso+trackIso [GeV]");
   
       hIsoMixSig[ieta][ipt]->SetLineColor(kBlack);
       hIsoMixSig[ieta][ipt]->SetMarkerColor(kBlack);
       hIsoMixBkg[ieta][ipt]->SetMarkerColor(kRed);
       hIsoMixBkg[ieta][ipt]->SetLineColor(kRed);
       hIsoMixSBData[ieta][ipt]->SetMarkerColor(kRed);
       hIsoMixSBData[ieta][ipt]->SetLineColor(kRed);
       hIsoMixSBMC[ieta][ipt]->SetMarkerColor(kRed);
       hIsoMixSBMC[ieta][ipt]->SetLineColor(kRed);
       hIsoMixData[ieta][ipt]->SetMarkerColor(kBlue);
       hIsoMixData[ieta][ipt]->SetLineColor(kBlue);

       //====================== Making histograms ====================== //
       
       cout << "making histograms from mixed MC signal samples" << endl;     
       TCut allCut = basicCut + sigCut;
       makePlot(mixTree,mixWeight,mixPtHatLo,mixPtHatHi,Form("%s",var.data()),allCut,hIsoMixSig[ieta][ipt],normalize);

       cout << "making histograms from mixed MC background samples" << endl;     
       allCut = basicCut + bkgCut;
       makePlot(mixTree,mixWeight,mixPtHatLo,mixPtHatHi,Form("%s",var.data()),allCut,hIsoMixBkg[ieta][ipt],normalize);

       cout << "making histogram for mixed sideband MC background samples" << endl;
       allCut = basicSBCut + bkgCut;
       makePlot(mixTree,mixWeight,mixPtHatLo,mixPtHatHi,Form("%s",var.data()),allCut,hIsoMixSBMC[ieta][ipt],normalize);
     

       cout << "hIsoMixSig[" << ieta << "][" << ipt << "]->Integral()  = " << hIsoMixSig[ieta][ipt]->Integral(0,MAXBIN_ALLOWED) << endl;
       cout << "hIsoMixBkg[" << ieta << "][" << ipt << "]->Integral()  = " << hIsoMixBkg[ieta][ipt]->Integral(0,MAXBIN_ALLOWED) << endl;

       if(dataTree.size()>0)
	 {
	   cout << "making histograms for data samples" << endl;
	   allCut = basicCut;
	   makePlot(dataTree,dataWeight,dataPtHatLo,dataPtHatHi,Form("%s",var.data()),allCut,hIsoMixData[ieta][ipt],normalize);
	   cout << "making histogram for sideband background data" << endl;
	   allCut = basicSBCut;
	   makePlot(dataTree,dataWeight,dataPtHatLo,dataPtHatHi,Form("%s",var.data()),allCut,hIsoMixSBData[ieta][ipt],normalize);

	   cout << "making histogram for all sideband background data" << endl;
	   hIsoAllSBData[ieta]->Add(hIsoMixSBData[ieta][ipt],1.0);	   
	 }
       else
	 {
	   hIsoMixData[ieta][ipt]->Reset();
	   hIsoMixData[ieta][ipt]->Sumw2();
	   hIsoMixData[ieta][ipt]->Add(hIsoMixSig[ieta][ipt],hIsoMixBkg[ieta][ipt],1.0,1.0);
	   cout << "hIsoMixData[" << ieta << "][" << ipt << "]->Integral()  = " << hIsoMixData[ieta][ipt]->Integral(0,MAXBIN_ALLOWED) << endl;
	 }
     }// end of looping over pt bins
   } // end of looping over eta bins



   // settings for purity TGraphAsymmetryErrors
   double fBinsPtMidPoint[nPtBin]={0};
   double fBinsPtError[nPtBin]={0};
   
   for(int ipt=0; ipt < nPtBin; ipt++)
     {
       fBinsPtMidPoint[ipt] = 0.5*(fBinsPt[ipt+1]+fBinsPt[ipt]);
       fBinsPtError[ipt] = 0.5*(fBinsPt[ipt+1]-fBinsPt[ipt]);
     }

   // purity
   TGraphAsymmErrors *hTruthPurity[nEtaBin];
   TGraphAsymmErrors *hTruthPurity5GeV[nEtaBin];


   int histoMaxBin = hTemplate->GetNbinsX();
   int histoMaxBin5GeV = hTemplate->FindBin(5.0)-1;

   cout << "histoMaxBin = " << histoMaxBin << " and histoMaxBin5GeV = " << histoMaxBin5GeV << endl;

   // calculate MC truth purity
   for(int ieta = 0; ieta < nEtaBin; ieta++){

     double fYPurity[nPtBin]={0};
     double fYPurityErr[nPtBin]={0};
     double fYPurity5GeV[nPtBin]={0};
     double fYPurityErr5GeV[nPtBin]={0};

     for(int ipt=0; ipt < nPtBin; ipt++){

       double nsig     = 0;
       double nbkg     = 0;       
       double nsig_err = 0;
       double nbkg_err = 0;

       histoError(hIsoMixSig[ieta][ipt],histoMaxBin,nsig,nsig_err);
       histoError(hIsoMixBkg[ieta][ipt],histoMaxBin,nbkg,nbkg_err);
 
       cout << "nsig = " << nsig << " +- " << nsig_err << endl;
       cout << "nbkg = " << nbkg << " +- " << nbkg_err << endl;

       if((nsig+nbkg) < 1e-6)continue;
       double purity = 0;
       double purityErr = 0;

       ratioErr(nsig,nsig_err,nbkg,nbkg_err,purity,purityErr);
       fYPurity[ipt]    = purity;
       fYPurityErr[ipt] = purityErr;
       cout << "purity = " << purity << " +- " << purityErr << endl;


       double nsig5GeV = 0;
       double nbkg5GeV = 0;       
       double nsig_err5GeV = 0;
       double nbkg_err5GeV = 0;

       histoError(hIsoMixSig[ieta][ipt],histoMaxBin5GeV,nsig5GeV,nsig_err5GeV);
       histoError(hIsoMixBkg[ieta][ipt],histoMaxBin5GeV,nbkg5GeV,nbkg_err5GeV);
 
       cout << "nsig = " << nsig5GeV << " +- " << nsig_err5GeV << endl;
       cout << "nbkg = " << nbkg5GeV << " +- " << nbkg_err5GeV << endl;


       if((nsig5GeV+nbkg5GeV) < 1e-6)continue;
       double purity5GeV = 0;
       double purityErr5GeV = 0;

       ratioErr(nsig5GeV,nsig_err5GeV,nbkg5GeV,nbkg_err5GeV,purity5GeV,purityErr5GeV);
       fYPurity5GeV[ipt]    = purity5GeV;
       fYPurityErr5GeV[ipt] = purityErr5GeV;
       cout << "purity for isolation < 5 GeV = " << purity5GeV << " +- " << purityErr5GeV << endl;

     } // end of loop over pt bins

     
   // now define purity histogram
     hTruthPurity[ieta] = new TGraphAsymmErrors(nPtBin,fBinsPtMidPoint,fYPurity,
						fBinsPtError,fBinsPtError,
						fYPurityErr,fYPurityErr);
     sprintf(tmp,"hTruthPurity_Eta_%.2f_%.2f",
	     fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
     hTruthPurity[ieta] -> SetName(tmp);
     sprintf(tmp,"%.2f < |#eta(#gamma)| < %.2f",
	     fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
     hTruthPurity[ieta]->SetTitle(tmp);
     hTruthPurity[ieta]->GetXaxis()->SetTitle("p_{T}(#gamma) [GeV/c]");
     hTruthPurity[ieta]->GetYaxis()->SetTitle("MC Truth Purity");


   // now define purity histogram for isolation < 5 GeV
     hTruthPurity5GeV[ieta] = new TGraphAsymmErrors(nPtBin,fBinsPtMidPoint,fYPurity5GeV,
						    fBinsPtError,fBinsPtError,
						    fYPurityErr5GeV,fYPurityErr5GeV);
     sprintf(tmp,"hTruthPurity5GeV_Eta_%.2f_%.2f",
	     fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
     hTruthPurity5GeV[ieta] -> SetName(tmp);
     sprintf(tmp,"%.2f < |#eta(#gamma)| < %.2f for isolation < 5 GeV",
	     fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
     hTruthPurity5GeV[ieta]->SetTitle(tmp);
     hTruthPurity5GeV[ieta]->GetXaxis()->SetTitle("p_{T}(#gamma) [GeV/c]");
     hTruthPurity5GeV[ieta]->GetYaxis()->SetTitle("MC Truth Purity for Iso < 5 GeV");

   } // end of loop over eta bins


   
   
   // display the comparison

   TLegend* leg = new TLegend(0.609,0.722,0.889,0.896);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->SetTextSize(0.040);
   leg->SetBorderSize(0);

   TCanvas* c1 = new TCanvas("c1","",500,500);
   //   c1->SetLogy(1);
   for(int ieta = 0; ieta < nEtaBin; ieta++){
     for(int ipt=0; ipt < nPtBin; ipt++){

       hIsoMixData[ieta][ipt]->Draw("histe");
       hIsoMixSig[ieta][ipt]->Draw("histe,same");
       hIsoMixBkg[ieta][ipt]->Draw("histe,same");

       leg->Clear();
       leg->AddEntry(hIsoMixData[ieta][ipt],"MC Data");
       leg->AddEntry(hIsoMixSig[ieta][ipt],"Signal");
       leg->AddEntry(hIsoMixBkg[ieta][ipt],"Background");
       leg->Draw("same");
    
       sprintf(tmp,"figures/%s_Template_Eta_%.2f_%.2f_Et_%d_%d",outputName.data(),
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1],
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1]);

       std::string histoName = tmp;
       std::string canvasName = histoName + ".eps"; 
       c1->Print(canvasName.data());
       canvasName = histoName + ".gif";
       c1->Print(canvasName.data());
     } // end of loop over pt bins
   } // end of loop over eta bins



   // dump histogram to a root file
   std::string histoFile = outputName + "_histo.root";

   TFile* outFile = new TFile(histoFile.data(),"recreate");

   hIsoSpikeData->Write();

   for(int ieta = 0; ieta < nEtaBin; ieta++){
     hTruthPurity[ieta]->Write();
     hTruthPurity5GeV[ieta]->Write();
     
     hIsoAllSBData[ieta]->Write();

     for(int ipt=0; ipt < nPtBin; ipt++){
       hIsoMixSig[ieta][ipt]->Write();
       hIsoMixBkg[ieta][ipt]->Write();
       hIsoMixData[ieta][ipt]->Write();
       hIsoMixSBData[ieta][ipt]->Write();
       hIsoMixSBMC[ieta][ipt]->Write();

     } // end of looping over pt bins
   } // end of looping over eta bins
   
   outFile->Close();
   
}
