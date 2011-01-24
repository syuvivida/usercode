/*==============================================================================================     
                  
Produce 2-D histograms for combined isolation
std::string outputName = "FileName";// prefix for output file name
std::string xvar;                   // variables for making histograms in xaxis   
int xnbin;                         // total number of bins in x axis
double xmin, double xmax;          // minimum and maximum of histograms in xaxis
std::string yvar;                  // variables for making histograms in yaxis
int ynbin;                         // total number of bins in y axis
double ymin, double ymax;          // minimum and maximum of histograms in yaxis
double split;                     // if 1.0 for unsplit and original MC samples, 2.0 for split MC samples
bool normalize;                   // store the histogram by normalizing the integrated area to 1.0 
      
make2DHistos(std::string outputName="";
std::string xvar="(ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04)",
int xnbin=120, double xmin=-1.0, double xmax=11.0, 
std::string yvar="pt",
int ynbin=500, double ymin=0.0, double ymax=500.0,
double split = 1.0,bool normalize=false)

==============================================================================================*/

#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <vector>
#include <string>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include <TEntryList.h>
#include <string>
#include "tightSelection.h"

using namespace std;
const int MAXBIN_ALLOWED=1000;
const double lumi = 36.0; // 1.2/pb
// the pt and eta binning
//const float fBinsEta[]={0,1.4442,1.566,2.5};
const float fBinsEta[]={0,0.9,0.9,1.4442,1.566,2.1,2.1,2.5};

const int nEtaBin = (sizeof(fBinsEta)/sizeof(fBinsEta[0]))/2;


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
	      vector<TH2F*> DHisto,
	      std::string xvar, std::string yvar, 
	      TCut cut,TH2F* h,bool norm)
{
  TH2F *hRes = (TH2F*)h->Clone("hRes");
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
      TH2F *htmp = (TH2F*)h->Clone("htmp");
      htmp->Reset();
      cout << htmp->GetName() << endl;
      sigTree[i]->Draw(Form("%s:(%s)>>htmp",yvar.data(),xvar.data()),allCut);
      htmp->Sumw2();
      htmp->Scale(sigWeight[i]);

      // printing out information
      cout << "scale = " << sigWeight[i] << endl;
      cout << "After scaling htmp -> entries() " << htmp->GetEntries() << endl;
      cout << "After scaling htmp -> Integral() " << htmp->Integral(0,MAXBIN_ALLOWED) << endl;
      cout << "After scaling htmp -> GetMean()  " << htmp->GetMean() << endl;
      cout << "After scaling htmp -> GetRMS()  " << htmp->GetRMS() << endl;
      hRes->Add(htmp);
      
      if(DHisto.size()>0)DHisto[i]->Add(htmp);
      delete htmp;

    }

  h->Sumw2();
  h->Add(hRes);
  if(norm)h->Scale(1.0/(double)h->Integral(0,MAXBIN_ALLOWED,0,MAXBIN_ALLOWED));
  cout << "After scaling h-> entries() " << h->GetEntries() << endl;
  cout << "After scaling h-> Integral() " << h->Integral(0,MAXBIN_ALLOWED) << endl;
  cout << "After scaling h -> GetMean()  " << h->GetMean() << endl;
  cout << "After scaling h -> GetRMS()  " << h->GetRMS() << endl;
   
  delete hRes;
}


// calling functions
void make2DHistos(std::string outputName="", 
		  std::string xvar="(ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04)",
		  int xnbin=300, double xmin=-1.0, double xmax=29.0, 
		  std::string yvar="pt",
		  int ynbin=500, double ymin=0.0, double ymax=500.0,
		  double split = 1.0,
		  bool normalize=false)
{
  TH2F *hTemplate = new TH2F("hTemplate","",xnbin,xmin,xmax,ynbin,ymin,ymax);

  std::string xTitle = "ISO (GeV)";
  std::string yTitle = "p_{T}(#gamma) (GeV/c)";

  char tmp[1000];

  if(outputName=="")outputName = xvar;
  cout << "output file prefix = " << outputName << endl;

  vector<TH2F*> dummyHisto; dummyHisto.clear();
  vector<TH2F*> sigHisto; sigHisto.clear();
  vector<TH2F*> bkgHisto; bkgHisto.clear();
 
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
  int countHistos=0;
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
	countHistos++;
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

   
  TH2F *hIsoMixSig[nEtaBin];
  TH2F *hIsoMixBkg[nEtaBin];
  TH2F *hIsoMixData[nEtaBin];
  TH2F *hIsoMixSBMC[nEtaBin];
  TH2F *hIsoMixSBData[nEtaBin];
  TH2F *hIsoAllSBData[nEtaBin];

  // making templates for spikes
  TH2F *hIsoSpikeData;
  hIsoSpikeData = (TH2F*)hTemplate->Clone("hIsoSpikeData");
  hIsoSpikeData->Reset();
  hIsoSpikeData->SetXTitle(xTitle.data());
  hIsoSpikeData->SetYTitle(yTitle.data());
  cout << "making histograms for spikes " << endl;
  TCut basicSpikeCut = eventCut  + trigCut + rsCutEB + spikeCut;

  makePlot(dataTree,dataWeight,dataPtHatLo,dataPtHatHi,dummyHisto,
	   xvar,yvar,basicSpikeCut,hIsoSpikeData,normalize);

  
  // first looping over eta bins
  for(int ieta = 0; ieta < nEtaBin; ieta++){

    sprintf(tmp,"abs(scEta)>%.5lf && abs(scEta)<%.5lf",
	    fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
    TCut etaCut = tmp;
     
    int etaBinStart = fBinsEta[ieta*2]*10000;
    int etaBinEnd   = fBinsEta[ieta*2+1]*10000;

    hIsoAllSBData[ieta] = (TH2F*)hTemplate->Clone();
    sprintf(tmp,"hIsoAllSBData_Eta_%d_%d", etaBinStart, etaBinEnd);
    hIsoAllSBData[ieta]->SetName(tmp);
    hIsoAllSBData[ieta]->Reset();
    hIsoAllSBData[ieta]->Sumw2();
    hIsoAllSBData[ieta]->SetXTitle(xTitle.data());
     
    //====================== Defining cuts ====================== //

    TCut IDCut    =  fBinsEta[ieta*2+1]*10000 < 15000? rsCutEB: rsCutEE;
    TCut IDSidebandCut = fBinsEta[ieta*2+1]*10000 < 15000? 
      rsCutSidebandEB : rsCutSidebandEE;
    TCut basicCut   = eventCut+ etaCut + IDCut + removeSpikeCut;
    TCut basicSBCut = eventCut+ etaCut + IDSidebandCut + removeSpikeCut;


    //====================== booking histograms ====================== //

    hIsoMixSig[ieta] = (TH2F*)hTemplate->Clone();
    sprintf(tmp,"TemplateS_Eta_%d_%d",
	    etaBinStart,etaBinEnd); 

    hIsoMixSig[ieta]->SetName(tmp);
    hIsoMixSig[ieta]->Reset();
    sprintf(tmp,"%.2f < |#eta(#gamma)| < %.2f",
	    fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
    hIsoMixSig[ieta]->SetTitle(tmp);
    hIsoMixSig[ieta]->SetXTitle(xTitle.data());

    hIsoMixBkg[ieta] = (TH2F*)hTemplate->Clone(); 
    sprintf(tmp,"TemplateB_Eta_%d_%d",
	    etaBinStart,etaBinEnd);
    hIsoMixBkg[ieta]->SetName(tmp);
    hIsoMixBkg[ieta]->Reset();
    sprintf(tmp,"%.2f < |#eta(#gamma)| < %.2f",
	    fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
    hIsoMixBkg[ieta]->SetTitle(tmp);
    hIsoMixBkg[ieta]->SetXTitle(xTitle.data());

    hIsoMixSBData[ieta] = (TH2F*)hTemplate->Clone(); 
    sprintf(tmp,"TemplateSBData_Eta_%d_%d",
	    etaBinStart,etaBinEnd);
    hIsoMixSBData[ieta]->SetName(tmp);
    hIsoMixSBData[ieta]->Reset();
    sprintf(tmp,"%.2f < |#eta(#gamma)| < %.2f",
	    fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
    hIsoMixSBData[ieta]->SetTitle(tmp);
    hIsoMixSBData[ieta]->SetXTitle(xTitle.data());

    hIsoMixSBMC[ieta] = (TH2F*)hTemplate->Clone(); 
    sprintf(tmp,"TemplateSBMC_Eta_%d_%d",
	    etaBinStart,etaBinEnd);
    hIsoMixSBMC[ieta]->SetName(tmp);
    hIsoMixSBMC[ieta]->Reset();
    sprintf(tmp,"%.2f < |#eta(#gamma)| < %.2f",
	    fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
    hIsoMixSBMC[ieta]->SetTitle(tmp);
    hIsoMixSBMC[ieta]->SetXTitle(xTitle.data());


    hIsoMixData[ieta] = (TH2F*)hTemplate->Clone();
    sprintf(tmp,"Template_Eta_%d_%d",
	    etaBinStart,etaBinEnd);
    hIsoMixData[ieta]->SetName(tmp);
    hIsoMixData[ieta]->Reset();
    sprintf(tmp,"%.2f < |#eta(#gamma)| < %.2f",
	    fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
    hIsoMixData[ieta]->SetTitle(tmp);
    hIsoMixData[ieta]->SetXTitle(xTitle.data());
   
    hIsoMixSig[ieta]->SetLineColor(kBlack);
    hIsoMixSig[ieta]->SetMarkerColor(kBlack);
    hIsoMixBkg[ieta]->SetMarkerColor(kRed);
    hIsoMixBkg[ieta]->SetLineColor(kRed);
    hIsoMixSBData[ieta]->SetMarkerColor(kRed);
    hIsoMixSBData[ieta]->SetLineColor(kRed);
    hIsoMixSBMC[ieta]->SetMarkerColor(kRed);
    hIsoMixSBMC[ieta]->SetLineColor(kRed);
    hIsoMixData[ieta]->SetMarkerColor(kBlue);
    hIsoMixData[ieta]->SetLineColor(kBlue);

    //====================== Making histograms ====================== //
       
    cout << "making histograms from mixed MC signal samples" << endl;     
    TCut allCut = basicCut + sigCut;
    // adding distributions of variable from each signal MC file
    sigHisto.clear();
    const int totalSize = mixTree.size();
    TH2F* hIsoSigMCFile[totalSize];
    for(int iv=0; iv < totalSize; iv++){
      hIsoSigMCFile[iv] = (TH2F*)hTemplate->Clone();
      hIsoSigMCFile[iv]->SetName(Form("%s_%02i",hIsoMixSig[ieta]->GetName(),iv));
      hIsoSigMCFile[iv]->SetTitle(mixFile[iv].data());
      hIsoSigMCFile[iv]->Reset();
      sigHisto.push_back(hIsoSigMCFile[iv]);
    }
    makePlot(mixTree,mixWeight,mixPtHatLo,mixPtHatHi,sigHisto,
	     xvar,yvar,allCut,hIsoMixSig[ieta],normalize);

    cout << "making histograms from mixed MC background samples" << endl;     
    allCut = basicCut + bkgCut;
    // adding distributions of variable from each background MC file
    bkgHisto.clear();
    TH2F* hIsoBkgMCFile[totalSize];
    for(int iv=0; iv < totalSize; iv++){
      hIsoBkgMCFile[iv] = (TH2F*)hTemplate->Clone();
      hIsoBkgMCFile[iv]->SetName(Form("%s_%02i",hIsoMixBkg[ieta]->GetName(),iv));
      hIsoBkgMCFile[iv]->SetTitle(mixFile[iv].data());
      hIsoBkgMCFile[iv]->Reset();
      bkgHisto.push_back(hIsoBkgMCFile[iv]);
    }
    makePlot(mixTree,mixWeight,mixPtHatLo,mixPtHatHi,bkgHisto,
	     xvar,yvar,allCut,hIsoMixBkg[ieta],normalize);

    cout << "making histogram for mixed sideband MC background samples" << endl;
    allCut = basicSBCut;
    makePlot(mixTree,mixWeight,mixPtHatLo,mixPtHatHi,dummyHisto,xvar,yvar,allCut,hIsoMixSBMC[ieta],normalize);
     
    cout << "hIsoMixSig[" << ieta << "]->Integral()  = " << hIsoMixSig[ieta]->Integral(0,MAXBIN_ALLOWED,0,MAXBIN_ALLOWED) << endl;
    cout << "hIsoMixBkg[" << ieta << "]->Integral()  = " << hIsoMixBkg[ieta]->Integral(0,MAXBIN_ALLOWED,0,MAXBIN_ALLOWED) << endl;

    if(dataTree.size()>0)
      {
	cout << "making histograms for data samples" << endl;
	allCut = basicCut + trigCut;
	makePlot(dataTree,dataWeight,dataPtHatLo,dataPtHatHi,dummyHisto,xvar,yvar,allCut,hIsoMixData[ieta],normalize);
 	cout << "making histogram for sideband background data" << endl;
 	allCut = basicSBCut + trigCut;
 	makePlot(dataTree,dataWeight,dataPtHatLo,dataPtHatHi,dummyHisto,xvar,yvar,allCut,hIsoMixSBData[ieta],normalize);
	cout << "making histogram for all sideband background data" << endl;
	hIsoAllSBData[ieta]->Add(hIsoMixSBData[ieta],1.0);	   
      }
    else
      {
	hIsoMixData[ieta]->Reset();
	hIsoMixData[ieta]->Sumw2();
	hIsoMixData[ieta]->Add(hIsoMixSig[ieta],hIsoMixBkg[ieta],1.0,1.0);
	cout << "hIsoMixData[" << ieta << "]->Integral()  = " << hIsoMixData[ieta]->Integral(0,MAXBIN_ALLOWED,0,MAXBIN_ALLOWED) << endl;
      }
  } // end of looping over eta bins



  // dump histogram to a root file
  std::string histoFile = outputName + "_histo.root";

  TFile* outFile = new TFile(histoFile.data(),"recreate");
  hIsoSpikeData->Write();
  cout << "hey" << endl;

  cout << "sigHisto.size() = " << sigHisto.size() << endl;
  for(int iv=0; iv < sigHisto.size(); iv++)
    sigHisto[iv]->Write();

  cout << "bkgHisto.size() = " << bkgHisto.size() << endl;
  for(int iv=0; iv < bkgHisto.size(); iv++)
    bkgHisto[iv]->Write();


  for(int ieta = 0; ieta < nEtaBin; ieta++){
     
    hIsoAllSBData[ieta]->Write();

    hIsoMixSig[ieta]->Write();
    hIsoMixBkg[ieta]->Write();
    hIsoMixData[ieta]->Write();
    hIsoMixSBData[ieta]->Write();
    hIsoMixSBMC[ieta]->Write();

  } // end of looping over eta bins

   
  outFile->Close();
   
}
