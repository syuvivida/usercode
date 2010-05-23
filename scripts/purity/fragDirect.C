#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <vector>
#include <string>
#include <TCanvas.h>
#include <TLegend.h>
#include <string>
#include "tightSelection.h"

using namespace std;
const int MAXBIN_ALLOWED=1000;
const double lumi = 0.1; // 100/nb
   // the pt and eta binning
const double fBinsEta[]={0,1.45,1.7,2.5};
const double fBinsPt[]={15,20,30,50,80,120};
const int nPtBin = sizeof(fBinsPt)/sizeof(fBinsPt[0])-1;
const int nEtaBin = (sizeof(fBinsEta)/sizeof(fBinsEta[0]))/2;
   


// fill trees, open a list of files, then get trees, and save the trees in the vector
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
       
       sprintf(tmp, "ptHat >= %d && ptHat <= %d",ptHatLo[i], ptHatHi[i]);
       TCut ptHatCut = tmp;
       TCut allCut = cut + ptHatCut;
       cout << "Current cut = " << allCut.GetTitle() << endl;
       TH1F *htmp = (TH1F*)h->Clone();
       htmp->SetName("htmp");
       htmp->Reset(); 
       htmp->Sumw2();
       sigTree[i]->Draw(Form("%s>>htmp",var.data()),allCut);
       htmp->Scale(sigWeight[i]);
       cout << "scale = " << sigWeight[i] << endl;
       cout << "After scaling htmp -> entries() " << htmp->GetEntries() << endl;
       cout << "After scaling htmp -> Integral() " << htmp->Integral() << endl;
       cout << "After scaling htmp -> GetMean()  " << htmp->GetMean() << endl;
       cout << "After scaling htmp -> GetRMS()  " << htmp->GetRMS() << endl;
       hRes->Add(htmp);
       delete htmp;
   }
   h->Sumw2();
   h->Add(hRes);
   if(norm)h->Scale(1.0/(double)h->Integral(0,MAXBIN_ALLOWED));
   cout << "After scaling h-> entries() " << h->GetEntries() << endl;
   cout << "After scaling h-> Integral() " << h->Integral() << endl;
   cout << "After scaling h -> GetMean()  " << h->GetMean() << endl;
   cout << "After scaling h -> GetRMS()  " << h->GetRMS() << endl;
   
   delete hRes;
}


// calling functions, comparing the component of "direct" and "fragmentation" photons
void fragDirect(std::string outputName="", std::string var=
"(ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04)",
	       int nbin=24,
	       double xmin=-1.0, double xmax=11.0, int genIsoCutValue=5,
	       bool normalize=true )
{
  
  // the template histogram that determines the binning, xmin, and xmax for all histograms
  TH1F *hTemplate = new TH1F("hTemplate","",nbin,xmin,xmax);
  char tmp[1000];
  sprintf(tmp, "genCalIsoDR04 < %d",genIsoCutValue);
  TCut genIsoCut = tmp;


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

   // now reading the input file which describes the cross section, number of input events
   // pt hat limits

   FILE *fTable = fopen("inputFile.txt","r");
   
   int flag=1;   
   int nfile=0;

   char filename[100];
   while (flag!=-1){
     // first reading input file
     flag=fscanf(fTable,"%s",filename);
     std::string tempFile = filename;

     bool isPhotonJet = false;
     if(tempFile.find("PhotonJet") != std::string::npos)isPhotonJet=true;
     // read in x-section
     flag=fscanf(fTable,"%s",tmp);
     double cross=atof(tmp);
     // read in number of events
     flag=fscanf(fTable,"%s",tmp);
     double nevt=atof(tmp);
     double scale =lumi*cross/nevt;

     flag=fscanf(fTable,"%s",tmp);
     int ptHatLo=atof(tmp);
     
     flag=fscanf(fTable,"%s",tmp);
     int ptHatHi=atof(tmp);

    if (flag!=-1) {
      cout <<filename<<" "<<cross<<" "<<nevt<< " " << ptHatLo << " " << ptHatHi << endl;
      if(isPhotonJet)
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
      
      cout << "filling mixture" << endl;
      mixFile.push_back(tempFile);
      mixWeight.push_back(scale);
      mixPtHatLo.push_back(ptHatLo);
      mixPtHatHi.push_back(ptHatHi);      

      nfile++; 
    }	 
	
   } // finish reading inputFile.txt

   // now start filling trees according to the information in the inputFile.txt
   std::string treeName = "Analysis";

   // first photon trees
   fillTrees(phoFile,phoTree,treeName);
   const unsigned int nSize_pho = phoFile.size();
   if(phoTree.size()!= nSize_pho){cout << "pho error 1"<< endl; return;}
   if(phoWeight.size()!= nSize_pho){cout << "pho error 2"<< endl; return;}
   if(phoPtHatLo.size()!= nSize_pho){cout << "pho error 3"<< endl; return;}
   if(phoPtHatHi.size()!= nSize_pho){cout << "pho error 4"<< endl; return;}

   // second fill dijet MC trees
   fillTrees(jetFile,jetTree,treeName);
   const unsigned int nSize_jet = jetFile.size();
   if(jetTree.size()!= nSize_jet){cout << "jet error 1"<< endl; return;}
   if(jetWeight.size()!= nSize_jet){cout << "jet error 2"<< endl; return;}
   if(jetPtHatLo.size()!= nSize_jet){cout << "jet error 3"<< endl; return;}
   if(jetPtHatHi.size()!= nSize_jet){cout << "jet error 4"<< endl; return;}


   // last fill mix trees
   fillTrees(mixFile,mixTree,treeName);
   const unsigned int nSize_mix = mixFile.size();
   if(mixTree.size()!= nSize_mix){cout << "mix error 1"<< endl; return;}
   if(mixWeight.size()!= nSize_mix){cout << "mix error 2"<< endl; return;}
   if(mixPtHatLo.size()!= nSize_mix){cout << "mix error 3"<< endl; return;}
   if(mixPtHatHi.size()!= nSize_mix){cout << "mix error 4"<< endl; return;}


   // for comparison between direct and fragmentation
   TH1F *hIsoPho_1[nEtaBin][nPtBin];
   TH1F *hIsoPho_2[nEtaBin][nPtBin];


   // first looping over eta bins
   for(int ieta = 0; ieta < nEtaBin; ieta++){

     sprintf(tmp,"abs(eta)>%.2lf && abs(eta)<%.2lf",fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
     TCut etaCut = tmp;
     
     // second, looping over pt bins
     for(int ipt=0; ipt < nPtBin; ipt++){

       sprintf(tmp, "pt >= %d && pt < %d",(int)fBinsPt[ipt], (int)fBinsPt[ipt+1]);
       TCut ptCut    = tmp;     
       TCut kineCut  = ptCut + etaCut;
       TCut basicCut = eventCut+ ptCut + etaCut + rsCut + removeSpikeCut + genIsoCut;

       hIsoPho_1[ieta][ipt] = (TH1F*)hTemplate->Clone();
       sprintf(tmp,"Direct_Template_Et_%d_%d_Eta_%.2f_%.2f",(int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
       hIsoPho_1[ieta][ipt]->SetName(tmp);
       hIsoPho_1[ieta][ipt]->SetTitle(kineCut.GetTitle());
       sprintf(tmp,"%d < p_{T}(#gamma) < %d GeV, %.2f < |#eta(#gamma)| < %.2f",
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
       hIsoPho_1[ieta][ipt]->SetTitle(tmp);
       hIsoPho_1[ieta][ipt]->SetXTitle("ecalIso+hcalIso+trackIso [GeV]");

       hIsoPho_2[ieta][ipt] = (TH1F*)hTemplate->Clone();
       sprintf(tmp,"Fragmentation_Template_Et_%d_%d_Eta_%.2f_%.2f",(int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
       hIsoPho_2[ieta][ipt]->SetName(tmp);
       hIsoPho_2[ieta][ipt]->SetTitle(kineCut.GetTitle());
       sprintf(tmp,"%d < p_{T}(#gamma) < %d GeV, %.2f < |#eta(#gamma)| < %.2f",
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
       hIsoPho_2[ieta][ipt]->SetTitle(tmp);
       hIsoPho_2[ieta][ipt]->SetXTitle("ecalIso+hcalIso+trackIso [GeV]");
       

       cout << "making direct histograms from mixed MC samples" << endl;
       TCut allCut   = basicCut + hardScatterCut;
       makePlot(mixTree,mixWeight,mixPtHatLo,mixPtHatHi,Form("%s",var.data()),allCut,hIsoPho_1[ieta][ipt],normalize);

       cout << "making fragmentation histograms from mixed MC samples" << endl;
       allCut = basicCut + fragCut;
       makePlot(mixTree,mixWeight,mixPtHatLo,mixPtHatHi,Form("%s",var.data()),allCut,hIsoPho_2[ieta][ipt],normalize);
     
     }// end of looping over pt bins
   } // end of looping over eta bins



   // display the comparison

   TLegend* leg = new TLegend(0.52 ,0.706,0.72 ,0.905);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->SetTextSize(0.04);
   leg->SetBorderSize(0);

   TCanvas* c1 = new TCanvas("c1","",500,500);
   //   c1->SetLogy(1);
   for(int ieta = 0; ieta < nEtaBin; ieta++){
     for(int ipt=0; ipt < nPtBin; ipt++){

       hIsoPho_1[ieta][ipt]->GetYaxis()->SetDecimals();
       hIsoPho_1[ieta][ipt]->SetMarkerColor(kBlack);
       hIsoPho_1[ieta][ipt]->SetLineColor(kBlack);

       hIsoPho_2[ieta][ipt]->GetYaxis()->SetDecimals();
       hIsoPho_2[ieta][ipt]->SetMarkerColor(kRed);
       hIsoPho_2[ieta][ipt]->SetLineColor(kRed);
       hIsoPho_2[ieta][ipt]->SetLineStyle(2);


       Float_t max1 = hIsoPho_1[ieta][ipt]->GetMaximum() + 
	 hIsoPho_1[ieta][ipt]->GetBinError(hIsoPho_1[ieta][ipt]->GetMaximumBin());
       Float_t max2 = hIsoPho_2[ieta][ipt]->GetMaximum() + 
	 hIsoPho_2[ieta][ipt]->GetBinError(hIsoPho_2[ieta][ipt]->GetMaximumBin());

       if(max1>max2)
	 {
	   hIsoPho_1[ieta][ipt]->Draw("histe");
	   hIsoPho_2[ieta][ipt]->Draw("histesame");
	 }
       else
	 {
	   hIsoPho_2[ieta][ipt]->Draw("histe");
	   hIsoPho_1[ieta][ipt]->Draw("histesame");
	 }

       leg->Clear();
       sprintf(tmp,"#gamma + jet and Dijet MC (genIso < %d)",genIsoCutValue);
       //       leg->SetHeader(tmp);
       leg->AddEntry(hIsoPho_1[ieta][ipt],"Direct");
       leg->AddEntry(hIsoPho_2[ieta][ipt],"Fragmentation");
       leg->Draw("same");
    
       sprintf(tmp,"/mc/QCD_mess/figures/%s_Template_Et_%d_%d_Eta_%.2f_%.2f_genIso%d",outputName.data(),
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1],genIsoCutValue);

       std::string histoName = tmp;
       std::string canvasName = histoName + ".eps"; 
       c1->Print(canvasName.data()); 
       canvasName = histoName + ".gif";
       c1->Print(canvasName.data());
//        canvasName = histoName + ".C";
//        c1->Print(canvasName.data());
     } // end of loop over pt bins
   } // end of loop over eta bins



   // dump histogram to a root file
   std::string histoFile = outputName + "_histo.root";

   TFile* outFile = new TFile(histoFile.data(),"recreate");

   for(int ieta = 0; ieta < nEtaBin; ieta++){
        for(int ipt=0; ipt < nPtBin; ipt++){

       hIsoPho_1[ieta][ipt]->Write();
       hIsoPho_2[ieta][ipt]->Write();

     } // end of looping over pt bins
   } // end of looping over eta bins

   outFile->Close();

}
