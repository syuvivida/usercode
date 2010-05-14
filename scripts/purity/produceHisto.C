#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <vector>
#include <string>
#include <TCanvas.h>
#include <TLegend.h>
#include <string>
#include "tightSelection.h"

using namespace std;

void fillTrees(vector<string> &fileName,vector<TTree*> &trees,string treeName)
{
   for (unsigned int i=0; i<fileName.size(); i++){
      TFile *f = new TFile(fileName[i].data());
      TTree *t = (TTree*)f->FindObjectAny(treeName.data());
      if (t==0) cout <<"Error!"<<endl;
      trees.push_back(t);
   }

}

void makePlot(vector<TTree*> sigTree,vector<double> sigWeight,
	      vector<int> ptHatLo, vector<int> ptHatHi,
              std::string var,TCut cut,TH1F* h,bool norm)
{
   TH1F *hRes = (TH1F*)h->Clone();
   hRes->SetName("hRes");
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
       sigTree[i]->Draw(Form("%s>>htmp",var.data()),allCut);
       htmp->Sumw2();
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
   if(norm)h->Scale(1.0/(double)h->Integral(0,1000));
   cout << "After scaling h-> entries() " << h->GetEntries() << endl;
   cout << "After scaling h-> Integral() " << h->Integral() << endl;
   cout << "After scaling h -> GetMean()  " << h->GetMean() << endl;
   cout << "After scaling h -> GetRMS()  " << h->GetRMS() << endl;
   
   delete hRes;
}


void produceHisto(std::string outputName="", std::string var="0.692-9.240*ecalIso/pt-11.117*hcalIso/pt",double genIsoCut=5.0,bool normalize=false )
{
   if(outputName=="")outputName = var;
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


   double lumi = 1.0;
   FILE *fTable = fopen("inputFile.txt","r");
   
   int flag=1;   
   int nfile=0;
   while (flag!=-1){
     // first reading input file
     char filename[100];
     flag=fscanf(fTable,"%s",filename);
     std::string tempFile = filename;

     bool isPhotonJet = false;
     if(tempFile.find("PhotonJet") != std::string::npos)isPhotonJet=true;
     char tmp[1000];
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
	
   } 

   std::string treeName = "Analysis";

   // first photon trees
   fillTrees(phoFile,phoTree,treeName);
   const unsigned int nSize_pho = phoFile.size();
   if(phoTree.size()!= nSize_pho){cout << "error 1"<< endl; return;}
   if(phoWeight.size()!= nSize_pho){cout << "error 2"<< endl; return;}
   if(phoPtHatLo.size()!= nSize_pho){cout << "error 3"<< endl; return;}
   if(phoPtHatHi.size()!= nSize_pho){cout << "error 4"<< endl; return;}

   // second fill background trees
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
   
   TH1F *hTemplate = new TH1F("hTemplate","",70,-5,2);
   TH1F* hEcalIsoPho = (TH1F*)hTemplate->Clone();
   std::string histoName = outputName + "Pho";
   hEcalIsoPho->SetName(histoName.data());

   TH1F* hEcalIsoJet = (TH1F*)hTemplate->Clone();
   histoName = outputName + "Jet";
   hEcalIsoJet->SetName(histoName.data());

   TH1F* hEcalIsoMixSig = (TH1F*)hTemplate->Clone();
   histoName = outputName + "MixSig";
   hEcalIsoMixSig->SetName(histoName.data());

   TH1F* hEcalIsoMixBkg = (TH1F*)hTemplate->Clone();
   histoName = outputName + "MixBkg";
   hEcalIsoMixBkg->SetName(histoName.data());


   TH1F* hEcalIsoMixData = (TH1F*)hTemplate->Clone();
   histoName = outputName + "MixData";
   hEcalIsoMixData->SetName(histoName.data());
   
   cout << "making histograms from photon+jet MC samples" << endl;
   TCut allCut = basicCut + sigCut;
//    makePlot(phoTree,phoWeight,Form("%s",var.data()),hardScatterCut,hEcalIsoPho,normalize);
   makePlot(phoTree,phoWeight,phoPtHatLo,phoPtHatHi,Form("%s",var.data()),allCut,hEcalIsoPho,normalize);

   cout << "making histograms from dijet MC samples" << endl;
   allCut = basicCut + decayCut;
   makePlot(jetTree,jetWeight,jetPtHatLo,jetPtHatHi,Form("%s",var.data()),allCut,hEcalIsoJet,normalize);

   cout << "making histograms from mixed MC signal samples" << endl;     
   allCut = basicCut + sigCut;
   makePlot(mixTree,mixWeight,mixPtHatLo,mixPtHatHi,Form("%s",var.data()),allCut,hEcalIsoMixSig,normalize);

   cout << "making histograms from mixed MC background samples" << endl;     
   allCut = basicCut + bkgCut;
   makePlot(mixTree,mixWeight,mixPtHatLo,mixPtHatHi,Form("%s",var.data()),allCut,hEcalIsoMixBkg,normalize);
   
   hEcalIsoPho->SetMarkerColor(2);
   hEcalIsoPho->SetLineColor(2);
   hEcalIsoJet->SetMarkerColor(1);
   hEcalIsoJet->SetLineColor(1);
   hEcalIsoMixSig->SetMarkerColor(5);
   hEcalIsoMixSig->SetLineColor(5);
   hEcalIsoMixBkg->SetMarkerColor(4);
   hEcalIsoMixBkg->SetLineColor(4);


   cout << "hEcalIsoPho->Integral()  = " << hEcalIsoPho->Integral() << endl;
   cout << "hEcalIsoJet->Integral()  = " << hEcalIsoJet->Integral() << endl;
   cout << "hEcalIsoMixSig->Integral()  = " << hEcalIsoMixSig->Integral() << endl;
   cout << "hEcalIsoMixBkg->Integral()  = " << hEcalIsoMixBkg->Integral() << endl;

   hEcalIsoMixData->Reset();
   hEcalIsoMixData->Sumw2();
   hEcalIsoMixData->Add(hEcalIsoMixSig,hEcalIsoMixBkg,1.0,1.0);
   cout << "hEcalIsoMixData->Integral()  = " << hEcalIsoMixData->Integral() << endl;
   

   hEcalIsoPho->SetXTitle(var.data());
   hEcalIsoPho->Draw("hist");
   hEcalIsoJet->Draw("histesame");
   hEcalIsoMixSig->Draw("histesame");
   hEcalIsoMixBkg->Draw("histesame");
   TLegend* leg = new TLegend(0.5,0.6,0.7,0.9);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->SetTextSize(0.03);
   leg->SetBorderSize(0);
   leg->AddEntry(hEcalIsoPho,"#gamma+jet MC");
   leg->AddEntry(hEcalIsoJet,"Dijet MC");
   leg->AddEntry(hEcalIsoMixSig,"Mixed MC: signal");
   leg->AddEntry(hEcalIsoMixBkg,"Mixed MC: background");
   leg->Draw("same");

   // dump histogram to a root file
   std::string histoFile = outputName + "_histo.root";

   TFile* outFile = new TFile(histoFile.data(),"recreate");
   hEcalIsoPho->Write();
   hEcalIsoJet->Write();
   hEcalIsoMixSig->Write();
   hEcalIsoMixBkg->Write();
   hEcalIsoMixData->Write();
   outFile->Close();


}
