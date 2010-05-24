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
 

// calculate errors from weighted histograms
void histoError(TH1F* h, double& total_value, double& total_err)
{

  total_value = total_err = 0;
//   for(int i=0; i<= h->GetNbinsX()+1; i++)
  for(int i=1; i<= h->GetNbinsX(); i++)
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
void allHistos(std::string outputName="", std::string var="(ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04)",
	       int nbin=20,
	       double xmin=-1.0, double xmax=11.0,
	       bool isData=false,bool normalize=false)
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
	
   } 

   std::string treeName = "Analysis";

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
   

   // purity
   TH1F *hTruthPurity[nEtaBin];
  
   // first looping over eta bins
   for(int ieta = 0; ieta < nEtaBin; ieta++){

     sprintf(tmp,"abs(eta)>%.2lf && abs(eta)<%.2lf",
	     fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
     TCut etaCut = tmp;

     sprintf(tmp,"hTruthPurity_Eta_%.2f_%.2f",
	     fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
     hTruthPurity[ieta] = new TH1F(tmp,"",nPtBin, fBinsPt); 

     sprintf(tmp,"%.2f < |#eta(#gamma)| < %.2f",
	     fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
     hTruthPurity[ieta]->SetTitle(tmp);
     hTruthPurity[ieta]->SetXTitle("p_{T}(#gamma) [GeV/c]");
     
     // second, looping over pt bins
     for(int ipt=0; ipt < nPtBin; ipt++){

       sprintf(tmp, "pt >= %d && pt <  %d",(int)fBinsPt[ipt], (int)fBinsPt[ipt+1]);
       TCut ptCut    = tmp;     
       TCut kineCut  = ptCut + etaCut;
       TCut basicCut = eventCut+ ptCut + etaCut + rsCut + removeSpikeCut;

       hIsoMixSig[ieta][ipt] = (TH1F*)hTemplate->Clone();
       sprintf(tmp,"TemplateS_Et_%d_%d_Eta_%.2f_%.2f",
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
       hIsoMixSig[ieta][ipt]->SetName(tmp);
       sprintf(tmp,"%d < p_{T}(#gamma) < %d GeV, %.2f < |#eta(#gamma)| < %.2f",
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
       hIsoMixSig[ieta][ipt]->SetTitle(tmp);
       hIsoMixSig[ieta][ipt]->SetXTitle("ecalIso+hcalIso+trackIso [GeV]");

       hIsoMixBkg[ieta][ipt] = (TH1F*)hTemplate->Clone(); 
       sprintf(tmp,"TemplateB_Et_%d_%d_Eta_%.2f_%.2f",
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
       hIsoMixBkg[ieta][ipt]->SetName(tmp);
       sprintf(tmp,"%d < p_{T}(#gamma) < %d GeV, %.2f < |#eta(#gamma)| < %.2f",
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
       hIsoMixBkg[ieta][ipt]->SetTitle(tmp);
       hIsoMixBkg[ieta][ipt]->SetXTitle("ecalIso+hcalIso+trackIso [GeV]");


       hIsoMixData[ieta][ipt] = (TH1F*)hTemplate->Clone();
       sprintf(tmp,"Template_Et_%d_%d_Eta_%.2f_%.2f",
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
       hIsoMixData[ieta][ipt]->SetName(tmp);
       sprintf(tmp,"%d < p_{T}(#gamma) < %d GeV, %.2f < |#eta(#gamma)| < %.2f",
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
       hIsoMixData[ieta][ipt]->SetTitle(tmp);
       hIsoMixData[ieta][ipt]->SetXTitle("ecalIso+hcalIso+trackIso [GeV]");
   
       hIsoMixSig[ieta][ipt]->SetLineColor(kBlack);
       hIsoMixSig[ieta][ipt]->SetMarkerColor(kBlack);
       hIsoMixBkg[ieta][ipt]->SetMarkerColor(kRed);
       hIsoMixBkg[ieta][ipt]->SetLineColor(kRed);
       hIsoMixData[ieta][ipt]->SetMarkerColor(kBlue);
       hIsoMixData[ieta][ipt]->SetLineColor(kBlue);

       cout << "making histograms from mixed MC signal samples" << endl;     
       TCut allCut = basicCut + sigCut;
       makePlot(mixTree,mixWeight,mixPtHatLo,mixPtHatHi,Form("%s",var.data()),allCut,hIsoMixSig[ieta][ipt],normalize);

       cout << "making histograms from mixed MC background samples" << endl;     
       allCut = basicCut + bkgCut;
       makePlot(mixTree,mixWeight,mixPtHatLo,mixPtHatHi,Form("%s",var.data()),allCut,hIsoMixBkg[ieta][ipt],normalize);
   

       cout << "hIsoMixSig[" << ieta << "][" << ipt << "]->Integral()  = " << hIsoMixSig[ieta][ipt]->Integral(0,MAXBIN_ALLOWED) << endl;
       cout << "hIsoMixBkg[" << ieta << "][" << ipt << "]->Integral()  = " << hIsoMixBkg[ieta][ipt]->Integral(0,MAXBIN_ALLOWED) << endl;

       hIsoMixData[ieta][ipt]->Reset();
       hIsoMixData[ieta][ipt]->Sumw2();
       hIsoMixData[ieta][ipt]->Add(hIsoMixSig[ieta][ipt],hIsoMixBkg[ieta][ipt],1.0,1.0);
       cout << "hIsoMixData[" << ieta << "][" << ipt << "]->Integral()  = " << hIsoMixData[ieta][ipt]->Integral(0,MAXBIN_ALLOWED) << endl;

     }// end of looping over pt bins
   } // end of looping over eta bins


   // calculate MC truth purity
   for(int ieta = 0; ieta < nEtaBin; ieta++){
     for(int ipt=0; ipt < nPtBin; ipt++){

       double nsig     = 0;
       double nbkg     = 0;       
       double nsig_err = 0;
       double nbkg_err = 0;

       histoError(hIsoMixSig[ieta][ipt],nsig,nsig_err);
       histoError(hIsoMixBkg[ieta][ipt],nbkg,nbkg_err);
 
       cout << "nsig = " << nsig << " +- " << nsig_err << endl;
       cout << "nbkg = " << nbkg << " +- " << nbkg_err << endl;

       if((nsig+nbkg) < 1e-6)continue;
       double purity = 0;
       double purityErr = 0;

       ratioErr(nsig,nsig_err,nbkg,nbkg_err,purity,purityErr);

       cout << "purity = " << purity << " +- " << purityErr << endl;

       hTruthPurity[ieta]->SetBinContent(ipt+1,purity);
       hTruthPurity[ieta]->SetBinError(ipt+1,purityErr);
       
     }
   }

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
    
       sprintf(tmp,"/mc/QCD_mess/figures/%s_Template_Et_%d_%d_Eta_%.2f_%.2f",outputName.data(),
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);

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

   for(int ieta = 0; ieta < nEtaBin; ieta++){
     hTruthPurity[ieta]->Write();
     for(int ipt=0; ipt < nPtBin; ipt++){
       hIsoMixSig[ieta][ipt]->Write();
       hIsoMixBkg[ieta][ipt]->Write();
       hIsoMixData[ieta][ipt]->Write();

     } // end of looping over pt bins
   } // end of looping over eta bins
   
   outFile->Close();
   
}
