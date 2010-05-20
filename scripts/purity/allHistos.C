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
const int MAXBIN_ALLOWED=1000;

// calculate errors from weighted histograms
void histoError(TH1F* h, double& total_value, double& total_err)
{

  total_value = total_err = 0;
  for(int i=0; i<= h->GetNbinsX()+1; i++)
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


// calling functions
void allHistos(std::string outputName="", std::string var="(ecalRecHitSumEtConeDR04+hcalTowerSumEtConeDR04+trkSumPtHollowConeDR04)",
	       int nbin=40,
	       double xmin=-1.0, double xmax=39.0,
	       bool isData=false,bool normalize=true )
{
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


   double lumi = 1.0; // 1.0/pb
   FILE *fTable = fopen("inputFile.txt","r");
   
   int flag=1;   
   int nfile=0;
   char tmp[1000];
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

   double fBinsEta[]={0,1.45,1.7,2.5};
   double fBinsPt[]={15,20,30,50,80,120,5000};
   const int nPtBin = sizeof(fBinsPt)/sizeof(fBinsPt[0])-1;
   const int nEtaBin = 2; // 0~1.45, 1.55~2.5
   
   TH1F *hTemplate = new TH1F("hTemplate","",nbin,xmin,xmax);

   TH1F *hIsoPho[nEtaBin][nPtBin];
   TH1F *hIsoJet[nEtaBin][nPtBin];
   TH1F *hIsoMixSig[nEtaBin][nPtBin];
   TH1F *hIsoMixBkg[nEtaBin][nPtBin];
   TH1F *hIsoMixData[nEtaBin][nPtBin];
   
   // for comparison
   TH1F *hIsoPho_1[nEtaBin][nPtBin];
   TH1F *hIsoPho_2[nEtaBin][nPtBin];
   TH1F *hIsoPho_3[nEtaBin][nPtBin];
   TH1F *hIsoPho_4[nEtaBin][nPtBin];
   TH1F *hIsoPho_5[nEtaBin][nPtBin];


   // purity
   TH1F *hTruthPurity[nEtaBin];
  
   // first looping over eta bins
   for(int ieta = 0; ieta < nEtaBin; ieta++){

     TCut etaCut   = (ieta == 0)?  "abs(eta)<1.45": "abs(eta)>1.55 && abs(eta)<2.5";
     sprintf(tmp,"hTruthPurity%02i",ieta);
     hTruthPurity[ieta] = new TH1F(tmp,"",nPtBin, fBinsPt); 
     hTruthPurity[ieta]->SetTitle(etaCut.GetTitle());
     hTruthPurity[ieta]->SetXTitle("p_{T}(#gamma) [GeV/c]");
     
     // second, looping over pt bins
     for(int ipt=0; ipt < nPtBin; ipt++){

       sprintf(tmp, "pt >= %d && pt <= %d",(int)fBinsPt[ipt], (int)fBinsPt[ipt+1]);
       TCut ptCut    = tmp;     
       TCut kineCut  = ptCut + etaCut;
       TCut basicCut = ptCut + etaCut + rsCut;

       hIsoPho[ieta][ipt] = (TH1F*)hTemplate->Clone();
       std::string histoName = "isoPho";
       sprintf(tmp,"%s%02i%02i",histoName.data(),ieta,ipt);
       hIsoPho[ieta][ipt]->SetName(tmp);
       hIsoPho[ieta][ipt]->SetTitle(kineCut.GetTitle());
       hIsoPho[ieta][ipt]->SetXTitle(var.data());

       hIsoPho_1[ieta][ipt] = (TH1F*)hTemplate->Clone();
       histoName = "isoPho_1";
       sprintf(tmp,"%s%02i%02i",histoName.data(),ieta,ipt);
       hIsoPho_1[ieta][ipt]->SetName(tmp);
       hIsoPho_1[ieta][ipt]->SetTitle(kineCut.GetTitle());
       hIsoPho_1[ieta][ipt]->SetXTitle(var.data());

       hIsoPho_2[ieta][ipt] = (TH1F*)hTemplate->Clone();
       histoName = "isoPho_2";
       sprintf(tmp,"%s%02i%02i",histoName.data(),ieta,ipt);
       hIsoPho_2[ieta][ipt]->SetName(tmp);
       hIsoPho_2[ieta][ipt]->SetTitle(kineCut.GetTitle());
       hIsoPho_2[ieta][ipt]->SetXTitle(var.data());

       hIsoPho_3[ieta][ipt] = (TH1F*)hTemplate->Clone();
       histoName = "isoPho_3";
       sprintf(tmp,"%s%02i%02i",histoName.data(),ieta,ipt);
       hIsoPho_3[ieta][ipt]->SetName(tmp);
       hIsoPho_3[ieta][ipt]->SetTitle(kineCut.GetTitle());
       hIsoPho_3[ieta][ipt]->SetXTitle(var.data());

       hIsoPho_4[ieta][ipt] = (TH1F*)hTemplate->Clone();
       histoName = "isoPho_4";
       sprintf(tmp,"%s%02i%02i",histoName.data(),ieta,ipt);
       hIsoPho_4[ieta][ipt]->SetName(tmp);
       hIsoPho_4[ieta][ipt]->SetTitle(kineCut.GetTitle());
       hIsoPho_4[ieta][ipt]->SetXTitle(var.data());

       hIsoPho_5[ieta][ipt] = (TH1F*)hTemplate->Clone();
       histoName = "isoPho_5";
       sprintf(tmp,"%s%02i%02i",histoName.data(),ieta,ipt);
       hIsoPho_5[ieta][ipt]->SetName(tmp);
       hIsoPho_5[ieta][ipt]->SetTitle(kineCut.GetTitle());
       hIsoPho_5[ieta][ipt]->SetXTitle(var.data());


       hIsoJet[ieta][ipt] = (TH1F*)hTemplate->Clone();
       histoName = "isoJet";
       sprintf(tmp,"%s%02i%02i",histoName.data(),ieta,ipt);
       hIsoJet[ieta][ipt]->SetName(tmp);
       hIsoJet[ieta][ipt]->SetTitle(kineCut.GetTitle());
       hIsoJet[ieta][ipt]->SetXTitle(var.data());

       hIsoMixSig[ieta][ipt] = (TH1F*)hTemplate->Clone();
       sprintf(tmp,"TemplateS_Et_%d_%d_Eta_%.1f_%.1f",(int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
       hIsoMixSig[ieta][ipt]->SetName(tmp);
       hIsoMixSig[ieta][ipt]->SetTitle(kineCut.GetTitle());
       hIsoMixSig[ieta][ipt]->SetXTitle(var.data());

       hIsoMixBkg[ieta][ipt] = (TH1F*)hTemplate->Clone(); 
       sprintf(tmp,"TemplateB_Et_%d_%d_Eta_%.1f_%.1f",(int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
       sprintf(tmp,"%s%02i%02i",histoName.data(),ieta,ipt);
       hIsoMixBkg[ieta][ipt]->SetName(tmp);
       hIsoMixBkg[ieta][ipt]->SetTitle(kineCut.GetTitle());
       hIsoMixBkg[ieta][ipt]->SetXTitle(var.data());


       hIsoMixData[ieta][ipt] = (TH1F*)hTemplate->Clone();
       sprintf(tmp,"Template_Et_%d_%d_Eta_%.1f_%.1f",(int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
       hIsoMixData[ieta][ipt]->SetName(tmp);
       hIsoMixData[ieta][ipt]->SetTitle(kineCut.GetTitle());
       hIsoMixData[ieta][ipt]->SetXTitle(var.data());
   

       hIsoPho[ieta][ipt]->SetMarkerColor(2);
       hIsoPho[ieta][ipt]->SetLineColor(2);
       hIsoJet[ieta][ipt]->SetMarkerColor(1);
       hIsoJet[ieta][ipt]->SetLineColor(1);
       hIsoMixSig[ieta][ipt]->SetMarkerColor(5);
       hIsoMixSig[ieta][ipt]->SetLineColor(5);
       hIsoMixBkg[ieta][ipt]->SetMarkerColor(4);
       hIsoMixBkg[ieta][ipt]->SetLineColor(4);

       cout << "making histograms from photon+jet MC samples" << endl;
       TCut allCut   = basicCut + sigCut; 
       makePlot(phoTree,phoWeight,phoPtHatLo,phoPtHatHi,Form("%s",var.data()),allCut,hIsoPho[ieta][ipt],normalize);
      
       allCut = basicCut + sigCut1;
       makePlot(mixTree,mixWeight,mixPtHatLo,mixPtHatHi,Form("%s",var.data()),allCut,hIsoPho_1[ieta][ipt],normalize);

       allCut = basicCut + sigCut2;
       makePlot(mixTree,mixWeight,mixPtHatLo,mixPtHatHi,Form("%s",var.data()),allCut,hIsoPho_2[ieta][ipt],normalize);

       allCut = basicCut + sigCut3;
       makePlot(mixTree,mixWeight,mixPtHatLo,mixPtHatHi,Form("%s",var.data()),allCut,hIsoPho_3[ieta][ipt],normalize);

       allCut = basicCut + sigCut4;
       makePlot(mixTree,mixWeight,mixPtHatLo,mixPtHatHi,Form("%s",var.data()),allCut,hIsoPho_4[ieta][ipt],normalize);

       allCut = basicCut + sigCut5;
       makePlot(mixTree,mixWeight,mixPtHatLo,mixPtHatHi,Form("%s",var.data()),allCut,hIsoPho_5[ieta][ipt],normalize);

       cout << "making histograms from dijet MC samples" << endl;
       allCut = basicCut + decayCut;
       makePlot(jetTree,jetWeight,jetPtHatLo,jetPtHatHi,Form("%s",var.data()),allCut,hIsoJet[ieta][ipt],normalize);

       cout << "making histograms from mixed MC signal samples" << endl;     
       allCut = basicCut + sigCut;
       makePlot(mixTree,mixWeight,mixPtHatLo,mixPtHatHi,Form("%s",var.data()),allCut,hIsoMixSig[ieta][ipt],normalize);

       cout << "making histograms from mixed MC background samples" << endl;     
       allCut = basicCut + bkgCut;
       makePlot(mixTree,mixWeight,mixPtHatLo,mixPtHatHi,Form("%s",var.data()),allCut,hIsoMixBkg[ieta][ipt],normalize);
   


       cout << "hIsoPho[" << ieta << "][" << ipt << "]->Integral()  = " << hIsoPho[ieta][ipt]->Integral() << endl;
       cout << "hIsoJet[" << ieta << "][" << ipt << "]->Integral()  = " << hIsoJet[ieta][ipt]->Integral() << endl;
       cout << "hIsoMixSig[" << ieta << "][" << ipt << "]->Integral()  = " << hIsoMixSig[ieta][ipt]->Integral() << endl;
       cout << "hIsoMixBkg[" << ieta << "][" << ipt << "]->Integral()  = " << hIsoMixBkg[ieta][ipt]->Integral() << endl;

       hIsoMixData[ieta][ipt]->Reset();
       hIsoMixData[ieta][ipt]->Sumw2();
       hIsoMixData[ieta][ipt]->Add(hIsoMixSig[ieta][ipt],hIsoMixBkg[ieta][ipt],1.0,1.0);
       cout << "hIsoMixData[" << ieta << "][" << ipt << "]->Integral()  = " << hIsoMixData[ieta][ipt]->Integral() << endl;

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

   TLegend* leg = new TLegend(0.421,0.606,0.621,0.905);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->SetTextSize(0.03);
   leg->SetBorderSize(0);

   TCanvas* c1 = new TCanvas("c1","",500,500);
   //   c1->SetLogy(1);
   for(int ieta = 0; ieta < nEtaBin; ieta++){
     for(int ipt=0; ipt < nPtBin; ipt++){

       hIsoPho_1[ieta][ipt]->Draw("histe");
       hIsoPho_2[ieta][ipt]->Draw("histesame");
       hIsoPho_3[ieta][ipt]->Draw("histesame");
       hIsoPho_4[ieta][ipt]->Draw("histesame");
       hIsoPho_5[ieta][ipt]->Draw("histesame");

       hIsoPho_1[ieta][ipt]->SetMarkerColor(kBlue);
       hIsoPho_1[ieta][ipt]->SetLineColor(kBlue);

       hIsoPho_2[ieta][ipt]->SetMarkerColor(kRed);
       hIsoPho_2[ieta][ipt]->SetLineColor(kRed);
       hIsoPho_2[ieta][ipt]->SetLineStyle(2);

       hIsoPho_3[ieta][ipt]->SetMarkerColor(kGreen+2);
       hIsoPho_3[ieta][ipt]->SetLineColor(kGreen+2);

       hIsoPho_4[ieta][ipt]->SetMarkerColor(kBlack);
       hIsoPho_4[ieta][ipt]->SetLineColor(kBlack);

       hIsoPho_5[ieta][ipt]->SetMarkerColor(kViolet-5);
       hIsoPho_5[ieta][ipt]->SetLineColor(kViolet-5);
       hIsoPho_5[ieta][ipt]->SetLineStyle(2);


       leg->Clear();
       leg->SetHeader("#gamma + jet and Dijet MC");
       leg->AddEntry(hIsoPho_1[ieta][ipt],"genMomId==22");
       leg->AddEntry(hIsoPho_2[ieta][ipt],"genMomId<=22");
       leg->AddEntry(hIsoPho_3[ieta][ipt],"genMomId==22 && genCalIsoDR04<5");
       leg->AddEntry(hIsoPho_4[ieta][ipt],"genMomId<=22 && genCalIsoDR04<5");
       leg->AddEntry(hIsoPho_5[ieta][ipt],"genCalIsoDR04<5");
       leg->Draw("same");
    
       sprintf(tmp,"/mc/QCD_mess/figures/%s_Template_Et_%d_%d_Eta_%.1f_%.1f",outputName.data(),
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
       hIsoPho[ieta][ipt]->Write();
       hIsoPho_1[ieta][ipt]->Write();
       hIsoPho_2[ieta][ipt]->Write();
       hIsoPho_3[ieta][ipt]->Write();
       hIsoPho_4[ieta][ipt]->Write();
       hIsoPho_5[ieta][ipt]->Write();

       hIsoJet[ieta][ipt]->Write();
       hIsoMixSig[ieta][ipt]->Write();
       hIsoMixBkg[ieta][ipt]->Write();
       hIsoMixData[ieta][ipt]->Write();

     } // end of looping over pt bins
   } // end of looping over eta bins

   outFile->Close();

}
