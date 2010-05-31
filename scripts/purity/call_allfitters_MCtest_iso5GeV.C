#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TH1.h>
#include <string>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h> 
#include <TCanvas.h>
#include <TLegend.h>
#include "purity_twobin.C"
#include "ext-mL_fit.C"
#include "ext-mL_fit_unbinned.C"
#include "setTDRStyle.C"

using namespace std;

// the pt and eta binning
const double fBinsEta[]={0,1.45,1.7,2.5};
const double fBinsPt[]={15,20,30,50,80,120};
// const int nPtBin = sizeof(fBinsPt)/sizeof(fBinsPt[0])-1;
// const int nEtaBin = (sizeof(fBinsEta)/sizeof(fBinsEta[0]))/2;
const int nPtBin = 5;
const int nEtaBin = 2;

void ratioErr(Double_t n1, Double_t n1err, Double_t n2, Double_t n2err,
	      Double_t& ratio, Double_t& err)
{
  if(fabs(n1+n2)<1e-6){ratio=0;err=0;return;}
  ratio = (n1)/(n1+n2);
  
  err= pow(1/(n1+n2) - n1/(n1+n2)/(n1+n2),2)*n1err*n1err+
    pow(n1/(n1+n2)/(n1+n2),2)*n2err*n2err; 

  err = sqrt(err);

}


void call_allfitters_MCtest_iso5GeV(){

  setTDRStyle();
     // settings for purity TGraphAsymmetryErrors
  double fBinsPtMidPoint[nPtBin]={0};
  double fBinsPtError[nPtBin]={0};
  
  for(int ipt=0; ipt < nPtBin; ipt++)
    {
      fBinsPtMidPoint[ipt] = 0.5*(fBinsPt[ipt+1]+fBinsPt[ipt]);
      fBinsPtError[ipt] = 0.5*(fBinsPt[ipt+1]-fBinsPt[ipt]);
    }


  double mcpurity[nEtaBin][nPtBin]={{0}};
  double nsig_mc[nEtaBin][nPtBin]={{0}};
  double nbkg_mc[nEtaBin][nPtBin]={{0}};


  double purity_func[nEtaBin][nPtBin]={{0}};
  double purity_err_func[nEtaBin][nPtBin]={{0}};
  double nsig_func[nEtaBin][nPtBin]={{0}};
  double nbkg_func[nEtaBin][nPtBin]={{0}};
  double nsig_err_func[nEtaBin][nPtBin]={{0}};
  double nbkg_err_func[nEtaBin][nPtBin]={{0}};

  double purity_template[nEtaBin][nPtBin]={{0}};
  double purity_err_template[nEtaBin][nPtBin]={{0}};
  double nsig_template[nEtaBin][nPtBin]={{0}};
  double nbkg_template[nEtaBin][nPtBin]={{0}};
  double nsig_err_template[nEtaBin][nPtBin]={{0}};
  double nbkg_err_template[nEtaBin][nPtBin]={{0}};
 
  double purity_2bin[nEtaBin][nPtBin]={{0}};
  double purity_err_2bin[nEtaBin][nPtBin]={{0}};

  // first try two bin results
  TH1F* hTemplate_S[nEtaBin][nPtBin];
  TH1F* hTemplate_B[nEtaBin][nPtBin];
  TH1F* hTemplate_data[nEtaBin][nPtBin];
  TH1F* hdata_S[nEtaBin][nPtBin];
  TH1F* hdata_B[nEtaBin][nPtBin];
  TH1F* hdata_data[nEtaBin][nPtBin];


  TFile *inf_data = new TFile("/mc/QCD_mess/pcncu18_figures/fakedata_sumIso_histo.root");

  //   TFile *inf_template = new TFile("/mc/QCD_mess/pcncu18_figures/template_sumIso_histo.root");
  TFile *inf_template = new TFile("/mc/QCD_mess/pcncu18_figures/fakedata_sumIso_histo.root");


  char tmp[1000];
  for(int ieta = 0; ieta < nEtaBin; ieta++){
     for(int ipt=0; ipt < nPtBin; ipt++){

       // getting histograms from data root file
       sprintf(tmp,"Template_Et_%d_%d_Eta_%.2f_%.2f",
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
       cout << "looking for histogram " << tmp << " in file " << 
	 inf_data->GetName() << endl;
       hdata_data[ieta][ipt] = (TH1F*)inf_data->FindObjectAny(tmp);

       sprintf(tmp,"TemplateS_Et_%d_%d_Eta_%.2f_%.2f",
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
       cout << "looking for histogram " << tmp << " in file " << 
	 inf_data->GetName() << endl;
       hdata_S[ieta][ipt] = (TH1F*)inf_data->FindObjectAny(tmp);

       sprintf(tmp,"TemplateB_Et_%d_%d_Eta_%.2f_%.2f",
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
       cout << "looking for histogram " << tmp << " in file " << 
	 inf_data->GetName() << endl;
       hdata_B[ieta][ipt] = (TH1F*)inf_data->FindObjectAny(tmp);


       // getting histogram for template root file
       sprintf(tmp,"Template_Et_%d_%d_Eta_%.2f_%.2f",
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
       cout << "looking for histogram " << tmp << " in file " << 
	 inf_template->GetName() << endl;
       hTemplate_data[ieta][ipt] = (TH1F*)inf_template->FindObjectAny(tmp);

       sprintf(tmp,"TemplateS_Et_%d_%d_Eta_%.2f_%.2f",
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
       cout << "looking for histogram " << tmp << " in file " << 
	 inf_template->GetName() << endl;
       hTemplate_S[ieta][ipt] = (TH1F*)inf_template->FindObjectAny(tmp);

       sprintf(tmp,"TemplateB_Et_%d_%d_Eta_%.2f_%.2f",
	       (int)fBinsPt[ipt], (int)fBinsPt[ipt+1],
	       fBinsEta[ieta*2],fBinsEta[ieta*2+1]);
       cout << "looking for histogram " << tmp << " in file " << 
	 inf_template->GetName() << endl;
       hTemplate_B[ieta][ipt] = (TH1F*)inf_template->FindObjectAny(tmp);

     }
  }



  for(int ieta = 0; ieta < nEtaBin; ieta++){
    for(int ipt=0; ipt < nPtBin; ipt++){
  

      // 1st, get MC truth
      nsig_mc[ieta][ipt]  = hdata_S[ieta][ipt]->Integral(1,10);
      nbkg_mc[ieta][ipt]  = hdata_B[ieta][ipt]->Integral(1,10);
      if(nsig_mc[ieta][ipt]+nbkg_mc[ieta][ipt] < 1e-6)continue;
      mcpurity[ieta][ipt] = nsig_mc[ieta][ipt]/(nsig_mc[ieta][ipt]+nbkg_mc[ieta][ipt]);

      // 2nd, get unbinned fit
      Double_t* FuncFitResult;
      FuncFitResult = Ifit(hdata_data[ieta][ipt],
			   hTemplate_S[ieta][ipt],
			   hTemplate_B[ieta][ipt]);

 
      Double_t nsigfunc    = FuncFitResult[4];
      Double_t errnsigfunc = FuncFitResult[5];
      Double_t nbkgfunc    = FuncFitResult[6];
      Double_t errnbkgfunc = FuncFitResult[7];

      Double_t purityfunc = 0;
      Double_t errpurityfunc = 0;
 
      ratioErr(nsigfunc, errnsigfunc, nbkgfunc, errnbkgfunc,
	       purityfunc, errpurityfunc);

      purity_func[ieta][ipt]    = purityfunc;
      purity_err_func[ieta][ipt]= errpurityfunc;
      nsig_func[ieta][ipt]      = nsigfunc;
      nsig_err_func[ieta][ipt]  = errnsigfunc;
      nbkg_func[ieta][ipt]      = nbkgfunc;
      nbkg_err_func[ieta][ipt]  = errnbkgfunc;


      // 3rd, get template fit result
      Double_t* TemplateFitResult;
      TemplateFitResult = IfitBin(hdata_data[ieta][ipt],
				  hTemplate_S[ieta][ipt],
				  hTemplate_B[ieta][ipt]);

 
      Double_t nsigtemplate    = TemplateFitResult[4];
      Double_t errnsigtemplate = TemplateFitResult[5];
      Double_t nbkgtemplate    = TemplateFitResult[6];
      Double_t errnbkgtemplate = TemplateFitResult[7];

      Double_t puritytemplate = 0;
      Double_t errpuritytemplate = 0;
 
      ratioErr(nsigtemplate, errnsigtemplate, nbkgtemplate, errnbkgtemplate,
	       puritytemplate, errpuritytemplate);

      purity_template[ieta][ipt]    = puritytemplate;
      purity_err_template[ieta][ipt]= errpuritytemplate;
      nsig_template[ieta][ipt]      = nsigtemplate;
      nsig_err_template[ieta][ipt]  = errnsigtemplate;
      nbkg_template[ieta][ipt]      = nbkgtemplate;
      nbkg_err_template[ieta][ipt]  = errnbkgtemplate;



      // 4th, get two bin result
      Double_t* TwoBinResult;
      TwoBinResult = purity_twobin(hdata_data[ieta][ipt],
				   hTemplate_S[ieta][ipt],
				   hTemplate_B[ieta][ipt],10);

      purity_2bin[ieta][ipt]     = TwoBinResult[2];
      purity_err_2bin[ieta][ipt] = TwoBinResult[3];

    

    } // end of loop over pt bins
  } // end of loop over eta bins

  
  // printing out the results
  for(int ieta=0; ieta < nEtaBin; ieta++){
    for(int ipt=0; ipt < nPtBin; ipt++){
      cout << endl;
      cout << " ===== for ieta=" << ieta << ", ipt=" << ipt << " ======== " << endl;
      cout << " MC truth : " << endl;
      cout << " purity = " << mcpurity[ieta][ipt] << endl;
      cout << " nsig = " << nsig_mc[ieta][ipt] << endl;
      cout << " nbkg = " << nbkg_mc[ieta][ipt] << endl;
      cout << endl;
      cout << " Unbinned fit result : " << endl;
      cout << " purity = " << purity_func[ieta][ipt] << " +- " << purity_err_func[ieta][ipt] << endl;
      cout << " nsig = " <<   nsig_func[ieta][ipt] << " +- " << nsig_err_func[ieta][ipt] << endl;
      cout << " nbkg = " <<   nbkg_func[ieta][ipt] << " +- " << nbkg_err_func[ieta][ipt] << endl;
      cout << endl;
      cout << " Template fit result : " << endl;
      cout << " purity = " << purity_template[ieta][ipt] << " +- " << purity_err_template[ieta][ipt] << endl;
      cout << " nsig = " <<   nsig_template[ieta][ipt] << " +- " << nsig_err_template[ieta][ipt] << endl;
      cout << " nbkg = " <<   nbkg_template[ieta][ipt] << " +- " << nbkg_err_template[ieta][ipt] << endl;
      cout << endl;
      cout << " Two-bin result : " << endl;
      cout << " purity = " << purity_2bin[ieta][ipt] << " +- " << purity_err_2bin[ieta][ipt] << endl;

 
    }
  }
 
  
  TCanvas *c1 = new TCanvas("c1","",600,600);
  c1->Draw();
  TH2F *h2 = new TH2F("h2","",100, 0., 140, 100, 0, 1);
  h2->SetNdivisions(505,"XY");
  h2->SetXTitle("p_{T} (GeV)");
  h2->SetYTitle("purity");
  h2->Draw();

  TGraph *tg_mc_purity_EB = new TGraph(nPtBin, fBinsPtMidPoint, mcpurity[0]);
  tg_mc_purity_EB->SetName("tg_mc_purity_EB");
  tg_mc_purity_EB->SetMarkerSize(1.);
  tg_mc_purity_EB->Draw("p same");
  TGraphAsymmErrors *tgrs_purity_EB = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, 
							    purity_func[0], 
							    fBinsPtError, fBinsPtError, 
							    purity_err_func[0], purity_err_func[0]);
  tgrs_purity_EB->SetName("tgrs_purity_EB");
  tgrs_purity_EB->SetMarkerColor(4);
  tgrs_purity_EB->SetMarkerStyle(25);
  tgrs_purity_EB->SetLineColor(4);
  tgrs_purity_EB->SetLineWidth(2);
  tgrs_purity_EB->Draw("p e same");

  TGraphAsymmErrors *tgrs_purity_template_EB = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, 
								     purity_template[0], 
								     fBinsPtError,fBinsPtError, 
								     purity_err_template[0],
								     purity_err_template[0]);
  tgrs_purity_template_EB->SetName("tgrs_purity_template_EB");
  tgrs_purity_template_EB->SetMarkerColor(6);
  tgrs_purity_template_EB->SetMarkerStyle(26);
  tgrs_purity_template_EB->SetLineColor(6);
  tgrs_purity_template_EB->SetLineWidth(2);
  tgrs_purity_template_EB->Draw("p e same");

  TGraphAsymmErrors *tgrs_purity_twobin_EB = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, 
								   purity_2bin[0], 
								   fBinsPtError,fBinsPtError, 
								   purity_err_2bin[0],
								   purity_err_2bin[0]);
  tgrs_purity_twobin_EB->SetName("tgrs_purity_twobin_EB");
  tgrs_purity_twobin_EB->SetMarkerColor(kGreen);
  tgrs_purity_twobin_EB->SetMarkerStyle(3 );
  tgrs_purity_twobin_EB->SetLineColor(kGreen);
  tgrs_purity_twobin_EB->SetLineWidth(2);
  tgrs_purity_twobin_EB->Draw("p e same");

  tgrs_purity_template_EB->Draw("p e same");
  tg_mc_purity_EB->Draw("p same");

  TLegend *tleg = new TLegend(0.51, 0.17, 0.91, 0.44);
  char text[50];  
  tleg->SetHeader("EB Photons");
  tleg->SetFillColor(0);
  tleg->SetShadowColor(0);
  tleg->SetBorderSize(0);
  tleg->AddEntry(tg_mc_purity_EB,"MC truth","p");
  tleg->AddEntry(tgrs_purity_EB,"unbinned fit","pl");
  tleg->AddEntry(tgrs_purity_template_EB,"template fit","pl");
  tleg->AddEntry(tgrs_purity_twobin_EB,"two bin","pl");
  tleg->Draw();

  c1->SaveAs("plots/allfitter_purity_iso5GeV_EB.eps");
  c1->SaveAs("plots/allfitter_purity_iso5GeV_EB.gif"); 
  c1->SaveAs("plots/allfitter_purity_iso5GeV_EB.pdf");
  c1->SaveAs("plots/allfitter_purity_iso5GeV_EB.C");
  

  TCanvas *c2 = new TCanvas("c2","",600,600);
  c2->Draw();
  h2->Draw();
  TGraph *tg_mc_purity_EE = new TGraph(nPtBin, fBinsPtMidPoint, mcpurity[1]);
  tg_mc_purity_EE->SetName("tg_mc_purity_EE");
  tg_mc_purity_EE->SetMarkerSize(1.);
  tg_mc_purity_EE->Draw("p same");

  TGraphAsymmErrors *tgrs_purity_EE = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, 
							    purity_func[1], 
							    fBinsPtError, fBinsPtError, 
							    purity_err_func[1],purity_err_func[1]);
  tgrs_purity_EE->SetName("tgrs_purity_EE");
  tgrs_purity_EE->SetMarkerColor(4);
  tgrs_purity_EE->SetMarkerStyle(25);
  tgrs_purity_EE->SetLineColor(4);
  tgrs_purity_EE->SetLineWidth(2);
  tgrs_purity_EE->Draw("p e same");

  TGraphAsymmErrors *tgrs_purity_template_EE = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, 
								     purity_template[1], 
								     fBinsPtError,fBinsPtError, 
								     purity_err_template[1],
								     purity_err_template[1]);
  tgrs_purity_template_EE->SetName("tgrs_purity_template_EE");
  tgrs_purity_template_EE->SetMarkerColor(6);
  tgrs_purity_template_EE->SetMarkerStyle(26);
  tgrs_purity_template_EE->SetLineColor(6);
  tgrs_purity_template_EE->SetLineWidth(2);
  tgrs_purity_template_EE->Draw("p e same");

  TGraphAsymmErrors *tgrs_purity_twobin_EE = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, 
								   purity_2bin[1], 
								   fBinsPtError, fBinsPtError,
								   purity_err_2bin[1], 
								   purity_err_2bin[1]);
  tgrs_purity_twobin_EE->SetName("tgrs_purity_twobin_EE");
  tgrs_purity_twobin_EE->SetMarkerColor(kGreen);
  tgrs_purity_twobin_EE->SetMarkerStyle(3);
  tgrs_purity_twobin_EE->SetLineColor(kGreen);
  tgrs_purity_twobin_EE->SetLineWidth(2);
  tgrs_purity_twobin_EE->Draw("p e same");

  tgrs_purity_template_EE->Draw("p e same");
  tg_mc_purity_EE->Draw("p same");

  TLegend *tleg2 = new TLegend(0.15, 0.65, 0.55, 0.92);
  tleg2->SetHeader("EE Photons");
  tleg2->SetFillColor(0);
  tleg2->SetShadowColor(0);
  tleg2->SetBorderSize(0);
  tleg2->AddEntry(tg_mc_purity_EE,"MC truth","p");
  tleg2->AddEntry(tgrs_purity_EE,"unbinned fit","pl");
  tleg2->AddEntry(tgrs_purity_template_EE,"template fit","pl");
  tleg2->AddEntry(tgrs_purity_twobin_EE,"two bin","pl");
  tleg2->Draw();

  c2->SaveAs("plots/allfitter_purity_iso5GeV_EE.eps");
  c2->SaveAs("plots/allfitter_purity_iso5GeV_EE.gif");
  c2->SaveAs("plots/allfitter_purity_iso5GeV_EE.pdf");
  c2->SaveAs("plots/allfitter_purity_iso5GeV_EE.C");


  TCanvas *c3 = new TCanvas("c3","",600,600);
  c3->Draw();
  gPad->SetLogy();
  TH2F *h3 = new TH2F("h3","",100, 0., 140, 100, 0.1, 10000);  
  h3->SetNdivisions(505,"XY");
  h3->SetXTitle("p_{T} (GeV)");
  h3->SetYTitle("Signal yield / bin");
  h3->Draw();
  TGraph *tg_mc_sig_EB = new TGraph(nPtBin, fBinsPtMidPoint, nsig_mc[0]);
  tg_mc_sig_EB->SetName("tg_mc_sig_EB");
  tg_mc_sig_EB->SetMarkerSize(1);
  tg_mc_sig_EB->Draw("p same");

  TGraphAsymmErrors *tgrs_sig_EB = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, nsig_func[0], 
							 fBinsPtError, fBinsPtError, 
							 nsig_err_func[0], nsig_err_func[0]);
  tgrs_sig_EB->SetName("tgrs_sig_EB");
  tgrs_sig_EB->SetMarkerColor(4);
  tgrs_sig_EB->SetMarkerStyle(25);
  tgrs_sig_EB->SetLineColor(4);
  tgrs_sig_EB->SetLineWidth(2);
  tgrs_sig_EB->Draw("p e same");

  TGraphAsymmErrors *tgrs_sig_template_EB = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, 
								  nsig_template[0], 
								  fBinsPtError,fBinsPtError, 
								  nsig_err_template[0],
								  nsig_err_template[0]);
  tgrs_sig_template_EB->SetName("tgrs_sig_template_EB");
  tgrs_sig_template_EB->SetMarkerColor(6);
  tgrs_sig_template_EB->SetMarkerStyle(26);
  tgrs_sig_template_EB->SetLineColor(6);
  tgrs_sig_template_EB->SetLineWidth(2);
  tgrs_sig_template_EB->Draw("p e same");


  TLegend *tleg3 = new TLegend(0.193, 0.168, 0.642, 0.437);
  tleg3->SetHeader("EB Photons");
  tleg3->SetFillColor(0);
  tleg3->SetShadowColor(0);
  tleg3->SetBorderSize(0);
  tleg3->AddEntry(tg_mc_sig_EB,"MC truth","p");
  tleg3->AddEntry(tgrs_sig_EB,"unbinned fit","pl");
  tleg3->AddEntry(tgrs_sig_template_EB,"template fit","pl");
  tleg3->Draw();

  c3->SaveAs("plots/allfitter_yield_iso5GeV_EB.eps");
  c3->SaveAs("plots/allfitter_yield_iso5GeV_EB.gif");
  c3->SaveAs("plots/allfitter_yield_iso5GeV_EB.pdf");
  c3->SaveAs("plots/allfitter_yield_iso5GeV_EB.C");

  TCanvas *c4 = new TCanvas("c4","",600,600);
  c4->Draw();
  gPad->SetLogy();

  h3->Draw();
  TGraph *tg_mc_sig_EE = new TGraph(nPtBin, fBinsPtMidPoint, nsig_mc[1]);
  tg_mc_sig_EE->SetName("tg_mc_sig_EE");
  tg_mc_sig_EE->SetMarkerSize(1.);
  tg_mc_sig_EE->Draw("p same");

  TGraphAsymmErrors *tgrs_sig_EE = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, nsig_func[1], 
							 fBinsPtError,fBinsPtError, 
							 nsig_err_func[1],nsig_err_func[1]);
  tgrs_sig_EE->SetName("tgrs_sig_EE");
  tgrs_sig_EE->SetMarkerColor(4);
  tgrs_sig_EE->SetMarkerStyle(25);
  tgrs_sig_EE->SetLineColor(4);
  tgrs_sig_EE->SetLineWidth(2);
  tgrs_sig_EE->Draw("p e same");

  TGraphAsymmErrors *tgrs_sig_template_EE = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, 
								  nsig_template[1], 
								  fBinsPtError,fBinsPtError, 
								  nsig_err_template[1],
								  nsig_err_template[1]);
  tgrs_sig_template_EE->SetName("tgrs_sig_template_EE");
  tgrs_sig_template_EE->SetMarkerColor(6);
  tgrs_sig_template_EE->SetMarkerStyle(26);
  tgrs_sig_template_EE->SetLineColor(6);
  tgrs_sig_template_EE->SetLineWidth(2);
  tgrs_sig_template_EE->Draw("p e same");



  TLegend *tleg4 = new TLegend(0.193, 0.168, 0.642, 0.437);
  tleg4->SetHeader("EE Photons");
  tleg4->SetFillColor(0);
  tleg4->SetShadowColor(0);
  tleg4->SetBorderSize(0);
  tleg4->AddEntry(tg_mc_sig_EE,"MC truth","p");
  tleg4->AddEntry(tgrs_sig_EE,"unbinned fit","pl");
  tleg4->AddEntry(tgrs_sig_template_EE,"template fit","pl");
  tleg4->Draw();


  c4->SaveAs("plots/allfitter_yield_iso5GeV_EE.eps");
  c4->SaveAs("plots/allfitter_yield_iso5GeV_EE.gif");
  c4->SaveAs("plots/allfitter_yield_iso5GeV_EE.pdf");
  c4->SaveAs("plots/allfitter_yield_iso5GeV_EE.C");

  TFile* outFile = new TFile("isotemplate_yield_iso5GeV.root","recreate"); 
  tg_mc_purity_EB->Write();
  tgrs_purity_EB->Write();
  tgrs_purity_template_EB->Write();
  tgrs_purity_twobin_EB->Write();

  tg_mc_purity_EE->Write();
  tgrs_purity_EE->Write();
  tgrs_purity_template_EE->Write();
  tgrs_purity_twobin_EE->Write();

  tg_mc_sig_EB->Write();
  tgrs_sig_EB->Write();
  tgrs_sig_template_EB->Write();

  tg_mc_sig_EE->Write();
  tgrs_sig_EE->Write();
  tgrs_sig_template_EE->Write();

  outFile->Close();



}
