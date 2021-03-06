/*==============================================================================================
  Use Rong-Shyang's style histograms to call template binned likelihood fitter, unbined fitter, 
 two-bin fitter and compare results with MC truth. In addition, signal yields are saved in 
 txt and tex formats.
Get purity for isolation < 5 GeV only
  .L call_allfitters_rs_iso5GeV.C++

  (1) if fit data, and get corrected yield:

  call_allfitters_rs_iso5GeV(true, true, luminosity in nb)

  // template MC samples were produced normalized to 100/nb

  (2) if fit fake data from MC and get corrected yield:
 
  call_allfitters_rs_iso5GeV(false, true)


 Shin-Shan Eiko Yu, 2010.06.01

  ==============================================================================================*/



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
const int nPtBin = sizeof(fBinsPt)/sizeof(fBinsPt[0])-1;
const int nEtaBin = (sizeof(fBinsEta)/sizeof(fBinsEta[0]))/2;

void ratioErr(Double_t n1, Double_t n1err, Double_t n2, Double_t n2err,
	      Double_t& ratio, Double_t& err)
{
  if(fabs(n1+n2)<1e-6){ratio=0;err=0;return;}
  ratio = (n1)/(n1+n2);
  
  err= pow(1/(n1+n2) - n1/(n1+n2)/(n1+n2),2)*n1err*n1err+
    pow(n1/(n1+n2)/(n1+n2),2)*n2err*n2err; 

  err = sqrt(err);

}


void call_allfitters_rs_iso5GeV(bool fitData=false, bool doEffCorr=false, double lumi=8.02){
	
  cout << "do efficiency correction = " << doEffCorr << endl;
  cout << "fit data = " << fitData << endl;
  std::string isDataName = fitData? "data": "MC";


  setTDRStyle();
  // settings for purity TGraphAsymmetryErrors
  double fBinsPtMidPoint[nPtBin]={0};
  double fBinsPtError[nPtBin]={0};
  
  for(int ipt=0; ipt < nPtBin; ipt++)
    {
      fBinsPtMidPoint[ipt] = 0.5*(fBinsPt[ipt+1]+fBinsPt[ipt]);
      fBinsPtError[ipt] = 0.5*(fBinsPt[ipt+1]-fBinsPt[ipt]);
    }

  // first get efficiency
  double eff[nEtaBin][nPtBin]={{0}};
  double efferr[nEtaBin][nPtBin]={{0}};

  TFile* inf_eff = new TFile("Efficiency_allMC.root");
  TH1F* phoEff[nEtaBin];
  char tmp[1000];
  
  for(int i=0;i<nEtaBin;i++)
    {
      sprintf(tmp,"phoEff_Pt%d",i);
      phoEff[i] = (TH1F*)inf_eff->FindObjectAny(tmp);
      for(int ipt=0; ipt<nPtBin;ipt++)
	{
	  eff[i][ipt]    =phoEff[i]->GetBinContent(ipt+1);
	  efferr[i][ipt] =phoEff[i]->GetBinError(ipt+1); 
	  cout << "Efficiency in ieta=" << i << ", ipt = " << ipt << " bin = " << eff[i][ipt] << 
	    " +- " << efferr[i][ipt] << endl;
	}
   
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
  double nsig_2bin[nEtaBin][nPtBin]={{0}};
  double nsig_err_2bin[nEtaBin][nPtBin]={{0}};


  TH1F* hTemplate_S[nEtaBin][nPtBin];
  TH1F* hTemplate_B[nEtaBin][nPtBin];
  TH1F* hTemplate_data[nEtaBin][nPtBin];
  TH1F* hdata_S[nEtaBin][nPtBin];
  TH1F* hdata_B[nEtaBin][nPtBin];
  TH1F* hdata_data[nEtaBin][nPtBin];


  std::string dataFile     = fitData? "template_comb3Iso_template_0531.root":"template_comb3Iso_test.root";
  std::string templateFile = fitData? "template_comb3Iso_template_0531.root":"template_comb3Iso_template.root";

  TFile* inf_data = new TFile(dataFile.data());
  TFile* inf_template = new TFile(templateFile.data());
  
  char* dec[2] = {"EB","EE"};
  for(int ieta = 0; ieta < nEtaBin; ieta++){
     for(int ipt=0; ipt < nPtBin; ipt++){
       // getting histograms from data root file
       sprintf(tmp,"h_%s_comb3Iso_EGdata_pt%d",dec[ieta],(int)fBinsPt[ipt]);
       cout << "looking for histogram " << tmp << " in file " << 
	 inf_data->GetName() << endl;
       hdata_data[ieta][ipt] = (TH1F*)inf_data->FindObjectAny(tmp);
       hdata_data[ieta][ipt]->Rebin(6);

       sprintf(tmp,"h_%s_comb3Iso_sig_pt%d",dec[ieta],(int)fBinsPt[ipt]);
       cout << "looking for histogram " << tmp << " in file " << 
	 inf_data->GetName() << endl;
       hdata_S[ieta][ipt] = (TH1F*)inf_data->FindObjectAny(tmp);
       hdata_S[ieta][ipt]->Rebin(6);

       sprintf(tmp,"h_%s_comb3Iso_bkg_pt%d",dec[ieta],(int)fBinsPt[ipt]);
       cout << "looking for histogram " << tmp << " in file " << 
	 inf_data->GetName() << endl;
       hdata_B[ieta][ipt] = (TH1F*)inf_data->FindObjectAny(tmp);
       hdata_B[ieta][ipt]->Rebin(6);


       // getting histogram for template root file
       sprintf(tmp,"h_%s_comb3Iso_EGdata_pt%d",dec[ieta],(int)fBinsPt[ipt]);
       cout << "looking for histogram " << tmp << " in file " << 
	 inf_template->GetName() << endl;
       hTemplate_data[ieta][ipt] = (TH1F*)inf_template->FindObjectAny(tmp);
       hTemplate_data[ieta][ipt]->Rebin(6);

       sprintf(tmp,"h_%s_comb3Iso_sig_pt%d",dec[ieta],(int)fBinsPt[ipt]);
       cout << "looking for histogram " << tmp << " in file " << 
	 inf_template->GetName() << endl;
       hTemplate_S[ieta][ipt] = (TH1F*)inf_template->FindObjectAny(tmp);
       hTemplate_S[ieta][ipt]->Rebin(6);

       sprintf(tmp,"h_%s_comb3Iso_bkg_pt%d",dec[ieta],(int)fBinsPt[ipt]);
       cout << "looking for histogram " << tmp << " in file " << 
	 inf_template->GetName() << endl;
       hTemplate_B[ieta][ipt] = (TH1F*)inf_template->FindObjectAny(tmp);
       hTemplate_B[ieta][ipt]->Rebin(6);

     }
  }



  for(int ieta = 0; ieta < nEtaBin; ieta++){
    for(int ipt=0; ipt < nPtBin; ipt++){

      double scaleEff = doEffCorr? eff[ieta][ipt]: 1.0;

      // 1st, get MC truth
      nsig_mc[ieta][ipt]  = fitData? 
	hTemplate_S[ieta][ipt]->Integral(1,10)/scaleEff * (lumi/100.0):
	hdata_S[ieta][ipt]->Integral(1,10)/scaleEff;

      nbkg_mc[ieta][ipt]  = fitData? 
	hTemplate_B[ieta][ipt]->Integral(1,10)/scaleEff * (lumi/100.0):
	hdata_B[ieta][ipt]->Integral(1,10)/scaleEff;

      if(nsig_mc[ieta][ipt]+nbkg_mc[ieta][ipt] < 1e-6)continue;
      mcpurity[ieta][ipt] = nsig_mc[ieta][ipt]/(nsig_mc[ieta][ipt]+nbkg_mc[ieta][ipt]);

      // 2nd, get unbinned fit
      Double_t* FuncFitResult;
      if(fitData)
	FuncFitResult = Ifit(hdata_data[ieta][ipt],
			     hTemplate_S[ieta][ipt],
			     hTemplate_B[ieta][ipt],1,"EGdata_comb3Iso_et_0531.dat",
			     fBinsEta[ieta*2],fBinsEta[ieta*2+1],
			     fBinsPt[ipt],fBinsPt[ipt+1]);	
      else
	FuncFitResult = Ifit(hdata_data[ieta][ipt],
			     hTemplate_S[ieta][ipt],
			     hTemplate_B[ieta][ipt],0,"EGdata_comb3Iso_et_0531.dat",
			     fBinsEta[ieta*2],fBinsEta[ieta*2+1],
			     fBinsPt[ipt],fBinsPt[ipt+1]);	
 
      Double_t nsigfunc    = FuncFitResult[4]/scaleEff;
      Double_t errnsigfunc = FuncFitResult[5]/scaleEff;
      Double_t nbkgfunc    = FuncFitResult[6]/scaleEff;
      Double_t errnbkgfunc = FuncFitResult[7]/scaleEff;

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

 
      Double_t nsigtemplate    = TemplateFitResult[4]/scaleEff;
      Double_t errnsigtemplate = TemplateFitResult[5]/scaleEff;
      Double_t nbkgtemplate    = TemplateFitResult[6]/scaleEff;
      Double_t errnbkgtemplate = TemplateFitResult[7]/scaleEff;

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

      Double_t ndata = hdata_data[ieta][ipt]->Integral(1,10);
      purity_2bin[ieta][ipt]     = TwoBinResult[2];
      purity_err_2bin[ieta][ipt] = TwoBinResult[3];
      nsig_2bin[ieta][ipt]       = purity_2bin[ieta][ipt]*ndata/scaleEff;
      nsig_err_2bin[ieta][ipt]   = purity_err_2bin[ieta][ipt]*ndata/scaleEff;

    } // end of loop over pt bins
  } // end of loop over eta bins


  // printing out the results
  ofstream fout;
  fout.open("yield_iso5GeV.txt");

  ofstream fout2;
  fout2.open("purity_iso5GeV.txt");
 
  ofstream texout;
  texout.open("yield_iso5GeV.tex");

  ofstream texout2;
  texout2.open("purity_iso5GeV.tex");
 
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

      fout << Form("%d:%d",(int)fBinsPt[ipt], (int)fBinsPt[ipt+1]) << " ";
      fout << Form("& %.1f +- %.1f", nsig_mc[ieta][ipt],0.) << " ";

      fout2 << Form("%d:%d",(int)fBinsPt[ipt], (int)fBinsPt[ipt+1]) << " ";
      fout2 << Form("& %.3f +- %.3f", mcpurity[ieta][ipt],0.) << " ";

      texout << Form("%d:%d",(int)fBinsPt[ipt], (int)fBinsPt[ipt+1]) << " ";
      texout << Form("& %.1f $\\pm$ %.1f", nsig_mc[ieta][ipt],0.) << " ";

      texout2 << Form("%d:%d",(int)fBinsPt[ipt], (int)fBinsPt[ipt+1]) << " ";
      texout2 << Form("& %.3f $\\pm$ %.3f", mcpurity[ieta][ipt],0.) << " ";

      cout << " Unbinned fit result : " << endl;
      cout << " purity = " << purity_func[ieta][ipt] << " +- " << purity_err_func[ieta][ipt] << endl;
      cout << " nsig = " <<   nsig_func[ieta][ipt] << " +- " << nsig_err_func[ieta][ipt] << endl;
      cout << " nbkg = " <<   nbkg_func[ieta][ipt] << " +- " << nbkg_err_func[ieta][ipt] << endl;
      cout << endl;

      fout  << Form("& %.1f +- %.1f", nsig_func[ieta][ipt], nsig_err_func[ieta][ipt]) << " ";
      fout2 << Form("& %.3f +- %.3f", purity_func[ieta][ipt], purity_err_func[ieta][ipt]) << " ";

      texout << Form("& %.1f $\\pm$ %.1f", nsig_func[ieta][ipt], nsig_err_func[ieta][ipt]) << " ";
      texout2 << Form("& %.3f $\\pm$ %.3f", purity_func[ieta][ipt], purity_err_func[ieta][ipt]) << " ";

      cout << " Template fit result : " << endl;
      cout << " purity = " << purity_template[ieta][ipt] << " +- " << purity_err_template[ieta][ipt] << endl;
      cout << " nsig = " <<   nsig_template[ieta][ipt] << " +- " << nsig_err_template[ieta][ipt] << endl;
      cout << " nbkg = " <<   nbkg_template[ieta][ipt] << " +- " << nbkg_err_template[ieta][ipt] << endl;
      cout << endl;

      fout << Form("& %.1f +- %.1f", nsig_template[ieta][ipt], nsig_err_template[ieta][ipt]) << " ";
      fout2 << Form("& %.3f +- %.3f", purity_template[ieta][ipt], purity_err_template[ieta][ipt]) << " ";

      texout << Form("& %.1f $\\pm$ %.1f", nsig_template[ieta][ipt], nsig_err_template[ieta][ipt]) << " ";
      texout2 << Form("& %.3f $\\pm$ %.3f", purity_template[ieta][ipt], purity_err_template[ieta][ipt]) << " ";

      cout << " Two-bin result : " << endl;
      cout << " purity = " << purity_2bin[ieta][ipt] << " +- " << purity_err_2bin[ieta][ipt] << endl;
      cout << " nsig = " <<   nsig_2bin[ieta][ipt] << " +- " << nsig_err_2bin[ieta][ipt] << endl;

      fout << Form("& %.1f +- %.1f", nsig_2bin[ieta][ipt], nsig_err_2bin[ieta][ipt]) << " ";
      fout2 << Form("& %.3f +- %.3f", purity_2bin[ieta][ipt], purity_err_2bin[ieta][ipt]) << " ";

      texout << "\\\\" << endl;
      fout << endl;

      texout2 << "\\\\" << endl;
      fout2 << endl;

    }
  }

  texout.close();
  fout.close();

  texout2.close();
  fout2.close();

  
  TCanvas *c1 = new TCanvas("c1","",600,600);
  c1->Draw();
  TH2F *h2 = new TH2F("h2","",100, 0., 140, 100, 0, 1);
  h2->SetNdivisions(505,"XY");
  h2->SetXTitle("p_{T} (GeV)");
  h2->SetYTitle("purity for iso < 5 GeV");
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

  TLegend *tleg = new TLegend(0.15, 0.65, 0.55, 0.92);
  tleg->SetHeader("EB Photons");
  tleg->SetFillColor(0);
  tleg->SetFillStyle(0);
  tleg->SetShadowColor(0);
  tleg->SetBorderSize(0);
  tleg->AddEntry(tg_mc_purity_EB,"MC truth","p");
  tleg->AddEntry(tgrs_purity_EB,"unbinned fit","pl");
  tleg->AddEntry(tgrs_purity_template_EB,"template fit","pl");
  tleg->AddEntry(tgrs_purity_twobin_EB,"two bin","pl");
  tleg->Draw();

  c1->SaveAs(Form("plots/%s_allfitter_purity_EB_iso5GeV.eps",isDataName.data()));
  c1->SaveAs(Form("plots/%s_allfitter_purity_EB_iso5GeV.gif",isDataName.data()));
  c1->SaveAs(Form("plots/%s_allfitter_purity_EB_iso5GeV.C",isDataName.data()));

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
  tleg2->SetFillStyle(0);
  tleg2->SetShadowColor(0);
  tleg2->SetBorderSize(0);
  tleg2->AddEntry(tg_mc_purity_EE,"MC truth","p");
  tleg2->AddEntry(tgrs_purity_EE,"unbinned fit","pl");
  tleg2->AddEntry(tgrs_purity_template_EE,"template fit","pl");
  tleg2->AddEntry(tgrs_purity_twobin_EE,"two bin","pl");
  tleg2->Draw();

  c2->SaveAs(Form("plots/%s_allfitter_purity_EE_iso5GeV.eps",isDataName.data()));
  c2->SaveAs(Form("plots/%s_allfitter_purity_EE_iso5GeV.gif",isDataName.data()));
  c2->SaveAs(Form("plots/%s_allfitter_purity_EE_iso5GeV.C",isDataName.data()));


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

  TGraphAsymmErrors *tgrs_sig_twobin_EB = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, nsig_2bin[0], 
								fBinsPtError, fBinsPtError, 
								nsig_err_2bin[0],nsig_err_2bin[0]);
  tgrs_sig_twobin_EB->SetName("tgrs_sig_twobin_EB");
  tgrs_sig_twobin_EB->SetMarkerColor(kGreen);
  tgrs_sig_twobin_EB->SetMarkerStyle(3);
  tgrs_sig_twobin_EB->SetLineColor(kGreen);
  tgrs_sig_twobin_EB->SetLineWidth(2);
  tgrs_sig_twobin_EB->Draw("p e same");

  TLegend *tleg3 = new TLegend(0.193, 0.168, 0.642, 0.437);
  tleg3->SetHeader("EB Photons");
  tleg3->SetFillColor(0);
  tleg3->SetFillStyle(0);
  tleg3->SetShadowColor(0);
  tleg3->SetBorderSize(0);
  tleg3->AddEntry(tg_mc_sig_EB,"MC truth","p");
  tleg3->AddEntry(tgrs_sig_EB,"unbinned fit","pl");
  tleg3->AddEntry(tgrs_sig_template_EB,"template fit","pl");
  tleg3->AddEntry(tgrs_sig_twobin_EB,"two bin","pl");
  tleg3->Draw();

  c3->SaveAs(Form("plots/%s_allfitter_yield_EB_iso5GeV.eps",isDataName.data()));
  c3->SaveAs(Form("plots/%s_allfitter_yield_EB_iso5GeV.gif",isDataName.data()));
  c3->SaveAs(Form("plots/%s_allfitter_yield_EB_iso5GeV.C",isDataName.data()));

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

  TGraphAsymmErrors *tgrs_sig_twobin_EE = new TGraphAsymmErrors(nPtBin, fBinsPtMidPoint, nsig_2bin[1], 
							   fBinsPtError, fBinsPtError, 
							   nsig_err_2bin[1],nsig_err_2bin[1]);
  tgrs_sig_twobin_EE->SetName("tgrs_sig_twobin_EE");
  tgrs_sig_twobin_EE->SetMarkerColor(kGreen);
  tgrs_sig_twobin_EE->SetMarkerStyle(3);
  tgrs_sig_twobin_EE->SetLineColor(kGreen);
  tgrs_sig_twobin_EE->SetLineWidth(2);
  tgrs_sig_twobin_EE->Draw("p e same");


  TLegend *tleg4 = new TLegend(0.193, 0.168, 0.642, 0.437);
  tleg4->SetHeader("EE Photons");
  tleg4->SetFillColor(0);
  tleg4->SetFillStyle(0);
  tleg4->SetShadowColor(0);
  tleg4->SetBorderSize(0);
  tleg4->AddEntry(tg_mc_sig_EE,"MC truth","p");
  tleg4->AddEntry(tgrs_sig_EE,"unbinned fit","pl");
  tleg4->AddEntry(tgrs_sig_template_EE,"template fit","pl");
  tleg4->AddEntry(tgrs_sig_twobin_EE,"two bin","pl");
  tleg4->Draw();

  c4->SaveAs(Form("plots/%s_allfitter_yield_EE_iso5GeV.eps",isDataName.data()));
  c4->SaveAs(Form("plots/%s_allfitter_yield_EE_iso5GeV.gif",isDataName.data()));
  c4->SaveAs(Form("plots/%s_allfitter_yield_EE_iso5GeV.C",isDataName.data()));


  TFile* outFile = new TFile("isotemplate_yield_rongshyang_iso5GeV.root","recreate"); 
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
