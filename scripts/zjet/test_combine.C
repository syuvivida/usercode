#include <TH1.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TFile.h>
#include <iostream>
#include <string>

using namespace std;

void test_combine(std::string eikoName="h_jety",double correlation=1.0, 
		  bool logScale=false)
{

  TH1D* h_e;
  TH1D* h_ejesup;
  TH1D* h_ejesdn;

  TH1F* h_mu;
  TH1F* h_mujesup;
  TH1F* h_mujesdn;
  TH1F* h_mujes;

  TH1D* h_combine;

  // electron channel
  TFile* f_e = TFile::Open("cts_CorrectedPlotsZCut_Jes0_NoBkgSub.root");
  h_e  = (TH1D*)(f_e->Get(eikoName.data()));
  h_e  -> SetName("h_e");
  h_e  -> SetTitle("");
  h_e  -> Sumw2();
  h_e  -> Scale(1.0/h_e->Integral());
  h_e  -> SetYTitle("Arbitrary Unit");
  h_e  -> SetTitleOffset(2.0,"Y");
  h_e  -> GetYaxis()->SetDecimals();
  h_e  -> GetXaxis()->SetDecimals();
  h_e  -> SetLineColor(kBlue-7);
  h_e  -> SetMarkerColor(kBlue-7);
  h_e  -> SetMarkerSize(1);
  h_e  -> SetMarkerStyle(24);

  cout << "h_e integral = " << h_e->Integral() << endl;


  TFile* f_ejesup = TFile::Open("cts_CorrectedPlotsZCut_JesUp_NoBkgSub.root");
  h_ejesup  = (TH1D*)(f_ejesup->Get(eikoName.data()));
  h_ejesup  -> SetName("h_ejesup");
  h_ejesup  -> SetTitle("");
  h_ejesup  -> Sumw2();
  h_ejesup  -> Scale(1.0/h_ejesup->Integral());
  h_ejesup  -> SetYTitle("Arbitrary Unit");
  h_ejesup  -> SetTitleOffset(2.0,"Y");
  h_ejesup  -> GetYaxis()->SetDecimals();
  h_ejesup  -> GetXaxis()->SetDecimals();
  h_ejesup  -> SetLineColor(kBlue-7);
  h_ejesup  -> SetMarkerColor(kBlue-7);
  h_ejesup  -> SetMarkerSize(1);
  h_ejesup  -> SetMarkerStyle(24);


  cout << "h_ejesup integral = " << h_ejesup->Integral() << endl;


  TFile* f_ejesdn = TFile::Open("cts_CorrectedPlotsZCut_JesDn_NoBkgSub.root");
  h_ejesdn  = (TH1D*)(f_ejesdn->Get(eikoName.data()));
  h_ejesdn  -> SetName("h_ejesdn");
  h_ejesdn  -> SetTitle("");
  h_ejesdn  -> Sumw2();
  h_ejesdn  -> Scale(1.0/h_ejesdn->Integral());
  h_ejesdn  -> SetYTitle("Arbitrary Unit");
  h_ejesdn  -> SetTitleOffset(2.0,"Y");
  h_ejesdn  -> GetYaxis()->SetDecimals();
  h_ejesdn  -> GetXaxis()->SetDecimals();
  h_ejesdn  -> SetLineColor(kBlue-7);
  h_ejesdn  -> SetMarkerColor(kBlue-7);
  h_ejesdn  -> SetMarkerSize(1);
  h_ejesdn  -> SetMarkerStyle(24);


  cout << "h_ejesdn integral = " << h_ejesdn->Integral() << endl;

  h_combine= (TH1D*)h_e->Clone("h_combine");
  h_combine  -> Reset();
  h_combine  -> SetTitle("");
  h_combine  -> SetLineColor(1);
  h_combine  -> SetMarkerColor(1);
  h_combine  -> SetMarkerSize(1);
  h_combine  -> SetMarkerStyle(8);

  // muon channel
  TFile* f_mu = TFile::Open("DoubleMu2011_EffCorr_ZpT40_absY_051412.root");

  std::string kengName = "Z1jets_1jeta_BE";
  std::string xtitle   = "|Y(jet)|";

  if(eikoName=="h_ystar")
    {
      kengName = "DEta_per2_Z1jets_BE";
      xtitle = "0.5|Y_{Z}-Y_{jet}|";
    }

  else if(eikoName=="h_yB")
    {
      kengName = "SumEta_per2_Z1jets_BE";
      xtitle = "0.5|Y_{Z}+Y_{jet}|";
    }
  else if(eikoName=="h_jetpt")
    {
      kengName = "Z1jets_1jpt_BE";
      xtitle = "p_{T}(jet) [GeV/c]";
    }
  else if(eikoName=="h_jety")
    {
      kengName = "Z1jets_1jeta_BE";
      xtitle = "|Y(jet)|";
    }
  else if(eikoName=="h_zpt")
    {
      kengName = "dimuonpt1jet_BE";
      xtitle = "p_{T}(Z) [GeV/c]";
    }
  else if(eikoName=="h_zy")
    {
      kengName ="dimuoneta1jet_BE";
      xtitle = "|Y(Z)|";
    }
  

  h_mu = (TH1F*)(f_mu->Get(kengName.data()));
  h_mu -> Sumw2();
  h_mu  -> SetTitle("");
  h_mu -> Scale(1.0/h_mu->Integral());
  h_mu -> SetLineColor(kRed-7);
  h_mu -> SetMarkerColor(kRed-7);
  h_mu -> SetMarkerSize(1);
  h_mu -> SetMarkerStyle(21);

  cout << "h_mu integral = " << h_mu->Integral() << endl;

  TFile* f_jetsys_mu = TFile::Open("DoubleMu2011_JESuncertainty_JetY_061712.root");

  h_mujes = (TH1F*)(f_jetsys_mu->Get(kengName.data()));
  h_mujes -> Sumw2();
  h_mujes -> Scale(1.0/h_mujes->Integral());

  cout << "h_mujes integral = " << h_mujes->Integral() << endl;

  h_mujesup = (TH1F*)(f_jetsys_mu->Get(Form("%sUp",kengName.data())));
  h_mujesup -> Sumw2();
  h_mujesup -> Scale(1.0/h_mujesup->Integral());

  cout << "h_mujesup integral = " << h_mujesup->Integral() << endl;

  h_mujesdn = (TH1F*)(f_jetsys_mu->Get(Form("%sDn",kengName.data())));
  h_mujesdn -> Sumw2();
  h_mujesdn -> Scale(1.0/h_mujesdn->Integral());

  cout << "h_mujesdn integral = " << h_mujesdn->Integral() << endl;


  int nbin = h_e->GetNbinsX(); // electron file has smaller statistical 
                               // uncertainty

  // now loop over bins to compute chi2
  for(int i=1; i<=nbin; i++){

    // electron channel
    double value_e = h_e->GetBinContent(i);
    if(value_e < 1e-10)continue;
    double stat_e  = h_e->GetBinError(i);

    double syse_up = fabs(h_ejesup->GetBinContent(i) - 
			  h_e->GetBinContent(i));
    double syse_dn = fabs(h_ejesdn->GetBinContent(i) - 
			  h_e->GetBinContent(i));
    double sys_e = syse_up > syse_dn? syse_up: syse_dn;

    double total_e_2 = stat_e*stat_e+ sys_e*sys_e;

    // muon channel
    double value_m = h_mu->GetBinContent(i);
    if(value_m < 1e-10)continue;
    double stat_m  = h_mu->GetBinError(i);

    double value_m_jes = h_mujes->GetBinContent(i);
    if(value_m_jes <1e-10)continue;

    double rel_sysm_up = fabs(h_mujesup->GetBinContent(i) -
			     value_m_jes)/value_m_jes;
    
    double rel_sysm_dn = fabs(h_mujesdn->GetBinContent(i) -
			     value_m_jes)/value_m_jes;
   
    double rel_sysm = rel_sysm_up > rel_sysm_dn?
      rel_sysm_up: rel_sysm_dn;

    double sys_m = value_m*rel_sysm;

    double total_m_2 = stat_m*stat_m+ sys_m*sys_m;

    cout << "Bin " << i << ": rel_syse= " << sys_e/value_e<< "\t" 
	 << "rel_sysm = " << rel_sysm<< "\t" 
  	 << "rel_state= " << stat_e/value_e << "\t rel_statmu = " <<
      stat_m/value_m << endl;
    

    double alpha = (total_m_2 - correlation*sys_e*sys_m)/
      (total_e_2 + total_m_2 - 2*correlation*sys_e*sys_m);

    cout << "Bin " << i << " alpha = " << alpha << endl;

    double combined_value = alpha * value_e + 
      (1-alpha)*value_m;

    h_combine->SetBinContent(i,combined_value);

    double combined_error = 
      alpha*alpha*(total_e_2) 
      + (1-alpha)*(1-alpha)*(total_m_2) + 
      2*alpha*(1-alpha)*correlation *sys_e*sys_m;


    if(combined_error>1e-10)combined_error = sqrt(combined_error);
    else combined_error=1e-10;

    h_combine->SetBinError(i,combined_error);

  }

  cout << "h_combine integral = " << h_combine->Integral() << endl;

  for(int i=1; i <= h_combine->GetNbinsX(); i++){

    if(h_combine->GetBinContent(i)>1e-10)
      cout << "Bin " << i << ": " << h_combine->GetBinContent(i) 
	   << " +- " << h_combine->GetBinError(i) << "\t"
	   << "rele = " << 
	(h_combine->GetBinError(i)/h_combine->GetBinContent(i)) 
	   << endl;

  }


  h_e      ->SetXTitle(xtitle.data());
  h_mu     ->SetXTitle(xtitle.data());
  h_combine->SetXTitle(xtitle.data());


  TCanvas* c1 = new TCanvas("c1","",500,500);
  if(logScale)
    c1->SetLogy(1);
  float x1NDC = 0.67;
  float y1NDC = 0.764831;
  float x2NDC = 0.830;
  float y2NDC = 0.908898;

  h_e->Draw("e");
  h_mu->Draw("esame");
  h_combine->Draw("e1same");

  TLegend* leg = new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);  
  
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->AddEntry(h_combine, "combined");
  leg->AddEntry(h_e, "e");
  leg->AddEntry(h_mu, "#mu");
  leg->Draw("same");

  std::string remword  ="h_";
  std::string remword2 ="h";
  size_t pos  = eikoName.find(remword);
  if(pos!= std::string::npos)
    eikoName.replace(pos,remword2.length(),"");

  string dirName = "fig";
  gSystem->mkdir(dirName.data());

  string fileName = dirName + "/" + "combine_emu" + eikoName;
  c1->Print(Form("%s.eps",fileName.data()));
  c1->Print(Form("%s.gif",fileName.data()));
  c1->Print(Form("%s.pdf",fileName.data()));
}
