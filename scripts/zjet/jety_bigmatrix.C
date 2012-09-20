#include <TH1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TSystem.h>
#include <TLegend.h>
#include <iostream>
#include <string>

using namespace std;
//==========================================================================
//
//  Dump element content of TMatrixD and TVectorD
//
//==========================================================================
void dumpElements(TMatrixD& a)
{
  cout << endl << endl;
  const int nrows = a.GetNrows();
  const int ncols = a.GetNcols();
  if(nrows==ncols)
    cout << "determinent = " << a.Determinant() << endl;
  a.Print();
  cout << endl << endl;
  return;
}


void dumpElements(TVectorD& a)
{
  cout << endl << endl;
  a.Print();
  cout << endl << endl;
  return;
}

//==========================================================================
//
//  Main function
//
//==========================================================================

void jety_bigmatrix(
		    bool update=false,
		    bool acceptanceCorr=true,
		    double correlation=1.0, 
		    bool logScale=false 
		    )
{
  // declare histograms
  TH1D* h_e;
  TH1D* h_ejes;
  TH1D* h_ejesup;
  TH1D* h_ejesdn;
  TH1D* heff_e;

  TH1F* h_mu;
  TH1F* h_mujes;
  TH1F* h_mujesup;
  TH1F* h_mujesdn;
  TH1D* heff_mu;

  TH1D* h_combine; // for combining electron and muon channels  
  TH1D* h_corr; // for correction of acceptance


  TMatrixD covEle(15,15);
  TMatrixD covMuo(15,15);

  TMatrixD covEle_obj(15,15);
  TMatrixD covMuo_obj(15,15);

  TFile f_cove("mainCore/electron_covMatrix_Bayes.root");
  if (f_cove.IsZombie()) {
    cout << endl << "Error opening file" << f_cove.GetName() << endl;
    return;
  }
  else
    cout << endl << "Opened " << f_cove.GetName() << endl << endl;
 
  covEle =*((TMatrixD*)(f_cove.Get("covMatrix_Bayes_Yjet")));
  covEle.Print();
  
  TFile f_covm("mainCore/muon_covMatrix_jety.root");
  if (f_covm.IsZombie()) {
    cout << endl << "Error opening file" << f_covm.GetName() << endl;
    return;
  }
  else
    cout << endl << "Opened " << f_covm.GetName() << endl << endl;
 
  covMuo =*((TMatrixD*)(f_covm.Get("Bayes_YDiff")));
  covMuo.Print();


  //===========================================================
  // electron channel
  //===========================================================


  std::string eikoName="h_jety";
  std::string remword3  ="h_";
  std::string corrName = eikoName;
  size_t pos3  = corrName.find(remword3);
  if(pos3!= std::string::npos)
    corrName.replace(pos3,remword3.length(),"");


  // acceptance correction for electrons
  TFile f_crack("mainCore/ave_sherpamadgraph.root");
  if (f_crack.IsZombie()) {
    cout << endl << "Error opening file" << f_crack.GetName() << endl;
    return;
  }
  else
    cout << endl << "Opened " << f_crack.GetName() << endl << endl;

  h_corr = (TH1D*)(f_crack.Get(Form("have_%s",corrName.data())));

  // for debugging
  if( !acceptanceCorr )
    h_corr->Reset();
  for(int i=1; i<= h_corr->GetNbinsX(); i++)
    {
      if( !acceptanceCorr )
	{
	  h_corr->SetBinContent(i,1.0);
	  h_corr->SetBinError(i,1e-6);
	}
      cout << "Correction for bin " << i << " = " 
	   << h_corr->GetBinContent(i) << " +- " << h_corr->GetBinError(i) 
	   << endl;
    }


  //================================================================
  // efficiency errors
  //================================================================
 
  std::string xtitle;
  std::string muName;
  std::string kengName;
  std::string effElectronName;
  std::string effMuonName;

  
  effElectronName = "EfficienyVsYjet";
  effMuonName = "EffCorr_Yjet";


  TFile f_eff_electron("mainCore/EffiZCut_Jer0.root");
  if(f_eff_electron.IsZombie()){
    cout << endl << "Error opening file" << f_eff_electron.GetName() << endl << endl;
    return;
  }
  else
    cout << endl << "Opened " << f_eff_electron.GetName() << endl << endl;

  heff_e = (TH1D*)(f_eff_electron.Get(effElectronName.data()));
  heff_e -> SetName(Form("heff_e_%s",corrName.data()));


  TFile f_eff_muon("mainCore/EffCorr_root_091912.root");
  if(f_eff_muon.IsZombie()){
    cout << endl << "Error opening file" << f_eff_muon.GetName() << endl << endl;
    return;
  }
  else
    cout << endl << "Opened " << f_eff_muon.GetName() << endl << endl;

  heff_mu = (TH1D*)(f_eff_muon.Get(effMuonName.data()));
  heff_mu -> SetName(Form("heff_mu_%s",corrName.data()));


  xtitle = "|Y(jet)|";


  //===========================================================
  // electron channel
  //===========================================================
  TFile f_e("mainCore/cts_CorrectedPlotsZCut.root");
  if (f_e.IsZombie()) {
    cout << endl << "Error opening file" << f_e.GetName() << endl << endl;
    return;
  }
  else
    cout << endl << "Opened " << f_e.GetName() << endl << endl;

  h_e  = (TH1D*)(f_e.Get(eikoName.data()));
  h_e  -> SetName(Form("h_e_%s",corrName.data()));
  h_e  -> SetTitle("");
  h_e  -> Sumw2();

  h_e  -> Reset();

  double darko[12][2]={
    {13218.9 , 155.4},
    {13038.2 , 159.7},
    {12338.0 , 155.2},
    {12004.6 , 153.5},
    {11173.1 , 150.4},
    { 9959.5 , 142.6},
    { 9266.5 , 137.1},
    { 7725.4 , 126.6},
    { 6897.7 , 120.1},
    { 5876.9 , 110.5},
    { 4814.4 , 101.1},
    { 3277.8 , 80.8}
  };

  for(int i=1;i<=12;i++){
    h_e->SetBinContent(i, darko[i-1][0]);
    h_e->SetBinError(i, darko[i-1][1]);
  }
  //-------------------------------------------------------------------
  // use electron channel to count the number of non-zero bins
  // so to avoid singular error matrix
  int nbins = 0; 
  const double threshold = 1.0;
  for(int ie=1; ie<= h_e->GetNbinsX(); ie++)
    if(h_e->GetBinContent(ie)> threshold)nbins++;
  cout << "There are " << nbins << " bins with non-zero content" << endl;  
  //-------------------------------------------------------------------

  h_e  -> Divide(h_corr);
  double scale_electron = 1.0/h_e->Integral();

  covEle_obj = covEle * pow(scale_electron,2);
  covEle_obj.Print();

  h_e  -> Scale(scale_electron);
  h_e  -> SetYTitle("Arbitrary Unit");
  h_e  -> SetTitleOffset(2.0,"Y");
  h_e  -> GetYaxis()->SetDecimals();
  h_e  -> GetXaxis()->SetDecimals();
  h_e  -> SetLineColor(kBlue-7);
  h_e  -> SetMarkerColor(kBlue-7);
  h_e  -> SetMarkerSize(1);
  h_e  -> SetMarkerStyle(24);

  cout << "h_e integral = " << h_e->Integral() << endl;

  h_ejes= (TH1D*)(f_e.Get(eikoName.data()));
  h_ejes    -> SetName("h_ejes");
  h_ejes    -> Sumw2();
  h_ejes    -> Scale(1.0/h_ejes->Integral());
  cout << "h_ejes integral = " << h_ejes->Integral() << endl;

  TFile f_ejesup("mainCore/cts_CorrectedPlotsZCut_JesUp_NoBkgSub.root");
  if (f_ejesup.IsZombie()) {
    cout << endl << "Error opening file" << f_ejesup.GetName() << endl << endl;
    return;
  }
  else
    cout << endl << "Opened " << f_ejesup.GetName() << endl << endl;

  h_ejesup  = (TH1D*)(f_ejesup.Get(eikoName.data()));
  h_ejesup  -> SetName("h_ejesup");
  h_ejesup  -> SetTitle("");
  h_ejesup  -> Sumw2();
  h_ejesup  -> Scale(1.0/h_ejesup->Integral());


  cout << "h_ejesup integral = " << h_ejesup->Integral() << endl;

  TFile f_ejesdn("mainCore/cts_CorrectedPlotsZCut_JesDn_NoBkgSub.root");
  if (f_ejesdn.IsZombie()) {
    cout << endl << "Error opening file" << f_ejesdn.GetName() << endl << endl;
    return;
  }
  else
    cout << endl << "Opened " << f_ejesdn.GetName() << endl << endl;

  h_ejesdn  = (TH1D*)(f_ejesdn.Get(eikoName.data()));
  h_ejesdn  -> SetName("h_ejesdn");
  h_ejesdn  -> SetTitle("");
  h_ejesdn  -> Sumw2();
  h_ejesdn  -> Scale(1.0/h_ejesdn->Integral());

  cout << "h_ejesdn integral = " << h_ejesdn->Integral() << endl;

  h_combine= (TH1D*)h_e->Clone(Form("h_combine_%s",corrName.data()));
  h_combine  -> Reset();
  h_combine  -> SetTitle("");
  h_combine  -> SetLineColor(1);
  h_combine  -> SetMarkerColor(1);
  h_combine  -> SetMarkerSize(1);
  h_combine  -> SetMarkerStyle(8);


  //===========================================================
  // muon channel
  //===========================================================
  muName = "Bayes_UnfoldedYjet";
  kengName = "Z1jets_1jeta_BE";

  TFile f_mu("mainCore/Unfolded_Yjet_091912.root");
  if (f_mu.IsZombie()) {
    cout << endl << "Error opening file" << f_mu.GetName() << endl << endl;
    return;
  }
  else
    cout << endl << "Opened " << f_mu.GetName() << endl << endl;

  h_mu = (TH1F*)(f_mu.Get(muName.data()));
  h_mu -> Sumw2();
  h_mu -> SetName(Form("h_mu_%s",corrName.data()));
  h_mu  -> SetTitle("");

  double scale_muon = 1.0/h_mu->Integral();

  covMuo_obj = covMuo * pow(scale_muon,2);
  covMuo_obj.Print();

  h_mu -> Scale(scale_muon);
  h_mu -> SetLineColor(kRed-7);
  h_mu -> SetMarkerColor(kRed-7);
  h_mu -> SetMarkerSize(1);
  h_mu -> SetMarkerStyle(21);

  cout << "h_mu integral = " << h_mu->Integral() << endl;

  TFile f_jetsys_mu("mainCore/DoubleMu2011_JESuncertainty_JetY_061712.root");
  if (f_jetsys_mu.IsZombie()) {
    cout << endl << "Error opening file" << f_jetsys_mu.GetName() << endl << endl;
    return;
  }
  else
    cout << endl << "Opened " << f_jetsys_mu.GetName() << endl << endl;

  h_mujes = (TH1F*)(f_jetsys_mu.Get(kengName.data()));
  h_mujes -> Sumw2();
  h_mujes -> Scale(1.0/h_mujes->Integral());

  cout << "h_mujes integral = " << h_mujes->Integral() << endl;

  h_mujesup = (TH1F*)(f_jetsys_mu.Get(Form("%sUp",kengName.data())));
  h_mujesup -> Sumw2();
  h_mujesup -> Scale(1.0/h_mujesup->Integral());

  cout << "h_mujesup integral = " << h_mujesup->Integral() << endl;

  h_mujesdn = (TH1F*)(f_jetsys_mu.Get(Form("%sDn",kengName.data())));
  h_mujesdn -> Sumw2();
  h_mujesdn -> Scale(1.0/h_mujesdn->Integral());

  cout << "h_mujesdn integral = " << h_mujesdn->Integral() << endl;


  // declare the big matrix
  const int NELE=2*nbins;
  TMatrixD errorM(NELE,NELE);

  TMatrixD U(NELE,nbins);
  for(int irow=0; irow < NELE; irow++)
    for(int icol=0; icol < nbins; icol++)
      U(irow,icol) = ( (irow== icol) || (irow== icol+nbins) )? 1:0;
  // debug
  dumpElements(U);

  TMatrixD transposeU(nbins,NELE);
  transposeU.Transpose(U);
  dumpElements(transposeU);

  
  TVectorD measurement(NELE);

  // now loop over bins to set the matrix content
  for(int ibin=0; ibin<nbins; ibin++){

    // electron channel
    double value_e = h_e->GetBinContent(ibin+1);
    if(value_e < 1e-10)continue;

    double stat_mc_eff_e = heff_e->GetBinContent(ibin+1)<1e-6? 0: 
      value_e*heff_e->GetBinError(ibin+1)/heff_e->GetBinContent(ibin+1);

    
    double stat_e  = sqrt( covEle_obj(ibin,ibin) + 
			   pow(stat_mc_eff_e,2));

    double value_e_jes = h_ejes->GetBinContent(ibin+1);
    if(value_e_jes <1e-10)continue;

    double rel_syse_up = fabs(h_ejesup->GetBinContent(ibin+1) -
			      value_e_jes)/value_e_jes;
    
    double rel_syse_dn = fabs(h_ejesdn->GetBinContent(ibin+1) -
			      value_e_jes)/value_e_jes;
   
    double rel_syse = rel_syse_up > rel_syse_dn?
      rel_syse_up: rel_syse_dn;

    double sys_e = value_e*rel_syse;

    double total_e_2 = stat_e*stat_e+ sys_e*sys_e;


    // muon channel
    double value_m = h_mu->GetBinContent(ibin+1);
    if(value_m < 1e-10)continue;

    double stat_mc_eff_m = heff_mu->GetBinContent(ibin+1)<1e-6? 0: 
      value_m*heff_mu->GetBinError(ibin+1)/heff_mu->GetBinContent(ibin+1);

    double stat_m  = sqrt( covMuo_obj(ibin,ibin) + 
			   pow(stat_mc_eff_m,2));
    
    double value_m_jes = h_mujes->GetBinContent(ibin+1);
    if(value_m_jes <1e-10)continue;

    double rel_sysm_up = fabs(h_mujesup->GetBinContent(ibin+1) -
			      value_m_jes)/value_m_jes;
    
    double rel_sysm_dn = fabs(h_mujesdn->GetBinContent(ibin+1) -
			      value_m_jes)/value_m_jes;
   
    double rel_sysm = rel_sysm_up > rel_sysm_dn?
      rel_sysm_up: rel_sysm_dn;

    double sys_m = value_m*rel_sysm;

    double total_m_2 = stat_m*stat_m+ sys_m*sys_m;

    measurement(ibin) = value_e;
    measurement(nbins+ibin) = value_m;

    // first put electron error matrix component
    errorM(ibin,ibin) = total_e_2; 

    // then muon
    errorM(nbins+ibin,nbins+ibin) = total_m_2; 

    errorM(ibin,nbins+ibin) = errorM(nbins+ibin,ibin) = 
      correlation*sys_e*sys_m;


    // for now 
    // the correlation between different bins within the same lepton 
    // channel is zero

    for(int jbin=0; jbin<nbins; jbin++)
      {
	if(jbin!=ibin){
	  errorM(ibin,jbin)=covEle_obj(ibin,jbin);
	  errorM(nbins+ibin,nbins+jbin)=covMuo_obj(ibin,jbin);
	}	
      }

  } // loop over the number of bins
  

  // debug
  // print out measurement value  
  dumpElements(measurement);

  //   // print out error matrix component
  dumpElements(errorM);


  TMatrixD errorInverse = errorM;
  double* det;
  errorInverse.Invert(det);
  dumpElements(errorInverse);

  //   TMatrixD Unit = errorInverse*errorM;
  // dumpElements(Unit);
  

  TMatrixD matrixRight(nbins,NELE);  
  matrixRight = transposeU*errorInverse;
  dumpElements(matrixRight);

  TMatrixD matrixLeft(nbins,nbins);
  matrixLeft = transposeU*(errorInverse*U);
  dumpElements(matrixLeft);

  TMatrixD matrixLeftInverse = matrixLeft;
  double* det2;
  matrixLeftInverse.Invert(det2);
  dumpElements(matrixLeftInverse);

  TMatrixD lambda(nbins,NELE);
  lambda= matrixLeftInverse*matrixRight;
  dumpElements(lambda);
  
  TMatrixD transposeLambda(NELE,nbins);
  transposeLambda.Transpose(lambda);
  dumpElements(transposeLambda);

  dumpElements(lambda);

  TVectorD combined_value(nbins);
  combined_value = lambda*measurement;
  dumpElements(combined_value);

  TMatrixD combined_error(nbins,nbins);
  combined_error = lambda*(errorM*transposeLambda);
  

  // after all matrix operation, now set the histogram
  for(int i=0; i<nbins;i++){

    h_combine->SetBinContent(i+1, combined_value(i));
  
    double error = combined_error(i,i);
  
    if(error>1e-10)error = sqrt(error);
    else error=1e-10;

    h_combine->SetBinError(i+1,error);

  }
  
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

  std::string command = "recreate";
  if(update)command="update";
  TFile* outFile = new TFile("jety_bigmatrix.root", command.data());
  h_e      ->Write();
  h_mu     ->Write();
  h_combine->Write();
  heff_e   ->Write();
  heff_mu  ->Write();
  outFile->Close();


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

  string fileName = dirName + "/" + "final_emu" + eikoName;
  c1->Print(Form("%s.eps",fileName.data()));
  c1->Print(Form("%s.gif",fileName.data()));
  c1->Print(Form("%s.pdf",fileName.data()));
}
