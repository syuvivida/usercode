
void errmc(float nsig,float ntotal,float& eff, float& err)
{
  eff = nsig/ntotal;
  err = sqrt( (1-eff)*eff/ntotal);
//   cout << "efficiency = " << eff << " +- " << err << endl;
}

void errdataMC(float nsig, float nerr, float nsig1, float nerr1, float& ratio, float& err)
{

  ratio = nsig/nsig1;
  err = ratio*sqrt(nerr*nerr/nsig/nsig + nerr1*nerr1/nsig1/nsig1);
  cout << "Double ratio = " << ratio << " +- " << err << endl;
}


void countPhotonID()
{
  const double fBinsEta[]={0,1.45,1.7,2.5};
  const double fBinsPt[]={15,20,30,50,80,120};
  const int nPtBin = 5;
  const int nEtaBin = 2;

  float data[nEtaBin][nPtBin][2]={{{0.}}}; // the last index is before and after cuts
  float mc[nEtaBin][nPtBin][2]={{{0.}}};  // the last index is before and after cuts

  
  TFile* inf_central = new TFile("TightEESBData_template.root");
  TFile* inf_vary    = new TFile("vary_showershape.root");

  char* dec[2] = {"EB","EE"};

  for(int ieta=0; ieta < nEtaBin; ieta++)
    {
      for(int ipt=0; ipt < nPtBin; ipt++)
	{
	  data[ieta][ipt][0]= data[ieta][ipt][1] = 
	    mc[ieta][ipt][0]= mc[ieta][ipt][1] = 0.;
	  
	  TH1F* htemp = (TH1F*)inf_central->FindObjectAny(
	  Form("h_%s_comb3Iso_EGdata_pt%d", dec[ieta], (int)fBinsPt[ipt]));	
	  cout << "Find " << htemp->GetName() << " in " << inf_central->GetName() << endl;
	  data[ieta][ipt][0] = htemp->Integral();

	  TH1F* htemp1 = (TH1F*)inf_vary->FindObjectAny(
	  Form("h_%s_comb3Iso_EGdata_pt%d", dec[ieta], (int)fBinsPt[ipt]));	  
	  cout << "Find " << htemp1->GetName() << " in " << inf_vary->GetName() << endl;
	  data[ieta][ipt][1] = htemp1->Integral();

	  float dataRatio=0, dataRatioErr=0;
	  errmc(data[ieta][ipt][1],data[ieta][ipt][0],dataRatio,dataRatioErr);
	  
	  TH1F* htemp2 = (TH1F*)inf_central->FindObjectAny(
	   Form("h_%s_comb3Iso_sig_pt%d", dec[ieta], (int)fBinsPt[ipt]));
	  cout << "Find " << htemp2->GetName() << " in " << inf_central->GetName() << endl;

	  TH1F* htemp2_1 = (TH1F*)inf_central->FindObjectAny(
	   Form("h_%s_comb3Iso_bkg_pt%d", dec[ieta], (int)fBinsPt[ipt]));
	  cout << "Find " << htemp2_1->GetName() << " in " << inf_central->GetName() << endl;

	  mc[ieta][ipt][0] = htemp2->Integral() + htemp2_1->Integral();
	  

	  TH1F* htemp3 = (TH1F*)inf_vary->FindObjectAny(
	   Form("h_%s_comb3Iso_sig_pt%d", dec[ieta], (int)fBinsPt[ipt]));
	  cout << "Find " << htemp3->GetName() << " in " << inf_vary->GetName() << endl;

	  TH1F* htemp3_1 = (TH1F*)inf_vary->FindObjectAny(
	   Form("h_%s_comb3Iso_bkg_pt%d", dec[ieta], (int)fBinsPt[ipt]));
	  cout << "Find " << htemp3_1->GetName() << " in " << inf_vary->GetName() << endl;

	  mc[ieta][ipt][1] = htemp3->Integral() + htemp3_1->Integral();

	  
	  float mcRatio=0, mcRatioErr=0;
	  errmc(mc[ieta][ipt][1],mc[ieta][ipt][0],mcRatio,mcRatioErr);
	  
	  cout << "For " << dec[ieta] << " and pt = " << fBinsPt[ipt] << " -- " << fBinsPt[ipt+1] << endl; 
	  float doubleRatio=0, doubleRatioErr=0;
	  errdataMC(dataRatio,dataRatioErr,mcRatio,mcRatioErr,
		    doubleRatio, doubleRatioErr);
	  
	} // end of loop over pt bins
    } // end of loop over eta bins
  
  


}
