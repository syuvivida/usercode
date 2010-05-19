#include <TH1.h>
#include <string>
#include <iostream>

using namespace std;

const int MAXBIN_ALLOWED=1000;

void ratioErr(Double_t n1, Double_t n1err, Double_t n2, Double_t n2err,
	      Double_t& ratio, Double_t& err)
{
  if(fabs(n1+n2)<1e-6){ratio=0;err=0;return;}
  ratio = (n1)/(n1+n2);
  
  err= pow(1/(n1+n2) - n1/(n1+n2)/(n1+n2),2)*n1err*n1err+
    pow(n1/(n1+n2)/(n1+n2),2)*n2err*n2err; 

  err = sqrt(err);

}


void sigFracErr(Double_t rd, Double_t rderr, 
		Double_t rs, Double_t rserr,
		Double_t rb, Double_t rberr,
		Double_t& sigfrac, Double_t& err)
{

  if(fabs(rs-rb)<1e-6){sigfrac=0; err=0; return;}
  sigfrac = (rd-rb)/(rs-rb);

  err = pow(1/(rs-rb),2)*rderr*rderr + 
    pow((rd-rb)/(rs-rb)/(rs-rb),2)*rserr*rserr + 
    pow(-1.0/(rs-rb)+(rd-rb)/(rs-rb)/(rs-rb),2)*rberr*rberr;

  err= sqrt(err);

}



void sigFracErr_inOneBin(int bin, Double_t rd, Double_t rderr, 
			 Double_t rs, Double_t rserr,
			 Double_t rb, Double_t rberr,
			 Double_t& sigfrac, Double_t& err)
{

  if(fabs(rs-rb)<1e-6){sigfrac=0; err=0; return;}
  if(bin==1)
    {
      if(fabs(rd)<1e-6){sigfrac=0; err=0; return;}
      sigfrac = (rd-rb)*rs/(rs-rb)/rd;
    }
  else
    {
      if(fabs(1-rd)<1e-6){sigfrac=0; err=0; return;}
      sigfrac = (rd-rb)*(1-rs)/(rs-rb)/(1-rd);
    }


  if(bin==1)
    err = pow(rs/(rs-rb)/rd -(rd-rb)*rs/(rs-rb)/rd/rd,2)*rderr*rderr + 
      pow((rd-rb)/(rs-rb)/rd - (rd-rb)*rs/(rs-rb)/(rs-rb)/rd,2)*rserr*rserr + 
      pow(-rs/(rs-rb)/rd+(rd-rb)*rs/(rs-rb)/(rs-rb)/rd,2)*rberr*rberr;

  else
    err = pow((1-rs)/(rs-rb)/(1-rd) +(rd-rb)*(1-rs)/(rs-rb)/(1-rd)/(1-rd),2)*rderr*rderr + 
      pow(-(rd-rb)/(rs-rb)/(1-rd) - (rd-rb)*(1-rs)/(rs-rb)/(rs-rb)/(1-rd),2)*rserr*rserr + 
      pow(-(1-rs)/(rs-rb)/(1-rd)+(rd-rb)*(1-rs)/(rs-rb)/(rs-rb)/(1-rd),2)*rberr*rberr;

  err= sqrt(err);

}



void histoError(int MAX, TH1D* h, Double_t& bin1, Double_t& binerr1, 
		Double_t& bin2, Double_t&binerr2)
{

  bin1 = binerr1 = bin2 = binerr2 = 0;
  for(int i=0; i<= h->GetNbinsX()+1; i++)
    {
      if(i < MAX+1)
	{
	  bin1    += h->GetBinContent(i);
	  binerr1 += pow(h->GetBinError(i),2);
	}
      else
	{
	  bin2     += h->GetBinContent(i);
	  binerr2  += pow(h->GetBinError(i),2);

	}


    }

  binerr1   = sqrt(binerr1);
  binerr2   = sqrt(binerr2);
  return;

}

void purity_twobin(TH1D* dataInput, TH1D* sigTemplate, TH1D* bkgTemplate, 
		   std::string outputName, Double_t& purity_central, Double_t& purity_central_err)
{
  
  purity_central = 0;
  purity_central_err = 0;

  TH1D* signal_pos;
  TH1D* background_pos;
  TH1D* data;

  // 0 is the signal fraction from the total sample,
  // 1 is the signal fraction in the first bin,
  // 2 is the signal fraction in the second bin
  TH1D* result = new TH1D("result","",3,-0.5,2.5);

  Double_t scale=1.;

  data = (TH1D*)dataInput->Clone();
  data->SetName("data");
  data->SetLineColor(1);
  data->SetMarkerColor(1);
  data->SetXTitle("Fisher's isolation [GeV]");
  data->Sumw2();
  scale = 1.0/(Double_t)data->Integral(0,MAXBIN_ALLOWED); 
  cout << "scale for data = " << scale << endl;
  data->Scale(scale);


  signal_pos = (TH1D*)sigTemplate->Clone();
  signal_pos->SetName("signal_pos");
  signal_pos->SetLineColor(2);
  signal_pos->SetMarkerColor(2);
  signal_pos->SetFillColor(2);
  signal_pos->SetXTitle("Fisher's isolation [GeV]");
  signal_pos->Sumw2();
  scale = 1.0/(Double_t)signal_pos->Integral(0,MAXBIN_ALLOWED); 
  cout << "scale for signal template = " << scale << endl;
  signal_pos->Scale(scale);


  background_pos = (TH1D*)bkgTemplate->Clone();
  background_pos->SetName("background_pos");
  background_pos->SetLineColor(4);
  background_pos->SetMarkerColor(4);
  background_pos->SetFillColor(4);
  background_pos->SetXTitle("Fisher's isolation [GeV]");
  background_pos->Sumw2();
  scale = 1.0/(Double_t)background_pos->Integral(0,MAXBIN_ALLOWED); 
  cout << "scale for background template = " << scale << endl;
  background_pos->Scale(scale);


  // now determine where to cut

  TH1D* integral_signal = (TH1D*)sigTemplate->Clone();
  integral_signal->SetName("integral_signal");
  integral_signal->Reset();

  TH1D* integral_background = (TH1D*)bkgTemplate->Clone();
  integral_background->SetName("integral_background");
  integral_background->Reset();

  TH1D* integral_diff = (TH1D*)bkgTemplate->Clone();
  integral_diff->SetName("integral_diff");
  integral_diff->Reset();

  for(int i=1; i<= integral_signal->GetNbinsX(); i++){
     
    Double_t sigEff = signal_pos->Integral(i,MAXBIN_ALLOWED);
    Double_t bkgEff = background_pos->Integral(i,MAXBIN_ALLOWED);
    integral_signal->SetBinContent(i,sigEff);
    integral_background->SetBinContent(i,bkgEff);
    
    Double_t diff = fabs(sigEff - bkgEff);
    
    integral_diff->SetBinContent(i,diff);

  }
  

  int maxbin = integral_diff->GetMaximumBin();
  cout << "Maximum efficiency difference is at " << 
    integral_diff->GetBinLowEdge(maxbin) 
       << endl;

  cout << "Maximum efficiency difference Bin is at " << maxbin
     << endl;

  Double_t ndata1, nsig1, nbkg1; // number of events in the first bin
  Double_t ndata2, nsig2, nbkg2; // number of events in the second bin

  Double_t ndata1err, nsig1err, nbkg1err; // number of events in the first bin
  Double_t ndata2err, nsig2err, nbkg2err; // number of events in the second bin

  ndata1 = nsig1 = nbkg1 = ndata2 = nsig2 = nbkg2 = 0;
  ndata1err = nsig1err = nbkg1err = ndata2err = nsig2err = nbkg2err = 0;

  ndata1 = dataInput->Integral(0, maxbin);
  nsig1  = sigTemplate->Integral(0, maxbin);
  nbkg1  = bkgTemplate->Integral(0, maxbin);

  ndata2 = dataInput->Integral(maxbin+1, MAXBIN_ALLOWED);
  nsig2  = sigTemplate->Integral(maxbin+1, MAXBIN_ALLOWED);
  nbkg2  = bkgTemplate->Integral(maxbin+1, MAXBIN_ALLOWED);

  cout << ndata1 << "\t" << nsig1 << "\t" << nbkg1 << endl;
  cout << ndata2 << "\t" << nsig2 << "\t" << nbkg2 << endl;

  ndata1 = nsig1 = nbkg1 = ndata2 = nsig2 = nbkg2 = 0;
  ndata1err = nsig1err = nbkg1err = ndata2err = nsig2err = nbkg2err = 0;

  histoError(maxbin,dataInput,ndata1,ndata1err, ndata2,ndata2err);
  histoError(maxbin,sigTemplate,nsig1,nsig1err, nsig2,nsig2err);
  histoError(maxbin,bkgTemplate,nbkg1,nbkg1err, nbkg2,nbkg2err);


  cout << ndata1 << "\t" << nsig1 << "\t" << nbkg1 << endl;
  cout << ndata2 << "\t" << nsig2 << "\t" << nbkg2 << endl;

  cout << ndata1err << "\t" << nsig1err << "\t" << nbkg1err << endl;
  cout << ndata2err << "\t" << nsig2err << "\t" << nbkg2err << endl;

  
  Double_t ratdat;
  Double_t ratsig;
  Double_t ratbkg;
  
  
  Double_t signal_fraction = -1;
  Double_t errfrac = -1;

  Double_t errratdat;
  Double_t errratsig;
  Double_t errratbkg;  

  ratioErr(ndata1,ndata1err,ndata2,ndata2err,ratdat,errratdat);
  ratioErr(nsig1,nsig1err,nsig2,nsig2err,ratsig,errratsig);
  ratioErr(nbkg1,nbkg1err,nbkg2,nbkg2err,ratbkg,errratbkg);

  sigFracErr(ratdat, errratdat,
	     ratsig, errratsig,
	     ratbkg, errratbkg,
	     signal_fraction, errfrac);

  result->SetBinContent(1, signal_fraction);
  result->SetBinError(1, errfrac);

  cout << "Total signal fraction = " << signal_fraction << " +- " << 
    errfrac << endl;

  purity_central = signal_fraction;
  purity_central_err = errfrac;

  
  // signal_fraction within the first bin
  Double_t signal_fraction_iso1 = -1;
  Double_t errfrac_iso1 = -1;
  sigFracErr_inOneBin(
		      1, ratdat, errratdat,
		      ratsig, errratsig,
		      ratbkg, errratbkg,
		      signal_fraction_iso1, errfrac_iso1);

  Double_t signal_fraction_iso2 = -1;
  Double_t errfrac_iso2 = -1;

  sigFracErr_inOneBin(
		      2, ratdat, errratdat,
		      ratsig, errratsig,
		      ratbkg, errratbkg,
		      signal_fraction_iso2, errfrac_iso2);


  cout << "Signal fraction in the first bin = " << signal_fraction_iso1 << 
    " +- " << 
    errfrac_iso1 << endl;

  cout << "Signal fraction in the second bin = " << signal_fraction_iso2 << 
    " +- " << 
    errfrac_iso2 << endl;


  result->SetBinContent(2, signal_fraction_iso1);
  result->SetBinError(2, errfrac_iso1   );

  result->SetBinContent(3, signal_fraction_iso2);
  result->SetBinError(3, errfrac_iso2   );


  // dump histogram to a root file
  std::string histoFile = outputName + "_2bin.root";

  TFile* outFile = new TFile(histoFile.data(),"recreate");
  result->Write();
  dataInput->Write();
  data->Write();
  sigTemplate->Write();
  signal_pos->Write();
  bkgTemplate->Write();
  background_pos->Write();

  integral_signal->Write();
  integral_background->Write();
  integral_diff->Write();
  outFile->Close();



}
