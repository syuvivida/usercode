
#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"

void plotJESCorr(std::string filename)
{

  TF1* func1 = new TF1("func1","[0]+[1]*x");
  TFile* fin = TFile::Open(filename.data());
  
  const double ymin = 0.9;
  const double ymax = 1.1;
  const double offset=2.00;

  std::string ytitle = "JES correction";

  TCanvas* c1 = new TCanvas("c1","",1200,1000);
  c1->Divide(2,2);
  /////////////////////////////////////////////////////////////////
  c1->cd(1);
  int nbins = r_zy_corr->GetNbinsX();
  cout << "nbins = " << nbins << endl;
  TLine* a = new TLine(r_zy_corr->GetBinLowEdge(1),1.0,
		       r_zy_corr->GetBinLowEdge(nbins+1),1.0);
  a->SetLineWidth(2);
  a->SetLineStyle(2);
  a->SetLineColor(4);
  r_zy_corr->SetMinimum(ymin);
  r_zy_corr->SetMaximum(ymax);
  r_zy_corr->SetYTitle(ytitle.data());
  r_zy_corr->SetTitleOffset(offset,"Y");
  r_zy_corr->Draw("e1");
//   r_zy_corr->Fit("func1");
//   cout << "r_zy_corr correction range: " << func1->GetParameter(1)*200.0 << endl;
  a->Draw("same");

  /////////////////////////////////////////////////////////////////
  c1->cd(2);
  TLine* b = new TLine(r_jety_corr->GetBinLowEdge(1),1.0,
		       r_jety_corr->GetBinLowEdge(nbins+1),1.0);
  b->SetLineWidth(2);
  b->SetLineStyle(2);
  b->SetLineColor(4);
  r_jety_corr->SetMinimum(ymin);
  r_jety_corr->SetMaximum(ymax);
  r_jety_corr->SetYTitle(ytitle.data());
  r_jety_corr->SetTitleOffset(offset,"Y");
  r_jety_corr->Draw("e1");
//   r_jety_corr->Fit("func1");
//   cout << "r_jety_corr correction range: " << func1->GetParameter(1)*200.0 << endl;
  b->Draw("Same");


  /////////////////////////////////////////////////////////////////
  c1->cd(3);
  TLine* c = new TLine(r_ystar_corr->GetBinLowEdge(1),1.0,
		       r_ystar_corr->GetBinLowEdge(nbins+1),1.0);
  c->SetLineWidth(2);
  c->SetLineStyle(2);
  c->SetLineColor(4);
  r_ystar_corr->SetMinimum(ymin);
  r_ystar_corr->SetMaximum(ymax);
  r_ystar_corr->SetYTitle(ytitle.data());
  r_ystar_corr->SetTitleOffset(offset,"Y");
  r_ystar_corr->Draw("e1");
//   r_ystar_corr->Fit("func1");
//   cout << "r_ystar_corr correction range: " << func1->GetParameter(1)*200.0 << endl;
  c->Draw("same");

  /////////////////////////////////////////////////////////////////
  c1->cd(4);
  TLine* d = new TLine(r_yB_corr->GetBinLowEdge(1),1.0,
		       r_yB_corr->GetBinLowEdge(nbins+1),1.0);
  d->SetLineWidth(2);
  d->SetLineStyle(2);
  d->SetLineColor(4);
  r_yB_corr->SetMinimum(ymin);
  r_yB_corr->SetMaximum(ymax);
  r_yB_corr->SetYTitle(ytitle.data());
  r_yB_corr->SetTitleOffset(offset,"Y");
  r_yB_corr->Draw("e1");
//   r_yB_corr->Fit("func1");
//   cout << "r_yB_corr correction range: " << func1->GetParameter(1)*200.0 << endl; 
  d->Draw("same");
  
  gSystem->mkdir("kengJES");
  c1->Print("kengJES/JEScorr.eps");
  c1->Print("kengJES/JEScorr.gif");
  c1->Print("kengJES/JEScorr.pdf");

}
