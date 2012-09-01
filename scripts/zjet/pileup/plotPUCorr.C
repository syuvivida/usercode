#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"

void plotPUCorr(std::string filename)
{
  TFile* fin = TFile::Open(filename.data());
  
  const double ymin = 0.95;
  const double ymax = 1.05;

  TCanvas* c1 = new TCanvas("c1","",1200,1000);
  c1->Divide(2,2);
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
  r_zy_corr->Draw("hist");
  a->Draw("same");

  c1->cd(2);
  TLine* b = new TLine(r_jety_corr->GetBinLowEdge(1),1.0,
		       r_jety_corr->GetBinLowEdge(nbins+1),1.0);
  b->SetLineWidth(2);
  b->SetLineStyle(2);
  b->SetLineColor(4);
  r_jety_corr->SetMinimum(ymin);
  r_jety_corr->SetMaximum(ymax);
  r_jety_corr->Draw("hist");
  b->Draw("Same");

  c1->cd(3);
  TLine* c = new TLine(r_ystar_corr->GetBinLowEdge(1),1.0,
		       r_ystar_corr->GetBinLowEdge(nbins+1),1.0);
  c->SetLineWidth(2);
  c->SetLineStyle(2);
  c->SetLineColor(4);
  r_ystar_corr->SetMinimum(ymin);
  r_ystar_corr->SetMaximum(ymax);
  r_ystar_corr->Draw("hist");
  c->Draw("same");

  c1->cd(4);
  TLine* d = new TLine(r_yB_corr->GetBinLowEdge(1),1.0,
		       r_yB_corr->GetBinLowEdge(nbins+1),1.0);
  d->SetLineWidth(2);
  d->SetLineStyle(2);
  d->SetLineColor(4);
  r_yB_corr->SetMinimum(ymin);
  r_yB_corr->SetMaximum(ymax);
  r_yB_corr->Draw("hist");
  d->Draw("same");
  
  c1->Print("pileup/PUcorr.eps");
  c1->Print("pileup/PUcorr.gif");
  c1->Print("pileup/PUcorr.pdf");

}
