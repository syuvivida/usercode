#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"

void plotPUSys(std::string filename)
{
  TFile* fin = TFile::Open(filename.data());
  
  const double ymin = 0.95;
  const double ymax = 1.05;
  const double offset = 2.0;
  std::string ytitle = "PU Systematics/Central";

  TCanvas* c1 = new TCanvas("c1","",1200,1000);
  c1->Divide(2,2);
  c1->cd(1);
  int nbins = r_zy_corr->GetNbinsX();
  cout << "nbins = " << nbins << endl;
  r_zy_up->SetMinimum(ymin);
  r_zy_up->SetMaximum(ymax);
  r_zy_up->SetLineColor(2);
  r_zy_up->SetMarkerColor(2);
  r_zy_up->SetYTitle(ytitle.data());
  r_zy_up->SetTitleOffset(offset,"Y");
  r_zy_up->Draw("hist");

  r_zy_down->SetMinimum(ymin);
  r_zy_down->SetMaximum(ymax);
  r_zy_down->SetLineColor(4);
  r_zy_down->SetMarkerColor(4);
  r_zy_down->Draw("histsame");
  TLine* a = new TLine(r_zy_corr->GetBinLowEdge(1),1.0,
		       r_zy_corr->GetBinLowEdge(nbins+1),1.0);
  a->SetLineWidth(2);
  a->SetLineStyle(2);
  a->SetLineColor(1);
  a->Draw("same");
  float x1NDC = 0.551;
  float y1NDC = 0.713;
  float x2NDC = 0.689;
  float y2NDC = 0.905;

  TLegend* leg = new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  leg->AddEntry(r_zy_up, "#sigma_{minBias}=71.4 mb");
  leg->AddEntry(r_zy_down,"#sigma_{minBias}=64.6 mb");
  leg->Draw("same");



  /////////////////////////////////////////////////////////////////
  c1->cd(2);
  r_jety_up->SetMinimum(ymin);
  r_jety_up->SetMaximum(ymax);
  r_jety_up->SetLineColor(2);
  r_jety_up->SetMarkerColor(2);
  r_jety_up->SetYTitle(ytitle.data());
  r_jety_up->SetTitleOffset(offset,"Y");
  r_jety_up->Draw("hist");

  r_jety_down->SetMinimum(ymin);
  r_jety_down->SetMaximum(ymax);
  r_jety_down->SetLineColor(4);
  r_jety_down->SetMarkerColor(4);
  r_jety_down->Draw("histsame");

  TLine* b = new TLine(r_jety_corr->GetBinLowEdge(1),1.0,
		       r_jety_corr->GetBinLowEdge(nbins+1),1.0);
  b->SetLineWidth(2);
  b->SetLineStyle(2);
  b->SetLineColor(1);

  b->Draw("Same");
  leg->Draw("same");
  /////////////////////////////////////////////////////////////////
  c1->cd(3);
  r_ystar_up->SetMinimum(ymin);
  r_ystar_up->SetMaximum(ymax);
  r_ystar_up->SetLineColor(2);
  r_ystar_up->SetMarkerColor(2);
  r_ystar_up->SetYTitle(ytitle.data());
  r_ystar_up->SetTitleOffset(offset,"Y");
  r_ystar_up->Draw("hist");

  r_ystar_down->SetMinimum(ymin);
  r_ystar_down->SetMaximum(ymax);
  r_ystar_down->SetLineColor(4);
  r_ystar_down->SetMarkerColor(4);
  r_ystar_down->Draw("histsame");

  TLine* c = new TLine(r_ystar_corr->GetBinLowEdge(1),1.0,
		       r_ystar_corr->GetBinLowEdge(nbins+1),1.0);
  c->SetLineWidth(2);
  c->SetLineStyle(2);
  c->SetLineColor(1);

  c->Draw("same");
  leg->Draw("same");

  /////////////////////////////////////////////////////////////////
  c1->cd(4);
  r_yB_up->SetMinimum(ymin);
  r_yB_up->SetMaximum(ymax);
  r_yB_up->SetLineColor(2);
  r_yB_up->SetMarkerColor(2);
  r_yB_up->SetYTitle(ytitle.data());
  r_yB_up->SetTitleOffset(offset,"Y");
  r_yB_up->Draw("hist");

  r_yB_down->SetMinimum(ymin);
  r_yB_down->SetMaximum(ymax);
  r_yB_down->SetLineColor(4);
  r_yB_down->SetMarkerColor(4);
  r_yB_down->Draw("histsame");

  TLine* d = new TLine(r_yB_corr->GetBinLowEdge(1),1.0,
		       r_yB_corr->GetBinLowEdge(nbins+1),1.0);
  d->SetLineWidth(2);
  d->SetLineStyle(2);
  d->SetLineColor(1);

  d->Draw("same");
  leg->Draw("same");
  c1->Print("pileup/PUsys.eps");
  c1->Print("pileup/PUsys.gif");
  c1->Print("pileup/PUsys.pdf");

}
