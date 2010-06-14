#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"

void plotNPho()
{
  setTDRStyle();

  const int nbin=3;
  TFile *f1 = TFile::Open("npho_proj_comb.root");
  TH1F* hnpho_EB = (TH1F*)f1->FindObjectAny("h_EB_npho_EGdata_SIG");
  TH1F* hnpho_EE = (TH1F*)f1->FindObjectAny("h_EE_npho_EGdata_SIG");

  float npho_EB[10]={0};
  float nphoerr_EB[10]={0};
  float npho_EE[10]={0};
  float nphoerr_EE[10]={0};

  float x[10]={0,1,2};
  float dx[10]={0.5,0.5,0.5};

  for(int i=1;i<=nbin; i++)
    {
      npho_EB[i-1] = hnpho_EB->GetBinContent(i);
      nphoerr_EB[i-1] = hnpho_EB->GetBinError(i);
      cout << npho_EB[i-1] << "\t" << nphoerr_EB[i-1] << endl;
    }
  for(int i=1;i<=nbin; i++)
    {
      npho_EE[i-1] = hnpho_EE->GetBinContent(i);
      nphoerr_EE[i-1] = hnpho_EE->GetBinError(i);
      cout << npho_EE[i-1] << "\t" << nphoerr_EE[i-1] << endl;
    }

   
  TCanvas* c1 = new TCanvas("c1","",500,500);
  TH2F* h2 = new TH2F("h2","",3,-0.5,2.5,100,1,1000000);
  h2->SetXTitle("Number of photon candidates per event");
  h2->SetYTitle("Entries");
  h2->GetXaxis()->SetNdivisions(5);
  c1->SetLogy(1);
  h2->Draw();
//   hnpho_EB->SetXTitle("Number of photon candidates per event");
//   hnpho_EB->SetYTitle("Entries");
//   hnpho_EB->SetTitle("");
//   hnpho_EB->GetXaxis()->SetRangeUser(-0.5,2.5);
//   hnpho_EB->GetXaxis()->SetNdivisions(5);
  TGraphErrors *tg_npho_EB = new TGraphErrors(nbin, x, 
					      npho_EB, dx, nphoerr_EB);
  tg_npho_EB->SetMarkerColor(1);
  tg_npho_EB->SetMarkerStyle(8);
  tg_npho_EB->SetLineColor(1);
  tg_npho_EB->SetLineWidth(2);
  tg_npho_EB->Draw("p e same");

  TLegend *tleg = new TLegend(0.70, 0.75, 0.95, 0.92);
  tleg->SetHeader("Barrel");
  tleg->SetFillColor(0);
  tleg->SetShadowColor(0);
  tleg->SetBorderSize(0);
  tleg->Draw();

  c1->Print("npho_aftershowershapeID_EB.eps");
  c1->Print("npho_aftershowershapeID_EB.gif");
  c1->Print("npho_aftershowershapeID_EB.pdf");


  TCanvas* c2 = new TCanvas("c2","",500,500);
  h2->Draw();
  c2->SetLogy(1);
//   hnpho_EE->SetXTitle("Number of photon candidates per event");
//   hnpho_EE->SetYTitle("Entries");
//   hnpho_EE->SetTitle("");
//   hnpho_EE->Draw("e1");
//   hnpho_EE->GetXaxis()->SetRangeUser(-0.5,2.5);
//   hnpho_EE->GetXaxis()->SetNdivisions(5);
  TGraphErrors *tg_npho_EE = new TGraphErrors(nbin, x, 
					      npho_EE, dx, nphoerr_EE);
  tg_npho_EE->SetMarkerColor(1);
  tg_npho_EE->SetMarkerStyle(8);
  tg_npho_EE->SetLineColor(1);
  tg_npho_EE->SetLineWidth(2);
  tg_npho_EE->Draw("p e same");

  TLegend *tleg2 = new TLegend(0.70, 0.75, 0.95, 0.92);
  tleg2->SetHeader("Endcap");
  tleg2->SetFillColor(0);
  tleg2->SetShadowColor(0);
  tleg2->SetBorderSize(0);
  tleg2->Draw();


  c2->Print("npho_aftershowershapeID_EE.eps");
  c2->Print("npho_aftershowershapeID_EE.gif");
  c2->Print("npho_aftershowershapeID_EE.pdf");
  





}
