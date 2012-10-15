#include "/afs/cern.ch/user/s/syu/scripts/setTDRStyle.C"

void plotAllSys(std::string var)
{
  setTDRStyle();
  gStyle->SetOptStat(0);

  const double ymax=1.2;
  const double ymin=0.9;

  TH1D* hjes_stat;
  TH1D* hjes_syst;
  TH1D* hpu_stat;
  TH1D* hpu_syst;

  TH1D* heff;
  TH1D* hbkg;
  TH1D* hstat;

  std::string xtitle;
  double xmax=2.39999;
  double xmin=0.0;
  if(var=="h_ystar")
    {
      xtitle = "0.5|Y_{Z}-Y_{jet}|";	    
      xmax=1.799999;
    }
  else if(var=="h_yB")
    {
      xtitle = "0.5|Y_{Z}+Y_{jet}|";
      xmax=2.1999999;
    }
  else if(var=="h_jety")
    {
      xtitle = "|Y_{jet}|";
    }
  else if(var=="h_zy")
    {
      xtitle = "|Y_{Z}|";
      xmax = 2.199999;
    }
 
  // first get the histogram files
  TFile *fin  = TFile::Open("20121015_allsys/electron_forSteve.root");  
  cout << "Opening " << fin->GetName() << endl;

  std::string var1=var;
  std::string remword1 ="h_";
  size_t pos1  = var1.find(remword1);
  if(pos1!= std::string::npos)
    var1.replace(pos1,remword1.length(),"");
  

  hjes_stat = (TH1D*)(fin->Get(Form("jes_%s_corr_stat",var1.data())));
  hjes_syst = (TH1D*)(fin->Get(Form("jes_%s_corr_sys",var1.data())));

  hpu_stat = (TH1D*)(fin->Get(Form("pu_%s_corr_stat",var1.data())));
  hpu_syst = (TH1D*)(fin->Get(Form("pu_%s_corr_sys",var1.data())));

  heff     = (TH1D*)(fin->Get(Form("eff_%s_corr",var1.data())));
  hbkg     = (TH1D*)(fin->Get(Form("bkg_%s_corr",var1.data())));
  hstat    = (TH1D*)(fin->Get(Form("stat_%s_corr",var1.data())));

  heff->SetMinimum(ymin);
  heff->SetMaximum(ymax);

  hstat->SetMinimum(ymin);
  hstat->SetMaximum(ymax);

  hbkg->SetMinimum(ymin);
  hbkg->SetMaximum(ymax);

  TCanvas* c1 = new TCanvas("c1","",600,500);

  hstat->SetYTitle("Corrections");
  hstat->SetFillStyle(3006);
  hstat->SetFillColor(kViolet);
  hstat->SetLineColor(kViolet);
  hstat->SetMarkerColor(kViolet);
  hstat->SetMarkerSize(1e-4);
  hstat->GetXaxis()->SetRangeUser(xmin,xmax);
  hstat->GetXaxis()->SetNdivisions(5);
  hstat->GetXaxis()->SetDecimals();
  hstat->GetYaxis()->SetDecimals();
  hstat->Draw("ce3");

  heff->SetFillStyle(3350);
  heff->SetFillColor(5);
  heff->SetLineColor(5);
  heff->SetMarkerColor(5);
  heff->SetMarkerSize(1e-4);
  heff->GetXaxis()->SetRangeUser(xmin,xmax);
  heff->GetXaxis()->SetNdivisions(5);
  heff->Draw("ce3same");


  hjes_syst->SetFillStyle(3325);
  hjes_syst->SetFillColor(2);
  hjes_syst->SetLineColor(2);
  hjes_syst->SetMarkerColor(2);
  hjes_syst->SetMarkerSize(1e-4);
  hjes_syst->GetXaxis()->SetRangeUser(xmin,xmax);
  hjes_syst->GetXaxis()->SetNdivisions(5);
  hjes_syst->Draw("ce3same");

  hpu_syst->SetFillStyle(3395);
  hpu_syst->SetFillColor(4);
  hpu_syst->SetLineColor(4);
  hpu_syst->SetMarkerColor(4);
  hpu_syst->SetMarkerSize(1e-4);
  hpu_syst->GetXaxis()->SetRangeUser(xmin,xmax);
  hpu_syst->GetXaxis()->SetNdivisions(5);
  hpu_syst->Draw("ce3same");


  hbkg->SetFillStyle(3544);
  hbkg->SetFillColor(3);
  hbkg->SetLineColor(3);
  hbkg->SetMarkerColor(3);
  hbkg->SetMarkerSize(1e-4);
  hbkg->GetXaxis()->SetRangeUser(xmin,xmax);
  hbkg->GetXaxis()->SetNdivisions(5);
  hbkg->Draw("ce3same");

  float x1NDC = 0.208054;
  float y1NDC = 0.588983;
  float x2NDC = 0.419463;
  float y2NDC = 0.900424;

  TLegend* leg = new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  leg->SetBorderSize(0);
  leg->AddEntry(hstat,"Stat","f");
  leg->AddEntry(heff, "Eff","f");
  leg->AddEntry(hjes_syst, "JES","f");
  leg->AddEntry(hpu_syst, "PU","f");
  leg->AddEntry(hbkg, "Bkg","f");
  leg->Draw("same"); 

  c1->Print(Form("AllEleSys_%s.eps",var1.data()));
  c1->Print(Form("AllEleSys_%s.gif",var1.data()));
  c1->Print(Form("AllEleSys_%s.pdf",var1.data()));

}
		     
