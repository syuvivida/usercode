#include <TH2.h>
#include <TH1.h>
#include <TFile.h>
#include <iostream>
#include <string>
#include <TCanvas.h>

using namespace std;

void get_template_all(int mctest=0,const char* filename="proj_comb.root"){


  TFile *f = new TFile(filename);   
  TFile *f1 = new TFile(filename);
  cout << "projecting " << f->GetName() << endl;
  cout << "projecting " << f1->GetName() << endl;

  TH2F *h_EB_sig = (TH2F*)f->Get("h_EB_comb3Iso_et_sig_sum_SIG");
  TH2F *h_EB_bkg = (TH2F*)f->Get("h_EB_comb3Iso_et_bkg_sum_BKG");
  TH2F *h_EE_sig = (TH2F*)f->Get("h_EE_comb3Iso_et_sig_sum_SIG");
  TH2F *h_EE_bkg = (TH2F*)f->Get("h_EE_comb3Iso_et_bkg_sum_BKG");
  //for data
  TH2F *h_EB_EGdata = (TH2F*)f1->Get("h_EB_comb3Iso_et_EGdata_SIG");
  TH2F *h_EE_EGdata = (TH2F*)f1->Get("h_EE_comb3Iso_et_EGdata_SIG");
  TH2F *h_EB_EGdataSB = (TH2F*)f1->Get("h_EB_comb3IsoSB_et_EGdata_SIG");
  TH2F *h_EE_EGdataSB = (TH2F*)f1->Get("h_EE_comb3IsoSB_et_EGdata_SIG");
  if ( mctest==1 ) {
    TH2F *h_EB_EGdata = (TH2F*)f1->Get("h_EB_comb3Iso_et_sig_sum_SIG");
    TH2F *h_EE_EGdata = (TH2F*)f1->Get("h_EE_comb3Iso_et_sig_sum_SIG");
    TH2F *h_EB_EGdataSB = (TH2F*)f1->Get("h_EB_comb3IsoSB_et_bkg_sum_BKG");
    TH2F *h_EE_EGdataSB = (TH2F*)f1->Get("h_EE_comb3IsoSB_et_bkg_sum_BKG");
  }

  //pt15-20 bin4
  h_EB_sig->ProjectionX("h_EB_comb3Iso_sig_pt15", 4,4);
  h_EB_bkg->ProjectionX("h_EB_comb3Iso_bkg_pt15", 4,4);
  h_EB_EGdata->ProjectionX("h_EB_comb3Iso_EGdata_pt15", 4,4);  
  h_EB_EGdataSB->ProjectionX("h_EB_comb3IsoSB_EGdata_pt15", 4,60);  
  TH1F *h_EB_comb3Iso_sig_pt15 = (TH1F*)h_EB_comb3Iso_sig_pt15->Clone();
  TH1F *h_EB_comb3Iso_bkg_pt15 = (TH1F*)h_EB_comb3Iso_bkg_pt15->Clone();
  TH1F *h_EB_comb3Iso_EGdata_pt15 = (TH1F*)h_EB_comb3Iso_EGdata_pt15->Clone();
  TH1F *h_EB_comb3IsoSB_EGdata_pt15 = (TH1F*)h_EB_comb3IsoSB_EGdata_pt15->Clone();


  //pt20-30 bin5-6
  h_EB_sig->ProjectionX("h_EB_comb3Iso_sig_pt20", 5,6);
  h_EB_bkg->ProjectionX("h_EB_comb3Iso_bkg_pt20", 5,6);
  h_EB_EGdata->ProjectionX("h_EB_comb3Iso_EGdata_pt20", 5,6);  
  h_EB_EGdataSB->ProjectionX("h_EB_comb3IsoSB_EGdata_pt20", 5,60);  
  TH1F *h_EB_comb3Iso_sig_pt20 = (TH1F*)h_EB_comb3Iso_sig_pt20->Clone();
  TH1F *h_EB_comb3Iso_bkg_pt20 = (TH1F*)h_EB_comb3Iso_bkg_pt20->Clone();
  TH1F *h_EB_comb3Iso_EGdata_pt20 = (TH1F*)h_EB_comb3Iso_EGdata_pt20->Clone();
  TH1F *h_EB_comb3IsoSB_EGdata_pt20 = (TH1F*)h_EB_comb3IsoSB_EGdata_pt20->Clone();

  //pt30-50 bin7-10
  h_EB_sig->ProjectionX("h_EB_comb3Iso_sig_pt30", 7,10);
  h_EB_bkg->ProjectionX("h_EB_comb3Iso_bkg_pt30", 7,10);
  h_EB_EGdata->ProjectionX("h_EB_comb3Iso_EGdata_pt30", 7,10);  
  h_EB_EGdataSB->ProjectionX("h_EB_comb3IsoSB_EGdata_pt30", 7,60);  
  TH1F *h_EB_comb3Iso_sig_pt30 = (TH1F*)h_EB_comb3Iso_sig_pt30->Clone();
  TH1F *h_EB_comb3Iso_bkg_pt30 = (TH1F*)h_EB_comb3Iso_bkg_pt30->Clone();
  TH1F *h_EB_comb3Iso_EGdata_pt30 = (TH1F*)h_EB_comb3Iso_EGdata_pt30->Clone();
  TH1F *h_EB_comb3IsoSB_EGdata_pt30 = (TH1F*)h_EB_comb3IsoSB_EGdata_pt30->Clone();

  //pt50-80 bin11-16
  h_EB_sig->ProjectionX("h_EB_comb3Iso_sig_pt50", 11,16);
  h_EB_bkg->ProjectionX("h_EB_comb3Iso_bkg_pt50", 11,16);
  h_EB_EGdata->ProjectionX("h_EB_comb3Iso_EGdata_pt50", 11,16);  
  h_EB_EGdataSB->ProjectionX("h_EB_comb3IsoSB_EGdata_pt50", 11,60);  
  TH1F *h_EB_comb3Iso_sig_pt50 = (TH1F*)h_EB_comb3Iso_sig_pt50->Clone();
  TH1F *h_EB_comb3Iso_bkg_pt50 = (TH1F*)h_EB_comb3Iso_bkg_pt50->Clone();
  TH1F *h_EB_comb3Iso_EGdata_pt50 = (TH1F*)h_EB_comb3Iso_EGdata_pt50->Clone();
  TH1F *h_EB_comb3IsoSB_EGdata_pt50 = (TH1F*)h_EB_comb3IsoSB_EGdata_pt50->Clone();


  //pt80-120 bin17-24
  h_EB_sig->ProjectionX("h_EB_comb3Iso_sig_pt80", 17,24);
  h_EB_bkg->ProjectionX("h_EB_comb3Iso_bkg_pt80", 17,24);
  h_EB_EGdata->ProjectionX("h_EB_comb3Iso_EGdata_pt80", 17,24);  
  h_EB_EGdataSB->ProjectionX("h_EB_comb3IsoSB_EGdata_pt80", 17,60);  
  TH1F *h_EB_comb3Iso_sig_pt80 = (TH1F*)h_EB_comb3Iso_sig_pt80->Clone();
  TH1F *h_EB_comb3Iso_bkg_pt80 = (TH1F*)h_EB_comb3Iso_bkg_pt80->Clone();
  TH1F *h_EB_comb3Iso_EGdata_pt80 = (TH1F*)h_EB_comb3Iso_EGdata_pt80->Clone();
  TH1F *h_EB_comb3IsoSB_EGdata_pt80 = (TH1F*)h_EB_comb3IsoSB_EGdata_pt80->Clone();

  //pt120-170 bin25-34
  h_EB_sig->ProjectionX("h_EB_comb3Iso_sig_pt120", 25,34);
  h_EB_bkg->ProjectionX("h_EB_comb3Iso_bkg_pt120", 25,34);
  h_EB_EGdata->ProjectionX("h_EB_comb3Iso_EGdata_pt120", 25,34);  
  h_EB_EGdataSB->ProjectionX("h_EB_comb3IsoSB_EGdata_pt120", 25,60);  
  TH1F *h_EB_comb3Iso_sig_pt120 = (TH1F*)h_EB_comb3Iso_sig_pt120->Clone();
  TH1F *h_EB_comb3Iso_bkg_pt120 = (TH1F*)h_EB_comb3Iso_bkg_pt120->Clone();
  TH1F *h_EB_comb3Iso_EGdata_pt120 = (TH1F*)h_EB_comb3Iso_EGdata_pt120->Clone();
  TH1F *h_EB_comb3IsoSB_EGdata_pt120 = (TH1F*)h_EB_comb3IsoSB_EGdata_pt120->Clone();


  //pt170-230 bin35-46
  h_EB_sig->ProjectionX("h_EB_comb3Iso_sig_pt170", 35,46);
  h_EB_bkg->ProjectionX("h_EB_comb3Iso_bkg_pt170", 35,46);
  h_EB_EGdata->ProjectionX("h_EB_comb3Iso_EGdata_pt170", 35,46);  
  h_EB_EGdataSB->ProjectionX("h_EB_comb3IsoSB_EGdata_pt170", 35,60);  
  TH1F *h_EB_comb3Iso_sig_pt170 = (TH1F*)h_EB_comb3Iso_sig_pt170->Clone();
  TH1F *h_EB_comb3Iso_bkg_pt170 = (TH1F*)h_EB_comb3Iso_bkg_pt170->Clone();
  TH1F *h_EB_comb3Iso_EGdata_pt170 = (TH1F*)h_EB_comb3Iso_EGdata_pt170->Clone();
  TH1F *h_EB_comb3IsoSB_EGdata_pt170 = (TH1F*)h_EB_comb3IsoSB_EGdata_pt170->Clone();


  //pt170-300 bin47-60
  h_EB_sig->ProjectionX("h_EB_comb3Iso_sig_pt230", 47,60);
  h_EB_bkg->ProjectionX("h_EB_comb3Iso_bkg_pt230", 47,60);
  h_EB_EGdata->ProjectionX("h_EB_comb3Iso_EGdata_pt230", 47,60);  
  h_EB_EGdataSB->ProjectionX("h_EB_comb3IsoSB_EGdata_pt230", 47,60);  
  TH1F *h_EB_comb3Iso_sig_pt230 = (TH1F*)h_EB_comb3Iso_sig_pt230->Clone();
  TH1F *h_EB_comb3Iso_bkg_pt230 = (TH1F*)h_EB_comb3Iso_bkg_pt230->Clone();
  TH1F *h_EB_comb3Iso_EGdata_pt230 = (TH1F*)h_EB_comb3Iso_EGdata_pt230->Clone();
  TH1F *h_EB_comb3IsoSB_EGdata_pt230 = (TH1F*)h_EB_comb3IsoSB_EGdata_pt230->Clone();

  //-------------------------------------------------------------------------------------



  //pt15-20 bin4
  h_EE_sig->ProjectionX("h_EE_comb3Iso_sig_pt15", 4,4);
  h_EE_bkg->ProjectionX("h_EE_comb3Iso_bkg_pt15", 4,4);
  h_EE_EGdata->ProjectionX("h_EE_comb3Iso_EGdata_pt15", 4,4);  
  h_EE_EGdataSB->ProjectionX("h_EE_comb3IsoSB_EGdata_pt15", 4, 4);  
  TH1F *h_EE_comb3Iso_sig_pt15 = (TH1F*)h_EE_comb3Iso_sig_pt15->Clone();
  TH1F *h_EE_comb3Iso_bkg_pt15 = (TH1F*)h_EE_comb3Iso_bkg_pt15->Clone();
  TH1F *h_EE_comb3Iso_EGdata_pt15 = (TH1F*)h_EE_comb3Iso_EGdata_pt15->Clone();
  TH1F *h_EE_comb3IsoSB_EGdata_pt15 = (TH1F*)h_EE_comb3IsoSB_EGdata_pt15->Clone();

  //   h_EE_sig->RebinY(2);
  //   h_EE_bkg->RebinY(2);
  //   h_EE_EGdata->RebinY(2);

  //pt20-30 bin5-6
  h_EE_sig->ProjectionX("h_EE_comb3Iso_sig_pt20", 5,6);
  h_EE_bkg->ProjectionX("h_EE_comb3Iso_bkg_pt20", 5,6);
  h_EE_EGdata->ProjectionX("h_EE_comb3Iso_EGdata_pt20", 5,6);  
  h_EE_EGdataSB->ProjectionX("h_EE_comb3IsoSB_EGdata_pt20", 5,6 );  
  TH1F *h_EE_comb3Iso_sig_pt20 = (TH1F*)h_EE_comb3Iso_sig_pt20->Clone();
  TH1F *h_EE_comb3Iso_bkg_pt20 = (TH1F*)h_EE_comb3Iso_bkg_pt20->Clone();
  TH1F *h_EE_comb3Iso_EGdata_pt20 = (TH1F*)h_EE_comb3Iso_EGdata_pt20->Clone();
  TH1F *h_EE_comb3IsoSB_EGdata_pt20 = (TH1F*)h_EE_comb3IsoSB_EGdata_pt20->Clone();

  //pt30-50 bin7-10
  h_EE_sig->ProjectionX("h_EE_comb3Iso_sig_pt30", 7,10);
  h_EE_bkg->ProjectionX("h_EE_comb3Iso_bkg_pt30", 7,10);
  h_EE_EGdata->ProjectionX("h_EE_comb3Iso_EGdata_pt30", 7,10);  
  h_EE_EGdataSB->ProjectionX("h_EE_comb3IsoSB_EGdata_pt30", 7,10);  
  TH1F *h_EE_comb3Iso_sig_pt30 = (TH1F*)h_EE_comb3Iso_sig_pt30->Clone();
  TH1F *h_EE_comb3Iso_bkg_pt30 = (TH1F*)h_EE_comb3Iso_bkg_pt30->Clone();
  TH1F *h_EE_comb3Iso_EGdata_pt30 = (TH1F*)h_EE_comb3Iso_EGdata_pt30->Clone();
  TH1F *h_EE_comb3IsoSB_EGdata_pt30 = (TH1F*)h_EE_comb3IsoSB_EGdata_pt30->Clone();

  //pt50-80 bin11-16
  h_EE_sig->ProjectionX("h_EE_comb3Iso_sig_pt50", 11,16);
  h_EE_bkg->ProjectionX("h_EE_comb3Iso_bkg_pt50", 11,16);
  h_EE_EGdata->ProjectionX("h_EE_comb3Iso_EGdata_pt50", 11,16);  
  h_EE_EGdataSB->ProjectionX("h_EE_comb3IsoSB_EGdata_pt50", 11,16);  
  TH1F *h_EE_comb3Iso_sig_pt50 = (TH1F*)h_EE_comb3Iso_sig_pt50->Clone();
  TH1F *h_EE_comb3Iso_bkg_pt50 = (TH1F*)h_EE_comb3Iso_bkg_pt50->Clone();
  TH1F *h_EE_comb3Iso_EGdata_pt50 = (TH1F*)h_EE_comb3Iso_EGdata_pt50->Clone();
  TH1F *h_EE_comb3IsoSB_EGdata_pt50 = (TH1F*)h_EE_comb3IsoSB_EGdata_pt50->Clone();


  //pt80-120 bin17-24
  h_EE_sig->ProjectionX("h_EE_comb3Iso_sig_pt80", 17,24);
  h_EE_bkg->ProjectionX("h_EE_comb3Iso_bkg_pt80", 17,24);
  h_EE_EGdata->ProjectionX("h_EE_comb3Iso_EGdata_pt80", 17,24);  
  h_EE_EGdataSB->ProjectionX("h_EE_comb3IsoSB_EGdata_pt80", 17,24);  
  TH1F *h_EE_comb3Iso_sig_pt80 = (TH1F*)h_EE_comb3Iso_sig_pt80->Clone();
  TH1F *h_EE_comb3Iso_bkg_pt80 = (TH1F*)h_EE_comb3Iso_bkg_pt80->Clone();
  TH1F *h_EE_comb3Iso_EGdata_pt80 = (TH1F*)h_EE_comb3Iso_EGdata_pt80->Clone();
  TH1F *h_EE_comb3IsoSB_EGdata_pt80 = (TH1F*)h_EE_comb3IsoSB_EGdata_pt80->Clone();

  //pt120-170 bin25-34
  h_EE_sig->ProjectionX("h_EE_comb3Iso_sig_pt120", 25,34);
  h_EE_bkg->ProjectionX("h_EE_comb3Iso_bkg_pt120", 25,34);
  h_EE_EGdata->ProjectionX("h_EE_comb3Iso_EGdata_pt120", 25,34);  
  h_EE_EGdataSB->ProjectionX("h_EE_comb3IsoSB_EGdata_pt120", 25,34);  
  TH1F *h_EE_comb3Iso_sig_pt120 = (TH1F*)h_EE_comb3Iso_sig_pt120->Clone();
  TH1F *h_EE_comb3Iso_bkg_pt120 = (TH1F*)h_EE_comb3Iso_bkg_pt120->Clone();
  TH1F *h_EE_comb3Iso_EGdata_pt120 = (TH1F*)h_EE_comb3Iso_EGdata_pt120->Clone();
  TH1F *h_EE_comb3IsoSB_EGdata_pt120 = (TH1F*)h_EE_comb3IsoSB_EGdata_pt120->Clone();


  //pt170-230 bin35-46
  h_EE_sig->ProjectionX("h_EE_comb3Iso_sig_pt170", 35,46);
  h_EE_bkg->ProjectionX("h_EE_comb3Iso_bkg_pt170", 35,46);
  h_EE_EGdata->ProjectionX("h_EE_comb3Iso_EGdata_pt170", 35,46);  
  h_EE_EGdataSB->ProjectionX("h_EE_comb3IsoSB_EGdata_pt170", 35,46);  
  TH1F *h_EE_comb3Iso_sig_pt170 = (TH1F*)h_EE_comb3Iso_sig_pt170->Clone();
  TH1F *h_EE_comb3Iso_bkg_pt170 = (TH1F*)h_EE_comb3Iso_bkg_pt170->Clone();
  TH1F *h_EE_comb3Iso_EGdata_pt170 = (TH1F*)h_EE_comb3Iso_EGdata_pt170->Clone();
  TH1F *h_EE_comb3IsoSB_EGdata_pt170 = (TH1F*)h_EE_comb3IsoSB_EGdata_pt170->Clone();


  //pt170-300 bin47-60
  h_EE_sig->ProjectionX("h_EE_comb3Iso_sig_pt230", 47,60);
  h_EE_bkg->ProjectionX("h_EE_comb3Iso_bkg_pt230", 47,60);
  h_EE_EGdata->ProjectionX("h_EE_comb3Iso_EGdata_pt230", 47,60);  
  h_EE_EGdataSB->ProjectionX("h_EE_comb3IsoSB_EGdata_pt230", 47,60);  
  TH1F *h_EE_comb3Iso_sig_pt230 = (TH1F*)h_EE_comb3Iso_sig_pt230->Clone();
  TH1F *h_EE_comb3Iso_bkg_pt230 = (TH1F*)h_EE_comb3Iso_bkg_pt230->Clone();
  TH1F *h_EE_comb3Iso_EGdata_pt230 = (TH1F*)h_EE_comb3Iso_EGdata_pt230->Clone();
  TH1F *h_EE_comb3IsoSB_EGdata_pt230 = (TH1F*)h_EE_comb3IsoSB_EGdata_pt230->Clone();


  if(mctest==1) {
    h_EB_comb3Iso_EGdata_pt15->Add(h_EB_comb3Iso_bkg_pt15);
    h_EB_comb3Iso_EGdata_pt20->Add(h_EB_comb3Iso_bkg_pt20);
    h_EB_comb3Iso_EGdata_pt30->Add(h_EB_comb3Iso_bkg_pt30);
    h_EB_comb3Iso_EGdata_pt50->Add(h_EB_comb3Iso_bkg_pt50);
    h_EB_comb3Iso_EGdata_pt80->Add(h_EB_comb3Iso_bkg_pt80);
    h_EB_comb3Iso_EGdata_pt120->Add(h_EB_comb3Iso_bkg_pt120);
    h_EB_comb3Iso_EGdata_pt170->Add(h_EB_comb3Iso_bkg_pt170);
    h_EB_comb3Iso_EGdata_pt230->Add(h_EB_comb3Iso_bkg_pt230);

    h_EE_comb3Iso_EGdata_pt15->Add(h_EE_comb3Iso_bkg_pt15);
    h_EE_comb3Iso_EGdata_pt20->Add(h_EE_comb3Iso_bkg_pt20);
    h_EE_comb3Iso_EGdata_pt30->Add(h_EE_comb3Iso_bkg_pt30);
    h_EE_comb3Iso_EGdata_pt50->Add(h_EE_comb3Iso_bkg_pt50);
    h_EE_comb3Iso_EGdata_pt80->Add(h_EE_comb3Iso_bkg_pt80);
    h_EE_comb3Iso_EGdata_pt120->Add(h_EE_comb3Iso_bkg_pt120);
    h_EE_comb3Iso_EGdata_pt170->Add(h_EE_comb3Iso_bkg_pt170);
    h_EE_comb3Iso_EGdata_pt230->Add(h_EE_comb3Iso_bkg_pt230);
  }



  TCanvas *c1 = new TCanvas("c1","",900,600);
  c1->Divide(3,2);
  c1->Draw();

  float ymax=0;

  c1->cd(1);
  ymax = h_EB_comb3Iso_sig_pt15->GetMaximum();
  if(h_EB_comb3Iso_bkg_pt15->GetMaximum() > ymax ) {
    ymax = h_EB_comb3Iso_bkg_pt15->GetMaximum();
  }
  h_EB_comb3Iso_bkg_pt15->SetMaximum(ymax*1.2);
  h_EB_comb3Iso_bkg_pt15->Draw("hist e");
  h_EB_comb3Iso_sig_pt15->SetLineColor(4);
  h_EB_comb3Iso_sig_pt15->Draw("hist e same");


  c1->cd(2);
  ymax = h_EB_comb3Iso_sig_pt20->GetMaximum();
  if(h_EB_comb3Iso_bkg_pt20->GetMaximum() > ymax ) {
    ymax = h_EB_comb3Iso_bkg_pt20->GetMaximum();
  }
  h_EB_comb3Iso_bkg_pt20->SetMaximum(ymax*1.2);
  h_EB_comb3Iso_bkg_pt20->Draw("hist e");
  h_EB_comb3Iso_sig_pt20->SetLineColor(4);
  h_EB_comb3Iso_sig_pt20->Draw("hist e same");

  c1->cd(3);
  ymax = h_EB_comb3Iso_sig_pt30->GetMaximum();
  if(h_EB_comb3Iso_bkg_pt30->GetMaximum() > ymax ) {
    ymax = h_EB_comb3Iso_bkg_pt30->GetMaximum();
  }
  h_EB_comb3Iso_bkg_pt30->SetMaximum(ymax*1.2);
  h_EB_comb3Iso_bkg_pt30->Draw("hist e");
  h_EB_comb3Iso_sig_pt30->SetLineColor(4);
  h_EB_comb3Iso_sig_pt30->Draw("hist e same");

  c1->cd(4);
  ymax = h_EB_comb3Iso_sig_pt50->GetMaximum();
  if(h_EB_comb3Iso_bkg_pt50->GetMaximum() > ymax ) {
    ymax = h_EB_comb3Iso_bkg_pt50->GetMaximum();
  }
  h_EB_comb3Iso_bkg_pt50->SetMaximum(ymax*1.2);
  h_EB_comb3Iso_bkg_pt50->Draw("hist e");
  h_EB_comb3Iso_sig_pt50->SetLineColor(4);
  h_EB_comb3Iso_sig_pt50->Draw("hist e same");

  c1->cd(5);
  ymax = h_EB_comb3Iso_sig_pt80->GetMaximum();
  if(h_EB_comb3Iso_bkg_pt80->GetMaximum() > ymax ) {
    ymax = h_EB_comb3Iso_bkg_pt80->GetMaximum();
  }
  h_EB_comb3Iso_bkg_pt80->SetMaximum(ymax*1.2);
  h_EB_comb3Iso_bkg_pt80->Draw("hist e");
  h_EB_comb3Iso_sig_pt80->SetLineColor(4);
  h_EB_comb3Iso_sig_pt80->Draw("hist e same");

  c1->cd(6);
  ymax = h_EB_comb3Iso_sig_pt120->GetMaximum();
  if(h_EB_comb3Iso_bkg_pt120->GetMaximum() > ymax ) {
    ymax = h_EB_comb3Iso_bkg_pt120->GetMaximum();
  }
  h_EB_comb3Iso_bkg_pt120->SetMaximum(ymax*1.2);
  h_EB_comb3Iso_bkg_pt120->Draw("hist e");
  h_EB_comb3Iso_sig_pt120->SetLineColor(4);
  h_EB_comb3Iso_sig_pt120->Draw("hist e same");

  TCanvas *c2 = new TCanvas("c2","",900,800);
  c2->Divide(3,2);
  c2->Draw();

  c2->cd(1);
  ymax = h_EE_comb3Iso_sig_pt15->GetMaximum();
  if(h_EE_comb3Iso_bkg_pt15->GetMaximum() > ymax ) {
    ymax = h_EE_comb3Iso_bkg_pt15->GetMaximum();
  }
  h_EE_comb3Iso_bkg_pt15->SetMaximum(ymax*1.2);
  h_EE_comb3Iso_bkg_pt15->Draw("hist e");
  h_EE_comb3Iso_sig_pt15->SetLineColor(4);
  h_EE_comb3Iso_sig_pt15->Draw("hist e same");


  c2->cd(2);
  ymax = h_EE_comb3Iso_sig_pt20->GetMaximum();
  if(h_EE_comb3Iso_bkg_pt20->GetMaximum() > ymax ) {
    ymax = h_EE_comb3Iso_bkg_pt20->GetMaximum();
  }
  h_EE_comb3Iso_bkg_pt20->SetMaximum(ymax*1.2);
  h_EE_comb3Iso_bkg_pt20->Draw("hist e");
  h_EE_comb3Iso_sig_pt20->SetLineColor(4);
  h_EE_comb3Iso_sig_pt20->Draw("hist e same");

  c2->cd(3);
  ymax = h_EE_comb3Iso_sig_pt30->GetMaximum();
  if(h_EE_comb3Iso_bkg_pt30->GetMaximum() > ymax ) {
    ymax = h_EE_comb3Iso_bkg_pt30->GetMaximum();
  }
  h_EE_comb3Iso_bkg_pt30->SetMaximum(ymax*1.2);
  h_EE_comb3Iso_bkg_pt30->Draw("hist e");
  h_EE_comb3Iso_sig_pt30->SetLineColor(4);
  h_EE_comb3Iso_sig_pt30->Draw("hist e same");

  c2->cd(4);
  ymax = h_EE_comb3Iso_sig_pt50->GetMaximum();
  if(h_EE_comb3Iso_bkg_pt50->GetMaximum() > ymax ) {
    ymax = h_EE_comb3Iso_bkg_pt50->GetMaximum();
  }
  h_EE_comb3Iso_bkg_pt50->SetMaximum(ymax*1.2);
  h_EE_comb3Iso_bkg_pt50->Draw("hist e");
  h_EE_comb3Iso_sig_pt50->SetLineColor(4);
  h_EE_comb3Iso_sig_pt50->Draw("hist e same");

  c2->cd(5);
  ymax = h_EE_comb3Iso_sig_pt80->GetMaximum();
  if(h_EE_comb3Iso_bkg_pt80->GetMaximum() > ymax ) {
    ymax = h_EE_comb3Iso_bkg_pt80->GetMaximum();
  }
  h_EE_comb3Iso_bkg_pt80->SetMaximum(ymax*1.2);
  h_EE_comb3Iso_bkg_pt80->Draw("hist e");
  h_EE_comb3Iso_sig_pt80->SetLineColor(4);
  h_EE_comb3Iso_sig_pt80->Draw("hist e same");

  c2->cd(6);
  ymax = h_EE_comb3Iso_sig_pt120->GetMaximum();
  if(h_EE_comb3Iso_bkg_pt120->GetMaximum() > ymax ) {
    ymax = h_EE_comb3Iso_bkg_pt120->GetMaximum();
  }
  h_EE_comb3Iso_bkg_pt120->SetMaximum(ymax*1.2);
  h_EE_comb3Iso_bkg_pt120->Draw("hist e");
  h_EE_comb3Iso_sig_pt120->SetLineColor(4);
  h_EE_comb3Iso_sig_pt120->Draw("hist e same");

  std::string dataMC = mctest==0? "Integral_SBData" : "Integral_SBMC";
  TFile *fout = new TFile(Form("%sTemplate_%s",dataMC.data(),filename),"recreate");
  cout << "producing " << fout->GetName() << endl;

  h_EB_comb3Iso_sig_pt15->Write();
  h_EB_comb3Iso_bkg_pt15->Write();
  h_EB_comb3Iso_EGdata_pt15 ->Write();
  h_EB_comb3IsoSB_EGdata_pt15 ->Write();
                          

  h_EB_comb3Iso_sig_pt20->Write();
  h_EB_comb3Iso_bkg_pt20->Write();
  h_EB_comb3Iso_EGdata_pt20 ->Write();
  h_EB_comb3IsoSB_EGdata_pt20 ->Write();
                          

  h_EB_comb3Iso_sig_pt30->Write();
  h_EB_comb3Iso_bkg_pt30->Write();
  h_EB_comb3Iso_EGdata_pt30 ->Write();
  h_EB_comb3IsoSB_EGdata_pt30 ->Write();
                          
  h_EB_comb3Iso_sig_pt50->Write();
  h_EB_comb3Iso_bkg_pt50->Write();
  h_EB_comb3Iso_EGdata_pt50->Write();
  h_EB_comb3IsoSB_EGdata_pt50->Write();
                          
  h_EB_comb3Iso_sig_pt80->Write();
  h_EB_comb3Iso_bkg_pt80->Write();
  h_EB_comb3Iso_EGdata_pt80->Write();
  h_EB_comb3IsoSB_EGdata_pt80->Write();
                         
  h_EB_comb3Iso_sig_pt120->Write();
  h_EB_comb3Iso_bkg_pt120->Write();
  h_EB_comb3Iso_EGdata_pt120->Write();
  h_EB_comb3IsoSB_EGdata_pt120->Write();

  h_EB_comb3Iso_sig_pt170->Write();
  h_EB_comb3Iso_bkg_pt170->Write();
  h_EB_comb3Iso_EGdata_pt170->Write();
  h_EB_comb3IsoSB_EGdata_pt170->Write();

  h_EB_comb3Iso_sig_pt230->Write();
  h_EB_comb3Iso_bkg_pt230->Write();
  h_EB_comb3Iso_EGdata_pt230->Write();
  h_EB_comb3IsoSB_EGdata_pt230->Write();


  h_EE_comb3Iso_sig_pt15->Write();
  h_EE_comb3Iso_bkg_pt15->Write();
  h_EE_comb3Iso_EGdata_pt15->Write();
  h_EE_comb3IsoSB_EGdata_pt15->Write();

  h_EE_comb3Iso_sig_pt20->Write();
  h_EE_comb3Iso_bkg_pt20->Write();
  h_EE_comb3Iso_EGdata_pt20->Write();
  h_EE_comb3IsoSB_EGdata_pt20->Write();

  h_EE_comb3Iso_sig_pt30->Write();
  h_EE_comb3Iso_bkg_pt30->Write();
  h_EE_comb3Iso_EGdata_pt30->Write();
  h_EE_comb3IsoSB_EGdata_pt30->Write();
                          
  h_EE_comb3Iso_sig_pt50->Write();
  h_EE_comb3Iso_bkg_pt50->Write();
  h_EE_comb3Iso_EGdata_pt50->Write();
  h_EE_comb3IsoSB_EGdata_pt50->Write();
                          
  h_EE_comb3Iso_sig_pt80->Write();
  h_EE_comb3Iso_bkg_pt80->Write();
  h_EE_comb3Iso_EGdata_pt80->Write();
  h_EE_comb3IsoSB_EGdata_pt80->Write();
                          
  h_EE_comb3Iso_sig_pt120->Write();
  h_EE_comb3Iso_bkg_pt120->Write();
  h_EE_comb3Iso_EGdata_pt120->Write();
  h_EE_comb3IsoSB_EGdata_pt120->Write();

  h_EE_comb3Iso_sig_pt170->Write();
  h_EE_comb3Iso_bkg_pt170->Write();
  h_EE_comb3Iso_EGdata_pt170->Write();
  h_EE_comb3IsoSB_EGdata_pt170->Write();

  h_EE_comb3Iso_sig_pt230->Write();
  h_EE_comb3Iso_bkg_pt230->Write();
  h_EE_comb3Iso_EGdata_pt230->Write();
  h_EE_comb3IsoSB_EGdata_pt230->Write();

  fout->Close();



}
