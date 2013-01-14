#include "setTDRStyle.C"
#include "cutvalues.h"

void temp3Profiles(std::string inputfile,
		   std::string var="pf_dR_matchProb_rec"
		   )
{
  
  setTDRStyle();
  gStyle->SetOptStat(0);

  const int NHISTOS=6;
  TProfile* h[NHISTOS];

  std::string massname = inputfile;
  std::string remword  ="debugmjj_M";
  size_t pos  = massname.find(remword);
  if(pos!= std::string::npos)
    massname.replace(pos,remword.length(),"");
  std::string remword2  =".root";
  size_t pos2  = massname.find(remword2);
  if(pos2!= std::string::npos)
    massname.replace(pos2,remword2.length(),"");
  
  cout << massname << endl;

  char tempName[300];
  TFile *fmc;

  
  std::string varName[NHISTOS]={
    Form("%s0",var.data()),
    Form("%s_ass0",var.data()),
    Form("%s_either0",var.data()),
    Form("%s1",var.data()),
    Form("%s_ass1",var.data()),
    Form("%s_either1",var.data())
  };

  std::string mcName[NHISTOS]={
    "Leading: quarks from Z",
    "Leading: associated parton",
    "Leading: either",
    "Subleading: quarks from Z",
    "Subleading: associated parton",
    "Subleading: either"
  };


  fmc = TFile::Open(inputfile.data());

  for(int i=0; i<NHISTOS; i++){
    h[i]   = (TProfile*)(fmc->Get(varName[i].data()));
  }


  int COLOR[NHISTOS]={4,1,2, kGreen, kOrange, 6};
  int MARKERSTYLE[NHISTOS]={20,21,29,  4, 25,30};

  for(int i=0; i < NHISTOS; i++){

//     h[i]->SetTitle("");
    h[i]->GetXaxis()->SetNdivisions(5);
    h[i]->GetYaxis()->SetNdivisions(5);
//     h[i]->GetYaxis()->SetDecimals();
    h[i]->SetLineColor(COLOR[i]);
    h[i]->SetMarkerColor(COLOR[i]);
    h[i]->SetMarkerSize(1);
    h[i]->SetMarkerStyle(MARKERSTYLE[i]);
    h[i]->SetTitleOffset(1.0,"Y");

  }

  // if normalizing to the same area, set the scale 

  vector<float> maxArray;
  maxArray.clear();

  vector<float> minArray;
  minArray.clear();

  for(int ih=0; ih < NHISTOS; ih++){

    float max_this  =  h[ih]->GetBinError(h[ih]->GetMaximumBin()) + h[ih]->GetMaximum();
    maxArray.push_back(max_this);

    float min_this  = -h[ih]->GetBinError(h[ih]->GetMinimumBin()) + h[ih]->GetMinimum();
    minArray.push_back(min_this);

  }


  float max = *max_element(maxArray.begin(),maxArray.end());
  cout << "Max = " << max << endl;

  float min = *min_element(minArray.begin(),minArray.end());
  cout << "Min = " << min << endl;


  for(int ih=0; ih < NHISTOS; ih++){
    h[ih]->SetYTitle("Matching probability");
    h[ih]->SetTitle("");
    h[ih]->SetMaximum(1.8);
    h[ih]->SetMinimum(0.0);
  }

  cout << "here" << endl;
  TCanvas* c1 = new TCanvas("c1","",600,500);  
  gStyle->SetTitleX(0.1);
  h[0]->Draw("e");
  for(int ih=0; ih < NHISTOS-1; ih++)
    h[ih]->Draw("esame");

  h[2]->Draw("histesame");
  h[5]->Draw("histesame");


  float x1NDC = 0.449664;
  float y1NDC = 0.618644;
  float x2NDC = 0.669463;
  float y2NDC = 0.913136;

  TLegend* leg = new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  leg->SetHeader(Form("m_{H} = %s GeV",massname.data()));
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  for(int i=0; i < NHISTOS; i++)
    leg->AddEntry(h[i], mcName[i].data());
  leg->Draw("same");


  string dirName = "matchProb";
  gSystem->mkdir(dirName.data());



  std::string filename;
  std::string psname = dirName + "/mH" + massname + "_" + var;
  
  filename = psname + ".eps";
  c1->Print(filename.data());
  filename = psname + ".gif";
  c1->Print(filename.data());
  filename = psname + ".pdf";
  c1->Print(filename.data());
  //   c1->Close();
}
		     
