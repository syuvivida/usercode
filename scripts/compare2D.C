#include <TCut.h>
#include <selection.h>

void compare2D(TChain* t1, TChain* t2, 
	       std::string var1, std::string var2,
	       int decCode, 
	       std::string psname,
	       int xnbin, float xmin, float xmax, 
	       int ynbin, float ymin, float ymax, 
	       std::string xtitle="", std::string ytitle="",
	       TCut cut1="")
{

  // for cuts the same for barrel and endcap
  TCut allCut = cut1 + looseCut;

//   // only for data

//   // only for data
  TCut dataCut;
  if(psname.find("2360") != std::string::npos)
    dataCut= Cut_2360GeV;
  else
    dataCut= Cut_900GeV;

  if(decCode==1)
    allCut += barrelCut;
  else if(decCode==2)
    allCut += endcapCut;

//   cout << "The cuts are " << allCut << endl;
  
  gStyle->SetOptStat(0);


  float xbinwidth = (xmax-xmin)/(float)xnbin;
  float ybinwidth = (ymax-ymin)/(float)ynbin;

  char tempName[300];
  sprintf(tempName, "Candidates per %1.3lf", xbinwidth);

  
  TH2F* h1 = new TH2F("h1","",xnbin,xmin,xmax,ynbin,ymin,ymax);
  h1->SetXTitle(xtitle.data());  
  h1->SetYTitle(ytitle.data());  

  TH2F* h2 = new TH2F("h2","",xnbin,xmin,xmax,ynbin,ymin,ymax);
  h2->SetXTitle(xtitle.data());
  h2->SetYTitle(ytitle.data());  


  std::string temp = var2 + ":" + var1 + ">>h1";

  t1->Draw(temp.data(), allCut && dataCut,"colz");

  temp = var2 + ":" + var1 + ">>h2";
  t2->Draw(temp.data(), allCut,"colz");

  TCanvas* c1 = new TCanvas("c1","",1000,500);
  c1->Divide(2,1);
  c1->cd(1);
  h1->Draw("colz");
  TLatex *   tex = new TLatex(-2.778001,3.200524,"Data");
  tex->SetTextFont(42);
  tex->SetTextSize(0.06918021);
  tex->SetLineWidth(2);
  tex->Draw("same");
  c1->cd(2);
  h2->Draw("colz");
  TLatex *   tex2 = new TLatex(-2.778001,3.200524,"MC");
  tex2->SetTextFont(42);
  tex2->SetTextSize(0.06918021);
  tex2->SetLineWidth(2);
  tex2->Draw("same");

  cout << "Data entries = " << h1->GetEntries() << " MC entries = " << 
    h2->GetEntries() << endl;

  std::string dirname = 
    "/home/syu/testYenJie/CMSSW_3_3_6_patch3/src/CRAB/preprod_figures/";

  if(decCode==1) psname += "_barrel";
  else if(decCode==2) psname += "_endcap";
  else if(decCode==0) psname += "_all";

  std::string filename = dirname + psname + ".eps";
  c1->Print(filename.data());
  filename = dirname + psname + ".gif";  
  c1->Print(filename.data());  

 

}
		     
