void compareTwoTrees(std::string file1, std::string file2, 
		     std::string tree1, std::string var1,
		     int decCode, 
		     std::string psname,
		     int nbin, float xmin, float xmax, 
		     std::string xtitle="", std::string ytitle="",
		     std::string cut1="",
		     std::string tree2="", std::string var2="",
		     std::string cut2="")
{
  if(tree2 == ""  )tree2=tree1;
  if(var2 ==  "" )var2=var1;
  if(cut2 == "")cut2=cut1;
  
  gStyle->SetOptStat(0);
  TChain* t1 = new TChain(tree1.data());
  t1->Add(file1.data());

  TChain* t2 = new TChain(tree2.data());
  t2->Add(file2.data());

  float binwidth = (xmax-xmin)/(float)nbin;

  char tempName[300];
  sprintf(tempName, "Candidates per %1.3lf", binwidth);
  if(ytitle == "")ytitle = tempName;

  
  TH1F* h1 = new TH1F("h1","",nbin,xmin,xmax);
  h1->SetXTitle(xtitle.data());  
  h1->SetYTitle(ytitle.data());  

  TH1F* h2 = new TH1F("h2","",nbin,xmin,xmax);
  h2->SetXTitle(xtitle.data());
  h2->SetYTitle(ytitle.data());  

  TH1F* hscale = new TH1F("hscale","",nbin,xmin,xmax);
  hscale->SetXTitle(xtitle.data());
  hscale->SetYTitle("Data/MC");

  std::string temp = var1 + ">>h1";

  t1->Draw(temp.data(), cut1.data());
  h1->SetLineColor(4);
  h1->Draw();
  cout << "h1 mean = " << h1->GetMean() << " and width = " << h1->GetRMS() << 
    " and entries = " << h1->GetEntries() << endl;

  temp = var2 + ">>h2";
  t2->Draw(temp.data(), cut2.data());
  h2->SetLineColor(2);
  h2->Draw();
  cout << "h2 mean = " << h2->GetMean() << " and width = " << h2->GetRMS() << 
    " and entries = " << h2->GetEntries() << endl;

  
  gROOT->ProcessLine(".L computeChi2New.C");
  computeChi2New(h1,h2,hscale,psname,decCode,false);
  computeChi2New(h1,h2,hscale,psname,decCode,true);


}
		     
