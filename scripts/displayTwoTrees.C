void displayTwoTrees(std::string file1, std::string file2, 
		     std::string tree1, std::string var1,
		     std::string cut1="",
		     std::string tree2="", std::string var2="",
		     std::string cut2="")
{
  gROOT->Reset();
  gStyle->SetOptStat(111111);
  if(tree2 == ""  )tree2=tree1;
  if(var2 ==  "" )var2=var1;
  if(cut2 == "")cut2=cut1;
  
  TChain* t1 = new TChain(tree1.data());
  t1->Add(file1.data());

  TChain* t2 = new TChain(tree2.data());
  t2->Add(file2.data());
  
  TCanvas* c1 = new TCanvas("c1", file1.data(), 0, 0, 500,500);

  std::string temp = var1 + ">>h1";

  t1->Draw(temp.data(), cut1.data());
  h1->SetLineColor(4);
  h1->SetTitle(file1.data());


  cout << "h1 mean = " << h1->GetMean() << " and width = " << h1->GetRMS() << 
    " and entries = " << h1->GetEntries() << endl;

  TCanvas* c2 = new TCanvas("c2", file2.data(),500,0,500,500);
  temp = var2 + ">>h2";
  t2->Draw(temp.data(), cut2.data());
  h2->SetLineColor(2);
  h2->SetTitle(file2.data());
  h2->Draw();
  cout << "h2 mean = " << h2->GetMean() << " and width = " << h2->GetRMS() << 
    " and entries = " << h2->GetEntries() << endl;

  h1->SetMarkerSize(0);
  h2->SetMarkerSize(0);
  

  TCanvas* c3 = new TCanvas("c3","",0,0,1000,1000);
  c3->Divide(2,2);
  c3->cd(1);
  h1->Draw();
  c3->cd(2);
  h2->Draw();
  c3->cd(3);

  TLegend* leg = new TLegend(0.2385,0.6593,0.38940,0.8761) ; 
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(h1,file1.data());
  leg->AddEntry(h2,file2.data());
  leg->Draw("same");

  std::string outputfile = "test.gif";
  c3->Print(outputfile.data());

  gSystem->Sleep(2000);

  delete h1;
  delete h2;


}
		     
