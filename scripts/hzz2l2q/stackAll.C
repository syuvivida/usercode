{
  TStopwatch* myWatch = new TStopwatch();
  myWatch->Start();

  std::string inputFile = "hzz2l2qMC_muon.txt";

  gROOT->ProcessLine(".L stackHisto.C++");

  gROOT->ProcessLine(Form("stackHisto(\"h_nvtx0\",\"%s\",false)",
			  inputFile.data()));
  gROOT->ProcessLine(Form("stackHisto(\"h_nvtx1\",\"%s\",true)",
			  inputFile.data()));
  gROOT->ProcessLine(Form("stackHisto(\"h_nvtx2\",\"%s\",true)",
			  inputFile.data()));

  gROOT->ProcessLine(Form("stackHisto(\"h_eleRho0\",\"%s\",true)",
			  inputFile.data()));
  gROOT->ProcessLine(Form("stackHisto(\"h_eleRho1\",\"%s\",true)",
			  inputFile.data()));
  gROOT->ProcessLine(Form("stackHisto(\"h_eleRho2\",\"%s\",true)",
			  inputFile.data()));


  gROOT->ProcessLine(Form("stackHisto(\"h_muoRho0\",\"%s\",true)",
			  inputFile.data()));
  gROOT->ProcessLine(Form("stackHisto(\"h_muoRho1\",\"%s\",true)",
			  inputFile.data()));
  gROOT->ProcessLine(Form("stackHisto(\"h_muoRho2\",\"%s\",true)",
			  inputFile.data()));

  


  myWatch->Stop();
  cout << myWatch->RealTime() << " seconds has passed! " << endl; 

}
