{
  TStopwatch* myWatch = new TStopwatch();
  myWatch->Start();

  std::string inputFile = "hzz2l2qMC_muon.txt";

  gROOT->ProcessLine(".L compareDataStackMC.C++");

  gROOT->ProcessLine(Form("compareDataStackMC(\"h_nvtx0\",\"%s\",false,0,40)",
			  inputFile.data()));
  gROOT->ProcessLine(Form("compareDataStackMC(\"h_nvtx1\",\"%s\",true,0,40)",
			  inputFile.data()));
  gROOT->ProcessLine(Form("compareDataStackMC(\"h_nvtx2\",\"%s\",true,0,40)",
			  inputFile.data()));
  gROOT->ProcessLine(Form("compareDataStackMC(\"h_nvtx3\",\"%s\",true,0,40)",
			  inputFile.data()));

  gROOT->ProcessLine(Form("compareDataStackMC(\"h_eleRho0\",\"%s\",true,0,40)",
			  inputFile.data()));
  gROOT->ProcessLine(Form("compareDataStackMC(\"h_eleRho1\",\"%s\",true,0,40)",
			  inputFile.data()));
  gROOT->ProcessLine(Form("compareDataStackMC(\"h_eleRho2\",\"%s\",true,0,40)",
			  inputFile.data()));
  gROOT->ProcessLine(Form("compareDataStackMC(\"h_eleRho3\",\"%s\",true,0,40)",
			  inputFile.data()));


  gROOT->ProcessLine(Form("compareDataStackMC(\"h_muoRho0\",\"%s\",true,0,6)",
			  inputFile.data()));
  gROOT->ProcessLine(Form("compareDataStackMC(\"h_muoRho1\",\"%s\",true,0,6)",
			  inputFile.data()));
  gROOT->ProcessLine(Form("compareDataStackMC(\"h_muoRho2\",\"%s\",true,0,6)",
			  inputFile.data()));
  gROOT->ProcessLine(Form("compareDataStackMC(\"h_muoRho3\",\"%s\",true,0,6)",
			  inputFile.data()));
  


  myWatch->Stop();
  cout << myWatch->RealTime() << " seconds has passed! " << endl; 

}
