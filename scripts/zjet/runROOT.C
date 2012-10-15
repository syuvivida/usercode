void runROOT(std::string className)

{
  TStopwatch* myWatch = new TStopwatch();
  myWatch->Start();
  gROOT->ProcessLine(Form(".L %s.C",className.data()));

  gROOT->ProcessLine(Form("%s(\"h_zy\",false);",className.data()));
  gROOT->ProcessLine(Form("%s(\"h_jety\",true);",className.data()));
  gROOT->ProcessLine(Form("%s(\"h_ystar\",true);",className.data()));
  gROOT->ProcessLine(Form("%s(\"h_yB\",true);",className.data()));


  myWatch->Stop();
  cout << myWatch->RealTime() << " seconds has passed! " << endl; 

}
