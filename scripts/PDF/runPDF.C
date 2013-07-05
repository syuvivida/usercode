void runPDF(std::string className, std::string fileName, int PID)

{
  TStopwatch* myWatch = new TStopwatch();
  myWatch->Start();
  gROOT->ProcessLine(".L /afs/cern.ch/work/s/syu/jetsub/LHAPDF/lib/libLHAPDF.so");

  //  gROOT->ProcessLine(Form(".L %s.C++",className.data()));
  gROOT->ProcessLine(Form(".L %s_C.so",className.data()));

  gROOT->ProcessLine(Form("%s a(\"%s\");",className.data(),fileName.data()));

  gROOT->ProcessLine(Form("a.Loop(%d);",PID));

  myWatch->Stop();
  cout << myWatch->RealTime() << " seconds has passed! " << endl; 

}
