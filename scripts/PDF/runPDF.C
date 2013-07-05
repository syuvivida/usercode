void runPDF(std::string className)

{
  TStopwatch* myWatch = new TStopwatch();
  myWatch->Start();
  gROOT->ProcessLine(".L /afs/cern.ch/work/s/syu/jetsub/LHAPDF/lib/libLHAPDF.so");

  gROOT->ProcessLine(Form(".L %s.C++",className.data()));

  gROOT->ProcessLine(Form("%s a;",className.data()));

  gROOT->ProcessLine("a.Loop(13);");

  myWatch->Stop();
  cout << myWatch->RealTime() << " seconds has passed! " << endl; 

}
