{
  TStopwatch* myWatch = new TStopwatch();
  myWatch->Start();
  gROOT->ProcessLine(".L combineMC_eff.C++");

  gROOT->ProcessLine(
"combineMC_eff\(\"h_yB_EB_85_95\",\"y^{boost}\",2,-2.4,2.4)");

  myWatch->Stop();
  cout << myWatch->RealTime() << " seconds has passed! " << endl; 

}
