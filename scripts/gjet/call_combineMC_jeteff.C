{
  TStopwatch* myWatch = new TStopwatch();
  myWatch->Start();
  gROOT->ProcessLine(".L combineMC_eff.C++");

  gROOT->ProcessLine(
"combineMC_eff(\"h_genjetpt_alljets\",\"p_{T}^{GEN}(all jets)\",5,30,200)"
);

  gROOT->ProcessLine(
"combineMC_eff(\"h_genjeteta_alljets\",\"#eta^{GEN}(all jets)\",2,-3.0,3.0)"
);

  gROOT->ProcessLine(
"combineMC_eff(\"h_recjetpt_alljets\",\"p_{T}^{REC}(all jets)\",5,30,200)"
);

  gROOT->ProcessLine(
"combineMC_eff(\"h_recjeteta_alljets\",\"#eta^{REC}(all jets)\",2,-2.4,2.4)"
);

  gROOT->ProcessLine(
"combineMC_eff(\"h_genjetpt_1stjet\",\"p_{T}^{GEN}(leading jet)\",5,30,200)"
);

  gROOT->ProcessLine(
"combineMC_eff(\"h_genjeteta_1stjet\",\"#eta^{GEN}(leading jet)\",2,-3.0,3.0)"
);

  gROOT->ProcessLine(
"combineMC_eff(\"h_recjetpt_1stjet\",\"p_{T}^{REC}(leading jet)\",5,30,200)"
);

  gROOT->ProcessLine(
"combineMC_eff(\"h_recjeteta_1stjet\",\"#eta^{REC}(leading jet)\",2,-2.4,2.4)"
);


  myWatch->Stop();
  cout << myWatch->RealTime() << " seconds has passed! " << endl; 

}
