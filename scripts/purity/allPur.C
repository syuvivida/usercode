
void allPur() {
  gROOT->ProcessLine(".L Laurent_IsoPur_7TeVMC.C++");
  gROOT->ProcessLine("Laurent_IsoPur_7TeVMC iso");
  gROOT->ProcessLine("iso.Loop(1)");
  gROOT->ProcessLine("iso.Loop(2)");
  gROOT->ProcessLine("iso.Loop(3)");
  //gROOT->ProcessLine(".! kate isoOut.log &");
  //gROOT->ProcessLine(".q");
}










