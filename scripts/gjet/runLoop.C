
void runLoop(std::string className, std::string fileName)
{
  TStopwatch* myWatch = new TStopwatch();
  myWatch->Start();
  std::string command = ".L " + className + ".C++";
  gROOT->ProcessLine(command.data());

  command = className + " a(\"" + fileName + "\");";
  gROOT->ProcessLine(command.data());

  command = "a.Loop(true)";
  gROOT->ProcessLine(command.data());
  myWatch->Stop();
  cout << myWatch->RealTime() << " seconds has passed! " << endl; 

}
