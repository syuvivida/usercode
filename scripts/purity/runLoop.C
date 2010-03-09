
void runLoop(std::string className, std::string fileName, 
	     std::string lead, std::string match)
{

  std::string command = ".L " + className + ".C++";
  gROOT->ProcessLine(command.data());

  command = className + " a(\"" + fileName + "\");";
  gROOT->ProcessLine(command.data());

  command = "a.SetUseLeadPhoton(" + lead + ");";
  gROOT->ProcessLine(command.data());

  command = "a.SetMatching(" + match + ");";
  gROOT->ProcessLine(command.data());

  command = "a.Loop();";
  gROOT->ProcessLine(command.data());


}
