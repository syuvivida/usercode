void exec_Dir(std::string dataDir, std::string filestring, 
	       bool isBarrel, bool applyR9)
{

  TChain* pho = new TChain("RECOTrigger/root");

  TSystemDirectory *base = new TSystemDirectory("root","root");
  std::string outerDir = "/home/syu/Photon_Code/CMSSW_3_1_2/src/syu/MyObjectCounter/test/" + dataDir + "/";
  std::string filename = outerDir+ "/res";

  base->SetDirectory(filename.data());
  TList *listOfFiles = base->GetListOfFiles();
  TIter fileIt(listOfFiles);
  TFile *fileH = new TFile();
  int nfile=0;
  while(fileH = ((TFile*)fileIt())) {
    std::string fileN = fileH->GetName();
    std::string baseString = "root";
    if( fileH->IsFolder())  continue;
    if(fileN.find(baseString) == std::string::npos)continue;
    cout << fileN.data() << endl;
    nfile++;
    pho->Add(fileN.data());
  }

  std::cout << "Opened " << nfile << " files" << std::endl;


  gROOT->ProcessLine(".L /home/syu/Photon_Code/CMSSW_3_1_2/src/syu/MyObjectCounter/test/TrigEff_withL3Match.C++");  
  TrigEff_withL3Match(pho,filestring, isBarrel, applyR9);
  gROOT->ProcessLine(".L /home/syu/Photon_Code/CMSSW_3_1_2/src/syu/MyObjectCounter/test/TrigEff_AsymmetryErrors.C++");
  gROOT->ProcessLine(".L /home/syu/Photon_Code/CMSSW_3_1_2/src/syu/MyObjectCounter/test/easyMatchTrigEff.C++");  
  easyMatchTrigEff(filestring, isBarrel, applyR9);

  TChain* gen = new TChain("GenTrig/root");

  TIter fileIt2(listOfFiles);
  TFile *fileH2 = new TFile();
  int nfile2=0;
  while(fileH2 = ((TFile*)fileIt2())) {
    std::string fileN = fileH2->GetName();
    std::string baseString = "root";
    if( fileH2->IsFolder())  continue;
    if(fileN.find(baseString) == std::string::npos)continue;
    cout << fileN.data() << endl;
    nfile2++;
    gen->Add(fileN.data());
  }

  std::cout << "Opened " << nfile2 << " files" << std::endl;

  filestring = "gen" + filestring;
  gROOT->ProcessLine(".L /home/syu/Photon_Code/CMSSW_3_1_2/src/syu/MyObjectCounter/test/GenTrigEff_withL3Match.C++");  
  GenTrigEff_withL3Match(gen,filestring, isBarrel);
  easyMatchTrigEff(filestring, isBarrel, false);


}
