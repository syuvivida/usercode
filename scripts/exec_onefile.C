void exec_onefile(std::string inputfile, std::string filestring, 
	       bool isBarrel, bool applyR9)
{

  TChain* pho = new TChain("RECOTrigger/root");
  pho->Add(inputfile.data());
  gROOT->ProcessLine(".L /home/syu/Photon_Code/CMSSW_3_1_2/src/syu/MyObjectCounter/test/TrigEff_withL3Match.C++");  
  TrigEff_withL3Match(pho,filestring, isBarrel, applyR9);
  gROOT->ProcessLine(".L /home/syu/Photon_Code/CMSSW_3_1_2/src/syu/MyObjectCounter/test/TrigEff_AsymmetryErrors.C++");
  gROOT->ProcessLine(".L /home/syu/Photon_Code/CMSSW_3_1_2/src/syu/MyObjectCounter/test/easyMatchTrigEff.C++");  
  easyMatchTrigEff(filestring, isBarrel, applyR9);


  //generator-level matching
  TChain* gen = new TChain("GenTrig/root");
  gen->Add(inputfile.data());
  filestring = "gen" + filestring;
  gROOT->ProcessLine(".L /home/syu/Photon_Code/CMSSW_3_1_2/src/syu/MyObjectCounter/test/GenTrigEff_withL3Match.C++");  
  GenTrigEff_withL3Match(gen,filestring, isBarrel);
  easyMatchTrigEff(filestring, isBarrel, false);
}
