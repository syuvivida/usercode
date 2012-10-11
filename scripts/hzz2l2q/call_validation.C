#include "displayTwoTreesRatio.C"
void call_validation()
{
  TStopwatch* myWatch = new TStopwatch();
  myWatch->Start();

//   std::string inFile = "DoubleElectron_Run2012A_filteredtree.root";
  
//   std::string inFile = "GluGluToHToZZTo2L2Q_M-400_8TeV-powheg-pythia6_filteredtree.root";
  std::string inFile = "DoubleMu_Run2012A_filteredtree.root";

  /*
  displayTwoTreesRatio(inFile.data(),"genParM_",100,0,1000,false,"test",
  		       "genParId_==25");
   displayTwoTreesRatio(inFile.data(),"PU_nTrueInt",60,0,60);

  displayTwoTreesRatio(inFile.data(),"EvtInfo_VertexX_[0]",100,0.23,0.26,false,"EvtInfo_VertexX_");
  displayTwoTreesRatio(inFile.data(),"EvtInfo_VertexY_[0]",100,0.38,0.41,false,"EvtInfo_VertexY_");

//   displayTwoTreesRatio(inFile.data(),"EvtInfo_VertexX_[0]",100,0.05,0.10,false,
// 		       "EvtInfo_VertexX_");
//   displayTwoTreesRatio(inFile.data(),"EvtInfo_VertexY_[0]",100,0.04,0.08,false,
// 		       "EvtInfo_VertexY_");
  
  displayTwoTreesRatio(inFile.data(),"EvtInfo_EventNum",100,0,100000000);
  displayTwoTreesRatio(inFile.data(),"EvtInfo_RunNum",100,150000,300000);
  displayTwoTreesRatio(inFile.data(),"EvtInfo_VertexZ_[0]",100,-30,30,false,
  "EvtInfo_VertexZ_");
  */
  displayTwoTreesRatio(inFile.data(),"EvtInfo_NumVtx",30,0,30);

  /*
  displayTwoTreesRatio(inFile.data(),"patJetPfAk05Pt_",50,0, 500);
  displayTwoTreesRatio(inFile.data(),"patJetPfAk05Eta_",50,-3.0,3.0);
  displayTwoTreesRatio(inFile.data(),"patJetPfAk05Phi_",50,0,TMath::Pi());
  displayTwoTreesRatio(inFile.data(),"patJetPfAk05CharMulti_",100,0,100);
  displayTwoTreesRatio(inFile.data(),"eleRho",100,0,70);
  displayTwoTreesRatio(inFile.data(),"muoRho",100,0,10);

  displayTwoTreesRatio(inFile.data(),"higgsM",50,0,1500);
  displayTwoTreesRatio(inFile.data(),"higgsPt",50, 0,500);
  displayTwoTreesRatio(inFile.data(),"higgsEta",50, -6,6);
  displayTwoTreesRatio(inFile.data(),"higgsPhi",50, 0,TMath::Pi());

  displayTwoTreesRatio(inFile.data(),"zllPt",50, 0,500);
  displayTwoTreesRatio(inFile.data(),"zllEta",50, -6,6);
  displayTwoTreesRatio(inFile.data(),"zllPhi",50, 0,TMath::Pi());

  displayTwoTreesRatio(inFile.data(),"zjjPt",50, 0,500);
  displayTwoTreesRatio(inFile.data(),"zjjEta",50, -6,6);
  displayTwoTreesRatio(inFile.data(),"zjjPhi",50, 0,TMath::Pi());

  displayTwoTreesRatio(inFile.data(),"zllM",100,60,120);
  displayTwoTreesRatio(inFile.data(),"zjjM",100, 0,300);

  displayTwoTreesRatio(inFile.data(),"nBTags",3,-0.5,2.5);
  displayTwoTreesRatio(inFile.data(),"lepType",2,-0.5,1.5);
  displayTwoTreesRatio(inFile.data(),"heliLD",100,0.0,1.0);
  */
  displayTwoTreesRatio(inFile.data(),"metSig", 50,0, 50,true);
  /*
  displayTwoTreesRatio(inFile.data(),"patJetPfAk05NeutEmEFr_",100,0,1,true);
  displayTwoTreesRatio(inFile.data(),"patJetPfAk05NeutHadEFr_",100,0,1,true);
  displayTwoTreesRatio(inFile.data(),"patJetPfAk05CharEmEFr_",100,0,1,true);
  displayTwoTreesRatio(inFile.data(),"patJetPfAk05CharHadEFr_",100,0,1,true);

  displayTwoTreesRatio(inFile.data(),"jetPt",50,0, 500);
  displayTwoTreesRatio(inFile.data(),"jetEta",50,-3.0,3.0);
  displayTwoTreesRatio(inFile.data(),"jetPhi",50,0,TMath::Pi());
  */
  myWatch->Stop();
  cout << myWatch->RealTime() << " seconds has passed! " << endl; 

}
