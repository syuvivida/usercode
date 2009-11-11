#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TChain.h>
#include <TBranch.h>
#include <iostream.h>
#include <fstream.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystemDirectory.h>
#include <TProfile.h>					
#include <TList.h>
#include <TF1.h>
#include <vector>
#include <TLorentzVector.h>
#include "TrigEff_AsymmetryErrors.C"
#include <TGraphAsymmErrors.h>

using namespace SY_NT;

void easyTrigEff(std::string prefix, std::string dirsuffix = "triggerdata")
{

//   const int NHISTOS = 18;
//   const int NHISTOS = 22;

  const int NHISTOS = 28;
  const int NEFF = NHISTOS/2;

  TH1F* horiginal[NHISTOS];

  TFile* histoFile;
 
  std::string outerDir = "/home/syu/Photon_Code/CMSSW_3_1_2/src/syu/MyObjectCounter/test/" + dirsuffix + "/";

  std::string fileName = outerDir + prefix + ".root";
 
  cout << "opening " << fileName << endl;
  histoFile = TFile::Open(fileName.data());
  
  char name[300];
  std::string title[NEFF];
  
  title[0]="HLT_Photon15_L1R/HLT_L1SingleEG5 (all reconstructed photons)";
  title[1]="HLT_Photon15_L1R/HLT_L1SingleEG5 (photons from hard scattering)";
  title[2]="HLT_Photon15_L1R/HLT_L1SingleEG5 (photons from jets)";
  title[3]="HLT_Photon15_L1R/HLT_L1SingleEG5 (photons from quark radiation)";
  title[4]="HLT_Photon15_L1R/HLT_L1SingleEG5 (photons from gluon)";
  title[5]="HLT_Photon15_L1R/HLT_L1SingleEG5 (generator photon)";
  title[6]="HLT_Photon15_L1R/HLT_Mu5 (all reconstructed photon)";
  title[7]="HLT_Photon25_L1R/HLT_L1SingleEG5 (all reconstructed photon)";
  title[8]="HLT_Photon20_Iso/HLT_L1SingleEG5 (all reconstructed photon)";
  title[9]="HLT_Photon15_L1R/HLT_L1SingleEG5 (photons from hard scattering)";
  title[10]="HLT_L1SingleEG5/HLT_Mu5 (all reconstructed photons)";
  title[11]="HLT_Photon10_L1R/HLT_L1SingleEG5 (photons from hard scattering)";
  title[12]="HLT_L1SingleEG5/reconstructed (photons from hard scattering)";
  title[13]="HLT_L1SingleEG5/reconstructed (photons from jets)";
  
  horiginal[0] = (TH1F*)(histoFile->Get("h_recgetdeno"));
  horiginal[1] = (TH1F*)(histoFile->Get("h_recgetnumr"));
  horiginal[2] = (TH1F*)(histoFile->Get("h_getdeno"));
  horiginal[3] = (TH1F*)(histoFile->Get("h_getnumr"));
  horiginal[4] = (TH1F*)(histoFile->Get("h_jetgetdeno"));
  horiginal[5] = (TH1F*)(histoFile->Get("h_jetgetnumr"));
  horiginal[6] = (TH1F*)(histoFile->Get("h_qgetdeno"));
  horiginal[7] = (TH1F*)(histoFile->Get("h_qgetnumr"));
  horiginal[8] = (TH1F*)(histoFile->Get("h_ggetdeno"));
  horiginal[9] = (TH1F*)(histoFile->Get("h_ggetnumr"));
  horiginal[10] = (TH1F*)(histoFile->Get("h_gengetdeno"));
  horiginal[11] = (TH1F*)(histoFile->Get("h_gengetnumr"));
  horiginal[12] = (TH1F*)(histoFile->Get("h_mugetdeno"));
  horiginal[13] = (TH1F*)(histoFile->Get("h_mugetnumr"));
  horiginal[14] = (TH1F*)(histoFile->Get("h_25allgetdeno"));
  horiginal[15] = (TH1F*)(histoFile->Get("h_25allgetnumr"));
  horiginal[16] = (TH1F*)(histoFile->Get("h_isoallgetdeno"));
  horiginal[17] = (TH1F*)(histoFile->Get("h_isoallgetnumr"));
  horiginal[18] = (TH1F*)(histoFile->Get("h_etadeno"));
  horiginal[19] = (TH1F*)(histoFile->Get("h_etadeno"));
  horiginal[20] = (TH1F*)(histoFile->Get("h_l1mugetdeno"));
  horiginal[21] = (TH1F*)(histoFile->Get("h_l1mugetnumr"));

  horiginal[22] = (TH1F*)(histoFile->Get("h_get10deno"));
  horiginal[23] = (TH1F*)(histoFile->Get("h_get10numr"));

  horiginal[24] = (TH1F*)(histoFile->Get("h_get5deno"));
  horiginal[25] = (TH1F*)(histoFile->Get("h_get5numr"));

  horiginal[26] = (TH1F*)(histoFile->Get("h_jetget5deno"));
  horiginal[27] = (TH1F*)(histoFile->Get("h_jetget5numr"));


  TH1F* hdeno[NEFF];
  TH1F* hnumr[NEFF];

  for(int k=0; k < NEFF   ; k++)
    { 
      int deno = k*2;
      int numr = k*2+1;
      hdeno[k] = (TH1F*)horiginal[deno]->Clone();
      sprintf(name,"hdeno%02i",k);
      hdeno[k]->SetName(name);

      hnumr[k] = (TH1F*)horiginal[numr]->Clone();
      sprintf(name,"hnumr%02i",k);
      hnumr[k]->SetName(name);
    }
   
  
  TGraphAsymmErrors *heff[NEFF];
  for(int i=0; i< NEFF;i++)
    {
      heff[i] = MyDivide(hdeno[i],hnumr[i]);  
      sprintf(name,"heff%02i",i);
      heff[i]->SetName(name);
      heff[i]->SetTitle(title[i].data());
      heff[i]->GetXaxis()->SetTitle("Photon E_{T} [GeV]");
      heff[i]->GetYaxis()->SetTitle("Trigger efficiency");                     
      heff[i]->Draw("ap");
    }

  std::string outFileName = outerDir + prefix + "_trigeff.root";
  TFile* outFile = new TFile(outFileName.data(),"recreate");
  for(int i=0;i< NEFF;i++)
    heff[i]->Write();
  
  outFile->Close();
  
}


