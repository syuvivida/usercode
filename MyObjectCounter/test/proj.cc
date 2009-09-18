#include "TTree.h"
#include "TH1F.h"
#include "format.h"

void proj()
{
	TH1F *h_pt  = new TH1F("h_pt","Electron p_{T} (GeV/c)",60.,0.,300);
        
        TFile *f1 = new TFile("results_eiko.root");
        TTree *root = (TTree*)f1->Get("MyObjectCounter/root");
	

        LepInfoBranches LepInfo;
        LepInfo.Register(root);
	
        for(int entry=0;entry<root->GetEntries();entry++) {
                root->GetEntry(entry);
                for(int i=0;i<LepInfo.Size;i++)  
		  h_pt->Fill(LepInfo.Pt[i]);
        }
	h_pt->Draw();
}
