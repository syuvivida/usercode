#include <TH1.h>
#include <TRandom2.h>
#include <TCanvas.h>

void toy_example(){

  TH1F* h1 = new TH1F("h1","",30,-3,3);
  TRandom2* gndm = new TRandom2();

  for(int i=0;i<1000;i++)
    {
      float x = -3 + 6*gndm->Rndm();
      h1->Fill(x);
    }

  TCanvas* c1 = new TCanvas("c1","",500,500);
  h1->SetMinimum(0);
  h1->Draw();
  c1->Print("h1.pdf");

}
