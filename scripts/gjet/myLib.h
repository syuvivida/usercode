#include <TLorentzVector.h>
#include <TMath.h>

Double_t deltaPhi(TLorentzVector l1, TLorentzVector l2) { 
  Double_t result = l1.Phi()-l2.Phi();
  while (result > TMath::Pi()) result -= 2*TMath::Pi();
  while (result <= -TMath::Pi()) result += 2*TMath::Pi();
  return result;
}

Double_t cosThetaStar(TLorentzVector l1, TLorentzVector l2) { 

  Double_t ystar = 0.5*fabs(l1.Rapidity()-l2.Rapidity());
  Double_t result = TMath::TanH(ystar);
  return result;

}

Double_t chiPair(TLorentzVector l1, TLorentzVector l2) { 

  Double_t ystar = 0.5*fabs(l1.Rapidity()-l2.Rapidity());
  Double_t result = TMath::Exp(2*ystar);
  return result;
}
