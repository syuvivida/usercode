#include <TLorentzVector.h>
#include <TMath.h>

Bool_t separated(TLorentzVector l1, TLorentzVector l2);
Double_t zgamma(TLorentzVector l1, TLorentzVector l2);
Double_t deltaPhi(TLorentzVector l1, TLorentzVector l2);
Double_t cosThetaStar(TLorentzVector l1, TLorentzVector l2);
Double_t chiPair(TLorentzVector l1, TLorentzVector l2);
Double_t pstar(TLorentzVector l1, TLorentzVector l2);
Double_t yB(TLorentzVector l1, TLorentzVector l2);


Bool_t separated(TLorentzVector l1, TLorentzVector l2) { 
  Bool_t result = false;
  if(l1.DeltaR(l2)>0.5)result=true; 
  return result;
}

// variable proposed by JETPHOX authors
Double_t zgamma(TLorentzVector l1, TLorentzVector l2) { 
  Double_t dphi = l1.Phi()-l2.Phi();
  while (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
  while (dphi <= -TMath::Pi()) dphi += 2*TMath::Pi();

  Double_t result = fabs(l2.Pt())>1e-6? - (l1.Pt()/l2.Pt())*TMath::Cos(dphi) : -999.0;

  return result;
}

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


// some CM frame variables
Double_t pstar(TLorentzVector l1, TLorentzVector l2) { 

  Double_t ystar = 0.5*fabs(l1.Rapidity()-l2.Rapidity());  
  Double_t result = l1.Pt()*TMath::CosH(ystar);
  return result;

}


Double_t yB(TLorentzVector l1, TLorentzVector l2) { 

  Double_t result = 0.5*(l1.Rapidity()+l2.Rapidity());  
  return result;

}


