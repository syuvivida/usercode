#include <TLorentzVector.h>
#include <TMath.h>

namespace eiko {

  Bool_t separated(TLorentzVector l1, TLorentzVector l2);
  Double_t zgamma(TLorentzVector l1, TLorentzVector l2);
  Double_t deltaPhi(TLorentzVector l1, TLorentzVector l2);
  Double_t deltaR(double eta1, double phi1, double eta2, double phi2);
  Double_t cosThetaStar(TLorentzVector l1, TLorentzVector l2);
  Double_t cosThetaStar_BoostToCM(TLorentzVector l1, TLorentzVector l2);
  Double_t cosThetaStar_ZBoostToCM(TLorentzVector l1, TLorentzVector l2);
  Double_t chiPair(TLorentzVector l1, TLorentzVector l2);
  Double_t pstar(TLorentzVector l1, TLorentzVector l2);
  Double_t pstar_BoostToCM(TLorentzVector l1, TLorentzVector l2);
  Double_t pstar_ZBoostToCM(TLorentzVector l1, TLorentzVector l2);
  Double_t yB(TLorentzVector l1, TLorentzVector l2);
  Double_t ystar(TLorentzVector l1, TLorentzVector l2);
  Double_t ystar_BoostToCM(TLorentzVector l1, TLorentzVector l2);
  Double_t ystar_ZBoostToCM(TLorentzVector l1, TLorentzVector l2);

  void BoostToCM_3D(TLorentzVector& l1, TLorentzVector& l2);
  void BoostToCM_Z(TLorentzVector& l1, TLorentzVector& l2);

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

  Double_t deltaR(double eta1, double phi1, double eta2, double phi2)
  {
    
    Double_t deta = eta1-eta2;
    double dphi = phi1-phi2;
    while (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
    while (dphi <= -TMath::Pi()) dphi += 2*TMath::Pi();
    double dR = sqrt(deta*deta+dphi*dphi);
    return dR;

  }

  Double_t deltaPhi(TLorentzVector l1, TLorentzVector l2) { 
    Double_t result = l1.Phi()-l2.Phi();
    while (result > TMath::Pi()) result -= 2*TMath::Pi();
    while (result <= -TMath::Pi()) result += 2*TMath::Pi();
    return fabs(result);
  }

  Double_t cosThetaStar(TLorentzVector l1, TLorentzVector l2) { 

    Double_t ystar = 0.5*fabs(l1.Rapidity()-l2.Rapidity());
    Double_t result = TMath::TanH(ystar);
    return result;

  }

  // boost along the sum of two momenta
  Double_t cosThetaStar_BoostToCM(TLorentzVector l1, TLorentzVector l2) { 

    TLorentzVector l1_copy(l1);
    TLorentzVector l2_copy(l2);

    BoostToCM_3D(l1_copy, l2_copy);

    Double_t result = fabs(l1_copy.CosTheta());

    return result;

  }


  // boost along the z direction of the sum of two momenta
  Double_t cosThetaStar_ZBoostToCM(TLorentzVector l1, TLorentzVector l2) { 

    TLorentzVector l1_copy(l1);
    TLorentzVector l2_copy(l2);

    BoostToCM_Z(l1_copy, l2_copy);

    Double_t result = fabs(l1_copy.CosTheta());

    return result;

  }


  Double_t chiPair(TLorentzVector l1, TLorentzVector l2) { 

    Double_t ystar = 0.5*fabs(l1.Rapidity()-l2.Rapidity());
    Double_t result = TMath::Exp(2*ystar);
    return result;
  }


  // some CM frame variables, pstar is momentum magnitude of photon in the 
  // center of mass frame pt/sinTheta*

  Double_t pstar(TLorentzVector l1, TLorentzVector l2) { 

    Double_t ystar = 0.5*fabs(l1.Rapidity()-l2.Rapidity());  
    Double_t result = l1.Pt()*TMath::CosH(ystar);
    return result;

  }

  Double_t pstar_BoostToCM(TLorentzVector l1, TLorentzVector l2){

    TLorentzVector l1_copy(l1);
    TLorentzVector l2_copy(l2);

    BoostToCM_3D(l1_copy, l2_copy);

    Double_t result = l1_copy.Rho();

    return result;


  }

  Double_t pstar_ZBoostToCM(TLorentzVector l1, TLorentzVector l2){

    TLorentzVector l1_copy(l1);
    TLorentzVector l2_copy(l2);

    BoostToCM_Z(l1_copy, l2_copy);

    Double_t result = l1_copy.Rho();

    return result;


  }


  Double_t yB(TLorentzVector l1, TLorentzVector l2) { 

    Double_t result = 0.5*(l1.Rapidity()+l2.Rapidity());  
    return result;

  }

  Double_t ystar(TLorentzVector l1, TLorentzVector l2){

    Double_t result =0.5*(l1.Rapidity()-l2.Rapidity());
    return result;

  }

  Double_t ystar_BoostToCM(TLorentzVector l1, TLorentzVector l2){

    TLorentzVector l1_copy(l1);
    TLorentzVector l2_copy(l2);

    BoostToCM_3D(l1_copy, l2_copy);

    Double_t result = l1_copy.Rapidity();

    return result;

  }

  Double_t ystar_ZBoostToCM(TLorentzVector l1, TLorentzVector l2){

    TLorentzVector l1_copy(l1);
    TLorentzVector l2_copy(l2);

    BoostToCM_Z(l1_copy, l2_copy);

    Double_t result = l1_copy.Rapidity();

    return result;

  }


  void BoostToCM_3D(TLorentzVector& l1, TLorentzVector& l2){

    TLorentzVector sum_p4 = (l1+l2);

    TVector3 boostVector = -sum_p4.BoostVector();

    l1.Boost(boostVector); 
    l2.Boost(boostVector); 

    return;
  }

  void BoostToCM_Z(TLorentzVector& l1, TLorentzVector& l2){

    TLorentzVector sum_p4 = (l1+l2);
    sum_p4.SetPx(0.0);
    sum_p4.SetPy(0.0);

    TVector3 boostVector = -sum_p4.BoostVector();

    l1.Boost(boostVector); 
    l2.Boost(boostVector); 

    return;
  }


} // end of namespace eiko
