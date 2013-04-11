#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Math/Vector4.hh"
#include "Rivet/Projections/Thrust.hh"

namespace Rivet {

  class CMS_SMP_12_004_Z : public Analysis {
  public:

    //Constructor:
    CMS_SMP_12_004_Z()
      : Analysis("CMS_SMP_12_004_Z")
    {
      setBeams(PROTON, PROTON);
      setNeedsCrossSection(true);
    }


    void init()
    {

      //full final state
      const FinalState fs(-5.0,5.0);
      addProjection(fs, "FS");
      //Z finders for electrons and muons
      const ZFinder zfe(fs, -2.1, 2.1, 20*GeV, 11, 76*GeV, 106.*GeV, 0.1, true, false);
      const ZFinder zfm(fs, -2.1, 2.1, 20*GeV, 13, 76*GeV, 106.*GeV, 0.1, true, false);
      addProjection(zfe, "ZFE");
      addProjection(zfm, "ZFM");
      //jets  
      const FastJets jets(fs, FastJets::ANTIKT, 0.5);
      addProjection(jets, "JETS");

      //Histograms with data
      _histYZ             = bookHistogram1D(1, 1, 1);
      _histYJet           = bookHistogram1D(2, 1, 1); 
      _histYSum           = bookHistogram1D(3, 1, 1); 
      _histYDif           = bookHistogram1D(4, 1, 1); 
 
    }

    void analyze(const Event& event){
      const double weight = event.weight();
      //apply the Z finders
      const ZFinder& zfe = applyProjection<ZFinder>(event, "ZFE");
      const ZFinder& zfm = applyProjection<ZFinder>(event, "ZFM");

      //if no Z found, veto
      if (zfe.empty() && zfm.empty())
        vetoEvent;

      //Choose the Z candidate
      const ParticleVector& z = !zfm.empty() ? zfm.bosons() : zfe.bosons();
      const ParticleVector& clusteredConstituents = !zfm.empty() ? zfm.constituents() : zfe.constituents();
      
      //determine whether we are in boosted regime
      if (z[0].momentum().pT()<40*GeV)
        vetoEvent;

      //build the jets
      const FastJets& jetfs = applyProjection<FastJets>(event, "JETS");
      Jets jets = jetfs.jetsByPt(30.*GeV, MAXDOUBLE, -2.4, 2.4);

      //clean the jets against the lepton candidates, as in the paper, with a DeltaR cut of 0.4 against the clustered leptons
      std::vector<const Jet*> cleanedJets;
      for (unsigned int i = 0; i < jets.size(); ++i){
        bool isolated = true; 
        for (unsigned j = 0; j < clusteredConstituents.size(); ++j){
          if (deltaR(clusteredConstituents[j].momentum().vector3(), jets[i].momentum().vector3()) < 0.5){
            isolated=false;
            break;
          }
        }
        if (isolated)
          cleanedJets.push_back(&jets[i]);
      }

      unsigned int Njets = cleanedJets.size();
      //require at least 1 jet
      if (Njets !=1)
        vetoEvent;

      double yz   = z[0].momentum().rapidity();
      double yjet = cleanedJets[0]->momentum().rapidity();
      double ysum = 0.5*(yz+yjet);
      double ydif = 0.5*(yz-yjet);

      _histYZ     ->fill(fabs(yz));
      _histYJet   ->fill(fabs(yjet));
      _histYSum   ->fill(fabs(ysum));
      _histYDif   ->fill(fabs(ydif));      

    }

    void normalizeNoOverflows(AIDA::IHistogram1D* plot, double integral){
      double factor=1.;
      if (plot->sumBinHeights()>0 && plot->sumAllBinHeights()>0)
        factor = plot->sumAllBinHeights()/plot->sumBinHeights();
      normalize(plot, factor*integral);  
    }

    void finalize() 
    {
      normalizeNoOverflows(_histYZ         ,1.);
      normalizeNoOverflows(_histYJet       ,1.);
      normalizeNoOverflows(_histYSum         ,1.);
      normalizeNoOverflows(_histYDif         ,1.);
    }
    
  private:

    AIDA::IHistogram1D* _histYZ;             
    AIDA::IHistogram1D* _histYJet;           
    AIDA::IHistogram1D* _histYSum;
    AIDA::IHistogram1D* _histYDif;
  };
  
  AnalysisBuilder<CMS_SMP_12_004_Z> plugin_CMS_SMP_12_004_Z;
  
}



   
