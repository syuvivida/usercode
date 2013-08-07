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

#include "Rivet/Projections/LeadingParticlesFinalState.hh"

namespace Rivet {


  class CMS_SMP_12_004_PhotonJet : public Analysis {
  public:

    //Constructor:
    CMS_SMP_12_004_PhotonJet()  	
    : Analysis("CMS_SMP_12_004_PhotonJet")
    {
      setBeams(PROTON, PROTON);
      setNeedsCrossSection(true);
    }


    void init()
    {

      //full final state
      const FinalState fs(-5.0,5.0);
      addProjection(fs, "FS");

      //try to get Photon
      LeadingParticlesFinalState photonfs(FinalState(-2.5, 2.5, 40.0*GeV));
      photonfs.addParticleId(PHOTON);
      addProjection(photonfs, "LeadingPhoton");

      //jets  

      const FastJets jets(fs, FastJets::ANTIKT, 0.5);
      addProjection(jets, "JETS");

      //Histograms with data
      _histYPhoton         = bookHistogram1D(1, 1, 1);
      _histYJet           = bookHistogram1D(2, 1, 1); 
      _histYSum           = bookHistogram1D(3, 1, 1); 
      _histYDif           = bookHistogram1D(4, 1, 1); 
 
    }

    void analyze(const Event& event){
      const double weight = event.weight();


        // Get the photon

        const FinalState& photonfs = applyProjection<FinalState>(event, "LeadingPhoton");
        if (photonfs.particles().size() < 1) {
          vetoEvent;
        }
        const FourMomentum photon = photonfs.particles().front().momentum();
        double pt_P = photon.pT()/GeV;
        double eta_P = photon.eta();
        double phi_P = photon.phi();

	if(pt_P < 40 ){
            vetoEvent;
        }
        if(fabs(eta_P) > 1.4442 ){
            vetoEvent;
        }

      //build the jets
      const FastJets& jetfs = applyProjection<FastJets>(event, "JETS");
      Jets jets = jetfs.jetsByPt(30.*GeV, MAXDOUBLE, -2.4, 2.4);

        if (jets.size() < 1) {
          vetoEvent;
        }

      //clean the jets against the lepton candidates, as in the paper, with a DeltaR cut of 0.4 against the clustered leptons
      std::vector<const Jet*> cleanedJets;
      for (unsigned int i = 0; i < jets.size(); ++i){
        bool isolated = true; 
	
	if (deltaR(photon.vector3(),jets[i].momentum().vector3())<0.5){
	  isolated=false;
  	}

        if (isolated)
          cleanedJets.push_back(&jets[i]);
      }

      unsigned int Njets = cleanedJets.size();
      //require at least 1 jet
      if (Njets !=1)
        vetoEvent;
                            
      double yPhoton   = photon.rapidity();
      double yjet = cleanedJets[0]->momentum().rapidity();
      double ysum = 0.5*(yPhoton+yjet);
      double ydif = 0.5*(yPhoton-yjet);

      _histYPhoton  ->fill(fabs(yPhoton));
      _histYJet     ->fill(fabs(yjet));
      _histYSum     ->fill(fabs(ysum));
      _histYDif     ->fill(fabs(ydif));      

    }

    void normalizeNoOverflows(AIDA::IHistogram1D* plot, double integral){
      double factor=1.;
      if (plot->sumBinHeights()>0 && plot->sumAllBinHeights()>0)
        factor = plot->sumAllBinHeights()/plot->sumBinHeights();
      normalize(plot, factor*integral);  
    }

    void finalize() 
    {
      normalizeNoOverflows(_histYPhoton         ,1.);
      normalizeNoOverflows(_histYJet       ,1.);
      normalizeNoOverflows(_histYSum         ,1.);
      normalizeNoOverflows(_histYDif         ,1.);
    }
    
  private:

    AIDA::IHistogram1D* _histYPhoton;             
    AIDA::IHistogram1D* _histYJet;           
    AIDA::IHistogram1D* _histYSum;
    AIDA::IHistogram1D* _histYDif;
  };
  
  AnalysisBuilder<CMS_SMP_12_004_PhotonJet> plugin_CMS_SMP_12_004_PhotonJet;
  
}



   
