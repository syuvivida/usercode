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

  class CMS_SMP_12_004 : public Analysis {
  public:

    //Constructor:
    CMS_SMP_12_004()
      : Analysis("CMS_SMP_12_004")
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

      //try to get Photon
      LeadingParticlesFinalState photonfs(FinalState(-2.5, 2.5, 40.0*GeV));
      photonfs.addParticleId(PHOTON);
      addProjection(photonfs, "LeadingPhoton");

      //jets  
      const FastJets jets(fs, FastJets::ANTIKT, 0.5);
      addProjection(jets, "JETS");

      //Histograms with data
      _hist1YZ             = bookHistogram1D(1, 1, 1);
      _hist1YJet           = bookHistogram1D(2, 1, 1); 
      _hist1YSum           = bookHistogram1D(3, 1, 1); 
      _hist1YDif           = bookHistogram1D(4, 1, 1); 

      _hist2YPhoton        = bookHistogram1D(5, 1, 1);
      _hist2YJet           = bookHistogram1D(6, 1, 1);
      _hist2YSum           = bookHistogram1D(7, 1, 1);
      _hist2YDif           = bookHistogram1D(8, 1, 1);

    }

    void makeZCut(const Event& event){
      

      const double weight = event.weight();
      //apply the Z finders
      const ZFinder& zfe = applyProjection<ZFinder>(event, "ZFE");
      const ZFinder& zfm = applyProjection<ZFinder>(event, "ZFM");

      //if no Z found, veto
      if (zfe.empty() && zfm.empty())
	return;

      //Choose the Z candidate
      const ParticleVector& z = !zfm.empty() ? zfm.bosons() : zfe.bosons();
      const ParticleVector& clusteredConstituents = !zfm.empty() ? zfm.constituents() : zfe.constituents();
      
      //determine whether we are in boosted regime
      if (z[0].momentum().pT()<40*GeV)
	return;
	
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
	return;
      
      double yz   = z[0].momentum().rapidity();
      double yjet = cleanedJets[0]->momentum().rapidity();
      double ysum = 0.5*(yz+yjet);
      double ydif = 0.5*(yz-yjet);

      _hist1YZ     ->fill(fabs(yz));
      _hist1YJet   ->fill(fabs(yjet));
      _hist1YSum   ->fill(fabs(ysum));
      _hist1YDif   ->fill(fabs(ydif));      

      return ;
    }


    void makePhotonCut(const Event& event){
      const double weight = event.weight();


        // Get the photon

        const FinalState& photonfs = applyProjection<FinalState>(event, "LeadingPhoton");
        if (photonfs.particles().size() < 1) {
          return;
        }
        const FourMomentum photon = photonfs.particles().front().momentum();
        double pt_P = photon.pT()/GeV;
        double eta_P = photon.eta();
        double phi_P = photon.phi();

        if(pt_P < 40 ){
            return;
        }
        if(fabs(eta_P) > 1.4442 ){
            return;
        }

      //build the jets
      const FastJets& jetfs = applyProjection<FastJets>(event, "JETS");
      Jets jets = jetfs.jetsByPt(30.*GeV, MAXDOUBLE, -2.4, 2.4);

        if (jets.size() < 1) {
          return;
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
        return;

      double yPhoton   = photon.rapidity();
      double yjet = cleanedJets[0]->momentum().rapidity();
      double ysum = 0.5*(yPhoton+yjet);
      double ydif = 0.5*(yPhoton-yjet);

      _hist2YPhoton  ->fill(fabs(yPhoton));
      _hist2YJet     ->fill(fabs(yjet));
      _hist2YSum     ->fill(fabs(ysum));
      _hist2YDif     ->fill(fabs(ydif));

      return;	
    }


    void analyze(const Event& event){
      makeZCut(event);
      makePhotonCut(event);

    }


void normalizeByContents(AIDA::IHistogram1D* plot){
      double binwidth = plot->axis().binWidth(1);
      normalize(plot, binwidth);  
    }
/* 
the function normalize(HistoID histo,double norm) is to normalize the histogram histo to area = norm . And our normalization defination in data histogram is sum over all bin-content = 1. Thus we require our MC histogram to normalize to area = area*binwidth.
*/

    void finalize() 
    {

     normalizeByContents(_hist1YZ );
     normalizeByContents(_hist1YJet );
     normalizeByContents(_hist1YSum );
     normalizeByContents(_hist1YDif );
     normalizeByContents(_hist2YPhoton );
     normalizeByContents(_hist2YJet );
     normalizeByContents(_hist2YSum );
     normalizeByContents(_hist2YDif );


    }
    
  private:
    AIDA::IHistogram1D* _hist1YZ;             
    AIDA::IHistogram1D* _hist1YJet;           
    AIDA::IHistogram1D* _hist1YSum;
    AIDA::IHistogram1D* _hist1YDif;

    AIDA::IHistogram1D* _hist2YPhoton;
    AIDA::IHistogram1D* _hist2YJet;
    AIDA::IHistogram1D* _hist2YSum;
    AIDA::IHistogram1D* _hist2YDif;

  };
  AnalysisBuilder<CMS_SMP_12_004> plugin_CMS_SMP_12_004;
  
}



   
