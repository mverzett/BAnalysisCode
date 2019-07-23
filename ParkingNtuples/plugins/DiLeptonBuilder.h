/*

  Class used to build and select di-lepton pairs to be later used for the 
  construction of B candidates

  original author: Mauro Verzetti (CERN)

 */

#ifndef BAnalysisCode_ParkingNtuples_DiLeptonBuilder
#define BAnalysisCode_ParkingNtuples_DiLeptonBuilder

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "helpers.h"

template<typename Lepton, typename Fitter>
class DiLeptonBuilder {
private:
  const StringCutObjectSelector<Lepton> l1_selection_; // cut on leading lepton
  const StringCutObjectSelector<Lepton> l2_selection_; // cut on sub-leading lepton
  const StringCutObjectSelector<pat::CompositeCandidate> dilepton_pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> dilepton_post_vtx_selection_; // cut on the di-lepton after the SV fit

public:
  typedef ChachedObject<Lepton> CachedLepton;
  typedef std::vector< CachedLepton > CachedLeptonCollection;
  typedef std::map< std::pair<size_t, size_t>, pat::CompositeCandidate > DiLeptonCache;

  struct DiLeptonInfo { // all pointer are not owning
    CachedLepton *l1 = 0;
    CachedLepton *l2 = 0;
    pat::CompositeCandidate *dilepton = 0;
  };
  typedef std::vector<DiLeptonInfo> DiLeptonInfoCollection;

  DiLeptonBuilder(const edm::ParameterSet& cfg):
    l1_selection_{cfg.getParameter<std::string>("l1Selection")},
    l2_selection_{cfg.getParameter<std::string>("l2Selection")},
    dilepton_pre_vtx_selection_{cfg.getParameter<std::string>("diLeptonPreVtxSelection")},
    dilepton_post_vtx_selection_{cfg.getParameter<std::string>("diLeptonPostVtxSelection")} {}

  ~DiLeptonBuilder() {};

  std::unique_ptr<DiLeptonBuilder<Lepton, Fitter>::DiLeptonInfoCollection> build(
    CachedLeptonCollection& leptons,
    DiLeptonCache& cache // It's important that is pass by ref and is not const, we are going to modify it!
    ) const {

    auto ret_value = std::make_unique<DiLeptonInfoCollection>();

    for(size_t il1 = 0; il1 < leptons.size(); ++il1) {
      Lepton* l1 = leptons[il1].obj;
      bool l1_as_l1 = l1_selection_(*l1);
      bool l1_as_l2 = l2_selection_(*l1);
      if( !(l1_as_l1 || l1_as_l2) ) continue; 
    
      for(size_t il2 = il1+1; il2 < leptons.size(); ++il2) {
        Lepton* l2 = leptons[il2].obj;
        bool l2_as_l1 = l1_selection_(*l2);
        bool l2_as_l2 = l2_selection_(*l2);
        if( !(l2_as_l1 || l2_as_l2) ) continue;

        bool l1_l2 = l1_as_l1 && l2_as_l2;
        bool l2_l1 = l2_as_l1 && l1_as_l2;
        if( !(l1_l2 || l2_l1) ) continue; // make sure both leptons pass opposite selections

        // check cache for the pair
        auto key = std::make_pair(il1, il2);
        auto iter = cache.find(key);
      
        // no luck, create the candidate and store the relevant values
        if(iter == cache.end()) { 
          cache[key]; // Does this work?
          iter = cache.find(key); // re-create the iterator to point to the correct item
          iter->second.addDaughter( *l1 );
          iter->second.addDaughter( *l2 );
          iter->second.addUserFloat("deltaR", reco::deltaR(*l1, *l2));
          // add here more stuff if needed
        }

        pat::CompositeCandidate &lepton_pair = iter->second;
        lepton_pair.setP4(l1->p4() + l2->p4());
        lepton_pair.setCharge(l1->charge() + l2->charge());
        //std::cout << "(" << il1 << ", " << il2 << "): " << lepton_pair.charge() << " / " << lepton_pair.mass() << std::endl;
        if( !dilepton_pre_vtx_selection_(lepton_pair) ) continue; // before making the SV, cut on the info we have
        //std::cout << "PASS!" << std::endl;

        // check if we have the SV already
        if( !lepton_pair.hasUserFloat("sv_chi2") ) {
          Fitter fitter(
            {*leptons[il1].transient_track, *leptons[il2].transient_track},
            {l1->mass(), l2->mass()},
            {LEP_SIGMA, LEP_SIGMA} //some small sigma for the particle mass
            );
          lepton_pair.addUserFloat("sv_chi2", fitter.chi2());
          lepton_pair.addUserFloat("sv_ndof", fitter.dof()); // float??
          lepton_pair.addUserFloat("sv_prob", fitter.prob());
          // if needed, add here more stuff
        }

        // cut on the SV info
        if( !dilepton_post_vtx_selection_(lepton_pair) ) continue;

        // make the info object to append to the vector
        DiLeptonInfo info;
        info.l1 = &leptons[il1];
        info.l2 = &leptons[il2];
        info.dilepton = &lepton_pair;
        ret_value->push_back(info);
      
      } // loop over l2
    } // loop over l1

    return std::move(ret_value);
  }
};

#endif
