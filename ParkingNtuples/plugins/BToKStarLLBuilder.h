/*

  Class used to build and select B->Kll candidates

  original author: Mauro Verzetti (CERN)

 */

#ifndef BAnalysisCode_ParkingNtuples_BToKStarLLBuilder
#define BAnalysisCode_ParkingNtuples_BToKStarLLBuilder

#include "DiLeptonBuilder.h"
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
#include <algorithm>

template<typename Lepton, typename Fitter>
class BToKStarLLBuilder {
public:
  typedef std::vector<Lepton> LeptonCollection;
  typedef std::map< std::pair<size_t, size_t>, pat::CompositeCandidate > DiLeptonCache;

  BToKStarLLBuilder(const edm::ParameterSet&);
  ~BToKStarLLBuilder() {}

  // Nothing is const here as we modify everything
  std::unique_ptr< pat::CompositeCandidateCollection > build(LeptonCollection&, CahcedTrackCollection&, DiLeptonCache &) const;

private:
  const DiLeptonBuilder<Lepton, Fitter> ll_builder_;
  const StringCutObjectSelector<reco::Track> k_selection_; // cut on sub-leading lepton
  const StringCutObjectSelector<reco::Track> pi_selection_; // cut on sub-leading lepton
  const StringCutObjectSelector<pat::CompositeCandidate> candidate_pre_vtx_selection_; // cut on the candidate before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> candidate_post_vtx_selection_; // cut on the candidate after the SV fit
};

template<typename Lepton, typename Fitter>
BToKStarLLBuilder::BToKStarLLBuilder(const edm::ParameterSet& cfg):
  ll_builder_{cfg},
  k_selection_{cfg.getParameter<std::string>("kSelection")},
  k_selection_{cfg.getParameter<std::string>("piSelection")},
  candidate_pre_vtx_selection_{cfg.getParameter<std::string>("candidatePreVtxSelection")},
  candidate_post_vtx_selection_{cfg.getParameter<std::string>("candidatePostVtxSelection")} {}

template<typename Lepton, typename Fitter>
std::unique_ptr< pat::CompositeCandidateCollection >
BToKStarLLBuilder::build(LeptonCollection& leptons, CahcedTrackCollection& tracks, DiLeptonCache &cache) const {
  auto ret_val = std::make_unique<pat::CompositeCandidateCollection>();
  // get dilepton pairs
  auto lepton_pairs = ll_builder_.build(leptons, cache);

  // get number of stored electrons and tracks
  int nleps = -1;
  for(const auto& l : leptons) std::max(nleps, l.idx);
  
  int ntrks = -1;
  for(const auto& t : tracks) std::max(ntrks, t.idx);

  for(auto &lepton_pair : lepton_pairs) {
    for(auto &k_track : tracks) {
      if( !k_selection_(*k_track.obj) ) continue;

      for(auto &pi_track : tracks) {
        if(&pi_track == &k_track) continue;
        if( !pi_selection_(*pi_track.obj) ) continue;
      
        pat::CompositeCandidate cand;
        cand.addDaughter( *lepton_pair.l1->obj );
        cand.addDaughter( *lepton_pair.l2->obj );
        cand.addDaughter( *k_track.obj );
        cand.addDaughter( *pi_track.obj );
      
        // propagate UserFloats from the dilepton pair
        for(const auto& float_name : lepton_pair.dilepton->userFloatNames()) {
          cand.addUserFloat(
            "ll_"+float_name, 
            lepton_pair.dilepton->userFloat(float_name)
            );
        }
        // TODO add meaningful variables
      
        if( !candidate_pre_vtx_selection_(cand) ) continue;

        Fitter fitter(
          {*lepton_pair.l1->transient_track, *lepton_pair.l1->transient_track, k_track.transient_track, pi_track.transient_track},
          {lepton_pair.l1->obj->mass(), lepton_pair.l2->obj->mass(), K_MASS, PI_MASS},
          {LEP_SIGMA, LEP_SIGMA, K_SIGMA, PI_SIGMA}, //some small sigma for the lepton mass
          );
        cand.addUserFloat("sv_chi2", fitter.chi());
        cand.addUserFloat("sv_ndof", fitter.dof());
        cand.addUserFloat("sv_prob", ChiSquaredProbability(
                          fitter.chi(), fitter.dof()) );

        if( !candidate_post_vtx_selection_(cand) ) continue;

        // store the candidate and mark the indices
        int l1_idx = (lepton_pair.l1->idx == -1) ? nleps++ : lepton_pair.l1->idx;
        int l2_idx = (lepton_pair.l2->idx == -1) ? nleps++ : lepton_pair.l2->idx;
        int k_idx = (k_track.idx == -1) ? ntrks++ : k_track.idx;
        int pi_idx = (pi_track.idx == -1) ? ntrks++ : pi_track.idx;
        lepton_pair.l1->idx = l1_idx;
        lepton_pair.l2->idx = l2_idx;
        k_track.idx = trk_idx;
        pi_track.idx = pi_idx;
        cand.addUserInt("l1_idx", l1_idx);
        cand.addUserInt("l2_idx", l2_idx);
        cand.addUserInt("k_idx", trk_idx);
        cand.addUserInt("pi_idx", trk_idx);        

        ret_val->push_back(cand);
      } // loop over pi
    } // loop over k
  } // loop over lepton pairs

  return std::move(ret_val);
}

#endif  
