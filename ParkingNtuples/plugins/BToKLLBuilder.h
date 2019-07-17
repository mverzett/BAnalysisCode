/*

  Class used to build and select B->Kll candidates

  original author: Mauro Verzetti (CERN)

 */

#ifndef BAnalysisCode_ParkingNtuples_BToKLLBuilder
#define BAnalysisCode_ParkingNtuples_BToKLLBuilder

#include "DiLeptonBuilder.h"
#include <vector>
#include <memory>
#include <map>
#include <string>
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "helpers.h"
#include <algorithm>

template<typename Lepton, typename Fitter>
class BToKLLBuilder {
public:
  typedef std::vector< ChachedObject<Lepton> > LeptonCollection;
  typedef std::map< std::pair<size_t, size_t>, pat::CompositeCandidate > DiLeptonCache;

  BToKLLBuilder(const edm::ParameterSet&);
  ~BToKLLBuilder() {}

  // Nothing is const here as we modify everything
  std::unique_ptr< pat::CompositeCandidateCollection > build(LeptonCollection&, CachedCandidateCollection&, DiLeptonCache &) const;

private:
  const bool active_ = false;
  const DiLeptonBuilder<Lepton, Fitter> ll_builder_;
  const StringCutObjectSelector<pat::PackedCandidate> k_selection_; // cut on sub-leading lepton
  const StringCutObjectSelector<pat::CompositeCandidate> candidate_pre_vtx_selection_; // cut on the candidate before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> candidate_post_vtx_selection_; // cut on the candidate after the SV fit
};

template<typename Lepton, typename Fitter>
BToKLLBuilder<Lepton, Fitter>::BToKLLBuilder(const edm::ParameterSet& cfg):
  active_{cfg.getParameter<bool>("active")},
  ll_builder_{cfg},
  k_selection_{cfg.getParameter<std::string>("kSelection")},
  candidate_pre_vtx_selection_{cfg.getParameter<std::string>("candidatePreVtxSelection")},
  candidate_post_vtx_selection_{cfg.getParameter<std::string>("candidatePostVtxSelection")} {}

template<typename Lepton, typename Fitter>
std::unique_ptr< pat::CompositeCandidateCollection >
BToKLLBuilder<Lepton, Fitter>::build(LeptonCollection& leptons, CachedCandidateCollection& tracks, DiLeptonCache &cache) const {
  auto ret_val = std::make_unique<pat::CompositeCandidateCollection>();
  if(!active_) return std::move(ret_val);
  
  // get dilepton pairs
  auto lepton_pairs = ll_builder_.build(leptons, cache);

  // get number of stored electrons and tracks
  int nleps = -1;
  for(const auto& l : leptons) std::max(nleps, l.idx);
  
  int ntrks = -1;
  for(const auto& t : tracks) std::max(ntrks, t.idx);

  for(auto &lepton_pair : *lepton_pairs) {
    for(auto &track : tracks) {
      if( !k_selection_(*track.obj) ) continue;
      
      pat::CompositeCandidate cand;
      cand.addDaughter( *lepton_pair.l1->obj );
      cand.addDaughter( *lepton_pair.l2->obj );
      cand.addDaughter( *track.obj );

      auto dr_info = min_max_dr(cand);
      cand.addUserFloat("min_dr", dr_info.first);
      cand.addUserFloat("max_dr", dr_info.second);
      
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
        {*lepton_pair.l1->transient_track, *lepton_pair.l2->transient_track, *track.transient_track},
        {lepton_pair.l1->obj->mass(), lepton_pair.l2->obj->mass(), K_MASS},
        {LEP_SIGMA, LEP_SIGMA, K_SIGMA} //some small sigma for the lepton mass
        );
      cand.addUserFloat("sv_chi2", fitter.chi2());
      cand.addUserFloat("sv_ndof", fitter.dof()); // float??
      cand.addUserFloat("sv_prob", fitter.prob());

      if( !candidate_post_vtx_selection_(cand) ) continue;

      // store the candidate and mark the indices
      int l1_idx = (lepton_pair.l1->idx == -1) ? nleps++ : lepton_pair.l1->idx;
      int l2_idx = (lepton_pair.l2->idx == -1) ? nleps++ : lepton_pair.l2->idx;
      int trk_idx = (track.idx == -1) ? ntrks++ : track.idx;
      lepton_pair.l1->idx = l1_idx;
      lepton_pair.l2->idx = l2_idx;
      track.idx = trk_idx;
      cand.addUserInt("l1_idx", l1_idx);
      cand.addUserInt("l2_idx", l2_idx);
      cand.addUserInt("k_idx", trk_idx);

      ret_val->push_back(cand);
    }
  }

  return std::move(ret_val);
}

#endif  
