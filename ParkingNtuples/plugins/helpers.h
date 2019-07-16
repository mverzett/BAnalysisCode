/*
  Simple set of constants, functions and classes that do not belong anywhere specifically

  original author: Mauro Verzetti (CERN)
 */

#ifndef BAnalysisCode_ParkingNtuples_helpers
#define BAnalysisCode_ParkingNtuples_helpers

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h" //just a typedef to double -.-'
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include <vector>
#include <algorithm>
#include <limits>

template<typename T>
struct ChachedObject {
  // Convenience struct to store the object, its transient track, and the
  // index the object will have in the final collection (-1 if not assigned)
  // all pointers are non-owning, hence just pointers
  ChachedObject() {};
  ChachedObject(T* i, reco::TransientTrack* t, int index = -1):
    transient_track{t},
    obj{i},
    idx{index} {};
  reco::TransientTrack* transient_track = 0; 
  T* obj = 0;
  int idx = -1;
};

typedef ChachedObject<pat::Muon> CachedMuon;
typedef ChachedObject<pat::Electron> CachedElectron;
typedef ChachedObject<pat::PackedCandidate> CachedCandidate;

typedef std::vector<CachedMuon> CachedMuonCollection;
typedef std::vector<CachedElectron> CachedElectronCollection;
typedef std::vector<CachedCandidate> CachedCandidateCollection;

constexpr double K_MASS = 0.493677;
constexpr double PI_MASS = 0.139571;
constexpr float LEP_SIGMA = 0.0000001;
constexpr float K_SIGMA = 0.000016;
constexpr float PI_SIGMA = 0.000016;

inline std::pair<float, float> min_max_dr(const pat::CompositeCandidate & cand) {
  float min_dr = std::numeric_limits<float>::max();
  float max_dr = 0.;
  for(size_t i = 0; i < cand.numberOfDaughters(); ++i) {
    for(size_t j = i+1; j < cand.numberOfDaughters(); ++j) {
      float dr = reco::deltaR(*cand.daughter(i), *cand.daughter(j));
      min_dr = std::min(min_dr, dr);
      max_dr = std::max(max_dr, dr);
    }
  }
  return make_pair(min_dr, max_dr);
}

#endif
