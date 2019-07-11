/*
  Simple set of constants, functions and classes that do not belong anywhere specifically

  original author: Mauro Verzetti (CERN)
 */

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h" //just a typedef to double -.-'
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

template<typename T>
struct ChachedObject {
  // Convenience struct to store the object, its transient track, and the
  // index the object will have in the final collection (-1 if not assigned)
  // all pointers are non-owning, hence just pointers
  TransientTrack* transient_track = 0; 
  T* obj = 0;
  int idx = -1;
};
typedef ChachedObject<pat::Muon> CachedMuon;
typedef ChachedObject<pat::Electron> CachedElectron;
typedef ChachedObject<reco::Track> CachedTrack;

constexpr double LEP_SIGMA = 0.0000001;
constexpr double K_MASS = 0.493677;
constexpr double K_SIGMA = 0.000016;
constexpr double PI_MASS = 0.139571;
constexpr double PI_SIGMA = 0.000016;


