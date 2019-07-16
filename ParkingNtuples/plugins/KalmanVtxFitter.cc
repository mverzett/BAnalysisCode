#include "KalmanVtxFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"


KalmanVtxFitter::KalmanVtxFitter(const std::vector<reco::TransientTrack> tracks, 
                           const std::vector<double> masses, 
                           const std::vector<double> sigmas):
  n_particles_{masses.size()} {
  
  KalmanVertexFitter vtxFitter(true);
  fitted_vtx_ = vtxFitter.vertex(tracks); 
  if(!fitted_vtx_.isValid()) {
    success_ = false;
    return;
  }
  success_ = true;
  if(fitted_vtx_.hasRefittedTracks()) {
    refitted_tracks_ = fitted_vtx_.refittedTracks();
  }
}
