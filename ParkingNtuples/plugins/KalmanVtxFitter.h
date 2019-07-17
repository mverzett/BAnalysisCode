#ifndef BAnalysisCode_ParkingNtuples_KalmanVtxFitter
#define BAnalysisCode_ParkingNtuples_KalmanVtxFitter

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include <vector>

class KalmanVtxFitter {
public: 
  KalmanVtxFitter():
    fitted_vtx_{}, 
    refitted_tracks_{} {};

  KalmanVtxFitter(const std::vector<reco::TransientTrack> tracks, 
               const std::vector<double> masses, 
               const std::vector<double> sigmas);

  ~KalmanVtxFitter() {};

  bool success() const {return success_;}
  float chi2() const {return success_ ? fitted_vtx_.totalChiSquared()  : 999;}
  float dof() const  {return success_ ? fitted_vtx_.degreesOfFreedom() : -1;}
  float prob() const {
    return success_ ? ChiSquaredProbability(chi2(), dof()) : 0.;
  }
  
  GlobalPoint fitted_vtx() const {
    return fitted_vtx_.position();
  }
  GlobalError fitted_vtx_uncertainty() const {
    return fitted_vtx_.positionError();
  }

private:
  float kin_chi2_ = 0.;
  float kin_ndof_ = 0.;
  size_t n_particles_ = 0;
  bool success_ = false;

  TransientVertex fitted_vtx_;
  std::vector<reco::TransientTrack> refitted_tracks_;
};
#endif
