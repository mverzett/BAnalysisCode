#ifndef KALMANVTXFITTER_H
#define KALMANVTXFITTER_H

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TLorentzVector.h"
#include <vector>

class KalmanVtxFitter{
public: 

  KalmanVtxFitter(std::vector<reco::TransientTrack> & ttks, bool refit,std::vector<float> &mass);   
  virtual ~KalmanVtxFitter();

  bool success() {return m_success;}
  GlobalVector Mother_Momentum() {
     return GlobalVector(vB.Px(),vB.Py(),vB.Pz()); }  
  ParticleMass Mother_Mass() { return  vB.M();}
  float Mother_Charge() { return tch; }
  float chi() {return myVertex.totalChiSquared();}
  float dof() {return myVertex.degreesOfFreedom();}
  GlobalVector Daughter_Momentum(unsigned int idaughter);
  ParticleMass Daughter_Mass(unsigned int idaughter);
  float Daughter_Charge(unsigned int idaughter);
  GlobalPoint Mother_XYZ(){ return myVertex.position(); }
  GlobalError Mother_XYZError(){ return myVertex.positionError(); }
  double Mother_PtError(){ return -999; }
  double Mother_EtaError(){ return -999; }
  double Mother_PhiError(){ return -999; }
  //KVF does not compute errors on refited Pt


private:
  bool m_success; TransientVertex myVertex; float tch=0;
  std::vector<reco::TransientTrack> KFrefitted; TLorentzVector vB;
  std::vector<float> mass_;
};
#endif
