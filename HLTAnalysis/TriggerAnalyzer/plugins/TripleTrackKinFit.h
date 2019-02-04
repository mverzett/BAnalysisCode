#ifndef TRIPLETRACKKINFIT_H
#define TRIPLETRACKKINFIT_H



#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/CombinedKinematicConstraint.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackPointingKinematicConstraint.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <vector>

class TripleTrackKinFit{
public: 
//  TripleTrackKinFit();
  TripleTrackKinFit(std::vector<RefCountedKinematicParticle>& allParticles);
  TripleTrackKinFit(std::vector<RefCountedKinematicParticle>& allParticles,ParticleMass Kstar_m);    
  virtual ~TripleTrackKinFit();
  bool success();

  GlobalVector Mother_Momentum(bool refit) {
    if (refit)
      return bs_state.globalMomentum();
    else
      return UnfittedMotherMomentum();    
        }  
  ParticleMass Mother_Mass(bool refit) {
    if (refit)
      return bs_state.mass();
    else
      return UnfittedMotherMass();
}
  float Mother_Charge() {return bs_state.particleCharge();}
  float chi() {return b_dec_vertex->chiSquared();}
  float dof() {return b_dec_vertex->degreesOfFreedom();}
  GlobalVector Daughter_Momentum(unsigned int idaughter,bool refit);
  ParticleMass Daughter_Mass(unsigned int idaughter,bool refit);
  GlobalVector UnfittedMotherMomentum();
  ParticleMass UnfittedMotherMass();
  float Daughter_Charge(unsigned int idaughter,bool refit);
  GlobalPoint Mother_XYZ(){ return b_dec_vertex->position(); }
  GlobalError Mother_XYZError(){ return b_dec_vertex->error(); }
  double Mother_PtError(){ return bs_track.track().ptError(); }
  double Mother_EtaError(){ return bs_track.track().etaError(); }
  double Mother_PhiError(){ return bs_track.track().phiError(); }


private:
  bool m_success;
  RefCountedKinematicTree bsTree;
  RefCountedKinematicVertex b_dec_vertex;
  KinematicState bs_state;
  RefCountedKinematicParticle b_s;
  std::vector< RefCountedKinematicParticle > bs_children;
  reco::TransientTrack bs_track;
  unsigned int m_npart;
};
#endif
