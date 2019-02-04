#include "TripleTrackKinFit.h"
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


//TripleTrackKinFit::TripleTrackKinFit():
//  m_success(false) {}


TripleTrackKinFit::TripleTrackKinFit(std::vector<RefCountedKinematicParticle>& allParticles){ 
 // std::cout<<"part "<<allParticles.size()<<std::endl;
  KinematicConstrainedVertexFitter kcvFitter;    
  bsTree = kcvFitter.fit(allParticles);
  m_npart=allParticles.size();

}

TripleTrackKinFit::TripleTrackKinFit(std::vector<RefCountedKinematicParticle>& allParticles,ParticleMass Kstar_m){
    KinematicConstrainedVertexFitter kcvFitter;
    MultiTrackKinematicConstraint * Kstar_c = new TwoTrackMassKinematicConstraint(Kstar_m);
    bsTree = kcvFitter.fit(allParticles,Kstar_c);
    m_npart=allParticles.size();
   }


bool TripleTrackKinFit::success(){
  if (bsTree->isEmpty()) { m_success=false; return false;}
  if(!bsTree->isValid()) { m_success=false; return false;}
  if (!bsTree->isConsistent()) {m_success=false;  return false;}
  bsTree->movePointerToTheTop(); b_s=bsTree->currentParticle();
  b_dec_vertex=bsTree->currentDecayVertex();
  if (!b_s->currentState().isValid()){ m_success=false; return false;}
  bs_state=b_s->currentState();
  //bs_state=b_s->initialState();
  if(!b_dec_vertex->vertexIsValid()){m_success=false; return false;}
  bs_children = bsTree->finalStateParticles();
  if(bs_children.size()!=m_npart){ m_success=false; return false;}
  bs_track=b_s->refittedTransientTrack();
  m_success=true; return true;  
}

GlobalVector TripleTrackKinFit::Daughter_Momentum(unsigned int idaughter,bool refit){
  KinematicState child_state;
  if(refit) child_state=bs_children.at(idaughter)->currentState();
  else child_state=bs_children.at(idaughter)->initialState();
  return child_state.globalMomentum();
}

ParticleMass TripleTrackKinFit::Daughter_Mass(unsigned int idaughter,bool refit){
  KinematicState child_state;
  if(refit) child_state=bs_children.at(idaughter)->currentState();
  else child_state=bs_children.at(idaughter)->initialState();
  return child_state.mass();
}

float TripleTrackKinFit::Daughter_Charge(unsigned int idaughter,bool refit){
  KinematicState child_state;
  if(refit) child_state=bs_children.at(idaughter)->currentState();
  else child_state=bs_children.at(idaughter)->initialState();
  return child_state.particleCharge();
}

GlobalVector TripleTrackKinFit::UnfittedMotherMomentum(){
  TLorentzVector vMother; TLorentzVector v1;
  for (unsigned int i=0; i<m_npart; i++){
     KinematicState child_state=bs_children.at(i)->initialState();
     v1.SetPtEtaPhiM( child_state.globalMomentum().perp(), child_state.globalMomentum().eta(), child_state.globalMomentum().phi(), child_state.mass());
    vMother+=v1;
  }
  return GlobalVector(vMother.Px(),vMother.Py(),vMother.Pz());
}

ParticleMass TripleTrackKinFit::UnfittedMotherMass(){
  TLorentzVector vMother; TLorentzVector v1;
  for (unsigned int i=0; i<m_npart; i++){
    KinematicState child_state=bs_children.at(i)->initialState();
    v1.SetPtEtaPhiM( child_state.globalMomentum().perp(), child_state.globalMomentum().eta(), child_state.globalMomentum().phi(), child_state.mass());
    vMother+=v1;
  }
  return vMother.M();
}

TripleTrackKinFit::~TripleTrackKinFit() {}
