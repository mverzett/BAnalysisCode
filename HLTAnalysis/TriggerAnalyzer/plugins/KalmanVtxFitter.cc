#include "KalmanVtxFitter.h"


KalmanVtxFitter::KalmanVtxFitter(std::vector<reco::TransientTrack> & ttks,bool refit,std::vector<float> &mass){
  mass_=mass;
  KalmanVertexFitter vtxFitter(true);
  myVertex=vtxFitter.vertex(ttks); 
  if (!myVertex.isValid() ) {m_success=false; return;}
  m_success=true;
  if (myVertex.hasRefittedTracks() && refit)
     KFrefitted=myVertex.refittedTracks();
  else
    KFrefitted=ttks;  
  for (reco::TransientTrack & trk: KFrefitted){
    tch+=trk.charge(); TLorentzVector vtemp;
    vtemp.SetPtEtaPhiM(trk.track().pt(),trk.track().eta(),trk.track().phi(),mass[&trk-&KFrefitted[0]]);
    vB+=vtemp;
  }
 }

KalmanVtxFitter::~KalmanVtxFitter() {}

GlobalVector KalmanVtxFitter::Daughter_Momentum(unsigned int idaughter){
  return GlobalVector(KFrefitted[idaughter].track().px(),KFrefitted[idaughter].track().py(),KFrefitted[idaughter].track().pz());
}

float KalmanVtxFitter::Daughter_Charge(unsigned int idaughter){
  return KFrefitted[idaughter].charge();
}

ParticleMass KalmanVtxFitter::Daughter_Mass(unsigned int idaughter){
  return mass_[idaughter];
}




