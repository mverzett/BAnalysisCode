#ifndef BKLLDECAY_H
#define BKLLDECAY_H

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/CombinedKinematicConstraint.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "NtupleContent.h"
#include "TripleTrackKinFit.h"
#include <vector>
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"



class BKlldecay{
   
public:
  BKlldecay(unsigned int &nmupairs_, std::vector<std::pair<unsigned int,unsigned int>>& used_lep_tracks, std::vector<reco::TransientTrack> & vmuTrack1,std::vector<reco::TransientTrack> & vmuTrack2, std::vector<reco::TransientTrack> & vKTrack, float & beam_xz_,float & beam_yz_ ,bool RefitMuTracksOnly_);
  virtual ~BKlldecay();
  void Fill(NtupleContent & nt);
  void ProbCut(float & Pchi2BMuMuK) { probcut_=Pchi2BMuMuK; }
  void CosCut(float & CosThetaCut) { coscut_=CosThetaCut; }
  void MassCuts(float &min, float &max) { massmin_=min; massmax_=max;  }

private:
  ParticleMass part_mass=0.1056583; float part_sigma=0.0000001;
  ParticleMass kaon_mass=0.493677; float kaon_sigma=0.000016;
  float chi=0,ndf=0; float probcut_=-1,coscut_=-2,massmin_=0,massmax_=100; 
  math::XYZVector vperp; math::XYZVector pperp; std::vector<float> temp;  
  unsigned int &nmupairs;
  std::vector<std::pair<unsigned int,unsigned int>>& used_muTracks_index;
  std::vector<reco::TransientTrack> & muTrack1; std::vector<reco::TransientTrack> & muTrack2; 
  std::vector<reco::TransientTrack> & KTrack; float & beam_xz; float & beam_yz; 
  bool &  RefitMuTracksOnly;
  std::vector<float> tempK,tempBpt,tempBx,tempBex,tempBept;
  //reco::TrackBase::Point vertex_point; 
  //  math::XYZPoint vertex_point;
};

#endif 

