#ifndef LOWPTELOBJECTS_H
#define LOWPTELOBJECTS_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "TLorentzVector.h"
#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "NtupleContent.h"


class LowPtElObjects{
 public:
   LowPtElObjects(edm::EDGetToken & electrons_, edm::EDGetTokenT<edm::ValueMap<float>> & mva_unbiased_ , const edm::Event& iEvent );
   LowPtElObjects(edm::EDGetToken & electrons_,edm::EDGetToken & gsfTracks_, edm::EDGetTokenT<edm::ValueMap<float>> & mva_unbiased_ , const edm::Event& iEvent );
  ~LowPtElObjects();
   void AddElectrons(NtupleContent & nt);
   void AddGsfTracks(NtupleContent & nt);
  

 private:
   edm::Handle<std::vector<pat::Electron>> electrons;  
   edm::Handle<std::vector<reco::GsfTrack>> gsfTracks;   
   edm::Handle<edm::ValueMap<float>> mva_unbiased;

};
#endif
