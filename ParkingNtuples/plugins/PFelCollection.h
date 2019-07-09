#ifndef PFELCOLLECTION_H
#define PFELCOLLECTION_H

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


class PFelCollection{
 public:
   PFelCollection(edm::EDGetToken & electrons_,edm::EDGetTokenT<edm::ValueMap<bool>> & soft_id_,edm::EDGetTokenT<edm::ValueMap<bool>> & medium_id_ ,edm::EDGetTokenT<edm::ValueMap<bool>> & tight_id_, const edm::Event& iEvent );
  ~PFelCollection();
   void AddElectrons(NtupleContent & nt);
  
  

 private:
   edm::Handle<std::vector<pat::Electron>> electrons;
   edm::Handle<edm::ValueMap<bool>> soft_id;
   edm::Handle<edm::ValueMap<bool>> medium_id;
   edm::Handle<edm::ValueMap<bool>> tight_id;


};
#endif
