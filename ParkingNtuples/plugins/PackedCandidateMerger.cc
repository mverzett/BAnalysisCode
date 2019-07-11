// Merges the PF and LowPT collections, sets the isPF and isLowPt 
// UserInt's accordingly

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include <limits>
#include <algorithm>

class PackedCandidateMerger : public edm::global::EDProducer<> {
public:
  explicit PackedCandidateMerger(const edm::ParameterSet &cfg):
    packed_{ consumes<pat::PackedCandidateCollection>( cfg.getParameter<edm::InputTag>("packedSrc") )},
    lost_{ consumes<pat::PackedCandidateCollection>( cfg.getParameter<edm::InputTag>("lostSrc") )} {
      produces<pat::PackedCandidateCollection>();
  }

  ~PackedCandidateMerger() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const edm::EDGetTokenT<pat::PackedCandidateCollection> packed_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> lost_;
};

void PackedCandidateMerger::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  edm::Handle<pat::PackedCandidateCollection> packed;
  evt.getByToken(packed_, packed);

  edm::Handle<pat::PackedCandidateCollection> lost;
  evt.getByToken(lost_, lost);
  
  std::unique_ptr<pat::PackedCandidateCollection> out(new pat::PackedCandidateCollection());

  std::copy(packed->begin(), packed->end(), std::back_inserter(*out));
  std::copy(lost->begin(), lost->end(), std::back_inserter(*out));

  evt.put(std::move(out));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PackedCandidateMerger);
