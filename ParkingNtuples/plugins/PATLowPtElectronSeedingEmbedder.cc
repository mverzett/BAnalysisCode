// Unfortunately the (un)biased seeding BDT output is stored in a ValueMap
// with the GSFTrackRef as key, therefore we need a dedicated plugin 
// to embed that into the pat::Electron

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/PatCandidates/interface/Electron.h"

class PATLowPtElectronSeedingEmbedder : public edm::global::EDProducer<> {
public:
  explicit PATLowPtElectronSeedingEmbedder(const edm::ParameterSet &cfg):
    src_{ consumes<pat::ElectronCollection>( cfg.getParameter<edm::InputTag>("src") )},
    biased_seed_{ consumes<edm::ValueMap<float> >( cfg.getParameter<edm::InputTag>("biasedSeeding") )},
    unbiased_seed_{ consumes<edm::ValueMap<float> >( cfg.getParameter<edm::InputTag>("unbiasedSeeding") )} {
      produces<pat::ElectronCollection>();
  }

  ~PATLowPtElectronSeedingEmbedder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const edm::EDGetTokenT<pat::ElectronCollection> src_;
  const edm::EDGetTokenT<edm::ValueMap<float> > biased_seed_;
  const edm::EDGetTokenT<edm::ValueMap<float> > unbiased_seed_; 
};

void PATLowPtElectronSeedingEmbedder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  edm::Handle<pat::ElectronCollection> electrons;
  evt.getByToken(src_, electrons);

  edm::Handle<edm::ValueMap<float> > biased;
  evt.getByToken(biased_seed_, biased);
  
  edm::Handle<edm::ValueMap<float> > unbiased;
  evt.getByToken(unbiased_seed_, unbiased);
  
  std::unique_ptr<pat::ElectronCollection> out(new pat::ElectronCollection(*electrons));

  for(auto &ele : *out) {
    auto available_ids = ele.electronIDs();
    available_ids.push_back(
      std::make_pair(
        "biased_seed",
        (*biased)[ele.gsfTrack()]
        )
      );
    available_ids.push_back(
      std::make_pair(
        "unbiased_seed", 
        (*unbiased)[ele.gsfTrack()]
        )
      );
    
    //for some reason this overwrites the existing IDs, hence all the above trick
    ele.setElectronIDs(available_ids);
  }

  evt.put(std::move(out));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATLowPtElectronSeedingEmbedder);
