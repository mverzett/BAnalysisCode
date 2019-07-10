// Merges the PF and LowPT collections, sets the isPF and isLowPt 
// UserInt's accordingly

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include <limits>
#include <algorithm>

class ElectronMerger : public edm::global::EDProducer<> {
public:
  explicit ElectronMerger(const edm::ParameterSet &cfg):
    lowpt_src_{ consumes<pat::ElectronCollection>( cfg.getParameter<edm::InputTag>("lowptSrc") )},
    pf_src_{ consumes<pat::ElectronCollection>( cfg.getParameter<edm::InputTag>("pfSrc") )},
    dr_cleaning_{cfg.getParameter<double>("drForCleaning")},
    dz_cleaning_{cfg.getParameter<double>("dzForCleaning")} {
      produces<pat::ElectronCollection>();
  }

  ~ElectronMerger() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const edm::EDGetTokenT<pat::ElectronCollection> lowpt_src_;
  const edm::EDGetTokenT<pat::ElectronCollection> pf_src_;
  const double dr_cleaning_;
  const double dz_cleaning_;
};

void ElectronMerger::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  edm::Handle<pat::ElectronCollection> lowpt;
  evt.getByToken(lowpt_src_, lowpt);

  edm::Handle<pat::ElectronCollection> pf;
  evt.getByToken(pf_src_, pf);
  
  std::unique_ptr<pat::ElectronCollection> out(new pat::ElectronCollection);

  for(auto ele : *pf) {
    ele.addUserInt("isPF", 1);
    ele.addUserInt("isLowPt", 0);
    out->push_back(ele);
  }

  for(auto ele : *lowpt) {
    double min_dr = std::numeric_limits<double>::max();
    for(const auto& pfele : *pf) {
      if(fabs(pfele.vz() - ele.vz()) > dz_cleaning_) continue;
      min_dr = std::min(min_dr, reco::deltaR(ele, pfele));
    }
    if(min_dr < dr_cleaning_) continue;
    ele.addUserInt("isPF", 0);
    ele.addUserInt("isLowPt", 1);
    out->push_back(ele);
  }

  evt.put(std::move(out));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ElectronMerger);
