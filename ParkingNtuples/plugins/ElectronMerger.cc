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
    dz_cleaning_{cfg.getParameter<double>("dzForCleaning")},
    use_gsf_mode_for_p4_{cfg.getParameter<bool>("useGsfModeForP4")} {
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
  const bool use_gsf_mode_for_p4_;
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
    bool clean_out = false;
    for(const auto& pfele : *pf) {
      clean_out |= (
        fabs(pfele.vz() - ele.vz()) < dz_cleaning_ &&
        reco::deltaR(ele, pfele) < dr_cleaning_
        );
    }
    if(clean_out) continue;
    ele.addUserInt("isPF", 0);
    ele.addUserInt("isLowPt", 1);
    if(use_gsf_mode_for_p4_) {
      reco::Candidate::PolarLorentzVector p4(
        ele.gsfTrack()->ptMode(),
        ele.gsfTrack()->etaMode(),
        ele.gsfTrack()->phiMode(),
        ele.mass()
        );
      ele.setP4(p4);
    }
    out->push_back(ele);
  }

  std::sort(
    out->begin(), out->end(), 
    [] (pat::Electron e1, pat::Electron e2) -> bool {return e1.pt() > e2.pt();}
    );
  evt.put(std::move(out));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ElectronMerger);
