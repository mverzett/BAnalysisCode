// Package:    HLTAnalysis/TriggerAnalyzer
// Class:      TriggerAnalyzer
// 
/**\class TriggerAnalyzer TriggerAnalyzer.cc HLTAnalysis/TriggerAnalyzer/plugins/TriggerAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:
//                George Karathanasis (georgios.karathanasis@cern.ch)
//         Created:  Thu, 23 Mar 2017 17:40:23 GMT
//
//


// system include files
#include <memory>
#include <string>
#include <iostream>
#include <vector>
// user include files
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDFilter.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1TUtmTriggerMenu.h"

#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "TTree.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "TripleTrackKinFit.h"
#include "GeneratorBTree.h"
#include "NtupleContent.h"
#include "LowPtElObjects.h"
#include "PFelCollection.h"
#include "HLTL1tree.h"
#include "BKlldecay.h"
#include "BKstarlldecay.h"

#include "helpers.h"
#include "DiLeptonBuilder.h"
#include "BToKLLBuilder.h"
#include "BToKStarLLBuilder.h"
#include "KinVtxFitter.h"

using namespace std;


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.
class ParkingNtupleMaker : public edm::one::EDFilter<edm::one::WatchRuns>  {

public:
  explicit ParkingNtupleMaker(const edm::ParameterSet&);
  ~ParkingNtupleMaker();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  void         beginRun(edm::Run const& iEvent, edm::EventSetup const&) ;
  virtual bool filter(edm::Event&, edm::EventSetup const&) override;
  void         endRun(edm::Run const& iEvent, edm::EventSetup const&){};
  virtual void endJob() override;
//   virtual void beginRun(const edm::Run &, const edm::EventSetup &);
  inline void fill_empty(edm::Event& evt) {
    // utility function to fill all empty collections
    auto electrons_out = make_unique<pat::ElectronCollection>();
    auto muons_out = make_unique<pat::MuonCollection>();
    auto cands_out = make_unique<pat::PackedCandidateCollection>();
    auto b_kmumu = make_unique<pat::CompositeCandidateCollection>();
    auto b_kee = make_unique<pat::CompositeCandidateCollection>();

    evt.put(std::move(electrons_out), "electrons");
    evt.put(std::move(muons_out), "muons");
    evt.put(std::move(cands_out), "tracks");
    evt.put(std::move(b_kmumu), "BToKMuMu");
    evt.put(std::move(b_kee), "BToKEE");
  
  }

  const edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  const edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken_;
  const edm::EDGetTokenT<pat::ElectronCollection> electronsToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> tracksToken_;
  const edm::EDGetTokenT<pat::MuonCollection> muonsToken_;
//  edm::EDGetToken photonToken_;
  const edm::EDGetTokenT<GlobalAlgBlkBxCollection> l1resultToken_;
  const edm::EDGetTokenT<l1t::MuonBxCollection> l1MuonsToken_;
  const vector<string> seed_;
  const edm::EDGetTokenT<edm::TriggerResults> trgresultsToken_; 
  const edm::EDGetTokenT<vector<pat::TriggerObjectStandAlone>> trigobjectsToken_;
  const vector<string> HLTPath_;
  HLTPrescaleProvider hltPrescaleProvider_;
  const edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
  const edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;

  // Candidate builders
  const BToKLLBuilder<pat::Muon, KinVtxFitter> b_to_kmumu_builder_;
  const BToKLLBuilder<pat::Electron, KinVtxFitter> b_to_kee_builder_;
  // const BToKStarLLBuilder<CachedMuon, KinVtxFitter> b_to_kstarmumu_builder_;
  // const BToKStarLLBuilder<CachedElectron, KinVtxFitter> b_to_kstaree_builder_;

  ///options
  const double dr_trig_mu_;
  const double electron_dz_cut_;
  const double mu_trk_min_dr_;

  // bool UseDirectlyGenBeeK=false; 
  // double DRgenCone=10;
  // int KIdToMatch=-1; 
  // int LepIdToMatch=-1; 
  // int BpdgIdToMatch=-1;
  // bool IsResonantDecayToMatch=false; 

  TString * algoBitToName = new TString[512];
  // ----------member data ---------------------------
};


// constructors and destructor
//
ParkingNtupleMaker::ParkingNtupleMaker(const edm::ParameterSet& iConfig): 
  beamSpotToken_{consumes<reco::BeamSpot>(iConfig.getParameter <edm::InputTag>("beamSpot"))},
  vtxToken_{consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))},
  electronsToken_{consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>  ("electrons"))},
  tracksToken_{consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("tracks"))},
  muonsToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  l1resultToken_(consumes<GlobalAlgBlkBxCollection>(iConfig.getParameter<edm::InputTag>("l1seed"))),
  l1MuonsToken_(consumes<l1t::MuonBxCollection>(iConfig.getParameter<edm::InputTag>("l1muons"))),
  seed_(iConfig.getParameter<vector<string> >("seed")),
  trgresultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag> ("triggerresults"))),
  trigobjectsToken_(consumes<vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag> ("triggerobjects"))),
  HLTPath_(iConfig.getParameter<vector<string> >("HLTPath")),
  hltPrescaleProvider_ (iConfig, consumesCollector(), *this),
  prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
  packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
  b_to_kmumu_builder_{iConfig.getParameter<edm::ParameterSet>("BToKMuMu")},
  b_to_kee_builder_{iConfig.getParameter<edm::ParameterSet>("BToKEE")},
  dr_trig_mu_{iConfig.getParameter<double>("DRForTrigMu")},
  electron_dz_cut_{iConfig.getParameter<double>("DzMaxFromTrigMu")},
  mu_trk_min_dr_{iConfig.getParameter<double>("minDRMuToTrk")} {

  produces<pat::ElectronCollection>("electrons");
  produces<pat::MuonCollection>("muons");
  produces<pat::PackedCandidateCollection>("tracks");
  produces<pat::CompositeCandidateCollection>("BToKMuMu");
  produces<pat::CompositeCandidateCollection>("BToKEE");
  
  
}

ParkingNtupleMaker::~ParkingNtupleMaker()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
bool
ParkingNtupleMaker::filter(edm::Event& iEvent, edm::EventSetup const& iSetup)
{
  using namespace std;
  //Get a few collections to apply basic electron ID
  //Get data
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(beamSpotToken_,theBeamSpot); 
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  //continue if there are no vertices
  if (vertices->size()==0) {    
    fill_empty(iEvent);
    return false;
  }
  edm::Handle<std::vector<pat::Electron>> electrons;
  iEvent.getByToken(electronsToken_, electrons); 

  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonsToken_,muons);

  edm::Handle<pat::PackedCandidateCollection> tracks;
  iEvent.getByToken(tracksToken_, tracks);

  edm::ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

  edm::Handle<GlobalAlgBlkBxCollection> l1result;
  iEvent.getByToken(l1resultToken_,l1result);
  if (count==0){
    edm::ESHandle<L1TUtmTriggerMenu> menu;
    iSetup.get<L1TUtmTriggerMenuRcd>().get(menu);
    if (l1result.isValid()) {
      for (auto const & keyval: menu->getAlgorithmMap()) {
        std::string const & trigName  = keyval.second.getName();
        unsigned int index = keyval.second.getIndex();  
        int itrig = index;
        algoBitToName[itrig] = TString( trigName );
      }
      count++;
    }
  }
  
  // per-event quantities
  reco::Vertex& firstGoodVertex = vertices->front();
  for (const reco::Vertex &vtx : *vertices) {
    if ( vtx.isFake() || !vtx.isValid() ) continue;
    firstGoodVertex = vtx;
    break;
  }
  auto vertex_point = firstGoodVertex.position();
  
  // std::pair<float,float> EtaPhiE1(-10,-10),EtaPhiE2(-10,-10),EtaPhiK(-10,-10);
  // if(! iEvent.isRealData() ){
  //   if (UseDirectlyGenBeeK){
  //     gen.GenKLepLepEtaPhi(IsResonantDecayToMatch,BpdgIdToMatch,LepIdToMatch,KIdToMatch);
  //     EtaPhiK=gen.EtaPhiK();
  //     EtaPhiE1=gen.EtaPhiL();
  //     EtaPhiE2=gen.EtaPhiaL();
  //   }
  // }
  // if (UseDirectlyGenBeeK && (EtaPhiK.second==-10 || EtaPhiE1.second==-10 || EtaPhiE2.second==-10) ) {    
  //   fill_empty(iEvent);
  //   return false;
  // }
    
  //save trigger related information  
  std::pair<float,float> TrgMu_EtaPhi(-100,-100);
  HLTL1tree trigger(
    iEvent, l1resultToken_, l1MuonsToken_, 
    trgresultsToken_, trigobjectsToken_);
  trigger.L1trigger(algoBitToName, seed_); 
  trigger.HLTtrigger(HLTPath_, hltPrescaleProvider_); 
  trigger.HLTobjects(HLTPath_);
  if(!trigger.HLTPathFire()) {
    fill_empty(iEvent);
    return false;
  }
  TrgMu_EtaPhi = std::make_pair(
    trigger.GetHighestPtHLTObject()[1],
    trigger.GetHighestPtHLTObject()[2]
    );
  
  // Stores to keep the objects we care about
  // We are not going to access them directly from here often
  std::vector<pat::Muon> muon_store;
  std::vector<pat::Electron> electron_store;
  std::vector<pat::PackedCandidate> candidate_store;
  std::vector<reco::TransientTrack> muon_ttracks;
  std::vector<reco::TransientTrack> electron_ttracks;
  std::vector<reco::TransientTrack> candidate_ttracks;
  
  //  save offline muons (all or saveOnlyHLTFires)
  //std::cout << "Muons" << std::endl;
  bool found_trig_muon = false;
  double trig_mu_z = 0.;
  for (const pat::Muon &mu : *muons){    
    const reco::MuonPFIsolation& isol = mu.pfIsolationR04();
    double mu_iso=(isol.sumChargedHadronPt+max(0.,isol.sumNeutralHadronEt+isol.sumPhotonEt-0.5*isol.sumPUPt))/mu.pt();
    
    // Store info in stores
    muon_store.push_back(mu);    
    auto &new_mu = muon_store.back();

    // Add variables we care about
    bool is_triggering = deltaR(TrgMu_EtaPhi.first, TrgMu_EtaPhi.second, new_mu.eta(), new_mu.phi()) < dr_trig_mu_;
    if(is_triggering) trig_mu_z = new_mu.vz();
    found_trig_muon |= is_triggering;
    new_mu.addUserInt("triggering", is_triggering);
    new_mu.addUserInt("soft_id", new_mu.isSoftMuon(firstGoodVertex));
    new_mu.addUserInt("tight_id", new_mu.isTightMuon(firstGoodVertex));
    new_mu.addUserFloat("dz", new_mu.bestTrack()->dz(vertex_point));
    new_mu.addUserFloat("dxy", new_mu.bestTrack()->dxy(vertex_point));
    new_mu.addUserFloat("mu_iso", mu_iso);

    muon_ttracks.emplace_back(reco::TransientTrack(*mu.bestTrack(), &(*bFieldHandle))); // Better way with the TTrack builder?
  }

  if (found_trig_muon) {
    fill_empty(iEvent);
    return false;
  }
  
  //electrons  
  //std::cout << "Electrons" << std::endl;
  for(const pat::Electron &el : *electrons) {
    if(fabs(trig_mu_z - el.vz()) > electron_dz_cut_) continue;
    double iso = el.pfIsolationVariables().sumChargedHadronPt + max(
      0.0, 
      el.pfIsolationVariables().sumNeutralHadronEt + el.pfIsolationVariables().sumPhotonEt - 0.5 * el.pfIsolationVariables().sumPUPt
      ) / el.pt();

    // Store info in stores
    electron_store.push_back(el);
    electron_ttracks.emplace_back(reco::TransientTrack(*el.bestTrack(), &(*bFieldHandle))); // Better way with the TTrack builder?
    auto &new_el = electron_store.back();

    // Add variables we care about
    new_el.addUserFloat("iso", iso);
    new_el.addUserFloat("dz" , el.bestTrack()->dz(vertex_point) );
    new_el.addUserFloat("dxy", el.bestTrack()->dxy(vertex_point));
  }
  
  //std::cout << "Tracks:" << std::endl;
  for (const pat::PackedCandidate &trk : *tracks){      
    if (fabs(trig_mu_z - trk.vz()) > electron_dz_cut_) continue;
    // TODO
    // if (UseDirectlyGenBeeK) {
    //   if (deltaR(EtaPhiK.first,EtaPhiK.second,trk.eta(),trk.phi())>DRgenCone)  continue;
    // } 
    bool mu_overlap = false;
    for (const pat::Muon & mu :*muons) {
      mu_overlap |= reco::deltaR(mu, trk) < mu_trk_min_dr_;
    }
    if(mu_overlap) continue;

    // TODO?
    // for (const pat::Electron & el : *electrons) {     
    //   if (deltaR(el.eta(),el.phi(),trk.eta(),trk.phi())<LepTrkExclusionCone)
    //     isE=true;
    // }
     
    // Store info in stores
    candidate_store.push_back(trk);
    candidate_ttracks.emplace_back(reco::TransientTrack(trk.pseudoTrack(), &(*bFieldHandle))); // Better way with the TTrack builder?
  }

  // Support collections for the objects, they contain pointers
  // to the vectors above, in an orderly fashion
  // cached vectors need to be built here AFTER the original vectors
  // are done, otherwise the vector self-growing feature will
  // mess up the pointer location and lead to segfaults
  CachedMuonCollection cached_muons = vec_to_cached(muon_store, muon_ttracks, true); // true: keep all muons no matter what
  CachedElectronCollection cached_electrons = vec_to_cached(electron_store, electron_ttracks, false); // false: keep only those electron that make a B candidate
  CachedCandidateCollection cached_candidates = vec_to_cached(candidate_store, candidate_ttracks, false);
  
  // Di-lepton caches (to be used later when building the B candidates)
  // the key is the index of the two leptons
  std::map< std::pair<size_t, size_t>, pat::CompositeCandidate > dimu_cache;
  std::map< std::pair<size_t, size_t>, pat::CompositeCandidate > diel_cache;

  // Create candidate collections
  auto b_kmumu = b_to_kmumu_builder_.build(*theBeamSpot, cached_muons, cached_candidates, dimu_cache);
  auto b_kee   = b_to_kee_builder_.build(*theBeamSpot, cached_electrons, cached_candidates, diel_cache);

  // Check if all the output candidates are empty 
  if(b_kmumu->empty() && b_kee->empty()) { 
    fill_empty(iEvent);
    return false;
  }  

  // Copy only candidates that have an index, meaning they formed a candidate somewhere
  CachedMuonCollection      used_muons;
  CachedElectronCollection  used_electrons;
  CachedCandidateCollection used_candidates;
  
  std::copy_if(
    cached_muons.begin(), cached_muons.end(), std::back_inserter(used_muons),
    [] (const CachedMuon &m) -> bool {return m.idx >= 0;}
    );
  std::copy_if(
    cached_electrons.begin(), cached_electrons.end(), std::back_inserter(used_electrons),
    [] (const CachedElectron &m) -> bool {return m.idx >= 0;}
    );
  std::copy_if(
    cached_candidates.begin(), cached_candidates.end(), std::back_inserter(used_candidates),
    [] (const CachedCandidate &m) -> bool {return m.idx >= 0;}
    );

  // Now, sort the used collection according to the indices
  std::sort(
    used_muons.begin(), used_muons.end(), 
    [] (const CachedMuon &m1, const CachedMuon &m2) -> bool {return m1.idx < m2.idx;}
    );
  std::sort(
    used_electrons.begin(), used_electrons.end(), 
    [] (const CachedElectron &m1, const CachedElectron &m2) -> bool {return m1.idx < m2.idx;}
    );
  std::sort(
    used_candidates.begin(), used_candidates.end(), 
    [] (const CachedCandidate &m1, const CachedCandidate &m2) -> bool {return m1.idx < m2.idx;}
    );

  // Finally, create the output collections from the pointer contents
  auto electrons_out = cached_to_vec(used_electrons);
  auto muons_out     = cached_to_vec(used_muons);
  auto cands_out     = cached_to_vec(used_candidates);

  iEvent.put(std::move(electrons_out), "electrons");
  iEvent.put(std::move(muons_out), "muons");
  iEvent.put(std::move(cands_out), "tracks");
  iEvent.put(std::move(b_kmumu), "BToKMuMu");
  iEvent.put(std::move(b_kee), "BToKEE");

  return true;
}


// ------------ method called once each job just before starting event loop  ------------
void 
ParkingNtupleMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ParkingNtupleMaker::endJob() 
{
}

// -------------------------------------
void ParkingNtupleMaker::beginRun( edm::Run const& run,  edm::EventSetup const& iSetup) {

  bool changed(true);
  if(hltPrescaleProvider_.init(run, iSetup, "HLT", changed)) {
    // if init returns TRUE, initialisation has succeeded!
    if(changed) {
      // The HLT config has actually changed wrt the previous Run
      std::cout << "Initalizing HLTConfigProvider"  << std::endl;
    }
  } 
  else {
    std::cout << " HLT config extraction failure with process name HLT" << std::endl;
  }
}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ParkingNtupleMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ParkingNtupleMaker);
