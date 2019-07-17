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
class ParkingNtupleMaker : public edm::one::EDFilter<edm::one::SharedResources,edm::one::WatchRuns>  {

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

  std::vector<std::vector<float>> track_DCA(std::vector<reco::TransientTrack> ttks);
  std::vector<GlobalVector>refit_tracks(TransientVertex myVertex,std::vector<reco::TransientTrack> tracks);

  const edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  const edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken_;
  const edm::EDGetTokenT<pat::ElectronCollection> electronsToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> tracksToken_;
  edm::EDGetToken muonsToken_;
//  edm::EDGetToken photonToken_;
  edm::EDGetTokenT<GlobalAlgBlkBxCollection> l1resultToken_;
  edm::EDGetToken l1MuonsToken_;
  vector<string> Seed_;
  edm::EDGetTokenT<edm::TriggerResults> trgresultsToken_; 
  edm::EDGetTokenT<vector<pat::TriggerObjectStandAlone>> trigobjectsToken_;
  vector<string> HLTPath_;
  HLTPrescaleProvider hltPrescaleProvider_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;

  // Candidate builders
  const BToKLLBuilder<pat::Muon, KinVtxFitter> b_to_kmumu_builder_;
  const BToKLLBuilder<pat::Electron, KinVtxFitter> b_to_kee_builder_;
  // const BToKStarLLBuilder<CachedMuon, KinVtxFitter> b_to_kstarmumu_builder_;
  // const BToKStarLLBuilder<CachedElectron, KinVtxFitter> b_to_kstaree_builder_;

  edm::Service<TFileService> fs;
  TTree * t1; 
  NtupleContent nt;
  int nevts=0;
  ///options
  bool data=true; 
  bool saveTracks=true; 
  bool saveHLT=true; 
  bool saveL1=true;
  bool saveOnlyHLTFires=false; 
  double track_pt_cut_forB=0;
  double min_muon_pt_cut_forB=0; 
  bool reconstructBMuMuK=true;
  double max_muon_pt_cut_forB=0;
  bool Pointing_constraint=false; 
  bool reconstructBMuMuKstar=true;
  double Pchi2BMuMuK=-1; 
  double MLLmax_Cut=100; 
  double MLLmin_Cut=-1;
  bool SkipEventWithNoBToMuMuK=false; 
  bool UseBeamspot=false; 
  bool AddeeK=false;
  double MBmax_Cut=100; 
  double MBmin_Cut=-1; 
  
  double  LepTrkExclusionCone=-1; 
  double EtaTrk_Cut=5; 
  bool AddLostTracks=true;
  double MKstarMin_Cut=0.5; 
  double MKstarMax_Cut=1.5; 
  
  std::string RefitTracks="none"; 
  bool UsePFeForCos=true; 
  bool OnlyKee=false;
  bool SkipEventWithNoBToMuMuKstar=false; 
  double Electron1PtCut=0;
  double Electron2PtCut=0; 
  double ElectronDzCut=0; 
  double TrgConeCut=-1;
  double MVAEl1Cut=-20; 
  double MVAEl2Cut=-20;
  double CosThetaCut=-1; 
  bool UseDirectlyGenBeeK=false; 
  double DRgenCone=10;
  int KIdToMatch=-1; 
  int LepIdToMatch=-1; 
  int BpdgIdToMatch=-1;
  bool IsResonantDecayToMatch=false; 
  std::string NtupleOutputClasses="auto"; 
  bool RetrieveMuFromTrk=false; 
  double maxPtTrk=0;
  double DzeeMaxCut=1000; 
  double PtBminCut=0;
  //internal
  std::vector<reco::TransientTrack> KTrack;
  std::vector<unsigned int> KTrack_index;
  unsigned int nmupairs=0; 
  float TrgmuDz=0,DRtrgMu=100; 
  int count=0;
  std::vector<reco::TransientTrack> muttks,ettks;
  std::vector<reco::TransientTrack> muTrack1,muTrack2,muPfTrack1,muPfTrack2;
  std::vector<std::pair<unsigned int,unsigned int>> used_muTrack_index,used_eTrack_index,used_muTrack_pfTrack_index;
  std::vector<reco::CandidatePtr> footprint; 
  int nmupfpairs=0;
  reco::TrackBase::Point vertex_point;

  TString * algoBitToName = new TString[512];
  l1t::L1TGlobalUtil *fGtUtil;

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
// photonToken_(consumes<std::vector<pat::Photon>>(iConfig.getParameter<edm::InputTag>("photons"))),
// Tracks_(consumes<std::vector<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracks"))),
  l1resultToken_(consumes<GlobalAlgBlkBxCollection>(iConfig.getParameter<edm::InputTag>("l1seed"))),
  l1MuonsToken_(consumes<l1t::MuonBxCollection>(iConfig.getParameter<edm::InputTag>("l1muons"))),
  Seed_(iConfig.getParameter<vector<string> >("Seed")),
  trgresultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag> ("triggerresults"))),
  trigobjectsToken_(consumes<vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag> ("triggerobjects"))),
  HLTPath_(iConfig.getParameter<vector<string> >("HLTPath")),
  hltPrescaleProvider_ (iConfig, consumesCollector(), *this),
  prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
  packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
  b_to_kmumu_builder_{iConfig.getParameter<edm::ParameterSet>("BToKMuMu")},
  b_to_kee_builder_{iConfig.getParameter<edm::ParameterSet>("BToKEE")} {

  edm::ParameterSet runParameters=iConfig.getParameter<edm::ParameterSet>("RunParameters");
  data=runParameters.getParameter<bool>("Data");
  saveTracks=runParameters.getParameter<bool>("SaveTracks");
  saveHLT=runParameters.getParameter<bool>("SaveHLT");
  saveL1=runParameters.getParameter<bool>("SaveL1");
  saveOnlyHLTFires=runParameters.getParameter<bool>("SaveResultsOnlyIfAPathFired");
  reconstructBMuMuK=runParameters.getParameter<bool>("ReconstructBMuMuK");
  reconstructBMuMuKstar=runParameters.getParameter<bool>("ReconstructBMuMuKstar");

  min_muon_pt_cut_forB=runParameters.getParameter<double>("MuonMinPtCut");
  max_muon_pt_cut_forB=runParameters.getParameter<double>("MuonMaxPtCut");
  track_pt_cut_forB=runParameters.getParameter<double>("TrackPtCutForB"); 
  Pchi2BMuMuK=runParameters.getParameter<double>("ProbBMuMuKcut");
  SkipEventWithNoBToMuMuK=runParameters.getParameter<bool>("SkipEventWithNoBToMuMuK");
  SkipEventWithNoBToMuMuKstar=runParameters.getParameter<bool>("SkipEventWithNoBToMuMuKstar");
  UseBeamspot=runParameters.getParameter<bool>("UseBeamspot");
  AddeeK=runParameters.getParameter<bool>("AddeeK");
  MLLmax_Cut=runParameters.getParameter<double>("MLLmax_Cut");
  MLLmin_Cut=runParameters.getParameter<double>("MLLmin_Cut");
  MBmax_Cut=runParameters.getParameter<double>("MBmax_Cut");
  MBmin_Cut=runParameters.getParameter<double>("MBmin_Cut");
  LepTrkExclusionCone=runParameters.getParameter<double>("LepTrkExclusionCone");
  EtaTrk_Cut=runParameters.getParameter<double>("EtaTrk_Cut");
  AddLostTracks=runParameters.getParameter<bool>("AddLostTracks");
  RefitTracks=runParameters.getParameter<std::string>("RefitTracks");
  UsePFeForCos=runParameters.getParameter<bool>("UsePFeForCos");
  OnlyKee=runParameters.getParameter<bool>("OnlyKee");
  ElectronDzCut=runParameters.getParameter<double>("ElectronDzCut");
  Electron1PtCut=runParameters.getParameter<double>("Electron1PtCut");
  Electron2PtCut=runParameters.getParameter<double>("Electron2PtCut");
  TrgConeCut=runParameters.getParameter<double>("TrgConeCut");
  MVAEl1Cut=runParameters.getParameter<double>("MVAEl1Cut");
  MVAEl2Cut=runParameters.getParameter<double>("MVAEl2Cut");
  CosThetaCut=runParameters.getParameter<double>("CosThetaCut");
  UseDirectlyGenBeeK=runParameters.getParameter<bool>("UseDirectlyGenBeeK");
  DRgenCone=runParameters.getParameter<double>("DRgenCone");
  KIdToMatch= runParameters.getParameter<int>("KIdToMatch");
  LepIdToMatch= runParameters.getParameter<int>("LepIdToMatch");
  BpdgIdToMatch=runParameters.getParameter<int>("BpdgIdToMatch");
  IsResonantDecayToMatch=runParameters.getParameter<bool>("IsResonantDecayToMatch");
  NtupleOutputClasses=runParameters.getParameter<std::string>("NtupleOutputClasses");
  maxPtTrk=runParameters.getParameter<double>("maxPtTrk");
  RetrieveMuFromTrk=runParameters.getParameter<bool>("RetrieveMuFromTrk");
  DzeeMaxCut=runParameters.getParameter<double>("DzeeMaxCut");
  PtBminCut=runParameters.getParameter<double>("PtBminCut");


  produces<pat::ElectronCollection>("electrons");
  produces<pat::MuonCollection>("muons");
  produces<pat::PackedCandidateCollection>("tracks");
  produces<pat::CompositeCandidateCollection>("BToKMuMu");
  produces<pat::CompositeCandidateCollection>("BToKEE");
  
  
//   fGtUtil = new l1t::L1TGlobalUtil(iConfig, 
//                                 consumesCollector(), 
//                                 *this, 
//                                 iConfig.getParameter<edm::InputTag>("l1results"), 
//                                 iConfig.getParameter<edm::InputTag>("l1results")
//                                 );
//   
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

std::vector<std::vector<float>> ParkingNtupleMaker::track_DCA(std::vector<reco::TransientTrack> ttks) {
  std::vector<std::vector<float>> dca;
  std::vector<float> def;
  def.push_back(-9999999);
  if(ttks.size()<2) {
    dca.push_back(def);
    return dca;
  }
  for(unsigned int tk1=0; tk1<ttks.size(); tk1++){
    TrajectoryStateClosestToPoint mu1TS = ttks[tk1].impactPointTSCP();
    std::vector<float> temp1;
    for(unsigned int tk2=tk1+1; tk2<ttks.size(); tk2++){
      TrajectoryStateClosestToPoint mu2TS = ttks[tk2].impactPointTSCP();      
      if (mu1TS.isValid() && mu2TS.isValid()) {
        ClosestApproachInRPhi cdca;
      cdca.calculate(mu1TS.theState(), mu2TS.theState());
        if (cdca.status()) temp1.push_back(cdca.distance());
        else temp1.push_back(-999999999);
      }
      else temp1.push_back(-99999999);
    }
    dca.push_back(temp1);         
  }  
  return dca;
}

std::vector<GlobalVector>
ParkingNtupleMaker::refit_tracks(TransientVertex myVertex,std::vector<reco::TransientTrack> tracks){
  std::auto_ptr<TrajectoryStateClosestToPoint> traj1;
  std::auto_ptr<TrajectoryStateClosestToPoint> traj2;
  GlobalPoint vtxPos(myVertex.position().x(), myVertex.position().y(), myVertex.position().z());
  GlobalVector gvmu1,gvmu2;
  if(myVertex.hasRefittedTracks()){
    std::vector<reco::TransientTrack> refited;
    refited=myVertex.refittedTracks();
    reco::TransientTrack* Track1 = &refited[0];
    reco::TransientTrack* Track2= &refited[1];
    traj1.reset(new TrajectoryStateClosestToPoint(Track1->trajectoryStateClosestToPoint(vtxPos)));
    traj2.reset(new TrajectoryStateClosestToPoint(Track2->trajectoryStateClosestToPoint(vtxPos)));    
    gvmu1=traj1->momentum();
    gvmu2=traj2->momentum();
  }         
  else {
    traj1.reset(new TrajectoryStateClosestToPoint(tracks[0].trajectoryStateClosestToPoint(vtxPos)));
    traj2.reset(new TrajectoryStateClosestToPoint(tracks[1].trajectoryStateClosestToPoint(vtxPos)));
    gvmu1=traj1->momentum();
    gvmu2=traj2->momentum();
  }                       
  std::vector<GlobalVector> gvmu;
  gvmu.push_back(gvmu1);
  gvmu.push_back(gvmu2);
  return gvmu;
}

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

  //clear
  footprint.clear();
  nt.ClearVariables();
  ettks.clear();
  muttks.clear();
  KTrack.clear();
  KTrack_index.clear();
  muTrack1.clear();
  muTrack2.clear();
  used_muTrack_index.clear();
  used_eTrack_index.clear();
  muPfTrack1.clear();
  muPfTrack2.clear();
  used_muTrack_pfTrack_index.clear();
  nmupairs=0;
  TrgmuDz=0;
  DRtrgMu=100;
  
  // per-event quantities
  nt.event=iEvent.id().event();
  nt.run_number=iEvent.id().run();
  nt.ls=iEvent.luminosityBlock();
  nt.beam_x=theBeamSpot->x0();
  nt.beam_y=theBeamSpot->y0();
  nt.beam_z=theBeamSpot->z0();
  nt.beam_ex= theBeamSpot->x0Error();
  nt.beam_ey= theBeamSpot->y0Error();
  nt.beam_ez= theBeamSpot->z0Error();
  float beam_xz=theBeamSpot->dxdz();
  float beam_yz =theBeamSpot->dydz();
  const  reco::Vertex firstGoodVertex=vertices->front();
  int firstGoodVertexIdx = 0;
  for (const reco::Vertex &vtx : *vertices) {
    bool isFake = vtx.isFake();
    if ( isFake || !vtx.isValid() ) continue;
    if (firstGoodVertexIdx==0){
      firstGoodVertexIdx=1; 
      nt.pvertex_x=vtx.x();
      nt.pvertex_y=vtx.y();
      nt.pvertex_z=vtx.z();
      nt.pvertex_ex=vtx.xError();
      nt.pvertex_ey=vtx.yError();
      nt.pvertex_ez=vtx.zError();
    }
    nt.vertex_x.push_back(vtx.x());
    nt.vertex_y.push_back(vtx.y());
    nt.vertex_z.push_back(vtx.z());
    nt.vertex_ex.push_back(vtx.xError());
    nt.vertex_ey.push_back(vtx.yError());
    nt.vertex_ez.push_back(vtx.zError());
    nt.vertex_chi.push_back(vtx.chi2());
    nt.vertex_ndof.push_back(vtx.ndof());
  }

  vertex_point.SetCoordinates(nt.pvertex_x,nt.pvertex_y,nt.pvertex_z);
   //K0 fit kalman
   //   KalmanVertexFitter theKalmanFitter(false);
   // TransientVertex K0vertex;
//   const pat::MET &theMet = met->front();
//   nt.ptmet=theMet.et(); nt.phimet=theMet.phi();
// 

   std::pair<float,float> EtaPhiE1(-10,-10),EtaPhiE2(-10,-10),EtaPhiK(-10,-10);
   if(! iEvent.isRealData() ){
     GeneratorBTree gen(prunedGenToken_,packedGenToken_,iEvent);
     gen.Fill(nt);
     if (UseDirectlyGenBeeK){
       gen.GenKLepLepEtaPhi(IsResonantDecayToMatch,BpdgIdToMatch,LepIdToMatch,KIdToMatch);
       EtaPhiK=gen.EtaPhiK();
       EtaPhiE1=gen.EtaPhiL();
       EtaPhiE2=gen.EtaPhiaL();
     }
   }
   if (UseDirectlyGenBeeK && (EtaPhiK.second==-10 || EtaPhiE1.second==-10 || EtaPhiE2.second==-10) ) {    
     fill_empty(iEvent);
     return false;
   }
  
  
  //save trigger related information
  
  //find max
  std::pair<float,float> TrgMu_EtaPhi(-100,-100);
  if (saveHLT || saveL1){
    HLTL1tree trigger(iEvent, l1resultToken_, l1MuonsToken_, trgresultsToken_, trigobjectsToken_);
    if (saveL1){ 
      trigger.L1objects(nt);
      trigger.L1trigger(algoBitToName,Seed_); 
      trigger.FillL1(nt);
    }
    if (saveHLT){
      trigger.HLTtrigger(HLTPath_, hltPrescaleProvider_); 
      trigger.HLTobjects(HLTPath_);
      if (saveOnlyHLTFires && !trigger.HLTPathFire()) {
        fill_empty(iEvent);
        return false;
      }
      trigger.FillHLT(nt); 
      trigger.FillObj(nt);
      if (trigger.HLTPathFire())
        TrgMu_EtaPhi=std::make_pair(trigger.GetHighestPtHLTObject()[1],trigger.GetHighestPtHLTObject()[2]);
    }
  } 
//    cout<<"here trg"<<endl;
  nevts++;
  
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
  for (const pat::Muon &mu : *muons){    
    // for (unsigned int i=0, n=mu.numberOfSourceCandidatePtrs(); i<n; ++i){
     //   footprint.push_back(mu.sourceCandidatePtr(i));
   // }
    nt.muon_pt.push_back(mu.pt());
    nt.muon_phi.push_back(mu.phi());
    nt.muon_eta.push_back(mu.eta());
    nt.muon_charge.push_back(mu.charge());
    nt.muon_global_flag.push_back(mu.isGlobalMuon());
    nt.muon_standalone_flag.push_back(mu.isStandAloneMuon());
    nt.muon_tracker_flag.push_back(mu.isTrackerMuon());
    nt.muon_vx.push_back(mu.vx()); 
    nt.muon_vy.push_back(mu.vy());
    nt.muon_vz.push_back(mu.vz()); 
    nt.muon_edz.push_back(mu.dzError());
    nt.muon_dz.emplace_back(mu.bestTrack()->dz(vertex_point));
    nt.muon_dxy.emplace_back(mu.bestTrack()->dxy(vertex_point));
    nt.muon_edxy.emplace_back(mu.dxyError());
    nt.muon_d0.emplace_back(mu.bestTrack()->d0());
    nt.muon_ed0.emplace_back(mu.bestTrack()->d0Error());
    nt.muon_trkpt.emplace_back(mu.bestTrack()->pt()); 
    nt.muon_trketa.emplace_back(mu.bestTrack()->eta());
    nt.muon_trkphi.emplace_back(mu.bestTrack()->phi());
    nt.muon_medium.push_back(mu.isMediumMuon());
    nt.muon_loose.push_back(mu.isLooseMuon());
    nt.muon_tight.push_back(mu.isTightMuon(firstGoodVertex));
    nt.muon_soft.push_back(mu.isSoftMuon(firstGoodVertex));
    if (deltaR(TrgMu_EtaPhi.first,TrgMu_EtaPhi.second,mu.eta(),mu.phi())<DRtrgMu){
      DRtrgMu=deltaR(TrgMu_EtaPhi.first,TrgMu_EtaPhi.second,mu.eta(),mu.phi());
      TrgmuDz=mu.vz();  
      nt.muon_trgIndex=nt.nmuons;
    }
    const reco::MuonPFIsolation&  isol=mu.pfIsolationR04();
    double mu_iso=(isol.sumChargedHadronPt+max(0.,isol.sumNeutralHadronEt+isol.sumPhotonEt-0.5*isol.sumPUPt))/mu.pt();
    nt.muon_iso.push_back(mu_iso);
    muttks.emplace_back(reco::TransientTrack(*mu.bestTrack(),&(*bFieldHandle)));
    nt.nmuons++;
    
    // Store info in stores
    muon_store.push_back(mu);    
    muon_ttracks.emplace_back(reco::TransientTrack(*mu.bestTrack(), &(*bFieldHandle))); // Better way with the TTrack builder?
  }
  if (DRtrgMu>TrgConeCut && TrgConeCut>0 && saveOnlyHLTFires) {
    fill_empty(iEvent);
    return false;
  }
    
  //electrons  
  //std::cout << "Electrons" << std::endl;
  trigger::size_type eindex=-1; 
  for(const pat::Electron &el : *electrons){
    eindex++;
    if(fabs(TrgmuDz - el.vz()) > ElectronDzCut) continue;
    bool is_lowpt = el.userInt("isLowPt");
    float pt = is_lowpt ? el.gsfTrack()->ptMode()  : el.pt();
    nt.el_pt.push_back(pt);
    nt.el_eta.push_back(is_lowpt ? el.gsfTrack()->etaMode() : el.eta());
    nt.el_phi.push_back(is_lowpt ? el.gsfTrack()->phiMode() : el.phi());
    nt.el_mva_biased.push_back(is_lowpt ? el.electronID("biased_seed") : -10);
    nt.el_mva_unbiased.push_back(is_lowpt ? el.electronID("unbiased_seed") : -10);
    nt.el_veto.push_back(  is_lowpt ? -10 : el.electronID("cutBasedElectronID-Fall17-94X-V2-veto"));
    nt.el_soft.push_back(  is_lowpt ? -10 : el.electronID("mvaEleID-Fall17-noIso-V2-wpLoose"));
    nt.el_medium.push_back(is_lowpt ? -10 : el.electronID("mvaEleID-Fall17-noIso-V2-wp90"));
    nt.el_tight.push_back( is_lowpt ? -10 : el.electronID("mvaEleID-Fall17-noIso-V2-wp80"));
    nt.el_mva_map_value.push_back(-10);
    nt.el_charge.push_back(el.charge());
    nt.el_dz.emplace_back(el.bestTrack()->dz(vertex_point));
    nt.el_dxy.emplace_back(el.bestTrack()->dxy(vertex_point));
    nt.el_edxy.push_back(el.dxyError());
    nt.el_edz.push_back(el.dzError());
    nt.el_vx.push_back(el.vx());
    nt.el_vy.push_back(el.vy());
    nt.el_vz.push_back(el.vz());
    double iso= el.pfIsolationVariables().sumChargedHadronPt+max(0.0,el.pfIsolationVariables().sumNeutralHadronEt+el.pfIsolationVariables().sumPhotonEt-0.5*el.pfIsolationVariables().sumPUPt)/pt;
    nt.el_iso.push_back(iso);
    nt.el_islowpt.push_back(is_lowpt);
    nt.el_trkpt.push_back(el.bestTrack()->pt());
    nt.el_trketa.push_back(el.bestTrack()->eta());
    nt.el_trkphi.push_back(el.bestTrack()->phi()); 
    ettks.emplace_back(reco::TransientTrack(*el.bestTrack(),&(*bFieldHandle)));
    nt.nel++;

    // Store info in stores
    electron_store.push_back(el);
    electron_ttracks.emplace_back(reco::TransientTrack(*el.bestTrack(), &(*bFieldHandle))); // Better way with the TTrack builder?
  }
  
  //std::cout << "Tracks:" << std::endl;
  for (const pat::PackedCandidate &trk : *tracks){      
    if (fabs(TrgmuDz - trk.vz()) > ElectronDzCut) continue;
    if (UseDirectlyGenBeeK){
      if (deltaR(EtaPhiK.first,EtaPhiK.second,trk.eta(),trk.phi())>DRgenCone)  continue;
    } 
    bool isMu=false;
    bool isE=false;
    for (const pat::Muon & mu :*muons){
      if (deltaR(mu.eta(),mu.phi(),trk.eta(),trk.phi())<LepTrkExclusionCone) isMu=true; 
    }
    for (const pat::Electron & el : *electrons) {     
      if (deltaR(el.eta(),el.phi(),trk.eta(),trk.phi())<LepTrkExclusionCone)
        isE=true;
    }
     
    if(isMu || isE ) continue;
    nt.track_pt.push_back(trk.pt());
    nt.track_eta.push_back(trk.eta());
    nt.track_phi.push_back(trk.phi());
    nt.track_charge.push_back(trk.charge());
     
    nt.track_highPurity.push_back(trk.trackHighPurity());
    nt.track_norm_chi2.emplace_back(trk.pseudoTrack().normalizedChi2());
    nt.track_dxy.push_back(trk.dxy(vertex_point));
    nt.track_dz.push_back(trk.dz(vertex_point));
    nt.track_validhits.push_back(trk.numberOfHits());
    nt.track_losthits.push_back(trk.lostInnerHits());
    nt.track_fromPV.push_back(trk.fromPV());
    nt.ntracks++; 

    // Store info in stores
    candidate_store.push_back(trk);
    candidate_ttracks.emplace_back(reco::TransientTrack(trk.pseudoTrack(), &(*bFieldHandle))); // Better way with the TTrack builder?
  }

  if(!reconstructBMuMuK && !reconstructBMuMuKstar){ 
    t1->Fill(); 
    fill_empty(iEvent);
    return false;
  }  

  // Support collections for the objects, they contain pointers
  // to the vectors above, in an orderly fashion
  // cached vectors need to be built here AFTER the original vectors
  // are done, otherwise the vector self-growing feature will
  // mess up the pointer location and lead to segfaults
  CachedMuonCollection cached_muons = vec_to_cached(muon_store, muon_ttracks, true);
  CachedElectronCollection cached_electrons = vec_to_cached(electron_store, electron_ttracks, true);
  CachedCandidateCollection cached_candidates = vec_to_cached(candidate_store, candidate_ttracks, true);
  
  // Di-lepton caches (to be used later when building the B candidates)
  // the key is the index of the two leptons
  std::map< std::pair<size_t, size_t>, pat::CompositeCandidate > dimu_cache;
  std::map< std::pair<size_t, size_t>, pat::CompositeCandidate > diel_cache;

  // Create candidate collections
  //std::cout << "B->Kuu" << std::endl;
  auto b_kmumu = b_to_kmumu_builder_.build(cached_muons, cached_candidates, dimu_cache);
  //std::cout << "B->Kee" << std::endl;
  auto b_kee   = b_to_kee_builder_.build(cached_electrons, cached_candidates, diel_cache);

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

  //muon pairs 
  if(!OnlyKee){
    for (int imu=0; imu<nt.nmuons; ++imu){
      if (!nt.muon_soft[imu]) continue;
      if(nt.muon_pt[imu]< min_muon_pt_cut_forB) continue;
      for (int imu2=imu+1; imu2<nt.nmuons; ++imu2){
        if (!nt.muon_soft[imu2]) continue; 
        if(nt.muon_pt[imu2]< min_muon_pt_cut_forB) continue;
        if (nt.muon_charge[imu]==nt.muon_charge[imu2]) continue;
        if (nt.muon_pt[imu]< max_muon_pt_cut_forB && nt.muon_pt[imu2]< max_muon_pt_cut_forB) continue;
        TLorentzVector vmu1,vmu2; 
        vmu1.SetPtEtaPhiM(nt.muon_pt[imu], nt.muon_eta[imu] ,nt.muon_phi[imu] , 0.105);
        vmu2.SetPtEtaPhiM(nt.muon_pt[imu2],nt.muon_eta[imu2],nt.muon_phi[imu2], 0.105);
        if((vmu1+vmu2).M()<MLLmin_Cut || (vmu1+vmu2).M()>MLLmax_Cut) continue;
        muTrack1.push_back(muttks[imu]); 
        muTrack2.push_back(muttks[imu2]);
        used_muTrack_index.emplace_back(imu,imu2);
        nmupairs++;
      }
    }
  } 
 
  //e pairs
  if(AddeeK || OnlyKee ){
    for(int iel=0; iel<nt.nel; ++iel){
      if (UseDirectlyGenBeeK && !iEvent.isRealData() ){
        if ( deltaR(EtaPhiE1.first, EtaPhiE1.second, nt.el_eta[iel], nt.el_phi[iel]) > DRgenCone && \
             deltaR(EtaPhiE2.first, EtaPhiE2.second, nt.el_eta[iel], nt.el_phi[iel]) > DRgenCone) 
             continue;
      }
      for (int iel2=iel+1; iel2<nt.nel; ++iel2){
        if (nt.el_charge[iel]==nt.el_charge[iel2]) continue;
        if (nt.el_pt[iel2] < Electron1PtCut && \
            nt.el_pt[iel]  < Electron1PtCut) 
            continue;
        if (UseDirectlyGenBeeK && !iEvent.isRealData() ){
          if (deltaR(EtaPhiE1.first, EtaPhiE1.second, nt.el_eta[iel2], nt.el_phi[iel2]) > DRgenCone && \
              deltaR(EtaPhiE2.first, EtaPhiE2.second, nt.el_eta[iel2], nt.el_phi[iel2]) > DRgenCone) 
              continue;
        }
        if (nt.el_islowpt[iel] && nt.el_islowpt[iel2]){
          if (nt.el_mva_unbiased[iel]  < MVAEl1Cut && \
              nt.el_mva_unbiased[iel2] < MVAEl1Cut) 
              continue;
        }
        if ( fabs(nt.el_vz[iel]-nt.el_vz[iel2]) > DzeeMaxCut) continue;
        TLorentzVector vel1,vel2; 
        vel1.SetPtEtaPhiM(nt.el_pt[iel] , nt.el_eta[iel] , nt.el_phi[iel] , 0.000511);
        vel2.SetPtEtaPhiM(nt.el_pt[iel2], nt.el_eta[iel2], nt.el_phi[iel2], 0.000511);
        if( (vel1+vel2).M() < MLLmin_Cut || \
            (vel1+vel2).M() > MLLmax_Cut) 
            continue; 
        muTrack1.push_back(ettks[iel]); 
        muTrack2.push_back(ettks[iel2]);
        used_eTrack_index.emplace_back(iel,iel2);
      }
    }
  }

  if (used_muTrack_index.size()==0 && used_eTrack_index.size()==0){
    if      (reconstructBMuMuK     && SkipEventWithNoBToMuMuK     && !reconstructBMuMuKstar) {
      fill_empty(iEvent);
      return false;
    }
    else if (reconstructBMuMuKstar && SkipEventWithNoBToMuMuKstar && !reconstructBMuMuK    ) {
      fill_empty(iEvent);
      return false;
    }
    else if (reconstructBMuMuK     && SkipEventWithNoBToMuMuK     && SkipEventWithNoBToMuMuKstar && reconstructBMuMuKstar) {
      fill_empty(iEvent);
      return false;
    }
  }
 
  if (RetrieveMuFromTrk){
    nmupfpairs=nmupairs;
    for (const pat::Muon & mu : *muons){
      if (fabs(TrgmuDz-mu.vz())>ElectronDzCut) continue;
      if (mu.pt()<min_muon_pt_cut_forB) continue;
      if (fabs(mu.eta())>2.5) continue;
      if (!mu.isSoftMuon(firstGoodVertex)) continue;
      for(const pat::PackedCandidate &trk : *tracks){
        if(mu.charge() == trk.charge()) continue;
        if(trk.charge()==0) continue;
        if (trk.pt()< min_muon_pt_cut_forB) continue;
        if (fabs(trk.eta())>EtaTrk_Cut) continue;
        if (fabs(TrgmuDz-trk.vz())>ElectronDzCut) continue;
        if(trk.pt()>maxPtTrk) continue;
        TLorentzVector vmu1,vmu2;
        vmu1.SetPtEtaPhiM(mu.pt() , mu.eta() , mu.phi() , 0.105);
        vmu1.SetPtEtaPhiM(trk.pt(), trk.eta(), trk.phi(), 0.105);
        if ((vmu1+vmu2).M()<MLLmin_Cut || \
            (vmu1+vmu2).M()>MLLmax_Cut ) 
            continue;
        muTrack1.emplace_back(reco::TransientTrack(*mu.bestTrack(),   &(*bFieldHandle) ) );
        muTrack2.emplace_back(reco::TransientTrack(trk.pseudoTrack(), &(*bFieldHandle) ) );
        used_muTrack_pfTrack_index.emplace_back(std::make_pair( &mu-&(muons->at(0)), &trk-&(*tracks)[0]));      
        nmupfpairs++;
      } 
    }
  }
  
  int index=-1;
  if ( (used_muTrack_index.size()>0 || used_eTrack_index.size()>0) && \
       (reconstructBMuMuK || reconstructBMuMuKstar) ){
    for(const pat::PackedCandidate &trk : *tracks){
      reco::Track ptrk = trk.pseudoTrack();
      index++;
      if (ptrk.pt() < track_pt_cut_forB) continue;
      if (fabs(ptrk.eta()) > EtaTrk_Cut) continue;
      if (UseDirectlyGenBeeK){
        if (deltaR(EtaPhiK.first,EtaPhiK.second,ptrk.eta(),ptrk.phi())>DRgenCone)
   	      continue;
      }
      if (fabs(TrgmuDz-ptrk.vz())>ElectronDzCut) continue;
      bool isMu = false; 
      bool isE  = false;
      for (const pat::Muon & mu : *muons){
        if (deltaR(mu.eta(),mu.phi(),ptrk.eta(),ptrk.phi())<LepTrkExclusionCone) 
        isMu=true; 
      }
      // sara: don't understand
      for (const pat::Electron & el: *electrons) {     
        if (deltaR(el.eta(),el.phi(),ptrk.eta(),ptrk.phi())<LepTrkExclusionCone) 
          isE=true;
      } 
      if(isMu || isE ) continue;
      KTrack.emplace_back(reco::TransientTrack(ptrk,&(*bFieldHandle)));
      KTrack_index.push_back(&trk-&(*tracks)[0]);  
    }
  }



  //add pf
  if (used_muTrack_pfTrack_index.size()>0)
    used_muTrack_index.insert(used_muTrack_index.end(), used_muTrack_pfTrack_index.begin(), used_muTrack_pfTrack_index.end());
  //add ee chanel 
  if (used_eTrack_index.size()>0)
    used_muTrack_index.insert(used_muTrack_index.end(), used_eTrack_index.begin(), used_eTrack_index.end());
 
  ///building B 
  if (used_muTrack_index.size()>0 && KTrack_index.size()>0 && reconstructBMuMuK){
    BKlldecay Kll(nmupairs,
                  used_muTrack_index,
                  muTrack1,
                  muTrack2,
                  KTrack,
                  beam_xz,
                  beam_yz,
                  RefitTracks
                  );
    Kll.ProbCut(Pchi2BMuMuK); 
    Kll.CosCut(CosThetaCut); 
    Kll.MassCuts(MBmin_Cut,MBmax_Cut); 
    Kll.PtBCut(PtBminCut); 
    Kll.SetNumTrkForMu(nmupfpairs); 
    Kll.Fill(nt);  
  } 

  if (used_muTrack_index.size()>0 && KTrack_index.size()>0 && reconstructBMuMuKstar){
    BKstarlldecay Kstarll(nmupairs,
                          used_muTrack_index,
                          muTrack1,
                          muTrack2,
                          KTrack,
                          beam_xz,
                          beam_yz,
                          RefitTracks
                          );
    Kstarll.CombineTracks(0.600,1.100);
    Kstarll.Fill(nt);
 } 

  if (used_muTrack_index.size()>0 && KTrack_index.size()>0 && reconstructBMuMuKstar){
    BKstarlldecay phill(nmupairs,
                        used_muTrack_index,
                        muTrack1,
                        muTrack2,
                        KTrack,
                        beam_xz,
                        beam_yz,
                        RefitTracks
                        );
    phill.Runphill(true); 
    phill.CombineTracks(0.600,1.100);
    phill.FillPhill(nt);
  }

  //skip useless events
  if (nt.NRb_mass.size()==0 && nt.NRbks_mass.size()==0){
    if      (reconstructBMuMuK     && SkipEventWithNoBToMuMuK     && !reconstructBMuMuKstar) {
      fill_empty(iEvent);
      return false;
    }
    else if (reconstructBMuMuKstar && SkipEventWithNoBToMuMuKstar && !reconstructBMuMuK    ) {
      fill_empty(iEvent);
      return false;
    }
    else if (reconstructBMuMuK     && SkipEventWithNoBToMuMuK     && SkipEventWithNoBToMuMuKstar && reconstructBMuMuKstar) {
      fill_empty(iEvent);
      return false;
    }
  }  
  
  if (NtupleOutputClasses=="flat") nt.Flatting();
  t1->Fill();

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
  t1=fs->make<TTree>("mytree","mytree");
  
  nt.SetTree(t1);
  TString ToSave="";
  if (NtupleOutputClasses=="auto" || NtupleOutputClasses=="flat" || NtupleOutputClasses=="lite" ){
    if (NtupleOutputClasses=="flat") ToSave+="Flat_";
    if (NtupleOutputClasses=="lite") ToSave+="Lite_";
    if (reconstructBMuMuK) ToSave+="KLL_";
    if (reconstructBMuMuKstar) ToSave+="KstarLL_";
    if (!data) ToSave+="GEN_";
    if (saveTracks) ToSave+="TRK_";
  } else ToSave=NtupleOutputClasses;

  nt.SetNtupleVariables(ToSave);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
ParkingNtupleMaker::endJob() 
{
  std::cout<<"Processed "<<nevts<<std::endl;
  if (NtupleOutputClasses=="flat"){
    for (unsigned int ib=0; ib<nt.NRb_pt_eta_phi.size(); ib++){
      nt.NRb_pt.emplace_back(nt.NRb_pt_eta_phi[ib][0]);
      cout<<nt.NRb_pt_eta_phi[ib][0]<<endl;
      nt.NRb_eta.emplace_back(nt.NRb_pt_eta_phi[ib][1]);
      nt.NRb_phi.emplace_back(nt.NRb_pt_eta_phi[ib][2]);
      nt.NRb_Kpt.emplace_back(nt.NRb_Kpt_eta_phi[ib][0]);
      nt.NRb_Keta.emplace_back(nt.NRb_Kpt_eta_phi[ib][1]);
      nt.NRb_Kphi.emplace_back(nt.NRb_Kpt_eta_phi[ib][2]);
      nt.NRb_l1pt.emplace_back(nt.NRb_l1pt_eta_phi[ib][0]);
      nt.NRb_l1eta.emplace_back(nt.NRb_l1pt_eta_phi[ib][1]);
      nt.NRb_l1phi.emplace_back(nt.NRb_l1pt_eta_phi[ib][2]);
      nt.NRb_l2pt.emplace_back(nt.NRb_l2pt_eta_phi[ib][0]);
      nt.NRb_l2eta.emplace_back(nt.NRb_l2pt_eta_phi[ib][1]);
      nt.NRb_l2phi.emplace_back(nt.NRb_l2pt_eta_phi[ib][2]);
    }
   for (unsigned int ib=0; ib<nt.NRbks_pt_eta_phi.size(); ib++){
      nt.NRbks_pt.emplace_back(nt.NRbks_pt_eta_phi[ib][0]);
      nt.NRbks_eta.emplace_back(nt.NRbks_pt_eta_phi[ib][1]);
      nt.NRbks_phi.emplace_back(nt.NRbks_pt_eta_phi[ib][2]);
      nt.NRbks_Kpt.emplace_back(nt.NRbks_Kpt_eta_phi[ib][0]);
      nt.NRbks_Keta.emplace_back(nt.NRbks_Kpt_eta_phi[ib][1]);
      nt.NRbks_Kphi.emplace_back(nt.NRbks_Kpt_eta_phi[ib][2]);
      nt.NRbks_Pipt.emplace_back(nt.NRbks_Pipt_eta_phi[ib][0]);
      nt.NRbks_Pieta.emplace_back(nt.NRbks_Pipt_eta_phi[ib][1]);
      nt.NRbks_Piphi.emplace_back(nt.NRbks_Pipt_eta_phi[ib][2]);
      nt.NRbks_l1pt.emplace_back(nt.NRbks_l1pt_eta_phi[ib][0]);
      nt.NRbks_l1eta.emplace_back(nt.NRbks_l1pt_eta_phi[ib][1]);
      nt.NRbks_l1phi.emplace_back(nt.NRbks_l1pt_eta_phi[ib][2]);
      nt.NRbks_l2pt.emplace_back(nt.NRbks_l2pt_eta_phi[ib][0]);
      nt.NRbks_l2eta.emplace_back(nt.NRbks_l2pt_eta_phi[ib][1]);
      nt.NRbks_l2phi.emplace_back(nt.NRbks_l2pt_eta_phi[ib][2]);
    }
  }
}

// -------------------------------------
void ParkingNtupleMaker::beginRun( edm::Run const& run,  edm::EventSetup const& iSetup) {

  bool changed(true);
  if (hltPrescaleProvider_.init(run,iSetup,"HLT",changed)) {
    // if init returns TRUE, initialisation has succeeded!
    if (changed) {
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
