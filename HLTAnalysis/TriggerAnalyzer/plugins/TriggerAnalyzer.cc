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
#include <iostream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include <vector>
#include "TTree.h"
#include "TMath.h"
#include <string>
#include <iostream>
#include "DataFormats/Common/interface/Ref.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TLorentzVector.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1TUtmTriggerMenu.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "TripleTrackKinFit.h"
#include "GeneratorBTree.h"
#include "NtupleContent.h"
#include "LowPtElObjects.h"
#include "PFelCollection.h"
#include "HLTL1tree.h"
#include "BKlldecay.h"
#include "BKstarlldecay.h"
#include "helper.h"

using namespace std;


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.
template<typename T1>
class TriggerAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {

  typedef std::vector<T1> T1Collection;
  typedef edm::Ref<T1Collection> T1Ref;
  typedef edm::AssociationMap<edm::OneToValue<std::vector<T1>, float > > T1IsolationMap;

public:
  explicit TriggerAnalyzer(const edm::ParameterSet&);
  ~TriggerAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  std::vector<std::vector<float>> track_DCA(std::vector<reco::TransientTrack> ttks);
  std::vector<GlobalVector>refit_tracks(TransientVertex myVertex,std::vector<reco::TransientTrack> tracks);
  float Dphi(float phi1,float phi2);
  float DR(float eta1,float phi1,float eta2, float phi2);

  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken_;
  edm::EDGetToken electronsToken_;
  edm::EDGetToken lowPtElectronsToken_;
  edm::EDGetToken lowPtGsfTracksToken_;
  edm::EDGetToken pfElectronsToken_;
  edm::EDGetToken muonsToken_;
  edm::EDGetToken jetsToken_;
  edm::EDGetToken metToken_;
//  edm::EDGetToken photonToken_;
  edm::EDGetToken PFCands_;
  edm::EDGetToken LostTracks_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapVetoToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapSoftToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapMediumToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapTightToken_;
  edm::EDGetTokenT<edm::ValueMap<int> > elIdMapValueToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > eleBWPToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > eleUnBWPToken_;
  edm::EDGetTokenT<GlobalAlgBlkBxCollection> l1resultToken_;
  edm::EDGetToken l1MuonsToken_;
  edm::EDGetToken l1JetsToken_;
  edm::EDGetToken l1MetToken_;
  vector<string> Seed_;
  edm::EDGetTokenT<edm::TriggerResults> trgresultsToken_; 
  edm::EDGetTokenT<vector<pat::TriggerObjectStandAlone>> trigobjectsToken_;
  vector<string> HLTPath_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;

  edm::Service<TFileService> fs;
  TTree * t1; 
  NtupleContent nt;
  int nevts=0;
  ///options
   bool data=true; bool saveTracks=true; bool saveHLT=true; bool saveL1=true;
   bool saveOnlyHLTFires=false; double track_pt_cut_forB=0;
   double min_muon_pt_cut_forB=0; bool reconstructBMuMuK=true;
   double max_muon_pt_cut_forB=0;
   bool Pointing_constraint=false; bool reconstructBMuMuKstar=true;
  double Pchi2BMuMuK=-1; double MLLmax_Cut=100; double MLLmin_Cut=-1;
  bool SkipEventWithNoBToMuMuK=false; bool UseBeamspot=false; bool AddeeK=false;
  double MBmax_Cut=100; double MBmin_Cut=-1; 
  double  LepTrkExclusionCone=-1; double EtaTrk_Cut=5; bool AddLostTracks=true;
  double MKstarMin_Cut=0.5; double MKstarMax_Cut=1.5; 
  std::string RefitTracks="none"; bool UsePFeForCos=true; bool OnlyKee=false;
  bool SkipEventWithNoBToMuMuKstar=false; double Electron1PtCut=0;
  double Electron2PtCut=0; double ElectronDzCut=0; double TrgConeCut=-1;
  bool IsLowpTE=false; double MVAEl1Cut=-20; double MVAEl2Cut=-20;
  double CosThetaCut=-1; bool UseDirectlyGenBeeK=false; double DRgenCone=10;
  int KIdToMatch=-1; int LepIdToMatch=-1; int BpdgIdToMatch=-1;
  bool IsResonantDecayToMatch=false; bool AddLowPtElAsCol=false;
  bool AddLowPtGsfTrkAsCol=false; bool AddPFElAsCol=false;
  std::string NtupleOutputClasses="auto"; bool CombineElCol=false;
  double CombineCone=0; bool RetrieveMuFromTrk=false; double maxPtTrk=0;
  double DzeeMaxCut=1000; double PtBminCut=0;
  //internal
  std::vector<std::pair<float,float>> PFe_EtaPhi;
  std::vector<reco::TransientTrack> KTrack;
  std::vector<unsigned int> KTrack_index;
  unsigned int nmupairs=0; float TrgmuDz=0,DRtrgMu=100; int count=0;
  std::vector<reco::TransientTrack> muttks,ettks;
  std::vector<reco::TransientTrack> muTrack1,muTrack2,muPfTrack1,muPfTrack2;
  std::vector<std::pair<unsigned int,unsigned int>> used_muTrack_index,used_eTrack_index,used_muTrack_pfTrack_index;
  std::vector<reco::CandidatePtr> footprint; std::vector<pat::PackedCandidate> tracks;
  int nmupfpairs=0;
  reco::TrackBase::Point vertex_point;

 TString * algoBitToName = new TString[512];
      // ----------member data ---------------------------
};


// constructors and destructor
//
template<typename T1>
TriggerAnalyzer<T1>::TriggerAnalyzer(const edm::ParameterSet& iConfig): 
 beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter <edm::InputTag>("beamSpot"))),
 vtxToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
 electronsToken_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>  ("electrons"))),
 lowPtElectronsToken_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>  ("lowptElectrons"))),
 lowPtGsfTracksToken_(consumes<vector<reco::GsfTrack>>(iConfig.getParameter<edm::InputTag>  ("lowptGsftracks"))),
  pfElectronsToken_(consumes<vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>  ("pfElectrons"))),
 muonsToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
 jetsToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>  ("jets"))),
 metToken_(consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("met"))),
// photonToken_(consumes<std::vector<pat::Photon>>(iConfig.getParameter<edm::InputTag>("photons"))),
 PFCands_(consumes<std::vector<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("PFCands"))),
 LostTracks_(consumes<std::vector<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("losttracks"))),
// Tracks_(consumes<std::vector<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracks"))),
 eleIdMapVetoToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleIdMapVeto"))),
 eleIdMapSoftToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleIdMapSoft"))),
 eleIdMapMediumToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleIdMapMedium"))),
 eleIdMapTightToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleIdMapTight"))),
 elIdMapValueToken_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("eleIdMapValue"))),
 eleBWPToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("eleBiasedWP"))),
  eleUnBWPToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("eleUnbiasedWP"))),
 l1resultToken_(consumes<GlobalAlgBlkBxCollection>(iConfig.getParameter<edm::InputTag>("l1seed"))),
 l1MuonsToken_(consumes<l1t::MuonBxCollection>(iConfig.getParameter<edm::InputTag>("l1muons"))),
 l1JetsToken_(consumes<l1t::JetBxCollection>(iConfig.getParameter<edm::InputTag>("l1jets"))),
 l1MetToken_(consumes<BXVector<l1t::EtSum> >(iConfig.getParameter<edm::InputTag>("l1met"))),
 Seed_(iConfig.getParameter<vector<string> >("Seed")),
 trgresultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag> ("triggerresults"))),
 trigobjectsToken_(consumes<vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag> ("triggerobjects"))),
 HLTPath_(iConfig.getParameter<vector<string> >("HLTPath")),
 prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed")))
{

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
 IsLowpTE=runParameters.getParameter<bool>("IsLowpTE");
 MVAEl1Cut=runParameters.getParameter<double>("MVAEl1Cut");
 MVAEl2Cut=runParameters.getParameter<double>("MVAEl2Cut");
 CosThetaCut=runParameters.getParameter<double>("CosThetaCut");
 UseDirectlyGenBeeK=runParameters.getParameter<bool>("UseDirectlyGenBeeK");
 DRgenCone=runParameters.getParameter<double>("DRgenCone");
 KIdToMatch= runParameters.getParameter<int>("KIdToMatch");
 LepIdToMatch= runParameters.getParameter<int>("LepIdToMatch");
 BpdgIdToMatch=runParameters.getParameter<int>("BpdgIdToMatch");
 IsResonantDecayToMatch=runParameters.getParameter<bool>("IsResonantDecayToMatch");
 AddLowPtElAsCol=runParameters.getParameter<bool>("AddLowPtElAsCol");
 AddLowPtGsfTrkAsCol=runParameters.getParameter<bool>("AddLowGsfTrkAsCol");
 AddPFElAsCol=runParameters.getParameter<bool>("AddPFElAsCol");
 NtupleOutputClasses=runParameters.getParameter<std::string>("NtupleOutputClasses");
 CombineElCol=runParameters.getParameter<bool>("CombineElCol");
 CombineCone=runParameters.getParameter<double>("CombineCone");
 maxPtTrk=runParameters.getParameter<double>("maxPtTrk");
 RetrieveMuFromTrk=runParameters.getParameter<bool>("RetrieveMuFromTrk");
 DzeeMaxCut=runParameters.getParameter<double>("DzeeMaxCut");
 PtBminCut=runParameters.getParameter<double>("PtBminCut");
}

template<typename T1>
TriggerAnalyzer<T1>::~TriggerAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------

template<typename T1> 
float TriggerAnalyzer<T1>::Dphi(float phi1,float phi2){
float result = phi1 - phi2;
 while (result > float(M_PI)) result -= float(2*M_PI);
    while (result <= -float(M_PI)) result += float(2*M_PI);
return result;

}

template<typename T1> 
float TriggerAnalyzer<T1>::DR(float eta1,float phi1,float eta2, float phi2){
  return TMath::Sqrt((eta1-eta2)*(eta1-eta2)+Dphi(phi1,phi2)*Dphi(phi1,phi2));
}



template<typename T1> 
std::vector<std::vector<float>> TriggerAnalyzer<T1>::track_DCA(std::vector<reco::TransientTrack> ttks) {
 std::vector<std::vector<float>> dca;
 std::vector<float> def;
 def.push_back(-9999999);
 if(ttks.size()<2) {dca.push_back(def); return dca; }
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

template<typename T1>
std::vector<GlobalVector>
TriggerAnalyzer<T1>::refit_tracks(TransientVertex myVertex,std::vector<reco::TransientTrack> tracks){
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
        gvmu1=traj1->momentum();  gvmu2=traj2->momentum(); 
     }         
     else {
        traj1.reset(new TrajectoryStateClosestToPoint(tracks[0].trajectoryStateClosestToPoint(vtxPos)));
        traj2.reset(new TrajectoryStateClosestToPoint(tracks[1].trajectoryStateClosestToPoint(vtxPos)));
        gvmu1=traj1->momentum();  gvmu2=traj2->momentum();
     }                       
   std::vector<GlobalVector> gvmu; gvmu.push_back(gvmu1); gvmu.push_back(gvmu2);
   return gvmu;
}

template<typename T1>
void
TriggerAnalyzer<T1>::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace std;
  //Get a few collections to apply basic electron ID
  //Get data
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(beamSpotToken_,theBeamSpot); 
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  //continue if there are no vertices
  if (vertices->size()==0) return;
  edm::Handle<std::vector<pat::Electron>> electrons;
  iEvent.getByToken(electronsToken_, electrons); 
  edm::Handle<std::vector<pat::Electron>> lowpte;
  iEvent.getByToken(lowPtElectronsToken_,lowpte);
  edm::Handle<std::vector<pat::Jet>> jets;
  iEvent.getByToken(jetsToken_, jets);
  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonsToken_,muons);
  edm::Handle<std::vector<pat::MET>> met;
  iEvent.getByToken(metToken_,met);
  edm::Handle<vector<pat::PackedCandidate>> tracks1;
  iEvent.getByToken(PFCands_, tracks1);
  edm::Handle<vector<pat::PackedCandidate>> tracks2;
  iEvent.getByToken(LostTracks_, tracks2);
  edm::Handle<edm::ValueMap<bool> > ele_veto_id;
  iEvent.getByToken(eleIdMapVetoToken_ ,ele_veto_id);
  edm::Handle<edm::ValueMap<bool> > ele_soft_id;
  iEvent.getByToken(eleIdMapSoftToken_ ,ele_soft_id);
  edm::Handle<edm::ValueMap<bool> > ele_medium_id;
  iEvent.getByToken(eleIdMapMediumToken_ ,ele_medium_id);
  edm::Handle<edm::ValueMap<bool> > ele_tight_id;
  iEvent.getByToken(eleIdMapTightToken_ ,ele_tight_id);
  edm::Handle<edm::ValueMap<int> > ele_mva_id_value;
  iEvent.getByToken( elIdMapValueToken_ ,ele_mva_id_value);
  edm::Handle<edm::ValueMap<float> > ele_mva_wp_biased;
  iEvent.getByToken( eleBWPToken_ ,ele_mva_wp_biased);
  edm::Handle<edm::ValueMap<float> > ele_mva_wp_unbiased;
  iEvent.getByToken( eleUnBWPToken_ ,ele_mva_wp_unbiased);
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
        int itrig = index; algoBitToName[itrig] = TString( trigName );
      }
      count++;
    }
  }
  //clear
  footprint.clear(); nt.ClearVariables(); PFe_EtaPhi.clear(); ettks.clear();
  muttks.clear(); KTrack.clear(); KTrack_index.clear(); tracks.clear();
  muTrack1.clear(); muTrack2.clear(); used_muTrack_index.clear();
  used_eTrack_index.clear(); muPfTrack1.clear(); muPfTrack2.clear();
  used_muTrack_pfTrack_index.clear();
  nmupairs=0; TrgmuDz=0; DRtrgMu=100;
  
  nt.event=iEvent.id().event();
  nt.run_number=iEvent.id().run();
  nt.ls=iEvent.luminosityBlock();
  nt.beam_x=theBeamSpot->x0(); nt.beam_y=theBeamSpot->y0(); 
  nt.beam_z=theBeamSpot->z0(); nt.beam_ex= theBeamSpot->x0Error(); 
  nt.beam_ey= theBeamSpot->y0Error(); nt.beam_ez= theBeamSpot->z0Error();
  float beam_xz=theBeamSpot->dxdz(); float beam_yz =theBeamSpot->dydz();
  const  reco::Vertex firstGoodVertex=vertices->front();
  int firstGoodVertexIdx = 0;
  for (const reco::Vertex &vtx : *vertices) {
    bool isFake = vtx.isFake();
    if ( isFake || !vtx.isValid() ) continue;
    if (firstGoodVertexIdx==0){
      firstGoodVertexIdx=1; 
      nt.pvertex_x=vtx.x(); nt.pvertex_y=vtx.y(); nt.pvertex_z=vtx.z();
      nt.pvertex_ex=vtx.xError(); nt.pvertex_ey=vtx.yError(); nt.pvertex_ez=vtx.zError();
    }
    nt.vertex_x.push_back(vtx.x()); nt.vertex_y.push_back(vtx.y()); 
    nt.vertex_z.push_back(vtx.z()); nt.vertex_ex.push_back(vtx.xError());
    nt.vertex_ey.push_back(vtx.yError()); nt.vertex_ez.push_back(vtx.zError());
    nt.vertex_chi.push_back(vtx.chi2()); nt.vertex_ndof.push_back(vtx.ndof());
  }

   vertex_point.SetCoordinates(nt.pvertex_x,nt.pvertex_y,nt.pvertex_z);
   //K0 fit kalman
   //   KalmanVertexFitter theKalmanFitter(false);
   // TransientVertex K0vertex;
  const pat::MET &theMet = met->front();
  nt.ptmet=theMet.et(); nt.phimet=theMet.phi();

  for (const pat::Jet &jet : *jets){
    if (fabs(jet.eta())>2.5) continue;
    nt.jet_pt.push_back(jet.pt()); nt.jet_eta.push_back(jet.eta());
    nt.jet_phi.push_back(jet.phi());
    nt.jet_cEmEF.push_back(jet.chargedEmEnergyFraction());
    nt.jet_cHEF.push_back(jet.chargedHadronEnergyFraction());
    nt.jet_cHMult.push_back(jet.chargedHadronMultiplicity());
    nt.jet_cMuEF.push_back(jet.chargedMuEnergyFraction());
    nt.jet_cMult.push_back(jet.chargedMultiplicity());
    nt.jet_MuEF.push_back(jet.muonEnergyFraction());
    nt.jet_eEF.push_back(jet.electronEnergyFraction());
    nt.jet_nEmEF.push_back(jet.neutralEmEnergyFraction());
    nt.jet_nHEF.push_back(jet.neutralHadronEnergyFraction());
    nt.jet_nMult.push_back(jet.neutralMultiplicity());
    nt.jet_pEF.push_back(jet.photonEnergyFraction());
    nt.njets++;
  }

   std::pair<float,float> EtaPhiE1(-10,-10),EtaPhiE2(-10,-10),EtaPhiK(-10,-10);
   if(!data){
     GeneratorBTree gen(prunedGenToken_,packedGenToken_,iEvent);
     gen.Fill(nt);
     if (UseDirectlyGenBeeK){
       gen.GenKLepLepEtaPhi(IsResonantDecayToMatch,BpdgIdToMatch,LepIdToMatch,KIdToMatch);
       EtaPhiK=gen.EtaPhiK(); EtaPhiE1=gen.EtaPhiL(); EtaPhiE2=gen.EtaPhiaL();
     }
   }
   if (UseDirectlyGenBeeK && (EtaPhiK.second==-10 || EtaPhiE1.second==-10 || EtaPhiE2.second==-10) ) return;
  //find max
  std::pair<float,float> TrgMu_EtaPhi(-100,-100);
  if (saveHLT || saveL1){
    HLTL1tree trigger(iEvent, l1resultToken_, l1MuonsToken_, l1JetsToken_, l1MetToken_, trgresultsToken_, trigobjectsToken_);
    if (saveL1){ 
      trigger.L1objetcs(nt);
      trigger.L1trigger(algoBitToName,Seed_); trigger.FillL1(nt);
    }
    if (saveHLT){
      trigger.HLTtrigger(HLTPath_); trigger.HLTobjects(HLTPath_);
      if (saveOnlyHLTFires && !trigger.HLTPathFire()) return;
      trigger.FillHLT(nt); trigger.FillObj(nt);
      if (trigger.HLTPathFire())
        TrgMu_EtaPhi=std::make_pair(trigger.GetHighestPtHLTObject()[1],trigger.GetHighestPtHLTObject()[2]);
    }
   } 
  //  cout<<"here trg"<<endl;
  nevts++;
  for (const pat::Muon &mu : *muons){    
    // for (unsigned int i=0, n=mu.numberOfSourceCandidatePtrs(); i<n; ++i){
     //   footprint.push_back(mu.sourceCandidatePtr(i));
   // }
    nt.muon_pt.push_back(mu.pt()); nt.muon_phi.push_back(mu.phi());
    nt.muon_eta.push_back(mu.eta()); nt.muon_charge.push_back(mu.charge());
    nt.muon_global_flag.push_back(mu.isGlobalMuon());
    nt.muon_standalone_flag.push_back(mu.isStandAloneMuon());
    nt.muon_tracker_flag.push_back(mu.isTrackerMuon());
    nt.muon_vx.push_back(mu.vx()); nt.muon_vy.push_back(mu.vy());
    nt.muon_vz.push_back(mu.vz()); nt.muon_edz.push_back(mu.dzError());
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
    if (DR(TrgMu_EtaPhi.first,TrgMu_EtaPhi.second,mu.eta(),mu.phi())<DRtrgMu){
      DRtrgMu=DR(TrgMu_EtaPhi.first,TrgMu_EtaPhi.second,mu.eta(),mu.phi());
      TrgmuDz=mu.vz();  nt.muon_trgIndex=nt.nmuons;
    }
    const reco::MuonPFIsolation&  isol=mu.pfIsolationR04();
    double mu_iso=(isol.sumChargedHadronPt+max(0.,isol.sumNeutralHadronEt+isol.sumPhotonEt-0.5*isol.sumPUPt))/mu.pt();
    nt.muon_iso.push_back(mu_iso);
    muttks.emplace_back(reco::TransientTrack(*mu.bestTrack(),&(*bFieldHandle)));
    nt.nmuons++;
    }
    if (DRtrgMu>TrgConeCut && TrgConeCut>0 && saveOnlyHLTFires) return;
    
    //electrons  
    trigger::size_type eindex=-1; 
    for(const pat::Electron &el : *electrons){
      eindex++;
   //   for (unsigned int i=0, n=el.numberOfSourceCandidatePtrs(); i<n; ++i){
  //      footprint.push_back(el.sourceCandidatePtr(i));
//      }
      float pt=0,eta=-30,phi=-30;
      if (fabs(TrgmuDz-el.vz())>ElectronDzCut) continue;
      if (!el.passConversionVeto()) continue;
      if (IsLowpTE){
        reco::GsfTrackRef seed = el.gsfTrack();
        if ( seed.isNull() ) continue;
        pt=el.gsfTrack()->ptMode(); eta=el.gsfTrack()->etaMode();
        phi=el.gsfTrack()->phiMode();
        if (pt<Electron2PtCut) continue;
        if ( (*ele_mva_wp_unbiased)[seed]<MVAEl2Cut) continue;
        nt.el_mva_biased.push_back((*ele_mva_wp_biased)[seed]);
        nt.el_mva_unbiased.push_back((*ele_mva_wp_unbiased)[seed]);   
      } else{
        pt=el.pt(); eta=el.eta(); phi=el.phi();
        if (pt<Electron2PtCut) continue;
        nt.el_mva_out.push_back(el.mva_e_pi());
        nt.el_mva_iso.push_back(el.mva_Isolated());
        const edm::Ptr<pat::Electron> elePtr(electrons,eindex);
        nt.el_veto.push_back((*ele_veto_id)[elePtr]);
        nt.el_soft.push_back((*ele_soft_id)[elePtr]);
        nt.el_medium.push_back((*ele_medium_id)[elePtr]);
        nt.el_tight.push_back((*ele_tight_id)[elePtr]);
        nt.el_mva_map_value.push_back((*ele_mva_id_value)[elePtr]);
        nt.el_mva_biased.push_back(-10); nt.el_mva_unbiased.push_back(-10); 
      }
      nt.el_pt.push_back(pt); nt.el_eta.push_back(eta);
      nt.el_phi.push_back(phi); nt.el_charge.push_back(el.charge());
      nt.el_dz.emplace_back(el.bestTrack()->dz(vertex_point));
      nt.el_dxy.emplace_back(el.bestTrack()->dxy(vertex_point));
      nt.el_edxy.push_back(el.dxyError()); nt.el_edz.push_back(el.dzError());
      nt.el_vx.push_back(el.vx()); nt.el_vy.push_back(el.vy());
      nt.el_vz.push_back(el.vz());
      double iso= el.pfIsolationVariables().sumChargedHadronPt+max(0.0,el.pfIsolationVariables().sumNeutralHadronEt+el.pfIsolationVariables().sumPhotonEt-0.5*el.pfIsolationVariables().sumPUPt)/pt;
      nt.el_iso.push_back(iso); nt.el_islowpt.push_back(IsLowpTE);
      nt.el_trkpt.push_back(el.bestTrack()->pt()); nt.el_trketa.push_back(el.bestTrack()->eta());
      nt.el_trkphi.push_back(el.bestTrack()->phi()); 
      ettks.emplace_back(reco::TransientTrack(*el.bestTrack(),&(*bFieldHandle)));
      nt.nel++;
      if (CombineElCol)  PFe_EtaPhi.emplace_back(eta,phi);
    }
   
    //combine with low pT
    if (CombineElCol){
     for(const pat::Electron &el : *lowpte){
       float pt=0,eta=-30,phi=-30;
       if (fabs(TrgmuDz-el.vz())>ElectronDzCut) continue;
       if (!el.passConversionVeto()) continue;
       reco::GsfTrackRef seed = el.gsfTrack();
       if ( seed.isNull() ) continue;
       pt=el.gsfTrack()->ptMode(); eta=el.gsfTrack()->etaMode();
       phi=el.gsfTrack()->phiMode();
       if (pt<Electron2PtCut) continue;
       if ( (*ele_mva_wp_unbiased)[seed]<MVAEl2Cut) continue;
       bool Elsaved=false;
       for (std::pair<float,float> & pfe : PFe_EtaPhi){
         if (DR(pfe.first,pfe.second,eta,phi)<CombineCone)
           Elsaved=true;
       }
       if (Elsaved) continue;
       nt.el_mva_biased.emplace_back((*ele_mva_wp_biased)[seed]);
       nt.el_mva_unbiased.emplace_back((*ele_mva_wp_unbiased)[seed]);         
       nt.el_pt.push_back(pt); nt.el_eta.push_back(eta);
       nt.el_phi.push_back(phi); nt.el_charge.push_back(el.charge());
       nt.el_dz.emplace_back(el.bestTrack()->dz(vertex_point));
       nt.el_dxy.emplace_back(el.bestTrack()->dxy(vertex_point));
       nt.el_edxy.push_back(el.dxyError()); nt.el_edz.push_back(el.dzError());
       nt.el_vx.push_back(el.vx()); nt.el_vy.push_back(el.vy());
       nt.el_vz.push_back(el.vz());
       double iso= el.pfIsolationVariables().sumChargedHadronPt+max(0.0,el.pfIsolationVariables().sumNeutralHadronEt+el.pfIsolationVariables().sumPhotonEt-0.5*el.pfIsolationVariables().sumPUPt)/pt;
       nt.el_iso.push_back(iso); nt.el_islowpt.push_back(true);
       nt.el_trkpt.emplace_back(el.bestTrack()->pt());
       nt.el_trketa.emplace_back(el.bestTrack()->eta());
       nt.el_trkphi.emplace_back(el.bestTrack()->phi()); 
       ettks.emplace_back(reco::TransientTrack(*el.bestTrack(),&(*bFieldHandle)));
       nt.nel++;
      }
    }    
    //additional low pt el for combination
    if (AddLowPtElAsCol && AddLowPtGsfTrkAsCol){
      LowPtElObjects lowelpt(lowPtElectronsToken_,lowPtGsfTracksToken_,eleUnBWPToken_,iEvent); 
      lowelpt.AddElectrons(nt);  lowelpt.AddGsfTracks(nt);
    } else if (AddLowPtElAsCol){
      LowPtElObjects lowelpt(lowPtElectronsToken_,eleUnBWPToken_,iEvent); 
      lowelpt.AddElectrons(nt);
    } else if (AddLowPtGsfTrkAsCol){
      LowPtElObjects lowelpt(lowPtElectronsToken_,lowPtGsfTracksToken_,eleUnBWPToken_,iEvent); 
      lowelpt.AddGsfTracks(nt);
    }
    //add PF as collection
   if (AddPFElAsCol){
      PFelCollection PFelectrons(pfElectronsToken_,eleIdMapSoftToken_,eleIdMapMediumToken_,eleIdMapTightToken_,iEvent);
      PFelectrons.AddElectrons(nt);
   }
    /*  std::vector<std::vector<float>> mu_DCA=track_DCA(muttks);
  if (mu_DCA.size()==0) { std::vector<float> d1; d1.push_back(-99); muon_DCA.push_back(d1); }
  else muon_DCA=mu_DCA;*/
   //save tracks
  for (const pat::PackedCandidate &trk : *tracks1) {
   if (!trk.trackHighPurity()) continue;
   tracks.push_back(trk); }
  
  if (AddLostTracks){
   for (const pat::PackedCandidate &trk : *tracks2) {
    tracks.push_back(trk); }}

 if (saveTracks){
  for (const pat::PackedCandidate &trk : tracks){      
   if(trk.charge()==0) continue;
   if(fabs(trk.pdgId())!=211) continue;
   if(!trk.hasTrackDetails())continue;
   if (trk.pt()< track_pt_cut_forB) continue;
   if (fabs(trk.eta())>EtaTrk_Cut) continue;
   if (fabs(TrgmuDz-trk.vz())>ElectronDzCut) continue;
   if (UseDirectlyGenBeeK){
     if (DR(EtaPhiK.first,EtaPhiK.second,trk.eta(),trk.phi())>DRgenCone)
       continue;
    } 
   bool isMu=false; bool isE=false;
   for (const pat::Muon & mu :*muons){
     if (DR(mu.eta(),mu.phi(),trk.eta(),trk.phi())<LepTrkExclusionCone) isMu=true; 
   }
   if (CombineElCol || !IsLowpTE){
    for (const pat::Electron & el : *electrons) {     
      if (DR(el.eta(),el.phi(),trk.eta(),trk.phi())<LepTrkExclusionCone)
          isE=true;
     }}
   if(isMu || isE ) continue;
   nt.track_pt.push_back(trk.pt()); nt.track_eta.push_back(trk.eta());
   nt.track_phi.push_back(trk.phi()); nt.track_charge.push_back(trk.charge());
   if(trk.trackHighPurity()) nt.track_highPurity.push_back(1);
   else nt.track_highPurity.push_back(0);
   nt.track_norm_chi2.emplace_back(trk.pseudoTrack().normalizedChi2());
   nt.track_dxy.push_back(trk.dxy(vertex_point));
   nt.track_dz.push_back(trk.dz(vertex_point));
   nt.track_validhits.push_back(trk.numberOfHits());
   nt.track_losthits.push_back(trk.lostInnerHits());
   nt.track_fromPV.push_back(trk.fromPV()); nt.ntracks++;  
  }
 }//save tracks

 if(!reconstructBMuMuK && !reconstructBMuMuKstar){ t1->Fill(); return; }  
 //muon pairs 
 if(!OnlyKee){
  for (int imu=0; imu<nt.nmuons; ++imu){
    if (!nt.muon_soft[imu]) continue;
    if(nt.muon_pt[imu]< min_muon_pt_cut_forB) continue;
    for (int imu2=imu+1; imu2<nt.nmuons; ++imu2){
     if (!nt.muon_soft[imu2]) continue; 
     if(nt.muon_pt[imu2]< min_muon_pt_cut_forB) continue;
     if (nt.muon_charge[imu]==nt.muon_charge[imu2]) continue;
     if (nt.muon_pt[imu]< max_muon_pt_cut_forB && nt.muon_pt[imu2]< max_muon_pt_cut_forB)
         continue;
     TLorentzVector vmu1,vmu2; 
     vmu1.SetPtEtaPhiM(nt.muon_pt[imu],nt.muon_eta[imu],nt.muon_phi[imu],0.105);
     vmu2.SetPtEtaPhiM(nt.muon_pt[imu2],nt.muon_eta[imu2],nt.muon_phi[imu2],0.105);
     if((vmu1+vmu2).M()<MLLmin_Cut || (vmu1+vmu2).M()>MLLmax_Cut) continue;
     muTrack1.push_back(muttks[imu]); muTrack2.push_back(muttks[imu2]);
     used_muTrack_index.emplace_back(imu,imu2);
     nmupairs++;
    }
  }
 } 
 
 //e pairs
 if(AddeeK || OnlyKee ){
  for(int iel=0; iel<nt.nel; ++iel){
    if (UseDirectlyGenBeeK && !data){
      if (DR(EtaPhiE1.first,EtaPhiE1.second,nt.el_eta[iel],nt.el_phi[iel])>DRgenCone && DR(EtaPhiE2.first,EtaPhiE2.second,nt.el_eta[iel],nt.el_phi[iel])>DRgenCone) continue;
    }
    for (int iel2=iel+1; iel2<nt.nel; ++iel2){
      if (nt.el_charge[iel]==nt.el_charge[iel2]) continue;
      if (nt.el_pt[iel2]<Electron1PtCut && nt.el_pt[iel]<Electron1PtCut) continue;
      if (UseDirectlyGenBeeK && !data){
        if (DR(EtaPhiE1.first,EtaPhiE1.second,nt.el_eta[iel2],nt.el_phi[iel2])>DRgenCone && DR(EtaPhiE2.first,EtaPhiE2.second,nt.el_eta[iel2],nt.el_phi[iel2])>DRgenCone) continue;
      }
       if (nt.el_islowpt[iel] && nt.el_islowpt[iel2]){
         if (nt.el_mva_unbiased[iel]<MVAEl1Cut && nt.el_mva_unbiased[iel2]<MVAEl1Cut) continue;
       }
       if(fabs(nt.el_vz[iel]-nt.el_vz[iel2])>DzeeMaxCut) continue;
       TLorentzVector vel1,vel2; 
       vel1.SetPtEtaPhiM(nt.el_pt[iel],nt.el_eta[iel],nt.el_phi[iel],0.000511);
       vel2.SetPtEtaPhiM(nt.el_pt[iel2],nt.el_eta[iel2],nt.el_phi[iel2],0.000511);
       if((vel1+vel2).M()<MLLmin_Cut || (vel1+vel2).M()>MLLmax_Cut) continue; 
       muTrack1.push_back(ettks[iel]); muTrack2.push_back(ettks[iel2]);
       used_eTrack_index.emplace_back(iel,iel2);
    }
  }
 }
 /*cout<<"size "<<used_eTrack_index.size()<<endl;
 for (unsigned int i=0; i<used_eTrack_index.size(); i++)
 cout<<"fr "<<used_eTrack_index[i].first<<" sec "<<used_eTrack_index[i].second<<endl;*/
 if (used_muTrack_index.size()==0 && used_eTrack_index.size()==0){
   if (reconstructBMuMuK && SkipEventWithNoBToMuMuK && !reconstructBMuMuKstar) return;
   else if (reconstructBMuMuKstar && SkipEventWithNoBToMuMuKstar && !reconstructBMuMuK) return;
   else if (reconstructBMuMuK && SkipEventWithNoBToMuMuK && SkipEventWithNoBToMuMuKstar && reconstructBMuMuKstar) return;
}
 
if (RetrieveMuFromTrk){
   nmupfpairs=nmupairs;
   for (const pat::Muon & mu : *muons){
     if (fabs(TrgmuDz-mu.vz())>ElectronDzCut) continue;
     if (mu.pt()<min_muon_pt_cut_forB) continue;
     if (fabs(mu.eta())>2.5) continue;
     if (!mu.isSoftMuon(firstGoodVertex)) continue;
     for(const pat::PackedCandidate &trk: tracks){
       if(mu.charge()==trk.charge()) continue;
       if(trk.charge()==0) continue;
       if(fabs(trk.pdgId())!=211) continue;
       if(!trk.hasTrackDetails())continue;
       if (trk.pt()< min_muon_pt_cut_forB) continue;
       if (fabs(trk.eta())>EtaTrk_Cut) continue;
       if (fabs(TrgmuDz-trk.vz())>ElectronDzCut) continue;
       if(trk.pt()>maxPtTrk) continue;
       TLorentzVector vmu1,vmu2;
       vmu1.SetPtEtaPhiM(mu.pt(),mu.eta(),mu.phi(),0.105);
       vmu1.SetPtEtaPhiM(trk.pt(),trk.eta(),trk.phi(),0.105);
       if ((vmu1+vmu2).M()<MLLmin_Cut || (vmu1+vmu2).M()>MLLmax_Cut) continue;
       muTrack1.emplace_back(reco::TransientTrack(*mu.bestTrack(),&(*bFieldHandle)));
       muTrack2.emplace_back(reco::TransientTrack(trk.pseudoTrack(),&(*bFieldHandle)));
       used_muTrack_pfTrack_index.emplace_back(std::make_pair(&mu-&(muons->at(0)),&trk-&tracks[0]));      
       nmupfpairs++;
     }
    }
}
int index=-1;
if((used_muTrack_index.size()>0 || used_eTrack_index.size()>0) && (reconstructBMuMuK || reconstructBMuMuKstar)){
 for (const pat::PackedCandidate & trk: tracks){
   index++;
   if(trk.charge()==0) continue;
   if(fabs(trk.pdgId())!=211) continue;
   if(!trk.hasTrackDetails())continue;
   reco::Track ptrk=trk.pseudoTrack();
   if (ptrk.pt()< track_pt_cut_forB) continue;
   if (fabs(ptrk.eta())>EtaTrk_Cut) continue;
   if (UseDirectlyGenBeeK){
      if (DR(EtaPhiK.first,EtaPhiK.second,ptrk.eta(),ptrk.phi())>DRgenCone)
	continue;
    }
   if (fabs(TrgmuDz-trk.vz())>ElectronDzCut) continue;
   bool isMu=false; bool isE=false;
   for (const pat::Muon & mu : *muons){
     if (DR(mu.eta(),mu.phi(),ptrk.eta(),ptrk.phi())<LepTrkExclusionCone) 
       isMu=true; 
   }
   if (!IsLowpTE || CombineElCol){
     for (const pat::Electron & el: *electrons) {     
        if (DR(el.eta(),el.phi(),ptrk.eta(),ptrk.phi())<LepTrkExclusionCone) 
           isE=true;
   } }
   if(isMu || isE ) continue;
   KTrack.emplace_back(reco::TransientTrack(ptrk,&(*bFieldHandle)));
   KTrack_index.push_back(&trk-&tracks[0]);  
   }
 }



//add pf
  if (used_muTrack_pfTrack_index.size()>0)
    used_muTrack_index.insert(used_muTrack_index.end(),used_muTrack_pfTrack_index.begin(),used_muTrack_pfTrack_index.end());
 //add ee chanel 
 if (used_eTrack_index.size()>0)
   used_muTrack_index.insert(used_muTrack_index.end(),used_eTrack_index.begin(),used_eTrack_index.end());
 
 ///building B 
 if (used_muTrack_index.size()>0 && KTrack_index.size()>0 && reconstructBMuMuK){

    BKlldecay Kll(nmupairs,used_muTrack_index,muTrack1,muTrack2,KTrack,beam_xz,beam_yz,RefitTracks);
    Kll.ProbCut(Pchi2BMuMuK); Kll.CosCut(CosThetaCut); 
    Kll.MassCuts(MBmin_Cut,MBmax_Cut); Kll.PtBCut(PtBminCut); 
    Kll.SetNumTrkForMu(nmupfpairs); 
    Kll.Fill(nt);  
 } 

 if (used_muTrack_index.size()>0 && KTrack_index.size()>0 && reconstructBMuMuKstar){
    BKstarlldecay Kstarll(nmupairs,used_muTrack_index,muTrack1,muTrack2,KTrack,beam_xz,beam_yz,RefitTracks);
    Kstarll.CombineTracks(0.600,1.100);
    Kstarll.Fill(nt);
 } 

 if (used_muTrack_index.size()>0 && KTrack_index.size()>0 && reconstructBMuMuKstar){
   BKstarlldecay phill(nmupairs,used_muTrack_index,muTrack1,muTrack2,KTrack,beam_xz,beam_yz,RefitTracks);
    phill.Runphill(true); phill.CombineTracks(0.600,1.100);
    phill.FillPhill(nt);
  }

  //skip useless events
  if (nt.NRb_mass.size()==0 && nt.NRbks_mass.size()==0){
   if (reconstructBMuMuK && SkipEventWithNoBToMuMuK && !reconstructBMuMuKstar) return;
   else if (reconstructBMuMuKstar && SkipEventWithNoBToMuMuKstar && !reconstructBMuMuK) return;
   else if (reconstructBMuMuK && SkipEventWithNoBToMuMuK && SkipEventWithNoBToMuMuKstar && reconstructBMuMuKstar) return;
}  
if (NtupleOutputClasses=="flat") nt.Flatting();
  t1->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
template<typename T1>
void 
TriggerAnalyzer<T1>::beginJob()
{
  t1=fs->make<TTree>("mytree","mytree");
  
  nt.SetTree(t1);
  TString ToSave="";
  if (NtupleOutputClasses=="auto" || NtupleOutputClasses=="flat" || NtupleOutputClasses=="lite" ){
    if (NtupleOutputClasses=="flat") ToSave+="Flat_";
    if (NtupleOutputClasses=="lite") ToSave+="Lite_";
    if (reconstructBMuMuK) ToSave+="KLL_";
    if (reconstructBMuMuKstar) ToSave+="KstarLL_";
    if (AddLowPtElAsCol) ToSave+="LowPtEl_";
    if (AddLowPtGsfTrkAsCol) ToSave+="LowPtGsf_";
    if (!data) ToSave+="GEN_";
    if (saveTracks) ToSave+="TRK_";
  } else ToSave=NtupleOutputClasses;

  nt.SetNtupleVariables(ToSave);

}

// ------------ method called once each job just after ending the event loop  ------------
template<typename T1>
void 
TriggerAnalyzer<T1>::endJob() 
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

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template<typename T1>
void
TriggerAnalyzer<T1>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
typedef TriggerAnalyzer<reco::RecoEcalCandidate> TriggerAnalyzerb;
DEFINE_FWK_MODULE(TriggerAnalyzerb);
