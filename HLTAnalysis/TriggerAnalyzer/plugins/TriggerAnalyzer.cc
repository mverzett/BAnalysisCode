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
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "HLTrigger/HLTcore/interface/defaultModuleLabel.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/L1Trigger/interface/EtSumHelper.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateIsolation.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "L1Trigger/L1TNtuples/interface/MuonID.h"
#include <vector>
#include "TTree.h"
#include <string>
#include <iostream>
#include "TMath.h"
#include <cmath>
#include "DataFormats/Common/interface/Ref.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/CombinedKinematicConstraint.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackPointingKinematicConstraint.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1TUtmAlgorithm.h"
#include "CondFormats/L1TObjects/interface/L1TUtmTriggerMenu.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1TGlobal/interface/GlobalExtBlk.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "TripleTrackKinFit.h"
#include "GeneratorBTree.h"
#include "NtupleContent.h"
#include "LowPtElObjects.h"
#include "PFelCollection.h"
#include "HLTL1tree.h"
/*namespace edm {
  class ConfigurationDescriptions;
  }*/



using namespace std;
//using namespace edm;


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
  std::pair<std::vector<float>,std::vector<std::vector<float>>> L1Analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  std::vector<std::pair<string,int>> L1SeedAnalyze(const edm::Event& iEvent,TString * algoBitToName, std::vector<string> Seed);
  std::pair<std::vector<float>,std::vector<std::vector<std::vector<float>>>> HLTAnalyze(const edm::Event& iEvent, const edm::EventSetup& iSetup,std::vector<string> HLTPath);
  std::vector<std::vector<float>> track_DCA(std::vector<reco::TransientTrack> ttks);
  std::vector<GlobalVector>refit_tracks(TransientVertex myVertex,std::vector<reco::TransientTrack> tracks);
//  static bool findmax(const vector<float> &a,const vector<float> &b);
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
  edm::EDGetToken photonToken_;
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

 
//  edm::ParameterSet const& conf;
 

  edm::Service<TFileService> fs;
  TTree * t1;
  NtupleContent nt;
  int good_vertex=0;
  int totb1=0,totb2=0,totb3=0;
  int trg_counter=0,fire_counter=0,l3_counter=0,ntracks=0,nkaons=0;
  std::vector<std::vector<int>>muon_vtx_IsValid,el_vtx_IsValid;
  int npv=0;
 
  ///options
   bool data=true; bool saveTracks=true; bool saveKshort=true;
  bool saveHLT=true; bool saveL1=true; bool saveOnlyHLTFires=false;
  double track_pt_cut_forB=0; double muon_pt_cut_forB=0;
  bool reconstructBMuMuK=true; bool Pointing_constraint=false;
  bool reconstructBMuMuKstar=true;
  double Pchi2BMuMuK=-1; double MLLmax_Cut=100; double MLLmin_Cut=-1;
  bool SkipEventWithNoBToMuMuK=false; bool UseBeamspot=false; bool AddeeK=false;
  double MBmax_Cut=100; double MBmin_Cut=-1;
  double  LepTrkExclusionCone=-1; double EtaTrk_Cut=5; bool AddLostTracks=true;
  double MKstarMin_Cut=0.5; double MKstarMax_Cut=1.5; bool RefitTracks=true;
  bool RefitMuTracksOnly=false; bool UsePFeForCos=true; bool OnlyKee=false;
  bool SkipEventWithNoBToMuMuKstar=false; double Electron1PtCut=0;
  double Electron2PtCut=0; double ElectronDzCut=0; double TrgConeCut=-1;
  bool IsLowpTE=false; double MVAEl1Cut=-20; double MVAEl2Cut=-20;
  double CosThetaCut=-1; bool UseDirectlyGenBeeK=false; double DRgenCone=10;
  int KIdToMatch=-1; int LepIdToMatch=-1; int BpdgIdToMatch=-1;
  bool IsResonantDecayToMatch=false; bool AddLowPtElAsCol=false;
  bool AddLowPtGsfTrkAsCol=false; bool AddPFElAsCol=false;
  std::string NtupleOutputClasses="auto";
 //internal
 ParticleMass part_mass = 0.1056583; float part_sigma = 0.0000001;
 ParticleMass kaon_mass = 0.493677; float kaon_sigma = 0.000016;
 float chi = 0.; float ndf = 0.;
 std::pair<std::vector<float>,std::vector<std::vector<float>>> l1objects;
 std::vector<std::pair<string,int>> l1seeds;
 std::vector<std::shared_ptr<reco::TransientTrack>> KTrack;
 std::vector<unsigned int> KTrack_index;
 std::vector<std::pair<std::shared_ptr<reco::TransientTrack>,std::shared_ptr<reco::TransientTrack>>> KstarTrack;
 std::vector<std::pair<unsigned int,unsigned int>> KstarTrack_index;
 unsigned int nmupairs=0;
std::vector<std::shared_ptr<reco::TransientTrack>> muTrack1,muTrack2;//,eTrack1,eTrack2;
std::vector<std::pair<unsigned int,unsigned int>> used_muTrack_index,used_eTrack_index;
 std::vector<reco::CandidatePtr> footprint;
 std::vector<pat::PackedCandidate> tracks;


 TString * algoBitToName = new TString[512];
 int count=0;
      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
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
 photonToken_(consumes<std::vector<pat::Photon>>(iConfig.getParameter<edm::InputTag>("photons"))),
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

 muon_pt_cut_forB=runParameters.getParameter<double>("MuonPtCutForB");
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
 RefitTracks=runParameters.getParameter<bool>("RefitTracks");
 RefitMuTracksOnly=runParameters.getParameter<bool>("RefitMuTracksOnly");
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
}

template<typename T1>
TriggerAnalyzer<T1>::~TriggerAnalyzer()
{
  // cout<<"total "<<trg_counter<<" fires "<<fire_counter<<" l3 "<<l3_counter<<endl;
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  cout<<"total trg2(Mu17) "<<totb1<<" trg3(Mu20) "<<totb2<<" trg5(Mu27) "<<totb3<<endl;
  
  
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
  using namespace edm;
  using namespace reco;
  using namespace trigger;


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
  edm::Handle<std::vector<pat::Jet>> jets;
  iEvent.getByToken(jetsToken_, jets);
  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonsToken_,muons);
  edm::Handle<std::vector<pat::MET>> met;
  iEvent.getByToken(metToken_,met);
  edm::Handle<std::vector<pat::Photon>> photons;
  iEvent.getByToken(photonToken_,photons);
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

  footprint.clear(); nt.ClearVariables();
 
  nt.event=iEvent.id().event();
  nt.run_number=iEvent.id().run();
  nt.ls=iEvent.luminosityBlock();
  nt.beam_x=theBeamSpot->x0(); nt.beam_y=theBeamSpot->y0(); 
  nt.beam_z=theBeamSpot->z0(); nt.beam_ex= theBeamSpot->x0Error(); 
  nt.beam_ey= theBeamSpot->y0Error(); nt.beam_ez= theBeamSpot->z0Error();
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

   reco::TrackBase::Point  vertex_point; vertex_point.SetCoordinates(nt.pvertex_x,nt.pvertex_y,nt.pvertex_z);
   //K0 fit kalman
   KalmanVertexFitter theKalmanFitter(false);
   TransientVertex K0vertex;
   std::pair<float,float> EtaPhiE1(-10,-10),EtaPhiE2(-10,-10),EtaPhiK(-10,-10);
   if(!data){
     GeneratorBTree gen(prunedGenToken_,packedGenToken_,iEvent);
     nt.ngenB=gen.BMother_Pt().size();
     nt.genB_pt=gen.BMother_Pt(); nt.genB_phi=gen.BMother_Phi();
     nt.genB_eta=gen.BMother_Eta(); nt.genB_pdgId=gen.BMother_PdgId(); 
     nt.genB_Bindex=gen.BMother_Bid();
     nt.genB_daughter_pt=gen.BDaughter_Pt(); nt.genB_daughter_eta=gen.BDaughter_Eta();
     nt.genB_daughter_phi=gen.BDaughter_Phi(); nt.genB_daughter_pdgId=gen.BDaughter_PdgId(); 
     nt.genB_daughter_Bindex=gen.BDaughter_Bid(); nt.genB_daughter_Dindex=gen.BDaughter_Did(); 
     nt.genB_granddaughter_pt=gen.BGDaughter_Pt(); nt.genB_granddaughter_eta=gen.BGDaughter_Eta();
     nt.genB_granddaughter_phi=gen.BGDaughter_Phi(); nt.genB_granddaughter_pdgId=gen.BGDaughter_PdgId();
     nt.genB_granddaughter_Bindex=gen.BGDaughter_Bid();
     nt.genB_granddaughter_Dindex=gen.BGDaughter_Did();
     nt.ngenLep=gen.Lep_Pt().size();     
     nt.genLep_pt=gen.Lep_Pt(); nt.genLep_phi=gen.Lep_Phi();
     nt.genLep_eta=gen.Lep_Eta(); nt.genLep_pdgId=gen.Lep_pdgId(); 
     nt.genLep_mom=gen.Lep_Mother();
     if (UseDirectlyGenBeeK){
      gen.GenKLepLepEtaPhi(IsResonantDecayToMatch,BpdgIdToMatch,LepIdToMatch,KIdToMatch);
      EtaPhiK=gen.EtaPhiK(); EtaPhiE1=gen.EtaPhiL(); EtaPhiE2=gen.EtaPhiaL();
     }
   }

  //find max
  std::pair<float,float> TrgMu_EtaPhi(-100,-100);
  if (saveHLT || saveL1){
    HLTL1tree trigger(iEvent, l1resultToken_, l1MuonsToken_, l1JetsToken_, l1MetToken_, trgresultsToken_, trigobjectsToken_);
    if (saveL1){ 
      trigger.L1objetcs(nt);
      trigger.L1trigger(algoBitToName,Seed_);
      nt.l1_seed1=(trigger.GetL1Decision()[0]) ? 1 :0;
      nt.l1_seed2=(trigger.GetL1Decision()[1]) ? 1 :0;
      nt.l1_seed3=(trigger.GetL1Decision()[2]) ? 1 :0;
      nt.l1_seed4=(trigger.GetL1Decision()[3]) ? 1 :0;
      nt.l1_seed5=(trigger.GetL1Decision()[4]) ? 1 :0;
      nt.l1_seed6=(trigger.GetL1Decision()[5]) ? 1 :0;
    }
    if (saveHLT){
      trigger.HLTtrigger(HLTPath_); trigger.HLTobjects(HLTPath_);
      if (saveOnlyHLTFires && !trigger.HLTPathFire()) return;
      nt.trigger1=(trigger.GetHLTDecision()[0]) ? 1:0;
      nt.trigger2=(trigger.GetHLTDecision()[1]) ? 1:0;
      nt.trigger3=(trigger.GetHLTDecision()[2]) ? 1:0;
      nt.trigger4=(trigger.GetHLTDecision()[3]) ? 1:0;
      nt.trigger5=(trigger.GetHLTDecision()[4]) ? 1:0;
      nt.trigger6=(trigger.GetHLTDecision()[5]) ? 1:0;
      nt.trigger7=(trigger.GetHLTDecision()[6]) ? 1:0;
      nt.trigger8=(trigger.GetHLTDecision()[7]) ? 1:0;
      if (nt.trigger1==1)
        nt.tr1_obj_pt_eta_phi=trigger.GetHLTObjects(0);
      if (nt.trigger2==1)
        nt.tr2_obj_pt_eta_phi=trigger.GetHLTObjects(1);
      if (nt.trigger3==1)
        nt.tr3_obj_pt_eta_phi=trigger.GetHLTObjects(2);
      if (nt.trigger4==1)
        nt.tr4_obj_pt_eta_phi=trigger.GetHLTObjects(3);
      if (nt.trigger5==1)
        nt.tr5_obj_pt_eta_phi=trigger.GetHLTObjects(4);
      if (nt.trigger6==1)
        nt.tr6_obj_pt_eta_phi=trigger.GetHLTObjects(5);
      if (nt.trigger7==1)
        nt.tr7_obj_pt_eta_phi=trigger.GetHLTObjects(6);
      if (nt.trigger8==1)
        nt.tr8_obj_pt_eta_phi=trigger.GetHLTObjects(7);
      if (trigger.HLTPathFire())
        TrgMu_EtaPhi=std::make_pair(trigger.GetHighestPtHLTObject()[1],trigger.GetHighestPtHLTObject()[2]);
    }
   }
   // totb1+=nt.trigger2; totb2+=nt.trigger3; totb3+=nt.trigger5;
   // ef (TrgMu_EtaPhi.second<-10 && saveHLT) return;

 
  std::vector<reco::TransientTrack> muttks,ettks,jpsimuttks,jpsiettks;
  std::vector<int> muindex,elindex,mucharge,elcharge;
  float TrgmuDz=0,DRtrgMu=100;
  for (const pat::Muon &mu : *muons){    
     for (unsigned int i=0, n=mu.numberOfSourceCandidatePtrs(); i<n; ++i){
        footprint.push_back(mu.sourceCandidatePtr(i));
    }
    nt.muon_pt.push_back(mu.pt()); nt.muon_phi.push_back(mu.phi());
    nt.muon_eta.push_back(mu.eta()); nt.muon_charge.push_back(mu.charge());
    const Track * mutrack= mu.bestTrack();
    nt.muon_dz.push_back(mutrack->dz(vertex_point));
    nt.muon_dxy.push_back(mutrack->dxy(vertex_point));
    nt.muon_global_flag.push_back(mu.isGlobalMuon());
    nt.muon_standalone_flag.push_back(mu.isStandAloneMuon());
    nt.muon_tracker_flag.push_back(mu.isTrackerMuon());
    nt.muon_vx.push_back(mu.vx()); nt.muon_vy.push_back(mu.vy());
    nt.muon_vz.push_back(mu.vz()); nt.muon_edz.push_back(mu.dzError());
    nt.muon_edxy.push_back(mu.dxyError()); nt.muon_d0.push_back(mutrack->d0());
    nt.muon_ed0.push_back(mutrack->d0Error());
    nt.muon_medium.push_back(mu.isMediumMuon());
    nt.muon_loose.push_back(mu.isLooseMuon());
    nt.muon_trkpt.push_back(mutrack->pt()); nt.muon_trketa.push_back(mutrack->eta());
    nt.muon_trkphi.push_back(mutrack->phi());
    nt.muon_tight.push_back(mu.isTightMuon(firstGoodVertex));
    nt.muon_soft.push_back(mu.isSoftMuon(firstGoodVertex));
    if (DR(TrgMu_EtaPhi.first,TrgMu_EtaPhi.second,mu.eta(),mu.phi())<DRtrgMu){
      DRtrgMu=DR(TrgMu_EtaPhi.first,TrgMu_EtaPhi.second,mu.eta(),mu.phi());
      TrgmuDz=mu.vz();  nt.muon_trgIndex=nt.nmuons;
    }
   // if (mu.isSoftMuon(firstGoodVertex))snmu++;
    const MuonPFIsolation&  isol=mu.pfIsolationR04();
    double mu_iso=(isol.sumChargedHadronPt+max(0.,isol.sumNeutralHadronEt+isol.sumPhotonEt-0.5*isol.sumPUPt))/mu.pt();
    nt.muon_iso.push_back(mu_iso);
    muttks.push_back(reco::TransientTrack(*mutrack,&(*bFieldHandle)));   
    nt.nmuons++;
    }
    if (DRtrgMu>TrgConeCut && TrgConeCut>0 && saveOnlyHLTFires) return;
    
      
    size_type eindex=-1;
    for(const pat::Electron &el : *electrons){
      eindex++;
      for (unsigned int i=0, n=el.numberOfSourceCandidatePtrs(); i<n; ++i){
        footprint.push_back(el.sourceCandidatePtr(i));
      }
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
      }
      nt.el_pt.push_back(pt); nt.el_eta.push_back(eta);
      nt.el_phi.push_back(phi); nt.el_charge.push_back(el.charge());
      const Track * eltrack= el.bestTrack();
      nt.el_dz.push_back(eltrack->dz(vertex_point));
      nt.el_dxy.push_back(eltrack->dxy(vertex_point));
      nt.el_edxy.push_back(el.dxyError()); nt.el_edz.push_back(el.dzError());
      nt.el_vx.push_back(el.vx()); nt.el_vy.push_back(el.vy());
      nt.el_vz.push_back(el.vz());
      double iso= el.pfIsolationVariables().sumChargedHadronPt+max(0.0,el.pfIsolationVariables().sumNeutralHadronEt+el.pfIsolationVariables().sumPhotonEt-0.5*el.pfIsolationVariables().sumPUPt)/el.pt();
      nt.el_iso.push_back(iso); nt.el_islowpt.push_back(IsLowpTE);
      nt.el_trkpt.push_back(eltrack->pt()); nt.el_trketa.push_back(eltrack->eta());
      nt.el_trkphi.push_back(eltrack->phi()); 
      ettks.push_back(reco::TransientTrack(*eltrack,&(*bFieldHandle)));
      nt.nel++;
    }
    //additional low pt el
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
    //  if ( nt.el_pt.size()<2 && OnlyKee){ return;}
  tracks.clear();  
  for (const pat::PackedCandidate &trk : *tracks1) {
   if (!trk.trackHighPurity()) continue;
   tracks.push_back(trk);  }
 
if (AddLostTracks){
   for (const pat::PackedCandidate &trk : *tracks2) {
    tracks.push_back(trk); }}

 if (saveTracks){
  for (const pat::PackedCandidate &trk : tracks){      
   if (trk.pt()<0.2) continue;
   if(trk.charge()==0) continue;
   if(fabs(trk.pdgId())==11 || fabs(trk.pdgId())==13 || fabs(trk.pdgId())==22 || fabs(trk.pdgId())==130) continue;
   bool skip=false;
   for (unsigned int i=0; i<nt.muon_eta.size(); i++){
     if (DR(nt.muon_eta[i],nt.muon_phi[i],trk.eta(),trk.phi())<0.02) skip=true;
   }
   if (!IsLowpTE){
     for (unsigned int i=0; i<nt.el_eta.size(); i++){
       if (DR(nt.el_eta[i],nt.el_phi[i],trk.eta(),trk.phi())<0.02) skip=true;
     }
   }
   if (skip) continue;
   nt.track_pt.push_back(trk.pt()); nt.track_eta.push_back(trk.eta());
   nt.track_phi.push_back(trk.phi()); nt.track_charge.push_back(trk.charge());
   if(trk.trackHighPurity()) nt.track_highPurity.push_back(1);
   else nt.track_highPurity.push_back(0);
   if(trk.hasTrackDetails()){  
     const reco::Track ptrk= trk.pseudoTrack();
     nt.track_norm_chi2.push_back(ptrk.normalizedChi2());
   }
   else nt.track_norm_chi2.push_back(-1);
   nt.track_dxy.push_back(trk.dxy(vertex_point));
   nt.track_dz.push_back(trk.dz(vertex_point));
   nt.track_validhits.push_back(trk.numberOfHits());
   nt.track_losthits.push_back(trk.lostInnerHits());
   nt.track_fromPV.push_back(trk.fromPV()); nt.ntracks++;  
  }
 }//save tracks
 
 nmupairs=0;
 muTrack1.clear(); muTrack2.clear();//,eTrack1,eTrack2;
 used_muTrack_index.clear(); used_eTrack_index.clear();
 if( (reconstructBMuMuK || reconstructBMuMuKstar) && !OnlyKee){
  for (std::vector<pat::Muon>::const_iterator mu=muons->begin(); mu!=muons->end(); ++mu){
    if (!mu->isSoftMuon(firstGoodVertex)) continue;
    if(mu->pt()< muon_pt_cut_forB) continue;
    auto i=mu-muons->begin();
    for (std::vector<pat::Muon>::const_iterator mu2=mu+1; mu2!=muons->end(); ++mu2){
     if (!mu2->isSoftMuon(firstGoodVertex)) continue; 
     if(mu2->pt()< muon_pt_cut_forB) continue;
     if (mu->charge()==mu2->charge()) continue;
     auto i2=mu2-muons->begin();
     TLorentzVector vmu1,vmu2; 
     vmu1.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),0.105);
     vmu2.SetPtEtaPhiM(mu2->pt(),mu2->eta(),mu2->phi(),0.105);
     if((vmu1+vmu2).M()<MLLmin_Cut || (vmu1+vmu2).M()>MLLmax_Cut) continue;
     const Track * mutrack1= mu->bestTrack();
     const Track * mutrack2= mu2->bestTrack();      
     auto mt1=std::make_shared<reco::TransientTrack> (reco::TransientTrack(*mutrack1,&(*bFieldHandle)));
     auto mt2=std::make_shared<reco::TransientTrack> (reco::TransientTrack(*mutrack2,&(*bFieldHandle)));
     muTrack1.push_back(mt1); muTrack2.push_back(mt2);
     used_muTrack_index.emplace_back(i,i2);
     nmupairs++;
    }
  }
 }

int cel=-1;
if((AddeeK || OnlyKee )  && (reconstructBMuMuK || reconstructBMuMuKstar) ){
  for(std::vector<pat::Electron>::const_iterator el=electrons->begin(); el!=electrons->end(); ++el){
   float pt1=0,eta1=-30,phi1=-30;
   if (!el->passConversionVeto()) continue;
   if (fabs(TrgmuDz-el->vz())>ElectronDzCut) continue;
   if (IsLowpTE){
     reco::GsfTrackRef seed = el->gsfTrack();
     if ( seed.isNull() ) continue;
     pt1=seed->ptMode(); eta1=seed->etaMode(); phi1=seed->phiMode();
     if ( (*ele_mva_wp_unbiased)[seed]<MVAEl2Cut) continue;
   } else{
     pt1=el->pt(); eta1=el->eta(); phi1=el->phi();
   }
   if (pt1<Electron2PtCut) continue;
   cel++;
   int cel2=cel;
   if (UseDirectlyGenBeeK && !data){
    if (DR(EtaPhiE1.first,EtaPhiE1.second,eta1,phi1)>DRgenCone && DR(EtaPhiE2.first,EtaPhiE2.second,eta1,phi1)>DRgenCone)
      continue;
    }
    for (std::vector<pat::Electron>::const_iterator el2=el+1; el2!=electrons->end(); ++el2){
      float pt2=0,eta2=-30,phi2=-30;
      if (!el2->passConversionVeto()) continue; 
      if (fabs(TrgmuDz-el2->vz())>ElectronDzCut) continue;
      if (IsLowpTE){
        reco::GsfTrackRef seed2 = el2->gsfTrack();
        if ( seed2.isNull() ) continue;
        pt2=seed2->ptMode(); eta2=seed2->etaMode(); phi2=seed2->phiMode();
        if (pt2<Electron2PtCut) continue;
        if ( (*ele_mva_wp_unbiased)[seed2]<MVAEl2Cut) continue;
      } else{
        pt2=el2->pt(); eta2=el2->eta(); phi2=el2->phi();
        if (pt2<Electron2PtCut) continue;
      } 
      cel2++;
      if (UseDirectlyGenBeeK && !data){
        if (DR(EtaPhiE1.first,EtaPhiE1.second,eta2,phi2)>DRgenCone && DR(EtaPhiE2.first,EtaPhiE2.second,eta2,phi2)>DRgenCone)
	  continue;
      }
      if (el->charge()==el2->charge()) continue;
      if (pt2<Electron1PtCut && pt1<Electron1PtCut) continue;
      if (IsLowpTE){
        reco::GsfTrackRef seed = el->gsfTrack();
        reco::GsfTrackRef seed2 = el2->gsfTrack();
        if ( (*ele_mva_wp_unbiased)[seed2]<MVAEl1Cut &&  (*ele_mva_wp_unbiased)[seed]<MVAEl1Cut) continue;
      }
      TLorentzVector vel1,vel2; 
      vel1.SetPtEtaPhiM(pt1,eta1,phi1,0.000511);
      vel2.SetPtEtaPhiM(pt2,eta2,phi2,0.000511);
      if((vel1+vel2).M()<MLLmin_Cut || (vel1+vel2).M()>MLLmax_Cut) continue; 
      const Track * etrack1= el->bestTrack();
      const Track * etrack2= el2->bestTrack();      
      auto mt1=std::make_shared<reco::TransientTrack> (reco::TransientTrack(*etrack1,&(*bFieldHandle)));
      auto mt2=std::make_shared<reco::TransientTrack> (reco::TransientTrack(*etrack2,&(*bFieldHandle)));
      muTrack1.push_back(mt1); muTrack2.push_back(mt2);
      used_eTrack_index.emplace_back(cel,cel2);
    }
}
}

if(OnlyKee && used_eTrack_index.size()==0 && ( ( reconstructBMuMuK && SkipEventWithNoBToMuMuK) || (SkipEventWithNoBToMuMuKstar && reconstructBMuMuKstar) ) ) return;
KTrack.clear(); KTrack_index.clear();
int index=-1;
if((used_muTrack_index.size()>0 || used_eTrack_index.size()>0) &&( reconstructBMuMuK || reconstructBMuMuKstar)){
for ( const pat::PackedCandidate & trk: tracks){
  index++;
  if(trk.charge()==0) continue;
   if(fabs(trk.pdgId())!=211) continue;
   if(!trk.hasTrackDetails())continue;
   if (trk.pt()< track_pt_cut_forB) continue;
   if (fabs(trk.eta())>EtaTrk_Cut) continue;
   if (UseDirectlyGenBeeK){
      if (DR(EtaPhiK.first,EtaPhiK.second,trk.eta(),trk.phi())>DRgenCone)
	continue;
    }
   if (fabs(TrgmuDz-trk.vz())>ElectronDzCut) continue;
   bool isMu=false; bool isE=false;
   for (const pat::Muon & mu : *muons)
       if (DR(mu.eta(),mu.phi(),trk.eta(),trk.phi())<LepTrkExclusionCone) isMu=true;
   if (!IsLowpTE){
     for (const pat::Electron & el : *electrons) 
       if (DR(el.eta(),el.phi(),trk.eta(),trk.phi())<LepTrkExclusionCone) isE=true;}
   if(isMu || isE ) { continue; std::cout<<"mu or e "<<std::endl;}
   const reco::Track & ptrk=trk.pseudoTrack();
   auto Ktrk=std::make_shared<reco::TransientTrack> (reco::TransientTrack(ptrk,&(*bFieldHandle)));
   KTrack.push_back(Ktrk); KTrack_index.push_back(&trk-&tracks[0]);
   }
 }

KstarTrack.clear(); KstarTrack_index.clear();
if (reconstructBMuMuKstar){
  for(unsigned int iks=0; iks<KTrack_index.size(); iks++){
   unsigned int iks1=KTrack_index[iks];
   const pat::PackedCandidate & trk1=tracks[iks1];  
   for(unsigned int iks2=iks+1; iks2<KTrack_index.size(); iks2++){
     unsigned int iks2b=KTrack_index[iks2];
     const pat::PackedCandidate & trk2= tracks[iks2b];
     if (trk1.charge()==trk2.charge()) continue;
     TLorentzVector vK,vPi;
     vK.SetPtEtaPhiM(trk1.pt(),trk1.eta(),trk1.phi(),0.493);
     vPi.SetPtEtaPhiM(trk2.pt(),trk2.eta(),trk2.phi(),0.139);
     if ( (vK+vPi).M()>MKstarMin_Cut-0.1 && (vK+vPi).M()<MKstarMax_Cut+0.1){   
       KinematicParticleFactoryFromTransientTrack pFactory;
       ParticleMass kaon_mass = 0.493677; float kaon_sigma = 0.000016;
       ParticleMass pion_mass = 0.139; float pion_sigma = 0.000016;
       float chi = 0.; float ndf = 0.;
       vector<RefCountedKinematicParticle> allParticles;
       allParticles.push_back(pFactory.particle(*KTrack.at(iks),kaon_mass,chi,ndf,kaon_sigma));
       allParticles.push_back(pFactory.particle(*KTrack.at(iks2),pion_mass,chi,ndf,pion_sigma));
       TripleTrackKinFit fitter(allParticles);
       if (fitter.success()) {        
	 if (fitter.Mother_Mass(true)>MKstarMin_Cut && fitter.Mother_Mass(true)<MKstarMax_Cut){
           KstarTrack.emplace_back(std::make_pair(KTrack[iks],KTrack[iks2]));
           KstarTrack_index.emplace_back(std::make_pair(iks1,iks2b));
	  }
        }
      }
     vK.SetPtEtaPhiM(trk1.pt(),trk1.eta(),trk1.phi(),0.139);
     vPi.SetPtEtaPhiM(trk2.pt(),trk2.eta(),trk2.phi(),0.493);
     if ((vK+vPi).M()>MKstarMin_Cut-0.1 && (vK+vPi).M()<MKstarMax_Cut+0.1){
       std::vector<reco::TransientTrack> tempTracks;
       tempTracks.push_back(*KTrack.at(iks)); 
       tempTracks.push_back(*KTrack.at(iks2)); 
       KinematicParticleFactoryFromTransientTrack pFactory;
       ParticleMass kaon_mass = 0.493677; float kaon_sigma = 0.000016;
       ParticleMass pion_mass = 0.139; float pion_sigma = 0.000016;
       float chi = 0.; float ndf = 0.;
       vector<RefCountedKinematicParticle> allParticles;
       allParticles.push_back(pFactory.particle(*KTrack.at(iks),pion_mass,chi,ndf,pion_sigma));
       allParticles.push_back(pFactory.particle(*KTrack.at(iks2),kaon_mass,chi,ndf,kaon_sigma));
       TripleTrackKinFit fitter(allParticles);
       if (fitter.success()) {        
	 if (fitter.Mother_Mass(true)>MKstarMin_Cut && fitter.Mother_Mass(true)<MKstarMax_Cut){     
             KstarTrack.emplace_back(std::make_pair(KTrack[iks2],KTrack[iks]));
             KstarTrack_index.emplace_back(std::make_pair(iks2b,iks1));
	     }
	  }
     }  
   }//trk2
}//trk1
}
  if (SkipEventWithNoBToMuMuKstar && KstarTrack.size()==0) return;
//add ee chanel
 
 if (used_muTrack_index.size()>0 && used_eTrack_index.size()>0 ){
    used_muTrack_index.insert(used_muTrack_index.end(),used_eTrack_index.begin(),used_eTrack_index.end());
 }
else if (used_muTrack_index.size()==0 && used_eTrack_index.size()>0 )
    used_muTrack_index=used_eTrack_index;
 
//cout<<"event "<<event<<endl;
 ///building B 
if (used_muTrack_index.size()>0 && KTrack_index.size()>0 && reconstructBMuMuK){
for(unsigned int imu=0; imu<used_muTrack_index.size(); imu++){
  bool IsE=false; //cout<<"new pair "<<imu<<endl;
   
   if ((imu>nmupairs ||imu==nmupairs )  && AddeeK){ IsE=true; }
  for(unsigned int ik=0; ik<KTrack_index.size(); ik++){
    unsigned int imu1=used_muTrack_index.at(imu).first;
    unsigned int imu2=used_muTrack_index.at(imu).second;
    unsigned int ik1=KTrack_index.at(ik);
    TLorentzVector vk,vmu1,vmu2;
    float m=0.105;
    if (IsE)  m=0.000511;
    if (RefitMuTracksOnly && IsE) RefitTracks=false;
    if (RefitMuTracksOnly && !IsE) RefitTracks=true;
    if (!IsE){
       vmu1.SetPtEtaPhiM(nt.muon_pt.at(imu1),nt.muon_eta.at(imu1),nt.muon_phi.at(imu1),m);
       vmu2.SetPtEtaPhiM(nt.muon_pt.at(imu2),nt.muon_eta.at(imu2),nt.muon_phi.at(imu2),m);
    }
    else {
      vmu1.SetPtEtaPhiM(nt.el_pt.at(imu1),nt.el_eta.at(imu1),nt.el_phi.at(imu1),m);
      vmu2.SetPtEtaPhiM(nt.el_pt.at(imu2),nt.el_eta.at(imu2),nt.el_phi.at(imu2),m);
    }
    //cout<<"B e1 pt "<<vmu1.Pt()<<" index "<<imu1<<" e2 pt "<<vmu2.Pt()<<" index "<<imu2<<endl;
    typename std::vector<pat::PackedCandidate>::const_iterator trk=tracks.begin();
     std::advance(trk,ik1);         
     vk.SetPtEtaPhiM(trk->pt(),trk->eta(),trk->phi(),0.495);
     if (DR(vmu1.Eta(),vmu1.Phi(),vk.Eta(),vk.Phi())<0.005) continue;
     if (DR(vmu2.Eta(),vmu2.Phi(),vk.Eta(),vk.Phi())<0.005) continue;
     if ((vmu1+vmu2+vk).M()<MBmin_Cut || (vmu1+vmu2+vk).M()>MBmax_Cut) continue;
      KinematicParticleFactoryFromTransientTrack pFactory;
      part_mass = 0.1056583; part_sigma = 0.0000001;
      if (IsE) part_mass=0.000511;
      ParticleMass kaon_mass = 0.493677; float kaon_sigma = 0.000016;
      float chi = 0.; float ndf = 0.;
      vector<RefCountedKinematicParticle> allParticles;
      allParticles.push_back(pFactory.particle(*muTrack1.at(imu),part_mass,chi,ndf,part_sigma));
      allParticles.push_back(pFactory.particle(*muTrack2.at(imu),part_mass,chi,ndf,part_sigma));
      allParticles.push_back(pFactory.particle(*KTrack.at(ik),kaon_mass,chi,ndf,kaon_sigma));
      TripleTrackKinFit fitter(allParticles);
      if (!fitter.success()) continue;
      float ChiProb= ChiSquaredProbability(fitter.chi(),fitter.dof());
      if (ChiProb<Pchi2BMuMuK) continue;
      GlobalPoint Dispbeamspot(-1*((theBeamSpot->x0()-fitter.Mother_XYZ().x())+(fitter.Mother_XYZ().z()-theBeamSpot->z0()) * theBeamSpot->dxdz()),-1*((theBeamSpot->y0()-fitter.Mother_XYZ().y())+ (fitter.Mother_XYZ().z()-theBeamSpot->z0()) * theBeamSpot->dydz()), 0);     
     

      math::XYZVector vperp(Dispbeamspot.x(),Dispbeamspot.y(),0.);     
      math::XYZVector pperp;
      if (IsE){
           pperp.SetXYZ((vmu1+vmu2+vk).Px(),(vmu1+vmu2+vk).Py(),0);}
      else{
         pperp.SetXYZ(fitter.Mother_Momentum(RefitTracks).x(),fitter.Mother_Momentum(RefitTracks).y(),0);}
      if (vperp.Dot(pperp)/(vperp.R()*pperp.R())<CosThetaCut) continue;   
      // if ( fitter.Mother_Mass(RefitTracks)<MBmin_Cut || fitter.Mother_Mass(RefitTracks)>MBmax_Cut) continue;
      nt.NRb_cosTheta2D.push_back(vperp.Dot(pperp)/(vperp.R()*pperp.R()));
      nt.NRb_bspot_lxy.push_back( Dispbeamspot.perp());
      nt.NRb_bspot_elxy.push_back(fitter.Mother_XYZError().rerr(Dispbeamspot));
      nt.NRb_trk_sdxy.push_back(KTrack.at(ik)->track().dxy(vertex_point)/KTrack.at(ik)->track().dxyError());      
      nt.NRb_mass.push_back(fitter.Mother_Mass(RefitTracks)); 
      std::vector<float> tempK; tempK.push_back(vk.Pt());
      tempK.push_back(vk.Eta()); tempK.push_back(vk.Phi());
      nt.NRb_KUNFITpt_eta_phi.push_back(tempK);
      std::vector<float> tempBpt; 
      tempBpt.push_back(fitter.Mother_Momentum(RefitTracks).perp());
      tempBpt.push_back(fitter.Mother_Momentum(RefitTracks).eta());
      tempBpt.push_back(fitter.Mother_Momentum(RefitTracks).phi());
      nt.NRb_pt_eta_phi.push_back(tempBpt);
      nt.NRb_charge.push_back(fitter.Mother_Charge());
      std::vector<float> tempBx; 
      tempBx.push_back(fitter.Mother_XYZ().x());
      tempBx.push_back(fitter.Mother_XYZ().y());
      tempBx.push_back(fitter.Mother_XYZ().z());
      nt.NRb_x_y_z.push_back(tempBx);
      std::vector<float> tempBex;
      tempBex.push_back(fitter.Mother_XYZError().cxx());
      tempBex.push_back(fitter.Mother_XYZError().cyy());
      tempBex.push_back(fitter.Mother_XYZError().czz());
      nt.NRb_ex_ey_ez.push_back(tempBex); 
      std::vector<float> tempBept;
      tempBept.push_back(fitter.Mother_PtError());
      tempBept.push_back(fitter.Mother_EtaError());
      tempBept.push_back(fitter.Mother_PhiError());
      nt.NRb_ept_eeta_ephi.push_back(tempBept);
      nt.NRb_chi_prob.push_back(ChiProb); 
      if(IsE) nt.NRb_mudecay.push_back(0);
      else  nt.NRb_mudecay.push_back(1);
	 
      nt.NRb_lep1Id.push_back(imu1); nt.NRb_lep2Id.push_back(imu2);
      
      TLorentzVector refl1,refl2;
     for(unsigned int ichild=0; ichild<allParticles.size(); ichild++){
       std::vector<float> temp;
       temp.push_back(fitter.Daughter_Momentum(ichild,RefitTracks).perp());
       temp.push_back(fitter.Daughter_Momentum(ichild,RefitTracks).eta());
       temp.push_back(fitter.Daughter_Momentum(ichild,RefitTracks).phi());
       temp.push_back(fitter.Daughter_Charge(ichild,RefitTracks));
        
       if(ichild==2)
           nt.NRb_Kpt_eta_phi.push_back(temp);
       else if (ichild==0){
          nt.NRb_l1pt_eta_phi.push_back(temp);
          refl1.SetPtEtaPhiM(temp[0],temp[1],temp[2],m); }
       else{
          nt.NRb_l2pt_eta_phi.push_back(temp);    
          refl2.SetPtEtaPhiM(temp[0],temp[1],temp[2],m); }      
     }
     nt.NRb_mll.push_back((refl1+refl2).M()); nt.NRb_vtx_index.push_back(-1);
    //isolation
    double temp_iso04=0,temp_iso08=0;
    for(unsigned int ikIso=0; ikIso<KTrack.size(); ikIso++){
        Track temptrk=KTrack[ikIso]->track();
        if (DR(tempBpt[1],tempBpt[2],temptrk.eta(),temptrk.phi())<0.4){
           temp_iso04+=temptrk.pt(); 
           }
        if (DR(tempBpt[1],tempBpt[2],temptrk.eta(),temptrk.phi())<0.8){
            temp_iso08+=temptrk.pt(); 
           }
       }
     nt.NRb_iso04.push_back( temp_iso04); nt.NRb_iso08.push_back( temp_iso08); 
  
   }
 }
} 

if (used_muTrack_index.size()>0 && KTrack_index.size()>0 && reconstructBMuMuKstar){
for(unsigned int imu=0; imu<used_muTrack_index.size(); imu++){
  unsigned int imu1=used_muTrack_index.at(imu).first;
  unsigned int imu2=used_muTrack_index.at(imu).second;
  bool IsE=false;
  if ((imu>nmupairs ||imu==nmupairs )  && AddeeK) IsE=true; 
  for(unsigned int ik=0; ik<KstarTrack_index.size(); ik++){
    unsigned int ik1=KstarTrack_index.at(ik).first;
    unsigned int ipi1=KstarTrack_index.at(ik).second;
    TLorentzVector vk,vpi,vmu1,vmu2;
    float m=0.105;
    if (IsE)  m=0.000511;
    if (RefitMuTracksOnly && IsE) RefitTracks=false;
    if (RefitMuTracksOnly && !IsE) RefitTracks=true;
    if (!IsE){
      vmu1.SetPtEtaPhiM(nt.muon_pt.at(imu1),nt.muon_eta.at(imu1),nt.muon_phi.at(imu1),m);
      vmu2.SetPtEtaPhiM(nt.muon_pt.at(imu2),nt.muon_eta.at(imu2),nt.muon_phi.at(imu2),m);
      }
   else {
      vmu1.SetPtEtaPhiM(nt.el_pt.at(imu1),nt.el_eta.at(imu1),nt.el_phi.at(imu1),m);
      vmu2.SetPtEtaPhiM(nt.el_pt.at(imu2),nt.el_eta.at(imu2),nt.el_phi.at(imu2),m);
      }
   
    typename std::vector<pat::PackedCandidate>::const_iterator trk=tracks.begin();
    std::advance(trk,ik1);  
    typename std::vector<pat::PackedCandidate>::const_iterator trk2=tracks.begin();
    std::advance(trk2,ipi1);
    vk.SetPtEtaPhiM(trk->pt(),trk->eta(),trk->phi(),0.493);
    vpi.SetPtEtaPhiM(trk2->pt(),trk2->eta(),trk2->phi(),0.139);

     if ((vmu1+vmu2+vk+vpi).M()<MBmin_Cut || (vmu1+vmu2+vk+vpi).M()>MBmax_Cut) continue;
      
      KinematicParticleFactoryFromTransientTrack pFactory;
      part_mass = 0.1056583; part_sigma = 0.0000001;
      if (IsE) part_mass=0.000511;
      ParticleMass kaon_mass = 0.493677; float kaon_sigma = 0.000016;
      ParticleMass pion_mass = 0.139; float pion_sigma = 0.000016;
      float chi = 0.; float ndf = 0.;
      vector<RefCountedKinematicParticle> allParticles;
      allParticles.push_back(pFactory.particle(*muTrack1.at(imu),part_mass,chi,ndf,part_sigma));
      allParticles.push_back(pFactory.particle(*muTrack2.at(imu),part_mass,chi,ndf,part_sigma));
      allParticles.push_back(pFactory.particle(*KstarTrack.at(ik).first,kaon_mass,chi,ndf,kaon_sigma));
      allParticles.push_back(pFactory.particle(*KstarTrack.at(ik).second,pion_mass,chi,ndf,pion_sigma));
      //      ParticleMass Kstar_m = 0.896;
      TripleTrackKinFit fitter(allParticles);
      if (!fitter.success()) continue;

      float ChiProb= ChiSquaredProbability(fitter.chi(),fitter.dof());
      if (ChiProb<Pchi2BMuMuK) continue;

      nt.NRbks_k_sdxy.push_back(KstarTrack.at(ik).first->track().dxy(vertex_point)/KstarTrack.at(ik).first->track().dxyError());      
      nt.NRbks_pi_sdxy.push_back(KstarTrack.at(ik).second->track().dxy(vertex_point)/KstarTrack.at(ik).second->track().dxyError());
      nt.NRbks_mass.push_back(fitter.Mother_Mass(RefitTracks)); 

      std::vector<float> tempBpt; 
      tempBpt.push_back(fitter.Mother_Momentum(RefitTracks).perp());
      tempBpt.push_back(fitter.Mother_Momentum(RefitTracks).eta());
      tempBpt.push_back(fitter.Mother_Momentum(RefitTracks).phi());
      nt.NRbks_pt_eta_phi.push_back(tempBpt);
      nt.NRbks_charge.push_back(fitter.Mother_Charge());
      std::vector<float> tempBx; 
      tempBx.push_back(fitter.Mother_XYZ().x());
      tempBx.push_back(fitter.Mother_XYZ().y());
      tempBx.push_back(fitter.Mother_XYZ().z());
      nt.NRbks_x_y_z.push_back(tempBx);
      std::vector<float> tempBex;
      tempBex.push_back(fitter.Mother_XYZError().cxx());
      tempBex.push_back(fitter.Mother_XYZError().cyy());
      tempBex.push_back(fitter.Mother_XYZError().czz());
      nt.NRbks_ex_ey_ez.push_back(tempBex); 
      std::vector<float> tempBept;
      tempBept.push_back(fitter.Mother_PtError());
      tempBept.push_back(fitter.Mother_EtaError());
      tempBept.push_back(fitter.Mother_PhiError());
      nt.NRbks_ept_eeta_ephi.push_back(tempBept);
      nt.NRbks_chi_prob.push_back(ChiProb); 
      if(IsE) nt.NRbks_mudecay.push_back(0);
      else  nt.NRbks_mudecay.push_back(1);
	 
      nt.NRbks_lep1Id.push_back(imu1); nt.NRbks_lep2Id.push_back(imu2);
      GlobalPoint Dispbeamspot(-1*((theBeamSpot->x0()-fitter.Mother_XYZ().x())+(fitter.Mother_XYZ().z()-theBeamSpot->z0()) * theBeamSpot->dxdz()),-1*((theBeamSpot->y0()-fitter.Mother_XYZ().y())+ (fitter.Mother_XYZ().z()-theBeamSpot->z0()) * theBeamSpot->dydz()), 0);     
      nt.NRbks_bspot_lxy.push_back( Dispbeamspot.perp());
      nt.NRbks_bspot_elxy.push_back(fitter.Mother_XYZError().rerr(Dispbeamspot));
      math::XYZVector pperp;
      if (IsE && UsePFeForCos){
           pperp.SetXYZ((vmu1+vmu2+vk).Px(),(vmu1+vmu2+vk).Py(),0);}
      else{
         pperp.SetXYZ(fitter.Mother_Momentum(RefitTracks).x(),fitter.Mother_Momentum(RefitTracks).y(),0);}
      math::XYZVector vperp(Dispbeamspot.x(),Dispbeamspot.y(),0.);
      nt.NRbks_cosTheta2D.push_back(vperp.Dot(pperp)/(vperp.R()*pperp.R()));
      
      TLorentzVector refl1,refl2,refK,refPi;
      for(unsigned int ichild=0; ichild<allParticles.size(); ichild++){
       std::vector<float> temp;
       temp.push_back(fitter.Daughter_Momentum(ichild,RefitTracks).perp());
       temp.push_back(fitter.Daughter_Momentum(ichild,RefitTracks).eta());
       temp.push_back(fitter.Daughter_Momentum(ichild,RefitTracks).phi());
       temp.push_back(fitter.Daughter_Charge(ichild,RefitTracks));
   
       if(ichild==3){
          nt.NRbks_Pipt_eta_phi.push_back(temp);
          refPi.SetPtEtaPhiM(temp[0],temp[1],temp[2],0.139);}
       else if(ichild==2){
          nt.NRbks_Kpt_eta_phi.push_back(temp);
          refK.SetPtEtaPhiM(temp[0],temp[1],temp[2],0.493); }
       else if (ichild==0){
          nt.NRbks_l1pt_eta_phi.push_back(temp);
          refl1.SetPtEtaPhiM(temp[0],temp[1],temp[2],0.105); }
       else{
          nt.NRbks_l2pt_eta_phi.push_back(temp);    
          refl2.SetPtEtaPhiM(temp[0],temp[1],temp[2],0.105); }     
     }
     nt.NRbks_mll.push_back((refl1+refl2).M());  nt.NRbks_ksmass.push_back((refK+refPi).M()); 
//   cout<<"M(K,pi) "<<(refK+refPi).M()<<"  M(l1,l2) "<<(refl1+refl2).M()<<endl;
   }
 }
} 
 
  if (nt.NRb_mass.size()==0 && SkipEventWithNoBToMuMuK) return ;
  if (nt.NRbks_mass.size()==0 && SkipEventWithNoBToMuMuKstar) return ;

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
  if (NtupleOutputClasses=="auto"){
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
