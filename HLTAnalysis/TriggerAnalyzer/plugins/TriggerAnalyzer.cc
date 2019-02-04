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
//                george karathanasis
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
#include "DataFormats/Common/interface/ValueMap.h"


#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateIsolation.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "HLTrigger/Egamma/plugins/HLTGenericFilter.h"
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

  float Dphi(float phi1,float phi2);
  float DR(float eta1,float phi1,float eta2, float phi2);

  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken_;
  edm::EDGetToken electronsToken_;
  edm::EDGetToken muonsToken_;
  edm::EDGetToken jetsToken_;
  edm::EDGetToken metToken_;
  edm::EDGetToken photonToken_;
  edm::EDGetToken PFCands_;
  edm::EDGetToken LostTracks_;
//  edm::EDGetToken Tracks_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapVetoToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapSoftToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapMediumToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapTightToken_;
  edm::EDGetTokenT<edm::ValueMap<int> > elIdMapValueToken_;
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
 
 int totb1=0,totb2=0,totb3=0;
  int trigger1=0,trigger2=0,trigger3=0,trigger4=0,trigger5=0,trigger6=0,trigger7=0,trigger8=0;
  int l1_seed1=0,l1_seed2=0,l1_seed3=0,l1_seed4=0,l1_seed5=0,l1_seed6=0;
  int event=0,nmuons=0,njets=0,total_triggers=0,good_vertex=1,nel=0,njpsi_mumu=0,njpsi_ee=0;
  float beam_x=0,beam_y=0,beam_z=0,beam_ex=0,beam_ey=0,beam_ez=0;
  double ptmet=0,phimet=-999,l1met=0,l1met_phi=-999,l1ht=0,l1hf_met=0,l1hf_met_phi=-999;
  int trg_counter=0,fire_counter=0,l3_counter=0,ntracks=0,nkaons=0;
  std::vector<float> muon_pt,muon_eta,muon_phi,muon_qual,muon_charge,muon_global_flag,muon_tracker_flag,muon_standalone_flag,muon_dxy,muon_dz,muon_edxy,muon_edz,muon_d0,muon_ed0,muon_vx,muon_vz,muon_vy,muon_iso,l1muon_eta,l1muon_iso,l1muon_pt,l1muon_phi,l1muon_qual,l1jet_pt,l1jet_phi,l1jet_eta,l1jet_iso,l1jet_qual,muon_trkpt,muon_trketa,muon_trkphi;
  std::vector<std::vector<float>> tr1_obj_pt_eta_phi,tr2_obj_pt_eta_phi,tr3_obj_pt_eta_phi,tr4_obj_pt_eta_phi,tr5_obj_pt_eta_phi,tr6_obj_pt_eta_phi,tr7_obj_pt_eta_phi,tr8_obj_pt_eta_phi;
  std::vector<float> el_pt,el_eta,el_phi,el_charge,el_vx,el_vy,el_vz,el_dxy,el_dz,el_mva_out,el_mva_iso,el_iso,el_trkpt,el_trketa,el_trkphi,el_edxy,el_edz;
  float pvertex_x,pvertex_y,pvertex_z,pvertex_ex,pvertex_ey,pvertex_ez;
  std::vector<float> vertex_x,vertex_y,vertex_z,vertex_ex,vertex_ey,vertex_ez,vertex_chi,vertex_ndof;
  std::vector<bool>el_veto,el_soft,el_medium,el_tight;
  std::vector<float>el_mva_map_value;
  std::vector<float> jet_pt,jet_eta,jet_phi,jet_cEmEF,jet_cHEF,jet_cHMult,jet_cMuEF,jet_cMult,jet_MuEF,jet_nEmEF,jet_nHEF,jet_nMult,jet_pEF,jet_eEF;
  std::vector<bool> muon_medium,muon_loose,muon_tight,muon_soft;
  std::vector<float> charged_pf_pt,charged_pf_eta,charged_pf_phi;
  std::vector<float> track_pt,track_eta,track_phi,track_norm_chi2,track_charge,track_dxy,track_dz,track_validhits,track_losthits,track_fromPV,track_highPurity;
 
  std::vector<std::vector<int>>muon_vtx_IsValid,el_vtx_IsValid;
int run_number=0,ls=0;
  int npv=0;

;
  std::vector<float> NRb_mass,NRb_chi_prob,NRb_charge,NRb_mll,NRb_ll_prob,NRb_trk_sdxy,NRb_trk_chi_norm,NRb_vtx_index,NRb_bspot_lxy,NRb_bspot_elxy,NRb_llsvbsv,NRb_cosTheta2D,NRb_iso04,NRb_iso08,NRb_biso02,NRb_biso04,NRb_biso06,NRb_biso1p2;

  std::vector<std::vector<float>> NRb_Kpt_eta_phi,NRb_KUNFITpt_eta_phi,NRb_pt_eta_phi,NRb_l1pt_eta_phi,NRb_l2pt_eta_phi,NRb_x_y_z,NRb_ex_ey_ez,NRb_ept_eeta_ephi;
 std::vector<unsigned int> NRb_mudecay,NRb_lep1Id,NRb_lep2Id;

  float  ngenB;
  std::vector<float> genB_pt,genB_phi,genB_eta,genB_pdgId,genB_Bindex,
  genB_daughter_pt,genB_daughter_eta,genB_daughter_phi,genB_daughter_pdgId,
  genB_daughter_Bindex,genB_daughter_Dindex,genB_granddaughter_pt,
  genB_granddaughter_eta,genB_granddaughter_phi,genB_granddaughter_pdgId,
  genB_granddaughter_Bindex,genB_granddaughter_Dindex;
  float ngenLep;
  std::vector<float> genLep_pt,genLep_phi,genLep_eta,genLep_pdgId,genLep_mom;
  std::vector<float> genMu_pt,genMu_eta,genMu_phi,genMu_ch,genMu_motherId,genMu_gmotherId,genMu_dxy,genMu_edxy,genMu_dir,genMu_indir;
  std::vector<float> NRbks_k_sdxy,NRbks_pi_sdxy,NRbks_mass,NRbks_charge,NRbks_chi_prob,NRbks_bspot_lxy, NRbks_bspot_elxy,NRbks_cosTheta2D; 
  std::vector<std::vector<float>> NRbks_pt_eta_phi,NRbks_x_y_z,NRbks_ex_ey_ez,NRbks_ept_eeta_ephi; 
  std::vector<unsigned int> NRbks_mudecay,NRbks_lep1Id, NRbks_lep2Id;
 
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
std::pair<std::vector<float>,std::vector<std::vector<float>>>  TriggerAnalyzer<T1>::L1Analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

 
  edm::Handle<GlobalAlgBlkBxCollection> l1result;
  iEvent.getByToken(l1resultToken_,l1result);
  edm::Handle<l1t::MuonBxCollection> l1Muons;
  iEvent.getByToken(l1MuonsToken_,l1Muons);
  edm::Handle<l1t::JetBxCollection> l1Jets;
  iEvent.getByToken(l1JetsToken_,l1Jets);
  edm::Handle<BXVector<l1t::EtSum> > l1Met;
  iEvent.getByToken(l1MetToken_,l1Met);
 

  std::vector<float> l1met;
  std::vector< std::vector<float>> l1event;
  std::vector<float> l1muon_pt,l1muon_eta,l1muon_phi,l1muon_iso,l1muon_qual,l1jet_pt,l1jet_eta,l1jet_phi,l1jet_iso,l1jet_qual;
  for(typename std::vector< l1t::Muon >::const_iterator mu=l1Muons->begin(0); mu !=l1Muons->end(0); mu++){
    l1muon_pt.push_back(mu->et());
    l1muon_eta.push_back(mu->eta());
    l1muon_phi.push_back(mu->phi());
    l1muon_iso.push_back(mu->hwIso());
    l1muon_qual.push_back(mu->hwQual());
}
  l1event.push_back(l1muon_pt); l1event.push_back(l1muon_eta); l1event.push_back(l1muon_phi);  l1event.push_back(l1muon_iso);  l1event.push_back(l1muon_qual); 
   for(typename std::vector< l1t::Jet >::const_iterator jet=l1Jets->begin(0); jet!=l1Jets->end(0); jet++){
      l1jet_pt.push_back(jet->et());
      l1jet_eta.push_back(jet->eta());
      l1jet_phi.push_back(jet->phi());
      l1jet_iso.push_back(jet->hwIso());
      l1jet_qual.push_back(jet->hwQual()); 
}
  l1event.push_back(l1jet_pt); l1event.push_back(l1jet_eta); l1event.push_back(l1jet_phi);  l1event.push_back(l1jet_iso);  l1event.push_back(l1jet_qual); 
  for(typename std::vector< l1t::EtSum >::const_iterator et=l1Met->begin(0); et !=l1Met->end(0); et++){
  if (et->getType()==l1t::EtSum::kMissingEt) {l1met.push_back(et->et()); l1met.push_back(et->phi()); }
  if (et->getType()==l1t::EtSum::kTotalHt) l1met.push_back(et->et()); 
  if (et->getType()==l1t::EtSum::kMissingEtHF){ l1met.push_back(et->et()); l1met.push_back(et->phi()); }
    }
  
      return std::make_pair(l1met,l1event);
}
template<typename T1>
std::vector<std::pair<string,int>>  TriggerAnalyzer<T1>::L1SeedAnalyze(const edm::Event& iEvent, TString * algoBitToName,std::vector< string> Seed){
   using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;
   edm::Handle<GlobalAlgBlkBxCollection> l1result;
  iEvent.getByToken(l1resultToken_,l1result); 
  std::vector<std::pair<string,int>> res;
if (l1result.isValid()) {
GlobalAlgBlk const &result = l1result->at(0, 0);
for (unsigned int iseed=0; iseed<Seed.size(); ++iseed){
 bool sfire=false;
 for (unsigned int itrg=0; itrg<result.maxPhysicsTriggers; ++itrg){
   if (result.getAlgoDecisionFinal(itrg)!=1) continue;
   std::string l1trigName = static_cast<const char *>(algoBitToName[itrg]);
   if (l1trigName!=Seed[iseed]) continue;
   sfire=true;
   }
 if (sfire) {
   std::pair<string,int> temp(Seed[iseed],1);
   res.push_back(temp);
  }  
 else {
    std::pair<string,int> temp(Seed[iseed],0);
    res.push_back(temp);
  }
 //  std::cout<<"bit "<<itrg<<" algo "<<algoBitToName[itrg]<<" reult "<<result.getAlgoDecisionFinal(itrg)<<endl;   
}
}
return res;
}



template<typename T1>
std::pair<std::vector<float>,std::vector<std::vector<std::vector<float>>>>  TriggerAnalyzer<T1>::HLTAnalyze(const edm::Event& iEvent, const edm::EventSetup& iSetup,std::vector< string> HLTPath ){
   using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;
 
  edm::Handle<vector<pat::TriggerObjectStandAlone>> triggerObjects;
  iEvent.getByToken(trigobjectsToken_ ,triggerObjects);
  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trgresultsToken_, trigResults);
  std::vector<float> fires;
  std::vector<std::vector<std::vector<float>>> trg_event; 
  const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);
  
for (unsigned int ifilter=0; ifilter<HLTPath.size(); ifilter++){    
    std::vector<std::vector<float>> tot_tr_obj_pt_eta_phi;
     for(pat::TriggerObjectStandAlone itrg :*triggerObjects){
      if(!itrg.id(83)) continue;
       bool save=false;
       itrg.unpackPathNames(trigName);
       std::vector<std::string> const & pathnames = itrg.pathNames(true,true);
       for(std::string name : pathnames){
          if (name.find(HLTPath[ifilter]) != std::string::npos) save=true;
	   }
        std::vector<float> tr_obj_pt_eta_phi;
	if(!save) continue;
        tr_obj_pt_eta_phi.push_back(itrg.pt());
        tr_obj_pt_eta_phi.push_back(itrg.eta());
        tr_obj_pt_eta_phi.push_back(itrg.phi());
        tr_obj_pt_eta_phi.push_back(itrg.charge());
        tot_tr_obj_pt_eta_phi.push_back( tr_obj_pt_eta_phi);
	}
   if(tot_tr_obj_pt_eta_phi.size()>0){
     trg_event.push_back(tot_tr_obj_pt_eta_phi);}            
   else{
     std::vector<float> def; def.push_back(-1); 
     std::vector<std::vector<float>> def2; def2.push_back(def);
     trg_event.push_back(def2);
   }
  }
  //paths
  float fire0=0,fire1=0,fire2=0,fire3=0,fire4=0,fire5=0,fire6=0,fire7=0;
    
    if( !trigResults.failedToGet() ) {
    int N_Triggers = trigResults->size();
    for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
      TString TrigPath =trigName.triggerName(i_Trig);
      if (!trigResults->accept(i_Trig)) continue;
      if (TrigPath.Contains(HLTPath[0])) fire0=1;
      if (TrigPath.Contains(HLTPath[1])) fire1=1;       
      if (TrigPath.Contains(HLTPath[2])) fire2=1; 
      if (TrigPath.Contains(HLTPath[3])) fire3=1;
      if (TrigPath.Contains(HLTPath[4])) fire4=1;
      if (TrigPath.Contains(HLTPath[5])) fire5=1;
      if (TrigPath.Contains(HLTPath[6])) fire6=1;
      if (TrigPath.Contains(HLTPath[7])) fire7=1;
       } 
      }
      							   	
    fires.push_back(fire0); fires.push_back(fire1); fires.push_back(fire2);
    fires.push_back(fire3); fires.push_back(fire4); fires.push_back(fire5);
    fires.push_back(fire6); fires.push_back(fire7);
    //  return std::make_pair(fires,tot_tr_obj_pt_eta_phi);
     return std::make_pair(fires,trg_event);
}


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
    std::auto_ptr<TrajectoryStateClosestToPoint> traj3;
    GlobalPoint vtxPos(myVertex.position().x(), myVertex.position().y(), myVertex.position().z());
     GlobalVector gvmu1,gvmu2,gvmu3;
     if(myVertex.hasRefittedTracks()){
        std::vector<reco::TransientTrack> refited;
        refited=myVertex.refittedTracks();
        reco::TransientTrack* Track1 = &refited[0];
        reco::TransientTrack* Track2= &refited[1];
        traj1.reset(new TrajectoryStateClosestToPoint(Track1->trajectoryStateClosestToPoint(vtxPos)));
        traj2.reset(new TrajectoryStateClosestToPoint(Track2->trajectoryStateClosestToPoint(vtxPos)));    
        if (tracks.size()==3) {
           reco::TransientTrack* Track3= &refited[2];
           traj3.reset(new TrajectoryStateClosestToPoint(Track3->trajectoryStateClosestToPoint(vtxPos)));  
         }
         if(Track1->charge()!=1) cout<<"SKATA!!!!!!!!!!!!!!!!!!!"<<endl;
         if(Track2->charge()!=-1) cout<<"SKATA!!!!!!!!!!!!!!!!!!!"<<endl;
         gvmu1=traj1->momentum();  gvmu2=traj2->momentum(); 
         if (tracks.size()==3) gvmu3=traj3->momentum();
         
                  }         
               else {
                   traj1.reset(new TrajectoryStateClosestToPoint(tracks[0].trajectoryStateClosestToPoint(vtxPos)));
                   traj2.reset(new TrajectoryStateClosestToPoint(tracks[1].trajectoryStateClosestToPoint(vtxPos)));
                 //  traj3.reset(new TrajectoryStateClosestToPoint(tracks[2].trajectoryStateClosestToPoint(vtxPos)));

                    gvmu1=traj1->momentum();  gvmu2=traj2->momentum();
                    if (tracks.size()==3){
                       traj3.reset(new TrajectoryStateClosestToPoint(tracks[2].trajectoryStateClosestToPoint(vtxPos)));
                        gvmu3=traj3->momentum();}
               }               

   std::vector<GlobalVector> gvmu; gvmu.push_back(gvmu1); gvmu.push_back(gvmu2);
   if (tracks.size()==3) gvmu.push_back(gvmu3);
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
//  edm::Handle<vector<reco::Track>> tracks3;
//  iEvent.getByToken(Tracks_, tracks3);
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
  edm::ESHandle<MagneticField> bFieldHandle;
 iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
// GeneratorBTree gen(prunedGenToken_,packedGenToken_,iEvent);
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



 vertex_x.clear();  vertex_y.clear();  vertex_z.clear();  vertex_ex.clear();  vertex_ey.clear();  vertex_ez.clear();  vertex_chi.clear();  vertex_ndof.clear(); npv=0;
  pvertex_x=-9999; pvertex_y=-9999; pvertex_z=-9999; pvertex_ex=-9999; pvertex_ey=-9999; pvertex_ez=-9999,beam_x=-9999,beam_y=-9999,beam_z=-99999,beam_ex=-9999,beam_ey=-9999,beam_ez=-99999;
  footprint.clear(); 
  //continue or not
/*  bool kill_event_all=false,kill_event_1_2_3_5=false;
  if ( (saveOnlyHLTFires || saveOnlyHLT_1_2_3_5_Fires)  && saveHLT){
    std::pair<std::vector<float>,std::vector<std::vector<std::vector<float>>>> trgresult_check=HLTAnalyze(iEvent,iSetup,HLTPath_,HLTFilter_);     
    if (trgresult_check.first[0]==0 && trgresult_check.first[1]==0 && trgresult_check.first[2]==0 &&  trgresult_check.first[3]==0 && trgresult_check.first[4]==0 && trgresult_check.first[5]==0) kill_event_all=true;
     if (trgresult_check.first[0]==0 && trgresult_check.first[1]==0 && trgresult_check.first[2]==0  && trgresult_check.first[4]==0) kill_event_1_2_3_5=true;
   
   }
  if(saveOnlyHLTFires && kill_event_all) { return;}
  if( saveOnlyHLT_1_2_3_5_Fires &&  kill_event_1_2_3_5) { return;}*/
  //Loop 
  // std::cout<<" ev="<<event<<std::endl;
  trigger1=0; trigger2=0; trigger3=0; trigger4=0; trigger5=0; trigger6=0;
   trigger7=0; trigger8=0;
  event=iEvent.id().event(); nmuons=0; njets=0; nel=0; ntracks=0;
  //muons
  muon_pt.clear(); muon_eta.clear(); muon_phi.clear(); muon_qual.clear();
  muon_charge.clear(); muon_dz.clear(); muon_dxy.clear(); muon_global_flag.clear();
  muon_standalone_flag.clear(); muon_tracker_flag.clear();  muon_edz.clear();
  muon_edxy.clear(); muon_d0.clear(); muon_ed0.clear(); muon_vx.clear();
  muon_vy.clear(); muon_vz.clear(); muon_iso.clear(); muon_soft.clear();
  muon_loose.clear(); muon_medium.clear(); muon_tight.clear();
  muon_trkpt.clear(); muon_trketa.clear(); muon_trkphi.clear();
  //e
  el_pt.clear(); el_eta.clear(); el_phi.clear(); el_charge.clear();
  el_vx.clear(); el_vy.clear(); el_vz.clear(); el_dxy.clear(); el_dz.clear();
  el_edxy.clear(); el_edz.clear(); el_mva_out.clear(); el_mva_iso.clear();
  el_iso.clear(); el_veto.clear(); el_soft.clear(); el_medium.clear();
  el_tight.clear();
  //hlt
  tr1_obj_pt_eta_phi.clear(); tr2_obj_pt_eta_phi.clear(); tr3_obj_pt_eta_phi.clear();
  tr4_obj_pt_eta_phi.clear(); tr5_obj_pt_eta_phi.clear(); tr6_obj_pt_eta_phi.clear();
   tr7_obj_pt_eta_phi.clear(); tr8_obj_pt_eta_phi.clear();
  //l1
  l1muon_pt.clear(); l1muon_eta.clear(); l1muon_phi.clear(); l1muon_iso.clear();
  l1muon_qual.clear(); l1jet_pt.clear(); l1jet_eta.clear(); l1jet_phi.clear();
  l1jet_iso.clear(); l1jet_qual.clear();


  ptmet=-1; phimet=-99; jet_pt.clear(); jet_phi.clear(); jet_eta.clear();
  l1met=-1; l1met_phi=-999; l1ht=-1; l1hf_met=-1; l1hf_met_phi=-999;
  

jet_cEmEF.clear(); jet_cHEF.clear(); jet_cHMult.clear(); jet_cMuEF.clear();
jet_cMult.clear(); jet_MuEF.clear(); jet_nEmEF.clear(); jet_nHEF.clear();
jet_nMult.clear(); jet_pEF.clear(); jet_eEF.clear();
 track_pt.clear(); track_eta.clear(); track_phi.clear(); track_norm_chi2.clear(); track_charge.clear();  track_dxy.clear(); track_dz.clear(); track_validhits.clear(); track_losthits.clear(); track_fromPV.clear(); track_highPurity.clear();
 
  el_trkpt.clear(); el_trketa.clear(); el_trkphi.clear();

 

  ngenB=0;
  genB_pt.clear(); genB_phi.clear(); genB_eta.clear(); genB_pdgId.clear(); 
  genB_Bindex.clear(); genB_daughter_pt.clear(); genB_daughter_eta.clear(); 
  genB_daughter_phi.clear(); genB_daughter_pdgId.clear(); genB_daughter_Bindex.clear();
  genB_daughter_Dindex.clear(); genB_granddaughter_pt.clear(); 
  genB_granddaughter_eta.clear(); genB_granddaughter_phi.clear(); 
  genB_granddaughter_pdgId.clear(); genB_granddaughter_Bindex.clear(); 
  genB_granddaughter_Dindex.clear();
 ngenLep=0; 
  genLep_pt.clear(); genLep_phi.clear(); genLep_eta.clear(); 
  genLep_pdgId.clear(); genLep_mom.clear();
 el_mva_map_value.clear();
 NRb_mass.clear(); NRb_chi_prob.clear(); NRb_x_y_z.clear(); NRb_lep1Id.clear(); NRb_lep2Id.clear(); NRb_charge.clear();
 NRb_Kpt_eta_phi.clear(); NRb_pt_eta_phi.clear(); NRb_l1pt_eta_phi.clear();
 NRb_l2pt_eta_phi.clear(); NRb_mll.clear(); NRb_KUNFITpt_eta_phi.clear();
l1_seed1=0; l1_seed2=0; l1_seed3=0; l1_seed4=0; l1_seed5=0; l1_seed6=0;
NRb_ex_ey_ez.clear(); NRb_ept_eeta_ephi.clear(); NRb_mudecay.clear();
 NRb_ll_prob.clear(); NRb_trk_sdxy.clear(); NRb_vtx_index.clear();
 NRb_trk_chi_norm.clear(); NRb_vtx_index.clear(); NRb_bspot_lxy.clear(); 
 NRb_bspot_elxy.clear(); NRb_llsvbsv.clear(); NRb_cosTheta2D.clear();
 NRb_iso04.clear(); NRb_iso08.clear(); NRb_biso02.clear(); NRb_biso04.clear(); 
 NRb_biso06.clear(); NRb_biso1p2.clear();

 genMu_dir.clear(); genMu_indir.clear();

 NRbks_k_sdxy.clear();  NRbks_pi_sdxy.clear(); NRbks_mass.clear(); 
 NRbks_pt_eta_phi.clear(); NRbks_charge.clear(); NRbks_x_y_z.clear();
 NRbks_ex_ey_ez.clear(); NRbks_ept_eeta_ephi.clear(); NRbks_chi_prob.clear(); 
 NRbks_mudecay.clear(); NRbks_lep1Id.clear(); NRbks_lep2Id.clear();
 NRbks_bspot_lxy.clear(); NRbks_bspot_elxy.clear(); NRbks_cosTheta2D.clear();
 

  run_number=iEvent.id().run();
  ls=iEvent.luminosityBlock();

  beam_x= theBeamSpot->x0();  beam_y= theBeamSpot->y0(); beam_z= theBeamSpot->z0();
  beam_ex= theBeamSpot->x0Error(); beam_ey= theBeamSpot->y0Error(); beam_ez= theBeamSpot->z0Error();



   const  reco::Vertex firstGoodVertex=vertices->front();
  int firstGoodVertexIdx = 0;
  float pvx=-99,pvy=-99,pvz=-99,pevx=-99,pevy=-99,pevz=-99;
   for (const reco::Vertex &vtx : *vertices) {
    bool isFake = vtx.isFake();
     if ( isFake || !vtx.isValid () ) continue;
     if (firstGoodVertexIdx==0){
      firstGoodVertexIdx=1; 
      pvx=vtx.x(); pvy=vtx.y(); pvz=vtx.z(); pevx=vtx.xError(); pevy=vtx.yError(); pevz=vtx.zError();    }
      vertex_x.push_back(vtx.x());  vertex_y.push_back(vtx.y());  vertex_z.push_back(vtx.z()); vertex_ex.push_back(vtx.xError());  vertex_ey.push_back(vtx.yError());  vertex_ez.push_back(vtx.zError()); vertex_chi.push_back(vtx.chi2()); vertex_ndof.push_back(vtx.ndof());

  }

   pvertex_x=pvx; pvertex_y=pvy; pvertex_z=pvz; 
   pvertex_ex=pevx; pvertex_ey=pevy; pvertex_ez=pevz;
   reco::TrackBase::Point  vertex_point; vertex_point.SetCoordinates(pvx,pvy,pvz);

if(!data){
     GeneratorBTree gen(prunedGenToken_,packedGenToken_,iEvent);
     ngenB=gen.BMother_Pt().size();
     genB_pt=gen.BMother_Pt(); genB_phi=gen.BMother_Phi();
     genB_eta=gen.BMother_Eta(); genB_pdgId=gen.BMother_PdgId(); 
     genB_Bindex=gen.BMother_Bid();
     genB_daughter_pt=gen.BDaughter_Pt(); genB_daughter_eta=gen.BDaughter_Eta();
     genB_daughter_phi=gen.BDaughter_Phi(); genB_daughter_pdgId=gen.BDaughter_PdgId(); 
     genB_daughter_Bindex=gen.BDaughter_Bid(); genB_daughter_Dindex=gen.BDaughter_Did(); 
     genB_granddaughter_pt=gen.BGDaughter_Pt(); genB_granddaughter_eta=gen.BGDaughter_Eta();
     genB_granddaughter_phi=gen.BGDaughter_Phi(); genB_granddaughter_pdgId=gen.BGDaughter_PdgId();
     genB_granddaughter_Bindex=gen.BGDaughter_Bid();
     genB_granddaughter_Dindex=gen.BGDaughter_Did();
     ngenLep=gen.genLep_Pt().size();     
     genLep_pt=gen.genLep_Pt(); genLep_phi=gen.genLep_Phi();
     genLep_eta=gen.genLep_Eta(); genLep_pdgId=gen.genLep_pdgId(); 
     genLep_mom=gen.genLep_Mother();
 }

 
 //std::vector<std::vector<float>> muskatatest;
 ///////////////////
 //  int snmu=0;
   //    OAEParametrizedMagneticField * paramField_ = new OAEParametrizedMagneticField("3_8T");
    std::vector<reco::TransientTrack> muttks,ettks,jpsimuttks,jpsiettks;
     std::vector<int> muindex,elindex,mucharge,elcharge;
  for (const pat::Muon &mu : *muons){    
     for (unsigned int i=0, n=mu.numberOfSourceCandidatePtrs(); i<n; ++i){
        footprint.push_back(mu.sourceCandidatePtr(i));
    }
    muon_pt.push_back(mu.pt());  muon_phi.push_back(mu.phi());
    muon_eta.push_back(mu.eta()); muon_charge.push_back(mu.charge());
    const Track * mutrack= mu.bestTrack();
    muon_dz.push_back(mutrack->dz(vertex_point));
    muon_dxy.push_back(mutrack->dxy(vertex_point));
    muon_global_flag.push_back(mu.isGlobalMuon());
    muon_standalone_flag.push_back(mu.isStandAloneMuon());
    muon_tracker_flag.push_back(mu.isTrackerMuon());
    muon_vx.push_back(mu.vx());  muon_vy.push_back(mu.vy());
    muon_vz.push_back(mu.vz());  muon_edz.push_back(mu.dzError());
    muon_edxy.push_back(mu.dxyError());  muon_d0.push_back(mutrack->d0());
    muon_ed0.push_back(mutrack->d0Error());
    muon_medium.push_back(mu.isMediumMuon());
    muon_loose.push_back(mu.isLooseMuon());
    muon_trkpt.push_back(mutrack->pt()); muon_trketa.push_back(mutrack->eta());
    muon_trkphi.push_back(mutrack->phi());
    muon_tight.push_back(mu.isTightMuon(firstGoodVertex));
    muon_soft.push_back(mu.isSoftMuon(firstGoodVertex));
   // if (mu.isSoftMuon(firstGoodVertex))snmu++;
    const MuonPFIsolation&  isol=mu.pfIsolationR04();
    double mu_iso=(isol.sumChargedHadronPt+max(0.,isol.sumNeutralHadronEt+isol.sumPhotonEt-0.5*isol.sumPUPt))/mu.pt();
    muon_iso.push_back(mu_iso);
    muttks.push_back(reco::TransientTrack(*mutrack,&(*bFieldHandle)));   
    nmuons++;
    }
    size_type eindex=-1;
    for(const pat::Electron &el : *electrons){
      eindex++;
      for (unsigned int i=0, n=el.numberOfSourceCandidatePtrs(); i<n; ++i){
        footprint.push_back(el.sourceCandidatePtr(i));
      }
      if (!el.passConversionVeto()) continue;
      const edm::Ptr<pat::Electron> elePtr(electrons,eindex);
      el_pt.push_back(el.pt()); el_eta.push_back(el.eta());
      el_phi.push_back(el.phi()); el_charge.push_back(el.charge());
      const Track * eltrack= el.bestTrack();
      el_dz.push_back(eltrack->dz(vertex_point));
      el_dxy.push_back(eltrack->dxy(vertex_point));
      el_edxy.push_back(el.dxyError());  el_edz.push_back(el.dzError());
      el_vx.push_back(el.vx()); el_vy.push_back(el.vy());
      el_vz.push_back(el.vz()); el_mva_out.push_back(el.mva_e_pi());
      el_mva_iso.push_back(el.mva_Isolated());
      double iso= el.pfIsolationVariables().sumChargedHadronPt+max(0.0,el.pfIsolationVariables().sumNeutralHadronEt+el.pfIsolationVariables().sumPhotonEt-0.5*el.pfIsolationVariables().sumPUPt)/el.pt();
      el_iso.push_back(iso); el_veto.push_back((*ele_veto_id)[elePtr]);
      el_soft.push_back((*ele_soft_id)[elePtr]);
      el_medium.push_back((*ele_medium_id)[elePtr]);
      el_tight.push_back((*ele_tight_id)[elePtr]);
      el_mva_map_value.push_back((*ele_mva_id_value)[elePtr]);
      el_trkpt.push_back(eltrack->pt());  el_trketa.push_back(eltrack->eta());
      el_trkphi.push_back(eltrack->phi()); 
      ettks.push_back(reco::TransientTrack(*eltrack,&(*bFieldHandle)));
        
      nel++;
    }
    /*  std::vector<std::vector<float>> mu_DCA=track_DCA(muttks);
  if (mu_DCA.size()==0) { std::vector<float> d1; d1.push_back(-99); muon_DCA.push_back(d1); }
  else muon_DCA=mu_DCA;*/
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
   for (unsigned int i=0; i<muon_eta.size(); i++){
     if (DR(muon_eta[i],muon_phi[i],trk.eta(),trk.phi())<0.02) skip=true;
   }
   for (unsigned int i=0; i<el_eta.size(); i++){
     if (DR(el_trketa[i],el_trkphi[i],trk.eta(),trk.phi())<0.02) skip=true;
   }
   if (skip) continue;
   track_pt.push_back(trk.pt()); track_eta.push_back(trk.eta());
   track_phi.push_back(trk.phi()); track_charge.push_back(trk.charge());
   if(trk.trackHighPurity()) track_highPurity.push_back(1);
   else track_highPurity.push_back(0);
   if(trk.hasTrackDetails()){  
     const reco::Track ptrk= trk.pseudoTrack();
     track_norm_chi2.push_back(ptrk.normalizedChi2());
   }
   else track_norm_chi2.push_back(-1);
   track_dxy.push_back(trk.dxy(vertex_point));
   track_dz.push_back(trk.dz(vertex_point));
   track_validhits.push_back(trk.numberOfHits());
   track_losthits.push_back(trk.lostInnerHits());
   track_fromPV.push_back(trk.fromPV()); ntracks++;  
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
    if (!el->passConversionVeto()) continue;
    cel++;
    if(el->pt()< muon_pt_cut_forB) continue;
    int cel2=cel;
    for (std::vector<pat::Electron>::const_iterator el2=el+1; el2!=electrons->end(); ++el2){
      if (!el2->passConversionVeto()) continue; 
      cel2++;
      if(el2->pt()< muon_pt_cut_forB) continue;
      if (el->charge()==el2->charge()) continue;
      TLorentzVector vel1,vel2; 
      vel1.SetPtEtaPhiM(el->pt(),el->eta(),el->phi(),0.000511);
      vel2.SetPtEtaPhiM(el2->pt(),el2->eta(),el2->phi(),0.000511);
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
if(OnlyKee && used_eTrack_index.size()==0 && ( ( reconstructBMuMuK && SkipEventWithNoBToMuMuK) ) ) return;
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
   if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(&tracks,&trk-&tracks[0])) != footprint.end())    continue;  
   bool isMu=false; bool isE=false;
   for (const pat::Muon & mu : *muons)
       if (DR(mu.eta(),mu.phi(),trk.eta(),trk.phi())<LepTrkExclusionCone) isMu=true;
  for (const pat::Electron & el : *electrons) 
       if (DR(el.eta(),el.phi(),trk.eta(),trk.phi())<LepTrkExclusionCone) isE=true;
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
      if ((vK+vPi).M()>MKstarMin_Cut && (vK+vPi).M()<MKstarMax_Cut){
           KstarTrack.emplace_back(std::make_pair(KTrack[iks],KTrack[iks2]));
           KstarTrack_index.emplace_back(std::make_pair(iks1,iks2b));
       }
      vK.SetPtEtaPhiM(trk1.pt(),trk1.eta(),trk1.phi(),0.139);
      vPi.SetPtEtaPhiM(trk2.pt(),trk2.eta(),trk2.phi(),0.493);
      if ((vK+vPi).M()>MKstarMin_Cut && (vK+vPi).M()<MKstarMax_Cut){
           KstarTrack.emplace_back(std::make_pair(KTrack[iks2],KTrack[iks]));
           KstarTrack_index.emplace_back(std::make_pair(iks2,iks));
       }
   }  
  }
}
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
    // else  RefitTracks=true;
    if (!IsE){
       vmu1.SetPtEtaPhiM(muon_pt.at(imu1),muon_eta.at(imu1),muon_phi.at(imu1),m);
       vmu2.SetPtEtaPhiM(muon_pt.at(imu2),muon_eta.at(imu2),muon_phi.at(imu2),m);
     }
    else {
      vmu1.SetPtEtaPhiM(el_pt.at(imu1),el_eta.at(imu1),el_phi.at(imu1),m);
      vmu2.SetPtEtaPhiM(el_pt.at(imu2),el_eta.at(imu2),el_phi.at(imu2),m);
          }
    //   if(IsE) cout<<" pT "<<vmu1.Pt()<<" pT2 "<<muTrack1.at(imu)->track().pt()<<endl;
    typename std::vector<pat::PackedCandidate>::const_iterator trk=tracks.begin();
     std::advance(trk,ik1);  
        
     vk.SetPtEtaPhiM(trk->pt(),trk->eta(),trk->phi(),0.495);
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
      NRb_trk_sdxy.push_back(KTrack.at(ik)->track().dxy(vertex_point)/KTrack.at(ik)->track().dxyError());      
      NRb_mass.push_back(fitter.Mother_Mass(RefitTracks)); 
      std::vector<float> tempK; tempK.push_back(vk.Pt());
      tempK.push_back(vk.Eta()); tempK.push_back(vk.Phi());
      NRb_KUNFITpt_eta_phi.push_back(tempK);
      std::vector<float> tempBpt; 
      tempBpt.push_back(fitter.Mother_Momentum(RefitTracks).perp());
      tempBpt.push_back(fitter.Mother_Momentum(RefitTracks).eta());
      tempBpt.push_back(fitter.Mother_Momentum(RefitTracks).phi());
      NRb_pt_eta_phi.push_back(tempBpt);
      NRb_charge.push_back(fitter.Mother_Charge());
      std::vector<float> tempBx; 
      tempBx.push_back(fitter.Mother_XYZ().x());
      tempBx.push_back(fitter.Mother_XYZ().y());
      tempBx.push_back(fitter.Mother_XYZ().z());
      NRb_x_y_z.push_back(tempBx);
      std::vector<float> tempBex;
      tempBex.push_back(fitter.Mother_XYZError().cxx());
      tempBex.push_back(fitter.Mother_XYZError().cyy());
      tempBex.push_back(fitter.Mother_XYZError().czz());
      NRb_ex_ey_ez.push_back(tempBex); 
      std::vector<float> tempBept;
      tempBept.push_back(fitter.Mother_PtError());
      tempBept.push_back(fitter.Mother_EtaError());
      tempBept.push_back(fitter.Mother_PhiError());
      NRb_ept_eeta_ephi.push_back(tempBept);
      NRb_chi_prob.push_back(ChiProb); 
      if(IsE) NRb_mudecay.push_back(0);
      else  NRb_mudecay.push_back(1);
	 
      NRb_lep1Id.push_back(imu1); NRb_lep2Id.push_back(imu2);
      GlobalPoint Dispbeamspot(-1*((theBeamSpot->x0()-fitter.Mother_XYZ().x())+(fitter.Mother_XYZ().z()-theBeamSpot->z0()) * theBeamSpot->dxdz()),-1*((theBeamSpot->y0()-fitter.Mother_XYZ().y())+ (fitter.Mother_XYZ().z()-theBeamSpot->z0()) * theBeamSpot->dydz()), 0);     
      NRb_bspot_lxy.push_back( Dispbeamspot.perp());
      NRb_bspot_elxy.push_back(fitter.Mother_XYZError().rerr(Dispbeamspot));

      math::XYZVector vperp(Dispbeamspot.x(),Dispbeamspot.y(),0.);     
      math::XYZVector pperp;
      if (IsE && UsePFeForCos){
           pperp.SetXYZ((vmu1+vmu2+vk).Px(),(vmu1+vmu2+vk).Py(),0);}
      else{
         pperp.SetXYZ(fitter.Mother_Momentum(RefitTracks).x(),fitter.Mother_Momentum(RefitTracks).y(),0);}
      NRb_cosTheta2D.push_back(vperp.Dot(pperp)/(vperp.R()*pperp.R()));
      TLorentzVector refl1,refl2;
     for(unsigned int ichild=0; ichild<allParticles.size(); ichild++){
       std::vector<float> temp;
       temp.push_back(fitter.Daughter_Momentum(ichild,RefitTracks).perp());
       temp.push_back(fitter.Daughter_Momentum(ichild,RefitTracks).eta());
       temp.push_back(fitter.Daughter_Momentum(ichild,RefitTracks).phi());
       temp.push_back(fitter.Daughter_Charge(ichild,RefitTracks));
        
       if(ichild==2)
           NRb_Kpt_eta_phi.push_back(temp);
       else if (ichild==0){
          NRb_l1pt_eta_phi.push_back(temp);
          refl1.SetPtEtaPhiM(temp[0],temp[1],temp[2],m);     }
       else{
          NRb_l2pt_eta_phi.push_back(temp);    
          refl2.SetPtEtaPhiM(temp[0],temp[1],temp[2],m);     }
       
     }
     NRb_mll.push_back((refl1+refl2).M()); NRb_vtx_index.push_back(-1);
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
     NRb_iso04.push_back( temp_iso04); NRb_iso08.push_back( temp_iso08); 
  
   }
 }
} 

if (used_muTrack_index.size()>0 && KTrack_index.size()>0 && reconstructBMuMuKstar){
for(unsigned int imu=0; imu<used_muTrack_index.size(); imu++){
  for(unsigned int ik=0; ik<KstarTrack_index.size(); ik++){
    unsigned int imu1=used_muTrack_index.at(imu).first;
    unsigned int imu2=used_muTrack_index.at(imu).second;
    unsigned int ik1=KstarTrack_index.at(ik).first;
    unsigned int ipi1=KstarTrack_index.at(ik).second;
    TLorentzVector vk,vpi,vmu1,vmu2;
    bool IsE=false;
    if ((imu>nmupairs ||imu==nmupairs )  && AddeeK) IsE=true;
    float m=0.105;
    if (IsE)  m=0.000511;
    if (RefitMuTracksOnly && IsE) RefitTracks=false;
    else  RefitTracks=true;
    if (!IsE){
       auto  mu1=muons->at(imu1);  auto  mu2=muons->at(imu2);
       vmu1.SetPtEtaPhiM(mu1.pt(),mu1.eta(),mu1.phi(),m);
       vmu2.SetPtEtaPhiM(mu2.pt(),mu2.eta(),mu2.phi(),m);
     }
   else {
      auto  el1=electrons->at(imu1);  auto  el2=electrons->at(imu2);
      vmu1.SetPtEtaPhiM(el1.pt(),el1.eta(),el1.phi(),m);
      vmu2.SetPtEtaPhiM(el2.pt(),el2.eta(),el2.phi(),m);
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
      ParticleMass Kstar_m = 0.896;
      TripleTrackKinFit fitter(allParticles,Kstar_m);
      if (!fitter.success()) continue;

      float ChiProb= ChiSquaredProbability(fitter.chi(),fitter.dof());
      if (ChiProb<Pchi2BMuMuK) continue;

      NRbks_k_sdxy.push_back(KstarTrack.at(ik).first->track().dxy(vertex_point)/KstarTrack.at(ik).first->track().dxyError());      
      NRbks_pi_sdxy.push_back(KstarTrack.at(ik).second->track().dxy(vertex_point)/KstarTrack.at(ik).second->track().dxyError());
      NRbks_mass.push_back(fitter.Mother_Mass(RefitTracks)); 

      std::vector<float> tempBpt; 
      tempBpt.push_back(fitter.Mother_Momentum(RefitTracks).perp());
      tempBpt.push_back(fitter.Mother_Momentum(RefitTracks).eta());
      tempBpt.push_back(fitter.Mother_Momentum(RefitTracks).phi());
      NRbks_pt_eta_phi.push_back(tempBpt);
      NRbks_charge.push_back(fitter.Mother_Charge());
      std::vector<float> tempBx; 
      tempBx.push_back(fitter.Mother_XYZ().x());
      tempBx.push_back(fitter.Mother_XYZ().y());
      tempBx.push_back(fitter.Mother_XYZ().z());
      NRbks_x_y_z.push_back(tempBx);
      std::vector<float> tempBex;
      tempBex.push_back(fitter.Mother_XYZError().cxx());
      tempBex.push_back(fitter.Mother_XYZError().cyy());
      tempBex.push_back(fitter.Mother_XYZError().czz());
      NRbks_ex_ey_ez.push_back(tempBex); 
      std::vector<float> tempBept;
      tempBept.push_back(fitter.Mother_PtError());
      tempBept.push_back(fitter.Mother_EtaError());
      tempBept.push_back(fitter.Mother_PhiError());
      NRbks_ept_eeta_ephi.push_back(tempBept);
      NRbks_chi_prob.push_back(ChiProb); 
      if(IsE) NRbks_mudecay.push_back(0);
      else  NRbks_mudecay.push_back(1);
	 
      NRbks_lep1Id.push_back(imu1); NRbks_lep2Id.push_back(imu2);
      GlobalPoint Dispbeamspot(-1*((theBeamSpot->x0()-fitter.Mother_XYZ().x())+(fitter.Mother_XYZ().z()-theBeamSpot->z0()) * theBeamSpot->dxdz()),-1*((theBeamSpot->y0()-fitter.Mother_XYZ().y())+ (fitter.Mother_XYZ().z()-theBeamSpot->z0()) * theBeamSpot->dydz()), 0);     
      NRbks_bspot_lxy.push_back( Dispbeamspot.perp());
      NRbks_bspot_elxy.push_back(fitter.Mother_XYZError().rerr(Dispbeamspot));
     
      math::XYZVector pperp(fitter.Mother_Momentum(true).x(),fitter.Mother_Momentum(true).y(),0);
      math::XYZVector vperp(Dispbeamspot.x(),Dispbeamspot.y(),0.);
      NRbks_cosTheta2D.push_back(vperp.Dot(pperp)/(vperp.R()*pperp.R()));
//      TLorentzVector refl1,refl2;
    /* for(unsigned int ichild=0; ichild<allParticles.size(); ichild++){
       std::vector<float> temp;
       temp.push_back(fitter.Daughter_Momentum(ichild).perp());
       temp.push_back(fitter.Daughter_Momentum(ichild).eta());
       temp.push_back(fitter.Daughter_Momentum(ichild).phi());
       temp.push_back(fitter.Daughter_Charge(ichild));
   
       if(ichild==2)
           NRbks_Kpt_eta_phi.push_back(temp);
       else if (ichild==0){
          NRbks_l1pt_eta_phi.push_back(temp);
          refl1.SetPtEtaPhiM(temp[0],temp[1],temp[2],m);     }
       else{
          NRbks_l2pt_eta_phi.push_back(temp);    
          refl2.SetPtEtaPhiM(temp[0],temp[1],temp[2],m);     }
       
     }*/
//     NRbks_mll.push_back((refl1+refl2).M()); NRbks_vtx_index.push_back(-1);
    //isolation
  /*  double temp_iso04=0,temp_iso08=0;
    for(unsigned int ikIso=0; ikIso<KTrack.size(); ikIso++){
        Track temptrk=KTrack[ikIso]->track();
        if (DR(tempBpt[1],tempBpt[2],temptrk.eta(),temptrk.phi())<0.4){
           temp_iso04+=temptrk.pt(); 
           }
        if (DR(tempBpt[1],tempBpt[2],temptrk.eta(),temptrk.phi())<0.8){
            temp_iso08+=temptrk.pt(); 
           }
       }
     NRb_iso04.push_back( temp_iso04); NRb_iso08.push_back( temp_iso08); 
  */
   }


 }
} 
 
 
if (NRb_mass.size()==0 && SkipEventWithNoBToMuMuK) return ;

const pat::MET &theMet = met->front();
    ptmet=theMet.et();
    phimet=theMet.phi();

  for (const pat::Jet &jet : *jets){

    if (fabs(jet.eta())>2.5) continue;
    jet_pt.push_back(jet.pt()); jet_eta.push_back(jet.eta());
    jet_phi.push_back(jet.phi());
    jet_cEmEF.push_back(jet.chargedEmEnergyFraction());
    jet_cHEF.push_back(jet.chargedHadronEnergyFraction());
    jet_cHMult.push_back(jet.chargedHadronMultiplicity());
    jet_cMuEF.push_back(jet.chargedMuEnergyFraction());
    jet_cMult.push_back(jet.chargedMultiplicity());
    jet_MuEF.push_back(jet.muonEnergyFraction());
    jet_eEF.push_back(jet.electronEnergyFraction());
    jet_nEmEF.push_back(jet.neutralEmEnergyFraction());
    jet_nHEF.push_back(jet.neutralHadronEnergyFraction());
    jet_nMult.push_back(jet.neutralMultiplicity());
    jet_pEF.push_back(jet.photonEnergyFraction());
    njets++;
  }

  if(saveL1){
      l1objects=L1Analyze(iEvent, iSetup);
      l1muon_pt=l1objects.second[0];    l1muon_eta=l1objects.second[1];
      l1muon_phi=l1objects.second[2];   l1muon_iso=l1objects.second[3];
      l1muon_qual=l1objects.second[4];  l1jet_pt=l1objects.second[5];
      l1jet_eta=l1objects.second[6];    l1jet_phi=l1objects.second[7];
      l1jet_iso=l1objects.second[8];    l1jet_qual=l1objects.second[9];
      l1met=l1objects.first[0];         l1met_phi=l1objects.first[1]; 
      l1ht=l1objects.first[2];          l1hf_met=l1objects.first[3];
      l1hf_met_phi=l1objects.first[4];
      l1seeds=L1SeedAnalyze(iEvent,algoBitToName,Seed_);
      l1_seed1=l1seeds[0].second; l1_seed2=l1seeds[1].second; 
      l1_seed3=l1seeds[2].second; l1_seed4=l1seeds[3].second; 
      l1_seed5=l1seeds[4].second; l1_seed6=l1seeds[5].second;
   }
  if (saveHLT){
  
    std::pair<std::vector<float>,std::vector<std::vector<std::vector<float>>>> trgresult=HLTAnalyze(iEvent,iSetup,HLTPath_);     
     trigger1=trgresult.first[0]; trigger2=trgresult.first[1]; trigger3=trgresult.first[2];  
     trigger4=trgresult.first[3]; trigger5=trgresult.first[4]; trigger6=trgresult.first[5];   
     trigger7=trgresult.first[6]; trigger8=trgresult.first[7];
      if(trigger1==1) tr1_obj_pt_eta_phi=trgresult.second[0];  
      if(trigger2==1) tr2_obj_pt_eta_phi=trgresult.second[1];  
      if(trigger3==1) tr3_obj_pt_eta_phi=trgresult.second[2];   
      if(trigger4==1) tr4_obj_pt_eta_phi=trgresult.second[3]; 
      if(trigger5==1) tr5_obj_pt_eta_phi=trgresult.second[4]; 
      if(trigger6==1) tr6_obj_pt_eta_phi=trgresult.second[5];  
      if(trigger7==1) tr7_obj_pt_eta_phi=trgresult.second[6];
      if(trigger8==1) tr8_obj_pt_eta_phi=trgresult.second[7];
    }
    totb1+=trigger2; totb2+=trigger3; totb3+=trigger5;
    //cout<<"B trg1 "<<trigger1<<" trg2 "<<trigger2<<" trg3 "<<trigger3<<endl;
  t1->Fill();


}


// ------------ method called once each job just before starting event loop  ------------
template<typename T1>
void 
TriggerAnalyzer<T1>::beginJob()
{
  t1=fs->make<TTree>("mytree","mytree");
  t1->Branch("event",&event); t1->Branch("run_number",&run_number);
  t1->Branch("ls",&ls);
  t1->Branch("HLT_path1",&trigger1); t1->Branch("HLT_path2",&trigger2);
  t1->Branch("HLT_path3",&trigger3); t1->Branch("HLT_path4",&trigger4);
  t1->Branch("HLT_path5",&trigger5); t1->Branch("HLT_path6",&trigger6);
  t1->Branch("HLT_path7",&trigger5); t1->Branch("HLT_path8",&trigger6);
  t1->Branch("L1_seed1",&l1_seed1); t1->Branch("L1_seed2",&l1_seed2);
  t1->Branch("L1_seed3",&l1_seed3); t1->Branch("L1_seed4",&l1_seed4);
  t1->Branch("L1_seed5",&l1_seed5); t1->Branch("L1_seed6",&l1_seed6);
  t1->Branch("met_pt",&ptmet); t1->Branch("met_phi",&phimet);
  t1->Branch("nmuon",&nmuons); t1->Branch("muon_pt",&muon_pt);
  t1->Branch("muon_eta",&muon_eta); t1->Branch("muon_phi",&muon_phi);
  t1->Branch("muon_charge",&muon_charge); t1->Branch("muon_qual",&muon_qual);
  t1->Branch("muon_dxy",&muon_dxy); t1->Branch("muon_dz",&muon_dz);
  t1->Branch("muon_edxy",&muon_edxy); t1->Branch("muon_edz",&muon_edz);
  t1->Branch("muon_d0",&muon_d0); t1->Branch("muon_ed0",&muon_ed0);
  t1->Branch("muon_vx",&muon_vx); t1->Branch("muon_vy",&muon_vy);
  t1->Branch("muon_vz",&muon_vz); t1->Branch("muon_iso",&muon_iso);
  t1->Branch("muon_global_flag",&muon_global_flag);
  t1->Branch("muon_tracker_flag",&muon_tracker_flag);
  t1->Branch("muon_standalone_flag",&muon_standalone_flag);
  t1->Branch("muon_soft",&muon_soft); t1->Branch("muon_loose",&muon_loose);
  t1->Branch("muon_medium",&muon_medium);  t1->Branch("muon_tight",&muon_tight);
  t1->Branch("jet_pt",&jet_pt); t1->Branch("jet_eta",&jet_eta);
  t1->Branch("jet_phi",&jet_phi); t1->Branch("jet_cEmEF",&jet_cEmEF);
  t1->Branch("jet_cHEF",&jet_cHEF); t1->Branch("jet_cHMult",&jet_cHMult);
  t1->Branch("jet_cMuEF",&jet_cMuEF); t1->Branch("jet_cMult",&jet_cMult);
  t1->Branch("jet_MuEF",&jet_MuEF); t1->Branch("jet_nEmEF",&jet_nEmEF);
  t1->Branch("jet_nHEF",&jet_nHEF); t1->Branch("jet_nMult",&jet_nMult);
  t1->Branch("jet_eEF",&jet_eEF); t1->Branch("jet_pEF",&jet_pEF);
  t1->Branch("njets",&njets); t1->Branch("nelectron",&nel);
  t1->Branch("el_pt",&el_pt); t1->Branch("el_eta",&el_eta);
  t1->Branch("el_phi",&el_phi); t1->Branch("el_charge",&el_charge);
  t1->Branch("el_dxy",&el_dxy); t1->Branch("el_dz",&el_dz);
  t1->Branch("el_edxy",&el_edxy); t1->Branch("el_edz",&el_edz);
  t1->Branch("el_vx",&el_vx); t1->Branch("el_vy",&el_vy);
  t1->Branch("el_vz",&el_vz); t1->Branch("el_mva_out",&el_mva_out);
  t1->Branch("el_mva_iso",&el_mva_iso); t1->Branch("el_iso",&el_iso);
  t1->Branch("el_veto",&el_veto); t1->Branch("el_soft",&el_soft);
  t1->Branch("el_medium",&el_medium); t1->Branch("el_tight",&el_tight);
  t1->Branch("el_mva_map_value",&el_mva_map_value);
  t1->Branch("HLT1Obj_pt_eta_phi_charge",&tr1_obj_pt_eta_phi);
  t1->Branch("HLT2Obj_pt_eta_phi_charge",&tr2_obj_pt_eta_phi);
  t1->Branch("HLT3Obj_pt_eta_phi_charge",&tr3_obj_pt_eta_phi);
  t1->Branch("HLT4Obj_pt_eta_phi_charge",&tr4_obj_pt_eta_phi);
  t1->Branch("HLT5Obj_pt_eta_phi_charge",&tr5_obj_pt_eta_phi);
  t1->Branch("HLT6Obj_pt_eta_phi_charge",&tr6_obj_pt_eta_phi);
  t1->Branch("HLT7Obj_pt_eta_phi_charge",&tr7_obj_pt_eta_phi);
  t1->Branch("HLT8Obj_pt_eta_phi_charge",&tr8_obj_pt_eta_phi);
  t1->Branch("l1muon_pt",&l1muon_pt); t1->Branch("l1muon_eta",&l1muon_eta);
  t1->Branch("l1muon_phi",&l1muon_phi); t1->Branch("l1muon_iso",&l1muon_iso);
  t1->Branch("l1muon_qual",&l1muon_qual);
  t1->Branch("l1jet_pt",&l1jet_pt); t1->Branch("l1jet_eta",&l1jet_eta);
  t1->Branch("l1jet_phi",&l1jet_phi); t1->Branch("l1jet_iso",&l1jet_iso);
  t1->Branch("l1jet_qual",&l1jet_qual); t1->Branch("l1met",&l1met);
  t1->Branch("l1met_phi",&l1met_phi); t1->Branch("l1ht",&l1ht);
  t1->Branch("l1hf_met",&l1hf_met); t1->Branch("l1hf_met_phi",&l1hf_met_phi);
  t1->Branch("ngenB",&ngenB); t1->Branch("genB_pt",&genB_pt);
  t1->Branch("genB_phi",&genB_phi); t1->Branch("genB_eta",&genB_eta);
  t1->Branch("genB_pdgId",&genB_pdgId); t1->Branch("genB_Bindex",&genB_Bindex);
  t1->Branch("genB_daughter_pt",&genB_daughter_pt);
  t1->Branch("genB_daughter_eta",&genB_daughter_eta);
  t1->Branch("genB_daughter_phi",&genB_daughter_phi);
  t1->Branch("genB_daughter_pdgId",&genB_daughter_pdgId);
  t1->Branch("genB_daughter_Bindex",&genB_daughter_Bindex);
  t1->Branch("genB_daughter_Dindex",&genB_daughter_Dindex);
  t1->Branch("genB_granddaughter_pt",&genB_granddaughter_pt);
  t1->Branch("genB_granddaughter_eta",&genB_granddaughter_eta);
  t1->Branch("genB_granddaughter_phi",&genB_granddaughter_phi);
  t1->Branch("genB_granddaughter_pdgId",&genB_granddaughter_pdgId);
  t1->Branch("genB_granddaughter_Bindex",&genB_granddaughter_Bindex);
  t1->Branch("genB_granddaughter_Dindex",&genB_granddaughter_Dindex);
  t1->Branch("ngenLep",&ngenLep); t1->Branch("genLep_pt",&genLep_pt);
  t1->Branch("genLep_phi",&genLep_phi); t1->Branch("genLep_eta",&genLep_eta);
  t1->Branch("genLep_pdgId",&genLep_pdgId); t1->Branch("genLep_mom",&genLep_mom);


  //  t1->Branch("muon_DCA",&muon_DCA);
t1->Branch("muon_trkpt",&muon_trkpt);
t1->Branch("muon_trketa",&muon_trketa);
t1->Branch("muon_trkphi",&muon_trkphi);
t1->Branch("el_trkpt",&el_trkpt);
t1->Branch("el_trketa",&el_trketa);
t1->Branch("el_trkphi",&el_trkphi);
  
  t1->Branch("pvertex_x",&pvertex_x);
  t1->Branch("pvertex_y",&pvertex_y);
  t1->Branch("pvertex_z",&pvertex_z);
  t1->Branch("pvertex_ex",&pvertex_ex);
  t1->Branch("pvertex_ey",&pvertex_ey);
  t1->Branch("pvertex_ez",&pvertex_ez);
  t1->Branch("ntracks",&ntracks);
  t1->Branch("track_pt",&track_pt);
  t1->Branch("track_eta",&track_eta);
  t1->Branch("track_phi",&track_phi);
  t1->Branch("track_norm_chi2",&track_norm_chi2);
  t1->Branch("track_charge",&track_charge);
  t1->Branch("track_dxy",&track_dxy);
  t1->Branch("track_dz",&track_dz);
  t1->Branch("track_validhits",&track_validhits);
  t1->Branch("track_losthits",&track_losthits);
  t1->Branch("track_fromPV",&track_fromPV);
  t1->Branch("track_highPurity",&track_highPurity);
t1->Branch("NRb_pt_eta_phi",&NRb_pt_eta_phi);
t1->Branch("NRb_mass",&NRb_mass);
t1->Branch("NRb_chi_prob",&NRb_chi_prob);
t1->Branch("NRb_x_y_z",&NRb_x_y_z);
t1->Branch("NRb_ex_ey_ez",&NRb_ex_ey_ez);
t1->Branch("NRb_ept_eeta_ephi",&NRb_ept_eeta_ephi);
t1->Branch("NRb_mudecay",&NRb_mudecay);
t1->Branch("NRb_Kpt_eta_phi_charge",&NRb_Kpt_eta_phi);
t1->Branch("NRb_KUNFITpt_eta_phi",&NRb_KUNFITpt_eta_phi);
t1->Branch("NRb_charge",&NRb_charge);
t1->Branch("NRb_l1pt_eta_phi_charge",&NRb_l1pt_eta_phi);
t1->Branch("NRb_l2pt_eta_phi_charge",&NRb_l2pt_eta_phi);
t1->Branch("NRb_lep1Id",&NRb_lep1Id);
t1->Branch("NRb_lep2Id",&NRb_lep2Id);
t1->Branch("NRb_mll",&NRb_mll);
t1->Branch("NRb_ll_prob",&NRb_ll_prob);
t1->Branch("NRb_trk_sdxy",&NRb_trk_sdxy);
t1->Branch("NRb_vtx_index",&NRb_vtx_index);
t1->Branch("NRb_trk_chi_norm",&NRb_trk_chi_norm);
 t1->Branch("NRb_bspot_lxy",&NRb_bspot_lxy);
 t1->Branch("NRb_bspot_elxy",&NRb_bspot_elxy);
 t1->Branch("NRb_llsvbsv",&NRb_llsvbsv);
 t1->Branch("NRb_cosTheta2D",&NRb_cosTheta2D);
t1->Branch("NRb_iso04",&NRb_iso04);
t1->Branch("NRb_iso08",&NRb_iso08);
t1->Branch("NRb_biso02",&NRb_biso02);
t1->Branch("NRb_biso04",&NRb_biso04);
t1->Branch("NRb_biso06",&NRb_biso06);
t1->Branch("NRb_biso1p2",&NRb_biso1p2);
t1->Branch("NRbks_k_sdxy",&NRbks_k_sdxy);
t1->Branch("NRbks_pi_sdxy",&NRbks_pi_sdxy);
t1->Branch("NRbks_mass",&NRbks_mass);
t1->Branch("NRbks_charge",&NRbks_charge);
t1->Branch("NRbks_chi_prob",&NRbks_chi_prob);
t1->Branch("NRbks_bspot_lxy",&NRbks_bspot_lxy);
t1->Branch("NRbks_bspot_elxy",&NRbks_bspot_elxy);
t1->Branch("NRbks_cosTheta2D",&NRbks_cosTheta2D);
t1->Branch("NRbks_pt_eta_phi",&NRbks_pt_eta_phi);
t1->Branch("NRbks_x_y_z",&NRbks_x_y_z);
t1->Branch("NRbks_ept_eeta_ephi",&NRbks_ept_eeta_ephi);
t1->Branch("NRbks_ex_ey_ez",&NRbks_ex_ey_ez);
t1->Branch("NRbks_mudecay",&NRbks_mudecay);
t1->Branch("NRbks_lep1Id",&NRbks_lep1Id);
t1->Branch("NRbks_lep2Id",&NRbks_lep2Id);

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
