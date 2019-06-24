#ifndef HLTL1TREE_H
#define HLTL1TREE_H
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "HLTrigger/HLTcore/interface/defaultModuleLabel.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include <vector>
#include "TString.h"
#include "NtupleContent.h"

class HLTL1tree{

 public:
  HLTL1tree(const edm::Event& iEvent,edm::EDGetTokenT<GlobalAlgBlkBxCollection> & l1resultToken_, edm::EDGetToken & l1MuonsToken_,edm::EDGetToken & l1JetsToken_,edm::EDGetToken & l1MetToken_,edm::EDGetTokenT<edm::TriggerResults> & trgresultsToken_,  edm::EDGetTokenT<vector<pat::TriggerObjectStandAlone>> & trigobjectsToken_);
  ~HLTL1tree();
  void L1objetcs(NtupleContent& nt);
  void L1trigger(TString * algoBitToName,std::vector<std::string> &seed);
  void HLTtrigger(std::vector<std::string>& HLTPaths);
  std::vector<bool> GetL1Decision() {return l1seeds;}
  std::vector<bool> GetHLTDecision() {return hltpaths;}
  void HLTobjects(std::vector<std::string>& HLTPath);
  std::vector<std::vector<float>> GetHLTObjects(unsigned int i){ 
    return HLTObjects[i];
  }
  std::vector<float> GetHighestPtHLTObject();
  bool HLTPathFire() {return evtFire;}
  void FillL1(NtupleContent& nt){
	  nt.l1_seed1=(l1seeds[0]) ? 1:0;
      nt.l1_seed2=(l1seeds[1]) ? 1:0;
      nt.l1_seed3=(l1seeds[2]) ? 1:0;
      nt.l1_seed4=(l1seeds[3]) ? 1:0;
      nt.l1_seed5=(l1seeds[4]) ? 1:0;
      nt.l1_seed6=(l1seeds[5]) ? 1:0;
      }
  void FillHLT(NtupleContent& nt){
	  nt.trigger1=(hltpaths[0]) ? 1:0;
      nt.trigger2=(hltpaths[1]) ? 1:0;
      nt.trigger3=(hltpaths[2]) ? 1:0;
      nt.trigger4=(hltpaths[3]) ? 1:0;
      nt.trigger5=(hltpaths[4]) ? 1:0;
      nt.trigger6=(hltpaths[5]) ? 1:0;
      nt.trigger7=(hltpaths[6]) ? 1:0;
      nt.trigger8=(hltpaths[7]) ? 1:0;
      }
   void FillObj(NtupleContent & nt){
	  if (hltpaths[0]) nt.tr1_obj_pt_eta_phi=GetHLTObjects(0);
      if (hltpaths[1]) nt.tr2_obj_pt_eta_phi=GetHLTObjects(1);
      if (hltpaths[2]) nt.tr3_obj_pt_eta_phi=GetHLTObjects(2);
      if (hltpaths[3]) nt.tr4_obj_pt_eta_phi=GetHLTObjects(3);
      if (hltpaths[4]) nt.tr5_obj_pt_eta_phi=GetHLTObjects(4);
      if (hltpaths[5]) nt.tr6_obj_pt_eta_phi=GetHLTObjects(5);
      if (hltpaths[6]) nt.tr7_obj_pt_eta_phi=GetHLTObjects(6);
      if (hltpaths[7]) nt.tr8_obj_pt_eta_phi=GetHLTObjects(7);
	   }
   

 private:
  edm::Handle<GlobalAlgBlkBxCollection> l1result;
  edm::Handle<l1t::MuonBxCollection> l1Muons;
  edm::Handle<l1t::JetBxCollection> l1Jets;
  edm::Handle<BXVector<l1t::EtSum> > l1Met;
  edm::Handle<vector<pat::TriggerObjectStandAlone>> triggerObjects;
  edm::Handle<edm::TriggerResults> trigResults;
  bool l1valid; std::vector<bool> l1seeds; std::vector<bool> hltpaths;
  edm::TriggerNames trigName; std::vector<std::vector<std::vector<float>>> HLTObjects;
  std::vector<std::vector<float>> ActiveTrgObject;
  bool evtFire=false;
};



#endif
