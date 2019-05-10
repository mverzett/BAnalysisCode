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
