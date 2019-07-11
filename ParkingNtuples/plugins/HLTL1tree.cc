#include "HLTL1tree.h"

HLTL1tree::HLTL1tree(const edm::Event& iEvent,edm::EDGetTokenT<GlobalAlgBlkBxCollection> & l1resultToken_, 
                           edm::EDGetToken & l1MuonsToken_,
                           edm::EDGetTokenT<edm::TriggerResults> & trgresultsToken_,
                           edm::EDGetTokenT<vector<pat::TriggerObjectStandAlone>> & trigobjectsToken_
                    )
{
  iEvent.getByToken(l1resultToken_,l1result);
  iEvent.getByToken(l1MuonsToken_,l1Muons);
  iEvent.getByToken(trigobjectsToken_ ,triggerObjects);
  iEvent.getByToken(trgresultsToken_, trigResults);

  trigName = iEvent.triggerNames(*trigResults);
  l1valid  = l1result.isValid();
}

HLTL1tree::~HLTL1tree(){}

void HLTL1tree::L1objects(NtupleContent& nt){
  for(typename std::vector< l1t::Muon >::const_iterator mu=l1Muons->begin(0); mu !=l1Muons->end(0); mu++){
    nt.l1muon_pt  .push_back(mu->et()); 
    nt.l1muon_eta .push_back(mu->eta());
    nt.l1muon_phi .push_back(mu->phi()); 
    nt.l1muon_qual.push_back(mu->hwQual());
  }
}

void HLTL1tree::L1trigger(TString * algoBitToName,std::vector<std::string> & seed){

  if (!l1valid) {
    for (unsigned int iseed = 0; iseed < seed.size(); iseed++)
      l1seeds.push_back(false);
    return;
  }

  GlobalAlgBlk const &result = l1result->at(0, 0);
  for (unsigned int iseed = 0; iseed < seed.size(); iseed++){
    bool sfire=false;
    for (unsigned int itrg = 0; itrg < result.maxPhysicsTriggers; ++itrg){
      if (result.getAlgoDecisionFinal(itrg)!=1) continue;
      std::string l1trigName = static_cast<const char *>(algoBitToName[itrg]);
      if (l1trigName!=seed[iseed]) continue;
      sfire=true;
    }
    if (sfire) l1seeds.push_back(true);
    else       l1seeds.push_back(false);
  }

}


void HLTL1tree::HLTtrigger(std::vector<std::string>& HLTPaths,
                           HLTPrescaleProvider & hltPrescaleProvider_
                          ){
  if( trigResults.failedToGet() ){
    for (unsigned int trg=0; trg < HLTPaths.size(); trg++)
      hltpaths.push_back(false);
    return;
  }
  int N_Triggers = trigResults->size();
 
  bool fire = false;
  for (unsigned int itrg=0; itrg < HLTPaths.size(); itrg++){
    fire = false;
    for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
      TString TrigPath =trigName.triggerName(i_Trig);
      if (TrigPath.Contains(HLTPaths[itrg])  &&                                   // trigger path is included in our list
          trigResults->accept(i_Trig)        &&                                   // trigger path is fired
          hltPrescaleProvider_.rejectedByHLTPrescaler(*trigResults, i_Trig) == 0  // trigger path is not rejected by applying prescale
         ){ 
         fire    = true; 
         evtFire = true;
      }
    }
    if (fire) hltpaths.push_back(true);
    else      hltpaths.push_back(false);
  }
}


void HLTL1tree::HLTobjects(std::vector<std::string>& HLTPath){
  for (unsigned int ifilter=0; ifilter<HLTPath.size(); ifilter++){
    std::vector<std::vector<float>> tot_trgobj_PtEtaPhi;
    for(pat::TriggerObjectStandAlone itrg :*triggerObjects){
      if(!itrg.id(83)) continue;
      bool save=false;
      itrg.unpackPathNames(trigName);
      std::vector<std::string> const & pathnames = itrg.pathNames(true,true);
      for(std::string name : pathnames){
        if (name.find(HLTPath[ifilter]) != std::string::npos) save=true;
      }
      std::vector<float> trgobj_PtEtaPhi;
      if(!save) continue;
      trgobj_PtEtaPhi.push_back(itrg.pt());
      trgobj_PtEtaPhi.push_back(itrg.eta());
      trgobj_PtEtaPhi.push_back(itrg.phi());
      trgobj_PtEtaPhi.push_back(itrg.charge());
      tot_trgobj_PtEtaPhi.push_back(trgobj_PtEtaPhi);
      ActiveTrgObject.push_back(trgobj_PtEtaPhi);
   }
   HLTObjects.push_back(tot_trgobj_PtEtaPhi);
  }
}





std::vector<float> HLTL1tree::GetHighestPtHLTObject(){
  if (ActiveTrgObject.size()==0){
   std::vector<float> def; def.push_back(-1);
   ActiveTrgObject.push_back(def);
  } else
    std::sort(ActiveTrgObject.begin(),ActiveTrgObject.end(),[](const std::vector<float>& p1, const std::vector<float>& p2) { return p1[0] > p2[0]; });

  return ActiveTrgObject[0];
}
