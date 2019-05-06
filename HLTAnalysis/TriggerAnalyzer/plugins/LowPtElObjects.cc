#include "LowPtElObjects.h"


LowPtElObjects::LowPtElObjects(edm::EDGetToken & electrons_,edm::EDGetTokenT<edm::ValueMap<float>> & mva_unbiased_, const edm::Event& iEvent ){
  iEvent.getByToken(electrons_,electrons);
  iEvent.getByToken(mva_unbiased_,mva_unbiased);
}

LowPtElObjects::LowPtElObjects(edm::EDGetToken & electrons_,edm::EDGetToken & gsfTracks_, edm::EDGetTokenT<edm::ValueMap<float>> & mva_unbiased_, const edm::Event& iEvent ){
  iEvent.getByToken(electrons_,electrons);
  iEvent.getByToken(gsfTracks_,gsfTracks);
  iEvent.getByToken(mva_unbiased_,mva_unbiased);
}

LowPtElObjects::~LowPtElObjects(){}

void LowPtElObjects::AddElectrons(NtupleContent & nt){
  for (const pat::Electron & el: *electrons){
    if (!el.passConversionVeto()) continue;
    reco::GsfTrackRef seed = el.gsfTrack();
    if (seed.isNull()) continue;
    nt.lowel_pt.push_back(el.pt()); nt.lowel_eta.push_back(el.eta());  
    nt.lowel_phi.push_back(el.phi()); nt.lowel_charge.push_back(el.charge()); 
    nt.lowel_vx.push_back(el.vx()); nt.lowel_vy.push_back(el.vy());
    nt.lowel_vz.push_back(el.vz());  const reco::Track * trk=el.bestTrack();
    nt.lowel_trkpt.push_back(trk->pt()); nt.lowel_trketa.push_back(trk->eta()); 
    nt.lowel_trkphi.push_back(trk->phi());
    nt.lowel_mva_unbiased.push_back((*mva_unbiased)[seed]);
  }

}
 
void LowPtElObjects::AddGsfTracks(NtupleContent & nt){
  for (const reco::GsfTrack & gsf: *gsfTracks){
    nt.lowgsf_pt.push_back(gsf.pt()); nt.lowgsf_eta.push_back(gsf.eta());  
    nt.lowgsf_phi.push_back(gsf.phi()); nt.lowgsf_charge.push_back(gsf.charge()); 
    nt.lowgsf_vx.push_back(gsf.vx()); nt.lowgsf_vy.push_back(gsf.vy());
    nt.lowgsf_vz.push_back(gsf.vz()); nt.lowgsf_modept.push_back(gsf.ptMode());
    nt.lowgsf_modeeta.push_back(gsf.etaMode()); 
    nt.lowgsf_modephi.push_back(gsf.phiMode());
    const reco::GsfTrackRef seed=edm::Ref(gsfTracks,&gsf-&(gsfTracks->at(0)));
    nt.lowgsf_mva_unbiased.push_back((*mva_unbiased)[seed]);
  }

}
  
