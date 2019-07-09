#include "PFelCollection.h"


PFelCollection::PFelCollection(edm::EDGetToken & electrons_,edm::EDGetTokenT<edm::ValueMap<bool>> & soft_id_,edm::EDGetTokenT<edm::ValueMap<bool>> & medium_id_ ,edm::EDGetTokenT<edm::ValueMap<bool>> & tight_id_, const edm::Event& iEvent){
   iEvent.getByToken(electrons_,electrons); iEvent.getByToken(soft_id_,soft_id);
   iEvent.getByToken(medium_id_,medium_id); iEvent.getByToken(tight_id_,tight_id);
}

PFelCollection::~PFelCollection(){}

void PFelCollection::AddElectrons(NtupleContent & nt){
  trigger::size_type eindex=-1;
  for (const pat::Electron & el : *electrons){
    eindex++;
    if (!el.passConversionVeto()) continue;
    nt.pfel_pt.push_back(el.pt());  nt.pfel_eta.push_back(el.eta());
    nt.pfel_phi.push_back(el.phi()); nt.pfel_vx.push_back(el.vx()); 
    nt.pfel_vy.push_back(el.vy()); nt.pfel_vz.push_back(el.vz());
    nt.pfel_charge.push_back(el.charge());
    const edm::Ptr<pat::Electron> elePtr(electrons,eindex);
    nt.pfel_soft.push_back((*soft_id)[elePtr]);
    nt.pfel_medium.push_back((*medium_id)[elePtr]);
    nt.pfel_tight.push_back((*tight_id)[elePtr]);
  }
}
