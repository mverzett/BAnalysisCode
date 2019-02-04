#include "GeneratorBTree.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"


GeneratorBTree::GeneratorBTree(edm::EDGetTokenT<edm::View<reco::GenParticle> > & prunedGenToken_,edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > & packedGenToken_,const edm::Event& iEvent){
  
 
   iEvent.getByToken(prunedGenToken_,pruned);
   iEvent.getByToken(packedGenToken_,packed);
//    std::cout<<"new ev "<<std::endl;

   finalLeptons( pruned, packed);
   float IB=0;
   for(const reco::GenParticle & gen: *pruned){
     if((gen.pdgId()%1000)/100!=5 && (gen.pdgId()%10000)/1000!=5 && (gen.pdgId()%1000)/100!=-5 && (gen.pdgId()%10000)/1000!=-5) continue;
  //   std::cout<<" id "<<gen.pdgId()<<" pt "<<gen.pt()<<" eta "<<gen.eta()<<" phi "<<gen.phi()<<" ndaught "<<gen.numberOfDaughters()<<std::endl;         
     if (isExcitedB(gen)) continue;
     
      BMotherPt.push_back(gen.pt()); BMotherEta.push_back(gen.eta()); 
      BMotherPhi.push_back(gen.phi()); BMotherPdgId.push_back(gen.pdgId()); 
      BMotherBid.push_back(IB);
      float ID=0;
      for (unsigned int idaughter=0; idaughter<gen.numberOfDaughters(); idaughter++){
        const reco::Candidate & daught=*(gen.daughter(idaughter));
        BDaughterPt.push_back(daught.pt());  BDaughterEta.push_back(daught.eta());
        BDaughterPhi.push_back(daught.phi()); BDaughterPdgId.push_back(daught.pdgId());
        BDaughterBid.push_back(IB); BDaughterDid.push_back(ID);
//        std::cout<<"daugther id "<<daught.pdgId()<<" pt "<<daught.pt()<<std::endl;
        
        for (unsigned int igdaughter=0; igdaughter<daught.numberOfDaughters(); igdaughter++){ 
            const reco::Candidate & gdaught=*(daught.daughter(igdaughter));
            BGDaughterPt.push_back(gdaught.pt());  BGDaughterEta.push_back(gdaught.eta());
            BGDaughterPhi.push_back(gdaught.phi()); BGDaughterPdgId.push_back(gdaught.pdgId());
            BGDaughterBid.push_back(IB); BGDaughterDid.push_back(ID);
       }
       ID++;
      }
      IB++;     
   }
  
}


GeneratorBTree::~GeneratorBTree(){}

bool GeneratorBTree::isExcitedB(const reco::GenParticle & gen){
  bool IsB=false;
  for (unsigned int i=0; i<gen.numberOfDaughters(); i++){
      const reco::Candidate & temp=*(gen.daughter(i));
       if( (temp.pdgId()%1000)/100==5 || (temp.pdgId()%10000)/1000==5 || (temp.pdgId()%1000)/100==-5 || (temp.pdgId()%10000)/1000==-5)
           IsB=true;
      }
return IsB;
}

void GeneratorBTree::finalLeptons(edm::Handle<edm::View<reco::GenParticle> > & pruned, edm::Handle<edm::View<pat::PackedGenParticle> > & packed){
  for(const pat::PackedGenParticle & gen: *packed){
    if (gen.pdgId()!=11 && gen.pdgId()!=-11 && gen.pdgId()!=13 && gen.pdgId()!=-13) continue;
    genLepPt.push_back(gen.pt());  genLepEta.push_back(gen.eta()); 
    genLepPhi.push_back(gen.phi()); genLepPdgId.push_back(gen.pdgId());
 //   std::cout<<"lep pt "<<gen.pt()<<" id "<<gen.pdgId()<<std::endl;
    const reco::Candidate * motherInPrunedCollection =gen.mother(0);
    if (motherInPrunedCollection == nullptr){  genLepMother.push_back(-std::numeric_limits<float>::max()); /*std::cout<<"no mom"<<std::endl;*/}
    else{ genLepMother.push_back(motherInPrunedCollection->pdgId());
          /* std::cout<<" mom "<<motherInPrunedCollection->pdgId()<<std::endl;*/}

  }




}




