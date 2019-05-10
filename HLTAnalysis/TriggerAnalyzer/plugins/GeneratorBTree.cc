#include "GeneratorBTree.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"


GeneratorBTree::GeneratorBTree(edm::EDGetTokenT<edm::View<reco::GenParticle> > & prunedGenToken_,edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > & packedGenToken_,const edm::Event& iEvent){  
   iEvent.getByToken(prunedGenToken_,pruned);
   iEvent.getByToken(packedGenToken_,packed);
  
   finalLeptonsTracks( pruned, packed); float IB=0;
   for(const reco::GenParticle & gen: *pruned){
     if((gen.pdgId()%1000)/100!=5 && (gen.pdgId()%10000)/1000!=5 && (gen.pdgId()%1000)/100!=-5 && (gen.pdgId()%10000)/1000!=-5) continue;
     if (isExcitedB(gen)) continue;
      TLorentzVector vel1,vel2,vk;
      BMotherPt.push_back(gen.pt()); BMotherEta.push_back(gen.eta()); 
      BMotherPhi.push_back(gen.phi()); BMotherPdgId.push_back(gen.pdgId()); 
      BMotherBid.push_back(IB); BMotherM.push_back(gen.mass());
      float ID=0;
      for (unsigned int idaughter=0; idaughter<gen.numberOfDaughters(); idaughter++){
        const reco::Candidate & daught=*(gen.daughter(idaughter));
        BDaughterPt.push_back(daught.pt());  BDaughterEta.push_back(daught.eta());
        BDaughterPhi.push_back(daught.phi()); BDaughterPdgId.push_back(daught.pdgId());
        BDaughterBid.push_back(IB); BDaughterDid.push_back(ID);
       
       if (fabs(daught.pdgId())==313 && fabs(gen.pdgId())==511 && daught.numberOfDaughters()==0 && genKPt.size()==genPiPt.size() ){
          for (unsigned int iks=0; iks<genKstarPt.size(); iks++){
            if (daught.pdgId()!=genKstarId[iks] || daught.pt()!=genKstarPt[iks] || daught.eta()!=genKstarEta[iks] || daught.phi()!=genKstarPhi[iks] || daught.mass()!=genKstarM[iks]) continue;
            BGDaughterPt.push_back(genKPt[iks]); BGDaughterPt.push_back(genPiPt[iks]);
            BGDaughterEta.push_back(genKEta[iks]); BGDaughterEta.push_back(genPiEta[iks]);
            BGDaughterPhi.push_back(genKPhi[iks]); BGDaughterPhi.push_back(genPiPhi[iks]); 
            BGDaughterPdgId.push_back(genKId[iks]); BGDaughterPdgId.push_back(genPiId[iks]);
            BGDaughterBid.push_back(IB); BGDaughterBid.push_back(IB);
            BGDaughterDid.push_back(ID); BGDaughterDid.push_back(ID);
         }
       }
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


void GeneratorBTree::finalLeptonsTracks(edm::Handle<edm::View<reco::GenParticle> > & pruned, edm::Handle<edm::View<pat::PackedGenParticle> > & packed){
  for(const pat::PackedGenParticle & gen: *packed){
     if (fabs(gen.pdgId())==321 ){
         const reco::Candidate * motherInPrunedCollection =gen.mother(0);
         if(fabs(motherInPrunedCollection->pdgId())==313 && fabs(motherInPrunedCollection->mother()->pdgId())==511){   
             genKstarId.push_back(motherInPrunedCollection->pdgId());
             genKstarM.push_back(motherInPrunedCollection->mass()); 
             genKstarPt.push_back(motherInPrunedCollection->pt());
             genKstarEta.push_back(motherInPrunedCollection->eta());
             genKstarPhi.push_back(motherInPrunedCollection->phi());
             genKId.push_back(gen.pdgId()); genKM.push_back(gen.mass()); 
             genKPt.push_back(gen.pt()); genKEta.push_back(gen.eta());
             genKPhi.push_back(gen.phi());
          }
      }
      if (fabs(gen.pdgId())==211){
        const reco::Candidate * motherInPrunedCollection =gen.mother(0);
        if(fabs(motherInPrunedCollection->pdgId())==313 && fabs(motherInPrunedCollection->mother()->pdgId())==511){ 
          genPiId.push_back(gen.pdgId()); genPiM.push_back(gen.mass());
          genPiPt.push_back(gen.pt()); genPiEta.push_back(gen.eta());
          genPiPhi.push_back(gen.phi());
     }
   }
      if (gen.pdgId()!=11 && gen.pdgId()!=-11 && gen.pdgId()!=13 && gen.pdgId()!=-13) continue;
      genLepPt.push_back(gen.pt());  genLepEta.push_back(gen.eta()); 
      genLepPhi.push_back(gen.phi()); genLepPdgId.push_back(gen.pdgId());
      const reco::Candidate * motherInPrunedCollection =gen.mother(0);
      if (motherInPrunedCollection == nullptr)
         genLepMother.push_back(-std::numeric_limits<float>::max()); 
      else 
         genLepMother.push_back(motherInPrunedCollection->pdgId());
 }
}

void GeneratorBTree::GenKLepLepEtaPhi(bool ResonantDecay,int BpdgId,int DaughterLepId,int DaughterHadronId){
 bool found=false;
 for (unsigned int ib=0; ib<BMotherPdgId.size(); ib++){
   if (fabs(BMotherPdgId.at(ib))!=BpdgId) continue;
   unsigned int intmax=std::numeric_limits<unsigned int>::max();
   unsigned int Kindex=intmax,Lepindex=intmax,aLepindex=intmax;
   for (unsigned int idaughter=0; idaughter<BDaughterPt.size(); idaughter++){
     if (BDaughterBid.at(idaughter)!=BMotherBid.at(ib)) continue;
     if (fabs(BDaughterPdgId.at(idaughter))==DaughterHadronId)
        Kindex=idaughter;
     if (!ResonantDecay){
        if (BDaughterPdgId.at(idaughter)==DaughterLepId)
           Lepindex=idaughter;
        if (BDaughterPdgId.at(idaughter)==-1*DaughterLepId)
           aLepindex=idaughter;
     }else if (ResonantDecay && fabs(BDaughterPdgId.at(idaughter))==443){
        for(unsigned int igdaughter=0; igdaughter<BGDaughterPdgId.size(); igdaughter++){
          if (BGDaughterBid.at(igdaughter)!=BMotherBid.at(ib)) continue;
          if (BGDaughterDid.at(igdaughter)!=BDaughterDid.at(idaughter)) continue;
          if (BGDaughterPdgId.at(igdaughter)==DaughterLepId)
             Lepindex=igdaughter;
          if (BGDaughterPdgId.at(igdaughter)==-1*DaughterLepId)
             aLepindex=igdaughter;
        }
     }    
   }//dauthers loop
   if (Kindex==intmax || Lepindex==intmax || aLepindex==intmax)
      continue;
    found=true;
    EtaPhiHadron=std::make_pair(BDaughterEta.at(Kindex),BDaughterPhi.at(Kindex));
    if (!ResonantDecay){
      EtaPhiLep=std::make_pair(BDaughterEta.at(Lepindex),BDaughterPhi.at(Lepindex));
      EtaPhiaLep=std::make_pair(BDaughterEta.at(aLepindex),BDaughterPhi.at(aLepindex));
    } else{
      EtaPhiLep=std::make_pair(BGDaughterEta.at(Lepindex),BGDaughterPhi.at(Lepindex));
      EtaPhiaLep=std::make_pair(BGDaughterEta.at(aLepindex),BGDaughterPhi.at(aLepindex));
    }
    break;
    }//B loop
 
 if (!found){
   EtaPhiHadron=std::make_pair(-10,-10); EtaPhiLep=std::make_pair(-10,-10);
   EtaPhiaLep=std::make_pair(-10,-10);
 }
}//end void
