#ifndef GENERATORBTREE_H
#define GENERATORBTREE_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "TLorentzVector.h"
#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include <limits>
#include <cstddef>

class GeneratorBTree{

public:
   GeneratorBTree(edm::EDGetTokenT<edm::View<reco::GenParticle> > & prunedGenToken_,edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > & packedGenToken_, const edm::Event& iEvent);
   virtual ~GeneratorBTree();
   bool isExcitedB(const reco::GenParticle & gen);
   void finalLeptonsTracks(edm::Handle<edm::View<reco::GenParticle> > & pruned, edm::Handle<edm::View<pat::PackedGenParticle> > & packed);
   void GenKLepLepEtaPhi(bool ResonantDecay,int BpdgId,int DaughterLepId,int DaughterHadronId);
   std::vector<float> BMother_Pt(){ return BMotherPt;}
   std::vector<float> BMother_Eta(){ return BMotherEta;}
   std::vector<float> BMother_Phi(){ return BMotherPhi;}
   std::vector<float> BMother_PdgId(){ return BMotherPdgId;}
   std::vector<float> BMother_Bid(){ return BMotherBid;}
   std::vector<float> BDaughter_Pt(){ return BDaughterPt;}
   std::vector<float> BDaughter_Eta(){ return BDaughterEta;}
   std::vector<float> BDaughter_Phi(){ return BDaughterPhi;}
   std::vector<float> BDaughter_PdgId(){ return BDaughterPdgId;}
   std::vector<float> BDaughter_Bid(){ return BDaughterBid;}
   std::vector<float> BDaughter_Did(){ return BDaughterDid;}
   std::vector<float> BGDaughter_Pt(){ return BGDaughterPt;}
   std::vector<float> BGDaughter_Eta(){ return BGDaughterEta;}
   std::vector<float> BGDaughter_Phi(){ return BGDaughterPhi;}
   std::vector<float> BGDaughter_PdgId(){ return BGDaughterPdgId;}
   std::vector<float> BGDaughter_Bid(){ return BGDaughterBid;}
   std::vector<float> BGDaughter_Did(){ return BGDaughterDid;}
   std::vector<float> Lep_Pt(){return genLepPt;}
   std::vector<float> Lep_Eta(){return genLepEta;}
   std::vector<float> Lep_Phi(){return genLepPhi;}
   std::vector<float> Lep_pdgId(){return genLepPdgId;}
   std::vector<float> Lep_Mother(){return genLepMother;}
   std::pair<float,float> EtaPhiK(){return EtaPhiHadron;}
   std::pair<float,float> EtaPhiL(){return EtaPhiLep;}
   std::pair<float,float> EtaPhiaL(){return EtaPhiaLep;}
  

private:
  edm::Handle<edm::View<reco::GenParticle> > pruned;
  edm::Handle<edm::View<pat::PackedGenParticle> > packed;
  std::vector<float> BMotherPt,BMotherEta,BMotherPhi,BMotherPdgId,BMotherBid;
  std::vector<float> BDaughterPt,BDaughterEta,BDaughterPhi,BDaughterPdgId,BDaughterBid,BDaughterDid;
  std::vector<float> BGDaughterPt,BGDaughterEta,BGDaughterPhi,BGDaughterPdgId,BGDaughterBid,BGDaughterDid;
  std::vector<float> genLepPt,genLepEta,genLepPhi,genLepPdgId,genLepMother;
  std::vector<float> genKstarId,genKstarM,genKstarPt,genKstarEta,genKstarPhi;
  std::vector<float> genKId,genKM,genKPt,genKEta,genKPhi;
  std::vector<float> genPiId,genPiM,genPiPt,genPiEta,genPiPhi;
  std::pair<float,float> EtaPhiLep,EtaPhiaLep,EtaPhiHadron;
//std::array<std::vector<float>,4> BDaughters;
};

#endif
