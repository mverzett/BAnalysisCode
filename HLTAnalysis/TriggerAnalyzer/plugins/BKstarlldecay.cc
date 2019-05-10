#include "BKstarlldecay.h"

BKstarlldecay::BKstarlldecay(unsigned int &nmupairs_, std::vector<std::pair<unsigned int,unsigned int>>& used_lep_tracks, std::vector<reco::TransientTrack> & vmuTrack1,std::vector<reco::TransientTrack> & vmuTrack2, std::vector<reco::TransientTrack> & vKTrack, float & beam_xz_,float & beam_yz_ ,bool & RefitMuTracksOnly_) : nmupairs(nmupairs_), used_muTracks_index(used_lep_tracks), muTrack1(vmuTrack1), muTrack2(vmuTrack2), KTrack(vKTrack), beam_xz(beam_xz_),beam_yz(beam_yz_),RefitMuTracksOnly(RefitMuTracksOnly_)
  {
    temp.reserve(4); tempK.reserve(3); tempBpt.reserve(3); tempBx.reserve(3);
    tempBex.reserve(3); tempBept.reserve(3); indexpair.clear();
  }

BKstarlldecay::~BKstarlldecay(){};

void BKstarlldecay::CombineTracks(float & MPairmin,float & MPairmax){
  MPairmin_=MPairmin; MPairmax_=MPairmax; int idx=-1;
  for(std::vector<reco::TransientTrack>::iterator trk1=KTrack.begin(); trk1!=KTrack.end(); ++trk1){
   for(std::vector<reco::TransientTrack>::iterator trk2=trk1+1; trk2!=KTrack.end(); ++trk2){
     if (trk1->charge()==trk2->charge()) continue;  
     idx++;
     vK.SetPtEtaPhiM(trk1->track().pt(),trk1->track().eta(),trk1->track().phi(),0.493);
     vPi.SetPtEtaPhiM(trk2->track().pt(),trk2->track().eta(),trk2->track().phi(),0.139);
     if ( (vK+vPi).M()>MPairmin-0.1 && (vK+vPi).M()<MPairmax+0.1){    
       if (FormKstar(*trk1,*trk2)){
         KstarTracks.emplace_back(std::make_pair(*trk1,*trk2));
         indexpair.push_back(idx);
       }
     }
     vK.SetPtEtaPhiM(trk2->track().pt(),trk2->track().eta(),trk2->track().phi(),0.493);
     vPi.SetPtEtaPhiM(trk1->track().pt(),trk1->track().eta(),trk1->track().phi(),0.139);
     if ( (vK+vPi).M()>MPairmin-0.1 && (vK+vPi).M()<MPairmax+0.1){    
       if (FormKstar(*trk2,*trk1)){
         KstarTracks.emplace_back(std::make_pair(*trk2,*trk1));
         indexpair.push_back(idx);
       }
     }  
   }//trk2
}//trk1
}

bool BKstarlldecay::FormKstar(reco::TransientTrack & trkK,reco::TransientTrack & trkPi ){
    std::vector<RefCountedKinematicParticle> allParticles;
    allParticles.push_back(pFactory.particle(trkK,kaon_mass,chi,ndf,kaon_sigma));
    allParticles.push_back(pFactory.particle(trkPi,pion_mass,chi,ndf,pion_sigma));
    TripleTrackKinFit fitter(allParticles);
    if (fitter.success()) {        
       if (fitter.Mother_Mass(true)>MPairmin_ && fitter.Mother_Mass(true)<MPairmax_)  return true;   
    }
    return false;   
}

void BKstarlldecay::Fill(NtupleContent &nt){
  vperp.SetXYZ(nt.beam_x,nt.beam_y,nt.beam_z);
  for(const std::pair<unsigned int, unsigned int> &imu : used_muTracks_index){
    bool IsE(false); float m(0.105); 
    unsigned int mindex(&imu-&used_muTracks_index[0]);
    if ( mindex>nmupairs || mindex==nmupairs ){ 
      IsE=true; m=0.000511; part_mass=0.000511;
    }
    unsigned int imu1(imu.first); unsigned int imu2(imu.second);
    TLorentzVector vmu1,vmu2;
    if (!IsE){
      vmu1.SetPtEtaPhiM(nt.muon_pt.at(imu1),nt.muon_eta.at(imu1),nt.muon_phi.at(imu1),m);
      vmu2.SetPtEtaPhiM(nt.muon_pt.at(imu2),nt.muon_eta.at(imu2),nt.muon_phi.at(imu2),m); }
    else {
      vmu1.SetPtEtaPhiM(nt.el_pt.at(imu1),nt.el_eta.at(imu1),nt.el_phi.at(imu1),m);
      vmu2.SetPtEtaPhiM(nt.el_pt.at(imu2),nt.el_eta.at(imu2),nt.el_phi.at(imu2),m); }
    for(const std::pair<reco::TransientTrack,reco::TransientTrack> & pairtrk : KstarTracks){
       tempBpt.clear(); tempBx.clear(); tempBex.clear(); tempBept.clear();
       bool RefitTracks=false;  const reco::TransientTrack & trkK=pairtrk.first;
       const reco::TransientTrack & trkPi=pairtrk.second;
       if (RefitMuTracksOnly && !IsE) RefitTracks=true;
       vK.SetPtEtaPhiM(trkK.track().pt(),trkK.track().eta(),trkK.track().phi(),0.495);
       vPi.SetPtEtaPhiM(trkPi.track().pt(),trkPi.track().eta(),trkPi.track().phi(),0.139);
       // if (deltaR(vmu1.Eta(),vmu1.Phi(),vK.Eta(),vK.Phi())<0.005) continue;
       // if (deltaR(vmu2.Eta(),vmu2.Phi(),vK.Eta(),vK.Phi())<0.005) continue;
       // if (deltaR(vmu1.Eta(),vmu1.Phi(),vPi.Eta(),vPi.Phi())<0.005) continue;
       //  if (deltaR(vmu2.Eta(),vmu2.Phi(),vPi.Eta(),vPi.Phi())<0.005) continue;
       if ((vmu1+vmu2+vK+vPi).M()<massmin_ || (vmu1+vmu2+vK+vPi).M()>massmax_) continue;
       KinematicParticleFactoryFromTransientTrack pFactory;
       std::vector<RefCountedKinematicParticle> allParticles;
       allParticles.push_back(pFactory.particle(muTrack1.at(mindex),part_mass,chi,ndf,part_sigma));
       allParticles.push_back(pFactory.particle(muTrack2.at(mindex),part_mass,chi,ndf,part_sigma));
       allParticles.push_back(pFactory.particle(trkK,kaon_mass,chi,ndf,kaon_sigma));
       allParticles.push_back(pFactory.particle(trkPi,pion_mass,chi,ndf,pion_sigma));
      TripleTrackKinFit fitter(allParticles);
      if (!fitter.success()) continue;
      float ChiProb(ChiSquaredProbability(fitter.chi(),fitter.dof()));
      if (ChiProb<probcut_) continue;
      GlobalPoint Dispbeamspot(-1*((nt.beam_x-fitter.Mother_XYZ().x())+(fitter.Mother_XYZ().z()-nt.beam_z)* beam_xz),-1*((nt.beam_y-fitter.Mother_XYZ().y())+ (fitter.Mother_XYZ().z()-nt.beam_z) * beam_yz), 0);              
      if (IsE)
         pperp.SetXYZ((vmu1+vmu2+vK+vPi).Px(),(vmu1+vmu2+vK+vPi).Py(),0);
      else
         pperp.SetXYZ(fitter.Mother_Momentum(RefitTracks).x(),fitter.Mother_Momentum(RefitTracks).y(),0);
      float cos(vperp.Dot(pperp)/(vperp.R()*pperp.R()));
      if (cos<coscut_) continue;   
      nt.NRbks_cosTheta2D.push_back(cos);
      nt.NRbks_bspot_lxy.push_back(Dispbeamspot.perp());
      nt.NRbks_bspot_elxy.push_back(fitter.Mother_XYZError().rerr(Dispbeamspot));         
      nt.NRbks_mass.push_back(fitter.Mother_Mass(RefitTracks));
      nt.NRbks_chi_prob.push_back(ChiProb); 
      nt.NRbks_charge.push_back(fitter.Mother_Charge());

      for ( auto & x : {vK.Pt(),vK.Eta(),vK.Phi()}) tempK.push_back(x);
      nt.NRbks_KUNFITpt_eta_phi.push_back(tempK);
    
      for ( auto & x : {vPi.Pt(),vPi.Eta(),vPi.Phi()}) tempK.push_back(x);
      nt.NRbks_PiUNFITpt_eta_phi.push_back(tempK);

      for (auto & x : {fitter.Mother_Momentum(RefitTracks).perp(),fitter.Mother_Momentum(RefitTracks).eta()} ) tempBpt.push_back(x);
      tempBpt.push_back(fitter.Mother_Momentum(RefitTracks).phi());
      nt.NRbks_pt_eta_phi.push_back(tempBpt);

      for (auto & x : {fitter.Mother_XYZ().x(),fitter.Mother_XYZ().y(),fitter.Mother_XYZ().z()} ) tempBx.push_back(x);
      nt.NRbks_x_y_z.push_back(tempBx);
     
      for (auto & x : {fitter.Mother_XYZError().cxx(),fitter.Mother_XYZError().cyy(),fitter.Mother_XYZError().czz()} ) tempBex.push_back(x);
      nt.NRbks_ex_ey_ez.push_back(tempBex); 
      
      for (auto & x : {fitter.Mother_PtError(),fitter.Mother_EtaError(),fitter.Mother_PhiError()}) tempBept.push_back(x);
      nt.NRbks_ept_eeta_ephi.push_back(tempBept);

      if(IsE) nt.NRbks_mudecay.push_back(0);
      else  nt.NRbks_mudecay.push_back(1);

      nt.NRbks_lep1Id.push_back(imu1); nt.NRbks_lep2Id.push_back(imu2);
      nt.NRbks_mll.push_back((vmu1+vmu2).M()); 
      nt.NRbks_trkpair_index.push_back(indexpair[&pairtrk-&KstarTracks[0]]);
      nt.NRbks_k_sdxy.push_back(trkK.track().dxy(math::XYZPoint(nt.pvertex_x,nt.pvertex_y,nt.pvertex_z))/trkK.track().dxyError());      
      nt.NRbks_pi_sdxy.push_back(trkPi.track().dxy(math::XYZPoint(nt.pvertex_x,nt.pvertex_y,nt.pvertex_z))/trkPi.track().dxyError());
      for(unsigned int ichild=0; ichild<allParticles.size(); ichild++){
        temp.clear();      
        temp.push_back(fitter.Daughter_Momentum(ichild,RefitTracks).perp());
        temp.push_back(fitter.Daughter_Momentum(ichild,RefitTracks).eta());
        temp.push_back(fitter.Daughter_Momentum(ichild,RefitTracks).phi());
        temp.push_back(fitter.Daughter_Charge(ichild,RefitTracks)); 
        if (fitter.Daughter_Mass(ichild,RefitTracks)==0.493){
          nt.NRbks_Kpt_eta_phi.push_back(temp);
          vK.SetPtEtaPhiM(temp[0],temp[1],temp[2],0.493); }
        else if (fitter.Daughter_Mass(ichild,RefitTracks)==0.139) {
          nt.NRbks_Pipt_eta_phi.push_back(temp);
          vPi.SetPtEtaPhiM(temp[0],temp[1],temp[2],0.139); }
        else if (ichild==0 && fitter.Daughter_Mass(ichild,RefitTracks)==0.105)
          nt.NRbks_l1pt_eta_phi.push_back(temp); 
        else if (ichild==0 && fitter.Daughter_Mass(ichild,RefitTracks)==0.105) 
          nt.NRbks_l2pt_eta_phi.push_back(temp);    
      }
      nt.NRbks_ksmass.push_back((vK+vPi).M());
      }
     }
}




