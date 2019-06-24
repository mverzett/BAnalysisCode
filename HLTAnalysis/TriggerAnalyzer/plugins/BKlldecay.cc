#include "BKlldecay.h"

BKlldecay::BKlldecay(unsigned int &nmupairs_, std::vector<std::pair<unsigned int,unsigned int>>& used_lep_tracks, std::vector<reco::TransientTrack> & vmuTrack1,std::vector<reco::TransientTrack> & vmuTrack2, std::vector<reco::TransientTrack> & vKTrack, float beam_xz_,float beam_yz_ ,std::string RefitTracks_) : nmupairs(nmupairs_), used_muTracks_index(used_lep_tracks), muTrack1(vmuTrack1), muTrack2(vmuTrack2), KTrack(vKTrack), beam_xz(beam_xz_),beam_yz(beam_yz_),RefitBTracks(RefitTracks_)
  {
    temp.reserve(4); tempK.reserve(3); tempBpt.reserve(3); tempBx.reserve(3);
    tempBex.reserve(3); tempBept.reserve(3); 
    if (RefitBTracks=="mu"){
      RefitMu=true; RefitEl=false; 
    }else if (RefitBTracks=="both"){
      RefitMu=true; RefitEl=true; 
    }else {
      RefitMu=false; RefitEl=false;}
  
  }

BKlldecay::~BKlldecay(){};

void BKlldecay::Fill(NtupleContent &nt){
    
    for(const std::pair<unsigned int, unsigned int> &imu : used_muTracks_index){
      bool IsE(false); float m(0.105); 
      unsigned int mindex(&imu-&used_muTracks_index[0]);
      if ( mindex>nmupairs || mindex==nmupairs ){ 
        IsE=true; m=0.000511; part_mass=0.000511;
      }
      unsigned int imu1(imu.first); unsigned int imu2(imu.second);
      TLorentzVector vk,vmu1,vmu2;
      // std::cout<<imu1<<"  "<<imu2<<std::endl;
      if (!IsE){
         vmu1.SetPtEtaPhiM(nt.muon_pt[imu1],nt.muon_eta[imu1],nt.muon_phi[imu1],m);
         if (mindex<Mu2IsTrk && Mu2IsTrk>0)
          vmu2.SetPtEtaPhiM(nt.muon_pt[imu2],nt.muon_eta[imu2],nt.muon_phi[imu2],m); 
         else if ( (mindex==Mu2IsTrk || mindex>Mu2IsTrk) && Mu2IsTrk>0)
          vmu2.SetPtEtaPhiM(muTrack2[mindex].track().pt(),muTrack2[mindex].track().eta(),muTrack2[mindex].track().phi(),m); 
         else
          vmu2.SetPtEtaPhiM(nt.muon_pt[imu2],nt.muon_eta[imu2],nt.muon_phi[imu2],m);
         }
      else {
         vmu1.SetPtEtaPhiM(nt.el_pt[imu1],nt.el_eta[imu1],nt.el_phi[imu1],m);
         vmu2.SetPtEtaPhiM(nt.el_pt[imu2],nt.el_eta[imu2],nt.el_phi[imu2],m); }
     RefitTracks=false; 
     if (RefitMu && RefitEl) RefitTracks=true; 
     else if (RefitMu && !IsE) RefitTracks=true;
     for(const reco::TransientTrack & trk : KTrack){
       //  std::cout<<"K "<<trk.track().pt()<<" eta "<<trk.track().eta()<<" phi "<<trk.track().phi()<<std::endl;
      tempBpt.clear(); tempBx.clear(); tempBex.clear(); tempBept.clear();
      tempK.clear();
      vk.SetPtEtaPhiM(trk.track().pt(),trk.track().eta(),trk.track().phi(),0.495);
      if (deltaR(vmu1.Eta(),vmu1.Phi(),vk.Eta(),vk.Phi())<0.001) continue;
      if (deltaR(vmu2.Eta(),vmu2.Phi(),vk.Eta(),vk.Phi())<0.001) continue;
      if ((vmu1+vmu2+vk).M()<massmin_ || (vmu1+vmu2+vk).M()>massmax_) continue;
      if ((vmu1+vmu2+vk).Pt()<ptbmin_) continue;
      //std::cout<<"passed 1"<<std::endl;
        KinematicParticleFactoryFromTransientTrack pFactory;
        std::vector<RefCountedKinematicParticle> allParticles;
        allParticles.emplace_back(pFactory.particle(muTrack1[mindex],part_mass,chi,ndf,part_sigma));
        allParticles.emplace_back(pFactory.particle(muTrack2[mindex],part_mass,chi,ndf,part_sigma));
        allParticles.emplace_back(pFactory.particle(trk,kaon_mass,chi,ndf,kaon_sigma));
        TripleTrackKinFit fitter(allParticles);
      if (!fitter.success()) continue;
      float ChiProb(ChiSquaredProbability(fitter.chi(),fitter.dof()));
      if (ChiProb<probcut_) continue;
            
      GlobalPoint Dispbeamspot(-1*((nt.beam_x-fitter.Mother_XYZ().x())+(fitter.Mother_XYZ().z()-nt.beam_z)* beam_xz),-1*((nt.beam_y-fitter.Mother_XYZ().y())+ (fitter.Mother_XYZ().z()-nt.beam_z) * beam_yz), 0);              
      vperp.SetXYZ(Dispbeamspot.x(),Dispbeamspot.y(),0.);
      if (!RefitTracks)
         pperp.SetXYZ((vmu1+vmu2+vk).Px(),(vmu1+vmu2+vk).Py(),0);
      else
         pperp.SetXYZ(fitter.Mother_Momentum(RefitTracks).x(),fitter.Mother_Momentum(RefitTracks).y(),0);
      float cos(vperp.Dot(pperp)/(vperp.R()*pperp.R()));
      if (cos<coscut_) continue; 
      ///std::cout<<"passed 2"<<std::endl;  
      nt.NRb_cosTheta2D.push_back(cos);
      nt.NRb_bspot_lxy.push_back(Dispbeamspot.perp());
      nt.NRb_bspot_elxy.push_back(fitter.Mother_XYZError().rerr(Dispbeamspot));
      //  nt.NRb_trk_sdxy.push_back(trk.track().dxy(vertex_point)/trk.track().dxyError());      
      nt.NRb_trk_sdxy.emplace_back(trk.track().dxy(math::XYZPoint(nt.pvertex_x,nt.pvertex_y,nt.pvertex_z)));
      if (!RefitTracks)
        nt.NRb_mass.push_back((vmu1+vmu2+vk).M());
      else
        nt.NRb_mass.push_back(fitter.Mother_Mass(RefitTracks));
      nt.NRb_charge.push_back(fitter.Mother_Charge());

      for ( auto & x : {vk.Pt(),vk.Eta(),vk.Phi()}) tempK.push_back(x);
      nt.NRb_KUNFITpt_eta_phi.push_back(tempK);

      if (!RefitTracks){
	for (auto & x : {(vmu1+vmu2+vk).Pt(),(vmu1+vmu2+vk).Eta(),(vmu1+vmu2+vk).Phi()} ) tempBpt.push_back(x);      
      } else{
         for (auto & x : {fitter.Mother_Momentum(RefitTracks).perp(),fitter.Mother_Momentum(RefitTracks).eta()} ) tempBpt.push_back(x);
         tempBpt.push_back(fitter.Mother_Momentum(RefitTracks).phi());
      }
      nt.NRb_pt_eta_phi.push_back(tempBpt);

      for (auto & x : {fitter.Mother_XYZ().x(),fitter.Mother_XYZ().y(),fitter.Mother_XYZ().z()} ) tempBx.push_back(x);
      nt.NRb_x_y_z.push_back(tempBx);
     
      for (auto & x : {fitter.Mother_XYZError().cxx(),fitter.Mother_XYZError().cyy(),fitter.Mother_XYZError().czz()} ) tempBex.push_back(x);
      nt.NRb_ex_ey_ez.push_back(tempBex); 
      
      for (auto & x : {fitter.Mother_PtError(),fitter.Mother_EtaError(),fitter.Mother_PhiError()}) tempBept.push_back(x);
      nt.NRb_ept_eeta_ephi.push_back(tempBept);

      nt.NRb_chi_prob.push_back(ChiProb); 
      if(IsE) nt.NRb_mudecay.push_back(0);
      else  nt.NRb_mudecay.push_back(1);
      nt.NRb_lep1Id.push_back(imu1); nt.NRb_lep2Id.push_back(imu2);
      nt.NRb_mll.push_back((vmu1+vmu2).M()); nt.NRb_vtx_index.push_back(-1);
      for(unsigned int ichild=0; ichild<allParticles.size(); ichild++){
        temp.clear();      
        temp.push_back(fitter.Daughter_Momentum(ichild,true).perp());
        temp.push_back(fitter.Daughter_Momentum(ichild,true).eta());
        temp.push_back(fitter.Daughter_Momentum(ichild,true).phi());
        temp.push_back(fitter.Daughter_Charge(ichild,true)); 
        switch(ichild){
         case 0:{ nt.NRb_l1pt_eta_phi.push_back(temp); continue;}
         case 1:{ nt.NRb_l2pt_eta_phi.push_back(temp); continue;}   
         case 2:{ nt.NRb_Kpt_eta_phi.push_back(temp); continue;}
	}
      }
     }
    }
}

