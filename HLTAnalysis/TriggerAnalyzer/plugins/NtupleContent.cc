#include "NtupleContent.h"

NtupleContent::NtupleContent(){}
NtupleContent::~NtupleContent(){}

void NtupleContent::SetTree(TTree *mytree){
  t1=mytree; }

void NtupleContent::ClearVariables(){
 //trg
 trigger1=0,trigger2=0,trigger3=0,trigger4=0,trigger5=0,trigger6=0,trigger7=0,trigger8=0;
 l1_seed1=0,l1_seed2=0,l1_seed3=0,l1_seed4=0,l1_seed5=0,l1_seed6=0;
 l1met=0,l1met_phi=-99,l1ht=0,l1hf_met=0,l1hf_met_phi=-99;
 l1muon_eta.clear(); l1muon_pt.clear(); l1muon_phi.clear(); l1muon_qual.clear();
 l1jet_pt.clear(); l1jet_phi.clear(); l1jet_eta.clear(); l1jet_iso.clear();
 l1jet_qual.clear();
 tr1_obj_pt_eta_phi.clear(); tr2_obj_pt_eta_phi.clear(); 
 tr3_obj_pt_eta_phi.clear(); tr4_obj_pt_eta_phi.clear(); 
 tr5_obj_pt_eta_phi.clear(); tr6_obj_pt_eta_phi.clear();
 tr7_obj_pt_eta_phi.clear(); tr8_obj_pt_eta_phi.clear();
 //general
  event=0,run_number=0,ls=0;
  beam_x=-999,beam_y=-999,beam_z=-999,beam_ex=-999,beam_ey=-999,beam_ez=-999;
  pvertex_x=-999,pvertex_y=-999,pvertex_z=-999,pvertex_ex=-999,pvertex_ey=-999,
  pvertex_ez=-999;
  vertex_x.clear(); vertex_y.clear(); vertex_z.clear(); vertex_ex.clear(); 
  vertex_ey.clear(); vertex_ez.clear(); vertex_chi.clear(); vertex_ndof.clear();
  nmuons=0,nel=0,njets=0,ntracks=0;
  //muons
  muon_pt.clear(); muon_eta.clear(); muon_phi.clear(); muon_qual.clear();
  muon_charge.clear();  muon_global_flag.clear(); muon_tracker_flag.clear(); 
  muon_standalone_flag.clear(); muon_dxy.clear(); muon_dz.clear();
  muon_edxy.clear(); muon_edz.clear(); muon_d0.clear(); muon_ed0.clear(); 
  muon_vx.clear(); muon_vz.clear(); muon_vy.clear(); muon_iso.clear();
  muon_trkpt.clear(); muon_trketa.clear(); muon_trkphi.clear(); 
  muon_medium.clear(); muon_loose.clear(); muon_tight.clear(); muon_soft.clear();
  muon_trgIndex=0;
  //electrons
  el_pt.clear(); el_eta.clear(); el_phi.clear(); el_charge.clear();
  el_vx.clear(); el_vy.clear(); el_vz.clear(); el_dxy.clear(); el_dz.clear();
  el_mva_out.clear(); el_mva_iso.clear(); el_iso.clear(); el_trkpt.clear(); 
  el_trketa.clear(); el_trkphi.clear(); el_edxy.clear(); el_edz.clear(); 
  el_mva_map_value.clear(); el_veto.clear(); el_soft.clear(); el_medium.clear();
  el_tight.clear(); el_mva_biased.clear(); el_mva_unbiased.clear();
  el_islowpt.clear();
  //lowpt as additional class
  lowel_pt.clear(); lowel_eta.clear(); lowel_phi.clear(); lowel_charge.clear();
  lowel_vx.clear(); lowel_vy.clear(); lowel_vz.clear(); lowel_mva_unbiased.clear(); 
  lowel_trkpt.clear(); lowel_trketa.clear(); lowel_trkphi.clear();
  //lowgsf as additional class
   lowgsf_pt.clear(); lowgsf_eta.clear(); lowgsf_phi.clear(); lowgsf_charge.clear();
  lowgsf_vx.clear(); lowgsf_vy.clear(); lowgsf_vz.clear(); lowgsf_mva_unbiased.clear();
  lowgsf_modept.clear(); lowgsf_modeeta.clear(); lowgsf_modephi.clear();
  //PFel as additional class
  pfel_pt.clear(); pfel_eta.clear(); pfel_phi.clear(); pfel_charge.clear();
  pfel_vx.clear(); pfel_vy.clear(); pfel_vz.clear(); pfel_soft.clear();
  pfel_medium.clear(); pfel_tight.clear();
  //tracks
  track_pt.clear(); track_eta.clear(); track_phi.clear(); track_charge.clear();
  track_norm_chi2.clear(); track_dxy.clear(); track_dz.clear(); 
  track_validhits.clear(); track_losthits.clear(); track_fromPV.clear(); 
  track_highPurity.clear();
  //met
  ptmet=0,phimet=-999;
  //jets
  jet_pt.clear(); jet_eta.clear(); jet_phi.clear(); jet_cEmEF.clear();
  jet_cHEF.clear(); jet_cHMult.clear(); jet_cMuEF.clear(); jet_cMult.clear();
  jet_MuEF.clear(); jet_nEmEF.clear(); jet_nHEF.clear(); jet_nMult.clear();
  jet_pEF.clear(); jet_eEF.clear();
  //B->Kll
   NRb_mass.clear(); NRb_chi_prob.clear(); NRb_charge.clear();
   NRb_mll.clear(); NRb_ll_prob.clear(); NRb_trk_sdxy.clear(); 
   NRb_trk_chi_norm.clear(); NRb_vtx_index.clear(); 
   NRb_bspot_lxy.clear(); NRb_bspot_elxy.clear();
   NRb_llsvbsv.clear(); NRb_cosTheta2D.clear(); NRb_iso04.clear(); 
   NRb_iso08.clear(); NRb_biso02.clear(); NRb_biso04.clear(); 
   NRb_biso06.clear();  NRb_biso1p2.clear();
   NRb_Kpt_eta_phi.clear(); NRb_KUNFITpt_eta_phi.clear();
   NRb_pt_eta_phi.clear(); NRb_l1pt_eta_phi.clear(); 
   NRb_l2pt_eta_phi.clear(); NRb_x_y_z.clear(); NRb_ex_ey_ez.clear(); 
   NRb_ept_eeta_ephi.clear(); NRb_trkId.clear();
   NRb_mudecay.clear(); NRb_lep1Id.clear(); NRb_lep2Id.clear();
   //B->K*ll
   NRbks_k_sdxy.clear(); NRbks_pi_sdxy.clear(); NRbks_mass.clear();
   NRbks_charge.clear(); NRbks_chi_prob.clear(); NRbks_bspot_lxy.clear(); 
   NRbks_bspot_elxy.clear(); NRbks_cosTheta2D.clear(); 
   NRbks_mll.clear(); NRbks_ksmass.clear(); 
   NRbks_pt_eta_phi.clear(); NRbks_x_y_z.clear();
   NRbks_ex_ey_ez.clear(); NRbks_ept_eeta_ephi.clear(); 
   NRbks_l2pt_eta_phi.clear(); NRbks_l1pt_eta_phi.clear(); 
   NRbks_Kpt_eta_phi.clear(); NRbks_Pipt_eta_phi.clear(); 
   NRbks_KUNFITpt_eta_phi.clear(); NRbks_PiUNFITpt_eta_phi.clear(); 
   NRbks_mudecay.clear(); NRbks_lep1Id.clear(); NRbks_lep2Id.clear();
   NRbks_trkpair_index.clear();
   //B->phill
   NRbphi_k1_sdxy.clear(); NRbphi_k2_sdxy.clear(); NRbphi_mass.clear();
   NRbphi_charge.clear(); NRbphi_chi_prob.clear(); NRbphi_bspot_lxy.clear(); 
   NRbphi_bspot_elxy.clear(); NRbphi_cosTheta2D.clear(); 
   NRbphi_mll.clear(); NRbphi_phimass.clear(); NRbphi_pt_eta_phi.clear();
   NRbphi_x_y_z.clear(); NRbphi_ex_ey_ez.clear(); NRbphi_ept_eeta_ephi.clear(); 
   NRbphi_l2pt_eta_phi.clear(); NRbphi_l1pt_eta_phi.clear(); 
   NRbphi_K1pt_eta_phi.clear(); NRbphi_K2pt_eta_phi.clear(); 
   NRbphi_K1UNFITpt_eta_phi.clear(); NRbphi_K2UNFITpt_eta_phi.clear(); 
   NRbphi_mudecay.clear(); NRbphi_lep1Id.clear(); NRbphi_lep2Id.clear();
   //gen
  ngenB=0,ngenLep=0;
  genB_pt.clear(); genB_phi.clear(); genB_eta.clear(); genB_pdgId.clear();
  genB_Bindex.clear(); genB_daughter_pt.clear(); genB_daughter_eta.clear();
  genB_daughter_phi.clear(); genB_daughter_pdgId.clear();
  genB_daughter_Bindex.clear(); genB_daughter_Dindex.clear(); 
  genB_granddaughter_pt.clear(); genB_granddaughter_eta.clear(); 
  genB_granddaughter_phi.clear(); genB_granddaughter_pdgId.clear();
  genB_granddaughter_Bindex.clear(); genB_granddaughter_Dindex.clear();
  genLep_pt.clear(); genLep_phi.clear(); genLep_eta.clear();
  genLep_pdgId.clear(); genLep_mom.clear(); 

}

void NtupleContent::SetNtupleVariables(TString Vars){
  bool writeGen=false,writeTrk=false,writeKll=false,writeKstarll=false,writeLowPtEl=false,writeLowPtGsf=false,writePFel=false,Lite=false,Flat=false,writePhill=false; 
  if (Vars.Contains("ALL")){
    writeGen=true,writeTrk=true,writeKll=true,writeKstarll=true,writeLowPtEl=true,writeLowPtGsf=true,writePFel=true;}
  if (Vars.Contains("GEN"))
      writeGen=true;
  if (Vars.Contains("TRK"))
      writeTrk=true;
  if (Vars.Contains("KLL"))
      writeKll=true;
  if (Vars.Contains("KstarLL"))
      writeKstarll=true;
  if(Vars.Contains("LowPtEl"))
     writeLowPtEl=true;
  if(Vars.Contains("LowPtGsf"))
     writeLowPtGsf=true;
  if(Vars.Contains("PFel"))
     writePFel=true;
  if (Vars.Contains("Lite"))
     Lite=true;
  if (Vars.Contains("Flat"))
     Flat=true;

  t1->Branch("event",&event); t1->Branch("run_number",&run_number);
  t1->Branch("ls",&ls);
  t1->Branch("beam_x",&beam_x); t1->Branch("beam_y",&beam_y);
  t1->Branch("beam_z",&beam_z); 
  
  if (Lite || Flat){
   t1->Branch("vertex_x",&vertex_x); t1->Branch("vertex_y",&vertex_y);
   t1->Branch("vertex_z",&vertex_z);
   t1->Branch("HLT_path1",&trigger1); t1->Branch("HLT_path2",&trigger2);
   t1->Branch("HLT_path3",&trigger3); t1->Branch("HLT_path4",&trigger4);
   t1->Branch("HLT_path5",&trigger5); t1->Branch("HLT_path6",&trigger6);
   t1->Branch("HLT_path7",&trigger7); t1->Branch("HLT_path8",&trigger8);
   t1->Branch("L1_seed1",&l1_seed1); t1->Branch("L1_seed2",&l1_seed2);
   t1->Branch("L1_seed3",&l1_seed3); t1->Branch("L1_seed4",&l1_seed4);
   t1->Branch("L1_seed5",&l1_seed5); t1->Branch("L1_seed6",&l1_seed6);
   t1->Branch("nmuon",&nmuons); t1->Branch("muon_pt",&muon_pt);
   t1->Branch("muon_eta",&muon_eta); t1->Branch("muon_phi",&muon_phi);
   t1->Branch("muon_charge",&muon_charge); t1->Branch("muon_dxy",&muon_dxy);
   t1->Branch("muon_edxy",&muon_edxy); t1->Branch("muon_soft",&muon_soft);
   t1->Branch("muon_trgIndex",&muon_trgIndex); t1->Branch("nelectron",&nel);
   t1->Branch("el_pt",&el_pt); t1->Branch("el_eta",&el_eta);
   t1->Branch("el_phi",&el_phi); t1->Branch("el_charge",&el_charge);
   t1->Branch("el_mva_unbiased",&el_mva_unbiased);
   t1->Branch("el_islowpt",&el_islowpt);
   if (writeKll){
     if (Flat){
       t1->Branch("NRb_pt",&NRb_pt); t1->Branch("NRb_eta",&NRb_eta);
       t1->Branch("NRb_phi",&NRb_phi); t1->Branch("NRb_Kpt",&NRb_Kpt);
       t1->Branch("NRb_Keta",&NRb_Keta); t1->Branch("NRb_Kphi",&NRb_Kphi);
       t1->Branch("NRb_l1pt",&NRb_l1pt); t1->Branch("NRb_eta",&NRb_l1eta);
       t1->Branch("NRb_l1phi",&NRb_l1phi); t1->Branch("NRb_l2pt",&NRb_l2pt); 
       t1->Branch("NRb_l2eta",&NRb_l2eta);  t1->Branch("NRb_l2phi",&NRb_l2phi);
     } else {
       t1->Branch("NRb_pt_eta_phi",&NRb_pt_eta_phi); 
       t1->Branch("NRb_Kpt_eta_phi_charge",&NRb_Kpt_eta_phi);
       t1->Branch("NRb_l1pt_eta_phi_charge",&NRb_l1pt_eta_phi);
       t1->Branch("NRb_l2pt_eta_phi_charge",&NRb_l2pt_eta_phi);
     }

    t1->Branch("NRb_mass",&NRb_mass);
    t1->Branch("NRb_chi_prob",&NRb_chi_prob); t1->Branch("NRb_mudecay",&NRb_mudecay);
     t1->Branch("NRb_charge",&NRb_charge);
     t1->Branch("NRb_lep1Id",&NRb_lep1Id); t1->Branch("NRb_lep2Id",&NRb_lep2Id);
     t1->Branch("NRb_bspot_lxy",&NRb_bspot_lxy); t1->Branch("NRb_bspot_elxy",&NRb_bspot_elxy);
     t1->Branch("NRb_cosTheta2D",&NRb_cosTheta2D); t1->Branch("NRb_mll",&NRb_mll);}
   if (writeKstarll){
    if (Flat){
     t1->Branch("NRbks_pt",&NRbks_pt); t1->Branch("NRbks_eta",&NRbks_eta);
     t1->Branch("NRbks_phi",&NRbks_phi); t1->Branch("NRbks_Kpt",&NRbks_Kpt);
     t1->Branch("NRbks_Keta",&NRbks_Keta); t1->Branch("NRbks_Kphi",&NRbks_Kphi);
     t1->Branch("NRbks_l1pt",&NRbks_l1pt); t1->Branch("NRbks_eta",&NRbks_l1eta);
     t1->Branch("NRbks_l1phi",&NRbks_l1phi); t1->Branch("NRbks_l2pt",&NRbks_l2pt);
     t1->Branch("NRbks_l2eta",&NRbks_l2eta); t1->Branch("NRbks_l2phi",&NRbks_l2phi);
     t1->Branch("NRbks_Pipt",&NRbks_Pipt);
     t1->Branch("NRbks_Pieta",&NRbks_Pieta); t1->Branch("NRbks_Piphi",&NRbks_Piphi);
    } else {
       t1->Branch("NRbks_pt_eta_phi",&NRbks_pt_eta_phi); 
       t1->Branch("NRbks_Kpt_eta_phi_charge",&NRbks_Kpt_eta_phi);
       t1->Branch("NRbks_l1pt_eta_phi_charge",&NRbks_l1pt_eta_phi);
       t1->Branch("NRbks_l2pt_eta_phi_charge",&NRbks_l2pt_eta_phi);
       t1->Branch("NRbks_Pipt_eta_phi",&NRbks_Pipt_eta_phi);
     }
    t1->Branch("NRbks_mass",&NRbks_mass); t1->Branch("NRbks_charge",&NRbks_charge);
    t1->Branch("NRbks_chi_prob",&NRbks_chi_prob); t1->Branch("NRbks_bspot_lxy",&NRbks_bspot_lxy);
    t1->Branch("NRbks_bspot_elxy",&NRbks_bspot_elxy); t1->Branch("NRbks_cosTheta2D",&NRbks_cosTheta2D);
    t1->Branch("NRbks_mll",&NRbks_mll); t1->Branch("NRbks_ksmass",&NRbks_ksmass);
    t1->Branch("NRbks_trkpair_index",&NRbks_trkpair_index);
   }
   return;
  }


  t1->Branch("beam_ex",&beam_ex);
  t1->Branch("beam_ey",&beam_ey); t1->Branch("beam_ez",&beam_ez);
  t1->Branch("pvertex_x",&pvertex_x); t1->Branch("pvertex_y",&pvertex_y);
  t1->Branch("pvertex_z",&pvertex_z); t1->Branch("pvertex_ex",&pvertex_ex);
  t1->Branch("pvertex_ey",&pvertex_ey); t1->Branch("pvertex_ez",&pvertex_ez);
  t1->Branch("vertex_x",&vertex_x); t1->Branch("vertex_y",&vertex_y);
  t1->Branch("vertex_z",&vertex_z); t1->Branch("vertex_ex",&vertex_ex);
  t1->Branch("vertex_ey",&vertex_ey); t1->Branch("vertex_ez",&vertex_ez);
  t1->Branch("vertex_chi",&vertex_chi); t1->Branch("vertex_ndof",&vertex_ndof);
  //triggers
  t1->Branch("HLT_path1",&trigger1); t1->Branch("HLT_path2",&trigger2);
  t1->Branch("HLT_path3",&trigger3); t1->Branch("HLT_path4",&trigger4);
  t1->Branch("HLT_path5",&trigger5); t1->Branch("HLT_path6",&trigger6);
  t1->Branch("HLT_path7",&trigger7); t1->Branch("HLT_path8",&trigger8);
  t1->Branch("L1_seed1",&l1_seed1); t1->Branch("L1_seed2",&l1_seed2);
  t1->Branch("L1_seed3",&l1_seed3); t1->Branch("L1_seed4",&l1_seed4);
  t1->Branch("L1_seed5",&l1_seed5); t1->Branch("L1_seed6",&l1_seed6);
  t1->Branch("HLT1Obj_pt_eta_phi_charge",&tr1_obj_pt_eta_phi);
  t1->Branch("HLT2Obj_pt_eta_phi_charge",&tr2_obj_pt_eta_phi);
  t1->Branch("HLT3Obj_pt_eta_phi_charge",&tr3_obj_pt_eta_phi);
  t1->Branch("HLT4Obj_pt_eta_phi_charge",&tr4_obj_pt_eta_phi);
  t1->Branch("HLT5Obj_pt_eta_phi_charge",&tr5_obj_pt_eta_phi);
  t1->Branch("HLT6Obj_pt_eta_phi_charge",&tr6_obj_pt_eta_phi);
  t1->Branch("HLT7Obj_pt_eta_phi_charge",&tr7_obj_pt_eta_phi);
  t1->Branch("HLT8Obj_pt_eta_phi_charge",&tr8_obj_pt_eta_phi);
  t1->Branch("l1muon_pt",&l1muon_pt); t1->Branch("l1muon_eta",&l1muon_eta);
  t1->Branch("l1muon_phi",&l1muon_phi); t1->Branch("l1muon_qual",&l1muon_qual);
  t1->Branch("l1jet_pt",&l1jet_pt); t1->Branch("l1jet_eta",&l1jet_eta);
  t1->Branch("l1jet_phi",&l1jet_phi); t1->Branch("l1jet_iso",&l1jet_iso);
  t1->Branch("l1jet_qual",&l1jet_qual); t1->Branch("l1met",&l1met);
  t1->Branch("l1met_phi",&l1met_phi); t1->Branch("l1ht",&l1ht);
  t1->Branch("l1hf_met",&l1hf_met); t1->Branch("l1hf_met_phi",&l1hf_met_phi);
  
  //gen
  if (writeGen){
   t1->Branch("ngenB",&ngenB); t1->Branch("genB_pt",&genB_pt);
   t1->Branch("genB_phi",&genB_phi); t1->Branch("genB_eta",&genB_eta);
   t1->Branch("genB_pdgId",&genB_pdgId); t1->Branch("genB_Bindex",&genB_Bindex);
   t1->Branch("genB_daughter_pt",&genB_daughter_pt);
   t1->Branch("genB_daughter_eta",&genB_daughter_eta);
   t1->Branch("genB_daughter_phi",&genB_daughter_phi);
   t1->Branch("genB_daughter_pdgId",&genB_daughter_pdgId);
   t1->Branch("genB_daughter_Bindex",&genB_daughter_Bindex);
   t1->Branch("genB_daughter_Dindex",&genB_daughter_Dindex);
   t1->Branch("genB_granddaughter_pt",&genB_granddaughter_pt);
   t1->Branch("genB_granddaughter_eta",&genB_granddaughter_eta);
   t1->Branch("genB_granddaughter_phi",&genB_granddaughter_phi);
   t1->Branch("genB_granddaughter_pdgId",&genB_granddaughter_pdgId);
   t1->Branch("genB_granddaughter_Bindex",&genB_granddaughter_Bindex);
   t1->Branch("genB_granddaughter_Dindex",&genB_granddaughter_Dindex);
   t1->Branch("ngenLep",&ngenLep); t1->Branch("genLep_pt",&genLep_pt);
   t1->Branch("genLep_phi",&genLep_phi); t1->Branch("genLep_eta",&genLep_eta);
   t1->Branch("genLep_pdgId",&genLep_pdgId); t1->Branch("genLep_mom",&genLep_mom);
   }
  //muons
  t1->Branch("nmuon",&nmuons); t1->Branch("muon_pt",&muon_pt);
  t1->Branch("muon_eta",&muon_eta); t1->Branch("muon_phi",&muon_phi);
  t1->Branch("muon_charge",&muon_charge); t1->Branch("muon_qual",&muon_qual);
  t1->Branch("muon_dxy",&muon_dxy); t1->Branch("muon_dz",&muon_dz);
  t1->Branch("muon_edxy",&muon_edxy); t1->Branch("muon_edz",&muon_edz);
  t1->Branch("muon_d0",&muon_d0); t1->Branch("muon_ed0",&muon_ed0);
  t1->Branch("muon_vx",&muon_vx); t1->Branch("muon_vy",&muon_vy);
  t1->Branch("muon_vz",&muon_vz); t1->Branch("muon_iso",&muon_iso);
  t1->Branch("muon_global_flag",&muon_global_flag);
  t1->Branch("muon_tracker_flag",&muon_tracker_flag);
  t1->Branch("muon_standalone_flag",&muon_standalone_flag);
  t1->Branch("muon_soft",&muon_soft); t1->Branch("muon_loose",&muon_loose);
  t1->Branch("muon_medium",&muon_medium);  t1->Branch("muon_tight",&muon_tight);
  t1->Branch("muon_trkpt",&muon_trkpt); t1->Branch("muon_trketa",&muon_trketa);
  t1->Branch("muon_trkphi",&muon_trkphi); t1->Branch("muon_trgIndex",&muon_trgIndex);
  //electrons
  t1->Branch("nelectron",&nel); t1->Branch("el_pt",&el_pt);
  t1->Branch("el_eta",&el_eta); t1->Branch("el_phi",&el_phi);
  t1->Branch("el_charge",&el_charge); t1->Branch("el_dxy",&el_dxy);
  t1->Branch("el_dz",&el_dz); t1->Branch("el_edxy",&el_edxy);
  t1->Branch("el_edz",&el_edz); t1->Branch("el_vx",&el_vx);
  t1->Branch("el_vy",&el_vy); t1->Branch("el_vz",&el_vz);
  t1->Branch("el_mva_out",&el_mva_out); t1->Branch("el_mva_iso",&el_mva_iso);
  t1->Branch("el_iso",&el_iso); t1->Branch("el_veto",&el_veto);
  t1->Branch("el_soft",&el_soft); t1->Branch("el_medium",&el_medium);
  t1->Branch("el_tight",&el_tight); t1->Branch("el_trkpt",&el_trkpt);
  t1->Branch("el_trketa",&el_trketa); t1->Branch("el_trkphi",&el_trkphi);
  t1->Branch("el_mva_map_value",&el_mva_map_value); 
  t1->Branch("el_mva_biased",&el_mva_biased); 
  t1->Branch("el_mva_unbiased",&el_mva_unbiased); 
  t1->Branch("el_islowpt",&el_islowpt);
  if (writeLowPtEl){
   t1->Branch("lowel_pt",&lowel_pt); t1->Branch("lowel_eta",&lowel_eta);
   t1->Branch("lowel_phi",&lowel_phi); t1->Branch("lowel_charge",&lowel_charge);
   t1->Branch("lowel_vx",&lowel_vx); t1->Branch("lowel_vy",&lowel_vy);
   t1->Branch("lowel_vz",&lowel_vz); t1->Branch("lowel_mva_unbiased",&lowel_mva_unbiased);
   t1->Branch("lowel_trkpt",&lowel_trkpt); t1->Branch("lowel_trketa",&lowel_trketa);
   t1->Branch("lowel_trkphi",&lowel_trkphi);
  }
  if (writeLowPtGsf){
   t1->Branch("lowgsf_pt",&lowgsf_pt); t1->Branch("lowgsf_eta",&lowgsf_eta);
   t1->Branch("lowgsf_phi",&lowgsf_phi); t1->Branch("lowgsf_charge",&lowgsf_charge);
   t1->Branch("lowgsf_vx",&lowgsf_vx); t1->Branch("lowgsf_vy",&lowgsf_vy);
   t1->Branch("lowgsf_vz",&lowgsf_vz); t1->Branch("lowgsf_mva_unbiased",&lowgsf_mva_unbiased);
   t1->Branch("lowgsf_modept",&lowgsf_modept); t1->Branch("lowgsf_modeeta",&lowgsf_modeeta);
   t1->Branch("lowgsf_modephi",&lowgsf_modephi);
  }
 if (writePFel){
  t1->Branch("pfel_pt",&pfel_pt); t1->Branch("pfel_eta",&pfel_eta);
  t1->Branch("pfel_phi",&pfel_phi);
  t1->Branch("pfel_pfcharge",&pfel_charge);
  t1->Branch("pfel_vx",&pfel_vx); t1->Branch("pfel_vy",&pfel_vy);
  t1->Branch("pfel_vz",&pfel_vz); t1->Branch("pfel_soft",&pfel_soft);
  t1->Branch("pfel_medium",&el_medium); t1->Branch("pfel_tight",&pfel_tight);
 }
  //tracks
  if (writeTrk){
    t1->Branch("ntracks",&ntracks); t1->Branch("track_pt",&track_pt);
    t1->Branch("track_eta",&track_eta); t1->Branch("track_phi",&track_phi);
    t1->Branch("track_norm_chi2",&track_norm_chi2); t1->Branch("track_charge",&track_charge);
    t1->Branch("track_dxy",&track_dxy); t1->Branch("track_dz",&track_dz);
    t1->Branch("track_validhits",&track_validhits); t1->Branch("track_losthits",&track_losthits);
    t1->Branch("track_fromPV",&track_fromPV); t1->Branch("track_highPurity",&track_highPurity);
  }
  //jets
  t1->Branch("jet_pt",&jet_pt); t1->Branch("jet_eta",&jet_eta);
  t1->Branch("jet_phi",&jet_phi); t1->Branch("jet_cEmEF",&jet_cEmEF);
  t1->Branch("jet_cHEF",&jet_cHEF); t1->Branch("jet_cHMult",&jet_cHMult);
  t1->Branch("jet_cMuEF",&jet_cMuEF); t1->Branch("jet_cMult",&jet_cMult);
  t1->Branch("jet_MuEF",&jet_MuEF); t1->Branch("jet_nEmEF",&jet_nEmEF);
  t1->Branch("jet_nHEF",&jet_nHEF); t1->Branch("jet_nMult",&jet_nMult);
  t1->Branch("jet_eEF",&jet_eEF); t1->Branch("jet_pEF",&jet_pEF);
  //met
  t1->Branch("met_pt",&ptmet); t1->Branch("met_phi",&phimet);
   //B->Kll
 if (writeKll){
   t1->Branch("NRb_pt_eta_phi",&NRb_pt_eta_phi); t1->Branch("NRb_mass",&NRb_mass);
   t1->Branch("NRb_chi_prob",&NRb_chi_prob); t1->Branch("NRb_x_y_z",&NRb_x_y_z);
   t1->Branch("NRb_ex_ey_ez",&NRb_ex_ey_ez);
   t1->Branch("NRb_ept_eeta_ephi",&NRb_ept_eeta_ephi);
   t1->Branch("NRb_mudecay",&NRb_mudecay);
   t1->Branch("NRb_Kpt_eta_phi_charge",&NRb_Kpt_eta_phi);
   t1->Branch("NRb_KUNFITpt_eta_phi",&NRb_KUNFITpt_eta_phi);
   t1->Branch("NRb_charge",&NRb_charge);
   t1->Branch("NRb_l1pt_eta_phi_charge",&NRb_l1pt_eta_phi);
   t1->Branch("NRb_l2pt_eta_phi_charge",&NRb_l2pt_eta_phi);
   t1->Branch("NRb_lep1Id",&NRb_lep1Id); t1->Branch("NRb_lep2Id",&NRb_lep2Id);
   t1->Branch("NRb_mll",&NRb_mll); t1->Branch("NRb_ll_prob",&NRb_ll_prob);
   t1->Branch("NRb_trk_sdxy",&NRb_trk_sdxy); t1->Branch("NRb_vtx_index",&NRb_vtx_index);
   t1->Branch("NRb_trk_chi_norm",&NRb_trk_chi_norm);
   if (writeTrk) t1->Branch("NRb_trkId",&NRb_trkId);
   t1->Branch("NRb_bspot_lxy",&NRb_bspot_lxy);
   t1->Branch("NRb_bspot_elxy",&NRb_bspot_elxy);
   t1->Branch("NRb_llsvbsv",&NRb_llsvbsv);
   t1->Branch("NRb_cosTheta2D",&NRb_cosTheta2D);
   t1->Branch("NRb_iso04",&NRb_iso04); t1->Branch("NRb_iso08",&NRb_iso08);
   t1->Branch("NRb_biso02",&NRb_biso02); t1->Branch("NRb_biso04",&NRb_biso04);
   t1->Branch("NRb_biso06",&NRb_biso06); t1->Branch("NRb_biso1p2",&NRb_biso1p2);
} 
//B->K*ll
if (writeKstarll){
  
  t1->Branch("NRbks_mass",&NRbks_mass); t1->Branch("NRbks_charge",&NRbks_charge);
  t1->Branch("NRbks_chi_prob",&NRbks_chi_prob);
  t1->Branch("NRbks_bspot_lxy",&NRbks_bspot_lxy);
  t1->Branch("NRbks_bspot_elxy",&NRbks_bspot_elxy);
  t1->Branch("NRbks_cosTheta2D",&NRbks_cosTheta2D);
  t1->Branch("NRbks_pt_eta_phi",&NRbks_pt_eta_phi);
  t1->Branch("NRbks_x_y_z",&NRbks_x_y_z);
  t1->Branch("NRbks_ept_eeta_ephi",&NRbks_ept_eeta_ephi);
  t1->Branch("NRbks_ex_ey_ez",&NRbks_ex_ey_ez);
  t1->Branch("NRbks_mudecay",&NRbks_mudecay);
  t1->Branch("NRbks_lep1Id",&NRbks_lep1Id); t1->Branch("NRbks_lep2Id",&NRbks_lep2Id);
  t1->Branch("NRbks_Kpt_eta_phi",&NRbks_Kpt_eta_phi);
  t1->Branch("NRbks_Pipt_eta_phi",&NRbks_Pipt_eta_phi);
  t1->Branch("NRbks_l1pt_eta_phi",&NRbks_l1pt_eta_phi);
  t1->Branch("NRbks_l2pt_eta_phi",&NRbks_l2pt_eta_phi);
  t1->Branch("NRbks_mll",&NRbks_mll); t1->Branch("NRbks_ksmass",&NRbks_ksmass);
  t1->Branch("NRbks_KUNFITpt_eta_phi",&NRbks_KUNFITpt_eta_phi);
  t1->Branch("NRbks_PiUNFITpt_eta_phi",&NRbks_PiUNFITpt_eta_phi);
  t1->Branch("NRbks_trkpair_index",&NRbks_trkpair_index);
  t1->Branch("NRbks_k_sdxy",&NRbks_k_sdxy); t1->Branch("NRbks_pi_sdxy",&NRbks_pi_sdxy);
}
//B->phill
if (writePhill){
  
  t1->Branch("NRbphi_mass",&NRbphi_mass); t1->Branch("NRbphi_charge",&NRbphi_charge);
  t1->Branch("NRbphi_chi_prob",&NRbphi_chi_prob);
  t1->Branch("NRbphi_bspot_lxy",&NRbphi_bspot_lxy);
  t1->Branch("NRbphi_bspot_elxy",&NRbphi_bspot_elxy);
  t1->Branch("NRbphi_cosTheta2D",&NRbphi_cosTheta2D);
  t1->Branch("NRbphi_pt_eta_phi",&NRbphi_pt_eta_phi);
  t1->Branch("NRbphi_x_y_z",&NRbphi_x_y_z);
  t1->Branch("NRbphi_ept_eeta_ephi",&NRbphi_ept_eeta_ephi);
  t1->Branch("NRbphi_ex_ey_ez",&NRbphi_ex_ey_ez);
  t1->Branch("NRbphi_mudecay",&NRbphi_mudecay);
  t1->Branch("NRbphi_lep1Id",&NRbphi_lep1Id); t1->Branch("NRbphi_lep2Id",&NRbphi_lep2Id);
  t1->Branch("NRbphi_K1pt_eta_phi",&NRbphi_K1pt_eta_phi);
  t1->Branch("NRbphi_K2pt_eta_phi",&NRbphi_K2pt_eta_phi);
  t1->Branch("NRbphi_l1pt_eta_phi",&NRbphi_l1pt_eta_phi);
  t1->Branch("NRbphi_l2pt_eta_phi",&NRbphi_l2pt_eta_phi);
  t1->Branch("NRbphi_mll",&NRbphi_mll); t1->Branch("NRbphi_phimass",&NRbphi_phimass);
  t1->Branch("NRbphi_K1UNFITpt_eta_phi",&NRbphi_K1UNFITpt_eta_phi);
  t1->Branch("NRbphi_K2UNFITpt_eta_phi",&NRbphi_K2UNFITpt_eta_phi);
  t1->Branch("NRbphi_k1_sdxy",&NRbphi_k1_sdxy); t1->Branch("NRbphi_k2_sdxy",&NRbphi_k2_sdxy);
}

 
}
