#ifndef NTUPLECONTENT_H
#define NTUPLECONTENT_H
#include <vector>
#include "TTree.h"
#include <boost/range/combine.hpp>
#include <boost/tuple/tuple.hpp>

class NtupleContent{

public:
   NtupleContent();
   virtual ~NtupleContent();
   void SetTree(TTree *mytree);
   void SetNtupleVariables(TString Vars);
   void DefineVariables();
   void ClearVariables();
   //general
   int event,run_number,ls;
   float beam_x,beam_y,beam_z,beam_ex,beam_ey,beam_ez;
   float pvertex_x,pvertex_y,pvertex_z,pvertex_ex,pvertex_ey,pvertex_ez;
   std::vector<float> vertex_x,vertex_y,vertex_z,vertex_ex,vertex_ey,
   vertex_ez,vertex_chi,vertex_ndof;
   // trigger vars
   int trigger1,trigger2,trigger3,trigger4,trigger5,trigger6,trigger7,trigger8;
   int l1_seed1,l1_seed2,l1_seed3,l1_seed4,l1_seed5,l1_seed6;
   std::vector<std::vector<float>> tr1_obj_pt_eta_phi,tr2_obj_pt_eta_phi,
   tr3_obj_pt_eta_phi,tr4_obj_pt_eta_phi,tr5_obj_pt_eta_phi,tr6_obj_pt_eta_phi,
   tr7_obj_pt_eta_phi,tr8_obj_pt_eta_phi;
   float l1met,l1met_phi,l1ht,l1hf_met,l1hf_met_phi;
   std::vector<float> l1muon_eta,l1muon_iso,l1muon_pt,l1muon_phi,l1muon_qual,l1jet_pt,l1jet_phi,l1jet_eta,l1jet_iso,l1jet_qual;
   int nmuons,nel,ntracks,njets,muon_trgIndex;
   //muons
   std::vector<float> muon_pt,muon_eta,muon_phi,muon_qual,muon_charge,
   muon_global_flag,muon_tracker_flag,muon_standalone_flag,muon_dxy,muon_dz,
   muon_edxy,muon_edz,muon_d0,muon_ed0,muon_vx,muon_vz,muon_vy,muon_iso,
   muon_trkpt,muon_trketa,muon_trkphi;
   std::vector<bool> muon_medium,muon_loose,muon_tight,muon_soft;
   //electrons
   std::vector<float> el_pt,el_eta,el_phi,el_charge,el_vx,el_vy,el_vz,el_dxy,
   el_dz,el_mva_out,el_mva_iso,el_iso,el_trkpt,el_trketa,el_trkphi,el_edxy,el_edz,el_mva_map_value,el_mva_biased,el_mva_unbiased;
   std::vector<bool>el_veto,el_soft,el_medium,el_tight,el_islowpt;
   //lowpt as additional class
   std::vector<float> lowel_pt,lowel_eta,lowel_phi,lowel_charge,lowel_vx,lowel_vy,lowel_vz,lowel_mva_unbiased,lowel_trkpt,lowel_trketa,lowel_trkphi;
  //lowgsf as additional class
   std::vector<float> lowgsf_pt,lowgsf_eta,lowgsf_phi,lowgsf_charge,lowgsf_vx,lowgsf_vy,lowgsf_vz,lowgsf_mva_unbiased,lowgsf_modept,lowgsf_modeeta,lowgsf_modephi;
   //PF e as aditional class
   std::vector<float> pfel_pt,pfel_eta,pfel_phi,pfel_charge,pfel_vx,pfel_vy,pfel_vz;
   std::vector<bool> pfel_soft,pfel_medium,pfel_tight;
   //tracks
   std::vector<float> track_pt,track_eta,track_phi,track_norm_chi2,track_charge,
   track_dxy,track_dz,track_validhits,track_losthits,track_fromPV,track_highPurity;
   //met
   float ptmet,phimet;
   //jets
   std::vector<float> jet_pt,jet_eta,jet_phi,jet_cEmEF,jet_cHEF,jet_cHMult,
   jet_cMuEF,jet_cMult,jet_MuEF,jet_nEmEF,jet_nHEF,jet_nMult,jet_pEF,jet_eEF;
   //Reconstruction B->Kll
   std::vector<float> NRb_mass,NRb_chi_prob,NRb_charge,NRb_mll,
   NRb_ll_prob,NRb_trk_sdxy,NRb_trk_chi_norm,NRb_vtx_index,NRb_bspot_lxy
   ,NRb_bspot_elxy,NRb_llsvbsv,NRb_cosTheta2D,NRb_iso04,NRb_iso08,
   NRb_biso02,NRb_biso04,NRb_biso06,NRb_biso1p2;
   std::vector<std::vector<float>> NRb_Kpt_eta_phi,NRb_KUNFITpt_eta_phi,
   NRb_pt_eta_phi,NRb_l1pt_eta_phi,NRb_l2pt_eta_phi,
   NRb_x_y_z,NRb_ex_ey_ez,NRb_ept_eeta_ephi;
   std::vector<unsigned int> NRb_mudecay,NRb_lep1Id,NRb_lep2Id,NRb_trkId;
   //Reconstruction B->K*ll
   std::vector<float> NRbks_k_sdxy,NRbks_pi_sdxy,NRbks_mass,NRbks_charge,
   NRbks_chi_prob,NRbks_bspot_lxy, NRbks_bspot_elxy,
   NRbks_cosTheta2D,NRbks_mll,NRbks_ksmass; 
   std::vector<std::vector<float>> NRbks_pt_eta_phi,NRbks_x_y_z,
   NRbks_ex_ey_ez,NRbks_ept_eeta_ephi,NRbks_l2pt_eta_phi,
   NRbks_l1pt_eta_phi,NRbks_Kpt_eta_phi,NRbks_Pipt_eta_phi,
   NRbks_KUNFITpt_eta_phi,NRbks_PiUNFITpt_eta_phi; 
   std::vector<unsigned int> NRbks_mudecay,NRbks_lep1Id, NRbks_lep2Id;
   std::vector<int> NRbks_trkpair_index;
   //Reconstruction B->phill
   std::vector<float> NRbphi_k1_sdxy,NRbphi_k2_sdxy,NRbphi_mass,NRbphi_charge,
   NRbphi_chi_prob,NRbphi_bspot_lxy, NRbphi_bspot_elxy,
   NRbphi_cosTheta2D,NRbphi_mll,NRbphi_phimass; 
   std::vector<std::vector<float>> NRbphi_pt_eta_phi,NRbphi_x_y_z,
   NRbphi_ex_ey_ez,NRbphi_ept_eeta_ephi,NRbphi_l1pt_eta_phi,
   NRbphi_l2pt_eta_phi,NRbphi_K1pt_eta_phi,NRbphi_K2pt_eta_phi,
   NRbphi_K1UNFITpt_eta_phi,NRbphi_K2UNFITpt_eta_phi; 
   std::vector<unsigned int> NRbphi_mudecay,NRbphi_lep1Id, NRbphi_lep2Id;
   //gen
   float ngenB,ngenLep;
  std::vector<float> genB_pt,genB_phi,genB_eta,genB_pdgId,genB_Bindex,
  genB_daughter_pt,genB_daughter_eta,genB_daughter_phi,genB_daughter_pdgId,
  genB_daughter_Bindex,genB_daughter_Dindex,genB_granddaughter_pt,
  genB_granddaughter_eta,genB_granddaughter_phi,genB_granddaughter_pdgId,
  genB_granddaughter_Bindex,genB_granddaughter_Dindex;
  std::vector<float> genLep_pt,genLep_phi,genLep_eta,genLep_pdgId,genLep_mom; 
  //flat
  std::vector<float> NRb_pt,NRb_eta,NRb_phi,NRb_Kpt,NRb_Keta,NRb_Kphi,NRb_l1pt,
  NRb_l1eta,NRb_l1phi,NRb_l2pt,NRb_l2eta,NRb_l2phi,
  NRbks_pt,NRbks_eta,NRbks_phi,NRbks_Kpt,NRbks_Keta,NRbks_Kphi,NRbks_l1pt,
  NRbks_l1eta,NRbks_l1phi,NRbks_l2pt,NRbks_l2eta,NRbks_l2phi,NRbks_Pipt,
  NRbks_Pieta,NRbks_Piphi;
    
private:
   TTree * t1;

};
#endif
