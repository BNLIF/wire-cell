// cz: code modified from tutorials/tmva/TMVAClassification.C

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "tagger.h"

int main( int argc, char** argv )
{
  TaggerInfo tagger;
  
  tagger.pio_2_v_dis2 = new std::vector<float>;
  tagger.pio_2_v_angle2 = new std::vector<float>;
  tagger.pio_2_v_acc_length = new std::vector<float>;
  tagger.pio_2_v_flag = new std::vector<float>;
  tagger.sig_1_v_angle = new std::vector<float>;
  tagger.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger.sig_1_v_energy = new std::vector<float>;
  tagger.sig_1_v_energy_1 = new std::vector<float>;
  tagger.sig_1_v_flag = new std::vector<float>;
  tagger.sig_2_v_energy = new std::vector<float>;
  tagger.sig_2_v_shower_angle = new std::vector<float>;
  tagger.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_flag = new std::vector<float>;
  tagger.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_energy = new std::vector<float>;
  tagger.stw_2_v_angle = new std::vector<float>;
  tagger.stw_2_v_dir_length = new std::vector<float>;
  tagger.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_flag = new std::vector<float>;
  tagger.stw_3_v_angle = new std::vector<float>;
  tagger.stw_3_v_dir_length = new std::vector<float>;
  tagger.stw_3_v_energy = new std::vector<float>;
  tagger.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_3_v_flag = new std::vector<float>;
  tagger.stw_4_v_angle = new std::vector<float>;
  tagger.stw_4_v_dis = new std::vector<float>;
  tagger.stw_4_v_energy = new std::vector<float>;
  tagger.stw_4_v_flag = new std::vector<float>;
  tagger.br3_3_v_energy = new std::vector<float>;
  tagger.br3_3_v_angle = new std::vector<float>;
  tagger.br3_3_v_dir_length = new std::vector<float>;
  tagger.br3_3_v_length = new std::vector<float>;
  tagger.br3_3_v_flag = new std::vector<float>;
  tagger.br3_5_v_dir_length = new std::vector<float>;
  tagger.br3_5_v_total_length = new std::vector<float>;
  tagger.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger.br3_5_v_n_seg = new std::vector<float>;
  tagger.br3_5_v_angle = new std::vector<float>;
  tagger.br3_5_v_sg_length = new std::vector<float>;
  tagger.br3_5_v_energy = new std::vector<float>;
  tagger.br3_5_v_n_main_segs = new std::vector<float>;
  tagger.br3_5_v_n_segs = new std::vector<float>;
  tagger.br3_5_v_shower_main_length = new std::vector<float>;
  tagger.br3_5_v_shower_total_length = new std::vector<float>;
  tagger.br3_5_v_flag = new std::vector<float>;
  tagger.br3_6_v_angle = new std::vector<float>;
  tagger.br3_6_v_angle1 = new std::vector<float>;
  tagger.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger.br3_6_v_direct_length = new std::vector<float>;
  tagger.br3_6_v_length = new std::vector<float>;
  tagger.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger.br3_6_v_energy = new std::vector<float>;
  tagger.br3_6_v_flag = new std::vector<float>;
  tagger.tro_1_v_particle_type = new std::vector<float>;
  tagger.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger.tro_1_v_min_dis = new std::vector<float>;
  tagger.tro_1_v_sg1_length = new std::vector<float>;
  tagger.tro_1_v_shower_main_length = new std::vector<float>;
  tagger.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger.tro_1_v_tmp_length = new std::vector<float>;
  tagger.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger.tro_1_v_flag = new std::vector<float>;
  tagger.tro_2_v_energy = new std::vector<float>;
  tagger.tro_2_v_stem_length = new std::vector<float>;
  tagger.tro_2_v_iso_angle = new std::vector<float>;
  tagger.tro_2_v_max_length = new std::vector<float>;
  tagger.tro_2_v_angle = new std::vector<float>;
  tagger.tro_2_v_flag = new std::vector<float>;
  tagger.tro_4_v_dir2_mag = new std::vector<float>;
  tagger.tro_4_v_angle = new std::vector<float>;
  tagger.tro_4_v_angle1 = new std::vector<float>;
  tagger.tro_4_v_angle2 = new std::vector<float>;
  tagger.tro_4_v_length = new std::vector<float>;
  tagger.tro_4_v_length1 = new std::vector<float>;
  tagger.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_energy = new std::vector<float>;
  tagger.tro_4_v_shower_main_length = new std::vector<float>;
  tagger.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger.tro_4_v_flag = new std::vector<float>;
  tagger.tro_5_v_max_angle = new std::vector<float>;
  tagger.tro_5_v_min_angle = new std::vector<float>;
  tagger.tro_5_v_max_length = new std::vector<float>;
  tagger.tro_5_v_iso_angle = new std::vector<float>;
  tagger.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger.tro_5_v_min_count = new std::vector<float>;
  tagger.tro_5_v_max_count = new std::vector<float>;
  tagger.tro_5_v_energy = new std::vector<float>;
  tagger.tro_5_v_flag = new std::vector<float>;
  tagger.lol_1_v_energy = new std::vector<float>;
  tagger.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_1_v_nseg = new std::vector<float>;
  tagger.lol_1_v_angle = new std::vector<float>;
  tagger.lol_1_v_flag = new std::vector<float>;
  tagger.lol_2_v_length = new std::vector<float>;
  tagger.lol_2_v_angle = new std::vector<float>;
  tagger.lol_2_v_type = new std::vector<float>;
  tagger.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_2_v_energy = new std::vector<float>;
  tagger.lol_2_v_shower_main_length = new std::vector<float>;
  tagger.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger.lol_2_v_flag = new std::vector<float>;
  tagger.cosmict_flag_10 = new std::vector<float>;
  tagger.cosmict_10_flag_inside = new std::vector<float>;
  tagger.cosmict_10_vtx_z = new std::vector<float>;
  tagger.cosmict_10_flag_shower = new std::vector<float>;
  tagger.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger.cosmict_10_angle_beam = new std::vector<float>;
  tagger.cosmict_10_length = new std::vector<float>;
  tagger.numu_cc_flag_1 = new std::vector<float>;
  tagger.numu_cc_1_particle_type = new std::vector<float>;
  tagger.numu_cc_1_length = new std::vector<float>;
  tagger.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger.numu_cc_1_direct_length = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger.numu_cc_flag_2 = new std::vector<float>;
  tagger.numu_cc_2_length = new std::vector<float>;
  tagger.numu_cc_2_total_length = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_all = new std::vector<float>;
  tagger.pio_2_v_dis2 = new std::vector<float>;
  tagger.pio_2_v_angle2 = new std::vector<float>;
  tagger.pio_2_v_acc_length = new std::vector<float>;
  tagger.pio_2_v_flag = new std::vector<float>;
  tagger.sig_1_v_angle = new std::vector<float>;
  tagger.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger.sig_1_v_energy = new std::vector<float>;
  tagger.sig_1_v_energy_1 = new std::vector<float>;
  tagger.sig_1_v_flag = new std::vector<float>;
  tagger.sig_2_v_energy = new std::vector<float>;
  tagger.sig_2_v_shower_angle = new std::vector<float>;
  tagger.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_flag = new std::vector<float>;
  tagger.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_energy = new std::vector<float>;
  tagger.stw_2_v_angle = new std::vector<float>;
  tagger.stw_2_v_dir_length = new std::vector<float>;
  tagger.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_flag = new std::vector<float>;
  tagger.stw_3_v_angle = new std::vector<float>;
  tagger.stw_3_v_dir_length = new std::vector<float>;
  tagger.stw_3_v_energy = new std::vector<float>;
  tagger.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_3_v_flag = new std::vector<float>;
  tagger.stw_4_v_angle = new std::vector<float>;
  tagger.stw_4_v_dis = new std::vector<float>;
  tagger.stw_4_v_energy = new std::vector<float>;
  tagger.stw_4_v_flag = new std::vector<float>;
  tagger.br3_3_v_energy = new std::vector<float>;
  tagger.br3_3_v_angle = new std::vector<float>;
  tagger.br3_3_v_dir_length = new std::vector<float>;
  tagger.br3_3_v_length = new std::vector<float>;
  tagger.br3_3_v_flag = new std::vector<float>;
  tagger.br3_5_v_dir_length = new std::vector<float>;
  tagger.br3_5_v_total_length = new std::vector<float>;
  tagger.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger.br3_5_v_n_seg = new std::vector<float>;
  tagger.br3_5_v_angle = new std::vector<float>;
  tagger.br3_5_v_sg_length = new std::vector<float>;
  tagger.br3_5_v_energy = new std::vector<float>;
  tagger.br3_5_v_n_main_segs = new std::vector<float>;
  tagger.br3_5_v_n_segs = new std::vector<float>;
  tagger.br3_5_v_shower_main_length = new std::vector<float>;
  tagger.br3_5_v_shower_total_length = new std::vector<float>;
  tagger.br3_5_v_flag = new std::vector<float>;
  tagger.br3_6_v_angle = new std::vector<float>;
  tagger.br3_6_v_angle1 = new std::vector<float>;
  tagger.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger.br3_6_v_direct_length = new std::vector<float>;
  tagger.br3_6_v_length = new std::vector<float>;
  tagger.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger.br3_6_v_energy = new std::vector<float>;
  tagger.br3_6_v_flag = new std::vector<float>;
  tagger.tro_1_v_particle_type = new std::vector<float>;
  tagger.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger.tro_1_v_min_dis = new std::vector<float>;
  tagger.tro_1_v_sg1_length = new std::vector<float>;
  tagger.tro_1_v_shower_main_length = new std::vector<float>;
  tagger.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger.tro_1_v_tmp_length = new std::vector<float>;
  tagger.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger.tro_1_v_flag = new std::vector<float>;
  tagger.tro_2_v_energy = new std::vector<float>;
  tagger.tro_2_v_stem_length = new std::vector<float>;
  tagger.tro_2_v_iso_angle = new std::vector<float>;
  tagger.tro_2_v_max_length = new std::vector<float>;
  tagger.tro_2_v_angle = new std::vector<float>;
  tagger.tro_2_v_flag = new std::vector<float>;
  tagger.tro_4_v_dir2_mag = new std::vector<float>;
  tagger.tro_4_v_angle = new std::vector<float>;
  tagger.tro_4_v_angle1 = new std::vector<float>;
  tagger.tro_4_v_angle2 = new std::vector<float>;
  tagger.tro_4_v_length = new std::vector<float>;
  tagger.tro_4_v_length1 = new std::vector<float>;
  tagger.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_energy = new std::vector<float>;
  tagger.tro_4_v_shower_main_length = new std::vector<float>;
  tagger.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger.tro_4_v_flag = new std::vector<float>;
  tagger.tro_5_v_max_angle = new std::vector<float>;
  tagger.tro_5_v_min_angle = new std::vector<float>;
  tagger.tro_5_v_max_length = new std::vector<float>;
  tagger.tro_5_v_iso_angle = new std::vector<float>;
  tagger.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger.tro_5_v_min_count = new std::vector<float>;
  tagger.tro_5_v_max_count = new std::vector<float>;
  tagger.tro_5_v_energy = new std::vector<float>;
  tagger.tro_5_v_flag = new std::vector<float>;
  tagger.lol_1_v_energy = new std::vector<float>;
  tagger.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_1_v_nseg = new std::vector<float>;
  tagger.lol_1_v_angle = new std::vector<float>;
  tagger.lol_1_v_flag = new std::vector<float>;
  tagger.lol_2_v_length = new std::vector<float>;
  tagger.lol_2_v_angle = new std::vector<float>;
  tagger.lol_2_v_type = new std::vector<float>;
  tagger.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_2_v_energy = new std::vector<float>;
  tagger.lol_2_v_shower_main_length = new std::vector<float>;
  tagger.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger.lol_2_v_flag = new std::vector<float>;
  tagger.cosmict_flag_10 = new std::vector<float>;
  tagger.cosmict_10_flag_inside = new std::vector<float>;
  tagger.cosmict_10_vtx_z = new std::vector<float>;
  tagger.cosmict_10_flag_shower = new std::vector<float>;
  tagger.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger.cosmict_10_angle_beam = new std::vector<float>;
  tagger.cosmict_10_length = new std::vector<float>;
  tagger.numu_cc_flag_1 = new std::vector<float>;
  tagger.numu_cc_1_particle_type = new std::vector<float>;
  tagger.numu_cc_1_length = new std::vector<float>;
  tagger.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger.numu_cc_1_direct_length = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger.numu_cc_flag_2 = new std::vector<float>;
  tagger.numu_cc_2_length = new std::vector<float>;
  tagger.numu_cc_2_total_length = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_all = new std::vector<float>;
  
  // intrinsic nueCC
  TFile *file1 = new TFile("run1_intrinsic_nue_POT1.2E23.root");
  TTree *t1 = (TTree*)file1->Get("bdt");
  set_tree_address(t1, tagger);
  double pot_1 = 1.2e23/2.;

  TFile *file2 = new TFile("run3_intrinsic_nue_POT8.3E22.root");
  TTree *t2 = (TTree*)file2->Get("bdt");
  set_tree_address(t2, tagger);
  double pot_2 = 8.3e22/2.;

  // intrinsic nueCC lowE patch
  TFile *file3 = new TFile("run1_intrinsic_nue_LowE_POT6.1E23.root");
  TTree *t3 = (TTree*)file3->Get("bdt");
  set_tree_address(t3, tagger);
  double pot_3 = 6.1e23/2.;

  TFile *file4 = new TFile("run3_intrinsic_nue_LowE_POT6.0E23.root");
  TTree *t4 = (TTree*)file4->Get("bdt");
  set_tree_address(t4, tagger);
  double pot_4 = 6.0e23/2.;

  // nu overlay
  TFile *file5 = new TFile("run1_bnb_nu_POT1.2E21.root");
  TTree *t5 = (TTree*)file5->Get("bdt");
  set_tree_address(t5, tagger);
  double pot_5 = 1.2e21/2.;

  TFile *file6 = new TFile("run3_bnb_nu_POT1.2E21.root");
  TTree *t6 = (TTree*)file6->Get("bdt");
  set_tree_address(t6, tagger);
  double pot_6 = 1.21e21/2.;
  
  // EXTBNB
  TFile *file7 = new TFile("run1_ext_bnb_C1_gt10_wcp_v00_14_00_POT1.2E20.root");
  TTree *t7 = (TTree*)file7->Get("bdt");
  set_tree_address(t7, tagger);
  double pot_7 = 1.2e20/2.;

  TFile *file8 = new TFile("run3_ext_bnb_F_G1_POT1.9E20.root");
  TTree *t8 = (TTree*)file8->Get("bdt");
  set_tree_address(t8, tagger);
  double pot_8 = 1.92e20/2.;


  // nu overlay low E patch ...
  TFile *file9 = new TFile("run1_bnb_nu_LowE_POT1.6E21.root");
  TTree *t9 = (TTree*)file9->Get("bdt");
  set_tree_address(t9, tagger);
  double pot_9 = 1.6e21/2.;

  TFile *file10 = new TFile("run3_bnb_nu_LowE_POT1.5E21.root");
  TTree *t10 = (TTree*)file10->Get("bdt");
  set_tree_address(t10, tagger);
  double pot_10 = 1.5e21/2.;

  
  TFile *new_file = new TFile("bdt.root","RECREATE");
  TTree* Tsig = new TTree("sig","sig");
  TTree* Tbkg = new TTree("bkg","bkg");
  Tsig->SetDirectory(new_file);
  Tbkg->SetDirectory(new_file);
  put_tree_address(Tsig, tagger);
  put_tree_address(Tbkg, tagger);

  // signal, intrinsic nue ...
  for (Int_t i=0;i!=t1->GetEntries();i++){
    t1->GetEntry(i);

    // weight corrupted ...
    if (std::isnan(tagger.weight_cv) || std::isnan(tagger.weight_spline) || std::isinf(tagger.weight_cv) || std::isinf(tagger.weight_spline)) continue;
    // no EM showers ...
    if (tagger.br_filled!=1) continue;
    // odd sub-run number ...
    if (tagger.subrun %2 == 1) continue;

    
    
    tagger.weight = tagger.weight_spline * tagger.weight_cv * (1 + tagger.weight_lee);
    tagger.lowEweight = (1. + 5 * tagger.weight_lee)/(1. + tagger.weight_lee);

    // scale intrinsic nue higher ...
    if (tagger.truth_nuEnergy > 400) tagger.weight *= (pot_1+pot_2+pot_3+pot_4)/(pot_1+pot_2);
    
    if (tagger.nuvtx_diff < 1 && tagger.showervtx_diff < 1 && tagger.truth_isCC == 1 &&
	abs(tagger.truth_nuPdg)==12 && tagger.truth_vtxInside == 1){
      Tsig->Fill();
    }
  }
  
  for (Int_t i=0;i!=t2->GetEntries();i++){
    t2->GetEntry(i);

    // weight corrupted ...
    if (std::isnan(tagger.weight_cv) || std::isnan(tagger.weight_spline) || std::isinf(tagger.weight_cv) || std::isinf(tagger.weight_spline)) continue;
    // no EM showers ...
    if (tagger.br_filled!=1) continue;
    // odd sub-run number ...
    if (tagger.subrun %2 == 1) continue;

    
    
    tagger.weight = tagger.weight_spline * tagger.weight_cv * (1 + tagger.weight_lee);
    tagger.lowEweight = (1. + 5 * tagger.weight_lee)/(1. + tagger.weight_lee);

    // scale intrinsic nue higher ...
    if (tagger.truth_nuEnergy > 400) tagger.weight *= (pot_1+pot_2+pot_3+pot_4)/(pot_1+pot_2);
    
    if (tagger.nuvtx_diff < 1 && tagger.showervtx_diff < 1 && tagger.truth_isCC == 1 &&
	abs(tagger.truth_nuPdg)==12 && tagger.truth_vtxInside == 1){
      Tsig->Fill();
    }
  }

  // low energy patch ...
  for (Int_t i=0;i!=t3->GetEntries();i++){
    t3->GetEntry(i);

    // weight corrupted ...
    if (std::isnan(tagger.weight_cv) || std::isnan(tagger.weight_spline) || std::isinf(tagger.weight_cv) || std::isinf(tagger.weight_spline)) continue;
    // no EM showers ...
    if (tagger.br_filled!=1) continue;
    // odd sub-run number ...
    if (tagger.subrun %2 == 1) continue;

    
    
    tagger.weight = tagger.weight_spline * tagger.weight_cv * (1 + tagger.weight_lee);
    tagger.lowEweight = (1. + 5 * tagger.weight_lee)/(1. + tagger.weight_lee);

            
    if (tagger.nuvtx_diff < 1 && tagger.showervtx_diff < 1 && tagger.truth_isCC == 1 &&
	abs(tagger.truth_nuPdg)==12 && tagger.truth_vtxInside == 1){
      Tsig->Fill();
    }
  }

  for (Int_t i=0;i!=t4->GetEntries();i++){
    t4->GetEntry(i);

    // weight corrupted ...
    if (std::isnan(tagger.weight_cv) || std::isnan(tagger.weight_spline) || std::isinf(tagger.weight_cv) || std::isinf(tagger.weight_spline)) continue;
    // no EM showers ...
    if (tagger.br_filled!=1) continue;
    // odd sub-run number ...
    if (tagger.subrun %2 == 1) continue;

    
    
    tagger.weight = tagger.weight_spline * tagger.weight_cv * (1 + tagger.weight_lee);
    tagger.lowEweight = (1. + 5 * tagger.weight_lee)/(1. + tagger.weight_lee);
   
    
    if (tagger.nuvtx_diff < 1 && tagger.showervtx_diff < 1 && tagger.truth_isCC == 1 &&
	abs(tagger.truth_nuPdg)==12 && tagger.truth_vtxInside == 1){
      Tsig->Fill();
    }
  }
  
  
  
  

  // background  numu overlay ...
  for (Int_t i=0;i!=t5->GetEntries();i++){
    t5->GetEntry(i);

    if (std::isnan(tagger.weight_cv) || std::isnan(tagger.weight_spline) || std::isinf(tagger.weight_cv) || std::isinf(tagger.weight_spline)) continue;
    // require EM shower ...
    if (tagger.br_filled!=1) continue;
    // only use even subrun number to train ...
    if (tagger.subrun %2 == 1) continue;
    
    // also need to exclude the NC nu-electron elastic scattering 
    if (tagger.truth_isCC==1 && abs(tagger.truth_nuPdg)==12 || tagger.truth_nuIntType==1098 ) continue;

    if (tagger.truth_nuEnergy > 400)
      tagger.weight = tagger.weight_spline * tagger.weight_cv * (pot_1+pot_2+pot_3+pot_4)/(pot_5+pot_6 );
    else
      tagger.weight = tagger.weight_spline * tagger.weight_cv * (pot_1+pot_2+pot_3+pot_4)/(pot_5+pot_6 +pot_9 + pot_10);
    
    tagger.lowEweight = 1;
    
    Tbkg->Fill();
  }

  for (Int_t i=0;i!=t6->GetEntries();i++){
    t6->GetEntry(i);

    if (std::isnan(tagger.weight_cv) || std::isnan(tagger.weight_spline) || std::isinf(tagger.weight_cv) || std::isinf(tagger.weight_spline)) continue;
    // require EM shower ...
    if (tagger.br_filled!=1) continue;
    // only use even subrun number to train ...
    if (tagger.subrun %2 == 1) continue;
    
    // also need to exclude the NC nu-electron elastic scattering 
    if (tagger.truth_isCC==1 && abs(tagger.truth_nuPdg)==12 || tagger.truth_nuIntType==1098 ) continue;

    if (tagger.truth_nuEnergy > 400)
      tagger.weight = tagger.weight_spline * tagger.weight_cv * (pot_1+pot_2+pot_3+pot_4)/(pot_5+pot_6);
    else
      tagger.weight = tagger.weight_spline * tagger.weight_cv * (pot_1+pot_2+pot_3+pot_4)/(pot_5+pot_6+pot_9+pot_10);
    
    tagger.lowEweight = 1;
    
    Tbkg->Fill();
  }

  // low energy patck nu overlay ...
  for (Int_t i=0;i!=t9->GetEntries();i++){
    t9->GetEntry(i);

    if (std::isnan(tagger.weight_cv) || std::isnan(tagger.weight_spline) || std::isinf(tagger.weight_cv) || std::isinf(tagger.weight_spline)) continue;
    // require EM shower ...
    if (tagger.br_filled!=1) continue;
    // only use even subrun number to train ...
    if (tagger.subrun %2 == 1) continue;
    
    // also need to exclude the NC nu-electron elastic scattering 
    if (tagger.truth_isCC==1 && abs(tagger.truth_nuPdg)==12 || tagger.truth_nuIntType==1098 ) continue;

    tagger.weight = tagger.weight_spline * tagger.weight_cv * (pot_1+pot_2+pot_3+pot_4)/(pot_5+pot_6 +pot_9 + pot_10);
    
    tagger.lowEweight = 1;
    
    Tbkg->Fill();
  }
  

  for (Int_t i=0;i!=t10->GetEntries();i++){
    t10->GetEntry(i);

    if (std::isnan(tagger.weight_cv) || std::isnan(tagger.weight_spline) || std::isinf(tagger.weight_cv) || std::isinf(tagger.weight_spline)) continue;
    // require EM shower ...
    if (tagger.br_filled!=1) continue;
    // only use even subrun number to train ...
    if (tagger.subrun %2 == 1) continue;
    
    // also need to exclude the NC nu-electron elastic scattering 
    if (tagger.truth_isCC==1 && abs(tagger.truth_nuPdg)==12 || tagger.truth_nuIntType==1098 ) continue;

    tagger.weight = tagger.weight_spline * tagger.weight_cv * (pot_1+pot_2+pot_3+pot_4)/(pot_5+pot_6 +pot_9 + pot_10);
    
    tagger.lowEweight = 1;
    
    Tbkg->Fill();
  }
  
  


  // EXT background 
  for (Int_t i=0;i!=t7->GetEntries();i++){
    t7->GetEntry(i);

    //    if (std::isnan(tagger.weight_cv) || std::isnan(tagger.weight_spline) || std::isinf(tagger.weight_cv) || std::isinf(tagger.weight_spline)) continue;
    // require EM shower ...
    if (tagger.br_filled!=1) continue;
    // only use even subrun number to train ...
    if (tagger.subrun %2 == 1) continue;
    
    // also need to exclude the NC nu-electron elastic scattering 
    //if (tagger.truth_isCC==1 && abs(tagger.truth_nuPdg)==12 || tagger.truth_nuIntType==1098 ) continue;
    
    tagger.weight = (pot_1+pot_2+pot_3+pot_4)/(pot_7+pot_8);
    tagger.lowEweight = 1;
    
    Tbkg->Fill();
  }

  for (Int_t i=0;i!=t8->GetEntries();i++){
    t8->GetEntry(i);

    //    if (std::isnan(tagger.weight_cv) || std::isnan(tagger.weight_spline) || std::isinf(tagger.weight_cv) || std::isinf(tagger.weight_spline)) continue;
    // require EM shower ...
    if (tagger.br_filled!=1) continue;
    // only use even subrun number to train ...
    if (tagger.subrun %2 == 1) continue;
    
    // also need to exclude the NC nu-electron elastic scattering 
    //if (tagger.truth_isCC==1 && abs(tagger.truth_nuPdg)==12 || tagger.truth_nuIntType==1098 ) continue;
    
    tagger.weight = (pot_1+pot_2+pot_3+pot_4)/(pot_7+pot_8);
    tagger.lowEweight = 1;
    
    Tbkg->Fill();
  }
  
  
  new_file->Write();
  new_file->Close();
  
  return 0;
}


