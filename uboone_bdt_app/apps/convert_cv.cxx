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
#include "TMath.h"

#include "WCPLEEANA/tagger.h"

#include "WCPLEEANA/eval.h"

using namespace std;
using namespace LEEana;

//#include "WCPLEEANA/bdt.h"
#include "WCPLEEANA/pot.h"
#include "WCPLEEANA/pfeval.h"
#include "WCPLEEANA/kine.h"

int main( int argc, char** argv )
{
  if (argc < 3) {
    std::cout << "merge_det #input_file_cv #output_file " << std::endl;
    return -1;
  }
  TString input_file_cv = argv[1];
  TString out_file = argv[2];

  float fail_percentage = 0.2;
   for (Int_t i=1;i!=argc;i++){
    switch(argv[i][1]){
    case 't':
      fail_percentage = atof(&argv[i][2]);
      break;
    }
   }
   bool flag_data = true;

  TFile *file1 = new TFile(input_file_cv);
  TTree *T_BDTvars_cv = (TTree*)file1->Get("wcpselection/T_BDTvars");
  TTree *T_eval_cv = (TTree*)file1->Get("wcpselection/T_eval");
  TTree *T_pot_cv = (TTree*)file1->Get("wcpselection/T_pot");
  TTree *T_PFeval_cv = (TTree*)file1->Get("wcpselection/T_PFeval");
  TTree *T_KINEvars_cv = (TTree*)file1->Get("wcpselection/T_KINEvars");

  if (T_eval_cv->GetBranch("weight_cv")) flag_data =false;
  
  

  EvalInfo eval_cv;
  eval_cv.file_type = new std::string();
  POTInfo pot_cv;
  TaggerInfo tagger_cv;
  PFevalInfo pfeval_cv;
  KineInfo kine_cv;
  kine_cv.kine_energy_particle = new std::vector<float>;
  kine_cv.kine_energy_info = new std::vector<int>;
  kine_cv.kine_particle_type = new std::vector<int>;
  kine_cv.kine_energy_included = new std::vector<int>;
  
  
  tagger_cv.pio_2_v_dis2 = new std::vector<float>;
  tagger_cv.pio_2_v_angle2 = new std::vector<float>;
  tagger_cv.pio_2_v_acc_length = new std::vector<float>;
  tagger_cv.pio_2_v_flag = new std::vector<float>;
  tagger_cv.sig_1_v_angle = new std::vector<float>;
  tagger_cv.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger_cv.sig_1_v_energy = new std::vector<float>;
  tagger_cv.sig_1_v_energy_1 = new std::vector<float>;
  tagger_cv.sig_1_v_flag = new std::vector<float>;
  tagger_cv.sig_2_v_energy = new std::vector<float>;
  tagger_cv.sig_2_v_shower_angle = new std::vector<float>;
  tagger_cv.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger_cv.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger_cv.sig_2_v_flag = new std::vector<float>;
  tagger_cv.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.stw_2_v_energy = new std::vector<float>;
  tagger_cv.stw_2_v_angle = new std::vector<float>;
  tagger_cv.stw_2_v_dir_length = new std::vector<float>;
  tagger_cv.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger_cv.stw_2_v_flag = new std::vector<float>;
  tagger_cv.stw_3_v_angle = new std::vector<float>;
  tagger_cv.stw_3_v_dir_length = new std::vector<float>;
  tagger_cv.stw_3_v_energy = new std::vector<float>;
  tagger_cv.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.stw_3_v_flag = new std::vector<float>;
  tagger_cv.stw_4_v_angle = new std::vector<float>;
  tagger_cv.stw_4_v_dis = new std::vector<float>;
  tagger_cv.stw_4_v_energy = new std::vector<float>;
  tagger_cv.stw_4_v_flag = new std::vector<float>;
  tagger_cv.br3_3_v_energy = new std::vector<float>;
  tagger_cv.br3_3_v_angle = new std::vector<float>;
  tagger_cv.br3_3_v_dir_length = new std::vector<float>;
  tagger_cv.br3_3_v_length = new std::vector<float>;
  tagger_cv.br3_3_v_flag = new std::vector<float>;
  tagger_cv.br3_5_v_dir_length = new std::vector<float>;
  tagger_cv.br3_5_v_total_length = new std::vector<float>;
  tagger_cv.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger_cv.br3_5_v_n_seg = new std::vector<float>;
  tagger_cv.br3_5_v_angle = new std::vector<float>;
  tagger_cv.br3_5_v_sg_length = new std::vector<float>;
  tagger_cv.br3_5_v_energy = new std::vector<float>;
  tagger_cv.br3_5_v_n_main_segs = new std::vector<float>;
  tagger_cv.br3_5_v_n_segs = new std::vector<float>;
  tagger_cv.br3_5_v_shower_main_length = new std::vector<float>;
  tagger_cv.br3_5_v_shower_total_length = new std::vector<float>;
  tagger_cv.br3_5_v_flag = new std::vector<float>;
  tagger_cv.br3_6_v_angle = new std::vector<float>;
  tagger_cv.br3_6_v_angle1 = new std::vector<float>;
  tagger_cv.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger_cv.br3_6_v_direct_length = new std::vector<float>;
  tagger_cv.br3_6_v_length = new std::vector<float>;
  tagger_cv.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger_cv.br3_6_v_energy = new std::vector<float>;
  tagger_cv.br3_6_v_flag = new std::vector<float>;
  tagger_cv.tro_1_v_particle_type = new std::vector<float>;
  tagger_cv.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger_cv.tro_1_v_min_dis = new std::vector<float>;
  tagger_cv.tro_1_v_sg1_length = new std::vector<float>;
  tagger_cv.tro_1_v_shower_main_length = new std::vector<float>;
  tagger_cv.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger_cv.tro_1_v_tmp_length = new std::vector<float>;
  tagger_cv.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger_cv.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger_cv.tro_1_v_flag = new std::vector<float>;
  tagger_cv.tro_2_v_energy = new std::vector<float>;
  tagger_cv.tro_2_v_stem_length = new std::vector<float>;
  tagger_cv.tro_2_v_iso_angle = new std::vector<float>;
  tagger_cv.tro_2_v_max_length = new std::vector<float>;
  tagger_cv.tro_2_v_angle = new std::vector<float>;
  tagger_cv.tro_2_v_flag = new std::vector<float>;
  tagger_cv.tro_4_v_dir2_mag = new std::vector<float>;
  tagger_cv.tro_4_v_angle = new std::vector<float>;
  tagger_cv.tro_4_v_angle1 = new std::vector<float>;
  tagger_cv.tro_4_v_angle2 = new std::vector<float>;
  tagger_cv.tro_4_v_length = new std::vector<float>;
  tagger_cv.tro_4_v_length1 = new std::vector<float>;
  tagger_cv.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger_cv.tro_4_v_energy = new std::vector<float>;
  tagger_cv.tro_4_v_shower_main_length = new std::vector<float>;
  tagger_cv.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger_cv.tro_4_v_flag = new std::vector<float>;
  tagger_cv.tro_5_v_max_angle = new std::vector<float>;
  tagger_cv.tro_5_v_min_angle = new std::vector<float>;
  tagger_cv.tro_5_v_max_length = new std::vector<float>;
  tagger_cv.tro_5_v_iso_angle = new std::vector<float>;
  tagger_cv.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger_cv.tro_5_v_min_count = new std::vector<float>;
  tagger_cv.tro_5_v_max_count = new std::vector<float>;
  tagger_cv.tro_5_v_energy = new std::vector<float>;
  tagger_cv.tro_5_v_flag = new std::vector<float>;
  tagger_cv.lol_1_v_energy = new std::vector<float>;
  tagger_cv.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger_cv.lol_1_v_nseg = new std::vector<float>;
  tagger_cv.lol_1_v_angle = new std::vector<float>;
  tagger_cv.lol_1_v_flag = new std::vector<float>;
  tagger_cv.lol_2_v_length = new std::vector<float>;
  tagger_cv.lol_2_v_angle = new std::vector<float>;
  tagger_cv.lol_2_v_type = new std::vector<float>;
  tagger_cv.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger_cv.lol_2_v_energy = new std::vector<float>;
  tagger_cv.lol_2_v_shower_main_length = new std::vector<float>;
  tagger_cv.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger_cv.lol_2_v_flag = new std::vector<float>;
  tagger_cv.cosmict_flag_10 = new std::vector<float>;
  tagger_cv.cosmict_10_flag_inside = new std::vector<float>;
  tagger_cv.cosmict_10_vtx_z = new std::vector<float>;
  tagger_cv.cosmict_10_flag_shower = new std::vector<float>;
  tagger_cv.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger_cv.cosmict_10_angle_beam = new std::vector<float>;
  tagger_cv.cosmict_10_length = new std::vector<float>;
  tagger_cv.numu_cc_flag_1 = new std::vector<float>;
  tagger_cv.numu_cc_1_particle_type = new std::vector<float>;
  tagger_cv.numu_cc_1_length = new std::vector<float>;
  tagger_cv.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger_cv.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger_cv.numu_cc_1_direct_length = new std::vector<float>;
  tagger_cv.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger_cv.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger_cv.numu_cc_flag_2 = new std::vector<float>;
  tagger_cv.numu_cc_2_length = new std::vector<float>;
  tagger_cv.numu_cc_2_total_length = new std::vector<float>;
  tagger_cv.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger_cv.numu_cc_2_n_daughter_all = new std::vector<float>;
  tagger_cv.pio_2_v_dis2 = new std::vector<float>;
  tagger_cv.pio_2_v_angle2 = new std::vector<float>;
  tagger_cv.pio_2_v_acc_length = new std::vector<float>;
  tagger_cv.pio_2_v_flag = new std::vector<float>;
  tagger_cv.sig_1_v_angle = new std::vector<float>;
  tagger_cv.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger_cv.sig_1_v_energy = new std::vector<float>;
  tagger_cv.sig_1_v_energy_1 = new std::vector<float>;
  tagger_cv.sig_1_v_flag = new std::vector<float>;
  tagger_cv.sig_2_v_energy = new std::vector<float>;
  tagger_cv.sig_2_v_shower_angle = new std::vector<float>;
  tagger_cv.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger_cv.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger_cv.sig_2_v_flag = new std::vector<float>;
  tagger_cv.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.stw_2_v_energy = new std::vector<float>;
  tagger_cv.stw_2_v_angle = new std::vector<float>;
  tagger_cv.stw_2_v_dir_length = new std::vector<float>;
  tagger_cv.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger_cv.stw_2_v_flag = new std::vector<float>;
  tagger_cv.stw_3_v_angle = new std::vector<float>;
  tagger_cv.stw_3_v_dir_length = new std::vector<float>;
  tagger_cv.stw_3_v_energy = new std::vector<float>;
  tagger_cv.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.stw_3_v_flag = new std::vector<float>;
  tagger_cv.stw_4_v_angle = new std::vector<float>;
  tagger_cv.stw_4_v_dis = new std::vector<float>;
  tagger_cv.stw_4_v_energy = new std::vector<float>;
  tagger_cv.stw_4_v_flag = new std::vector<float>;
  tagger_cv.br3_3_v_energy = new std::vector<float>;
  tagger_cv.br3_3_v_angle = new std::vector<float>;
  tagger_cv.br3_3_v_dir_length = new std::vector<float>;
  tagger_cv.br3_3_v_length = new std::vector<float>;
  tagger_cv.br3_3_v_flag = new std::vector<float>;
  tagger_cv.br3_5_v_dir_length = new std::vector<float>;
  tagger_cv.br3_5_v_total_length = new std::vector<float>;
  tagger_cv.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger_cv.br3_5_v_n_seg = new std::vector<float>;
  tagger_cv.br3_5_v_angle = new std::vector<float>;
  tagger_cv.br3_5_v_sg_length = new std::vector<float>;
  tagger_cv.br3_5_v_energy = new std::vector<float>;
  tagger_cv.br3_5_v_n_main_segs = new std::vector<float>;
  tagger_cv.br3_5_v_n_segs = new std::vector<float>;
  tagger_cv.br3_5_v_shower_main_length = new std::vector<float>;
  tagger_cv.br3_5_v_shower_total_length = new std::vector<float>;
  tagger_cv.br3_5_v_flag = new std::vector<float>;
  tagger_cv.br3_6_v_angle = new std::vector<float>;
  tagger_cv.br3_6_v_angle1 = new std::vector<float>;
  tagger_cv.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger_cv.br3_6_v_direct_length = new std::vector<float>;
  tagger_cv.br3_6_v_length = new std::vector<float>;
  tagger_cv.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger_cv.br3_6_v_energy = new std::vector<float>;
  tagger_cv.br3_6_v_flag = new std::vector<float>;
  tagger_cv.tro_1_v_particle_type = new std::vector<float>;
  tagger_cv.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger_cv.tro_1_v_min_dis = new std::vector<float>;
  tagger_cv.tro_1_v_sg1_length = new std::vector<float>;
  tagger_cv.tro_1_v_shower_main_length = new std::vector<float>;
  tagger_cv.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger_cv.tro_1_v_tmp_length = new std::vector<float>;
  tagger_cv.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger_cv.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger_cv.tro_1_v_flag = new std::vector<float>;
  tagger_cv.tro_2_v_energy = new std::vector<float>;
  tagger_cv.tro_2_v_stem_length = new std::vector<float>;
  tagger_cv.tro_2_v_iso_angle = new std::vector<float>;
  tagger_cv.tro_2_v_max_length = new std::vector<float>;
  tagger_cv.tro_2_v_angle = new std::vector<float>;
  tagger_cv.tro_2_v_flag = new std::vector<float>;
  tagger_cv.tro_4_v_dir2_mag = new std::vector<float>;
  tagger_cv.tro_4_v_angle = new std::vector<float>;
  tagger_cv.tro_4_v_angle1 = new std::vector<float>;
  tagger_cv.tro_4_v_angle2 = new std::vector<float>;
  tagger_cv.tro_4_v_length = new std::vector<float>;
  tagger_cv.tro_4_v_length1 = new std::vector<float>;
  tagger_cv.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger_cv.tro_4_v_energy = new std::vector<float>;
  tagger_cv.tro_4_v_shower_main_length = new std::vector<float>;
  tagger_cv.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger_cv.tro_4_v_flag = new std::vector<float>;
  tagger_cv.tro_5_v_max_angle = new std::vector<float>;
  tagger_cv.tro_5_v_min_angle = new std::vector<float>;
  tagger_cv.tro_5_v_max_length = new std::vector<float>;
  tagger_cv.tro_5_v_iso_angle = new std::vector<float>;
  tagger_cv.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger_cv.tro_5_v_min_count = new std::vector<float>;
  tagger_cv.tro_5_v_max_count = new std::vector<float>;
  tagger_cv.tro_5_v_energy = new std::vector<float>;
  tagger_cv.tro_5_v_flag = new std::vector<float>;
  tagger_cv.lol_1_v_energy = new std::vector<float>;
  tagger_cv.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger_cv.lol_1_v_nseg = new std::vector<float>;
  tagger_cv.lol_1_v_angle = new std::vector<float>;
  tagger_cv.lol_1_v_flag = new std::vector<float>;
  tagger_cv.lol_2_v_length = new std::vector<float>;
  tagger_cv.lol_2_v_angle = new std::vector<float>;
  tagger_cv.lol_2_v_type = new std::vector<float>;
  tagger_cv.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger_cv.lol_2_v_energy = new std::vector<float>;
  tagger_cv.lol_2_v_shower_main_length = new std::vector<float>;
  tagger_cv.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger_cv.lol_2_v_flag = new std::vector<float>;
  tagger_cv.cosmict_flag_10 = new std::vector<float>;
  tagger_cv.cosmict_10_flag_inside = new std::vector<float>;
  tagger_cv.cosmict_10_vtx_z = new std::vector<float>;
  tagger_cv.cosmict_10_flag_shower = new std::vector<float>;
  tagger_cv.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger_cv.cosmict_10_angle_beam = new std::vector<float>;
  tagger_cv.cosmict_10_length = new std::vector<float>;
  tagger_cv.numu_cc_flag_1 = new std::vector<float>;
  tagger_cv.numu_cc_1_particle_type = new std::vector<float>;
  tagger_cv.numu_cc_1_length = new std::vector<float>;
  tagger_cv.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger_cv.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger_cv.numu_cc_1_direct_length = new std::vector<float>;
  tagger_cv.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger_cv.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger_cv.numu_cc_flag_2 = new std::vector<float>;
  tagger_cv.numu_cc_2_length = new std::vector<float>;
  tagger_cv.numu_cc_2_total_length = new std::vector<float>;
  tagger_cv.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger_cv.numu_cc_2_n_daughter_all = new std::vector<float>;

  set_tree_address(T_BDTvars_cv, tagger_cv,2 );
  
  if (flag_data){
    set_tree_address(T_eval_cv, eval_cv,2);
    
    set_tree_address(T_PFeval_cv, pfeval_cv,2);
  }else{
    set_tree_address(T_eval_cv, eval_cv);
    
    set_tree_address(T_PFeval_cv, pfeval_cv);
  }
  
  set_tree_address(T_pot_cv, pot_cv);
  set_tree_address(T_KINEvars_cv, kine_cv);
  
  T_BDTvars_cv->SetBranchStatus("*",0);
  T_BDTvars_cv->SetBranchStatus("numu_cc_flag",1);
  T_BDTvars_cv->SetBranchStatus("numu_score",1);
  T_BDTvars_cv->SetBranchStatus("nue_score",1);
  T_BDTvars_cv->SetBranchStatus("cosmict_flag",1);

  T_eval_cv->SetBranchStatus("*",0);
  T_eval_cv->SetBranchStatus("match_isFC",1);
  T_eval_cv->SetBranchStatus("match_found",1);
  if (T_eval_cv->GetBranch("match_found_asInt")) T_eval_cv->SetBranchStatus("match_found_asInt",1); 
  T_eval_cv->SetBranchStatus("stm_eventtype",1);
  T_eval_cv->SetBranchStatus("stm_lowenergy",1);
  T_eval_cv->SetBranchStatus("stm_LM",1);
  T_eval_cv->SetBranchStatus("stm_TGM",1);
  T_eval_cv->SetBranchStatus("stm_STM",1);
  T_eval_cv->SetBranchStatus("stm_FullDead",1);
  T_eval_cv->SetBranchStatus("stm_clusterlength",1);
  T_eval_cv->SetBranchStatus("match_energy",1);

  if (!flag_data){
    T_eval_cv->SetBranchStatus("weight_spline",1);
    T_eval_cv->SetBranchStatus("weight_cv",1);
    T_eval_cv->SetBranchStatus("weight_lee",1);
    T_eval_cv->SetBranchStatus("weight_change",1);
    // MC enable truth information ...
    T_eval_cv->SetBranchStatus("truth_isCC",1);
    T_eval_cv->SetBranchStatus("truth_nuPdg",1);
    T_eval_cv->SetBranchStatus("truth_vtxInside",1);
    T_eval_cv->SetBranchStatus("truth_nuEnergy",1);
    T_eval_cv->SetBranchStatus("truth_vtxX",1);
    T_eval_cv->SetBranchStatus("truth_vtxY",1);
    T_eval_cv->SetBranchStatus("truth_vtxZ",1);

    T_eval_cv->SetBranchStatus("match_completeness_energy",1);
    T_eval_cv->SetBranchStatus("truth_energyInside",1);
  }
  T_eval_cv->SetBranchStatus("run",1);
  T_eval_cv->SetBranchStatus("subrun",1);
  T_eval_cv->SetBranchStatus("event",1);

  
  T_KINEvars_cv->SetBranchStatus("*",0);
  T_KINEvars_cv->SetBranchStatus("kine_reco_Enu",1);
  T_KINEvars_cv->SetBranchStatus("kine_pio_mass",1);
  T_KINEvars_cv->SetBranchStatus("kine_pio_flag",1);
  T_KINEvars_cv->SetBranchStatus("kine_pio_vtx_dis",1);
  T_KINEvars_cv->SetBranchStatus("kine_pio_energy_1",1);
  T_KINEvars_cv->SetBranchStatus("kine_pio_theta_1",1);
  T_KINEvars_cv->SetBranchStatus("kine_pio_phi_1",1);
  T_KINEvars_cv->SetBranchStatus("kine_pio_dis_1",1);
  T_KINEvars_cv->SetBranchStatus("kine_pio_energy_2",1);
  T_KINEvars_cv->SetBranchStatus("kine_pio_theta_2",1);
  T_KINEvars_cv->SetBranchStatus("kine_pio_phi_2",1);
  T_KINEvars_cv->SetBranchStatus("kine_pio_dis_2",1);
  T_KINEvars_cv->SetBranchStatus("kine_pio_angle",1);

  T_PFeval_cv->SetBranchStatus("*",0);
  T_PFeval_cv->SetBranchStatus("reco_nuvtxX",1);
  T_PFeval_cv->SetBranchStatus("reco_nuvtxY",1);
  T_PFeval_cv->SetBranchStatus("reco_nuvtxZ",1);
  T_PFeval_cv->SetBranchStatus("reco_muonMomentum",1);
  T_PFeval_cv->SetBranchStatus("reco_showerKE",1);
  if (!flag_data){
    T_PFeval_cv->SetBranchStatus("nuvtx_diff",1);
    T_PFeval_cv->SetBranchStatus("showervtx_diff",1);
    T_PFeval_cv->SetBranchStatus("muonvtx_diff",1);
  }
  if (pfeval_cv.flag_NCDelta){
    T_PFeval_cv->SetBranchStatus("reco_protonMomentum",1);
    if (!flag_data){
      T_PFeval_cv->SetBranchStatus("truth_NCDelta",1);
      T_PFeval_cv->SetBranchStatus("truth_NprimPio",1);
    }
    
  }
  

  std::map<std::pair<int, int>, int> map_re_entry_cv;
  std::map<std::pair<int, int>, std::set<std::pair<int, int> > > map_rs_re_cv;

  bool flag_presel = false;

  int num_check = 0;
  
  for (int i=0;i!=T_eval_cv->GetEntries();i++){
    T_eval_cv->GetEntry(i);
    T_BDTvars_cv->GetEntry(i);
    
    map_rs_re_cv[std::make_pair(eval_cv.run, eval_cv.subrun)].insert(std::make_pair(eval_cv.run, eval_cv.event));
    
    int tmp_match_found = eval_cv.match_found;
    if (eval_cv.is_match_found_int) tmp_match_found = eval_cv.match_found_asInt;
    
    flag_presel = false;
    if (tmp_match_found != 0 && eval_cv.stm_eventtype != 0 && eval_cv.stm_lowenergy ==0 && eval_cv.stm_LM ==0 && eval_cv.stm_TGM ==0 && eval_cv.stm_STM==0 && eval_cv.stm_FullDead == 0 && eval_cv.stm_clusterlength >0) {
      flag_presel = true; // preselection ...
    }
    
    if (tmp_match_found == -1  || (tmp_match_found == 1 && eval_cv.stm_lowenergy == -1) || (flag_presel && tagger_cv.numu_cc_flag == -1)) {
      //num_check ++;
      continue;
    }

    //if (eval_cv.run == 7486 && eval_cv.event==3964) std::cout << flag_presel << " " << tagger_cv.numu_cc_flag << std::endl;
    
    map_re_entry_cv[std::make_pair(eval_cv.run, eval_cv.event)] = i;
  }

  //  std::cout << num_check << std::endl;
  
  
  
  
  std::map<std::pair<int, int>, std::pair<int, double> > map_rs_entry_pot_cv;
  for (Int_t i=0;i!=T_pot_cv->GetEntries();i++){
    T_pot_cv->GetEntry(i);
    map_rs_entry_pot_cv[std::make_pair(pot_cv.runNo,pot_cv.subRunNo)] = std::make_pair(i, pot_cv.pot_tor875);
  }


  
 
    
  // run, subrun ...
  std::set<std::pair<int,int> > remove_set;
  std::map<std::pair<int,int>, int> map_rs_failed;
  for (auto it =  map_rs_re_cv.begin(); it !=  map_rs_re_cv.end(); it++){
    int failed_num = 0;
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      if (map_re_entry_cv.find(*it1) == map_re_entry_cv.end()) failed_num++;
    }
    if (failed_num > it->second.size() * fail_percentage && failed_num !=1
	|| failed_num > it->second.size() * 0.33 && failed_num==1){
      remove_set.insert(it->first);
    }
    map_rs_failed[it->first] = failed_num;
  }
  //  std::cout << remove_set.size() << std::endl;

  

  TFile *file3 = new TFile(out_file,"RECREATE");
  file3->mkdir("wcpselection");

  file3->cd("wcpselection");
  
  TTree *t1_cv = new TTree("T_eval","T_eval");
  TTree *t2_cv = new TTree("T_pot","T_pot");
  TTree *t3_cv = new TTree("T_PFeval", "T_PFeval");
  TTree *t5_cv = new TTree("T_KINEvars", "T_KINEvars");
  TTree *t4_cv = new TTree("T_BDTvars","T_BDTvars");

  put_tree_address(t4_cv, tagger_cv,2);

  if (flag_data){
    put_tree_address(t1_cv, eval_cv,2);
    
    put_tree_address(t3_cv, pfeval_cv,2);
  }else{
    put_tree_address(t1_cv, eval_cv);
    
    put_tree_address(t3_cv, pfeval_cv);
  }
  
  put_tree_address(t2_cv, pot_cv);
  put_tree_address(t5_cv, kine_cv);
  

  for (auto it = map_re_entry_cv.begin(); it != map_re_entry_cv.end(); it++){
    T_BDTvars_cv->GetEntry(it->second);
    T_eval_cv->GetEntry(it->second);
    T_KINEvars_cv->GetEntry(it->second);
    T_PFeval_cv->GetEntry(it->second);

    if (remove_set.find(std::make_pair(eval_cv.run, eval_cv.subrun)) != remove_set.end()) continue;

    t1_cv->Fill();
    t3_cv->Fill();
    t4_cv->Fill();
    t5_cv->Fill();
  }

  double cv_pot=0;
  double cv1_pot = 0;
  float pass_ratio;
  t2_cv->Branch("pass_ratio",&pass_ratio,"pass_ratio/F");
  
  for (auto it = map_rs_entry_pot_cv.begin(); it != map_rs_entry_pot_cv.end(); it++){
    T_pot_cv->GetEntry(it->second.first);
    cv_pot += it->second.second;

    if(remove_set.find(it->first) != remove_set.end()) continue;
    if (map_rs_re_cv[it->first].size()==0) {
      //continue;
      pass_ratio = 1;
      //      std::cout << pot_cv.runNo << " " << pot_cv.subRunNo << " " << pot_cv.pot_tor875 << std::endl;
    }else{
      pass_ratio = 1 - map_rs_failed[it->first] * 1.0 / map_rs_re_cv[it->first].size();
    }
    cv1_pot += it->second.second * pass_ratio;

    pot_cv.pot_tor875 *= pass_ratio;
    pot_cv.pot_tor875good *= pass_ratio;
    
    t2_cv->Fill();
  }

  std::cout << out_file << std::endl;
  std::cout << "Events: " << t1_cv->GetEntries()<<"/"<<T_eval_cv->GetEntries() << std::endl;
  std::cout << "POT:    " << cv1_pot << " " << cv_pot << std::endl;
  
  file3->Write();
  file3->Close();
  
  
  return 0;
}
