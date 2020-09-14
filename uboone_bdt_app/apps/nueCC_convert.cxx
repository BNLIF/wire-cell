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

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

using namespace std;


float cal_br3_3_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& br3_3_v_energy,
		    float& br3_3_v_angle,
		    float& br3_3_v_dir_length,
		    float& br3_3_v_length);
float cal_br3_5_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& br3_5_v_dir_length,
		    float& br3_5_v_total_length,
		    float& br3_5_v_flag_avoid_muon_check,
		    float& br3_5_v_n_seg,
		    float& br3_5_v_angle,
		    float& br3_5_v_sg_length,
		    float& br3_5_v_energy,
		    float& br3_5_v_n_main_segs,
		    float& br3_5_v_n_segs,
		    float& br3_5_v_shower_main_length,
		    float& br3_5_v_shower_total_length);
float cal_br3_6_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& br3_6_v_angle,
		    float& br3_6_v_angle1,
		    float& br3_6_v_flag_shower_trajectory,
		    float& br3_6_v_direct_length,
		    float& br3_6_v_length,
		    float& br3_6_v_n_other_vtx_segs,
		    float& br3_6_v_energy);
float cal_pio_2_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& pio_2_v_dis2,
		    float& pio_2_v_angle2,
		    float& pio_2_v_acc_length);
float cal_stw_2_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& stw_2_v_medium_dQ_dx,
		    float& stw_2_v_energy,
		    float& stw_2_v_angle,
		    float& stw_2_v_dir_length,
		    float& stw_2_v_max_dQ_dx);
float cal_stw_3_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& stw_3_v_angle,
		    float& stw_3_v_dir_length,
		    float& stw_3_v_energy,
		    float& stw_3_v_medium_dQ_dx);
float cal_stw_4_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& stw_4_v_angle,
		    float& stw_4_v_dis,
		    float& stw_4_v_energy);
float cal_sig_1_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& sig_1_v_angle,
		    float& sig_1_v_flag_single_shower,
		    float& sig_1_v_energy,
		    float& sig_1_v_energy_1);
float cal_sig_2_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& sig_2_v_energy,
		    float& sig_2_v_shower_angle,
		    float& sig_2_v_flag_single_shower,
		    float& sig_2_v_medium_dQ_dx,
		    float& sig_2_v_start_dQ_dx);
float cal_lol_1_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& lol_1_v_energy,
		    float& lol_1_v_vtx_n_segs,
		    float& lol_1_v_nseg,
		    float& lol_1_v_angle);
float cal_lol_2_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& lol_2_v_length,
		    float& lol_2_v_angle,
		    float& lol_2_v_type,
		    float& lol_2_v_vtx_n_segs,
		    float& lol_2_v_energy,
		    float& lol_2_v_shower_main_length,
		    float& lol_2_v_flag_dir_weak);
float cal_tro_1_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& tro_1_v_particle_type,
		    float& tro_1_v_flag_dir_weak,
		    float& tro_1_v_min_dis,
		    float& tro_1_v_sg1_length,
		    float& tro_1_v_shower_main_length,
		    float& tro_1_v_max_n_vtx_segs,
		    float& tro_1_v_tmp_length,
		    float& tro_1_v_medium_dQ_dx,
		    float& tro_1_v_dQ_dx_cut,
		    float& tro_1_v_flag_shower_topology);
float cal_tro_2_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& tro_2_v_energy,
		    float& tro_2_v_stem_length,
		    float& tro_2_v_iso_angle,
		    float& tro_2_v_max_length,
		    float& tro_2_v_angle);
float cal_tro_4_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& tro_4_v_dir2_mag,
		    float& tro_4_v_angle,
		    float& tro_4_v_angle1,
		    float& tro_4_v_angle2,
		    float& tro_4_v_length,
		    float& tro_4_v_length1,
		    float& tro_4_v_medium_dQ_dx,
		    float& tro_4_v_end_dQ_dx,
		    float& tro_4_v_energy,
		    float& tro_4_v_shower_main_length,
		    float& tro_4_v_flag_shower_trajectory);
float cal_tro_5_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& tro_5_v_max_angle,
		    float& tro_5_v_min_angle,
		    float& tro_5_v_max_length,
		    float& tro_5_v_iso_angle,
		    float& tro_5_v_n_vtx_segs,
		    float& tro_5_v_min_count,
		    float& tro_5_v_max_count,
		    float& tro_5_v_energy);


float cal_br3_3_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& br3_3_v_energy,
		    float& br3_3_v_angle,
		    float& br3_3_v_dir_length,
		    float& br3_3_v_length){
  float val = default_val;
  if (tagger_info.br3_3_v_angle->size()>0){
    val = 1e9;
    for (size_t i = 0; i!= tagger_info.br3_3_v_energy->size(); i++){
      br3_3_v_energy = tagger_info.br3_3_v_energy->at(i);
      br3_3_v_angle = tagger_info.br3_3_v_angle->at(i);
      br3_3_v_dir_length = tagger_info.br3_3_v_dir_length->at(i);
      br3_3_v_length = tagger_info.br3_3_v_length->at(i);
      
      float tmp_val = reader.EvaluateMVA("MyBDT");
      if (tmp_val < val)     val = tmp_val;
    }
  }
  
  return val;
}
float cal_br3_5_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& br3_5_v_dir_length,
		    float& br3_5_v_total_length,
		    float& br3_5_v_flag_avoid_muon_check,
		    float& br3_5_v_n_seg,
		    float& br3_5_v_angle,
		    float& br3_5_v_sg_length,
		    float& br3_5_v_energy,
		    float& br3_5_v_n_main_segs,
		    float& br3_5_v_n_segs,
		    float& br3_5_v_shower_main_length,
		    float& br3_5_v_shower_total_length){
  float val = default_val;

  if (tagger_info.br3_5_v_dir_length->size()>0){
    val = 1e9;
    for (size_t i=0;i!=tagger_info.br3_5_v_dir_length->size();i++){
      br3_5_v_dir_length = tagger_info.br3_5_v_dir_length->at(i);
      br3_5_v_total_length = tagger_info.br3_5_v_total_length->at(i);
      br3_5_v_flag_avoid_muon_check = tagger_info.br3_5_v_flag_avoid_muon_check->at(i);
      br3_5_v_n_seg = tagger_info.br3_5_v_n_seg->at(i);
      br3_5_v_angle = tagger_info.br3_5_v_angle->at(i);
      br3_5_v_sg_length = tagger_info.br3_5_v_sg_length->at(i);
      br3_5_v_energy = tagger_info.br3_5_v_energy->at(i);
      br3_5_v_n_segs = tagger_info.br3_5_v_n_segs->at(i);
      br3_5_v_shower_main_length = tagger_info.br3_5_v_shower_main_length->at(i);
      br3_5_v_shower_total_length = tagger_info.br3_5_v_shower_total_length->at(i);      
      
      float tmp_val = reader.EvaluateMVA("MyBDT");
      if (tmp_val < val) val = tmp_val;
    }

  }
  
  return val;
}

float cal_br3_6_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& br3_6_v_angle,
		    float& br3_6_v_angle1,
		    float& br3_6_v_flag_shower_trajectory,
		    float& br3_6_v_direct_length,
		    float& br3_6_v_length,
		    float& br3_6_v_n_other_vtx_segs,
		    float& br3_6_v_energy){
  float val = default_val;
  if (tagger_info.br3_6_v_angle->size() > 0){
    val = 1e9;
    for (size_t i=0; i!= tagger_info.br3_6_v_angle->size();i++){
      br3_6_v_angle = tagger_info.br3_6_v_angle->at(i);
      br3_6_v_angle1 = tagger_info.br3_6_v_angle1->at(i);
      br3_6_v_flag_shower_trajectory = tagger_info.br3_6_v_flag_shower_trajectory->at(i);
      br3_6_v_direct_length = tagger_info.br3_6_v_direct_length->at(i);
      br3_6_v_length = tagger_info.br3_6_v_length->at(i);
      br3_6_v_n_other_vtx_segs = tagger_info.br3_6_v_n_other_vtx_segs->at(i);
      br3_6_v_energy = tagger_info.br3_6_v_energy->at(i);
      
      float tmp_val = reader.EvaluateMVA("MyBDT");
      if (tmp_val < val) val = tmp_val;
    }
  }
  
  return val;
}
float cal_pio_2_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& pio_2_v_dis2,
		    float& pio_2_v_angle2,
		    float& pio_2_v_acc_length){
  float val = default_val;
  if (tagger_info.pio_2_v_dis2->size()>0){
    val = 1e9;
    for (size_t i=0;i!=tagger_info.pio_2_v_dis2->size();i++){
      pio_2_v_dis2 = tagger_info.pio_2_v_dis2->at(i);
      pio_2_v_angle2 = tagger_info.pio_2_v_angle2->at(i);
      pio_2_v_acc_length = tagger_info.pio_2_v_acc_length->at(i);
      
      float tmp_val = reader.EvaluateMVA("MyBDT");
      if (tmp_val < val) val = tmp_val;
    }
  }
  
  return val;
}
float cal_stw_2_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& stw_2_v_medium_dQ_dx,
		    float& stw_2_v_energy,
		    float& stw_2_v_angle,
		    float& stw_2_v_dir_length,
		    float& stw_2_v_max_dQ_dx){
  float val = default_val;

  if (tagger_info.stw_2_v_energy->size()>0){
    val = 1e9;
    for (size_t i=0;i!=tagger_info.stw_2_v_medium_dQ_dx->size();i++){
      stw_2_v_medium_dQ_dx = tagger_info.stw_2_v_medium_dQ_dx->at(i);
      stw_2_v_energy = tagger_info.stw_2_v_energy->at(i);
      stw_2_v_angle = tagger_info.stw_2_v_angle->at(i);
      stw_2_v_dir_length = tagger_info.stw_2_v_dir_length->at(i);
      stw_2_v_max_dQ_dx = tagger_info.stw_2_v_max_dQ_dx->at(i);
      
      float tmp_val = reader.EvaluateMVA("MyBDT");
      if (tmp_val < val) val = tmp_val;
    }
  }
  
  return val;
}
float cal_stw_3_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& stw_3_v_angle,
		    float& stw_3_v_dir_length,
		    float& stw_3_v_energy,
		    float& stw_3_v_medium_dQ_dx){
  float val = default_val;

  if (tagger_info.stw_3_v_angle->size()>0){
    val = 1e9;
    for(size_t i=0;i!=tagger_info.stw_3_v_angle->size();i++){
      stw_3_v_angle = tagger_info.stw_3_v_angle->at(i);
      stw_3_v_dir_length = tagger_info.stw_3_v_dir_length->at(i);
      stw_3_v_energy = tagger_info.stw_3_v_energy->at(i);
      stw_3_v_medium_dQ_dx = tagger_info.stw_3_v_medium_dQ_dx->at(i);
      
      float tmp_val = reader.EvaluateMVA("MyBDT");
      if (tmp_val < val) val = tmp_val;
    }
  }
  
  return val;
}
float cal_stw_4_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& stw_4_v_angle,
		    float& stw_4_v_dis,
		    float& stw_4_v_energy){
  float val = default_val;
  if (tagger_info.stw_4_v_angle->size()>0){
    val = 1e9;
    for (size_t i=0;i!=tagger_info.stw_4_v_angle->size();i++){
      stw_4_v_angle = tagger_info.stw_4_v_angle->at(i);
      stw_4_v_dis = tagger_info.stw_4_v_dis->at(i);
      stw_4_v_energy = tagger_info.stw_4_v_energy->at(i);
      
      float tmp_val = reader.EvaluateMVA("MyBDT");
      if (tmp_val < val) val = tmp_val;
    }
  }
  return val;
}
float cal_sig_1_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& sig_1_v_angle,
		    float& sig_1_v_flag_single_shower,
		    float& sig_1_v_energy,
		    float& sig_1_v_energy_1){
  float val = default_val;
  if (tagger_info.sig_1_v_angle->size()>0){
    val = 1e9;
    for (size_t i=0;i!=tagger_info.sig_1_v_angle->size();i++){
      sig_1_v_angle = tagger_info.sig_1_v_angle->at(i);
      sig_1_v_flag_single_shower = tagger_info.sig_1_v_flag_single_shower->at(i);
      sig_1_v_energy = tagger_info.sig_1_v_energy->at(i);
      sig_1_v_energy_1 = tagger_info.sig_1_v_energy_1->at(i);
      
      float tmp_val = reader.EvaluateMVA("MyBDT");
      if (tmp_val < val) val = tmp_val;
    }
  }
  return val;
}
float cal_sig_2_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& sig_2_v_energy,
		    float& sig_2_v_shower_angle,
		    float& sig_2_v_flag_single_shower,
		    float& sig_2_v_medium_dQ_dx,
		    float& sig_2_v_start_dQ_dx){
  float val = default_val;
  if (tagger_info.sig_2_v_energy->size()>0){
    val = 1e9;
    for (size_t i=0; i!= tagger_info.sig_2_v_energy->size();i++){
      sig_2_v_energy = tagger_info.sig_2_v_energy->at(i);
      sig_2_v_shower_angle = tagger_info.sig_2_v_shower_angle->at(i);
      sig_2_v_flag_single_shower = tagger_info.sig_2_v_flag_single_shower->at(i);
      sig_2_v_medium_dQ_dx = tagger_info.sig_2_v_medium_dQ_dx->at(i);
      sig_2_v_start_dQ_dx = tagger_info.sig_2_v_start_dQ_dx->at(i);
      
      float tmp_val = reader.EvaluateMVA("MyBDT");
      if (tmp_val < val) val = tmp_val;
    }
  }
  
  return val;
}
float cal_lol_1_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& lol_1_v_energy,
		    float& lol_1_v_vtx_n_segs,
		    float& lol_1_v_nseg,
		    float& lol_1_v_angle){
  float val = default_val;
  if (tagger_info.lol_1_v_energy->size()>0){
    val = 1e9;
    for (size_t i=0;i!=tagger_info.lol_1_v_energy->size(); i++){
      lol_1_v_energy = tagger_info.lol_1_v_energy->at(i);
      lol_1_v_vtx_n_segs = tagger_info.lol_1_v_vtx_n_segs->at(i);
      lol_1_v_nseg = tagger_info.lol_1_v_nseg->at(i);
      lol_1_v_angle = tagger_info.lol_1_v_angle->at(i);
      
      float tmp_val = reader.EvaluateMVA("MyBDT");
      if (tmp_val < val) val = tmp_val;
    }
  }
  return val;
}
float cal_lol_2_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& lol_2_v_length,
		    float& lol_2_v_angle,
		    float& lol_2_v_type,
		    float& lol_2_v_vtx_n_segs,
		    float& lol_2_v_energy,
		    float& lol_2_v_shower_main_length,
		    float& lol_2_v_flag_dir_weak){
  float val = default_val;
  if (tagger_info.lol_2_v_length->size()>0){
    val = 1e9;
    for (size_t i=0;i!=tagger_info.lol_2_v_length->size();i++){
      lol_2_v_length = tagger_info.lol_2_v_length->at(i);
      lol_2_v_angle = tagger_info.lol_2_v_angle->at(i);
      lol_2_v_type = tagger_info.lol_2_v_type->at(i);
      lol_2_v_vtx_n_segs = tagger_info.lol_2_v_vtx_n_segs->at(i);
      lol_2_v_energy = tagger_info.lol_2_v_energy->at(i);
      lol_2_v_shower_main_length = tagger_info.lol_2_v_shower_main_length->at(i);
      lol_2_v_flag_dir_weak = tagger_info.lol_2_v_flag_dir_weak->at(i);
      
      float tmp_val = reader.EvaluateMVA("MyBDT");
      if (tmp_val < val) val = tmp_val;
    }
  }
  
  return val;
}
float cal_tro_1_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& tro_1_v_particle_type,
		    float& tro_1_v_flag_dir_weak,
		    float& tro_1_v_min_dis,
		    float& tro_1_v_sg1_length,
		    float& tro_1_v_shower_main_length,
		    float& tro_1_v_max_n_vtx_segs,
		    float& tro_1_v_tmp_length,
		    float& tro_1_v_medium_dQ_dx,
		    float& tro_1_v_dQ_dx_cut,
		    float& tro_1_v_flag_shower_topology){
  float val = default_val;
  if (tagger_info.tro_1_v_dQ_dx_cut->size()>0){
    val = 1e9;
    for (size_t i=0;i!=tagger_info.tro_1_v_particle_type->size();i++){
      tro_1_v_particle_type = tagger_info.tro_1_v_particle_type->at(i);
      tro_1_v_flag_dir_weak = tagger_info.tro_1_v_flag_dir_weak->at(i);
      tro_1_v_min_dis = tagger_info.tro_1_v_min_dis->at(i);
      tro_1_v_sg1_length = tagger_info.tro_1_v_sg1_length->at(i);
      tro_1_v_shower_main_length = tagger_info.tro_1_v_shower_main_length->at(i);
      tro_1_v_max_n_vtx_segs = tagger_info.tro_1_v_max_n_vtx_segs->at(i);
      tro_1_v_tmp_length = tagger_info.tro_1_v_tmp_length->at(i);
      tro_1_v_medium_dQ_dx = tagger_info.tro_1_v_medium_dQ_dx->at(i);
      tro_1_v_dQ_dx_cut = tagger_info.tro_1_v_dQ_dx_cut->at(i);
      tro_1_v_flag_shower_topology = tagger_info.tro_1_v_flag_shower_topology->at(i);
      
      float tmp_val = reader.EvaluateMVA("MyBDT");
      if (tmp_val < val) val = tmp_val;
    }
  }
  return val;
}
float cal_tro_2_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& tro_2_v_energy,
		    float& tro_2_v_stem_length,
		    float& tro_2_v_iso_angle,
		    float& tro_2_v_max_length,
		    float& tro_2_v_angle){
  float val = default_val;
  if (tagger_info.tro_2_v_energy->size()>0){
    val = 1e9;
    for (size_t i=0;i!=tagger_info.tro_2_v_energy->size();i++){

      tro_2_v_energy = tagger_info.tro_2_v_energy->at(i);
      tro_2_v_stem_length = tagger_info.tro_2_v_stem_length->at(i);
      tro_2_v_iso_angle = tagger_info.tro_2_v_iso_angle->at(i);
      tro_2_v_max_length = tagger_info.tro_2_v_max_length->at(i);
      tro_2_v_angle = tagger_info.tro_2_v_angle->at(i);
      
      float tmp_val = reader.EvaluateMVA("MyBDT");
      if (tmp_val < val) val = tmp_val;
    }
  }
  return val;
}
float cal_tro_4_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& tro_4_v_dir2_mag,
		    float& tro_4_v_angle,
		    float& tro_4_v_angle1,
		    float& tro_4_v_angle2,
		    float& tro_4_v_length,
		    float& tro_4_v_length1,
		    float& tro_4_v_medium_dQ_dx,
		    float& tro_4_v_end_dQ_dx,
		    float& tro_4_v_energy,
		    float& tro_4_v_shower_main_length,
		    float& tro_4_v_flag_shower_trajectory){
  float val = default_val;
  if (tagger_info.tro_4_v_angle->size()>0){
    val = 1e9;
     for (size_t i=0; i!= tagger_info.tro_4_v_dir2_mag->size(); i++){

       tro_4_v_dir2_mag = tagger_info.tro_4_v_dir2_mag->at(i);
       tro_4_v_angle = tagger_info.tro_4_v_angle->at(i);
       tro_4_v_angle1 = tagger_info.tro_4_v_angle1->at(i);
       tro_4_v_angle2 = tagger_info.tro_4_v_angle2->at(i);
       tro_4_v_length = tagger_info.tro_4_v_length->at(i);
       tro_4_v_length1 = tagger_info.tro_4_v_length1->at(i);
       tro_4_v_medium_dQ_dx = tagger_info.tro_4_v_medium_dQ_dx->at(i);
       tro_4_v_end_dQ_dx = tagger_info.tro_4_v_end_dQ_dx->at(i);
       tro_4_v_energy = tagger_info.tro_4_v_energy->at(i);
       tro_4_v_shower_main_length = tagger_info.tro_4_v_shower_main_length->at(i);
       tro_4_v_flag_shower_trajectory = tagger_info.tro_4_v_flag_shower_trajectory->at(i);
       
       float tmp_val = reader.EvaluateMVA("MyBDT");
       if (tmp_val < val) val = tmp_val;
     }
  }
  return val;
}
float cal_tro_5_bdt(float default_val , TaggerInfo& tagger_info, TMVA::Reader& reader,float& tro_5_v_max_angle,
		    float& tro_5_v_min_angle,
		    float& tro_5_v_max_length,
		    float& tro_5_v_iso_angle,
		    float& tro_5_v_n_vtx_segs,
		    float& tro_5_v_min_count,
		    float& tro_5_v_max_count,
		    float& tro_5_v_energy){
  float val = default_val;
  if (tagger_info.tro_5_v_energy->size()>0){
    val = 1e9;
    for (size_t i=0;i!=tagger_info.tro_5_v_max_angle->size();i++){

      tro_5_v_max_angle = tagger_info.tro_5_v_max_angle->at(i);
      tro_5_v_min_angle = tagger_info.tro_5_v_min_angle->at(i);
      tro_5_v_max_length = tagger_info.tro_5_v_max_length->at(i);
      tro_5_v_iso_angle = tagger_info.tro_5_v_iso_angle->at(i);
      tro_5_v_n_vtx_segs = tagger_info.tro_5_v_n_vtx_segs->at(i);
      tro_5_v_min_count = tagger_info.tro_5_v_min_count->at(i);
      tro_5_v_max_count = tagger_info.tro_5_v_max_count->at(i);
      tro_5_v_energy = tagger_info.tro_5_v_energy->at(i);
      
      float tmp_val = reader.EvaluateMVA("MyBDT");
      if (tmp_val < val) val = tmp_val;
    }
  }
  return val;
}

int main( int argc, char** argv )
{
  int process = 1;

  for (Int_t i=1;i!=argc;i++){
    switch(argv[i][1]){
    case 'p':
      process = atoi(&argv[i][2]);
      break;
    }
  }
  
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

  TString filename;
  if (process==1){
    filename = "bdt.root";
  }else{
    filename = "bdt_validation.root";
  }
  
  TFile *new_file = new TFile(filename,"RECREATE");
  TTree* Tsig = new TTree("sig","sig");
  TTree* Tbkg = new TTree("bkg","bkg");
  Tsig->SetDirectory(new_file);
  Tbkg->SetDirectory(new_file);
  put_tree_address(Tsig, tagger);
  put_tree_address(Tbkg, tagger);

  // BDT stuff 
  TMVA::Reader reader_br3_3;
  float br3_3_v_energy;
  float br3_3_v_angle;
  float br3_3_v_dir_length;
  float br3_3_v_length;
  reader_br3_3.AddVariable("br3_3_v_energy",&br3_3_v_energy);
  reader_br3_3.AddVariable("br3_3_v_angle",&br3_3_v_angle);
  reader_br3_3.AddVariable("br3_3_v_dir_length",&br3_3_v_dir_length);
  reader_br3_3.AddVariable("br3_3_v_length",&br3_3_v_length);
  reader_br3_3.BookMVA( "MyBDT", "weights/br3_3_BDT.weights.xml");

  TMVA::Reader reader_br3_5;
  float br3_5_v_dir_length;
  float br3_5_v_total_length;
  float br3_5_v_flag_avoid_muon_check;
  float br3_5_v_n_seg;
  float br3_5_v_angle;
  float br3_5_v_sg_length;
  float br3_5_v_energy;
  float br3_5_v_n_main_segs;
  float br3_5_v_n_segs;
  float br3_5_v_shower_main_length;
  float br3_5_v_shower_total_length;
  reader_br3_5.AddVariable("br3_5_v_dir_length",&br3_5_v_dir_length);
  reader_br3_5.AddVariable("br3_5_v_total_length",&br3_5_v_total_length);
  reader_br3_5.AddVariable("br3_5_v_flag_avoid_muon_check",&br3_5_v_flag_avoid_muon_check);
  reader_br3_5.AddVariable("br3_5_v_n_seg",&br3_5_v_n_seg);
  reader_br3_5.AddVariable("br3_5_v_angle",&br3_5_v_angle);
  reader_br3_5.AddVariable("br3_5_v_sg_length",&br3_5_v_sg_length);
  reader_br3_5.AddVariable("br3_5_v_energy",&br3_5_v_energy);
  reader_br3_5.AddVariable("br3_5_v_n_segs",&br3_5_v_n_segs);
  reader_br3_5.AddVariable("br3_5_v_shower_main_length",&br3_5_v_shower_main_length);
  reader_br3_5.AddVariable("br3_5_v_shower_total_length",&br3_5_v_shower_total_length);
  reader_br3_5.BookMVA( "MyBDT", "weights/br3_5_BDT.weights.xml");

  TMVA::Reader reader_br3_6;
  float br3_6_v_angle;
  float br3_6_v_angle1;
  float br3_6_v_flag_shower_trajectory;
  float br3_6_v_direct_length;
  float br3_6_v_length;
  float br3_6_v_n_other_vtx_segs;
  float br3_6_v_energy;
  reader_br3_6.AddVariable("br3_6_v_angle",&br3_6_v_angle);
  reader_br3_6.AddVariable("br3_6_v_angle1",&br3_6_v_angle1);
  reader_br3_6.AddVariable("br3_6_v_flag_shower_trajectory",&br3_6_v_flag_shower_trajectory);
  reader_br3_6.AddVariable("br3_6_v_direct_length",&br3_6_v_direct_length);
  reader_br3_6.AddVariable("br3_6_v_length",&br3_6_v_length);
  reader_br3_6.AddVariable("br3_6_v_n_other_vtx_segs",&br3_6_v_n_other_vtx_segs);
  reader_br3_6.AddVariable("br3_6_v_energy",&br3_6_v_energy);
  reader_br3_6.BookMVA( "MyBDT", "weights/br3_6_BDT.weights.xml");

  TMVA::Reader reader_pio_2;
  float pio_2_v_dis2;
  float pio_2_v_angle2;
  float pio_2_v_acc_length;
  reader_pio_2.AddVariable("pio_2_v_dis2",&pio_2_v_dis2);
  reader_pio_2.AddVariable("pio_2_v_angle2",&pio_2_v_angle2);
  reader_pio_2.AddVariable("pio_2_v_acc_length",&pio_2_v_acc_length);
  reader_pio_2.AddVariable("pio_mip_id",&tagger.pio_mip_id);
  reader_pio_2.BookMVA( "MyBDT", "weights/pio_2_BDT.weights.xml");

  TMVA::Reader reader_lol_1;
  float lol_1_v_energy;
  float lol_1_v_vtx_n_segs;
  float lol_1_v_nseg;
  float lol_1_v_angle;
  reader_lol_1.AddVariable("lol_1_v_energy",&lol_1_v_energy);
  reader_lol_1.AddVariable("lol_1_v_vtx_n_segs",&lol_1_v_vtx_n_segs);
  reader_lol_1.AddVariable("lol_1_v_nseg",&lol_1_v_nseg);
  reader_lol_1.AddVariable("lol_1_v_angle",&lol_1_v_angle);
  reader_lol_1.BookMVA( "MyBDT", "weights/lol_1_BDT.weights.xml");

  TMVA::Reader reader_lol_2;
  float lol_2_v_length;
  float lol_2_v_angle;
  float lol_2_v_type;
  float lol_2_v_vtx_n_segs;
  float lol_2_v_energy;
  float lol_2_v_shower_main_length;
  float lol_2_v_flag_dir_weak;
  reader_lol_2.AddVariable("lol_2_v_length",&lol_2_v_length);
  reader_lol_2.AddVariable("lol_2_v_angle",&lol_2_v_angle);
  reader_lol_2.AddVariable("lol_2_v_type",&lol_2_v_type);
  reader_lol_2.AddVariable("lol_2_v_vtx_n_segs",&lol_2_v_vtx_n_segs);
  reader_lol_2.AddVariable("lol_2_v_energy",&lol_2_v_energy);
  reader_lol_2.AddVariable("lol_2_v_shower_main_length",&lol_2_v_shower_main_length);
  reader_lol_2.AddVariable("lol_2_v_flag_dir_weak",&lol_2_v_flag_dir_weak);
  reader_lol_2.BookMVA( "MyBDT", "weights/lol_2_BDT.weights.xml");

  TMVA::Reader reader_tro_1;
  float tro_1_v_particle_type;
  float tro_1_v_flag_dir_weak;
  float tro_1_v_min_dis;
  float tro_1_v_sg1_length;
  float tro_1_v_shower_main_length;
  float tro_1_v_max_n_vtx_segs;
  float tro_1_v_tmp_length;
  float tro_1_v_medium_dQ_dx;
  float tro_1_v_dQ_dx_cut;
  float tro_1_v_flag_shower_topology;
  reader_tro_1.AddVariable("tro_1_v_particle_type",&tro_1_v_particle_type);
  reader_tro_1.AddVariable("tro_1_v_flag_dir_weak",&tro_1_v_flag_dir_weak);
  reader_tro_1.AddVariable("tro_1_v_min_dis",&tro_1_v_min_dis);
  reader_tro_1.AddVariable("tro_1_v_sg1_length",&tro_1_v_sg1_length);
  reader_tro_1.AddVariable("tro_1_v_shower_main_length",&tro_1_v_shower_main_length);
  reader_tro_1.AddVariable("tro_1_v_max_n_vtx_segs",&tro_1_v_max_n_vtx_segs);
  reader_tro_1.AddVariable("tro_1_v_tmp_length",&tro_1_v_tmp_length);
  reader_tro_1.AddVariable("tro_1_v_medium_dQ_dx",&tro_1_v_medium_dQ_dx);
  reader_tro_1.AddVariable("tro_1_v_dQ_dx_cut", &tro_1_v_dQ_dx_cut);
  reader_tro_1.AddVariable("tro_1_v_flag_shower_topology", &tro_1_v_flag_shower_topology);
  reader_tro_1.BookMVA( "MyBDT", "weights/tro_1_BDT.weights.xml");

  TMVA::Reader reader_tro_2;
  float tro_2_v_energy;
  float tro_2_v_stem_length;
  float tro_2_v_iso_angle;
  float tro_2_v_max_length;
  float tro_2_v_angle;
  reader_tro_2.AddVariable("tro_2_v_energy",&tro_2_v_energy);
  reader_tro_2.AddVariable("tro_2_v_stem_length",&tro_2_v_stem_length);
  reader_tro_2.AddVariable("tro_2_v_iso_angle",&tro_2_v_iso_angle);
  reader_tro_2.AddVariable("tro_2_v_max_length",&tro_2_v_max_length);
  reader_tro_2.AddVariable("tro_2_v_angle",&tro_2_v_angle);
  reader_tro_2.BookMVA( "MyBDT", "weights/tro_2_BDT.weights.xml");

  TMVA::Reader reader_tro_4;

  float tro_4_v_dir2_mag;
  float tro_4_v_angle;
  float tro_4_v_angle1;
  float tro_4_v_angle2;
  float tro_4_v_length;
  float tro_4_v_length1;
  float tro_4_v_medium_dQ_dx;
  float tro_4_v_end_dQ_dx;
  float tro_4_v_energy;
  float tro_4_v_shower_main_length;
  float tro_4_v_flag_shower_trajectory;
  
  reader_tro_4.AddVariable("tro_4_v_dir2_mag",&tro_4_v_dir2_mag);
  reader_tro_4.AddVariable("tro_4_v_angle",&tro_4_v_angle);
  reader_tro_4.AddVariable("tro_4_v_angle1",&tro_4_v_angle1);
  reader_tro_4.AddVariable("tro_4_v_angle2",&tro_4_v_angle2);
  reader_tro_4.AddVariable("tro_4_v_length",&tro_4_v_length);
  reader_tro_4.AddVariable("tro_4_v_length1",&tro_4_v_length1);
  reader_tro_4.AddVariable("tro_4_v_medium_dQ_dx",&tro_4_v_medium_dQ_dx);
  reader_tro_4.AddVariable("tro_4_v_end_dQ_dx",&tro_4_v_end_dQ_dx);
  reader_tro_4.AddVariable("tro_4_v_energy",&tro_4_v_energy);
  reader_tro_4.AddVariable("tro_4_v_shower_main_length",&tro_4_v_shower_main_length);
  reader_tro_4.AddVariable("tro_4_v_flag_shower_trajectory",&tro_4_v_flag_shower_trajectory);
  reader_tro_4.BookMVA( "MyBDT", "weights/tro_4_BDT.weights.xml");

  TMVA::Reader reader_tro_5;
  float tro_5_v_max_angle;
  float tro_5_v_min_angle;
  float tro_5_v_max_length;
  float tro_5_v_iso_angle;
  float tro_5_v_n_vtx_segs;
  float tro_5_v_min_count;
  float tro_5_v_max_count;
  float tro_5_v_energy;
  reader_tro_5.AddVariable("tro_5_v_max_angle",&tro_5_v_max_angle);
  reader_tro_5.AddVariable("tro_5_v_min_angle",&tro_5_v_min_angle);
  reader_tro_5.AddVariable("tro_5_v_max_length",&tro_5_v_max_length);
  reader_tro_5.AddVariable("tro_5_v_iso_angle",&tro_5_v_iso_angle);
  reader_tro_5.AddVariable("tro_5_v_n_vtx_segs",&tro_5_v_n_vtx_segs);
  reader_tro_5.AddVariable("tro_5_v_min_count",&tro_5_v_min_count);
  reader_tro_5.AddVariable("tro_5_v_max_count",&tro_5_v_max_count);
  reader_tro_5.AddVariable("tro_5_v_energy",&tro_5_v_energy);
  reader_tro_5.BookMVA( "MyBDT", "weights/tro_5_BDT.weights.xml");

  TMVA::Reader reader_sig_1;
  float sig_1_v_angle;
  float sig_1_v_flag_single_shower;
  float sig_1_v_energy;
  float sig_1_v_energy_1;
  reader_sig_1.AddVariable("sig_1_v_angle",&sig_1_v_angle);
  reader_sig_1.AddVariable("sig_1_v_flag_single_shower",&sig_1_v_flag_single_shower);
  reader_sig_1.AddVariable("sig_1_v_energy",&sig_1_v_energy);
  reader_sig_1.AddVariable("sig_1_v_energy_1",&sig_1_v_energy_1);
  reader_sig_1.BookMVA( "MyBDT", "weights/sig_1_BDT.weights.xml");

  TMVA::Reader reader_sig_2;
  float sig_2_v_energy;
  float sig_2_v_shower_angle;
  float sig_2_v_flag_single_shower;
  float sig_2_v_medium_dQ_dx;
  float sig_2_v_start_dQ_dx;
  
  reader_sig_2.AddVariable("sig_2_v_energy",&sig_2_v_energy);
  reader_sig_2.AddVariable("sig_2_v_shower_angle",&sig_2_v_shower_angle);
  reader_sig_2.AddVariable("sig_2_v_flag_single_shower",&sig_2_v_flag_single_shower);
  reader_sig_2.AddVariable("sig_2_v_medium_dQ_dx",&sig_2_v_medium_dQ_dx);  
  reader_sig_2.AddVariable("sig_2_v_start_dQ_dx",&sig_2_v_start_dQ_dx);

  reader_sig_2.BookMVA( "MyBDT", "weights/sig_2_BDT.weights.xml");

  TMVA::Reader reader_stw_2;
  
  float stw_2_v_medium_dQ_dx;
  float stw_2_v_energy;
  float stw_2_v_angle;
  float stw_2_v_dir_length;
  float stw_2_v_max_dQ_dx;
  
  reader_stw_2.AddVariable("stw_2_v_medium_dQ_dx",&stw_2_v_medium_dQ_dx);
  reader_stw_2.AddVariable("stw_2_v_energy",&stw_2_v_energy);
  reader_stw_2.AddVariable("stw_2_v_angle",&stw_2_v_angle);
  reader_stw_2.AddVariable("stw_2_v_dir_length",&stw_2_v_dir_length);
  reader_stw_2.AddVariable("stw_2_v_max_dQ_dx",&stw_2_v_max_dQ_dx);
  
  reader_stw_2.BookMVA( "MyBDT", "weights/stw_2_BDT.weights.xml");

  TMVA::Reader reader_stw_3;

  float stw_3_v_angle;
  float stw_3_v_dir_length;
  float stw_3_v_energy;
  float stw_3_v_medium_dQ_dx;
  
  reader_stw_3.AddVariable("stw_3_v_angle",&stw_3_v_angle);
  reader_stw_3.AddVariable("stw_3_v_dir_length",&stw_3_v_dir_length);
  reader_stw_3.AddVariable("stw_3_v_energy",&stw_3_v_energy);
  reader_stw_3.AddVariable("stw_3_v_medium_dQ_dx",&stw_3_v_medium_dQ_dx);
  
  reader_stw_3.BookMVA( "MyBDT", "weights/stw_3_BDT.weights.xml");

  TMVA::Reader reader_stw_4;

  float stw_4_v_angle;
  float stw_4_v_dis;
  float stw_4_v_energy;
  
  reader_stw_4.AddVariable("stw_4_v_angle",&stw_4_v_angle);
  reader_stw_4.AddVariable("stw_4_v_dis",&stw_4_v_dis);
  reader_stw_4.AddVariable("stw_4_v_energy",&stw_4_v_energy);
  reader_stw_4.BookMVA( "MyBDT", "weights/stw_4_BDT.weights.xml");


  
  
  // signal, intrinsic nue ...
  for (Int_t i=0;i!=t1->GetEntries();i++){
    t1->GetEntry(i);

    // weight corrupted ...
    if (std::isnan(tagger.weight_cv) || std::isnan(tagger.weight_spline) || std::isinf(tagger.weight_cv) || std::isinf(tagger.weight_spline)) continue;
    // no EM showers ...
    if (tagger.br_filled!=1) continue;
    // odd sub-run number ...
    if (tagger.subrun %2 == 1 && process==1 || tagger.subrun%2 == 0 && process != 1) continue;

    
    
    tagger.weight = tagger.weight_spline * tagger.weight_cv * (1 + tagger.weight_lee);
    tagger.lowEweight = (1. + 5 * tagger.weight_lee)/(1. + tagger.weight_lee);

    // scale intrinsic nue higher ...
    if (tagger.truth_nuEnergy > 400) tagger.weight *= (pot_1+pot_2+pot_3+pot_4)/(pot_1+pot_2);
    
    if ( (tagger.nuvtx_diff < 1 && tagger.showervtx_diff < 1 &&process == 1 || process != 1) &&
	 tagger.truth_isCC == 1 && abs(tagger.truth_nuPdg)==12 && tagger.truth_vtxInside == 1){

      tagger.br3_3_score     = cal_br3_3_bdt(0.3, tagger,  reader_br3_3, br3_3_v_energy,  br3_3_v_angle,  br3_3_v_dir_length, br3_3_v_length);
      tagger.br3_5_score     = cal_br3_5_bdt(0.42, tagger,  reader_br3_5, br3_5_v_dir_length, br3_5_v_total_length, br3_5_v_flag_avoid_muon_check, br3_5_v_n_seg, br3_5_v_angle, br3_5_v_sg_length, br3_5_v_energy, br3_5_v_n_main_segs, br3_5_v_n_segs, br3_5_v_shower_main_length, br3_5_v_shower_total_length);
      tagger.br3_6_score     = cal_br3_6_bdt(0.75, tagger, reader_br3_6, br3_6_v_angle, br3_6_v_angle1, br3_6_v_flag_shower_trajectory, br3_6_v_direct_length, br3_6_v_length, br3_6_v_n_other_vtx_segs, br3_6_v_energy);
      tagger.pio_2_score     = cal_pio_2_bdt(0.2,  tagger,  reader_pio_2, pio_2_v_dis2, pio_2_v_angle2, pio_2_v_acc_length);
      tagger.stw_2_score     = cal_stw_2_bdt(0.7, tagger, reader_stw_2, stw_2_v_medium_dQ_dx, stw_2_v_energy, stw_2_v_angle, stw_2_v_dir_length, stw_2_v_max_dQ_dx);
      tagger.stw_3_score     = cal_stw_3_bdt(0.5, tagger, reader_stw_3, stw_3_v_angle, stw_3_v_dir_length, stw_3_v_energy, stw_3_v_medium_dQ_dx);
      tagger.stw_4_score     = cal_stw_4_bdt(0.7, tagger, reader_stw_4, stw_4_v_angle, stw_4_v_dis, stw_4_v_energy);
      tagger.sig_1_score     = cal_sig_1_bdt(0.59, tagger,  reader_sig_1, sig_1_v_angle, sig_1_v_flag_single_shower, sig_1_v_energy, sig_1_v_energy_1);
      tagger.sig_2_score     = cal_sig_2_bdt(0.55, tagger, reader_sig_2, sig_2_v_energy, sig_2_v_shower_angle, sig_2_v_flag_single_shower, sig_2_v_medium_dQ_dx, sig_2_v_start_dQ_dx);
      tagger.lol_1_score     = cal_lol_1_bdt(0.85, tagger, reader_lol_1, lol_1_v_energy, lol_1_v_vtx_n_segs, lol_1_v_nseg, lol_1_v_angle);
      tagger.lol_2_score     = cal_lol_2_bdt(0.7, tagger, reader_lol_2,  lol_2_v_length, lol_2_v_angle, lol_2_v_type, lol_2_v_vtx_n_segs, lol_2_v_energy, lol_2_v_shower_main_length, lol_2_v_flag_dir_weak);
      tagger.tro_1_score     = cal_tro_1_bdt(0.28, tagger, reader_tro_1, tro_1_v_particle_type, tro_1_v_flag_dir_weak, tro_1_v_min_dis, tro_1_v_sg1_length,
  tro_1_v_shower_main_length, tro_1_v_max_n_vtx_segs, tro_1_v_tmp_length, tro_1_v_medium_dQ_dx, tro_1_v_dQ_dx_cut, tro_1_v_flag_shower_topology);
      tagger.tro_2_score     = cal_tro_2_bdt(0.35, tagger, reader_tro_2, tro_2_v_energy, tro_2_v_stem_length, tro_2_v_iso_angle, tro_2_v_max_length, tro_2_v_angle);
      tagger.tro_4_score     = cal_tro_4_bdt(0.33, tagger, reader_tro_4, tro_4_v_dir2_mag, tro_4_v_angle, tro_4_v_angle1, tro_4_v_angle2, tro_4_v_length, tro_4_v_length1, tro_4_v_medium_dQ_dx, tro_4_v_end_dQ_dx, tro_4_v_energy, tro_4_v_shower_main_length, tro_4_v_flag_shower_trajectory);
      tagger.tro_5_score     = cal_tro_5_bdt(0.5, tagger, reader_tro_5, tro_5_v_max_angle, tro_5_v_min_angle, tro_5_v_max_length, tro_5_v_iso_angle, tro_5_v_n_vtx_segs, tro_5_v_min_count, tro_5_v_max_count, tro_5_v_energy);

      
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
    if (tagger.subrun %2 == 1 && process==1 || tagger.subrun%2 == 0 && process != 1) continue;
    //    if (tagger.subrun %2 == 1) continue;

    
    
    tagger.weight = tagger.weight_spline * tagger.weight_cv * (1 + tagger.weight_lee);
    tagger.lowEweight = (1. + 5 * tagger.weight_lee)/(1. + tagger.weight_lee);

    // scale intrinsic nue higher ...
    if (tagger.truth_nuEnergy > 400) tagger.weight *= (pot_1+pot_2+pot_3+pot_4)/(pot_1+pot_2);
    
    if ( (tagger.nuvtx_diff < 1 && tagger.showervtx_diff < 1 &&process == 1 || process != 1 ) && tagger.truth_isCC == 1 &&
	abs(tagger.truth_nuPdg)==12 && tagger.truth_vtxInside == 1){

      tagger.br3_3_score     = cal_br3_3_bdt(0.3, tagger,  reader_br3_3, br3_3_v_energy,  br3_3_v_angle,  br3_3_v_dir_length, br3_3_v_length);
      tagger.br3_5_score     = cal_br3_5_bdt(0.42, tagger,  reader_br3_5, br3_5_v_dir_length, br3_5_v_total_length, br3_5_v_flag_avoid_muon_check, br3_5_v_n_seg, br3_5_v_angle, br3_5_v_sg_length, br3_5_v_energy, br3_5_v_n_main_segs, br3_5_v_n_segs, br3_5_v_shower_main_length, br3_5_v_shower_total_length);
      tagger.br3_6_score     = cal_br3_6_bdt(0.75, tagger, reader_br3_6, br3_6_v_angle, br3_6_v_angle1, br3_6_v_flag_shower_trajectory, br3_6_v_direct_length, br3_6_v_length, br3_6_v_n_other_vtx_segs, br3_6_v_energy);
      tagger.pio_2_score     = cal_pio_2_bdt(0.2,  tagger,  reader_pio_2, pio_2_v_dis2, pio_2_v_angle2, pio_2_v_acc_length);
      tagger.stw_2_score     = cal_stw_2_bdt(0.7, tagger, reader_stw_2, stw_2_v_medium_dQ_dx, stw_2_v_energy, stw_2_v_angle, stw_2_v_dir_length, stw_2_v_max_dQ_dx);
      tagger.stw_3_score     = cal_stw_3_bdt(0.5, tagger, reader_stw_3, stw_3_v_angle, stw_3_v_dir_length, stw_3_v_energy, stw_3_v_medium_dQ_dx);
      tagger.stw_4_score     = cal_stw_4_bdt(0.7, tagger, reader_stw_4, stw_4_v_angle, stw_4_v_dis, stw_4_v_energy);
      tagger.sig_1_score     = cal_sig_1_bdt(0.59, tagger,  reader_sig_1, sig_1_v_angle, sig_1_v_flag_single_shower, sig_1_v_energy, sig_1_v_energy_1);
      tagger.sig_2_score     = cal_sig_2_bdt(0.55, tagger, reader_sig_2, sig_2_v_energy, sig_2_v_shower_angle, sig_2_v_flag_single_shower, sig_2_v_medium_dQ_dx, sig_2_v_start_dQ_dx);
      tagger.lol_1_score     = cal_lol_1_bdt(0.85, tagger, reader_lol_1, lol_1_v_energy, lol_1_v_vtx_n_segs, lol_1_v_nseg, lol_1_v_angle);
      tagger.lol_2_score     = cal_lol_2_bdt(0.7, tagger, reader_lol_2,  lol_2_v_length, lol_2_v_angle, lol_2_v_type, lol_2_v_vtx_n_segs, lol_2_v_energy, lol_2_v_shower_main_length, lol_2_v_flag_dir_weak);
      tagger.tro_1_score     = cal_tro_1_bdt(0.28, tagger, reader_tro_1, tro_1_v_particle_type, tro_1_v_flag_dir_weak, tro_1_v_min_dis, tro_1_v_sg1_length,
  tro_1_v_shower_main_length, tro_1_v_max_n_vtx_segs, tro_1_v_tmp_length, tro_1_v_medium_dQ_dx, tro_1_v_dQ_dx_cut, tro_1_v_flag_shower_topology);
      tagger.tro_2_score     = cal_tro_2_bdt(0.35, tagger, reader_tro_2, tro_2_v_energy, tro_2_v_stem_length, tro_2_v_iso_angle, tro_2_v_max_length, tro_2_v_angle);
      tagger.tro_4_score     = cal_tro_4_bdt(0.33, tagger, reader_tro_4, tro_4_v_dir2_mag, tro_4_v_angle, tro_4_v_angle1, tro_4_v_angle2, tro_4_v_length, tro_4_v_length1, tro_4_v_medium_dQ_dx, tro_4_v_end_dQ_dx, tro_4_v_energy, tro_4_v_shower_main_length, tro_4_v_flag_shower_trajectory);
      tagger.tro_5_score     = cal_tro_5_bdt(0.5, tagger, reader_tro_5, tro_5_v_max_angle, tro_5_v_min_angle, tro_5_v_max_length, tro_5_v_iso_angle, tro_5_v_n_vtx_segs, tro_5_v_min_count, tro_5_v_max_count, tro_5_v_energy);
      
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
    if (tagger.subrun %2 == 1 && process==1 || tagger.subrun%2 == 0 && process != 1) continue;
    //    if (tagger.subrun %2 == 1) continue;

    
    
    tagger.weight = tagger.weight_spline * tagger.weight_cv * (1 + tagger.weight_lee);
    tagger.lowEweight = (1. + 5 * tagger.weight_lee)/(1. + tagger.weight_lee);

            
    if ( (tagger.nuvtx_diff < 1 && tagger.showervtx_diff < 1 &&process == 1 || process != 1) && tagger.truth_isCC == 1 &&
	abs(tagger.truth_nuPdg)==12 && tagger.truth_vtxInside == 1){

      tagger.br3_3_score     = cal_br3_3_bdt(0.3, tagger,  reader_br3_3, br3_3_v_energy,  br3_3_v_angle,  br3_3_v_dir_length, br3_3_v_length);
      tagger.br3_5_score     = cal_br3_5_bdt(0.42, tagger,  reader_br3_5, br3_5_v_dir_length, br3_5_v_total_length, br3_5_v_flag_avoid_muon_check, br3_5_v_n_seg, br3_5_v_angle, br3_5_v_sg_length, br3_5_v_energy, br3_5_v_n_main_segs, br3_5_v_n_segs, br3_5_v_shower_main_length, br3_5_v_shower_total_length);
      tagger.br3_6_score     = cal_br3_6_bdt(0.75, tagger, reader_br3_6, br3_6_v_angle, br3_6_v_angle1, br3_6_v_flag_shower_trajectory, br3_6_v_direct_length, br3_6_v_length, br3_6_v_n_other_vtx_segs, br3_6_v_energy);
      tagger.pio_2_score     = cal_pio_2_bdt(0.2,  tagger,  reader_pio_2, pio_2_v_dis2, pio_2_v_angle2, pio_2_v_acc_length);
      tagger.stw_2_score     = cal_stw_2_bdt(0.7, tagger, reader_stw_2, stw_2_v_medium_dQ_dx, stw_2_v_energy, stw_2_v_angle, stw_2_v_dir_length, stw_2_v_max_dQ_dx);
      tagger.stw_3_score     = cal_stw_3_bdt(0.5, tagger, reader_stw_3, stw_3_v_angle, stw_3_v_dir_length, stw_3_v_energy, stw_3_v_medium_dQ_dx);
      tagger.stw_4_score     = cal_stw_4_bdt(0.7, tagger, reader_stw_4, stw_4_v_angle, stw_4_v_dis, stw_4_v_energy);
      tagger.sig_1_score     = cal_sig_1_bdt(0.59, tagger,  reader_sig_1, sig_1_v_angle, sig_1_v_flag_single_shower, sig_1_v_energy, sig_1_v_energy_1);
      tagger.sig_2_score     = cal_sig_2_bdt(0.55, tagger, reader_sig_2, sig_2_v_energy, sig_2_v_shower_angle, sig_2_v_flag_single_shower, sig_2_v_medium_dQ_dx, sig_2_v_start_dQ_dx);
      tagger.lol_1_score     = cal_lol_1_bdt(0.85, tagger, reader_lol_1, lol_1_v_energy, lol_1_v_vtx_n_segs, lol_1_v_nseg, lol_1_v_angle);
      tagger.lol_2_score     = cal_lol_2_bdt(0.7, tagger, reader_lol_2,  lol_2_v_length, lol_2_v_angle, lol_2_v_type, lol_2_v_vtx_n_segs, lol_2_v_energy, lol_2_v_shower_main_length, lol_2_v_flag_dir_weak);
      tagger.tro_1_score     = cal_tro_1_bdt(0.28, tagger, reader_tro_1, tro_1_v_particle_type, tro_1_v_flag_dir_weak, tro_1_v_min_dis, tro_1_v_sg1_length,
  tro_1_v_shower_main_length, tro_1_v_max_n_vtx_segs, tro_1_v_tmp_length, tro_1_v_medium_dQ_dx, tro_1_v_dQ_dx_cut, tro_1_v_flag_shower_topology);
      tagger.tro_2_score     = cal_tro_2_bdt(0.35, tagger, reader_tro_2, tro_2_v_energy, tro_2_v_stem_length, tro_2_v_iso_angle, tro_2_v_max_length, tro_2_v_angle);
      tagger.tro_4_score     = cal_tro_4_bdt(0.33, tagger, reader_tro_4, tro_4_v_dir2_mag, tro_4_v_angle, tro_4_v_angle1, tro_4_v_angle2, tro_4_v_length, tro_4_v_length1, tro_4_v_medium_dQ_dx, tro_4_v_end_dQ_dx, tro_4_v_energy, tro_4_v_shower_main_length, tro_4_v_flag_shower_trajectory);
      tagger.tro_5_score     = cal_tro_5_bdt(0.5, tagger, reader_tro_5, tro_5_v_max_angle, tro_5_v_min_angle, tro_5_v_max_length, tro_5_v_iso_angle, tro_5_v_n_vtx_segs, tro_5_v_min_count, tro_5_v_max_count, tro_5_v_energy);
      
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
    if (tagger.subrun %2 == 1 && process==1 || tagger.subrun%2 == 0 && process != 1) continue;
    //    if (tagger.subrun %2 == 1) continue;

    
    
    tagger.weight = tagger.weight_spline * tagger.weight_cv * (1 + tagger.weight_lee);
    tagger.lowEweight = (1. + 5 * tagger.weight_lee)/(1. + tagger.weight_lee);
   
    
    if ( (tagger.nuvtx_diff < 1 && tagger.showervtx_diff < 1 &&process == 1 || process != 1) && tagger.truth_isCC == 1 &&
	abs(tagger.truth_nuPdg)==12 && tagger.truth_vtxInside == 1){

      tagger.br3_3_score     = cal_br3_3_bdt(0.3, tagger,  reader_br3_3, br3_3_v_energy,  br3_3_v_angle,  br3_3_v_dir_length, br3_3_v_length);
      tagger.br3_5_score     = cal_br3_5_bdt(0.42, tagger,  reader_br3_5, br3_5_v_dir_length, br3_5_v_total_length, br3_5_v_flag_avoid_muon_check, br3_5_v_n_seg, br3_5_v_angle, br3_5_v_sg_length, br3_5_v_energy, br3_5_v_n_main_segs, br3_5_v_n_segs, br3_5_v_shower_main_length, br3_5_v_shower_total_length);
      tagger.br3_6_score     = cal_br3_6_bdt(0.75, tagger, reader_br3_6, br3_6_v_angle, br3_6_v_angle1, br3_6_v_flag_shower_trajectory, br3_6_v_direct_length, br3_6_v_length, br3_6_v_n_other_vtx_segs, br3_6_v_energy);
      tagger.pio_2_score     = cal_pio_2_bdt(0.2,  tagger,  reader_pio_2, pio_2_v_dis2, pio_2_v_angle2, pio_2_v_acc_length);
      tagger.stw_2_score     = cal_stw_2_bdt(0.7, tagger, reader_stw_2, stw_2_v_medium_dQ_dx, stw_2_v_energy, stw_2_v_angle, stw_2_v_dir_length, stw_2_v_max_dQ_dx);
      tagger.stw_3_score     = cal_stw_3_bdt(0.5, tagger, reader_stw_3, stw_3_v_angle, stw_3_v_dir_length, stw_3_v_energy, stw_3_v_medium_dQ_dx);
      tagger.stw_4_score     = cal_stw_4_bdt(0.7, tagger, reader_stw_4, stw_4_v_angle, stw_4_v_dis, stw_4_v_energy);
      tagger.sig_1_score     = cal_sig_1_bdt(0.59, tagger,  reader_sig_1, sig_1_v_angle, sig_1_v_flag_single_shower, sig_1_v_energy, sig_1_v_energy_1);
      tagger.sig_2_score     = cal_sig_2_bdt(0.55, tagger, reader_sig_2, sig_2_v_energy, sig_2_v_shower_angle, sig_2_v_flag_single_shower, sig_2_v_medium_dQ_dx, sig_2_v_start_dQ_dx);
      tagger.lol_1_score     = cal_lol_1_bdt(0.85, tagger, reader_lol_1, lol_1_v_energy, lol_1_v_vtx_n_segs, lol_1_v_nseg, lol_1_v_angle);
      tagger.lol_2_score     = cal_lol_2_bdt(0.7, tagger, reader_lol_2,  lol_2_v_length, lol_2_v_angle, lol_2_v_type, lol_2_v_vtx_n_segs, lol_2_v_energy, lol_2_v_shower_main_length, lol_2_v_flag_dir_weak);
      tagger.tro_1_score     = cal_tro_1_bdt(0.28, tagger, reader_tro_1, tro_1_v_particle_type, tro_1_v_flag_dir_weak, tro_1_v_min_dis, tro_1_v_sg1_length,
  tro_1_v_shower_main_length, tro_1_v_max_n_vtx_segs, tro_1_v_tmp_length, tro_1_v_medium_dQ_dx, tro_1_v_dQ_dx_cut, tro_1_v_flag_shower_topology);
      tagger.tro_2_score     = cal_tro_2_bdt(0.35, tagger, reader_tro_2, tro_2_v_energy, tro_2_v_stem_length, tro_2_v_iso_angle, tro_2_v_max_length, tro_2_v_angle);
      tagger.tro_4_score     = cal_tro_4_bdt(0.33, tagger, reader_tro_4, tro_4_v_dir2_mag, tro_4_v_angle, tro_4_v_angle1, tro_4_v_angle2, tro_4_v_length, tro_4_v_length1, tro_4_v_medium_dQ_dx, tro_4_v_end_dQ_dx, tro_4_v_energy, tro_4_v_shower_main_length, tro_4_v_flag_shower_trajectory);
      tagger.tro_5_score     = cal_tro_5_bdt(0.5, tagger, reader_tro_5, tro_5_v_max_angle, tro_5_v_min_angle, tro_5_v_max_length, tro_5_v_iso_angle, tro_5_v_n_vtx_segs, tro_5_v_min_count, tro_5_v_max_count, tro_5_v_energy);
      
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
    if (tagger.subrun %2 == 1 && process==1 || tagger.subrun%2 == 0 && process != 1) continue;
    //    if (tagger.subrun %2 == 1) continue;
    
    // also need to exclude the NC nu-electron elastic scattering 
    if ((tagger.truth_isCC==1 && abs(tagger.truth_nuPdg)==12 || tagger.truth_nuIntType==1098) && process == 1 ) continue;

    if (tagger.truth_nuEnergy > 400)
      tagger.weight = tagger.weight_spline * tagger.weight_cv * (pot_1+pot_2+pot_3+pot_4)/(pot_5+pot_6 );
    else
      tagger.weight = tagger.weight_spline * tagger.weight_cv * (pot_1+pot_2+pot_3+pot_4)/(pot_5+pot_6 +pot_9 + pot_10);
    
    tagger.lowEweight = 1;


    tagger.br3_3_score     = cal_br3_3_bdt(0.3, tagger,  reader_br3_3, br3_3_v_energy,  br3_3_v_angle,  br3_3_v_dir_length, br3_3_v_length);
      tagger.br3_5_score     = cal_br3_5_bdt(0.42, tagger,  reader_br3_5, br3_5_v_dir_length, br3_5_v_total_length, br3_5_v_flag_avoid_muon_check, br3_5_v_n_seg, br3_5_v_angle, br3_5_v_sg_length, br3_5_v_energy, br3_5_v_n_main_segs, br3_5_v_n_segs, br3_5_v_shower_main_length, br3_5_v_shower_total_length);
      tagger.br3_6_score     = cal_br3_6_bdt(0.75, tagger, reader_br3_6, br3_6_v_angle, br3_6_v_angle1, br3_6_v_flag_shower_trajectory, br3_6_v_direct_length, br3_6_v_length, br3_6_v_n_other_vtx_segs, br3_6_v_energy);
      tagger.pio_2_score     = cal_pio_2_bdt(0.2,  tagger,  reader_pio_2, pio_2_v_dis2, pio_2_v_angle2, pio_2_v_acc_length);
      tagger.stw_2_score     = cal_stw_2_bdt(0.7, tagger, reader_stw_2, stw_2_v_medium_dQ_dx, stw_2_v_energy, stw_2_v_angle, stw_2_v_dir_length, stw_2_v_max_dQ_dx);
      tagger.stw_3_score     = cal_stw_3_bdt(0.5, tagger, reader_stw_3, stw_3_v_angle, stw_3_v_dir_length, stw_3_v_energy, stw_3_v_medium_dQ_dx);
      tagger.stw_4_score     = cal_stw_4_bdt(0.7, tagger, reader_stw_4, stw_4_v_angle, stw_4_v_dis, stw_4_v_energy);
      tagger.sig_1_score     = cal_sig_1_bdt(0.59, tagger,  reader_sig_1, sig_1_v_angle, sig_1_v_flag_single_shower, sig_1_v_energy, sig_1_v_energy_1);
      tagger.sig_2_score     = cal_sig_2_bdt(0.55, tagger, reader_sig_2, sig_2_v_energy, sig_2_v_shower_angle, sig_2_v_flag_single_shower, sig_2_v_medium_dQ_dx, sig_2_v_start_dQ_dx);
      tagger.lol_1_score     = cal_lol_1_bdt(0.85, tagger, reader_lol_1, lol_1_v_energy, lol_1_v_vtx_n_segs, lol_1_v_nseg, lol_1_v_angle);
      tagger.lol_2_score     = cal_lol_2_bdt(0.7, tagger, reader_lol_2,  lol_2_v_length, lol_2_v_angle, lol_2_v_type, lol_2_v_vtx_n_segs, lol_2_v_energy, lol_2_v_shower_main_length, lol_2_v_flag_dir_weak);
      tagger.tro_1_score     = cal_tro_1_bdt(0.28, tagger, reader_tro_1, tro_1_v_particle_type, tro_1_v_flag_dir_weak, tro_1_v_min_dis, tro_1_v_sg1_length,
  tro_1_v_shower_main_length, tro_1_v_max_n_vtx_segs, tro_1_v_tmp_length, tro_1_v_medium_dQ_dx, tro_1_v_dQ_dx_cut, tro_1_v_flag_shower_topology);
      tagger.tro_2_score     = cal_tro_2_bdt(0.35, tagger, reader_tro_2, tro_2_v_energy, tro_2_v_stem_length, tro_2_v_iso_angle, tro_2_v_max_length, tro_2_v_angle);
      tagger.tro_4_score     = cal_tro_4_bdt(0.33, tagger, reader_tro_4, tro_4_v_dir2_mag, tro_4_v_angle, tro_4_v_angle1, tro_4_v_angle2, tro_4_v_length, tro_4_v_length1, tro_4_v_medium_dQ_dx, tro_4_v_end_dQ_dx, tro_4_v_energy, tro_4_v_shower_main_length, tro_4_v_flag_shower_trajectory);
      tagger.tro_5_score     = cal_tro_5_bdt(0.5, tagger, reader_tro_5, tro_5_v_max_angle, tro_5_v_min_angle, tro_5_v_max_length, tro_5_v_iso_angle, tro_5_v_n_vtx_segs, tro_5_v_min_count, tro_5_v_max_count, tro_5_v_energy);


    
    Tbkg->Fill();
  }

  for (Int_t i=0;i!=t6->GetEntries();i++){
    t6->GetEntry(i);

    if (std::isnan(tagger.weight_cv) || std::isnan(tagger.weight_spline) || std::isinf(tagger.weight_cv) || std::isinf(tagger.weight_spline)) continue;
    // require EM shower ...
    if (tagger.br_filled!=1) continue;
    // only use even subrun number to train ...
    if (tagger.subrun %2 == 1 && process==1 || tagger.subrun%2 == 0 && process != 1) continue;
    //    if (tagger.subrun %2 == 1) continue;
    
    // also need to exclude the NC nu-electron elastic scattering 
    if ((tagger.truth_isCC==1 && abs(tagger.truth_nuPdg)==12 || tagger.truth_nuIntType==1098)&&process==1 ) continue;

    if (tagger.truth_nuEnergy > 400)
      tagger.weight = tagger.weight_spline * tagger.weight_cv * (pot_1+pot_2+pot_3+pot_4)/(pot_5+pot_6);
    else
      tagger.weight = tagger.weight_spline * tagger.weight_cv * (pot_1+pot_2+pot_3+pot_4)/(pot_5+pot_6+pot_9+pot_10);
    
    tagger.lowEweight = 1;

    tagger.br3_3_score     = cal_br3_3_bdt(0.3, tagger,  reader_br3_3, br3_3_v_energy,  br3_3_v_angle,  br3_3_v_dir_length, br3_3_v_length);
      tagger.br3_5_score     = cal_br3_5_bdt(0.42, tagger,  reader_br3_5, br3_5_v_dir_length, br3_5_v_total_length, br3_5_v_flag_avoid_muon_check, br3_5_v_n_seg, br3_5_v_angle, br3_5_v_sg_length, br3_5_v_energy, br3_5_v_n_main_segs, br3_5_v_n_segs, br3_5_v_shower_main_length, br3_5_v_shower_total_length);
      tagger.br3_6_score     = cal_br3_6_bdt(0.75, tagger, reader_br3_6, br3_6_v_angle, br3_6_v_angle1, br3_6_v_flag_shower_trajectory, br3_6_v_direct_length, br3_6_v_length, br3_6_v_n_other_vtx_segs, br3_6_v_energy);
      tagger.pio_2_score     = cal_pio_2_bdt(0.2,  tagger,  reader_pio_2, pio_2_v_dis2, pio_2_v_angle2, pio_2_v_acc_length);
      tagger.stw_2_score     = cal_stw_2_bdt(0.7, tagger, reader_stw_2, stw_2_v_medium_dQ_dx, stw_2_v_energy, stw_2_v_angle, stw_2_v_dir_length, stw_2_v_max_dQ_dx);
      tagger.stw_3_score     = cal_stw_3_bdt(0.5, tagger, reader_stw_3, stw_3_v_angle, stw_3_v_dir_length, stw_3_v_energy, stw_3_v_medium_dQ_dx);
      tagger.stw_4_score     = cal_stw_4_bdt(0.7, tagger, reader_stw_4, stw_4_v_angle, stw_4_v_dis, stw_4_v_energy);
      tagger.sig_1_score     = cal_sig_1_bdt(0.59, tagger,  reader_sig_1, sig_1_v_angle, sig_1_v_flag_single_shower, sig_1_v_energy, sig_1_v_energy_1);
      tagger.sig_2_score     = cal_sig_2_bdt(0.55, tagger, reader_sig_2, sig_2_v_energy, sig_2_v_shower_angle, sig_2_v_flag_single_shower, sig_2_v_medium_dQ_dx, sig_2_v_start_dQ_dx);
      tagger.lol_1_score     = cal_lol_1_bdt(0.85, tagger, reader_lol_1, lol_1_v_energy, lol_1_v_vtx_n_segs, lol_1_v_nseg, lol_1_v_angle);
      tagger.lol_2_score     = cal_lol_2_bdt(0.7, tagger, reader_lol_2,  lol_2_v_length, lol_2_v_angle, lol_2_v_type, lol_2_v_vtx_n_segs, lol_2_v_energy, lol_2_v_shower_main_length, lol_2_v_flag_dir_weak);
      tagger.tro_1_score     = cal_tro_1_bdt(0.28, tagger, reader_tro_1, tro_1_v_particle_type, tro_1_v_flag_dir_weak, tro_1_v_min_dis, tro_1_v_sg1_length,
  tro_1_v_shower_main_length, tro_1_v_max_n_vtx_segs, tro_1_v_tmp_length, tro_1_v_medium_dQ_dx, tro_1_v_dQ_dx_cut, tro_1_v_flag_shower_topology);
      tagger.tro_2_score     = cal_tro_2_bdt(0.35, tagger, reader_tro_2, tro_2_v_energy, tro_2_v_stem_length, tro_2_v_iso_angle, tro_2_v_max_length, tro_2_v_angle);
      tagger.tro_4_score     = cal_tro_4_bdt(0.33, tagger, reader_tro_4, tro_4_v_dir2_mag, tro_4_v_angle, tro_4_v_angle1, tro_4_v_angle2, tro_4_v_length, tro_4_v_length1, tro_4_v_medium_dQ_dx, tro_4_v_end_dQ_dx, tro_4_v_energy, tro_4_v_shower_main_length, tro_4_v_flag_shower_trajectory);
      tagger.tro_5_score     = cal_tro_5_bdt(0.5, tagger, reader_tro_5, tro_5_v_max_angle, tro_5_v_min_angle, tro_5_v_max_length, tro_5_v_iso_angle, tro_5_v_n_vtx_segs, tro_5_v_min_count, tro_5_v_max_count, tro_5_v_energy);

    
    Tbkg->Fill();
  }

  // low energy patck nu overlay ...
  for (Int_t i=0;i!=t9->GetEntries();i++){
    t9->GetEntry(i);

    if (std::isnan(tagger.weight_cv) || std::isnan(tagger.weight_spline) || std::isinf(tagger.weight_cv) || std::isinf(tagger.weight_spline)) continue;
    // require EM shower ...
    if (tagger.br_filled!=1) continue;
    // only use even subrun number to train ...
    if (tagger.subrun %2 == 1 && process==1 || tagger.subrun%2 == 0 && process != 1) continue;
    //    if (tagger.subrun %2 == 1) continue;
    
    // also need to exclude the NC nu-electron elastic scattering 
    if ((tagger.truth_isCC==1 && abs(tagger.truth_nuPdg)==12 || tagger.truth_nuIntType==1098) && process==1 ) continue;

    tagger.weight = tagger.weight_spline * tagger.weight_cv * (pot_1+pot_2+pot_3+pot_4)/(pot_5+pot_6 +pot_9 + pot_10);
    
    tagger.lowEweight = 1;

    tagger.br3_3_score     = cal_br3_3_bdt(0.3, tagger,  reader_br3_3, br3_3_v_energy,  br3_3_v_angle,  br3_3_v_dir_length, br3_3_v_length);
      tagger.br3_5_score     = cal_br3_5_bdt(0.42, tagger,  reader_br3_5, br3_5_v_dir_length, br3_5_v_total_length, br3_5_v_flag_avoid_muon_check, br3_5_v_n_seg, br3_5_v_angle, br3_5_v_sg_length, br3_5_v_energy, br3_5_v_n_main_segs, br3_5_v_n_segs, br3_5_v_shower_main_length, br3_5_v_shower_total_length);
      tagger.br3_6_score     = cal_br3_6_bdt(0.75, tagger, reader_br3_6, br3_6_v_angle, br3_6_v_angle1, br3_6_v_flag_shower_trajectory, br3_6_v_direct_length, br3_6_v_length, br3_6_v_n_other_vtx_segs, br3_6_v_energy);
      tagger.pio_2_score     = cal_pio_2_bdt(0.2,  tagger,  reader_pio_2, pio_2_v_dis2, pio_2_v_angle2, pio_2_v_acc_length);
      tagger.stw_2_score     = cal_stw_2_bdt(0.7, tagger, reader_stw_2, stw_2_v_medium_dQ_dx, stw_2_v_energy, stw_2_v_angle, stw_2_v_dir_length, stw_2_v_max_dQ_dx);
      tagger.stw_3_score     = cal_stw_3_bdt(0.5, tagger, reader_stw_3, stw_3_v_angle, stw_3_v_dir_length, stw_3_v_energy, stw_3_v_medium_dQ_dx);
      tagger.stw_4_score     = cal_stw_4_bdt(0.7, tagger, reader_stw_4, stw_4_v_angle, stw_4_v_dis, stw_4_v_energy);
      tagger.sig_1_score     = cal_sig_1_bdt(0.59, tagger,  reader_sig_1, sig_1_v_angle, sig_1_v_flag_single_shower, sig_1_v_energy, sig_1_v_energy_1);
      tagger.sig_2_score     = cal_sig_2_bdt(0.55, tagger, reader_sig_2, sig_2_v_energy, sig_2_v_shower_angle, sig_2_v_flag_single_shower, sig_2_v_medium_dQ_dx, sig_2_v_start_dQ_dx);
      tagger.lol_1_score     = cal_lol_1_bdt(0.85, tagger, reader_lol_1, lol_1_v_energy, lol_1_v_vtx_n_segs, lol_1_v_nseg, lol_1_v_angle);
      tagger.lol_2_score     = cal_lol_2_bdt(0.7, tagger, reader_lol_2,  lol_2_v_length, lol_2_v_angle, lol_2_v_type, lol_2_v_vtx_n_segs, lol_2_v_energy, lol_2_v_shower_main_length, lol_2_v_flag_dir_weak);
      tagger.tro_1_score     = cal_tro_1_bdt(0.28, tagger, reader_tro_1, tro_1_v_particle_type, tro_1_v_flag_dir_weak, tro_1_v_min_dis, tro_1_v_sg1_length,
  tro_1_v_shower_main_length, tro_1_v_max_n_vtx_segs, tro_1_v_tmp_length, tro_1_v_medium_dQ_dx, tro_1_v_dQ_dx_cut, tro_1_v_flag_shower_topology);
      tagger.tro_2_score     = cal_tro_2_bdt(0.35, tagger, reader_tro_2, tro_2_v_energy, tro_2_v_stem_length, tro_2_v_iso_angle, tro_2_v_max_length, tro_2_v_angle);
      tagger.tro_4_score     = cal_tro_4_bdt(0.33, tagger, reader_tro_4, tro_4_v_dir2_mag, tro_4_v_angle, tro_4_v_angle1, tro_4_v_angle2, tro_4_v_length, tro_4_v_length1, tro_4_v_medium_dQ_dx, tro_4_v_end_dQ_dx, tro_4_v_energy, tro_4_v_shower_main_length, tro_4_v_flag_shower_trajectory);
      tagger.tro_5_score     = cal_tro_5_bdt(0.5, tagger, reader_tro_5, tro_5_v_max_angle, tro_5_v_min_angle, tro_5_v_max_length, tro_5_v_iso_angle, tro_5_v_n_vtx_segs, tro_5_v_min_count, tro_5_v_max_count, tro_5_v_energy);

    
    Tbkg->Fill();
  }
  

  for (Int_t i=0;i!=t10->GetEntries();i++){
    t10->GetEntry(i);

    if (std::isnan(tagger.weight_cv) || std::isnan(tagger.weight_spline) || std::isinf(tagger.weight_cv) || std::isinf(tagger.weight_spline)) continue;
    // require EM shower ...
    if (tagger.br_filled!=1) continue;
    // only use even subrun number to train ...
    if (tagger.subrun %2 == 1 && process==1 || tagger.subrun%2 == 0 && process != 1) continue;
    //    if (tagger.subrun %2 == 1) continue;
    
    // also need to exclude the NC nu-electron elastic scattering 
    if ((tagger.truth_isCC==1 && abs(tagger.truth_nuPdg)==12 || tagger.truth_nuIntType==1098) && process==1 ) continue;

    tagger.weight = tagger.weight_spline * tagger.weight_cv * (pot_1+pot_2+pot_3+pot_4)/(pot_5+pot_6 +pot_9 + pot_10);
    
    tagger.lowEweight = 1;

    tagger.br3_3_score     = cal_br3_3_bdt(0.3, tagger,  reader_br3_3, br3_3_v_energy,  br3_3_v_angle,  br3_3_v_dir_length, br3_3_v_length);
      tagger.br3_5_score     = cal_br3_5_bdt(0.42, tagger,  reader_br3_5, br3_5_v_dir_length, br3_5_v_total_length, br3_5_v_flag_avoid_muon_check, br3_5_v_n_seg, br3_5_v_angle, br3_5_v_sg_length, br3_5_v_energy, br3_5_v_n_main_segs, br3_5_v_n_segs, br3_5_v_shower_main_length, br3_5_v_shower_total_length);
      tagger.br3_6_score     = cal_br3_6_bdt(0.75, tagger, reader_br3_6, br3_6_v_angle, br3_6_v_angle1, br3_6_v_flag_shower_trajectory, br3_6_v_direct_length, br3_6_v_length, br3_6_v_n_other_vtx_segs, br3_6_v_energy);
      tagger.pio_2_score     = cal_pio_2_bdt(0.2,  tagger,  reader_pio_2, pio_2_v_dis2, pio_2_v_angle2, pio_2_v_acc_length);
      tagger.stw_2_score     = cal_stw_2_bdt(0.7, tagger, reader_stw_2, stw_2_v_medium_dQ_dx, stw_2_v_energy, stw_2_v_angle, stw_2_v_dir_length, stw_2_v_max_dQ_dx);
      tagger.stw_3_score     = cal_stw_3_bdt(0.5, tagger, reader_stw_3, stw_3_v_angle, stw_3_v_dir_length, stw_3_v_energy, stw_3_v_medium_dQ_dx);
      tagger.stw_4_score     = cal_stw_4_bdt(0.7, tagger, reader_stw_4, stw_4_v_angle, stw_4_v_dis, stw_4_v_energy);
      tagger.sig_1_score     = cal_sig_1_bdt(0.59, tagger,  reader_sig_1, sig_1_v_angle, sig_1_v_flag_single_shower, sig_1_v_energy, sig_1_v_energy_1);
      tagger.sig_2_score     = cal_sig_2_bdt(0.55, tagger, reader_sig_2, sig_2_v_energy, sig_2_v_shower_angle, sig_2_v_flag_single_shower, sig_2_v_medium_dQ_dx, sig_2_v_start_dQ_dx);
      tagger.lol_1_score     = cal_lol_1_bdt(0.85, tagger, reader_lol_1, lol_1_v_energy, lol_1_v_vtx_n_segs, lol_1_v_nseg, lol_1_v_angle);
      tagger.lol_2_score     = cal_lol_2_bdt(0.7, tagger, reader_lol_2,  lol_2_v_length, lol_2_v_angle, lol_2_v_type, lol_2_v_vtx_n_segs, lol_2_v_energy, lol_2_v_shower_main_length, lol_2_v_flag_dir_weak);
      tagger.tro_1_score     = cal_tro_1_bdt(0.28, tagger, reader_tro_1, tro_1_v_particle_type, tro_1_v_flag_dir_weak, tro_1_v_min_dis, tro_1_v_sg1_length,
  tro_1_v_shower_main_length, tro_1_v_max_n_vtx_segs, tro_1_v_tmp_length, tro_1_v_medium_dQ_dx, tro_1_v_dQ_dx_cut, tro_1_v_flag_shower_topology);
      tagger.tro_2_score     = cal_tro_2_bdt(0.35, tagger, reader_tro_2, tro_2_v_energy, tro_2_v_stem_length, tro_2_v_iso_angle, tro_2_v_max_length, tro_2_v_angle);
      tagger.tro_4_score     = cal_tro_4_bdt(0.33, tagger, reader_tro_4, tro_4_v_dir2_mag, tro_4_v_angle, tro_4_v_angle1, tro_4_v_angle2, tro_4_v_length, tro_4_v_length1, tro_4_v_medium_dQ_dx, tro_4_v_end_dQ_dx, tro_4_v_energy, tro_4_v_shower_main_length, tro_4_v_flag_shower_trajectory);
      tagger.tro_5_score     = cal_tro_5_bdt(0.5, tagger, reader_tro_5, tro_5_v_max_angle, tro_5_v_min_angle, tro_5_v_max_length, tro_5_v_iso_angle, tro_5_v_n_vtx_segs, tro_5_v_min_count, tro_5_v_max_count, tro_5_v_energy);

    
    Tbkg->Fill();
  }
  
  


  // EXT background 
  for (Int_t i=0;i!=t7->GetEntries();i++){
    t7->GetEntry(i);

    //    if (std::isnan(tagger.weight_cv) || std::isnan(tagger.weight_spline) || std::isinf(tagger.weight_cv) || std::isinf(tagger.weight_spline)) continue;
    // require EM shower ...
    if (tagger.br_filled!=1) continue;
    // only use even subrun number to train ...
    if (tagger.subrun %2 == 1 && process==1 || tagger.subrun%2 == 0 && process != 1) continue;
    //    if (tagger.subrun %2 == 1) continue;
    
    // also need to exclude the NC nu-electron elastic scattering 
    //if (tagger.truth_isCC==1 && abs(tagger.truth_nuPdg)==12 || tagger.truth_nuIntType==1098 ) continue;
    
    tagger.weight = (pot_1+pot_2+pot_3+pot_4)/(pot_7+pot_8);
    tagger.lowEweight = 1;

    tagger.br3_3_score     = cal_br3_3_bdt(0.3, tagger,  reader_br3_3, br3_3_v_energy,  br3_3_v_angle,  br3_3_v_dir_length, br3_3_v_length);
      tagger.br3_5_score     = cal_br3_5_bdt(0.42, tagger,  reader_br3_5, br3_5_v_dir_length, br3_5_v_total_length, br3_5_v_flag_avoid_muon_check, br3_5_v_n_seg, br3_5_v_angle, br3_5_v_sg_length, br3_5_v_energy, br3_5_v_n_main_segs, br3_5_v_n_segs, br3_5_v_shower_main_length, br3_5_v_shower_total_length);
      tagger.br3_6_score     = cal_br3_6_bdt(0.75, tagger, reader_br3_6, br3_6_v_angle, br3_6_v_angle1, br3_6_v_flag_shower_trajectory, br3_6_v_direct_length, br3_6_v_length, br3_6_v_n_other_vtx_segs, br3_6_v_energy);
      tagger.pio_2_score     = cal_pio_2_bdt(0.2,  tagger,  reader_pio_2, pio_2_v_dis2, pio_2_v_angle2, pio_2_v_acc_length);
      tagger.stw_2_score     = cal_stw_2_bdt(0.7, tagger, reader_stw_2, stw_2_v_medium_dQ_dx, stw_2_v_energy, stw_2_v_angle, stw_2_v_dir_length, stw_2_v_max_dQ_dx);
      tagger.stw_3_score     = cal_stw_3_bdt(0.5, tagger, reader_stw_3, stw_3_v_angle, stw_3_v_dir_length, stw_3_v_energy, stw_3_v_medium_dQ_dx);
      tagger.stw_4_score     = cal_stw_4_bdt(0.7, tagger, reader_stw_4, stw_4_v_angle, stw_4_v_dis, stw_4_v_energy);
      tagger.sig_1_score     = cal_sig_1_bdt(0.59, tagger,  reader_sig_1, sig_1_v_angle, sig_1_v_flag_single_shower, sig_1_v_energy, sig_1_v_energy_1);
      tagger.sig_2_score     = cal_sig_2_bdt(0.55, tagger, reader_sig_2, sig_2_v_energy, sig_2_v_shower_angle, sig_2_v_flag_single_shower, sig_2_v_medium_dQ_dx, sig_2_v_start_dQ_dx);
      tagger.lol_1_score     = cal_lol_1_bdt(0.85, tagger, reader_lol_1, lol_1_v_energy, lol_1_v_vtx_n_segs, lol_1_v_nseg, lol_1_v_angle);
      tagger.lol_2_score     = cal_lol_2_bdt(0.7, tagger, reader_lol_2,  lol_2_v_length, lol_2_v_angle, lol_2_v_type, lol_2_v_vtx_n_segs, lol_2_v_energy, lol_2_v_shower_main_length, lol_2_v_flag_dir_weak);
      tagger.tro_1_score     = cal_tro_1_bdt(0.28, tagger, reader_tro_1, tro_1_v_particle_type, tro_1_v_flag_dir_weak, tro_1_v_min_dis, tro_1_v_sg1_length,
  tro_1_v_shower_main_length, tro_1_v_max_n_vtx_segs, tro_1_v_tmp_length, tro_1_v_medium_dQ_dx, tro_1_v_dQ_dx_cut, tro_1_v_flag_shower_topology);
      tagger.tro_2_score     = cal_tro_2_bdt(0.35, tagger, reader_tro_2, tro_2_v_energy, tro_2_v_stem_length, tro_2_v_iso_angle, tro_2_v_max_length, tro_2_v_angle);
      tagger.tro_4_score     = cal_tro_4_bdt(0.33, tagger, reader_tro_4, tro_4_v_dir2_mag, tro_4_v_angle, tro_4_v_angle1, tro_4_v_angle2, tro_4_v_length, tro_4_v_length1, tro_4_v_medium_dQ_dx, tro_4_v_end_dQ_dx, tro_4_v_energy, tro_4_v_shower_main_length, tro_4_v_flag_shower_trajectory);
      tagger.tro_5_score     = cal_tro_5_bdt(0.5, tagger, reader_tro_5, tro_5_v_max_angle, tro_5_v_min_angle, tro_5_v_max_length, tro_5_v_iso_angle, tro_5_v_n_vtx_segs, tro_5_v_min_count, tro_5_v_max_count, tro_5_v_energy);

    
    Tbkg->Fill();
  }

  for (Int_t i=0;i!=t8->GetEntries();i++){
    t8->GetEntry(i);

    //    if (std::isnan(tagger.weight_cv) || std::isnan(tagger.weight_spline) || std::isinf(tagger.weight_cv) || std::isinf(tagger.weight_spline)) continue;
    // require EM shower ...
    if (tagger.br_filled!=1) continue;
    // only use even subrun number to train ...
    if (tagger.subrun %2 == 1 && process==1 || tagger.subrun%2 == 0 && process != 1) continue;
    //    if (tagger.subrun %2 == 1) continue;
    
    // also need to exclude the NC nu-electron elastic scattering 
    //if (tagger.truth_isCC==1 && abs(tagger.truth_nuPdg)==12 || tagger.truth_nuIntType==1098 ) continue;
    
    tagger.weight = (pot_1+pot_2+pot_3+pot_4)/(pot_7+pot_8);
    tagger.lowEweight = 1;

    tagger.br3_3_score     = cal_br3_3_bdt(0.3, tagger,  reader_br3_3, br3_3_v_energy,  br3_3_v_angle,  br3_3_v_dir_length, br3_3_v_length);
      tagger.br3_5_score     = cal_br3_5_bdt(0.42, tagger,  reader_br3_5, br3_5_v_dir_length, br3_5_v_total_length, br3_5_v_flag_avoid_muon_check, br3_5_v_n_seg, br3_5_v_angle, br3_5_v_sg_length, br3_5_v_energy, br3_5_v_n_main_segs, br3_5_v_n_segs, br3_5_v_shower_main_length, br3_5_v_shower_total_length);
      tagger.br3_6_score     = cal_br3_6_bdt(0.75, tagger, reader_br3_6, br3_6_v_angle, br3_6_v_angle1, br3_6_v_flag_shower_trajectory, br3_6_v_direct_length, br3_6_v_length, br3_6_v_n_other_vtx_segs, br3_6_v_energy);
      tagger.pio_2_score     = cal_pio_2_bdt(0.2,  tagger,  reader_pio_2, pio_2_v_dis2, pio_2_v_angle2, pio_2_v_acc_length);
      tagger.stw_2_score     = cal_stw_2_bdt(0.7, tagger, reader_stw_2, stw_2_v_medium_dQ_dx, stw_2_v_energy, stw_2_v_angle, stw_2_v_dir_length, stw_2_v_max_dQ_dx);
      tagger.stw_3_score     = cal_stw_3_bdt(0.5, tagger, reader_stw_3, stw_3_v_angle, stw_3_v_dir_length, stw_3_v_energy, stw_3_v_medium_dQ_dx);
      tagger.stw_4_score     = cal_stw_4_bdt(0.7, tagger, reader_stw_4, stw_4_v_angle, stw_4_v_dis, stw_4_v_energy);
      tagger.sig_1_score     = cal_sig_1_bdt(0.59, tagger,  reader_sig_1, sig_1_v_angle, sig_1_v_flag_single_shower, sig_1_v_energy, sig_1_v_energy_1);
      tagger.sig_2_score     = cal_sig_2_bdt(0.55, tagger, reader_sig_2, sig_2_v_energy, sig_2_v_shower_angle, sig_2_v_flag_single_shower, sig_2_v_medium_dQ_dx, sig_2_v_start_dQ_dx);
      tagger.lol_1_score     = cal_lol_1_bdt(0.85, tagger, reader_lol_1, lol_1_v_energy, lol_1_v_vtx_n_segs, lol_1_v_nseg, lol_1_v_angle);
      tagger.lol_2_score     = cal_lol_2_bdt(0.7, tagger, reader_lol_2,  lol_2_v_length, lol_2_v_angle, lol_2_v_type, lol_2_v_vtx_n_segs, lol_2_v_energy, lol_2_v_shower_main_length, lol_2_v_flag_dir_weak);
      tagger.tro_1_score     = cal_tro_1_bdt(0.28, tagger, reader_tro_1, tro_1_v_particle_type, tro_1_v_flag_dir_weak, tro_1_v_min_dis, tro_1_v_sg1_length,
  tro_1_v_shower_main_length, tro_1_v_max_n_vtx_segs, tro_1_v_tmp_length, tro_1_v_medium_dQ_dx, tro_1_v_dQ_dx_cut, tro_1_v_flag_shower_topology);
      tagger.tro_2_score     = cal_tro_2_bdt(0.35, tagger, reader_tro_2, tro_2_v_energy, tro_2_v_stem_length, tro_2_v_iso_angle, tro_2_v_max_length, tro_2_v_angle);
      tagger.tro_4_score     = cal_tro_4_bdt(0.33, tagger, reader_tro_4, tro_4_v_dir2_mag, tro_4_v_angle, tro_4_v_angle1, tro_4_v_angle2, tro_4_v_length, tro_4_v_length1, tro_4_v_medium_dQ_dx, tro_4_v_end_dQ_dx, tro_4_v_energy, tro_4_v_shower_main_length, tro_4_v_flag_shower_trajectory);
      tagger.tro_5_score     = cal_tro_5_bdt(0.5, tagger, reader_tro_5, tro_5_v_max_angle, tro_5_v_min_angle, tro_5_v_max_length, tro_5_v_iso_angle, tro_5_v_n_vtx_segs, tro_5_v_min_count, tro_5_v_max_count, tro_5_v_energy);

    
    Tbkg->Fill();
  }
  
  
  new_file->Write();
  new_file->Close();
  
  return 0;
}


