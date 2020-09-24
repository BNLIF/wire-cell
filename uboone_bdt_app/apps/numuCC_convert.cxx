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
#include "TMath.h"

#include "tagger.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

using namespace std;

float cal_numu_bdts_xgboost(TaggerInfo& tagger_info, TMVA::Reader& reader);

float cal_cosmict_10_bdt(float default_val,TaggerInfo& tagger_info, TMVA::Reader& reader,
			 float& cosmict_10_vtx_z,
			 float& cosmict_10_flag_shower,
			 float& cosmict_10_flag_dir_weak,
			 float& cosmict_10_angle_beam,
			 float& cosmict_10_length);

float cal_numu_1_bdt(float default_val,TaggerInfo& tagger_info,TMVA::Reader& reader,
		     float& numu_cc_flag_1,
		     float& numu_cc_1_particle_type,
		     float& numu_cc_1_length,
		     float& numu_cc_1_medium_dQ_dx,
		     float& numu_cc_1_dQ_dx_cut,
		     float& numu_cc_1_direct_length,
		     float& numu_cc_1_n_daughter_tracks,
		     float& numu_cc_1_n_daughter_all);
float cal_numu_2_bdt(float default_val,TaggerInfo& tagger_info,TMVA::Reader& reader,
		     float& numu_cc_2_length,
		     float& numu_cc_2_total_length,
		     float& numu_cc_2_n_daughter_tracks,
		     float& numu_cc_2_n_daughter_all);

float cal_numu_bdts_xgboost(TaggerInfo& tagger_info, TMVA::Reader& reader)
{
  float val = -10;

  double val1 = reader.EvaluateMVA("MyBDT");
  
  val = TMath::Log10( (1+val1)/(1-val1) );
  
  return val;
}

float cal_cosmict_10_bdt(float default_val,TaggerInfo& tagger_info, TMVA::Reader& reader,
			 float& cosmict_10_vtx_z,
			 float& cosmict_10_flag_shower,
			 float& cosmict_10_flag_dir_weak,
			 float& cosmict_10_angle_beam,
			 float& cosmict_10_length){
  float val = default_val;
  
  if (tagger_info.cosmict_10_length->size()>0){
    val = 1e9;
    for (size_t i=0;i!=tagger_info.cosmict_10_length->size();i++){
      cosmict_10_vtx_z = tagger_info.cosmict_10_vtx_z->at(i);
      cosmict_10_flag_shower = tagger_info.cosmict_10_flag_shower->at(i);
      cosmict_10_flag_dir_weak = tagger_info.cosmict_10_flag_dir_weak->at(i);
      cosmict_10_angle_beam = tagger_info.cosmict_10_angle_beam->at(i);
      cosmict_10_length = tagger_info.cosmict_10_length->at(i);

      if (std::isnan(cosmict_10_angle_beam)) cosmict_10_angle_beam = 0;
      
      float tmp_bdt =  reader.EvaluateMVA("MyBDT");
      if (tmp_bdt < val) val = tmp_bdt;
    }
  }

  return val;
}

float cal_numu_1_bdt(float default_val,TaggerInfo& tagger_info,TMVA::Reader& reader,
		     float& numu_cc_flag_1,
		     float& numu_cc_1_particle_type,
		     float& numu_cc_1_length,
		     float& numu_cc_1_medium_dQ_dx,
		     float& numu_cc_1_dQ_dx_cut,
		     float& numu_cc_1_direct_length,
		     float& numu_cc_1_n_daughter_tracks,
		     float& numu_cc_1_n_daughter_all){
  float val = default_val;
  
  
  if (tagger_info.numu_cc_1_particle_type->size()>0){
    val = -1e9;
    for (size_t i=0;i!=tagger_info.numu_cc_1_particle_type->size();i++){
      numu_cc_flag_1 = tagger_info.numu_cc_flag_1->at(i);
      numu_cc_1_particle_type= tagger_info.numu_cc_1_particle_type->at(i);
      numu_cc_1_length= tagger_info.numu_cc_1_length->at(i);
      numu_cc_1_medium_dQ_dx= tagger_info.numu_cc_1_medium_dQ_dx->at(i);
      numu_cc_1_dQ_dx_cut= tagger_info.numu_cc_1_dQ_dx_cut->at(i);
      numu_cc_1_direct_length= tagger_info.numu_cc_1_direct_length->at(i);
      numu_cc_1_n_daughter_tracks= tagger_info.numu_cc_1_n_daughter_tracks->at(i);
      numu_cc_1_n_daughter_all= tagger_info.numu_cc_1_n_daughter_all->at(i);

      if (std::isinf(numu_cc_1_dQ_dx_cut))  numu_cc_1_dQ_dx_cut = 10;
      
      float tmp_bdt =  reader.EvaluateMVA("MyBDT");
      if (tmp_bdt > val) val = tmp_bdt;
    }
  }
  
  return val;
}
float cal_numu_2_bdt(float default_val,TaggerInfo& tagger_info,TMVA::Reader& reader,
		     float& numu_cc_2_length,
		     float& numu_cc_2_total_length,
		     float& numu_cc_2_n_daughter_tracks,
		     float& numu_cc_2_n_daughter_all){
  float val = default_val;

  if (tagger_info.numu_cc_2_length->size()>0){
    val = -1e9;
    for (size_t i=0;i!=tagger_info.numu_cc_2_length->size();i++){
      numu_cc_2_length = tagger_info.numu_cc_2_length->at(i);
      numu_cc_2_total_length = tagger_info.numu_cc_2_total_length->at(i);
      numu_cc_2_n_daughter_tracks = tagger_info.numu_cc_2_n_daughter_tracks->at(i);
      numu_cc_2_n_daughter_all = tagger_info.numu_cc_2_n_daughter_all->at(i);
	
      float tmp_bdt =  reader.EvaluateMVA("MyBDT");
      if (tmp_bdt > val) val = tmp_bdt;
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
  
  // Now read a file and set address ...
  TFile *file2 = new TFile("run1_bnb_nu_POT1.2E21.root");
  TTree *t2 = (TTree*)file2->Get("bdt");
  set_tree_address(t2, tagger);
  double pot_2 = 6e20;

  TFile *file3 = new TFile("run3_bnb_nu_POT1.2E21.root");
  TTree *t3 = (TTree*)file3->Get("bdt");
  set_tree_address(t3, tagger);
  double pot_3 = 6e20;

  TFile *file4 = new TFile("run1_ext_bnb_C1_gt10_wcp_v00_14_00_POT1.2E20.root");
  TTree *t4 = (TTree*)file4->Get("bdt");
  set_tree_address(t4, tagger);
  double pot_4 = 1.2e20/2.;

  TFile *file5 = new TFile("run3_ext_bnb_F_G1_POT1.9E20.root");
  TTree *t5 = (TTree*)file5->Get("bdt");
  set_tree_address(t5, tagger);
  double pot_5 = 1.9e20/2.;

  TMVA::Reader reader_cosmict_10;
  float cosmict_10_vtx_z;
  float cosmict_10_flag_shower;
  float cosmict_10_flag_dir_weak;
  float cosmict_10_angle_beam;
  float cosmict_10_length;
  
  reader_cosmict_10.AddVariable("cosmict_10_vtx_z",&cosmict_10_vtx_z);
  reader_cosmict_10.AddVariable("cosmict_10_flag_shower",&cosmict_10_flag_shower);
  reader_cosmict_10.AddVariable("cosmict_10_flag_dir_weak",&cosmict_10_flag_dir_weak);
  reader_cosmict_10.AddVariable("cosmict_10_angle_beam",&cosmict_10_angle_beam);
  reader_cosmict_10.AddVariable("cosmict_10_length",&cosmict_10_length);
      
  reader_cosmict_10.BookMVA( "MyBDT", "weights/cos_tagger_10.weights.xml");


  TMVA::Reader reader_numu_1;
  
  float numu_cc_flag_1;
  float numu_cc_1_particle_type;
  float numu_cc_1_length;
  float numu_cc_1_medium_dQ_dx;
  float numu_cc_1_dQ_dx_cut;
  float numu_cc_1_direct_length;
  float numu_cc_1_n_daughter_tracks;
  float numu_cc_1_n_daughter_all;
  
  reader_numu_1.AddVariable("numu_cc_1_particle_type",&numu_cc_1_particle_type);
  reader_numu_1.AddVariable("numu_cc_1_length",&numu_cc_1_length);
  reader_numu_1.AddVariable("numu_cc_1_medium_dQ_dx",&numu_cc_1_medium_dQ_dx);
  reader_numu_1.AddVariable("numu_cc_1_dQ_dx_cut",&numu_cc_1_dQ_dx_cut);
  reader_numu_1.AddVariable("numu_cc_1_direct_length",&numu_cc_1_direct_length);
  reader_numu_1.AddVariable("numu_cc_1_n_daughter_tracks",&numu_cc_1_n_daughter_tracks);
  reader_numu_1.AddVariable("numu_cc_1_n_daughter_all",&numu_cc_1_n_daughter_all);
      
  reader_numu_1.BookMVA( "MyBDT", "weights/numu_tagger1.weights.xml");


  TMVA::Reader reader_numu_2;
  float numu_cc_2_length;
  float numu_cc_2_total_length;
  float numu_cc_2_n_daughter_tracks;
  float numu_cc_2_n_daughter_all;
  
  reader_numu_2.AddVariable("numu_cc_2_length",&numu_cc_2_length);
  reader_numu_2.AddVariable("numu_cc_2_total_length",&numu_cc_2_total_length);
  reader_numu_2.AddVariable("numu_cc_2_n_daughter_tracks",&numu_cc_2_n_daughter_tracks);
  reader_numu_2.AddVariable("numu_cc_2_n_daughter_all",&numu_cc_2_n_daughter_all);
  
  reader_numu_2.BookMVA( "MyBDT", "weights/numu_tagger2.weights.xml");


  TMVA::Reader reader;
  
  reader.AddVariable("numu_cc_flag_3", &tagger.numu_cc_flag_3);
  reader.AddVariable("numu_cc_3_particle_type", &tagger.numu_cc_3_particle_type);
  reader.AddVariable("numu_cc_3_max_length", &tagger.numu_cc_3_max_length);
  reader.AddVariable("numu_cc_3_track_length",&tagger.numu_cc_3_acc_track_length);
  reader.AddVariable("numu_cc_3_max_length_all",&tagger.numu_cc_3_max_length_all);
  reader.AddVariable("numu_cc_3_max_muon_length",&tagger.numu_cc_3_max_muon_length);
  reader.AddVariable("numu_cc_3_n_daughter_tracks",&tagger.numu_cc_3_n_daughter_tracks);
  reader.AddVariable("numu_cc_3_n_daughter_all",&tagger.numu_cc_3_n_daughter_all);
  reader.AddVariable("cosmict_flag_2", &tagger.cosmict_flag_2);
  reader.AddVariable("cosmict_2_filled", &tagger.cosmict_2_filled);
  reader.AddVariable("cosmict_2_particle_type",&tagger.cosmict_2_particle_type);
  reader.AddVariable("cosmict_2_n_muon_tracks",&tagger.cosmict_2_n_muon_tracks);
  reader.AddVariable("cosmict_2_total_shower_length",&tagger.cosmict_2_total_shower_length);
  reader.AddVariable("cosmict_2_flag_inside",&tagger.cosmict_2_flag_inside);
  reader.AddVariable("cosmict_2_angle_beam",&tagger.cosmict_2_angle_beam);
  reader.AddVariable("cosmict_2_flag_dir_weak", &tagger.cosmict_2_flag_dir_weak);
  reader.AddVariable("cosmict_2_dQ_dx_end", &tagger.cosmict_2_dQ_dx_end);
  reader.AddVariable("cosmict_2_dQ_dx_front", &tagger.cosmict_2_dQ_dx_front);
  reader.AddVariable("cosmict_2_theta", &tagger.cosmict_2_theta);
  reader.AddVariable("cosmict_2_phi", &tagger.cosmict_2_phi);
  reader.AddVariable("cosmict_2_valid_tracks", &tagger.cosmict_2_valid_tracks);
  reader.AddVariable("cosmict_flag_4", &tagger.cosmict_flag_4);
  reader.AddVariable("cosmict_4_filled", &tagger.cosmict_4_filled);
  reader.AddVariable("cosmict_4_flag_inside", &tagger.cosmict_4_flag_inside);
  reader.AddVariable("cosmict_4_angle_beam", &tagger.cosmict_4_angle_beam);
  reader.AddVariable("cosmict_4_connected_showers", &tagger.cosmict_4_connected_showers);
  reader.AddVariable("cosmict_flag_3", &tagger.cosmict_flag_3);
  reader.AddVariable("cosmict_3_filled", &tagger.cosmict_3_filled);
  reader.AddVariable("cosmict_3_flag_inside", &tagger.cosmict_3_flag_inside);
  reader.AddVariable("cosmict_3_angle_beam", &tagger.cosmict_3_angle_beam);
  reader.AddVariable("cosmict_3_flag_dir_weak", &tagger.cosmict_3_flag_dir_weak);
  reader.AddVariable("cosmict_3_dQ_dx_end", &tagger.cosmict_3_dQ_dx_end);
  reader.AddVariable("cosmict_3_dQ_dx_front", &tagger.cosmict_3_dQ_dx_front);
  reader.AddVariable("cosmict_3_theta", &tagger.cosmict_3_theta);
  reader.AddVariable("cosmict_3_phi", &tagger.cosmict_3_phi);
  reader.AddVariable("cosmict_3_valid_tracks", &tagger.cosmict_3_valid_tracks);
  reader.AddVariable("cosmict_flag_5", &tagger.cosmict_flag_5);
  reader.AddVariable("cosmict_5_filled", &tagger.cosmict_5_filled);
  reader.AddVariable("cosmict_5_flag_inside", &tagger.cosmict_5_flag_inside);
  reader.AddVariable("cosmict_5_angle_beam", &tagger.cosmict_5_angle_beam);
  reader.AddVariable("cosmict_5_connected_showers", &tagger.cosmict_5_connected_showers);
  reader.AddVariable("cosmict_flag_6", &tagger.cosmict_flag_6);
  reader.AddVariable("cosmict_6_filled", &tagger.cosmict_6_filled);
  reader.AddVariable("cosmict_6_flag_dir_weak", &tagger.cosmict_6_flag_dir_weak);
  reader.AddVariable("cosmict_6_flag_inside", &tagger.cosmict_6_flag_inside);
  reader.AddVariable("cosmict_6_angle", &tagger.cosmict_6_angle);
  reader.AddVariable("cosmict_flag_7", &tagger.cosmict_flag_7);
  reader.AddVariable("cosmict_7_filled", &tagger.cosmict_7_filled);
  reader.AddVariable("cosmict_7_flag_sec", &tagger.cosmict_7_flag_sec);
  reader.AddVariable("cosmict_7_n_muon_tracks", &tagger.cosmict_7_n_muon_tracks);
  reader.AddVariable("cosmict_7_total_shower_length", &tagger.cosmict_7_total_shower_length);
  reader.AddVariable("cosmict_7_flag_inside", &tagger.cosmict_7_flag_inside);
  reader.AddVariable("cosmict_7_angle_beam", &tagger.cosmict_7_angle_beam);
  reader.AddVariable("cosmict_7_flag_dir_weak", &tagger.cosmict_7_flag_dir_weak);
  reader.AddVariable("cosmict_7_dQ_dx_end", &tagger.cosmict_7_dQ_dx_end);
  reader.AddVariable("cosmict_7_dQ_dx_front", &tagger.cosmict_7_dQ_dx_front);
  reader.AddVariable("cosmict_7_theta", &tagger.cosmict_7_theta);
  reader.AddVariable("cosmict_7_phi", &tagger.cosmict_7_phi);
  reader.AddVariable("cosmict_flag_8", &tagger.cosmict_flag_8);
  reader.AddVariable("cosmict_8_filled", &tagger.cosmict_8_filled);
  reader.AddVariable("cosmict_8_flag_out", &tagger.cosmict_8_flag_out);
  reader.AddVariable("cosmict_8_muon_length", &tagger.cosmict_8_muon_length);
  reader.AddVariable("cosmict_8_acc_length", &tagger.cosmict_8_acc_length);
  reader.AddVariable("cosmict_flag_9", &tagger.cosmict_flag_9);
  reader.AddVariable("cosmic_flag", &tagger.cosmic_flag);
  reader.AddVariable("cosmic_filled", &tagger.cosmic_filled);
  reader.AddVariable("cosmict_flag", &tagger.cosmict_flag);
  reader.AddVariable("numu_cc_flag", &tagger.numu_cc_flag);
  reader.AddVariable("cosmict_flag_1", &tagger.cosmict_flag_1);
  reader.AddVariable("kine_reco_Enu",&tagger.kine_reco_Enu);
  reader.AddVariable("match_isFC",&tagger.match_isFC);
  reader.AddVariable("cosmict_10_score", &tagger.cosmict_10_score);
  reader.AddVariable("numu_1_score", &tagger.numu_1_score);
  reader.AddVariable("numu_2_score", &tagger.numu_2_score);

  reader.BookMVA( "MyBDT", "weights/numu_scalars_scores_0923.xml");
  

  TString filename;
  if (process == 1){
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

  // add numu vector and then scalar BDTs ... 

  

  
  
  for (Int_t i=0;i!=t2->GetEntries();i++){
    t2->GetEntry(i);

    // weight not good ...
    if (std::isnan(tagger.weight_cv) || std::isnan(tagger.weight_spline) || std::isinf(tagger.weight_cv) || std::isinf(tagger.weight_spline)) continue;

    // effectively generic neutrino selection
    //if (tagger.kine_reco_Enu == 0) continue;

    // odd sub run number ...
    if (tagger.subrun %2 == 1 && process == 1 || tagger.subrun %2 == 0 && process != 1) continue;

    if (tagger.weight_spline * tagger.weight_cv <=0 || tagger.weight_spline * tagger.weight_cv > 1000){
      tagger.weight = 1;
    }else{
      tagger.weight = tagger.weight_spline * tagger.weight_cv ;
    }
    
    tagger.lowEweight = 1;


    // BDT calculations
    tagger.numu_1_score = cal_numu_1_bdt(-0.4,tagger, reader_numu_1, numu_cc_flag_1,
    					      numu_cc_1_particle_type,
    					      numu_cc_1_length,
    					      numu_cc_1_medium_dQ_dx,
    					      numu_cc_1_dQ_dx_cut,
    					      numu_cc_1_direct_length,
    					      numu_cc_1_n_daughter_tracks,
    					      numu_cc_1_n_daughter_all);
    tagger.numu_2_score = cal_numu_2_bdt(-0.1,tagger,reader_numu_2,
					      numu_cc_2_length,
					      numu_cc_2_total_length,
					      numu_cc_2_n_daughter_tracks,
					      numu_cc_2_n_daughter_all);
    tagger.cosmict_10_score = cal_cosmict_10_bdt(0.7, tagger, reader_cosmict_10,
						      cosmict_10_vtx_z,
						      cosmict_10_flag_shower,
						      cosmict_10_flag_dir_weak,
						      cosmict_10_angle_beam,
						      cosmict_10_length);
    if (std::isnan(tagger.cosmict_4_angle_beam)) tagger.cosmict_4_angle_beam = 0;
    tagger.numu_score = cal_numu_bdts_xgboost(tagger,reader);
    
    if (tagger.truth_isCC==1 && abs(tagger.truth_nuPdg)==14 && tagger.truth_vtxInside ==1  ) {
      Tsig->Fill(); 
    }else{
      Tbkg->Fill();
    }
  }


  for (Int_t i=0;i!=t3->GetEntries();i++){
    t3->GetEntry(i);

    // weight not good ...
    if (std::isnan(tagger.weight_cv) || std::isnan(tagger.weight_spline) || std::isinf(tagger.weight_cv) || std::isinf(tagger.weight_spline)) continue;

    // effectively generic neutrino selection
    // if (tagger.kine_reco_Enu == 0) continue;

    // odd sub run number skip ...
    if (tagger.subrun %2 == 1 && process == 1 || tagger.subrun %2 == 0 && process != 1) continue;
    //    if (tagger.subrun %2 == 1) continue;
    if (tagger.weight_spline * tagger.weight_cv <=0 || tagger.weight_spline * tagger.weight_cv > 1000){
      tagger.weight = 1;
    }else{
      tagger.weight = tagger.weight_spline * tagger.weight_cv ;
    }
    tagger.lowEweight = 1;


    // BDT calculations
    tagger.numu_1_score = cal_numu_1_bdt(-0.4,tagger, reader_numu_1, numu_cc_flag_1,
    					      numu_cc_1_particle_type,
    					      numu_cc_1_length,
    					      numu_cc_1_medium_dQ_dx,
    					      numu_cc_1_dQ_dx_cut,
    					      numu_cc_1_direct_length,
    					      numu_cc_1_n_daughter_tracks,
    					      numu_cc_1_n_daughter_all);
    tagger.numu_2_score = cal_numu_2_bdt(-0.1,tagger,reader_numu_2,
					      numu_cc_2_length,
					      numu_cc_2_total_length,
					      numu_cc_2_n_daughter_tracks,
					      numu_cc_2_n_daughter_all);
    tagger.cosmict_10_score = cal_cosmict_10_bdt(0.7, tagger, reader_cosmict_10,
						      cosmict_10_vtx_z,
						      cosmict_10_flag_shower,
						      cosmict_10_flag_dir_weak,
						      cosmict_10_angle_beam,
						      cosmict_10_length);

    if (std::isnan(tagger.cosmict_4_angle_beam)) tagger.cosmict_4_angle_beam = 0;
    tagger.numu_score = cal_numu_bdts_xgboost(tagger,reader);
    
    if (tagger.truth_isCC==1 && abs(tagger.truth_nuPdg)==14 && tagger.truth_vtxInside ==1  ) {
      Tsig->Fill(); 
    }else{
      Tbkg->Fill();
    }
  }

  // run 1 ext
  for (Int_t i=0;i!=t4->GetEntries();i++){
    t4->GetEntry(i);

    // effectively generic neutrino selection
    // if (tagger.kine_reco_Enu == 0) continue;

    // odd sub run number ...
    if (tagger.subrun %2 == 1 && process == 1 || tagger.subrun %2 == 0 && process != 1) continue;
    //    if (tagger.subrun %2 == 1) continue;
    
    tagger.weight = (pot_2+pot_3)/(pot_4+pot_5);
    tagger.lowEweight = 1;

    // BDT calculations
    tagger.numu_1_score = cal_numu_1_bdt(-0.4,tagger, reader_numu_1, numu_cc_flag_1,
    					      numu_cc_1_particle_type,
    					      numu_cc_1_length,
    					      numu_cc_1_medium_dQ_dx,
    					      numu_cc_1_dQ_dx_cut,
    					      numu_cc_1_direct_length,
    					      numu_cc_1_n_daughter_tracks,
    					      numu_cc_1_n_daughter_all);
    tagger.numu_2_score = cal_numu_2_bdt(-0.1,tagger,reader_numu_2,
					      numu_cc_2_length,
					      numu_cc_2_total_length,
					      numu_cc_2_n_daughter_tracks,
					      numu_cc_2_n_daughter_all);
    tagger.cosmict_10_score = cal_cosmict_10_bdt(0.7, tagger, reader_cosmict_10,
						      cosmict_10_vtx_z,
						      cosmict_10_flag_shower,
						      cosmict_10_flag_dir_weak,
						      cosmict_10_angle_beam,
						      cosmict_10_length);
    
    if (std::isnan(tagger.cosmict_4_angle_beam)) tagger.cosmict_4_angle_beam = 0;
    tagger.numu_score = cal_numu_bdts_xgboost(tagger,reader);
    Tbkg->Fill();
    
  }

  // Run3  ext
  for (Int_t i=0;i!=t5->GetEntries();i++){
    t5->GetEntry(i);

    // effectively generic neutrino selection
    //if (tagger.kine_reco_Enu == 0) continue;

    // odd sub run number ...
    if (tagger.subrun %2 == 1 && process == 1 || tagger.subrun %2 == 0 && process != 1) continue;
    //    if (tagger.subrun %2 == 1) continue;
    
    tagger.weight =  (pot_2+pot_3)/(pot_4+pot_5);
    tagger.lowEweight = 1;

    // BDT calculations
    tagger.numu_1_score = cal_numu_1_bdt(-0.4,tagger, reader_numu_1, numu_cc_flag_1,
    					      numu_cc_1_particle_type,
    					      numu_cc_1_length,
    					      numu_cc_1_medium_dQ_dx,
    					      numu_cc_1_dQ_dx_cut,
    					      numu_cc_1_direct_length,
    					      numu_cc_1_n_daughter_tracks,
    					      numu_cc_1_n_daughter_all);
    tagger.numu_2_score = cal_numu_2_bdt(-0.1,tagger,reader_numu_2,
					      numu_cc_2_length,
					      numu_cc_2_total_length,
					      numu_cc_2_n_daughter_tracks,
					      numu_cc_2_n_daughter_all);
    tagger.cosmict_10_score = cal_cosmict_10_bdt(0.7, tagger, reader_cosmict_10,
						      cosmict_10_vtx_z,
						      cosmict_10_flag_shower,
						      cosmict_10_flag_dir_weak,
						      cosmict_10_angle_beam,
						      cosmict_10_length);
    if (std::isnan(tagger.cosmict_4_angle_beam)) tagger.cosmict_4_angle_beam = 0;
    tagger.numu_score = cal_numu_bdts_xgboost(tagger,reader);
    Tbkg->Fill();
   
  }
  
  

  
  

  // later need to add EXTBNB ...
  
  new_file->Write();
  new_file->Close();
  
  return 0;
}


