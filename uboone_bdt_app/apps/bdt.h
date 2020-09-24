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

float cal_bdts_xgboost(TaggerInfo& tagger_info, TMVA::Reader& reader);

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


float cal_bdts_xgboost(TaggerInfo& tagger_info, TMVA::Reader& reader){
  float val = 0; // background like ...
  float default_val = -15;

  if (tagger_info.br_filled==1){
    // protection of variables ...
    if(tagger_info.mip_min_dis>1000) tagger_info.mip_min_dis = 1000.0;
    if(tagger_info.mip_quality_shortest_length>1000) tagger_info.mip_quality_shortest_length = 1000;
    if(std::isnan(tagger_info.mip_quality_shortest_angle)) tagger_info.mip_quality_shortest_angle = 0;
    if(std::isnan(tagger_info.stem_dir_ratio)) tagger_info.stem_dir_ratio = 1.0; 
    
    
    double val1 = reader.EvaluateMVA("MyBDT");
    val = TMath::Log10( (1+val1)/(1-val1) );
  }else{
    val = default_val;
  }

  return val;
}
