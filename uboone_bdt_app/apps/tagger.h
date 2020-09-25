

 //
  struct TaggerInfo
  {
    // cosmic tagger ones, one case of cosmics ...
    float cosmic_flag; // default is true ...
    float cosmic_filled;   // if things are filled
    // variable ...
    float cosmic_n_solid_tracks;    // used
    float cosmic_energy_main_showers; // used
    float cosmic_energy_indirect_showers; // used
    // not used ...
    float cosmic_energy_direct_showers;  
    float cosmic_n_direct_showers;
    float cosmic_n_indirect_showers;
    float cosmic_n_main_showers;


    // shower gap identification
    float gap_flag;
    float gap_filled;
    // main variables:
    float gap_n_bad;
    float gap_n_points;
    float gap_energy;
    float gap_flag_single_shower;
    float gap_flag_parallel;
    // not directly used
    float gap_flag_prolong_u;
    float gap_flag_prolong_v;
    float gap_flag_prolong_w;
    float gap_num_valid_tracks;

    
    // mip_quality
    float mip_quality_flag;
    float mip_quality_filled;
    // variables
    float mip_quality_energy;
    float mip_quality_overlap;
    float mip_quality_n_showers;
    float mip_quality_n_tracks;
    float mip_quality_flag_inside_pi0;
    float mip_quality_n_pi0_showers;
    float mip_quality_shortest_length;
    float mip_quality_acc_length;
    float mip_quality_shortest_angle;
    float mip_quality_flag_proton;
    

    // mip identification
    float mip_flag;
    float mip_filled;
    // variables
    float mip_energy; // checked
    float mip_n_end_reduction;     // checked
    float mip_n_first_mip; // checked
    float mip_n_first_non_mip; // checked
    float mip_n_first_non_mip_1; // checked
    float mip_n_first_non_mip_2; // checked
    float mip_vec_dQ_dx_0; // checked
    float mip_vec_dQ_dx_1; // checked
    float mip_max_dQ_dx_sample; // checked
    float mip_n_below_threshold; // checked
    float mip_n_below_zero;  // checked
    float mip_n_lowest;  // checked
    float mip_n_highest;  // checked
    float mip_lowest_dQ_dx;  // checked
    float mip_highest_dQ_dx;  // checked
    float mip_medium_dQ_dx;  // checked
    float mip_stem_length;  // checked
    float mip_length_main;   // checked
    float mip_length_total; // checked
    float mip_angle_beam;  // checked
    float mip_iso_angle;  // checked
    float mip_n_vertex;   // checked
    float mip_n_good_tracks;  // checked
    float mip_E_indirect_max_energy;  // checked
    float mip_flag_all_above;     // checked
    float mip_min_dQ_dx_5;  // checked
    float mip_n_other_vertex;   // checked
    float mip_n_stem_size;  // checked
    float mip_flag_stem_trajectory;  // checked
    float mip_min_dis;   // checked


    // extra
    float mip_vec_dQ_dx_2;
    float mip_vec_dQ_dx_3;
    float mip_vec_dQ_dx_4;
    float mip_vec_dQ_dx_5;
    float mip_vec_dQ_dx_6;
    float mip_vec_dQ_dx_7;
    float mip_vec_dQ_dx_8;
    float mip_vec_dQ_dx_9;
    float mip_vec_dQ_dx_10;
    float mip_vec_dQ_dx_11;
    float mip_vec_dQ_dx_12;
    float mip_vec_dQ_dx_13;
    float mip_vec_dQ_dx_14;
    float mip_vec_dQ_dx_15;
    float mip_vec_dQ_dx_16;
    float mip_vec_dQ_dx_17;
    float mip_vec_dQ_dx_18;
    float mip_vec_dQ_dx_19;
    
    // shower pi0 identification
    float pio_flag;
    float pio_mip_id;    // pre-condition to run this tagger
    float pio_filled;
    float pio_flag_pio;  // only when this is true, we run the first part of tagger ...
    // first part of tagger ...
    float pio_1_flag;
    float pio_1_mass;
    float pio_1_pio_type;
    float pio_1_energy_1;
    float pio_1_energy_2;
    float pio_1_dis_1;
    float pio_1_dis_2;
    // second part of tagger
    std::vector<float> *pio_2_v_dis2;
    std::vector<float> *pio_2_v_angle2;
    std::vector<float> *pio_2_v_acc_length;
    std::vector<float> *pio_2_v_flag;

    
    //single shower pi0 case ...
    std::vector<float> *sig_1_v_angle;
    std::vector<float> *sig_1_v_flag_single_shower;
    std::vector<float> *sig_1_v_energy;
    std::vector<float> *sig_1_v_energy_1;
    std::vector<float> *sig_1_v_flag;

    std::vector<float> *sig_2_v_energy;
    std::vector<float> *sig_2_v_shower_angle;
    std::vector<float> *sig_2_v_flag_single_shower;
    std::vector<float> *sig_2_v_medium_dQ_dx;
    std::vector<float> *sig_2_v_start_dQ_dx;
    std::vector<float> *sig_2_v_flag;

    float sig_flag;

     // multiple gamma I
    float mgo_energy;
    float mgo_max_energy;
    float mgo_total_energy;
    float mgo_n_showers;
    float mgo_max_energy_1;
    float mgo_max_energy_2;
    float mgo_total_other_energy;
    float mgo_n_total_showers;
    float mgo_total_other_energy_1;
    float mgo_flag;
    
    // multiple gamma II
    float mgt_flag_single_shower;
    float mgt_max_energy;
    float mgt_energy;
    float mgt_total_other_energy;
    float mgt_max_energy_1;
    float mgt_e_indirect_max_energy;
    float mgt_e_direct_max_energy;
    float mgt_n_direct_showers;
    float mgt_e_direct_total_energy;
    float mgt_flag_indirect_max_pio;
    float mgt_e_indirect_total_energy;
    float mgt_flag;

    // shower to wall
    float stw_1_energy;
    float stw_1_dis;
    float stw_1_dQ_dx;
    float stw_1_flag_single_shower;
    float stw_1_n_pi0;
    float stw_1_num_valid_tracks;
    float stw_1_flag;

    std::vector<float> *stw_2_v_medium_dQ_dx;
    std::vector<float> *stw_2_v_energy;
    std::vector<float> *stw_2_v_angle;
    std::vector<float> *stw_2_v_dir_length;
    std::vector<float> *stw_2_v_max_dQ_dx;
    std::vector<float> *stw_2_v_flag;

    std::vector<float> *stw_3_v_angle;
    std::vector<float> *stw_3_v_dir_length;
    std::vector<float> *stw_3_v_energy;
    std::vector<float> *stw_3_v_medium_dQ_dx;
    std::vector<float> *stw_3_v_flag;

    std::vector<float> *stw_4_v_angle;
    std::vector<float> *stw_4_v_dis;
    std::vector<float> *stw_4_v_energy;
    std::vector<float> *stw_4_v_flag;

    float stw_flag;

    // single photon  cases ...
    float spt_flag_single_shower;
    float spt_energy;
    float spt_shower_main_length;
    float spt_shower_total_length;
    float spt_angle_beam;
    float spt_angle_vertical;
    float spt_max_dQ_dx;
    float spt_angle_beam_1;
    float spt_angle_drift;
    float spt_angle_drift_1;
    float spt_num_valid_tracks;
    float spt_n_vtx_segs;
    float spt_max_length;
    float spt_flag;

    // stem length ...
    float stem_len_energy;
    float stem_len_length;
    float stem_len_flag_avoid_muon_check;
    float stem_len_num_daughters;
    float stem_len_daughter_length;
    float stem_len_flag;

    // low-energy michel
    float lem_shower_total_length;
    float lem_shower_main_length;
    float lem_n_3seg;
    float lem_e_charge;
    float lem_e_dQdx;
    float lem_shower_num_segs;
    float lem_shower_num_main_segs;
    float lem_flag;

     // broken muon ...
    float brm_n_mu_segs;
    float brm_Ep;
    float brm_energy;
    float brm_acc_length;
    float brm_shower_total_length;
    float brm_connected_length;
    float brm_n_size;
    float brm_acc_direct_length;
    float brm_n_shower_main_segs;
    float brm_n_mu_main;
    float brm_flag;


    // compare with muon
    float cme_mu_energy;
    float cme_energy;
    float cme_mu_length;
    float cme_length;
    float cme_angle_beam;
    float cme_flag;


    // angular cut ...
    float anc_energy;
    float anc_angle;
    float anc_max_angle;
    float anc_max_length;
    float anc_acc_forward_length;
    float anc_acc_backward_length;
    float anc_acc_forward_length1;
    float anc_shower_main_length;
    float anc_shower_total_length;
    float anc_flag_main_outside;
    float anc_flag;

        
    // stem direction
    float stem_dir_flag;
    float stem_dir_flag_single_shower;
    float stem_dir_filled;
    float stem_dir_angle;
    float stem_dir_energy;
    float stem_dir_angle1;
    float stem_dir_angle2;
    float stem_dir_angle3;
    float stem_dir_ratio;

     // vertex inside shower
    float vis_1_filled;
    float vis_1_n_vtx_segs;
    float vis_1_energy;
    float vis_1_num_good_tracks;
    float vis_1_max_angle;
    float vis_1_max_shower_angle;
    float vis_1_tmp_length1;
    float vis_1_tmp_length2;
    float vis_1_particle_type;
    float vis_1_flag;

    float vis_2_filled;
    float vis_2_n_vtx_segs;
    float vis_2_min_angle;
    float vis_2_min_weak_track;
    float vis_2_angle_beam;
    float vis_2_min_angle1;
    float vis_2_iso_angle1;
    float vis_2_min_medium_dQ_dx;
    float vis_2_min_length;
    float vis_2_sg_length;
    float vis_2_max_angle;
    float vis_2_max_weak_track;
    float vis_2_flag;

    float vis_flag;


    // bad reconstruction_1
    float br_filled;
    float br1_flag;
    
    //bad reconstruction 1_1
    float br1_1_flag;
    float br1_1_shower_type;
    float br1_1_vtx_n_segs;
    float br1_1_energy;
    float br1_1_n_segs;
    float br1_1_flag_sg_topology;
    float br1_1_flag_sg_trajectory;
    float br1_1_sg_length;

    // bad reconstruction 1_2
    float br1_2_flag;
    float br1_2_energy;
    float br1_2_n_connected;
    float br1_2_max_length;
    float br1_2_n_connected_1;
    float br1_2_vtx_n_segs;
    float br1_2_n_shower_segs;
    float br1_2_max_length_ratio;
    float br1_2_shower_length;
    
    // bad_reconstruction 1_3
    float br1_3_flag;
    float br1_3_energy;
    float br1_3_n_connected_p;
    float br1_3_max_length_p;
    float br1_3_n_shower_segs;
    float br1_3_flag_sg_topology;
    float br1_3_flag_sg_trajectory;
    float br1_3_n_shower_main_segs;
    float br1_3_sg_length;


    // bad reconstruction 2
    float br2_flag;
    float br2_flag_single_shower;
    float br2_num_valid_tracks;
    float br2_energy;
    float br2_angle1;
    float br2_angle2;
    float br2_angle;
    float br2_angle3;
    float br2_n_shower_main_segs;
    float br2_max_angle;
    float br2_sg_length;
    float br2_flag_sg_trajectory;


     //bad reconstruction 3
    float br3_1_energy;
    float br3_1_n_shower_segments;
    float br3_1_sg_flag_trajectory;
    float br3_1_sg_direct_length;
    float br3_1_sg_length;
    float br3_1_total_main_length;
    float br3_1_total_length;
    float br3_1_iso_angle;
    float br3_1_sg_flag_topology;
    float br3_1_flag;

    float br3_2_n_ele;
    float br3_2_n_other;
    float br3_2_energy;
    float br3_2_total_main_length;
    float br3_2_total_length;
    float br3_2_other_fid;
    float br3_2_flag;

    std::vector<float> *br3_3_v_energy;
    std::vector<float> *br3_3_v_angle;
    std::vector<float> *br3_3_v_dir_length;
    std::vector<float> *br3_3_v_length;
    std::vector<float> *br3_3_v_flag;

    float br3_4_acc_length;
    float br3_4_total_length;
    float br3_4_energy;
    float br3_4_flag;

    std::vector<float> *br3_5_v_dir_length;
    std::vector<float> *br3_5_v_total_length;
    std::vector<float> *br3_5_v_flag_avoid_muon_check;
    std::vector<float> *br3_5_v_n_seg;
    std::vector<float> *br3_5_v_angle;
    std::vector<float> *br3_5_v_sg_length;
    std::vector<float> *br3_5_v_energy;
    std::vector<float> *br3_5_v_n_main_segs;
    std::vector<float> *br3_5_v_n_segs;
    std::vector<float> *br3_5_v_shower_main_length;
    std::vector<float> *br3_5_v_shower_total_length;
    std::vector<float> *br3_5_v_flag;

    std::vector<float> *br3_6_v_angle;
    std::vector<float> *br3_6_v_angle1;
    std::vector<float> *br3_6_v_flag_shower_trajectory;
    std::vector<float> *br3_6_v_direct_length;
    std::vector<float> *br3_6_v_length;
    std::vector<float> *br3_6_v_n_other_vtx_segs;
    std::vector<float> *br3_6_v_energy;
    std::vector<float> *br3_6_v_flag;

    float br3_7_energy;
    float br3_7_min_angle;
    float br3_7_sg_length;
    float br3_7_shower_main_length;
    float br3_7_flag;

    float br3_8_max_dQ_dx;
    float br3_8_energy;
    float br3_8_n_main_segs;
    float br3_8_shower_main_length;
    float br3_8_shower_length;
    float br3_8_flag;
    
    float br3_flag;
    
    // BR 4
    float br4_1_shower_main_length;
    float br4_1_shower_total_length;
    float br4_1_min_dis;
    float br4_1_energy;
    float br4_1_flag_avoid_muon_check;
    float br4_1_n_vtx_segs;
    float br4_1_n_main_segs;
    float br4_1_flag;

    float br4_2_ratio_45;
    float br4_2_ratio_35;
    float br4_2_ratio_25;
    float br4_2_ratio_15;
    float br4_2_energy;
    float br4_2_ratio1_45;
    float br4_2_ratio1_35;
    float br4_2_ratio1_25;
    float br4_2_ratio1_15;
    float br4_2_iso_angle;
    float br4_2_iso_angle1;
    float br4_2_angle;
    float br4_2_flag;

    float br4_flag;

     // track overclustering ..
    std::vector<float> *tro_1_v_particle_type;
    std::vector<float> *tro_1_v_flag_dir_weak;
    std::vector<float> *tro_1_v_min_dis;
    std::vector<float> *tro_1_v_sg1_length;
    std::vector<float> *tro_1_v_shower_main_length;
    std::vector<float> *tro_1_v_max_n_vtx_segs;
    std::vector<float> *tro_1_v_tmp_length;
    std::vector<float> *tro_1_v_medium_dQ_dx;
    std::vector<float> *tro_1_v_dQ_dx_cut;
    std::vector<float> *tro_1_v_flag_shower_topology;
    std::vector<float> *tro_1_v_flag;

    std::vector<float> *tro_2_v_energy;
    std::vector<float> *tro_2_v_stem_length;
    std::vector<float> *tro_2_v_iso_angle;
    std::vector<float> *tro_2_v_max_length;
    std::vector<float> *tro_2_v_angle;
    std::vector<float> *tro_2_v_flag;

    float tro_3_stem_length;
    float tro_3_n_muon_segs;
    float tro_3_energy;
    float tro_3_flag;

    std::vector<float> *tro_4_v_dir2_mag;
    std::vector<float> *tro_4_v_angle;
    std::vector<float> *tro_4_v_angle1;
    std::vector<float> *tro_4_v_angle2;
    std::vector<float> *tro_4_v_length;
    std::vector<float> *tro_4_v_length1;
    std::vector<float> *tro_4_v_medium_dQ_dx;
    std::vector<float> *tro_4_v_end_dQ_dx;
    std::vector<float> *tro_4_v_energy;
    std::vector<float> *tro_4_v_shower_main_length;
    std::vector<float> *tro_4_v_flag_shower_trajectory;
    std::vector<float> *tro_4_v_flag;

    std::vector<float> *tro_5_v_max_angle;
    std::vector<float> *tro_5_v_min_angle;
    std::vector<float> *tro_5_v_max_length;
    std::vector<float> *tro_5_v_iso_angle;
    std::vector<float> *tro_5_v_n_vtx_segs;
    std::vector<float> *tro_5_v_min_count;
    std::vector<float> *tro_5_v_max_count;
    std::vector<float> *tro_5_v_energy;
    std::vector<float> *tro_5_v_flag;

    float tro_flag;

    // high energy overlap 
    float hol_1_n_valid_tracks;
    float hol_1_min_angle;
    float hol_1_energy;
    float hol_1_flag_all_shower;
    float hol_1_min_length;
    float hol_1_flag;

    float hol_2_min_angle;
    float hol_2_medium_dQ_dx;
    float hol_2_ncount;
    float hol_2_energy;
    float hol_2_flag;

    float hol_flag;
    
    // low-energy overlap ...
    float lol_flag;

    std::vector<float> *lol_1_v_energy;
    std::vector<float> *lol_1_v_vtx_n_segs;
    std::vector<float> *lol_1_v_nseg;
    std::vector<float> *lol_1_v_angle;
    std::vector<float> *lol_1_v_flag;

    std::vector<float> *lol_2_v_length;
    std::vector<float> *lol_2_v_angle;
    std::vector<float> *lol_2_v_type;
    std::vector<float> *lol_2_v_vtx_n_segs;
    std::vector<float> *lol_2_v_energy;
    std::vector<float> *lol_2_v_shower_main_length;
    std::vector<float> *lol_2_v_flag_dir_weak;
    std::vector<float> *lol_2_v_flag;

    float lol_3_angle_beam;
    float lol_3_n_valid_tracks;
    float lol_3_min_angle;
    float lol_3_vtx_n_segs;
    float lol_3_energy;
    float lol_3_shower_main_length;
    float lol_3_n_out;
    float lol_3_n_sum;    
    float lol_3_flag;


    // cosmic tagger ...
    float cosmict_flag_1; // fiducial volume vertex
    float cosmict_flag_2;  // single muon
    float cosmict_flag_3;  // single muon (long)
    float cosmict_flag_4;  // kinematics muon
    float cosmict_flag_5; // kinematics muon (long)
    float cosmict_flag_6; // special ...
    float cosmict_flag_7;  // muon+ michel
    float cosmict_flag_8;  // muon + michel + special
    float cosmict_flag_9;  // this tagger is relevant for nueCC, see "cosmic tagger ones, one case of cosmics ..." (frist one ...)
    std::vector<float> *cosmict_flag_10;  // front upstream (dirt)
    float cosmict_flag;

    // single muon
    float cosmict_2_filled;
    float cosmict_2_particle_type;
    float cosmict_2_n_muon_tracks;
    float cosmict_2_total_shower_length;
    float cosmict_2_flag_inside;
    float cosmict_2_angle_beam;
    float cosmict_2_flag_dir_weak;
    float cosmict_2_dQ_dx_end;
    float cosmict_2_dQ_dx_front;
    float cosmict_2_theta;
    float cosmict_2_phi;
    float cosmict_2_valid_tracks;

    // signel muon (long)
    float cosmict_3_filled;
    float cosmict_3_flag_inside;
    float cosmict_3_angle_beam;
    float cosmict_3_flag_dir_weak;
    float cosmict_3_dQ_dx_end;
    float cosmict_3_dQ_dx_front;
    float cosmict_3_theta;
    float cosmict_3_phi;
    float cosmict_3_valid_tracks;

    // kinematics muon
    float cosmict_4_filled;
    float cosmict_4_flag_inside;
    float cosmict_4_angle_beam;
    float cosmict_4_connected_showers;  // need to be careful about the nueCC ...
    
    // kinematics muon (long)
    float cosmict_5_filled;
    float cosmict_5_flag_inside;
    float cosmict_5_angle_beam;
    float cosmict_5_connected_showers;

    // special
    float cosmict_6_filled;
    float cosmict_6_flag_dir_weak;
    float cosmict_6_flag_inside;
    float cosmict_6_angle;

    // muon + michel
    float cosmict_7_filled;
    float cosmict_7_flag_sec;
    float cosmict_7_n_muon_tracks;
    float cosmict_7_total_shower_length;
    float cosmict_7_flag_inside;
    float cosmict_7_angle_beam;
    float cosmict_7_flag_dir_weak;
    float cosmict_7_dQ_dx_end;
    float cosmict_7_dQ_dx_front;
    float cosmict_7_theta;
    float cosmict_7_phi;
    
    // muon + michel + special
    float cosmict_8_filled;
    float cosmict_8_flag_out;
    float cosmict_8_muon_length;
    float cosmict_8_acc_length;
    
    // front upstream (dirt)
    std::vector<float> *cosmict_10_flag_inside;
    std::vector<float> *cosmict_10_vtx_z;
    std::vector<float> *cosmict_10_flag_shower;
    std::vector<float> *cosmict_10_flag_dir_weak;
    std::vector<float> *cosmict_10_angle_beam;
    std::vector<float> *cosmict_10_length;
    
    // numu vs. nc tagger
    float numu_cc_flag;

    // single muon connected to neutrino vertex
    std::vector<float> *numu_cc_flag_1;
    std::vector<float> *numu_cc_1_particle_type;
    std::vector<float> *numu_cc_1_length;
    std::vector<float> *numu_cc_1_medium_dQ_dx;
    std::vector<float> *numu_cc_1_dQ_dx_cut;
    std::vector<float> *numu_cc_1_direct_length;
    std::vector<float> *numu_cc_1_n_daughter_tracks;
    std::vector<float> *numu_cc_1_n_daughter_all;
    
    // long muon connected to neutrino vertex
    std::vector<float> *numu_cc_flag_2;
    std::vector<float> *numu_cc_2_length;
    std::vector<float> *numu_cc_2_total_length;
    std::vector<float> *numu_cc_2_n_daughter_tracks;
    std::vector<float> *numu_cc_2_n_daughter_all;
    
    // any muon ...
    float numu_cc_flag_3;
    float numu_cc_3_particle_type;
    float numu_cc_3_max_length;
    float numu_cc_3_acc_track_length;
    float numu_cc_3_max_length_all;
    float numu_cc_3_max_muon_length;
    float numu_cc_3_n_daughter_tracks;
    float numu_cc_3_n_daughter_all;

    // numu_BDTs
    float cosmict_2_4_score;
    float cosmict_3_5_score;
    float cosmict_6_score;
    float cosmict_7_score;
    float cosmict_8_score;
    // vector ...
    float cosmict_10_score;
    
    // vector
    float numu_1_score;
    float numu_2_score;
    // scalar
    float numu_3_score;

    // total one
    float cosmict_score;
    float numu_score;
    
    
    
    // nue BDTs
    float mipid_score;
    float gap_score;
    float hol_lol_score;
    float cme_anc_score;
    float mgo_mgt_score;
    float br1_score;
    float br3_score;
    float br3_3_score;
    float br3_5_score;
    float br3_6_score;
    float stemdir_br2_score;
    float trimuon_score;
    float br4_tro_score;
    float mipquality_score;
    float pio_1_score;
    float pio_2_score;
    float stw_spt_score;
    float vis_1_score;
    float vis_2_score;
    float stw_2_score;
    float stw_3_score;
    float stw_4_score;
    float sig_1_score;
    float sig_2_score;
    float lol_1_score;
    float lol_2_score;
    float tro_1_score;
    float tro_2_score;
    float tro_4_score;
    float tro_5_score;
    float nue_score;

    // additional variables ...
    Int_t run;
    Int_t subrun;
    Int_t event;
    

    float nuvtx_diff;
    float showervtx_diff;
    float muonvtx_diff;

    float match_isFC;
    float truth_isCC;
    float truth_vtxInside;
    float truth_nuPdg;
    float truth_nuEnergy;
    float truth_nuIntType;
    float truth_energyInside;
    float weight_spline;
    float weight_cv;
    float weight_lee;
    
    float kine_reco_Enu;
    float kine_reco_add_energy;
    float kine_pio_mass;
    float kine_pio_vtx_dis;
    float kine_pio_flag;
    float kine_pio_energy_1;
    float kine_pio_theta_1;
    float kine_pio_phi_1;
    float kine_pio_dis_1;
    float kine_pio_energy_2;
    float kine_pio_theta_2;
    float kine_pio_phi_2;
    float kine_pio_dis_2;
    float kine_pio_angle;
    
    
    float event_type;


    // original ...
    float weight; // standard weight 
    float lowEweight; // extra weight for the training ...
    
  };

void set_tree_address(TTree *tree0, TaggerInfo& tagger_info, int flag = 1);
void put_tree_address(TTree *Tsig, TaggerInfo& tagger_info, int flag = 1);


void set_tree_address(TTree *tree0, TaggerInfo& tagger_info, int flag){

  
    // cosmic tagger
  tree0->SetBranchAddress("cosmic_flag", &tagger_info.cosmic_flag);
  tree0->SetBranchAddress("cosmic_n_solid_tracks",&tagger_info.cosmic_n_solid_tracks);
  tree0->SetBranchAddress("cosmic_energy_main_showers",&tagger_info.cosmic_energy_main_showers);
  tree0->SetBranchAddress("cosmic_energy_direct_showers",&tagger_info.cosmic_energy_direct_showers);
  tree0->SetBranchAddress("cosmic_energy_indirect_showers",&tagger_info.cosmic_energy_indirect_showers);
  tree0->SetBranchAddress("cosmic_n_direct_showers",&tagger_info.cosmic_n_direct_showers);
  tree0->SetBranchAddress("cosmic_n_indirect_showers",&tagger_info.cosmic_n_indirect_showers);
  tree0->SetBranchAddress("cosmic_n_main_showers",&tagger_info.cosmic_n_main_showers);
  tree0->SetBranchAddress("cosmic_filled",&tagger_info.cosmic_filled);
    
    // gap tagger
  tree0->SetBranchAddress("gap_flag",&tagger_info.gap_flag);
  tree0->SetBranchAddress("gap_flag_prolong_u",&tagger_info.gap_flag_prolong_u);
  tree0->SetBranchAddress("gap_flag_prolong_v",&tagger_info.gap_flag_prolong_v);
  tree0->SetBranchAddress("gap_flag_prolong_w",&tagger_info.gap_flag_prolong_w);
  tree0->SetBranchAddress("gap_flag_parallel",&tagger_info.gap_flag_parallel);
  tree0->SetBranchAddress("gap_n_points",&tagger_info.gap_n_points);
  tree0->SetBranchAddress("gap_n_bad",&tagger_info.gap_n_bad);
  tree0->SetBranchAddress("gap_energy",&tagger_info.gap_energy);
  tree0->SetBranchAddress("gap_num_valid_tracks",&tagger_info.gap_num_valid_tracks);
  tree0->SetBranchAddress("gap_flag_single_shower",&tagger_info.gap_flag_single_shower);
  tree0->SetBranchAddress("gap_filled",&tagger_info.gap_filled);
 
  // mip quality
  tree0->SetBranchAddress("mip_quality_flag",&tagger_info.mip_quality_flag);
  tree0->SetBranchAddress("mip_quality_energy",&tagger_info.mip_quality_energy);
  tree0->SetBranchAddress("mip_quality_overlap",&tagger_info.mip_quality_overlap);
  tree0->SetBranchAddress("mip_quality_n_showers",&tagger_info.mip_quality_n_showers);
  tree0->SetBranchAddress("mip_quality_n_tracks",&tagger_info.mip_quality_n_tracks);
  tree0->SetBranchAddress("mip_quality_flag_inside_pi0",&tagger_info.mip_quality_flag_inside_pi0);
  tree0->SetBranchAddress("mip_quality_n_pi0_showers",&tagger_info.mip_quality_n_pi0_showers);
  tree0->SetBranchAddress("mip_quality_shortest_length",&tagger_info.mip_quality_shortest_length);
  tree0->SetBranchAddress("mip_quality_acc_length",&tagger_info.mip_quality_acc_length);
  tree0->SetBranchAddress("mip_quality_shortest_angle",&tagger_info.mip_quality_shortest_angle);
  tree0->SetBranchAddress("mip_quality_flag_proton",&tagger_info.mip_quality_flag_proton);
  tree0->SetBranchAddress("mip_quality_filled",&tagger_info.mip_quality_filled);

    // mip
  tree0->SetBranchAddress("mip_flag",&tagger_info.mip_flag);
  tree0->SetBranchAddress("mip_energy",&tagger_info.mip_energy);
  tree0->SetBranchAddress("mip_n_end_reduction",&tagger_info.mip_n_end_reduction);
  tree0->SetBranchAddress("mip_n_first_mip",&tagger_info.mip_n_first_mip);
  tree0->SetBranchAddress("mip_n_first_non_mip",&tagger_info.mip_n_first_non_mip);
  tree0->SetBranchAddress("mip_n_first_non_mip_1",&tagger_info.mip_n_first_non_mip_1);
  tree0->SetBranchAddress("mip_n_first_non_mip_2",&tagger_info.mip_n_first_non_mip_2);

  tree0->SetBranchAddress("mip_vec_dQ_dx_0",&tagger_info.mip_vec_dQ_dx_0);
  tree0->SetBranchAddress("mip_vec_dQ_dx_1",&tagger_info.mip_vec_dQ_dx_1);
  tree0->SetBranchAddress("mip_vec_dQ_dx_2",&tagger_info.mip_vec_dQ_dx_2);
  tree0->SetBranchAddress("mip_vec_dQ_dx_3",&tagger_info.mip_vec_dQ_dx_3);
  tree0->SetBranchAddress("mip_vec_dQ_dx_4",&tagger_info.mip_vec_dQ_dx_4);
  tree0->SetBranchAddress("mip_vec_dQ_dx_5",&tagger_info.mip_vec_dQ_dx_5);
  tree0->SetBranchAddress("mip_vec_dQ_dx_6",&tagger_info.mip_vec_dQ_dx_6);
  tree0->SetBranchAddress("mip_vec_dQ_dx_7",&tagger_info.mip_vec_dQ_dx_7);
  tree0->SetBranchAddress("mip_vec_dQ_dx_8",&tagger_info.mip_vec_dQ_dx_8);
  tree0->SetBranchAddress("mip_vec_dQ_dx_9",&tagger_info.mip_vec_dQ_dx_9);
  tree0->SetBranchAddress("mip_vec_dQ_dx_10",&tagger_info.mip_vec_dQ_dx_10);
  tree0->SetBranchAddress("mip_vec_dQ_dx_11",&tagger_info.mip_vec_dQ_dx_11);
  tree0->SetBranchAddress("mip_vec_dQ_dx_12",&tagger_info.mip_vec_dQ_dx_12);
  tree0->SetBranchAddress("mip_vec_dQ_dx_13",&tagger_info.mip_vec_dQ_dx_13);
  tree0->SetBranchAddress("mip_vec_dQ_dx_14",&tagger_info.mip_vec_dQ_dx_14);
  tree0->SetBranchAddress("mip_vec_dQ_dx_15",&tagger_info.mip_vec_dQ_dx_15);
  tree0->SetBranchAddress("mip_vec_dQ_dx_16",&tagger_info.mip_vec_dQ_dx_16);
  tree0->SetBranchAddress("mip_vec_dQ_dx_17",&tagger_info.mip_vec_dQ_dx_17);
  tree0->SetBranchAddress("mip_vec_dQ_dx_18",&tagger_info.mip_vec_dQ_dx_18);
  tree0->SetBranchAddress("mip_vec_dQ_dx_19",&tagger_info.mip_vec_dQ_dx_19);

  tree0->SetBranchAddress("mip_max_dQ_dx_sample",&tagger_info.mip_max_dQ_dx_sample);
  tree0->SetBranchAddress("mip_n_below_threshold",&tagger_info.mip_n_below_threshold);
  tree0->SetBranchAddress("mip_n_below_zero",&tagger_info.mip_n_below_zero);
  tree0->SetBranchAddress("mip_n_lowest",&tagger_info.mip_n_lowest);
  tree0->SetBranchAddress("mip_n_highest",&tagger_info.mip_n_highest);

  tree0->SetBranchAddress("mip_lowest_dQ_dx",&tagger_info.mip_lowest_dQ_dx);
  tree0->SetBranchAddress("mip_highest_dQ_dx",&tagger_info.mip_highest_dQ_dx);
  tree0->SetBranchAddress("mip_medium_dQ_dx",&tagger_info.mip_medium_dQ_dx);
  tree0->SetBranchAddress("mip_stem_length",&tagger_info.mip_stem_length);
  tree0->SetBranchAddress("mip_length_main",&tagger_info.mip_length_main);
  tree0->SetBranchAddress("mip_length_total",&tagger_info.mip_length_total);
  tree0->SetBranchAddress("mip_angle_beam",&tagger_info.mip_angle_beam);
  tree0->SetBranchAddress("mip_iso_angle",&tagger_info.mip_iso_angle);

  tree0->SetBranchAddress("mip_n_vertex",&tagger_info.mip_n_vertex);
  tree0->SetBranchAddress("mip_n_good_tracks",&tagger_info.mip_n_good_tracks);
  tree0->SetBranchAddress("mip_E_indirect_max_energy",&tagger_info.mip_E_indirect_max_energy);
  tree0->SetBranchAddress("mip_flag_all_above",&tagger_info.mip_flag_all_above);
  tree0->SetBranchAddress("mip_min_dQ_dx_5",&tagger_info.mip_min_dQ_dx_5);
  tree0->SetBranchAddress("mip_n_other_vertex",&tagger_info.mip_n_other_vertex);
  tree0->SetBranchAddress("mip_n_stem_size",&tagger_info.mip_n_stem_size);
  tree0->SetBranchAddress("mip_flag_stem_trajectory",&tagger_info.mip_flag_stem_trajectory);
  tree0->SetBranchAddress("mip_min_dis",&tagger_info.mip_min_dis);
  tree0->SetBranchAddress("mip_filled",&tagger_info.mip_filled);

  // pio ...
  tree0->SetBranchAddress("pio_flag",&tagger_info.pio_flag);
  tree0->SetBranchAddress("pio_mip_id",&tagger_info.pio_mip_id);
  tree0->SetBranchAddress("pio_filled",&tagger_info.pio_filled);
  tree0->SetBranchAddress("pio_flag_pio",&tagger_info.pio_flag_pio);
    
  tree0->SetBranchAddress("pio_1_flag",&tagger_info.pio_1_flag);
  tree0->SetBranchAddress("pio_1_mass",&tagger_info.pio_1_mass);
  tree0->SetBranchAddress("pio_1_pio_type",&tagger_info.pio_1_pio_type);
  tree0->SetBranchAddress("pio_1_energy_1",&tagger_info.pio_1_energy_1);
  tree0->SetBranchAddress("pio_1_energy_2",&tagger_info.pio_1_energy_2);
  tree0->SetBranchAddress("pio_1_dis_1",&tagger_info.pio_1_dis_1);
  tree0->SetBranchAddress("pio_1_dis_2",&tagger_info.pio_1_dis_2);
    
  tree0->SetBranchAddress("pio_2_v_dis2",&tagger_info.pio_2_v_dis2);
  tree0->SetBranchAddress("pio_2_v_angle2",&tagger_info.pio_2_v_angle2);
  tree0->SetBranchAddress("pio_2_v_acc_length",&tagger_info.pio_2_v_acc_length);
  tree0->SetBranchAddress("pio_2_v_flag",&tagger_info.pio_2_v_flag);
  
  // bad reconstruction ...
  tree0->SetBranchAddress("stem_dir_flag",&tagger_info.stem_dir_flag);
  tree0->SetBranchAddress("stem_dir_flag_single_shower",&tagger_info.stem_dir_flag_single_shower);
  tree0->SetBranchAddress("stem_dir_filled",&tagger_info.stem_dir_filled);
  tree0->SetBranchAddress("stem_dir_angle",&tagger_info.stem_dir_angle);
  tree0->SetBranchAddress("stem_dir_energy",&tagger_info.stem_dir_energy);
  tree0->SetBranchAddress("stem_dir_angle1",&tagger_info.stem_dir_angle1);
  tree0->SetBranchAddress("stem_dir_angle2",&tagger_info.stem_dir_angle2);
  tree0->SetBranchAddress("stem_dir_angle3",&tagger_info.stem_dir_angle3);
  tree0->SetBranchAddress("stem_dir_ratio",&tagger_info.stem_dir_ratio);

  tree0->SetBranchAddress("br_filled",&tagger_info.br_filled);
    
  tree0->SetBranchAddress("br1_flag",&tagger_info.br1_flag);
    
  tree0->SetBranchAddress("br1_1_flag",&tagger_info.br1_1_flag);
  tree0->SetBranchAddress("br1_1_shower_type",&tagger_info.br1_1_shower_type);
  tree0->SetBranchAddress("br1_1_vtx_n_segs",&tagger_info.br1_1_vtx_n_segs);
  tree0->SetBranchAddress("br1_1_energy",&tagger_info.br1_1_energy);
  tree0->SetBranchAddress("br1_1_n_segs",&tagger_info.br1_1_n_segs);
  tree0->SetBranchAddress("br1_1_flag_sg_topology",&tagger_info.br1_1_flag_sg_topology);
  tree0->SetBranchAddress("br1_1_flag_sg_trajectory",&tagger_info.br1_1_flag_sg_trajectory);
  tree0->SetBranchAddress("br1_1_sg_length",&tagger_info.br1_1_sg_length);

  tree0->SetBranchAddress("br1_2_flag",&tagger_info.br1_2_flag);
  tree0->SetBranchAddress("br1_2_energy",&tagger_info.br1_2_energy);
  tree0->SetBranchAddress("br1_2_n_connected",&tagger_info.br1_2_n_connected);
  tree0->SetBranchAddress("br1_2_max_length",&tagger_info.br1_2_max_length);
  tree0->SetBranchAddress("br1_2_n_connected_1",&tagger_info.br1_2_n_connected_1);
  tree0->SetBranchAddress("br1_2_vtx_n_segs",&tagger_info.br1_2_vtx_n_segs);
  tree0->SetBranchAddress("br1_2_n_shower_segs",&tagger_info.br1_2_n_shower_segs);
  tree0->SetBranchAddress("br1_2_max_length_ratio",&tagger_info.br1_2_max_length_ratio);
  tree0->SetBranchAddress("br1_2_shower_length",&tagger_info.br1_2_shower_length);

  tree0->SetBranchAddress("br1_3_flag",&tagger_info.br1_3_flag);
  tree0->SetBranchAddress("br1_3_energy",&tagger_info.br1_3_energy);
  tree0->SetBranchAddress("br1_3_n_connected_p",&tagger_info.br1_3_n_connected_p);
  tree0->SetBranchAddress("br1_3_max_length_p",&tagger_info.br1_3_max_length_p);
  tree0->SetBranchAddress("br1_3_n_shower_segs",&tagger_info.br1_3_n_shower_segs);
  tree0->SetBranchAddress("br1_3_flag_sg_topology",&tagger_info.br1_3_flag_sg_topology);
  tree0->SetBranchAddress("br1_3_flag_sg_trajectory",&tagger_info.br1_3_flag_sg_trajectory);
  tree0->SetBranchAddress("br1_3_n_shower_main_segs",&tagger_info.br1_3_n_shower_main_segs);
  tree0->SetBranchAddress("br1_3_sg_length",&tagger_info.br1_3_sg_length);
    
  tree0->SetBranchAddress("br2_flag",&tagger_info.br2_flag);
  tree0->SetBranchAddress("br2_flag_single_shower",&tagger_info.br2_flag_single_shower);
  tree0->SetBranchAddress("br2_num_valid_tracks",&tagger_info.br2_num_valid_tracks);
  tree0->SetBranchAddress("br2_energy",&tagger_info.br2_energy);
  tree0->SetBranchAddress("br2_angle1",&tagger_info.br2_angle1);
  tree0->SetBranchAddress("br2_angle2",&tagger_info.br2_angle2);
  tree0->SetBranchAddress("br2_angle",&tagger_info.br2_angle);
  tree0->SetBranchAddress("br2_angle3",&tagger_info.br2_angle3);
  tree0->SetBranchAddress("br2_n_shower_main_segs",&tagger_info.br2_n_shower_main_segs);
  tree0->SetBranchAddress("br2_max_angle",&tagger_info.br2_max_angle);
  tree0->SetBranchAddress("br2_sg_length",&tagger_info.br2_sg_length);
  tree0->SetBranchAddress("br2_flag_sg_trajectory",&tagger_info.br2_flag_sg_trajectory);


  tree0->SetBranchAddress("lol_flag",&tagger_info.lol_flag);

  tree0->SetBranchAddress("lol_1_v_energy",&tagger_info.lol_1_v_energy);
  tree0->SetBranchAddress("lol_1_v_vtx_n_segs",&tagger_info.lol_1_v_vtx_n_segs);
  tree0->SetBranchAddress("lol_1_v_nseg",&tagger_info.lol_1_v_nseg);
  tree0->SetBranchAddress("lol_1_v_angle",&tagger_info.lol_1_v_angle);
  tree0->SetBranchAddress("lol_1_v_flag",&tagger_info.lol_1_v_flag);
  
  tree0->SetBranchAddress("lol_2_v_length",&tagger_info.lol_2_v_length);
  tree0->SetBranchAddress("lol_2_v_angle",&tagger_info.lol_2_v_angle);
  tree0->SetBranchAddress("lol_2_v_type",&tagger_info.lol_2_v_type);
  tree0->SetBranchAddress("lol_2_v_vtx_n_segs",&tagger_info.lol_2_v_vtx_n_segs);
  tree0->SetBranchAddress("lol_2_v_energy",&tagger_info.lol_2_v_energy);
  tree0->SetBranchAddress("lol_2_v_shower_main_length",&tagger_info.lol_2_v_shower_main_length);
  tree0->SetBranchAddress("lol_2_v_flag_dir_weak",&tagger_info.lol_2_v_flag_dir_weak);
  tree0->SetBranchAddress("lol_2_v_flag",&tagger_info.lol_2_v_flag);
  
  tree0->SetBranchAddress("lol_3_angle_beam",&tagger_info.lol_3_angle_beam);
  tree0->SetBranchAddress("lol_3_n_valid_tracks",&tagger_info.lol_3_n_valid_tracks);
  tree0->SetBranchAddress("lol_3_min_angle",&tagger_info.lol_3_min_angle);
  tree0->SetBranchAddress("lol_3_vtx_n_segs",&tagger_info.lol_3_vtx_n_segs);
  tree0->SetBranchAddress("lol_3_energy",&tagger_info.lol_3_energy);
  tree0->SetBranchAddress("lol_3_shower_main_length",&tagger_info.lol_3_shower_main_length);
  tree0->SetBranchAddress("lol_3_n_out",&tagger_info.lol_3_n_out);
  tree0->SetBranchAddress("lol_3_n_sum",&tagger_info.lol_3_n_sum);
  tree0->SetBranchAddress("lol_3_flag",&tagger_info.lol_3_flag);
    
  tree0->SetBranchAddress("br3_1_energy",&tagger_info.br3_1_energy);
  tree0->SetBranchAddress("br3_1_n_shower_segments",&tagger_info.br3_1_n_shower_segments);
  tree0->SetBranchAddress("br3_1_sg_flag_trajectory",&tagger_info.br3_1_sg_flag_trajectory);
  tree0->SetBranchAddress("br3_1_sg_direct_length",&tagger_info.br3_1_sg_direct_length);
  tree0->SetBranchAddress("br3_1_sg_length",&tagger_info.br3_1_sg_length);
  tree0->SetBranchAddress("br3_1_total_main_length",&tagger_info.br3_1_total_main_length);
  tree0->SetBranchAddress("br3_1_total_length",&tagger_info.br3_1_total_length);
  tree0->SetBranchAddress("br3_1_iso_angle",&tagger_info.br3_1_iso_angle);
  tree0->SetBranchAddress("br3_1_sg_flag_topology",&tagger_info.br3_1_sg_flag_topology);
  tree0->SetBranchAddress("br3_1_flag",&tagger_info.br3_1_flag);

  tree0->SetBranchAddress("br3_2_n_ele",&tagger_info.br3_2_n_ele);
  tree0->SetBranchAddress("br3_2_n_other",&tagger_info.br3_2_n_other);
  tree0->SetBranchAddress("br3_2_energy",&tagger_info.br3_2_energy);
  tree0->SetBranchAddress("br3_2_total_main_length",&tagger_info.br3_2_total_main_length);
  tree0->SetBranchAddress("br3_2_total_length",&tagger_info.br3_2_total_length);
  tree0->SetBranchAddress("br3_2_other_fid",&tagger_info.br3_2_other_fid);
  tree0->SetBranchAddress("br3_2_flag",&tagger_info.br3_2_flag);

  tree0->SetBranchAddress("br3_3_v_energy",&tagger_info.br3_3_v_energy);
  tree0->SetBranchAddress("br3_3_v_angle",&tagger_info.br3_3_v_angle);
  tree0->SetBranchAddress("br3_3_v_dir_length",&tagger_info.br3_3_v_dir_length);
  tree0->SetBranchAddress("br3_3_v_length",&tagger_info.br3_3_v_length);
  tree0->SetBranchAddress("br3_3_v_flag",&tagger_info.br3_3_v_flag);
  
  tree0->SetBranchAddress("br3_4_acc_length", &tagger_info.br3_4_acc_length);
  tree0->SetBranchAddress("br3_4_total_length", &tagger_info.br3_4_total_length);
  tree0->SetBranchAddress("br3_4_energy", &tagger_info.br3_4_energy);
  tree0->SetBranchAddress("br3_4_flag", &tagger_info.br3_4_flag);
    
  tree0->SetBranchAddress("br3_5_v_dir_length", &tagger_info.br3_5_v_dir_length);
  tree0->SetBranchAddress("br3_5_v_total_length", &tagger_info.br3_5_v_total_length);
  tree0->SetBranchAddress("br3_5_v_flag_avoid_muon_check", &tagger_info.br3_5_v_flag_avoid_muon_check);
  tree0->SetBranchAddress("br3_5_v_n_seg", &tagger_info.br3_5_v_n_seg);
  tree0->SetBranchAddress("br3_5_v_angle", &tagger_info.br3_5_v_angle);
  tree0->SetBranchAddress("br3_5_v_sg_length", &tagger_info.br3_5_v_sg_length);
  tree0->SetBranchAddress("br3_5_v_energy", &tagger_info.br3_5_v_energy);
  tree0->SetBranchAddress("br3_5_v_n_main_segs", &tagger_info.br3_5_v_n_main_segs);
  tree0->SetBranchAddress("br3_5_v_n_segs", &tagger_info.br3_5_v_n_segs);
  tree0->SetBranchAddress("br3_5_v_shower_main_length", &tagger_info.br3_5_v_shower_main_length);
  tree0->SetBranchAddress("br3_5_v_shower_total_length", &tagger_info.br3_5_v_shower_total_length);
  tree0->SetBranchAddress("br3_5_v_flag", &tagger_info.br3_5_v_flag);
  
  tree0->SetBranchAddress("br3_6_v_angle",&tagger_info.br3_6_v_angle);
  tree0->SetBranchAddress("br3_6_v_angle1",&tagger_info.br3_6_v_angle1);
  tree0->SetBranchAddress("br3_6_v_flag_shower_trajectory",&tagger_info.br3_6_v_flag_shower_trajectory);
  tree0->SetBranchAddress("br3_6_v_direct_length",&tagger_info.br3_6_v_direct_length);
  tree0->SetBranchAddress("br3_6_v_length",&tagger_info.br3_6_v_length);
  tree0->SetBranchAddress("br3_6_v_n_other_vtx_segs",&tagger_info.br3_6_v_n_other_vtx_segs);
  tree0->SetBranchAddress("br3_6_v_energy",&tagger_info.br3_6_v_energy);
  tree0->SetBranchAddress("br3_6_v_flag",&tagger_info.br3_6_v_flag);
  
  tree0->SetBranchAddress("br3_7_energy",&tagger_info.br3_7_energy);
  tree0->SetBranchAddress("br3_7_min_angle",&tagger_info.br3_7_min_angle);
  tree0->SetBranchAddress("br3_7_sg_length",&tagger_info.br3_7_sg_length);
  tree0->SetBranchAddress("br3_7_main_length",&tagger_info.br3_7_shower_main_length);
  tree0->SetBranchAddress("br3_7_flag",&tagger_info.br3_7_flag);

  tree0->SetBranchAddress("br3_8_max_dQ_dx",&tagger_info.br3_8_max_dQ_dx);
  tree0->SetBranchAddress("br3_8_energy",&tagger_info.br3_8_energy);
  tree0->SetBranchAddress("br3_8_n_main_segs",&tagger_info.br3_8_n_main_segs);
  tree0->SetBranchAddress("br3_8_shower_main_length",&tagger_info.br3_8_shower_main_length);
  tree0->SetBranchAddress("br3_8_shower_length",&tagger_info.br3_8_shower_length);
  tree0->SetBranchAddress("br3_8_flag",&tagger_info.br3_8_flag);

  tree0->SetBranchAddress("br3_flag",&tagger_info.br3_flag);


  tree0->SetBranchAddress("br4_1_shower_main_length", &tagger_info.br4_1_shower_main_length);
  tree0->SetBranchAddress("br4_1_shower_total_length", &tagger_info.br4_1_shower_total_length);
  tree0->SetBranchAddress("br4_1_min_dis", &tagger_info.br4_1_min_dis);
  tree0->SetBranchAddress("br4_1_energy", &tagger_info.br4_1_energy);
  tree0->SetBranchAddress("br4_1_flag_avoid_muon_check", &tagger_info.br4_1_flag_avoid_muon_check);
  tree0->SetBranchAddress("br4_1_n_vtx_segs", &tagger_info.br4_1_n_vtx_segs);
  tree0->SetBranchAddress("br4_1_n_main_segs", &tagger_info.br4_1_n_main_segs);
  tree0->SetBranchAddress("br4_1_flag", &tagger_info.br4_1_flag);

  tree0->SetBranchAddress("br4_2_ratio_45", &tagger_info.br4_2_ratio_45);
  tree0->SetBranchAddress("br4_2_ratio_35", &tagger_info.br4_2_ratio_35);
  tree0->SetBranchAddress("br4_2_ratio_25", &tagger_info.br4_2_ratio_25);
  tree0->SetBranchAddress("br4_2_ratio_15", &tagger_info.br4_2_ratio_15);
  tree0->SetBranchAddress("br4_2_energy",   &tagger_info.br4_2_energy);
  tree0->SetBranchAddress("br4_2_ratio1_45", &tagger_info.br4_2_ratio1_45);
  tree0->SetBranchAddress("br4_2_ratio1_35", &tagger_info.br4_2_ratio1_35);
  tree0->SetBranchAddress("br4_2_ratio1_25", &tagger_info.br4_2_ratio1_25);
  tree0->SetBranchAddress("br4_2_ratio1_15", &tagger_info.br4_2_ratio1_15);
  tree0->SetBranchAddress("br4_2_iso_angle", &tagger_info.br4_2_iso_angle);
  tree0->SetBranchAddress("br4_2_iso_angle1", &tagger_info.br4_2_iso_angle1);
  tree0->SetBranchAddress("br4_2_angle", &tagger_info.br4_2_angle);
  tree0->SetBranchAddress("br4_2_flag", &tagger_info.br4_2_flag);

  tree0->SetBranchAddress("br4_flag", &tagger_info.br4_flag);
    

  tree0->SetBranchAddress("hol_1_n_valid_tracks", &tagger_info.hol_1_n_valid_tracks);
  tree0->SetBranchAddress("hol_1_min_angle", &tagger_info.hol_1_min_angle);
  tree0->SetBranchAddress("hol_1_energy", &tagger_info.hol_1_energy);
  tree0->SetBranchAddress("hol_1_flag_all_shower", &tagger_info.hol_1_flag_all_shower);
  tree0->SetBranchAddress("hol_1_min_length", &tagger_info.hol_1_min_length);
  tree0->SetBranchAddress("hol_1_flag", &tagger_info.hol_1_flag);

  tree0->SetBranchAddress("hol_2_min_angle", &tagger_info.hol_2_min_angle);
  tree0->SetBranchAddress("hol_2_medium_dQ_dx", &tagger_info.hol_2_medium_dQ_dx);
  tree0->SetBranchAddress("hol_2_ncount", &tagger_info.hol_2_ncount);
  tree0->SetBranchAddress("hol_2_energy", &tagger_info.hol_2_energy);
  tree0->SetBranchAddress("hol_2_flag", &tagger_info.hol_2_flag);

  tree0->SetBranchAddress("hol_flag", &tagger_info.hol_flag);
    

  tree0->SetBranchAddress("vis_1_filled",&tagger_info.vis_1_filled);
  tree0->SetBranchAddress("vis_1_n_vtx_segs",&tagger_info.vis_1_n_vtx_segs);
  tree0->SetBranchAddress("vis_1_energy",&tagger_info.vis_1_energy);
  tree0->SetBranchAddress("vis_1_num_good_tracks",&tagger_info.vis_1_num_good_tracks);
  tree0->SetBranchAddress("vis_1_max_angle",&tagger_info.vis_1_max_angle);
  tree0->SetBranchAddress("vis_1_max_shower_angle",&tagger_info.vis_1_max_shower_angle);
  tree0->SetBranchAddress("vis_1_tmp_length1",&tagger_info.vis_1_tmp_length1);
  tree0->SetBranchAddress("vis_1_tmp_length2",&tagger_info.vis_1_tmp_length2);
  tree0->SetBranchAddress("vis_1_particle_type",&tagger_info.vis_1_particle_type);
  tree0->SetBranchAddress("vis_1_flag",&tagger_info.vis_1_flag);

  tree0->SetBranchAddress("vis_2_filled",&tagger_info.vis_2_filled);
  tree0->SetBranchAddress("vis_2_n_vtx_segs",&tagger_info.vis_2_n_vtx_segs);
  tree0->SetBranchAddress("vis_2_min_angle",&tagger_info.vis_2_min_angle);
  tree0->SetBranchAddress("vis_2_min_weak_track",&tagger_info.vis_2_min_weak_track);
  tree0->SetBranchAddress("vis_2_angle_beam",&tagger_info.vis_2_angle_beam);
  tree0->SetBranchAddress("vis_2_min_angle1",&tagger_info.vis_2_min_angle1);
  tree0->SetBranchAddress("vis_2_iso_angle1",&tagger_info.vis_2_iso_angle1);
  tree0->SetBranchAddress("vis_2_min_medium_dQ_dx",&tagger_info.vis_2_min_medium_dQ_dx);
  tree0->SetBranchAddress("vis_2_min_length",&tagger_info.vis_2_min_length);
  tree0->SetBranchAddress("vis_2_sg_length",&tagger_info.vis_2_sg_length);
  tree0->SetBranchAddress("vis_2_max_angle",&tagger_info.vis_2_max_angle);
  tree0->SetBranchAddress("vis_2_max_weak_track",&tagger_info.vis_2_max_weak_track);
  tree0->SetBranchAddress("vis_2_flag",&tagger_info.vis_2_flag);
  
  tree0->SetBranchAddress("vis_flag",&tagger_info.vis_flag);
  

  tree0->SetBranchAddress("stem_len_energy", &tagger_info.stem_len_energy);
  tree0->SetBranchAddress("stem_len_length", &tagger_info.stem_len_length);
  tree0->SetBranchAddress("stem_len_flag_avoid_muon_check", &tagger_info.stem_len_flag_avoid_muon_check);
  tree0->SetBranchAddress("stem_len_num_daughters", &tagger_info.stem_len_num_daughters);
  tree0->SetBranchAddress("stem_len_daughter_length", &tagger_info.stem_len_daughter_length);
  tree0->SetBranchAddress("stem_len_flag", &tagger_info.stem_len_flag);

  tree0->SetBranchAddress("brm_n_mu_segs",&tagger_info.brm_n_mu_segs);
  tree0->SetBranchAddress("brm_Ep",&tagger_info.brm_Ep);
  tree0->SetBranchAddress("brm_energy",&tagger_info.brm_energy);
  tree0->SetBranchAddress("brm_acc_length",&tagger_info.brm_acc_length);
  tree0->SetBranchAddress("brm_shower_total_length",&tagger_info.brm_shower_total_length);
  tree0->SetBranchAddress("brm_connected_length",&tagger_info.brm_connected_length);
  tree0->SetBranchAddress("brm_n_size",&tagger_info.brm_n_size);
  tree0->SetBranchAddress("brm_acc_direct_length",&tagger_info.brm_acc_direct_length);
  tree0->SetBranchAddress("brm_n_shower_main_segs",&tagger_info.brm_n_shower_main_segs);
  tree0->SetBranchAddress("brm_n_mu_main",&tagger_info.brm_n_mu_main);
  tree0->SetBranchAddress("brm_flag",&tagger_info.brm_flag);

  tree0->SetBranchAddress("cme_mu_energy",&tagger_info.cme_mu_energy);
  tree0->SetBranchAddress("cme_energy",&tagger_info.cme_energy);
  tree0->SetBranchAddress("cme_mu_length",&tagger_info.cme_mu_length);
  tree0->SetBranchAddress("cme_length",&tagger_info.cme_length);
  tree0->SetBranchAddress("cme_angle_beam",&tagger_info.cme_angle_beam);
  tree0->SetBranchAddress("cme_flag",&tagger_info.cme_flag);
  
  tree0->SetBranchAddress("anc_energy",&tagger_info.anc_energy);
  tree0->SetBranchAddress("anc_angle",&tagger_info.anc_angle);
  tree0->SetBranchAddress("anc_max_angle",&tagger_info.anc_max_angle);
  tree0->SetBranchAddress("anc_max_length",&tagger_info.anc_max_length);
  tree0->SetBranchAddress("anc_acc_forward_length",&tagger_info.anc_acc_forward_length);
  tree0->SetBranchAddress("anc_acc_backward_length",&tagger_info.anc_acc_backward_length);
  tree0->SetBranchAddress("anc_acc_forward_length1",&tagger_info.anc_acc_forward_length1);
  tree0->SetBranchAddress("anc_shower_main_length",&tagger_info.anc_shower_main_length);
  tree0->SetBranchAddress("anc_shower_total_length",&tagger_info.anc_shower_total_length);
  tree0->SetBranchAddress("anc_flag_main_outside",&tagger_info.anc_flag_main_outside);
  tree0->SetBranchAddress("anc_flag",&tagger_info.anc_flag);

  tree0->SetBranchAddress("lem_shower_total_length",&tagger_info.lem_shower_total_length);
  tree0->SetBranchAddress("lem_shower_main_length",&tagger_info.lem_shower_main_length);
  tree0->SetBranchAddress("lem_n_3seg",&tagger_info.lem_n_3seg);
  tree0->SetBranchAddress("lem_e_charge",&tagger_info.lem_e_charge);
  tree0->SetBranchAddress("lem_e_dQdx",&tagger_info.lem_e_dQdx);
  tree0->SetBranchAddress("lem_shower_num_segs",&tagger_info.lem_shower_num_segs);
  tree0->SetBranchAddress("lem_shower_num_main_segs",&tagger_info.lem_shower_num_main_segs);
  tree0->SetBranchAddress("lem_flag",&tagger_info.lem_flag);

  tree0->SetBranchAddress("stw_1_energy",&tagger_info.stw_1_energy);
  tree0->SetBranchAddress("stw_1_dis",&tagger_info.stw_1_dis);
  tree0->SetBranchAddress("stw_1_dQ_dx",&tagger_info.stw_1_dQ_dx);
  tree0->SetBranchAddress("stw_1_flag_single_shower",&tagger_info.stw_1_flag_single_shower);
  tree0->SetBranchAddress("stw_1_n_pi0",&tagger_info.stw_1_n_pi0);
  tree0->SetBranchAddress("stw_1_num_valid_tracks",&tagger_info.stw_1_num_valid_tracks);
  tree0->SetBranchAddress("stw_1_flag",&tagger_info.stw_1_flag);
  
  tree0->SetBranchAddress("stw_2_v_medium_dQ_dx", &tagger_info.stw_2_v_medium_dQ_dx);
  tree0->SetBranchAddress("stw_2_v_energy", &tagger_info.stw_2_v_energy);
  tree0->SetBranchAddress("stw_2_v_angle", &tagger_info.stw_2_v_angle);
  tree0->SetBranchAddress("stw_2_v_dir_length", &tagger_info.stw_2_v_dir_length);
  tree0->SetBranchAddress("stw_2_v_max_dQ_dx", &tagger_info.stw_2_v_max_dQ_dx);
  tree0->SetBranchAddress("stw_2_v_flag", &tagger_info.stw_2_v_flag);
  
  tree0->SetBranchAddress("stw_3_v_angle",&tagger_info.stw_3_v_angle);
  tree0->SetBranchAddress("stw_3_v_dir_length",&tagger_info.stw_3_v_dir_length);
  tree0->SetBranchAddress("stw_3_v_energy",&tagger_info.stw_3_v_energy);
  tree0->SetBranchAddress("stw_3_v_medium_dQ_dx",&tagger_info.stw_3_v_medium_dQ_dx);
  tree0->SetBranchAddress("stw_3_v_flag",&tagger_info.stw_3_v_flag);
  
  tree0->SetBranchAddress("stw_4_v_angle",&tagger_info.stw_4_v_angle);
  tree0->SetBranchAddress("stw_4_v_dis",&tagger_info.stw_4_v_dis);
  tree0->SetBranchAddress("stw_4_v_energy",&tagger_info.stw_4_v_energy);
  tree0->SetBranchAddress("stw_4_v_flag",&tagger_info.stw_4_v_flag);
  
  tree0->SetBranchAddress("stw_flag", &tagger_info.stw_flag);

  tree0->SetBranchAddress("spt_flag_single_shower", &tagger_info.spt_flag_single_shower);
  tree0->SetBranchAddress("spt_energy", &tagger_info.spt_energy);
  tree0->SetBranchAddress("spt_shower_main_length", &tagger_info.spt_shower_main_length);
  tree0->SetBranchAddress("spt_shower_total_length", &tagger_info.spt_shower_total_length);
  tree0->SetBranchAddress("spt_angle_beam", &tagger_info.spt_angle_beam);
  tree0->SetBranchAddress("spt_angle_vertical", &tagger_info.spt_angle_vertical);
  tree0->SetBranchAddress("spt_max_dQ_dx", &tagger_info.spt_max_dQ_dx);
  tree0->SetBranchAddress("spt_angle_beam_1", &tagger_info.spt_angle_beam_1);
  tree0->SetBranchAddress("spt_angle_drift", &tagger_info.spt_angle_drift);
  tree0->SetBranchAddress("spt_angle_drift_1", &tagger_info.spt_angle_drift_1);
  tree0->SetBranchAddress("spt_num_valid_tracks", &tagger_info.spt_num_valid_tracks);
  tree0->SetBranchAddress("spt_n_vtx_segs", &tagger_info.spt_n_vtx_segs);
  tree0->SetBranchAddress("spt_max_length", &tagger_info.spt_max_length);
  tree0->SetBranchAddress("spt_flag", &tagger_info.spt_flag);

  tree0->SetBranchAddress("mgo_energy",&tagger_info.mgo_energy);
  tree0->SetBranchAddress("mgo_max_energy",&tagger_info.mgo_max_energy);
  tree0->SetBranchAddress("mgo_total_energy",&tagger_info.mgo_total_energy);
  tree0->SetBranchAddress("mgo_n_showers",&tagger_info.mgo_n_showers);
  tree0->SetBranchAddress("mgo_max_energy_1",&tagger_info.mgo_max_energy_1);
  tree0->SetBranchAddress("mgo_max_energy_2",&tagger_info.mgo_max_energy_2);
  tree0->SetBranchAddress("mgo_total_other_energy",&tagger_info.mgo_total_other_energy);
  tree0->SetBranchAddress("mgo_n_total_showers",&tagger_info.mgo_n_total_showers);
  tree0->SetBranchAddress("mgo_total_other_energy_1",&tagger_info.mgo_total_other_energy_1);
  tree0->SetBranchAddress("mgo_flag",&tagger_info.mgo_flag);
    
  tree0->SetBranchAddress("mgt_flag_single_shower",&tagger_info.mgt_flag_single_shower);
  tree0->SetBranchAddress("mgt_max_energy",&tagger_info.mgt_max_energy);
  tree0->SetBranchAddress("mgt_energy",&tagger_info.mgt_energy);
  tree0->SetBranchAddress("mgt_total_other_energy",&tagger_info.mgt_total_other_energy);
  tree0->SetBranchAddress("mgt_max_energy_1",&tagger_info.mgt_max_energy_1);
  tree0->SetBranchAddress("mgt_e_indirect_max_energy",&tagger_info.mgt_e_indirect_max_energy);
  tree0->SetBranchAddress("mgt_e_direct_max_energy",&tagger_info.mgt_e_direct_max_energy);
  tree0->SetBranchAddress("mgt_n_direct_showers",&tagger_info.mgt_n_direct_showers);
  tree0->SetBranchAddress("mgt_e_direct_total_energy",&tagger_info.mgt_e_direct_total_energy);
  tree0->SetBranchAddress("mgt_flag_indirect_max_pio",&tagger_info.mgt_flag_indirect_max_pio);
  tree0->SetBranchAddress("mgt_e_indirect_total_energy",&tagger_info.mgt_e_indirect_total_energy);
  tree0->SetBranchAddress("mgt_flag",&tagger_info.mgt_flag);
    
  tree0->SetBranchAddress("sig_1_v_angle",&tagger_info.sig_1_v_angle);
  tree0->SetBranchAddress("sig_1_v_flag_single_shower",&tagger_info.sig_1_v_flag_single_shower);
  tree0->SetBranchAddress("sig_1_v_energy",&tagger_info.sig_1_v_energy);
  tree0->SetBranchAddress("sig_1_v_energy_1",&tagger_info.sig_1_v_energy_1);
  tree0->SetBranchAddress("sig_1_v_flag",&tagger_info.sig_1_v_flag);
  
  tree0->SetBranchAddress("sig_2_v_energy",&tagger_info.sig_2_v_energy);
  tree0->SetBranchAddress("sig_2_v_shower_angle",&tagger_info.sig_2_v_shower_angle);
  tree0->SetBranchAddress("sig_2_v_flag_single_shower",&tagger_info.sig_2_v_flag_single_shower);
  tree0->SetBranchAddress("sig_2_v_medium_dQ_dx",&tagger_info.sig_2_v_medium_dQ_dx);
  tree0->SetBranchAddress("sig_2_v_start_dQ_dx",&tagger_info.sig_2_v_start_dQ_dx);
  tree0->SetBranchAddress("sig_2_v_flag",&tagger_info.sig_2_v_flag);
  
  tree0->SetBranchAddress("sig_flag",&tagger_info.sig_flag);

  tree0->SetBranchAddress("tro_1_v_particle_type",&tagger_info.tro_1_v_particle_type);
  tree0->SetBranchAddress("tro_1_v_flag_dir_weak",&tagger_info.tro_1_v_flag_dir_weak);
  tree0->SetBranchAddress("tro_1_v_min_dis",&tagger_info.tro_1_v_min_dis);
  tree0->SetBranchAddress("tro_1_v_sg1_length",&tagger_info.tro_1_v_sg1_length);
  tree0->SetBranchAddress("tro_1_v_shower_main_length",&tagger_info.tro_1_v_shower_main_length);
  tree0->SetBranchAddress("tro_1_v_max_n_vtx_segs",&tagger_info.tro_1_v_max_n_vtx_segs);
  tree0->SetBranchAddress("tro_1_v_tmp_length",&tagger_info.tro_1_v_tmp_length);
  tree0->SetBranchAddress("tro_1_v_medium_dQ_dx",&tagger_info.tro_1_v_medium_dQ_dx);
  tree0->SetBranchAddress("tro_1_v_dQ_dx_cut",&tagger_info.tro_1_v_dQ_dx_cut);
  tree0->SetBranchAddress("tro_1_v_flag_shower_topology",&tagger_info.tro_1_v_flag_shower_topology);
  tree0->SetBranchAddress("tro_1_v_flag",&tagger_info.tro_1_v_flag);
  
  tree0->SetBranchAddress("tro_2_v_energy",&tagger_info.tro_2_v_energy);
  tree0->SetBranchAddress("tro_2_v_stem_length",&tagger_info.tro_2_v_stem_length);
  tree0->SetBranchAddress("tro_2_v_iso_angle",&tagger_info.tro_2_v_iso_angle);
  tree0->SetBranchAddress("tro_2_v_max_length",&tagger_info.tro_2_v_max_length);
  tree0->SetBranchAddress("tro_2_v_angle",&tagger_info.tro_2_v_angle);
  tree0->SetBranchAddress("tro_2_v_flag",&tagger_info.tro_2_v_flag);
  
  tree0->SetBranchAddress("tro_3_stem_length",&tagger_info.tro_3_stem_length);
  tree0->SetBranchAddress("tro_3_n_muon_segs",&tagger_info.tro_3_n_muon_segs);
  tree0->SetBranchAddress("tro_3_energy",&tagger_info.tro_3_energy);
  tree0->SetBranchAddress("tro_3_flag",&tagger_info.tro_3_flag);
  
  tree0->SetBranchAddress("tro_4_v_dir2_mag",&tagger_info.tro_4_v_dir2_mag);
  tree0->SetBranchAddress("tro_4_v_angle",&tagger_info.tro_4_v_angle);
  tree0->SetBranchAddress("tro_4_v_angle1",&tagger_info.tro_4_v_angle1);
  tree0->SetBranchAddress("tro_4_v_angle2",&tagger_info.tro_4_v_angle2);
  tree0->SetBranchAddress("tro_4_v_length",&tagger_info.tro_4_v_length);
  tree0->SetBranchAddress("tro_4_v_length1",&tagger_info.tro_4_v_length1);
  tree0->SetBranchAddress("tro_4_v_medium_dQ_dx",&tagger_info.tro_4_v_medium_dQ_dx);
  tree0->SetBranchAddress("tro_4_v_end_dQ_dx",&tagger_info.tro_4_v_end_dQ_dx);
  tree0->SetBranchAddress("tro_4_v_energy",&tagger_info.tro_4_v_energy);
  tree0->SetBranchAddress("tro_4_v_shower_main_length",&tagger_info.tro_4_v_shower_main_length);
  tree0->SetBranchAddress("tro_4_v_flag_shower_trajectory",&tagger_info.tro_4_v_flag_shower_trajectory);
  tree0->SetBranchAddress("tro_4_v_flag",&tagger_info.tro_4_v_flag);
  
  tree0->SetBranchAddress("tro_5_v_max_angle",&tagger_info.tro_5_v_max_angle);
  tree0->SetBranchAddress("tro_5_v_min_angle",&tagger_info.tro_5_v_min_angle);
  tree0->SetBranchAddress("tro_5_v_max_length",&tagger_info.tro_5_v_max_length);
  tree0->SetBranchAddress("tro_5_v_iso_angle",&tagger_info.tro_5_v_iso_angle);
  tree0->SetBranchAddress("tro_5_v_n_vtx_segs",&tagger_info.tro_5_v_n_vtx_segs);
  tree0->SetBranchAddress("tro_5_v_min_count",&tagger_info.tro_5_v_min_count);
  tree0->SetBranchAddress("tro_5_v_max_count",&tagger_info.tro_5_v_max_count);
  tree0->SetBranchAddress("tro_5_v_energy",&tagger_info.tro_5_v_energy);
  tree0->SetBranchAddress("tro_5_v_flag",&tagger_info.tro_5_v_flag);
  
  tree0->SetBranchAddress("tro_flag",&tagger_info.tro_flag);


  // cosmic tagger ...
  
  tree0->SetBranchAddress("cosmict_flag_1",&tagger_info.cosmict_flag_1);
  tree0->SetBranchAddress("cosmict_flag_2",&tagger_info.cosmict_flag_2);
  tree0->SetBranchAddress("cosmict_flag_3",&tagger_info.cosmict_flag_3);
  tree0->SetBranchAddress("cosmict_flag_4",&tagger_info.cosmict_flag_4);
  tree0->SetBranchAddress("cosmict_flag_5",&tagger_info.cosmict_flag_5);
  tree0->SetBranchAddress("cosmict_flag_6",&tagger_info.cosmict_flag_6);
  tree0->SetBranchAddress("cosmict_flag_7",&tagger_info.cosmict_flag_7);
  tree0->SetBranchAddress("cosmict_flag_8",&tagger_info.cosmict_flag_8);
  tree0->SetBranchAddress("cosmict_flag_9",&tagger_info.cosmict_flag_9);
  tree0->SetBranchAddress("cosmict_flag_10",&tagger_info.cosmict_flag_10);
  tree0->SetBranchAddress("cosmict_flag",&tagger_info.cosmict_flag);

  tree0->SetBranchAddress("cosmict_2_filled",&tagger_info.cosmict_2_filled);
  tree0->SetBranchAddress("cosmict_2_particle_type",&tagger_info.cosmict_2_particle_type);
  tree0->SetBranchAddress("cosmict_2_n_muon_tracks",&tagger_info.cosmict_2_n_muon_tracks);
  tree0->SetBranchAddress("cosmict_2_total_shower_length",&tagger_info.cosmict_2_total_shower_length);
  tree0->SetBranchAddress("cosmict_2_flag_inside",&tagger_info.cosmict_2_flag_inside);
  tree0->SetBranchAddress("cosmict_2_angle_beam",&tagger_info.cosmict_2_angle_beam);
  tree0->SetBranchAddress("cosmict_2_flag_dir_weak",&tagger_info.cosmict_2_flag_dir_weak);
  tree0->SetBranchAddress("cosmict_2_dQ_dx_end",&tagger_info.cosmict_2_dQ_dx_end);
  tree0->SetBranchAddress("cosmict_2_dQ_dx_front",&tagger_info.cosmict_2_dQ_dx_front);
  tree0->SetBranchAddress("cosmict_2_theta",&tagger_info.cosmict_2_theta);
  tree0->SetBranchAddress("cosmict_2_phi",&tagger_info.cosmict_2_phi);
  tree0->SetBranchAddress("cosmict_2_valid_tracks",&tagger_info.cosmict_2_valid_tracks);

  tree0->SetBranchAddress("cosmict_3_filled",&tagger_info.cosmict_3_filled);
  tree0->SetBranchAddress("cosmict_3_flag_inside",&tagger_info.cosmict_3_flag_inside);
  tree0->SetBranchAddress("cosmict_3_angle_beam",&tagger_info.cosmict_3_angle_beam);
  tree0->SetBranchAddress("cosmict_3_flag_dir_weak",&tagger_info.cosmict_3_flag_dir_weak);
  tree0->SetBranchAddress("cosmict_3_dQ_dx_end",&tagger_info.cosmict_3_dQ_dx_end);
  tree0->SetBranchAddress("cosmict_3_dQ_dx_front",&tagger_info.cosmict_3_dQ_dx_front);
  tree0->SetBranchAddress("cosmict_3_theta",&tagger_info.cosmict_3_theta);
  tree0->SetBranchAddress("cosmict_3_phi",&tagger_info.cosmict_3_phi);
  tree0->SetBranchAddress("cosmict_3_valid_tracks",&tagger_info.cosmict_3_valid_tracks);
    
  tree0->SetBranchAddress("cosmict_4_filled",&tagger_info.cosmict_4_filled);
  tree0->SetBranchAddress("cosmict_4_flag_inside",&tagger_info.cosmict_4_flag_inside);
  tree0->SetBranchAddress("cosmict_4_angle_beam",&tagger_info.cosmict_4_angle_beam);
  tree0->SetBranchAddress("cosmict_4_connected_showers",&tagger_info.cosmict_4_connected_showers);

  tree0->SetBranchAddress("cosmict_5_filled",&tagger_info.cosmict_5_filled);
  tree0->SetBranchAddress("cosmict_5_flag_inside",&tagger_info.cosmict_5_flag_inside);
  tree0->SetBranchAddress("cosmict_5_angle_beam",&tagger_info.cosmict_5_angle_beam);
  tree0->SetBranchAddress("cosmict_5_connected_showers",&tagger_info.cosmict_5_connected_showers);

  tree0->SetBranchAddress("cosmict_6_filled",&tagger_info.cosmict_6_filled);
  tree0->SetBranchAddress("cosmict_6_flag_dir_weak",&tagger_info.cosmict_6_flag_dir_weak);
  tree0->SetBranchAddress("cosmict_6_flag_inside",&tagger_info.cosmict_6_flag_inside);
  tree0->SetBranchAddress("cosmict_6_angle",&tagger_info.cosmict_6_angle);
    

  tree0->SetBranchAddress("cosmict_7_filled",&tagger_info.cosmict_7_filled);
  tree0->SetBranchAddress("cosmict_7_flag_sec",&tagger_info.cosmict_7_flag_sec);
  tree0->SetBranchAddress("cosmict_7_n_muon_tracks",&tagger_info.cosmict_7_n_muon_tracks);
  tree0->SetBranchAddress("cosmict_7_total_shower_length",&tagger_info.cosmict_7_total_shower_length);
  tree0->SetBranchAddress("cosmict_7_flag_inside",&tagger_info.cosmict_7_flag_inside);
  tree0->SetBranchAddress("cosmict_7_angle_beam",&tagger_info.cosmict_7_angle_beam);
  tree0->SetBranchAddress("cosmict_7_flag_dir_weak",&tagger_info.cosmict_7_flag_dir_weak);
  tree0->SetBranchAddress("cosmict_7_dQ_dx_end",&tagger_info.cosmict_7_dQ_dx_end);
  tree0->SetBranchAddress("cosmict_7_dQ_dx_front",&tagger_info.cosmict_7_dQ_dx_front);
  tree0->SetBranchAddress("cosmict_7_theta",&tagger_info.cosmict_7_theta);
  tree0->SetBranchAddress("cosmict_7_phi",&tagger_info.cosmict_7_phi);

  tree0->SetBranchAddress("cosmict_8_filled",&tagger_info.cosmict_8_filled);
  tree0->SetBranchAddress("cosmict_8_flag_out",&tagger_info.cosmict_8_flag_out);
  tree0->SetBranchAddress("cosmict_8_muon_length",&tagger_info.cosmict_8_muon_length);
  tree0->SetBranchAddress("cosmict_8_acc_length",&tagger_info.cosmict_8_acc_length);

  tree0->SetBranchAddress("cosmict_10_flag_inside",&tagger_info.cosmict_10_flag_inside);
  tree0->SetBranchAddress("cosmict_10_vtx_z",&tagger_info.cosmict_10_vtx_z);
  tree0->SetBranchAddress("cosmict_10_flag_shower",&tagger_info.cosmict_10_flag_shower);
  tree0->SetBranchAddress("cosmict_10_flag_dir_weak",&tagger_info.cosmict_10_flag_dir_weak);
  tree0->SetBranchAddress("cosmict_10_angle_beam",&tagger_info.cosmict_10_angle_beam);
  tree0->SetBranchAddress("cosmict_10_length",&tagger_info.cosmict_10_length);
  
  tree0->SetBranchAddress("numu_cc_flag",&tagger_info.numu_cc_flag);
  
  tree0->SetBranchAddress("numu_cc_flag_1",&tagger_info.numu_cc_flag_1);
  tree0->SetBranchAddress("numu_cc_1_particle_type",&tagger_info.numu_cc_1_particle_type);
  tree0->SetBranchAddress("numu_cc_1_length",&tagger_info.numu_cc_1_length);
  tree0->SetBranchAddress("numu_cc_1_medium_dQ_dx",&tagger_info.numu_cc_1_medium_dQ_dx);
  tree0->SetBranchAddress("numu_cc_1_dQ_dx_cut",&tagger_info.numu_cc_1_dQ_dx_cut);
  tree0->SetBranchAddress("numu_cc_1_direct_length",&tagger_info.numu_cc_1_direct_length);
  tree0->SetBranchAddress("numu_cc_1_n_daughter_tracks",&tagger_info.numu_cc_1_n_daughter_tracks);
  tree0->SetBranchAddress("numu_cc_1_n_daughter_all",&tagger_info.numu_cc_1_n_daughter_all);
  
  tree0->SetBranchAddress("numu_cc_flag_2",&tagger_info.numu_cc_flag_2);
  tree0->SetBranchAddress("numu_cc_2_length",&tagger_info.numu_cc_2_length);
  tree0->SetBranchAddress("numu_cc_2_total_length",&tagger_info.numu_cc_2_total_length);
  tree0->SetBranchAddress("numu_cc_2_n_daughter_tracks",&tagger_info.numu_cc_2_n_daughter_tracks);
  tree0->SetBranchAddress("numu_cc_2_n_daughter_all",&tagger_info.numu_cc_2_n_daughter_all);
  
  tree0->SetBranchAddress("numu_cc_flag_3",&tagger_info.numu_cc_flag_3);
  tree0->SetBranchAddress("numu_cc_3_particle_type",&tagger_info.numu_cc_3_particle_type);
  tree0->SetBranchAddress("numu_cc_3_max_length",&tagger_info.numu_cc_3_max_length);
  tree0->SetBranchAddress("numu_cc_3_track_length",&tagger_info.numu_cc_3_acc_track_length);
  tree0->SetBranchAddress("numu_cc_3_max_length_all",&tagger_info.numu_cc_3_max_length_all);
  tree0->SetBranchAddress("numu_cc_3_max_muon_length",&tagger_info.numu_cc_3_max_muon_length);
  tree0->SetBranchAddress("numu_cc_3_n_daughter_tracks",&tagger_info.numu_cc_3_n_daughter_tracks);
  tree0->SetBranchAddress("numu_cc_3_n_daughter_all",&tagger_info.numu_cc_3_n_daughter_all);

  // numu BDTs
  tree0->SetBranchAddress("cosmict_2_4_score",&tagger_info.cosmict_2_4_score);
  tree0->SetBranchAddress("cosmict_3_5_score",&tagger_info.cosmict_3_5_score);
  tree0->SetBranchAddress("cosmict_6_score",&tagger_info.cosmict_6_score);
  tree0->SetBranchAddress("cosmict_7_score",&tagger_info.cosmict_7_score);
  tree0->SetBranchAddress("cosmict_8_score",&tagger_info.cosmict_8_score);
  tree0->SetBranchAddress("cosmict_10_score",&tagger_info.cosmict_10_score);

  tree0->SetBranchAddress("numu_1_score",&tagger_info.numu_1_score);
  tree0->SetBranchAddress("numu_2_score",&tagger_info.numu_2_score);
  tree0->SetBranchAddress("numu_3_score",&tagger_info.numu_3_score);

  tree0->SetBranchAddress("cosmict_score",&tagger_info.cosmict_score);
  tree0->SetBranchAddress("numu_score",&tagger_info.numu_score);
    
    
  // BDTs ...
  tree0->SetBranchAddress("mipid_score",&tagger_info.mipid_score);
  tree0->SetBranchAddress("gap_score",&tagger_info.gap_score);
  tree0->SetBranchAddress("hol_lol_score",&tagger_info.hol_lol_score);
  tree0->SetBranchAddress("cme_anc_score",&tagger_info.cme_anc_score);
  tree0->SetBranchAddress("mgo_mgt_score",&tagger_info.mgo_mgt_score);
  tree0->SetBranchAddress("br1_score",&tagger_info.br1_score);
  
  tree0->SetBranchAddress("br3_score",&tagger_info.br3_score);
  tree0->SetBranchAddress("br3_3_score",&tagger_info.br3_3_score);
  tree0->SetBranchAddress("br3_5_score",&tagger_info.br3_5_score);
  tree0->SetBranchAddress("br3_6_score",&tagger_info.br3_6_score);
  tree0->SetBranchAddress("stemdir_br2_score",&tagger_info.stemdir_br2_score);
  tree0->SetBranchAddress("trimuon_score",&tagger_info.trimuon_score);
  
  tree0->SetBranchAddress("br4_tro_score",&tagger_info.br4_tro_score);
  tree0->SetBranchAddress("mipquality_score",&tagger_info.mipquality_score);
  tree0->SetBranchAddress("pio_1_score",&tagger_info.pio_1_score);
  tree0->SetBranchAddress("pio_2_score",&tagger_info.pio_2_score);
  tree0->SetBranchAddress("stw_spt_score",&tagger_info.stw_spt_score);
  tree0->SetBranchAddress("vis_1_score",&tagger_info.vis_1_score);
  
  tree0->SetBranchAddress("vis_2_score",&tagger_info.vis_2_score);
  tree0->SetBranchAddress("stw_2_score",&tagger_info.stw_2_score);
  tree0->SetBranchAddress("stw_3_score",&tagger_info.stw_3_score);
  tree0->SetBranchAddress("stw_4_score",&tagger_info.stw_4_score);
  tree0->SetBranchAddress("sig_1_score",&tagger_info.sig_1_score);
  tree0->SetBranchAddress("sig_2_score",&tagger_info.sig_2_score);
  
  tree0->SetBranchAddress("lol_1_score",&tagger_info.lol_1_score);
  tree0->SetBranchAddress("lol_2_score",&tagger_info.lol_2_score);
  tree0->SetBranchAddress("tro_1_score",&tagger_info.tro_1_score);
  tree0->SetBranchAddress("tro_2_score",&tagger_info.tro_2_score);
  tree0->SetBranchAddress("tro_4_score",&tagger_info.tro_4_score);
  tree0->SetBranchAddress("tro_5_score",&tagger_info.tro_5_score);
  tree0->SetBranchAddress("nue_score",&tagger_info.nue_score);

  tree0->SetBranchAddress("run",&tagger_info.run);
  tree0->SetBranchAddress("subrun",&tagger_info.subrun);
  tree0->SetBranchAddress("event",&tagger_info.event);
  
  tree0->SetBranchAddress("nuvtx_diff",&tagger_info.nuvtx_diff);
  tree0->SetBranchAddress("showervtx_diff",&tagger_info.showervtx_diff);
  tree0->SetBranchAddress("muonvtx_diff",&tagger_info.muonvtx_diff);
  
  if (flag==1) tree0->SetBranchAddress("match_isFC",&tagger_info.match_isFC);
  if (flag==1) tree0->SetBranchAddress("truth_isCC",&tagger_info.truth_isCC);
  if (flag==1) tree0->SetBranchAddress("truth_vtxInside",&tagger_info.truth_vtxInside);
  if (flag==1) tree0->SetBranchAddress("truth_nuPdg",&tagger_info.truth_nuPdg);
  if (flag==1) tree0->SetBranchAddress("truth_nuEnergy",&tagger_info.truth_nuEnergy);
  if (flag==1) tree0->SetBranchAddress("truth_nuIntType",&tagger_info.truth_nuIntType);
  if (flag==1) tree0->SetBranchAddress("truth_energyInside",&tagger_info.truth_energyInside);
  if (flag==1) tree0->SetBranchAddress("weight_spline",&tagger_info.weight_spline);
  if (flag==1) tree0->SetBranchAddress("weight_cv",&tagger_info.weight_cv);
  if (flag==1) tree0->SetBranchAddress("weight_lee",&tagger_info.weight_lee);
  if (flag==1) tree0->SetBranchAddress("kine_reco_Enu",&tagger_info.kine_reco_Enu);
  
  if (flag==1) tree0->SetBranchAddress("kine_reco_add_energy",&tagger_info.kine_reco_add_energy);
  if (flag==1) tree0->SetBranchAddress("kine_pio_mass",&tagger_info.kine_pio_mass);
  if (flag==1) tree0->SetBranchAddress("kine_pio_vtx_dis",&tagger_info.kine_pio_vtx_dis);
  if (flag==1) tree0->SetBranchAddress("kine_pio_flag",&tagger_info.kine_pio_flag);
  if (flag==1) tree0->SetBranchAddress("kine_pio_energy_1",&tagger_info.kine_pio_energy_1);
  if (flag==1) tree0->SetBranchAddress("kine_pio_theta_1",&tagger_info.kine_pio_theta_1);
  if (flag==1) tree0->SetBranchAddress("kine_pio_phi_1",&tagger_info.kine_pio_phi_1);
  if (flag==1) tree0->SetBranchAddress("kine_pio_dis_1",&tagger_info.kine_pio_dis_1);
  if (flag==1) tree0->SetBranchAddress("kine_pio_energy_2",&tagger_info.kine_pio_energy_2);
  if (flag==1) tree0->SetBranchAddress("kine_pio_theta_2",&tagger_info.kine_pio_theta_2);
  if (flag==1) tree0->SetBranchAddress("kine_pio_phi_2",&tagger_info.kine_pio_phi_2);
  if (flag==1) tree0->SetBranchAddress("kine_pio_dis_2",&tagger_info.kine_pio_dis_2);
  if (flag==1) tree0->SetBranchAddress("kine_pio_angle",&tagger_info.kine_pio_angle);
  
  if (flag==1) tree0->SetBranchAddress("event_type",&tagger_info.event_type);
}

void put_tree_address(TTree *T_tagger, TaggerInfo& tagger_info, int flag){
   // cosmic tagger
    T_tagger->Branch("cosmic_flag", &tagger_info.cosmic_flag, "cosmic_flag/F");
    T_tagger->Branch("cosmic_n_solid_tracks",&tagger_info.cosmic_n_solid_tracks,"cosmic_n_solid_tracks/F");
    T_tagger->Branch("cosmic_energy_main_showers",&tagger_info.cosmic_energy_main_showers,"cosmic_energy_main_showers/F");
    T_tagger->Branch("cosmic_energy_direct_showers",&tagger_info.cosmic_energy_direct_showers,"cosmic_energy_direct_showers/F");
    T_tagger->Branch("cosmic_energy_indirect_showers",&tagger_info.cosmic_energy_indirect_showers,"cosmic_energy_indirect_showers/F");
    T_tagger->Branch("cosmic_n_direct_showers",&tagger_info.cosmic_n_direct_showers,"cosmic_n_direct_showers/F");
    T_tagger->Branch("cosmic_n_indirect_showers",&tagger_info.cosmic_n_indirect_showers,"cosmic_n_indirect_showers/F");
    T_tagger->Branch("cosmic_n_main_showers",&tagger_info.cosmic_n_main_showers,"cosmic_n_main_showers/F");
    T_tagger->Branch("cosmic_filled",&tagger_info.cosmic_filled,"cosmic_filled/F");
    
    // gap tagger
    T_tagger->Branch("gap_flag",&tagger_info.gap_flag,"gap_flag/F");
    T_tagger->Branch("gap_flag_prolong_u",&tagger_info.gap_flag_prolong_u,"gap_flag_prolong_u/F");
    T_tagger->Branch("gap_flag_prolong_v",&tagger_info.gap_flag_prolong_v,"gap_flag_prolong_v/F");
    T_tagger->Branch("gap_flag_prolong_w",&tagger_info.gap_flag_prolong_w,"gap_flag_prolong_w/F");
    T_tagger->Branch("gap_flag_parallel",&tagger_info.gap_flag_parallel,"gap_flag_parallel/F");
    T_tagger->Branch("gap_n_points",&tagger_info.gap_n_points,"gap_n_points/F");
    T_tagger->Branch("gap_n_bad",&tagger_info.gap_n_bad,"gap_n_bad/F");
    T_tagger->Branch("gap_energy",&tagger_info.gap_energy,"gap_energy/F");
    T_tagger->Branch("gap_num_valid_tracks",&tagger_info.gap_num_valid_tracks,"gap_num_valid_tracks/F");
    T_tagger->Branch("gap_flag_single_shower",&tagger_info.gap_flag_single_shower,"gap_flag_single_shower/F");
    T_tagger->Branch("gap_filled",&tagger_info.gap_filled,"gap_filled/F");

    // mip quality
    T_tagger->Branch("mip_quality_flag",&tagger_info.mip_quality_flag,"mip_quality_flag/F");
    T_tagger->Branch("mip_quality_energy",&tagger_info.mip_quality_energy,"mip_quality_energy/F");
    T_tagger->Branch("mip_quality_overlap",&tagger_info.mip_quality_overlap,"mip_quality_overlap/F");
    T_tagger->Branch("mip_quality_n_showers",&tagger_info.mip_quality_n_showers,"mip_quality_n_showers/F");
    T_tagger->Branch("mip_quality_n_tracks",&tagger_info.mip_quality_n_tracks,"mip_quality_n_tracks/F");
    T_tagger->Branch("mip_quality_flag_inside_pi0",&tagger_info.mip_quality_flag_inside_pi0,"mip_quality_flag_inside_pi0/F");
    T_tagger->Branch("mip_quality_n_pi0_showers",&tagger_info.mip_quality_n_pi0_showers,"mip_quality_n_pi0_showers/F");
    T_tagger->Branch("mip_quality_shortest_length",&tagger_info.mip_quality_shortest_length,"mip_quality_shortest_length/F");
    T_tagger->Branch("mip_quality_acc_length",&tagger_info.mip_quality_acc_length,"mip_quality_acc_length/F");
    T_tagger->Branch("mip_quality_shortest_angle",&tagger_info.mip_quality_shortest_angle,"mip_quality_shortest_angle/F");
    T_tagger->Branch("mip_quality_flag_proton",&tagger_info.mip_quality_flag_proton,"mip_quality_flag_proton/F");
    T_tagger->Branch("mip_quality_filled",&tagger_info.mip_quality_filled,"mip_quality_filled/F");

    // mip
    T_tagger->Branch("mip_flag",&tagger_info.mip_flag,"mip_flag/F");
    T_tagger->Branch("mip_energy",&tagger_info.mip_energy,"mip_energy/F");
    T_tagger->Branch("mip_n_end_reduction",&tagger_info.mip_n_end_reduction,"mip_n_end_reduction/F");
    T_tagger->Branch("mip_n_first_mip",&tagger_info.mip_n_first_mip,"mip_n_first_mip/F");
    T_tagger->Branch("mip_n_first_non_mip",&tagger_info.mip_n_first_non_mip,"mip_n_first_non_mip/F");
    T_tagger->Branch("mip_n_first_non_mip_1",&tagger_info.mip_n_first_non_mip_1,"mip_n_first_non_mip_1/F");
    T_tagger->Branch("mip_n_first_non_mip_2",&tagger_info.mip_n_first_non_mip_2,"mip_n_first_non_mip_2/F");

    T_tagger->Branch("mip_vec_dQ_dx_0",&tagger_info.mip_vec_dQ_dx_0,"mip_vec_dQ_dx_0/F");
    T_tagger->Branch("mip_vec_dQ_dx_1",&tagger_info.mip_vec_dQ_dx_1,"mip_vec_dQ_dx_1/F");
    T_tagger->Branch("mip_vec_dQ_dx_2",&tagger_info.mip_vec_dQ_dx_2,"mip_vec_dQ_dx_2/F");
    T_tagger->Branch("mip_vec_dQ_dx_3",&tagger_info.mip_vec_dQ_dx_3,"mip_vec_dQ_dx_3/F");
    T_tagger->Branch("mip_vec_dQ_dx_4",&tagger_info.mip_vec_dQ_dx_4,"mip_vec_dQ_dx_4/F");
    T_tagger->Branch("mip_vec_dQ_dx_5",&tagger_info.mip_vec_dQ_dx_5,"mip_vec_dQ_dx_5/F");
    T_tagger->Branch("mip_vec_dQ_dx_6",&tagger_info.mip_vec_dQ_dx_6,"mip_vec_dQ_dx_6/F");
    T_tagger->Branch("mip_vec_dQ_dx_7",&tagger_info.mip_vec_dQ_dx_7,"mip_vec_dQ_dx_7/F");
    T_tagger->Branch("mip_vec_dQ_dx_8",&tagger_info.mip_vec_dQ_dx_8,"mip_vec_dQ_dx_8/F");
    T_tagger->Branch("mip_vec_dQ_dx_9",&tagger_info.mip_vec_dQ_dx_9,"mip_vec_dQ_dx_9/F");
    T_tagger->Branch("mip_vec_dQ_dx_10",&tagger_info.mip_vec_dQ_dx_10,"mip_vec_dQ_dx_10/F");
    T_tagger->Branch("mip_vec_dQ_dx_11",&tagger_info.mip_vec_dQ_dx_11,"mip_vec_dQ_dx_11/F");
    T_tagger->Branch("mip_vec_dQ_dx_12",&tagger_info.mip_vec_dQ_dx_12,"mip_vec_dQ_dx_12/F");
    T_tagger->Branch("mip_vec_dQ_dx_13",&tagger_info.mip_vec_dQ_dx_13,"mip_vec_dQ_dx_13/F");
    T_tagger->Branch("mip_vec_dQ_dx_14",&tagger_info.mip_vec_dQ_dx_14,"mip_vec_dQ_dx_14/F");
    T_tagger->Branch("mip_vec_dQ_dx_15",&tagger_info.mip_vec_dQ_dx_15,"mip_vec_dQ_dx_15/F");
    T_tagger->Branch("mip_vec_dQ_dx_16",&tagger_info.mip_vec_dQ_dx_16,"mip_vec_dQ_dx_16/F");
    T_tagger->Branch("mip_vec_dQ_dx_17",&tagger_info.mip_vec_dQ_dx_17,"mip_vec_dQ_dx_17/F");
    T_tagger->Branch("mip_vec_dQ_dx_18",&tagger_info.mip_vec_dQ_dx_18,"mip_vec_dQ_dx_18/F");
    T_tagger->Branch("mip_vec_dQ_dx_19",&tagger_info.mip_vec_dQ_dx_19,"mip_vec_dQ_dx_19/F");

    T_tagger->Branch("mip_max_dQ_dx_sample",&tagger_info.mip_max_dQ_dx_sample,"mip_max_dQ_dx_sample/F");
    T_tagger->Branch("mip_n_below_threshold",&tagger_info.mip_n_below_threshold,"mip_n_below_threshold/F");
    T_tagger->Branch("mip_n_below_zero",&tagger_info.mip_n_below_zero,"mip_n_below_zero/F");
    T_tagger->Branch("mip_n_lowest",&tagger_info.mip_n_lowest,"mip_n_lowest/F");
    T_tagger->Branch("mip_n_highest",&tagger_info.mip_n_highest,"mip_n_highest/F");

    T_tagger->Branch("mip_lowest_dQ_dx",&tagger_info.mip_lowest_dQ_dx,"mip_lowest_dQ_dx/F");
    T_tagger->Branch("mip_highest_dQ_dx",&tagger_info.mip_highest_dQ_dx,"mip_highest_dQ_dx/F");
    T_tagger->Branch("mip_medium_dQ_dx",&tagger_info.mip_medium_dQ_dx,"mip_medium_dQ_dx/F");
    T_tagger->Branch("mip_stem_length",&tagger_info.mip_stem_length,"mip_stem_length/F");
    T_tagger->Branch("mip_length_main",&tagger_info.mip_length_main,"mip_length_main/F");
    T_tagger->Branch("mip_length_total",&tagger_info.mip_length_total,"mip_length_total/F");
    T_tagger->Branch("mip_angle_beam",&tagger_info.mip_angle_beam,"mip_angle_beam/F");
    T_tagger->Branch("mip_iso_angle",&tagger_info.mip_iso_angle,"mip_iso_angle/F");

    T_tagger->Branch("mip_n_vertex",&tagger_info.mip_n_vertex,"mip_n_vertex/F");
    T_tagger->Branch("mip_n_good_tracks",&tagger_info.mip_n_good_tracks,"mip_n_good_tracks/F");
    T_tagger->Branch("mip_E_indirect_max_energy",&tagger_info.mip_E_indirect_max_energy,"mip_E_indirect_max_energy/F");
    T_tagger->Branch("mip_flag_all_above",&tagger_info.mip_flag_all_above,"mip_flag_all_above/F");
    T_tagger->Branch("mip_min_dQ_dx_5",&tagger_info.mip_min_dQ_dx_5,"mip_min_dQ_dx_5/F");
    T_tagger->Branch("mip_n_other_vertex",&tagger_info.mip_n_other_vertex,"mip_n_other_vertex/F");
    T_tagger->Branch("mip_n_stem_size",&tagger_info.mip_n_stem_size,"mip_n_stem_size/F");
    T_tagger->Branch("mip_flag_stem_trajectory",&tagger_info.mip_flag_stem_trajectory,"mip_flag_stem_trajectory/F");
    T_tagger->Branch("mip_min_dis",&tagger_info.mip_min_dis,"mip_min_dis/F");
    T_tagger->Branch("mip_filled",&tagger_info.mip_filled,"mip_filled/F");

    // pio ...
    T_tagger->Branch("pio_flag",&tagger_info.pio_flag,"pio_flag/F");
    T_tagger->Branch("pio_mip_id",&tagger_info.pio_mip_id,"pio_mip_id/F");
    T_tagger->Branch("pio_filled",&tagger_info.pio_filled,"pio_filled/F");
    T_tagger->Branch("pio_flag_pio",&tagger_info.pio_flag_pio,"pio_flag_pio/F");
    
    T_tagger->Branch("pio_1_flag",&tagger_info.pio_1_flag,"pio_1_flag/F");
    T_tagger->Branch("pio_1_mass",&tagger_info.pio_1_mass,"pio_1_mass/F");
    T_tagger->Branch("pio_1_pio_type",&tagger_info.pio_1_pio_type,"pio_1_pio_type/F");
    T_tagger->Branch("pio_1_energy_1",&tagger_info.pio_1_energy_1,"pio_1_energy_1/F");
    T_tagger->Branch("pio_1_energy_2",&tagger_info.pio_1_energy_2,"pio_1_energy_2/F");
    T_tagger->Branch("pio_1_dis_1",&tagger_info.pio_1_dis_1,"pio_1_dis_1/F");
    T_tagger->Branch("pio_1_dis_2",&tagger_info.pio_1_dis_2,"pio_1_dis_2/F");
    
    T_tagger->Branch("pio_2_v_dis2",&tagger_info.pio_2_v_dis2);
    T_tagger->Branch("pio_2_v_angle2",&tagger_info.pio_2_v_angle2);
    T_tagger->Branch("pio_2_v_acc_length",&tagger_info.pio_2_v_acc_length);
    T_tagger->Branch("pio_2_v_flag",&tagger_info.pio_2_v_flag);

    
    
    
    // bad reconstruction ...
    T_tagger->Branch("stem_dir_flag",&tagger_info.stem_dir_flag,"stem_dir_flag/F");
    T_tagger->Branch("stem_dir_flag_single_shower",&tagger_info.stem_dir_flag_single_shower,"stem_dir_flag_single_shower/F");
    T_tagger->Branch("stem_dir_filled",&tagger_info.stem_dir_filled,"stem_dir_filled/F");
    T_tagger->Branch("stem_dir_angle",&tagger_info.stem_dir_angle,"stem_dir_angle/F");
    T_tagger->Branch("stem_dir_energy",&tagger_info.stem_dir_energy,"stem_dir_energy/F");
    T_tagger->Branch("stem_dir_angle1",&tagger_info.stem_dir_angle1,"stem_dir_angle1/F");
    T_tagger->Branch("stem_dir_angle2",&tagger_info.stem_dir_angle2,"stem_dir_angle2/F");
    T_tagger->Branch("stem_dir_angle3",&tagger_info.stem_dir_angle3,"stem_dir_angle3/F");
    T_tagger->Branch("stem_dir_ratio",&tagger_info.stem_dir_ratio,"stem_dir_ratio/F");

    T_tagger->Branch("br_filled",&tagger_info.br_filled,"br_filled/F");
    
    T_tagger->Branch("br1_flag",&tagger_info.br1_flag,"br1_flag/F");
    
    T_tagger->Branch("br1_1_flag",&tagger_info.br1_1_flag,"br1_1_flag/F");
    T_tagger->Branch("br1_1_shower_type",&tagger_info.br1_1_shower_type,"br1_1_shower_type/F");
    T_tagger->Branch("br1_1_vtx_n_segs",&tagger_info.br1_1_vtx_n_segs,"br1_1_vtx_n_segs/F");
    T_tagger->Branch("br1_1_energy",&tagger_info.br1_1_energy,"br1_1_energy/F");
    T_tagger->Branch("br1_1_n_segs",&tagger_info.br1_1_n_segs,"br1_1_n_segs/F");
    T_tagger->Branch("br1_1_flag_sg_topology",&tagger_info.br1_1_flag_sg_topology,"br1_1_flag_sg_topology/F");
    T_tagger->Branch("br1_1_flag_sg_trajectory",&tagger_info.br1_1_flag_sg_trajectory,"br1_1_flag_sg_trajectory/F");
    T_tagger->Branch("br1_1_sg_length",&tagger_info.br1_1_sg_length,"br1_1_sg_length/F");

    T_tagger->Branch("br1_2_flag",&tagger_info.br1_2_flag,"br1_2_flag/F");
    T_tagger->Branch("br1_2_energy",&tagger_info.br1_2_energy,"br1_2_energy/F");
    T_tagger->Branch("br1_2_n_connected",&tagger_info.br1_2_n_connected,"br1_2_n_connected/F");
    T_tagger->Branch("br1_2_max_length",&tagger_info.br1_2_max_length,"br1_2_max_length/F");
    T_tagger->Branch("br1_2_n_connected_1",&tagger_info.br1_2_n_connected_1,"br1_2_n_connected_1/F");
    T_tagger->Branch("br1_2_vtx_n_segs",&tagger_info.br1_2_vtx_n_segs,"br1_2_vtx_n_segs/F");
    T_tagger->Branch("br1_2_n_shower_segs",&tagger_info.br1_2_n_shower_segs,"br1_2_n_shower_segs/F");
    T_tagger->Branch("br1_2_max_length_ratio",&tagger_info.br1_2_max_length_ratio,"br1_2_max_length_ratio/F");
    T_tagger->Branch("br1_2_shower_length",&tagger_info.br1_2_shower_length,"br1_2_shower_length/F");

    T_tagger->Branch("br1_3_flag",&tagger_info.br1_3_flag,"br1_3_flag/F");
    T_tagger->Branch("br1_3_energy",&tagger_info.br1_3_energy,"br1_3_energy/F");
    T_tagger->Branch("br1_3_n_connected_p",&tagger_info.br1_3_n_connected_p,"br1_3_n_connected_p/F");
    T_tagger->Branch("br1_3_max_length_p",&tagger_info.br1_3_max_length_p,"br1_3_max_length_p/F");
    T_tagger->Branch("br1_3_n_shower_segs",&tagger_info.br1_3_n_shower_segs,"br1_3_n_shower_segs/F");
    T_tagger->Branch("br1_3_flag_sg_topology",&tagger_info.br1_3_flag_sg_topology,"br1_3_flag_sg_topology/F");
    T_tagger->Branch("br1_3_flag_sg_trajectory",&tagger_info.br1_3_flag_sg_trajectory,"br1_3_flag_sg_trajectory/F");
    T_tagger->Branch("br1_3_n_shower_main_segs",&tagger_info.br1_3_n_shower_main_segs,"br1_3_n_shower_main_segs/F");
    T_tagger->Branch("br1_3_sg_length",&tagger_info.br1_3_sg_length,"br1_3_sg_length/F");
    
    T_tagger->Branch("br2_flag",&tagger_info.br2_flag,"br2_flag/F");
    T_tagger->Branch("br2_flag_single_shower",&tagger_info.br2_flag_single_shower,"br2_flag_single_shower/F");
    T_tagger->Branch("br2_num_valid_tracks",&tagger_info.br2_num_valid_tracks,"br2_num_valid_tracks/F");
    T_tagger->Branch("br2_energy",&tagger_info.br2_energy,"br2_energy/F");
    T_tagger->Branch("br2_angle1",&tagger_info.br2_angle1,"br2_angle1/F");
    T_tagger->Branch("br2_angle2",&tagger_info.br2_angle2,"br2_angle2/F");
    T_tagger->Branch("br2_angle",&tagger_info.br2_angle,"br2_angle/F");
    T_tagger->Branch("br2_angle3",&tagger_info.br2_angle3,"br2_angle3/F");
    T_tagger->Branch("br2_n_shower_main_segs",&tagger_info.br2_n_shower_main_segs,"br2_n_shower_main_segs/F");
    T_tagger->Branch("br2_max_angle",&tagger_info.br2_max_angle,"br2_max_angle/F");
    T_tagger->Branch("br2_sg_length",&tagger_info.br2_sg_length,"br2_sg_length/F");
    T_tagger->Branch("br2_flag_sg_trajectory",&tagger_info.br2_flag_sg_trajectory,"br2_flag_sg_trajectory/F");


    T_tagger->Branch("lol_flag",&tagger_info.lol_flag,"lol_flag/F");

    T_tagger->Branch("lol_1_v_energy",&tagger_info.lol_1_v_energy);
    T_tagger->Branch("lol_1_v_vtx_n_segs",&tagger_info.lol_1_v_vtx_n_segs);
    T_tagger->Branch("lol_1_v_nseg",&tagger_info.lol_1_v_nseg);
    T_tagger->Branch("lol_1_v_angle",&tagger_info.lol_1_v_angle);
    T_tagger->Branch("lol_1_v_flag",&tagger_info.lol_1_v_flag);

    T_tagger->Branch("lol_2_v_length",&tagger_info.lol_2_v_length);
    T_tagger->Branch("lol_2_v_angle",&tagger_info.lol_2_v_angle);
    T_tagger->Branch("lol_2_v_type",&tagger_info.lol_2_v_type);
    T_tagger->Branch("lol_2_v_vtx_n_segs",&tagger_info.lol_2_v_vtx_n_segs);
    T_tagger->Branch("lol_2_v_energy",&tagger_info.lol_2_v_energy);
    T_tagger->Branch("lol_2_v_shower_main_length",&tagger_info.lol_2_v_shower_main_length);
    T_tagger->Branch("lol_2_v_flag_dir_weak",&tagger_info.lol_2_v_flag_dir_weak);
    T_tagger->Branch("lol_2_v_flag",&tagger_info.lol_2_v_flag);

    T_tagger->Branch("lol_3_angle_beam",&tagger_info.lol_3_angle_beam,"lol_3_angle_beam/F");
    T_tagger->Branch("lol_3_n_valid_tracks",&tagger_info.lol_3_n_valid_tracks,"lol_3_n_valid_tracks/F");
    T_tagger->Branch("lol_3_min_angle",&tagger_info.lol_3_min_angle,"lol_3_min_angle/F");
    T_tagger->Branch("lol_3_vtx_n_segs",&tagger_info.lol_3_vtx_n_segs,"lol_3_vtx_n_segs/F");
    T_tagger->Branch("lol_3_energy",&tagger_info.lol_3_energy,"lol_3_energy/F");
    T_tagger->Branch("lol_3_shower_main_length",&tagger_info.lol_3_shower_main_length,"lol_3_shower_main_length/F");
    T_tagger->Branch("lol_3_n_out",&tagger_info.lol_3_n_out,"lol_3_n_out/F");
    T_tagger->Branch("lol_3_n_sum",&tagger_info.lol_3_n_sum,"lol_3_n_sum/F");
    T_tagger->Branch("lol_3_flag",&tagger_info.lol_3_flag,"lol_3_flag/F");
    
    T_tagger->Branch("br3_1_energy",&tagger_info.br3_1_energy,"br3_1_energy/F");
    T_tagger->Branch("br3_1_n_shower_segments",&tagger_info.br3_1_n_shower_segments,"br3_1_n_shower_segments/F");
    T_tagger->Branch("br3_1_sg_flag_trajectory",&tagger_info.br3_1_sg_flag_trajectory,"br3_1_sg_flag_trajectory/F");
    T_tagger->Branch("br3_1_sg_direct_length",&tagger_info.br3_1_sg_direct_length,"br3_1_sg_direct_length/F");
    T_tagger->Branch("br3_1_sg_length",&tagger_info.br3_1_sg_length,"br3_1_sg_length/F");
    T_tagger->Branch("br3_1_total_main_length",&tagger_info.br3_1_total_main_length,"br3_1_total_main_length/F");
    T_tagger->Branch("br3_1_total_length",&tagger_info.br3_1_total_length,"br3_1_total_length/F");
    T_tagger->Branch("br3_1_iso_angle",&tagger_info.br3_1_iso_angle,"br3_1_iso_angle/F");
    T_tagger->Branch("br3_1_sg_flag_topology",&tagger_info.br3_1_sg_flag_topology,"br3_1_sg_flag_topology/F");
    T_tagger->Branch("br3_1_flag",&tagger_info.br3_1_flag,"br3_1_flag/F");

    T_tagger->Branch("br3_2_n_ele",&tagger_info.br3_2_n_ele,"br3_2_n_ele/F");
    T_tagger->Branch("br3_2_n_other",&tagger_info.br3_2_n_other,"br3_2_n_other/F");
    T_tagger->Branch("br3_2_energy",&tagger_info.br3_2_energy,"br3_2_energy/F");
    T_tagger->Branch("br3_2_total_main_length",&tagger_info.br3_2_total_main_length,"br3_2_total_main_length/F");
    T_tagger->Branch("br3_2_total_length",&tagger_info.br3_2_total_length,"br3_2_total_length/F");
    T_tagger->Branch("br3_2_other_fid",&tagger_info.br3_2_other_fid,"br3_2_other_fid/F");
    T_tagger->Branch("br3_2_flag",&tagger_info.br3_2_flag,"br3_2_flag/F");

    T_tagger->Branch("br3_3_v_energy",&tagger_info.br3_3_v_energy);
    T_tagger->Branch("br3_3_v_angle",&tagger_info.br3_3_v_angle);
    T_tagger->Branch("br3_3_v_dir_length",&tagger_info.br3_3_v_dir_length);
    T_tagger->Branch("br3_3_v_length",&tagger_info.br3_3_v_length);
    T_tagger->Branch("br3_3_v_flag",&tagger_info.br3_3_v_flag);

    T_tagger->Branch("br3_4_acc_length", &tagger_info.br3_4_acc_length, "br3_4_acc_length/F");
    T_tagger->Branch("br3_4_total_length", &tagger_info.br3_4_total_length, "br3_4_total_length/F");
    T_tagger->Branch("br3_4_energy", &tagger_info.br3_4_energy, "br3_4_energy/F");
    T_tagger->Branch("br3_4_flag", &tagger_info.br3_4_flag, "br3_4_flag/F");
    
    T_tagger->Branch("br3_5_v_dir_length", &tagger_info.br3_5_v_dir_length);
    T_tagger->Branch("br3_5_v_total_length", &tagger_info.br3_5_v_total_length);
    T_tagger->Branch("br3_5_v_flag_avoid_muon_check", &tagger_info.br3_5_v_flag_avoid_muon_check);
    T_tagger->Branch("br3_5_v_n_seg", &tagger_info.br3_5_v_n_seg);
    T_tagger->Branch("br3_5_v_angle", &tagger_info.br3_5_v_angle);
    T_tagger->Branch("br3_5_v_sg_length", &tagger_info.br3_5_v_sg_length);
    T_tagger->Branch("br3_5_v_energy", &tagger_info.br3_5_v_energy);
    T_tagger->Branch("br3_5_v_n_main_segs", &tagger_info.br3_5_v_n_main_segs);
    T_tagger->Branch("br3_5_v_n_segs", &tagger_info.br3_5_v_n_segs);
    T_tagger->Branch("br3_5_v_shower_main_length", &tagger_info.br3_5_v_shower_main_length);
    T_tagger->Branch("br3_5_v_shower_total_length", &tagger_info.br3_5_v_shower_total_length);
    T_tagger->Branch("br3_5_v_flag", &tagger_info.br3_5_v_flag);

    T_tagger->Branch("br3_6_v_angle",&tagger_info.br3_6_v_angle);
    T_tagger->Branch("br3_6_v_angle1",&tagger_info.br3_6_v_angle1);
    T_tagger->Branch("br3_6_v_flag_shower_trajectory",&tagger_info.br3_6_v_flag_shower_trajectory);
    T_tagger->Branch("br3_6_v_direct_length",&tagger_info.br3_6_v_direct_length);
    T_tagger->Branch("br3_6_v_length",&tagger_info.br3_6_v_length);
    T_tagger->Branch("br3_6_v_n_other_vtx_segs",&tagger_info.br3_6_v_n_other_vtx_segs);
    T_tagger->Branch("br3_6_v_energy",&tagger_info.br3_6_v_energy);
    T_tagger->Branch("br3_6_v_flag",&tagger_info.br3_6_v_flag);

    T_tagger->Branch("br3_7_energy",&tagger_info.br3_7_energy,"br3_7_energy/F");
    T_tagger->Branch("br3_7_min_angle",&tagger_info.br3_7_min_angle,"br3_7_min_angle/F");
    T_tagger->Branch("br3_7_sg_length",&tagger_info.br3_7_sg_length,"br3_7_sg_length/F");
    T_tagger->Branch("br3_7_main_length",&tagger_info.br3_7_shower_main_length,"br3_7_shower_main_length/F");
    T_tagger->Branch("br3_7_flag",&tagger_info.br3_7_flag,"br3_7_flag/F");

    T_tagger->Branch("br3_8_max_dQ_dx",&tagger_info.br3_8_max_dQ_dx,"br3_8_max_dQ_dx/F");
    T_tagger->Branch("br3_8_energy",&tagger_info.br3_8_energy,"br3_8_energy/F");
    T_tagger->Branch("br3_8_n_main_segs",&tagger_info.br3_8_n_main_segs,"br3_8_n_main_segs/F");
    T_tagger->Branch("br3_8_shower_main_length",&tagger_info.br3_8_shower_main_length,"br3_8_shower_main_length/F");
    T_tagger->Branch("br3_8_shower_length",&tagger_info.br3_8_shower_length,"br3_8_shower_length/F");
    T_tagger->Branch("br3_8_flag",&tagger_info.br3_8_flag,"br3_8_flag/F");

    T_tagger->Branch("br3_flag",&tagger_info.br3_flag,"br3_flag/F");


    T_tagger->Branch("br4_1_shower_main_length", &tagger_info.br4_1_shower_main_length, "br4_1_shower_main_length/F");
    T_tagger->Branch("br4_1_shower_total_length", &tagger_info.br4_1_shower_total_length, "br4_1_shower_total_length/F");
    T_tagger->Branch("br4_1_min_dis", &tagger_info.br4_1_min_dis, "br4_1_min_dis/F");
    T_tagger->Branch("br4_1_energy", &tagger_info.br4_1_energy, "br4_1_energy/F");
    T_tagger->Branch("br4_1_flag_avoid_muon_check", &tagger_info.br4_1_flag_avoid_muon_check, "br4_1_flag_avoid_muon_check/F");
    T_tagger->Branch("br4_1_n_vtx_segs", &tagger_info.br4_1_n_vtx_segs, "br4_1_n_vtx_segs/F");
    T_tagger->Branch("br4_1_n_main_segs", &tagger_info.br4_1_n_main_segs, "br4_1_n_main_segs/F");
    T_tagger->Branch("br4_1_flag", &tagger_info.br4_1_flag, "br4_1_flag/F");

    T_tagger->Branch("br4_2_ratio_45", &tagger_info.br4_2_ratio_45, "br4_2_ratio_45/F");
    T_tagger->Branch("br4_2_ratio_35", &tagger_info.br4_2_ratio_35, "br4_2_ratio_35/F");
    T_tagger->Branch("br4_2_ratio_25", &tagger_info.br4_2_ratio_25, "br4_2_ratio_25/F");
    T_tagger->Branch("br4_2_ratio_15", &tagger_info.br4_2_ratio_15, "br4_2_ratio_15/F");
    T_tagger->Branch("br4_2_energy",   &tagger_info.br4_2_energy, "br4_2_energy/F");
    T_tagger->Branch("br4_2_ratio1_45", &tagger_info.br4_2_ratio1_45, "br4_2_ratio1_45/F");
    T_tagger->Branch("br4_2_ratio1_35", &tagger_info.br4_2_ratio1_35, "br4_2_ratio1_35/F");
    T_tagger->Branch("br4_2_ratio1_25", &tagger_info.br4_2_ratio1_25, "br4_2_ratio1_25/F");
    T_tagger->Branch("br4_2_ratio1_15", &tagger_info.br4_2_ratio1_15, "br4_2_ratio1_15/F");
    T_tagger->Branch("br4_2_iso_angle", &tagger_info.br4_2_iso_angle, "br4_2_iso_angle/F");
    T_tagger->Branch("br4_2_iso_angle1", &tagger_info.br4_2_iso_angle1, "br4_2_iso_angle1/F");
    T_tagger->Branch("br4_2_angle", &tagger_info.br4_2_angle, "br4_2_angle/F");
    T_tagger->Branch("br4_2_flag", &tagger_info.br4_2_flag, "br4_2_flag/F");

    T_tagger->Branch("br4_flag", &tagger_info.br4_flag, "br4_flag/F");
    

    T_tagger->Branch("hol_1_n_valid_tracks", &tagger_info.hol_1_n_valid_tracks,"hol_1_n_valid_tracks/F");
    T_tagger->Branch("hol_1_min_angle", &tagger_info.hol_1_min_angle,"hol_1_min_angle/F");
    T_tagger->Branch("hol_1_energy", &tagger_info.hol_1_energy,"hol_1_energy/F");
    T_tagger->Branch("hol_1_flag_all_shower", &tagger_info.hol_1_flag_all_shower,"hol_1_flag_all_shower/F");
    T_tagger->Branch("hol_1_min_length", &tagger_info.hol_1_min_length,"hol_1_min_length/F");
    T_tagger->Branch("hol_1_flag", &tagger_info.hol_1_flag,"hol_1_flag/F");

    T_tagger->Branch("hol_2_min_angle", &tagger_info.hol_2_min_angle,"hol_2_min_angle/F");
    T_tagger->Branch("hol_2_medium_dQ_dx", &tagger_info.hol_2_medium_dQ_dx,"hol_2_medium_dQ_dx/F");
    T_tagger->Branch("hol_2_ncount", &tagger_info.hol_2_ncount,"hol_2_ncount/F");
    T_tagger->Branch("hol_2_energy", &tagger_info.hol_2_energy,"hol_2_energy/F");
    T_tagger->Branch("hol_2_flag", &tagger_info.hol_2_flag,"hol_2_flag/F");

    T_tagger->Branch("hol_flag", &tagger_info.hol_flag,"hol_flag/F");
    

    T_tagger->Branch("vis_1_filled",&tagger_info.vis_1_filled,"vis_1_filled/F");
    T_tagger->Branch("vis_1_n_vtx_segs",&tagger_info.vis_1_n_vtx_segs,"vis_1_n_vtx_segs/F");
    T_tagger->Branch("vis_1_energy",&tagger_info.vis_1_energy,"vis_1_energy/F");
    T_tagger->Branch("vis_1_num_good_tracks",&tagger_info.vis_1_num_good_tracks,"vis_1_num_good_tracks/F");
    T_tagger->Branch("vis_1_max_angle",&tagger_info.vis_1_max_angle,"vis_1_max_angle/F");
    T_tagger->Branch("vis_1_max_shower_angle",&tagger_info.vis_1_max_shower_angle,"vis_1_max_shower_angle/F");
    T_tagger->Branch("vis_1_tmp_length1",&tagger_info.vis_1_tmp_length1,"vis_1_tmp_length1/F");
    T_tagger->Branch("vis_1_tmp_length2",&tagger_info.vis_1_tmp_length2,"vis_1_tmp_length2/F");
    T_tagger->Branch("vis_1_particle_type",&tagger_info.vis_1_particle_type,"vis_1_particle_type/F");
    T_tagger->Branch("vis_1_flag",&tagger_info.vis_1_flag,"vis_1_flag/F");

    T_tagger->Branch("vis_2_filled",&tagger_info.vis_2_filled,"vis_2_filled/F");
    T_tagger->Branch("vis_2_n_vtx_segs",&tagger_info.vis_2_n_vtx_segs,"vis_2_n_vtx_segs/F");
    T_tagger->Branch("vis_2_min_angle",&tagger_info.vis_2_min_angle,"vis_2_min_angle/F");
    T_tagger->Branch("vis_2_min_weak_track",&tagger_info.vis_2_min_weak_track,"vis_2_min_weak_track/F");
    T_tagger->Branch("vis_2_angle_beam",&tagger_info.vis_2_angle_beam,"vis_2_angle_beam/F");
    T_tagger->Branch("vis_2_min_angle1",&tagger_info.vis_2_min_angle1,"vis_2_min_angle1/F");
    T_tagger->Branch("vis_2_iso_angle1",&tagger_info.vis_2_iso_angle1,"vis_2_iso_angle1/F");
    T_tagger->Branch("vis_2_min_medium_dQ_dx",&tagger_info.vis_2_min_medium_dQ_dx,"vis_2_min_medium_dQ_dx/F");
    T_tagger->Branch("vis_2_min_length",&tagger_info.vis_2_min_length,"vis_2_min_length/F");
    T_tagger->Branch("vis_2_sg_length",&tagger_info.vis_2_sg_length,"vis_2_sg_length/F");
    T_tagger->Branch("vis_2_max_angle",&tagger_info.vis_2_max_angle,"vis_2_max_angle/F");
    T_tagger->Branch("vis_2_max_weak_track",&tagger_info.vis_2_max_weak_track,"vis_2_max_weak_track/F");
    T_tagger->Branch("vis_2_flag",&tagger_info.vis_2_flag,"vis_2_flag/F");

    T_tagger->Branch("vis_flag",&tagger_info.vis_flag,"vis_flag/F");
    

    T_tagger->Branch("stem_len_energy", &tagger_info.stem_len_energy, "stem_len_energy/F");
    T_tagger->Branch("stem_len_length", &tagger_info.stem_len_length, "stem_len_length/F");
    T_tagger->Branch("stem_len_flag_avoid_muon_check", &tagger_info.stem_len_flag_avoid_muon_check, "stem_len_flag_avoid_muon_check/F");
    T_tagger->Branch("stem_len_num_daughters", &tagger_info.stem_len_num_daughters, "stem_len_num_daughters/F");
    T_tagger->Branch("stem_len_daughter_length", &tagger_info.stem_len_daughter_length, "stem_len_daughter_length/F");
    T_tagger->Branch("stem_len_flag", &tagger_info.stem_len_flag, "stem_len_flag/F");

    T_tagger->Branch("brm_n_mu_segs",&tagger_info.brm_n_mu_segs,"brm_n_mu_segs/F");
    T_tagger->Branch("brm_Ep",&tagger_info.brm_Ep,"brm_Ep/F");
    T_tagger->Branch("brm_energy",&tagger_info.brm_energy,"brm_energy/F");
    T_tagger->Branch("brm_acc_length",&tagger_info.brm_acc_length,"brm_acc_length/F");
    T_tagger->Branch("brm_shower_total_length",&tagger_info.brm_shower_total_length,"brm_shower_total_length/F");
    T_tagger->Branch("brm_connected_length",&tagger_info.brm_connected_length,"brm_connected_length/F");
    T_tagger->Branch("brm_n_size",&tagger_info.brm_n_size,"brm_n_size/F");
    T_tagger->Branch("brm_acc_direct_length",&tagger_info.brm_acc_direct_length,"brm_acc_direct_length/F");
    T_tagger->Branch("brm_n_shower_main_segs",&tagger_info.brm_n_shower_main_segs,"brm_n_shower_main_segs/F");
    T_tagger->Branch("brm_n_mu_main",&tagger_info.brm_n_mu_main,"brm_n_mu_main/F");
    T_tagger->Branch("brm_flag",&tagger_info.brm_flag,"brm_flag/F");

    T_tagger->Branch("cme_mu_energy",&tagger_info.cme_mu_energy,"cme_mu_energy/F");
    T_tagger->Branch("cme_energy",&tagger_info.cme_energy,"cme_energy/F");
    T_tagger->Branch("cme_mu_length",&tagger_info.cme_mu_length,"cme_mu_length/F");
    T_tagger->Branch("cme_length",&tagger_info.cme_length,"cme_length/F");
    T_tagger->Branch("cme_angle_beam",&tagger_info.cme_angle_beam,"cme_angle_beam/F");
    T_tagger->Branch("cme_flag",&tagger_info.cme_flag,"cme_flag/F");

    T_tagger->Branch("anc_energy",&tagger_info.anc_energy,"anc_energy/F");
    T_tagger->Branch("anc_angle",&tagger_info.anc_angle,"anc_angle/F");
    T_tagger->Branch("anc_max_angle",&tagger_info.anc_max_angle,"anc_max_angle/F");
    T_tagger->Branch("anc_max_length",&tagger_info.anc_max_length,"anc_max_length/F");
    T_tagger->Branch("anc_acc_forward_length",&tagger_info.anc_acc_forward_length,"anc_acc_forward_length/F");
    T_tagger->Branch("anc_acc_backward_length",&tagger_info.anc_acc_backward_length,"anc_acc_backward_length/F");
    T_tagger->Branch("anc_acc_forward_length1",&tagger_info.anc_acc_forward_length1,"anc_acc_forward_length1/F");
    T_tagger->Branch("anc_shower_main_length",&tagger_info.anc_shower_main_length,"anc_shower_main_length/F");
    T_tagger->Branch("anc_shower_total_length",&tagger_info.anc_shower_total_length,"anc_shower_total_length/F");
    T_tagger->Branch("anc_flag_main_outside",&tagger_info.anc_flag_main_outside,"anc_flag_main_outside/F");
    T_tagger->Branch("anc_flag",&tagger_info.anc_flag,"anc_flag/F");

    T_tagger->Branch("lem_shower_total_length",&tagger_info.lem_shower_total_length,"lem_shower_total_length/F");
    T_tagger->Branch("lem_shower_main_length",&tagger_info.lem_shower_main_length,"lem_shower_main_length/F");
    T_tagger->Branch("lem_n_3seg",&tagger_info.lem_n_3seg,"lem_n_3seg/F");
    T_tagger->Branch("lem_e_charge",&tagger_info.lem_e_charge,"lem_e_charge/F");
    T_tagger->Branch("lem_e_dQdx",&tagger_info.lem_e_dQdx,"lem_e_dQdx/F");
    T_tagger->Branch("lem_shower_num_segs",&tagger_info.lem_shower_num_segs,"lem_shower_num_segs/F");
    T_tagger->Branch("lem_shower_num_main_segs",&tagger_info.lem_shower_num_main_segs,"lem_shower_num_main_segs/F");
    T_tagger->Branch("lem_flag",&tagger_info.lem_flag,"lem_flag/F");

    T_tagger->Branch("stw_1_energy",&tagger_info.stw_1_energy,"stw_1_energy/F");
    T_tagger->Branch("stw_1_dis",&tagger_info.stw_1_dis,"stw_1_dis/F");
    T_tagger->Branch("stw_1_dQ_dx",&tagger_info.stw_1_dQ_dx,"stw_1_dQ_dx/F");
    T_tagger->Branch("stw_1_flag_single_shower",&tagger_info.stw_1_flag_single_shower,"stw_1_flag_single_shower/F");
    T_tagger->Branch("stw_1_n_pi0",&tagger_info.stw_1_n_pi0,"stw_1_n_pi0/F");
    T_tagger->Branch("stw_1_num_valid_tracks",&tagger_info.stw_1_num_valid_tracks,"stw_1_num_valid_tracks/F");
    T_tagger->Branch("stw_1_flag",&tagger_info.stw_1_flag,"stw_1_flag/F");
    
    T_tagger->Branch("stw_2_v_medium_dQ_dx", &tagger_info.stw_2_v_medium_dQ_dx);
    T_tagger->Branch("stw_2_v_energy", &tagger_info.stw_2_v_energy);
    T_tagger->Branch("stw_2_v_angle", &tagger_info.stw_2_v_angle);
    T_tagger->Branch("stw_2_v_dir_length", &tagger_info.stw_2_v_dir_length);
    T_tagger->Branch("stw_2_v_max_dQ_dx", &tagger_info.stw_2_v_max_dQ_dx);
    T_tagger->Branch("stw_2_v_flag", &tagger_info.stw_2_v_flag);

    T_tagger->Branch("stw_3_v_angle",&tagger_info.stw_3_v_angle);
    T_tagger->Branch("stw_3_v_dir_length",&tagger_info.stw_3_v_dir_length);
    T_tagger->Branch("stw_3_v_energy",&tagger_info.stw_3_v_energy);
    T_tagger->Branch("stw_3_v_medium_dQ_dx",&tagger_info.stw_3_v_medium_dQ_dx);
    T_tagger->Branch("stw_3_v_flag",&tagger_info.stw_3_v_flag);

    T_tagger->Branch("stw_4_v_angle",&tagger_info.stw_4_v_angle);
    T_tagger->Branch("stw_4_v_dis",&tagger_info.stw_4_v_dis);
    T_tagger->Branch("stw_4_v_energy",&tagger_info.stw_4_v_energy);
    T_tagger->Branch("stw_4_v_flag",&tagger_info.stw_4_v_flag);
    
    T_tagger->Branch("stw_flag", &tagger_info.stw_flag,"stw_flag/F");

    T_tagger->Branch("spt_flag_single_shower", &tagger_info.spt_flag_single_shower, "spt_flag_single_shower/F");
    T_tagger->Branch("spt_energy", &tagger_info.spt_energy, "spt_energy/F");
    T_tagger->Branch("spt_shower_main_length", &tagger_info.spt_shower_main_length, "spt_shower_main_length/F");
    T_tagger->Branch("spt_shower_total_length", &tagger_info.spt_shower_total_length, "spt_shower_total_length/F");
    T_tagger->Branch("spt_angle_beam", &tagger_info.spt_angle_beam, "spt_angle_beam/F");
    T_tagger->Branch("spt_angle_vertical", &tagger_info.spt_angle_vertical, "spt_angle_vertical/F");
    T_tagger->Branch("spt_max_dQ_dx", &tagger_info.spt_max_dQ_dx, "spt_max_dQ_dx/F");
    T_tagger->Branch("spt_angle_beam_1", &tagger_info.spt_angle_beam_1, "spt_angle_beam_1/F");
    T_tagger->Branch("spt_angle_drift", &tagger_info.spt_angle_drift, "spt_angle_drift/F");
    T_tagger->Branch("spt_angle_drift_1", &tagger_info.spt_angle_drift_1, "spt_angle_drift_1/F");
    T_tagger->Branch("spt_num_valid_tracks", &tagger_info.spt_num_valid_tracks, "spt_num_valid_tracks/F");
    T_tagger->Branch("spt_n_vtx_segs", &tagger_info.spt_n_vtx_segs, "spt_n_vtx_segs/F");
    T_tagger->Branch("spt_max_length", &tagger_info.spt_max_length, "spt_max_length/F");
    T_tagger->Branch("spt_flag", &tagger_info.spt_flag, "spt_flag/F");

    T_tagger->Branch("mgo_energy",&tagger_info.mgo_energy,"mgo_energy/F");
    T_tagger->Branch("mgo_max_energy",&tagger_info.mgo_max_energy,"mgo_max_energy/F");
    T_tagger->Branch("mgo_total_energy",&tagger_info.mgo_total_energy,"mgo_total_energy/F");
    T_tagger->Branch("mgo_n_showers",&tagger_info.mgo_n_showers,"mgo_n_showers/F");
    T_tagger->Branch("mgo_max_energy_1",&tagger_info.mgo_max_energy_1,"mgo_max_energy_1/F");
    T_tagger->Branch("mgo_max_energy_2",&tagger_info.mgo_max_energy_2,"mgo_max_energy_2/F");
    T_tagger->Branch("mgo_total_other_energy",&tagger_info.mgo_total_other_energy,"mgo_total_other_energy/F");
    T_tagger->Branch("mgo_n_total_showers",&tagger_info.mgo_n_total_showers,"mgo_n_total_showers/F");
    T_tagger->Branch("mgo_total_other_energy_1",&tagger_info.mgo_total_other_energy_1,"mgo_total_other_energy_1/F");
    T_tagger->Branch("mgo_flag",&tagger_info.mgo_flag,"mgo_flag/F");
    
    T_tagger->Branch("mgt_flag_single_shower",&tagger_info.mgt_flag_single_shower,"mgt_flag_single_shower/F");
    T_tagger->Branch("mgt_max_energy",&tagger_info.mgt_max_energy,"mgt_max_energy/F");
    T_tagger->Branch("mgt_energy",&tagger_info.mgt_energy,"mgt_energy/F");
    T_tagger->Branch("mgt_total_other_energy",&tagger_info.mgt_total_other_energy,"mgt_total_other_energy/F");
    T_tagger->Branch("mgt_max_energy_1",&tagger_info.mgt_max_energy_1,"mgt_max_energy_1/F");
    T_tagger->Branch("mgt_e_indirect_max_energy",&tagger_info.mgt_e_indirect_max_energy,"mgt_e_indirect_max_energy/F");
    T_tagger->Branch("mgt_e_direct_max_energy",&tagger_info.mgt_e_direct_max_energy,"mgt_e_direct_max_energy/F");
    T_tagger->Branch("mgt_n_direct_showers",&tagger_info.mgt_n_direct_showers,"mgt_n_direct_showers/F");
    T_tagger->Branch("mgt_e_direct_total_energy",&tagger_info.mgt_e_direct_total_energy,"mgt_e_direct_total_energy/F");
    T_tagger->Branch("mgt_flag_indirect_max_pio",&tagger_info.mgt_flag_indirect_max_pio,"mgt_flag_indirect_max_pio/F");
    T_tagger->Branch("mgt_e_indirect_total_energy",&tagger_info.mgt_e_indirect_total_energy,"mgt_e_indirect_total_energy/F");
    T_tagger->Branch("mgt_flag",&tagger_info.mgt_flag,"mgt_flag/F");
    
    T_tagger->Branch("sig_1_v_angle",&tagger_info.sig_1_v_angle);
    T_tagger->Branch("sig_1_v_flag_single_shower",&tagger_info.sig_1_v_flag_single_shower);
    T_tagger->Branch("sig_1_v_energy",&tagger_info.sig_1_v_energy);
    T_tagger->Branch("sig_1_v_energy_1",&tagger_info.sig_1_v_energy_1);
    T_tagger->Branch("sig_1_v_flag",&tagger_info.sig_1_v_flag);

    T_tagger->Branch("sig_2_v_energy",&tagger_info.sig_2_v_energy);
    T_tagger->Branch("sig_2_v_shower_angle",&tagger_info.sig_2_v_shower_angle);
    T_tagger->Branch("sig_2_v_flag_single_shower",&tagger_info.sig_2_v_flag_single_shower);
    T_tagger->Branch("sig_2_v_medium_dQ_dx",&tagger_info.sig_2_v_medium_dQ_dx);
    T_tagger->Branch("sig_2_v_start_dQ_dx",&tagger_info.sig_2_v_start_dQ_dx);
    T_tagger->Branch("sig_2_v_flag",&tagger_info.sig_2_v_flag);

    T_tagger->Branch("sig_flag",&tagger_info.sig_flag, "sig_flag/F");

    T_tagger->Branch("tro_1_v_particle_type",&tagger_info.tro_1_v_particle_type);
    T_tagger->Branch("tro_1_v_flag_dir_weak",&tagger_info.tro_1_v_flag_dir_weak);
    T_tagger->Branch("tro_1_v_min_dis",&tagger_info.tro_1_v_min_dis);
    T_tagger->Branch("tro_1_v_sg1_length",&tagger_info.tro_1_v_sg1_length);
    T_tagger->Branch("tro_1_v_shower_main_length",&tagger_info.tro_1_v_shower_main_length);
    T_tagger->Branch("tro_1_v_max_n_vtx_segs",&tagger_info.tro_1_v_max_n_vtx_segs);
    T_tagger->Branch("tro_1_v_tmp_length",&tagger_info.tro_1_v_tmp_length);
    T_tagger->Branch("tro_1_v_medium_dQ_dx",&tagger_info.tro_1_v_medium_dQ_dx);
    T_tagger->Branch("tro_1_v_dQ_dx_cut",&tagger_info.tro_1_v_dQ_dx_cut);
    T_tagger->Branch("tro_1_v_flag_shower_topology",&tagger_info.tro_1_v_flag_shower_topology);
    T_tagger->Branch("tro_1_v_flag",&tagger_info.tro_1_v_flag);

    T_tagger->Branch("tro_2_v_energy",&tagger_info.tro_2_v_energy);
    T_tagger->Branch("tro_2_v_stem_length",&tagger_info.tro_2_v_stem_length);
    T_tagger->Branch("tro_2_v_iso_angle",&tagger_info.tro_2_v_iso_angle);
    T_tagger->Branch("tro_2_v_max_length",&tagger_info.tro_2_v_max_length);
    T_tagger->Branch("tro_2_v_angle",&tagger_info.tro_2_v_angle);
    T_tagger->Branch("tro_2_v_flag",&tagger_info.tro_2_v_flag);

    T_tagger->Branch("tro_3_stem_length",&tagger_info.tro_3_stem_length,"tro_3_stem_length/F");
    T_tagger->Branch("tro_3_n_muon_segs",&tagger_info.tro_3_n_muon_segs,"tro_3_n_muon_segs/F");
    T_tagger->Branch("tro_3_energy",&tagger_info.tro_3_energy,"tro_3_energy/F");
    T_tagger->Branch("tro_3_flag",&tagger_info.tro_3_flag,"tro_3_flag/F");

    T_tagger->Branch("tro_4_v_dir2_mag",&tagger_info.tro_4_v_dir2_mag);
    T_tagger->Branch("tro_4_v_angle",&tagger_info.tro_4_v_angle);
    T_tagger->Branch("tro_4_v_angle1",&tagger_info.tro_4_v_angle1);
    T_tagger->Branch("tro_4_v_angle2",&tagger_info.tro_4_v_angle2);
    T_tagger->Branch("tro_4_v_length",&tagger_info.tro_4_v_length);
    T_tagger->Branch("tro_4_v_length1",&tagger_info.tro_4_v_length1);
    T_tagger->Branch("tro_4_v_medium_dQ_dx",&tagger_info.tro_4_v_medium_dQ_dx);
    T_tagger->Branch("tro_4_v_end_dQ_dx",&tagger_info.tro_4_v_end_dQ_dx);
    T_tagger->Branch("tro_4_v_energy",&tagger_info.tro_4_v_energy);
    T_tagger->Branch("tro_4_v_shower_main_length",&tagger_info.tro_4_v_shower_main_length);
    T_tagger->Branch("tro_4_v_flag_shower_trajectory",&tagger_info.tro_4_v_flag_shower_trajectory);
    T_tagger->Branch("tro_4_v_flag",&tagger_info.tro_4_v_flag);

    T_tagger->Branch("tro_5_v_max_angle",&tagger_info.tro_5_v_max_angle);
    T_tagger->Branch("tro_5_v_min_angle",&tagger_info.tro_5_v_min_angle);
    T_tagger->Branch("tro_5_v_max_length",&tagger_info.tro_5_v_max_length);
    T_tagger->Branch("tro_5_v_iso_angle",&tagger_info.tro_5_v_iso_angle);
    T_tagger->Branch("tro_5_v_n_vtx_segs",&tagger_info.tro_5_v_n_vtx_segs);
    T_tagger->Branch("tro_5_v_min_count",&tagger_info.tro_5_v_min_count);
    T_tagger->Branch("tro_5_v_max_count",&tagger_info.tro_5_v_max_count);
    T_tagger->Branch("tro_5_v_energy",&tagger_info.tro_5_v_energy);
    T_tagger->Branch("tro_5_v_flag",&tagger_info.tro_5_v_flag);

    T_tagger->Branch("tro_flag",&tagger_info.tro_flag,"tro_flag/F");


    // cosmic tagger ...
    T_tagger->Branch("cosmict_flag_1",&tagger_info.cosmict_flag_1,"cosmict_flag_1/F");
    T_tagger->Branch("cosmict_flag_2",&tagger_info.cosmict_flag_2,"cosmict_flag_2/F");
    T_tagger->Branch("cosmict_flag_3",&tagger_info.cosmict_flag_3,"cosmict_flag_3/F");
    T_tagger->Branch("cosmict_flag_4",&tagger_info.cosmict_flag_4,"cosmict_flag_4/F");
    T_tagger->Branch("cosmict_flag_5",&tagger_info.cosmict_flag_5,"cosmict_flag_5/F");
    T_tagger->Branch("cosmict_flag_6",&tagger_info.cosmict_flag_6,"cosmict_flag_6/F");
    T_tagger->Branch("cosmict_flag_7",&tagger_info.cosmict_flag_7,"cosmict_flag_7/F");
    T_tagger->Branch("cosmict_flag_8",&tagger_info.cosmict_flag_8,"cosmict_flag_8/F");
    T_tagger->Branch("cosmict_flag_9",&tagger_info.cosmict_flag_9,"cosmict_flag_9/F");
    T_tagger->Branch("cosmict_flag_10",&tagger_info.cosmict_flag_10);
    T_tagger->Branch("cosmict_flag",&tagger_info.cosmict_flag,"cosmict_flag/F");

    T_tagger->Branch("cosmict_2_filled",&tagger_info.cosmict_2_filled,"cosmict_2_filled/F");
    T_tagger->Branch("cosmict_2_particle_type",&tagger_info.cosmict_2_particle_type,"cosmict_2_particle_type/F");
    T_tagger->Branch("cosmict_2_n_muon_tracks",&tagger_info.cosmict_2_n_muon_tracks,"cosmict_2_n_muon_tracks/F");
    T_tagger->Branch("cosmict_2_total_shower_length",&tagger_info.cosmict_2_total_shower_length,"cosmict_2_total_shower_length/F");
    T_tagger->Branch("cosmict_2_flag_inside",&tagger_info.cosmict_2_flag_inside,"cosmict_2_flag_inside/F");
    T_tagger->Branch("cosmict_2_angle_beam",&tagger_info.cosmict_2_angle_beam,"cosmict_2_angle_beam/F");
    T_tagger->Branch("cosmict_2_flag_dir_weak",&tagger_info.cosmict_2_flag_dir_weak,"cosmict_2_flag_dir_weak/F");
    T_tagger->Branch("cosmict_2_dQ_dx_end",&tagger_info.cosmict_2_dQ_dx_end,"cosmict_2_dQ_dx_end/F");
    T_tagger->Branch("cosmict_2_dQ_dx_front",&tagger_info.cosmict_2_dQ_dx_front,"cosmict_2_dQ_dx_front/F");
    T_tagger->Branch("cosmict_2_theta",&tagger_info.cosmict_2_theta,"cosmict_2_theta/F");
    T_tagger->Branch("cosmict_2_phi",&tagger_info.cosmict_2_phi,"cosmict_2_phi/F");
    T_tagger->Branch("cosmict_2_valid_tracks",&tagger_info.cosmict_2_valid_tracks,"cosmict_2_valid_tracks/F");

    T_tagger->Branch("cosmict_3_filled",&tagger_info.cosmict_3_filled,"cosmict_3_filled/F");
    T_tagger->Branch("cosmict_3_flag_inside",&tagger_info.cosmict_3_flag_inside,"cosmict_3_flag_inside/F");
    T_tagger->Branch("cosmict_3_angle_beam",&tagger_info.cosmict_3_angle_beam,"cosmict_3_angle_beam/F");
    T_tagger->Branch("cosmict_3_flag_dir_weak",&tagger_info.cosmict_3_flag_dir_weak,"cosmict_3_flag_dir_weak/F");
    T_tagger->Branch("cosmict_3_dQ_dx_end",&tagger_info.cosmict_3_dQ_dx_end,"cosmict_3_dQ_dx_end/F");
    T_tagger->Branch("cosmict_3_dQ_dx_front",&tagger_info.cosmict_3_dQ_dx_front,"cosmict_3_dQ_dx_front/F");
    T_tagger->Branch("cosmict_3_theta",&tagger_info.cosmict_3_theta,"cosmict_3_theta/F");
    T_tagger->Branch("cosmict_3_phi",&tagger_info.cosmict_3_phi,"cosmict_3_phi/F");
    T_tagger->Branch("cosmict_3_valid_tracks",&tagger_info.cosmict_3_valid_tracks,"cosmict_3_valid_tracks/F");
    
    T_tagger->Branch("cosmict_4_filled",&tagger_info.cosmict_4_filled,"cosmict_4_filled/F");
    T_tagger->Branch("cosmict_4_flag_inside",&tagger_info.cosmict_4_flag_inside,"cosmict_4_flag_inside/F");
    T_tagger->Branch("cosmict_4_angle_beam",&tagger_info.cosmict_4_angle_beam,"cosmict_4_angle_beam/F");
    T_tagger->Branch("cosmict_4_connected_showers",&tagger_info.cosmict_4_connected_showers,"cosmict_4_connected_showers/F");

    T_tagger->Branch("cosmict_5_filled",&tagger_info.cosmict_5_filled,"cosmict_5_filled/F");
    T_tagger->Branch("cosmict_5_flag_inside",&tagger_info.cosmict_5_flag_inside,"cosmict_5_flag_inside/F");
    T_tagger->Branch("cosmict_5_angle_beam",&tagger_info.cosmict_5_angle_beam,"cosmict_5_angle_beam/F");
    T_tagger->Branch("cosmict_5_connected_showers",&tagger_info.cosmict_5_connected_showers,"cosmict_5_connected_showers/F");

    T_tagger->Branch("cosmict_6_filled",&tagger_info.cosmict_6_filled,"cosmict_6_filled/F");
    T_tagger->Branch("cosmict_6_flag_dir_weak",&tagger_info.cosmict_6_flag_dir_weak,"cosmict_6_flag_dir_weak/F");
    T_tagger->Branch("cosmict_6_flag_inside",&tagger_info.cosmict_6_flag_inside,"cosmict_6_flag_inside/F");
    T_tagger->Branch("cosmict_6_angle",&tagger_info.cosmict_6_angle,"cosmict_6_angle/F");
    

    T_tagger->Branch("cosmict_7_filled",&tagger_info.cosmict_7_filled,"cosmict_7_filled/F");
    T_tagger->Branch("cosmict_7_flag_sec",&tagger_info.cosmict_7_flag_sec,"cosmict_7_flag_sec/F");
    T_tagger->Branch("cosmict_7_n_muon_tracks",&tagger_info.cosmict_7_n_muon_tracks,"cosmict_7_n_muon_tracks/F");
    T_tagger->Branch("cosmict_7_total_shower_length",&tagger_info.cosmict_7_total_shower_length,"cosmict_7_total_shower_length/F");
    T_tagger->Branch("cosmict_7_flag_inside",&tagger_info.cosmict_7_flag_inside,"cosmict_7_flag_inside/F");
    T_tagger->Branch("cosmict_7_angle_beam",&tagger_info.cosmict_7_angle_beam,"cosmict_7_angle_beam/F");
    T_tagger->Branch("cosmict_7_flag_dir_weak",&tagger_info.cosmict_7_flag_dir_weak,"cosmict_7_flag_dir_weak/F");
    T_tagger->Branch("cosmict_7_dQ_dx_end",&tagger_info.cosmict_7_dQ_dx_end,"cosmict_7_dQ_dx_end/F");
    T_tagger->Branch("cosmict_7_dQ_dx_front",&tagger_info.cosmict_7_dQ_dx_front,"cosmict_7_dQ_dx_front/F");
    T_tagger->Branch("cosmict_7_theta",&tagger_info.cosmict_7_theta,"cosmict_7_theta/F");
    T_tagger->Branch("cosmict_7_phi",&tagger_info.cosmict_7_phi,"cosmict_7_phi/F");

    T_tagger->Branch("cosmict_8_filled",&tagger_info.cosmict_8_filled,"cosmict_8_filled/F");
    T_tagger->Branch("cosmict_8_flag_out",&tagger_info.cosmict_8_flag_out,"cosmict_8_flag_out/F");
    T_tagger->Branch("cosmict_8_muon_length",&tagger_info.cosmict_8_muon_length,"cosmict_8_muon_length/F");
    T_tagger->Branch("cosmict_8_acc_length",&tagger_info.cosmict_8_acc_length,"cosmict_8_acc_length/F");

    T_tagger->Branch("cosmict_10_flag_inside",&tagger_info.cosmict_10_flag_inside);
    T_tagger->Branch("cosmict_10_vtx_z",&tagger_info.cosmict_10_vtx_z);
    T_tagger->Branch("cosmict_10_flag_shower",&tagger_info.cosmict_10_flag_shower);
    T_tagger->Branch("cosmict_10_flag_dir_weak",&tagger_info.cosmict_10_flag_dir_weak);
    T_tagger->Branch("cosmict_10_angle_beam",&tagger_info.cosmict_10_angle_beam);
    T_tagger->Branch("cosmict_10_length",&tagger_info.cosmict_10_length);
    
    T_tagger->Branch("numu_cc_flag",&tagger_info.numu_cc_flag,"numu_cc_flag/F");

    T_tagger->Branch("numu_cc_flag_1",&tagger_info.numu_cc_flag_1);
    T_tagger->Branch("numu_cc_1_particle_type",&tagger_info.numu_cc_1_particle_type);
    T_tagger->Branch("numu_cc_1_length",&tagger_info.numu_cc_1_length);
    T_tagger->Branch("numu_cc_1_medium_dQ_dx",&tagger_info.numu_cc_1_medium_dQ_dx);
    T_tagger->Branch("numu_cc_1_dQ_dx_cut",&tagger_info.numu_cc_1_dQ_dx_cut);
    T_tagger->Branch("numu_cc_1_direct_length",&tagger_info.numu_cc_1_direct_length);
    T_tagger->Branch("numu_cc_1_n_daughter_tracks",&tagger_info.numu_cc_1_n_daughter_tracks);
    T_tagger->Branch("numu_cc_1_n_daughter_all",&tagger_info.numu_cc_1_n_daughter_all);

    T_tagger->Branch("numu_cc_flag_2",&tagger_info.numu_cc_flag_2);
    T_tagger->Branch("numu_cc_2_length",&tagger_info.numu_cc_2_length);
    T_tagger->Branch("numu_cc_2_total_length",&tagger_info.numu_cc_2_total_length);
    T_tagger->Branch("numu_cc_2_n_daughter_tracks",&tagger_info.numu_cc_2_n_daughter_tracks);
    T_tagger->Branch("numu_cc_2_n_daughter_all",&tagger_info.numu_cc_2_n_daughter_all);

    T_tagger->Branch("numu_cc_flag_3",&tagger_info.numu_cc_flag_3,"numu_cc_flag_3/F");
    T_tagger->Branch("numu_cc_3_particle_type",&tagger_info.numu_cc_3_particle_type,"numu_cc_3_particle_type/F");
    T_tagger->Branch("numu_cc_3_max_length",&tagger_info.numu_cc_3_max_length,"numu_cc_3_max_length/F");
    T_tagger->Branch("numu_cc_3_track_length",&tagger_info.numu_cc_3_acc_track_length,"numu_cc_3_acc_track_length/F");
    T_tagger->Branch("numu_cc_3_max_length_all",&tagger_info.numu_cc_3_max_length_all,"numu_cc_3_max_length_all/F");
    T_tagger->Branch("numu_cc_3_max_muon_length",&tagger_info.numu_cc_3_max_muon_length,"numu_cc_3_max_muon_length/F");
    T_tagger->Branch("numu_cc_3_n_daughter_tracks",&tagger_info.numu_cc_3_n_daughter_tracks,"numu_cc_3_n_daughter_tracks/F");
    T_tagger->Branch("numu_cc_3_n_daughter_all",&tagger_info.numu_cc_3_n_daughter_all,"numu_cc_3_n_daughter_all/F");

    // numu BDTs
    T_tagger->Branch("cosmict_2_4_score",&tagger_info.cosmict_2_4_score, "cosmict_2_4_score/F");
    T_tagger->Branch("cosmict_3_5_score",&tagger_info.cosmict_3_5_score, "cosmict_3_5_score/F");
    T_tagger->Branch("cosmict_6_score",&tagger_info.cosmict_6_score, "cosmict_6_score/F");
    T_tagger->Branch("cosmict_7_score",&tagger_info.cosmict_7_score, "cosmict_7_score/F");
    T_tagger->Branch("cosmict_8_score",&tagger_info.cosmict_8_score, "cosmict_8_score/F");
    T_tagger->Branch("cosmict_10_score",&tagger_info.cosmict_10_score, "cosmict_10_score/F");

    T_tagger->Branch("numu_1_score",&tagger_info.numu_1_score,"numu_1_score/F");
    T_tagger->Branch("numu_2_score",&tagger_info.numu_2_score,"numu_2_score/F");
    T_tagger->Branch("numu_3_score",&tagger_info.numu_3_score,"numu_3_score/F");

    T_tagger->Branch("cosmict_score",&tagger_info.cosmict_score,"cosmict_score/F");
    T_tagger->Branch("numu_score",&tagger_info.numu_score,"numu_score/F");
    
    
    // BDTs ...
    T_tagger->Branch("mipid_score",&tagger_info.mipid_score,"mipid_score/F");
    T_tagger->Branch("gap_score",&tagger_info.gap_score,"gap_score/F");
    T_tagger->Branch("hol_lol_score",&tagger_info.hol_lol_score,"hol_lol_score/F");
    T_tagger->Branch("cme_anc_score",&tagger_info.cme_anc_score,"cme_anc_score/F"); 
    T_tagger->Branch("mgo_mgt_score",&tagger_info.mgo_mgt_score,"mgo_mgt_score/F");
    T_tagger->Branch("br1_score",&tagger_info.br1_score,"br1_score/F");
    
    T_tagger->Branch("br3_score",&tagger_info.br3_score,"br3_score/F");
    T_tagger->Branch("br3_3_score",&tagger_info.br3_3_score,"br3_3_score/F");
    T_tagger->Branch("br3_5_score",&tagger_info.br3_5_score,"br3_5_score/F");
    T_tagger->Branch("br3_6_score",&tagger_info.br3_6_score,"br3_6_score/F");
    T_tagger->Branch("stemdir_br2_score",&tagger_info.stemdir_br2_score,"stemdir_br2_score/F");
    T_tagger->Branch("trimuon_score",&tagger_info.trimuon_score,"trimuon_score/F");
    
    T_tagger->Branch("br4_tro_score",&tagger_info.br4_tro_score,"br4_tro_score/F");
    T_tagger->Branch("mipquality_score",&tagger_info.mipquality_score,"mipquality_score/F");
    T_tagger->Branch("pio_1_score",&tagger_info.pio_1_score,"pio_1_score/F");
    T_tagger->Branch("pio_2_score",&tagger_info.pio_2_score,"pio_2_score/F");
    T_tagger->Branch("stw_spt_score",&tagger_info.stw_spt_score,"stw_spt_score/F");
    T_tagger->Branch("vis_1_score",&tagger_info.vis_1_score,"vis_1_score/F");
    
    T_tagger->Branch("vis_2_score",&tagger_info.vis_2_score,"vis_2_score/F");
    T_tagger->Branch("stw_2_score",&tagger_info.stw_2_score,"stw_2_score/F");
    T_tagger->Branch("stw_3_score",&tagger_info.stw_3_score,"stw_3_score/F");
    T_tagger->Branch("stw_4_score",&tagger_info.stw_4_score,"stw_4_score/F");
    T_tagger->Branch("sig_1_score",&tagger_info.sig_1_score,"sig_1_score/F");
    T_tagger->Branch("sig_2_score",&tagger_info.sig_2_score,"sig_2_score/F");

    T_tagger->Branch("lol_1_score",&tagger_info.lol_1_score,"lol_1_score/F");
    T_tagger->Branch("lol_2_score",&tagger_info.lol_2_score,"lol_2_score/F");
    T_tagger->Branch("tro_1_score",&tagger_info.tro_1_score,"tro_1_score/F");
    T_tagger->Branch("tro_2_score",&tagger_info.tro_2_score,"tro_2_score/F");
    T_tagger->Branch("tro_4_score",&tagger_info.tro_4_score,"tro_4_score/F");
    T_tagger->Branch("tro_5_score",&tagger_info.tro_5_score,"tro_5_score/F");
    T_tagger->Branch("nue_score",&tagger_info.nue_score,"nue_score/F");

    
    T_tagger->Branch("run",&tagger_info.run,"data/I");
    T_tagger->Branch("subrun",&tagger_info.subrun,"data/I");
    T_tagger->Branch("event",&tagger_info.event,"data/I");
    
    T_tagger->Branch("nuvtx_diff",&tagger_info.nuvtx_diff,"data/F");
    T_tagger->Branch("showervtx_diff",&tagger_info.showervtx_diff,"data/F");
    T_tagger->Branch("muonvtx_diff",&tagger_info.muonvtx_diff,"data/F");
    
    T_tagger->Branch("match_isFC",&tagger_info.match_isFC,"data/F");
    if (flag==1) T_tagger->Branch("truth_isCC",&tagger_info.truth_isCC,"data/F");
    if (flag==1) T_tagger->Branch("truth_vtxInside",&tagger_info.truth_vtxInside,"data/F");
    if (flag==1) T_tagger->Branch("truth_nuPdg",&tagger_info.truth_nuPdg,"data/F");
    if (flag==1) T_tagger->Branch("truth_nuEnergy",&tagger_info.truth_nuEnergy,"data/F");
    if (flag==1) T_tagger->Branch("truth_nuIntType",&tagger_info.truth_nuIntType,"data/F");
    if (flag==1) T_tagger->Branch("truth_energyInside",&tagger_info.truth_energyInside,"data/F");
    if (flag==1) T_tagger->Branch("weight_spline",&tagger_info.weight_spline,"data/F");
    if (flag==1) T_tagger->Branch("weight_cv",&tagger_info.weight_cv,"data/F");
    if (flag==1) T_tagger->Branch("weight_lee",&tagger_info.weight_lee,"data/F");
    T_tagger->Branch("kine_reco_Enu",&tagger_info.kine_reco_Enu,"data/F");
    
    if (flag==1) T_tagger->Branch("kine_reco_add_energy",&tagger_info.kine_reco_add_energy,"data/F");
    if (flag==1) T_tagger->Branch("kine_pio_mass",&tagger_info.kine_pio_mass,"data/F");
    if (flag==1) T_tagger->Branch("kine_pio_vtx_dis",&tagger_info.kine_pio_vtx_dis,"data/F");
    if (flag==1) T_tagger->Branch("kine_pio_flag",&tagger_info.kine_pio_flag,"data/F");
    if (flag==1) T_tagger->Branch("kine_pio_energy_1",&tagger_info.kine_pio_energy_1,"data/F");
    if (flag==1) T_tagger->Branch("kine_pio_theta_1",&tagger_info.kine_pio_theta_1,"data/F");
    if (flag==1) T_tagger->Branch("kine_pio_phi_1",&tagger_info.kine_pio_phi_1,"data/F");
    if (flag==1) T_tagger->Branch("kine_pio_dis_1",&tagger_info.kine_pio_dis_1,"data/F");
    if (flag==1) T_tagger->Branch("kine_pio_energy_2",&tagger_info.kine_pio_energy_2,"data/F");
    if (flag==1) T_tagger->Branch("kine_pio_theta_2",&tagger_info.kine_pio_theta_2,"data/F");
    if (flag==1) T_tagger->Branch("kine_pio_phi_2",&tagger_info.kine_pio_phi_2,"data/F");
    if (flag==1) T_tagger->Branch("kine_pio_dis_2",&tagger_info.kine_pio_dis_2,"data/F");
    if (flag==1) T_tagger->Branch("kine_pio_angle",&tagger_info.kine_pio_angle,"data/F");

    
    if (flag==1) T_tagger->Branch("event_type",&tagger_info.event_type,"data/F");

    if (flag==1) T_tagger->Branch("weight",&tagger_info.weight,"data/F");
    if (flag==1) T_tagger->Branch("lowEweight",&tagger_info.lowEweight,"data/F");
}
