#ifndef UBOONE_LEE_CUTS
#define UBOONE_LEE_CUTS

// define cuts here ...

#include "tagger.h"
#include "kine.h"
#include "eval.h"

// generic neutrino cuts
TCut generic_cut = "match_found == 1 && stm_eventtype != 0 &&stm_lowenergy ==0 && stm_LM ==0 && stm_TGM ==0 && stm_STM==0 && stm_FullDead == 0 && stm_cluster_length >15";
bool is_generic(EvalInfo& info);

// preselection cuts
TCut preselect_cut = "match_found == 1 && stm_eventtype != 0 &&stm_lowenergy ==0 && stm_LM ==0 && stm_TGM ==0 && stm_STM==0 && stm_FullDead == 0 && stm_cluster_length > 0";
bool is_preselection(EvalInfo& info); 

// nueCC cuts
TCut nueCC_cut = "numu_cc_flag >=0 && nue_score > 7.0";
bool is_nueCC(TaggerInfo& tagger_info);

// numuCC cuts
TCut numuCC_cut = "numu_cc_flag >=0 && numu_score > 0.9";
bool is_numuCC(TaggerInfo& tagger_info);

// pio cuts
TCut pi0_cut = "kine_pio_flag==1 && kine_pio_energy_1 > 15 && kine_pio_energy_2 > 15 && kine_pio_dis_1 < 80 && kine_pio_dis_2 < 80 && kine_pio_angle > 20 && kine_pio_vtx_dis < 1";
bool is_pi0();

// NC cuts
TCut NC_cut = "(!cosmict_flag) && numu_score < 0.0";
bool is_NC(TaggerInfo& tagger_info);



bool is_pi0(KineInfo& kine){
  bool flag = false;

  if (kine.kine_pio_flag==1 && kine.kine_pio_energy_1 > 15 && kine.kine_pio_energy_2 > 15 && kine.kine_pio_dis_1 < 80 && kine.kine_pio_dis_2 < 80 && kine.kine_pio_angle > 20 && kine.kine_pio_vtx_dis < 1)
    flag = true;

  
  return flag;
}



bool is_NC(TaggerInfo& tagger_info){
  bool flag = false;
  if ((!tagger_info.cosmict_flag) && tagger_info.numu_score < 0)
    flag = true;
  
  return flag;
}


bool is_numuCC{TaggerInfo& tagger_info){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9)
    flag = true;
  
  return flag;
}


bool is_nueCC(TaggerInfo& tagger_info){
  bool flag = false;

  if (tagger_info.numu_cc_flag >=0 && tagger_info.nue_score > 7.0)
    flag = true;
  
  return flag;
}


bool is_generic(EvalInfo& eval){
  // not very useful for the main analysis
  bool flag = is_preselection(eval);

  flag = flag && (eval.stm_cluster_length > 15);
  return flag;
}

bool is_preselection(EvalInfo& eval){
  bool flag = false;

  // match code ...
  int tmp_match_found = eval.match_found;
  if (eval.is_match_found_int){
    tmp_match_found = eval.match_found_asInt;
  }

  if (tmp_match_found == 1 && eval.stm_eventtype != 0 && eval.stm_lowenergy ==0 && eval.stm_LM ==0 && eval.stm_TGM ==0 && eval.stm_STM==0 && eval.stm_FullDead == 0 && eval.stm_cluster_length >0) flag = true;
  
  
  return flag;
}


#endif
