#ifndef UBOONE_LEE_CUTS
#define UBOONE_LEE_CUTS

// define cuts here ...
#include "TCut.h"

#include "tagger.h"
#include "kine.h"
#include "eval.h"

namespace LEEana{
  
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
TCut pi0_cut = "kine_pio_flag==1 && kine_pio_energy_1 > 40 && kine_pio_energy_2 > 25 && kine_pio_dis_1 < 110 && kine_pio_dis_2 < 120 && kine_pio_angle > 0 && kine_pio_vtx_dis < 9 && kine_pio_angle < 174 && kine_pio_mass > 22 && kine_pio_mass < 300";
bool is_pi0(KineInfo& kine);

// NC cuts
TCut NC_cut = "(!cosmict_flag) && numu_score < 0.0";
bool is_NC(TaggerInfo& tagger_info);


TCut FC_cut = "match_isFC==1";
TCut PC_cut = "match_isFC==0";

bool is_FC(EvalInfo& eval);


TCut truth_nueCC_inside = "abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1";
TCut truth_numuCC_inside = "abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1";
bool is_truth_nueCC_inside(EvalInfo& eval);
bool is_truth_numuCC_inside(EvalInfo& eval);
}


bool LEEana::is_truth_nueCC_inside(EvalInfo& eval){
  bool flag = false;

  if (fabs(eval.truth_nuPdg)==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1)
    flag = true;
  
  return flag;
}

bool LEEana::is_truth_numuCC_inside(EvalInfo& eval){
   bool flag = false;

  if (fabs(eval.truth_nuPdg)==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1)
    flag = true;
  
  return flag;
}



bool LEEana::is_FC(EvalInfo& eval){
  if (eval.match_isFC){
    return true;
  }else{
    return false;
  }
}


bool LEEana::is_pi0(KineInfo& kine){
  bool flag = false;

  if (kine.kine_pio_flag==1 && kine.kine_pio_energy_1 > 40 && kine.kine_pio_energy_2 > 25 && kine.kine_pio_dis_1 < 110 && kine.kine_pio_dis_2 < 120 && kine.kine_pio_angle > 0 && kine.kine_pio_angle < 174 && kine.kine_pio_vtx_dis < 9 && kine.kine_pio_mass > 22 && kine.kine_pio_mass < 300)
    flag = true;

  
  return flag;
}



bool LEEana::is_NC(TaggerInfo& tagger_info){
  bool flag = false;
  if ((!tagger_info.cosmict_flag) && tagger_info.numu_score < 0)
    flag = true;
  
  return flag;
}


bool LEEana::is_numuCC(TaggerInfo& tagger_info){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9)
    flag = true;
  
  return flag;
}



bool LEEana::is_nueCC(TaggerInfo& tagger_info){
  bool flag = false;

  if (tagger_info.numu_cc_flag >=0 && tagger_info.nue_score > 7.0)
    flag = true;
  
  return flag;
}


bool LEEana::is_generic(EvalInfo& eval){
  // not very useful for the main analysis
  bool flag = is_preselection(eval);

  flag = flag && (eval.stm_clusterlength > 15);
  return flag;
}

bool LEEana::is_preselection(EvalInfo& eval){
  bool flag = false;

  // match code ...
  int tmp_match_found = eval.match_found;
  if (eval.is_match_found_int){
    tmp_match_found = eval.match_found_asInt;
  }

  if (tmp_match_found == 1 && eval.stm_eventtype != 0 && eval.stm_lowenergy ==0 && eval.stm_LM ==0 && eval.stm_TGM ==0 && eval.stm_STM==0 && eval.stm_FullDead == 0 && eval.stm_clusterlength >0) flag = true;
  
  
  return flag;
}


#endif
