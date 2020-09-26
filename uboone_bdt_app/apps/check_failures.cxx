#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"

using namespace std;

int main( int argc, char** argv )
{
  if (argc < 2) {
    std::cout << "./check_passing_rate #filename -d[data 1/MC 0]" << std::endl;
    return -1;
  }

  int flag_data = 0; // MC (with truth information)
  int flag_verbose = 0;
  for (Int_t i=1;i!=argc;i++){
    switch(argv[i][1]){
    case 'd':
      flag_data = atoi(&argv[i][2]); // whether this is data or not ...
      break;
    case 'v':
      flag_verbose = atoi(&argv[i][2]); // whether this is data or not ...
      break;
    }
  }
  
  TString filename = argv[1];

  
  TFile *file = new TFile(filename);
  TTree *T_pot = (TTree*)file->Get("wcpselection/T_pot"); // POT tree ...
  bool flag_pot = 0;
  if (T_pot->GetEntries() != 0 ) flag_pot = 1;

  Int_t runNo, subRunNo;
  Double_t pot_tor875;
  T_pot->SetBranchAddress("runNo",&runNo);
  T_pot->SetBranchAddress("subRunNo",&subRunNo);
  T_pot->SetBranchAddress("pot_tor875",&pot_tor875);
  
  double total_pot = 0;
  std::map<std::pair<int,int>, double> map_rs_pot;
  std::map<std::pair<int, int>, int> map_rs_n;
  std::map<std::pair<int, int>, std::set<int> > map_rs_f1p5; // Reco 1.5
  std::map<std::pair<int, int>, std::set<int> > map_rs_f2stm; // Reco2 stm
  std::map<std::pair<int, int>, std::set<int> > map_rs_f2pr; // Reco2 Pattern recognition
  
  
  for (Int_t i=0;i!=T_pot->GetEntries();i++){
    T_pot->GetEntry(i);
    total_pot += pot_tor875;
    map_rs_pot[std::make_pair(runNo, subRunNo)] = pot_tor875;
    map_rs_n[std::make_pair(runNo, subRunNo)] = 0; // set 0
  }

  std::cout << filename << std::endl;
  std::cout << "Total POT: " << total_pot << std::endl;

  // check failure mode
  TTree *T_eval = (TTree*)file->Get("wcpselection/T_eval");
  T_eval->SetBranchStatus("*",0);
  Int_t stm_eventtype, stm_lowenergy, stm_LM, stm_TGM, stm_STM, stm_FullDead;
  Float_t stm_cluster_length;
  
  T_eval->SetBranchStatus("stm_eventtype",1); T_eval->SetBranchAddress("stm_eventtype",&stm_eventtype);
  T_eval->SetBranchStatus("stm_lowenergy",1); T_eval->SetBranchAddress("stm_lowenergy",&stm_lowenergy);
  T_eval->SetBranchStatus("stm_LM",1); T_eval->SetBranchAddress("stm_LM",&stm_LM);
  T_eval->SetBranchStatus("stm_TGM",1); T_eval->SetBranchAddress("stm_TGM",&stm_TGM);
  T_eval->SetBranchStatus("stm_STM",1); T_eval->SetBranchAddress("stm_STM",&stm_STM);
  T_eval->SetBranchStatus("stm_FullDead",1); T_eval->SetBranchAddress("stm_FullDead",&stm_FullDead);
  T_eval->SetBranchStatus("stm_clusterlength",1); T_eval->SetBranchAddress("stm_clusterlength",&stm_cluster_length);

  bool match_found;
  T_eval->SetBranchStatus("match_found",1); T_eval->SetBranchAddress("match_found",&match_found);
  Int_t run, subrun, event;
  T_eval->SetBranchStatus("run",1); T_eval->SetBranchAddress("run",&run);
  T_eval->SetBranchStatus("subrun",1); T_eval->SetBranchAddress("subrun",&subrun);
  T_eval->SetBranchStatus("event",1); T_eval->SetBranchAddress("event",&event);

  if (!T_eval->GetBranch("truth_isCC")) flag_data = 1;
  
  Int_t truth_nuPdg;
  Bool_t truth_isCC;
  Bool_t truth_isFC;
  Bool_t truth_vtxInside;
  Float_t weight_spline;
  Float_t weight_cv;
  Float_t weight_lee;
  if (flag_data==0){
    T_eval->SetBranchStatus("truth_nuPdg",1); T_eval->SetBranchAddress("truth_nuPdg",&truth_nuPdg);
    T_eval->SetBranchStatus("truth_isCC",1); T_eval->SetBranchAddress("truth_isCC",&truth_isCC);
    T_eval->SetBranchStatus("truth_isFC",1); T_eval->SetBranchAddress("truth_isFC",&truth_isFC);
    T_eval->SetBranchStatus("truth_vtxInside",1); T_eval->SetBranchAddress("truth_vtxInside",&truth_vtxInside);
    T_eval->SetBranchStatus("weight_spline",1); T_eval->SetBranchAddress("weight_spline",&weight_spline);
    T_eval->SetBranchStatus("weight_cv",1); T_eval->SetBranchAddress("weight_cv",&weight_cv);
    T_eval->SetBranchStatus("weight_lee",1); T_eval->SetBranchAddress("weight_lee",&weight_lee);
  }
  
  
  TTree *T_PFeval = (TTree*)file->Get("wcpselection/T_PFeval");
  T_PFeval->SetBranchStatus("*",0);
  Float_t reco_nuvtxX;
  T_PFeval->SetBranchAddress("reco_nuvtxX", &reco_nuvtxX);

  TTree *T_BDTvars = (TTree*)file->Get("wcpselection/T_BDTvars");
  T_BDTvars->SetBranchStatus("*",0);
  Float_t numu_cc_flag;
  T_BDTvars->SetBranchStatus("numu_cc_flag",1); T_BDTvars->SetBranchAddress("numu_cc_flag",&numu_cc_flag);

  Float_t numu_score;
  Float_t nue_score;
  Float_t cosmict_flag;
  T_BDTvars->SetBranchStatus("numu_score",1);  T_BDTvars->SetBranchAddress("numu_score",&numu_score);
  T_BDTvars->SetBranchStatus("nue_score",1);  T_BDTvars->SetBranchAddress("nue_score",&nue_score);
  T_BDTvars->SetBranchStatus("cosmict_flag",1);  T_BDTvars->SetBranchAddress("cosmict_flag",&cosmict_flag);

  
  TTree *T_kine = (TTree*)file->Get("wcpselection/T_KINEvars"); // Kine tree for pi0 ...
  T_kine->SetBranchStatus("*",0);
  Int_t kine_pio_flag;
  Float_t kine_pio_vtx_dis, kine_pio_angle;
  Float_t kine_pio_energy_1, kine_pio_theta_1, kine_pio_phi_1, kine_pio_dis_1;
  Float_t kine_pio_energy_2, kine_pio_theta_2, kine_pio_phi_2, kine_pio_dis_2;
  
  T_kine->SetBranchStatus("kine_pio_flag",1); T_kine->SetBranchAddress("kine_pio_flag",&kine_pio_flag);
  T_kine->SetBranchStatus("kine_pio_vtx_dis",1); T_kine->SetBranchAddress("kine_pio_vtx_dis",&kine_pio_vtx_dis);
  T_kine->SetBranchStatus("kine_pio_angle",1); T_kine->SetBranchAddress("kine_pio_angle",&kine_pio_angle);
  T_kine->SetBranchStatus("kine_pio_energy_1",1); T_kine->SetBranchAddress("kine_pio_energy_1",&kine_pio_energy_1);
  T_kine->SetBranchStatus("kine_pio_theta_1",1); T_kine->SetBranchAddress("kine_pio_theta_1",&kine_pio_theta_1);
  T_kine->SetBranchStatus("kine_pio_phi_1",1); T_kine->SetBranchAddress("kine_pio_phi_1",&kine_pio_phi_1);
  T_kine->SetBranchStatus("kine_pio_dis_1",1); T_kine->SetBranchAddress("kine_pio_dis_1",&kine_pio_dis_1);
  T_kine->SetBranchStatus("kine_pio_energy_2",1); T_kine->SetBranchAddress("kine_pio_energy_2",&kine_pio_energy_2);
  T_kine->SetBranchStatus("kine_pio_theta_2",1); T_kine->SetBranchAddress("kine_pio_theta_2",&kine_pio_theta_2);
  T_kine->SetBranchStatus("kine_pio_phi_2",1); T_kine->SetBranchAddress("kine_pio_phi_2",&kine_pio_phi_2);
  T_kine->SetBranchStatus("kine_pio_dis_2",1); T_kine->SetBranchAddress("kine_pio_dis_2",&kine_pio_dis_2);
  
  // generic neutrino's rate (without weight, and with weight for MC only ...)
  Float_t n_gen[2]={0,0}; // passing generic selection, the second one is with weight ...
  Float_t n_gen_err2[2]={0,0};
  Float_t n_all[2]={0,0}; // passing all, the second one is with weight 
  Float_t n_all_err2[2]={0,0};

  Float_t n_gen_numu[2]={0,0};
  Float_t n_gen_numu_err2[2]={0,0};
  Float_t n_all_numu[2]={0,0};
  Float_t n_all_numu_err2[2]={0,0}; 
  
  // numuCC passing rate after generic neutrino, also need to calculate efficiency

  // neCC passing rate after generic neutrino, also need to calculate efficiency

  // pi0 passing rate after generic neutrino
  
  
  bool flag_presel = false;
  bool flag_generic = false;
  bool flag_numu = false;
  bool flag_gen_numu = false;
  
  double weight = 1;
  for (Int_t i=0;i!=T_eval->GetEntries();i++){
    T_eval->GetEntry(i);
    T_PFeval->GetEntry(i);
    
    flag_presel = false;
    flag_generic = false;

    if (match_found != 0 && stm_eventtype != 0 && stm_lowenergy ==0 && stm_LM ==0 && stm_TGM ==0 && stm_STM==0 && stm_FullDead == 0 && stm_cluster_length >0) {
      flag_presel = true; // preselection ...
      if (stm_cluster_length > 15) flag_generic = true; // generic
    }
    
    if (flag_data ==0){
      flag_numu = false;
      flag_gen_numu = false;
      
      // deal with weight
      if (std::isnan(weight_cv) || std::isinf(weight_cv) || weight_cv <= 0 || weight_cv > 1000) weight_cv = 1;
      if (std::isnan(weight_spline) || std::isinf(weight_spline) || weight_spline <=0 || weight_spline > 1000) weight_spline = 1;
      weight = weight_cv * weight_spline;
      
      if (truth_vtxInside==1 && truth_nuPdg == 14 && truth_isCC == 1) flag_numu = true;
      if (flag_numu && flag_generic) flag_gen_numu = true;
    }
        
    
    
    n_all[0] ++; n_all[1] += weight;
    n_all_err2[0] ++; n_all_err2[1] += pow(weight,2);
    
    if (flag_generic){
      n_gen[0] ++; n_gen[1] += weight;
      n_gen_err2[0] ++; n_gen_err2[1] += pow(weight,2);
    }
    if (flag_numu){
      n_all_numu[0] ++; n_all_numu[1] += weight;
      n_all_numu_err2[0]++; n_all_numu_err2[1] += pow(weight,2);
    }
    if (flag_gen_numu){
      n_gen_numu[0] ++; n_gen_numu[1] += weight;
      n_gen_numu_err2[0]++; n_gen_numu_err2[1] += pow(weight,2);
    }
   

    
    // add events ...
    map_rs_n[std::make_pair(run, subrun)] ++;
    // Reco 1.5 failure ...
    if (match_found == -1) map_rs_f1p5[std::make_pair(run, subrun)].insert(event);
    if (match_found == 1 && stm_lowenergy == -1) map_rs_f2stm[std::make_pair(run, subrun)].insert(event);
    if (flag_presel && numu_cc_flag == -1) map_rs_f2pr[std::make_pair(run, subrun)].insert(event);
  }

  // Reco 1.5
  {
    double tmp_pot = 0;
    for (auto it = map_rs_f1p5.begin(); it != map_rs_f1p5.end(); it++){
      tmp_pot += map_rs_pot[it->first];
    }
    std::cout << "Reco 1.5 failure in subrun: " << map_rs_f1p5.size() << "/" << map_rs_n.size() << " ratio: " << map_rs_f1p5.size()*1.0 /map_rs_n.size() << "   POT: " << tmp_pot << " ratio: " << tmp_pot/total_pot << std::endl;

    tmp_pot = 0;
    int tmp_event = 0;
    for (auto it = map_rs_f1p5.begin(); it != map_rs_f1p5.end(); it++){
      if (flag_verbose) std::cout << "1p5: " << it->first.first << " " << it->first.second << " " << map_rs_f1p5[it->first].size() << "/" << map_rs_n[it->first] << " " << std::endl; 

      tmp_event += map_rs_f1p5[it->first].size();
      tmp_pot += map_rs_f1p5[it->first].size() * 1.0 / map_rs_n[it->first] * map_rs_pot[it->first];
    }
    std::cout << "Reco 1.5 failure in event: " << tmp_event << "/" << T_eval->GetEntries() << "  POT: " << tmp_pot << " ratio " << tmp_pot/total_pot << std::endl;
    
  }

  // Reco 2 stm
  {
    double tmp_pot = 0;
    for (auto it = map_rs_f2stm.begin(); it != map_rs_f2stm.end(); it++){
      //  bool a = map_rs_pot.find(it->first) == map_rs_pot.end();
      //std::cout << it->first.first << " " << it->first.second << " " << a << std::endl;
      tmp_pot += map_rs_pot[it->first];
    }
    std::cout << "Reco 2.0 stm failure in subrun: " << map_rs_f2stm.size() << "/" << map_rs_n.size() << " POT: " << tmp_pot << " ratio " << tmp_pot/total_pot << std::endl;

    tmp_pot = 0;
    int tmp_event = 0;
    for (auto it = map_rs_f2stm.begin(); it != map_rs_f2stm.end(); it++){
      if (flag_verbose) std::cout << "stm: " << it->first.first << " " << it->first.second << " " << map_rs_f2stm[it->first].size() << "/" << map_rs_n[it->first] << " " << std::endl; 

      tmp_event += map_rs_f2stm[it->first].size();
      tmp_pot += map_rs_f2stm[it->first].size() * 1.0 / map_rs_n[it->first] * map_rs_pot[it->first];
    }
    std::cout << "Reco 2.0 stm failure in event: " << tmp_event << "/" << T_eval->GetEntries() << "  POT: " << tmp_pot << " ratio " << tmp_pot/total_pot << std::endl;   
  }

  // Reco 2 pr
  {
    double tmp_pot = 0;
    for (auto it = map_rs_f2pr.begin(); it != map_rs_f2pr.end(); it++){
      tmp_pot += map_rs_pot[it->first];
    }
    std::cout << "Reco 2.0 PR failure in subrun: " << map_rs_f2pr.size() << "/" << map_rs_n.size() << " POT: " << tmp_pot << " ratio " << tmp_pot/total_pot << std::endl;

    tmp_pot = 0;
    int tmp_event = 0;
    for (auto it = map_rs_f2pr.begin(); it != map_rs_f2pr.end(); it++){
      if (flag_verbose) std::cout << "PR: " << it->first.first << " " << it->first.second << " " << map_rs_f2pr[it->first].size() << "/" << map_rs_n[it->first] << " " << std::endl; 

      tmp_event += map_rs_f2pr[it->first].size();
      tmp_pot += map_rs_f2pr[it->first].size() * 1.0 / map_rs_n[it->first] * map_rs_pot[it->first];
    }
    std::cout << "Reco 2.0 PR failure in event: " << tmp_event << "/" << T_eval->GetEntries() << "  POT: " << tmp_pot << " ratio " << tmp_pot/total_pot << std::endl;
  }
  
  // generic neutrino

  // Float_t n_gen_numu[2]={0,0};
  // Float_t n_gen_numu_err2[2]={0,0};
  // Float_t n_all_numu[2]={0,0};
  // Float_t n_all_numu_err2[2]={0,0}; 
  {
    double gen_ratio;
    double gen_ratio_err;
    double gen_ratio_weight;
    double gen_ratio_weight_err;

    gen_ratio = n_gen[0]/n_all[0];
    gen_ratio_err = sqrt(n_gen_err2[0] * pow(n_all[0]-n_gen[0],2) + pow(n_gen[0],2) * (n_all_err2[0] - n_gen_err2[0]))/pow(n_all[0],2);
    gen_ratio_weight = n_gen[1]/n_all[1];
    gen_ratio_weight_err = sqrt(n_gen_err2[1] * pow(n_all[1] - n_gen[1],2) + pow(n_gen[1],2) * (n_all_err2[1] - n_gen_err2[1]))/pow(n_all[1],2);

    std::cout << "Generic passing rate: " << gen_ratio << "+-" << gen_ratio_err << "  w. weight: " << gen_ratio_weight << "+-" << gen_ratio_weight_err << std::endl;

    if (flag_data == 0){
      gen_ratio = n_gen_numu[0]/n_all_numu[0];
      gen_ratio_err = sqrt(n_gen_numu_err2[0] * pow(n_all_numu[0]-n_gen_numu[0],2) + pow(n_gen_numu[0],2) * (n_all_numu_err2[0] - n_gen_numu_err2[0]) )/pow(n_all_numu[0],2);
      gen_ratio_weight = n_gen_numu[1]/n_all_numu[1];
      gen_ratio_weight_err = sqrt(n_gen_numu_err2[1] * pow(n_all_numu[1]-n_gen_numu[1],2) + pow(n_gen_numu[1],2) * (n_all_numu_err2[1] - n_gen_numu_err2[1]) )/pow(n_all_numu[1],2);
      
      std::cout << "Generic numuCC eff: " << gen_ratio << "+-" << gen_ratio_err << "  w. weight: " << gen_ratio_weight << "+-" << gen_ratio_weight_err << std::endl;
    }
  }

  

  
  
  
  
  
  //std::cout << filename << " " << POT << std::endl;

  return 0;
}
