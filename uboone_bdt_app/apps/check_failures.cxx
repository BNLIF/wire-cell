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
  double total_pot_even = 0;
  double total_pot_odd = 0;
  
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
    if (subRunNo %2 == 0) total_pot_even += pot_tor875;
    if (subRunNo %2 == 1) total_pot_odd += pot_tor875;
    
  }

  std::cout << filename << std::endl;
  std::cout << "Total POT: " << total_pot << "   " << total_pot_even << " (even)   " << total_pot_odd << " (odd)" << std::endl;

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
  T_PFeval->SetBranchStatus("reco_nuvtxX",1); T_PFeval->SetBranchAddress("reco_nuvtxX", &reco_nuvtxX);

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

  Float_t n_gen_nue[2]={0,0};
  Float_t n_gen_nue_err2[2]={0,0};
  Float_t n_all_nue[2]={0,0};
  Float_t n_all_nue_err2[2]={0,0};

  Float_t n_gen_nuelee[2]={0,0};
  Float_t n_gen_nuelee_err2[2]={0,0};
  Float_t n_all_nuelee[2]={0,0};
  Float_t n_all_nuelee_err2[2]={0,0};
  
  // numuCC passing rate after generic neutrino, also need to calculate efficiency  (even vs. odd)

  Float_t n_gen_pr_numu[4]={0,0,0,0};
  Float_t n_gen_pr_numu_err2[4]={0,0,0,0};
  Float_t n_gen_all[4]={0,0,0,0};
  Float_t n_gen_all_err2[4]={0,0,0,0};

  Float_t n_gen_prt_numu[4]={0,0,0,0};
  Float_t n_gen_prt_numu_err2[4]={0,0,0,0};
  Float_t n_gen_all_numu[4]={0,0,0,0};
  Float_t n_gen_all_numu_err2[4]={0,0,0,0};

  // pi0 passing rate after generic neutrino
  Float_t n_gen_pr_pio[2]={0,0};
  Float_t n_gen_pr_pio_err2[2]={0,0};
  
  // neCC passing rate after generic neutrino, also need to calculate efficiency (event vs. odd, intrinsic vs. LEE)
  Float_t n_gen_pr_nue[4]={0,0,0,0};
  Float_t n_gen_pr_nue_err2[4]={0,0,0,0};
  
  Float_t n_gen_prt_nue[4]={0,0,0,0};
  Float_t n_gen_prt_nue_err2[4]={0,0,0,0};
  Float_t n_gen_all_nue[4]={0,0,0,0};
  Float_t n_gen_all_nue_err2[4]={0,0,0,0};

  Float_t n_gen_pr_nuelee[4]={0,0,0,0};
  Float_t n_gen_pr_nuelee_err2[4]={0,0,0,0};
  
  Float_t n_gen_prt_nuelee[4]={0,0,0,0};
  Float_t n_gen_prt_nuelee_err2[4]={0,0,0,0};
  Float_t n_gen_all_nuelee[4]={0,0,0,0};
  Float_t n_gen_all_nuelee_err2[4]={0,0,0,0};
  
  
  bool flag_presel = false;
  bool flag_generic = false;
  bool flag_numu = false;
  bool flag_gen_numu = false;

  bool flag_nue = false;
  bool flag_gen_nue = false;

  bool flag_gen_pr_numu = false;
  bool flag_gen_pr_nue = false;
  bool flag_gen_pr_pio = false;
  
  
  double weight = 1;
  for (Int_t i=0;i!=T_eval->GetEntries();i++){
    T_eval->GetEntry(i);
    T_PFeval->GetEntry(i);
    T_BDTvars->GetEntry(i);
    T_kine->GetEntry(i);
    
    flag_presel = false;
    flag_generic = false;
    flag_gen_pr_numu = false;
    flag_gen_pr_nue = false;
    flag_gen_pr_pio = false;
    
    if (match_found != 0 &&
	stm_eventtype != 0 && stm_lowenergy ==0 && stm_LM ==0 && stm_TGM ==0 && stm_STM==0 && stm_FullDead == 0 && stm_cluster_length >0) {
      flag_presel = true; // preselection ...
      if (stm_cluster_length > 15) flag_generic = true; // generic
    }

    if (numu_score > 0.9 && flag_generic) flag_gen_pr_numu = true;
    if (nue_score > 7.0 && flag_generic) flag_gen_pr_nue = true;
    
    if (flag_generic && kine_pio_flag==1 && kine_pio_energy_1 > 15 && kine_pio_energy_2 > 15 && kine_pio_dis_1 < 80 && kine_pio_dis_2 < 80 && kine_pio_angle > 20 && kine_pio_vtx_dis < 1) flag_gen_pr_pio = true;
    
    if (flag_data ==0){
      flag_numu = false;
      flag_gen_numu = false;
      
      flag_nue = false;
      flag_gen_nue = false;
      
      // deal with weight
      if (std::isnan(weight_cv) || std::isinf(weight_cv) || weight_cv <= 0 || weight_cv > 1000) weight_cv = 1;
      if (std::isnan(weight_spline) || std::isinf(weight_spline) || weight_spline <=0 || weight_spline > 1000) weight_spline = 1;
      weight = weight_cv * weight_spline;
      
      if (truth_vtxInside==1 && truth_nuPdg == 14 && truth_isCC == 1) flag_numu = true;
      if (flag_numu && flag_generic) flag_gen_numu = true;

      if (truth_vtxInside==1 && abs(truth_nuPdg) == 12 && truth_isCC==1) flag_nue = true;
      if (flag_nue && flag_generic) flag_gen_nue = true;
    }
        
    
    // generic neutrino selections
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

    
    if (flag_nue){
      n_all_nue[0] ++; n_all_nue[1] += weight;
      n_all_nue_err2[0]++; n_all_nue_err2[1] += pow(weight,2);
    }
    if (flag_gen_nue){
      n_gen_nue[0] ++; n_gen_nue[1] += weight;
      n_gen_nue_err2[0]++; n_gen_nue_err2[1] += pow(weight,2);
    }

    if (weight_lee <0) weight_lee = 0;
    if (flag_nue){
      n_all_nuelee[0] ++; n_all_nuelee[1] += weight*weight_lee;
      n_all_nuelee_err2[0]++; n_all_nuelee_err2[1] += pow(weight*weight_lee,2);
    }
    if (flag_gen_nue){
      n_gen_nuelee[0] ++; n_gen_nuelee[1] += weight*weight_lee;
      n_gen_nuelee_err2[0]++; n_gen_nuelee_err2[1] += pow(weight*weight_lee,2);
    }

   
    // PR ...
    
  
    if (flag_gen_pr_numu){
      if (subrun%2==0){
	n_gen_pr_numu[0] ++; n_gen_pr_numu[1] += weight;
	n_gen_pr_numu_err2[0] ++; n_gen_pr_numu_err2[1] += pow(weight,2);
      }else{
	n_gen_pr_numu[2] ++; n_gen_pr_numu[3] += weight;
	n_gen_pr_numu_err2[2] ++; n_gen_pr_numu_err2[3] += pow(weight,2);
      }
      
      if (flag_numu){
	if (subrun%2==0){
	  n_gen_prt_numu[0] ++; n_gen_prt_numu[1] += weight;
	  n_gen_prt_numu_err2[0] ++; n_gen_prt_numu_err2[1] += pow(weight,2);
	}else{
	  n_gen_prt_numu[2] ++; n_gen_prt_numu[3] += weight;
	  n_gen_prt_numu_err2[2] ++; n_gen_prt_numu_err2[3] += pow(weight,2);
	}
      }      
    }

    //    if (flag_gen_pr_numu &&  (!flag_generic)) std::cout << "Something wrong! " << std::endl;
    if (flag_gen_pr_nue){
      if (subrun%2==0){
	n_gen_pr_nue[0] ++; n_gen_pr_nue[1] += weight;
	n_gen_pr_nue_err2[0] ++; n_gen_pr_nue_err2[1] += pow(weight,2);
      }else{
	n_gen_pr_nue[2] ++; n_gen_pr_nue[3] += weight;
	n_gen_pr_nue_err2[2] ++; n_gen_pr_nue_err2[3] += pow(weight,2);
      }
      
      if (flag_nue){
	if (subrun%2==0){
	  n_gen_prt_nue[0] ++; n_gen_prt_nue[1] += weight;
	  n_gen_prt_nue_err2[0] ++; n_gen_prt_nue_err2[1] += pow(weight,2);

	  n_gen_prt_nuelee[0] ++; n_gen_prt_nuelee[1] += weight*weight_lee;
	  n_gen_prt_nuelee_err2[0] ++; n_gen_prt_nuelee_err2[1] += pow(weight*weight_lee,2);
	}else{
	  n_gen_prt_nue[2] ++; n_gen_prt_nue[3] += weight;
	  n_gen_prt_nue_err2[2] ++; n_gen_prt_nue_err2[3] += pow(weight,2);

	  n_gen_prt_nuelee[2] ++; n_gen_prt_nuelee[3] += weight*weight_lee;
	  n_gen_prt_nuelee_err2[2] ++; n_gen_prt_nuelee_err2[3] += pow(weight*weight_lee,2);
	}
      }      
    }

    
    if (flag_generic){
      if (subrun%2==0){
	n_gen_all[0] ++; n_gen_all[1] += weight;
	n_gen_all_err2[0] ++; n_gen_all_err2[1] += pow(weight,2);
      }else{
	n_gen_all[2] ++; n_gen_all[3] += weight;
	n_gen_all_err2[2] ++; n_gen_all_err2[3] += pow(weight,2);
      }
      
      if (flag_numu){
	if (subrun%2==0){
	  n_gen_all_numu[0] ++; n_gen_all_numu[1] += weight;
	  n_gen_all_numu_err2[0] ++; n_gen_all_numu_err2[1] += pow(weight,2);
	}else{
	  n_gen_all_numu[2] ++; n_gen_all_numu[3] += weight;
	  n_gen_all_numu_err2[2] ++; n_gen_all_numu_err2[3] += pow(weight,2);
	}
      }

      if (flag_nue){
	if (subrun%2==0){
	  n_gen_all_nue[0] ++; n_gen_all_nue[1] += weight;
	  n_gen_all_nue_err2[0] ++; n_gen_all_nue_err2[1] += pow(weight,2);
	  
	  n_gen_all_nuelee[0] ++; n_gen_all_nuelee[1] += weight*weight_lee;
	  n_gen_all_nuelee_err2[0] ++; n_gen_all_nuelee_err2[1] += pow(weight*weight_lee,2);
	}else{
	  n_gen_all_nue[2] ++; n_gen_all_nue[3] += weight;
	  n_gen_all_nue_err2[2] ++; n_gen_all_nue_err2[3] += pow(weight,2);

	  n_gen_all_nuelee[2] ++; n_gen_all_nuelee[3] += weight*weight_lee;
	  n_gen_all_nuelee_err2[2] ++; n_gen_all_nuelee_err2[3] += pow(weight*weight_lee,2);
	}
      }
    }

    // pio
    if (flag_gen_pr_pio){
      n_gen_pr_pio[0]++; n_gen_pr_pio[1] += weight;
      n_gen_pr_pio_err2[0]++; n_gen_pr_pio_err2[1] += pow(weight,2);
    }
    
    
    // add events ...
    map_rs_n[std::make_pair(run, subrun)] ++;
    // Reco 1.5 failure ...
    if (match_found == -1) map_rs_f1p5[std::make_pair(run, subrun)].insert(event);
    if (match_found == 1 && stm_lowenergy == -1) map_rs_f2stm[std::make_pair(run, subrun)].insert(event);
    if (flag_presel && numu_cc_flag == -1) map_rs_f2pr[std::make_pair(run, subrun)].insert(event);
  }

  std::cout << endl;
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
  std::cout << endl;
  // generic neutrino

  {
    double gen_ratio;
    double gen_ratio_err;
    double gen_ratio_weight;
    double gen_ratio_weight_err;
    double gen_ratio_weight1;
    double gen_ratio_weight1_err;

    gen_ratio = n_gen[0]/n_all[0];
    gen_ratio_err = sqrt(n_gen_err2[0] * pow(n_all[0]-n_gen[0],2) + pow(n_gen[0],2) * (n_all_err2[0] - n_gen_err2[0]))/pow(n_all[0],2);
    gen_ratio_weight = n_gen[1]/n_all[1];
    gen_ratio_weight_err = sqrt(n_gen_err2[1] * pow(n_all[1] - n_gen[1],2) + pow(n_gen[1],2) * (n_all_err2[1] - n_gen_err2[1]))/pow(n_all[1],2);

    std::cout << "Generic passing rate: " << gen_ratio << "+-" << gen_ratio_err << "  w. weight: " << gen_ratio_weight << "+-" << gen_ratio_weight_err << std::endl;

    gen_ratio = n_gen[0]/total_pot * 5e19;
    gen_ratio_err = sqrt(n_gen_err2[0])/total_pot * 5e19;

    gen_ratio_weight = n_gen[1]/total_pot * 5e19;
    gen_ratio_weight_err = sqrt(n_gen_err2[1])/total_pot * 5e19;
    
    std::cout << "Generic rate @ 5e19POT: " << gen_ratio << "+-" << gen_ratio_err << "  w. weight: " << gen_ratio_weight << "+-" << gen_ratio_weight_err << std::endl;
    
    if (flag_data == 0){
      gen_ratio = n_gen_numu[0]/n_all_numu[0];
      gen_ratio_err = sqrt(n_gen_numu_err2[0] * pow(n_all_numu[0]-n_gen_numu[0],2) + pow(n_gen_numu[0],2) * (n_all_numu_err2[0] - n_gen_numu_err2[0]) )/pow(n_all_numu[0],2);
      gen_ratio_weight = n_gen_numu[1]/n_all_numu[1];
      gen_ratio_weight_err = sqrt(n_gen_numu_err2[1] * pow(n_all_numu[1]-n_gen_numu[1],2) + pow(n_gen_numu[1],2) * (n_all_numu_err2[1] - n_gen_numu_err2[1]) )/pow(n_all_numu[1],2);
      
      std::cout << "Generic numuCC eff: " << gen_ratio << "+-" << gen_ratio_err << "  w. weight: " << gen_ratio_weight << "+-" << gen_ratio_weight_err << std::endl;

      gen_ratio = n_gen_nue[0]/n_all_nue[0];
      gen_ratio_err = sqrt(n_gen_nue_err2[0] * pow(n_all_nue[0]-n_gen_nue[0],2) + pow(n_gen_nue[0],2) * (n_all_nue_err2[0] - n_gen_nue_err2[0]) )/pow(n_all_nue[0],2);
      gen_ratio_weight = n_gen_nue[1]/n_all_nue[1];
      gen_ratio_weight_err = sqrt(n_gen_nue_err2[1] * pow(n_all_nue[1]-n_gen_nue[1],2) + pow(n_gen_nue[1],2) * (n_all_nue_err2[1] - n_gen_nue_err2[1]) )/pow(n_all_nue[1],2);

      gen_ratio_weight1 = n_gen_nuelee[1]/n_all_nuelee[1];
      gen_ratio_weight1_err = sqrt(n_gen_nuelee_err2[1] * pow(n_all_nuelee[1]-n_gen_nuelee[1],2) + pow(n_gen_nuelee[1],2) * (n_all_nuelee_err2[1] - n_gen_nuelee_err2[1]) )/pow(n_all_nuelee[1],2);
      
      std::cout << "Generic nueCC eff: " << gen_ratio << "+-" << gen_ratio_err << "  w. weight: " << gen_ratio_weight << "+-" << gen_ratio_weight_err << " LEE:" << gen_ratio_weight1 << "+-" << gen_ratio_weight1_err << std::endl;

    }
  }
  std::cout << endl;

  // numu passing rate
   // Float_t n_gen_pr_numu[4]={0,0,0,0};
  // Float_t n_gen_pr_numu_err2[4]={0,0,0,0};
  // Float_t n_gen_all[4]={0,0,0,0};
  // Float_t n_gen_all_err2[4]={0,0,0,0};
  {
    double gen_ratio;
    double gen_ratio_err;
    double gen_ratio_weight;
    double gen_ratio_weight_err;
    double gen_ratio_weight1;
    double gen_ratio_weight1_err;
    
    double gen_ratio1;
    double gen_ratio1_err;
    double gen_ratio1_weight;
    double gen_ratio1_weight_err;
    double gen_ratio1_weight1;
    double gen_ratio1_weight1_err;

    gen_ratio = n_gen_pr_numu[0]/n_gen_all[0];
    gen_ratio_err = sqrt(n_gen_pr_numu_err2[0] * pow(n_gen_all[0]-n_gen_pr_numu[0],2) + pow(n_gen_pr_numu[0],2) * (n_gen_all_err2[0] - n_gen_pr_numu_err2[0]))/pow(n_gen_all[0],2);
    
    gen_ratio_weight = n_gen_pr_numu[1]/n_gen_all[1];
    gen_ratio_weight_err = sqrt(n_gen_pr_numu_err2[1] * pow(n_gen_all[1]-n_gen_pr_numu[1],2) + pow(n_gen_pr_numu[1],2) * (n_gen_all_err2[1] - n_gen_pr_numu_err2[1]))/pow(n_gen_all[1],2);

    gen_ratio1 = n_gen_pr_numu[2]/n_gen_all[2];
    gen_ratio1_err = sqrt(n_gen_pr_numu_err2[2] * pow(n_gen_all[2]-n_gen_pr_numu[2],2) + pow(n_gen_pr_numu[2],2) * (n_gen_all_err2[2] - n_gen_pr_numu_err2[2]))/pow(n_gen_all[2],2);
    
    gen_ratio1_weight = n_gen_pr_numu[3]/n_gen_all[3];
    gen_ratio1_weight_err = sqrt(n_gen_pr_numu_err2[3] * pow(n_gen_all[3]-n_gen_pr_numu[3],2) + pow(n_gen_pr_numu[3],2) * (n_gen_all_err2[3] - n_gen_pr_numu_err2[3]))/pow(n_gen_all[3],2);

    std::cout << "numuCC passing rate w.r.t. generic (even): " << gen_ratio << "+-" << gen_ratio_err << "  w. weight: " << gen_ratio_weight << "+-" << gen_ratio_weight_err << std::endl;
    std::cout << "numuCC passing rate w.r.t. generic (odd): "  << gen_ratio1 << "+-" << gen_ratio1_err << "  w. weight: " << gen_ratio1_weight << "+-" << gen_ratio1_weight_err << std::endl;


    gen_ratio = n_gen_pr_numu[0]/total_pot_even * 5e19;
    gen_ratio_err = sqrt(n_gen_pr_numu_err2[0])/total_pot_even * 5e19;

    gen_ratio_weight = n_gen_pr_numu[1]/total_pot_even * 5e19;
    gen_ratio_weight_err = sqrt(n_gen_pr_numu_err2[1])/total_pot_even * 5e19;

    gen_ratio1 = n_gen_pr_numu[2]/total_pot_odd * 5e19;
    gen_ratio1_err = sqrt(n_gen_pr_numu_err2[2])/total_pot_odd * 5e19;

    gen_ratio1_weight = n_gen_pr_numu[3]/total_pot_odd * 5e19;
    gen_ratio1_weight_err = sqrt(n_gen_pr_numu_err2[3])/total_pot_odd * 5e19;
    
    std::cout << "numuCC rate (even) @ 5e19POT: " << gen_ratio << "+-" << gen_ratio_err << "  w. weight: " << gen_ratio_weight << "+-" << gen_ratio_weight_err << std::endl;
    std::cout << "numuCC rate (odd) @ 5e19POT: " << gen_ratio1 << "+-" << gen_ratio1_err << "  w. weight: " << gen_ratio1_weight << "+-" << gen_ratio1_weight_err << std::endl;

    // Float_t n_gen_prt_numu[4]={0,0,0,0};
    // Float_t n_gen_prt_numu_err2[4]={0,0,0,0};
    // Float_t n_gen_all_numu[4]={0,0,0,0};
    // Float_t n_gen_all_numu_err2[4]={0,0,0,0};
    if (flag_data == 0){
      gen_ratio = n_gen_prt_numu[0]/n_gen_all_numu[0];
      gen_ratio_err = sqrt(n_gen_prt_numu_err2[0] * pow(n_gen_all_numu[0]-n_gen_prt_numu[0],2) + pow(n_gen_prt_numu[0],2) * (n_gen_all_numu_err2[0] - n_gen_prt_numu_err2[0]))/pow(n_gen_all_numu[0],2);
      
      gen_ratio_weight = n_gen_prt_numu[1]/n_gen_all_numu[1];
      //    std::cout << n_gen_all_numu[1] << " " << n_gen_prt_numu[1] << " " << n_gen_all_numu_err2[1] << " " << n_gen_prt_numu_err2[1] << " " <<   pow(n_gen_all_numu[1]-n_gen_prt_numu[1],2) << " " <<  pow(n_gen_prt_numu[1],2) << " " << (n_gen_all_numu_err2[1] - n_gen_prt_numu_err2[1]) << " " << std::endl;
      
      gen_ratio_weight_err = sqrt(n_gen_prt_numu_err2[1] * pow(n_gen_all_numu[1]-n_gen_prt_numu[1],2) + pow(n_gen_prt_numu[1],2) * (n_gen_all_numu_err2[1] - n_gen_prt_numu_err2[1]))/pow(n_gen_all_numu[1],2);
      
      gen_ratio1 = n_gen_prt_numu[2]/n_gen_all_numu[2];
      gen_ratio1_err = sqrt(n_gen_prt_numu_err2[2] * pow(n_gen_all_numu[2]-n_gen_prt_numu[2],2) + pow(n_gen_prt_numu[2],2) * (n_gen_all_numu_err2[2] - n_gen_prt_numu_err2[2]))/pow(n_gen_all_numu[2],2);
      
      gen_ratio1_weight = n_gen_prt_numu[3]/n_gen_all_numu[3];
      gen_ratio1_weight_err = sqrt(n_gen_prt_numu_err2[3] * pow(n_gen_all_numu[3]-n_gen_prt_numu[3],2) + pow(n_gen_prt_numu[3],2) * (n_gen_all_numu_err2[3] - n_gen_prt_numu_err2[3]))/pow(n_gen_all_numu[3],2);
      
      std::cout << "numuCC eff. w.r.t. generic (even): " << gen_ratio << "+-" << gen_ratio_err << "  w. weight: " << gen_ratio_weight << "+-" << gen_ratio_weight_err << std::endl;
      std::cout << "numuCC eff. w.r.t. generic (odd): "  << gen_ratio1 << "+-" << gen_ratio1_err << "  w. weight: " << gen_ratio1_weight << "+-" << gen_ratio1_weight_err << std::endl;
    }
    std::cout << std::endl;
    
    // pio ...
    gen_ratio = n_gen_pr_pio[0]/n_gen[0];
    gen_ratio_err = sqrt(n_gen_pr_pio_err2[0] * pow(n_gen[0]-n_gen_pr_pio[0],2) + pow(n_gen_pr_pio[0],2) * (n_gen_err2[0] - n_gen_pr_pio_err2[0]))/pow(n_gen[0],2);
    gen_ratio_weight = n_gen_pr_pio[1]/n_gen[1];
    gen_ratio_weight_err = sqrt(n_gen_pr_pio_err2[1] * pow(n_gen[1] - n_gen_pr_pio[1],2) + pow(n_gen_pr_pio[1],2) * (n_gen_err2[1] - n_gen_pr_pio_err2[1]))/pow(n_gen[1],2);

    std::cout << "Pi0 passing rate w.r.t. generic: " << gen_ratio << "+-" << gen_ratio_err << "  w. weight: " << gen_ratio_weight << "+-" << gen_ratio_weight_err << std::endl;

    gen_ratio = n_gen_pr_pio[0]/total_pot * 5e19;
    gen_ratio_err = sqrt(n_gen_pr_pio_err2[0])/total_pot * 5e19;

    gen_ratio_weight = n_gen_pr_pio[1]/total_pot * 5e19;
    gen_ratio_weight_err = sqrt(n_gen_pr_pio_err2[1])/total_pot * 5e19;
    
    std::cout << "Pi0 rate @ 5e19POT: " << gen_ratio << "+-" << gen_ratio_err << "  w. weight: " << gen_ratio_weight << "+-" << gen_ratio_weight_err << std::endl;

    std::cout << std::endl;

    // nueCC
    gen_ratio = n_gen_pr_nue[0]/n_gen_all[0];
    gen_ratio_err = sqrt(n_gen_pr_nue_err2[0] * pow(n_gen_all[0]-n_gen_pr_nue[0],2) + pow(n_gen_pr_nue[0],2) * (n_gen_all_err2[0] - n_gen_pr_nue_err2[0]))/pow(n_gen_all[0],2);
    
    gen_ratio_weight = n_gen_pr_nue[1]/n_gen_all[1];
    gen_ratio_weight_err = sqrt(n_gen_pr_nue_err2[1] * pow(n_gen_all[1]-n_gen_pr_nue[1],2) + pow(n_gen_pr_nue[1],2) * (n_gen_all_err2[1] - n_gen_pr_nue_err2[1]))/pow(n_gen_all[1],2);

    gen_ratio1 = n_gen_pr_nue[2]/n_gen_all[2];
    gen_ratio1_err = sqrt(n_gen_pr_nue_err2[2] * pow(n_gen_all[2]-n_gen_pr_nue[2],2) + pow(n_gen_pr_nue[2],2) * (n_gen_all_err2[2] - n_gen_pr_nue_err2[2]))/pow(n_gen_all[2],2);
    
    gen_ratio1_weight = n_gen_pr_nue[3]/n_gen_all[3];
    gen_ratio1_weight_err = sqrt(n_gen_pr_nue_err2[3] * pow(n_gen_all[3]-n_gen_pr_nue[3],2) + pow(n_gen_pr_nue[3],2) * (n_gen_all_err2[3] - n_gen_pr_nue_err2[3]))/pow(n_gen_all[3],2);

    std::cout << "nueCC passing rate w.r.t. generic (even): " << gen_ratio << "+-" << gen_ratio_err << "  w. weight: " << gen_ratio_weight << "+-" << gen_ratio_weight_err << std::endl;
    std::cout << "nueCC passing rate w.r.t. generic (odd): "  << gen_ratio1 << "+-" << gen_ratio1_err << "  w. weight: " << gen_ratio1_weight << "+-" << gen_ratio1_weight_err << std::endl;


    gen_ratio = n_gen_pr_nue[0]/total_pot_even * 5e19;
    gen_ratio_err = sqrt(n_gen_pr_nue_err2[0])/total_pot_even * 5e19;

    gen_ratio_weight = n_gen_pr_nue[1]/total_pot_even * 5e19;
    gen_ratio_weight_err = sqrt(n_gen_pr_nue_err2[1])/total_pot_even * 5e19;

    gen_ratio1 = n_gen_pr_nue[2]/total_pot_odd * 5e19;
    gen_ratio1_err = sqrt(n_gen_pr_nue_err2[2])/total_pot_odd * 5e19;

    gen_ratio1_weight = n_gen_pr_nue[3]/total_pot_odd * 5e19;
    gen_ratio1_weight_err = sqrt(n_gen_pr_nue_err2[3])/total_pot_odd * 5e19;
    
    std::cout << "nueCC rate (even) @ 5e19POT: " << gen_ratio << "+-" << gen_ratio_err << "  w. weight: " << gen_ratio_weight << "+-" << gen_ratio_weight_err << std::endl;
    std::cout << "nueCC rate (odd) @ 5e19POT: " << gen_ratio1 << "+-" << gen_ratio1_err << "  w. weight: " << gen_ratio1_weight << "+-" << gen_ratio1_weight_err << std::endl;

    // Float_t n_gen_prt_nue[4]={0,0,0,0};
    // Float_t n_gen_prt_nue_err2[4]={0,0,0,0};
    // Float_t n_gen_all_nue[4]={0,0,0,0};
    // Float_t n_gen_all_nue_err2[4]={0,0,0,0};
    if (flag_data == 0){
      gen_ratio = n_gen_prt_nue[0]/n_gen_all_nue[0];
      gen_ratio_err = sqrt(n_gen_prt_nue_err2[0] * pow(n_gen_all_nue[0]-n_gen_prt_nue[0],2) + pow(n_gen_prt_nue[0],2) * (n_gen_all_nue_err2[0] - n_gen_prt_nue_err2[0]))/pow(n_gen_all_nue[0],2);
      
      gen_ratio_weight = n_gen_prt_nue[1]/n_gen_all_nue[1];
      //    std::cout << n_gen_all_nue[1] << " " << n_gen_prt_nue[1] << " " << n_gen_all_nue_err2[1] << " " << n_gen_prt_nue_err2[1] << " " <<   pow(n_gen_all_nue[1]-n_gen_prt_nue[1],2) << " " <<  pow(n_gen_prt_nue[1],2) << " " << (n_gen_all_nue_err2[1] - n_gen_prt_nue_err2[1]) << " " << std::endl;
      
      gen_ratio_weight_err = sqrt(n_gen_prt_nue_err2[1] * pow(n_gen_all_nue[1]-n_gen_prt_nue[1],2) + pow(n_gen_prt_nue[1],2) * (n_gen_all_nue_err2[1] - n_gen_prt_nue_err2[1]))/pow(n_gen_all_nue[1],2);
      
      gen_ratio1 = n_gen_prt_nue[2]/n_gen_all_nue[2];
      gen_ratio1_err = sqrt(n_gen_prt_nue_err2[2] * pow(n_gen_all_nue[2]-n_gen_prt_nue[2],2) + pow(n_gen_prt_nue[2],2) * (n_gen_all_nue_err2[2] - n_gen_prt_nue_err2[2]))/pow(n_gen_all_nue[2],2);
      
      gen_ratio1_weight = n_gen_prt_nue[3]/n_gen_all_nue[3];
      gen_ratio1_weight_err = sqrt(n_gen_prt_nue_err2[3] * pow(n_gen_all_nue[3]-n_gen_prt_nue[3],2) + pow(n_gen_prt_nue[3],2) * (n_gen_all_nue_err2[3] - n_gen_prt_nue_err2[3]))/pow(n_gen_all_nue[3],2);
      
      std::cout << "nueCC eff. w.r.t. generic (even): " << gen_ratio << "+-" << gen_ratio_err << "  w. weight: " << gen_ratio_weight << "+-" << gen_ratio_weight_err << std::endl;
      std::cout << "nueCC eff. w.r.t. generic (odd): "  << gen_ratio1 << "+-" << gen_ratio1_err << "  w. weight: " << gen_ratio1_weight << "+-" << gen_ratio1_weight_err << std::endl;

      gen_ratio = n_gen_prt_nuelee[0]/n_gen_all_nuelee[0];
      gen_ratio_err = sqrt(n_gen_prt_nuelee_err2[0] * pow(n_gen_all_nuelee[0]-n_gen_prt_nuelee[0],2) + pow(n_gen_prt_nuelee[0],2) * (n_gen_all_nuelee_err2[0] - n_gen_prt_nuelee_err2[0]))/pow(n_gen_all_nuelee[0],2);
      
      gen_ratio_weight = n_gen_prt_nuelee[1]/n_gen_all_nuelee[1];
      //    std::cout << n_gen_all_nuelee[1] << " " << n_gen_prt_nuelee[1] << " " << n_gen_all_nuelee_err2[1] << " " << n_gen_prt_nuelee_err2[1] << " " <<   pow(n_gen_all_nuelee[1]-n_gen_prt_nuelee[1],2) << " " <<  pow(n_gen_prt_nuelee[1],2) << " " << (n_gen_all_nuelee_err2[1] - n_gen_prt_nuelee_err2[1]) << " " << std::endl;
      
      gen_ratio_weight_err = sqrt(n_gen_prt_nuelee_err2[1] * pow(n_gen_all_nuelee[1]-n_gen_prt_nuelee[1],2) + pow(n_gen_prt_nuelee[1],2) * (n_gen_all_nuelee_err2[1] - n_gen_prt_nuelee_err2[1]))/pow(n_gen_all_nuelee[1],2);
      
      gen_ratio1 = n_gen_prt_nuelee[2]/n_gen_all_nuelee[2];
      gen_ratio1_err = sqrt(n_gen_prt_nuelee_err2[2] * pow(n_gen_all_nuelee[2]-n_gen_prt_nuelee[2],2) + pow(n_gen_prt_nuelee[2],2) * (n_gen_all_nuelee_err2[2] - n_gen_prt_nuelee_err2[2]))/pow(n_gen_all_nuelee[2],2);
      
      gen_ratio1_weight = n_gen_prt_nuelee[3]/n_gen_all_nuelee[3];
      gen_ratio1_weight_err = sqrt(n_gen_prt_nuelee_err2[3] * pow(n_gen_all_nuelee[3]-n_gen_prt_nuelee[3],2) + pow(n_gen_prt_nuelee[3],2) * (n_gen_all_nuelee_err2[3] - n_gen_prt_nuelee_err2[3]))/pow(n_gen_all_nuelee[3],2);
      
      std::cout << "LEE eff. w.r.t. generic (even): " << gen_ratio << "+-" << gen_ratio_err << "  w. weight: " << gen_ratio_weight << "+-" << gen_ratio_weight_err << std::endl;
      std::cout << "LEE eff. w.r.t. generic (odd): "  << gen_ratio1 << "+-" << gen_ratio1_err << "  w. weight: " << gen_ratio1_weight << "+-" << gen_ratio1_weight_err << std::endl;

      
    }

    
  }

  

  
  
  
  std::cout << std::endl << std::endl << std::endl;
  
  //std::cout << filename << " " << POT << std::endl;

  return 0;
}
