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

#include "WCPLEEANA/tagger.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"


#include "WCPLEEANA/eval.h"

using namespace std;
using namespace LEEana;

#include "WCPLEEANA/bdt.h"
#include "WCPLEEANA/pot.h"
#include "WCPLEEANA/pfeval.h"
#include "WCPLEEANA/kine.h"

int main( int argc, char** argv )
{
  if (argc < 3) {
    std::cout << "bdt_convert #input_file #output_file -c[weight_cut_val] -l[traing_list] -g[global_file_type]" << std::endl;
    return -1;
  }

  TString input_file = argv[1];
  TString out_file = argv[2];

  float weight_cut_val = 1000;
  float fail_percentage = 0.15;

  TString training_list = "";
  string global_file_type = "";
  for (Int_t i=1;i!=argc;i++){
    switch(argv[i][1]){
    case 'c':
      weight_cut_val = atof(&argv[i][2]);
      break;
    case 't':
      fail_percentage = atof(&argv[i][2]);
      break;
    case 'l':
      training_list = &argv[i][2];
      break;
    case 'g':
      global_file_type = &argv[i][2];
      break;
    }
  }
  
  bool flag_check_run_subrun = false;
  bool flag_use_global_file_type = false;
  if (global_file_type != "") flag_use_global_file_type = true;
  
  
  std::map<string, std::set<std::pair<int, int> > > map_type_run_subrun;
  if (training_list != ""){
    flag_check_run_subrun = true;
    
    ifstream infile(training_list);
    string tmp_type;
    int run, subrun;
    while(!infile.eof()){
      infile >> tmp_type >> run >> subrun;
      map_type_run_subrun[tmp_type].insert(std::make_pair(run, subrun));
    }
    // std::cout << map_type_run_subrun.size() << std::endl;
    // return 0;
  }
  
  
  
  bool flag_data = true;
  //std::cout << input_file << " " << out_file << std::endl;


  TFile *file1 = new TFile(input_file);
  TTree *T_BDTvars = (TTree*)file1->Get("wcpselection/T_BDTvars");
  
  TTree *T_eval = (TTree*)file1->Get("wcpselection/T_eval");
  TTree *T_pot = (TTree*)file1->Get("wcpselection/T_pot");
  TTree *T_PFeval = (TTree*)file1->Get("wcpselection/T_PFeval");
  TTree *T_KINEvars = (TTree*)file1->Get("wcpselection/T_KINEvars");

  

  
  if (T_eval->GetBranch("weight_cv")) flag_data =false;
  if (T_eval->GetBranch("file_type")) flag_use_global_file_type = false;

  // std::cout << flag_use_global_file_type << " " << flag_check_run_subrun << std::endl;
  // return 0;

  
  //  std::cout << T_eval->GetEntries() << std::endl;
  TFile *file2 = new TFile(out_file,"RECREATE");
  file2->mkdir("wcpselection");
  file2->cd("wcpselection");
  TTree *t4 = new TTree("T_BDTvars","T_BDTvars");
  TTree *t1 = new TTree("T_eval","T_eval");
  TTree *t2 = new TTree("T_pot","T_pot");
  TTree *t3 = new TTree("T_PFeval", "T_PFeval");
  TTree *t5 = new TTree("T_KINEvars", "T_KINEvars");
  //TTree *t1 = T_eval->CloneTree(-1,"");
  //TTree *t2 = T_pot->CloneTree(-1,"");
  //TTree *t3 = T_PFeval->CloneTree(-1,"");
  //  TTree *t5 = T_KINEvars->CloneTree(-1,"");
  
  
  EvalInfo eval;
  eval.file_type = new std::string();
  
  POTInfo pot;
  TaggerInfo tagger;
  PFevalInfo pfeval;
  KineInfo kine;

  kine.kine_energy_particle = new std::vector<float>;
  kine.kine_energy_info = new std::vector<int>;
  kine.kine_particle_type = new std::vector<int>;
  kine.kine_energy_included = new std::vector<int>;
  
  
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

  set_tree_address(T_BDTvars, tagger,2 );
  put_tree_address(t4, tagger,2);

  if (flag_data){
    set_tree_address(T_eval, eval,2);
    put_tree_address(t1, eval,2);

    set_tree_address(T_PFeval, pfeval,2);
    put_tree_address(t3, pfeval,2);
  }else{
    set_tree_address(T_eval, eval);
    put_tree_address(t1, eval);

    set_tree_address(T_PFeval, pfeval);
    put_tree_address(t3, pfeval);
  }

  set_tree_address(T_pot, pot);
  put_tree_address(t2, pot);



  set_tree_address(T_KINEvars, kine);
  put_tree_address(t5, kine);
  
  
  //  bool match_isFC;
  //T_eval->SetBranchAddress("match_isFC",&match_isFC);
  //  T_KINEvars->SetBranchAddress("kine_reco_Enu",&tagger.kine_reco_Enu);


  // BDTs ...
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

  
  // total BDTs ... to be added ...
  TMVA::Reader reader;
  reader.AddVariable("match_isFC",&tagger.match_isFC);
  reader.AddVariable("kine_reco_Enu",&tagger.kine_reco_Enu);
  
  reader.AddVariable("cme_mu_energy",&tagger.cme_mu_energy);
  reader.AddVariable("cme_energy",&tagger.cme_energy);
  reader.AddVariable("cme_mu_length",&tagger.cme_mu_length);
  reader.AddVariable("cme_length",&tagger.cme_length);
  reader.AddVariable("cme_angle_beam",&tagger.cme_angle_beam);
  reader.AddVariable("anc_angle",&tagger.anc_angle);
  reader.AddVariable("anc_max_angle",&tagger.anc_max_angle);
  reader.AddVariable("anc_max_length",&tagger.anc_max_length);
  reader.AddVariable("anc_acc_forward_length",&tagger.anc_acc_forward_length);
  reader.AddVariable("anc_acc_backward_length",&tagger.anc_acc_backward_length);
  reader.AddVariable("anc_acc_forward_length1",&tagger.anc_acc_forward_length1);
  reader.AddVariable("anc_shower_main_length",&tagger.anc_shower_main_length);
  reader.AddVariable("anc_shower_total_length",&tagger.anc_shower_total_length);
  reader.AddVariable("anc_flag_main_outside",&tagger.anc_flag_main_outside);
  reader.AddVariable("gap_flag_prolong_u",&tagger.gap_flag_prolong_u);
  reader.AddVariable("gap_flag_prolong_v",&tagger.gap_flag_prolong_v);
  reader.AddVariable("gap_flag_prolong_w",&tagger.gap_flag_prolong_w);
  reader.AddVariable("gap_flag_parallel",&tagger.gap_flag_parallel);
  reader.AddVariable("gap_n_points",&tagger.gap_n_points);
  reader.AddVariable("gap_n_bad",&tagger.gap_n_bad);
  reader.AddVariable("gap_energy",&tagger.gap_energy);
  reader.AddVariable("gap_num_valid_tracks",&tagger.gap_num_valid_tracks);
  reader.AddVariable("gap_flag_single_shower",&tagger.gap_flag_single_shower);
  reader.AddVariable("hol_1_n_valid_tracks",&tagger.hol_1_n_valid_tracks);
  reader.AddVariable("hol_1_min_angle",&tagger.hol_1_min_angle);
  reader.AddVariable("hol_1_energy",&tagger.hol_1_energy);
  reader.AddVariable("hol_1_min_length",&tagger.hol_1_min_length);
  reader.AddVariable("hol_2_min_angle",&tagger.hol_2_min_angle);
  reader.AddVariable("hol_2_medium_dQ_dx",&tagger.hol_2_medium_dQ_dx);
  reader.AddVariable("hol_2_ncount",&tagger.hol_2_ncount);
  reader.AddVariable("lol_3_angle_beam",&tagger.lol_3_angle_beam);
  reader.AddVariable("lol_3_n_valid_tracks",&tagger.lol_3_n_valid_tracks);
  reader.AddVariable("lol_3_min_angle",&tagger.lol_3_min_angle);
  reader.AddVariable("lol_3_vtx_n_segs",&tagger.lol_3_vtx_n_segs);
  reader.AddVariable("lol_3_shower_main_length",&tagger.lol_3_shower_main_length);
  reader.AddVariable("lol_3_n_out",&tagger.lol_3_n_out);
  reader.AddVariable("lol_3_n_sum",&tagger.lol_3_n_sum);
  reader.AddVariable("hol_1_flag_all_shower",&tagger.hol_1_flag_all_shower); // naming issue
  reader.AddVariable("mgo_energy",&tagger.mgo_energy);
  reader.AddVariable("mgo_max_energy",&tagger.mgo_max_energy);
  reader.AddVariable("mgo_total_energy",&tagger.mgo_total_energy);
  reader.AddVariable("mgo_n_showers",&tagger.mgo_n_showers);
  reader.AddVariable("mgo_max_energy_1",&tagger.mgo_max_energy_1);
  reader.AddVariable("mgo_max_energy_2",&tagger.mgo_max_energy_2);
  reader.AddVariable("mgo_total_other_energy",&tagger.mgo_total_other_energy);
  reader.AddVariable("mgo_n_total_showers",&tagger.mgo_n_total_showers);
  reader.AddVariable("mgo_total_other_energy_1",&tagger.mgo_total_other_energy_1);
  reader.AddVariable("mgt_flag_single_shower",&tagger.mgt_flag_single_shower);
  reader.AddVariable("mgt_max_energy",&tagger.mgt_max_energy);
  reader.AddVariable("mgt_total_other_energy",&tagger.mgt_total_other_energy);
  reader.AddVariable("mgt_max_energy_1",&tagger.mgt_max_energy_1);
  reader.AddVariable("mgt_e_indirect_max_energy",&tagger.mgt_e_indirect_max_energy);
  reader.AddVariable("mgt_e_direct_max_energy",&tagger.mgt_e_direct_max_energy);
  reader.AddVariable("mgt_n_direct_showers",&tagger.mgt_n_direct_showers);
  reader.AddVariable("mgt_e_direct_total_energy",&tagger.mgt_e_direct_total_energy);
  reader.AddVariable("mgt_flag_indirect_max_pio",&tagger.mgt_flag_indirect_max_pio);
  reader.AddVariable("mgt_e_indirect_total_energy",&tagger.mgt_e_indirect_total_energy);
  reader.AddVariable("mip_quality_energy",&tagger.mip_quality_energy);
  reader.AddVariable("mip_quality_overlap",&tagger.mip_quality_overlap);
  reader.AddVariable("mip_quality_n_showers",&tagger.mip_quality_n_showers);
  reader.AddVariable("mip_quality_n_tracks",&tagger.mip_quality_n_tracks);
  reader.AddVariable("mip_quality_flag_inside_pi0",&tagger.mip_quality_flag_inside_pi0);
  reader.AddVariable("mip_quality_n_pi0_showers",&tagger.mip_quality_n_pi0_showers);
  reader.AddVariable("mip_quality_shortest_length",&tagger.mip_quality_shortest_length);
  reader.AddVariable("mip_quality_acc_length",&tagger.mip_quality_acc_length);
  reader.AddVariable("mip_quality_shortest_angle",&tagger.mip_quality_shortest_angle);
  reader.AddVariable("mip_quality_flag_proton",&tagger.mip_quality_flag_proton);
  reader.AddVariable("br1_1_shower_type",&tagger.br1_1_shower_type);
  reader.AddVariable("br1_1_vtx_n_segs",&tagger.br1_1_vtx_n_segs);
  reader.AddVariable("br1_1_energy",&tagger.br1_1_energy);
  reader.AddVariable("br1_1_n_segs",&tagger.br1_1_n_segs);
  reader.AddVariable("br1_1_flag_sg_topology",&tagger.br1_1_flag_sg_topology);
  reader.AddVariable("br1_1_flag_sg_trajectory",&tagger.br1_1_flag_sg_trajectory);
  reader.AddVariable("br1_1_sg_length",&tagger.br1_1_sg_length);
  reader.AddVariable("br1_2_n_connected",&tagger.br1_2_n_connected);
  reader.AddVariable("br1_2_max_length",&tagger.br1_2_max_length);
  reader.AddVariable("br1_2_n_connected_1",&tagger.br1_2_n_connected_1);
  reader.AddVariable("br1_2_n_shower_segs",&tagger.br1_2_n_shower_segs);
  reader.AddVariable("br1_2_max_length_ratio",&tagger.br1_2_max_length_ratio);
  reader.AddVariable("br1_2_shower_length",&tagger.br1_2_shower_length);
  reader.AddVariable("br1_3_n_connected_p",&tagger.br1_3_n_connected_p);
  reader.AddVariable("br1_3_max_length_p",&tagger.br1_3_max_length_p);
  reader.AddVariable("br1_3_n_shower_main_segs",&tagger.br1_3_n_shower_main_segs);
  reader.AddVariable("br3_1_energy",&tagger.br3_1_energy);
  reader.AddVariable("br3_1_n_shower_segments",&tagger.br3_1_n_shower_segments);
  reader.AddVariable("br3_1_sg_flag_trajectory",&tagger.br3_1_sg_flag_trajectory);
  reader.AddVariable("br3_1_sg_direct_length",&tagger.br3_1_sg_direct_length);
  reader.AddVariable("br3_1_sg_length",&tagger.br3_1_sg_length);
  reader.AddVariable("br3_1_total_main_length",&tagger.br3_1_total_main_length);
  reader.AddVariable("br3_1_total_length",&tagger.br3_1_total_length);
  reader.AddVariable("br3_1_iso_angle",&tagger.br3_1_iso_angle);
  reader.AddVariable("br3_1_sg_flag_topology",&tagger.br3_1_sg_flag_topology);
  reader.AddVariable("br3_2_n_ele",&tagger.br3_2_n_ele);
  reader.AddVariable("br3_2_n_other",&tagger.br3_2_n_other);
  reader.AddVariable("br3_2_other_fid",&tagger.br3_2_other_fid);
  reader.AddVariable("br3_4_acc_length",&tagger.br3_4_acc_length);
  reader.AddVariable("br3_4_total_length",&tagger.br3_4_total_length);
  reader.AddVariable("br3_7_min_angle",&tagger.br3_7_min_angle);
  reader.AddVariable("br3_8_max_dQ_dx",&tagger.br3_8_max_dQ_dx);
  reader.AddVariable("br3_8_n_main_segs",&tagger.br3_8_n_main_segs);
  reader.AddVariable("vis_1_n_vtx_segs",&tagger.vis_1_n_vtx_segs);
  reader.AddVariable("vis_1_energy",&tagger.vis_1_energy);
  reader.AddVariable("vis_1_num_good_tracks",&tagger.vis_1_num_good_tracks);
  reader.AddVariable("vis_1_max_angle",&tagger.vis_1_max_angle);
  reader.AddVariable("vis_1_max_shower_angle",&tagger.vis_1_max_shower_angle);
  reader.AddVariable("vis_1_tmp_length1",&tagger.vis_1_tmp_length1);
  reader.AddVariable("vis_1_tmp_length2",&tagger.vis_1_tmp_length2);
  reader.AddVariable("vis_2_n_vtx_segs",&tagger.vis_2_n_vtx_segs);
  reader.AddVariable("vis_2_min_angle",&tagger.vis_2_min_angle);
  reader.AddVariable("vis_2_min_weak_track",&tagger.vis_2_min_weak_track);
  reader.AddVariable("vis_2_angle_beam",&tagger.vis_2_angle_beam);
  reader.AddVariable("vis_2_min_angle1",&tagger.vis_2_min_angle1);
  reader.AddVariable("vis_2_iso_angle1",&tagger.vis_2_iso_angle1);
  reader.AddVariable("vis_2_min_medium_dQ_dx",&tagger.vis_2_min_medium_dQ_dx);
  reader.AddVariable("vis_2_min_length",&tagger.vis_2_min_length);
  reader.AddVariable("vis_2_sg_length",&tagger.vis_2_sg_length);
  reader.AddVariable("vis_2_max_angle",&tagger.vis_2_max_angle);
  reader.AddVariable("vis_2_max_weak_track",&tagger.vis_2_max_weak_track);
  reader.AddVariable("pio_1_mass",&tagger.pio_1_mass);
  reader.AddVariable("pio_1_pio_type",&tagger.pio_1_pio_type);
  reader.AddVariable("pio_1_energy_1",&tagger.pio_1_energy_1);
  reader.AddVariable("pio_1_energy_2",&tagger.pio_1_energy_2);
  reader.AddVariable("pio_1_dis_1",&tagger.pio_1_dis_1);
  reader.AddVariable("pio_1_dis_2",&tagger.pio_1_dis_2);
  reader.AddVariable("pio_mip_id",&tagger.pio_mip_id);
  reader.AddVariable("stem_dir_flag_single_shower",&tagger.stem_dir_flag_single_shower);
  reader.AddVariable("stem_dir_angle",&tagger.stem_dir_angle);
  reader.AddVariable("stem_dir_energy",&tagger.stem_dir_energy);
  reader.AddVariable("stem_dir_angle1",&tagger.stem_dir_angle1);
  reader.AddVariable("stem_dir_angle2",&tagger.stem_dir_angle2);
  reader.AddVariable("stem_dir_angle3",&tagger.stem_dir_angle3);
  reader.AddVariable("stem_dir_ratio",&tagger.stem_dir_ratio);
  reader.AddVariable("br2_num_valid_tracks",&tagger.br2_num_valid_tracks);
  reader.AddVariable("br2_n_shower_main_segs",&tagger.br2_n_shower_main_segs);
  reader.AddVariable("br2_max_angle",&tagger.br2_max_angle);
  reader.AddVariable("br2_sg_length",&tagger.br2_sg_length);
  reader.AddVariable("br2_flag_sg_trajectory",&tagger.br2_flag_sg_trajectory);
  reader.AddVariable("stem_len_energy",&tagger.stem_len_energy);
  reader.AddVariable("stem_len_length",&tagger.stem_len_length);
  reader.AddVariable("stem_len_flag_avoid_muon_check",&tagger.stem_len_flag_avoid_muon_check);
  reader.AddVariable("stem_len_num_daughters",&tagger.stem_len_num_daughters);
  reader.AddVariable("stem_len_daughter_length",&tagger.stem_len_daughter_length);
  reader.AddVariable("brm_n_mu_segs",&tagger.brm_n_mu_segs);
  reader.AddVariable("brm_Ep",&tagger.brm_Ep);
  reader.AddVariable("brm_acc_length",&tagger.brm_acc_length);
  reader.AddVariable("brm_shower_total_length",&tagger.brm_shower_total_length);
  reader.AddVariable("brm_connected_length",&tagger.brm_connected_length);
  reader.AddVariable("brm_n_size",&tagger.brm_n_size);
  reader.AddVariable("brm_n_shower_main_segs",&tagger.brm_n_shower_main_segs);
  reader.AddVariable("brm_n_mu_main",&tagger.brm_n_mu_main);
  reader.AddVariable("lem_shower_main_length",&tagger.lem_shower_main_length);
  reader.AddVariable("lem_n_3seg",&tagger.lem_n_3seg);
  reader.AddVariable("lem_e_charge",&tagger.lem_e_charge);
  reader.AddVariable("lem_e_dQdx",&tagger.lem_e_dQdx);
  reader.AddVariable("lem_shower_num_main_segs",&tagger.lem_shower_num_main_segs);
  reader.AddVariable("brm_acc_direct_length",&tagger.brm_acc_direct_length); // naming issue
  reader.AddVariable("stw_1_energy",&tagger.stw_1_energy);
  reader.AddVariable("stw_1_dis",&tagger.stw_1_dis);
  reader.AddVariable("stw_1_dQ_dx",&tagger.stw_1_dQ_dx);
  reader.AddVariable("stw_1_flag_single_shower",&tagger.stw_1_flag_single_shower);
  reader.AddVariable("stw_1_n_pi0",&tagger.stw_1_n_pi0);
  reader.AddVariable("stw_1_num_valid_tracks",&tagger.stw_1_num_valid_tracks);
  reader.AddVariable("spt_shower_main_length",&tagger.spt_shower_main_length);
  reader.AddVariable("spt_shower_total_length",&tagger.spt_shower_total_length);
  reader.AddVariable("spt_angle_beam",&tagger.spt_angle_beam);
  reader.AddVariable("spt_angle_vertical",&tagger.spt_angle_vertical);
  reader.AddVariable("spt_max_dQ_dx",&tagger.spt_max_dQ_dx);
  reader.AddVariable("spt_angle_beam_1",&tagger.spt_angle_beam_1);
  reader.AddVariable("spt_angle_drift",&tagger.spt_angle_drift);
  reader.AddVariable("spt_angle_drift_1",&tagger.spt_angle_drift_1);
  reader.AddVariable("spt_num_valid_tracks",&tagger.spt_num_valid_tracks);
  reader.AddVariable("spt_n_vtx_segs",&tagger.spt_n_vtx_segs);
  reader.AddVariable("spt_max_length",&tagger.spt_max_length);
  reader.AddVariable("mip_energy",&tagger.mip_energy);
  reader.AddVariable("mip_n_end_reduction",&tagger.mip_n_end_reduction);
  reader.AddVariable("mip_n_first_mip",&tagger.mip_n_first_mip);
  reader.AddVariable("mip_n_first_non_mip",&tagger.mip_n_first_non_mip);
  reader.AddVariable("mip_n_first_non_mip_1",&tagger.mip_n_first_non_mip_1);
  reader.AddVariable("mip_n_first_non_mip_2",&tagger.mip_n_first_non_mip_2);
  reader.AddVariable("mip_vec_dQ_dx_0",&tagger.mip_vec_dQ_dx_0);
  reader.AddVariable("mip_vec_dQ_dx_1",&tagger.mip_vec_dQ_dx_1);
  reader.AddVariable("mip_max_dQ_dx_sample",&tagger.mip_max_dQ_dx_sample);
  reader.AddVariable("mip_n_below_threshold",&tagger.mip_n_below_threshold);
  reader.AddVariable("mip_n_below_zero",&tagger.mip_n_below_zero);
  reader.AddVariable("mip_n_lowest",&tagger.mip_n_lowest);
  reader.AddVariable("mip_n_highest",&tagger.mip_n_highest);
  reader.AddVariable("mip_lowest_dQ_dx",&tagger.mip_lowest_dQ_dx);
  reader.AddVariable("mip_highest_dQ_dx",&tagger.mip_highest_dQ_dx);
  reader.AddVariable("mip_medium_dQ_dx",&tagger.mip_medium_dQ_dx);
  reader.AddVariable("mip_stem_length",&tagger.mip_stem_length);
  reader.AddVariable("mip_length_main",&tagger.mip_length_main);
  reader.AddVariable("mip_length_total",&tagger.mip_length_total);
  reader.AddVariable("mip_angle_beam",&tagger.mip_angle_beam);
  reader.AddVariable("mip_iso_angle",&tagger.mip_iso_angle);
  reader.AddVariable("mip_n_vertex",&tagger.mip_n_vertex);
  reader.AddVariable("mip_n_good_tracks",&tagger.mip_n_good_tracks);
  reader.AddVariable("mip_E_indirect_max_energy",&tagger.mip_E_indirect_max_energy);
  reader.AddVariable("mip_flag_all_above",&tagger.mip_flag_all_above);
  reader.AddVariable("mip_min_dQ_dx_5",&tagger.mip_min_dQ_dx_5);
  reader.AddVariable("mip_n_other_vertex",&tagger.mip_n_other_vertex);
  reader.AddVariable("mip_n_stem_size",&tagger.mip_n_stem_size);
  reader.AddVariable("mip_flag_stem_trajectory",&tagger.mip_flag_stem_trajectory);
  reader.AddVariable("mip_min_dis",&tagger.mip_min_dis);
  reader.AddVariable("mip_vec_dQ_dx_2",&tagger.mip_vec_dQ_dx_2);
  reader.AddVariable("mip_vec_dQ_dx_3",&tagger.mip_vec_dQ_dx_3);
  reader.AddVariable("mip_vec_dQ_dx_4",&tagger.mip_vec_dQ_dx_4);
  reader.AddVariable("mip_vec_dQ_dx_5",&tagger.mip_vec_dQ_dx_5);
  reader.AddVariable("mip_vec_dQ_dx_6",&tagger.mip_vec_dQ_dx_6);
  reader.AddVariable("mip_vec_dQ_dx_7",&tagger.mip_vec_dQ_dx_7);
  reader.AddVariable("mip_vec_dQ_dx_8",&tagger.mip_vec_dQ_dx_8);
  reader.AddVariable("mip_vec_dQ_dx_9",&tagger.mip_vec_dQ_dx_9);
  reader.AddVariable("mip_vec_dQ_dx_10",&tagger.mip_vec_dQ_dx_10);
  reader.AddVariable("mip_vec_dQ_dx_11",&tagger.mip_vec_dQ_dx_11);
  reader.AddVariable("mip_vec_dQ_dx_12",&tagger.mip_vec_dQ_dx_12);
  reader.AddVariable("mip_vec_dQ_dx_13",&tagger.mip_vec_dQ_dx_13);
  reader.AddVariable("mip_vec_dQ_dx_14",&tagger.mip_vec_dQ_dx_14);
  reader.AddVariable("mip_vec_dQ_dx_15",&tagger.mip_vec_dQ_dx_15);
  reader.AddVariable("mip_vec_dQ_dx_16",&tagger.mip_vec_dQ_dx_16);
  reader.AddVariable("mip_vec_dQ_dx_17",&tagger.mip_vec_dQ_dx_17);
  reader.AddVariable("mip_vec_dQ_dx_18",&tagger.mip_vec_dQ_dx_18);
  reader.AddVariable("mip_vec_dQ_dx_19",&tagger.mip_vec_dQ_dx_19);
  reader.AddVariable("br3_3_score",&tagger.br3_3_score);
  reader.AddVariable("br3_5_score",&tagger.br3_5_score);
  reader.AddVariable("br3_6_score",&tagger.br3_6_score);
  reader.AddVariable("pio_2_score",&tagger.pio_2_score);
  reader.AddVariable("stw_2_score",&tagger.stw_2_score);
  reader.AddVariable("stw_3_score",&tagger.stw_3_score);
  reader.AddVariable("stw_4_score",&tagger.stw_4_score);
  reader.AddVariable("sig_1_score",&tagger.sig_1_score);
  reader.AddVariable("sig_2_score",&tagger.sig_2_score);
  reader.AddVariable("lol_1_score",&tagger.lol_1_score);
  reader.AddVariable("lol_2_score",&tagger.lol_2_score);
  reader.AddVariable("tro_1_score",&tagger.tro_1_score);
  reader.AddVariable("tro_2_score",&tagger.tro_2_score);
  reader.AddVariable("tro_4_score",&tagger.tro_4_score);
  reader.AddVariable("tro_5_score",&tagger.tro_5_score);
  reader.AddVariable("br4_1_shower_main_length",&tagger.br4_1_shower_main_length);
  reader.AddVariable("br4_1_shower_total_length",&tagger.br4_1_shower_total_length);
  reader.AddVariable("br4_1_min_dis",&tagger.br4_1_min_dis);
  reader.AddVariable("br4_1_energy",&tagger.br4_1_energy);
  reader.AddVariable("br4_1_flag_avoid_muon_check",&tagger.br4_1_flag_avoid_muon_check);
  reader.AddVariable("br4_1_n_vtx_segs",&tagger.br4_1_n_vtx_segs);
  reader.AddVariable("br4_2_ratio_45",&tagger.br4_2_ratio_45);
  reader.AddVariable("br4_2_ratio_35",&tagger.br4_2_ratio_35);
  reader.AddVariable("br4_2_ratio_25",&tagger.br4_2_ratio_25);
  reader.AddVariable("br4_2_ratio_15",&tagger.br4_2_ratio_15);
  reader.AddVariable("br4_2_ratio1_45",&tagger.br4_2_ratio1_45);
  reader.AddVariable("br4_2_ratio1_35",&tagger.br4_2_ratio1_35);
  reader.AddVariable("br4_2_ratio1_25",&tagger.br4_2_ratio1_25);
  reader.AddVariable("br4_2_ratio1_15",&tagger.br4_2_ratio1_15);
  reader.AddVariable("br4_2_iso_angle",&tagger.br4_2_iso_angle);
  reader.AddVariable("br4_2_iso_angle1",&tagger.br4_2_iso_angle1);
  reader.AddVariable("br4_2_angle",&tagger.br4_2_angle);
  reader.AddVariable("tro_3_stem_length",&tagger.tro_3_stem_length);
  reader.AddVariable("tro_3_n_muon_segs",&tagger.tro_3_n_muon_segs);
  reader.AddVariable("br4_1_n_main_segs",&tagger.br4_1_n_main_segs); // naming issue

  reader.BookMVA( "MyBDT", "./weights/XGB_nue_seed2_0923.xml");
  //  reader.BookMVA( "MyBDT", "./weights/xgboost_set8seed7_kaicheng_0819.xml");

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


  TMVA::Reader reader_numu;
  
  reader_numu.AddVariable("numu_cc_flag_3", &tagger.numu_cc_flag_3);
  reader_numu.AddVariable("numu_cc_3_particle_type", &tagger.numu_cc_3_particle_type);
  reader_numu.AddVariable("numu_cc_3_max_length", &tagger.numu_cc_3_max_length);
  reader_numu.AddVariable("numu_cc_3_track_length",&tagger.numu_cc_3_acc_track_length);
  reader_numu.AddVariable("numu_cc_3_max_length_all",&tagger.numu_cc_3_max_length_all);
  reader_numu.AddVariable("numu_cc_3_max_muon_length",&tagger.numu_cc_3_max_muon_length);
  reader_numu.AddVariable("numu_cc_3_n_daughter_tracks",&tagger.numu_cc_3_n_daughter_tracks);
  reader_numu.AddVariable("numu_cc_3_n_daughter_all",&tagger.numu_cc_3_n_daughter_all);
  reader_numu.AddVariable("cosmict_flag_2", &tagger.cosmict_flag_2);
  reader_numu.AddVariable("cosmict_2_filled", &tagger.cosmict_2_filled);
  reader_numu.AddVariable("cosmict_2_particle_type",&tagger.cosmict_2_particle_type);
  reader_numu.AddVariable("cosmict_2_n_muon_tracks",&tagger.cosmict_2_n_muon_tracks);
  reader_numu.AddVariable("cosmict_2_total_shower_length",&tagger.cosmict_2_total_shower_length);
  reader_numu.AddVariable("cosmict_2_flag_inside",&tagger.cosmict_2_flag_inside);
  reader_numu.AddVariable("cosmict_2_angle_beam",&tagger.cosmict_2_angle_beam);
  reader_numu.AddVariable("cosmict_2_flag_dir_weak", &tagger.cosmict_2_flag_dir_weak);
  reader_numu.AddVariable("cosmict_2_dQ_dx_end", &tagger.cosmict_2_dQ_dx_end);
  reader_numu.AddVariable("cosmict_2_dQ_dx_front", &tagger.cosmict_2_dQ_dx_front);
  reader_numu.AddVariable("cosmict_2_theta", &tagger.cosmict_2_theta);
  reader_numu.AddVariable("cosmict_2_phi", &tagger.cosmict_2_phi);
  reader_numu.AddVariable("cosmict_2_valid_tracks", &tagger.cosmict_2_valid_tracks);
  reader_numu.AddVariable("cosmict_flag_4", &tagger.cosmict_flag_4);
  reader_numu.AddVariable("cosmict_4_filled", &tagger.cosmict_4_filled);
  reader_numu.AddVariable("cosmict_4_flag_inside", &tagger.cosmict_4_flag_inside);
  reader_numu.AddVariable("cosmict_4_angle_beam", &tagger.cosmict_4_angle_beam);
  reader_numu.AddVariable("cosmict_4_connected_showers", &tagger.cosmict_4_connected_showers);
  reader_numu.AddVariable("cosmict_flag_3", &tagger.cosmict_flag_3);
  reader_numu.AddVariable("cosmict_3_filled", &tagger.cosmict_3_filled);
  reader_numu.AddVariable("cosmict_3_flag_inside", &tagger.cosmict_3_flag_inside);
  reader_numu.AddVariable("cosmict_3_angle_beam", &tagger.cosmict_3_angle_beam);
  reader_numu.AddVariable("cosmict_3_flag_dir_weak", &tagger.cosmict_3_flag_dir_weak);
  reader_numu.AddVariable("cosmict_3_dQ_dx_end", &tagger.cosmict_3_dQ_dx_end);
  reader_numu.AddVariable("cosmict_3_dQ_dx_front", &tagger.cosmict_3_dQ_dx_front);
  reader_numu.AddVariable("cosmict_3_theta", &tagger.cosmict_3_theta);
  reader_numu.AddVariable("cosmict_3_phi", &tagger.cosmict_3_phi);
  reader_numu.AddVariable("cosmict_3_valid_tracks", &tagger.cosmict_3_valid_tracks);
  reader_numu.AddVariable("cosmict_flag_5", &tagger.cosmict_flag_5);
  reader_numu.AddVariable("cosmict_5_filled", &tagger.cosmict_5_filled);
  reader_numu.AddVariable("cosmict_5_flag_inside", &tagger.cosmict_5_flag_inside);
  reader_numu.AddVariable("cosmict_5_angle_beam", &tagger.cosmict_5_angle_beam);
  reader_numu.AddVariable("cosmict_5_connected_showers", &tagger.cosmict_5_connected_showers);
  reader_numu.AddVariable("cosmict_flag_6", &tagger.cosmict_flag_6);
  reader_numu.AddVariable("cosmict_6_filled", &tagger.cosmict_6_filled);
  reader_numu.AddVariable("cosmict_6_flag_dir_weak", &tagger.cosmict_6_flag_dir_weak);
  reader_numu.AddVariable("cosmict_6_flag_inside", &tagger.cosmict_6_flag_inside);
  reader_numu.AddVariable("cosmict_6_angle", &tagger.cosmict_6_angle);
  reader_numu.AddVariable("cosmict_flag_7", &tagger.cosmict_flag_7);
  reader_numu.AddVariable("cosmict_7_filled", &tagger.cosmict_7_filled);
  reader_numu.AddVariable("cosmict_7_flag_sec", &tagger.cosmict_7_flag_sec);
  reader_numu.AddVariable("cosmict_7_n_muon_tracks", &tagger.cosmict_7_n_muon_tracks);
  reader_numu.AddVariable("cosmict_7_total_shower_length", &tagger.cosmict_7_total_shower_length);
  reader_numu.AddVariable("cosmict_7_flag_inside", &tagger.cosmict_7_flag_inside);
  reader_numu.AddVariable("cosmict_7_angle_beam", &tagger.cosmict_7_angle_beam);
  reader_numu.AddVariable("cosmict_7_flag_dir_weak", &tagger.cosmict_7_flag_dir_weak);
  reader_numu.AddVariable("cosmict_7_dQ_dx_end", &tagger.cosmict_7_dQ_dx_end);
  reader_numu.AddVariable("cosmict_7_dQ_dx_front", &tagger.cosmict_7_dQ_dx_front);
  reader_numu.AddVariable("cosmict_7_theta", &tagger.cosmict_7_theta);
  reader_numu.AddVariable("cosmict_7_phi", &tagger.cosmict_7_phi);
  reader_numu.AddVariable("cosmict_flag_8", &tagger.cosmict_flag_8);
  reader_numu.AddVariable("cosmict_8_filled", &tagger.cosmict_8_filled);
  reader_numu.AddVariable("cosmict_8_flag_out", &tagger.cosmict_8_flag_out);
  reader_numu.AddVariable("cosmict_8_muon_length", &tagger.cosmict_8_muon_length);
  reader_numu.AddVariable("cosmict_8_acc_length", &tagger.cosmict_8_acc_length);
  reader_numu.AddVariable("cosmict_flag_9", &tagger.cosmict_flag_9);
  reader_numu.AddVariable("cosmic_flag", &tagger.cosmic_flag);
  reader_numu.AddVariable("cosmic_filled", &tagger.cosmic_filled);
  reader_numu.AddVariable("cosmict_flag", &tagger.cosmict_flag);
  reader_numu.AddVariable("numu_cc_flag", &tagger.numu_cc_flag);
  reader_numu.AddVariable("cosmict_flag_1", &tagger.cosmict_flag_1);
  reader_numu.AddVariable("kine_reco_Enu",&tagger.kine_reco_Enu);
  reader_numu.AddVariable("match_isFC",&tagger.match_isFC);
  reader_numu.AddVariable("cosmict_10_score", &tagger.cosmict_10_score);
  reader_numu.AddVariable("numu_1_score", &tagger.numu_1_score);
  reader_numu.AddVariable("numu_2_score", &tagger.numu_2_score);

  reader_numu.BookMVA( "MyBDT", "weights/numu_scalars_scores_0923.xml");


  std::map<std::pair<int, int>, int> map_rs_n;
  std::map<std::pair<int, int>, std::set<int> > map_rs_f1p5; // Reco 1.5
  std::map<std::pair<int, int>, std::set<int> > map_rs_f2stm; // Reco2 stm
  std::map<std::pair<int, int>, std::set<int> > map_rs_f2pr; // Reco2 Pattern recognition
  
  T_eval->SetBranchStatus("*",0);
  T_eval->SetBranchStatus("stm_eventtype",1); 
  T_eval->SetBranchStatus("stm_lowenergy",1); 
  T_eval->SetBranchStatus("stm_LM",1); 
  T_eval->SetBranchStatus("stm_TGM",1); 
  T_eval->SetBranchStatus("stm_STM",1); 
  T_eval->SetBranchStatus("stm_FullDead",1); 
  T_eval->SetBranchStatus("stm_clusterlength",1);
  T_eval->SetBranchStatus("match_found",1);
  T_eval->SetBranchStatus("run",1); 
  T_eval->SetBranchStatus("subrun",1); 
  T_eval->SetBranchStatus("event",1);
  T_eval->SetBranchStatus("file_type",1);

  if (T_eval->GetBranch("match_found_asInt")){
    T_eval->SetBranchStatus("match_found_asInt",1);
  }
  
  T_BDTvars->SetBranchStatus("*",0);
  T_BDTvars->SetBranchStatus("numu_cc_flag",1);

  std::set<std::pair<int,int> > remove_set;
  
  bool flag_presel = false;
  for (Int_t i=0;i!=T_eval->GetEntries();i++){
    
    T_eval->GetEntry(i);
    T_BDTvars->GetEntry(i);

    if (flag_check_run_subrun){
      if (flag_use_global_file_type){
	(*eval.file_type) = global_file_type;
      }
      auto it1 = map_type_run_subrun.find(*eval.file_type);
      
      if (it1 != map_type_run_subrun.end()){
	
	// hack for now ...
	if (it1->second.find(std::make_pair(eval.run, eval.subrun)) != it1->second.end()) {
	  remove_set.insert(std::make_pair(eval.run, eval.subrun));
	  continue;
	}
      }

      //std::cout << flag_use_global_file_type << " " << *eval.file_type  << " " << eval.run << " " << eval.subrun << " " << remove_set.size() << std::endl;
    }

    
    
    int tmp_match_found = eval.match_found;
    if (eval.is_match_found_int) tmp_match_found = eval.match_found_asInt;
    
    flag_presel = false;
    if (tmp_match_found != 0 && eval.stm_eventtype != 0 && eval.stm_lowenergy ==0 && eval.stm_LM ==0 && eval.stm_TGM ==0 && eval.stm_STM==0 && eval.stm_FullDead == 0 && eval.stm_clusterlength >0) {
      flag_presel = true; // preselection ...
    }
    
    map_rs_n[std::make_pair(eval.run, eval.subrun)] ++;
    if (tmp_match_found == -1) map_rs_f1p5[std::make_pair(eval.run, eval.subrun)].insert(eval.event);
    if (tmp_match_found == 1 && eval.stm_lowenergy == -1) map_rs_f2stm[std::make_pair(eval.run, eval.subrun)].insert(eval.event);
    if (flag_presel && tagger.numu_cc_flag == -1) map_rs_f2pr[std::make_pair(eval.run, eval.subrun)].insert(eval.event);
  }


  
  for (auto it = map_rs_f1p5.begin(); it!= map_rs_f1p5.end(); it++){
    if ( map_rs_f1p5[it->first].size()  > map_rs_n[it->first] * fail_percentage && map_rs_f1p5[it->first].size() != 1
	 || map_rs_f1p5[it->first].size()  > map_rs_n[it->first] * 0.33 && map_rs_f1p5[it->first].size() == 1)
      remove_set.insert(it->first);
  }
  for (auto it = map_rs_f2stm.begin(); it != map_rs_f2stm.end(); it++){
    if (map_rs_f2stm[it->first].size() > map_rs_n[it->first]* fail_percentage && map_rs_f2stm[it->first].size() != 1
	|| map_rs_f2stm[it->first].size() > map_rs_n[it->first]* 0.33 && map_rs_f2stm[it->first].size() == 1)
      remove_set.insert(it->first);
  }
  for (auto it = map_rs_f2pr.begin(); it != map_rs_f2pr.end(); it++){
    if (map_rs_f2pr[it->first].size() > map_rs_n[it->first] * fail_percentage && map_rs_f2pr[it->first].size() != 1 || map_rs_f2pr[it->first].size() > map_rs_n[it->first] * 0.33 && map_rs_f2pr[it->first].size() == 1
	)
      remove_set.insert(it->first);
  }

  //  std::cout << remove_set.size() << std::endl;
  

  T_eval->SetBranchStatus("*",1);
  T_BDTvars->SetBranchStatus("*",1);
  
  for (int i=0;i!=T_BDTvars->GetEntries();i++){
    eval.weight_change = false;
    T_BDTvars->GetEntry(i);
    T_eval->GetEntry(i); tagger.match_isFC = eval.match_isFC;
    T_KINEvars->GetEntry(i); tagger.kine_reco_Enu = kine.kine_reco_Enu;
    T_PFeval->GetEntry(i);

    if (remove_set.find(std::make_pair(eval.run, eval.subrun)) != remove_set.end()) continue;

    
    
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
    tagger.nue_score       = cal_bdts_xgboost( tagger,  reader);
    
    // BDT calculations
    tagger.numu_1_score = cal_numu_1_bdt(-0.4, tagger, reader_numu_1, numu_cc_flag_1,
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
    if (std::isnan(tagger.cosmict_7_angle_beam)) tagger.cosmict_7_angle_beam = 0;
    if (std::isnan(tagger.cosmict_7_theta)) tagger.cosmict_7_theta = 0;
    if (std::isnan(tagger.cosmict_7_phi)) tagger.cosmict_7_phi = 0;
    
    tagger.numu_score = cal_numu_bdts_xgboost(tagger,reader_numu); 

    // limit the cut val ...
    if (std::isnan(eval.weight_spline) || std::isinf(eval.weight_spline) ||
	std::isnan(eval.weight_cv) || std::isinf(eval.weight_cv) ||
	eval.weight_spline * eval.weight_cv <=0 ||
	eval.weight_spline * eval.weight_cv > weight_cut_val){
      eval.weight_spline = 1;
      eval.weight_cv = 1;
      eval.weight_change = true;
    }
    
    
    t4->Fill();
    t1->Fill();
    t3->Fill();
    t5->Fill();
  }

  for (Int_t i=0;i!=T_pot->GetEntries();i++){
    T_pot->GetEntry(i);
    
    if (remove_set.find(std::make_pair(pot.runNo, pot.subRunNo)) != remove_set.end()) continue;
    

    t2->Fill();
  }


  for (auto it = remove_set.begin(); it!= remove_set.end(); it++){
    std::cout <<"remove  run:" << it->first << " subrun:" << it->second << std::endl;
  }
  
  file2->Write();
  file2->Close();
  
  return 0;

  
  
}
