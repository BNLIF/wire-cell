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

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

using namespace std;

TFile *input = 0;
TFile *output = 0;
TMVA::Factory *factory = 0;
TMVA::DataLoader *dataloader = 0;

struct numuInfo
{
  // original cosmic 1 
  float cosmict_flag_1;
  
  // numu 1
  std::vector<float> *numu_cc_flag_1;// = new std::vector<float>;
  std::vector<float> *numu_cc_1_particle_type;//= new std::vector<float>;
  std::vector<float> *numu_cc_1_length;//= new std::vector<float>;
  std::vector<float> *numu_cc_1_medium_dQ_dx;//= new std::vector<float>;
  std::vector<float> *numu_cc_1_dQ_dx_cut;//= new std::vector<float>;
  std::vector<float> *numu_cc_1_direct_length;//= new std::vector<float>;
  std::vector<float> *numu_cc_1_n_daughter_tracks;//= new std::vector<float>;
  std::vector<float> *numu_cc_1_n_daughter_all;//= new std::vector<float>;

  //numu 2
  std::vector<float> *numu_cc_flag_2;// = new std::vector<float>;
  std::vector<float> *numu_cc_2_length;// = new std::vector<float>;
  std::vector<float> *numu_cc_2_total_length;// = new std::vector<float>;
  std::vector<float> *numu_cc_2_n_daughter_tracks;// = new std::vector<float>;
  std::vector<float> *numu_cc_2_n_daughter_all;// = new std::vector<float>;

  //numu 3
  float numu_cc_flag_3;
  float numu_cc_3_particle_type;
  float numu_cc_3_max_length;
  float numu_cc_3_acc_track_length;
  float numu_cc_3_max_length_all;
  float numu_cc_3_max_muon_length;
  float numu_cc_3_n_daughter_tracks;
  float numu_cc_3_n_daughter_all;

  // cosmict 2_4
  float cosmict_flag_2;
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
  
  float cosmict_flag_4;  // kinematics muon
  float cosmict_4_filled;
  float cosmict_4_flag_inside;
  float cosmict_4_angle_beam;
  float cosmict_4_connected_showers;  // need to be careful about the nueCC ...

  //cosmict 3_5
  float cosmict_flag_3;
  float cosmict_3_filled;
  float cosmict_3_flag_inside;
  float cosmict_3_angle_beam;
  float cosmict_3_flag_dir_weak;
  float cosmict_3_dQ_dx_end;
  float cosmict_3_dQ_dx_front;
  float cosmict_3_theta;
  float cosmict_3_phi;
  float cosmict_3_valid_tracks;
  
  float cosmict_flag_5;  // kinematics muon
  float cosmict_5_filled;
  float cosmict_5_flag_inside;
  float cosmict_5_angle_beam;
  float cosmict_5_connected_showers;
  
  //cosmict 6
  float cosmict_flag_6; // special ...
  float cosmict_6_filled;
  float cosmict_6_flag_dir_weak;
  float cosmict_6_flag_inside;
  float cosmict_6_angle;
  
  
  // cosmict 7
  float cosmict_flag_7;  // muon+ michel
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
  
  // cosmict 8
  float cosmict_flag_8;  // muon + michel + special
  float cosmict_8_filled;
  float cosmict_8_flag_out;
  float cosmict_8_muon_length;
  float cosmict_8_acc_length;
  
  // cosmict 10
  std::vector<float> *cosmict_flag_10;// = new std::vector<float>;  // front upstream (dirt)
  std::vector<float> *cosmict_10_flag_inside;//= new std::vector<float> ;
  std::vector<float> *cosmict_10_vtx_z;//= new std::vector<float>;
  std::vector<float> *cosmict_10_flag_shower;//= new std::vector<float>;
  std::vector<float> *cosmict_10_flag_dir_weak;//= new std::vector<float>;
  std::vector<float> *cosmict_10_angle_beam;//= new std::vector<float>;
  std::vector<float> *cosmict_10_length;//= new std::vector<float>;


  // cosmict 9
  float cosmict_flag_9;
  
  float cosmic_flag;
  float cosmic_filled;

  // overall
  float cosmict_flag;
  float numu_cc_flag;


  //scores
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
  float numu_score;

  float weight;
};

void set_tree_address(TTree *tree0, numuInfo& tagger_info);
void put_tree_address(TTree *Tsig, numuInfo& tagger_info);
void merge_files();
void convert_file(int level);

float cal_cosmict_2_4_bdt(float default_val,numuInfo& tagger_info);
float cal_cosmict_3_5_bdt(float default_val,numuInfo& tagger_info);
float cal_cosmict_6_bdt(float default_val,numuInfo& tagger_info);
float cal_cosmict_7_bdt(float default_val,numuInfo& tagger_info);
float cal_cosmict_8_bdt(float default_val,numuInfo& tagger_info);
float cal_cosmict_10_bdt(float default_val,numuInfo& tagger_info);

float cal_numu_1_bdt(float default_val,numuInfo& tagger_info,TMVA::Reader& reader,
		     float& numu_cc_flag_1,
		     float& numu_cc_1_particle_type,
		     float& numu_cc_1_length,
		     float& numu_cc_1_medium_dQ_dx,
		     float& numu_cc_1_dQ_dx_cut,
		     float& numu_cc_1_direct_length,
		     float& numu_cc_1_n_daughter_tracks,
		     float& numu_cc_1_n_daughter_all);
float cal_numu_2_bdt(float default_val,numuInfo& tagger_info,TMVA::Reader& reader,
		     float& numu_cc_2_length,
		     float& numu_cc_2_total_length,
		     float& numu_cc_2_n_daughter_tracks,
		     float& numu_cc_2_n_daughter_all);
float cal_numu_3_bdt(float default_val,numuInfo& tagger_info,TMVA::Reader& reader);


void set_tree_address(TTree *tree0, numuInfo& tagger_info){
  tree0->SetBranchAddress("numu_cc_flag",&tagger_info.numu_cc_flag);
  tree0->SetBranchAddress("cosmict_flag",&tagger_info.cosmict_flag);
  tree0->SetBranchAddress("cosmic_flag",&tagger_info.cosmic_flag);
  tree0->SetBranchAddress("cosmic_filled",&tagger_info.cosmic_filled);
  tree0->SetBranchAddress("cosmict_flag_9",&tagger_info.cosmict_flag_9);
  tree0->SetBranchAddress("cosmict_flag_1",&tagger_info.cosmict_flag_1);

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

  tree0->SetBranchAddress("cosmict_flag_2", &tagger_info.cosmict_flag_2);
  tree0->SetBranchAddress("cosmict_2_filled", &tagger_info.cosmict_2_filled);
  tree0->SetBranchAddress("cosmict_2_particle_type", &tagger_info.cosmict_2_particle_type);
  tree0->SetBranchAddress("cosmict_2_n_muon_tracks", &tagger_info.cosmict_2_n_muon_tracks);
  tree0->SetBranchAddress("cosmict_2_total_shower_length", &tagger_info.cosmict_2_total_shower_length);
  tree0->SetBranchAddress("cosmict_2_flag_inside", &tagger_info.cosmict_2_flag_inside);
  tree0->SetBranchAddress("cosmict_2_angle_beam", &tagger_info.cosmict_2_angle_beam);
  tree0->SetBranchAddress("cosmict_2_flag_dir_weak", &tagger_info.cosmict_2_flag_dir_weak);
  tree0->SetBranchAddress("cosmict_2_dQ_dx_end", &tagger_info.cosmict_2_dQ_dx_end);
  tree0->SetBranchAddress("cosmict_2_dQ_dx_front", &tagger_info.cosmict_2_dQ_dx_front);
  tree0->SetBranchAddress("cosmict_2_theta", &tagger_info.cosmict_2_theta);
  tree0->SetBranchAddress("cosmict_2_phi", &tagger_info.cosmict_2_phi);
  tree0->SetBranchAddress("cosmict_2_valid_tracks", &tagger_info.cosmict_2_valid_tracks);
  
  tree0->SetBranchAddress("cosmict_flag_4",&tagger_info.cosmict_flag_4);
  tree0->SetBranchAddress("cosmict_4_filled",&tagger_info.cosmict_4_filled);
  tree0->SetBranchAddress("cosmict_4_flag_inside",&tagger_info.cosmict_4_flag_inside);
  tree0->SetBranchAddress("cosmict_4_angle_beam",&tagger_info.cosmict_4_angle_beam);
  tree0->SetBranchAddress("cosmict_4_connected_showers",&tagger_info.cosmict_4_connected_showers);

  tree0->SetBranchAddress("cosmict_flag_3", &tagger_info.cosmict_flag_3);
  tree0->SetBranchAddress("cosmict_3_filled",&tagger_info.cosmict_3_filled);
  tree0->SetBranchAddress("cosmict_3_flag_inside",&tagger_info.cosmict_3_flag_inside);
  tree0->SetBranchAddress("cosmict_3_angle_beam",&tagger_info.cosmict_3_angle_beam);
  tree0->SetBranchAddress("cosmict_3_flag_dir_weak",&tagger_info.cosmict_3_flag_dir_weak);
  tree0->SetBranchAddress("cosmict_3_dQ_dx_end",&tagger_info.cosmict_3_dQ_dx_end);
  tree0->SetBranchAddress("cosmict_3_dQ_dx_front",&tagger_info.cosmict_3_dQ_dx_front);
  tree0->SetBranchAddress("cosmict_3_theta",&tagger_info.cosmict_3_theta);
  tree0->SetBranchAddress("cosmict_3_phi",&tagger_info.cosmict_3_phi);
  tree0->SetBranchAddress("cosmict_3_valid_tracks",&tagger_info.cosmict_3_valid_tracks);
  
  tree0->SetBranchAddress("cosmict_flag_5",&tagger_info.cosmict_flag_5);
  tree0->SetBranchAddress("cosmict_5_filled",&tagger_info.cosmict_5_filled);
  tree0->SetBranchAddress("cosmict_5_flag_inside",&tagger_info.cosmict_5_flag_inside);
  tree0->SetBranchAddress("cosmict_5_angle_beam",&tagger_info.cosmict_5_angle_beam);
  tree0->SetBranchAddress("cosmict_5_connected_showers",&tagger_info.cosmict_5_connected_showers);

  tree0->SetBranchAddress("cosmict_flag_6",&tagger_info.cosmict_flag_6);
  tree0->SetBranchAddress("cosmict_6_filled",&tagger_info.cosmict_6_filled);
  tree0->SetBranchAddress("cosmict_6_flag_dir_weak",&tagger_info.cosmict_6_flag_dir_weak);
  tree0->SetBranchAddress("cosmict_6_flag_inside",&tagger_info.cosmict_6_flag_inside);
  tree0->SetBranchAddress("cosmict_6_angle",&tagger_info.cosmict_6_angle);

  tree0->SetBranchAddress("cosmict_flag_7",&tagger_info.cosmict_flag_7);
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

  tree0->SetBranchAddress("cosmict_flag_8",&tagger_info.cosmict_flag_8);
  tree0->SetBranchAddress("cosmict_8_filled",&tagger_info.cosmict_8_filled);
  tree0->SetBranchAddress("cosmict_8_flag_out",&tagger_info.cosmict_8_flag_out);
  tree0->SetBranchAddress("cosmict_8_muon_length",&tagger_info.cosmict_8_muon_length);
  tree0->SetBranchAddress("cosmict_8_acc_length",&tagger_info.cosmict_8_acc_length);

  tree0->SetBranchAddress("cosmict_flag_1",&tagger_info.cosmict_flag_1);
  tree0->SetBranchAddress("cosmict_flag_10",&tagger_info.cosmict_flag_10);
  tree0->SetBranchAddress("cosmict_10_flag_inside",&tagger_info.cosmict_10_flag_inside);
  tree0->SetBranchAddress("cosmict_10_vtx_z",&tagger_info.cosmict_10_vtx_z);
  tree0->SetBranchAddress("cosmict_10_flag_shower",&tagger_info.cosmict_10_flag_shower);
  tree0->SetBranchAddress("cosmict_10_flag_dir_weak",&tagger_info.cosmict_10_flag_dir_weak);
  tree0->SetBranchAddress("cosmict_10_angle_beam",&tagger_info.cosmict_10_angle_beam);
  tree0->SetBranchAddress("cosmict_10_length",&tagger_info.cosmict_10_length);
}

void put_tree_address(TTree *Tsig, numuInfo& tagger_info){
  Tsig->Branch("weight",&tagger_info.weight,"data/F");
  Tsig->Branch("numu_cc_flag",&tagger_info.numu_cc_flag,"data/F");
  Tsig->Branch("cosmict_flag",&tagger_info.cosmict_flag,"data/F");
  Tsig->Branch("cosmic_flag",&tagger_info.cosmic_flag,"data/F");
  Tsig->Branch("cosmic_filled",&tagger_info.cosmic_filled,"data/F");
  Tsig->Branch("cosmict_flag_9",&tagger_info.cosmict_flag_9,"data/F");
  Tsig->Branch("cosmict_flag_1",&tagger_info.cosmict_flag_1,"data/F");

  Tsig->Branch("numu_cc_flag_1",&tagger_info.numu_cc_flag_1);
  Tsig->Branch("numu_cc_1_particle_type",&tagger_info.numu_cc_1_particle_type);
  Tsig->Branch("numu_cc_1_length",&tagger_info.numu_cc_1_length);
  Tsig->Branch("numu_cc_1_medium_dQ_dx",&tagger_info.numu_cc_1_medium_dQ_dx);
  Tsig->Branch("numu_cc_1_dQ_dx_cut",&tagger_info.numu_cc_1_dQ_dx_cut);
  Tsig->Branch("numu_cc_1_direct_length",&tagger_info.numu_cc_1_direct_length);
  Tsig->Branch("numu_cc_1_n_daughter_tracks",&tagger_info.numu_cc_1_n_daughter_tracks);
  Tsig->Branch("numu_cc_1_n_daughter_all",&tagger_info.numu_cc_1_n_daughter_all);

  Tsig->Branch("numu_cc_flag_2",&tagger_info.numu_cc_flag_2);
  Tsig->Branch("numu_cc_2_length",&tagger_info.numu_cc_2_length);
  Tsig->Branch("numu_cc_2_total_length",&tagger_info.numu_cc_2_total_length);
  Tsig->Branch("numu_cc_2_n_daughter_tracks",&tagger_info.numu_cc_2_n_daughter_tracks);
  Tsig->Branch("numu_cc_2_n_daughter_all",&tagger_info.numu_cc_2_n_daughter_all);

  
  Tsig->Branch("numu_cc_flag_3",&tagger_info.numu_cc_flag_3,"numu_cc_flag_3/F");
  Tsig->Branch("numu_cc_3_particle_type",&tagger_info.numu_cc_3_particle_type,"numu_cc_3_particle_type/F");
  Tsig->Branch("numu_cc_3_max_length",&tagger_info.numu_cc_3_max_length,"numu_cc_3_max_length/F");
  Tsig->Branch("numu_cc_3_track_length",&tagger_info.numu_cc_3_acc_track_length,"numu_cc_3_acc_track_length/F");
  Tsig->Branch("numu_cc_3_max_length_all",&tagger_info.numu_cc_3_max_length_all,"numu_cc_3_max_length_all/F");
  Tsig->Branch("numu_cc_3_max_muon_length",&tagger_info.numu_cc_3_max_muon_length,"numu_cc_3_max_muon_length/F");
  Tsig->Branch("numu_cc_3_n_daughter_tracks",&tagger_info.numu_cc_3_n_daughter_tracks,"numu_cc_3_n_daughter_tracks/F");
  Tsig->Branch("numu_cc_3_n_daughter_all",&tagger_info.numu_cc_3_n_daughter_all,"numu_cc_3_n_daughter_all/F");

  Tsig->Branch("cosmict_flag_2",&tagger_info.cosmict_flag_2,"cosmict_flag_2/F");
  Tsig->Branch("cosmict_2_filled",&tagger_info.cosmict_2_filled,"cosmict_2_filled/F");
  Tsig->Branch("cosmict_2_particle_type",&tagger_info.cosmict_2_particle_type,"cosmict_2_particle_type/F");
  Tsig->Branch("cosmict_2_n_muon_tracks",&tagger_info.cosmict_2_n_muon_tracks,"cosmict_2_n_muon_tracks/F");
  Tsig->Branch("cosmict_2_total_shower_length",&tagger_info.cosmict_2_total_shower_length,"cosmict_2_total_shower_length/F");
  Tsig->Branch("cosmict_2_flag_inside",&tagger_info.cosmict_2_flag_inside,"cosmict_2_flag_inside/F");
  Tsig->Branch("cosmict_2_angle_beam",&tagger_info.cosmict_2_angle_beam,"cosmict_2_angle_beam/F");
  Tsig->Branch("cosmict_2_flag_dir_weak",&tagger_info.cosmict_2_flag_dir_weak,"cosmict_2_flag_dir_weak/F");
  Tsig->Branch("cosmict_2_dQ_dx_end",&tagger_info.cosmict_2_dQ_dx_end,"cosmict_2_dQ_dx_end/F");
  Tsig->Branch("cosmict_2_dQ_dx_front",&tagger_info.cosmict_2_dQ_dx_front,"cosmict_2_dQ_dx_front/F");
  Tsig->Branch("cosmict_2_theta",&tagger_info.cosmict_2_theta,"cosmict_2_theta/F");
  Tsig->Branch("cosmict_2_phi",&tagger_info.cosmict_2_phi,"cosmict_2_phi/F");
  Tsig->Branch("cosmict_2_valid_tracks",&tagger_info.cosmict_2_valid_tracks,"cosmict_2_valid_tracks/F");
  
  Tsig->Branch("cosmict_flag_4",&tagger_info.cosmict_flag_4,"cosmict_flag_4/F");
  Tsig->Branch("cosmict_4_filled",&tagger_info.cosmict_4_filled,"cosmict_4_filled/F");
  Tsig->Branch("cosmict_4_flag_inside",&tagger_info.cosmict_4_flag_inside,"cosmict_4_flag_inside/F");
  Tsig->Branch("cosmict_4_angle_beam",&tagger_info.cosmict_4_angle_beam,"cosmict_4_angle_beam/F");
  Tsig->Branch("cosmict_4_connected_showers",&tagger_info.cosmict_4_connected_showers,"cosmict_4_connected_showers/F");

  
  Tsig->Branch("cosmict_flag_3",&tagger_info.cosmict_flag_3,"cosmict_flag_3/F");
  Tsig->Branch("cosmict_3_filled",&tagger_info.cosmict_3_filled,"cosmict_3_filled/F");
  Tsig->Branch("cosmict_3_flag_inside",&tagger_info.cosmict_3_flag_inside,"cosmict_3_flag_inside/F");
  Tsig->Branch("cosmict_3_angle_beam",&tagger_info.cosmict_3_angle_beam,"cosmict_3_angle_beam/F");
  Tsig->Branch("cosmict_3_flag_dir_weak",&tagger_info.cosmict_3_flag_dir_weak,"cosmict_3_flag_dir_weak/F");
  Tsig->Branch("cosmict_3_dQ_dx_end",&tagger_info.cosmict_3_dQ_dx_end,"cosmict_3_dQ_dx_end/F");
  Tsig->Branch("cosmict_3_dQ_dx_front",&tagger_info.cosmict_3_dQ_dx_front,"cosmict_3_dQ_dx_front/F");
  Tsig->Branch("cosmict_3_theta",&tagger_info.cosmict_3_theta,"cosmict_3_theta/F");
  Tsig->Branch("cosmict_3_phi",&tagger_info.cosmict_3_phi,"cosmict_3_phi/F");
  Tsig->Branch("cosmict_3_valid_tracks",&tagger_info.cosmict_3_valid_tracks,"cosmict_3_valid_tracks/F");
  
  Tsig->Branch("cosmict_flag_5",&tagger_info.cosmict_flag_5,"cosmict_flag_5/F");
  Tsig->Branch("cosmict_5_filled",&tagger_info.cosmict_5_filled,"cosmict_5_filled/F");
  Tsig->Branch("cosmict_5_flag_inside",&tagger_info.cosmict_5_flag_inside,"cosmict_5_flag_inside/F");
  Tsig->Branch("cosmict_5_angle_beam",&tagger_info.cosmict_5_angle_beam,"cosmict_5_angle_beam/F");
  Tsig->Branch("cosmict_5_connected_showers",&tagger_info.cosmict_5_connected_showers,"cosmict_5_connected_showers/F");

  Tsig->Branch("cosmict_flag_6",&tagger_info.cosmict_flag_6,"cosmict_flag_6/F");
  Tsig->Branch("cosmict_6_filled",&tagger_info.cosmict_6_filled,"cosmict_6_filled/F");
  Tsig->Branch("cosmict_6_flag_dir_weak",&tagger_info.cosmict_6_flag_dir_weak,"cosmict_6_flag_dir_weak/F");
  Tsig->Branch("cosmict_6_flag_inside",&tagger_info.cosmict_6_flag_inside,"cosmict_6_flag_inside/F");
  Tsig->Branch("cosmict_6_angle",&tagger_info.cosmict_6_angle,"cosmict_6_angle/F");

  Tsig->Branch("cosmict_flag_7",&tagger_info.cosmict_flag_7,"cosmict_flag_7/F");
  Tsig->Branch("cosmict_7_filled",&tagger_info.cosmict_7_filled,"cosmict_7_filled/F");
  Tsig->Branch("cosmict_7_flag_sec",&tagger_info.cosmict_7_flag_sec,"cosmict_7_flag_sec/F");
  Tsig->Branch("cosmict_7_n_muon_tracks",&tagger_info.cosmict_7_n_muon_tracks,"cosmict_7_n_muon_tracks/F");
  Tsig->Branch("cosmict_7_total_shower_length",&tagger_info.cosmict_7_total_shower_length,"cosmict_7_total_shower_length/F");
  Tsig->Branch("cosmict_7_flag_inside",&tagger_info.cosmict_7_flag_inside,"cosmict_7_flag_inside/F");
  Tsig->Branch("cosmict_7_angle_beam",&tagger_info.cosmict_7_angle_beam,"cosmict_7_angle_beam/F");
  Tsig->Branch("cosmict_7_flag_dir_weak",&tagger_info.cosmict_7_flag_dir_weak,"cosmict_7_flag_dir_weak/F");
  Tsig->Branch("cosmict_7_dQ_dx_end",&tagger_info.cosmict_7_dQ_dx_end,"cosmict_7_dQ_dx_end/F");
  Tsig->Branch("cosmict_7_dQ_dx_front",&tagger_info.cosmict_7_dQ_dx_front,"cosmict_7_dQ_dx_front/F");
  Tsig->Branch("cosmict_7_theta",&tagger_info.cosmict_7_theta,"cosmict_7_theta/F");
  Tsig->Branch("cosmict_7_phi",&tagger_info.cosmict_7_phi,"cosmict_7_phi/F");

  Tsig->Branch("cosmict_flag_8",&tagger_info.cosmict_flag_8,"cosmict_flag_8/F");
  Tsig->Branch("cosmict_8_filled",&tagger_info.cosmict_8_filled,"cosmict_8_filled/F");
  Tsig->Branch("cosmict_8_flag_out",&tagger_info.cosmict_8_flag_out,"cosmict_8_flag_out/F");
  Tsig->Branch("cosmict_8_muon_length",&tagger_info.cosmict_8_muon_length,"cosmict_8_muon_length/F");
  Tsig->Branch("cosmict_8_acc_length",&tagger_info.cosmict_8_acc_length,"cosmict_8_acc_length/F");

  Tsig->Branch("cosmict_flag_10",&tagger_info.cosmict_flag_10);
  Tsig->Branch("cosmict_10_flag_inside",&tagger_info.cosmict_10_flag_inside);
  Tsig->Branch("cosmict_10_vtx_z",&tagger_info.cosmict_10_vtx_z);
  Tsig->Branch("cosmict_10_flag_shower",&tagger_info.cosmict_10_flag_shower);
  Tsig->Branch("cosmict_10_flag_dir_weak",&tagger_info.cosmict_10_flag_dir_weak);
  Tsig->Branch("cosmict_10_angle_beam",&tagger_info.cosmict_10_angle_beam);
  Tsig->Branch("cosmict_10_length",&tagger_info.cosmict_10_length);
}

float cal_cosmict_2_4_bdt(float default_val,numuInfo& tagger_info){
  float val = default_val;
  TMVA::Reader reader;

  return val;
}
float cal_cosmict_3_5_bdt(float default_val,numuInfo& tagger_info){
  float val = default_val;
  TMVA::Reader reader;

  return val;
}
float cal_cosmict_6_bdt(float default_val,numuInfo& tagger_info){
  float val = default_val;
  TMVA::Reader reader;

  return val;
}
float cal_cosmict_7_bdt(float default_val,numuInfo& tagger_info){
  float val = default_val;
  TMVA::Reader reader;

  return val;
}
float cal_cosmict_8_bdt(float default_val,numuInfo& tagger_info){
  float val = default_val;
  TMVA::Reader reader;

  return val;
}
float cal_cosmict_10_bdt(float default_val,numuInfo& tagger_info){
  float val = default_val;
  TMVA::Reader reader;

  return val;
}

float cal_numu_1_bdt(float default_val,numuInfo& tagger_info,TMVA::Reader& reader,
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

      float tmp_bdt =  reader.EvaluateMVA("MyBDT");
      if (tmp_bdt > val) val = tmp_bdt;
    }
  }
  
  return val;
}
float cal_numu_2_bdt(float default_val,numuInfo& tagger_info,TMVA::Reader& reader,
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
float cal_numu_3_bdt(float default_val,numuInfo& tagger_info,TMVA::Reader& reader){
  float val = default_val;

  val = reader.EvaluateMVA("MyBDT");

  return val;
}

void merge_files(){
  TFile *file0 = new TFile("T_tagger_overlayBNB_signal_1-7.root"); // sig x1
  TFile *file1 = new TFile("T_tagger_overlayBNB_background.root"); // bkg x1
  TFile *file2 = new TFile("T_tagger_extBNB_0-9.root"); // bkg x0.45
  
  TTree *tree0 = (TTree*)file0->Get("T_tagger");
  TTree *tree1 = (TTree*)file1->Get("T_tagger");
  TTree *tree2 = (TTree*)file2->Get("T_tagger");

  numuInfo tagger_info;

  tagger_info.numu_cc_flag_1 = new std::vector<float>;
  tagger_info.numu_cc_1_particle_type= new std::vector<float>;
  tagger_info.numu_cc_1_length= new std::vector<float>;
  tagger_info.numu_cc_1_medium_dQ_dx= new std::vector<float>;
  tagger_info.numu_cc_1_dQ_dx_cut= new std::vector<float>;
  tagger_info.numu_cc_1_direct_length= new std::vector<float>;
  tagger_info.numu_cc_1_n_daughter_tracks= new std::vector<float>;
  tagger_info.numu_cc_1_n_daughter_all= new std::vector<float>;

  tagger_info.numu_cc_flag_2 = new std::vector<float>;
  tagger_info.numu_cc_2_length = new std::vector<float>;
  tagger_info.numu_cc_2_total_length = new std::vector<float>;
  tagger_info.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger_info.numu_cc_2_n_daughter_all = new std::vector<float>;

  tagger_info.cosmict_flag_10= new std::vector<float>;  // front upstream (dirt)
  tagger_info.cosmict_10_flag_inside= new std::vector<float> ;
  tagger_info.cosmict_10_vtx_z= new std::vector<float>;
  tagger_info.cosmict_10_flag_shower= new std::vector<float>;
  tagger_info.cosmict_10_flag_dir_weak= new std::vector<float>;
  tagger_info.cosmict_10_angle_beam= new std::vector<float>;
  tagger_info.cosmict_10_length= new std::vector<float>;

  
  set_tree_address(tree0, tagger_info);
  set_tree_address(tree1, tagger_info);
  set_tree_address(tree2, tagger_info);

  
  
  TFile *new_file = new TFile("merge.root","RECREATE");
  TTree *Tsig = new TTree("sig","sig");
  TTree *Tbkg = new TTree("bkg","bkg");
  Tsig->SetDirectory(new_file);
  Tbkg->SetDirectory(new_file);

 
  put_tree_address(Tsig, tagger_info);
  put_tree_address(Tbkg, tagger_info);

 

  
  for (Int_t i=0;i!=tree0->GetEntries();i++){
    tree0->GetEntry(i);
    
    tagger_info.weight = 1;
    Tsig->Fill();
  }
  
  for (Int_t i=0;i!=tree1->GetEntries();i++){
    tree1->GetEntry(i);
    
    tagger_info.weight = 1;
    Tbkg->Fill();
  }
  
  for (Int_t i=0;i!=tree2->GetEntries();i++){
    tree2->GetEntry(i);
    
    tagger_info.weight = 0.45;
    Tbkg->Fill();
  }

  cout << "signal tree entries: " << Tsig->GetEntries() << endl;
  cout << "background tree entries: " << Tbkg->GetEntries() << endl;
  
  new_file->Write();
  new_file->Close();
}

void convert_file(int level){
   numuInfo tagger_info;

  tagger_info.numu_cc_flag_1 = new std::vector<float>;
  tagger_info.numu_cc_1_particle_type= new std::vector<float>;
  tagger_info.numu_cc_1_length= new std::vector<float>;
  tagger_info.numu_cc_1_medium_dQ_dx= new std::vector<float>;
  tagger_info.numu_cc_1_dQ_dx_cut= new std::vector<float>;
  tagger_info.numu_cc_1_direct_length= new std::vector<float>;
  tagger_info.numu_cc_1_n_daughter_tracks= new std::vector<float>;
  tagger_info.numu_cc_1_n_daughter_all= new std::vector<float>;

  tagger_info.numu_cc_flag_2 = new std::vector<float>;
  tagger_info.numu_cc_2_length = new std::vector<float>;
  tagger_info.numu_cc_2_total_length = new std::vector<float>;
  tagger_info.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger_info.numu_cc_2_n_daughter_all = new std::vector<float>;

  tagger_info.cosmict_flag_10= new std::vector<float>;  // front upstream (dirt)
  tagger_info.cosmict_10_flag_inside= new std::vector<float> ;
  tagger_info.cosmict_10_vtx_z= new std::vector<float>;
  tagger_info.cosmict_10_flag_shower= new std::vector<float>;
  tagger_info.cosmict_10_flag_dir_weak= new std::vector<float>;
  tagger_info.cosmict_10_angle_beam= new std::vector<float>;
  tagger_info.cosmict_10_length= new std::vector<float>;

  TFile *file = new TFile("merge.root");
  TTree *sig = (TTree*)file->Get("sig");
  TTree *bkg = (TTree*)file->Get("bkg");

  set_tree_address(sig, tagger_info);
  set_tree_address(bkg, tagger_info);

  sig->SetBranchAddress("weight",&tagger_info.weight);
  bkg->SetBranchAddress("weight",&tagger_info.weight);


  TFile *new_file = new TFile(Form("level_%d.root",level),"RECREATE");
  TTree *Tsig = new TTree("sig","sig");
  TTree *Tbkg = new TTree("bkg","bkg");
  Tsig->SetDirectory(new_file);
  Tbkg->SetDirectory(new_file);
  
  put_tree_address(Tsig, tagger_info);
  put_tree_address(Tbkg, tagger_info);
  
  Tsig->Branch("cosmict_2_4_score",&tagger_info.cosmict_2_4_score,"data/F");
  Tsig->Branch("cosmict_3_5_score",&tagger_info.cosmict_3_5_score,"data/F");
  Tsig->Branch("cosmict_6_score",&tagger_info.cosmict_6_score,"data/F");
  Tsig->Branch("cosmict_7_score",&tagger_info.cosmict_7_score,"data/F");
  Tsig->Branch("cosmict_8_score",&tagger_info.cosmict_8_score,"data/F");
  Tsig->Branch("cosmict_10_score",&tagger_info.cosmict_10_score,"data/F");
  Tsig->Branch("numu_1_score",&tagger_info.numu_1_score,"data/F");
  Tsig->Branch("numu_2_score",&tagger_info.numu_2_score,"data/F");
  Tsig->Branch("numu_3_score",&tagger_info.numu_3_score,"data/F");
  if (level!=1) Tsig->Branch("numu_score",&tagger_info.numu_score,"data/F");

  Tbkg->Branch("cosmict_2_4_score",&tagger_info.cosmict_2_4_score,"data/F");
  Tbkg->Branch("cosmict_3_5_score",&tagger_info.cosmict_3_5_score,"data/F");
  Tbkg->Branch("cosmict_6_score",&tagger_info.cosmict_6_score,"data/F");
  Tbkg->Branch("cosmict_7_score",&tagger_info.cosmict_7_score,"data/F");
  Tbkg->Branch("cosmict_8_score",&tagger_info.cosmict_8_score,"data/F");
  Tbkg->Branch("cosmict_10_score",&tagger_info.cosmict_10_score,"data/F");
  Tbkg->Branch("numu_1_score",&tagger_info.numu_1_score,"data/F");
  Tbkg->Branch("numu_2_score",&tagger_info.numu_2_score,"data/F");
  Tbkg->Branch("numu_3_score",&tagger_info.numu_3_score,"data/F");
  if (level!=1) Tbkg->Branch("numu_score",&tagger_info.numu_score,"data/F");


   //numu reader
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


  TMVA::Reader reader_numu_3;
  reader_numu_3.AddVariable( "numu_cc_3_particle_type", &tagger_info.numu_cc_3_particle_type);
  reader_numu_3.AddVariable( "numu_cc_3_max_length", &tagger_info.numu_cc_3_max_length);
  reader_numu_3.AddVariable( "numu_cc_3_acc_track_length", &tagger_info.numu_cc_3_acc_track_length);
  reader_numu_3.AddVariable( "numu_cc_3_max_length_all", &tagger_info.numu_cc_3_max_length_all);
  reader_numu_3.AddVariable( "numu_cc_3_max_muon_length", &tagger_info.numu_cc_3_max_muon_length);
  reader_numu_3.AddVariable( "numu_cc_3_n_daughter_tracks", &tagger_info.numu_cc_3_n_daughter_tracks);
  reader_numu_3.AddVariable( "numu_cc_3_n_daughter_all", &tagger_info.numu_cc_3_n_daughter_all);
    
  reader_numu_3.BookMVA( "MyBDT", "weights/numu_tagger3.weights.xml");
  
  
  for (int i=0;i!=sig->GetEntries();i++){
    sig->GetEntry(i);
    tagger_info.numu_1_score = cal_numu_1_bdt(-0.4,tagger_info, reader_numu_1, numu_cc_flag_1,
    					      numu_cc_1_particle_type,
    					      numu_cc_1_length,
    					      numu_cc_1_medium_dQ_dx,
    					      numu_cc_1_dQ_dx_cut,
    					      numu_cc_1_direct_length,
    					      numu_cc_1_n_daughter_tracks,
    					      numu_cc_1_n_daughter_all);
    tagger_info.numu_2_score = cal_numu_2_bdt(-0.1,tagger_info,reader_numu_2,
					      numu_cc_2_length,
					      numu_cc_2_total_length,
					      numu_cc_2_n_daughter_tracks,
					      numu_cc_2_n_daughter_all);
    tagger_info.numu_3_score = cal_numu_3_bdt(-0.2,tagger_info, reader_numu_3);
    Tsig->Fill();
  }

  for (int i=0;i!=bkg->GetEntries();i++){
    bkg->GetEntry(i);
    tagger_info.numu_1_score = cal_numu_1_bdt(-0.4,tagger_info, reader_numu_1, numu_cc_flag_1,
    					      numu_cc_1_particle_type,
    					      numu_cc_1_length,
    					      numu_cc_1_medium_dQ_dx,
    					      numu_cc_1_dQ_dx_cut,
    					      numu_cc_1_direct_length,
    					      numu_cc_1_n_daughter_tracks,
    					      numu_cc_1_n_daughter_all);
    tagger_info.numu_2_score = cal_numu_2_bdt(-0.1,tagger_info,reader_numu_2,
					      numu_cc_2_length,
					      numu_cc_2_total_length,
					      numu_cc_2_n_daughter_tracks,
					      numu_cc_2_n_daughter_all);
    tagger_info.numu_3_score = cal_numu_3_bdt(-0.2,tagger_info, reader_numu_3);   
    
    Tbkg->Fill();
  }
  
  
  new_file->Write();
  new_file->Close();
  
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

  if (process==1){
    std::cout << "Merge the initial rootfile" << std::endl;
    merge_files(); // three weight values ...
    
  }else if (process == 2){
    std::cout << "Prepare the initial rootfile" << std::endl;
    convert_file(1); // calculate first round BDT scores ...
    
  }else if (process == 3){
    std::cout << "Training BDT in ROOT TMVA " << std::endl;
    //    Run_r1();    
  }else if (process == 4){
    std::cout << "Validating BDT in ROOT TMVA " << std::endl;
  }
  
  return 1;
    
}
