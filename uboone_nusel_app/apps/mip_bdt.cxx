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


void Run_r1();
void InitBDT_r1();
void TestEvaluate(TString filename = "round_1.root");

void Run_r2();
void InitBDT_r2();

void convert_file();


void convert_file(){
  TFile *file = new TFile("bdtfile_0718.root");
  TTree *sig = (TTree*)file->Get("sig");
  TTree *bkg = (TTree*)file->Get("bkg");


  float trueEdep;
  float weight;
  float lowEweight;
  Int_t nueTag;

  sig->SetBranchAddress("trueEdep",&trueEdep);
  sig->SetBranchAddress("weight",&weight);
  sig->SetBranchAddress("lowEweight",&lowEweight);
  sig->SetBranchAddress("nueTag",&nueTag);
  
  bkg->SetBranchAddress("trueEdep",&trueEdep);
  bkg->SetBranchAddress("weight",&weight);
  bkg->SetBranchAddress("lowEweight",&lowEweight);
  bkg->SetBranchAddress("nueTag",&nueTag);

   // mip identification
  int mip_flag;
  double mip_energy;
  int mip_n_end_reduction;    
  int mip_n_first_mip;
  int mip_n_first_non_mip;
  int mip_n_first_non_mip_1;
  int mip_n_first_non_mip_2;
  double mip_vec_dQ_dx_0;
  double mip_vec_dQ_dx_1;
  double mip_max_dQ_dx_sample;
  int mip_n_below_threshold;
  int mip_n_below_zero;
  int mip_n_lowest;
  int mip_n_highest;
  double mip_lowest_dQ_dx;
  double mip_highest_dQ_dx;
  double mip_medium_dQ_dx;
  double mip_stem_length;
  double mip_length_main;
  double mip_length_total;
  double mip_angle_beam;
  double mip_iso_angle;
  int mip_n_vertex;
  int mip_n_good_tracks;
  double mip_E_indirect_max_energy;
  int mip_flag_all_above;
  double mip_min_dQ_dx_5;
  int mip_n_other_vertex; 
  int mip_n_stem_size;
  int mip_flag_stem_trajectory;
  double mip_min_dis;

  sig->SetBranchAddress("mip_flag",&mip_flag);
  sig->SetBranchAddress("mip_energy",&mip_energy);
  sig->SetBranchAddress("mip_n_end_reduction",&mip_n_end_reduction);
  sig->SetBranchAddress("mip_n_first_mip",&mip_n_first_mip);
  sig->SetBranchAddress("mip_n_first_non_mip",&mip_n_first_non_mip);
  sig->SetBranchAddress("mip_n_first_non_mip_1",&mip_n_first_non_mip_1);
  sig->SetBranchAddress("mip_n_first_non_mip_2",&mip_n_first_non_mip_2);
  sig->SetBranchAddress("mip_vec_dQ_dx_0",&mip_vec_dQ_dx_0);
  sig->SetBranchAddress("mip_vec_dQ_dx_1",&mip_vec_dQ_dx_1);
  sig->SetBranchAddress("mip_max_dQ_dx_sample",&mip_max_dQ_dx_sample);
  sig->SetBranchAddress("mip_n_below_threshold",&mip_n_below_threshold);
  sig->SetBranchAddress("mip_n_below_zero",&mip_n_below_zero);
  sig->SetBranchAddress("mip_n_lowest",&mip_n_lowest);
  sig->SetBranchAddress("mip_n_highest",&mip_n_highest);
  sig->SetBranchAddress("mip_lowest_dQ_dx",&mip_lowest_dQ_dx);
  sig->SetBranchAddress("mip_highest_dQ_dx",&mip_highest_dQ_dx);
  sig->SetBranchAddress("mip_medium_dQ_dx",&mip_medium_dQ_dx);
  sig->SetBranchAddress("mip_stem_length",&mip_stem_length);
  sig->SetBranchAddress("mip_length_main",&mip_length_main);
  sig->SetBranchAddress("mip_length_total",&mip_length_total);
  sig->SetBranchAddress("mip_angle_beam",&mip_angle_beam);
  sig->SetBranchAddress("mip_iso_angle",&mip_iso_angle);
  sig->SetBranchAddress("mip_n_vertex",&mip_n_vertex);
  sig->SetBranchAddress("mip_n_good_tracks",&mip_n_good_tracks);
  sig->SetBranchAddress("mip_E_indirect_max_energy",&mip_E_indirect_max_energy);
  sig->SetBranchAddress("mip_flag_all_above",&mip_flag_all_above);
  sig->SetBranchAddress("mip_min_dQ_dx_5",&mip_min_dQ_dx_5);
  sig->SetBranchAddress("mip_n_other_vertex",&mip_n_other_vertex);
  sig->SetBranchAddress("mip_n_stem_size",&mip_n_stem_size);
  sig->SetBranchAddress("mip_flag_stem_trajectory",&mip_flag_stem_trajectory);
  sig->SetBranchAddress("mip_min_dis",&mip_min_dis);


  bkg->SetBranchAddress("mip_flag",&mip_flag);
  bkg->SetBranchAddress("mip_energy",&mip_energy);
  bkg->SetBranchAddress("mip_n_end_reduction",&mip_n_end_reduction);
  bkg->SetBranchAddress("mip_n_first_mip",&mip_n_first_mip);
  bkg->SetBranchAddress("mip_n_first_non_mip",&mip_n_first_non_mip);
  bkg->SetBranchAddress("mip_n_first_non_mip_1",&mip_n_first_non_mip_1);
  bkg->SetBranchAddress("mip_n_first_non_mip_2",&mip_n_first_non_mip_2);
  bkg->SetBranchAddress("mip_vec_dQ_dx_0",&mip_vec_dQ_dx_0);
  bkg->SetBranchAddress("mip_vec_dQ_dx_1",&mip_vec_dQ_dx_1);
  bkg->SetBranchAddress("mip_max_dQ_dx_sample",&mip_max_dQ_dx_sample);
  bkg->SetBranchAddress("mip_n_below_threshold",&mip_n_below_threshold);
  bkg->SetBranchAddress("mip_n_below_zero",&mip_n_below_zero);
  bkg->SetBranchAddress("mip_n_lowest",&mip_n_lowest);
  bkg->SetBranchAddress("mip_n_highest",&mip_n_highest);
  bkg->SetBranchAddress("mip_lowest_dQ_dx",&mip_lowest_dQ_dx);
  bkg->SetBranchAddress("mip_highest_dQ_dx",&mip_highest_dQ_dx);
  bkg->SetBranchAddress("mip_medium_dQ_dx",&mip_medium_dQ_dx);
  bkg->SetBranchAddress("mip_stem_length",&mip_stem_length);
  bkg->SetBranchAddress("mip_length_main",&mip_length_main);
  bkg->SetBranchAddress("mip_length_total",&mip_length_total);
  bkg->SetBranchAddress("mip_angle_beam",&mip_angle_beam);
  bkg->SetBranchAddress("mip_iso_angle",&mip_iso_angle);
  bkg->SetBranchAddress("mip_n_vertex",&mip_n_vertex);
  bkg->SetBranchAddress("mip_n_good_tracks",&mip_n_good_tracks);
  bkg->SetBranchAddress("mip_E_indirect_max_energy",&mip_E_indirect_max_energy);
  bkg->SetBranchAddress("mip_flag_all_above",&mip_flag_all_above);
  bkg->SetBranchAddress("mip_min_dQ_dx_5",&mip_min_dQ_dx_5);
  bkg->SetBranchAddress("mip_n_other_vertex",&mip_n_other_vertex);
  bkg->SetBranchAddress("mip_n_stem_size",&mip_n_stem_size);
  bkg->SetBranchAddress("mip_flag_stem_trajectory",&mip_flag_stem_trajectory);
  bkg->SetBranchAddress("mip_min_dis",&mip_min_dis);
  
  
  
  TFile *new_file = new TFile("reduced.root","RECREATE");
  TTree *Tsig = new TTree("sig","sig");
  TTree *Tbkg = new TTree("bkg","bkg");
  Tsig->SetDirectory(new_file);
  Tbkg->SetDirectory(new_file);
  
  Tsig->Branch("trueEdep",&trueEdep,"data/F");
  Tsig->Branch("weight",&weight,"data/F");
  Tsig->Branch("lowEweight",&lowEweight,"data/F");
  Tsig->Branch("nueTag",&nueTag,"data/I");
    
  Tbkg->Branch("trueEdep",&trueEdep,"data/F");
  Tbkg->Branch("weight",&weight,"data/F");
  Tbkg->Branch("lowEweight",&lowEweight,"data/F");
  Tbkg->Branch("nueTag",&nueTag,"data/I");
  
  float mip_flag_f;
  float mip_energy_f;
  float mip_n_end_reduction_f;    
  float mip_n_first_mip_f;
  float mip_n_first_non_mip_f;
  float mip_n_first_non_mip_1_f;
  float mip_n_first_non_mip_2_f;
  float mip_vec_dQ_dx_0_f;
  float mip_vec_dQ_dx_1_f;
  float mip_max_dQ_dx_sample_f;
  float mip_n_below_threshold_f;
  float mip_n_below_zero_f;
  float mip_n_lowest_f;
  float mip_n_highest_f;
  float mip_lowest_dQ_dx_f;
  float mip_highest_dQ_dx_f;
  float mip_medium_dQ_dx_f;
  float mip_stem_length_f;
  float mip_length_main_f;
  float mip_length_total_f;
  float mip_angle_beam_f;
  float mip_iso_angle_f;
  float mip_n_vertex_f;
  float mip_n_good_tracks_f;
  float mip_E_indirect_max_energy_f;
  float mip_flag_all_above_f;
  float mip_min_dQ_dx_5_f;
  float mip_n_other_vertex_f; 
  float mip_n_stem_size_f;
  float mip_flag_stem_trajectory_f;
  float mip_min_dis_f;

  
  Tsig->Branch("mip_flag",&mip_flag_f,"mip_flag/F");
  Tsig->Branch("mip_energy",&mip_energy_f,"mip_energy/F");
  Tsig->Branch("mip_n_end_reduction",&mip_n_end_reduction_f,"mip_n_end_reduction/F");
  Tsig->Branch("mip_n_first_mip",&mip_n_first_mip_f,"mip_n_first_mip/F");
  Tsig->Branch("mip_n_first_non_mip",&mip_n_first_non_mip_f,"mip_n_first_non_mip/F");
  Tsig->Branch("mip_n_first_non_mip_1",&mip_n_first_non_mip_1_f,"mip_n_first_non_mip_1/F");
  Tsig->Branch("mip_n_first_non_mip_2",&mip_n_first_non_mip_2_f,"mip_n_first_non_mip_2/F");
  Tsig->Branch("mip_vec_dQ_dx_0",&mip_vec_dQ_dx_0_f,"mip_vec_dQ_dx_0/F");
  Tsig->Branch("mip_vec_dQ_dx_1",&mip_vec_dQ_dx_1_f,"mip_vec_dQ_dx_1/F");
  Tsig->Branch("mip_max_dQ_dx_sample",&mip_max_dQ_dx_sample_f,"mip_max_dQ_dx_sample/F");
  Tsig->Branch("mip_n_below_threshold",&mip_n_below_threshold_f,"mip_n_below_threshold/F");
  Tsig->Branch("mip_n_below_zero",&mip_n_below_zero_f,"mip_n_below_zero/F");
  Tsig->Branch("mip_n_lowest",&mip_n_lowest_f,"mip_n_lowest/F");
  Tsig->Branch("mip_n_highest",&mip_n_highest_f,"mip_n_highest/F");
  Tsig->Branch("mip_lowest_dQ_dx",&mip_lowest_dQ_dx_f,"mip_lowest_dQ_dx/F");
  Tsig->Branch("mip_highest_dQ_dx",&mip_highest_dQ_dx_f,"mip_highest_dQ_dx/F");
  Tsig->Branch("mip_medium_dQ_dx",&mip_medium_dQ_dx_f,"mip_medium_dQ_dx/F");
  Tsig->Branch("mip_stem_length",&mip_stem_length_f,"mip_stem_length/F");
  Tsig->Branch("mip_length_main",&mip_length_main_f,"mip_length_main/F");
  Tsig->Branch("mip_length_total",&mip_length_total_f,"mip_length_total/F");
  Tsig->Branch("mip_angle_beam",&mip_angle_beam_f,"mip_angle_beam/F");
  Tsig->Branch("mip_iso_angle",&mip_iso_angle_f,"mip_iso_angle/F");
  Tsig->Branch("mip_n_vertex",&mip_n_vertex_f,"mip_n_vertex/F");
  Tsig->Branch("mip_n_good_tracks",&mip_n_good_tracks_f,"mip_n_good_tracks/F");
  Tsig->Branch("mip_E_indirect_max_energy",&mip_E_indirect_max_energy_f,"mip_E_indirect_max_energy/F");
  Tsig->Branch("mip_flag_all_above",&mip_flag_all_above_f,"mip_flag_all_above/F");
  Tsig->Branch("mip_min_dQ_dx_5",&mip_min_dQ_dx_5_f,"mip_min_dQ_dx_5/F");
  Tsig->Branch("mip_n_other_vertex",&mip_n_other_vertex_f,"mip_n_other_vertex/F");
  Tsig->Branch("mip_n_stem_size",&mip_n_stem_size_f,"mip_n_stem_size/F");
  Tsig->Branch("mip_flag_stem_trajectory",&mip_flag_stem_trajectory_f,"mip_flag_stem_trajectory/F");
  Tsig->Branch("mip_min_dis",&mip_min_dis_f,"mip_min_dis/F");

  
  Tbkg->Branch("mip_flag",&mip_flag_f,"mip_flag/F");
  Tbkg->Branch("mip_energy",&mip_energy_f,"mip_energy/F");
  Tbkg->Branch("mip_n_end_reduction",&mip_n_end_reduction_f,"mip_n_end_reduction/F");
  Tbkg->Branch("mip_n_first_mip",&mip_n_first_mip_f,"mip_n_first_mip/F");
  Tbkg->Branch("mip_n_first_non_mip",&mip_n_first_non_mip_f,"mip_n_first_non_mip/F");
  Tbkg->Branch("mip_n_first_non_mip_1",&mip_n_first_non_mip_1_f,"mip_n_first_non_mip_1/F");
  Tbkg->Branch("mip_n_first_non_mip_2",&mip_n_first_non_mip_2_f,"mip_n_first_non_mip_2/F");
  Tbkg->Branch("mip_vec_dQ_dx_0",&mip_vec_dQ_dx_0_f,"mip_vec_dQ_dx_0/F");
  Tbkg->Branch("mip_vec_dQ_dx_1",&mip_vec_dQ_dx_1_f,"mip_vec_dQ_dx_1/F");
  Tbkg->Branch("mip_max_dQ_dx_sample",&mip_max_dQ_dx_sample_f,"mip_max_dQ_dx_sample/F");
  Tbkg->Branch("mip_n_below_threshold",&mip_n_below_threshold_f,"mip_n_below_threshold/F");
  Tbkg->Branch("mip_n_below_zero",&mip_n_below_zero_f,"mip_n_below_zero/F");
  Tbkg->Branch("mip_n_lowest",&mip_n_lowest_f,"mip_n_lowest/F");
  Tbkg->Branch("mip_n_highest",&mip_n_highest_f,"mip_n_highest/F");
  Tbkg->Branch("mip_lowest_dQ_dx",&mip_lowest_dQ_dx_f,"mip_lowest_dQ_dx/F");
  Tbkg->Branch("mip_highest_dQ_dx",&mip_highest_dQ_dx_f,"mip_highest_dQ_dx/F");
  Tbkg->Branch("mip_medium_dQ_dx",&mip_medium_dQ_dx_f,"mip_medium_dQ_dx/F");
  Tbkg->Branch("mip_stem_length",&mip_stem_length_f,"mip_stem_length/F");
  Tbkg->Branch("mip_length_main",&mip_length_main_f,"mip_length_main/F");
  Tbkg->Branch("mip_length_total",&mip_length_total_f,"mip_length_total/F");
  Tbkg->Branch("mip_angle_beam",&mip_angle_beam_f,"mip_angle_beam/F");
  Tbkg->Branch("mip_iso_angle",&mip_iso_angle_f,"mip_iso_angle/F");
  Tbkg->Branch("mip_n_vertex",&mip_n_vertex_f,"mip_n_vertex/F");
  Tbkg->Branch("mip_n_good_tracks",&mip_n_good_tracks_f,"mip_n_good_tracks/F");
  Tbkg->Branch("mip_E_indirect_max_energy",&mip_E_indirect_max_energy_f,"mip_E_indirect_max_energy/F");
  Tbkg->Branch("mip_flag_all_above",&mip_flag_all_above_f,"mip_flag_all_above/F");
  Tbkg->Branch("mip_min_dQ_dx_5",&mip_min_dQ_dx_5_f,"mip_min_dQ_dx_5/F");
  Tbkg->Branch("mip_n_other_vertex",&mip_n_other_vertex_f,"mip_n_other_vertex/F");
  Tbkg->Branch("mip_n_stem_size",&mip_n_stem_size_f,"mip_n_stem_size/F");
  Tbkg->Branch("mip_flag_stem_trajectory",&mip_flag_stem_trajectory_f,"mip_flag_stem_trajectory/F");
  Tbkg->Branch("mip_min_dis",&mip_min_dis_f,"mip_min_dis/F");

  
  for (Int_t i=0;i!=sig->GetEntries();i++){
    sig->GetEntry(i);
    
    mip_flag_f = mip_flag;
    mip_energy_f = mip_energy;
    mip_n_end_reduction_f = mip_n_end_reduction;     
    mip_n_first_mip_f = mip_n_first_mip;
    mip_n_first_non_mip_f = mip_n_first_non_mip;
    mip_n_first_non_mip_1_f = mip_n_first_non_mip_1;
    mip_n_first_non_mip_2_f = mip_n_first_non_mip_2;
    mip_vec_dQ_dx_0_f = mip_vec_dQ_dx_0;
    mip_vec_dQ_dx_1_f = mip_vec_dQ_dx_1;
    mip_max_dQ_dx_sample_f = mip_max_dQ_dx_sample;
    mip_n_below_threshold_f = mip_n_below_threshold;
    mip_n_below_zero_f = mip_n_below_zero;
    mip_n_lowest_f = mip_n_lowest;
    mip_n_highest_f = mip_n_highest;
    mip_lowest_dQ_dx_f = mip_lowest_dQ_dx;
    mip_highest_dQ_dx_f = mip_highest_dQ_dx;
    mip_medium_dQ_dx_f = mip_medium_dQ_dx;
    mip_stem_length_f = mip_stem_length;
    mip_length_main_f = mip_length_main;
    mip_length_total_f = mip_length_total;
    mip_angle_beam_f = mip_angle_beam;
    mip_iso_angle_f = mip_iso_angle;
    mip_n_vertex_f = mip_n_vertex;
    mip_n_good_tracks_f = mip_n_good_tracks;
    mip_E_indirect_max_energy_f = mip_E_indirect_max_energy;
    mip_flag_all_above_f = mip_flag_all_above;
    mip_min_dQ_dx_5_f = mip_min_dQ_dx_5;
    mip_n_other_vertex_f = mip_n_other_vertex;  
    mip_n_stem_size_f = mip_n_stem_size;
    mip_flag_stem_trajectory_f = mip_flag_stem_trajectory;
    mip_min_dis_f = mip_min_dis;
    
    Tsig->Fill();
  }

   for (Int_t i=0;i!=bkg->GetEntries();i++){
    bkg->GetEntry(i);
    
    mip_flag_f = mip_flag;
    mip_energy_f = mip_energy;
    mip_n_end_reduction_f = mip_n_end_reduction;     
    mip_n_first_mip_f = mip_n_first_mip;
    mip_n_first_non_mip_f = mip_n_first_non_mip;
    mip_n_first_non_mip_1_f = mip_n_first_non_mip_1;
    mip_n_first_non_mip_2_f = mip_n_first_non_mip_2;
    mip_vec_dQ_dx_0_f = mip_vec_dQ_dx_0;
    mip_vec_dQ_dx_1_f = mip_vec_dQ_dx_1;
    mip_max_dQ_dx_sample_f = mip_max_dQ_dx_sample;
    mip_n_below_threshold_f = mip_n_below_threshold;
    mip_n_below_zero_f = mip_n_below_zero;
    mip_n_lowest_f = mip_n_lowest;
    mip_n_highest_f = mip_n_highest;
    mip_lowest_dQ_dx_f = mip_lowest_dQ_dx;
    mip_highest_dQ_dx_f = mip_highest_dQ_dx;
    mip_medium_dQ_dx_f = mip_medium_dQ_dx;
    mip_stem_length_f = mip_stem_length;
    mip_length_main_f = mip_length_main;
    mip_length_total_f = mip_length_total;
    mip_angle_beam_f = mip_angle_beam;
    mip_iso_angle_f = mip_iso_angle;
    mip_n_vertex_f = mip_n_vertex;
    mip_n_good_tracks_f = mip_n_good_tracks;
    mip_E_indirect_max_energy_f = mip_E_indirect_max_energy;
    mip_flag_all_above_f = mip_flag_all_above;
    mip_min_dQ_dx_5_f = mip_min_dQ_dx_5;
    mip_n_other_vertex_f = mip_n_other_vertex;  
    mip_n_stem_size_f = mip_n_stem_size;
    mip_flag_stem_trajectory_f = mip_flag_stem_trajectory;
    mip_min_dis_f = mip_min_dis;
    
    
    Tbkg->Fill();
  }
  
  
  new_file->Write();
  new_file->Close();
  
  
}

void Run_r2(){
  TMVA::Tools::Instance();

  TString fname = "./round_1.root";
  input = TFile::Open( fname ); // check if file in local directory exists
  cout << "Using input file: " << input->GetName() << std::endl;

  TString outfileName( "TMVA_r2.root" );
  output = TFile::Open( outfileName, "RECREATE" );

  InitBDT_r2();
  
  // Train MVAs using the set of training events
  factory->TrainAllMethods();
  
  // Evaluate all MVAs using the set of test events
  factory->TestAllMethods();
  
  // Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();
  
  output->Close();
  cout << "Results saved at: " << output->GetName() << std::endl;
  
  delete factory;
  delete dataloader;
  
  TestEvaluate("round_2.root");
  
  // Launch the GUI for the root macros
  if (!gROOT->IsBatch()) TMVA::TMVAGui( output->GetName() );
  
}

void Run_r1()
{
    TMVA::Tools::Instance();

    TString fname = "./reduced.root";
    input = TFile::Open( fname ); // check if file in local directory exists
    cout << "Using input file: " << input->GetName() << std::endl;

    TString outfileName( "TMVA.root" );
    output = TFile::Open( outfileName, "RECREATE" );


    InitBDT_r1();


    // Train MVAs using the set of training events
    factory->TrainAllMethods();

    // Evaluate all MVAs using the set of test events
    factory->TestAllMethods();

    // Evaluate and compare performance of all configured MVAs
    factory->EvaluateAllMethods();
    
    output->Close();
    cout << "Results saved at: " << output->GetName() << std::endl;

    delete factory;
    delete dataloader;

    TestEvaluate("round_1.root");

    // Launch the GUI for the root macros
    if (!gROOT->IsBatch()) TMVA::TMVAGui( output->GetName() );

}




void InitBDT_r1()
{
    factory = new TMVA::Factory( "Test", output,
        "!V:!Silent:Color:DrawProgressBar:"
				 //    "Transformations=I;D;P;G,D:"
				 "AnalysisType=Classification" );
    			 


    dataloader = new TMVA::DataLoader("dataset");
    // Define the input variables that shall be used for the MVA training
    // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
    // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
    
    
    dataloader->AddVariable("mip_energy","mip_energy","MeV",'F');
    dataloader->AddVariable("mip_n_end_reduction","mip_n_end_reduction","",'F');
    dataloader->AddVariable("mip_n_first_mip","mip_n_first_mip","",'F');
    dataloader->AddVariable("mip_n_first_non_mip","mip_n_first_non_mip","",'F');
    dataloader->AddVariable("mip_n_first_non_mip_1","mip_n_first_non_mip_1","",'F');
    dataloader->AddVariable("mip_n_first_non_mip_2","mip_n_first_non_mip_2","",'F');
    dataloader->AddVariable("mip_vec_dQ_dx_0","mip_vec_dQ_dx_0","MeV/cm",'F');
    dataloader->AddVariable("mip_vec_dQ_dx_1","mip_vec_dQ_dx_1","MeV/cm",'F');
    dataloader->AddVariable("mip_max_dQ_dx_sample","mip_max_dQ_dx_sample","",'F');
    dataloader->AddVariable("mip_n_below_threshold","mip_n_below_threshold","",'F');
    dataloader->AddVariable("mip_n_below_zero","mip_n_below_zero","",'F');
    dataloader->AddVariable("mip_n_lowest","mip_n_lowest","",'F');
    dataloader->AddVariable("mip_n_highest","mip_n_highest","",'F');
    dataloader->AddVariable("mip_lowest_dQ_dx","mip_lowest_dQ_dx","MeV/cm",'F');
    dataloader->AddVariable("mip_highest_dQ_dx","mip_highest_dQ_dx","MeV/cm",'F');
    dataloader->AddVariable("mip_medium_dQ_dx","mip_medium_dQ_dx","MeV/cm",'F');
    dataloader->AddVariable("mip_stem_length","mip_stem_length","cm",'F');
    dataloader->AddVariable("mip_length_main","mip_length_main","cm",'F');
    dataloader->AddVariable("mip_length_total","mip_length_total","cm",'F');
    dataloader->AddVariable("mip_angle_beam","mip_angle_beam","deg",'F');
    dataloader->AddVariable("mip_iso_angle","mip_iso_angle","deg",'F');
    dataloader->AddVariable("mip_n_vertex","mip_n_vertex","",'F');
    dataloader->AddVariable("mip_n_good_tracks","mip_n_good_tracks","",'F');
    dataloader->AddVariable("mip_E_indirect_max_energy","mip_E_indirect_max_energy","MeV",'F');
    dataloader->AddVariable("mip_flag_all_above","mip_flag_all_above","",'F');
    dataloader->AddVariable("mip_min_dQ_dx_5","mip_min_dQ_dx_5","MeV/cm",'F');
    dataloader->AddVariable("mip_n_other_vertex","mip_n_other_vertex","",'F');
    dataloader->AddVariable("mip_n_stem_size","mip_n_stem_size","",'F');
    dataloader->AddVariable("mip_flag_stem_trajectory","mip_flag_stem_trajectory","",'F');
    dataloader->AddVariable("mip_min_dis","mip_min_dis","",'F');
    
    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "1>0"; //  4911/44194
    TCut mycut_b = "mip_flag == 0"; //  8033/21070
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=30000:"
					    "nTrain_Background=7000:"
					    "nTest_Signal=10000:"
					    "nTest_Background=1033:"
					    "SplitMode=Random:"
					    "NormMode=NumEvents:"
					    "!V" );
  
  // variations of BDTs
  // BDT: uses Adaptive Boost
  // BDTG: uses Gradient Boost
  // BDTB: uses Bagging
  // BDTD: decorrelation + Adaptive Boost
  // BDTF: allow usage of fisher discriminant for node splitting
  factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT",
		      "!H:!V:"
		      "VarTransform=D,P,U,G,D,N:"
		      "NTrees=850:"
		      "MinNodeSize=2.5%:"
		      "MaxDepth=3:"
		      "BoostType=AdaBoost:"
		      "AdaBoostBeta=0.5:"
		      "UseBaggedBoost:"
		      "BaggedSampleFraction=0.5:"
		      "SeparationType=GiniIndex:"
		      "nCuts=20");
  
}

void InitBDT_r2()
{
  factory = new TMVA::Factory( "Test", output,
        "!V:!Silent:Color:DrawProgressBar:"
				 //    "Transformations=I;D;P;G,D:"
				 "AnalysisType=Classification" );
    			 


    dataloader = new TMVA::DataLoader("dataset");
    // Define the input variables that shall be used for the MVA training
    // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
    // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
    
    
    dataloader->AddVariable("mip_energy","mip_energy","MeV",'F');
    dataloader->AddVariable("mip_n_end_reduction","mip_n_end_reduction","",'F');
    dataloader->AddVariable("mip_n_first_mip","mip_n_first_mip","",'F');
    dataloader->AddVariable("mip_n_first_non_mip","mip_n_first_non_mip","",'F');
    dataloader->AddVariable("mip_n_first_non_mip_1","mip_n_first_non_mip_1","",'F');
    dataloader->AddVariable("mip_n_first_non_mip_2","mip_n_first_non_mip_2","",'F');
    dataloader->AddVariable("mip_vec_dQ_dx_0","mip_vec_dQ_dx_0","MeV/cm",'F');
    dataloader->AddVariable("mip_vec_dQ_dx_1","mip_vec_dQ_dx_1","MeV/cm",'F');
    dataloader->AddVariable("mip_max_dQ_dx_sample","mip_max_dQ_dx_sample","",'F');
    dataloader->AddVariable("mip_n_below_threshold","mip_n_below_threshold","",'F');
    dataloader->AddVariable("mip_n_below_zero","mip_n_below_zero","",'F');
    dataloader->AddVariable("mip_n_lowest","mip_n_lowest","",'F');
    dataloader->AddVariable("mip_n_highest","mip_n_highest","",'F');
    dataloader->AddVariable("mip_lowest_dQ_dx","mip_lowest_dQ_dx","MeV/cm",'F');
    dataloader->AddVariable("mip_highest_dQ_dx","mip_highest_dQ_dx","MeV/cm",'F');
    dataloader->AddVariable("mip_medium_dQ_dx","mip_medium_dQ_dx","MeV/cm",'F');
    dataloader->AddVariable("mip_stem_length","mip_stem_length","cm",'F');
    dataloader->AddVariable("mip_length_main","mip_length_main","cm",'F');
    dataloader->AddVariable("mip_length_total","mip_length_total","cm",'F');
    dataloader->AddVariable("mip_angle_beam","mip_angle_beam","deg",'F');
    dataloader->AddVariable("mip_iso_angle","mip_iso_angle","deg",'F');
    dataloader->AddVariable("mip_n_vertex","mip_n_vertex","",'F');
    dataloader->AddVariable("mip_n_good_tracks","mip_n_good_tracks","",'F');
    dataloader->AddVariable("mip_E_indirect_max_energy","mip_E_indirect_max_energy","MeV",'F');
    dataloader->AddVariable("mip_flag_all_above","mip_flag_all_above","",'F');
    dataloader->AddVariable("mip_min_dQ_dx_5","mip_min_dQ_dx_5","MeV/cm",'F');
    dataloader->AddVariable("mip_n_other_vertex","mip_n_other_vertex","",'F');
    dataloader->AddVariable("mip_n_stem_size","mip_n_stem_size","",'F');
    dataloader->AddVariable("mip_flag_stem_trajectory","mip_flag_stem_trajectory","",'F');
    dataloader->AddVariable("mip_min_dis","mip_min_dis","",'F');
    
    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "1>0"; // 
    TCut mycut_b = "mip_flag == 0 || mip_bdt < -0.05"; // 10467
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=30000:"
					    "nTrain_Background=9000:"
					    "nTest_Signal=10000:"
					    "nTest_Background=1467:"
					    "SplitMode=Random:"
					    "NormMode=NumEvents:"
					    "!V" );
  
  // variations of BDTs
  // BDT: uses Adaptive Boost
  // BDTG: uses Gradient Boost
  // BDTB: uses Bagging
  // BDTD: decorrelation + Adaptive Boost
  // BDTF: allow usage of fisher discriminant for node splitting
  factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT",
		      "!H:!V:"
		      "VarTransform=D,P,U,G,D,N:"
		      "NTrees=850:"
		      "MinNodeSize=2.5%:"
		      "MaxDepth=3:"
		      "BoostType=AdaBoost:"
		      "AdaBoostBeta=0.5:"
		      "UseBaggedBoost:"
		      "BaggedSampleFraction=0.5:"
		      "SeparationType=GiniIndex:"
		      "nCuts=20");
  
}


void TestEvaluate(TString filename)
{
  

  Float_t bdt_value = 0;
    
  TFile *file = new TFile("reduced.root");
  TTree *sig = (TTree*)file->Get("sig");
  TTree *bkg = (TTree*)file->Get("bkg");

  float trueEdep;
  float weight;
  float lowEweight;
  Int_t nueTag;
  
  sig->SetBranchAddress("trueEdep",&trueEdep);
  sig->SetBranchAddress("weight",&weight);
  sig->SetBranchAddress("lowEweight",&lowEweight);
  sig->SetBranchAddress("nueTag",&nueTag);
  
  bkg->SetBranchAddress("trueEdep",&trueEdep);
  bkg->SetBranchAddress("weight",&weight);
  bkg->SetBranchAddress("lowEweight",&lowEweight);
  bkg->SetBranchAddress("nueTag",&nueTag);

  float mip_flag;
  float mip_energy;
  float mip_n_end_reduction;    
  float mip_n_first_mip;
  float mip_n_first_non_mip;
  float mip_n_first_non_mip_1;
  float mip_n_first_non_mip_2;
  float mip_vec_dQ_dx_0;
  float mip_vec_dQ_dx_1;
  float mip_max_dQ_dx_sample;
  float mip_n_below_threshold;
  float mip_n_below_zero;
  float mip_n_lowest;
  float mip_n_highest;
  float mip_lowest_dQ_dx;
  float mip_highest_dQ_dx;
  float mip_medium_dQ_dx;
  float mip_stem_length;
  float mip_length_main;
  float mip_length_total;
  float mip_angle_beam;
  float mip_iso_angle;
  float mip_n_vertex;
  float mip_n_good_tracks;
  float mip_E_indirect_max_energy;
  float mip_flag_all_above;
  float mip_min_dQ_dx_5;
  float mip_n_other_vertex; 
  float mip_n_stem_size;
  float mip_flag_stem_trajectory;
  float mip_min_dis;
  
  sig->SetBranchAddress("mip_flag",&mip_flag);
  sig->SetBranchAddress("mip_energy",&mip_energy);
  sig->SetBranchAddress("mip_n_end_reduction",&mip_n_end_reduction);
  sig->SetBranchAddress("mip_n_first_mip",&mip_n_first_mip);
  sig->SetBranchAddress("mip_n_first_non_mip",&mip_n_first_non_mip);
  sig->SetBranchAddress("mip_n_first_non_mip_1",&mip_n_first_non_mip_1);
  sig->SetBranchAddress("mip_n_first_non_mip_2",&mip_n_first_non_mip_2);
  sig->SetBranchAddress("mip_vec_dQ_dx_0",&mip_vec_dQ_dx_0);
  sig->SetBranchAddress("mip_vec_dQ_dx_1",&mip_vec_dQ_dx_1);
  sig->SetBranchAddress("mip_max_dQ_dx_sample",&mip_max_dQ_dx_sample);
  sig->SetBranchAddress("mip_n_below_threshold",&mip_n_below_threshold);
  sig->SetBranchAddress("mip_n_below_zero",&mip_n_below_zero);
  sig->SetBranchAddress("mip_n_lowest",&mip_n_lowest);
  sig->SetBranchAddress("mip_n_highest",&mip_n_highest);
  sig->SetBranchAddress("mip_lowest_dQ_dx",&mip_lowest_dQ_dx);
  sig->SetBranchAddress("mip_highest_dQ_dx",&mip_highest_dQ_dx);
  sig->SetBranchAddress("mip_medium_dQ_dx",&mip_medium_dQ_dx);
  sig->SetBranchAddress("mip_stem_length",&mip_stem_length);
  sig->SetBranchAddress("mip_length_main",&mip_length_main);
  sig->SetBranchAddress("mip_length_total",&mip_length_total);
  sig->SetBranchAddress("mip_angle_beam",&mip_angle_beam);
  sig->SetBranchAddress("mip_iso_angle",&mip_iso_angle);
  sig->SetBranchAddress("mip_n_vertex",&mip_n_vertex);
  sig->SetBranchAddress("mip_n_good_tracks",&mip_n_good_tracks);
  sig->SetBranchAddress("mip_E_indirect_max_energy",&mip_E_indirect_max_energy);
  sig->SetBranchAddress("mip_flag_all_above",&mip_flag_all_above);
  sig->SetBranchAddress("mip_min_dQ_dx_5",&mip_min_dQ_dx_5);
  sig->SetBranchAddress("mip_n_other_vertex",&mip_n_other_vertex);
  sig->SetBranchAddress("mip_n_stem_size",&mip_n_stem_size);
  sig->SetBranchAddress("mip_flag_stem_trajectory",&mip_flag_stem_trajectory);
  sig->SetBranchAddress("mip_min_dis",&mip_min_dis);


  bkg->SetBranchAddress("mip_flag",&mip_flag);
  bkg->SetBranchAddress("mip_energy",&mip_energy);
  bkg->SetBranchAddress("mip_n_end_reduction",&mip_n_end_reduction);
  bkg->SetBranchAddress("mip_n_first_mip",&mip_n_first_mip);
  bkg->SetBranchAddress("mip_n_first_non_mip",&mip_n_first_non_mip);
  bkg->SetBranchAddress("mip_n_first_non_mip_1",&mip_n_first_non_mip_1);
  bkg->SetBranchAddress("mip_n_first_non_mip_2",&mip_n_first_non_mip_2);
  bkg->SetBranchAddress("mip_vec_dQ_dx_0",&mip_vec_dQ_dx_0);
  bkg->SetBranchAddress("mip_vec_dQ_dx_1",&mip_vec_dQ_dx_1);
  bkg->SetBranchAddress("mip_max_dQ_dx_sample",&mip_max_dQ_dx_sample);
  bkg->SetBranchAddress("mip_n_below_threshold",&mip_n_below_threshold);
  bkg->SetBranchAddress("mip_n_below_zero",&mip_n_below_zero);
  bkg->SetBranchAddress("mip_n_lowest",&mip_n_lowest);
  bkg->SetBranchAddress("mip_n_highest",&mip_n_highest);
  bkg->SetBranchAddress("mip_lowest_dQ_dx",&mip_lowest_dQ_dx);
  bkg->SetBranchAddress("mip_highest_dQ_dx",&mip_highest_dQ_dx);
  bkg->SetBranchAddress("mip_medium_dQ_dx",&mip_medium_dQ_dx);
  bkg->SetBranchAddress("mip_stem_length",&mip_stem_length);
  bkg->SetBranchAddress("mip_length_main",&mip_length_main);
  bkg->SetBranchAddress("mip_length_total",&mip_length_total);
  bkg->SetBranchAddress("mip_angle_beam",&mip_angle_beam);
  bkg->SetBranchAddress("mip_iso_angle",&mip_iso_angle);
  bkg->SetBranchAddress("mip_n_vertex",&mip_n_vertex);
  bkg->SetBranchAddress("mip_n_good_tracks",&mip_n_good_tracks);
  bkg->SetBranchAddress("mip_E_indirect_max_energy",&mip_E_indirect_max_energy);
  bkg->SetBranchAddress("mip_flag_all_above",&mip_flag_all_above);
  bkg->SetBranchAddress("mip_min_dQ_dx_5",&mip_min_dQ_dx_5);
  bkg->SetBranchAddress("mip_n_other_vertex",&mip_n_other_vertex);
  bkg->SetBranchAddress("mip_n_stem_size",&mip_n_stem_size);
  bkg->SetBranchAddress("mip_flag_stem_trajectory",&mip_flag_stem_trajectory);
  bkg->SetBranchAddress("mip_min_dis",&mip_min_dis);
 
    
  
  TFile *new_file = new TFile(filename,"RECREATE");
  TTree *Tsig = new TTree("sig","sig");
  TTree *Tbkg = new TTree("bkg","bkg");
  Tsig->SetDirectory(new_file);
  Tbkg->SetDirectory(new_file);

  Tsig->Branch("trueEdep",&trueEdep,"data/F");
  Tsig->Branch("weight",&weight,"data/F");
  Tsig->Branch("lowEweight",&lowEweight,"data/F");
  Tsig->Branch("nueTag",&nueTag,"data/I");

  Tbkg->Branch("trueEdep",&trueEdep,"data/F");
  Tbkg->Branch("weight",&weight,"data/F");
  Tbkg->Branch("lowEweight",&lowEweight,"data/F");
  Tbkg->Branch("nueTag",&nueTag,"data/I");
  
  
  Tsig->Branch("mip_flag",&mip_flag,"mip_flag/F");
  Tsig->Branch("mip_energy",&mip_energy,"mip_energy/F");
  Tsig->Branch("mip_n_end_reduction",&mip_n_end_reduction,"mip_n_end_reduction/F");
  Tsig->Branch("mip_n_first_mip",&mip_n_first_mip,"mip_n_first_mip/F");
  Tsig->Branch("mip_n_first_non_mip",&mip_n_first_non_mip,"mip_n_first_non_mip/F");
  Tsig->Branch("mip_n_first_non_mip_1",&mip_n_first_non_mip_1,"mip_n_first_non_mip_1/F");
  Tsig->Branch("mip_n_first_non_mip_2",&mip_n_first_non_mip_2,"mip_n_first_non_mip_2/F");
  Tsig->Branch("mip_vec_dQ_dx_0",&mip_vec_dQ_dx_0,"mip_vec_dQ_dx_0/F");
  Tsig->Branch("mip_vec_dQ_dx_1",&mip_vec_dQ_dx_1,"mip_vec_dQ_dx_1/F");
  Tsig->Branch("mip_max_dQ_dx_sample",&mip_max_dQ_dx_sample,"mip_max_dQ_dx_sample/F");
  Tsig->Branch("mip_n_below_threshold",&mip_n_below_threshold,"mip_n_below_threshold/F");
  Tsig->Branch("mip_n_below_zero",&mip_n_below_zero,"mip_n_below_zero/F");
  Tsig->Branch("mip_n_lowest",&mip_n_lowest,"mip_n_lowest/F");
  Tsig->Branch("mip_n_highest",&mip_n_highest,"mip_n_highest/F");
  Tsig->Branch("mip_lowest_dQ_dx",&mip_lowest_dQ_dx,"mip_lowest_dQ_dx/F");
  Tsig->Branch("mip_highest_dQ_dx",&mip_highest_dQ_dx,"mip_highest_dQ_dx/F");
  Tsig->Branch("mip_medium_dQ_dx",&mip_medium_dQ_dx,"mip_medium_dQ_dx/F");
  Tsig->Branch("mip_stem_length",&mip_stem_length,"mip_stem_length/F");
  Tsig->Branch("mip_length_main",&mip_length_main,"mip_length_main/F");
  Tsig->Branch("mip_length_total",&mip_length_total,"mip_length_total/F");
  Tsig->Branch("mip_angle_beam",&mip_angle_beam,"mip_angle_beam/F");
  Tsig->Branch("mip_iso_angle",&mip_iso_angle,"mip_iso_angle/F");
  Tsig->Branch("mip_n_vertex",&mip_n_vertex,"mip_n_vertex/F");
  Tsig->Branch("mip_n_good_tracks",&mip_n_good_tracks,"mip_n_good_tracks/F");
  Tsig->Branch("mip_E_indirect_max_energy",&mip_E_indirect_max_energy,"mip_E_indirect_max_energy/F");
  Tsig->Branch("mip_flag_all_above",&mip_flag_all_above,"mip_flag_all_above/F");
  Tsig->Branch("mip_min_dQ_dx_5",&mip_min_dQ_dx_5,"mip_min_dQ_dx_5/F");
  Tsig->Branch("mip_n_other_vertex",&mip_n_other_vertex,"mip_n_other_vertex/F");
  Tsig->Branch("mip_n_stem_size",&mip_n_stem_size,"mip_n_stem_size/F");
  Tsig->Branch("mip_flag_stem_trajectory",&mip_flag_stem_trajectory,"mip_flag_stem_trajectory/F");
  Tsig->Branch("mip_min_dis",&mip_min_dis,"mip_min_dis/F");

  
  Tbkg->Branch("mip_flag",&mip_flag,"mip_flag/F");
  Tbkg->Branch("mip_energy",&mip_energy,"mip_energy/F");
  Tbkg->Branch("mip_n_end_reduction",&mip_n_end_reduction,"mip_n_end_reduction/F");
  Tbkg->Branch("mip_n_first_mip",&mip_n_first_mip,"mip_n_first_mip/F");
  Tbkg->Branch("mip_n_first_non_mip",&mip_n_first_non_mip,"mip_n_first_non_mip/F");
  Tbkg->Branch("mip_n_first_non_mip_1",&mip_n_first_non_mip_1,"mip_n_first_non_mip_1/F");
  Tbkg->Branch("mip_n_first_non_mip_2",&mip_n_first_non_mip_2,"mip_n_first_non_mip_2/F");
  Tbkg->Branch("mip_vec_dQ_dx_0",&mip_vec_dQ_dx_0,"mip_vec_dQ_dx_0/F");
  Tbkg->Branch("mip_vec_dQ_dx_1",&mip_vec_dQ_dx_1,"mip_vec_dQ_dx_1/F");
  Tbkg->Branch("mip_max_dQ_dx_sample",&mip_max_dQ_dx_sample,"mip_max_dQ_dx_sample/F");
  Tbkg->Branch("mip_n_below_threshold",&mip_n_below_threshold,"mip_n_below_threshold/F");
  Tbkg->Branch("mip_n_below_zero",&mip_n_below_zero,"mip_n_below_zero/F");
  Tbkg->Branch("mip_n_lowest",&mip_n_lowest,"mip_n_lowest/F");
  Tbkg->Branch("mip_n_highest",&mip_n_highest,"mip_n_highest/F");
  Tbkg->Branch("mip_lowest_dQ_dx",&mip_lowest_dQ_dx,"mip_lowest_dQ_dx/F");
  Tbkg->Branch("mip_highest_dQ_dx",&mip_highest_dQ_dx,"mip_highest_dQ_dx/F");
  Tbkg->Branch("mip_medium_dQ_dx",&mip_medium_dQ_dx,"mip_medium_dQ_dx/F");
  Tbkg->Branch("mip_stem_length",&mip_stem_length,"mip_stem_length/F");
  Tbkg->Branch("mip_length_main",&mip_length_main,"mip_length_main/F");
  Tbkg->Branch("mip_length_total",&mip_length_total,"mip_length_total/F");
  Tbkg->Branch("mip_angle_beam",&mip_angle_beam,"mip_angle_beam/F");
  Tbkg->Branch("mip_iso_angle",&mip_iso_angle,"mip_iso_angle/F");
  Tbkg->Branch("mip_n_vertex",&mip_n_vertex,"mip_n_vertex/F");
  Tbkg->Branch("mip_n_good_tracks",&mip_n_good_tracks,"mip_n_good_tracks/F");
  Tbkg->Branch("mip_E_indirect_max_energy",&mip_E_indirect_max_energy,"mip_E_indirect_max_energy/F");
  Tbkg->Branch("mip_flag_all_above",&mip_flag_all_above,"mip_flag_all_above/F");
  Tbkg->Branch("mip_min_dQ_dx_5",&mip_min_dQ_dx_5,"mip_min_dQ_dx_5/F");
  Tbkg->Branch("mip_n_other_vertex",&mip_n_other_vertex,"mip_n_other_vertex/F");
  Tbkg->Branch("mip_n_stem_size",&mip_n_stem_size,"mip_n_stem_size/F");
  Tbkg->Branch("mip_flag_stem_trajectory",&mip_flag_stem_trajectory,"mip_flag_stem_trajectory/F");
  Tbkg->Branch("mip_min_dis",&mip_min_dis,"mip_min_dis/F");
  
  
  Tsig->Branch("mip_bdt",&bdt_value,"data/F");
  Tbkg->Branch("mip_bdt",&bdt_value,"data/F");
  
  TMVA::Reader *reader = new TMVA::Reader();
  reader->AddVariable("mip_energy",&mip_energy);
  reader->AddVariable("mip_n_end_reduction",&mip_n_end_reduction);
  reader->AddVariable("mip_n_first_mip",&mip_n_first_mip);
  reader->AddVariable("mip_n_first_non_mip",&mip_n_first_non_mip);
  reader->AddVariable("mip_n_first_non_mip_1",&mip_n_first_non_mip_1);
  reader->AddVariable("mip_n_first_non_mip_2",&mip_n_first_non_mip_2);
  reader->AddVariable("mip_vec_dQ_dx_0",&mip_vec_dQ_dx_0);
  reader->AddVariable("mip_vec_dQ_dx_1",&mip_vec_dQ_dx_1);
  reader->AddVariable("mip_max_dQ_dx_sample",&mip_max_dQ_dx_sample);
  reader->AddVariable("mip_n_below_threshold",&mip_n_below_threshold);
  reader->AddVariable("mip_n_below_zero",&mip_n_below_zero);
  reader->AddVariable("mip_n_lowest",&mip_n_lowest);
  reader->AddVariable("mip_n_highest",&mip_n_highest);
  reader->AddVariable("mip_lowest_dQ_dx",&mip_lowest_dQ_dx);
  reader->AddVariable("mip_highest_dQ_dx",&mip_highest_dQ_dx);
  reader->AddVariable("mip_medium_dQ_dx",&mip_medium_dQ_dx);
  reader->AddVariable("mip_stem_length",&mip_stem_length);
  reader->AddVariable("mip_length_main",&mip_length_main);
  reader->AddVariable("mip_length_total",&mip_length_total);
  reader->AddVariable("mip_angle_beam",&mip_angle_beam);
  reader->AddVariable("mip_iso_angle",&mip_iso_angle);
  reader->AddVariable("mip_n_vertex",&mip_n_vertex);
  reader->AddVariable("mip_n_good_tracks",&mip_n_good_tracks);
  reader->AddVariable("mip_E_indirect_max_energy",&mip_E_indirect_max_energy);
  reader->AddVariable("mip_flag_all_above",&mip_flag_all_above);
  reader->AddVariable("mip_min_dQ_dx_5",&mip_min_dQ_dx_5);
  reader->AddVariable("mip_n_other_vertex",&mip_n_other_vertex);
  reader->AddVariable("mip_n_stem_size",&mip_n_stem_size);
  reader->AddVariable("mip_flag_stem_trajectory",&mip_flag_stem_trajectory);
  reader->AddVariable("mip_min_dis",&mip_min_dis);
  
  reader->BookMVA( "MyBDT", "dataset/weights/Test_BDT.weights.xml");
  

  for (Int_t i=0;i!=sig->GetEntries();i++){
    sig->GetEntry(i);
    bdt_value = reader->EvaluateMVA("MyBDT");
    Tsig->Fill();
  }
  
  for (Int_t i=0;i!=bkg->GetEntries();i++){
    bkg->GetEntry(i);
    bdt_value = reader->EvaluateMVA("MyBDT");
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
    std::cout << "Prepare the initial rootfile" << std::endl;
    convert_file();
  }else if (process == 2){
    std::cout << "testing BDT in ROOT TMVA (1st round)" << std::endl;
    Run_r1();
  }else if (process == 3){
    std::cout << "testing BDT in ROOT TMVA (2nd round)" << std::endl;
    Run_r2();
  }
  
  return 1;
    
}

