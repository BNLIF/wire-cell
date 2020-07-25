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
void TestEvaluate(TString filename);


void Run_r2();
void InitBDT_r2();


void convert_file();


void convert_file(){
  TFile *file = new TFile("bdtfile_0718.root");
  TTree *sig = (TTree*)file->Get("sig");
  TTree *bkg = (TTree*)file->Get("bkg");

  
  int run, event;

  sig->SetBranchAddress("run",&run);
  sig->SetBranchAddress("event",&event);

  bkg->SetBranchAddress("run",&run);
  bkg->SetBranchAddress("event",&event);

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

   int truth_inFV;
  int truth_CC;
  int truth_nue;
  int truth_cosmic;

  sig->SetBranchAddress("truth_inFV",&truth_inFV);
  sig->SetBranchAddress("truth_CC",&truth_CC);
  sig->SetBranchAddress("truth_nue",&truth_nue);
  sig->SetBranchAddress("truth_cosmic",&truth_cosmic);

  bkg->SetBranchAddress("truth_inFV",&truth_inFV);
  bkg->SetBranchAddress("truth_CC",&truth_CC);
  bkg->SetBranchAddress("truth_nue",&truth_nue);
  bkg->SetBranchAddress("truth_cosmic",&truth_cosmic);
  
  // shower to wall
  double stw_1_energy;
  double stw_1_dis;
  double stw_1_dQ_dx;
  int stw_1_flag_single_shower;
  int stw_1_n_pi0;
  int stw_1_num_valid_tracks;
  int stw_1_flag;

  // int spt_flag_single_shower;
  //double spt_energy;
  double spt_shower_main_length;
  double spt_shower_total_length;
  double spt_angle_beam;
  double spt_angle_vertical;
  double spt_max_dQ_dx;
  double spt_angle_beam_1;
  double spt_angle_drift;
  double spt_angle_drift_1;
  int spt_num_valid_tracks;
  double spt_n_vtx_segs;
  double spt_max_length;
  int spt_flag;

  sig->SetBranchAddress("spt_shower_main_length", &spt_shower_main_length);
  sig->SetBranchAddress("spt_shower_total_length", &spt_shower_total_length);
  sig->SetBranchAddress("spt_angle_beam", &spt_angle_beam);
  sig->SetBranchAddress("spt_angle_vertical", &spt_angle_vertical);
  sig->SetBranchAddress("spt_max_dQ_dx", &spt_max_dQ_dx);
  sig->SetBranchAddress("spt_angle_beam_1", &spt_angle_beam_1);
  sig->SetBranchAddress("spt_angle_drift", &spt_angle_drift);
  sig->SetBranchAddress("spt_angle_drift_1", &spt_angle_drift_1);
  sig->SetBranchAddress("spt_num_valid_tracks", &spt_num_valid_tracks);
  sig->SetBranchAddress("spt_n_vtx_segs", &spt_n_vtx_segs);
  sig->SetBranchAddress("spt_max_length", &spt_max_length);
  sig->SetBranchAddress("spt_flag", &spt_flag);
  
  sig->SetBranchAddress("stw_1_energy",&stw_1_energy);
  sig->SetBranchAddress("stw_1_dis",&stw_1_dis);
  sig->SetBranchAddress("stw_1_dQ_dx",&stw_1_dQ_dx);
  sig->SetBranchAddress("stw_1_flag_single_shower",&stw_1_flag_single_shower);
  sig->SetBranchAddress("stw_1_n_pi0",&stw_1_n_pi0);
  sig->SetBranchAddress("stw_1_num_valid_tracks",&stw_1_num_valid_tracks);
  sig->SetBranchAddress("stw_1_flag",&stw_1_flag);
  

  bkg->SetBranchAddress("spt_shower_main_length", &spt_shower_main_length);
  bkg->SetBranchAddress("spt_shower_total_length", &spt_shower_total_length);
  bkg->SetBranchAddress("spt_angle_beam", &spt_angle_beam);
  bkg->SetBranchAddress("spt_angle_vertical", &spt_angle_vertical);
  bkg->SetBranchAddress("spt_max_dQ_dx", &spt_max_dQ_dx);
  bkg->SetBranchAddress("spt_angle_beam_1", &spt_angle_beam_1);
  bkg->SetBranchAddress("spt_angle_drift", &spt_angle_drift);
  bkg->SetBranchAddress("spt_angle_drift_1", &spt_angle_drift_1);
  bkg->SetBranchAddress("spt_num_valid_tracks", &spt_num_valid_tracks);
  bkg->SetBranchAddress("spt_n_vtx_segs", &spt_n_vtx_segs);
  bkg->SetBranchAddress("spt_max_length", &spt_max_length);
  bkg->SetBranchAddress("spt_flag", &spt_flag);
  
  bkg->SetBranchAddress("stw_1_energy",&stw_1_energy);
  bkg->SetBranchAddress("stw_1_dis",&stw_1_dis);
  bkg->SetBranchAddress("stw_1_dQ_dx",&stw_1_dQ_dx);
  bkg->SetBranchAddress("stw_1_flag_single_shower",&stw_1_flag_single_shower);
  bkg->SetBranchAddress("stw_1_n_pi0",&stw_1_n_pi0);
  bkg->SetBranchAddress("stw_1_num_valid_tracks",&stw_1_num_valid_tracks);
  bkg->SetBranchAddress("stw_1_flag",&stw_1_flag);
  
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

  Tsig->Branch("run",&run,"data/I");
  Tsig->Branch("event",&event,"data/I");

  Tbkg->Branch("run",&run,"data/I");
  Tbkg->Branch("event",&event,"data/I");

  Tsig->Branch("truth_inFV",&truth_inFV,"data/I");
  Tsig->Branch("truth_CC",&truth_CC,"data/I");
  Tsig->Branch("truth_nue",&truth_nue,"data/I");
  Tsig->Branch("truth_cosmic",&truth_cosmic,"data/I");

  Tbkg->Branch("truth_inFV",&truth_inFV,"data/I");
  Tbkg->Branch("truth_CC",&truth_CC,"data/I");
  Tbkg->Branch("truth_nue",&truth_nue,"data/I");
  Tbkg->Branch("truth_cosmic",&truth_cosmic,"data/I");
  
  
  float stw_1_energy_f;
  float stw_1_dis_f;
  float stw_1_dQ_dx_f;
  float stw_1_flag_single_shower_f;
  float stw_1_n_pi0_f;
  float stw_1_num_valid_tracks_f;
  float stw_1_flag_f;
  float spt_shower_main_length_f;
  float spt_shower_total_length_f;
  float spt_angle_beam_f;
  float spt_angle_vertical_f;
  float spt_max_dQ_dx_f;
  float spt_angle_beam_1_f;
  float spt_angle_drift_f;
  float spt_angle_drift_1_f;
  float spt_num_valid_tracks_f;
  float spt_n_vtx_segs_f;
  float spt_max_length_f;
  float spt_flag_f;

  Tsig->Branch("stw_1_energy",&stw_1_energy_f,"stw_1_energy/F");
  Tsig->Branch("stw_1_dis",&stw_1_dis_f,"stw_1_dis/F");
  Tsig->Branch("stw_1_dQ_dx",&stw_1_dQ_dx_f,"stw_1_dQ_dx/F");
  Tsig->Branch("stw_1_flag_single_shower",&stw_1_flag_single_shower_f,"stw_1_flag_single_shower/F");
  Tsig->Branch("stw_1_n_pi0",&stw_1_n_pi0_f,"stw_1_n_pi0/F");
  Tsig->Branch("stw_1_num_valid_tracks",&stw_1_num_valid_tracks_f,"stw_1_num_valid_tracks/F");
  Tsig->Branch("stw_1_flag",&stw_1_flag_f,"stw_1_flag/F");

  Tsig->Branch("spt_shower_main_length", &spt_shower_main_length_f, "spt_shower_main_length/F");
  Tsig->Branch("spt_shower_total_length", &spt_shower_total_length_f, "spt_shower_total_length/F");
  Tsig->Branch("spt_angle_beam", &spt_angle_beam_f, "spt_angle_beam/F");
  Tsig->Branch("spt_angle_vertical", &spt_angle_vertical_f, "spt_angle_vertical/F");
  Tsig->Branch("spt_max_dQ_dx", &spt_max_dQ_dx_f, "spt_max_dQ_dx/F");
  Tsig->Branch("spt_angle_beam_1", &spt_angle_beam_1_f, "spt_angle_beam_1/F");
  Tsig->Branch("spt_angle_drift", &spt_angle_drift_f, "spt_angle_drift/F");
  Tsig->Branch("spt_angle_drift_1", &spt_angle_drift_1_f, "spt_angle_drift_1/F");
  Tsig->Branch("spt_num_valid_tracks", &spt_num_valid_tracks_f, "spt_num_valid_tracks/F");
  Tsig->Branch("spt_n_vtx_segs", &spt_n_vtx_segs_f, "spt_n_vtx_segs/F");
  Tsig->Branch("spt_max_length", &spt_max_length_f, "spt_max_length/F");
  Tsig->Branch("spt_flag", &spt_flag_f, "spt_flag/F");
  
  Tbkg->Branch("stw_1_energy",&stw_1_energy_f,"stw_1_energy/F");
  Tbkg->Branch("stw_1_dis",&stw_1_dis_f,"stw_1_dis/F");
  Tbkg->Branch("stw_1_dQ_dx",&stw_1_dQ_dx_f,"stw_1_dQ_dx/F");
  Tbkg->Branch("stw_1_flag_single_shower",&stw_1_flag_single_shower_f,"stw_1_flag_single_shower/F");
  Tbkg->Branch("stw_1_n_pi0",&stw_1_n_pi0_f,"stw_1_n_pi0/F");
  Tbkg->Branch("stw_1_num_valid_tracks",&stw_1_num_valid_tracks_f,"stw_1_num_valid_tracks/F");
  Tbkg->Branch("stw_1_flag",&stw_1_flag_f,"stw_1_flag/F");

  Tbkg->Branch("spt_shower_main_length", &spt_shower_main_length_f, "spt_shower_main_length/F");
  Tbkg->Branch("spt_shower_total_length", &spt_shower_total_length_f, "spt_shower_total_length/F");
  Tbkg->Branch("spt_angle_beam", &spt_angle_beam_f, "spt_angle_beam/F");
  Tbkg->Branch("spt_angle_vertical", &spt_angle_vertical_f, "spt_angle_vertical/F");
  Tbkg->Branch("spt_max_dQ_dx", &spt_max_dQ_dx_f, "spt_max_dQ_dx/F");
  Tbkg->Branch("spt_angle_beam_1", &spt_angle_beam_1_f, "spt_angle_beam_1/F");
  Tbkg->Branch("spt_angle_drift", &spt_angle_drift_f, "spt_angle_drift/F");
  Tbkg->Branch("spt_angle_drift_1", &spt_angle_drift_1_f, "spt_angle_drift_1/F");
  Tbkg->Branch("spt_num_valid_tracks", &spt_num_valid_tracks_f, "spt_num_valid_tracks/F");
  Tbkg->Branch("spt_n_vtx_segs", &spt_n_vtx_segs_f, "spt_n_vtx_segs/F");
  Tbkg->Branch("spt_max_length", &spt_max_length_f, "spt_max_length/F");
  Tbkg->Branch("spt_flag", &spt_flag_f, "spt_flag/F");
  
  
  for (Int_t i=0;i!=sig->GetEntries();i++){
    sig->GetEntry(i);

    stw_1_energy_f = stw_1_energy;
    stw_1_dis_f = stw_1_dis;
    stw_1_dQ_dx_f = stw_1_dQ_dx;
    stw_1_flag_single_shower_f = stw_1_flag_single_shower;
    stw_1_n_pi0_f = stw_1_n_pi0;
    stw_1_num_valid_tracks_f = stw_1_num_valid_tracks;
    stw_1_flag_f = stw_1_flag;
    spt_shower_main_length_f = spt_shower_main_length;
    spt_shower_total_length_f = spt_shower_total_length;
    spt_angle_beam_f = spt_angle_beam;
    spt_angle_vertical_f = spt_angle_vertical;
    spt_max_dQ_dx_f = spt_max_dQ_dx;
    spt_angle_beam_1_f = spt_angle_beam_1;
    spt_angle_drift_f = spt_angle_drift;
    spt_angle_drift_1_f = spt_angle_drift_1;
    spt_num_valid_tracks_f = spt_num_valid_tracks;
    spt_n_vtx_segs_f = spt_n_vtx_segs;
    spt_max_length_f = spt_max_length;
    spt_flag_f = spt_flag;
    
    
    Tsig->Fill();
  }

   for (Int_t i=0;i!=bkg->GetEntries();i++){
    bkg->GetEntry(i);
    
    stw_1_energy_f = stw_1_energy;
    stw_1_dis_f = stw_1_dis;
    stw_1_dQ_dx_f = stw_1_dQ_dx;
    stw_1_flag_single_shower_f = stw_1_flag_single_shower;
    stw_1_n_pi0_f = stw_1_n_pi0;
    stw_1_num_valid_tracks_f = stw_1_num_valid_tracks;
    stw_1_flag_f = stw_1_flag;
    spt_shower_main_length_f = spt_shower_main_length;
    spt_shower_total_length_f = spt_shower_total_length;
    spt_angle_beam_f = spt_angle_beam;
    spt_angle_vertical_f = spt_angle_vertical;
    spt_max_dQ_dx_f = spt_max_dQ_dx;
    spt_angle_beam_1_f = spt_angle_beam_1;
    spt_angle_drift_f = spt_angle_drift;
    spt_angle_drift_1_f = spt_angle_drift_1;
    spt_num_valid_tracks_f = spt_num_valid_tracks;
    spt_n_vtx_segs_f = spt_n_vtx_segs;
    spt_max_length_f = spt_max_length;
    spt_flag_f = spt_flag;
    
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
    // float stw_1_energy_f;
    // float stw_1_dis_f;
    // float stw_1_dQ_dx_f;
    // float stw_1_flag_single_shower_f;
    // float stw_1_n_pi0_f;
    // float stw_1_num_valid_tracks_f;
    // float stw_1_flag_f;
    // float spt_shower_main_length_f;
    // float spt_shower_total_length_f;
    // float spt_angle_beam_f;
    // float spt_angle_vertical_f;
    // float spt_max_dQ_dx_f;
    // float spt_angle_beam_1_f;
    // float spt_angle_drift_f;
    // float spt_angle_drift_1_f;
    // float spt_num_valid_tracks_f;
    // float spt_n_vtx_segs_f;
    // float spt_max_length_f;
    // float spt_flag_f;
    
    dataloader->AddVariable("stw_1_energy","stw_1_energy","MeV",'F');
    dataloader->AddVariable("stw_1_dis","stw_1_dis","cm",'F');
    dataloader->AddVariable("stw_1_dQ_dx","stw_1_dQ_dx","MeV/cm",'F');
    dataloader->AddVariable("stw_1_flag_single_shower","stw_1_flag_single_shower","",'F');
    dataloader->AddVariable("stw_1_n_pi0","stw_1_n_pi0","",'F');
    dataloader->AddVariable("stw_1_num_valid_tracks","stw_1_num_valid_tracks","",'F');
    
    dataloader->AddVariable("spt_shower_main_length","spt_shower_main_length","cm",'F');
    dataloader->AddVariable("spt_shower_total_length","spt_shower_total_length","cm",'F');
    dataloader->AddVariable("spt_angle_beam","spt_angle_beam","",'F');
    dataloader->AddVariable("spt_angle_vertical","spt_angle_vertical","",'F');
    dataloader->AddVariable("spt_max_dQ_dx","spt_max_dQ_dx","MeV/cm",'F');
    dataloader->AddVariable("spt_angle_beam_1","spt_angle_beam_1","",'F');
    dataloader->AddVariable("spt_angle_drift","spt_angle_drift","",'F');
    dataloader->AddVariable("spt_angle_drift_1","spt_angle_drift_1","",'F');
    dataloader->AddVariable("spt_num_valid_tracks","spt_num_valid_tracks","",'F');
    dataloader->AddVariable("spt_n_vtx_segs","spt_n_vtx_segs","",'F');
    dataloader->AddVariable("spt_max_length","spt_max_length","cm",'F');

    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "1>0"; // 82 (stw) + 342 (spt) or 422 (total)/44194
    TCut mycut_b = "(spt_flag ==0 || stw_1_flag == 0) && (!(truth_nue==1 && truth_CC==1))"; // 404 (stw) + 1305 (spt) or 1663 (total)/21070
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=20000:"
					    "nTrain_Background=1400:"
					    "nTest_Signal=10000:"
					    "nTest_Background=211:"
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
   
    
    dataloader->AddVariable("stw_1_energy","stw_1_energy","MeV",'F');
    dataloader->AddVariable("stw_1_dis","stw_1_dis","cm",'F');
    dataloader->AddVariable("stw_1_dQ_dx","stw_1_dQ_dx","MeV/cm",'F');
    dataloader->AddVariable("stw_1_flag_single_shower","stw_1_flag_single_shower","",'F');
    dataloader->AddVariable("stw_1_n_pi0","stw_1_n_pi0","",'F');
    dataloader->AddVariable("stw_1_num_valid_tracks","stw_1_num_valid_tracks","",'F');
    
    dataloader->AddVariable("spt_shower_main_length","spt_shower_main_length","cm",'F');
    dataloader->AddVariable("spt_shower_total_length","spt_shower_total_length","cm",'F');
    dataloader->AddVariable("spt_angle_beam","spt_angle_beam","",'F');
    dataloader->AddVariable("spt_angle_vertical","spt_angle_vertical","",'F');
    dataloader->AddVariable("spt_max_dQ_dx","spt_max_dQ_dx","MeV/cm",'F');
    dataloader->AddVariable("spt_angle_beam_1","spt_angle_beam_1","",'F');
    dataloader->AddVariable("spt_angle_drift","spt_angle_drift","",'F');
    dataloader->AddVariable("spt_angle_drift_1","spt_angle_drift_1","",'F');
    dataloader->AddVariable("spt_num_valid_tracks","spt_num_valid_tracks","",'F');
    dataloader->AddVariable("spt_n_vtx_segs","spt_n_vtx_segs","",'F');
    dataloader->AddVariable("spt_max_length","spt_max_length","cm",'F');

    
    
    
    TTree *signalTree     = (TTree*)input->Get("sig");
  TTree *backgroundTree = (TTree*)input->Get("bkg");
  dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
  dataloader->AddBackgroundTree( backgroundTree, 1.0);
  // Set individual event weights (the variables must exist in the original TTree)
  dataloader->SetSignalWeightExpression( "weight * lowEweight" );
  dataloader->SetBackgroundWeightExpression( "weight " );
  
  // Apply additional cuts on the signal and background samples (can be different)
  //    TCut mycut_s = "pio_mip_id==0&&pio_filled==1&&pio_flag_pio==1"; // 831 events
  //    TCut mycut_b = "pio_mip_id==0&&pio_filled==1&&pio_flag_pio==1"; // 859 events
  
  TCut mycut_s = "1>0"; // 422/44194
  TCut mycut_b = "(stw_1_flag ==0 || spt_flag==0 || stw_spt_bdt <-0.1) && (!(truth_nue==1 && truth_CC==1))"; // 2181/21070
  
  dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					  "nTrain_Signal=20000:"
					  "nTrain_Background=1800:"
					  "nTest_Signal=10000:"
					  "nTest_Background=293:"
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

   int run, event;
  sig->SetBranchAddress("run",&run);
  sig->SetBranchAddress("event",&event);

  bkg->SetBranchAddress("run",&run);
  bkg->SetBranchAddress("event",&event);
  
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


  float stw_1_energy;
  float stw_1_dis;
  float stw_1_dQ_dx;
  float stw_1_flag_single_shower;
  float stw_1_n_pi0;
  float stw_1_num_valid_tracks;
  float stw_1_flag;
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

  sig->SetBranchAddress("spt_shower_main_length", &spt_shower_main_length);
  sig->SetBranchAddress("spt_shower_total_length", &spt_shower_total_length);
  sig->SetBranchAddress("spt_angle_beam", &spt_angle_beam);
  sig->SetBranchAddress("spt_angle_vertical", &spt_angle_vertical);
  sig->SetBranchAddress("spt_max_dQ_dx", &spt_max_dQ_dx);
  sig->SetBranchAddress("spt_angle_beam_1", &spt_angle_beam_1);
  sig->SetBranchAddress("spt_angle_drift", &spt_angle_drift);
  sig->SetBranchAddress("spt_angle_drift_1", &spt_angle_drift_1);
  sig->SetBranchAddress("spt_num_valid_tracks", &spt_num_valid_tracks);
  sig->SetBranchAddress("spt_n_vtx_segs", &spt_n_vtx_segs);
  sig->SetBranchAddress("spt_max_length", &spt_max_length);
  sig->SetBranchAddress("spt_flag", &spt_flag);
  
  sig->SetBranchAddress("stw_1_energy",&stw_1_energy);
  sig->SetBranchAddress("stw_1_dis",&stw_1_dis);
  sig->SetBranchAddress("stw_1_dQ_dx",&stw_1_dQ_dx);
  sig->SetBranchAddress("stw_1_flag_single_shower",&stw_1_flag_single_shower);
  sig->SetBranchAddress("stw_1_n_pi0",&stw_1_n_pi0);
  sig->SetBranchAddress("stw_1_num_valid_tracks",&stw_1_num_valid_tracks);
  sig->SetBranchAddress("stw_1_flag",&stw_1_flag);
  

  bkg->SetBranchAddress("spt_shower_main_length", &spt_shower_main_length);
  bkg->SetBranchAddress("spt_shower_total_length", &spt_shower_total_length);
  bkg->SetBranchAddress("spt_angle_beam", &spt_angle_beam);
  bkg->SetBranchAddress("spt_angle_vertical", &spt_angle_vertical);
  bkg->SetBranchAddress("spt_max_dQ_dx", &spt_max_dQ_dx);
  bkg->SetBranchAddress("spt_angle_beam_1", &spt_angle_beam_1);
  bkg->SetBranchAddress("spt_angle_drift", &spt_angle_drift);
  bkg->SetBranchAddress("spt_angle_drift_1", &spt_angle_drift_1);
  bkg->SetBranchAddress("spt_num_valid_tracks", &spt_num_valid_tracks);
  bkg->SetBranchAddress("spt_n_vtx_segs", &spt_n_vtx_segs);
  bkg->SetBranchAddress("spt_max_length", &spt_max_length);
  bkg->SetBranchAddress("spt_flag", &spt_flag);
  
  bkg->SetBranchAddress("stw_1_energy",&stw_1_energy);
  bkg->SetBranchAddress("stw_1_dis",&stw_1_dis);
  bkg->SetBranchAddress("stw_1_dQ_dx",&stw_1_dQ_dx);
  bkg->SetBranchAddress("stw_1_flag_single_shower",&stw_1_flag_single_shower);
  bkg->SetBranchAddress("stw_1_n_pi0",&stw_1_n_pi0);
  bkg->SetBranchAddress("stw_1_num_valid_tracks",&stw_1_num_valid_tracks);
  bkg->SetBranchAddress("stw_1_flag",&stw_1_flag);

   int truth_inFV;
  int truth_CC;
  int truth_nue;
  int truth_cosmic;

  sig->SetBranchAddress("truth_inFV",&truth_inFV);
  sig->SetBranchAddress("truth_CC",&truth_CC);
  sig->SetBranchAddress("truth_nue",&truth_nue);
  sig->SetBranchAddress("truth_cosmic",&truth_cosmic);

  bkg->SetBranchAddress("truth_inFV",&truth_inFV);
  bkg->SetBranchAddress("truth_CC",&truth_CC);
  bkg->SetBranchAddress("truth_nue",&truth_nue);
  bkg->SetBranchAddress("truth_cosmic",&truth_cosmic);
  
  TFile *new_file = new TFile(filename,"RECREATE");
  TTree *Tsig = new TTree("sig","sig");
  TTree *Tbkg = new TTree("bkg","bkg");
  Tsig->SetDirectory(new_file);
  Tbkg->SetDirectory(new_file);
  
  Tsig->Branch("run",&run,"data/I");
  Tsig->Branch("event",&event,"data/I");

  Tbkg->Branch("run",&run,"data/I");
  Tbkg->Branch("event",&event,"data/I");
  
  Tsig->Branch("trueEdep",&trueEdep,"data/F");
  Tsig->Branch("weight",&weight,"data/F");
  Tsig->Branch("lowEweight",&lowEweight,"data/F");
  Tsig->Branch("nueTag",&nueTag,"data/I");

  Tbkg->Branch("trueEdep",&trueEdep,"data/F");
  Tbkg->Branch("weight",&weight,"data/F");
  Tbkg->Branch("lowEweight",&lowEweight,"data/F");
  Tbkg->Branch("nueTag",&nueTag,"data/I");
  
  Tsig->Branch("truth_inFV",&truth_inFV,"data/I");
  Tsig->Branch("truth_CC",&truth_CC,"data/I");
  Tsig->Branch("truth_nue",&truth_nue,"data/I");
  Tsig->Branch("truth_cosmic",&truth_cosmic,"data/I");

  Tbkg->Branch("truth_inFV",&truth_inFV,"data/I");
  Tbkg->Branch("truth_CC",&truth_CC,"data/I");
  Tbkg->Branch("truth_nue",&truth_nue,"data/I");
  Tbkg->Branch("truth_cosmic",&truth_cosmic,"data/I");

  
  Tsig->Branch("stw_1_energy",&stw_1_energy,"stw_1_energy/F");
  Tsig->Branch("stw_1_dis",&stw_1_dis,"stw_1_dis/F");
  Tsig->Branch("stw_1_dQ_dx",&stw_1_dQ_dx,"stw_1_dQ_dx/F");
  Tsig->Branch("stw_1_flag_single_shower",&stw_1_flag_single_shower,"stw_1_flag_single_shower/F");
  Tsig->Branch("stw_1_n_pi0",&stw_1_n_pi0,"stw_1_n_pi0/F");
  Tsig->Branch("stw_1_num_valid_tracks",&stw_1_num_valid_tracks,"stw_1_num_valid_tracks/F");
  Tsig->Branch("stw_1_flag",&stw_1_flag,"stw_1_flag/F");

  Tsig->Branch("spt_shower_main_length", &spt_shower_main_length, "spt_shower_main_length/F");
  Tsig->Branch("spt_shower_total_length", &spt_shower_total_length, "spt_shower_total_length/F");
  Tsig->Branch("spt_angle_beam", &spt_angle_beam, "spt_angle_beam/F");
  Tsig->Branch("spt_angle_vertical", &spt_angle_vertical, "spt_angle_vertical/F");
  Tsig->Branch("spt_max_dQ_dx", &spt_max_dQ_dx, "spt_max_dQ_dx/F");
  Tsig->Branch("spt_angle_beam_1", &spt_angle_beam_1, "spt_angle_beam_1/F");
  Tsig->Branch("spt_angle_drift", &spt_angle_drift, "spt_angle_drift/F");
  Tsig->Branch("spt_angle_drift_1", &spt_angle_drift_1, "spt_angle_drift_1/F");
  Tsig->Branch("spt_num_valid_tracks", &spt_num_valid_tracks, "spt_num_valid_tracks/F");
  Tsig->Branch("spt_n_vtx_segs", &spt_n_vtx_segs, "spt_n_vtx_segs/F");
  Tsig->Branch("spt_max_length", &spt_max_length, "spt_max_length/F");
  Tsig->Branch("spt_flag", &spt_flag, "spt_flag/F");
  
  Tbkg->Branch("stw_1_energy",&stw_1_energy,"stw_1_energy/F");
  Tbkg->Branch("stw_1_dis",&stw_1_dis,"stw_1_dis/F");
  Tbkg->Branch("stw_1_dQ_dx",&stw_1_dQ_dx,"stw_1_dQ_dx/F");
  Tbkg->Branch("stw_1_flag_single_shower",&stw_1_flag_single_shower,"stw_1_flag_single_shower/F");
  Tbkg->Branch("stw_1_n_pi0",&stw_1_n_pi0,"stw_1_n_pi0/F");
  Tbkg->Branch("stw_1_num_valid_tracks",&stw_1_num_valid_tracks,"stw_1_num_valid_tracks/F");
  Tbkg->Branch("stw_1_flag",&stw_1_flag,"stw_1_flag/F");

  Tbkg->Branch("spt_shower_main_length", &spt_shower_main_length, "spt_shower_main_length/F");
  Tbkg->Branch("spt_shower_total_length", &spt_shower_total_length, "spt_shower_total_length/F");
  Tbkg->Branch("spt_angle_beam", &spt_angle_beam, "spt_angle_beam/F");
  Tbkg->Branch("spt_angle_vertical", &spt_angle_vertical, "spt_angle_vertical/F");
  Tbkg->Branch("spt_max_dQ_dx", &spt_max_dQ_dx, "spt_max_dQ_dx/F");
  Tbkg->Branch("spt_angle_beam_1", &spt_angle_beam_1, "spt_angle_beam_1/F");
  Tbkg->Branch("spt_angle_drift", &spt_angle_drift, "spt_angle_drift/F");
  Tbkg->Branch("spt_angle_drift_1", &spt_angle_drift_1, "spt_angle_drift_1/F");
  Tbkg->Branch("spt_num_valid_tracks", &spt_num_valid_tracks, "spt_num_valid_tracks/F");
  Tbkg->Branch("spt_n_vtx_segs", &spt_n_vtx_segs, "spt_n_vtx_segs/F");
  Tbkg->Branch("spt_max_length", &spt_max_length, "spt_max_length/F");
  Tbkg->Branch("spt_flag", &spt_flag, "spt_flag/F");

  Tsig->Branch("stw_spt_bdt", &bdt_value,"data/F");
  Tbkg->Branch("stw_spt_bdt", &bdt_value,"data/F");

  
  TMVA::Reader *reader = new TMVA::Reader();
  reader->AddVariable("stw_1_energy",&stw_1_energy);
  reader->AddVariable("stw_1_dis",&stw_1_dis);
  reader->AddVariable("stw_1_dQ_dx",&stw_1_dQ_dx);
  reader->AddVariable("stw_1_flag_single_shower",&stw_1_flag_single_shower);
  reader->AddVariable("stw_1_n_pi0",&stw_1_n_pi0);
  reader->AddVariable("stw_1_num_valid_tracks",&stw_1_num_valid_tracks);
  
  reader->AddVariable("spt_shower_main_length",&spt_shower_main_length);
  reader->AddVariable("spt_shower_total_length",&spt_shower_total_length);
  reader->AddVariable("spt_angle_beam",&spt_angle_beam);
  reader->AddVariable("spt_angle_vertical",&spt_angle_vertical);
  reader->AddVariable("spt_max_dQ_dx",&spt_max_dQ_dx);
  reader->AddVariable("spt_angle_beam_1",&spt_angle_beam_1);
  reader->AddVariable("spt_angle_drift",&spt_angle_drift);
  reader->AddVariable("spt_angle_drift_1",&spt_angle_drift_1);
  reader->AddVariable("spt_num_valid_tracks",&spt_num_valid_tracks);
  reader->AddVariable("spt_n_vtx_segs",&spt_n_vtx_segs);
  reader->AddVariable("spt_max_length",&spt_max_length);
  
  
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

