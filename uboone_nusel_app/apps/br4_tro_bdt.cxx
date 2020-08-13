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
  
  double br4_1_shower_main_length;
  double br4_1_shower_total_length;
  double br4_1_min_dis;
  double br4_1_energy;
  int br4_1_flag_avoid_muon_check;
  int br4_1_n_vtx_segs;
  int br4_1_n_main_segs;
    
  double br4_2_ratio_45;
  double br4_2_ratio_35;
  double br4_2_ratio_25;
  double br4_2_ratio_15;
  double br4_2_ratio1_45;
  double br4_2_ratio1_35;
  double br4_2_ratio1_25;
  double br4_2_ratio1_15;
  double br4_2_iso_angle;
  double br4_2_iso_angle1;
  double br4_2_angle;
    
  int br4_flag;

  double tro_3_stem_length;
  int tro_3_n_muon_segs;
  double tro_3_energy;
  int tro_3_flag;
  
  sig->SetBranchAddress("tro_3_stem_length",&tro_3_stem_length);
  sig->SetBranchAddress("tro_3_n_muon_segs",&tro_3_n_muon_segs);
  sig->SetBranchAddress("tro_3_energy",&tro_3_energy);
  sig->SetBranchAddress("tro_3_flag",&tro_3_flag);

  bkg->SetBranchAddress("tro_3_stem_length",&tro_3_stem_length);
  bkg->SetBranchAddress("tro_3_n_muon_segs",&tro_3_n_muon_segs);
  bkg->SetBranchAddress("tro_3_energy",&tro_3_energy);
  bkg->SetBranchAddress("tro_3_flag",&tro_3_flag);
  
  sig->SetBranchAddress("br4_1_shower_main_length", &br4_1_shower_main_length);
  sig->SetBranchAddress("br4_1_shower_total_length", &br4_1_shower_total_length);
  sig->SetBranchAddress("br4_1_min_dis", &br4_1_min_dis);
  sig->SetBranchAddress("br4_1_energy", &br4_1_energy);
  sig->SetBranchAddress("br4_1_flag_avoid_muon_check", &br4_1_flag_avoid_muon_check);
  sig->SetBranchAddress("br4_1_n_vtx_segs", &br4_1_n_vtx_segs);
  sig->SetBranchAddress("br4_1_br4_1_n_main_segs", &br4_1_n_main_segs);
    
  sig->SetBranchAddress("br4_2_ratio_45", &br4_2_ratio_45);
  sig->SetBranchAddress("br4_2_ratio_35", &br4_2_ratio_35);
  sig->SetBranchAddress("br4_2_ratio_25", &br4_2_ratio_25);
  sig->SetBranchAddress("br4_2_ratio_15", &br4_2_ratio_15);
  sig->SetBranchAddress("br4_2_ratio1_45", &br4_2_ratio1_45);
  sig->SetBranchAddress("br4_2_ratio1_35", &br4_2_ratio1_35);
  sig->SetBranchAddress("br4_2_ratio1_25", &br4_2_ratio1_25);
  sig->SetBranchAddress("br4_2_ratio1_15", &br4_2_ratio1_15);
  sig->SetBranchAddress("br4_2_iso_angle", &br4_2_iso_angle);
  sig->SetBranchAddress("br4_2_iso_angle1", &br4_2_iso_angle1);
  sig->SetBranchAddress("br4_2_angle", &br4_2_angle);
  
  sig->SetBranchAddress("br4_flag", &br4_flag);

  bkg->SetBranchAddress("br4_1_shower_main_length", &br4_1_shower_main_length);
  bkg->SetBranchAddress("br4_1_shower_total_length", &br4_1_shower_total_length);
  bkg->SetBranchAddress("br4_1_min_dis", &br4_1_min_dis);
  bkg->SetBranchAddress("br4_1_energy", &br4_1_energy);
  bkg->SetBranchAddress("br4_1_flag_avoid_muon_check", &br4_1_flag_avoid_muon_check);
  bkg->SetBranchAddress("br4_1_n_vtx_segs", &br4_1_n_vtx_segs);
  bkg->SetBranchAddress("br4_1_br4_1_n_main_segs", &br4_1_n_main_segs);
    
  bkg->SetBranchAddress("br4_2_ratio_45", &br4_2_ratio_45);
  bkg->SetBranchAddress("br4_2_ratio_35", &br4_2_ratio_35);
  bkg->SetBranchAddress("br4_2_ratio_25", &br4_2_ratio_25);
  bkg->SetBranchAddress("br4_2_ratio_15", &br4_2_ratio_15);
  bkg->SetBranchAddress("br4_2_ratio1_45", &br4_2_ratio1_45);
  bkg->SetBranchAddress("br4_2_ratio1_35", &br4_2_ratio1_35);
  bkg->SetBranchAddress("br4_2_ratio1_25", &br4_2_ratio1_25);
  bkg->SetBranchAddress("br4_2_ratio1_15", &br4_2_ratio1_15);
  bkg->SetBranchAddress("br4_2_iso_angle", &br4_2_iso_angle);
  bkg->SetBranchAddress("br4_2_iso_angle1", &br4_2_iso_angle1);
  bkg->SetBranchAddress("br4_2_angle", &br4_2_angle);
  
  bkg->SetBranchAddress("br4_flag", &br4_flag);
  
  
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

  Tsig->Branch("truth_inFV",&truth_inFV,"data/I");
  Tsig->Branch("truth_CC",&truth_CC,"data/I");
  Tsig->Branch("truth_nue",&truth_nue,"data/I");
  Tsig->Branch("truth_cosmic",&truth_cosmic,"data/I");

  Tbkg->Branch("truth_inFV",&truth_inFV,"data/I");
  Tbkg->Branch("truth_CC",&truth_CC,"data/I");
  Tbkg->Branch("truth_nue",&truth_nue,"data/I");
  Tbkg->Branch("truth_cosmic",&truth_cosmic,"data/I");
  
  float br4_1_shower_main_length_f;
  float br4_1_shower_total_length_f;
  float br4_1_min_dis_f;
  float br4_1_energy_f;
  float br4_1_flag_avoid_muon_check_f;
  float br4_1_n_vtx_segs_f;
  float br4_1_n_main_segs_f;
    
  float br4_2_ratio_45_f;
  float br4_2_ratio_35_f;
  float br4_2_ratio_25_f;
  float br4_2_ratio_15_f;
  float br4_2_ratio1_45_f;
  float br4_2_ratio1_35_f;
  float br4_2_ratio1_25_f;
  float br4_2_ratio1_15_f;
  float br4_2_iso_angle_f;
  float br4_2_iso_angle1_f;
  float br4_2_angle_f;
    
  float br4_flag_f;

  float tro_3_stem_length_f;
  float tro_3_n_muon_segs_f;
  float tro_3_energy_f;
  float tro_3_flag_f;


  Tsig->Branch("tro_3_stem_length",&tro_3_stem_length_f,"tro_3_stem_length/F");
  Tsig->Branch("tro_3_n_muon_segs",&tro_3_n_muon_segs_f,"tro_3_n_muon_segs/F");
  Tsig->Branch("tro_3_energy",&tro_3_energy_f,"tro_3_energy/F");
  Tsig->Branch("tro_3_flag",&tro_3_flag_f,"tro_3_flag/F");

  Tbkg->Branch("tro_3_stem_length",&tro_3_stem_length_f,"tro_3_stem_length/F");
  Tbkg->Branch("tro_3_n_muon_segs",&tro_3_n_muon_segs_f,"tro_3_n_muon_segs/F");
  Tbkg->Branch("tro_3_energy",&tro_3_energy_f,"tro_3_energy/F");
  Tbkg->Branch("tro_3_flag",&tro_3_flag_f,"tro_3_flag/F");
  

  
  Tsig->Branch("br4_1_shower_main_length", &br4_1_shower_main_length_f,"br4_1_shower_main_length/F");
  Tsig->Branch("br4_1_shower_total_length", &br4_1_shower_total_length_f,"br4_1_shower_total_length/F");
  Tsig->Branch("br4_1_min_dis", &br4_1_min_dis_f,"br4_1_min_dis/F");
  Tsig->Branch("br4_1_energy", &br4_1_energy_f,"br4_1_energy/F");
  Tsig->Branch("br4_1_flag_avoid_muon_check", &br4_1_flag_avoid_muon_check_f,"br4_1_flag_avoid_muon_check/F");
  Tsig->Branch("br4_1_n_vtx_segs", &br4_1_n_vtx_segs_f,"br4_1_n_vtx_segs/F");
  Tsig->Branch("br4_1_n_main_segs", &br4_1_n_main_segs_f,"br4_1_n_main_segs/F");
    
  Tsig->Branch("br4_2_ratio_45", &br4_2_ratio_45_f,"br4_2_ratio_45/F");
  Tsig->Branch("br4_2_ratio_35", &br4_2_ratio_35_f,"br4_2_ratio_35/F");
  Tsig->Branch("br4_2_ratio_25", &br4_2_ratio_25_f,"br4_2_ratio_25/F");
  Tsig->Branch("br4_2_ratio_15", &br4_2_ratio_15_f,"br4_2_ratio_15/F");
  Tsig->Branch("br4_2_ratio1_45", &br4_2_ratio1_45_f,"br4_2_ratio1_45/F");
  Tsig->Branch("br4_2_ratio1_35", &br4_2_ratio1_35_f,"br4_2_ratio1_35/F");
  Tsig->Branch("br4_2_ratio1_25", &br4_2_ratio1_25_f,"br4_2_ratio1_25/F");
  Tsig->Branch("br4_2_ratio1_15", &br4_2_ratio1_15_f,"br4_2_ratio1_15/F");
  Tsig->Branch("br4_2_iso_angle", &br4_2_iso_angle_f,"br4_2_iso_angle/F");
  Tsig->Branch("br4_2_iso_angle1", &br4_2_iso_angle1_f,"br4_2_iso_angle1/F");
  Tsig->Branch("br4_2_angle", &br4_2_angle_f,"br4_2_angle/F");
    
  Tsig->Branch("br4_flag", &br4_flag_f,"br4_flag/F");

  
  Tbkg->Branch("br4_1_shower_main_length", &br4_1_shower_main_length_f,"br4_1_shower_main_length/F");
  Tbkg->Branch("br4_1_shower_total_length", &br4_1_shower_total_length_f,"br4_1_shower_total_length/F");
  Tbkg->Branch("br4_1_min_dis", &br4_1_min_dis_f,"br4_1_min_dis/F");
  Tbkg->Branch("br4_1_energy", &br4_1_energy_f,"br4_1_energy/F");
  Tbkg->Branch("br4_1_flag_avoid_muon_check", &br4_1_flag_avoid_muon_check_f,"br4_1_flag_avoid_muon_check/F");
  Tbkg->Branch("br4_1_n_vtx_segs", &br4_1_n_vtx_segs_f,"br4_1_n_vtx_segs/F");
  Tbkg->Branch("br4_1_n_main_segs", &br4_1_n_main_segs_f,"br4_1_n_main_segs/F");
    
  Tbkg->Branch("br4_2_ratio_45", &br4_2_ratio_45_f,"br4_2_ratio_45/F");
  Tbkg->Branch("br4_2_ratio_35", &br4_2_ratio_35_f,"br4_2_ratio_35/F");
  Tbkg->Branch("br4_2_ratio_25", &br4_2_ratio_25_f,"br4_2_ratio_25/F");
  Tbkg->Branch("br4_2_ratio_15", &br4_2_ratio_15_f,"br4_2_ratio_15/F");
  Tbkg->Branch("br4_2_ratio1_45", &br4_2_ratio1_45_f,"br4_2_ratio1_45/F");
  Tbkg->Branch("br4_2_ratio1_35", &br4_2_ratio1_35_f,"br4_2_ratio1_35/F");
  Tbkg->Branch("br4_2_ratio1_25", &br4_2_ratio1_25_f,"br4_2_ratio1_25/F");
  Tbkg->Branch("br4_2_ratio1_15", &br4_2_ratio1_15_f,"br4_2_ratio1_15/F");
  Tbkg->Branch("br4_2_iso_angle", &br4_2_iso_angle_f,"br4_2_iso_angle/F");
  Tbkg->Branch("br4_2_iso_angle1", &br4_2_iso_angle1_f,"br4_2_iso_angle1/F");
  Tbkg->Branch("br4_2_angle", &br4_2_angle_f,"br4_2_angle/F");
    
  Tbkg->Branch("br4_flag", &br4_flag_f,"br4_flag/F");
  
  for (Int_t i=0;i!=sig->GetEntries();i++){
    sig->GetEntry(i);

    br4_1_shower_main_length_f = br4_1_shower_main_length;
    br4_1_shower_total_length_f = br4_1_shower_total_length;
    br4_1_min_dis_f = br4_1_min_dis;
    br4_1_energy_f = br4_1_energy;
    br4_1_flag_avoid_muon_check_f = br4_1_flag_avoid_muon_check;
    br4_1_n_vtx_segs_f = br4_1_n_vtx_segs;
    br4_1_n_main_segs_f = br4_1_n_main_segs;
    
    br4_2_ratio_45_f = br4_2_ratio_45;
    br4_2_ratio_35_f = br4_2_ratio_35;
    br4_2_ratio_25_f = br4_2_ratio_25;
    br4_2_ratio_15_f = br4_2_ratio_15;
    br4_2_ratio1_45_f = br4_2_ratio1_45;
    br4_2_ratio1_35_f = br4_2_ratio1_35;
    br4_2_ratio1_25_f = br4_2_ratio1_25;
    br4_2_ratio1_15_f = br4_2_ratio1_15;
    br4_2_iso_angle_f = br4_2_iso_angle;
    br4_2_iso_angle1_f = br4_2_iso_angle1;
    br4_2_angle_f = br4_2_angle;
    
    br4_flag_f = br4_flag;

    tro_3_stem_length_f = tro_3_stem_length;
    tro_3_n_muon_segs_f = tro_3_n_muon_segs;
    tro_3_energy_f = tro_3_energy;
    tro_3_flag_f = tro_3_flag;
    
    Tsig->Fill();
  }

   for (Int_t i=0;i!=bkg->GetEntries();i++){
    bkg->GetEntry(i);
    
    br4_1_shower_main_length_f = br4_1_shower_main_length;
    br4_1_shower_total_length_f = br4_1_shower_total_length;
    br4_1_min_dis_f = br4_1_min_dis;
    br4_1_energy_f = br4_1_energy;
    br4_1_flag_avoid_muon_check_f = br4_1_flag_avoid_muon_check;
    br4_1_n_vtx_segs_f = br4_1_n_vtx_segs;
    br4_1_n_main_segs_f = br4_1_n_main_segs;
    
    br4_2_ratio_45_f = br4_2_ratio_45;
    br4_2_ratio_35_f = br4_2_ratio_35;
    br4_2_ratio_25_f = br4_2_ratio_25;
    br4_2_ratio_15_f = br4_2_ratio_15;
    br4_2_ratio1_45_f = br4_2_ratio1_45;
    br4_2_ratio1_35_f = br4_2_ratio1_35;
    br4_2_ratio1_25_f = br4_2_ratio1_25;
    br4_2_ratio1_15_f = br4_2_ratio1_15;
    br4_2_iso_angle_f = br4_2_iso_angle;
    br4_2_iso_angle1_f = br4_2_iso_angle1;
    br4_2_angle_f = br4_2_angle;
    
    br4_flag_f = br4_flag;
    
    tro_3_stem_length_f = tro_3_stem_length;
    tro_3_n_muon_segs_f = tro_3_n_muon_segs;
    tro_3_energy_f = tro_3_energy;
    tro_3_flag_f = tro_3_flag;
    
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

   
    
    dataloader->AddVariable("br4_1_shower_main_length","br4_1_shower_main_length","cm",'F');
    dataloader->AddVariable("br4_1_shower_total_length","br4_1_shower_total_length","cm",'F');
    dataloader->AddVariable("br4_1_min_dis","br4_1_min_dis","cm",'F');
    dataloader->AddVariable("br4_1_energy","br4_1_energy","MeV",'F');
    dataloader->AddVariable("br4_1_flag_avoid_muon_check","br4_1_flag_avoid_muon_check","",'F');
    dataloader->AddVariable("br4_1_n_vtx_segs","br4_1_n_vtx_segs","",'F');
    dataloader->AddVariable("br4_1_n_main_segs","br4_1_n_main_segs","",'F');

    dataloader->AddVariable("br4_2_ratio_45","br4_2_ratio_45","",'F');
    dataloader->AddVariable("br4_2_ratio_35","br4_2_ratio_35","",'F');
    dataloader->AddVariable("br4_2_ratio_25","br4_2_ratio_25","",'F');
    dataloader->AddVariable("br4_2_ratio_15","br4_2_ratio_15","",'F');
    dataloader->AddVariable("br4_2_ratio1_45","br4_2_ratio1_45","",'F');
    dataloader->AddVariable("br4_2_ratio1_35","br4_2_ratio1_35","",'F');
    dataloader->AddVariable("br4_2_ratio1_25","br4_2_ratio1_25","",'F');
    dataloader->AddVariable("br4_2_ratio1_15","br4_2_ratio1_15","",'F');
    dataloader->AddVariable("br4_2_iso_angle","br4_2_iso_angle","deg",'F');
    dataloader->AddVariable("br4_2_iso_angle1","br4_2_iso_angle1","deg",'F');
    dataloader->AddVariable("br4_2_angle","br4_2_angle","deg",'F');
    
    dataloader->AddVariable("tro_3_stem_length","tro_3_stem_length","cm",'F');
    dataloader->AddVariable("tro_3_n_muon_segs","tro_3_n_muon_segs","",'F');
    
    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "1>0"; //  1903/44914
    TCut mycut_b = "(br4_flag==0 || tro_3_flag==0)&& (!(truth_nue==1 && truth_CC==1))"; //  5837/21070
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=20000:"
					    "nTrain_Background=5200:"
					    "nTest_Signal=10000:"
					    "nTest_Background=637:"
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

   
    
    dataloader->AddVariable("br4_1_shower_main_length","br4_1_shower_main_length","cm",'F');
    dataloader->AddVariable("br4_1_shower_total_length","br4_1_shower_total_length","cm",'F');
    dataloader->AddVariable("br4_1_min_dis","br4_1_min_dis","cm",'F');
    dataloader->AddVariable("br4_1_energy","br4_1_energy","MeV",'F');
    dataloader->AddVariable("br4_1_flag_avoid_muon_check","br4_1_flag_avoid_muon_check","",'F');
    dataloader->AddVariable("br4_1_n_vtx_segs","br4_1_n_vtx_segs","",'F');
    dataloader->AddVariable("br4_1_n_main_segs","br4_1_n_main_segs","",'F');

    dataloader->AddVariable("br4_2_ratio_45","br4_2_ratio_45","",'F');
    dataloader->AddVariable("br4_2_ratio_35","br4_2_ratio_35","",'F');
    dataloader->AddVariable("br4_2_ratio_25","br4_2_ratio_25","",'F');
    dataloader->AddVariable("br4_2_ratio_15","br4_2_ratio_15","",'F');
    dataloader->AddVariable("br4_2_ratio1_45","br4_2_ratio1_45","",'F');
    dataloader->AddVariable("br4_2_ratio1_35","br4_2_ratio1_35","",'F');
    dataloader->AddVariable("br4_2_ratio1_25","br4_2_ratio1_25","",'F');
    dataloader->AddVariable("br4_2_ratio1_15","br4_2_ratio1_15","",'F');
    dataloader->AddVariable("br4_2_iso_angle","br4_2_iso_angle","deg",'F');
    dataloader->AddVariable("br4_2_iso_angle1","br4_2_iso_angle1","deg",'F');
    dataloader->AddVariable("br4_2_angle","br4_2_angle","deg",'F');

    dataloader->AddVariable("tro_3_stem_length","tro_3_stem_length","cm",'F');
    dataloader->AddVariable("tro_3_n_muon_segs","tro_3_n_muon_segs","",'F');
    
    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "1>0"; //  
    TCut mycut_b = "(br4_flag==0 || tro_3_flag == 0 || br4_bdt < 0) && (!(truth_nue==1 && truth_CC==1))"; // 8161
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=20000:"
					    "nTrain_Background=7400:"
					    "nTest_Signal=10000:"
					    "nTest_Background=761:"
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
  
  float br4_1_shower_main_length;
  float br4_1_shower_total_length;
  float br4_1_min_dis;
  float br4_1_energy;
  float br4_1_flag_avoid_muon_check;
  float br4_1_n_vtx_segs;
  float br4_1_n_main_segs;
  
  float br4_2_ratio_45;
  float br4_2_ratio_35;
  float br4_2_ratio_25;
  float br4_2_ratio_15;
  float br4_2_ratio1_45;
  float br4_2_ratio1_35;
  float br4_2_ratio1_25;
  float br4_2_ratio1_15;
  float br4_2_iso_angle;
  float br4_2_iso_angle1;
  float br4_2_angle;
  
  float br4_flag;

  float tro_3_stem_length;
  float tro_3_n_muon_segs;
  float tro_3_energy;
  float tro_3_flag;
  
  sig->SetBranchAddress("tro_3_stem_length",&tro_3_stem_length);
  sig->SetBranchAddress("tro_3_n_muon_segs",&tro_3_n_muon_segs);
  sig->SetBranchAddress("tro_3_energy",&tro_3_energy);
  sig->SetBranchAddress("tro_3_flag",&tro_3_flag);

  bkg->SetBranchAddress("tro_3_stem_length",&tro_3_stem_length);
  bkg->SetBranchAddress("tro_3_n_muon_segs",&tro_3_n_muon_segs);
  bkg->SetBranchAddress("tro_3_energy",&tro_3_energy);
  bkg->SetBranchAddress("tro_3_flag",&tro_3_flag);

   sig->SetBranchAddress("br4_1_shower_main_length", &br4_1_shower_main_length);
  sig->SetBranchAddress("br4_1_shower_total_length", &br4_1_shower_total_length);
  sig->SetBranchAddress("br4_1_min_dis", &br4_1_min_dis);
  sig->SetBranchAddress("br4_1_energy", &br4_1_energy);
  sig->SetBranchAddress("br4_1_flag_avoid_muon_check", &br4_1_flag_avoid_muon_check);
  sig->SetBranchAddress("br4_1_n_vtx_segs", &br4_1_n_vtx_segs);
  sig->SetBranchAddress("br4_1_n_main_segs", &br4_1_n_main_segs);
    
  sig->SetBranchAddress("br4_2_ratio_45", &br4_2_ratio_45);
  sig->SetBranchAddress("br4_2_ratio_35", &br4_2_ratio_35);
  sig->SetBranchAddress("br4_2_ratio_25", &br4_2_ratio_25);
  sig->SetBranchAddress("br4_2_ratio_15", &br4_2_ratio_15);
  sig->SetBranchAddress("br4_2_ratio1_45", &br4_2_ratio1_45);
  sig->SetBranchAddress("br4_2_ratio1_35", &br4_2_ratio1_35);
  sig->SetBranchAddress("br4_2_ratio1_25", &br4_2_ratio1_25);
  sig->SetBranchAddress("br4_2_ratio1_15", &br4_2_ratio1_15);
  sig->SetBranchAddress("br4_2_iso_angle", &br4_2_iso_angle);
  sig->SetBranchAddress("br4_2_iso_angle1", &br4_2_iso_angle1);
  sig->SetBranchAddress("br4_2_angle", &br4_2_angle);
  
  sig->SetBranchAddress("br4_flag", &br4_flag);

  bkg->SetBranchAddress("br4_1_shower_main_length", &br4_1_shower_main_length);
  bkg->SetBranchAddress("br4_1_shower_total_length", &br4_1_shower_total_length);
  bkg->SetBranchAddress("br4_1_min_dis", &br4_1_min_dis);
  bkg->SetBranchAddress("br4_1_energy", &br4_1_energy);
  bkg->SetBranchAddress("br4_1_flag_avoid_muon_check", &br4_1_flag_avoid_muon_check);
  bkg->SetBranchAddress("br4_1_n_vtx_segs", &br4_1_n_vtx_segs);
  bkg->SetBranchAddress("br4_1_n_main_segs", &br4_1_n_main_segs);
    
  bkg->SetBranchAddress("br4_2_ratio_45", &br4_2_ratio_45);
  bkg->SetBranchAddress("br4_2_ratio_35", &br4_2_ratio_35);
  bkg->SetBranchAddress("br4_2_ratio_25", &br4_2_ratio_25);
  bkg->SetBranchAddress("br4_2_ratio_15", &br4_2_ratio_15);
  bkg->SetBranchAddress("br4_2_ratio1_45", &br4_2_ratio1_45);
  bkg->SetBranchAddress("br4_2_ratio1_35", &br4_2_ratio1_35);
  bkg->SetBranchAddress("br4_2_ratio1_25", &br4_2_ratio1_25);
  bkg->SetBranchAddress("br4_2_ratio1_15", &br4_2_ratio1_15);
  bkg->SetBranchAddress("br4_2_iso_angle", &br4_2_iso_angle);
  bkg->SetBranchAddress("br4_2_iso_angle1", &br4_2_iso_angle1);
  bkg->SetBranchAddress("br4_2_angle", &br4_2_angle);
  
  bkg->SetBranchAddress("br4_flag", &br4_flag);
   
  
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

  Tsig->Branch("truth_inFV",&truth_inFV,"data/I");
  Tsig->Branch("truth_CC",&truth_CC,"data/I");
  Tsig->Branch("truth_nue",&truth_nue,"data/I");
  Tsig->Branch("truth_cosmic",&truth_cosmic,"data/I");

  Tbkg->Branch("truth_inFV",&truth_inFV,"data/I");
  Tbkg->Branch("truth_CC",&truth_CC,"data/I");
  Tbkg->Branch("truth_nue",&truth_nue,"data/I");
  Tbkg->Branch("truth_cosmic",&truth_cosmic,"data/I");
  
  Tsig->Branch("tro_3_stem_length",&tro_3_stem_length,"tro_3_stem_length/F");
  Tsig->Branch("tro_3_n_muon_segs",&tro_3_n_muon_segs,"tro_3_n_muon_segs/F");
  Tsig->Branch("tro_3_energy",&tro_3_energy,"tro_3_energy/F");
  Tsig->Branch("tro_3_flag",&tro_3_flag,"tro_3_flag/F");

  Tbkg->Branch("tro_3_stem_length",&tro_3_stem_length,"tro_3_stem_length/F");
  Tbkg->Branch("tro_3_n_muon_segs",&tro_3_n_muon_segs,"tro_3_n_muon_segs/F");
  Tbkg->Branch("tro_3_energy",&tro_3_energy,"tro_3_energy/F");
  Tbkg->Branch("tro_3_flag",&tro_3_flag,"tro_3_flag/F");

  
  Tsig->Branch("br4_1_shower_main_length", &br4_1_shower_main_length,"br4_1_shower_main_length/F");
  Tsig->Branch("br4_1_shower_total_length", &br4_1_shower_total_length,"br4_1_shower_total_length/F");
  Tsig->Branch("br4_1_min_dis", &br4_1_min_dis,"br4_1_min_dis/F");
  Tsig->Branch("br4_1_energy", &br4_1_energy,"br4_1_energy/F");
  Tsig->Branch("br4_1_flag_avoid_muon_check", &br4_1_flag_avoid_muon_check,"br4_1_flag_avoid_muon_check/F");
  Tsig->Branch("br4_1_n_vtx_segs", &br4_1_n_vtx_segs,"br4_1_n_vtx_segs/F");
  Tsig->Branch("br4_1_n_main_segs", &br4_1_n_main_segs,"br4_1_n_main_segs/F");
    
  Tsig->Branch("br4_2_ratio_45", &br4_2_ratio_45,"br4_2_ratio_45/F");
  Tsig->Branch("br4_2_ratio_35", &br4_2_ratio_35,"br4_2_ratio_35/F");
  Tsig->Branch("br4_2_ratio_25", &br4_2_ratio_25,"br4_2_ratio_25/F");
  Tsig->Branch("br4_2_ratio_15", &br4_2_ratio_15,"br4_2_ratio_15/F");
  Tsig->Branch("br4_2_ratio1_45", &br4_2_ratio1_45,"br4_2_ratio1_45/F");
  Tsig->Branch("br4_2_ratio1_35", &br4_2_ratio1_35,"br4_2_ratio1_35/F");
  Tsig->Branch("br4_2_ratio1_25", &br4_2_ratio1_25,"br4_2_ratio1_25/F");
  Tsig->Branch("br4_2_ratio1_15", &br4_2_ratio1_15,"br4_2_ratio1_15/F");
  Tsig->Branch("br4_2_iso_angle", &br4_2_iso_angle,"br4_2_iso_angle/F");
  Tsig->Branch("br4_2_iso_angle1", &br4_2_iso_angle1,"br4_2_iso_angle1/F");
  Tsig->Branch("br4_2_angle", &br4_2_angle,"br4_2_angle/F");
    
  Tsig->Branch("br4_flag", &br4_flag,"br4_flag/F");

  
  Tbkg->Branch("br4_1_shower_main_length", &br4_1_shower_main_length,"br4_1_shower_main_length/F");
  Tbkg->Branch("br4_1_shower_total_length", &br4_1_shower_total_length,"br4_1_shower_total_length/F");
  Tbkg->Branch("br4_1_min_dis", &br4_1_min_dis,"br4_1_min_dis/F");
  Tbkg->Branch("br4_1_energy", &br4_1_energy,"br4_1_energy/F");
  Tbkg->Branch("br4_1_flag_avoid_muon_check", &br4_1_flag_avoid_muon_check,"br4_1_flag_avoid_muon_check/F");
  Tbkg->Branch("br4_1_n_vtx_segs", &br4_1_n_vtx_segs,"br4_1_n_vtx_segs/F");
  Tbkg->Branch("br4_1_n_main_segs", &br4_1_n_main_segs,"br4_1_n_main_segs/F");
    
  Tbkg->Branch("br4_2_ratio_45", &br4_2_ratio_45,"br4_2_ratio_45/F");
  Tbkg->Branch("br4_2_ratio_35", &br4_2_ratio_35,"br4_2_ratio_35/F");
  Tbkg->Branch("br4_2_ratio_25", &br4_2_ratio_25,"br4_2_ratio_25/F");
  Tbkg->Branch("br4_2_ratio_15", &br4_2_ratio_15,"br4_2_ratio_15/F");
  Tbkg->Branch("br4_2_ratio1_45", &br4_2_ratio1_45,"br4_2_ratio1_45/F");
  Tbkg->Branch("br4_2_ratio1_35", &br4_2_ratio1_35,"br4_2_ratio1_35/F");
  Tbkg->Branch("br4_2_ratio1_25", &br4_2_ratio1_25,"br4_2_ratio1_25/F");
  Tbkg->Branch("br4_2_ratio1_15", &br4_2_ratio1_15,"br4_2_ratio1_15/F");
  Tbkg->Branch("br4_2_iso_angle", &br4_2_iso_angle,"br4_2_iso_angle/F");
  Tbkg->Branch("br4_2_iso_angle1", &br4_2_iso_angle1,"br4_2_iso_angle1/F");
  Tbkg->Branch("br4_2_angle", &br4_2_angle,"br4_2_angle/F");
    
  Tbkg->Branch("br4_flag", &br4_flag,"br4_flag/F");
  
  
  Tsig->Branch("br4_bdt",&bdt_value,"data/F");
  Tbkg->Branch("br4_bdt",&bdt_value,"data/F");
  
  TMVA::Reader *reader = new TMVA::Reader();
 
  reader->AddVariable("br4_1_shower_main_length",&br4_1_shower_main_length);
  reader->AddVariable("br4_1_shower_total_length",&br4_1_shower_total_length);
  reader->AddVariable("br4_1_min_dis",&br4_1_min_dis);
  reader->AddVariable("br4_1_energy",&br4_1_energy);
  reader->AddVariable("br4_1_flag_avoid_muon_check",&br4_1_flag_avoid_muon_check);
  reader->AddVariable("br4_1_n_vtx_segs",&br4_1_n_vtx_segs);
  reader->AddVariable("br4_1_n_main_segs",&br4_1_n_main_segs);
  
  reader->AddVariable("br4_2_ratio_45",&br4_2_ratio_45);
  reader->AddVariable("br4_2_ratio_35",&br4_2_ratio_35);
  reader->AddVariable("br4_2_ratio_25",&br4_2_ratio_25);
  reader->AddVariable("br4_2_ratio_15",&br4_2_ratio_15);
  reader->AddVariable("br4_2_ratio1_45",&br4_2_ratio1_45);
  reader->AddVariable("br4_2_ratio1_35",&br4_2_ratio1_35);
  reader->AddVariable("br4_2_ratio1_25",&br4_2_ratio1_25);
  reader->AddVariable("br4_2_ratio1_15",&br4_2_ratio1_15);
  reader->AddVariable("br4_2_iso_angle",&br4_2_iso_angle);
  reader->AddVariable("br4_2_iso_angle1",&br4_2_iso_angle1);
  reader->AddVariable("br4_2_angle",&br4_2_angle);

  reader->AddVariable("tro_3_stem_length",&tro_3_stem_length);
  reader->AddVariable("tro_3_n_muon_segs",&tro_3_n_muon_segs);
  
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

