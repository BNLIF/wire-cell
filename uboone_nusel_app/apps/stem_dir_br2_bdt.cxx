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
  
  // stem direction
  int stem_dir_flag;
  int stem_dir_flag_single_shower;
  double stem_dir_angle;
  double stem_dir_energy;
  double stem_dir_angle1;
  double stem_dir_angle2;
  double stem_dir_angle3;
  double stem_dir_ratio;

  int br2_flag;
  int br2_num_valid_tracks;
  int br2_n_shower_main_segs;
  double br2_max_angle;
  double br2_sg_length;
  int br2_flag_sg_trajectory;
  
  sig->SetBranchAddress("stem_dir_flag",&stem_dir_flag);
  sig->SetBranchAddress("stem_dir_flag_single_shower",&stem_dir_flag_single_shower);
  sig->SetBranchAddress("stem_dir_angle",&stem_dir_angle);
  sig->SetBranchAddress("stem_dir_energy",&stem_dir_energy);
  sig->SetBranchAddress("stem_dir_angle1",&stem_dir_angle1);
  sig->SetBranchAddress("stem_dir_angle2",&stem_dir_angle2);
  sig->SetBranchAddress("stem_dir_angle3",&stem_dir_angle3);
  sig->SetBranchAddress("stem_dir_ratio",&stem_dir_ratio);

  sig->SetBranchAddress("br2_flag",&br2_flag);
  sig->SetBranchAddress("br2_num_valid_tracks",&br2_num_valid_tracks);
  sig->SetBranchAddress("br2_n_shower_main_segs",&br2_n_shower_main_segs);
  sig->SetBranchAddress("br2_max_angle",&br2_max_angle);
  sig->SetBranchAddress("br2_sg_length",&br2_sg_length);
  sig->SetBranchAddress("br2_flag_sg_trajectory",&br2_flag_sg_trajectory);

  
  bkg->SetBranchAddress("stem_dir_flag",&stem_dir_flag);
  bkg->SetBranchAddress("stem_dir_flag_single_shower",&stem_dir_flag_single_shower);
  bkg->SetBranchAddress("stem_dir_angle",&stem_dir_angle);
  bkg->SetBranchAddress("stem_dir_energy",&stem_dir_energy);
  bkg->SetBranchAddress("stem_dir_angle1",&stem_dir_angle1);
  bkg->SetBranchAddress("stem_dir_angle2",&stem_dir_angle2);
  bkg->SetBranchAddress("stem_dir_angle3",&stem_dir_angle3);
  bkg->SetBranchAddress("stem_dir_ratio",&stem_dir_ratio);

  bkg->SetBranchAddress("br2_flag",&br2_flag);
  bkg->SetBranchAddress("br2_num_valid_tracks",&br2_num_valid_tracks);
  bkg->SetBranchAddress("br2_n_shower_main_segs",&br2_n_shower_main_segs);
  bkg->SetBranchAddress("br2_max_angle",&br2_max_angle);
  bkg->SetBranchAddress("br2_sg_length",&br2_sg_length);
  bkg->SetBranchAddress("br2_flag_sg_trajectory",&br2_flag_sg_trajectory);
  
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
  
  // stem direction
  float stem_dir_flag_f;
  float stem_dir_flag_single_shower_f;
  float stem_dir_angle_f;
  float stem_dir_energy_f;
  float stem_dir_angle1_f;
  float stem_dir_angle2_f;
  float stem_dir_angle3_f;
  float stem_dir_ratio_f;

  float br2_flag_f;
  float br2_num_valid_tracks_f;
  float br2_n_shower_main_segs_f;
  float br2_max_angle_f;
  float br2_sg_length_f;
  float br2_flag_sg_trajectory_f;

  Tsig->Branch("stem_dir_flag",&stem_dir_flag_f,"stem_dir_flag/F");
  Tsig->Branch("stem_dir_flag_single_shower",&stem_dir_flag_single_shower_f,"stem_dir_flag_single_shower/F");
  Tsig->Branch("stem_dir_angle",&stem_dir_angle_f,"stem_dir_angle/F");
  Tsig->Branch("stem_dir_energy",&stem_dir_energy_f,"stem_dir_energy/F");
  Tsig->Branch("stem_dir_angle1",&stem_dir_angle1_f,"stem_dir_angle1/F");
  Tsig->Branch("stem_dir_angle2",&stem_dir_angle2_f,"stem_dir_angle2/F");
  Tsig->Branch("stem_dir_angle3",&stem_dir_angle3_f,"stem_dir_angle3/F");
  Tsig->Branch("stem_dir_ratio",&stem_dir_ratio_f,"stem_dir_ratio/F");
  
  Tsig->Branch("br2_flag",&br2_flag_f,"br2_flag/F");
  Tsig->Branch("br2_num_valid_tracks",&br2_num_valid_tracks_f,"br2_num_valid_tracks/F");
  Tsig->Branch("br2_n_shower_main_segs",&br2_n_shower_main_segs_f,"br2_n_shower_main_segs/F");
  Tsig->Branch("br2_max_angle",&br2_max_angle_f,"br2_max_angle/F");
  Tsig->Branch("br2_sg_length",&br2_sg_length_f,"br2_sg_length/F");
  Tsig->Branch("br2_flag_sg_trajectory",&br2_flag_sg_trajectory_f,"br2_flag_sg_trajectory/F");

  Tbkg->Branch("stem_dir_flag",&stem_dir_flag_f,"stem_dir_flag/F");
  Tbkg->Branch("stem_dir_flag_single_shower",&stem_dir_flag_single_shower_f,"stem_dir_flag_single_shower/F");
  Tbkg->Branch("stem_dir_angle",&stem_dir_angle_f,"stem_dir_angle/F");
  Tbkg->Branch("stem_dir_energy",&stem_dir_energy_f,"stem_dir_energy/F");
  Tbkg->Branch("stem_dir_angle1",&stem_dir_angle1_f,"stem_dir_angle1/F");
  Tbkg->Branch("stem_dir_angle2",&stem_dir_angle2_f,"stem_dir_angle2/F");
  Tbkg->Branch("stem_dir_angle3",&stem_dir_angle3_f,"stem_dir_angle3/F");
  Tbkg->Branch("stem_dir_ratio",&stem_dir_ratio_f,"stem_dir_ratio/F");
  
  Tbkg->Branch("br2_flag",&br2_flag_f,"br2_flag/F");
  Tbkg->Branch("br2_num_valid_tracks",&br2_num_valid_tracks_f,"br2_num_valid_tracks/F");
  Tbkg->Branch("br2_n_shower_main_segs",&br2_n_shower_main_segs_f,"br2_n_shower_main_segs/F");
  Tbkg->Branch("br2_max_angle",&br2_max_angle_f,"br2_max_angle/F");
  Tbkg->Branch("br2_sg_length",&br2_sg_length_f,"br2_sg_length/F");
  Tbkg->Branch("br2_flag_sg_trajectory",&br2_flag_sg_trajectory_f,"br2_flag_sg_trajectory/F");
  
  
  for (Int_t i=0;i!=sig->GetEntries();i++){
    sig->GetEntry(i);
    
    stem_dir_flag_f = stem_dir_flag;
    stem_dir_flag_single_shower_f = stem_dir_flag_single_shower;
    stem_dir_angle_f = stem_dir_angle;
    stem_dir_energy_f = stem_dir_energy;
    stem_dir_angle1_f = stem_dir_angle1;
    stem_dir_angle2_f = stem_dir_angle2;
    stem_dir_angle3_f = stem_dir_angle3 ;
    stem_dir_ratio_f = stem_dir_ratio;

    br2_flag_f = br2_flag;
    br2_num_valid_tracks_f = br2_num_valid_tracks;
    br2_n_shower_main_segs_f = br2_n_shower_main_segs;
    br2_max_angle_f = br2_max_angle;
    br2_sg_length_f = br2_sg_length;
    br2_flag_sg_trajectory_f = br2_flag_sg_trajectory;
    
    Tsig->Fill();
  }

   for (Int_t i=0;i!=bkg->GetEntries();i++){
    bkg->GetEntry(i);
    
    stem_dir_flag_f = stem_dir_flag;
    stem_dir_flag_single_shower_f = stem_dir_flag_single_shower;
    stem_dir_angle_f = stem_dir_angle;
    stem_dir_energy_f = stem_dir_energy;
    stem_dir_angle1_f = stem_dir_angle1;
    stem_dir_angle2_f = stem_dir_angle2;
    stem_dir_angle3_f = stem_dir_angle3 ;
    stem_dir_ratio_f = stem_dir_ratio;

    br2_flag_f = br2_flag;
    br2_num_valid_tracks_f = br2_num_valid_tracks;
    br2_n_shower_main_segs_f = br2_n_shower_main_segs;
    br2_max_angle_f = br2_max_angle;
    br2_sg_length_f = br2_sg_length;
    br2_flag_sg_trajectory_f = br2_flag_sg_trajectory;
    
    
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
   
    
    dataloader->AddVariable("stem_dir_flag_single_shower","stem_dir_flag_single_shower","",'F');
    dataloader->AddVariable("stem_dir_angle","stem_dir_angle","deg",'F');
    dataloader->AddVariable("stem_dir_energy","stem_dir_energy","MeV",'F');
    dataloader->AddVariable("stem_dir_angle1","stem_dir_angle1","deg",'F');
    dataloader->AddVariable("stem_dir_angle2","stem_dir_angle2","deg",'F');
    dataloader->AddVariable("stem_dir_angle3","stem_dir_angle3","deg",'F');
    dataloader->AddVariable("stem_dir_ratio","stem_dir_ratio","",'F');

    dataloader->AddVariable("br2_num_valid_tracks","br2_num_valid_tracks","",'F');
    dataloader->AddVariable("br2_n_shower_main_segs","br2_n_shower_main_segs","",'F');
    dataloader->AddVariable("br2_max_angle","br2_max_angle","deg",'F');
    dataloader->AddVariable("br2_sg_length","br2_sg_length","cm",'F');
    dataloader->AddVariable("br2_flag_sg_trajectory","br2_flag_sg_trajectory","",'F');
    
    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "1>0"; //  482 (br2) 1487 ( stem_dir)  1487 (both)/44194
    TCut mycut_b = "(stem_dir_flag==0 || br2_flag == 0) && (!(truth_nue==1 && truth_CC==1))"; //  3062 (br2) 6558 (stem_direction) 6660 (both)
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=20000:"
					    "nTrain_Background=5500:"
					    "nTest_Signal=10000:"
					    "nTest_Background=1068:"
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
   
    
    dataloader->AddVariable("stem_dir_flag_single_shower","stem_dir_flag_single_shower","",'F');
    dataloader->AddVariable("stem_dir_angle","stem_dir_angle","deg",'F');
    dataloader->AddVariable("stem_dir_energy","stem_dir_energy","MeV",'F');
    dataloader->AddVariable("stem_dir_angle1","stem_dir_angle1","deg",'F');
    dataloader->AddVariable("stem_dir_angle2","stem_dir_angle2","deg",'F');
    dataloader->AddVariable("stem_dir_angle3","stem_dir_angle3","deg",'F');
    dataloader->AddVariable("stem_dir_ratio","stem_dir_ratio","",'F');

    dataloader->AddVariable("br2_num_valid_tracks","br2_num_valid_tracks","",'F');
    dataloader->AddVariable("br2_n_shower_main_segs","br2_n_shower_main_segs","",'F');
    dataloader->AddVariable("br2_max_angle","br2_max_angle","deg",'F');
    dataloader->AddVariable("br2_sg_length","br2_sg_length","cm",'F');
    dataloader->AddVariable("br2_flag_sg_trajectory","br2_flag_sg_trajectory","",'F');
    
    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "1>0"; //  
    TCut mycut_b = "(stem_dir_flag==0 || br2_flag == 0 || stem_dir_br2_bdt < -0.06) && (!(truth_nue==1 && truth_CC==1))"; //  7573
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=20000:"
					    "nTrain_Background=6400:"
					    "nTest_Signal=10000:"
					    "nTest_Background=1001:"
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

  float stem_dir_flag;
  float stem_dir_flag_single_shower;
  float stem_dir_angle;
  float stem_dir_energy;
  float stem_dir_angle1;
  float stem_dir_angle2;
  float stem_dir_angle3;
  float stem_dir_ratio;
  
  float br2_flag;
  float br2_num_valid_tracks;
  float br2_n_shower_main_segs;
  float br2_max_angle;
  float br2_sg_length;
  float br2_flag_sg_trajectory;

  
   sig->SetBranchAddress("stem_dir_flag",&stem_dir_flag);
  sig->SetBranchAddress("stem_dir_flag_single_shower",&stem_dir_flag_single_shower);
  sig->SetBranchAddress("stem_dir_angle",&stem_dir_angle);
  sig->SetBranchAddress("stem_dir_energy",&stem_dir_energy);
  sig->SetBranchAddress("stem_dir_angle1",&stem_dir_angle1);
  sig->SetBranchAddress("stem_dir_angle2",&stem_dir_angle2);
  sig->SetBranchAddress("stem_dir_angle3",&stem_dir_angle3);
  sig->SetBranchAddress("stem_dir_ratio",&stem_dir_ratio);

  sig->SetBranchAddress("br2_flag",&br2_flag);
  sig->SetBranchAddress("br2_num_valid_tracks",&br2_num_valid_tracks);
  sig->SetBranchAddress("br2_n_shower_main_segs",&br2_n_shower_main_segs);
  sig->SetBranchAddress("br2_max_angle",&br2_max_angle);
  sig->SetBranchAddress("br2_sg_length",&br2_sg_length);
  sig->SetBranchAddress("br2_flag_sg_trajectory",&br2_flag_sg_trajectory);

  
  bkg->SetBranchAddress("stem_dir_flag",&stem_dir_flag);
  bkg->SetBranchAddress("stem_dir_flag_single_shower",&stem_dir_flag_single_shower);
  bkg->SetBranchAddress("stem_dir_angle",&stem_dir_angle);
  bkg->SetBranchAddress("stem_dir_energy",&stem_dir_energy);
  bkg->SetBranchAddress("stem_dir_angle1",&stem_dir_angle1);
  bkg->SetBranchAddress("stem_dir_angle2",&stem_dir_angle2);
  bkg->SetBranchAddress("stem_dir_angle3",&stem_dir_angle3);
  bkg->SetBranchAddress("stem_dir_ratio",&stem_dir_ratio);

  bkg->SetBranchAddress("br2_flag",&br2_flag);
  bkg->SetBranchAddress("br2_num_valid_tracks",&br2_num_valid_tracks);
  bkg->SetBranchAddress("br2_n_shower_main_segs",&br2_n_shower_main_segs);
  bkg->SetBranchAddress("br2_max_angle",&br2_max_angle);
  bkg->SetBranchAddress("br2_sg_length",&br2_sg_length);
  bkg->SetBranchAddress("br2_flag_sg_trajectory",&br2_flag_sg_trajectory);
   
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
  
  Tsig->Branch("stem_dir_flag",&stem_dir_flag,"stem_dir_flag/F");
  Tsig->Branch("stem_dir_flag_single_shower",&stem_dir_flag_single_shower,"stem_dir_flag_single_shower/F");
  Tsig->Branch("stem_dir_angle",&stem_dir_angle,"stem_dir_angle/F");
  Tsig->Branch("stem_dir_energy",&stem_dir_energy,"stem_dir_energy/F");
  Tsig->Branch("stem_dir_angle1",&stem_dir_angle1,"stem_dir_angle1/F");
  Tsig->Branch("stem_dir_angle2",&stem_dir_angle2,"stem_dir_angle2/F");
  Tsig->Branch("stem_dir_angle3",&stem_dir_angle3,"stem_dir_angle3/F");
  Tsig->Branch("stem_dir_ratio",&stem_dir_ratio,"stem_dir_ratio/F");
  
  Tsig->Branch("br2_flag",&br2_flag,"br2_flag/F");
  Tsig->Branch("br2_num_valid_tracks",&br2_num_valid_tracks,"br2_num_valid_tracks/F");
  Tsig->Branch("br2_n_shower_main_segs",&br2_n_shower_main_segs,"br2_n_shower_main_segs/F");
  Tsig->Branch("br2_max_angle",&br2_max_angle,"br2_max_angle/F");
  Tsig->Branch("br2_sg_length",&br2_sg_length,"br2_sg_length/F");
  Tsig->Branch("br2_flag_sg_trajectory",&br2_flag_sg_trajectory,"br2_flag_sg_trajectory/F");

  Tbkg->Branch("stem_dir_flag",&stem_dir_flag,"stem_dir_flag/F");
  Tbkg->Branch("stem_dir_flag_single_shower",&stem_dir_flag_single_shower,"stem_dir_flag_single_shower/F");
  Tbkg->Branch("stem_dir_angle",&stem_dir_angle,"stem_dir_angle/F");
  Tbkg->Branch("stem_dir_energy",&stem_dir_energy,"stem_dir_energy/F");
  Tbkg->Branch("stem_dir_angle1",&stem_dir_angle1,"stem_dir_angle1/F");
  Tbkg->Branch("stem_dir_angle2",&stem_dir_angle2,"stem_dir_angle2/F");
  Tbkg->Branch("stem_dir_angle3",&stem_dir_angle3,"stem_dir_angle3/F");
  Tbkg->Branch("stem_dir_ratio",&stem_dir_ratio,"stem_dir_ratio/F");
  
  Tbkg->Branch("br2_flag",&br2_flag,"br2_flag/F");
  Tbkg->Branch("br2_num_valid_tracks",&br2_num_valid_tracks,"br2_num_valid_tracks/F");
  Tbkg->Branch("br2_n_shower_main_segs",&br2_n_shower_main_segs,"br2_n_shower_main_segs/F");
  Tbkg->Branch("br2_max_angle",&br2_max_angle,"br2_max_angle/F");
  Tbkg->Branch("br2_sg_length",&br2_sg_length,"br2_sg_length/F");
  Tbkg->Branch("br2_flag_sg_trajectory",&br2_flag_sg_trajectory,"br2_flag_sg_trajectory/F");
  
  
  Tsig->Branch("stem_dir_br2_bdt",&bdt_value,"data/F");
  Tbkg->Branch("stem_dir_br2_bdt",&bdt_value,"data/F");
  
  TMVA::Reader *reader = new TMVA::Reader();
  reader->AddVariable("stem_dir_flag_single_shower",&stem_dir_flag_single_shower);
  reader->AddVariable("stem_dir_angle",&stem_dir_angle);
  reader->AddVariable("stem_dir_energy",&stem_dir_energy);
  reader->AddVariable("stem_dir_angle1",&stem_dir_angle1);
  reader->AddVariable("stem_dir_angle2",&stem_dir_angle2);
  reader->AddVariable("stem_dir_angle3",&stem_dir_angle3);
  reader->AddVariable("stem_dir_ratio",&stem_dir_ratio);
  
  reader->AddVariable("br2_num_valid_tracks",&br2_num_valid_tracks);
  reader->AddVariable("br2_n_shower_main_segs",&br2_n_shower_main_segs);
  reader->AddVariable("br2_max_angle",&br2_max_angle);
  reader->AddVariable("br2_sg_length",&br2_sg_length);
  reader->AddVariable("br2_flag_sg_trajectory",&br2_flag_sg_trajectory);
  
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

