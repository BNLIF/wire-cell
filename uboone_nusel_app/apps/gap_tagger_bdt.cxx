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

  int gap_flag;
  int gap_flag_prolong_u;
  int gap_flag_prolong_v;
  int gap_flag_prolong_w;
  int gap_flag_parallel;
  int gap_n_points;
  int gap_n_bad;
  double gap_energy;
  int gap_num_valid_tracks;
  int gap_flag_single_shower;
  int gap_filled;
  
  sig->SetBranchAddress("gap_flag",&gap_flag);
  sig->SetBranchAddress("gap_flag_prolong_u",&gap_flag_prolong_u);
  sig->SetBranchAddress("gap_flag_prolong_v",&gap_flag_prolong_v);
  sig->SetBranchAddress("gap_flag_prolong_w",&gap_flag_prolong_w);
  sig->SetBranchAddress("gap_flag_parallel",&gap_flag_parallel);
  sig->SetBranchAddress("gap_n_points",&gap_n_points);
  sig->SetBranchAddress("gap_n_bad",&gap_n_bad);
  sig->SetBranchAddress("gap_energy",&gap_energy);
  sig->SetBranchAddress("gap_num_valid_tracks",&gap_num_valid_tracks);
  sig->SetBranchAddress("gap_flag_single_shower",&gap_flag_single_shower);
  sig->SetBranchAddress("gap_filled",&gap_filled);

  bkg->SetBranchAddress("gap_flag",&gap_flag);
  bkg->SetBranchAddress("gap_flag_prolong_u",&gap_flag_prolong_u);
  bkg->SetBranchAddress("gap_flag_prolong_v",&gap_flag_prolong_v);
  bkg->SetBranchAddress("gap_flag_prolong_w",&gap_flag_prolong_w);
  bkg->SetBranchAddress("gap_flag_parallel",&gap_flag_parallel);
  bkg->SetBranchAddress("gap_n_points",&gap_n_points);
  bkg->SetBranchAddress("gap_n_bad",&gap_n_bad);
  bkg->SetBranchAddress("gap_energy",&gap_energy);
  bkg->SetBranchAddress("gap_num_valid_tracks",&gap_num_valid_tracks);
  bkg->SetBranchAddress("gap_flag_single_shower",&gap_flag_single_shower);
  bkg->SetBranchAddress("gap_filled",&gap_filled);

  
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

  
  float gap_flag_f;
  float gap_flag_prolong_u_f;
  float gap_flag_prolong_v_f;
  float gap_flag_prolong_w_f;
  float gap_flag_parallel_f;
  float gap_n_points_f;
  float gap_n_bad_f;
  float gap_energy_f;
  float gap_num_valid_tracks_f;
  float gap_flag_single_shower_f;
  float gap_filled_f;

  Tsig->Branch("gap_flag",&gap_flag_f,"gap_flag/F");
  Tsig->Branch("gap_flag_prolong_u",&gap_flag_prolong_u_f,"gap_flag_prolong_u/F");
  Tsig->Branch("gap_flag_prolong_v",&gap_flag_prolong_v_f,"gap_flag_prolong_v/F");
  Tsig->Branch("gap_flag_prolong_w",&gap_flag_prolong_w_f,"gap_flag_prolong_w/F");
  Tsig->Branch("gap_flag_parallel",&gap_flag_parallel_f,"gap_flag_parallel/F");
  Tsig->Branch("gap_n_points",&gap_n_points_f,"gap_n_points/F");
  Tsig->Branch("gap_n_bad",&gap_n_bad_f,"gap_n_bad/F");
  Tsig->Branch("gap_energy",&gap_energy_f,"gap_energy/F");
  Tsig->Branch("gap_num_valid_tracks",&gap_num_valid_tracks_f,"gap_num_valid_tracks/F");
  Tsig->Branch("gap_flag_single_shower",&gap_flag_single_shower_f,"gap_flag_single_shower/F");
  Tsig->Branch("gap_filled",&gap_filled_f,"gap_filled/F");
  
  Tbkg->Branch("gap_flag",&gap_flag_f,"gap_flag/F");
  Tbkg->Branch("gap_flag_prolong_u",&gap_flag_prolong_u_f,"gap_flag_prolong_u/F");
  Tbkg->Branch("gap_flag_prolong_v",&gap_flag_prolong_v_f,"gap_flag_prolong_v/F");
  Tbkg->Branch("gap_flag_prolong_w",&gap_flag_prolong_w_f,"gap_flag_prolong_w/F");
  Tbkg->Branch("gap_flag_parallel",&gap_flag_parallel_f,"gap_flag_parallel/F");
  Tbkg->Branch("gap_n_points",&gap_n_points_f,"gap_n_points/F");
  Tbkg->Branch("gap_n_bad",&gap_n_bad_f,"gap_n_bad/F");
  Tbkg->Branch("gap_energy",&gap_energy_f,"gap_energy/F");
  Tbkg->Branch("gap_num_valid_tracks",&gap_num_valid_tracks_f,"gap_num_valid_tracks/F");
  Tbkg->Branch("gap_flag_single_shower",&gap_flag_single_shower_f,"gap_flag_single_shower/F");
  Tbkg->Branch("gap_filled",&gap_filled_f,"gap_filled/F");
  
  for (Int_t i=0;i!=sig->GetEntries();i++){
    sig->GetEntry(i);

    gap_flag_f = gap_flag;
    gap_flag_prolong_u_f = gap_flag_prolong_u;
    gap_flag_prolong_v_f = gap_flag_prolong_v;
    gap_flag_prolong_w_f = gap_flag_prolong_w;
    gap_flag_parallel_f = gap_flag_parallel;
    gap_n_points_f = gap_n_points;
    gap_n_bad_f = gap_n_bad;
    gap_energy_f = gap_energy;
    gap_num_valid_tracks_f = gap_num_valid_tracks;
    gap_flag_single_shower_f = gap_flag_single_shower;
    gap_filled_f = gap_filled;
    
    Tsig->Fill();
  }

   for (Int_t i=0;i!=bkg->GetEntries();i++){
    bkg->GetEntry(i);

    gap_flag_f = gap_flag;
    gap_flag_prolong_u_f = gap_flag_prolong_u;
    gap_flag_prolong_v_f = gap_flag_prolong_v;
    gap_flag_prolong_w_f = gap_flag_prolong_w;
    gap_flag_parallel_f = gap_flag_parallel;
    gap_n_points_f = gap_n_points;
    gap_n_bad_f = gap_n_bad;
    gap_energy_f = gap_energy;
    gap_num_valid_tracks_f = gap_num_valid_tracks;
    gap_flag_single_shower_f = gap_flag_single_shower;
    gap_filled_f = gap_filled;
    
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
   
    
     dataloader->AddVariable("gap_flag_prolong_u","gap_flag_prolong_u","",'F');
     dataloader->AddVariable("gap_flag_prolong_v","gap_flag_prolong_v","",'F');
     dataloader->AddVariable("gap_flag_prolong_w","gap_flag_prolong_w","",'F');
     dataloader->AddVariable("gap_flag_parallel","gap_flag_parallel","",'F');
     dataloader->AddVariable("gap_n_points","gap_n_points","",'F');
     dataloader->AddVariable("gap_n_bad","gap_n_bad","",'F');
     dataloader->AddVariable("gap_energy","gap_energy","",'F');
     dataloader->AddVariable("gap_num_valid_tracks","gap_num_valid_tracks","",'F');
     dataloader->AddVariable("gap_flag_single_shower","gap_flag_single_shower","",'F');
     
    

    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "gap_filled==1"; // 962/44194
    TCut mycut_b = "gap_filled==1 && gap_flag==0"; // 2770/21070
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=20000:"
					    "nTrain_Background=2400:"
					    "nTest_Signal=10000:"
					    "nTest_Background=370:"
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
   
    
     dataloader->AddVariable("gap_flag_prolong_u","gap_flag_prolong_u","",'F');
     dataloader->AddVariable("gap_flag_prolong_v","gap_flag_prolong_v","",'F');
     dataloader->AddVariable("gap_flag_prolong_w","gap_flag_prolong_w","",'F');
     dataloader->AddVariable("gap_flag_parallel","gap_flag_parallel","",'F');
     dataloader->AddVariable("gap_n_points","gap_n_points","",'F');
     dataloader->AddVariable("gap_n_bad","gap_n_bad","",'F');
     dataloader->AddVariable("gap_energy","gap_energy","",'F');
     dataloader->AddVariable("gap_num_valid_tracks","gap_num_valid_tracks","",'F');
     dataloader->AddVariable("gap_flag_single_shower","gap_flag_single_shower","",'F');
     
    

    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "gap_filled==1"; // 
    TCut mycut_b = "gap_filled==1 && (gap_flag==0 || gap_bdt < 0)"; // 2862
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=20000:"
					    "nTrain_Background=2500:"
					    "nTest_Signal=10000:"
					    "nTest_Background=362:"
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

  
  float gap_flag;
  float gap_flag_prolong_u;
  float gap_flag_prolong_v;
  float gap_flag_prolong_w;
  float gap_flag_parallel;
  float gap_n_points;
  float gap_n_bad;
  float gap_energy;
  float gap_num_valid_tracks;
  float gap_flag_single_shower;
  float gap_filled;


  sig->SetBranchAddress("gap_flag",&gap_flag);
  sig->SetBranchAddress("gap_flag_prolong_u",&gap_flag_prolong_u);
  sig->SetBranchAddress("gap_flag_prolong_v",&gap_flag_prolong_v);
  sig->SetBranchAddress("gap_flag_prolong_w",&gap_flag_prolong_w);
  sig->SetBranchAddress("gap_flag_parallel",&gap_flag_parallel);
  sig->SetBranchAddress("gap_n_points",&gap_n_points);
  sig->SetBranchAddress("gap_n_bad",&gap_n_bad);
  sig->SetBranchAddress("gap_energy",&gap_energy);
  sig->SetBranchAddress("gap_num_valid_tracks",&gap_num_valid_tracks);
  sig->SetBranchAddress("gap_flag_single_shower",&gap_flag_single_shower);
  sig->SetBranchAddress("gap_filled",&gap_filled);

  bkg->SetBranchAddress("gap_flag",&gap_flag);
  bkg->SetBranchAddress("gap_flag_prolong_u",&gap_flag_prolong_u);
  bkg->SetBranchAddress("gap_flag_prolong_v",&gap_flag_prolong_v);
  bkg->SetBranchAddress("gap_flag_prolong_w",&gap_flag_prolong_w);
  bkg->SetBranchAddress("gap_flag_parallel",&gap_flag_parallel);
  bkg->SetBranchAddress("gap_n_points",&gap_n_points);
  bkg->SetBranchAddress("gap_n_bad",&gap_n_bad);
  bkg->SetBranchAddress("gap_energy",&gap_energy);
  bkg->SetBranchAddress("gap_num_valid_tracks",&gap_num_valid_tracks);
  bkg->SetBranchAddress("gap_flag_single_shower",&gap_flag_single_shower);
  bkg->SetBranchAddress("gap_filled",&gap_filled);

  
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
  
  Tsig->Branch("gap_flag",&gap_flag,"gap_flag/F");
  Tsig->Branch("gap_flag_prolong_u",&gap_flag_prolong_u,"gap_flag_prolong_u/F");
  Tsig->Branch("gap_flag_prolong_v",&gap_flag_prolong_v,"gap_flag_prolong_v/F");
  Tsig->Branch("gap_flag_prolong_w",&gap_flag_prolong_w,"gap_flag_prolong_w/F");
  Tsig->Branch("gap_flag_parallel",&gap_flag_parallel,"gap_flag_parallel/F");
  Tsig->Branch("gap_n_points",&gap_n_points,"gap_n_points/F");
  Tsig->Branch("gap_n_bad",&gap_n_bad,"gap_n_bad/F");
  Tsig->Branch("gap_energy",&gap_energy,"gap_energy/F");
  Tsig->Branch("gap_num_valid_tracks",&gap_num_valid_tracks,"gap_num_valid_tracks/F");
  Tsig->Branch("gap_flag_single_shower",&gap_flag_single_shower,"gap_flag_single_shower/F");
  Tsig->Branch("gap_filled",&gap_filled,"gap_filled/F");
  
  Tbkg->Branch("gap_flag",&gap_flag,"gap_flag/F");
  Tbkg->Branch("gap_flag_prolong_u",&gap_flag_prolong_u,"gap_flag_prolong_u/F");
  Tbkg->Branch("gap_flag_prolong_v",&gap_flag_prolong_v,"gap_flag_prolong_v/F");
  Tbkg->Branch("gap_flag_prolong_w",&gap_flag_prolong_w,"gap_flag_prolong_w/F");
  Tbkg->Branch("gap_flag_parallel",&gap_flag_parallel,"gap_flag_parallel/F");
  Tbkg->Branch("gap_n_points",&gap_n_points,"gap_n_points/F");
  Tbkg->Branch("gap_n_bad",&gap_n_bad,"gap_n_bad/F");
  Tbkg->Branch("gap_energy",&gap_energy,"gap_energy/F");
  Tbkg->Branch("gap_num_valid_tracks",&gap_num_valid_tracks,"gap_num_valid_tracks/F");
  Tbkg->Branch("gap_flag_single_shower",&gap_flag_single_shower,"gap_flag_single_shower/F");
  Tbkg->Branch("gap_filled",&gap_filled,"gap_filled/F");

  Tsig->Branch("gap_bdt",&bdt_value,"data/F");
  Tbkg->Branch("gap_bdt",&bdt_value,"data/F");
  
  TMVA::Reader *reader = new TMVA::Reader();

  reader->AddVariable("gap_flag_prolong_u",&gap_flag_prolong_u);
  reader->AddVariable("gap_flag_prolong_v",&gap_flag_prolong_v);
  reader->AddVariable("gap_flag_prolong_w",&gap_flag_prolong_w);
  reader->AddVariable("gap_flag_parallel",&gap_flag_parallel);
  reader->AddVariable("gap_n_points",&gap_n_points);
  reader->AddVariable("gap_n_bad",&gap_n_bad);
  reader->AddVariable("gap_energy",&gap_energy);
  reader->AddVariable("gap_num_valid_tracks",&gap_num_valid_tracks);
  reader->AddVariable("gap_flag_single_shower",&gap_flag_single_shower);
  
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

