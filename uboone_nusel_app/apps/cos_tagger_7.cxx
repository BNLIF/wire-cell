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
 TFile *file0 = new TFile("bdt.root"); // sig x1
 
 TTree *sig = (TTree*)file0->Get("sig");
 TTree *bkg = (TTree*)file0->Get("bkg");


 float cosmict_flag_1;
 
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

 sig->SetBranchAddress("cosmict_flag_1",&cosmict_flag_1);

 sig->SetBranchAddress("cosmict_flag_7",&cosmict_flag_7);
 sig->SetBranchAddress("cosmict_7_filled",&cosmict_7_filled);
 sig->SetBranchAddress("cosmict_7_flag_sec",&cosmict_7_flag_sec);
 sig->SetBranchAddress("cosmict_7_n_muon_tracks",&cosmict_7_n_muon_tracks);
 sig->SetBranchAddress("cosmict_7_total_shower_length",&cosmict_7_total_shower_length);
 sig->SetBranchAddress("cosmict_7_flag_inside",&cosmict_7_flag_inside);
 sig->SetBranchAddress("cosmict_7_angle_beam",&cosmict_7_angle_beam);
 sig->SetBranchAddress("cosmict_7_flag_dir_weak",&cosmict_7_flag_dir_weak);
 sig->SetBranchAddress("cosmict_7_dQ_dx_end",&cosmict_7_dQ_dx_end);
 sig->SetBranchAddress("cosmict_7_dQ_dx_front",&cosmict_7_dQ_dx_front);
 sig->SetBranchAddress("cosmict_7_theta",&cosmict_7_theta);
 sig->SetBranchAddress("cosmict_7_phi",&cosmict_7_phi);



 bkg->SetBranchAddress("cosmict_flag_1",&cosmict_flag_1);

 bkg->SetBranchAddress("cosmict_flag_7",&cosmict_flag_7);
 bkg->SetBranchAddress("cosmict_7_filled",&cosmict_7_filled);
 bkg->SetBranchAddress("cosmict_7_flag_sec",&cosmict_7_flag_sec);
 bkg->SetBranchAddress("cosmict_7_n_muon_tracks",&cosmict_7_n_muon_tracks);
 bkg->SetBranchAddress("cosmict_7_total_shower_length",&cosmict_7_total_shower_length);
 bkg->SetBranchAddress("cosmict_7_flag_inside",&cosmict_7_flag_inside);
 bkg->SetBranchAddress("cosmict_7_angle_beam",&cosmict_7_angle_beam);
 bkg->SetBranchAddress("cosmict_7_flag_dir_weak",&cosmict_7_flag_dir_weak);
 bkg->SetBranchAddress("cosmict_7_dQ_dx_end",&cosmict_7_dQ_dx_end);
 bkg->SetBranchAddress("cosmict_7_dQ_dx_front",&cosmict_7_dQ_dx_front);
 bkg->SetBranchAddress("cosmict_7_theta",&cosmict_7_theta);
 bkg->SetBranchAddress("cosmict_7_phi",&cosmict_7_phi);

 
 

 float weight;
 sig->SetBranchAddress("weight",&weight);
 bkg->SetBranchAddress("weight",&weight);
 
 TFile *new_file = new TFile("round_0.root","RECREATE");
 TTree *Stree = new TTree("TreeS","signal tree");
 TTree *Btree = new TTree("TreeB","background tree");
 Stree->SetDirectory(new_file);
 Btree->SetDirectory(new_file);



 Stree->Branch("weight",&weight,"weight/F");
 Btree->Branch("weight",&weight,"weight/F");
 
 Stree->Branch("cosmict_flag_1",&cosmict_flag_1,"cosmict_flag_1/F");
 Btree->Branch("cosmict_flag_1",&cosmict_flag_1,"cosmict_flag_1/F");

 Stree->Branch("cosmict_flag_7",&cosmict_flag_7,"cosmict_flag_7/F");
 Stree->Branch("cosmict_7_filled",&cosmict_7_filled,"cosmict_7_filled/F");
 Stree->Branch("cosmict_7_flag_sec",&cosmict_7_flag_sec,"cosmict_7_flag_sec/F");
 Stree->Branch("cosmict_7_n_muon_tracks",&cosmict_7_n_muon_tracks,"cosmict_7_n_muon_tracks/F");
 Stree->Branch("cosmict_7_total_shower_length",&cosmict_7_total_shower_length,"cosmict_7_total_shower_length/F");
 Stree->Branch("cosmict_7_flag_inside",&cosmict_7_flag_inside,"cosmict_7_flag_inside/F");
 Stree->Branch("cosmict_7_angle_beam",&cosmict_7_angle_beam,"cosmict_7_angle_beam/F");
 Stree->Branch("cosmict_7_flag_dir_weak",&cosmict_7_flag_dir_weak,"cosmict_7_flag_dir_weak/F");
 Stree->Branch("cosmict_7_dQ_dx_end",&cosmict_7_dQ_dx_end,"cosmict_7_dQ_dx_end/F");
 Stree->Branch("cosmict_7_dQ_dx_front",&cosmict_7_dQ_dx_front,"cosmict_7_dQ_dx_front/F");
 Stree->Branch("cosmict_7_theta",&cosmict_7_theta,"cosmict_7_theta/F");
 Stree->Branch("cosmict_7_phi",&cosmict_7_phi,"cosmict_7_phi/F");

 
 Btree->Branch("cosmict_flag_7",&cosmict_flag_7,"cosmict_flag_7/F");
 Btree->Branch("cosmict_7_filled",&cosmict_7_filled,"cosmict_7_filled/F");
 Btree->Branch("cosmict_7_flag_sec",&cosmict_7_flag_sec,"cosmict_7_flag_sec/F");
 Btree->Branch("cosmict_7_n_muon_tracks",&cosmict_7_n_muon_tracks,"cosmict_7_n_muon_tracks/F");
 Btree->Branch("cosmict_7_total_shower_length",&cosmict_7_total_shower_length,"cosmict_7_total_shower_length/F");
 Btree->Branch("cosmict_7_flag_inside",&cosmict_7_flag_inside,"cosmict_7_flag_inside/F");
 Btree->Branch("cosmict_7_angle_beam",&cosmict_7_angle_beam,"cosmict_7_angle_beam/F");
 Btree->Branch("cosmict_7_flag_dir_weak",&cosmict_7_flag_dir_weak,"cosmict_7_flag_dir_weak/F");
 Btree->Branch("cosmict_7_dQ_dx_end",&cosmict_7_dQ_dx_end,"cosmict_7_dQ_dx_end/F");
 Btree->Branch("cosmict_7_dQ_dx_front",&cosmict_7_dQ_dx_front,"cosmict_7_dQ_dx_front/F");
 Btree->Branch("cosmict_7_theta",&cosmict_7_theta,"cosmict_7_theta/F");
 Btree->Branch("cosmict_7_phi",&cosmict_7_phi,"cosmict_7_phi/F");
 
 for (Int_t i=0;i!=sig->GetEntries();i++){
    sig->GetEntry(i);
     
    Stree->Fill();
 }

 for (Int_t i=0;i!=bkg->GetEntries();i++){
   bkg->GetEntry(i);
    
   Btree->Fill();
 }



 cout << "signal tree entries: " << Stree->GetEntries() << " / " << Stree->GetEntries("cosmict_7_filled>0 && cosmict_flag_1==0")<< endl;
 cout << "background tree entries: " << Btree->GetEntries() << " / " << Btree->GetEntries("cosmict_7_filled>0 && (cosmict_flag_7>0)&& cosmict_flag_1==0")<< endl;
  
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
  
  TString fname = "./round_0.root";
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
 

  dataloader->AddVariable( "cosmict_7_flag_sec", "cosmict_7_flag_sec", "", 'F' );
  dataloader->AddVariable( "cosmict_7_n_muon_tracks", "cosmict_7_n_muon_tracks", "", 'F' );
  dataloader->AddVariable( "cosmict_7_total_shower_length", "cosmict_7_total_shower_length", "cm", 'F' );
  dataloader->AddVariable( "cosmict_7_flag_inside", "cosmict_7_flag_inside", "", 'F' );
  dataloader->AddVariable( "cosmict_7_angle_beam", "cosmict_7_angle_beam", "deg", 'F' );
  dataloader->AddVariable( "cosmict_7_flag_dir_weak", "cosmict_7_flag_dir_weak", "", 'F' );
  dataloader->AddVariable( "cosmict_7_dQ_dx_end", "cosmict_7_dQ_dx_end", "MeV/cm", 'F' );
  dataloader->AddVariable( "cosmict_7_dQ_dx_front", "cosmict_7_dQ_dx_front", "MeV/cm", 'F' );
  dataloader->AddVariable( "cosmict_7_theta", "cosmict_7_theta", "deg", 'F' );
  dataloader->AddVariable( "cosmict_7_phi", "cosmict_7_phi", "deg", 'F' );
 
 

  
 
  TTree *signalTree     = (TTree*)input->Get("TreeS");
  TTree *backgroundTree = (TTree*)input->Get("TreeB");
  dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
  dataloader->AddBackgroundTree( backgroundTree, 1.0);
  // Set individual event weights (the variables must exist in the original TTree)
  dataloader->SetSignalWeightExpression( "weight" );
  dataloader->SetBackgroundWeightExpression( "weight" );

  // Apply additional cuts on the signal and background samples (can be different)
  TCut mycut_s = "cosmict_7_filled>0 && cosmict_flag_1==0"; // 34340
  TCut mycut_b = "cosmict_7_filled>0 && (cosmict_flag_7>0 ) && cosmict_flag_1==0"; //  2642
    
  dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					  "nTrain_Signal=20000:"
					  "nTrain_Background=2200:"
					  "nTest_Signal=4000:"
					  "nTest_Background=442:"
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
     dataloader->AddVariable( "cosmict_7_flag_sec", "cosmict_7_flag_sec", "", 'F' );
  dataloader->AddVariable( "cosmict_7_n_muon_tracks", "cosmict_7_n_muon_tracks", "", 'F' );
  dataloader->AddVariable( "cosmict_7_total_shower_length", "cosmict_7_total_shower_length", "cm", 'F' );
  dataloader->AddVariable( "cosmict_7_flag_inside", "cosmict_7_flag_inside", "", 'F' );
  dataloader->AddVariable( "cosmict_7_angle_beam", "cosmict_7_angle_beam", "deg", 'F' );
  dataloader->AddVariable( "cosmict_7_flag_dir_weak", "cosmict_7_flag_dir_weak", "", 'F' );
  dataloader->AddVariable( "cosmict_7_dQ_dx_end", "cosmict_7_dQ_dx_end", "MeV/cm", 'F' );
  dataloader->AddVariable( "cosmict_7_dQ_dx_front", "cosmict_7_dQ_dx_front", "MeV/cm", 'F' );
  dataloader->AddVariable( "cosmict_7_theta", "cosmict_7_theta", "deg", 'F' );
  dataloader->AddVariable( "cosmict_7_phi", "cosmict_7_phi", "deg", 'F' );
 
    
   
    
    TTree *signalTree     = (TTree*)input->Get("TreeS");
    TTree *backgroundTree = (TTree*)input->Get("TreeB");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight" );
    dataloader->SetBackgroundWeightExpression( "weight" );

    // Apply additional cuts on the signal and background samples (can be different)
    TCut mycut_s = "cosmict_7_filled>0 && cosmict_flag_1==0"; // 
    TCut mycut_b = "cosmict_7_filled>0 && (cosmict_flag_7>0  || cosmict_7_bdt < 0.1) && cosmict_flag_1==0"; // 3382

      dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					  "nTrain_Signal=20000:"
					  "nTrain_Background=2800:"
					  "nTest_Signal=4000:"
					  "nTest_Background=582:"
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
  float cosmict_flag_1;

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

  
  
  
  
  
  float weight;
  
  TFile *file = new TFile("round_0.root");
  TTree *sig = (TTree*)file->Get("TreeS");
  TTree *bkg = (TTree*)file->Get("TreeB");

  sig->SetBranchAddress("cosmict_flag_1",&cosmict_flag_1);
  
  sig->SetBranchAddress("cosmict_flag_7",&cosmict_flag_7);
  sig->SetBranchAddress("cosmict_7_filled",&cosmict_7_filled);
  sig->SetBranchAddress("cosmict_7_flag_sec",&cosmict_7_flag_sec);
  sig->SetBranchAddress("cosmict_7_n_muon_tracks",&cosmict_7_n_muon_tracks);
  sig->SetBranchAddress("cosmict_7_total_shower_length",&cosmict_7_total_shower_length);
  sig->SetBranchAddress("cosmict_7_flag_inside",&cosmict_7_flag_inside);
  sig->SetBranchAddress("cosmict_7_angle_beam",&cosmict_7_angle_beam);
  sig->SetBranchAddress("cosmict_7_flag_dir_weak",&cosmict_7_flag_dir_weak);
  sig->SetBranchAddress("cosmict_7_dQ_dx_end",&cosmict_7_dQ_dx_end);
  sig->SetBranchAddress("cosmict_7_dQ_dx_front",&cosmict_7_dQ_dx_front);
  sig->SetBranchAddress("cosmict_7_theta",&cosmict_7_theta);
  sig->SetBranchAddress("cosmict_7_phi",&cosmict_7_phi);
  
  
  
  bkg->SetBranchAddress("cosmict_flag_1",&cosmict_flag_1);
  
  bkg->SetBranchAddress("cosmict_flag_7",&cosmict_flag_7);
  bkg->SetBranchAddress("cosmict_7_filled",&cosmict_7_filled);
  bkg->SetBranchAddress("cosmict_7_flag_sec",&cosmict_7_flag_sec);
  bkg->SetBranchAddress("cosmict_7_n_muon_tracks",&cosmict_7_n_muon_tracks);
  bkg->SetBranchAddress("cosmict_7_total_shower_length",&cosmict_7_total_shower_length);
  bkg->SetBranchAddress("cosmict_7_flag_inside",&cosmict_7_flag_inside);
  bkg->SetBranchAddress("cosmict_7_angle_beam",&cosmict_7_angle_beam);
  bkg->SetBranchAddress("cosmict_7_flag_dir_weak",&cosmict_7_flag_dir_weak);
  bkg->SetBranchAddress("cosmict_7_dQ_dx_end",&cosmict_7_dQ_dx_end);
  bkg->SetBranchAddress("cosmict_7_dQ_dx_front",&cosmict_7_dQ_dx_front);
  bkg->SetBranchAddress("cosmict_7_theta",&cosmict_7_theta);
  bkg->SetBranchAddress("cosmict_7_phi",&cosmict_7_phi);
  
  
  TFile *new_file = new TFile(filename,"RECREATE");
  TTree *Tsig = new TTree("TreeS","TreeS");
  TTree *Tbkg = new TTree("TreeB","TreeB");
  Tsig->SetDirectory(new_file);
  Tbkg->SetDirectory(new_file);

  
  Tsig->Branch("weight",&weight,"weight/F");
  Tsig->Branch("cosmict_flag_1",&cosmict_flag_1,"cosmict_flag_1/F");
 

  Tbkg->Branch("weight",&weight,"weight/F");
  Tbkg->Branch("cosmict_flag_1",&cosmict_flag_1,"cosmict_flag_1/F");
  
  Tsig->Branch("cosmict_flag_7",&cosmict_flag_7,"cosmict_flag_7/F");
  Tsig->Branch("cosmict_7_filled",&cosmict_7_filled,"cosmict_7_filled/F");
  Tsig->Branch("cosmict_7_flag_sec",&cosmict_7_flag_sec,"cosmict_7_flag_sec/F");
  Tsig->Branch("cosmict_7_n_muon_tracks",&cosmict_7_n_muon_tracks,"cosmict_7_n_muon_tracks/F");
  Tsig->Branch("cosmict_7_total_shower_length",&cosmict_7_total_shower_length,"cosmict_7_total_shower_length/F");
  Tsig->Branch("cosmict_7_flag_inside",&cosmict_7_flag_inside,"cosmict_7_flag_inside/F");
  Tsig->Branch("cosmict_7_angle_beam",&cosmict_7_angle_beam,"cosmict_7_angle_beam/F");
  Tsig->Branch("cosmict_7_flag_dir_weak",&cosmict_7_flag_dir_weak,"cosmict_7_flag_dir_weak/F");
  Tsig->Branch("cosmict_7_dQ_dx_end",&cosmict_7_dQ_dx_end,"cosmict_7_dQ_dx_end/F");
  Tsig->Branch("cosmict_7_dQ_dx_front",&cosmict_7_dQ_dx_front,"cosmict_7_dQ_dx_front/F");
  Tsig->Branch("cosmict_7_theta",&cosmict_7_theta,"cosmict_7_theta/F");
  Tsig->Branch("cosmict_7_phi",&cosmict_7_phi,"cosmict_7_phi/F");


  Tbkg->Branch("cosmict_flag_7",&cosmict_flag_7,"cosmict_flag_7/F");
  Tbkg->Branch("cosmict_7_filled",&cosmict_7_filled,"cosmict_7_filled/F");
  Tbkg->Branch("cosmict_7_flag_sec",&cosmict_7_flag_sec,"cosmict_7_flag_sec/F");
  Tbkg->Branch("cosmict_7_n_muon_tracks",&cosmict_7_n_muon_tracks,"cosmict_7_n_muon_tracks/F");
  Tbkg->Branch("cosmict_7_total_shower_length",&cosmict_7_total_shower_length,"cosmict_7_total_shower_length/F");
  Tbkg->Branch("cosmict_7_flag_inside",&cosmict_7_flag_inside,"cosmict_7_flag_inside/F");
  Tbkg->Branch("cosmict_7_angle_beam",&cosmict_7_angle_beam,"cosmict_7_angle_beam/F");
  Tbkg->Branch("cosmict_7_flag_dir_weak",&cosmict_7_flag_dir_weak,"cosmict_7_flag_dir_weak/F");
  Tbkg->Branch("cosmict_7_dQ_dx_end",&cosmict_7_dQ_dx_end,"cosmict_7_dQ_dx_end/F");
  Tbkg->Branch("cosmict_7_dQ_dx_front",&cosmict_7_dQ_dx_front,"cosmict_7_dQ_dx_front/F");
  Tbkg->Branch("cosmict_7_theta",&cosmict_7_theta,"cosmict_7_theta/F");
  Tbkg->Branch("cosmict_7_phi",&cosmict_7_phi,"cosmict_7_phi/F");
  
  
  sig->SetBranchAddress("weight",&weight);
  bkg->SetBranchAddress("weight",&weight);
  
  float bdt_value;
  Tsig->Branch("cosmict_7_bdt", &bdt_value, "cosmict_7_bdt/F");
  Tbkg->Branch("cosmict_7_bdt", &bdt_value, "cosmict_7_bdt/F");

  TMVA::Reader *reader = new TMVA::Reader();
  reader->AddVariable( "cosmict_7_flag_sec", &cosmict_7_flag_sec);
  reader->AddVariable( "cosmict_7_n_muon_tracks", &cosmict_7_n_muon_tracks);
  reader->AddVariable( "cosmict_7_total_shower_length", &cosmict_7_total_shower_length);
  reader->AddVariable( "cosmict_7_flag_inside", &cosmict_7_flag_inside);
  reader->AddVariable( "cosmict_7_angle_beam", &cosmict_7_angle_beam);
  reader->AddVariable( "cosmict_7_flag_dir_weak", &cosmict_7_flag_dir_weak);
  reader->AddVariable( "cosmict_7_dQ_dx_end", &cosmict_7_dQ_dx_end);
  reader->AddVariable( "cosmict_7_dQ_dx_front", &cosmict_7_dQ_dx_front);
  reader->AddVariable( "cosmict_7_theta", &cosmict_7_theta);
  reader->AddVariable( "cosmict_7_phi", &cosmict_7_phi);
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

