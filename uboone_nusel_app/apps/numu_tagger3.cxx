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

 float numu_cc_flag_3;
 float numu_cc_3_particle_type;
 float numu_cc_3_max_length;
 float numu_cc_3_acc_track_length;
 float numu_cc_3_max_length_all;
 float numu_cc_3_max_muon_length;
 float numu_cc_3_n_daughter_tracks;
 float numu_cc_3_n_daughter_all;
 
 sig->SetBranchAddress("cosmict_flag_1",&cosmict_flag_1);
 
 sig->SetBranchAddress("numu_cc_flag_3",&numu_cc_flag_3);
 sig->SetBranchAddress("numu_cc_3_particle_type",&numu_cc_3_particle_type);
 sig->SetBranchAddress("numu_cc_3_max_length",&numu_cc_3_max_length);
 sig->SetBranchAddress("numu_cc_3_track_length",&numu_cc_3_acc_track_length);
 sig->SetBranchAddress("numu_cc_3_max_length_all",&numu_cc_3_max_length_all);
 sig->SetBranchAddress("numu_cc_3_max_muon_length",&numu_cc_3_max_muon_length);
 sig->SetBranchAddress("numu_cc_3_n_daughter_tracks",&numu_cc_3_n_daughter_tracks);
 sig->SetBranchAddress("numu_cc_3_n_daughter_all",&numu_cc_3_n_daughter_all);

 bkg->SetBranchAddress("cosmict_flag_1",&cosmict_flag_1);
 
 bkg->SetBranchAddress("numu_cc_flag_3",&numu_cc_flag_3);
 bkg->SetBranchAddress("numu_cc_3_particle_type",&numu_cc_3_particle_type);
 bkg->SetBranchAddress("numu_cc_3_max_length",&numu_cc_3_max_length);
 bkg->SetBranchAddress("numu_cc_3_track_length",&numu_cc_3_acc_track_length);
 bkg->SetBranchAddress("numu_cc_3_max_length_all",&numu_cc_3_max_length_all);
 bkg->SetBranchAddress("numu_cc_3_max_muon_length",&numu_cc_3_max_muon_length);
 bkg->SetBranchAddress("numu_cc_3_n_daughter_tracks",&numu_cc_3_n_daughter_tracks);
 bkg->SetBranchAddress("numu_cc_3_n_daughter_all",&numu_cc_3_n_daughter_all);



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

 
 Stree->Branch("numu_cc_flag_3",&numu_cc_flag_3,"numu_cc_flag_3/F");
 Stree->Branch("numu_cc_3_particle_type",&numu_cc_3_particle_type,"numu_cc_3_particle_type/F");
 Stree->Branch("numu_cc_3_max_length",&numu_cc_3_max_length,"numu_cc_3_max_length/F");
 Stree->Branch("numu_cc_3_track_length",&numu_cc_3_acc_track_length,"numu_cc_3_acc_track_length/F");
 Stree->Branch("numu_cc_3_max_length_all",&numu_cc_3_max_length_all,"numu_cc_3_max_length_all/F");
 Stree->Branch("numu_cc_3_max_muon_length",&numu_cc_3_max_muon_length,"numu_cc_3_max_muon_length/F");
 Stree->Branch("numu_cc_3_n_daughter_tracks",&numu_cc_3_n_daughter_tracks,"numu_cc_3_n_daughter_tracks/F");
 Stree->Branch("numu_cc_3_n_daughter_all",&numu_cc_3_n_daughter_all,"numu_cc_3_n_daughter_all/F");
 

 Btree->Branch("numu_cc_flag_3",&numu_cc_flag_3,"numu_cc_flag_3/F");
 Btree->Branch("numu_cc_3_particle_type",&numu_cc_3_particle_type,"numu_cc_3_particle_type/F");
 Btree->Branch("numu_cc_3_max_length",&numu_cc_3_max_length,"numu_cc_3_max_length/F");
 Btree->Branch("numu_cc_3_track_length",&numu_cc_3_acc_track_length,"numu_cc_3_acc_track_length/F");
 Btree->Branch("numu_cc_3_max_length_all",&numu_cc_3_max_length_all,"numu_cc_3_max_length_all/F");
 Btree->Branch("numu_cc_3_max_muon_length",&numu_cc_3_max_muon_length,"numu_cc_3_max_muon_length/F");
 Btree->Branch("numu_cc_3_n_daughter_tracks",&numu_cc_3_n_daughter_tracks,"numu_cc_3_n_daughter_tracks/F");
 Btree->Branch("numu_cc_3_n_daughter_all",&numu_cc_3_n_daughter_all,"numu_cc_3_n_daughter_all/F");

 
 for (Int_t i=0;i!=sig->GetEntries();i++){
    sig->GetEntry(i);


     
    Stree->Fill();
 }

 for (Int_t i=0;i!=bkg->GetEntries();i++){
   bkg->GetEntry(i);

   
    
   Btree->Fill();
 }

 

 cout << "signal tree entries: " << Stree->GetEntries() << " / " << Stree->GetEntries("cosmict_flag_1==0")<< endl;
 cout << "background tree entries: " << Btree->GetEntries() << " / " << Btree->GetEntries(" (numu_cc_flag_3==0 )&& cosmict_flag_1==0")<< endl;
  
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
 

 dataloader->AddVariable( "numu_cc_3_particle_type", "numu_cc_3_particle_type", "", 'F' );
 dataloader->AddVariable( "numu_cc_3_max_length", "numu_cc_3_max_length", "cm", 'F' );
 dataloader->AddVariable( "numu_cc_3_acc_track_length", "numu_cc_3_acc_track_length", "cm", 'F' );
 dataloader->AddVariable( "numu_cc_3_max_length_all", "numu_cc_3_max_length_all", "cm", 'F' );
 dataloader->AddVariable( "numu_cc_3_max_muon_length", "numu_cc_3_max_muon_length", "cm", 'F' );
 dataloader->AddVariable( "numu_cc_3_n_daughter_tracks", "numu_cc_3_n_daughter_tracks", "", 'F' );
 dataloader->AddVariable( "numu_cc_3_n_daughter_all", "numu_cc_3_n_daughter_all", "", 'F' );


  
 
  TTree *signalTree     = (TTree*)input->Get("TreeS");
  TTree *backgroundTree = (TTree*)input->Get("TreeB");
  dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
  dataloader->AddBackgroundTree( backgroundTree, 1.0);
  // Set individual event weights (the variables must exist in the original TTree)
  dataloader->SetSignalWeightExpression( "weight" );
  dataloader->SetBackgroundWeightExpression( "weight" );

  // Apply additional cuts on the signal and background samples (can be different)
  TCut mycut_s = " cosmict_flag_1==0"; // 233853
  TCut mycut_b = " cosmict_flag_1==0"; // 70965
    
  dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					  "nTrain_Signal=80000:"
					  "nTrain_Background=60000:"
					  "nTest_Signal=20000:"
					  "nTest_Background=10000:"
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
    
 dataloader->AddVariable( "numu_cc_3_particle_type", "numu_cc_3_particle_type", "", 'F' );
 dataloader->AddVariable( "numu_cc_3_max_length", "numu_cc_3_max_length", "cm", 'F' );
 dataloader->AddVariable( "numu_cc_3_acc_track_length", "numu_cc_3_acc_track_length", "cm", 'F' );
 dataloader->AddVariable( "numu_cc_3_max_length_all", "numu_cc_3_max_length_all", "cm", 'F' );
 dataloader->AddVariable( "numu_cc_3_max_muon_length", "numu_cc_3_max_muon_length", "cm", 'F' );
 dataloader->AddVariable( "numu_cc_3_n_daughter_tracks", "numu_cc_3_n_daughter_tracks", "", 'F' );
 dataloader->AddVariable( "numu_cc_3_n_daughter_all", "numu_cc_3_n_daughter_all", "", 'F' );

    
   
    
    
    TTree *signalTree     = (TTree*)input->Get("TreeS");
    TTree *backgroundTree = (TTree*)input->Get("TreeB");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight" );
    dataloader->SetBackgroundWeightExpression( "weight" );

     // Apply additional cuts on the signal and background samples (can be different)
  TCut mycut_s = " cosmict_flag_1==0"; // 
  TCut mycut_b = " (numu_cc_flag_3==0 || numu_cc_3_bdt < 0) && cosmict_flag_1==0"; // 
    
  dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					  "nTrain_Signal=20000:"
					  "nTrain_Background=4500:"
					  "nTest_Signal=11000:"
					  "nTest_Background=728:"
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

  float numu_cc_flag_3;
  float numu_cc_3_particle_type;
  float numu_cc_3_max_length;
  float numu_cc_3_acc_track_length;
  float numu_cc_3_max_length_all;
  float numu_cc_3_max_muon_length;
  float numu_cc_3_n_daughter_tracks;
  float numu_cc_3_n_daughter_all;
  
   
  
  float weight;
  
  TFile *file = new TFile("round_0.root");
  TTree *sig = (TTree*)file->Get("TreeS");
  TTree *bkg = (TTree*)file->Get("TreeB");

  sig->SetBranchAddress("cosmict_flag_1",&cosmict_flag_1);
  
 sig->SetBranchAddress("numu_cc_flag_3",&numu_cc_flag_3);
 sig->SetBranchAddress("numu_cc_3_particle_type",&numu_cc_3_particle_type);
 sig->SetBranchAddress("numu_cc_3_max_length",&numu_cc_3_max_length);
 sig->SetBranchAddress("numu_cc_3_track_length",&numu_cc_3_acc_track_length);
 sig->SetBranchAddress("numu_cc_3_max_length_all",&numu_cc_3_max_length_all);
 sig->SetBranchAddress("numu_cc_3_max_muon_length",&numu_cc_3_max_muon_length);
 sig->SetBranchAddress("numu_cc_3_n_daughter_tracks",&numu_cc_3_n_daughter_tracks);
 sig->SetBranchAddress("numu_cc_3_n_daughter_all",&numu_cc_3_n_daughter_all);

 
  
  bkg->SetBranchAddress("cosmict_flag_1",&cosmict_flag_1);

  bkg->SetBranchAddress("numu_cc_flag_3",&numu_cc_flag_3);
 bkg->SetBranchAddress("numu_cc_3_particle_type",&numu_cc_3_particle_type);
 bkg->SetBranchAddress("numu_cc_3_max_length",&numu_cc_3_max_length);
 bkg->SetBranchAddress("numu_cc_3_track_length",&numu_cc_3_acc_track_length);
 bkg->SetBranchAddress("numu_cc_3_max_length_all",&numu_cc_3_max_length_all);
 bkg->SetBranchAddress("numu_cc_3_max_muon_length",&numu_cc_3_max_muon_length);
 bkg->SetBranchAddress("numu_cc_3_n_daughter_tracks",&numu_cc_3_n_daughter_tracks);
 bkg->SetBranchAddress("numu_cc_3_n_daughter_all",&numu_cc_3_n_daughter_all);
  
  
  
  TFile *new_file = new TFile(filename,"RECREATE");
  TTree *Tsig = new TTree("TreeS","TreeS");
  TTree *Tbkg = new TTree("TreeB","TreeB");
  Tsig->SetDirectory(new_file);
  Tbkg->SetDirectory(new_file);

  
  Tsig->Branch("weight",&weight,"weight/F");
  Tsig->Branch("cosmict_flag_1",&cosmict_flag_1,"cosmict_flag_1/F");
 

  Tbkg->Branch("weight",&weight,"weight/F");
  Tbkg->Branch("cosmict_flag_1",&cosmict_flag_1,"cosmict_flag_1/F");
  
 
  Tsig->Branch("numu_cc_flag_3",&numu_cc_flag_3,"numu_cc_flag_3/F");
 Tsig->Branch("numu_cc_3_particle_type",&numu_cc_3_particle_type,"numu_cc_3_particle_type/F");
 Tsig->Branch("numu_cc_3_max_length",&numu_cc_3_max_length,"numu_cc_3_max_length/F");
 Tsig->Branch("numu_cc_3_track_length",&numu_cc_3_acc_track_length,"numu_cc_3_acc_track_length/F");
 Tsig->Branch("numu_cc_3_max_length_all",&numu_cc_3_max_length_all,"numu_cc_3_max_length_all/F");
 Tsig->Branch("numu_cc_3_max_muon_length",&numu_cc_3_max_muon_length,"numu_cc_3_max_muon_length/F");
 Tsig->Branch("numu_cc_3_n_daughter_tracks",&numu_cc_3_n_daughter_tracks,"numu_cc_3_n_daughter_tracks/F");
 Tsig->Branch("numu_cc_3_n_daughter_all",&numu_cc_3_n_daughter_all,"numu_cc_3_n_daughter_all/F");


 Tbkg->Branch("numu_cc_flag_3",&numu_cc_flag_3,"numu_cc_flag_3/F");
 Tbkg->Branch("numu_cc_3_particle_type",&numu_cc_3_particle_type,"numu_cc_3_particle_type/F");
 Tbkg->Branch("numu_cc_3_max_length",&numu_cc_3_max_length,"numu_cc_3_max_length/F");
 Tbkg->Branch("numu_cc_3_track_length",&numu_cc_3_acc_track_length,"numu_cc_3_acc_track_length/F");
 Tbkg->Branch("numu_cc_3_max_length_all",&numu_cc_3_max_length_all,"numu_cc_3_max_length_all/F");
 Tbkg->Branch("numu_cc_3_max_muon_length",&numu_cc_3_max_muon_length,"numu_cc_3_max_muon_length/F");
 Tbkg->Branch("numu_cc_3_n_daughter_tracks",&numu_cc_3_n_daughter_tracks,"numu_cc_3_n_daughter_tracks/F");
 Tbkg->Branch("numu_cc_3_n_daughter_all",&numu_cc_3_n_daughter_all,"numu_cc_3_n_daughter_all/F");
 


 
  
  sig->SetBranchAddress("weight",&weight);
  bkg->SetBranchAddress("weight",&weight);
  
  float bdt_value;
  Tsig->Branch("numu_cc_3_bdt", &bdt_value, "numu_cc_3_bdt/F");
  Tbkg->Branch("numu_cc_3_bdt", &bdt_value, "numu_cc_3_bdt/F");
  
  TMVA::Reader *reader = new TMVA::Reader();
  reader->AddVariable( "numu_cc_3_particle_type", &numu_cc_3_particle_type);
  reader->AddVariable( "numu_cc_3_max_length", &numu_cc_3_max_length);
  reader->AddVariable( "numu_cc_3_acc_track_length", &numu_cc_3_acc_track_length);
  reader->AddVariable( "numu_cc_3_max_length_all", &numu_cc_3_max_length_all);
  reader->AddVariable( "numu_cc_3_max_muon_length", &numu_cc_3_max_muon_length);
  reader->AddVariable( "numu_cc_3_n_daughter_tracks", &numu_cc_3_n_daughter_tracks);
  reader->AddVariable( "numu_cc_3_n_daughter_all", &numu_cc_3_n_daughter_all);

    
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

