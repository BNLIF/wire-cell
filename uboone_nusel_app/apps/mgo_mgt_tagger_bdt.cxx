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

  
  double mgo_energy;
  double mgo_max_energy;
  double mgo_total_energy;
  int mgo_n_showers;
  double mgo_max_energy_1;
  double mgo_max_energy_2;
  double mgo_total_other_energy;
  int mgo_n_total_showers;
  double mgo_total_other_energy_1;
  int mgo_flag;

  int mgt_flag_single_shower;
  double mgt_max_energy;
  double mgt_total_other_energy;
  double mgt_max_energy_1;
  double mgt_e_indirect_max_energy;
  double mgt_e_direct_max_energy;
  int mgt_n_direct_showers;
  double mgt_e_direct_total_energy;
  int mgt_flag_indirect_max_pio;
  double mgt_e_indirect_total_energy;
  int mgt_flag;
  
  
  sig->SetBranchAddress("mgo_energy",&mgo_energy);
  sig->SetBranchAddress("mgo_max_energy",&mgo_max_energy);
  sig->SetBranchAddress("mgo_total_energy",&mgo_total_energy);
  sig->SetBranchAddress("mgo_n_showers",&mgo_n_showers);
  sig->SetBranchAddress("mgo_max_energy_1",&mgo_max_energy_1);
  sig->SetBranchAddress("mgo_max_energy_2",&mgo_max_energy_2);
  sig->SetBranchAddress("mgo_total_other_energy",&mgo_total_other_energy);
  sig->SetBranchAddress("mgo_n_total_showers",&mgo_n_total_showers);
  sig->SetBranchAddress("mgo_total_other_energy_1",&mgo_total_other_energy_1);
  sig->SetBranchAddress("mgo_flag",&mgo_flag);

  sig->SetBranchAddress("mgt_flag_single_shower",&mgt_flag_single_shower);
  sig->SetBranchAddress("mgt_max_energy",&mgt_max_energy);
  sig->SetBranchAddress("mgt_total_other_energy",&mgt_total_other_energy);
  sig->SetBranchAddress("mgt_max_energy_1",&mgt_max_energy_1);
  sig->SetBranchAddress("mgt_e_indirect_max_energy",&mgt_e_indirect_max_energy);
  sig->SetBranchAddress("mgt_e_direct_max_energy",&mgt_e_direct_max_energy);
  sig->SetBranchAddress("mgt_n_direct_showers",&mgt_n_direct_showers);
  sig->SetBranchAddress("mgt_e_direct_total_energy",&mgt_e_direct_total_energy);
  sig->SetBranchAddress("mgt_flag_indirect_max_pio",&mgt_flag_indirect_max_pio);
  sig->SetBranchAddress("mgt_e_indirect_total_energy",&mgt_e_indirect_total_energy);
  sig->SetBranchAddress("mgt_flag",&mgt_flag);

  bkg->SetBranchAddress("mgo_energy",&mgo_energy);
  bkg->SetBranchAddress("mgo_max_energy",&mgo_max_energy);
  bkg->SetBranchAddress("mgo_total_energy",&mgo_total_energy);
  bkg->SetBranchAddress("mgo_n_showers",&mgo_n_showers);
  bkg->SetBranchAddress("mgo_max_energy_1",&mgo_max_energy_1);
  bkg->SetBranchAddress("mgo_max_energy_2",&mgo_max_energy_2);
  bkg->SetBranchAddress("mgo_total_other_energy",&mgo_total_other_energy);
  bkg->SetBranchAddress("mgo_n_total_showers",&mgo_n_total_showers);
  bkg->SetBranchAddress("mgo_total_other_energy_1",&mgo_total_other_energy_1);
  bkg->SetBranchAddress("mgo_flag",&mgo_flag);

  bkg->SetBranchAddress("mgt_flag_single_shower",&mgt_flag_single_shower);
  bkg->SetBranchAddress("mgt_max_energy",&mgt_max_energy);
  bkg->SetBranchAddress("mgt_total_other_energy",&mgt_total_other_energy);
  bkg->SetBranchAddress("mgt_max_energy_1",&mgt_max_energy_1);
  bkg->SetBranchAddress("mgt_e_indirect_max_energy",&mgt_e_indirect_max_energy);
  bkg->SetBranchAddress("mgt_e_direct_max_energy",&mgt_e_direct_max_energy);
  bkg->SetBranchAddress("mgt_n_direct_showers",&mgt_n_direct_showers);
  bkg->SetBranchAddress("mgt_e_direct_total_energy",&mgt_e_direct_total_energy);
  bkg->SetBranchAddress("mgt_flag_indirect_max_pio",&mgt_flag_indirect_max_pio);
  bkg->SetBranchAddress("mgt_e_indirect_total_energy",&mgt_e_indirect_total_energy);
  bkg->SetBranchAddress("mgt_flag",&mgt_flag);

  
  
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
  
  float mgo_energy_f;
  float mgo_max_energy_f;
  float mgo_total_energy_f;
  float mgo_n_showers_f;
  float mgo_max_energy_1_f;
  float mgo_max_energy_2_f;
  float mgo_total_other_energy_f;
  float mgo_n_total_showers_f;
  float mgo_total_other_energy_1_f;
  float mgo_flag_f;

  float mgt_flag_single_shower_f;
  float mgt_max_energy_f;
  float mgt_total_other_energy_f;
  float mgt_max_energy_1_f;
  float mgt_e_indirect_max_energy_f;
  float mgt_e_direct_max_energy_f;
  float mgt_n_direct_showers_f;
  float mgt_e_direct_total_energy_f;
  float mgt_flag_indirect_max_pio_f;
  float mgt_e_indirect_total_energy_f;
  float mgt_flag_f;

  
  Tsig->Branch("mgo_energy",&mgo_energy_f,"mgo_energy/F");
  Tsig->Branch("mgo_max_energy",&mgo_max_energy_f,"mgo_max_energy/F");
  Tsig->Branch("mgo_total_energy",&mgo_total_energy_f,"mgo_total_energy/F");
  Tsig->Branch("mgo_n_showers",&mgo_n_showers_f,"mgo_n_showers/F");
  Tsig->Branch("mgo_max_energy_1",&mgo_max_energy_1_f,"mgo_max_energy_1/F");
  Tsig->Branch("mgo_max_energy_2",&mgo_max_energy_2_f,"mgo_max_energy_2/F");
  Tsig->Branch("mgo_total_other_energy",&mgo_total_other_energy_f,"mgo_total_other_energy/F");
  Tsig->Branch("mgo_n_total_showers",&mgo_n_total_showers_f,"mgo_n_total_showers/F");
  Tsig->Branch("mgo_total_other_energy_1",&mgo_total_other_energy_1_f,"mgo_total_other_energy_1/F");
  Tsig->Branch("mgo_flag",&mgo_flag_f,"mgo_flag/F");
  
  Tsig->Branch("mgt_flag_single_shower",&mgt_flag_single_shower_f,"mgt_flag_single_shower/F");
  Tsig->Branch("mgt_max_energy",&mgt_max_energy_f,"mgt_max_energy/F");
  Tsig->Branch("mgt_total_other_energy",&mgt_total_other_energy_f,"mgt_total_other_energy/F");
  Tsig->Branch("mgt_max_energy_1",&mgt_max_energy_1_f,"mgt_max_energy_1/F");
  Tsig->Branch("mgt_e_indirect_max_energy",&mgt_e_indirect_max_energy_f,"mgt_e_indirect_max_energy/F");
  Tsig->Branch("mgt_e_direct_max_energy",&mgt_e_direct_max_energy_f,"mgt_e_direct_max_energy/F");
  Tsig->Branch("mgt_n_direct_showers",&mgt_n_direct_showers_f,"mgt_n_direct_showers/F");
  Tsig->Branch("mgt_e_direct_total_energy",&mgt_e_direct_total_energy_f,"mgt_e_direct_total_energy/F");
  Tsig->Branch("mgt_flag_indirect_max_pio",&mgt_flag_indirect_max_pio_f,"mgt_flag_indirect_max_pio/F");
  Tsig->Branch("mgt_e_indirect_total_energy",&mgt_e_indirect_total_energy_f,"mgt_e_indirect_total_energy/F");
  Tsig->Branch("mgt_flag",&mgt_flag_f,"mgt_flag/F");

  
  Tbkg->Branch("mgo_energy",&mgo_energy_f,"mgo_energy/F");
  Tbkg->Branch("mgo_max_energy",&mgo_max_energy_f,"mgo_max_energy/F");
  Tbkg->Branch("mgo_total_energy",&mgo_total_energy_f,"mgo_total_energy/F");
  Tbkg->Branch("mgo_n_showers",&mgo_n_showers_f,"mgo_n_showers/F");
  Tbkg->Branch("mgo_max_energy_1",&mgo_max_energy_1_f,"mgo_max_energy_1/F");
  Tbkg->Branch("mgo_max_energy_2",&mgo_max_energy_2_f,"mgo_max_energy_2/F");
  Tbkg->Branch("mgo_total_other_energy",&mgo_total_other_energy_f,"mgo_total_other_energy/F");
  Tbkg->Branch("mgo_n_total_showers",&mgo_n_total_showers_f,"mgo_n_total_showers/F");
  Tbkg->Branch("mgo_total_other_energy_1",&mgo_total_other_energy_1_f,"mgo_total_other_energy_1/F");
  Tbkg->Branch("mgo_flag",&mgo_flag_f,"mgo_flag/F");
  
  Tbkg->Branch("mgt_flag_single_shower",&mgt_flag_single_shower_f,"mgt_flag_single_shower/F");
  Tbkg->Branch("mgt_max_energy",&mgt_max_energy_f,"mgt_max_energy/F");
  Tbkg->Branch("mgt_total_other_energy",&mgt_total_other_energy_f,"mgt_total_other_energy/F");
  Tbkg->Branch("mgt_max_energy_1",&mgt_max_energy_1_f,"mgt_max_energy_1/F");
  Tbkg->Branch("mgt_e_indirect_max_energy",&mgt_e_indirect_max_energy_f,"mgt_e_indirect_max_energy/F");
  Tbkg->Branch("mgt_e_direct_max_energy",&mgt_e_direct_max_energy_f,"mgt_e_direct_max_energy/F");
  Tbkg->Branch("mgt_n_direct_showers",&mgt_n_direct_showers_f,"mgt_n_direct_showers/F");
  Tbkg->Branch("mgt_e_direct_total_energy",&mgt_e_direct_total_energy_f,"mgt_e_direct_total_energy/F");
  Tbkg->Branch("mgt_flag_indirect_max_pio",&mgt_flag_indirect_max_pio_f,"mgt_flag_indirect_max_pio/F");
  Tbkg->Branch("mgt_e_indirect_total_energy",&mgt_e_indirect_total_energy_f,"mgt_e_indirect_total_energy/F");
  Tbkg->Branch("mgt_flag",&mgt_flag_f,"mgt_flag/F");
  
  
  for (Int_t i=0;i!=sig->GetEntries();i++){
    sig->GetEntry(i);
    
    mgo_energy_f = mgo_energy;
    mgo_max_energy_f = mgo_max_energy;
    mgo_total_energy_f = mgo_total_energy;
    mgo_n_showers_f = mgo_n_showers;
    mgo_max_energy_1_f = mgo_max_energy_1;
    mgo_max_energy_2_f = mgo_max_energy_2;
    mgo_total_other_energy_f = mgo_total_other_energy;
    mgo_n_total_showers_f = mgo_n_total_showers;
    mgo_total_other_energy_1_f = mgo_total_other_energy_1;
    mgo_flag_f = mgo_flag;
    
    mgt_flag_single_shower_f = mgt_flag_single_shower;
    mgt_max_energy_f = mgt_max_energy;
    mgt_total_other_energy_f = mgt_total_other_energy;
    mgt_max_energy_1_f = mgt_max_energy_1;
    mgt_e_indirect_max_energy_f = mgt_e_indirect_max_energy;
    mgt_e_direct_max_energy_f = mgt_e_direct_max_energy;
    mgt_n_direct_showers_f = mgt_n_direct_showers;
    mgt_e_direct_total_energy_f = mgt_e_direct_total_energy;
    mgt_flag_indirect_max_pio_f = mgt_flag_indirect_max_pio;
    mgt_e_indirect_total_energy_f = mgt_e_indirect_total_energy;
    mgt_flag_f = mgt_flag;
    
    Tsig->Fill();
  }

   for (Int_t i=0;i!=bkg->GetEntries();i++){
    bkg->GetEntry(i);
    
     mgo_energy_f = mgo_energy;
    mgo_max_energy_f = mgo_max_energy;
    mgo_total_energy_f = mgo_total_energy;
    mgo_n_showers_f = mgo_n_showers;
    mgo_max_energy_1_f = mgo_max_energy_1;
    mgo_max_energy_2_f = mgo_max_energy_2;
    mgo_total_other_energy_f = mgo_total_other_energy;
    mgo_n_total_showers_f = mgo_n_total_showers;
    mgo_total_other_energy_1_f = mgo_total_other_energy_1;
    mgo_flag_f = mgo_flag;
    
    mgt_flag_single_shower_f = mgt_flag_single_shower;
    mgt_max_energy_f = mgt_max_energy;
    mgt_total_other_energy_f = mgt_total_other_energy;
    mgt_max_energy_1_f = mgt_max_energy_1;
    mgt_e_indirect_max_energy_f = mgt_e_indirect_max_energy;
    mgt_e_direct_max_energy_f = mgt_e_direct_max_energy;
    mgt_n_direct_showers_f = mgt_n_direct_showers;
    mgt_e_direct_total_energy_f = mgt_e_direct_total_energy;
    mgt_flag_indirect_max_pio_f = mgt_flag_indirect_max_pio;
    mgt_e_indirect_total_energy_f = mgt_e_indirect_total_energy;
    mgt_flag_f = mgt_flag;
    
    
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
    
    
    dataloader->AddVariable("mgo_energy","mgo_energy","MeV",'F');
    dataloader->AddVariable("mgo_max_energy","mgo_max_energy","MeV",'F');
    dataloader->AddVariable("mgo_total_energy","mgo_total_energy","MeV",'F');
    dataloader->AddVariable("mgo_n_showers","mgo_n_showers","",'F');
    dataloader->AddVariable("mgo_max_energy_1","mgo_max_energy_1","MeV",'F');
    dataloader->AddVariable("mgo_max_energy_2","mgo_max_energy_2","MeV",'F');
    dataloader->AddVariable("mgo_total_other_energy","mgo_total_other_energy","MeV",'F');
    dataloader->AddVariable("mgo_n_total_showers","mgo_n_total_showers","",'F');
    dataloader->AddVariable("mgo_total_other_energy_1","mgo_total_other_energy_1","MeV",'F');
    dataloader->AddVariable("mgt_flag_single_shower","mgt_flag_single_shower","",'F');
    dataloader->AddVariable("mgt_max_energy","mgt_max_energy","MeV",'F');
    dataloader->AddVariable("mgt_total_other_energy","mgt_total_other_energy","MeV",'F');
    dataloader->AddVariable("mgt_max_energy_1","mgt_max_energy_1","MeV",'F');
    dataloader->AddVariable("mgt_e_indirect_max_energy","mgt_e_indirect_max_energy","MeV",'F');
    dataloader->AddVariable("mgt_e_direct_max_energy","mgt_e_direct_max_energy","MeV",'F');
    dataloader->AddVariable("mgt_n_direct_showers","mgt_n_direct_showers","",'F');
    dataloader->AddVariable("mgt_e_direct_total_energy","mgt_e_direct_total_energy","MeV",'F');
    dataloader->AddVariable("mgt_flag_indirect_max_pio","mgt_flag_indirect_max_pio","",'F');
    dataloader->AddVariable("mgt_e_indirect_total_energy","mgt_e_indirect_total_energy","MeV",'F');
    
        TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "1>0"; //  1345 (mgo) + 2558 (mgt) or 2753 (total)/ 44194
    TCut mycut_b = "mgo_flag==0 || mgt_flag==0"; // 2109 (mgo) + 3400 (mgt) + 3817 (total) / 21070
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=20000:"
					    "nTrain_Background=3400:"
					    "nTest_Signal=10000:"
					    "nTest_Background=417:"
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
    
    
    dataloader->AddVariable("mgo_energy","mgo_energy","MeV",'F');
    dataloader->AddVariable("mgo_max_energy","mgo_max_energy","MeV",'F');
    dataloader->AddVariable("mgo_total_energy","mgo_total_energy","MeV",'F');
    dataloader->AddVariable("mgo_n_showers","mgo_n_showers","",'F');
    dataloader->AddVariable("mgo_max_energy_1","mgo_max_energy_1","MeV",'F');
    dataloader->AddVariable("mgo_max_energy_2","mgo_max_energy_2","MeV",'F');
    dataloader->AddVariable("mgo_total_other_energy","mgo_total_other_energy","MeV",'F');
    dataloader->AddVariable("mgo_n_total_showers","mgo_n_total_showers","",'F');
    dataloader->AddVariable("mgo_total_other_energy_1","mgo_total_other_energy_1","MeV",'F');

    dataloader->AddVariable("mgt_flag_single_shower","mgt_flag_single_shower","",'F');
    dataloader->AddVariable("mgt_max_energy","mgt_max_energy","MeV",'F');
    dataloader->AddVariable("mgt_total_other_energy","mgt_total_other_energy","MeV",'F');
    dataloader->AddVariable("mgt_max_energy_1","mgt_max_energy_1","MeV",'F');
    dataloader->AddVariable("mgt_e_indirect_max_energy","mgt_e_indirect_max_energy","MeV",'F');
    dataloader->AddVariable("mgt_e_direct_max_energy","mgt_e_direct_max_energy","MeV",'F');
    dataloader->AddVariable("mgt_n_direct_showers","mgt_n_direct_showers","",'F');
    dataloader->AddVariable("mgt_e_direct_total_energy","mgt_e_direct_total_energy","MeV",'F');
    dataloader->AddVariable("mgt_flag_indirect_max_pio","mgt_flag_indirect_max_pio","",'F');
    dataloader->AddVariable("mgt_e_indirect_total_energy","mgt_e_indirect_total_energy","MeV",'F');
    
    

    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "1>0"; // 
    TCut mycut_b = "mgo_flag==0 || mgt_flag==0 || mgo_mgt_bdt < 0"; // 5729
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=20000:"
					    "nTrain_Background=5000:"
					    "nTest_Signal=10000:"
					    "nTest_Background=729:"
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

  
  float mgo_energy;
  float mgo_max_energy;
  float mgo_total_energy;
  float mgo_n_showers;
  float mgo_max_energy_1;
  float mgo_max_energy_2;
  float mgo_total_other_energy;
  float mgo_n_total_showers;
  float mgo_total_other_energy_1;
  float mgo_flag;
  
  float mgt_flag_single_shower;
  float mgt_max_energy;
  float mgt_total_other_energy;
  float mgt_max_energy_1;
  float mgt_e_indirect_max_energy;
  float mgt_e_direct_max_energy;
  float mgt_n_direct_showers;
  float mgt_e_direct_total_energy;
  float mgt_flag_indirect_max_pio;
  float mgt_e_indirect_total_energy;
  float mgt_flag;
  
 
  sig->SetBranchAddress("mgo_energy",&mgo_energy);
  sig->SetBranchAddress("mgo_max_energy",&mgo_max_energy);
  sig->SetBranchAddress("mgo_total_energy",&mgo_total_energy);
  sig->SetBranchAddress("mgo_n_showers",&mgo_n_showers);
  sig->SetBranchAddress("mgo_max_energy_1",&mgo_max_energy_1);
  sig->SetBranchAddress("mgo_max_energy_2",&mgo_max_energy_2);
  sig->SetBranchAddress("mgo_total_other_energy",&mgo_total_other_energy);
  sig->SetBranchAddress("mgo_n_total_showers",&mgo_n_total_showers);
  sig->SetBranchAddress("mgo_total_other_energy_1",&mgo_total_other_energy_1);
  sig->SetBranchAddress("mgo_flag",&mgo_flag);

  sig->SetBranchAddress("mgt_flag_single_shower",&mgt_flag_single_shower);
  sig->SetBranchAddress("mgt_max_energy",&mgt_max_energy);
  sig->SetBranchAddress("mgt_total_other_energy",&mgt_total_other_energy);
  sig->SetBranchAddress("mgt_max_energy_1",&mgt_max_energy_1);
  sig->SetBranchAddress("mgt_e_indirect_max_energy",&mgt_e_indirect_max_energy);
  sig->SetBranchAddress("mgt_e_direct_max_energy",&mgt_e_direct_max_energy);
  sig->SetBranchAddress("mgt_n_direct_showers",&mgt_n_direct_showers);
  sig->SetBranchAddress("mgt_e_direct_total_energy",&mgt_e_direct_total_energy);
  sig->SetBranchAddress("mgt_flag_indirect_max_pio",&mgt_flag_indirect_max_pio);
  sig->SetBranchAddress("mgt_e_indirect_total_energy",&mgt_e_indirect_total_energy);
  sig->SetBranchAddress("mgt_flag",&mgt_flag);

  bkg->SetBranchAddress("mgo_energy",&mgo_energy);
  bkg->SetBranchAddress("mgo_max_energy",&mgo_max_energy);
  bkg->SetBranchAddress("mgo_total_energy",&mgo_total_energy);
  bkg->SetBranchAddress("mgo_n_showers",&mgo_n_showers);
  bkg->SetBranchAddress("mgo_max_energy_1",&mgo_max_energy_1);
  bkg->SetBranchAddress("mgo_max_energy_2",&mgo_max_energy_2);
  bkg->SetBranchAddress("mgo_total_other_energy",&mgo_total_other_energy);
  bkg->SetBranchAddress("mgo_n_total_showers",&mgo_n_total_showers);
  bkg->SetBranchAddress("mgo_total_other_energy_1",&mgo_total_other_energy_1);
  bkg->SetBranchAddress("mgo_flag",&mgo_flag);

  bkg->SetBranchAddress("mgt_flag_single_shower",&mgt_flag_single_shower);
  bkg->SetBranchAddress("mgt_max_energy",&mgt_max_energy);
  bkg->SetBranchAddress("mgt_total_other_energy",&mgt_total_other_energy);
  bkg->SetBranchAddress("mgt_max_energy_1",&mgt_max_energy_1);
  bkg->SetBranchAddress("mgt_e_indirect_max_energy",&mgt_e_indirect_max_energy);
  bkg->SetBranchAddress("mgt_e_direct_max_energy",&mgt_e_direct_max_energy);
  bkg->SetBranchAddress("mgt_n_direct_showers",&mgt_n_direct_showers);
  bkg->SetBranchAddress("mgt_e_direct_total_energy",&mgt_e_direct_total_energy);
  bkg->SetBranchAddress("mgt_flag_indirect_max_pio",&mgt_flag_indirect_max_pio);
  bkg->SetBranchAddress("mgt_e_indirect_total_energy",&mgt_e_indirect_total_energy);
  bkg->SetBranchAddress("mgt_flag",&mgt_flag);

   
  
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
  
  
  Tsig->Branch("mgo_energy",&mgo_energy,"mgo_energy/F");
  Tsig->Branch("mgo_max_energy",&mgo_max_energy,"mgo_max_energy/F");
  Tsig->Branch("mgo_total_energy",&mgo_total_energy,"mgo_total_energy/F");
  Tsig->Branch("mgo_n_showers",&mgo_n_showers,"mgo_n_showers/F");
  Tsig->Branch("mgo_max_energy_1",&mgo_max_energy_1,"mgo_max_energy_1/F");
  Tsig->Branch("mgo_max_energy_2",&mgo_max_energy_2,"mgo_max_energy_2/F");
  Tsig->Branch("mgo_total_other_energy",&mgo_total_other_energy,"mgo_total_other_energy/F");
  Tsig->Branch("mgo_n_total_showers",&mgo_n_total_showers,"mgo_n_total_showers/F");
  Tsig->Branch("mgo_total_other_energy_1",&mgo_total_other_energy_1,"mgo_total_other_energy_1/F");
  Tsig->Branch("mgo_flag",&mgo_flag,"mgo_flag/F");
  
  Tsig->Branch("mgt_flag_single_shower",&mgt_flag_single_shower,"mgt_flag_single_shower/F");
  Tsig->Branch("mgt_max_energy",&mgt_max_energy,"mgt_max_energy/F");
  Tsig->Branch("mgt_total_other_energy",&mgt_total_other_energy,"mgt_total_other_energy/F");
  Tsig->Branch("mgt_max_energy_1",&mgt_max_energy_1,"mgt_max_energy_1/F");
  Tsig->Branch("mgt_e_indirect_max_energy",&mgt_e_indirect_max_energy,"mgt_e_indirect_max_energy/F");
  Tsig->Branch("mgt_e_direct_max_energy",&mgt_e_direct_max_energy,"mgt_e_direct_max_energy/F");
  Tsig->Branch("mgt_n_direct_showers",&mgt_n_direct_showers,"mgt_n_direct_showers/F");
  Tsig->Branch("mgt_e_direct_total_energy",&mgt_e_direct_total_energy,"mgt_e_direct_total_energy/F");
  Tsig->Branch("mgt_flag_indirect_max_pio",&mgt_flag_indirect_max_pio,"mgt_flag_indirect_max_pio/F");
  Tsig->Branch("mgt_e_indirect_total_energy",&mgt_e_indirect_total_energy,"mgt_e_indirect_total_energy/F");
  Tsig->Branch("mgt_flag",&mgt_flag,"mgt_flag/F");

  
  Tbkg->Branch("mgo_energy",&mgo_energy,"mgo_energy/F");
  Tbkg->Branch("mgo_max_energy",&mgo_max_energy,"mgo_max_energy/F");
  Tbkg->Branch("mgo_total_energy",&mgo_total_energy,"mgo_total_energy/F");
  Tbkg->Branch("mgo_n_showers",&mgo_n_showers,"mgo_n_showers/F");
  Tbkg->Branch("mgo_max_energy_1",&mgo_max_energy_1,"mgo_max_energy_1/F");
  Tbkg->Branch("mgo_max_energy_2",&mgo_max_energy_2,"mgo_max_energy_2/F");
  Tbkg->Branch("mgo_total_other_energy",&mgo_total_other_energy,"mgo_total_other_energy/F");
  Tbkg->Branch("mgo_n_total_showers",&mgo_n_total_showers,"mgo_n_total_showers/F");
  Tbkg->Branch("mgo_total_other_energy_1",&mgo_total_other_energy_1,"mgo_total_other_energy_1/F");
  Tbkg->Branch("mgo_flag",&mgo_flag,"mgo_flag/F");
  
  Tbkg->Branch("mgt_flag_single_shower",&mgt_flag_single_shower,"mgt_flag_single_shower/F");
  Tbkg->Branch("mgt_max_energy",&mgt_max_energy,"mgt_max_energy/F");
  Tbkg->Branch("mgt_total_other_energy",&mgt_total_other_energy,"mgt_total_other_energy/F");
  Tbkg->Branch("mgt_max_energy_1",&mgt_max_energy_1,"mgt_max_energy_1/F");
  Tbkg->Branch("mgt_e_indirect_max_energy",&mgt_e_indirect_max_energy,"mgt_e_indirect_max_energy/F");
  Tbkg->Branch("mgt_e_direct_max_energy",&mgt_e_direct_max_energy,"mgt_e_direct_max_energy/F");
  Tbkg->Branch("mgt_n_direct_showers",&mgt_n_direct_showers,"mgt_n_direct_showers/F");
  Tbkg->Branch("mgt_e_direct_total_energy",&mgt_e_direct_total_energy,"mgt_e_direct_total_energy/F");
  Tbkg->Branch("mgt_flag_indirect_max_pio",&mgt_flag_indirect_max_pio,"mgt_flag_indirect_max_pio/F");
  Tbkg->Branch("mgt_e_indirect_total_energy",&mgt_e_indirect_total_energy,"mgt_e_indirect_total_energy/F");
  Tbkg->Branch("mgt_flag",&mgt_flag,"mgt_flag/F");
  
  
  Tsig->Branch("mgo_mgt_bdt",&bdt_value,"data/F");
  Tbkg->Branch("mgo_mgt_bdt",&bdt_value,"data/F");
  
  TMVA::Reader *reader = new TMVA::Reader();
  reader->AddVariable("mgo_energy",&mgo_energy);
  reader->AddVariable("mgo_max_energy",&mgo_max_energy);
  reader->AddVariable("mgo_total_energy",&mgo_total_energy);
  reader->AddVariable("mgo_n_showers",&mgo_n_showers);
  reader->AddVariable("mgo_max_energy_1",&mgo_max_energy_1);
  reader->AddVariable("mgo_max_energy_2",&mgo_max_energy_2);
  reader->AddVariable("mgo_total_other_energy",&mgo_total_other_energy);
  reader->AddVariable("mgo_n_total_showers",&mgo_n_total_showers);
  reader->AddVariable("mgo_total_other_energy_1",&mgo_total_other_energy_1);
  
  reader->AddVariable("mgt_flag_single_shower",&mgt_flag_single_shower);
  reader->AddVariable("mgt_max_energy",&mgt_max_energy);
  reader->AddVariable("mgt_total_other_energy",&mgt_total_other_energy);
  reader->AddVariable("mgt_max_energy_1",&mgt_max_energy_1);
  reader->AddVariable("mgt_e_indirect_max_energy",&mgt_e_indirect_max_energy);
  reader->AddVariable("mgt_e_direct_max_energy",&mgt_e_direct_max_energy);
  reader->AddVariable("mgt_n_direct_showers",&mgt_n_direct_showers);
  reader->AddVariable("mgt_e_direct_total_energy",&mgt_e_direct_total_energy);
  reader->AddVariable("mgt_flag_indirect_max_pio",&mgt_flag_indirect_max_pio);
  reader->AddVariable("mgt_e_indirect_total_energy",&mgt_e_indirect_total_energy);
  
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

