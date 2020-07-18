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
void TestEvaluate_r1();


void Run_r2();
void InitBDT_r2();
void TestEvaluate_r2();


void convert_file();


void convert_file(){
  TFile *file = new TFile("bdtfile_0718.root");
  TTree *sig = (TTree*)file->Get("sig");
  TTree *bkg = (TTree*)file->Get("bkg");

  int br_filled;
  int mip_quality_flag;
  double mip_quality_energy;
  int mip_quality_overlap;
  int mip_quality_n_showers;
  int mip_quality_n_tracks;
  int mip_quality_flag_inside_pi0;
  int mip_quality_n_pi0_showers;
  double mip_quality_shortest_length;
  double mip_quality_acc_length;
  double mip_quality_shortest_angle;
  int mip_quality_flag_proton;
  int mip_quality_filled;

  float trueEdep;
  float weight;
  float lowEweight;
  Int_t nueTag;
  
  sig->SetBranchAddress("trueEdep",&trueEdep);
  sig->SetBranchAddress("weight",&weight);
  sig->SetBranchAddress("lowEweight",&lowEweight);
  sig->SetBranchAddress("nueTag",&nueTag);
  
    
  sig->SetBranchAddress("br_filled",&br_filled);
  sig->SetBranchAddress("mip_quality_flag",&mip_quality_flag);
  sig->SetBranchAddress("mip_quality_energy",&mip_quality_energy);
  sig->SetBranchAddress("mip_quality_overlap",&mip_quality_overlap);
  sig->SetBranchAddress("mip_quality_n_showers",&mip_quality_n_showers);
  sig->SetBranchAddress("mip_quality_n_tracks",&mip_quality_n_tracks);
  sig->SetBranchAddress("mip_quality_flag_inside_pi0",&mip_quality_flag_inside_pi0);
  sig->SetBranchAddress("mip_quality_n_pi0_showers",&mip_quality_n_pi0_showers);
  sig->SetBranchAddress("mip_quality_shortest_length",&mip_quality_shortest_length);
  sig->SetBranchAddress("mip_quality_acc_length",&mip_quality_acc_length);
  sig->SetBranchAddress("mip_quality_shortest_angle",&mip_quality_shortest_angle);
  sig->SetBranchAddress("mip_quality_flag_proton",&mip_quality_flag_proton);
  sig->SetBranchAddress("mip_quality_filled",&mip_quality_filled);
 

  bkg->SetBranchAddress("trueEdep",&trueEdep);
  bkg->SetBranchAddress("weight",&weight);
  bkg->SetBranchAddress("lowEweight",&lowEweight);
  bkg->SetBranchAddress("nueTag",&nueTag);
  
  bkg->SetBranchAddress("br_filled",&br_filled);
  bkg->SetBranchAddress("mip_quality_flag",&mip_quality_flag);
  bkg->SetBranchAddress("mip_quality_energy",&mip_quality_energy);
  bkg->SetBranchAddress("mip_quality_overlap",&mip_quality_overlap);
  bkg->SetBranchAddress("mip_quality_n_showers",&mip_quality_n_showers);
  bkg->SetBranchAddress("mip_quality_n_tracks",&mip_quality_n_tracks);
  bkg->SetBranchAddress("mip_quality_flag_inside_pi0",&mip_quality_flag_inside_pi0);
  bkg->SetBranchAddress("mip_quality_n_pi0_showers",&mip_quality_n_pi0_showers);
  bkg->SetBranchAddress("mip_quality_shortest_length",&mip_quality_shortest_length);
  bkg->SetBranchAddress("mip_quality_acc_length",&mip_quality_acc_length);
  bkg->SetBranchAddress("mip_quality_shortest_angle",&mip_quality_shortest_angle);
  bkg->SetBranchAddress("mip_quality_flag_proton",&mip_quality_flag_proton);
  bkg->SetBranchAddress("mip_quality_filled",&mip_quality_filled);

  TFile *new_file = new TFile("reduced.root","RECREATE");
  TTree *Tsig = new TTree("sig","sig");
  TTree *Tbkg = new TTree("bkg","bkg");
  Tsig->SetDirectory(new_file);
  Tbkg->SetDirectory(new_file);
  
  float mip_quality_flag_f;
  float mip_quality_energy_f;
  float mip_quality_overlap_f;
  float mip_quality_n_showers_f;
  float mip_quality_n_tracks_f;
  float mip_quality_flag_inside_pi0_f;
  float mip_quality_n_pi0_showers_f;
  float mip_quality_shortest_length_f;
  float mip_quality_acc_length_f;
  float mip_quality_shortest_angle_f;
  float mip_quality_flag_proton_f;
  float mip_quality_filled_f;

  Tsig->Branch("mip_quality_flag",&mip_quality_flag_f,"data/F");
  Tsig->Branch("mip_quality_energy",&mip_quality_energy_f,"data/F");
  Tsig->Branch("mip_quality_overlap",&mip_quality_overlap_f,"data/F");
  Tsig->Branch("mip_quality_n_showers",&mip_quality_n_showers_f,"data/F");
  Tsig->Branch("mip_quality_n_tracks",&mip_quality_n_tracks_f,"data/F");
  Tsig->Branch("mip_quality_flag_inside_pi0",&mip_quality_flag_inside_pi0_f,"data/F");
  Tsig->Branch("mip_quality_n_pi0_showers",&mip_quality_n_pi0_showers_f,"data/F");
  Tsig->Branch("mip_quality_shortest_length",&mip_quality_shortest_length_f,"data/F");
  Tsig->Branch("mip_quality_acc_length",&mip_quality_acc_length_f,"data/F");
  Tsig->Branch("mip_quality_shortest_angle",&mip_quality_shortest_angle_f,"data/F");
  Tsig->Branch("mip_quality_flag_proton",&mip_quality_flag_proton_f,"data/F");
  Tsig->Branch("mip_quality_filled",&mip_quality_filled_f,"data/F");

  Tsig->Branch("trueEdep",&trueEdep,"data/F");
  Tsig->Branch("weight",&weight,"data/F");
  Tsig->Branch("lowEweight",&lowEweight,"data/F");
  Tsig->Branch("nueTag",&nueTag,"data/I");
  
  Tbkg->Branch("mip_quality_flag",&mip_quality_flag_f,"data/F");
  Tbkg->Branch("mip_quality_energy",&mip_quality_energy_f,"data/F");
  Tbkg->Branch("mip_quality_overlap",&mip_quality_overlap_f,"data/F");
  Tbkg->Branch("mip_quality_n_showers",&mip_quality_n_showers_f,"data/F");
  Tbkg->Branch("mip_quality_n_tracks",&mip_quality_n_tracks_f,"data/F");
  Tbkg->Branch("mip_quality_flag_inside_pi0",&mip_quality_flag_inside_pi0_f,"data/F");
  Tbkg->Branch("mip_quality_n_pi0_showers",&mip_quality_n_pi0_showers_f,"data/F");
  Tbkg->Branch("mip_quality_shortest_length",&mip_quality_shortest_length_f,"data/F");
  Tbkg->Branch("mip_quality_acc_length",&mip_quality_acc_length_f,"data/F");
  Tbkg->Branch("mip_quality_shortest_angle",&mip_quality_shortest_angle_f,"data/F");
  Tbkg->Branch("mip_quality_flag_proton",&mip_quality_flag_proton_f,"data/F");
  Tbkg->Branch("mip_quality_filled",&mip_quality_filled_f,"data/F");

  Tbkg->Branch("trueEdep",&trueEdep,"data/F");
  Tbkg->Branch("weight",&weight,"data/F");
  Tbkg->Branch("lowEweight",&lowEweight,"data/F");
  Tbkg->Branch("nueTag",&nueTag,"data/I");
  
  for (Int_t i=0;i!=sig->GetEntries();i++){
    sig->GetEntry(i);

    mip_quality_flag_f = mip_quality_flag;
    mip_quality_energy_f = mip_quality_energy;
    mip_quality_overlap_f= mip_quality_overlap;
    mip_quality_n_showers_f= mip_quality_n_showers;
    mip_quality_n_tracks_f= mip_quality_n_tracks;
    mip_quality_flag_inside_pi0_f= mip_quality_flag_inside_pi0;
    mip_quality_n_pi0_showers_f= mip_quality_n_pi0_showers;
    mip_quality_shortest_length_f= mip_quality_shortest_length;
    mip_quality_acc_length_f= mip_quality_acc_length;
    mip_quality_shortest_angle_f= mip_quality_shortest_angle;
    mip_quality_flag_proton_f= mip_quality_flag_proton;
    mip_quality_filled_f= mip_quality_filled;
    

    // fixed incorrect filled flag ...
    mip_quality_filled_f = br_filled;
    
    // fix NaN problem in this variable 
    if (std::isnan(mip_quality_shortest_angle_f)) mip_quality_shortest_angle_f = 0;

    // limit the highest value ...
    if (mip_quality_shortest_length_f > 1000) mip_quality_shortest_length_f = 1000;
    
    
    Tsig->Fill();
  }

   for (Int_t i=0;i!=bkg->GetEntries();i++){
    bkg->GetEntry(i);

    mip_quality_flag_f = mip_quality_flag;
    mip_quality_energy_f = mip_quality_energy;
    mip_quality_overlap_f= mip_quality_overlap;
    mip_quality_n_showers_f= mip_quality_n_showers;
    mip_quality_n_tracks_f= mip_quality_n_tracks;
    mip_quality_flag_inside_pi0_f= mip_quality_flag_inside_pi0;
    mip_quality_n_pi0_showers_f= mip_quality_n_pi0_showers;
    mip_quality_shortest_length_f= mip_quality_shortest_length;
    mip_quality_acc_length_f= mip_quality_acc_length;
    mip_quality_shortest_angle_f= mip_quality_shortest_angle;
    mip_quality_flag_proton_f= mip_quality_flag_proton;
    mip_quality_filled_f= mip_quality_filled;
    

    // fixed incorrect filled flag ...
    mip_quality_filled_f = br_filled;
    
    // fix NaN problem in this variable 
    if (std::isnan(mip_quality_shortest_angle_f)) mip_quality_shortest_angle_f = 0;

    // limit the highest value ...
    if (mip_quality_shortest_length_f > 1000) mip_quality_shortest_length_f = 1000;
    
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

    TestEvaluate_r2();

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

    TestEvaluate_r1();

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
    // dataloader->AddVariable( "sum := var1+var2", "v1+v2", "", 'F' );
    // dataloader->AddVariable( "diff := var1-var2", "v1-v2", "", 'F' );
    // dataloader->AddVariable( "var3", "v3", "m", 'F' );
    // dataloader->AddVariable( "var4", "v4", "MeV", 'F' );

    dataloader->AddVariable("mip_quality_energy","mip_quality_energy","MeV",'D');
    dataloader->AddVariable("mip_quality_overlap","mip_quality_overlap","",'I');
    dataloader->AddVariable("mip_quality_n_showers","mip_quality_n_showers","",'I');
    dataloader->AddVariable("mip_quality_n_tracks","mip_quality_n_tracks","",'I');
    dataloader->AddVariable("mip_quality_flag_inside_pi0","mip_quality_flag_inside_pi0","",'I');
    dataloader->AddVariable("mip_quality_n_pi0_showers","mip_quality_n_pi0_showers","",'I');
    dataloader->AddVariable("mip_quality_shortest_length","mip_quality_shortest_length","cm",'D');
    dataloader->AddVariable("mip_quality_acc_length","mip_quality_acc_length","cm",'D');
    dataloader->AddVariable("mip_quality_shortest_angle","mip_quality_shortest_angle","",'D');
    dataloader->AddVariable("mip_quality_flag_proton","mip_quality_flag_proton","",'I');
    
    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight" );
    dataloader->SetBackgroundWeightExpression( "weight" );

    // Apply additional cuts on the signal and background samples (can be different)
    TCut mycut_s = "mip_quality_filled==1"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
    TCut mycut_b = "mip_quality_filled==1 && mip_quality_flag==0"; // for example: TCut mycutb = "abs(var1)<0.5";
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
        "nTrain_Signal=20000:"
        "nTrain_Background=851:"
	"nTest_Signal=5000:"
        "nTest_Background=130:"
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

    dataloader->AddVariable("mip_quality_energy","mip_quality_energy","MeV",'D');
    dataloader->AddVariable("mip_quality_overlap","mip_quality_overlap","",'I');
    dataloader->AddVariable("mip_quality_n_showers","mip_quality_n_showers","",'I');
    dataloader->AddVariable("mip_quality_n_tracks","mip_quality_n_tracks","",'I');
    dataloader->AddVariable("mip_quality_flag_inside_pi0","mip_quality_flag_inside_pi0","",'I');
    dataloader->AddVariable("mip_quality_n_pi0_showers","mip_quality_n_pi0_showers","",'I');
    dataloader->AddVariable("mip_quality_shortest_length","mip_quality_shortest_length","cm",'D');
    dataloader->AddVariable("mip_quality_acc_length","mip_quality_acc_length","cm",'D');
    dataloader->AddVariable("mip_quality_shortest_angle","mip_quality_shortest_angle","",'D');
    dataloader->AddVariable("mip_quality_flag_proton","mip_quality_flag_proton","",'I');
    
    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight" );
    dataloader->SetBackgroundWeightExpression( "weight" );

    // Apply additional cuts on the signal and background samples (can be different)
    TCut mycut_s = "mip_quality_filled==1"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
    TCut mycut_b = "mip_quality_filled==1 && (mip_quality_flag==0 || mip_quality_bdt < 0.1)"; // for example: TCut mycutb = "abs(var1)<0.5";
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
        "nTrain_Signal=20000:"
        "nTrain_Background=1200:"
	"nTest_Signal=5000:"
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


void TestEvaluate_r1()
{
  Float_t mip_quality_energy = 200;
  Float_t mip_quality_overlap = 0;
  Float_t mip_quality_n_showers = 0;
  Float_t mip_quality_n_tracks = 0;
  Float_t mip_quality_flag_inside_pi0= 0;
  Float_t mip_quality_n_pi0_showers = 0;
  Float_t mip_quality_shortest_length = 0;
  Float_t mip_quality_shortest_angle = 0;
  Float_t mip_quality_acc_length = 0;
  Float_t mip_quality_flag_proton = 0;
  Float_t bdt_value = 0;
  Float_t mip_quality_filled;
  Float_t mip_quality_flag;
  
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
  

  
  sig->SetBranchAddress("mip_quality_flag",&mip_quality_flag);
  sig->SetBranchAddress("mip_quality_energy",&mip_quality_energy);
  sig->SetBranchAddress("mip_quality_overlap",&mip_quality_overlap);
  sig->SetBranchAddress("mip_quality_n_showers",&mip_quality_n_showers);
  sig->SetBranchAddress("mip_quality_n_tracks",&mip_quality_n_tracks);
  sig->SetBranchAddress("mip_quality_flag_inside_pi0",&mip_quality_flag_inside_pi0);
  sig->SetBranchAddress("mip_quality_n_pi0_showers",&mip_quality_n_pi0_showers);
  sig->SetBranchAddress("mip_quality_shortest_length",&mip_quality_shortest_length);
  sig->SetBranchAddress("mip_quality_acc_length",&mip_quality_acc_length);
  sig->SetBranchAddress("mip_quality_shortest_angle",&mip_quality_shortest_angle);
  sig->SetBranchAddress("mip_quality_flag_proton",&mip_quality_flag_proton);
  sig->SetBranchAddress("mip_quality_filled",&mip_quality_filled);

  bkg->SetBranchAddress("mip_quality_flag",&mip_quality_flag);
  bkg->SetBranchAddress("mip_quality_energy",&mip_quality_energy);
  bkg->SetBranchAddress("mip_quality_overlap",&mip_quality_overlap);
  bkg->SetBranchAddress("mip_quality_n_showers",&mip_quality_n_showers);
  bkg->SetBranchAddress("mip_quality_n_tracks",&mip_quality_n_tracks);
  bkg->SetBranchAddress("mip_quality_flag_inside_pi0",&mip_quality_flag_inside_pi0);
  bkg->SetBranchAddress("mip_quality_n_pi0_showers",&mip_quality_n_pi0_showers);
  bkg->SetBranchAddress("mip_quality_shortest_length",&mip_quality_shortest_length);
  bkg->SetBranchAddress("mip_quality_acc_length",&mip_quality_acc_length);
  bkg->SetBranchAddress("mip_quality_shortest_angle",&mip_quality_shortest_angle);
  bkg->SetBranchAddress("mip_quality_flag_proton",&mip_quality_flag_proton);
  bkg->SetBranchAddress("mip_quality_filled",&mip_quality_filled);

  bkg->SetBranchAddress("trueEdep",&trueEdep);
  bkg->SetBranchAddress("weight",&weight);
  bkg->SetBranchAddress("lowEweight",&lowEweight);
  bkg->SetBranchAddress("nueTag",&nueTag);

  
  TFile *new_file = new TFile("round_1.root","RECREATE");
  TTree *Tsig = new TTree("sig","sig");
  TTree *Tbkg = new TTree("bkg","bkg");
  Tsig->SetDirectory(new_file);
  Tbkg->SetDirectory(new_file);
  
  Tsig->Branch("mip_quality_flag",&mip_quality_flag,"data/F");
  Tsig->Branch("mip_quality_energy",&mip_quality_energy,"data/F");
  Tsig->Branch("mip_quality_overlap",&mip_quality_overlap,"data/F");
  Tsig->Branch("mip_quality_n_showers",&mip_quality_n_showers,"data/F");
  Tsig->Branch("mip_quality_n_tracks",&mip_quality_n_tracks,"data/F");
  Tsig->Branch("mip_quality_flag_inside_pi0",&mip_quality_flag_inside_pi0,"data/F");
  Tsig->Branch("mip_quality_n_pi0_showers",&mip_quality_n_pi0_showers,"data/F");
  Tsig->Branch("mip_quality_shortest_length",&mip_quality_shortest_length,"data/F");
  Tsig->Branch("mip_quality_acc_length",&mip_quality_acc_length,"data/F");
  Tsig->Branch("mip_quality_shortest_angle",&mip_quality_shortest_angle,"data/F");
  Tsig->Branch("mip_quality_flag_proton",&mip_quality_flag_proton,"data/F");
  Tsig->Branch("mip_quality_filled",&mip_quality_filled,"data/F");
  Tsig->Branch("mip_quality_bdt",&bdt_value,"data/F");

  Tsig->Branch("trueEdep",&trueEdep,"data/F");
  Tsig->Branch("weight",&weight,"data/F");
  Tsig->Branch("lowEweight",&lowEweight,"data/F");
  Tsig->Branch("nueTag",&nueTag,"data/I");

  
  Tbkg->Branch("mip_quality_flag",&mip_quality_flag,"data/F");
  Tbkg->Branch("mip_quality_energy",&mip_quality_energy,"data/F");
  Tbkg->Branch("mip_quality_overlap",&mip_quality_overlap,"data/F");
  Tbkg->Branch("mip_quality_n_showers",&mip_quality_n_showers,"data/F");
  Tbkg->Branch("mip_quality_n_tracks",&mip_quality_n_tracks,"data/F");
  Tbkg->Branch("mip_quality_flag_inside_pi0",&mip_quality_flag_inside_pi0,"data/F");
  Tbkg->Branch("mip_quality_n_pi0_showers",&mip_quality_n_pi0_showers,"data/F");
  Tbkg->Branch("mip_quality_shortest_length",&mip_quality_shortest_length,"data/F");
  Tbkg->Branch("mip_quality_acc_length",&mip_quality_acc_length,"data/F");
  Tbkg->Branch("mip_quality_shortest_angle",&mip_quality_shortest_angle,"data/F");
  Tbkg->Branch("mip_quality_flag_proton",&mip_quality_flag_proton,"data/F");
  Tbkg->Branch("mip_quality_filled",&mip_quality_filled,"data/F");
  Tbkg->Branch("mip_quality_bdt",&bdt_value,"data/F");

  Tbkg->Branch("trueEdep",&trueEdep,"data/F");
  Tbkg->Branch("weight",&weight,"data/F");
  Tbkg->Branch("lowEweight",&lowEweight,"data/F");
  Tbkg->Branch("nueTag",&nueTag,"data/I");
  
  
  TMVA::Reader *reader = new TMVA::Reader();
  reader->AddVariable("mip_quality_energy",&mip_quality_energy);
  reader->AddVariable("mip_quality_overlap",&mip_quality_overlap);
  reader->AddVariable("mip_quality_n_showers",&mip_quality_n_showers);
  reader->AddVariable("mip_quality_n_tracks",&mip_quality_n_tracks);
  reader->AddVariable("mip_quality_flag_inside_pi0",&mip_quality_flag_inside_pi0);
  reader->AddVariable("mip_quality_n_pi0_showers",&mip_quality_n_pi0_showers);
  reader->AddVariable("mip_quality_shortest_length",&mip_quality_shortest_length);
  reader->AddVariable("mip_quality_acc_length",&mip_quality_acc_length);
  reader->AddVariable("mip_quality_shortest_angle",&mip_quality_shortest_angle);
  reader->AddVariable("mip_quality_flag_proton",&mip_quality_flag_proton);
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



void TestEvaluate_r2()
{
  Float_t mip_quality_energy = 200;
  Float_t mip_quality_overlap = 0;
  Float_t mip_quality_n_showers = 0;
  Float_t mip_quality_n_tracks = 0;
  Float_t mip_quality_flag_inside_pi0= 0;
  Float_t mip_quality_n_pi0_showers = 0;
  Float_t mip_quality_shortest_length = 0;
  Float_t mip_quality_shortest_angle = 0;
  Float_t mip_quality_acc_length = 0;
  Float_t mip_quality_flag_proton = 0;
  Float_t bdt_value = 0;
  Float_t mip_quality_filled;
  Float_t mip_quality_flag;
  
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
  
  sig->SetBranchAddress("mip_quality_flag",&mip_quality_flag);
  sig->SetBranchAddress("mip_quality_energy",&mip_quality_energy);
  sig->SetBranchAddress("mip_quality_overlap",&mip_quality_overlap);
  sig->SetBranchAddress("mip_quality_n_showers",&mip_quality_n_showers);
  sig->SetBranchAddress("mip_quality_n_tracks",&mip_quality_n_tracks);
  sig->SetBranchAddress("mip_quality_flag_inside_pi0",&mip_quality_flag_inside_pi0);
  sig->SetBranchAddress("mip_quality_n_pi0_showers",&mip_quality_n_pi0_showers);
  sig->SetBranchAddress("mip_quality_shortest_length",&mip_quality_shortest_length);
  sig->SetBranchAddress("mip_quality_acc_length",&mip_quality_acc_length);
  sig->SetBranchAddress("mip_quality_shortest_angle",&mip_quality_shortest_angle);
  sig->SetBranchAddress("mip_quality_flag_proton",&mip_quality_flag_proton);
  sig->SetBranchAddress("mip_quality_filled",&mip_quality_filled);

  bkg->SetBranchAddress("mip_quality_flag",&mip_quality_flag);
  bkg->SetBranchAddress("mip_quality_energy",&mip_quality_energy);
  bkg->SetBranchAddress("mip_quality_overlap",&mip_quality_overlap);
  bkg->SetBranchAddress("mip_quality_n_showers",&mip_quality_n_showers);
  bkg->SetBranchAddress("mip_quality_n_tracks",&mip_quality_n_tracks);
  bkg->SetBranchAddress("mip_quality_flag_inside_pi0",&mip_quality_flag_inside_pi0);
  bkg->SetBranchAddress("mip_quality_n_pi0_showers",&mip_quality_n_pi0_showers);
  bkg->SetBranchAddress("mip_quality_shortest_length",&mip_quality_shortest_length);
  bkg->SetBranchAddress("mip_quality_acc_length",&mip_quality_acc_length);
  bkg->SetBranchAddress("mip_quality_shortest_angle",&mip_quality_shortest_angle);
  bkg->SetBranchAddress("mip_quality_flag_proton",&mip_quality_flag_proton);
  bkg->SetBranchAddress("mip_quality_filled",&mip_quality_filled);

  bkg->SetBranchAddress("trueEdep",&trueEdep);
  bkg->SetBranchAddress("weight",&weight);
  bkg->SetBranchAddress("lowEweight",&lowEweight);
  bkg->SetBranchAddress("nueTag",&nueTag);

  
  TFile *new_file = new TFile("round_2.root","RECREATE");
  TTree *Tsig = new TTree("sig","sig");
  TTree *Tbkg = new TTree("bkg","bkg");
  Tsig->SetDirectory(new_file);
  Tbkg->SetDirectory(new_file);
  
  Tsig->Branch("mip_quality_flag",&mip_quality_flag,"data/F");
  Tsig->Branch("mip_quality_energy",&mip_quality_energy,"data/F");
  Tsig->Branch("mip_quality_overlap",&mip_quality_overlap,"data/F");
  Tsig->Branch("mip_quality_n_showers",&mip_quality_n_showers,"data/F");
  Tsig->Branch("mip_quality_n_tracks",&mip_quality_n_tracks,"data/F");
  Tsig->Branch("mip_quality_flag_inside_pi0",&mip_quality_flag_inside_pi0,"data/F");
  Tsig->Branch("mip_quality_n_pi0_showers",&mip_quality_n_pi0_showers,"data/F");
  Tsig->Branch("mip_quality_shortest_length",&mip_quality_shortest_length,"data/F");
  Tsig->Branch("mip_quality_acc_length",&mip_quality_acc_length,"data/F");
  Tsig->Branch("mip_quality_shortest_angle",&mip_quality_shortest_angle,"data/F");
  Tsig->Branch("mip_quality_flag_proton",&mip_quality_flag_proton,"data/F");
  Tsig->Branch("mip_quality_filled",&mip_quality_filled,"data/F");
  Tsig->Branch("mip_quality_bdt",&bdt_value,"data/F");

  Tsig->Branch("trueEdep",&trueEdep,"data/F");
  Tsig->Branch("weight",&weight,"data/F");
  Tsig->Branch("lowEweight",&lowEweight,"data/F");
  Tsig->Branch("nueTag",&nueTag,"data/I");

  
  Tbkg->Branch("mip_quality_flag",&mip_quality_flag,"data/F");
  Tbkg->Branch("mip_quality_energy",&mip_quality_energy,"data/F");
  Tbkg->Branch("mip_quality_overlap",&mip_quality_overlap,"data/F");
  Tbkg->Branch("mip_quality_n_showers",&mip_quality_n_showers,"data/F");
  Tbkg->Branch("mip_quality_n_tracks",&mip_quality_n_tracks,"data/F");
  Tbkg->Branch("mip_quality_flag_inside_pi0",&mip_quality_flag_inside_pi0,"data/F");
  Tbkg->Branch("mip_quality_n_pi0_showers",&mip_quality_n_pi0_showers,"data/F");
  Tbkg->Branch("mip_quality_shortest_length",&mip_quality_shortest_length,"data/F");
  Tbkg->Branch("mip_quality_acc_length",&mip_quality_acc_length,"data/F");
  Tbkg->Branch("mip_quality_shortest_angle",&mip_quality_shortest_angle,"data/F");
  Tbkg->Branch("mip_quality_flag_proton",&mip_quality_flag_proton,"data/F");
  Tbkg->Branch("mip_quality_filled",&mip_quality_filled,"data/F");
  Tbkg->Branch("mip_quality_bdt",&bdt_value,"data/F");

  Tbkg->Branch("trueEdep",&trueEdep,"data/F");
  Tbkg->Branch("weight",&weight,"data/F");
  Tbkg->Branch("lowEweight",&lowEweight,"data/F");
  Tbkg->Branch("nueTag",&nueTag,"data/I");
  
  
  TMVA::Reader *reader = new TMVA::Reader();
  reader->AddVariable("mip_quality_energy",&mip_quality_energy);
  reader->AddVariable("mip_quality_overlap",&mip_quality_overlap);
  reader->AddVariable("mip_quality_n_showers",&mip_quality_n_showers);
  reader->AddVariable("mip_quality_n_tracks",&mip_quality_n_tracks);
  reader->AddVariable("mip_quality_flag_inside_pi0",&mip_quality_flag_inside_pi0);
  reader->AddVariable("mip_quality_n_pi0_showers",&mip_quality_n_pi0_showers);
  reader->AddVariable("mip_quality_shortest_length",&mip_quality_shortest_length);
  reader->AddVariable("mip_quality_acc_length",&mip_quality_acc_length);
  reader->AddVariable("mip_quality_shortest_angle",&mip_quality_shortest_angle);
  reader->AddVariable("mip_quality_flag_proton",&mip_quality_flag_proton);
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

