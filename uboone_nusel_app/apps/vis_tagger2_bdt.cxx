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

  int vis_2_filled;
  int vis_2_n_vtx_segs;
  double vis_2_min_angle;
  int vis_2_min_weak_track;
  double vis_2_angle_beam;
  double vis_2_min_angle1;
  double vis_2_iso_angle1;
  double vis_2_min_medium_dQ_dx;
  double vis_2_min_length;
  double vis_2_sg_length;
  double vis_2_max_angle;
  int vis_2_max_weak_track;
  int vis_2_flag;

  sig->SetBranchAddress("vis_2_filled",&vis_2_filled);
  sig->SetBranchAddress("vis_2_n_vtx_segs",&vis_2_n_vtx_segs);
  sig->SetBranchAddress("vis_2_min_angle",&vis_2_min_angle);
  sig->SetBranchAddress("vis_2_min_weak_track",&vis_2_min_weak_track);
  sig->SetBranchAddress("vis_2_angle_beam",&vis_2_angle_beam);
  sig->SetBranchAddress("vis_2_min_angle1",&vis_2_min_angle1);
  sig->SetBranchAddress("vis_2_iso_angle1",&vis_2_iso_angle1);
  sig->SetBranchAddress("vis_2_min_medium_dQ_dx",&vis_2_min_medium_dQ_dx);
  sig->SetBranchAddress("vis_2_min_length",&vis_2_min_length);
  sig->SetBranchAddress("vis_2_sg_length",&vis_2_sg_length);
  sig->SetBranchAddress("vis_2_max_angle",&vis_2_max_angle);
  sig->SetBranchAddress("vis_2_max_weak_track",&vis_2_max_weak_track);
  sig->SetBranchAddress("vis_2_flab",&vis_2_flag);

  bkg->SetBranchAddress("vis_2_filled",&vis_2_filled);
  bkg->SetBranchAddress("vis_2_n_vtx_segs",&vis_2_n_vtx_segs);
  bkg->SetBranchAddress("vis_2_min_angle",&vis_2_min_angle);
  bkg->SetBranchAddress("vis_2_min_weak_track",&vis_2_min_weak_track);
  bkg->SetBranchAddress("vis_2_angle_beam",&vis_2_angle_beam);
  bkg->SetBranchAddress("vis_2_min_angle1",&vis_2_min_angle1);
  bkg->SetBranchAddress("vis_2_iso_angle1",&vis_2_iso_angle1);
  bkg->SetBranchAddress("vis_2_min_medium_dQ_dx",&vis_2_min_medium_dQ_dx);
  bkg->SetBranchAddress("vis_2_min_length",&vis_2_min_length);
  bkg->SetBranchAddress("vis_2_sg_length",&vis_2_sg_length);
  bkg->SetBranchAddress("vis_2_max_angle",&vis_2_max_angle);
  bkg->SetBranchAddress("vis_2_max_weak_track",&vis_2_max_weak_track);
  bkg->SetBranchAddress("vis_2_flab",&vis_2_flag);
  
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
  
  
  float vis_2_filled_f;
  float vis_2_n_vtx_segs_f;
  float vis_2_min_angle_f;
  float vis_2_min_weak_track_f;
  float vis_2_angle_beam_f;
  float vis_2_min_angle1_f;
  float vis_2_iso_angle1_f;
  float vis_2_min_medium_dQ_dx_f;
  float vis_2_min_length_f;
  float vis_2_sg_length_f;
  float vis_2_max_angle_f;
  float vis_2_max_weak_track_f;
  float vis_2_flag_f;

  
  Tsig->Branch("vis_2_filled",&vis_2_filled_f,"vis_2_filled/F");
  Tsig->Branch("vis_2_n_vtx_segs",&vis_2_n_vtx_segs_f,"vis_2_n_vtx_segs/F");
  Tsig->Branch("vis_2_min_angle",&vis_2_min_angle_f,"vis_2_min_angle/F");
  Tsig->Branch("vis_2_min_weak_track",&vis_2_min_weak_track_f,"vis_2_min_weak_track/F");
  Tsig->Branch("vis_2_angle_beam",&vis_2_angle_beam_f,"vis_2_angle_beam/F");
  Tsig->Branch("vis_2_min_angle1",&vis_2_min_angle1_f,"vis_2_min_angle1/F");
  Tsig->Branch("vis_2_iso_angle1",&vis_2_iso_angle1_f,"vis_2_iso_angle1/F");
  Tsig->Branch("vis_2_min_medium_dQ_dx",&vis_2_min_medium_dQ_dx_f,"vis_2_min_medium_dQ_dx/F");
  Tsig->Branch("vis_2_min_length",&vis_2_min_length_f,"vis_2_min_length/F");
  Tsig->Branch("vis_2_sg_length",&vis_2_sg_length_f,"vis_2_sg_length/F");
  Tsig->Branch("vis_2_max_angle",&vis_2_max_angle_f,"vis_2_max_angle/F");
  Tsig->Branch("vis_2_max_weak_track",&vis_2_max_weak_track_f,"vis_2_max_weak_track/F");
  Tsig->Branch("vis_2_flag",&vis_2_flag_f,"vis_2_flag/F");

  Tbkg->Branch("vis_2_filled",&vis_2_filled_f,"vis_2_filled/F");
  Tbkg->Branch("vis_2_n_vtx_segs",&vis_2_n_vtx_segs_f,"vis_2_n_vtx_segs/F");
  Tbkg->Branch("vis_2_min_angle",&vis_2_min_angle_f,"vis_2_min_angle/F");
  Tbkg->Branch("vis_2_min_weak_track",&vis_2_min_weak_track_f,"vis_2_min_weak_track/F");
  Tbkg->Branch("vis_2_angle_beam",&vis_2_angle_beam_f,"vis_2_angle_beam/F");
  Tbkg->Branch("vis_2_min_angle1",&vis_2_min_angle1_f,"vis_2_min_angle1/F");
  Tbkg->Branch("vis_2_iso_angle1",&vis_2_iso_angle1_f,"vis_2_iso_angle1/F");
  Tbkg->Branch("vis_2_min_medium_dQ_dx",&vis_2_min_medium_dQ_dx_f,"vis_2_min_medium_dQ_dx/F");
  Tbkg->Branch("vis_2_min_length",&vis_2_min_length_f,"vis_2_min_length/F");
  Tbkg->Branch("vis_2_sg_length",&vis_2_sg_length_f,"vis_2_sg_length/F");
  Tbkg->Branch("vis_2_max_angle",&vis_2_max_angle_f,"vis_2_max_angle/F");
  Tbkg->Branch("vis_2_max_weak_track",&vis_2_max_weak_track_f,"vis_2_max_weak_track/F");
  Tbkg->Branch("vis_2_flag",&vis_2_flag_f,"vis_2_flag/F");

  
  
  for (Int_t i=0;i!=sig->GetEntries();i++){
    sig->GetEntry(i);

    vis_2_filled_f = vis_2_filled;
    vis_2_n_vtx_segs_f = vis_2_n_vtx_segs;
    vis_2_min_angle_f = vis_2_min_angle;
    vis_2_min_weak_track_f = vis_2_min_weak_track;
    vis_2_angle_beam_f = vis_2_angle_beam;
    vis_2_min_angle1_f = vis_2_min_angle1;
    vis_2_iso_angle1_f = vis_2_iso_angle1;
    vis_2_min_medium_dQ_dx_f = vis_2_min_medium_dQ_dx;
    vis_2_min_length_f = vis_2_min_length;
    vis_2_sg_length_f = vis_2_sg_length;
    vis_2_max_angle_f = vis_2_max_angle; 
    vis_2_max_weak_track_f = vis_2_max_weak_track;
    vis_2_flag_f = vis_2_flag;
    
    
    Tsig->Fill();
  }

   for (Int_t i=0;i!=bkg->GetEntries();i++){
    bkg->GetEntry(i);
    
    vis_2_filled_f = vis_2_filled;
    vis_2_n_vtx_segs_f = vis_2_n_vtx_segs;
    vis_2_min_angle_f = vis_2_min_angle;
    vis_2_min_weak_track_f = vis_2_min_weak_track;
    vis_2_angle_beam_f = vis_2_angle_beam;
    vis_2_min_angle1_f = vis_2_min_angle1;
    vis_2_iso_angle1_f = vis_2_iso_angle1;
    vis_2_min_medium_dQ_dx_f = vis_2_min_medium_dQ_dx;
    vis_2_min_length_f = vis_2_min_length;
    vis_2_sg_length_f = vis_2_sg_length;
    vis_2_max_angle_f = vis_2_max_angle; 
    vis_2_max_weak_track_f = vis_2_max_weak_track;
    vis_2_flag_f = vis_2_flag;
    
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
    
    
    dataloader->AddVariable("vis_2_n_vtx_segs","vis_2_n_vtx_segs","",'F');
    dataloader->AddVariable("vis_2_min_angle","vis_2_min_angle","deg",'F');
    dataloader->AddVariable("vis_2_min_weak_track","vis_2_min_weak_track","deg",'F');
    dataloader->AddVariable("vis_2_angle_beam","vis_2_angle_beam","deg",'F');
    dataloader->AddVariable("vis_2_min_angle1","vis_2_min_angle1","deg",'F');
    dataloader->AddVariable("vis_2_iso_angle1","vis_2_iso_angle1","MeV/cm",'F');
    dataloader->AddVariable("vis_2_min_medium_dQ_dx","vis_2_min_medium_dQ_dx","cm",'F');
    dataloader->AddVariable("vis_2_min_length","vis_2_min_length","cm",'F');
    dataloader->AddVariable("vis_2_sg_length","vis_2_sg_length","cm",'F');
    dataloader->AddVariable("vis_2_max_angle","vis_2_max_angle","deg",'F');
    dataloader->AddVariable("vis_2_max_weak_track","vis_2_max_weak_track","",'F');
    

    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "vis_2_filled==1"; // 33/44914
    TCut mycut_b = "vis_2_filled==1 && (vis_2_flag==0)"; // 420/16639
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=20000:"
					    "nTrain_Background=380:"
					    "nTest_Signal=10000:"
					    "nTest_Background=40:"
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
    
    
    dataloader->AddVariable("vis_2_n_vtx_segs","vis_2_n_vtx_segs","",'F');
    dataloader->AddVariable("vis_2_min_angle","vis_2_min_angle","deg",'F');
    dataloader->AddVariable("vis_2_min_weak_track","vis_2_min_weak_track","deg",'F');
    dataloader->AddVariable("vis_2_angle_beam","vis_2_angle_beam","deg",'F');
    dataloader->AddVariable("vis_2_min_angle1","vis_2_min_angle1","deg",'F');
    dataloader->AddVariable("vis_2_iso_angle1","vis_2_iso_angle1","MeV/cm",'F');
    dataloader->AddVariable("vis_2_min_medium_dQ_dx","vis_2_min_medium_dQ_dx","cm",'F');
    dataloader->AddVariable("vis_2_min_length","vis_2_min_length","cm",'F');
    dataloader->AddVariable("vis_2_sg_length","vis_2_sg_length","cm",'F');
    dataloader->AddVariable("vis_2_max_angle","vis_2_max_angle","deg",'F');
    dataloader->AddVariable("vis_2_max_weak_track","vis_2_max_weak_track","",'F');
    

    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "vis_2_filled==1"; // 33/44914
    TCut mycut_b = "vis_2_filled==1 && (vis_2_flag==0 || vis_2_bdt < 0 )"; // 1219/16639
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=20000:"
					    "nTrain_Background=1000:"
					    "nTest_Signal=10000:"
					    "nTest_Background=219:"
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

  float vis_2_filled;
  float vis_2_n_vtx_segs;
  float vis_2_min_angle;
  float vis_2_min_weak_track;
  float vis_2_angle_beam;
  float vis_2_min_angle1;
  float vis_2_iso_angle1;
  float vis_2_min_medium_dQ_dx;
  float vis_2_min_length;
  float vis_2_sg_length;
  float vis_2_max_angle;
  float vis_2_max_weak_track;
  float vis_2_flag;
  
  
sig->SetBranchAddress("vis_2_filled",&vis_2_filled);
  sig->SetBranchAddress("vis_2_n_vtx_segs",&vis_2_n_vtx_segs);
  sig->SetBranchAddress("vis_2_min_angle",&vis_2_min_angle);
  sig->SetBranchAddress("vis_2_min_weak_track",&vis_2_min_weak_track);
  sig->SetBranchAddress("vis_2_angle_beam",&vis_2_angle_beam);
  sig->SetBranchAddress("vis_2_min_angle1",&vis_2_min_angle1);
  sig->SetBranchAddress("vis_2_iso_angle1",&vis_2_iso_angle1);
  sig->SetBranchAddress("vis_2_min_medium_dQ_dx",&vis_2_min_medium_dQ_dx);
  sig->SetBranchAddress("vis_2_min_length",&vis_2_min_length);
  sig->SetBranchAddress("vis_2_sg_length",&vis_2_sg_length);
  sig->SetBranchAddress("vis_2_max_angle",&vis_2_max_angle);
  sig->SetBranchAddress("vis_2_max_weak_track",&vis_2_max_weak_track);
  sig->SetBranchAddress("vis_2_flag",&vis_2_flag);

  bkg->SetBranchAddress("vis_2_filled",&vis_2_filled);
  bkg->SetBranchAddress("vis_2_n_vtx_segs",&vis_2_n_vtx_segs);
  bkg->SetBranchAddress("vis_2_min_angle",&vis_2_min_angle);
  bkg->SetBranchAddress("vis_2_min_weak_track",&vis_2_min_weak_track);
  bkg->SetBranchAddress("vis_2_angle_beam",&vis_2_angle_beam);
  bkg->SetBranchAddress("vis_2_min_angle1",&vis_2_min_angle1);
  bkg->SetBranchAddress("vis_2_iso_angle1",&vis_2_iso_angle1);
  bkg->SetBranchAddress("vis_2_min_medium_dQ_dx",&vis_2_min_medium_dQ_dx);
  bkg->SetBranchAddress("vis_2_min_length",&vis_2_min_length);
  bkg->SetBranchAddress("vis_2_sg_length",&vis_2_sg_length);
  bkg->SetBranchAddress("vis_2_max_angle",&vis_2_max_angle);
  bkg->SetBranchAddress("vis_2_max_weak_track",&vis_2_max_weak_track);
  bkg->SetBranchAddress("vis_2_flag",&vis_2_flag);
  
  
  
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
  

  Tsig->Branch("vis_2_filled",&vis_2_filled,"vis_2_filled/F");
  Tsig->Branch("vis_2_n_vtx_segs",&vis_2_n_vtx_segs,"vis_2_n_vtx_segs/F");
  Tsig->Branch("vis_2_min_angle",&vis_2_min_angle,"vis_2_min_angle/F");
  Tsig->Branch("vis_2_min_weak_track",&vis_2_min_weak_track,"vis_2_min_weak_track/F");
  Tsig->Branch("vis_2_angle_beam",&vis_2_angle_beam,"vis_2_angle_beam/F");
  Tsig->Branch("vis_2_min_angle1",&vis_2_min_angle1,"vis_2_min_angle1/F");
  Tsig->Branch("vis_2_iso_angle1",&vis_2_iso_angle1,"vis_2_iso_angle1/F");
  Tsig->Branch("vis_2_min_medium_dQ_dx",&vis_2_min_medium_dQ_dx,"vis_2_min_medium_dQ_dx/F");
  Tsig->Branch("vis_2_min_length",&vis_2_min_length,"vis_2_min_length/F");
  Tsig->Branch("vis_2_sg_length",&vis_2_sg_length,"vis_2_sg_length/F");
  Tsig->Branch("vis_2_max_angle",&vis_2_max_angle,"vis_2_max_angle/F");
  Tsig->Branch("vis_2_max_weak_track",&vis_2_max_weak_track,"vis_2_max_weak_track/F");
  Tsig->Branch("vis_2_flag",&vis_2_flag,"vis_2_flag/F");

  Tbkg->Branch("vis_2_filled",&vis_2_filled,"vis_2_filled/F");
  Tbkg->Branch("vis_2_n_vtx_segs",&vis_2_n_vtx_segs,"vis_2_n_vtx_segs/F");
  Tbkg->Branch("vis_2_min_angle",&vis_2_min_angle,"vis_2_min_angle/F");
  Tbkg->Branch("vis_2_min_weak_track",&vis_2_min_weak_track,"vis_2_min_weak_track/F");
  Tbkg->Branch("vis_2_angle_beam",&vis_2_angle_beam,"vis_2_angle_beam/F");
  Tbkg->Branch("vis_2_min_angle1",&vis_2_min_angle1,"vis_2_min_angle1/F");
  Tbkg->Branch("vis_2_iso_angle1",&vis_2_iso_angle1,"vis_2_iso_angle1/F");
  Tbkg->Branch("vis_2_min_medium_dQ_dx",&vis_2_min_medium_dQ_dx,"vis_2_min_medium_dQ_dx/F");
  Tbkg->Branch("vis_2_min_length",&vis_2_min_length,"vis_2_min_length/F");
  Tbkg->Branch("vis_2_sg_length",&vis_2_sg_length,"vis_2_sg_length/F");
  Tbkg->Branch("vis_2_max_angle",&vis_2_max_angle,"vis_2_max_angle/F");
  Tbkg->Branch("vis_2_max_weak_track",&vis_2_max_weak_track,"vis_2_max_weak_track/F");
  Tbkg->Branch("vis_2_flag",&vis_2_flag,"vis_2_flag/F");
  
  Tsig->Branch("vis_2_bdt", &bdt_value,"data/F");
  Tbkg->Branch("vis_2_bdt", &bdt_value,"data/F");

  
  TMVA::Reader *reader = new TMVA::Reader();
  
  reader->AddVariable("vis_2_n_vtx_segs",&vis_2_n_vtx_segs);
  reader->AddVariable("vis_2_min_angle",&vis_2_min_angle);
  reader->AddVariable("vis_2_min_weak_track",&vis_2_min_weak_track);
  reader->AddVariable("vis_2_angle_beam", &vis_2_angle_beam);
  reader->AddVariable("vis_2_min_angle1",&vis_2_min_angle1);
  reader->AddVariable("vis_2_iso_angle1",&vis_2_iso_angle1);
  reader->AddVariable("vis_2_min_medium_dQ_dx",&vis_2_min_medium_dQ_dx);
  reader->AddVariable("vis_2_min_length",&vis_2_min_length);
  reader->AddVariable("vis_2_sg_length",&vis_2_sg_length);
  reader->AddVariable("vis_2_max_angle",&vis_2_max_angle);
  reader->AddVariable("vis_2_max_weak_track",&vis_2_max_weak_track);
  
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

