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

  double br3_1_energy;
  int br3_1_n_shower_segments;
  int br3_1_sg_flag_trajectory;
  double br3_1_sg_direct_length;
  double br3_1_sg_length;
  double br3_1_total_main_length;
  double br3_1_total_length;
  double br3_1_iso_angle;
  int br3_1_sg_flag_topology;
  int br3_1_flag;
  
  int br3_2_n_ele;
  int br3_2_n_other;
  int br3_2_other_fid;
  int br3_2_flag;

  double br3_4_acc_length;
  double br3_4_total_length;
  int br3_4_flag;

  double br3_7_min_angle;
  int br3_7_flag;
  
  double br3_8_max_dQ_dx;
  int br3_8_n_main_segs;
  int br3_8_flag;

  sig->SetBranchAddress("br3_1_energy",&br3_1_energy);
  sig->SetBranchAddress("br3_1_n_shower_segments",&br3_1_n_shower_segments);
  sig->SetBranchAddress("br3_1_sg_flag_trajectory",&br3_1_sg_flag_trajectory);
  sig->SetBranchAddress("br3_1_sg_direct_length",&br3_1_sg_direct_length);
  sig->SetBranchAddress("br3_1_sg_length",&br3_1_sg_length);
  sig->SetBranchAddress("br3_1_total_main_length",&br3_1_total_main_length);
  sig->SetBranchAddress("br3_1_total_length",&br3_1_total_length);
  sig->SetBranchAddress("br3_1_iso_angle",&br3_1_iso_angle);
  sig->SetBranchAddress("br3_1_sg_flag_topology",&br3_1_sg_flag_topology);
  sig->SetBranchAddress("br3_1_flag",&br3_1_flag);
  
  sig->SetBranchAddress("br3_2_n_ele",&br3_2_n_ele);
  sig->SetBranchAddress("br3_2_n_other",&br3_2_n_other);
  sig->SetBranchAddress("br3_2_other_fid",&br3_2_other_fid);
  sig->SetBranchAddress("br3_2_flag",&br3_2_flag);
  
  sig->SetBranchAddress("br3_4_acc_length", &br3_4_acc_length);
  sig->SetBranchAddress("br3_4_total_length", &br3_4_total_length);
  sig->SetBranchAddress("br3_4_flag", &br3_4_flag);

  sig->SetBranchAddress("br3_7_min_angle",&br3_7_min_angle);
  sig->SetBranchAddress("br3_7_flag",&br3_7_flag);
  
  sig->SetBranchAddress("br3_8_max_dQ_dx",&br3_8_max_dQ_dx);
  sig->SetBranchAddress("br3_8_n_main_segs",&br3_8_n_main_segs);
  sig->SetBranchAddress("br3_8_flag",&br3_8_flag);

  bkg->SetBranchAddress("br3_1_energy",&br3_1_energy);
  bkg->SetBranchAddress("br3_1_n_shower_segments",&br3_1_n_shower_segments);
  bkg->SetBranchAddress("br3_1_sg_flag_trajectory",&br3_1_sg_flag_trajectory);
  bkg->SetBranchAddress("br3_1_sg_direct_length",&br3_1_sg_direct_length);
  bkg->SetBranchAddress("br3_1_sg_length",&br3_1_sg_length);
  bkg->SetBranchAddress("br3_1_total_main_length",&br3_1_total_main_length);
  bkg->SetBranchAddress("br3_1_total_length",&br3_1_total_length);
  bkg->SetBranchAddress("br3_1_iso_angle",&br3_1_iso_angle);
  bkg->SetBranchAddress("br3_1_sg_flag_topology",&br3_1_sg_flag_topology);
  bkg->SetBranchAddress("br3_1_flag",&br3_1_flag);
  
  bkg->SetBranchAddress("br3_2_n_ele",&br3_2_n_ele);
  bkg->SetBranchAddress("br3_2_n_other",&br3_2_n_other);
  bkg->SetBranchAddress("br3_2_other_fid",&br3_2_other_fid);
  bkg->SetBranchAddress("br3_2_flag",&br3_2_flag);
  
  bkg->SetBranchAddress("br3_4_acc_length", &br3_4_acc_length);
  bkg->SetBranchAddress("br3_4_total_length", &br3_4_total_length);
  bkg->SetBranchAddress("br3_4_flag", &br3_4_flag);

  bkg->SetBranchAddress("br3_7_min_angle",&br3_7_min_angle);
  bkg->SetBranchAddress("br3_7_flag",&br3_7_flag);
  
  bkg->SetBranchAddress("br3_8_max_dQ_dx",&br3_8_max_dQ_dx);
  bkg->SetBranchAddress("br3_8_n_main_segs",&br3_8_n_main_segs);
  bkg->SetBranchAddress("br3_8_flag",&br3_8_flag);

  
  
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


  float br3_1_energy_f;
  float br3_1_n_shower_segments_f;
  float br3_1_sg_flag_trajectory_f;
  float br3_1_sg_direct_length_f;
  float br3_1_sg_length_f;
  float br3_1_total_main_length_f;
  float br3_1_total_length_f;
  float br3_1_iso_angle_f;
  float br3_1_sg_flag_topology_f;
  float br3_1_flag_f;
  
  float br3_2_n_ele_f;
  float br3_2_n_other_f;
  float br3_2_other_fid_f;
  float br3_2_flag_f;

  float br3_4_acc_length_f;
  float br3_4_total_length_f;
  float br3_4_flag_f;

  float br3_7_min_angle_f;
  float br3_7_flag_f;
  
  float br3_8_max_dQ_dx_f;
  float br3_8_n_main_segs_f;
  float br3_8_flag_f;
  
  Tsig->Branch("br3_1_energy",&br3_1_energy_f,"br3_1_energy/F");
  Tsig->Branch("br3_1_n_shower_segments",&br3_1_n_shower_segments_f,"br3_1_n_shower_segments/F");
  Tsig->Branch("br3_1_sg_flag_trajectory",&br3_1_sg_flag_trajectory_f,"br3_1_sg_flag_trajectory/F");
  Tsig->Branch("br3_1_sg_direct_length",&br3_1_sg_direct_length_f,"br3_1_sg_direct_length/F");
  Tsig->Branch("br3_1_sg_length",&br3_1_sg_length_f,"br3_1_sg_length/F");
  Tsig->Branch("br3_1_total_main_length",&br3_1_total_main_length_f,"br3_1_total_main_length/F");
  Tsig->Branch("br3_1_total_length",&br3_1_total_length_f,"br3_1_total_length/F");
  Tsig->Branch("br3_1_iso_angle",&br3_1_iso_angle_f,"br3_1_iso_angle/F");
  Tsig->Branch("br3_1_sg_flag_topology",&br3_1_sg_flag_topology_f,"br3_1_sg_flag_topology/F");
  Tsig->Branch("br3_1_flag",&br3_1_flag_f,"br3_1_flag/F");
  
  Tsig->Branch("br3_2_n_ele",&br3_2_n_ele_f,"br3_2_n_ele/F");
  Tsig->Branch("br3_2_n_other",&br3_2_n_other_f,"br3_2_n_other/F");
  Tsig->Branch("br3_2_other_fid",&br3_2_other_fid_f,"br3_2_other_fid/F");
  Tsig->Branch("br3_2_flag",&br3_2_flag_f,"br3_2_flag/F");
  
  Tsig->Branch("br3_4_acc_length", &br3_4_acc_length_f, "br3_4_acc_length/F");
  Tsig->Branch("br3_4_total_length", &br3_4_total_length_f, "br3_4_total_length/F");
  Tsig->Branch("br3_4_flag", &br3_4_flag_f, "br3_4_flag/F");
  
  Tsig->Branch("br3_7_min_angle",&br3_7_min_angle_f,"br3_7_min_angle/F");
  Tsig->Branch("br3_7_flag",&br3_7_flag_f,"br3_7_flag/F");
  
  Tsig->Branch("br3_8_max_dQ_dx",&br3_8_max_dQ_dx_f,"br3_8_max_dQ_dx/F");
  Tsig->Branch("br3_8_n_main_segs",&br3_8_n_main_segs_f,"br3_8_n_main_segs/F");
  Tsig->Branch("br3_8_flag",&br3_8_flag_f,"br3_8_flag/F");


  Tbkg->Branch("br3_1_energy",&br3_1_energy_f,"br3_1_energy/F");
  Tbkg->Branch("br3_1_n_shower_segments",&br3_1_n_shower_segments_f,"br3_1_n_shower_segments/F");
  Tbkg->Branch("br3_1_sg_flag_trajectory",&br3_1_sg_flag_trajectory_f,"br3_1_sg_flag_trajectory/F");
  Tbkg->Branch("br3_1_sg_direct_length",&br3_1_sg_direct_length_f,"br3_1_sg_direct_length/F");
  Tbkg->Branch("br3_1_sg_length",&br3_1_sg_length_f,"br3_1_sg_length/F");
  Tbkg->Branch("br3_1_total_main_length",&br3_1_total_main_length_f,"br3_1_total_main_length/F");
  Tbkg->Branch("br3_1_total_length",&br3_1_total_length_f,"br3_1_total_length/F");
  Tbkg->Branch("br3_1_iso_angle",&br3_1_iso_angle_f,"br3_1_iso_angle/F");
  Tbkg->Branch("br3_1_sg_flag_topology",&br3_1_sg_flag_topology_f,"br3_1_sg_flag_topology/F");
  Tbkg->Branch("br3_1_flag",&br3_1_flag_f,"br3_1_flag/F");
  
  Tbkg->Branch("br3_2_n_ele",&br3_2_n_ele_f,"br3_2_n_ele/F");
  Tbkg->Branch("br3_2_n_other",&br3_2_n_other_f,"br3_2_n_other/F");
  Tbkg->Branch("br3_2_other_fid",&br3_2_other_fid_f,"br3_2_other_fid/F");
  Tbkg->Branch("br3_2_flag",&br3_2_flag_f,"br3_2_flag/F");
  
  Tbkg->Branch("br3_4_acc_length", &br3_4_acc_length_f, "br3_4_acc_length/F");
  Tbkg->Branch("br3_4_total_length", &br3_4_total_length_f, "br3_4_total_length/F");
  Tbkg->Branch("br3_4_flag", &br3_4_flag_f, "br3_4_flag/F");
  
  Tbkg->Branch("br3_7_min_angle",&br3_7_min_angle_f,"br3_7_min_angle/F");
  Tbkg->Branch("br3_7_flag",&br3_7_flag_f,"br3_7_flag/F");
  
  Tbkg->Branch("br3_8_max_dQ_dx",&br3_8_max_dQ_dx_f,"br3_8_max_dQ_dx/F");
  Tbkg->Branch("br3_8_n_main_segs",&br3_8_n_main_segs_f,"br3_8_n_main_segs/F");
  Tbkg->Branch("br3_8_flag",&br3_8_flag_f,"br3_8_flag/F");
  
  for (Int_t i=0;i!=sig->GetEntries();i++){
    sig->GetEntry(i);
    
    br3_1_energy_f = br3_1_energy;
    br3_1_n_shower_segments_f = br3_1_n_shower_segments;
    br3_1_sg_flag_trajectory_f = br3_1_sg_flag_trajectory;
    br3_1_sg_direct_length_f = br3_1_sg_direct_length;
    br3_1_sg_length_f = br3_1_sg_length;
    br3_1_total_main_length_f = br3_1_total_main_length;
    br3_1_total_length_f = br3_1_total_length;
    br3_1_iso_angle_f = br3_1_iso_angle;
    br3_1_sg_flag_topology_f = br3_1_sg_flag_topology;
    br3_1_flag_f = br3_1_flag;
    
    br3_2_n_ele_f = br3_2_n_ele;
    br3_2_n_other_f = br3_2_n_other;
    br3_2_other_fid_f = br3_2_other_fid;
    br3_2_flag_f = br3_2_flag;
    
    br3_4_acc_length_f = br3_4_acc_length;
    br3_4_total_length_f = br3_4_total_length;
    br3_4_flag_f = br3_4_flag;
    
    br3_7_min_angle_f = br3_7_min_angle;
    br3_7_flag_f = br3_7_flag; 
    
    br3_8_max_dQ_dx_f = br3_8_max_dQ_dx;
    br3_8_n_main_segs_f = br3_8_n_main_segs;
    br3_8_flag_f = br3_8_flag;
    
    Tsig->Fill();
  }

   for (Int_t i=0;i!=bkg->GetEntries();i++){
    bkg->GetEntry(i);
    
    br3_1_energy_f = br3_1_energy;
    br3_1_n_shower_segments_f = br3_1_n_shower_segments;
    br3_1_sg_flag_trajectory_f = br3_1_sg_flag_trajectory;
    br3_1_sg_direct_length_f = br3_1_sg_direct_length;
    br3_1_sg_length_f = br3_1_sg_length;
    br3_1_total_main_length_f = br3_1_total_main_length;
    br3_1_total_length_f = br3_1_total_length;
    br3_1_iso_angle_f = br3_1_iso_angle;
    br3_1_sg_flag_topology_f = br3_1_sg_flag_topology;
    br3_1_flag_f = br3_1_flag;
    
    br3_2_n_ele_f = br3_2_n_ele;
    br3_2_n_other_f = br3_2_n_other;
    br3_2_other_fid_f = br3_2_other_fid;
    br3_2_flag_f = br3_2_flag;
    
    br3_4_acc_length_f = br3_4_acc_length;
    br3_4_total_length_f = br3_4_total_length;
    br3_4_flag_f = br3_4_flag;
    
    br3_7_min_angle_f = br3_7_min_angle;
    br3_7_flag_f = br3_7_flag; 
    
    br3_8_max_dQ_dx_f = br3_8_max_dQ_dx;
    br3_8_n_main_segs_f = br3_8_n_main_segs;
    br3_8_flag_f = br3_8_flag;
    
    
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
    
    
    dataloader->AddVariable("br3_1_energy","br3_1_energy","MeV",'F');
    dataloader->AddVariable("br3_1_n_shower_segments","br3_1_n_shower_segments","",'F');
    dataloader->AddVariable("br3_1_sg_flag_trajectory","br3_1_sg_flag_trajectory","",'F');
    dataloader->AddVariable("br3_1_sg_direct_length","br3_1_sg_direct_length","cm",'F');
    dataloader->AddVariable("br3_1_sg_length","br3_1_sg_length","cm",'F');
    dataloader->AddVariable("br3_1_total_main_length","br3_1_total_main_length","cm",'F');
    dataloader->AddVariable("br3_1_total_length","br3_1_total_length","cm",'F');
    dataloader->AddVariable("br3_1_iso_angle","br3_1_iso_angle","deg",'F');
    dataloader->AddVariable("br3_1_sg_flag_topology","br3_1_sg_flag_topology","",'F');

    dataloader->AddVariable("br3_2_n_ele","br3_2_n_ele","",'F');
    dataloader->AddVariable("br3_2_n_other","br3_2_n_other","",'F');
    dataloader->AddVariable("br3_2_other_fid","br3_2_other_fid","",'F');

    dataloader->AddVariable("br3_4_acc_length","br3_4_acc_length","cm",'F');
    dataloader->AddVariable("br3_4_total_length","br3_4_total_length","cm",'F');

    dataloader->AddVariable("br3_7_min_angle","br3_7_min_angle","deg",'F');

    dataloader->AddVariable("br3_8_max_dQ_dx","br3_8_max_dQ_dx","MeV/cm",'F');
    dataloader->AddVariable("br3_8_n_main_segs","br3_8_n_main_segs","",'F');
    
    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "1>0"; //    682/44194
    TCut mycut_b = "br3_1_flag == 0 || br3_2_flag == 0 || br3_4_flag ==0 || br3_7_flag == 0 || br3_8_flag == 0"; // 4225/21070

    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=20000:"
					    "nTrain_Background=3500:"
					    "nTest_Signal=10000:"
					    "nTest_Background=725:"
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
    
    
    dataloader->AddVariable("br3_1_energy","br3_1_energy","MeV",'F');
    dataloader->AddVariable("br3_1_n_shower_segments","br3_1_n_shower_segments","",'F');
    dataloader->AddVariable("br3_1_sg_flag_trajectory","br3_1_sg_flag_trajectory","",'F');
    dataloader->AddVariable("br3_1_sg_direct_length","br3_1_sg_direct_length","cm",'F');
    dataloader->AddVariable("br3_1_sg_length","br3_1_sg_length","cm",'F');
    dataloader->AddVariable("br3_1_total_main_length","br3_1_total_main_length","cm",'F');
    dataloader->AddVariable("br3_1_total_length","br3_1_total_length","cm",'F');
    dataloader->AddVariable("br3_1_iso_angle","br3_1_iso_angle","deg",'F');
    dataloader->AddVariable("br3_1_sg_flag_topology","br3_1_sg_flag_topology","",'F');

    dataloader->AddVariable("br3_2_n_ele","br3_2_n_ele","",'F');
    dataloader->AddVariable("br3_2_n_other","br3_2_n_other","",'F');
    dataloader->AddVariable("br3_2_other_fid","br3_2_other_fid","",'F');

    dataloader->AddVariable("br3_4_acc_length","br3_4_acc_length","cm",'F');
    dataloader->AddVariable("br3_4_total_length","br3_4_total_length","cm",'F');

    dataloader->AddVariable("br3_7_min_angle","br3_7_min_angle","deg",'F');

    dataloader->AddVariable("br3_8_max_dQ_dx","br3_8_max_dQ_dx","MeV/cm",'F');
    dataloader->AddVariable("br3_8_n_main_segs","br3_8_n_main_segs","",'F');
    
    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "1>0"; //  
    TCut mycut_b = "br3_1_flag == 0 || br3_2_flag == 0 || br3_4_flag ==0 || br3_7_flag == 0 || br3_8_flag == 0 || br3_bdt < -0.06"; // 6242
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=20000:"
					    "nTrain_Background=5000:"
					    "nTest_Signal=10000:"
					    "nTest_Background=1242:"
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
  
  float br3_1_energy;
  float br3_1_n_shower_segments;
  float br3_1_sg_flag_trajectory;
  float br3_1_sg_direct_length;
  float br3_1_sg_length;
  float br3_1_total_main_length;
  float br3_1_total_length;
  float br3_1_iso_angle;
  float br3_1_sg_flag_topology;
  float br3_1_flag;
  
  float br3_2_n_ele;
  float br3_2_n_other;
  float br3_2_other_fid;
  float br3_2_flag;
  
  float br3_4_acc_length;
  float br3_4_total_length;
  float br3_4_flag;
  
  float br3_7_min_angle;
  float br3_7_flag;
  
  float br3_8_max_dQ_dx;
  float br3_8_n_main_segs;
  float br3_8_flag;
  
  
  sig->SetBranchAddress("br3_1_energy",&br3_1_energy);
  sig->SetBranchAddress("br3_1_n_shower_segments",&br3_1_n_shower_segments);
  sig->SetBranchAddress("br3_1_sg_flag_trajectory",&br3_1_sg_flag_trajectory);
  sig->SetBranchAddress("br3_1_sg_direct_length",&br3_1_sg_direct_length);
  sig->SetBranchAddress("br3_1_sg_length",&br3_1_sg_length);
  sig->SetBranchAddress("br3_1_total_main_length",&br3_1_total_main_length);
  sig->SetBranchAddress("br3_1_total_length",&br3_1_total_length);
  sig->SetBranchAddress("br3_1_iso_angle",&br3_1_iso_angle);
  sig->SetBranchAddress("br3_1_sg_flag_topology",&br3_1_sg_flag_topology);
  sig->SetBranchAddress("br3_1_flag",&br3_1_flag);
  
  sig->SetBranchAddress("br3_2_n_ele",&br3_2_n_ele);
  sig->SetBranchAddress("br3_2_n_other",&br3_2_n_other);
  sig->SetBranchAddress("br3_2_other_fid",&br3_2_other_fid);
  sig->SetBranchAddress("br3_2_flag",&br3_2_flag);
  
  sig->SetBranchAddress("br3_4_acc_length", &br3_4_acc_length);
  sig->SetBranchAddress("br3_4_total_length", &br3_4_total_length);
  sig->SetBranchAddress("br3_4_flag", &br3_4_flag);

  sig->SetBranchAddress("br3_7_min_angle",&br3_7_min_angle);
  sig->SetBranchAddress("br3_7_flag",&br3_7_flag);
  
  sig->SetBranchAddress("br3_8_max_dQ_dx",&br3_8_max_dQ_dx);
  sig->SetBranchAddress("br3_8_n_main_segs",&br3_8_n_main_segs);
  sig->SetBranchAddress("br3_8_flag",&br3_8_flag);

  bkg->SetBranchAddress("br3_1_energy",&br3_1_energy);
  bkg->SetBranchAddress("br3_1_n_shower_segments",&br3_1_n_shower_segments);
  bkg->SetBranchAddress("br3_1_sg_flag_trajectory",&br3_1_sg_flag_trajectory);
  bkg->SetBranchAddress("br3_1_sg_direct_length",&br3_1_sg_direct_length);
  bkg->SetBranchAddress("br3_1_sg_length",&br3_1_sg_length);
  bkg->SetBranchAddress("br3_1_total_main_length",&br3_1_total_main_length);
  bkg->SetBranchAddress("br3_1_total_length",&br3_1_total_length);
  bkg->SetBranchAddress("br3_1_iso_angle",&br3_1_iso_angle);
  bkg->SetBranchAddress("br3_1_sg_flag_topology",&br3_1_sg_flag_topology);
  bkg->SetBranchAddress("br3_1_flag",&br3_1_flag);
  
  bkg->SetBranchAddress("br3_2_n_ele",&br3_2_n_ele);
  bkg->SetBranchAddress("br3_2_n_other",&br3_2_n_other);
  bkg->SetBranchAddress("br3_2_other_fid",&br3_2_other_fid);
  bkg->SetBranchAddress("br3_2_flag",&br3_2_flag);
  
  bkg->SetBranchAddress("br3_4_acc_length", &br3_4_acc_length);
  bkg->SetBranchAddress("br3_4_total_length", &br3_4_total_length);
  bkg->SetBranchAddress("br3_4_flag", &br3_4_flag);

  bkg->SetBranchAddress("br3_7_min_angle",&br3_7_min_angle);
  bkg->SetBranchAddress("br3_7_flag",&br3_7_flag);
  
  bkg->SetBranchAddress("br3_8_max_dQ_dx",&br3_8_max_dQ_dx);
  bkg->SetBranchAddress("br3_8_n_main_segs",&br3_8_n_main_segs);
  bkg->SetBranchAddress("br3_8_flag",&br3_8_flag);

   
  
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
  
  
  Tsig->Branch("br3_1_energy",&br3_1_energy,"br3_1_energy/F");
  Tsig->Branch("br3_1_n_shower_segments",&br3_1_n_shower_segments,"br3_1_n_shower_segments/F");
  Tsig->Branch("br3_1_sg_flag_trajectory",&br3_1_sg_flag_trajectory,"br3_1_sg_flag_trajectory/F");
  Tsig->Branch("br3_1_sg_direct_length",&br3_1_sg_direct_length,"br3_1_sg_direct_length/F");
  Tsig->Branch("br3_1_sg_length",&br3_1_sg_length,"br3_1_sg_length/F");
  Tsig->Branch("br3_1_total_main_length",&br3_1_total_main_length,"br3_1_total_main_length/F");
  Tsig->Branch("br3_1_total_length",&br3_1_total_length,"br3_1_total_length/F");
  Tsig->Branch("br3_1_iso_angle",&br3_1_iso_angle,"br3_1_iso_angle/F");
  Tsig->Branch("br3_1_sg_flag_topology",&br3_1_sg_flag_topology,"br3_1_sg_flag_topology/F");
  Tsig->Branch("br3_1_flag",&br3_1_flag,"br3_1_flag/F");
  
  Tsig->Branch("br3_2_n_ele",&br3_2_n_ele,"br3_2_n_ele/F");
  Tsig->Branch("br3_2_n_other",&br3_2_n_other,"br3_2_n_other/F");
  Tsig->Branch("br3_2_other_fid",&br3_2_other_fid,"br3_2_other_fid/F");
  Tsig->Branch("br3_2_flag",&br3_2_flag,"br3_2_flag/F");
  
  Tsig->Branch("br3_4_acc_length", &br3_4_acc_length, "br3_4_acc_length/F");
  Tsig->Branch("br3_4_total_length", &br3_4_total_length, "br3_4_total_length/F");
  Tsig->Branch("br3_4_flag", &br3_4_flag, "br3_4_flag/F");
  
  Tsig->Branch("br3_7_min_angle",&br3_7_min_angle,"br3_7_min_angle/F");
  Tsig->Branch("br3_7_flag",&br3_7_flag,"br3_7_flag/F");
  
  Tsig->Branch("br3_8_max_dQ_dx",&br3_8_max_dQ_dx,"br3_8_max_dQ_dx/F");
  Tsig->Branch("br3_8_n_main_segs",&br3_8_n_main_segs,"br3_8_n_main_segs/F");
  Tsig->Branch("br3_8_flag",&br3_8_flag,"br3_8_flag/F");


  Tbkg->Branch("br3_1_energy",&br3_1_energy,"br3_1_energy/F");
  Tbkg->Branch("br3_1_n_shower_segments",&br3_1_n_shower_segments,"br3_1_n_shower_segments/F");
  Tbkg->Branch("br3_1_sg_flag_trajectory",&br3_1_sg_flag_trajectory,"br3_1_sg_flag_trajectory/F");
  Tbkg->Branch("br3_1_sg_direct_length",&br3_1_sg_direct_length,"br3_1_sg_direct_length/F");
  Tbkg->Branch("br3_1_sg_length",&br3_1_sg_length,"br3_1_sg_length/F");
  Tbkg->Branch("br3_1_total_main_length",&br3_1_total_main_length,"br3_1_total_main_length/F");
  Tbkg->Branch("br3_1_total_length",&br3_1_total_length,"br3_1_total_length/F");
  Tbkg->Branch("br3_1_iso_angle",&br3_1_iso_angle,"br3_1_iso_angle/F");
  Tbkg->Branch("br3_1_sg_flag_topology",&br3_1_sg_flag_topology,"br3_1_sg_flag_topology/F");
  Tbkg->Branch("br3_1_flag",&br3_1_flag,"br3_1_flag/F");
  
  Tbkg->Branch("br3_2_n_ele",&br3_2_n_ele,"br3_2_n_ele/F");
  Tbkg->Branch("br3_2_n_other",&br3_2_n_other,"br3_2_n_other/F");
  Tbkg->Branch("br3_2_other_fid",&br3_2_other_fid,"br3_2_other_fid/F");
  Tbkg->Branch("br3_2_flag",&br3_2_flag,"br3_2_flag/F");
  
  Tbkg->Branch("br3_4_acc_length", &br3_4_acc_length, "br3_4_acc_length/F");
  Tbkg->Branch("br3_4_total_length", &br3_4_total_length, "br3_4_total_length/F");
  Tbkg->Branch("br3_4_flag", &br3_4_flag, "br3_4_flag/F");
  
  Tbkg->Branch("br3_7_min_angle",&br3_7_min_angle,"br3_7_min_angle/F");
  Tbkg->Branch("br3_7_flag",&br3_7_flag,"br3_7_flag/F");
  
  Tbkg->Branch("br3_8_max_dQ_dx",&br3_8_max_dQ_dx,"br3_8_max_dQ_dx/F");
  Tbkg->Branch("br3_8_n_main_segs",&br3_8_n_main_segs,"br3_8_n_main_segs/F");
  Tbkg->Branch("br3_8_flag",&br3_8_flag,"br3_8_flag/F");
  
  
  Tsig->Branch("br3_bdt",&bdt_value,"data/F");
  Tbkg->Branch("br3_bdt",&bdt_value,"data/F");
  
  TMVA::Reader *reader = new TMVA::Reader();
  reader->AddVariable("br3_1_energy",&br3_1_energy);
  reader->AddVariable("br3_1_n_shower_segments",&br3_1_n_shower_segments);
  reader->AddVariable("br3_1_sg_flag_trajectory",&br3_1_sg_flag_trajectory);
  reader->AddVariable("br3_1_sg_direct_length",&br3_1_sg_direct_length);
  reader->AddVariable("br3_1_sg_length",&br3_1_sg_length);
  reader->AddVariable("br3_1_total_main_length",&br3_1_total_main_length);
  reader->AddVariable("br3_1_total_length",&br3_1_total_length);
  reader->AddVariable("br3_1_iso_angle",&br3_1_iso_angle);
  reader->AddVariable("br3_1_sg_flag_topology",&br3_1_sg_flag_topology);
  reader->AddVariable("br3_2_n_ele",&br3_2_n_ele);
  reader->AddVariable("br3_2_n_other",&br3_2_n_other);
  reader->AddVariable("br3_2_other_fid",&br3_2_other_fid);
  reader->AddVariable("br3_4_acc_length",&br3_4_acc_length);
  reader->AddVariable("br3_4_total_length",&br3_4_total_length);
  reader->AddVariable("br3_7_min_angle",&br3_7_min_angle);
  reader->AddVariable("br3_8_max_dQ_dx",&br3_8_max_dQ_dx);
  reader->AddVariable("br3_8_n_main_segs",&br3_8_n_main_segs);
  
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

