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
  
  int br1_flag;
  
  //bad reconstruction 1_1
  int br1_1_shower_type;
  int br1_1_vtx_n_segs;
  double br1_1_energy;
  int br1_1_n_segs;
  int br1_1_flag_sg_topology;
  int br1_1_flag_sg_trajectory;
  double br1_1_sg_length;
  
  // bad reconstruction 1_2
  int br1_2_n_connected;
  double br1_2_max_length;
  int br1_2_n_connected_1;
  int br1_2_n_shower_segs;
  double br1_2_max_length_ratio;
  double br1_2_shower_length;
  
  // bad_reconstruction 1_3
  int br1_3_n_connected_p;
  double br1_3_max_length_p;
  int br1_3_n_shower_main_segs;


  sig->SetBranchAddress("br1_flag",&br1_flag);

  sig->SetBranchAddress("br1_1_shower_type",&br1_1_shower_type);
  sig->SetBranchAddress("br1_1_vtx_n_segs",&br1_1_vtx_n_segs);
  sig->SetBranchAddress("br1_1_energy",&br1_1_energy);
  sig->SetBranchAddress("br1_1_n_segs",&br1_1_n_segs);
  sig->SetBranchAddress("br1_1_flag_sg_topology",&br1_1_flag_sg_topology);
  sig->SetBranchAddress("br1_1_flag_sg_trajectory",&br1_1_flag_sg_trajectory);
  sig->SetBranchAddress("br1_1_sg_length",&br1_1_sg_length);
  
  sig->SetBranchAddress("br1_2_n_connected",&br1_2_n_connected);
  sig->SetBranchAddress("br1_2_max_length",&br1_2_max_length);
  sig->SetBranchAddress("br1_2_n_connected_1",&br1_2_n_connected_1);
  sig->SetBranchAddress("br1_2_n_shower_segs",&br1_2_n_shower_segs);
  sig->SetBranchAddress("br1_2_max_length_ratio",&br1_2_max_length_ratio);
  sig->SetBranchAddress("br1_2_shower_length",&br1_2_shower_length);
  
  sig->SetBranchAddress("br1_3_n_connected_p",&br1_3_n_connected_p);
  sig->SetBranchAddress("br1_3_max_length_p",&br1_3_max_length_p);
  sig->SetBranchAddress("br1_3_n_shower_main_segs",&br1_3_n_shower_main_segs);

  
  bkg->SetBranchAddress("br1_flag",&br1_flag);

  bkg->SetBranchAddress("br1_1_shower_type",&br1_1_shower_type);
  bkg->SetBranchAddress("br1_1_vtx_n_segs",&br1_1_vtx_n_segs);
  bkg->SetBranchAddress("br1_1_energy",&br1_1_energy);
  bkg->SetBranchAddress("br1_1_n_segs",&br1_1_n_segs);
  bkg->SetBranchAddress("br1_1_flag_sg_topology",&br1_1_flag_sg_topology);
  bkg->SetBranchAddress("br1_1_flag_sg_trajectory",&br1_1_flag_sg_trajectory);
  bkg->SetBranchAddress("br1_1_sg_length",&br1_1_sg_length);
  
  bkg->SetBranchAddress("br1_2_n_connected",&br1_2_n_connected);
  bkg->SetBranchAddress("br1_2_max_length",&br1_2_max_length);
  bkg->SetBranchAddress("br1_2_n_connected_1",&br1_2_n_connected_1);
  bkg->SetBranchAddress("br1_2_n_shower_segs",&br1_2_n_shower_segs);
  bkg->SetBranchAddress("br1_2_max_length_ratio",&br1_2_max_length_ratio);
  bkg->SetBranchAddress("br1_2_shower_length",&br1_2_shower_length);
  
  bkg->SetBranchAddress("br1_3_n_connected_p",&br1_3_n_connected_p);
  bkg->SetBranchAddress("br1_3_max_length_p",&br1_3_max_length_p);
  bkg->SetBranchAddress("br1_3_n_shower_main_segs",&br1_3_n_shower_main_segs);
  
  
  
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


  
  float br1_flag_f;
  
  //bad reconstruction 1_1
  float br1_1_shower_type_f;
  float br1_1_vtx_n_segs_f;
  float br1_1_energy_f;
  float br1_1_n_segs_f;
  float br1_1_flag_sg_topology_f;
  float br1_1_flag_sg_trajectory_f;
  float br1_1_sg_length_f;
  
  // bad reconstruction 1_2
  float br1_2_n_connected_f;
  float br1_2_max_length_f;
  float br1_2_n_connected_1_f;
  float br1_2_n_shower_segs_f;
  float br1_2_max_length_ratio_f;
  float br1_2_shower_length_f;
  
  // bad_reconstruction 1_3
  float br1_3_n_connected_p_f;
  float br1_3_max_length_p_f;
  float br1_3_n_shower_main_segs_f;

  
   Tsig->Branch("br1_flag",&br1_flag_f,"br1_flag/F");
   Tsig->Branch("br1_1_shower_type",&br1_1_shower_type_f,"br1_1_shower_type/F");
   Tsig->Branch("br1_1_vtx_n_segs",&br1_1_vtx_n_segs_f,"br1_1_vtx_n_segs/F");
   Tsig->Branch("br1_1_energy",&br1_1_energy_f,"br1_1_energy/F");
   Tsig->Branch("br1_1_n_segs",&br1_1_n_segs_f,"br1_1_n_segs/F");
   Tsig->Branch("br1_1_flag_sg_topology",&br1_1_flag_sg_topology_f,"br1_1_flag_sg_topology/F");
   Tsig->Branch("br1_1_flag_sg_trajectory",&br1_1_flag_sg_trajectory_f,"br1_1_flag_sg_trajectory/F");
   Tsig->Branch("br1_1_sg_length",&br1_1_sg_length_f,"br1_1_sg_length/F");
   
   Tsig->Branch("br1_2_n_connected",&br1_2_n_connected_f,"br1_2_n_connected/F");
   Tsig->Branch("br1_2_max_length",&br1_2_max_length_f,"br1_2_max_length/F");
   Tsig->Branch("br1_2_n_connected_1",&br1_2_n_connected_1_f,"br1_2_n_connected_1/F");
   Tsig->Branch("br1_2_n_shower_segs",&br1_2_n_shower_segs_f,"br1_2_n_shower_segs/F");
   Tsig->Branch("br1_2_max_length_ratio",&br1_2_max_length_ratio_f,"br1_2_max_length_ratio/F");
   Tsig->Branch("br1_2_shower_length",&br1_2_shower_length_f,"br1_2_shower_length/F");

   Tsig->Branch("br1_3_n_connected_p",&br1_3_n_connected_p_f,"br1_3_n_connected_p/F");
   Tsig->Branch("br1_3_max_length_p",&br1_3_max_length_p_f,"br1_3_max_length_p/F");
   Tsig->Branch("br1_3_n_shower_main_segs",&br1_3_n_shower_main_segs_f,"br1_3_n_shower_main_segs/F");
   
   Tbkg->Branch("br1_flag",&br1_flag_f,"br1_flag/F");
   Tbkg->Branch("br1_1_shower_type",&br1_1_shower_type_f,"br1_1_shower_type/F");
   Tbkg->Branch("br1_1_vtx_n_segs",&br1_1_vtx_n_segs_f,"br1_1_vtx_n_segs/F");
   Tbkg->Branch("br1_1_energy",&br1_1_energy_f,"br1_1_energy/F");
   Tbkg->Branch("br1_1_n_segs",&br1_1_n_segs_f,"br1_1_n_segs/F");
   Tbkg->Branch("br1_1_flag_sg_topology",&br1_1_flag_sg_topology_f,"br1_1_flag_sg_topology/F");
   Tbkg->Branch("br1_1_flag_sg_trajectory",&br1_1_flag_sg_trajectory_f,"br1_1_flag_sg_trajectory/F");
   Tbkg->Branch("br1_1_sg_length",&br1_1_sg_length_f,"br1_1_sg_length/F");
   
   Tbkg->Branch("br1_2_n_connected",&br1_2_n_connected_f,"br1_2_n_connected/F");
   Tbkg->Branch("br1_2_max_length",&br1_2_max_length_f,"br1_2_max_length/F");
   Tbkg->Branch("br1_2_n_connected_1",&br1_2_n_connected_1_f,"br1_2_n_connected_1/F");
   Tbkg->Branch("br1_2_n_shower_segs",&br1_2_n_shower_segs_f,"br1_2_n_shower_segs/F");
   Tbkg->Branch("br1_2_max_length_ratio",&br1_2_max_length_ratio_f,"br1_2_max_length_ratio/F");
   Tbkg->Branch("br1_2_shower_length",&br1_2_shower_length_f,"br1_2_shower_length/F");

   Tbkg->Branch("br1_3_n_connected_p",&br1_3_n_connected_p_f,"br1_3_n_connected_p/F");
   Tbkg->Branch("br1_3_max_length_p",&br1_3_max_length_p_f,"br1_3_max_length_p/F");
   Tbkg->Branch("br1_3_n_shower_main_segs",&br1_3_n_shower_main_segs_f,"br1_3_n_shower_main_segs/F");
  
  for (Int_t i=0;i!=sig->GetEntries();i++){
    sig->GetEntry(i);
    
    br1_flag_f = br1_flag;
    
    //bad reconstruction 1_1
    br1_1_shower_type_f = br1_1_shower_type;
    br1_1_vtx_n_segs_f = br1_1_vtx_n_segs;
    br1_1_energy_f = br1_1_energy;
    br1_1_n_segs_f = br1_1_n_segs;
    br1_1_flag_sg_topology_f = br1_1_flag_sg_topology;
    br1_1_flag_sg_trajectory_f = br1_1_flag_sg_trajectory;
    br1_1_sg_length_f = br1_1_sg_length;
    
    // bad reconstruction 1_2
    br1_2_n_connected_f = br1_2_n_connected;
    br1_2_max_length_f = br1_2_max_length;
    br1_2_n_connected_1_f = br1_2_n_connected_1;
    br1_2_n_shower_segs_f = br1_2_n_shower_segs;
    br1_2_max_length_ratio_f = br1_2_max_length_ratio;
    br1_2_shower_length_f = br1_2_shower_length;
    
    // bad_reconstruction 1_3
    br1_3_n_connected_p_f = br1_3_n_connected_p;
    br1_3_max_length_p_f = br1_3_max_length_p; 
    br1_3_n_shower_main_segs_f = br1_3_n_shower_main_segs;
    
    Tsig->Fill();
  }

   for (Int_t i=0;i!=bkg->GetEntries();i++){
    bkg->GetEntry(i);
    
     br1_flag_f = br1_flag;
    
    //bad reconstruction 1_1
    br1_1_shower_type_f = br1_1_shower_type;
    br1_1_vtx_n_segs_f = br1_1_vtx_n_segs;
    br1_1_energy_f = br1_1_energy;
    br1_1_n_segs_f = br1_1_n_segs;
    br1_1_flag_sg_topology_f = br1_1_flag_sg_topology;
    br1_1_flag_sg_trajectory_f = br1_1_flag_sg_trajectory;
    br1_1_sg_length_f = br1_1_sg_length;
    
    // bad reconstruction 1_2
    br1_2_n_connected_f = br1_2_n_connected;
    br1_2_max_length_f = br1_2_max_length;
    br1_2_n_connected_1_f = br1_2_n_connected_1;
    br1_2_n_shower_segs_f = br1_2_n_shower_segs;
    br1_2_max_length_ratio_f = br1_2_max_length_ratio;
    br1_2_shower_length_f = br1_2_shower_length;
    
    // bad_reconstruction 1_3
    br1_3_n_connected_p_f = br1_3_n_connected_p;
    br1_3_max_length_p_f = br1_3_max_length_p; 
    br1_3_n_shower_main_segs_f = br1_3_n_shower_main_segs;
    
    
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
   
    
    dataloader->AddVariable("br1_1_shower_type","br1_1_shower_type","",'F');
    dataloader->AddVariable("br1_1_vtx_n_segs","br1_1_vtx_n_segs","",'F');
    dataloader->AddVariable("br1_1_energy","br1_1_energy","MeV",'F');
    dataloader->AddVariable("br1_1_n_segs","br1_1_n_segs","",'F');
    dataloader->AddVariable("br1_1_flag_sg_topology","br1_1_flag_sg_topology","",'F');
    dataloader->AddVariable("br1_1_flag_sg_trajectory","br1_1_flag_sg_trajectory","",'F');
    dataloader->AddVariable("br1_1_sg_length","br1_1_sg_length","cm",'F');

    dataloader->AddVariable("br1_2_n_connected","br1_2_n_connected","",'F');
    dataloader->AddVariable("br1_2_max_length","br1_2_max_length","cm",'F');
    dataloader->AddVariable("br1_2_n_connected_1","br1_2_n_connected_1","",'F');
    dataloader->AddVariable("br1_2_n_shower_segs","br1_2_n_shower_segs","",'F');
    dataloader->AddVariable("br1_2_max_length_ratio","br1_2_max_length_ratio","",'F');
    dataloader->AddVariable("br1_2_shower_length","br1_2_shower_length","cm",'F');

    dataloader->AddVariable("br1_3_n_connected_p","br1_3_n_connected_p","",'F');
    dataloader->AddVariable("br1_3_max_length_p","br1_3_max_length_p","cm",'F');
    dataloader->AddVariable("br1_3_n_shower_main_segs","br1_3_n_shower_main_segs","",'F');
    
    
    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "1>0"; //  1134/44194
    TCut mycut_b = "br1_flag == 0 && (!(truth_nue==1 && truth_CC==1))"; // 5819/21070
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=20000:"
					    "nTrain_Background=5200:"
					    "nTest_Signal=10000:"
					    "nTest_Background=619:"
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
   
    
    dataloader->AddVariable("br1_1_shower_type","br1_1_shower_type","",'F');
    dataloader->AddVariable("br1_1_vtx_n_segs","br1_1_vtx_n_segs","",'F');
    dataloader->AddVariable("br1_1_energy","br1_1_energy","MeV",'F');
    dataloader->AddVariable("br1_1_n_segs","br1_1_n_segs","",'F');
    dataloader->AddVariable("br1_1_flag_sg_topology","br1_1_flag_sg_topology","",'F');
    dataloader->AddVariable("br1_1_flag_sg_trajectory","br1_1_flag_sg_trajectory","",'F');
    dataloader->AddVariable("br1_1_sg_length","br1_1_sg_length","cm",'F');

    dataloader->AddVariable("br1_2_n_connected","br1_2_n_connected","",'F');
    dataloader->AddVariable("br1_2_max_length","br1_2_max_length","cm",'F');
    dataloader->AddVariable("br1_2_n_connected_1","br1_2_n_connected_1","",'F');
    dataloader->AddVariable("br1_2_n_shower_segs","br1_2_n_shower_segs","",'F');
    dataloader->AddVariable("br1_2_max_length_ratio","br1_2_max_length_ratio","",'F');
    dataloader->AddVariable("br1_2_shower_length","br1_2_shower_length","cm",'F');

    dataloader->AddVariable("br1_3_n_connected_p","br1_3_n_connected_p","",'F');
    dataloader->AddVariable("br1_3_max_length_p","br1_3_max_length_p","cm",'F');
    dataloader->AddVariable("br1_3_n_shower_main_segs","br1_3_n_shower_main_segs","",'F');
    
    
    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "1>0"; //  
    TCut mycut_b = "(br1_flag == 0 || br1_bdt < -0.05)&& (!(truth_nue==1 && truth_CC==1))"; //  6621
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=20000:"
					    "nTrain_Background=5410:"
					    "nTest_Signal=10000:"
					    "nTest_Background=1211:"
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

  
  float br1_flag;
  
  //bad reconstruction 1_1
  float br1_1_shower_type;
  float br1_1_vtx_n_segs;
  float br1_1_energy;
  float br1_1_n_segs;
  float br1_1_flag_sg_topology;
  float br1_1_flag_sg_trajectory;
  float br1_1_sg_length;
  
  // bad reconstruction 1_2
  float br1_2_n_connected;
  float br1_2_max_length;
  float br1_2_n_connected_1;
  float br1_2_n_shower_segs;
  float br1_2_max_length_ratio;
  float br1_2_shower_length;
  
  // bad_reconstruction 1_3
  float br1_3_n_connected_p;
  float br1_3_max_length_p;
  float br1_3_n_shower_main_segs;
  
   sig->SetBranchAddress("br1_flag",&br1_flag);

  sig->SetBranchAddress("br1_1_shower_type",&br1_1_shower_type);
  sig->SetBranchAddress("br1_1_vtx_n_segs",&br1_1_vtx_n_segs);
  sig->SetBranchAddress("br1_1_energy",&br1_1_energy);
  sig->SetBranchAddress("br1_1_n_segs",&br1_1_n_segs);
  sig->SetBranchAddress("br1_1_flag_sg_topology",&br1_1_flag_sg_topology);
  sig->SetBranchAddress("br1_1_flag_sg_trajectory",&br1_1_flag_sg_trajectory);
  sig->SetBranchAddress("br1_1_sg_length",&br1_1_sg_length);
  
  sig->SetBranchAddress("br1_2_n_connected",&br1_2_n_connected);
  sig->SetBranchAddress("br1_2_max_length",&br1_2_max_length);
  sig->SetBranchAddress("br1_2_n_connected_1",&br1_2_n_connected_1);
  sig->SetBranchAddress("br1_2_n_shower_segs",&br1_2_n_shower_segs);
  sig->SetBranchAddress("br1_2_max_length_ratio",&br1_2_max_length_ratio);
  sig->SetBranchAddress("br1_2_shower_length",&br1_2_shower_length);
  
  sig->SetBranchAddress("br1_3_n_connected_p",&br1_3_n_connected_p);
  sig->SetBranchAddress("br1_3_max_length_p",&br1_3_max_length_p);
  sig->SetBranchAddress("br1_3_n_shower_main_segs",&br1_3_n_shower_main_segs);

  
  bkg->SetBranchAddress("br1_flag",&br1_flag);

  bkg->SetBranchAddress("br1_1_shower_type",&br1_1_shower_type);
  bkg->SetBranchAddress("br1_1_vtx_n_segs",&br1_1_vtx_n_segs);
  bkg->SetBranchAddress("br1_1_energy",&br1_1_energy);
  bkg->SetBranchAddress("br1_1_n_segs",&br1_1_n_segs);
  bkg->SetBranchAddress("br1_1_flag_sg_topology",&br1_1_flag_sg_topology);
  bkg->SetBranchAddress("br1_1_flag_sg_trajectory",&br1_1_flag_sg_trajectory);
  bkg->SetBranchAddress("br1_1_sg_length",&br1_1_sg_length);
  
  bkg->SetBranchAddress("br1_2_n_connected",&br1_2_n_connected);
  bkg->SetBranchAddress("br1_2_max_length",&br1_2_max_length);
  bkg->SetBranchAddress("br1_2_n_connected_1",&br1_2_n_connected_1);
  bkg->SetBranchAddress("br1_2_n_shower_segs",&br1_2_n_shower_segs);
  bkg->SetBranchAddress("br1_2_max_length_ratio",&br1_2_max_length_ratio);
  bkg->SetBranchAddress("br1_2_shower_length",&br1_2_shower_length);
  
  bkg->SetBranchAddress("br1_3_n_connected_p",&br1_3_n_connected_p);
  bkg->SetBranchAddress("br1_3_max_length_p",&br1_3_max_length_p);
  bkg->SetBranchAddress("br1_3_n_shower_main_segs",&br1_3_n_shower_main_segs);
    
  
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

  
  Tsig->Branch("br1_flag",&br1_flag,"br1_flag/F");
   Tsig->Branch("br1_1_shower_type",&br1_1_shower_type,"br1_1_shower_type/F");
   Tsig->Branch("br1_1_vtx_n_segs",&br1_1_vtx_n_segs,"br1_1_vtx_n_segs/F");
   Tsig->Branch("br1_1_energy",&br1_1_energy,"br1_1_energy/F");
   Tsig->Branch("br1_1_n_segs",&br1_1_n_segs,"br1_1_n_segs/F");
   Tsig->Branch("br1_1_flag_sg_topology",&br1_1_flag_sg_topology,"br1_1_flag_sg_topology/F");
   Tsig->Branch("br1_1_flag_sg_trajectory",&br1_1_flag_sg_trajectory,"br1_1_flag_sg_trajectory/F");
   Tsig->Branch("br1_1_sg_length",&br1_1_sg_length,"br1_1_sg_length/F");
   
   Tsig->Branch("br1_2_n_connected",&br1_2_n_connected,"br1_2_n_connected/F");
   Tsig->Branch("br1_2_max_length",&br1_2_max_length,"br1_2_max_length/F");
   Tsig->Branch("br1_2_n_connected_1",&br1_2_n_connected_1,"br1_2_n_connected_1/F");
   Tsig->Branch("br1_2_n_shower_segs",&br1_2_n_shower_segs,"br1_2_n_shower_segs/F");
   Tsig->Branch("br1_2_max_length_ratio",&br1_2_max_length_ratio,"br1_2_max_length_ratio/F");
   Tsig->Branch("br1_2_shower_length",&br1_2_shower_length,"br1_2_shower_length/F");

   Tsig->Branch("br1_3_n_connected_p",&br1_3_n_connected_p,"br1_3_n_connected_p/F");
   Tsig->Branch("br1_3_max_length_p",&br1_3_max_length_p,"br1_3_max_length_p/F");
   Tsig->Branch("br1_3_n_shower_main_segs",&br1_3_n_shower_main_segs,"br1_3_n_shower_main_segs/F");
   
   Tbkg->Branch("br1_flag",&br1_flag,"br1_flag/F");
   Tbkg->Branch("br1_1_shower_type",&br1_1_shower_type,"br1_1_shower_type/F");
   Tbkg->Branch("br1_1_vtx_n_segs",&br1_1_vtx_n_segs,"br1_1_vtx_n_segs/F");
   Tbkg->Branch("br1_1_energy",&br1_1_energy,"br1_1_energy/F");
   Tbkg->Branch("br1_1_n_segs",&br1_1_n_segs,"br1_1_n_segs/F");
   Tbkg->Branch("br1_1_flag_sg_topology",&br1_1_flag_sg_topology,"br1_1_flag_sg_topology/F");
   Tbkg->Branch("br1_1_flag_sg_trajectory",&br1_1_flag_sg_trajectory,"br1_1_flag_sg_trajectory/F");
   Tbkg->Branch("br1_1_sg_length",&br1_1_sg_length,"br1_1_sg_length/F");
   
   Tbkg->Branch("br1_2_n_connected",&br1_2_n_connected,"br1_2_n_connected/F");
   Tbkg->Branch("br1_2_max_length",&br1_2_max_length,"br1_2_max_length/F");
   Tbkg->Branch("br1_2_n_connected_1",&br1_2_n_connected_1,"br1_2_n_connected_1/F");
   Tbkg->Branch("br1_2_n_shower_segs",&br1_2_n_shower_segs,"br1_2_n_shower_segs/F");
   Tbkg->Branch("br1_2_max_length_ratio",&br1_2_max_length_ratio,"br1_2_max_length_ratio/F");
   Tbkg->Branch("br1_2_shower_length",&br1_2_shower_length,"br1_2_shower_length/F");

   Tbkg->Branch("br1_3_n_connected_p",&br1_3_n_connected_p,"br1_3_n_connected_p/F");
   Tbkg->Branch("br1_3_max_length_p",&br1_3_max_length_p,"br1_3_max_length_p/F");
   Tbkg->Branch("br1_3_n_shower_main_segs",&br1_3_n_shower_main_segs,"br1_3_n_shower_main_segs/F");
 
  
  
  Tsig->Branch("br1_bdt",&bdt_value,"data/F");
  Tbkg->Branch("br1_bdt",&bdt_value,"data/F");
  
  TMVA::Reader *reader = new TMVA::Reader();
  reader->AddVariable("br1_1_shower_type",&br1_1_shower_type);
  reader->AddVariable("br1_1_vtx_n_segs",&br1_1_vtx_n_segs);
  reader->AddVariable("br1_1_energy",&br1_1_energy);
  reader->AddVariable("br1_1_n_segs",&br1_1_n_segs);
  reader->AddVariable("br1_1_flag_sg_topology",&br1_1_flag_sg_topology);
  reader->AddVariable("br1_1_flag_sg_trajectory",&br1_1_flag_sg_trajectory);
  reader->AddVariable("br1_1_sg_length",&br1_1_sg_length);

  reader->AddVariable("br1_2_n_connected",&br1_2_n_connected);
  reader->AddVariable("br1_2_max_length",&br1_2_max_length);
  reader->AddVariable("br1_2_n_connected_1",&br1_2_n_connected_1);
  reader->AddVariable("br1_2_n_shower_segs",&br1_2_n_shower_segs);
  reader->AddVariable("br1_2_max_length_ratio",&br1_2_max_length_ratio);
  reader->AddVariable("br1_2_shower_length",&br1_2_shower_length);

  reader->AddVariable("br1_3_n_connected_p",&br1_3_n_connected_p);
  reader->AddVariable("br1_3_max_length_p",&br1_3_max_length_p);
  reader->AddVariable("br1_3_n_shower_main_segs",&br1_3_n_shower_main_segs);
  
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

