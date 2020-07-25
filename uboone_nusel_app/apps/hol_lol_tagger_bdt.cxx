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
  
  int hol_1_n_valid_tracks;
  double hol_1_min_angle;
  double hol_1_energy;
  int hol_1_flag_all_shower;
  double hol_1_min_length;
  int hol_1_flag;
  
  double hol_2_min_angle;
  double hol_2_medium_dQ_dx;
  int hol_2_ncount;
  int hol_2_flag;

  double lol_3_angle_beam;
  int lol_3_n_valid_tracks;
  double lol_3_min_angle;
  int lol_3_vtx_n_segs;
  double lol_3_shower_main_length;
  int lol_3_n_out;
  int lol_3_n_sum;    
  int lol_3_flag;
  
  sig->SetBranchAddress("hol_1_n_valid_tracks", &hol_1_n_valid_tracks);
  sig->SetBranchAddress("hol_1_min_angle", &hol_1_min_angle);
  sig->SetBranchAddress("hol_1_energy", &hol_1_energy);
  sig->SetBranchAddress("hol_1_all_shower", &hol_1_flag_all_shower);
  sig->SetBranchAddress("hol_1_min_length", &hol_1_min_length);
  sig->SetBranchAddress("hol_1_flag", &hol_1_flag);
  
  sig->SetBranchAddress("hol_2_min_angle", &hol_2_min_angle);
  sig->SetBranchAddress("hol_2_medium_dQ_dx", &hol_2_medium_dQ_dx);
  sig->SetBranchAddress("hol_2_ncount", &hol_2_ncount);
  sig->SetBranchAddress("hol_2_flag", &hol_2_flag);

  sig->SetBranchAddress("lol_3_angle_beam",&lol_3_angle_beam);
  sig->SetBranchAddress("lol_3_n_valid_tracks",&lol_3_n_valid_tracks);
  sig->SetBranchAddress("lol_3_min_angle",&lol_3_min_angle);
  sig->SetBranchAddress("lol_3_vtx_n_segs",&lol_3_vtx_n_segs);
  sig->SetBranchAddress("lol_3_shower_main_length",&lol_3_shower_main_length);
  sig->SetBranchAddress("lol_3_n_out",&lol_3_n_out);
  sig->SetBranchAddress("lol_3_n_sum",&lol_3_n_sum);
  sig->SetBranchAddress("lol_3_flag",&lol_3_flag);

  bkg->SetBranchAddress("hol_1_n_valid_tracks", &hol_1_n_valid_tracks);
  bkg->SetBranchAddress("hol_1_min_angle", &hol_1_min_angle);
  bkg->SetBranchAddress("hol_1_energy", &hol_1_energy);
  bkg->SetBranchAddress("hol_1_all_shower", &hol_1_flag_all_shower);
  bkg->SetBranchAddress("hol_1_min_length", &hol_1_min_length);
  bkg->SetBranchAddress("hol_1_flag", &hol_1_flag);
  
  bkg->SetBranchAddress("hol_2_min_angle", &hol_2_min_angle);
  bkg->SetBranchAddress("hol_2_medium_dQ_dx", &hol_2_medium_dQ_dx);
  bkg->SetBranchAddress("hol_2_ncount", &hol_2_ncount);
  bkg->SetBranchAddress("hol_2_flag", &hol_2_flag);

  bkg->SetBranchAddress("lol_3_angle_beam",&lol_3_angle_beam);
  bkg->SetBranchAddress("lol_3_n_valid_tracks",&lol_3_n_valid_tracks);
  bkg->SetBranchAddress("lol_3_min_angle",&lol_3_min_angle);
  bkg->SetBranchAddress("lol_3_vtx_n_segs",&lol_3_vtx_n_segs);
  bkg->SetBranchAddress("lol_3_shower_main_length",&lol_3_shower_main_length);
  bkg->SetBranchAddress("lol_3_n_out",&lol_3_n_out);
  bkg->SetBranchAddress("lol_3_n_sum",&lol_3_n_sum);
  bkg->SetBranchAddress("lol_3_flag",&lol_3_flag);
  
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
  
  float hol_1_n_valid_tracks_f;
  float hol_1_min_angle_f;
  float hol_1_energy_f;
  float hol_1_flag_all_shower_f;
  float hol_1_min_length_f;
  float hol_1_flag_f;
  float hol_2_min_angle_f;
  float hol_2_medium_dQ_dx_f;
  float hol_2_ncount_f;
  float hol_2_flag_f;
  float lol_3_angle_beam_f;
  float lol_3_n_valid_tracks_f;
  float lol_3_min_angle_f;
  float lol_3_vtx_n_segs_f;
  float lol_3_shower_main_length_f;
  float lol_3_n_out_f;
  float lol_3_n_sum_f;    
  float lol_3_flag_f;

  Tsig->Branch("hol_1_n_valid_tracks", &hol_1_n_valid_tracks_f,"hol_1_n_valid_tracks/F");
  Tsig->Branch("hol_1_min_angle", &hol_1_min_angle_f,"hol_1_min_angle/F");
  Tsig->Branch("hol_1_energy", &hol_1_energy_f,"hol_1_energy/F");
  Tsig->Branch("hol_1_flag_all_shower", &hol_1_flag_all_shower_f,"hol_1_flag_all_shower/F");
  Tsig->Branch("hol_1_min_length", &hol_1_min_length_f,"hol_1_min_length/F");
  Tsig->Branch("hol_1_flag", &hol_1_flag_f,"hol_1_flag/F");
  
  Tsig->Branch("hol_2_min_angle", &hol_2_min_angle_f,"hol_2_min_angle/F");
  Tsig->Branch("hol_2_medium_dQ_dx", &hol_2_medium_dQ_dx_f,"hol_2_medium_dQ_dx/F");
  Tsig->Branch("hol_2_ncount", &hol_2_ncount_f,"hol_2_ncount/F");
  Tsig->Branch("hol_2_flag", &hol_2_flag_f,"hol_2_flag/F");

  Tsig->Branch("lol_3_angle_beam",&lol_3_angle_beam_f,"lol_3_angle_beam/F");
  Tsig->Branch("lol_3_n_valid_tracks",&lol_3_n_valid_tracks_f,"lol_3_n_valid_tracks/F");
  Tsig->Branch("lol_3_min_angle",&lol_3_min_angle_f,"lol_3_min_angle/F");
  Tsig->Branch("lol_3_vtx_n_segs",&lol_3_vtx_n_segs_f,"lol_3_vtx_n_segs/F");
  Tsig->Branch("lol_3_shower_main_length",&lol_3_shower_main_length_f,"lol_3_shower_main_length/F");
  Tsig->Branch("lol_3_n_out",&lol_3_n_out_f,"lol_3_n_out/F");
  Tsig->Branch("lol_3_n_sum",&lol_3_n_sum_f,"lol_3_n_sum/F");
  Tsig->Branch("lol_3_flag",&lol_3_flag_f,"lol_3_flag/F");
  
  Tbkg->Branch("hol_1_n_valid_tracks", &hol_1_n_valid_tracks_f,"hol_1_n_valid_tracks/F");
  Tbkg->Branch("hol_1_min_angle", &hol_1_min_angle_f,"hol_1_min_angle/F");
  Tbkg->Branch("hol_1_energy", &hol_1_energy_f,"hol_1_energy/F");
  Tbkg->Branch("hol_1_flag_all_shower", &hol_1_flag_all_shower_f,"hol_1_flag_all_shower/F");
  Tbkg->Branch("hol_1_min_length", &hol_1_min_length_f,"hol_1_min_length/F");
  Tbkg->Branch("hol_1_flag", &hol_1_flag_f,"hol_1_flag/F");
  
  Tbkg->Branch("hol_2_min_angle", &hol_2_min_angle_f,"hol_2_min_angle/F");
  Tbkg->Branch("hol_2_medium_dQ_dx", &hol_2_medium_dQ_dx_f,"hol_2_medium_dQ_dx/F");
  Tbkg->Branch("hol_2_ncount", &hol_2_ncount_f,"hol_2_ncount/F");
  Tbkg->Branch("hol_2_flag", &hol_2_flag_f,"hol_2_flag/F");

  Tbkg->Branch("lol_3_angle_beam",&lol_3_angle_beam_f,"lol_3_angle_beam/F");
  Tbkg->Branch("lol_3_n_valid_tracks",&lol_3_n_valid_tracks_f,"lol_3_n_valid_tracks/F");
  Tbkg->Branch("lol_3_min_angle",&lol_3_min_angle_f,"lol_3_min_angle/F");
  Tbkg->Branch("lol_3_vtx_n_segs",&lol_3_vtx_n_segs_f,"lol_3_vtx_n_segs/F");
  Tbkg->Branch("lol_3_shower_main_length",&lol_3_shower_main_length_f,"lol_3_shower_main_length/F");
  Tbkg->Branch("lol_3_n_out",&lol_3_n_out_f,"lol_3_n_out/F");
  Tbkg->Branch("lol_3_n_sum",&lol_3_n_sum_f,"lol_3_n_sum/F");
  Tbkg->Branch("lol_3_flag",&lol_3_flag_f,"lol_3_flag/F");
  
  for (Int_t i=0;i!=sig->GetEntries();i++){
    sig->GetEntry(i);
    
    hol_1_n_valid_tracks_f = hol_1_n_valid_tracks;
    hol_1_min_angle_f = hol_1_min_angle;
    hol_1_energy_f = hol_1_energy;
    hol_1_flag_all_shower_f = hol_1_flag_all_shower;
    hol_1_min_length_f = hol_1_min_length;
    hol_1_flag_f = hol_1_flag;
    hol_2_min_angle_f = hol_2_min_angle;
    hol_2_medium_dQ_dx_f = hol_2_medium_dQ_dx;
    hol_2_ncount_f = hol_2_ncount;
    hol_2_flag_f = hol_2_flag;
    lol_3_angle_beam_f = lol_3_angle_beam;
    lol_3_n_valid_tracks_f = lol_3_n_valid_tracks;
    lol_3_min_angle_f = lol_3_min_angle;
    lol_3_vtx_n_segs_f = lol_3_vtx_n_segs;
    lol_3_shower_main_length_f = lol_3_shower_main_length;
    lol_3_n_out_f = lol_3_n_out;
    lol_3_n_sum_f = lol_3_n_sum;    
    lol_3_flag_f = lol_3_flag;
    
    Tsig->Fill();
  }

   for (Int_t i=0;i!=bkg->GetEntries();i++){
    bkg->GetEntry(i);

    hol_1_n_valid_tracks_f = hol_1_n_valid_tracks;
    hol_1_min_angle_f = hol_1_min_angle;
    hol_1_energy_f = hol_1_energy;
    hol_1_flag_all_shower_f = hol_1_flag_all_shower;
    hol_1_min_length_f = hol_1_min_length;
    hol_1_flag_f = hol_1_flag;
    hol_2_min_angle_f = hol_2_min_angle;
    hol_2_medium_dQ_dx_f = hol_2_medium_dQ_dx;
    hol_2_ncount_f = hol_2_ncount;
    hol_2_flag_f = hol_2_flag;
    lol_3_angle_beam_f = lol_3_angle_beam;
    lol_3_n_valid_tracks_f = lol_3_n_valid_tracks;
    lol_3_min_angle_f = lol_3_min_angle;
    lol_3_vtx_n_segs_f = lol_3_vtx_n_segs;
    lol_3_shower_main_length_f = lol_3_shower_main_length;
    lol_3_n_out_f = lol_3_n_out;
    lol_3_n_sum_f = lol_3_n_sum;    
    lol_3_flag_f = lol_3_flag;
    
    
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
    
    
    dataloader->AddVariable("hol_1_n_valid_tracks","hol_1_n_valid_tracks","",'F');
    dataloader->AddVariable("hol_1_min_angle","hol_1_min_angle","deg",'F');
    dataloader->AddVariable("hol_1_energy","hol_1_energy","MeV",'F');
    dataloader->AddVariable("hol_1_flag_all_shower","hol_1_flag_all_shower","",'F');
    dataloader->AddVariable("hol_1_min_length","hol_1_min_length","cm",'F');

    dataloader->AddVariable("hol_2_min_angle","hol_2_min_angle","deg",'F');
    dataloader->AddVariable("hol_2_medium_dQ_dx","hol_2_medium_dQ_dx","MeV/cm",'F');
    dataloader->AddVariable("hol_2_ncount","hol_2_ncount","",'F');

    dataloader->AddVariable("lol_3_angle_beam","lol_3_angle_beam","deg",'F');
    dataloader->AddVariable("lol_3_n_valid_tracks","lol_3_n_valid_tracks","",'F');
    dataloader->AddVariable("lol_3_min_angle","lol_3_min_angle","deg",'F');
    dataloader->AddVariable("lol_3_vtx_n_segs","lol_3_vtx_n_segs","",'F');
    dataloader->AddVariable("lol_3_shower_main_length","lol_3_shower_main_length","cm",'F');
    dataloader->AddVariable("lol_3_n_out","lol_3_n_out","",'F');
    dataloader->AddVariable("lol_3_n_sum","lol_3_n_sum","",'F');
    
    

    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "1>0"; // 359 (stl) + 227 (brm)  21 (lem) or 588 (total)/44194
    TCut mycut_b = "(hol_1_flag==0 || hol_2_flag ==0 || lol_3_flag ==0) && (!(truth_nue==1 && truth_CC==1))"; // 750 (stl) + 274 (brm) 198 (lem) or 1087 (total)/21070
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=20000:"
					    "nTrain_Background=950:"
					    "nTest_Signal=10000:"
					    "nTest_Background=137:"
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
    
    
    dataloader->AddVariable("hol_1_n_valid_tracks","hol_1_n_valid_tracks","",'F');
    dataloader->AddVariable("hol_1_min_angle","hol_1_min_angle","deg",'F');
    dataloader->AddVariable("hol_1_energy","hol_1_energy","MeV",'F');
    dataloader->AddVariable("hol_1_flag_all_shower","hol_1_flag_all_shower","",'F');
    dataloader->AddVariable("hol_1_min_length","hol_1_min_length","cm",'F');

    dataloader->AddVariable("hol_2_min_angle","hol_2_min_angle","deg",'F');
    dataloader->AddVariable("hol_2_medium_dQ_dx","hol_2_medium_dQ_dx","MeV/cm",'F');
    dataloader->AddVariable("hol_2_ncount","hol_2_ncount","",'F');

    dataloader->AddVariable("lol_3_angle_beam","lol_3_angle_beam","deg",'F');
    dataloader->AddVariable("lol_3_n_valid_tracks","lol_3_n_valid_tracks","",'F');
    dataloader->AddVariable("lol_3_min_angle","lol_3_min_angle","deg",'F');
    dataloader->AddVariable("lol_3_vtx_n_segs","lol_3_vtx_n_segs","",'F');
    dataloader->AddVariable("lol_3_shower_main_length","lol_3_shower_main_length","cm",'F');
    dataloader->AddVariable("lol_3_n_out","lol_3_n_out","",'F');
    dataloader->AddVariable("lol_3_n_sum","lol_3_n_sum","",'F');
    
    

    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "1>0"; // 359 (stl) + 227 (brm)  21 (lem) or 588 (total)/44194
    TCut mycut_b = "(hol_1_flag==0 || hol_2_flag ==0 || lol_3_flag ==0 || hol_lol_bdt < 0.0)&& (!(truth_nue==1 && truth_CC==1))"; // 750 (stl) + 274 (brm) 198 (lem) or 1872 (total)/21070
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=20000:"
					    "nTrain_Background=1600:"
					    "nTest_Signal=10000:"
					    "nTest_Background=272:"
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
  
  float hol_1_n_valid_tracks;
  float hol_1_min_angle;
  float hol_1_energy;
  float hol_1_flag_all_shower;
  float hol_1_min_length;
  float hol_1_flag;
  float hol_2_min_angle;
  float hol_2_medium_dQ_dx;
  float hol_2_ncount;
  float hol_2_flag;
  float lol_3_angle_beam;
  float lol_3_n_valid_tracks;
  float lol_3_min_angle;
  float lol_3_vtx_n_segs;
  float lol_3_shower_main_length;
  float lol_3_n_out;
  float lol_3_n_sum;    
  float lol_3_flag;

  
  sig->SetBranchAddress("hol_1_n_valid_tracks", &hol_1_n_valid_tracks);
  sig->SetBranchAddress("hol_1_min_angle", &hol_1_min_angle);
  sig->SetBranchAddress("hol_1_energy", &hol_1_energy);
  sig->SetBranchAddress("hol_1_flag_all_shower", &hol_1_flag_all_shower);
  sig->SetBranchAddress("hol_1_min_length", &hol_1_min_length);
  sig->SetBranchAddress("hol_1_flag", &hol_1_flag);
  
  sig->SetBranchAddress("hol_2_min_angle", &hol_2_min_angle);
  sig->SetBranchAddress("hol_2_medium_dQ_dx", &hol_2_medium_dQ_dx);
  sig->SetBranchAddress("hol_2_ncount", &hol_2_ncount);
  sig->SetBranchAddress("hol_2_flag", &hol_2_flag);

  sig->SetBranchAddress("lol_3_angle_beam",&lol_3_angle_beam);
  sig->SetBranchAddress("lol_3_n_valid_tracks",&lol_3_n_valid_tracks);
  sig->SetBranchAddress("lol_3_min_angle",&lol_3_min_angle);
  sig->SetBranchAddress("lol_3_vtx_n_segs",&lol_3_vtx_n_segs);
  sig->SetBranchAddress("lol_3_shower_main_length",&lol_3_shower_main_length);
  sig->SetBranchAddress("lol_3_n_out",&lol_3_n_out);
  sig->SetBranchAddress("lol_3_n_sum",&lol_3_n_sum);
  sig->SetBranchAddress("lol_3_flag",&lol_3_flag);

  bkg->SetBranchAddress("hol_1_n_valid_tracks", &hol_1_n_valid_tracks);
  bkg->SetBranchAddress("hol_1_min_angle", &hol_1_min_angle);
  bkg->SetBranchAddress("hol_1_energy", &hol_1_energy);
  bkg->SetBranchAddress("hol_1_flag_all_shower", &hol_1_flag_all_shower);
  bkg->SetBranchAddress("hol_1_min_length", &hol_1_min_length);
  bkg->SetBranchAddress("hol_1_flag", &hol_1_flag);
  
  bkg->SetBranchAddress("hol_2_min_angle", &hol_2_min_angle);
  bkg->SetBranchAddress("hol_2_medium_dQ_dx", &hol_2_medium_dQ_dx);
  bkg->SetBranchAddress("hol_2_ncount", &hol_2_ncount);
  bkg->SetBranchAddress("hol_2_flag", &hol_2_flag);

  bkg->SetBranchAddress("lol_3_angle_beam",&lol_3_angle_beam);
  bkg->SetBranchAddress("lol_3_n_valid_tracks",&lol_3_n_valid_tracks);
  bkg->SetBranchAddress("lol_3_min_angle",&lol_3_min_angle);
  bkg->SetBranchAddress("lol_3_vtx_n_segs",&lol_3_vtx_n_segs);
  bkg->SetBranchAddress("lol_3_shower_main_length",&lol_3_shower_main_length);
  bkg->SetBranchAddress("lol_3_n_out",&lol_3_n_out);
  bkg->SetBranchAddress("lol_3_n_sum",&lol_3_n_sum);
  bkg->SetBranchAddress("lol_3_flag",&lol_3_flag);
   
  
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

  
  Tsig->Branch("hol_1_n_valid_tracks", &hol_1_n_valid_tracks,"hol_1_n_valid_tracks/F");
  Tsig->Branch("hol_1_min_angle", &hol_1_min_angle,"hol_1_min_angle/F");
  Tsig->Branch("hol_1_energy", &hol_1_energy,"hol_1_energy/F");
  Tsig->Branch("hol_1_flag_all_shower", &hol_1_flag_all_shower,"hol_1_flag_all_shower/F");
  Tsig->Branch("hol_1_min_length", &hol_1_min_length,"hol_1_min_length/F");
  Tsig->Branch("hol_1_flag", &hol_1_flag,"hol_1_flag/F");
  
  Tsig->Branch("hol_2_min_angle", &hol_2_min_angle,"hol_2_min_angle/F");
  Tsig->Branch("hol_2_medium_dQ_dx", &hol_2_medium_dQ_dx,"hol_2_medium_dQ_dx/F");
  Tsig->Branch("hol_2_ncount", &hol_2_ncount,"hol_2_ncount/F");
  Tsig->Branch("hol_2_flag", &hol_2_flag,"hol_2_flag/F");

  Tsig->Branch("lol_3_angle_beam",&lol_3_angle_beam,"lol_3_angle_beam/F");
  Tsig->Branch("lol_3_n_valid_tracks",&lol_3_n_valid_tracks,"lol_3_n_valid_tracks/F");
  Tsig->Branch("lol_3_min_angle",&lol_3_min_angle,"lol_3_min_angle/F");
  Tsig->Branch("lol_3_vtx_n_segs",&lol_3_vtx_n_segs,"lol_3_vtx_n_segs/F");
  Tsig->Branch("lol_3_shower_main_length",&lol_3_shower_main_length,"lol_3_shower_main_length/F");
  Tsig->Branch("lol_3_n_out",&lol_3_n_out,"lol_3_n_out/F");
  Tsig->Branch("lol_3_n_sum",&lol_3_n_sum,"lol_3_n_sum/F");
  Tsig->Branch("lol_3_flag",&lol_3_flag,"lol_3_flag/F");
  
  Tbkg->Branch("hol_1_n_valid_tracks", &hol_1_n_valid_tracks,"hol_1_n_valid_tracks/F");
  Tbkg->Branch("hol_1_min_angle", &hol_1_min_angle,"hol_1_min_angle/F");
  Tbkg->Branch("hol_1_energy", &hol_1_energy,"hol_1_energy/F");
  Tbkg->Branch("hol_1_flag_all_shower", &hol_1_flag_all_shower,"hol_1_flag_all_shower/F");
  Tbkg->Branch("hol_1_min_length", &hol_1_min_length,"hol_1_min_length/F");
  Tbkg->Branch("hol_1_flag", &hol_1_flag,"hol_1_flag/F");
  
  Tbkg->Branch("hol_2_min_angle", &hol_2_min_angle,"hol_2_min_angle/F");
  Tbkg->Branch("hol_2_medium_dQ_dx", &hol_2_medium_dQ_dx,"hol_2_medium_dQ_dx/F");
  Tbkg->Branch("hol_2_ncount", &hol_2_ncount,"hol_2_ncount/F");
  Tbkg->Branch("hol_2_flag", &hol_2_flag,"hol_2_flag/F");

  Tbkg->Branch("lol_3_angle_beam",&lol_3_angle_beam,"lol_3_angle_beam/F");
  Tbkg->Branch("lol_3_n_valid_tracks",&lol_3_n_valid_tracks,"lol_3_n_valid_tracks/F");
  Tbkg->Branch("lol_3_min_angle",&lol_3_min_angle,"lol_3_min_angle/F");
  Tbkg->Branch("lol_3_vtx_n_segs",&lol_3_vtx_n_segs,"lol_3_vtx_n_segs/F");
  Tbkg->Branch("lol_3_shower_main_length",&lol_3_shower_main_length,"lol_3_shower_main_length/F");
  Tbkg->Branch("lol_3_n_out",&lol_3_n_out,"lol_3_n_out/F");
  Tbkg->Branch("lol_3_n_sum",&lol_3_n_sum,"lol_3_n_sum/F");
  Tbkg->Branch("lol_3_flag",&lol_3_flag,"lol_3_flag/F");

  Tsig->Branch("hol_lol_bdt",&bdt_value,"data/F");
  Tbkg->Branch("hol_lol_bdt",&bdt_value,"data/F");
  
  TMVA::Reader *reader = new TMVA::Reader();

  reader->AddVariable("hol_1_n_valid_tracks",&hol_1_n_valid_tracks);
  reader->AddVariable("hol_1_min_angle",&hol_1_min_angle);
  reader->AddVariable("hol_1_energy",&hol_1_energy);
  reader->AddVariable("hol_1_flag_all_shower",&hol_1_flag_all_shower);
  reader->AddVariable("hol_1_min_length",&hol_1_min_length);

  reader->AddVariable("hol_2_min_angle",&hol_2_min_angle);
  reader->AddVariable("hol_2_medium_dQ_dx",&hol_2_medium_dQ_dx);
  reader->AddVariable("hol_2_ncount",&hol_2_ncount);
  
  reader->AddVariable("lol_3_angle_beam",&lol_3_angle_beam);
  reader->AddVariable("lol_3_n_valid_tracks",&lol_3_n_valid_tracks);
  reader->AddVariable("lol_3_min_angle",&lol_3_min_angle);
  reader->AddVariable("lol_3_vtx_n_segs",&lol_3_vtx_n_segs);
  reader->AddVariable("lol_3_shower_main_length",&lol_3_shower_main_length);
  reader->AddVariable("lol_3_n_out",&lol_3_n_out);
  reader->AddVariable("lol_3_n_sum",&lol_3_n_sum);
  
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

