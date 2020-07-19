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

  double stem_len_energy;
  double stem_len_length;
  int stem_len_flag_avoid_muon_check;
  int stem_len_num_daughters;
  double stem_len_daughter_length;
  int stem_len_flag;
  
  int brm_n_mu_segs;
  double brm_Ep;
  double brm_acc_length;
  double brm_shower_total_length;
  double brm_connected_length;
  int brm_n_size;
  double brm_acc_direct_length;
  int brm_n_shower_main_segs;
  int brm_n_mu_main;
  int brm_flag;

   // low-energy michel
  double lem_shower_main_length;
  int lem_n_3seg;
  double lem_e_charge;
  double lem_e_dQdx;
  int lem_shower_num_main_segs;
  int lem_flag;
  
  sig->SetBranchAddress("stem_len_energy", &stem_len_energy);
  sig->SetBranchAddress("stem_len_length", &stem_len_length);
  sig->SetBranchAddress("stem_len_flag_avoid_muon_check", &stem_len_flag_avoid_muon_check);
  sig->SetBranchAddress("stem_len_num_daughters", &stem_len_num_daughters);
  sig->SetBranchAddress("stem_len_daughter_length", &stem_len_daughter_length);
  sig->SetBranchAddress("stem_len_flag", &stem_len_flag);

  sig->SetBranchAddress("brm_n_mu_segs",&brm_n_mu_segs);
  sig->SetBranchAddress("brm_Ep",&brm_Ep);
  sig->SetBranchAddress("brm_acc_length",&brm_acc_length);
  sig->SetBranchAddress("brm_shower_total_length",&brm_shower_total_length);
  sig->SetBranchAddress("brm_connected_length",&brm_connected_length);
  sig->SetBranchAddress("brm_n_size",&brm_n_size);
  sig->SetBranchAddress("brm_nacc_direct_length",&brm_acc_direct_length);
  sig->SetBranchAddress("brm_n_shower_main_segs",&brm_n_shower_main_segs);
  sig->SetBranchAddress("brm_n_mu_main",&brm_n_mu_main);
  sig->SetBranchAddress("brm_n_flag",&brm_flag);
  
  sig->SetBranchAddress("lem_shower_main_length",&lem_shower_main_length);
  sig->SetBranchAddress("lem_n_3seg",&lem_n_3seg);
  sig->SetBranchAddress("lem_e_charge",&lem_e_charge);
  sig->SetBranchAddress("lem_e_dQdx",&lem_e_dQdx);
  sig->SetBranchAddress("lem_shower_num_main_segs",&lem_shower_num_main_segs);
  sig->SetBranchAddress("lem_flag",&lem_flag);
  
  
  bkg->SetBranchAddress("stem_len_energy", &stem_len_energy);
  bkg->SetBranchAddress("stem_len_length", &stem_len_length);
  bkg->SetBranchAddress("stem_len_flag_avoid_muon_check", &stem_len_flag_avoid_muon_check);
  bkg->SetBranchAddress("stem_len_num_daughters", &stem_len_num_daughters);
  bkg->SetBranchAddress("stem_len_daughter_length", &stem_len_daughter_length);
  bkg->SetBranchAddress("stem_len_flag", &stem_len_flag);

  bkg->SetBranchAddress("brm_n_mu_segs",&brm_n_mu_segs);
  bkg->SetBranchAddress("brm_Ep",&brm_Ep);
  bkg->SetBranchAddress("brm_acc_length",&brm_acc_length);
  bkg->SetBranchAddress("brm_shower_total_length",&brm_shower_total_length);
  bkg->SetBranchAddress("brm_connected_length",&brm_connected_length);
  bkg->SetBranchAddress("brm_n_size",&brm_n_size);
  bkg->SetBranchAddress("brm_nacc_direct_length",&brm_acc_direct_length);
  bkg->SetBranchAddress("brm_n_shower_main_segs",&brm_n_shower_main_segs);
  bkg->SetBranchAddress("brm_n_mu_main",&brm_n_mu_main);
  bkg->SetBranchAddress("brm_n_flag",&brm_flag);
  
  bkg->SetBranchAddress("lem_shower_main_length",&lem_shower_main_length);
  bkg->SetBranchAddress("lem_n_3seg",&lem_n_3seg);
  bkg->SetBranchAddress("lem_e_charge",&lem_e_charge);
  bkg->SetBranchAddress("lem_e_dQdx",&lem_e_dQdx);
  bkg->SetBranchAddress("lem_shower_num_main_segs",&lem_shower_num_main_segs);
  bkg->SetBranchAddress("lem_flag",&lem_flag);
 

  
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
  
 
  float stem_len_energy_f;
  float stem_len_length_f;
  float stem_len_flag_avoid_muon_check_f;
  float stem_len_num_daughters_f;
  float stem_len_daughter_length_f;
  float stem_len_flag_f;
  
  float brm_n_mu_segs_f;
  float brm_Ep_f;
  float brm_acc_length_f;
  float brm_shower_total_length_f;
  float brm_connected_length_f;
  float brm_n_size_f;
  float brm_acc_direct_length_f;
  float brm_n_shower_main_segs_f;
  float brm_n_mu_main_f;
  float brm_flag_f;

   // low-energy michel
  float lem_shower_main_length_f;
  float lem_n_3seg_f;
  float lem_e_charge_f;
  float lem_e_dQdx_f;
  float lem_shower_num_main_segs_f;
  float lem_flag_f;

  Tsig->Branch("stem_len_energy", &stem_len_energy_f, "stem_len_energy/F");
  Tsig->Branch("stem_len_length", &stem_len_length_f, "stem_len_length/F");
  Tsig->Branch("stem_len_flag_avoid_muon_check", &stem_len_flag_avoid_muon_check_f, "stem_len_flag_avoid_muon_check/F");
  Tsig->Branch("stem_len_num_daughters", &stem_len_num_daughters_f, "stem_len_num_daughters/F");
  Tsig->Branch("stem_len_daughter_length", &stem_len_daughter_length_f, "stem_len_daughter_length/F");
  Tsig->Branch("stem_len_flag", &stem_len_flag_f, "stem_len_flag/F");
  
  Tsig->Branch("brm_n_mu_segs",&brm_n_mu_segs_f,"brm_n_mu_segs/F");
  Tsig->Branch("brm_Ep",&brm_Ep_f,"brm_Ep/F");
  Tsig->Branch("brm_acc_length",&brm_acc_length_f,"brm_acc_length/F");
  Tsig->Branch("brm_shower_total_length",&brm_shower_total_length_f,"brm_shower_total_length/F");
  Tsig->Branch("brm_connected_length",&brm_connected_length_f,"brm_connected_length/F");
  Tsig->Branch("brm_n_size",&brm_n_size_f,"brm_n_size/F");
  Tsig->Branch("brm_acc_direct_length",&brm_acc_direct_length_f,"brm_acc_direct_length/F");
  Tsig->Branch("brm_n_shower_main_segs",&brm_n_shower_main_segs_f,"brm_n_shower_main_segs/F");
  Tsig->Branch("brm_n_mu_main",&brm_n_mu_main_f,"brm_n_mu_main/F");
  Tsig->Branch("brm_flag",&brm_flag_f,"brm_flag/F");
  
  Tsig->Branch("lem_shower_main_length",&lem_shower_main_length_f,"lem_shower_main_length/F");
  Tsig->Branch("lem_n_3seg",&lem_n_3seg_f,"lem_n_3seg/F");
  Tsig->Branch("lem_e_charge",&lem_e_charge_f,"lem_e_charge/F");
  Tsig->Branch("lem_e_dQdx",&lem_e_dQdx_f,"lem_e_dQdx/F");
  Tsig->Branch("lem_shower_num_main_segs",&lem_shower_num_main_segs_f,"lem_shower_num_main_segs/F");
  Tsig->Branch("lem_flag",&lem_flag_f,"lem_flag/F");

  Tbkg->Branch("stem_len_energy", &stem_len_energy_f, "stem_len_energy/F");
  Tbkg->Branch("stem_len_length", &stem_len_length_f, "stem_len_length/F");
  Tbkg->Branch("stem_len_flag_avoid_muon_check", &stem_len_flag_avoid_muon_check_f, "stem_len_flag_avoid_muon_check/F");
  Tbkg->Branch("stem_len_num_daughters", &stem_len_num_daughters_f, "stem_len_num_daughters/F");
  Tbkg->Branch("stem_len_daughter_length", &stem_len_daughter_length_f, "stem_len_daughter_length/F");
  Tbkg->Branch("stem_len_flag", &stem_len_flag_f, "stem_len_flag/F");
  
  Tbkg->Branch("brm_n_mu_segs",&brm_n_mu_segs_f,"brm_n_mu_segs/F");
  Tbkg->Branch("brm_Ep",&brm_Ep_f,"brm_Ep/F");
  Tbkg->Branch("brm_acc_length",&brm_acc_length_f,"brm_acc_length/F");
  Tbkg->Branch("brm_shower_total_length",&brm_shower_total_length_f,"brm_shower_total_length/F");
  Tbkg->Branch("brm_connected_length",&brm_connected_length_f,"brm_connected_length/F");
  Tbkg->Branch("brm_n_size",&brm_n_size_f,"brm_n_size/F");
  Tbkg->Branch("brm_acc_direct_length",&brm_acc_direct_length_f,"brm_acc_direct_length/F");
  Tbkg->Branch("brm_n_shower_main_segs",&brm_n_shower_main_segs_f,"brm_n_shower_main_segs/F");
  Tbkg->Branch("brm_n_mu_main",&brm_n_mu_main_f,"brm_n_mu_main/F");
  Tbkg->Branch("brm_flag",&brm_flag_f,"brm_flag/F");
  
  Tbkg->Branch("lem_shower_main_length",&lem_shower_main_length_f,"lem_shower_main_length/F");
  Tbkg->Branch("lem_n_3seg",&lem_n_3seg_f,"lem_n_3seg/F");
  Tbkg->Branch("lem_e_charge",&lem_e_charge_f,"lem_e_charge/F");
  Tbkg->Branch("lem_e_dQdx",&lem_e_dQdx_f,"lem_e_dQdx/F");
  Tbkg->Branch("lem_shower_num_main_segs",&lem_shower_num_main_segs_f,"lem_shower_num_main_segs/F");
  Tbkg->Branch("lem_flag",&lem_flag_f,"lem_flag/F");
  
  for (Int_t i=0;i!=sig->GetEntries();i++){
    sig->GetEntry(i);

    stem_len_energy_f = stem_len_energy;
    stem_len_length_f = stem_len_length;
    stem_len_flag_avoid_muon_check_f = stem_len_flag_avoid_muon_check;
    stem_len_num_daughters_f = stem_len_num_daughters;
    stem_len_daughter_length_f = stem_len_daughter_length;
    stem_len_flag_f = stem_len_flag;
    
    brm_n_mu_segs_f = brm_n_mu_segs;
    brm_Ep_f = brm_Ep;
    brm_acc_length_f = brm_acc_length;
    brm_shower_total_length_f = brm_shower_total_length;
    brm_connected_length_f = brm_connected_length;
    brm_n_size_f = brm_n_size;
    brm_acc_direct_length_f = brm_acc_direct_length;
    brm_n_shower_main_segs_f = brm_n_shower_main_segs;
    brm_n_mu_main_f = brm_n_mu_main;
    brm_flag_f = brm_flag;
    
    lem_shower_main_length_f = lem_shower_main_length;
    lem_n_3seg_f = lem_n_3seg;
    lem_e_charge_f = lem_e_charge;
    lem_e_dQdx_f = lem_e_dQdx;
    lem_shower_num_main_segs_f = lem_shower_num_main_segs;
    lem_flag_f = lem_flag;
    
    
    Tsig->Fill();
  }

   for (Int_t i=0;i!=bkg->GetEntries();i++){
    bkg->GetEntry(i);

    stem_len_energy_f = stem_len_energy;
    stem_len_length_f = stem_len_length;
    stem_len_flag_avoid_muon_check_f = stem_len_flag_avoid_muon_check;
    stem_len_num_daughters_f = stem_len_num_daughters;
    stem_len_daughter_length_f = stem_len_daughter_length;
    stem_len_flag_f = stem_len_flag;
    
    brm_n_mu_segs_f = brm_n_mu_segs;
    brm_Ep_f = brm_Ep;
    brm_acc_length_f = brm_acc_length;
    brm_shower_total_length_f = brm_shower_total_length;
    brm_connected_length_f = brm_connected_length;
    brm_n_size_f = brm_n_size;
    brm_acc_direct_length_f = brm_acc_direct_length;
    brm_n_shower_main_segs_f = brm_n_shower_main_segs;
    brm_n_mu_main_f = brm_n_mu_main;
    brm_flag_f = brm_flag;
    
    lem_shower_main_length_f = lem_shower_main_length;
    lem_n_3seg_f = lem_n_3seg;
    lem_e_charge_f = lem_e_charge;
    lem_e_dQdx_f = lem_e_dQdx;
    lem_shower_num_main_segs_f = lem_shower_num_main_segs;
    lem_flag_f = lem_flag;
    
    
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
    
    
    dataloader->AddVariable("stem_len_energy","stem_len_energy","MeV",'F');
    dataloader->AddVariable("stem_len_length","stem_len_length","cm",'F');
    dataloader->AddVariable("stem_len_flag_avoid_muon_check","stem_len_flag_avoid_muon_check","",'F');
    dataloader->AddVariable("stem_len_num_daughters","stem_len_num_daughters","",'F');
    dataloader->AddVariable("stem_len_daughter_length","stem_len_daughter_length","cm",'F');

    dataloader->AddVariable("brm_n_mu_segs","brm_n_mu_segs","",'F');
    dataloader->AddVariable("brm_Ep","brm_Ep","MeV",'F');
    dataloader->AddVariable("brm_acc_length","brm_acc_length","cm",'F');
    dataloader->AddVariable("brm_shower_total_length","brm_shower_total_length","cm",'F');
    dataloader->AddVariable("brm_connected_length","brm_connected_length","cm",'F');
    dataloader->AddVariable("brm_n_size","brm_n_size","",'F');
    dataloader->AddVariable("brm_acc_direct_length","brm_acc_direct_length","cm",'F');
    dataloader->AddVariable("brm_n_shower_main_segs","brm_n_shower_main_segs","",'F');
    dataloader->AddVariable("brm_n_mu_main","brm_n_mu_main","",'F');

    dataloader->AddVariable("lem_shower_main_length","lem_shower_main_length","cm",'F');
    dataloader->AddVariable("lem_n_3seg","lem_n_3seg","",'F');
    dataloader->AddVariable("lem_e_charge","lem_e_charge","MeV",'F');
    dataloader->AddVariable("lem_e_dQdx","lem_e_dQdx","MeV",'F');
    dataloader->AddVariable("lem_shower_num_main_segs","lem_shower_num_main_segs","",'F');

    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "1>0"; // 72 (stl) + 315 (brm)  261 (lem) or 637 (total)/44194
    TCut mycut_b = "stem_len_flag ==0 || brm_flag == 0 || lem_flag == 0"; // 429 (stl) + 1635 (brm) 1708 (lem) or 3739 (total)/21070
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=20000:"
					    "nTrain_Background=3300:"
					    "nTest_Signal=10000:"
					    "nTest_Background=439:"
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
    
    
    dataloader->AddVariable("stem_len_energy","stem_len_energy","MeV",'F');
    dataloader->AddVariable("stem_len_length","stem_len_length","cm",'F');
    dataloader->AddVariable("stem_len_flag_avoid_muon_check","stem_len_flag_avoid_muon_check","",'F');
    dataloader->AddVariable("stem_len_num_daughters","stem_len_num_daughters","",'F');
    dataloader->AddVariable("stem_len_daughter_length","stem_len_daughter_length","cm",'F');

    dataloader->AddVariable("brm_n_mu_segs","brm_n_mu_segs","",'F');
    dataloader->AddVariable("brm_Ep","brm_Ep","MeV",'F');
    dataloader->AddVariable("brm_acc_length","brm_acc_length","cm",'F');
    dataloader->AddVariable("brm_shower_total_length","brm_shower_total_length","cm",'F');
    dataloader->AddVariable("brm_connected_length","brm_connected_length","cm",'F');
    dataloader->AddVariable("brm_n_size","brm_n_size","",'F');
    dataloader->AddVariable("brm_acc_direct_length","brm_acc_direct_length","cm",'F');
    dataloader->AddVariable("brm_n_shower_main_segs","brm_n_shower_main_segs","",'F');
    dataloader->AddVariable("brm_n_mu_main","brm_n_mu_main","",'F');

    dataloader->AddVariable("lem_shower_main_length","lem_shower_main_length","cm",'F');
    dataloader->AddVariable("lem_n_3seg","lem_n_3seg","",'F');
    dataloader->AddVariable("lem_e_charge","lem_e_charge","MeV",'F');
    dataloader->AddVariable("lem_e_dQdx","lem_e_dQdx","MeV",'F');
    dataloader->AddVariable("lem_shower_num_main_segs","lem_shower_num_main_segs","",'F');

    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "1>0"; // 72 (stl) + 315 (brm)  261 (lem) or 637 (total)/44194
    TCut mycut_b = "stem_len_flag ==0 || brm_flag == 0 || lem_flag == 0 || stl_lem_brm_bdt < -0.2"; // 429 (stl) + 1635 (brm) 1708 (lem) or 3918 (total)/21070
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=20000:"
					    "nTrain_Background=3500:"
					    "nTest_Signal=10000:"
					    "nTest_Background=418:"
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

  float stem_len_energy;
  float stem_len_length;
  float stem_len_flag_avoid_muon_check;
  float stem_len_num_daughters;
  float stem_len_daughter_length;
  float stem_len_flag;
  
  float brm_n_mu_segs;
  float brm_Ep;
  float brm_acc_length;
  float brm_shower_total_length;
  float brm_connected_length;
  float brm_n_size;
  float brm_acc_direct_length;
  float brm_n_shower_main_segs;
  float brm_n_mu_main;
  float brm_flag;
  
  float lem_shower_main_length;
  float lem_n_3seg;
  float lem_e_charge;
  float lem_e_dQdx;
  float lem_shower_num_main_segs;
  float lem_flag;

  sig->SetBranchAddress("stem_len_energy", &stem_len_energy);
  sig->SetBranchAddress("stem_len_length", &stem_len_length);
  sig->SetBranchAddress("stem_len_flag_avoid_muon_check", &stem_len_flag_avoid_muon_check);
  sig->SetBranchAddress("stem_len_num_daughters", &stem_len_num_daughters);
  sig->SetBranchAddress("stem_len_daughter_length", &stem_len_daughter_length);
  sig->SetBranchAddress("stem_len_flag", &stem_len_flag);

  sig->SetBranchAddress("brm_n_mu_segs",&brm_n_mu_segs);
  sig->SetBranchAddress("brm_Ep",&brm_Ep);
  sig->SetBranchAddress("brm_acc_length",&brm_acc_length);
  sig->SetBranchAddress("brm_shower_total_length",&brm_shower_total_length);
  sig->SetBranchAddress("brm_connected_length",&brm_connected_length);
  sig->SetBranchAddress("brm_n_size",&brm_n_size);
  sig->SetBranchAddress("brm_acc_direct_length",&brm_acc_direct_length);
  sig->SetBranchAddress("brm_n_shower_main_segs",&brm_n_shower_main_segs);
  sig->SetBranchAddress("brm_n_mu_main",&brm_n_mu_main);
  sig->SetBranchAddress("brm_flag",&brm_flag);
  
  sig->SetBranchAddress("lem_shower_main_length",&lem_shower_main_length);
  sig->SetBranchAddress("lem_n_3seg",&lem_n_3seg);
  sig->SetBranchAddress("lem_e_charge",&lem_e_charge);
  sig->SetBranchAddress("lem_e_dQdx",&lem_e_dQdx);
  sig->SetBranchAddress("lem_shower_num_main_segs",&lem_shower_num_main_segs);
  sig->SetBranchAddress("lem_flag",&lem_flag);
  
  
  bkg->SetBranchAddress("stem_len_energy", &stem_len_energy);
  bkg->SetBranchAddress("stem_len_length", &stem_len_length);
  bkg->SetBranchAddress("stem_len_flag_avoid_muon_check", &stem_len_flag_avoid_muon_check);
  bkg->SetBranchAddress("stem_len_num_daughters", &stem_len_num_daughters);
  bkg->SetBranchAddress("stem_len_daughter_length", &stem_len_daughter_length);
  bkg->SetBranchAddress("stem_len_flag", &stem_len_flag);

  bkg->SetBranchAddress("brm_n_mu_segs",&brm_n_mu_segs);
  bkg->SetBranchAddress("brm_Ep",&brm_Ep);
  bkg->SetBranchAddress("brm_acc_length",&brm_acc_length);
  bkg->SetBranchAddress("brm_shower_total_length",&brm_shower_total_length);
  bkg->SetBranchAddress("brm_connected_length",&brm_connected_length);
  bkg->SetBranchAddress("brm_n_size",&brm_n_size);
  bkg->SetBranchAddress("brm_acc_direct_length",&brm_acc_direct_length);
  bkg->SetBranchAddress("brm_n_shower_main_segs",&brm_n_shower_main_segs);
  bkg->SetBranchAddress("brm_n_mu_main",&brm_n_mu_main);
  bkg->SetBranchAddress("brm_flag",&brm_flag);
  
  bkg->SetBranchAddress("lem_shower_main_length",&lem_shower_main_length);
  bkg->SetBranchAddress("lem_n_3seg",&lem_n_3seg);
  bkg->SetBranchAddress("lem_e_charge",&lem_e_charge);
  bkg->SetBranchAddress("lem_e_dQdx",&lem_e_dQdx);
  bkg->SetBranchAddress("lem_shower_num_main_segs",&lem_shower_num_main_segs);
  bkg->SetBranchAddress("lem_flag",&lem_flag);
 

  
  
  
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
  
  
  Tsig->Branch("stem_len_energy", &stem_len_energy, "stem_len_energy/F");
  Tsig->Branch("stem_len_length", &stem_len_length, "stem_len_length/F");
  Tsig->Branch("stem_len_flag_avoid_muon_check", &stem_len_flag_avoid_muon_check, "stem_len_flag_avoid_muon_check/F");
  Tsig->Branch("stem_len_num_daughters", &stem_len_num_daughters, "stem_len_num_daughters/F");
  Tsig->Branch("stem_len_daughter_length", &stem_len_daughter_length, "stem_len_daughter_length/F");
  Tsig->Branch("stem_len_flag", &stem_len_flag, "stem_len_flag/F");
  
  Tsig->Branch("brm_n_mu_segs",&brm_n_mu_segs,"brm_n_mu_segs/F");
  Tsig->Branch("brm_Ep",&brm_Ep,"brm_Ep/F");
  Tsig->Branch("brm_acc_length",&brm_acc_length,"brm_acc_length/F");
  Tsig->Branch("brm_shower_total_length",&brm_shower_total_length,"brm_shower_total_length/F");
  Tsig->Branch("brm_connected_length",&brm_connected_length,"brm_connected_length/F");
  Tsig->Branch("brm_n_size",&brm_n_size,"brm_n_size/F");
  Tsig->Branch("brm_acc_direct_length",&brm_acc_direct_length,"brm_acc_direct_length/F");
  Tsig->Branch("brm_n_shower_main_segs",&brm_n_shower_main_segs,"brm_n_shower_main_segs/F");
  Tsig->Branch("brm_n_mu_main",&brm_n_mu_main,"brm_n_mu_main/F");
  Tsig->Branch("brm_flag",&brm_flag,"brm_flag/F");
  
  Tsig->Branch("lem_shower_main_length",&lem_shower_main_length,"lem_shower_main_length/F");
  Tsig->Branch("lem_n_3seg",&lem_n_3seg,"lem_n_3seg/F");
  Tsig->Branch("lem_e_charge",&lem_e_charge,"lem_e_charge/F");
  Tsig->Branch("lem_e_dQdx",&lem_e_dQdx,"lem_e_dQdx/F");
  Tsig->Branch("lem_shower_num_main_segs",&lem_shower_num_main_segs,"lem_shower_num_main_segs/F");
  Tsig->Branch("lem_flag",&lem_flag,"lem_flag/F");

  Tbkg->Branch("stem_len_energy", &stem_len_energy, "stem_len_energy/F");
  Tbkg->Branch("stem_len_length", &stem_len_length, "stem_len_length/F");
  Tbkg->Branch("stem_len_flag_avoid_muon_check", &stem_len_flag_avoid_muon_check, "stem_len_flag_avoid_muon_check/F");
  Tbkg->Branch("stem_len_num_daughters", &stem_len_num_daughters, "stem_len_num_daughters/F");
  Tbkg->Branch("stem_len_daughter_length", &stem_len_daughter_length, "stem_len_daughter_length/F");
  Tbkg->Branch("stem_len_flag", &stem_len_flag, "stem_len_flag/F");
  
  Tbkg->Branch("brm_n_mu_segs",&brm_n_mu_segs,"brm_n_mu_segs/F");
  Tbkg->Branch("brm_Ep",&brm_Ep,"brm_Ep/F");
  Tbkg->Branch("brm_acc_length",&brm_acc_length,"brm_acc_length/F");
  Tbkg->Branch("brm_shower_total_length",&brm_shower_total_length,"brm_shower_total_length/F");
  Tbkg->Branch("brm_connected_length",&brm_connected_length,"brm_connected_length/F");
  Tbkg->Branch("brm_n_size",&brm_n_size,"brm_n_size/F");
  Tbkg->Branch("brm_acc_direct_length",&brm_acc_direct_length,"brm_acc_direct_length/F");
  Tbkg->Branch("brm_n_shower_main_segs",&brm_n_shower_main_segs,"brm_n_shower_main_segs/F");
  Tbkg->Branch("brm_n_mu_main",&brm_n_mu_main,"brm_n_mu_main/F");
  Tbkg->Branch("brm_flag",&brm_flag,"brm_flag/F");
  
  Tbkg->Branch("lem_shower_main_length",&lem_shower_main_length,"lem_shower_main_length/F");
  Tbkg->Branch("lem_n_3seg",&lem_n_3seg,"lem_n_3seg/F");
  Tbkg->Branch("lem_e_charge",&lem_e_charge,"lem_e_charge/F");
  Tbkg->Branch("lem_e_dQdx",&lem_e_dQdx,"lem_e_dQdx/F");
  Tbkg->Branch("lem_shower_num_main_segs",&lem_shower_num_main_segs,"lem_shower_num_main_segs/F");
  Tbkg->Branch("lem_flag",&lem_flag,"lem_flag/F");

  Tsig->Branch("stl_lem_brm_bdt", &bdt_value,"data/F");
  Tbkg->Branch("stl_lem_brm_bdt", &bdt_value,"data/F");

  
  TMVA::Reader *reader = new TMVA::Reader();

  reader->AddVariable("stem_len_energy",&stem_len_energy);
  reader->AddVariable("stem_len_length",&stem_len_length);
  reader->AddVariable("stem_len_flag_avoid_muon_check",&stem_len_flag_avoid_muon_check);
  reader->AddVariable("stem_len_num_daughters",&stem_len_num_daughters);
  reader->AddVariable("stem_len_daughter_length",&stem_len_daughter_length);
  
  reader->AddVariable("brm_n_mu_segs",&brm_n_mu_segs);
  reader->AddVariable("brm_Ep",&brm_Ep);
  reader->AddVariable("brm_acc_length",&brm_acc_length);
  reader->AddVariable("brm_shower_total_length",&brm_shower_total_length);
  reader->AddVariable("brm_connected_length",&brm_connected_length);
  reader->AddVariable("brm_n_size",&brm_n_size);
  reader->AddVariable("brm_acc_direct_length",&brm_acc_direct_length);
  reader->AddVariable("brm_n_shower_main_segs",&brm_n_shower_main_segs);
  reader->AddVariable("brm_n_mu_main",&brm_n_mu_main);
  
  reader->AddVariable("lem_shower_main_length",&lem_shower_main_length);
  reader->AddVariable("lem_n_3seg",&lem_n_3seg);
  reader->AddVariable("lem_e_charge",&lem_e_charge);
  reader->AddVariable("lem_e_dQdx",&lem_e_dQdx);
  reader->AddVariable("lem_shower_num_main_segs",&lem_shower_num_main_segs);
  
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

