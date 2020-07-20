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

  
   // compare with muon
  double cme_mu_energy;
  double cme_energy;
  double cme_mu_length;
  double cme_length;
  double cme_angle_beam;
  int cme_flag;

  double anc_angle;
  double anc_max_angle;
  double anc_max_length;
  double anc_acc_forward_length;
  double anc_acc_backward_length;
  double anc_acc_forward_length1;
  double anc_shower_main_length;
  double anc_shower_total_length;
  int anc_flag_main_outside;
  int anc_flag;
  
  sig->SetBranchAddress("cme_mu_energy",&cme_mu_energy);
  sig->SetBranchAddress("cme_energy",&cme_energy);
  sig->SetBranchAddress("cme_mu_length",&cme_mu_length);
  sig->SetBranchAddress("cme_length",&cme_length);
  sig->SetBranchAddress("cme_angle_beam",&cme_angle_beam);
  sig->SetBranchAddress("cme_flag",&cme_flag);

  sig->SetBranchAddress("anc_angle",&anc_angle);
  sig->SetBranchAddress("anc_max_angle",&anc_max_angle);
  sig->SetBranchAddress("anc_max_length",&anc_max_length);
  sig->SetBranchAddress("anc_acc_forward_length",&anc_acc_forward_length);
  sig->SetBranchAddress("anc_acc_backward_length",&anc_acc_backward_length);
  sig->SetBranchAddress("anc_acc_forward_length1",&anc_acc_forward_length1);
  sig->SetBranchAddress("anc_shower_main_length",&anc_shower_main_length);
  sig->SetBranchAddress("anc_shower_total_length",&anc_shower_total_length);
  sig->SetBranchAddress("anc_flag_main_outside",&anc_flag_main_outside);
  sig->SetBranchAddress("anc_flag",&anc_flag);

   bkg->SetBranchAddress("cme_mu_energy",&cme_mu_energy);
  bkg->SetBranchAddress("cme_energy",&cme_energy);
  bkg->SetBranchAddress("cme_mu_length",&cme_mu_length);
  bkg->SetBranchAddress("cme_length",&cme_length);
  bkg->SetBranchAddress("cme_angle_beam",&cme_angle_beam);
  bkg->SetBranchAddress("cme_flag",&cme_flag);

  bkg->SetBranchAddress("anc_angle",&anc_angle);
  bkg->SetBranchAddress("anc_max_angle",&anc_max_angle);
  bkg->SetBranchAddress("anc_max_length",&anc_max_length);
  bkg->SetBranchAddress("anc_acc_forward_length",&anc_acc_forward_length);
  bkg->SetBranchAddress("anc_acc_backward_length",&anc_acc_backward_length);
  bkg->SetBranchAddress("anc_acc_forward_length1",&anc_acc_forward_length1);
  bkg->SetBranchAddress("anc_shower_main_length",&anc_shower_main_length);
  bkg->SetBranchAddress("anc_shower_total_length",&anc_shower_total_length);
  bkg->SetBranchAddress("anc_flag_main_outside",&anc_flag_main_outside);
  bkg->SetBranchAddress("anc_flag",&anc_flag);
  
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
  
  float cme_mu_energy_f;
  float cme_energy_f;
  float cme_mu_length_f;
  float cme_length_f;
  float cme_angle_beam_f;
  float cme_flag_f;

  float anc_angle_f;
  float anc_max_angle_f;
  float anc_max_length_f;
  float anc_acc_forward_length_f;
  float anc_acc_backward_length_f;
  float anc_acc_forward_length1_f;
  float anc_shower_main_length_f;
  float anc_shower_total_length_f;
  float anc_flag_main_outside_f;
  float anc_flag_f;

  
  Tsig->Branch("cme_mu_energy",&cme_mu_energy_f,"cme_mu_energy/F");
  Tsig->Branch("cme_energy",&cme_energy_f,"cme_energy/F");
  Tsig->Branch("cme_mu_length",&cme_mu_length_f,"cme_mu_length/F");
  Tsig->Branch("cme_length",&cme_length_f,"cme_length/F");
  Tsig->Branch("cme_angle_beam",&cme_angle_beam_f,"cme_angle_beam/F");
  Tsig->Branch("cme_flag",&cme_flag_f,"cme_flag/F");
    
  Tsig->Branch("anc_angle",&anc_angle_f,"anc_angle/F");
  Tsig->Branch("anc_max_angle",&anc_max_angle_f,"anc_max_angle/F");
  Tsig->Branch("anc_max_length",&anc_max_length_f,"anc_max_length/F");
  Tsig->Branch("anc_acc_forward_length",&anc_acc_forward_length_f,"anc_acc_forward_length/F");
  Tsig->Branch("anc_acc_backward_length",&anc_acc_backward_length_f,"anc_acc_backward_length/F");
  Tsig->Branch("anc_acc_forward_length1",&anc_acc_forward_length1_f,"anc_acc_forward_length1/F");
  Tsig->Branch("anc_shower_main_length",&anc_shower_main_length_f,"anc_shower_main_length/F");
  Tsig->Branch("anc_shower_total_length",&anc_shower_total_length_f,"anc_shower_total_length/F");
  Tsig->Branch("anc_flag_main_outside",&anc_flag_main_outside_f,"anc_flag_main_outside/F");
  Tsig->Branch("anc_flag",&anc_flag_f,"anc_flag/F");

  Tbkg->Branch("cme_mu_energy",&cme_mu_energy_f,"cme_mu_energy/F");
  Tbkg->Branch("cme_energy",&cme_energy_f,"cme_energy/F");
  Tbkg->Branch("cme_mu_length",&cme_mu_length_f,"cme_mu_length/F");
  Tbkg->Branch("cme_length",&cme_length_f,"cme_length/F");
  Tbkg->Branch("cme_angle_beam",&cme_angle_beam_f,"cme_angle_beam/F");
  Tbkg->Branch("cme_flag",&cme_flag_f,"cme_flag/F");
    
  Tbkg->Branch("anc_angle",&anc_angle_f,"anc_angle/F");
  Tbkg->Branch("anc_max_angle",&anc_max_angle_f,"anc_max_angle/F");
  Tbkg->Branch("anc_max_length",&anc_max_length_f,"anc_max_length/F");
  Tbkg->Branch("anc_acc_forward_length",&anc_acc_forward_length_f,"anc_acc_forward_length/F");
  Tbkg->Branch("anc_acc_backward_length",&anc_acc_backward_length_f,"anc_acc_backward_length/F");
  Tbkg->Branch("anc_acc_forward_length1",&anc_acc_forward_length1_f,"anc_acc_forward_length1/F");
  Tbkg->Branch("anc_shower_main_length",&anc_shower_main_length_f,"anc_shower_main_length/F");
  Tbkg->Branch("anc_shower_total_length",&anc_shower_total_length_f,"anc_shower_total_length/F");
  Tbkg->Branch("anc_flag_main_outside",&anc_flag_main_outside_f,"anc_flag_main_outside/F");
  Tbkg->Branch("anc_flag",&anc_flag_f,"anc_flag/F");
  
  
  for (Int_t i=0;i!=sig->GetEntries();i++){
    sig->GetEntry(i);
    
    cme_mu_energy_f = cme_mu_energy;
    cme_energy_f = cme_energy;
    cme_mu_length_f = cme_mu_length;
    cme_length_f = cme_length;
    cme_angle_beam_f = cme_angle_beam;
    cme_flag_f = cme_flag;
    
    anc_angle_f = anc_angle;
    anc_max_angle_f = anc_max_angle;
    anc_max_length_f = anc_max_length;
    anc_acc_forward_length_f = anc_acc_forward_length;
    anc_acc_backward_length_f = anc_acc_backward_length;
    anc_acc_forward_length1_f = anc_acc_forward_length1;
    anc_shower_main_length_f = anc_shower_main_length;
    anc_shower_total_length_f = anc_shower_total_length;
    anc_flag_main_outside_f = anc_flag_main_outside;
    anc_flag_f = anc_flag;
    
    Tsig->Fill();
  }

   for (Int_t i=0;i!=bkg->GetEntries();i++){
    bkg->GetEntry(i);
    
    cme_mu_energy_f = cme_mu_energy;
    cme_energy_f = cme_energy;
    cme_mu_length_f = cme_mu_length;
    cme_length_f = cme_length;
    cme_angle_beam_f = cme_angle_beam;
    cme_flag_f = cme_flag;
    
    anc_angle_f = anc_angle;
    anc_max_angle_f = anc_max_angle;
    anc_max_length_f = anc_max_length;
    anc_acc_forward_length_f = anc_acc_forward_length;
    anc_acc_backward_length_f = anc_acc_backward_length;
    anc_acc_forward_length1_f = anc_acc_forward_length1;
    anc_shower_main_length_f = anc_shower_main_length;
    anc_shower_total_length_f = anc_shower_total_length;
    anc_flag_main_outside_f = anc_flag_main_outside;
    anc_flag_f = anc_flag;
    
    
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
   
    
    dataloader->AddVariable("cme_mu_energy","cme_mu_energy","MeV",'F');
    dataloader->AddVariable("cme_energy","cme_energy","MeV",'F');
    dataloader->AddVariable("cme_mu_length","cme_mu_length","cm",'F');
    dataloader->AddVariable("cme_length","cme_length","cm",'F');
    dataloader->AddVariable("cme_angle_beam","cme_angle_beam","deg",'F');

    dataloader->AddVariable("anc_angle","anc_angle","deg",'F');
    dataloader->AddVariable("anc_max_angle","anc_max_angle","deg",'F');
    dataloader->AddVariable("anc_max_length","anc_max_length","cm",'F');
    dataloader->AddVariable("anc_acc_forward_length","anc_acc_forward_length","cm",'F');
    dataloader->AddVariable("anc_acc_backward_length","anc_acc_backward_length","cm",'F');
    dataloader->AddVariable("anc_acc_forward_length1","anc_acc_forward_length1","cm",'F');
    dataloader->AddVariable("anc_shower_main_length","anc_shower_main_length","cm",'F');
    dataloader->AddVariable("anc_shower_total_length","anc_shower_total_length","cm",'F');
    dataloader->AddVariable("anc_flag_main_outside","anc_flag_main_outside","",'F');
    
    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "1>0"; //  3394 (cme)  2231 (anc) or 5314 (total)/44194
    TCut mycut_b = "cme_flag==0 || anc_flag==0"; //  8985 (cme) 4279 (anc) or 11677 (total) / 21070
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=20000:"
					    "nTrain_Background=10000:"
					    "nTest_Signal=10000:"
					    "nTest_Background=1677:"
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
   
    
    dataloader->AddVariable("cme_mu_energy","cme_mu_energy","MeV",'F');
    dataloader->AddVariable("cme_energy","cme_energy","MeV",'F');
    dataloader->AddVariable("cme_mu_length","cme_mu_length","cm",'F');
    dataloader->AddVariable("cme_length","cme_length","cm",'F');
    dataloader->AddVariable("cme_angle_beam","cme_angle_beam","deg",'F');

    dataloader->AddVariable("anc_angle","anc_angle","deg",'F');
    dataloader->AddVariable("anc_max_angle","anc_max_angle","deg",'F');
    dataloader->AddVariable("anc_max_length","anc_max_length","cm",'F');
    dataloader->AddVariable("anc_acc_forward_length","anc_acc_forward_length","cm",'F');
    dataloader->AddVariable("anc_acc_backward_length","anc_acc_backward_length","cm",'F');
    dataloader->AddVariable("anc_acc_forward_length1","anc_acc_forward_length1","cm",'F');
    dataloader->AddVariable("anc_shower_main_length","anc_shower_main_length","cm",'F');
    dataloader->AddVariable("anc_shower_total_length","anc_shower_total_length","cm",'F');
    dataloader->AddVariable("anc_flag_main_outside","anc_flag_main_outside","",'F');
    
    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "1>0"; //  
    TCut mycut_b = "cme_flag==0 || anc_flag==0 || cme_anc_bdt < 0"; // 12981
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=20000:"
					    "nTrain_Background=12000:"
					    "nTest_Signal=10000:"
					    "nTest_Background=981:"
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

  float cme_mu_energy;
  float cme_energy;
  float cme_mu_length;
  float cme_length;
  float cme_angle_beam;
  float cme_flag;

  float anc_angle;
  float anc_max_angle;
  float anc_max_length;
  float anc_acc_forward_length;
  float anc_acc_backward_length;
  float anc_acc_forward_length1;
  float anc_shower_main_length;
  float anc_shower_total_length;
  float anc_flag_main_outside;
  float anc_flag;

  sig->SetBranchAddress("cme_mu_energy",&cme_mu_energy);
  sig->SetBranchAddress("cme_energy",&cme_energy);
  sig->SetBranchAddress("cme_mu_length",&cme_mu_length);
  sig->SetBranchAddress("cme_length",&cme_length);
  sig->SetBranchAddress("cme_angle_beam",&cme_angle_beam);
  sig->SetBranchAddress("cme_flag",&cme_flag);

  sig->SetBranchAddress("anc_angle",&anc_angle);
  sig->SetBranchAddress("anc_max_angle",&anc_max_angle);
  sig->SetBranchAddress("anc_max_length",&anc_max_length);
  sig->SetBranchAddress("anc_acc_forward_length",&anc_acc_forward_length);
  sig->SetBranchAddress("anc_acc_backward_length",&anc_acc_backward_length);
  sig->SetBranchAddress("anc_acc_forward_length1",&anc_acc_forward_length1);
  sig->SetBranchAddress("anc_shower_main_length",&anc_shower_main_length);
  sig->SetBranchAddress("anc_shower_total_length",&anc_shower_total_length);
  sig->SetBranchAddress("anc_flag_main_outside",&anc_flag_main_outside);
  sig->SetBranchAddress("anc_flag",&anc_flag);

   bkg->SetBranchAddress("cme_mu_energy",&cme_mu_energy);
  bkg->SetBranchAddress("cme_energy",&cme_energy);
  bkg->SetBranchAddress("cme_mu_length",&cme_mu_length);
  bkg->SetBranchAddress("cme_length",&cme_length);
  bkg->SetBranchAddress("cme_angle_beam",&cme_angle_beam);
  bkg->SetBranchAddress("cme_flag",&cme_flag);

  bkg->SetBranchAddress("anc_angle",&anc_angle);
  bkg->SetBranchAddress("anc_max_angle",&anc_max_angle);
  bkg->SetBranchAddress("anc_max_length",&anc_max_length);
  bkg->SetBranchAddress("anc_acc_forward_length",&anc_acc_forward_length);
  bkg->SetBranchAddress("anc_acc_backward_length",&anc_acc_backward_length);
  bkg->SetBranchAddress("anc_acc_forward_length1",&anc_acc_forward_length1);
  bkg->SetBranchAddress("anc_shower_main_length",&anc_shower_main_length);
  bkg->SetBranchAddress("anc_shower_total_length",&anc_shower_total_length);
  bkg->SetBranchAddress("anc_flag_main_outside",&anc_flag_main_outside);
  bkg->SetBranchAddress("anc_flag",&anc_flag);
   
  
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
  
  
  Tsig->Branch("cme_mu_energy",&cme_mu_energy,"cme_mu_energy/F");
  Tsig->Branch("cme_energy",&cme_energy,"cme_energy/F");
  Tsig->Branch("cme_mu_length",&cme_mu_length,"cme_mu_length/F");
  Tsig->Branch("cme_length",&cme_length,"cme_length/F");
  Tsig->Branch("cme_angle_beam",&cme_angle_beam,"cme_angle_beam/F");
  Tsig->Branch("cme_flag",&cme_flag,"cme_flag/F");
    
  Tsig->Branch("anc_angle",&anc_angle,"anc_angle/F");
  Tsig->Branch("anc_max_angle",&anc_max_angle,"anc_max_angle/F");
  Tsig->Branch("anc_max_length",&anc_max_length,"anc_max_length/F");
  Tsig->Branch("anc_acc_forward_length",&anc_acc_forward_length,"anc_acc_forward_length/F");
  Tsig->Branch("anc_acc_backward_length",&anc_acc_backward_length,"anc_acc_backward_length/F");
  Tsig->Branch("anc_acc_forward_length1",&anc_acc_forward_length1,"anc_acc_forward_length1/F");
  Tsig->Branch("anc_shower_main_length",&anc_shower_main_length,"anc_shower_main_length/F");
  Tsig->Branch("anc_shower_total_length",&anc_shower_total_length,"anc_shower_total_length/F");
  Tsig->Branch("anc_flag_main_outside",&anc_flag_main_outside,"anc_flag_main_outside/F");
  Tsig->Branch("anc_flag",&anc_flag,"anc_flag/F");

  Tbkg->Branch("cme_mu_energy",&cme_mu_energy,"cme_mu_energy/F");
  Tbkg->Branch("cme_energy",&cme_energy,"cme_energy/F");
  Tbkg->Branch("cme_mu_length",&cme_mu_length,"cme_mu_length/F");
  Tbkg->Branch("cme_length",&cme_length,"cme_length/F");
  Tbkg->Branch("cme_angle_beam",&cme_angle_beam,"cme_angle_beam/F");
  Tbkg->Branch("cme_flag",&cme_flag,"cme_flag/F");
    
  Tbkg->Branch("anc_angle",&anc_angle,"anc_angle/F");
  Tbkg->Branch("anc_max_angle",&anc_max_angle,"anc_max_angle/F");
  Tbkg->Branch("anc_max_length",&anc_max_length,"anc_max_length/F");
  Tbkg->Branch("anc_acc_forward_length",&anc_acc_forward_length,"anc_acc_forward_length/F");
  Tbkg->Branch("anc_acc_backward_length",&anc_acc_backward_length,"anc_acc_backward_length/F");
  Tbkg->Branch("anc_acc_forward_length1",&anc_acc_forward_length1,"anc_acc_forward_length1/F");
  Tbkg->Branch("anc_shower_main_length",&anc_shower_main_length,"anc_shower_main_length/F");
  Tbkg->Branch("anc_shower_total_length",&anc_shower_total_length,"anc_shower_total_length/F");
  Tbkg->Branch("anc_flag_main_outside",&anc_flag_main_outside,"anc_flag_main_outside/F");
  Tbkg->Branch("anc_flag",&anc_flag,"anc_flag/F");
  
  
  Tsig->Branch("cme_anc_bdt",&bdt_value,"data/F");
  Tbkg->Branch("cme_anc_bdt",&bdt_value,"data/F");
  
  TMVA::Reader *reader = new TMVA::Reader();
  reader->AddVariable("cme_mu_energy",&cme_mu_energy);
  reader->AddVariable("cme_energy",&cme_energy);
  reader->AddVariable("cme_mu_length",&cme_mu_length);
  reader->AddVariable("cme_length",&cme_length);
  reader->AddVariable("cme_angle_beam",&cme_angle_beam);

  reader->AddVariable("anc_angle",&anc_angle);
  reader->AddVariable("anc_max_angle",&anc_max_angle);
  reader->AddVariable("anc_max_length",&anc_max_length);
  reader->AddVariable("anc_acc_forward_length",&anc_acc_forward_length);
  reader->AddVariable("anc_acc_backward_length",&anc_acc_backward_length);
  reader->AddVariable("anc_acc_forward_length1",&anc_acc_forward_length1);
  reader->AddVariable("anc_shower_main_length",&anc_shower_main_length);
  reader->AddVariable("anc_shower_total_length",&anc_shower_total_length);
  reader->AddVariable("anc_flag_main_outside",&anc_flag_main_outside);
  
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

