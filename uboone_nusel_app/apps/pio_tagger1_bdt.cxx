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
  TFile *file = new TFile("bdtfile_0718.root");
  TTree *sig = (TTree*)file->Get("sig");
  TTree *bkg = (TTree*)file->Get("bkg");

  int pio_flag;
  int pio_mip_id;
  int pio_filled;
  int pio_flag_pio;

  int pio_1_flag;
  double pio_1_mass;
  int pio_1_pio_type;
  double pio_1_energy_1;
  double pio_1_energy_2;
  double pio_1_dis_1;
  double pio_1_dis_2;

  sig->SetBranchAddress("pio_flag",&pio_flag);
  sig->SetBranchAddress("pio_mip_id",&pio_mip_id);
  sig->SetBranchAddress("pio_filled",&pio_filled);
  sig->SetBranchAddress("pio_flag_pio",&pio_flag_pio);
  
  sig->SetBranchAddress("pio_1_flag",&pio_1_flag);
  sig->SetBranchAddress("pio_1_mass",&pio_1_mass);
  sig->SetBranchAddress("pio_1_pio_type",&pio_1_pio_type);
  sig->SetBranchAddress("pio_1_energy_1",&pio_1_energy_1);
  sig->SetBranchAddress("pio_1_energy_2",&pio_1_energy_2);
  sig->SetBranchAddress("pio_1_dis_1",&pio_1_dis_1);
  sig->SetBranchAddress("pio_1_dis_2",&pio_1_dis_2);
  
  float trueEdep;
  float weight;
  float lowEweight;
  Int_t nueTag;

  sig->SetBranchAddress("trueEdep",&trueEdep);
  sig->SetBranchAddress("weight",&weight);
  sig->SetBranchAddress("lowEweight",&lowEweight);
  sig->SetBranchAddress("nueTag",&nueTag);
  
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
  

  bkg->SetBranchAddress("trueEdep",&trueEdep);
  bkg->SetBranchAddress("weight",&weight);
  bkg->SetBranchAddress("lowEweight",&lowEweight);
  bkg->SetBranchAddress("nueTag",&nueTag);
  
  bkg->SetBranchAddress("pio_flag",&pio_flag);
  bkg->SetBranchAddress("pio_mip_id",&pio_mip_id);
  bkg->SetBranchAddress("pio_filled",&pio_filled);
  bkg->SetBranchAddress("pio_flag_pio",&pio_flag_pio);
  
  bkg->SetBranchAddress("pio_1_flag",&pio_1_flag);
  bkg->SetBranchAddress("pio_1_mass",&pio_1_mass);
  bkg->SetBranchAddress("pio_1_pio_type",&pio_1_pio_type);
  bkg->SetBranchAddress("pio_1_energy_1",&pio_1_energy_1);
  bkg->SetBranchAddress("pio_1_energy_2",&pio_1_energy_2);
  bkg->SetBranchAddress("pio_1_dis_1",&pio_1_dis_1);
  bkg->SetBranchAddress("pio_1_dis_2",&pio_1_dis_2);

  TFile *new_file = new TFile("reduced.root","RECREATE");
  TTree *Tsig = new TTree("sig","sig");
  TTree *Tbkg = new TTree("bkg","bkg");
  Tsig->SetDirectory(new_file);
  Tbkg->SetDirectory(new_file);
  
  float pio_flag_f;
  float pio_mip_id_f;
  float pio_filled_f;
  float pio_flag_pio_f;

  float pio_1_flag_f;
  float pio_1_mass_f;
  float pio_1_pio_type_f;
  float pio_1_energy_1_f;
  float pio_1_energy_2_f;
  float pio_1_dis_1_f;
  float pio_1_dis_2_f;

  Tsig->Branch("trueEdep",&trueEdep,"data/F");
  Tsig->Branch("weight",&weight,"data/F");
  Tsig->Branch("lowEweight",&lowEweight,"data/F");
  Tsig->Branch("nueTag",&nueTag,"data/I");
  
  Tsig->Branch("pio_flag",&pio_flag_f,"pio_flag/F");
  Tsig->Branch("pio_mip_id",&pio_mip_id_f,"pio_mip_id/F");
  Tsig->Branch("pio_filled",&pio_filled_f,"pio_filled/F");
  Tsig->Branch("pio_flag_pio",&pio_flag_pio_f,"pio_flag_pio/F");

  Tsig->Branch("truth_inFV",&truth_inFV,"data/I");
  Tsig->Branch("truth_CC",&truth_CC,"data/I");
  Tsig->Branch("truth_nue",&truth_nue,"data/I");
  Tsig->Branch("truth_cosmic",&truth_cosmic,"data/I");

  Tbkg->Branch("truth_inFV",&truth_inFV,"data/I");
  Tbkg->Branch("truth_CC",&truth_CC,"data/I");
  Tbkg->Branch("truth_nue",&truth_nue,"data/I");
  Tbkg->Branch("truth_cosmic",&truth_cosmic,"data/I");

  
  Tsig->Branch("pio_1_flag",&pio_1_flag_f,"pio_1_flag/F");
  Tsig->Branch("pio_1_mass",&pio_1_mass_f,"pio_1_mass/F");
  Tsig->Branch("pio_1_pio_type",&pio_1_pio_type_f,"pio_1_pio_type/F");
  Tsig->Branch("pio_1_energy_1",&pio_1_energy_1_f,"pio_1_energy_1/F");
  Tsig->Branch("pio_1_energy_2",&pio_1_energy_2_f,"pio_1_energy_2/F");
  Tsig->Branch("pio_1_dis_1",&pio_1_dis_1_f,"pio_1_dis_1/F");
  Tsig->Branch("pio_1_dis_2",&pio_1_dis_2_f,"pio_1_dis_2/F");
  
  Tbkg->Branch("trueEdep",&trueEdep,"data/F");
  Tbkg->Branch("weight",&weight,"data/F");
  Tbkg->Branch("lowEweight",&lowEweight,"data/F");
  Tbkg->Branch("nueTag",&nueTag,"data/I");

  Tbkg->Branch("pio_flag",&pio_flag_f,"pio_flag/F");
  Tbkg->Branch("pio_mip_id",&pio_mip_id_f,"pio_mip_id/F");
  Tbkg->Branch("pio_filled",&pio_filled_f,"pio_filled/F");
  Tbkg->Branch("pio_flag_pio",&pio_flag_pio_f,"pio_flag_pio/F");
  
  Tbkg->Branch("pio_1_flag",&pio_1_flag_f,"pio_1_flag/F");
  Tbkg->Branch("pio_1_mass",&pio_1_mass_f,"pio_1_mass/F");
  Tbkg->Branch("pio_1_pio_type",&pio_1_pio_type_f,"pio_1_pio_type/F");
  Tbkg->Branch("pio_1_energy_1",&pio_1_energy_1_f,"pio_1_energy_1/F");
  Tbkg->Branch("pio_1_energy_2",&pio_1_energy_2_f,"pio_1_energy_2/F");
  Tbkg->Branch("pio_1_dis_1",&pio_1_dis_1_f,"pio_1_dis_1/F");
  Tbkg->Branch("pio_1_dis_2",&pio_1_dis_2_f,"pio_1_dis_2/F");
  
  for (Int_t i=0;i!=sig->GetEntries();i++){
    sig->GetEntry(i);

    
    pio_flag_f = pio_flag;
    pio_mip_id_f = pio_mip_id;
    pio_filled_f = pio_filled;
    pio_flag_pio_f = pio_flag_pio;

    pio_1_flag_f = pio_1_flag;
    pio_1_mass_f = pio_1_mass;
    pio_1_pio_type_f = pio_1_pio_type;
    pio_1_energy_1_f = pio_1_energy_1;
    pio_1_energy_2_f = pio_1_energy_2;
    pio_1_dis_1_f = pio_1_dis_1;
    pio_1_dis_2_f = pio_1_dis_2;
    
    
    Tsig->Fill();
  }

   for (Int_t i=0;i!=bkg->GetEntries();i++){
    bkg->GetEntry(i);
    
     pio_flag_f = pio_flag;
    pio_mip_id_f = pio_mip_id;
    pio_filled_f = pio_filled;
    pio_flag_pio_f = pio_flag_pio;

    pio_1_flag_f = pio_1_flag;
    pio_1_mass_f = pio_1_mass;
    pio_1_pio_type_f = pio_1_pio_type;
    pio_1_energy_1_f = pio_1_energy_1;
    pio_1_energy_2_f = pio_1_energy_2;
    pio_1_dis_1_f = pio_1_dis_1;
    pio_1_dis_2_f = pio_1_dis_2;
    
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

    dataloader->AddVariable("pio_1_mass","pio_1_mass","MeV2",'F');
    dataloader->AddVariable("pio_1_pio_type","pio_1_pio_type","",'F');
    dataloader->AddVariable("pio_1_energy_1","pio_1_energy_1","MeV",'F');
    dataloader->AddVariable("pio_1_energy_2","pio_1_energy_2","MeV",'F');
    dataloader->AddVariable("pio_1_dis_1","pio_1_dis_1","cm",'F');
    dataloader->AddVariable("pio_1_dis_2","pio_1_dis_2","cm",'F');
    dataloader->AddVariable("pio_mip_id","pio_mip_id","",'F');

    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );

    // Apply additional cuts on the signal and background samples (can be different)
    //    TCut mycut_s = "pio_mip_id==0&&pio_filled==1&&pio_flag_pio==1"; // 831 events
    //    TCut mycut_b = "pio_mip_id==0&&pio_filled==1&&pio_flag_pio==1"; // 859 events

    TCut mycut_s = "pio_filled==1&&pio_flag_pio==1"; // 6141 events  rejected 1143
    TCut mycut_b = "pio_filled==1&&pio_flag_pio==1 && (!(truth_nue==1 && truth_CC==1))"; // 3960 events   rejected 2892
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
        "nTrain_Signal=5000:"
        "nTrain_Background=3400:"
	"nTest_Signal=1141:"
        "nTest_Background=560:"
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


void TestEvaluate(TString filename)
{
  float pio_flag;
  float pio_mip_id;
  float pio_filled;
  float pio_flag_pio;

  float pio_1_flag;
  float pio_1_mass;
  float pio_1_pio_type;
  float pio_1_energy_1;
  float pio_1_energy_2;
  float pio_1_dis_1;
  float pio_1_dis_2;

  Float_t bdt_value = 0;
    
  TFile *file = new TFile("reduced.root");
  TTree *sig = (TTree*)file->Get("sig");
  TTree *bkg = (TTree*)file->Get("bkg");


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

  sig->SetBranchAddress("pio_flag",&pio_flag);
  sig->SetBranchAddress("pio_mip_id",&pio_mip_id);
  sig->SetBranchAddress("pio_filled",&pio_filled);
  sig->SetBranchAddress("pio_flag_pio",&pio_flag_pio);
  
  sig->SetBranchAddress("pio_1_flag",&pio_1_flag);
  sig->SetBranchAddress("pio_1_mass",&pio_1_mass);
  sig->SetBranchAddress("pio_1_pio_type",&pio_1_pio_type);
  sig->SetBranchAddress("pio_1_energy_1",&pio_1_energy_1);
  sig->SetBranchAddress("pio_1_energy_2",&pio_1_energy_2);
  sig->SetBranchAddress("pio_1_dis_1",&pio_1_dis_1);
  sig->SetBranchAddress("pio_1_dis_2",&pio_1_dis_2);
  
   
  bkg->SetBranchAddress("pio_flag",&pio_flag);
  bkg->SetBranchAddress("pio_mip_id",&pio_mip_id);
  bkg->SetBranchAddress("pio_filled",&pio_filled);
  bkg->SetBranchAddress("pio_flag_pio",&pio_flag_pio);
  
  bkg->SetBranchAddress("pio_1_flag",&pio_1_flag);
  bkg->SetBranchAddress("pio_1_mass",&pio_1_mass);
  bkg->SetBranchAddress("pio_1_pio_type",&pio_1_pio_type);
  bkg->SetBranchAddress("pio_1_energy_1",&pio_1_energy_1);
  bkg->SetBranchAddress("pio_1_energy_2",&pio_1_energy_2);
  bkg->SetBranchAddress("pio_1_dis_1",&pio_1_dis_1);
  bkg->SetBranchAddress("pio_1_dis_2",&pio_1_dis_2);
  
  
  TFile *new_file = new TFile("round_1.root","RECREATE");
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

  
  Tsig->Branch("pio_flag",&pio_flag,"pio_flag/F");
  Tsig->Branch("pio_mip_id",&pio_mip_id,"pio_mip_id/F");
  Tsig->Branch("pio_filled",&pio_filled,"pio_filled/F");
  Tsig->Branch("pio_flag_pio",&pio_flag_pio,"pio_flag_pio/F");
  
  Tsig->Branch("pio_1_flag",&pio_1_flag,"pio_1_flag/F");
  Tsig->Branch("pio_1_mass",&pio_1_mass,"pio_1_mass/F");
  Tsig->Branch("pio_1_pio_type",&pio_1_pio_type,"pio_1_pio_type/F");
  Tsig->Branch("pio_1_energy_1",&pio_1_energy_1,"pio_1_energy_1/F");
  Tsig->Branch("pio_1_energy_2",&pio_1_energy_2,"pio_1_energy_2/F");
  Tsig->Branch("pio_1_dis_1",&pio_1_dis_1,"pio_1_dis_1/F");
  Tsig->Branch("pio_1_dis_2",&pio_1_dis_2,"pio_1_dis_2/F");
  Tsig->Branch("pio_1_bdt",&bdt_value,"pio_1_bdt/F");
  

  Tbkg->Branch("pio_flag",&pio_flag,"pio_flag/F");
  Tbkg->Branch("pio_mip_id",&pio_mip_id,"pio_mip_id/F");
  Tbkg->Branch("pio_filled",&pio_filled,"pio_filled/F");
  Tbkg->Branch("pio_flag_pio",&pio_flag_pio,"pio_flag_pio/F");
  
  Tbkg->Branch("pio_1_flag",&pio_1_flag,"pio_1_flag/F");
  Tbkg->Branch("pio_1_mass",&pio_1_mass,"pio_1_mass/F");
  Tbkg->Branch("pio_1_pio_type",&pio_1_pio_type,"pio_1_pio_type/F");
  Tbkg->Branch("pio_1_energy_1",&pio_1_energy_1,"pio_1_energy_1/F");
  Tbkg->Branch("pio_1_energy_2",&pio_1_energy_2,"pio_1_energy_2/F");
  Tbkg->Branch("pio_1_dis_1",&pio_1_dis_1,"pio_1_dis_1/F");
  Tbkg->Branch("pio_1_dis_2",&pio_1_dis_2,"pio_1_dis_2/F");
  Tbkg->Branch("pio_1_bdt",&bdt_value,"pio_1_bdt/F");

  
  TMVA::Reader *reader = new TMVA::Reader();
  reader->AddVariable("pio_1_mass",&pio_1_mass);
  reader->AddVariable("pio_1_pio_type",&pio_1_pio_type);
  reader->AddVariable("pio_1_energy_1",&pio_1_energy_1);
  reader->AddVariable("pio_1_energy_2",&pio_1_energy_2);
  reader->AddVariable("pio_1_dis_1",&pio_1_dis_1);
  reader->AddVariable("pio_1_dis_2",&pio_1_dis_2);
  reader->AddVariable("pio_mip_id",&pio_mip_id);
  
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

