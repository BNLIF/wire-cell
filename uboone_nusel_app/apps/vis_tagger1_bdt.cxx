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

  int vis_1_filled;
  int vis_1_n_vtx_segs;
  double vis_1_energy;
  int vis_1_num_good_tracks;
  double vis_1_max_angle;
  double vis_1_max_shower_angle;
  double vis_1_tmp_length1;
  double vis_1_tmp_length2;
  double vis_1_particle_type;
  int vis_1_flag;
  
  sig->SetBranchAddress("vis_1_filled",&vis_1_filled);
  sig->SetBranchAddress("vis_1_n_vtx_segs",&vis_1_n_vtx_segs);
  sig->SetBranchAddress("vis_1_energy",&vis_1_energy);
  sig->SetBranchAddress("vis_1_num_good_tracks",&vis_1_num_good_tracks);
  sig->SetBranchAddress("vis_1_max_angle",&vis_1_max_angle);
  sig->SetBranchAddress("vis_1_max_shower_angle",&vis_1_max_shower_angle);
  sig->SetBranchAddress("vis_1_tmp_length1",&vis_1_tmp_length1);
  sig->SetBranchAddress("vis_1_tmp_length2",&vis_1_tmp_length2);
  sig->SetBranchAddress("vis_1_particle_type",&vis_1_particle_type);
  sig->SetBranchAddress("vis_1_flag",&vis_1_flag);

  bkg->SetBranchAddress("vis_1_filled",&vis_1_filled);
  bkg->SetBranchAddress("vis_1_n_vtx_segs",&vis_1_n_vtx_segs);
  bkg->SetBranchAddress("vis_1_energy",&vis_1_energy);
  bkg->SetBranchAddress("vis_1_num_good_tracks",&vis_1_num_good_tracks);
  bkg->SetBranchAddress("vis_1_max_angle",&vis_1_max_angle);
  bkg->SetBranchAddress("vis_1_max_shower_angle",&vis_1_max_shower_angle);
  bkg->SetBranchAddress("vis_1_tmp_length1",&vis_1_tmp_length1);
  bkg->SetBranchAddress("vis_1_tmp_length2",&vis_1_tmp_length2);
  bkg->SetBranchAddress("vis_1_particle_type",&vis_1_particle_type);
  bkg->SetBranchAddress("vis_1_flag",&vis_1_flag);

  
  TFile *new_file = new TFile("reduced.root","RECREATE");
  TTree *Tsig = new TTree("sig","sig");
  TTree *Tbkg = new TTree("bkg","bkg");
  Tsig->SetDirectory(new_file);
  Tbkg->SetDirectory(new_file);

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


Tsig->Branch("truth_inFV",&truth_inFV,"data/I");
  Tsig->Branch("truth_CC",&truth_CC,"data/I");
  Tsig->Branch("truth_nue",&truth_nue,"data/I");
  Tsig->Branch("truth_cosmic",&truth_cosmic,"data/I");

  Tbkg->Branch("truth_inFV",&truth_inFV,"data/I");
  Tbkg->Branch("truth_CC",&truth_CC,"data/I");
  Tbkg->Branch("truth_nue",&truth_nue,"data/I");
  Tbkg->Branch("truth_cosmic",&truth_cosmic,"data/I");
  
  Tsig->Branch("trueEdep",&trueEdep,"data/F");
  Tsig->Branch("weight",&weight,"data/F");
  Tsig->Branch("lowEweight",&lowEweight,"data/F");
  Tsig->Branch("nueTag",&nueTag,"data/I");
    
  Tbkg->Branch("trueEdep",&trueEdep,"data/F");
  Tbkg->Branch("weight",&weight,"data/F");
  Tbkg->Branch("lowEweight",&lowEweight,"data/F");
  Tbkg->Branch("nueTag",&nueTag,"data/I");
  

 
  float vis_1_filled_f;
  float vis_1_n_vtx_segs_f;
  float vis_1_energy_f;
  float vis_1_num_good_tracks_f;
  float vis_1_max_angle_f;
  float vis_1_max_shower_angle_f;
  float vis_1_tmp_length1_f;
  float vis_1_tmp_length2_f;
  float vis_1_particle_type_f;
  float vis_1_flag_f;

  Tsig->Branch("vis_1_filled",&vis_1_filled_f,"vis_1_filled/F");
  Tsig->Branch("vis_1_n_vtx_segs",&vis_1_n_vtx_segs_f,"vis_1_n_vtx_segs/F");
  Tsig->Branch("vis_1_energy",&vis_1_energy_f,"vis_1_energy/F");
  Tsig->Branch("vis_1_num_good_tracks",&vis_1_num_good_tracks_f,"vis_1_num_good_tracks/F");
  Tsig->Branch("vis_1_max_angle",&vis_1_max_angle_f,"vis_1_max_angle/F");
  Tsig->Branch("vis_1_max_shower_angle",&vis_1_max_shower_angle_f,"vis_1_max_shower_angle/F");
  Tsig->Branch("vis_1_tmp_length1",&vis_1_tmp_length1_f,"vis_1_tmp_length1/F");
  Tsig->Branch("vis_1_tmp_length2",&vis_1_tmp_length2_f,"vis_1_tmp_length2/F");
  Tsig->Branch("vis_1_particle_type",&vis_1_particle_type_f,"vis_1_particle_type/F");
  Tsig->Branch("vis_1_flag",&vis_1_flag_f,"vis_1_flag/F");

  Tbkg->Branch("vis_1_filled",&vis_1_filled_f,"vis_1_filled/F");
  Tbkg->Branch("vis_1_n_vtx_segs",&vis_1_n_vtx_segs_f,"vis_1_n_vtx_segs/F");
  Tbkg->Branch("vis_1_energy",&vis_1_energy_f,"vis_1_energy/F");
  Tbkg->Branch("vis_1_num_good_tracks",&vis_1_num_good_tracks_f,"vis_1_num_good_tracks/F");
  Tbkg->Branch("vis_1_max_angle",&vis_1_max_angle_f,"vis_1_max_angle/F");
  Tbkg->Branch("vis_1_max_shower_angle",&vis_1_max_shower_angle_f,"vis_1_max_shower_angle/F");
  Tbkg->Branch("vis_1_tmp_length1",&vis_1_tmp_length1_f,"vis_1_tmp_length1/F");
  Tbkg->Branch("vis_1_tmp_length2",&vis_1_tmp_length2_f,"vis_1_tmp_length2/F");
  Tbkg->Branch("vis_1_particle_type",&vis_1_particle_type_f,"vis_1_particle_type/F");
  Tbkg->Branch("vis_1_flag",&vis_1_flag_f,"vis_1_flag/F");
  
  
  for (Int_t i=0;i!=sig->GetEntries();i++){
    sig->GetEntry(i);

   vis_1_filled_f = vis_1_filled;
   vis_1_n_vtx_segs_f = vis_1_n_vtx_segs;
   vis_1_energy_f = vis_1_energy;
   vis_1_num_good_tracks_f = vis_1_num_good_tracks;
   vis_1_max_angle_f = vis_1_max_angle;
   vis_1_max_shower_angle_f = vis_1_max_shower_angle;
   vis_1_tmp_length1_f = vis_1_tmp_length1;
   vis_1_tmp_length2_f = vis_1_tmp_length2;
   vis_1_particle_type_f = vis_1_particle_type;
   vis_1_flag_f = vis_1_flag;
    
    
    Tsig->Fill();
  }

   for (Int_t i=0;i!=bkg->GetEntries();i++){
    bkg->GetEntry(i);
    
    vis_1_filled_f = vis_1_filled;
    vis_1_n_vtx_segs_f = vis_1_n_vtx_segs;
    vis_1_energy_f = vis_1_energy;
    vis_1_num_good_tracks_f = vis_1_num_good_tracks;
    vis_1_max_angle_f = vis_1_max_angle;
    vis_1_max_shower_angle_f = vis_1_max_shower_angle;
    vis_1_tmp_length1_f = vis_1_tmp_length1;
    vis_1_tmp_length2_f = vis_1_tmp_length2;
    vis_1_particle_type_f = vis_1_particle_type;
    vis_1_flag_f = vis_1_flag;
    
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
    
    
    dataloader->AddVariable("vis_1_n_vtx_segs","vis_1_n_vtx_segs","",'F');
    dataloader->AddVariable("vis_1_energy","vis_1_energy","MeV",'F');
    dataloader->AddVariable("vis_1_num_good_tracks","vis_1_num_good_tracks","",'F');
    dataloader->AddVariable("vis_1_max_angle","vis_1_max_angle","deg",'F');
    dataloader->AddVariable("vis_1_max_shower_angle","vis_1_max_shower_angle","deg",'F');
    dataloader->AddVariable("vis_1_tmp_length1","vis_1_tmp_length1","cm",'F');
    dataloader->AddVariable("vis_1_tmp_length2","vis_1_tmp_length2","cm",'F');
    //    dataloader->AddVariable("vis_1_particle_type","vis_1_particle_type","",'F');
    

    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "vis_1_filled==1"; // 106/37627
    TCut mycut_b = "vis_1_filled==1 && (vis_1_flag==0)&& (!(truth_nue==1 && truth_CC==1))"; // 911/16540
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=20000:"
					    "nTrain_Background=800:"
					    "nTest_Signal=10000:"
					    "nTest_Background=104:"
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
    
    
    dataloader->AddVariable("vis_1_n_vtx_segs","vis_1_n_vtx_segs","",'F');
    dataloader->AddVariable("vis_1_energy","vis_1_energy","MeV",'F');
    dataloader->AddVariable("vis_1_num_good_tracks","vis_1_num_good_tracks","",'F');
    dataloader->AddVariable("vis_1_max_angle","vis_1_max_angle","deg",'F');
    dataloader->AddVariable("vis_1_max_shower_angle","vis_1_max_shower_angle","deg",'F');
    dataloader->AddVariable("vis_1_tmp_length1","vis_1_tmp_length1","cm",'F');
    dataloader->AddVariable("vis_1_tmp_length2","vis_1_tmp_length2","cm",'F');
    //    dataloader->AddVariable("vis_1_particle_type","vis_1_particle_type","",'F');
    

    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight * lowEweight" );
    dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
      
    TCut mycut_s = "vis_1_filled==1"; // 106/37627
    TCut mycut_b = "vis_1_filled==1 && (vis_1_flag==0 || vis_1_bdt < 0)&& (!(truth_nue==1 && truth_CC==1))"; // 1364/16540
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
					    "nTrain_Signal=20000:"
					    "nTrain_Background=1100:"
					    "nTest_Signal=10000:"
					    "nTest_Background=218:"
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


  float vis_1_filled;
  float vis_1_n_vtx_segs;
  float vis_1_energy;
  float vis_1_num_good_tracks;
  float vis_1_max_angle;
  float vis_1_max_shower_angle;
  float vis_1_tmp_length1;
  float vis_1_tmp_length2;
  float vis_1_particle_type;
  float vis_1_flag;
  
  
  sig->SetBranchAddress("vis_1_filled",&vis_1_filled);
  sig->SetBranchAddress("vis_1_n_vtx_segs",&vis_1_n_vtx_segs);
  sig->SetBranchAddress("vis_1_energy",&vis_1_energy);
  sig->SetBranchAddress("vis_1_num_good_tracks",&vis_1_num_good_tracks);
  sig->SetBranchAddress("vis_1_max_angle",&vis_1_max_angle);
  sig->SetBranchAddress("vis_1_max_shower_angle",&vis_1_max_shower_angle);
  sig->SetBranchAddress("vis_1_tmp_length1",&vis_1_tmp_length1);
  sig->SetBranchAddress("vis_1_tmp_length2",&vis_1_tmp_length2);
  sig->SetBranchAddress("vis_1_particle_type",&vis_1_particle_type);
  sig->SetBranchAddress("vis_1_flag",&vis_1_flag);

  bkg->SetBranchAddress("vis_1_filled",&vis_1_filled);
  bkg->SetBranchAddress("vis_1_n_vtx_segs",&vis_1_n_vtx_segs);
  bkg->SetBranchAddress("vis_1_energy",&vis_1_energy);
  bkg->SetBranchAddress("vis_1_num_good_tracks",&vis_1_num_good_tracks);
  bkg->SetBranchAddress("vis_1_max_angle",&vis_1_max_angle);
  bkg->SetBranchAddress("vis_1_max_shower_angle",&vis_1_max_shower_angle);
  bkg->SetBranchAddress("vis_1_tmp_length1",&vis_1_tmp_length1);
  bkg->SetBranchAddress("vis_1_tmp_length2",&vis_1_tmp_length2);
  bkg->SetBranchAddress("vis_1_particle_type",&vis_1_particle_type);
  bkg->SetBranchAddress("vis_1_flag",&vis_1_flag);

  
  
  
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
  
  Tsig->Branch("vis_1_filled",&vis_1_filled,"vis_1_filled/F");
  Tsig->Branch("vis_1_n_vtx_segs",&vis_1_n_vtx_segs,"vis_1_n_vtx_segs/F");
  Tsig->Branch("vis_1_energy",&vis_1_energy,"vis_1_energy/F");
  Tsig->Branch("vis_1_num_good_tracks",&vis_1_num_good_tracks,"vis_1_num_good_tracks/F");
  Tsig->Branch("vis_1_max_angle",&vis_1_max_angle,"vis_1_max_angle/F");
  Tsig->Branch("vis_1_max_shower_angle",&vis_1_max_shower_angle,"vis_1_max_shower_angle/F");
  Tsig->Branch("vis_1_tmp_length1",&vis_1_tmp_length1,"vis_1_tmp_length1/F");
  Tsig->Branch("vis_1_tmp_length2",&vis_1_tmp_length2,"vis_1_tmp_length2/F");
  Tsig->Branch("vis_1_particle_type",&vis_1_particle_type,"vis_1_particle_type/F");
  Tsig->Branch("vis_1_flag",&vis_1_flag,"vis_1_flag/F");

  Tbkg->Branch("vis_1_filled",&vis_1_filled,"vis_1_filled/F");
  Tbkg->Branch("vis_1_n_vtx_segs",&vis_1_n_vtx_segs,"vis_1_n_vtx_segs/F");
  Tbkg->Branch("vis_1_energy",&vis_1_energy,"vis_1_energy/F");
  Tbkg->Branch("vis_1_num_good_tracks",&vis_1_num_good_tracks,"vis_1_num_good_tracks/F");
  Tbkg->Branch("vis_1_max_angle",&vis_1_max_angle,"vis_1_max_angle/F");
  Tbkg->Branch("vis_1_max_shower_angle",&vis_1_max_shower_angle,"vis_1_max_shower_angle/F");
  Tbkg->Branch("vis_1_tmp_length1",&vis_1_tmp_length1,"vis_1_tmp_length1/F");
  Tbkg->Branch("vis_1_tmp_length2",&vis_1_tmp_length2,"vis_1_tmp_length2/F");
  Tbkg->Branch("vis_1_particle_type",&vis_1_particle_type,"vis_1_particle_type/F");
  Tbkg->Branch("vis_1_flag",&vis_1_flag,"vis_1_flag/F");



  
  Tsig->Branch("vis_1_bdt", &bdt_value,"data/F");
  Tbkg->Branch("vis_1_bdt", &bdt_value,"data/F");

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


Tsig->Branch("truth_inFV",&truth_inFV,"data/I");
  Tsig->Branch("truth_CC",&truth_CC,"data/I");
  Tsig->Branch("truth_nue",&truth_nue,"data/I");
  Tsig->Branch("truth_cosmic",&truth_cosmic,"data/I");

  Tbkg->Branch("truth_inFV",&truth_inFV,"data/I");
  Tbkg->Branch("truth_CC",&truth_CC,"data/I");
  Tbkg->Branch("truth_nue",&truth_nue,"data/I");
  Tbkg->Branch("truth_cosmic",&truth_cosmic,"data/I");
  
  
  TMVA::Reader *reader = new TMVA::Reader();
  
  reader->AddVariable("vis_1_n_vtx_segs",&vis_1_n_vtx_segs);
  reader->AddVariable("vis_1_energy",&vis_1_energy);
  reader->AddVariable("vis_1_num_good_tracks",&vis_1_num_good_tracks);
  reader->AddVariable("vis_1_max_angle",&vis_1_max_angle);
  reader->AddVariable("vis_1_max_shower_angle",&vis_1_max_shower_angle);
  reader->AddVariable("vis_1_tmp_length1",&vis_1_tmp_length1);
  reader->AddVariable("vis_1_tmp_length2",&vis_1_tmp_length2);
  //  reader->AddVariable("vis_1_particle_type",&vis_1_particle_type);
    
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

