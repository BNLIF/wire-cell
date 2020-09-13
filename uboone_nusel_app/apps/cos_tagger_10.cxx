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
  
  TFile *file0 = new TFile("bdt.root"); // sig x1
  
  TTree *sig = (TTree*)file0->Get("sig");
  TTree *bkg = (TTree*)file0->Get("bkg");


  float cosmict_flag_1;

  float event_no = 0;
  
  std::vector<float> *cosmict_flag_10 = new std::vector<float>;  // front upstream (dirt)
  
  std::vector<float> *cosmict_10_flag_inside= new std::vector<float> ;
  std::vector<float> *cosmict_10_vtx_z= new std::vector<float>;
  std::vector<float> *cosmict_10_flag_shower= new std::vector<float>;
  std::vector<float> *cosmict_10_flag_dir_weak= new std::vector<float>;
  std::vector<float> *cosmict_10_angle_beam= new std::vector<float>;
  std::vector<float> *cosmict_10_length= new std::vector<float>;
  
  sig->SetBranchAddress("cosmict_flag_1",&cosmict_flag_1);
  sig->SetBranchAddress("cosmict_flag_10",&cosmict_flag_10);
  sig->SetBranchAddress("cosmict_10_flag_inside",&cosmict_10_flag_inside);
  sig->SetBranchAddress("cosmict_10_vtx_z",&cosmict_10_vtx_z);
  sig->SetBranchAddress("cosmict_10_flag_shower",&cosmict_10_flag_shower);
  sig->SetBranchAddress("cosmict_10_flag_dir_weak",&cosmict_10_flag_dir_weak);
  sig->SetBranchAddress("cosmict_10_angle_beam",&cosmict_10_angle_beam);
  sig->SetBranchAddress("cosmict_10_length",&cosmict_10_length);

  bkg->SetBranchAddress("cosmict_flag_1",&cosmict_flag_1);
  bkg->SetBranchAddress("cosmict_flag_10",&cosmict_flag_10);
  bkg->SetBranchAddress("cosmict_10_flag_inside",&cosmict_10_flag_inside);
  bkg->SetBranchAddress("cosmict_10_vtx_z",&cosmict_10_vtx_z);
  bkg->SetBranchAddress("cosmict_10_flag_shower",&cosmict_10_flag_shower);
  bkg->SetBranchAddress("cosmict_10_flag_dir_weak",&cosmict_10_flag_dir_weak);
  bkg->SetBranchAddress("cosmict_10_angle_beam",&cosmict_10_angle_beam);
  bkg->SetBranchAddress("cosmict_10_length",&cosmict_10_length);

  float weight;
  sig->SetBranchAddress("weight",&weight);
  bkg->SetBranchAddress("weight",&weight);
  
  TFile *new_file = new TFile("round_0.root","RECREATE");
  TTree *Tsig = new TTree("sig","sig");
  TTree *Tbkg = new TTree("bkg","bkg");
  Tsig->SetDirectory(new_file);
  Tbkg->SetDirectory(new_file);

  Tsig->Branch("event",&event_no,"data/F");
  Tbkg->Branch("event",&event_no,"data/F");

  Tsig->Branch("cosmict_flag_1",&cosmict_flag_1,"data/F");
  Tbkg->Branch("cosmict_flag_1",&cosmict_flag_1,"data/F");
  
  float cosmict_flag_10_f;
  
  float cosmict_10_flag_inside_f;
  float cosmict_10_vtx_z_f;
  float cosmict_10_flag_shower_f;
  float cosmict_10_flag_dir_weak_f;
  float cosmict_10_angle_beam_f;
  float cosmict_10_length_f;

  Tsig->Branch("cosmict_flag_10",&cosmict_flag_10_f,"data/F");
  Tsig->Branch("cosmict_10_flag_inside",&cosmict_10_flag_inside_f,"data/F");
  Tsig->Branch("cosmict_10_vtx_z",&cosmict_10_vtx_z_f,"data/F");
  Tsig->Branch("cosmict_10_flag_shower",&cosmict_10_flag_shower_f,"data/F");
  Tsig->Branch("cosmict_10_flag_dir_weak",&cosmict_10_flag_dir_weak_f,"data/F");
  Tsig->Branch("cosmict_10_angle_beam",&cosmict_10_angle_beam_f,"data/F");
  Tsig->Branch("cosmict_10_length",&cosmict_10_length_f,"data/F");
  
  Tbkg->Branch("cosmict_flag_10",&cosmict_flag_10_f,"data/F");
  Tbkg->Branch("cosmict_10_flag_inside",&cosmict_10_flag_inside_f,"data/F");
  Tbkg->Branch("cosmict_10_vtx_z",&cosmict_10_vtx_z_f,"data/F");
  Tbkg->Branch("cosmict_10_flag_shower",&cosmict_10_flag_shower_f,"data/F");
  Tbkg->Branch("cosmict_10_flag_dir_weak",&cosmict_10_flag_dir_weak_f,"data/F");
  Tbkg->Branch("cosmict_10_angle_beam",&cosmict_10_angle_beam_f,"data/F");
  Tbkg->Branch("cosmict_10_length",&cosmict_10_length_f,"data/F");

 
  Tsig->Branch("weight",&weight,"data/F");
  Tbkg->Branch("weight",&weight,"data/F");
  
  for (Int_t i=0;i!=sig->GetEntries();i++){
    sig->GetEntry(i);
    
    for (size_t j=0;j!=cosmict_flag_10->size();j++){

      cosmict_flag_10_f = cosmict_flag_10->at(j);
      cosmict_10_flag_inside_f = cosmict_10_flag_inside->at(j);
      cosmict_10_vtx_z_f = cosmict_10_vtx_z->at(j);
      cosmict_10_flag_shower_f = cosmict_10_flag_shower->at(j);
      cosmict_10_flag_dir_weak_f = cosmict_10_flag_dir_weak->at(j);
      cosmict_10_angle_beam_f = cosmict_10_angle_beam->at(j);
      cosmict_10_length_f = cosmict_10_length->at(j);

      if (std::isnan(cosmict_10_angle_beam_f)) cosmict_10_angle_beam_f = 0;
      
      Tsig->Fill();
    }
    
    event_no ++;
    
  }

  for (Int_t i=0;i!=bkg->GetEntries();i++){
    bkg->GetEntry(i);
    
    for (size_t j=0;j!=cosmict_flag_10->size();j++){

      cosmict_flag_10_f = cosmict_flag_10->at(j);
      cosmict_10_flag_inside_f = cosmict_10_flag_inside->at(j);
      cosmict_10_vtx_z_f = cosmict_10_vtx_z->at(j);
      cosmict_10_flag_shower_f = cosmict_10_flag_shower->at(j);
      cosmict_10_flag_dir_weak_f = cosmict_10_flag_dir_weak->at(j);
      cosmict_10_angle_beam_f = cosmict_10_angle_beam->at(j);
      cosmict_10_length_f = cosmict_10_length->at(j);

      if (std::isnan(cosmict_10_angle_beam_f)) cosmict_10_angle_beam_f = 0;
      
      Tbkg->Fill();
    }
    
    event_no ++;
    
  }

  
  
  cout << "signal tree entries: " << Tsig->GetEntries() << " / " << Tsig->GetEntries("cosmict_flag_1==0")<< endl;
  cout << "background tree entries: " << Tbkg->GetEntries() << " / " << Tbkg->GetEntries(" (cosmict_flag_10==1)&& cosmict_flag_1==0")<< endl;
  
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

    TString fname = "./round_0.root";
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
    
   

    //    dataloader->AddVariable("cosmict_10_flag_inside","cosmict_10_flag_inside","",'F');
    dataloader->AddVariable("cosmict_10_vtx_z","cosmict_10_vtx_z","cm",'F');
    dataloader->AddVariable("cosmict_10_flag_shower","cosmict_10_flag_shower","",'F');
    dataloader->AddVariable("cosmict_10_flag_dir_weak","cosmict_10_flag_dir_weak","",'F');
    dataloader->AddVariable("cosmict_10_angle_beam","cosmict_10_angle_beam","deg",'F');
    dataloader->AddVariable("cosmict_10_length","cosmict_10_length","cm",'F');
        
    
    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight " );
    dataloader->SetBackgroundWeightExpression( "weight " );

    // Apply additional cuts on the signal and background samples (can be different)
    TCut mycut_s = "cosmict_flag_1==0"; //  995
    TCut mycut_b = "cosmict_flag_1==0 && cosmict_flag_10==1"; //  76
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
        "nTrain_Signal=800:"
        "nTrain_Background=60:"
	"nTest_Signal=195:"
        "nTest_Background=16:"
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
        "NTrees=400:"
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
    
    //    dataloader->AddVariable("cosmict_10_flag_inside","cosmict_10_flag_inside","",'F');
    dataloader->AddVariable("cosmict_10_vtx_z","cosmict_10_vtx_z","cm",'F');
    dataloader->AddVariable("cosmict_10_flag_shower","cosmict_10_flag_shower","",'F');
    dataloader->AddVariable("cosmict_10_flag_dir_weak","cosmict_10_flag_dir_weak","",'F');
    dataloader->AddVariable("cosmict_10_angle_beam","cosmict_10_angle_beam","deg",'F');
    dataloader->AddVariable("cosmict_10_length","cosmict_10_length","cm",'F');
        
    
   
   
    
    
    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight " );
    dataloader->SetBackgroundWeightExpression( "weight " );

    // Apply additional cuts on the signal and background samples (can be different)
    TCut mycut_s = "cosmict_flag_1==0"; // 
    TCut mycut_b = "cosmict_flag_1==0 && (cosmict_flag_10==1 || cos_10_bdt < 0.2)"; // 
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
        "nTrain_Signal=800:"
        "nTrain_Background=90:"
	"nTest_Signal=195:"
        "nTest_Background=22:"
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
        "NTrees=400:"
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

  float event;
  float weight;
  float cosmict_flag_1;
  float cosmict_flag_10;
  
  float cosmict_10_flag_inside;
  float cosmict_10_vtx_z;
  float cosmict_10_flag_shower;
  float cosmict_10_flag_dir_weak;
  float cosmict_10_angle_beam;
  float cosmict_10_length;

  TFile *file = new TFile("round_0.root");
  TTree *sig = (TTree*)file->Get("sig");
  TTree *bkg = (TTree*)file->Get("bkg");

  sig->SetBranchAddress("event",&event);
  bkg->SetBranchAddress("event",&event);

  sig->SetBranchAddress("weight",&weight);
  bkg->SetBranchAddress("weight",&weight);
  
  sig->SetBranchAddress("cosmict_flag_1",&cosmict_flag_1);
  sig->SetBranchAddress("cosmict_flag_10",&cosmict_flag_10);
  sig->SetBranchAddress("cosmict_10_flag_inside",&cosmict_10_flag_inside);
  sig->SetBranchAddress("cosmict_10_vtx_z",&cosmict_10_vtx_z);
  sig->SetBranchAddress("cosmict_10_flag_shower",&cosmict_10_flag_shower);
  sig->SetBranchAddress("cosmict_10_flag_dir_weak",&cosmict_10_flag_dir_weak);
  sig->SetBranchAddress("cosmict_10_angle_beam",&cosmict_10_angle_beam);
  sig->SetBranchAddress("cosmict_10_length",&cosmict_10_length);

  bkg->SetBranchAddress("cosmict_flag_1",&cosmict_flag_1);
  bkg->SetBranchAddress("cosmict_flag_10",&cosmict_flag_10);
  bkg->SetBranchAddress("cosmict_10_flag_inside",&cosmict_10_flag_inside);
  bkg->SetBranchAddress("cosmict_10_vtx_z",&cosmict_10_vtx_z);
  bkg->SetBranchAddress("cosmict_10_flag_shower",&cosmict_10_flag_shower);
  bkg->SetBranchAddress("cosmict_10_flag_dir_weak",&cosmict_10_flag_dir_weak);
  bkg->SetBranchAddress("cosmict_10_angle_beam",&cosmict_10_angle_beam);
  bkg->SetBranchAddress("cosmict_10_length",&cosmict_10_length);

  
 
  TFile *new_file = new TFile(filename,"RECREATE");
  TTree *Tsig = new TTree("sig","sig");
  TTree *Tbkg = new TTree("bkg","bkg");
  Tsig->SetDirectory(new_file);
  Tbkg->SetDirectory(new_file);

 
  Tsig->Branch("event",&event,"data/F");
  Tbkg->Branch("event",&event,"data/F");

  Tsig->Branch("cos_10_bdt",&bdt_value,"data/F");
  Tbkg->Branch("cos_10_bdt",&bdt_value,"data/F");

  Tsig->Branch("cosmict_flag_1",&cosmict_flag_1,"data/F");
  Tbkg->Branch("cosmict_flag_1",&cosmict_flag_1,"data/F");
  
  Tbkg->Branch("cosmict_flag_10",&cosmict_flag_10,"data/F");
  Tbkg->Branch("cosmict_10_flag_inside",&cosmict_10_flag_inside,"data/F");
  Tbkg->Branch("cosmict_10_vtx_z",&cosmict_10_vtx_z,"data/F");
  Tbkg->Branch("cosmict_10_flag_shower",&cosmict_10_flag_shower,"data/F");
  Tbkg->Branch("cosmict_10_flag_dir_weak",&cosmict_10_flag_dir_weak,"data/F");
  Tbkg->Branch("cosmict_10_angle_beam",&cosmict_10_angle_beam,"data/F");
  Tbkg->Branch("cosmict_10_length",&cosmict_10_length,"data/F");
  Tbkg->Branch("weight",&weight,"data/F");

  Tsig->Branch("cosmict_flag_10",&cosmict_flag_10,"data/F");
  Tsig->Branch("cosmict_10_flag_inside",&cosmict_10_flag_inside,"data/F");
  Tsig->Branch("cosmict_10_vtx_z",&cosmict_10_vtx_z,"data/F");
  Tsig->Branch("cosmict_10_flag_shower",&cosmict_10_flag_shower,"data/F");
  Tsig->Branch("cosmict_10_flag_dir_weak",&cosmict_10_flag_dir_weak,"data/F");
  Tsig->Branch("cosmict_10_angle_beam",&cosmict_10_angle_beam,"data/F");
  Tsig->Branch("cosmict_10_length",&cosmict_10_length,"data/F");
  Tsig->Branch("weight",&weight,"data/F");
  
  TMVA::Reader *reader = new TMVA::Reader();
  reader->AddVariable("cosmict_10_vtx_z",&cosmict_10_vtx_z);
  reader->AddVariable("cosmict_10_flag_shower",&cosmict_10_flag_shower);
  reader->AddVariable("cosmict_10_flag_dir_weak",&cosmict_10_flag_dir_weak);
  reader->AddVariable("cosmict_10_angle_beam",&cosmict_10_angle_beam);
  reader->AddVariable("cosmict_10_length",&cosmict_10_length);
      
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

