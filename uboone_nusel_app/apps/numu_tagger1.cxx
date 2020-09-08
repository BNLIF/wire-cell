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
  
  std::vector<float> *numu_cc_flag_1 = new std::vector<float>;
  std::vector<float> *numu_cc_1_particle_type= new std::vector<float>;
  std::vector<float> *numu_cc_1_length= new std::vector<float>;
  std::vector<float> *numu_cc_1_medium_dQ_dx= new std::vector<float>;
  std::vector<float> *numu_cc_1_dQ_dx_cut= new std::vector<float>;
  std::vector<float> *numu_cc_1_direct_length= new std::vector<float>;
  std::vector<float> *numu_cc_1_n_daughter_tracks= new std::vector<float>;
  std::vector<float> *numu_cc_1_n_daughter_all= new std::vector<float>;
  
  sig->SetBranchAddress("cosmict_flag_1",&cosmict_flag_1);
  sig->SetBranchAddress("numu_cc_flag_1",&numu_cc_flag_1);
  sig->SetBranchAddress("numu_cc_1_particle_type",&numu_cc_1_particle_type);
  sig->SetBranchAddress("numu_cc_1_length",&numu_cc_1_length);
  sig->SetBranchAddress("numu_cc_1_medium_dQ_dx",&numu_cc_1_medium_dQ_dx);
  sig->SetBranchAddress("numu_cc_1_dQ_dx_cut",&numu_cc_1_dQ_dx_cut);
  sig->SetBranchAddress("numu_cc_1_direct_length",&numu_cc_1_direct_length);
  sig->SetBranchAddress("numu_cc_1_n_daughter_tracks",&numu_cc_1_n_daughter_tracks);
  sig->SetBranchAddress("numu_cc_1_n_daughter_all",&numu_cc_1_n_daughter_all);

  
  bkg->SetBranchAddress("cosmict_flag_1",&cosmict_flag_1);
  bkg->SetBranchAddress("numu_cc_flag_1",&numu_cc_flag_1);
  bkg->SetBranchAddress("numu_cc_1_particle_type",&numu_cc_1_particle_type);
  bkg->SetBranchAddress("numu_cc_1_length",&numu_cc_1_length);
  bkg->SetBranchAddress("numu_cc_1_medium_dQ_dx",&numu_cc_1_medium_dQ_dx);
  bkg->SetBranchAddress("numu_cc_1_dQ_dx_cut",&numu_cc_1_dQ_dx_cut);
  bkg->SetBranchAddress("numu_cc_1_direct_length",&numu_cc_1_direct_length);
  bkg->SetBranchAddress("numu_cc_1_n_daughter_tracks",&numu_cc_1_n_daughter_tracks);
  bkg->SetBranchAddress("numu_cc_1_n_daughter_all",&numu_cc_1_n_daughter_all);

  
 
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
  
  float numu_cc_flag_1_f;
  float numu_cc_1_particle_type_f;
  float numu_cc_1_length_f;
  float numu_cc_1_medium_dQ_dx_f;
  float numu_cc_1_dQ_dx_cut_f;
  float numu_cc_1_direct_length_f;
  float numu_cc_1_n_daughter_tracks_f;
  float numu_cc_1_n_daughter_all_f;
  
  Tsig->Branch("numu_cc_flag_1",&numu_cc_flag_1_f,"data/F");
  Tsig->Branch("numu_cc_1_particle_type",&numu_cc_1_particle_type_f,"data/F");
  Tsig->Branch("numu_cc_1_length",&numu_cc_1_length_f,"data/F");
  Tsig->Branch("numu_cc_1_medium_dQ_dx",&numu_cc_1_medium_dQ_dx_f,"data/F");
  Tsig->Branch("numu_cc_1_dQ_dx_cut",&numu_cc_1_dQ_dx_cut_f,"data/F");
  Tsig->Branch("numu_cc_1_direct_length",&numu_cc_1_direct_length_f,"data/F");
  Tsig->Branch("numu_cc_1_n_daughter_tracks",&numu_cc_1_n_daughter_tracks_f,"data/F");
  Tsig->Branch("numu_cc_1_n_daughter_all",&numu_cc_1_n_daughter_all_f,"data/F");

  Tbkg->Branch("numu_cc_flag_1",&numu_cc_flag_1_f,"data/F");
  Tbkg->Branch("numu_cc_1_particle_type",&numu_cc_1_particle_type_f,"data/F");
  Tbkg->Branch("numu_cc_1_length",&numu_cc_1_length_f,"data/F");
  Tbkg->Branch("numu_cc_1_medium_dQ_dx",&numu_cc_1_medium_dQ_dx_f,"data/F");
  Tbkg->Branch("numu_cc_1_dQ_dx_cut",&numu_cc_1_dQ_dx_cut_f,"data/F");
  Tbkg->Branch("numu_cc_1_direct_length",&numu_cc_1_direct_length_f,"data/F");
  Tbkg->Branch("numu_cc_1_n_daughter_tracks",&numu_cc_1_n_daughter_tracks_f,"data/F");
  Tbkg->Branch("numu_cc_1_n_daughter_all",&numu_cc_1_n_daughter_all_f,"data/F");

  

  Tsig->Branch("weight",&weight,"data/F");
  Tbkg->Branch("weight",&weight,"data/F");
  
  for (Int_t i=0;i!=sig->GetEntries();i++){
    sig->GetEntry(i);
    
    for (size_t j=0;j!=numu_cc_flag_1->size();j++){
      
      numu_cc_flag_1_f = numu_cc_flag_1->at(j);
      numu_cc_1_particle_type_f = numu_cc_1_particle_type->at(j);
      numu_cc_1_length_f = numu_cc_1_length->at(j);
      numu_cc_1_medium_dQ_dx_f = numu_cc_1_medium_dQ_dx->at(j);
      numu_cc_1_dQ_dx_cut_f = numu_cc_1_dQ_dx_cut->at(j);
      numu_cc_1_direct_length_f = numu_cc_1_direct_length->at(j);
      numu_cc_1_n_daughter_tracks_f = numu_cc_1_n_daughter_tracks->at(j);
      numu_cc_1_n_daughter_all_f = numu_cc_1_n_daughter_all->at(j);
      
      if (std::isinf(numu_cc_1_dQ_dx_cut_f))  numu_cc_1_dQ_dx_cut_f = 10;
      
      Tsig->Fill();
    }
    
    event_no ++;
    
  }

  for (Int_t i=0;i!=bkg->GetEntries();i++){
    bkg->GetEntry(i);
    

    for (size_t j=0;j!=numu_cc_flag_1->size();j++){
      
      numu_cc_flag_1_f = numu_cc_flag_1->at(j);
      numu_cc_1_particle_type_f = numu_cc_1_particle_type->at(j);
      numu_cc_1_length_f = numu_cc_1_length->at(j);
      numu_cc_1_medium_dQ_dx_f = numu_cc_1_medium_dQ_dx->at(j);
      numu_cc_1_dQ_dx_cut_f = numu_cc_1_dQ_dx_cut->at(j);
      numu_cc_1_direct_length_f = numu_cc_1_direct_length->at(j);
      numu_cc_1_n_daughter_tracks_f = numu_cc_1_n_daughter_tracks->at(j);
      numu_cc_1_n_daughter_all_f = numu_cc_1_n_daughter_all->at(j);

      if (std::isinf(numu_cc_1_dQ_dx_cut_f))  numu_cc_1_dQ_dx_cut_f = 10;
      
      Tbkg->Fill();
    }
    
    event_no ++;
    
  }


  
  cout << "signal tree entries: " << Tsig->GetEntries() << " / " << Tsig->GetEntries("(numu_cc_flag_1 == 1) && cosmict_flag_1==0")<< endl;
 cout << "background tree entries: " << Tbkg->GetEntries() << " / " << Tbkg->GetEntries(" cosmict_flag_1==0")<< endl;
  
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
    
   
  
  dataloader->AddVariable("numu_cc_1_particle_type","numu_cc_1_particle_type","",'F');
  dataloader->AddVariable("numu_cc_1_length","numu_cc_1_length","cm",'F');
  dataloader->AddVariable("numu_cc_1_medium_dQ_dx","numu_cc_1_medium_dQ_dx","MeV/cm",'F');
  dataloader->AddVariable("numu_cc_1_dQ_dx_cut","numu_cc_1_dQ_dx_cut","MeV/cm",'F');
  dataloader->AddVariable("numu_cc_1_direct_length","numu_cc_1_direct_length","cm",'F');
  dataloader->AddVariable("numu_cc_1_n_daughter_tracks","numu_cc_1_n_daughter_tracks","",'F');
  dataloader->AddVariable("numu_cc_1_n_daughter_all","numu_cc_1_n_daughter_all","",'F');
  
  
  TTree *signalTree     = (TTree*)input->Get("sig");
  TTree *backgroundTree = (TTree*)input->Get("bkg");
  dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
  dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
  dataloader->SetSignalWeightExpression( "weight " );
  dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
    TCut mycut_s = "cosmict_flag_1==0 && numu_cc_flag_1 == 1"; // 
    TCut mycut_b = "cosmict_flag_1==0"; // 
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
        "nTrain_Signal=20000:"
        "nTrain_Background=6000:"
	"nTest_Signal=2947:"
        "nTest_Background=1085:"
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

    
  
  dataloader->AddVariable("numu_cc_1_particle_type","numu_cc_1_particle_type","",'F');
  dataloader->AddVariable("numu_cc_1_length","numu_cc_1_length","cm",'F');
  dataloader->AddVariable("numu_cc_1_medium_dQ_dx","numu_cc_1_medium_dQ_dx","MeV/cm",'F');
  dataloader->AddVariable("numu_cc_1_dQ_dx_cut","numu_cc_1_dQ_dx_cut","MeV/cm",'F');
  dataloader->AddVariable("numu_cc_1_direct_length","numu_cc_1_direct_length","cm",'F');
  dataloader->AddVariable("numu_cc_1_n_daughter_tracks","numu_cc_1_n_daughter_tracks","",'F');
  dataloader->AddVariable("numu_cc_1_n_daughter_all","numu_cc_1_n_daughter_all","",'F');
  
  
  TTree *signalTree     = (TTree*)input->Get("sig");
  TTree *backgroundTree = (TTree*)input->Get("bkg");
  dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
  dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
  dataloader->SetSignalWeightExpression( "weight " );
  dataloader->SetBackgroundWeightExpression( "weight " );
    
    // Apply additional cuts on the signal and background samples (can be different)
    TCut mycut_s = "cosmict_flag_1==0 && (numu_cc_flag_1 == 1 || numu_cc_1_bdt > -0.05)"; // 
    TCut mycut_b = "cosmict_flag_1==0"; // 
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
        "nTrain_Signal=21000:"
        "nTrain_Background=6000:"
	"nTest_Signal=3310:"
        "nTest_Background=1085:"
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

  float event;
  float weight;
  float cosmict_flag_1;

  float numu_cc_flag_1;
  float numu_cc_1_particle_type;
  float numu_cc_1_length;
  float numu_cc_1_medium_dQ_dx;
  float numu_cc_1_dQ_dx_cut;
  float numu_cc_1_direct_length;
  float numu_cc_1_n_daughter_tracks;
  float numu_cc_1_n_daughter_all;
  
  TFile *file = new TFile("round_0.root");
  TTree *sig = (TTree*)file->Get("sig");
  TTree *bkg = (TTree*)file->Get("bkg");

  sig->SetBranchAddress("event",&event);
  bkg->SetBranchAddress("event",&event);

  sig->SetBranchAddress("weight",&weight);
  bkg->SetBranchAddress("weight",&weight);
  
  sig->SetBranchAddress("cosmict_flag_1",&cosmict_flag_1);
  bkg->SetBranchAddress("cosmict_flag_1",&cosmict_flag_1);

  
  sig->SetBranchAddress("numu_cc_flag_1",&numu_cc_flag_1);
  sig->SetBranchAddress("numu_cc_1_particle_type",&numu_cc_1_particle_type);
  sig->SetBranchAddress("numu_cc_1_length",&numu_cc_1_length);
  sig->SetBranchAddress("numu_cc_1_medium_dQ_dx",&numu_cc_1_medium_dQ_dx);
  sig->SetBranchAddress("numu_cc_1_dQ_dx_cut",&numu_cc_1_dQ_dx_cut);
  sig->SetBranchAddress("numu_cc_1_direct_length",&numu_cc_1_direct_length);
  sig->SetBranchAddress("numu_cc_1_n_daughter_tracks",&numu_cc_1_n_daughter_tracks);
  sig->SetBranchAddress("numu_cc_1_n_daughter_all",&numu_cc_1_n_daughter_all);

  
  bkg->SetBranchAddress("numu_cc_flag_1",&numu_cc_flag_1);
  bkg->SetBranchAddress("numu_cc_1_particle_type",&numu_cc_1_particle_type);
  bkg->SetBranchAddress("numu_cc_1_length",&numu_cc_1_length);
  bkg->SetBranchAddress("numu_cc_1_medium_dQ_dx",&numu_cc_1_medium_dQ_dx);
  bkg->SetBranchAddress("numu_cc_1_dQ_dx_cut",&numu_cc_1_dQ_dx_cut);
  bkg->SetBranchAddress("numu_cc_1_direct_length",&numu_cc_1_direct_length);
  bkg->SetBranchAddress("numu_cc_1_n_daughter_tracks",&numu_cc_1_n_daughter_tracks);
  bkg->SetBranchAddress("numu_cc_1_n_daughter_all",&numu_cc_1_n_daughter_all);
  
 
  TFile *new_file = new TFile(filename,"RECREATE");
  TTree *Tsig = new TTree("sig","sig");
  TTree *Tbkg = new TTree("bkg","bkg");
  Tsig->SetDirectory(new_file);
  Tbkg->SetDirectory(new_file);

 
  Tsig->Branch("event",&event,"data/F");
  Tbkg->Branch("event",&event,"data/F");
  
  Tsig->Branch("weight",&weight,"data/F");
  Tbkg->Branch("weight",&weight,"data/F");
  
  Tsig->Branch("numu_cc_1_bdt",&bdt_value,"data/F");
  Tbkg->Branch("numu_cc_1_bdt",&bdt_value,"data/F");

  Tsig->Branch("cosmict_flag_1",&cosmict_flag_1,"data/F");
  Tbkg->Branch("cosmict_flag_1",&cosmict_flag_1,"data/F");


  Tsig->Branch("numu_cc_flag_1",&numu_cc_flag_1,"data/F");
  Tsig->Branch("numu_cc_1_particle_type",&numu_cc_1_particle_type,"data/F");
  Tsig->Branch("numu_cc_1_length",&numu_cc_1_length,"data/F");
  Tsig->Branch("numu_cc_1_medium_dQ_dx",&numu_cc_1_medium_dQ_dx,"data/F");
  Tsig->Branch("numu_cc_1_dQ_dx_cut",&numu_cc_1_dQ_dx_cut,"data/F");
  Tsig->Branch("numu_cc_1_direct_length",&numu_cc_1_direct_length,"data/F");
  Tsig->Branch("numu_cc_1_n_daughter_tracks",&numu_cc_1_n_daughter_tracks,"data/F");
  Tsig->Branch("numu_cc_1_n_daughter_all",&numu_cc_1_n_daughter_all,"data/F");

  Tbkg->Branch("numu_cc_flag_1",&numu_cc_flag_1,"data/F");
  Tbkg->Branch("numu_cc_1_particle_type",&numu_cc_1_particle_type,"data/F");
  Tbkg->Branch("numu_cc_1_length",&numu_cc_1_length,"data/F");
  Tbkg->Branch("numu_cc_1_medium_dQ_dx",&numu_cc_1_medium_dQ_dx,"data/F");
  Tbkg->Branch("numu_cc_1_dQ_dx_cut",&numu_cc_1_dQ_dx_cut,"data/F");
  Tbkg->Branch("numu_cc_1_direct_length",&numu_cc_1_direct_length,"data/F");
  Tbkg->Branch("numu_cc_1_n_daughter_tracks",&numu_cc_1_n_daughter_tracks,"data/F");
  Tbkg->Branch("numu_cc_1_n_daughter_all",&numu_cc_1_n_daughter_all,"data/F");

  
 
  
  TMVA::Reader *reader = new TMVA::Reader();
  reader->AddVariable("numu_cc_1_particle_type",&numu_cc_1_particle_type);
  reader->AddVariable("numu_cc_1_length",&numu_cc_1_length);
  reader->AddVariable("numu_cc_1_medium_dQ_dx",&numu_cc_1_medium_dQ_dx);
  reader->AddVariable("numu_cc_1_dQ_dx_cut",&numu_cc_1_dQ_dx_cut);
  reader->AddVariable("numu_cc_1_direct_length",&numu_cc_1_direct_length);
  reader->AddVariable("numu_cc_1_n_daughter_tracks",&numu_cc_1_n_daughter_tracks);
  reader->AddVariable("numu_cc_1_n_daughter_all",&numu_cc_1_n_daughter_all);
      
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

