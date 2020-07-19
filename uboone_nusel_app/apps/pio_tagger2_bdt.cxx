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
void TestEvaluate_r1();


void Run_r2();
void InitBDT_r2();
void TestEvaluate_r2();


void convert_file();


void convert_file(){
  TFile *file = new TFile("bdtfile_0718.root");
  TTree *sig = (TTree*)file->Get("sig");
  TTree *bkg = (TTree*)file->Get("bkg");

  int run, event;

  sig->SetBranchAddress("run",&run);
  sig->SetBranchAddress("event",&event);

  bkg->SetBranchAddress("run",&run);
  bkg->SetBranchAddress("event",&event);
  
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

  
  int pio_flag;
  int pio_mip_id;
  int pio_filled;
  int pio_flag_pio;

  std::vector<double> *pio_2_v_dis2 = new std::vector<double>;
  std::vector<double> *pio_2_v_angle2 = new std::vector<double>;
  std::vector<double> *pio_2_v_acc_length = new std::vector<double>;
  std::vector<int> *pio_2_v_flag = new std::vector<int>;
  
  sig->SetBranchAddress("pio_flag",&pio_flag);
  sig->SetBranchAddress("pio_mip_id",&pio_mip_id);
  sig->SetBranchAddress("pio_filled",&pio_filled);
  sig->SetBranchAddress("pio_flag_pio",&pio_flag_pio);
  
  sig->SetBranchAddress("pio_2_v_dis2",&pio_2_v_dis2);
  sig->SetBranchAddress("pio_2_v_angle2",&pio_2_v_angle2);
  sig->SetBranchAddress("pio_2_v_acc_length",&pio_2_v_acc_length);
  sig->SetBranchAddress("pio_2_v_flag",&pio_2_v_flag);
      
  bkg->SetBranchAddress("pio_flag",&pio_flag);
  bkg->SetBranchAddress("pio_mip_id",&pio_mip_id);
  bkg->SetBranchAddress("pio_filled",&pio_filled);
  bkg->SetBranchAddress("pio_flag_pio",&pio_flag_pio);
  
  bkg->SetBranchAddress("pio_2_v_dis2",&pio_2_v_dis2);
  bkg->SetBranchAddress("pio_2_v_angle2",&pio_2_v_angle2);
  bkg->SetBranchAddress("pio_2_v_acc_length",&pio_2_v_acc_length);
  bkg->SetBranchAddress("pio_2_v_flag",&pio_2_v_flag);

  TFile *new_file = new TFile("reduced.root","RECREATE");
  TTree *Tsig = new TTree("sig","sig");
  TTree *Tbkg = new TTree("bkg","bkg");
  Tsig->SetDirectory(new_file);
  Tbkg->SetDirectory(new_file);
  
  float pio_flag_f;
  float pio_mip_id_f;
  float pio_filled_f;
  float pio_flag_pio_f;

  Tsig->Branch("run",&run,"data/I");
  Tsig->Branch("event",&event,"data/I");

  Tbkg->Branch("run",&run,"data/I");
  Tbkg->Branch("event",&event,"data/I");

  Tsig->Branch("trueEdep",&trueEdep,"data/F");
  Tsig->Branch("weight",&weight,"data/F");
  Tsig->Branch("lowEweight",&lowEweight,"data/F");
  Tsig->Branch("nueTag",&nueTag,"data/I");
  
  Tsig->Branch("pio_flag",&pio_flag_f,"pio_flag/F");
  Tsig->Branch("pio_mip_id",&pio_mip_id_f,"pio_mip_id/F");
  Tsig->Branch("pio_filled",&pio_filled_f,"pio_filled/F");
  Tsig->Branch("pio_flag_pio",&pio_flag_pio_f,"pio_flag_pio/F");
    
  
  Tbkg->Branch("trueEdep",&trueEdep,"data/F");
  Tbkg->Branch("weight",&weight,"data/F");
  Tbkg->Branch("lowEweight",&lowEweight,"data/F");
  Tbkg->Branch("nueTag",&nueTag,"data/I");

  Tbkg->Branch("pio_flag",&pio_flag_f,"pio_flag/F");
  Tbkg->Branch("pio_mip_id",&pio_mip_id_f,"pio_mip_id/F");
  Tbkg->Branch("pio_filled",&pio_filled_f,"pio_filled/F");
  Tbkg->Branch("pio_flag_pio",&pio_flag_pio_f,"pio_flag_pio/F");

  float pio_2_v_dis2_f;
  float pio_2_v_angle2_f;
  float pio_2_v_acc_length_f;
  float pio_2_v_flag_f;
  
  Tsig->Branch("pio_2_v_dis2",&pio_2_v_dis2_f,"pio_2_v_dis2/F");
  Tsig->Branch("pio_2_v_angle2",&pio_2_v_angle2_f,"pio_2_v_angle2/F");
  Tsig->Branch("pio_2_v_acc_length",&pio_2_v_acc_length_f,"pio_2_v_acc_length/F");
  Tsig->Branch("pio_2_v_flag",&pio_2_v_flag_f,"pio_2_v_flag/F");
  
  Tbkg->Branch("pio_2_v_dis2",&pio_2_v_dis2_f,"pio_2_v_dis2/F");
  Tbkg->Branch("pio_2_v_angle2",&pio_2_v_angle2_f,"pio_2_v_angle2/F");
  Tbkg->Branch("pio_2_v_acc_length",&pio_2_v_acc_length_f,"pio_2_v_acc_length/F");
  Tbkg->Branch("pio_2_v_flag",&pio_2_v_flag_f,"pio_2_v_flag/F");
  
  for (Int_t i=0;i!=sig->GetEntries();i++){
    sig->GetEntry(i);
    
    pio_flag_f = pio_flag;
    pio_mip_id_f = pio_mip_id;
    pio_filled_f = pio_filled;
    pio_flag_pio_f = pio_flag_pio;

    int ncount = 0;
    
    for (size_t j=0;j!=pio_2_v_dis2->size();j++){
      pio_2_v_dis2_f = pio_2_v_dis2->at(j);
      pio_2_v_angle2_f = pio_2_v_angle2->at(j);
      pio_2_v_acc_length_f = pio_2_v_acc_length->at(j);
      pio_2_v_flag_f = pio_2_v_flag->at(j);

      if (ncount <7 || pio_2_v_flag_f==0)
	Tsig->Fill();

      ncount ++;
    }
    
  }

   for (Int_t i=0;i!=bkg->GetEntries();i++){
    bkg->GetEntry(i);
    
    pio_flag_f = pio_flag;
    pio_mip_id_f = pio_mip_id;
    pio_filled_f = pio_filled;
    pio_flag_pio_f = pio_flag_pio;

    for (size_t j=0;j!=pio_2_v_dis2->size();j++){
      pio_2_v_dis2_f = pio_2_v_dis2->at(j);
      pio_2_v_angle2_f = pio_2_v_angle2->at(j);
      pio_2_v_acc_length_f = pio_2_v_acc_length->at(j);
      pio_2_v_flag_f = pio_2_v_flag->at(j);
      Tbkg->Fill();
    }
    
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

    TestEvaluate_r2();

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

    TestEvaluate_r1();

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

    dataloader->AddVariable("pio_2_v_dis2","pio_2_v_dis2","cm",'F');
    dataloader->AddVariable("pio_2_v_angle2","pio_2_v_angle2","",'F');
    dataloader->AddVariable("pio_2_v_acc_length","pio_2_v_acc_length","",'F');
    dataloader->AddVariable("pio_mip_id","pio_mip_id","",'F');
    
    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight " );
    dataloader->SetBackgroundWeightExpression( "weight " );

    // Apply additional cuts on the signal and background samples (can be different)
    //    TCut mycut_s = "pio_mip_id==0&&pio_filled==1&&pio_flag_pio==1"; // 831 events
    //    TCut mycut_b = "pio_mip_id==0&&pio_filled==1&&pio_flag_pio==1"; // 859 events

    TCut mycut_s = "pio_filled==1 && pio_flag_pio==0 "; // 308 / 265054, adding pi0 50
    TCut mycut_b = "pio_filled==1 && pio_flag_pio==0 && pio_2_v_flag==0 "; // 1147/471479, adding mip_id 340
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
        "nTrain_Signal=200000:"
        "nTrain_Background=1000:"
	"nTest_Signal=65054:"
        "nTest_Background=147:"
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

    dataloader->AddVariable("pio_2_v_dis2","pio_2_v_dis2","cm",'F');
    dataloader->AddVariable("pio_2_v_angle2","pio_2_v_angle2","",'F');
    dataloader->AddVariable("pio_2_v_acc_length","pio_2_v_acc_length","",'F');
    dataloader->AddVariable("pio_mip_id","pio_mip_id","",'F');
    
    TTree *signalTree     = (TTree*)input->Get("sig");
    TTree *backgroundTree = (TTree*)input->Get("bkg");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetSignalWeightExpression( "weight " );
    dataloader->SetBackgroundWeightExpression( "weight " );

    // Apply additional cuts on the signal and background samples (can be different)
    //    TCut mycut_s = "pio_mip_id==0&&pio_filled==1&&pio_flag_pio==1"; // 831 events
    //    TCut mycut_b = "pio_mip_id==0&&pio_filled==1&&pio_flag_pio==1"; // 859 events

    TCut mycut_s = "pio_filled==1 && pio_flag_pio==0 "; // 308 / 265054, adding pi0 50
    TCut mycut_b = "pio_filled==1 && pio_flag_pio==0 && (pio_2_v_flag==0 || pio_2_v_bdt<0)"; // 9828/471479, adding mip_id 340
    
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
        "nTrain_Signal=200000:"
        "nTrain_Background=7000:"
	"nTest_Signal=65054:"
        "nTest_Background=1828:"
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


void TestEvaluate_r1()
{
  float pio_flag;
  float pio_mip_id;
  float pio_filled;
  float pio_flag_pio;
    
 
  Float_t bdt_value = 0;
    
  TFile *file = new TFile("reduced.root");
  TTree *sig = (TTree*)file->Get("sig");
  TTree *bkg = (TTree*)file->Get("bkg");

  int run, event;
  sig->SetBranchAddress("run",&run);
  sig->SetBranchAddress("event",&event);

  bkg->SetBranchAddress("run",&run);
  bkg->SetBranchAddress("event",&event);
  
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
     
  bkg->SetBranchAddress("pio_flag",&pio_flag);
  bkg->SetBranchAddress("pio_mip_id",&pio_mip_id);
  bkg->SetBranchAddress("pio_filled",&pio_filled);
  bkg->SetBranchAddress("pio_flag_pio",&pio_flag_pio);

  float pio_2_v_dis2;
  float pio_2_v_angle2;
  float pio_2_v_acc_length;
  float pio_2_v_flag;
  
  sig->SetBranchAddress("pio_2_v_dis2",&pio_2_v_dis2);
  sig->SetBranchAddress("pio_2_v_angle2",&pio_2_v_angle2);
  sig->SetBranchAddress("pio_2_v_acc_length",&pio_2_v_acc_length);
  sig->SetBranchAddress("pio_2_v_flag",&pio_2_v_flag);

  bkg->SetBranchAddress("pio_2_v_dis2",&pio_2_v_dis2);
  bkg->SetBranchAddress("pio_2_v_angle2",&pio_2_v_angle2);
  bkg->SetBranchAddress("pio_2_v_acc_length",&pio_2_v_acc_length);
  bkg->SetBranchAddress("pio_2_v_flag",&pio_2_v_flag);
    
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
  
  Tsig->Branch("pio_flag",&pio_flag,"pio_flag/F");
  Tsig->Branch("pio_mip_id",&pio_mip_id,"pio_mip_id/F");
  Tsig->Branch("pio_filled",&pio_filled,"pio_filled/F");
  Tsig->Branch("pio_flag_pio",&pio_flag_pio,"pio_flag_pio/F");

  Tsig->Branch("run",&run,"data/I");
  Tsig->Branch("event",&event,"data/I");

  Tbkg->Branch("run",&run,"data/I");
  Tbkg->Branch("event",&event,"data/I");

  Tbkg->Branch("pio_flag",&pio_flag,"pio_flag/F");
  Tbkg->Branch("pio_mip_id",&pio_mip_id,"pio_mip_id/F");
  Tbkg->Branch("pio_filled",&pio_filled,"pio_filled/F");
  Tbkg->Branch("pio_flag_pio",&pio_flag_pio,"pio_flag_pio/F");
  
  Tsig->Branch("pio_2_v_dis2",&pio_2_v_dis2,"data/F");
  Tsig->Branch("pio_2_v_angle2",&pio_2_v_angle2,"data/F");
  Tsig->Branch("pio_2_v_acc_length",&pio_2_v_acc_length,"data/F");
  Tsig->Branch("pio_2_v_flag",&pio_2_v_flag,"data/F");
  Tsig->Branch("pio_2_v_bdt",&bdt_value,"data/F");
    
  Tbkg->Branch("pio_2_v_dis2",&pio_2_v_dis2,"data/F");
  Tbkg->Branch("pio_2_v_angle2",&pio_2_v_angle2,"data/F");
  Tbkg->Branch("pio_2_v_acc_length",&pio_2_v_acc_length,"data/F");
  Tbkg->Branch("pio_2_v_flag",&pio_2_v_flag,"data/F");
  Tbkg->Branch("pio_2_v_bdt",&bdt_value,"data/F");
  
  TMVA::Reader *reader = new TMVA::Reader();
  reader->AddVariable("pio_2_v_dis2",&pio_2_v_dis2);
  reader->AddVariable("pio_2_v_angle2",&pio_2_v_angle2);
  reader->AddVariable("pio_2_v_acc_length",&pio_2_v_acc_length);
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



void TestEvaluate_r2()
{
 
   float pio_flag;
  float pio_mip_id;
  float pio_filled;
  float pio_flag_pio;
    
 
  Float_t bdt_value = 0;
    
  TFile *file = new TFile("reduced.root");
  TTree *sig = (TTree*)file->Get("sig");
  TTree *bkg = (TTree*)file->Get("bkg");

  int run, event;
  sig->SetBranchAddress("run",&run);
  sig->SetBranchAddress("event",&event);

  bkg->SetBranchAddress("run",&run);
  bkg->SetBranchAddress("event",&event);
  
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
     
  bkg->SetBranchAddress("pio_flag",&pio_flag);
  bkg->SetBranchAddress("pio_mip_id",&pio_mip_id);
  bkg->SetBranchAddress("pio_filled",&pio_filled);
  bkg->SetBranchAddress("pio_flag_pio",&pio_flag_pio);

  float pio_2_v_dis2;
  float pio_2_v_angle2;
  float pio_2_v_acc_length;
  float pio_2_v_flag;
  
  sig->SetBranchAddress("pio_2_v_dis2",&pio_2_v_dis2);
  sig->SetBranchAddress("pio_2_v_angle2",&pio_2_v_angle2);
  sig->SetBranchAddress("pio_2_v_acc_length",&pio_2_v_acc_length);
  sig->SetBranchAddress("pio_2_v_flag",&pio_2_v_flag);

  bkg->SetBranchAddress("pio_2_v_dis2",&pio_2_v_dis2);
  bkg->SetBranchAddress("pio_2_v_angle2",&pio_2_v_angle2);
  bkg->SetBranchAddress("pio_2_v_acc_length",&pio_2_v_acc_length);
  bkg->SetBranchAddress("pio_2_v_flag",&pio_2_v_flag);
    
  TFile *new_file = new TFile("round_2.root","RECREATE");
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
  
  Tsig->Branch("pio_flag",&pio_flag,"pio_flag/F");
  Tsig->Branch("pio_mip_id",&pio_mip_id,"pio_mip_id/F");
  Tsig->Branch("pio_filled",&pio_filled,"pio_filled/F");
  Tsig->Branch("pio_flag_pio",&pio_flag_pio,"pio_flag_pio/F");

  Tsig->Branch("run",&run,"data/I");
  Tsig->Branch("event",&event,"data/I");

  Tbkg->Branch("run",&run,"data/I");
  Tbkg->Branch("event",&event,"data/I");

  Tbkg->Branch("pio_flag",&pio_flag,"pio_flag/F");
  Tbkg->Branch("pio_mip_id",&pio_mip_id,"pio_mip_id/F");
  Tbkg->Branch("pio_filled",&pio_filled,"pio_filled/F");
  Tbkg->Branch("pio_flag_pio",&pio_flag_pio,"pio_flag_pio/F");
  
  Tsig->Branch("pio_2_v_dis2",&pio_2_v_dis2,"data/F");
  Tsig->Branch("pio_2_v_angle2",&pio_2_v_angle2,"data/F");
  Tsig->Branch("pio_2_v_acc_length",&pio_2_v_acc_length,"data/F");
  Tsig->Branch("pio_2_v_flag",&pio_2_v_flag,"data/F");
  Tsig->Branch("pio_2_v_bdt",&bdt_value,"data/F");
    
  Tbkg->Branch("pio_2_v_dis2",&pio_2_v_dis2,"data/F");
  Tbkg->Branch("pio_2_v_angle2",&pio_2_v_angle2,"data/F");
  Tbkg->Branch("pio_2_v_acc_length",&pio_2_v_acc_length,"data/F");
  Tbkg->Branch("pio_2_v_flag",&pio_2_v_flag,"data/F");
  Tbkg->Branch("pio_2_v_bdt",&bdt_value,"data/F");
  
  TMVA::Reader *reader = new TMVA::Reader();
  reader->AddVariable("pio_2_v_dis2",&pio_2_v_dis2);
  reader->AddVariable("pio_2_v_angle2",&pio_2_v_angle2);
  reader->AddVariable("pio_2_v_acc_length",&pio_2_v_acc_length);
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

