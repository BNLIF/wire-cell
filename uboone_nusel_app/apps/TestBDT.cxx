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


void Run();
void InitInput();
void InitOutput();
void InitBDT();
void TestEvaluate();

void Run()
{
    TMVA::Tools::Instance();
    InitInput();
    InitOutput();
    InitBDT();


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

    TestEvaluate();

    // Launch the GUI for the root macros
    if (!gROOT->IsBatch()) TMVA::TMVAGui( output->GetName() );

}

void InitBDT()
{
    factory = new TMVA::Factory( "Test", output,
        "!V:!Silent:Color:DrawProgressBar:"
        "Transformations=I;D;P;G,D:"
        "AnalysisType=Classification" );


    dataloader = new TMVA::DataLoader("dataset");
    // Define the input variables that shall be used for the MVA training
    // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
    // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
    dataloader->AddVariable( "sum := var1+var2", "v1+v2", "", 'F' );
    dataloader->AddVariable( "diff := var1-var2", "v1-v2", "", 'F' );
    dataloader->AddVariable( "var3", "v3", "m", 'F' );
    dataloader->AddVariable( "var4", "v4", "MeV", 'F' );

    TTree *signalTree     = (TTree*)input->Get("TreeS");
    TTree *backgroundTree = (TTree*)input->Get("TreeB");
    dataloader->AddSignalTree(signalTree, 1.0); // can add the global event weight
    dataloader->AddBackgroundTree( backgroundTree, 1.0);
    // Set individual event weights (the variables must exist in the original TTree)
    dataloader->SetBackgroundWeightExpression( "weight" );

    // Apply additional cuts on the signal and background samples (can be different)
    TCut mycut_s = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
    TCut mycut_b = ""; // for example: TCut mycutb = "abs(var1)<0.5";
    dataloader->PrepareTrainingAndTestTree( mycut_s, mycut_b,
        "nTrain_Signal=4000:"
        "nTrain_Background=4000:"
        // "nTest_Signal=2000:"
        // "nTest_Background=2000:"
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


void TestEvaluate()
{
    Float_t sum, diff, var3, var4;
    
    TMVA::Reader *reader = new TMVA::Reader();
    reader->AddVariable( "var1+var2", &sum );
    reader->AddVariable( "var1-var2", &diff );
    reader->AddVariable( "var3", &var3 );
    reader->AddVariable( "var4", &var4 );
    reader->BookMVA( "MyBDT", "dataset/weights/Test_BDT.weights.xml");

    sum=1; diff=-1; var3=-2; var4=4;
    double bdt = reader->EvaluateMVA("MyBDT");
    cout << "Testing 1 entry: BDT score = "<< bdt << endl;

}

void InitInput()
{
    TString fname = "files/tmva_class_example.root";
    if (!gSystem->AccessPathName( fname )) {
        input = TFile::Open( fname ); // check if file in local directory exists
    }
    else {
        TFile::SetCacheFileDir(".");
        input = TFile::Open("http://root.cern.ch/files/tmva_class_example.root", "CACHEREAD");
    }
    if (!input) {
        std::cout << "ERROR: could not open data file" << std::endl;
        exit(1);
    }
    cout << "Using input file: " << input->GetName() << std::endl;
}

void InitOutput()
{
    TString outfileName( "TMVA.root" );
    output = TFile::Open( outfileName, "RECREATE" );
}


int main( int argc, char** argv )
{

    cout << "testing BDT in ROOT TMVA" << endl;
    Run();
    return 1;

}

