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
#include "TMath.h"

#include "tagger.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"


using namespace std;
using namespace LEEana;

int main( int argc, char** argv )
{
  int process = 1; // odd  (default training is odd ...)
  // 0 for even ...

  for (Int_t i=1;i!=argc;i++){
    switch(argv[i][1]){
    case 'p':
      process = atoi(&argv[i][2]);
      break;
    }
  }

  TString filename = argv[1];
  TString name = argv[2];
  
  std::set< std::pair<int, int> > set_run_subrun;
  TFile *file = new TFile(filename);
  TTree *T_bdt = (TTree*)file->Get("bdt");
  int run, subrun;
  T_bdt->SetBranchStatus("*",0);
  T_bdt->SetBranchStatus("run",1); T_bdt->SetBranchAddress("run", &run);
  T_bdt->SetBranchStatus("subrun",1); T_bdt->SetBranchAddress("subrun", &subrun);
  for (Int_t i = 0;i!=T_bdt->GetEntries();i++){
    T_bdt->GetEntry(i);
    if (subrun %2 == 1 && process==1 || subrun%2 == 0 && process != 1) continue;
    set_run_subrun.insert(std::make_pair(run,subrun));
  }

  for (auto it = set_run_subrun.begin(); it != set_run_subrun.end(); it++){
    std::cout << name << " " << it->first << " " << it->second << std::endl;
  }
  
  
  
  
}
