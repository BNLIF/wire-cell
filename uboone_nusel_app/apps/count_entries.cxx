#include "TFile.h"
#include "TTree.h"

#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
  if (argc < 2){
    cerr << "usage: wire-cell-uboone /path/to/match.root" << endl;
    return 1;
  }
  const char* root_file = argv[1];
  TFile *file1 = new TFile(root_file);
  TTree *T = (TTree*)file1->Get("Trun");
  Int_t run_no, subrun_no, event_no;
  T->SetBranchAddress("runNo",&run_no);
  T->SetBranchAddress("subRunNo",&subrun_no);
  T->SetBranchAddress("eventNo",&event_no);
  for (int i=0;i!=T->GetEntries();i++){
    T->GetEntry(i);
    std::cout << run_no << " " << subrun_no << " " << event_no << " " << root_file << " " << i << std::endl; 
  }
  
  //  std::cout << T->GetEntries() << std::endl;
}
