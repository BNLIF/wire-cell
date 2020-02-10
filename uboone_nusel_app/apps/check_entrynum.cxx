#include "TFile.h"
#include "TTree.h"

#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
  if (argc < 2){
    cerr << "usage: wire-cell-uboone /path/to/input.root run_subrun_event" << endl;
    return 1;
  }
  const char* root_file = argv[1];
  const char* event_check = argv[2];
  TFile *file1 = new TFile(root_file);
  TTree *T = (TTree*)file1->Get("Trun");
  Int_t run_no, subrun_no, event_no;
  T->SetBranchStatus("*", 0);
  T->SetBranchStatus("runNo", 1);
  T->SetBranchStatus("subRunNo", 1);
  T->SetBranchStatus("eventNo", 1);
  T->SetBranchAddress("runNo",&run_no);
  T->SetBranchAddress("subRunNo",&subrun_no);
  T->SetBranchAddress("eventNo",&event_no);
  for (int i=0;i!=T->GetEntries();i++){
    T->GetEntry(i);
    std::string event_runinfo = std::to_string(run_no)+"_"+std::to_string(subrun_no)+"_"+std::to_string(event_no);
    if( event_runinfo != std::string(event_check) ) continue;
    std::cout<<i<<endl;
    //std::cout << run_no << " " << subrun_no << " " << event_no << " " << root_file << " " << i << std::endl; 
  }
  //  std::cout << T->GetEntries() << std::endl;
}
