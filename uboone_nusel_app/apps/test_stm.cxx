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
  Int_t runNo, subRunNo, eventNo;
  T->SetBranchAddress("runNo",&runNo);
  T->SetBranchAddress("subRunNo",&subRunNo);
  T->SetBranchAddress("eventNo",&eventNo);
  unsigned int triggerbits;
  T->SetBranchAddress("triggerBits",&triggerbits);
  T->GetEntry(0);
  double lowerwindow = 0;
  double upperwindow = 0;
  
  //enlarge window ... 
  if(triggerbits==2048) { lowerwindow = 3.0; upperwindow = 5.0; }// bnb  
  if(triggerbits==512) { lowerwindow = 3.45; upperwindow = 5.45; } // extbnb

  TTree *T_match = (TTree*)file1->Get("T_match");
  Int_t tpc_cluster_id;
  Int_t event_type;
  Int_t flash_id;
  Double_t flash_time;
  T_match->SetBranchAddress("tpc_cluster_id",&tpc_cluster_id);
  T_match->SetBranchAddress("flash_id",&flash_id);
  T_match->SetBranchAddress("event_type",&event_type);
  T_match->SetBranchAddress("flash_time",&flash_time);

  for (int i =0;i!=T_match->GetEntries();i++){
      T_match->GetEntry(i);
      int flag_tgm = (event_type >> 3) & 1U;
      int flag_low_energy = (event_type >> 4) & 1U;
      int flag_lm = (event_type >> 1) & 1U;
      int flag_fully_contained = (event_type >> 2) & 1U;
      int flag_stm = (event_type >> 5) & 1U;
      int flag_full_detector_dead = (event_type >> 6) & 1U;
      
      std::cout << runNo << "_" << subRunNo << "_" << eventNo << " " << flash_id << " " << tpc_cluster_id << " " << flash_time << " " << event_type << " " << flag_low_energy << " " << flag_lm << " " << flag_tgm << " " << flag_fully_contained << " " << flag_stm << " " << flag_full_detector_dead << std::endl;
  }
}
