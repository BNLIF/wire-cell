#include "TFile.h"
#include "TTree.h"

#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
  if (argc < 2){
    cerr << "usage: test_neutrino /path/to/nue.root" << endl;
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
  if((triggerbits>>11) & 1U) { lowerwindow=3.1875; upperwindow=4.96876;} // bnb 
  if((triggerbits>>9) & 1U) { lowerwindow=3.5625; upperwindow=5.34376; } //extbnb 

  TTree *T_match = (TTree*)file1->Get("T_match");
  Int_t tpc_cluster_id;
  Int_t event_type;
  Int_t flash_id;
  Double_t flash_time;
  Double_t cluster_length;
  T_match->SetBranchAddress("tpc_cluster_id",&tpc_cluster_id);
  T_match->SetBranchAddress("flash_id",&flash_id);
  T_match->SetBranchAddress("neutrino_type",&event_type);
  T_match->SetBranchAddress("flash_time",&flash_time);
  for (int i =0;i!=T_match->GetEntries();i++){
    T_match->GetEntry(i);
    
    if (flash_time < lowerwindow || flash_time > upperwindow) continue;
    int flag_cosmic = (event_type >> 1) & 1U;
    int flag_numu = (event_type >> 2) & 1U;
    int flag_nc = (event_type >> 3) & 1U;
    int flag_long_muon = (event_type >> 4) & 1U;
    int flag_nue = (event_type >> 5) & 1U;

    
    
    std::cout << runNo << "_" << subRunNo << "_" << eventNo << " " << flash_id << " " << tpc_cluster_id << " " << flash_time << " " << event_type << " " << flag_cosmic << " " << flag_numu << " " << flag_nc << " " << flag_long_muon << " " << flag_nue <<  std::endl;
   
  }
}
