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

  bool image_fail=false;
  Float_t lm_cluster_length = 999999;
  TTree* T_eval = (TTree*)file1->Get("T_eval");
  if(T_eval->GetBranch("image_fail")){
    T_eval->SetBranchAddress("image_fail",&image_fail);
    T_eval->SetBranchAddress("lm_cluster_length",&lm_cluster_length);
  }

  unsigned int triggerbits;
  T->SetBranchAddress("triggerBits",&triggerbits);
  int time_offset;
  T->SetBranchAddress("time_offset",&time_offset);
  T->GetEntry(0);
  double lowerwindow = 0;
  double upperwindow = 0;
  
  //enlarge window ... 
  if((triggerbits>>11) & 1U) { lowerwindow=3.1875; upperwindow=4.96876;} // bnb 
  if ((triggerbits>>12) & 1U) { lowerwindow=4.9295; upperwindow=16.6483;} // NUMI
  if(((triggerbits>>9) & 1U)&& time_offset != 5) { lowerwindow=3.5625; upperwindow=5.34376; } //extbnb
  if (((triggerbits>>9) & 1U) && time_offset == 5) {lowerwindow=5.3045; upperwindow=17.0233;} // EXTNUMI

  TTree *T_match = (TTree*)file1->Get("T_match");
  Double_t flash_time;
  T_match->SetBranchAddress("flash_time",&flash_time);
  for (int i =0;i!=T_match->GetEntries();i++){
    T_match->GetEntry(i);
    T_eval->GetEntry(i); 
    if (image_fail==true) continue;
    if (flash_time <= lowerwindow || flash_time >= upperwindow) continue;

    std::cout << runNo << "_" << subRunNo << "_" << eventNo << " " << lm_cluster_length <<std::endl;


    
  }
}
