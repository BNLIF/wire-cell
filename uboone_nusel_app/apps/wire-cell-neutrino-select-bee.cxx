#include "TH1F.h"
#include "TH2F.h"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
  if (argc<3){
    cerr << "usage: wire-cell-neutrino-select-bee match.root type [0,1]" << std::endl;
    cerr << "0 for BNB [3,5]us and 1 for extBNB [3.45,5.45]us" << std::endl;
    
    return 0;
  }
  TString filename = argv[1];
  int trigger_type = atoi(argv[2]);
  
  TFile *file = new TFile(filename);
  TTree *T_flash = (TTree*)file->Get("T_flash");
  TTree *T_match = (TTree*)file->Get("T_match");
  TTree *Trun = (TTree*)file->Get("Trun");
  Double_t time;
  Int_t type;
  Int_t flash_id;
  T_flash->SetBranchAddress("time",&time);
  T_flash->SetBranchAddress("type",&type);
  T_flash->SetBranchAddress("flash_id",&flash_id);
  Int_t runNo, subRunNo, eventNo;
  Trun->SetBranchAddress("runNo",&runNo);
  Trun->SetBranchAddress("subRunNo",&subRunNo);
  Trun->SetBranchAddress("eventNo",&eventNo);
  Trun->GetEntry(0);
  Int_t tpc_cluster_id;
  T_match->SetBranchAddress("tpc_cluster_id",&tpc_cluster_id);
  T_match->SetBranchAddress("flash_id",&flash_id);
  

  Int_t saved_flash_id = -1;
  Double_t saved_time = 1e9;
  
  for (int i=0;i!=T_flash->GetEntries();i++){
    T_flash->GetEntry(i);
    if (type==2){
      if (trigger_type==0 && time >= 3 && time <=5){
	saved_flash_id = flash_id;
	saved_time = time;
      }else if (trigger_type==1 && time >=3.45 && time <= 5.45){
	saved_flash_id = flash_id;
	saved_time = time;
      }
    }
  }

  std::set<int> matched_flash_ids;
  for (int i=0;i!=T_match->GetEntries();i++){
    T_match->GetEntry(i);
    if (flash_id >= 0)
      matched_flash_ids.insert(flash_id);
  }
  
  if (matched_flash_ids.find(saved_flash_id)!=matched_flash_ids.end() ){
    std::cout << runNo << " " << subRunNo << " " << eventNo << " " << saved_time << std::endl;
    //std::cout << "Found flash " << saved_flash_id << " at " << saved_time << " us " << std::endl;
    return 1;
  }else{
    return 0;
    //    std::cout << "In-time beam flash not found" << std::endl;
  }

 
}
