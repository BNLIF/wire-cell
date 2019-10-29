#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TFile.h"
#include "TString.h"
#include "TGraph2D.h"
#include "TColor.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TTree.h"
#include <iostream>

#include <set>

using namespace std;

int main(int argc, char* argv[])
{
  if (argc < 2){
    cerr << "usage: wire-cell-uboone /path/to/match.root" << endl;
    return 1;
  }

  TString filename = argv[1];
  TFile *file = new TFile(filename);
  TTree *Trun = (TTree*)file->Get("Trun");
  
  Int_t runNo;
  Int_t eventNo;
  Int_t subRunNo;
  unsigned int triggerbits;
  Trun->SetBranchAddress("triggerBits",&triggerbits);
  Trun->SetBranchAddress("runNo",&runNo);
  Trun->SetBranchAddress("subRunNo",&subRunNo);
  Trun->SetBranchAddress("eventNo",&eventNo);
  Trun->GetEntry(0);

  double lowerwindow = 0;
  double upperwindow = 0;
  if((triggerbits>>11) & 1U) { lowerwindow = 3.0; upperwindow = 5.0;}//bnb
  if((triggerbits>>9) & 1U) { lowerwindow = 3.45; upperwindow = 5.45;} //extbnb

  
  // TTree *T_match = (TTree*)file->Get("T_match");
  // std::cout << T_match->GetEntries() << std::endl;

  TTree *T_flash = (TTree*)file->Get("T_flash");
  //std::cout << T_flash->GetEntries() << std::endl;
  Int_t flash_id;
  Double_t time;
  Int_t type;
  Double_t total_PE;
  T_flash->SetBranchAddress("flash_id",&flash_id);
  T_flash->SetBranchAddress("time",&time);
  T_flash->SetBranchAddress("type",&type);
  T_flash->SetBranchAddress("total_PE",&total_PE);
  
  bool flag_save = false;
  std::set<int> good_flash_ids;
  std::map<int,double> pe_map;
  std::map<int,double> time_map;
  for (Int_t i=0;i!=T_flash->GetEntries();i++){
    T_flash->GetEntry(i);
    if (time >=lowerwindow && time <= upperwindow){
      // std::cout << flash_id << " " << time << std::endl;
      good_flash_ids.insert(flash_id);
      pe_map[flash_id] = total_PE;
      time_map[flash_id] = time;
    }
  }
  if (good_flash_ids.size()>0) flag_save = true;

  if (flag_save){
    TTree *T_match = (TTree*)file->Get("T_match");
    Int_t temp_flash_id;
    Int_t tpc_cluster_id;
    Int_t event_type;

    Char_t flag_close_to_PMT;
    Char_t flag_at_x_boundary;

    bool flag_close_to_PMT1=true;
    bool flag_at_x_boundary1=true;
    
    Double_t ks_dis;
    Double_t chi2;
    Int_t ndf;
    Double_t pe[32];
    Double_t pe_pred[32];
    Double_t cluster_length;
    
    
    T_match->SetBranchAddress("flash_id",&temp_flash_id);
    T_match->SetBranchAddress("tpc_cluster_id",&tpc_cluster_id);
    T_match->SetBranchAddress("event_type",&event_type);
    
    T_match->SetBranchAddress("flag_close_to_PMT",&flag_close_to_PMT);
    T_match->SetBranchAddress("flag_at_x_boundary",&flag_at_x_boundary);

    T_match->SetBranchAddress("ks_dis",&ks_dis);
    T_match->SetBranchAddress("chi2",&chi2);
    T_match->SetBranchAddress("ndf",&ndf);
    T_match->SetBranchAddress("pe_pred",pe_pred);
    T_match->SetBranchAddress("pe_meas",pe);
    T_match->SetBranchAddress("cluster_length",&cluster_length);
    
    flag_save = false;
    double total_pred_pe;
    double max_flash_pe;
    std::set<int> matched_flash_ids;
    std::map<int,int> entry_map;
    
    for (Int_t i=0;i!=T_match->GetEntries();i++){
      T_match->GetEntry(i);

      if (good_flash_ids.find(temp_flash_id)!=good_flash_ids.end()){
	matched_flash_ids.insert(temp_flash_id);
	entry_map[temp_flash_id] = i;
      }
      
    }

    if (matched_flash_ids.size()>0){
      T_match->GetEntry(entry_map[(*matched_flash_ids.begin())]);
      flag_close_to_PMT1 = flag_close_to_PMT;
      flag_at_x_boundary1 = flag_at_x_boundary;
      //std::cout << flag_close_to_PMT1 << " " << flag_at_x_boundary1 << std::endl;
      if (good_flash_ids.find(temp_flash_id)!=good_flash_ids.end()){
	total_pred_pe = 0;
	max_flash_pe = 0;
	for (int j=0;j!=32;j++){
	  total_pred_pe += pe_pred[j];
	  if (pe[j] > max_flash_pe) max_flash_pe = pe[j];
	}
	
	
	flag_save = true;
      }
    }
    
    
    if (flag_save)
      std::cout << runNo << " " << subRunNo << " " << eventNo << " " <<
	temp_flash_id << " " << tpc_cluster_id << " " << time_map[temp_flash_id] << " " << pe_map[temp_flash_id] << " " << total_pred_pe << " " <<
	event_type << " " << flag_close_to_PMT1 << " " << flag_at_x_boundary1 << " " <<
	ks_dis << " " << chi2 << " " << ndf << " " << max_flash_pe << " " << cluster_length << std::endl;
  }
  
  //   
    
  

}
