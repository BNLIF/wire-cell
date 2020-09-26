#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"

using namespace std;

int main( int argc, char** argv )
{
  if (argc < 2) {
    std::cout << "./check_passing_rate #filename -d[data 1/MC 0]" << std::endl;
    return -1;
  }

  int flag_data = 0; // MC (with truth information)
  int flag_verbose = 0;
  for (Int_t i=1;i!=argc;i++){
    switch(argv[i][1]){
    case 'd':
      flag_data = atoi(&argv[i][2]); // whether this is data or not ...
      break;
    case 'v':
      flag_verbose = atoi(&argv[i][2]); // whether this is data or not ...
      break;
    }
  }
  
  TString filename = argv[1];

  
  TFile *file = new TFile(filename);
  TTree *T_pot = (TTree*)file->Get("wcpselection/T_pot"); // POT tree ...
  bool flag_pot = 0;
  if (T_pot->GetEntries() != 0 ) flag_pot = 1;

  Int_t runNo, subRunNo;
  Double_t pot_tor875;
  T_pot->SetBranchAddress("runNo",&runNo);
  T_pot->SetBranchAddress("subRunNo",&subRunNo);
  T_pot->SetBranchAddress("pot_tor875",&pot_tor875);
  
  double total_pot = 0;
  std::map<std::pair<int,int>, double> map_rs_pot;
  std::map<std::pair<int, int>, int> map_rs_n;
  std::map<std::pair<int, int>, std::set<int> > map_rs_f1p5; // Reco 1.5
  std::map<std::pair<int, int>, std::set<int> > map_rs_f2stm; // Reco2 stm
  std::map<std::pair<int, int>, std::set<int> > map_rs_f2pr; // Reco2 Pattern recognition
  
  
  for (Int_t i=0;i!=T_pot->GetEntries();i++){
    T_pot->GetEntry(i);
    total_pot += pot_tor875;
    map_rs_pot[std::make_pair(runNo, subRunNo)] = pot_tor875;
    map_rs_n[std::make_pair(runNo, subRunNo)] = 0; // set 0
  }

  std::cout << filename << std::endl;
  std::cout << "Total POT: " << total_pot << std::endl;

  // check failure mode
  TTree *T_eval = (TTree*)file->Get("wcpselection/T_eval");
  T_eval->SetBranchStatus("*",0);
  Int_t stm_eventtype, stm_lowenergy, stm_LM, stm_TGM, stm_STM, stm_FullDead;
  Float_t stm_cluster_length;
  
  T_eval->SetBranchStatus("stm_eventtype",1); T_eval->SetBranchAddress("stm_eventtype",&stm_eventtype);
  T_eval->SetBranchStatus("stm_lowenergy",1); T_eval->SetBranchAddress("stm_lowenergy",&stm_lowenergy);
  T_eval->SetBranchStatus("stm_LM",1); T_eval->SetBranchAddress("stm_LM",&stm_LM);
  T_eval->SetBranchStatus("stm_TGM",1); T_eval->SetBranchAddress("stm_TGM",&stm_TGM);
  T_eval->SetBranchStatus("stm_STM",1); T_eval->SetBranchAddress("stm_STM",&stm_STM);
  T_eval->SetBranchStatus("stm_FullDead",1); T_eval->SetBranchAddress("stm_FullDead",&stm_FullDead);
  T_eval->SetBranchStatus("stm_clusterlength",1); T_eval->SetBranchAddress("stm_clusterlength",&stm_cluster_length);

  bool match_found;
  T_eval->SetBranchStatus("match_found",1); T_eval->SetBranchAddress("match_found",&match_found);
  Int_t run, subrun, event;
  T_eval->SetBranchStatus("run",1); T_eval->SetBranchAddress("run",&run);
  T_eval->SetBranchStatus("subrun",1); T_eval->SetBranchAddress("subrun",&subrun);
  T_eval->SetBranchStatus("event",1); T_eval->SetBranchAddress("event",&event);

  TTree *T_PFeval = (TTree*)file->Get("wcpselection/T_PFeval");
  T_PFeval->SetBranchStatus("*",0);
  Float_t reco_nuvtxX;
  T_PFeval->SetBranchAddress("reco_nuvtxX", &reco_nuvtxX);

  
  bool flag_presel = false;
  bool flag_generic = false;
  for (Int_t i=0;i!=T_eval->GetEntries();i++){
    T_eval->GetEntry(i);
    T_PFeval->GetEntry(i);
    
    flag_presel = false;
    flag_generic = false;
    if (match_found != 0 && stm_eventtype != 0 && stm_lowenergy ==0 && stm_LM ==0 && stm_TGM ==0 && stm_STM==0 && stm_FullDead == 0 && stm_cluster_length >0) {
      flag_presel = true; // preselection ...
      if (stm_cluster_length > 15) flag_generic = true; // generic
      
    }

    // std::cout << match_found << " " << stm_lowenergy << std::endl;
    
    // add events ...
    map_rs_n[std::make_pair(run, subrun)] ++;
    // Reco 1.5 failure ...
    if (match_found == -1) map_rs_f1p5[std::make_pair(run, subrun)].insert(event);
    if (match_found == 1 && stm_lowenergy == -1) map_rs_f2stm[std::make_pair(run, subrun)].insert(event);
    if (flag_presel && reco_nuvtxX == -1) map_rs_f2pr[std::make_pair(run, subrun)].insert(event);
  }

  // Reco 1.5
  {
    double tmp_pot = 0;
    for (auto it = map_rs_f1p5.begin(); it != map_rs_f1p5.end(); it++){
      tmp_pot += map_rs_pot[it->first];
    }
    std::cout << "Reco 1.5 failure in subrun: " << map_rs_f1p5.size() << "/" << map_rs_n.size() << " ratio: " << map_rs_f1p5.size()*1.0 /map_rs_n.size() << "   POT: " << tmp_pot << " ratio: " << tmp_pot/total_pot << std::endl;

    tmp_pot = 0;
    int tmp_event = 0;
    for (auto it = map_rs_f1p5.begin(); it != map_rs_f1p5.end(); it++){
      if (flag_verbose) std::cout << "1p5: " << it->first.first << " " << it->first.second << " " << map_rs_f1p5[it->first].size() << "/" << map_rs_n[it->first] << " " << std::endl; 

      tmp_event += map_rs_f1p5[it->first].size();
      tmp_pot += map_rs_f1p5[it->first].size() * 1.0 / map_rs_n[it->first] * map_rs_pot[it->first];
    }
    std::cout << "Reco 1.5 failure in event: " << tmp_event << "/" << T_eval->GetEntries() << "  POT: " << tmp_pot << " ratio " << tmp_pot/total_pot << std::endl;
    
  }

  // Reco 2 stm
  {
    double tmp_pot = 0;
    for (auto it = map_rs_f2stm.begin(); it != map_rs_f2stm.end(); it++){
      //  bool a = map_rs_pot.find(it->first) == map_rs_pot.end();
      //std::cout << it->first.first << " " << it->first.second << " " << a << std::endl;
      tmp_pot += map_rs_pot[it->first];
    }
    std::cout << "Reco 2.0 stm failure in subrun: " << map_rs_f2stm.size() << "/" << map_rs_n.size() << " POT: " << tmp_pot << " ratio " << tmp_pot/total_pot << std::endl;

    tmp_pot = 0;
    int tmp_event = 0;
    for (auto it = map_rs_f2stm.begin(); it != map_rs_f2stm.end(); it++){
      if (flag_verbose) std::cout << "stm: " << it->first.first << " " << it->first.second << " " << map_rs_f2stm[it->first].size() << "/" << map_rs_n[it->first] << " " << std::endl; 

      tmp_event += map_rs_f2stm[it->first].size();
      tmp_pot += map_rs_f2stm[it->first].size() * 1.0 / map_rs_n[it->first] * map_rs_pot[it->first];
    }
    std::cout << "Reco 2.0 stm failure in event: " << tmp_event << "/" << T_eval->GetEntries() << "  POT: " << tmp_pot << " ratio " << tmp_pot/total_pot << std::endl;   
  }

  // Reco 2 pr
  {
    double tmp_pot = 0;
    for (auto it = map_rs_f2pr.begin(); it != map_rs_f2pr.end(); it++){
      tmp_pot += map_rs_pot[it->first];
    }
    std::cout << "Reco 2.0 PR failure in subrun: " << map_rs_f2pr.size() << "/" << map_rs_n.size() << " POT: " << tmp_pot << " ratio " << tmp_pot/total_pot << std::endl;

    tmp_pot = 0;
    int tmp_event = 0;
    for (auto it = map_rs_f2pr.begin(); it != map_rs_f2pr.end(); it++){
      if (flag_verbose) std::cout << "PR: " << it->first.first << " " << it->first.second << " " << map_rs_f2pr[it->first].size() << "/" << map_rs_n[it->first] << " " << std::endl; 

      tmp_event += map_rs_f2pr[it->first].size();
      tmp_pot += map_rs_f2pr[it->first].size() * 1.0 / map_rs_n[it->first] * map_rs_pot[it->first];
    }
    std::cout << "Reco 2.0 PR failure in event: " << tmp_event << "/" << T_eval->GetEntries() << "  POT: " << tmp_pot << " ratio " << tmp_pot/total_pot << std::endl;
  }
  
  
  //TTree *T_bdt= (TTree*)file->Get("T_BDTvars"); // BDT tree
  //TTree *T_kine = (TTree*)file->Get("T_KINEvars"); // Kine tree for pi0 ...
  //TTree *T_eval = (TTree*)file->Get("T_eval");  // evaluation tree
  

  
  
  
  
  
  //std::cout << filename << " " << POT << std::endl;

  return 0;
}
