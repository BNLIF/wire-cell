#include <cstdlib>
#include <iostream>
#include <fstream>
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

#include "WCPLEEANA/pot.h"

using namespace std;
using namespace LEEana;



int main(int argc, char** argv){

  if (argc < 2){
    std::cout << "pot_counting #bnb_file #extbnb_file -m[1 for BNB, 2 for both EXT and BNB] " << std::endl;
    return -1;
  }

  int mode = 1;
  for (Int_t i=1;i!=argc;i++){
    switch(argv[i][1]){
    case 'm':
      mode = atoi(&argv[i][2]);
      break;
    }
  }
  
  TString bnb_file = argv[1];
  TString extbnb_file = argv[2];

  int run, subrun;
  double trigger_no, pot;
  double total_bnb_trigger_no = 0, total_bnb_pot=0;
  double total_extbnb_trigger_no = 0;
  // calculate BNB number ...
  if (mode==1 || mode==2){
    std::map<std::pair<int, int>, std::pair<double, double> >  map_bnb_infos;
    ifstream infile("pot_bnb.txt");
    while(!infile.eof()){
      infile >> run >> subrun >> trigger_no >> pot;
      if (run >0)
	map_bnb_infos[std::make_pair(run,subrun)] = std::make_pair(trigger_no, pot);
    }
    TFile *file1 = new TFile(bnb_file);
    TTree *T_pot = (TTree*)file1->Get("wcpselection/T_pot");
    POTInfo pot;
    set_tree_address(T_pot, pot);
    float pass_ratio = 1;
    if (T_pot->GetBranch("pass_ratio")) T_pot->SetBranchAddress("pass_ratio", &pass_ratio);
   
    for (Int_t i=0;i!=T_pot->GetEntries();i++){
      T_pot->GetEntry(i);
      auto it = map_bnb_infos.find(std::make_pair(pot.runNo, pot.subRunNo));
      if (it == map_bnb_infos.end()){
	std::cout << "Run: " << pot.runNo << " subRun: " << pot.subRunNo << "not found!" << std::endl;
      }else{
	total_bnb_trigger_no += it->second.first * pass_ratio;
	total_bnb_pot += it->second.second * pass_ratio;
      }
    }

    std::cout << "BNB: total POT: " << total_bnb_pot << " total trigger: " << total_bnb_trigger_no << std::endl;
  }

  // calculate EXTBNB number ... 
  if (mode==2){
    std::map<std::pair<int, int>, double >  map_extbnb_infos;
    ifstream infile1("pot_extbnb.txt");
    while(!infile1.eof()){
      infile1 >> run >> subrun >> trigger_no;
      if (run >0)
	map_extbnb_infos[std::make_pair(run, subrun)] = trigger_no;
    }

    TFile *file2 = new TFile(extbnb_file);
    TTree *T_pot = (TTree*)file2->Get("wcpselection/T_pot");
    POTInfo pot;
    set_tree_address(T_pot, pot);
    float pass_ratio=1;
    if (T_pot->GetBranch("pass_ratio")) T_pot->SetBranchAddress("pass_ratio", &pass_ratio);
    
   
    for (Int_t i=0;i!=T_pot->GetEntries();i++){
      T_pot->GetEntry(i);
      auto it = map_extbnb_infos.find(std::make_pair(pot.runNo, pot.subRunNo));
      
      if (it == map_extbnb_infos.end()){
	std::cout << "Run: " << pot.runNo << " subRun: " << pot.subRunNo << "  not found!" << std::endl;
      }else{
	total_extbnb_trigger_no += it->second * pass_ratio;
      }
    }

    std::cout << "EXTBNB: total trigger: " << total_extbnb_trigger_no << std::endl;
    std::cout << "EXTBNB: total_pot: " << total_extbnb_trigger_no/total_bnb_trigger_no * total_bnb_pot << std::endl;
  }
  
  
}
