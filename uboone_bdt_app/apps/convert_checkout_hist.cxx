#include <iostream>

#include "WCPLEEANA/master_cov_matrix.h"

#include "TROOT.h"
#include "TMath.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"

#include "WCPLEEANA/cuts.h"
#include "WCPLEEANA/pot.h"
#include "WCPLEEANA/pfeval.h"

using namespace std;
using namespace LEEana;

int main( int argc, char** argv )
{
  CovMatrix cov;

  TString input_filename = argv[1];
  TString out_filename = argv[2];


  bool flag_data = true;
  TFile *file = new TFile(input_filename);

  TTree *T_BDTvars = (TTree*)file->Get("wcpselection/T_BDTvars");
  TTree *T_eval = (TTree*)file->Get("wcpselection/T_eval");
  TTree *T_pot = (TTree*)file->Get("wcpselection/T_pot");
  TTree *T_PFeval = (TTree*)file->Get("wcpselection/T_PFeval");
  TTree *T_KINEvars = (TTree*)file->Get("wcpselection/T_KINEvars");

  if (T_eval->GetBranch("weight_cv")) flag_data = false;

  EvalInfo eval;
  POTInfo pot;
  TaggerInfo tagger;
  PFevalInfo pfeval;
  KineInfo kine;

#include "init.txt"

  set_tree_address(T_BDTvars, tagger,2 );
  if (flag_data){
    set_tree_address(T_eval, eval,2);
    set_tree_address(T_PFeval, pfeval,2);
  }else{
    set_tree_address(T_eval, eval);
    set_tree_address(T_PFeval, pfeval);
  }
  set_tree_address(T_pot, pot);
  set_tree_address(T_KINEvars, kine);

  double total_pot = 0;
  for (Int_t i=0;i!=T_pot->GetEntries();i++){
    T_pot->GetEntry(i);
    total_pot += pot.pot_tor875;
  }
  double ext_pot = cov.get_ext_pot(input_filename);
  if (ext_pot != 0) total_pot = ext_pot;

  std::cout << "Total POT: " << total_pot << " external POT: " << ext_pot << std::endl;
  
  return 0;
}
