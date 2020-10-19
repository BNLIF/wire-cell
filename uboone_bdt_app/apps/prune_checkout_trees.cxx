#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"

using namespace std;
//using namespace LEEana;

int main( int argc, char** argv )
{
  if (argc < 2) {
    std::cout << "prune_checkout_trees #txtInputROOT #outputROOT" << std::endl;
    return -1;
  }

  string filelist(argv[1]);
  string outName(argv[2]);

  auto T_pot = new TChain("wcpselection/T_pot");
  auto T_eval = new TChain("wcpselection/T_eval");
  auto T_PFeval = new TChain("wcpselection/T_PFeval");
  auto T_BDTvars = new TChain("wcpselection/T_BDTvars");
  auto T_KINEvars = new TChain("wcpselection/T_KINEvars");


  std::ifstream in(filelist.c_str());
  std::string line;
  while (std::getline(in, line)) {
    T_pot->Add(line.c_str());
    T_eval->Add(line.c_str());
    T_PFeval->Add(line.c_str());
    T_BDTvars->Add(line.c_str());
    T_KINEvars->Add(line.c_str());
  }

  cout << T_pot->GetEntries() << endl;
  cout << T_eval->GetEntries() << endl;
  cout << T_PFeval->GetEntries() << endl;
  cout << T_BDTvars->GetEntries() << endl;
  cout << T_KINEvars->GetEntries() << endl;

  auto ofile = new TFile(outName.c_str(), "recreate");
  ofile->mkdir("wcpselection")->cd();
  T_pot->Merge(ofile, 0, "keep");
  T_eval->Merge(ofile, 0, "keep");
  T_PFeval->Merge(ofile, 0, "keep");
  T_BDTvars->Merge(ofile, 0, "keep");
  T_KINEvars->Merge(ofile, 0, "keep");
  ofile->Close();

  return 0;
}
