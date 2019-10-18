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
  std::cout << T->GetEntries() << std::endl;
}
