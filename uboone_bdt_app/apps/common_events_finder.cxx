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
  if (argc < 3) {
    std::cout << "./common_events_finder #file1 #file2" << std::endl;
    return -1;
  }

  TString filename1 = argv[1];
  TString filename2 = argv[2];
 
  std::cout << "Find common event lists for :" << filename1 << " " << filename2 << std::endl;


}
