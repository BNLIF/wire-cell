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
#include "TH1F.h"


#include "WCPLEEANA/master_cov_matrix.h"

using namespace std;
using namespace LEEana;

int main( int argc, char** argv )
{
  CovMatrix cov;
  // cov.print_ch_info();
  // cov.print_filetype_info();
  // cov.print_systematics();
  // cov.print_matrix();
  //  int ch = 2;
  // std::cout << cov.get_sys_xs_flux(ch) << " " << cov.get_sys_det(ch) << " " << cov.get_sys_add(ch).first << " " << cov.get_sys_add(ch).second << " " << cov.get_sys_mc_same(ch) << std::endl;
  // TFile *file = new TFile("temp.root","RECREATE");
  // TMatrixD* mat_collapse = cov.get_mat_collapse();
  // mat_collapse->Write("mat_collapse");
  // file->Write();
  // file->Close();
  //cov.print_cvfile_info();
  
  return 0;
}
