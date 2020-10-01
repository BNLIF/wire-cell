#ifndef LEEANAC_MASTER_COV_MATRIX
#define LEEANAC_MASTER_COV_MATRIX

#include "TString.h"
#include <map>

namespace LEEana{
  class CovMatrix{
  public:
    CovMatrix(TString filename = "./configurations/cov_input.txt");
    ~CovMatrix();

  private:
    // basic information about the channels
    std::map<int, TString> map_ch_name;
    std::map<int, TString> map_ch_var;
    std::map<int, std::tuple<int, double, double> > map_ch_hist;

    // information regarding ch and their filetype
    std::map<int, int > map_ch_filetype;
    std::map<int, std::vector<int> > map_filetype_chs;

    // information regarding systematics (xs_flux, det, additional, mc stat...)
    std::map<int, std::tuple<int, int, int, int> > map_ch_systematics;
    std::map<int, std::vector<int> > map_xfs_filetypes;
    std::map<int, std::vector<int> > map_det_filetypes;
    std::map<int, std::vector<int> > map_add_filetypes;

    // MC statistics ... same channels  (diagnal to start with ...)
    std::map<int, std::vector<int> > map_mcstat_same_covchs;
    
    // prepare covariance matrix structure ...
    std::map<int, int > map_ch_obsch;
    std::map<int, int > map_ch_covch;

    // covariance matrix in observation
    std::map<int, std::vector<int> > map_obsch_obscovbins;
    std::map<int, int> map_obscovbin_obsch;
  
    // covariance matrix internally ...
    std::map<int, std::vector<int> > map_covch_covbins;
    std::map<int, int> map_covbin_covch;

    // essentially the matrix ...
    std::map<int, int > map_covchbin_obschbin;

    
  };
}

#endif
