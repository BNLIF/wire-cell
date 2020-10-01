#ifndef LEEANAC_MASTER_COV_MATRIX
#define LEEANAC_MASTER_COV_MATRIX

#include "TString.h"
#include <map>
#include <set>

namespace LEEana{
  class CovMatrix{
  public:
    CovMatrix(TString filename = "./configurations/cov_input.txt");
    ~CovMatrix();

    void print_ch_info();
    void print_filetype_info();
    void print_systematics();
    void print_matrix();
    
    TString get_ch_name(int ch);
    TString get_ch_var(int ch);
    
    std::tuple<int, double, double> get_ch_hist(int ch);
    
  private:
    // basic information about the channels
    std::map<int, std::tuple<TString, TString, int, double, double> > map_ch_hist;
    // information regarding ch and their filetype
    std::map<int, int> map_ch_filetype;
    std::map<int, std::vector<int> > map_filetype_chs;

    // information regarding systematics (xs_flux, det, additional, mc stat...)
    std::map<int, std::tuple<int, int, float, int> > map_ch_systematics;
    std::set<int> xfs_filetypes;
    std::set<int> det_filetypes;
    std::set<int> add_filetypes;

    // MC statistics ... same channels  (diagnal to start with ...)
    std::map<int, std::vector<int> > map_mcstat_same_covchs;
    
    // prepare covariance matrix structure ...
    std::map<int, int> map_ch_obsch;
    std::map<int, int> map_ch_covch;

    std::map<int, int> map_obsch_nbin; // record the bin number + 1
    std::map<int, int> map_covch_nbin; // record the bin number + 1

    std::map<int, int> map_covch_obsch; // map ...

    // covariance matrix internally ...
    std::map<int, int> map_covch_startbin; 
    std::map<int, int> map_obsch_startbin;
    
    // covariance matrix in observation
    std::map<int, int> map_covchbin_obschbin;
    
  };
}

#endif
