#ifndef LEEANAC_MASTER_COV_MATRIX
#define LEEANAC_MASTER_COV_MATRIX

#include "TString.h"
#include "TMatrixD.h"
#include <map>
#include <set>

namespace LEEana{
  class CovMatrix{
  public:
    CovMatrix(TString cov_filename = "./configurations/cov_input.txt", TString cv_filename = "./configurations/cv_input.txt", TString file_filename = "./configurations/file_ch.txt");
    ~CovMatrix();

    void print_ch_info();
    void print_filetype_info();
    void print_systematics();
    void print_matrix();
    
    void print_cvfile_info();

    // histogram ...
    TString get_ch_name(int ch);
    TString get_ch_var(int ch);
    int get_ch(TString name);
    
    std::tuple<int, double, double> get_ch_hist(int ch);

    // ... filetype related ...
    int get_ch_filetype(int ch);
    std::vector<int> get_filetype_chs(int filetype);

    //
    std::set<int> get_xfs_filetypes(){return xfs_filetypes;};
    std::set<int> get_det_filetypes(){return det_filetypes;};
    std::set<int> get_add_filetypes(){return add_filetypes;};

    bool get_sys_xs_flux(int ch);
    bool get_sys_det(int ch);
    std::pair<bool, float> get_sys_add(int ch);
    int get_sys_mc_same(int ch);
    std::map<int, std::vector<int> > get_mcstat_same_covchs(){return map_mcstat_same_covchs;};

    int get_obsch(int ch);
    int get_covch(int ch);
    int get_obsch_fcov(int covch);
    std::pair<int, int> get_obsch_info(int obsch);
    std::pair<int, int> get_covch_info(int covch);

    TMatrixD* get_mat_collapse(){return mat_collapse;};

    double get_ext_pot(TString filename);
    
  private:
    TMatrixD* mat_collapse;
    
    // basic information about the channels
    std::map<int, std::tuple<TString, TString, int, double, double> > map_ch_hist;
    std::map<TString, int> map_name_ch;

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

    std::map<int, int> map_covch_obsch; // map ...
    
    std::map<int, int> map_obsch_nbin; // record the bin number + 1
    std::map<int, int> map_covch_nbin; // record the bin number + 1
    
    // covariance matrix internally ...
    std::map<int, int> map_covch_startbin; 
    std::map<int, int> map_obsch_startbin;
    
    // covariance matrix in observation
    std::map<int, int> map_covchbin_obschbin;

    // CV related input ...
    std::map<int, TString> map_filetype_name;
    std::map<int, std::vector<TString> > map_filetype_inputfiles;
    std::map<TString, std::tuple<int, int, TString, double> > map_inputfile_info;
    std::map<TString, std::vector<TString> > map_inputfile_cuts;
  };
}

#endif
