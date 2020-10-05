#ifndef LEEANAC_MASTER_COV_MATRIX
#define LEEANAC_MASTER_COV_MATRIX

#include "TString.h"
#include "TMatrixD.h"
#include "TH1F.h"
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
    
    std::tuple<int, float, float> get_ch_hist(int ch);

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
    std::map<int, std::vector<int> > get_mcstat_same_chs(){return map_mcstat_same_chs;};

    int get_obsch(int ch);
    int get_covch(int ch);
    int get_obsch_fcov(int covch);
    std::pair<int, int> get_obsch_info(int obsch);
    std::pair<int, int> get_covch_info(int covch);

    TMatrixD* get_mat_collapse(){return mat_collapse;};

    float get_ext_pot(TString filename);
    std::vector< std::tuple<TString, int, float, float, TString, TString, TString, TString > > get_histograms(TString filename, int flag = 0);
    std::map<TString, std::tuple<int, int, TString, float, int> > getp_map_inputfile_info(){return map_inputfile_info;};

    std::map<TString, std::pair<TString, int> > get_map_pred_histo_histo_err2_lee(){return map_pred_histo_histo_err2_lee;};

    // Now the cross uncertainty term
    std::map<std::pair<TString, TString>, std::pair<TString, int> > get_map_pair_hist_hists_cros(){return map_pair_histo_histos_cros;};

    std::map<int, std::set<std::set<std::pair<TString, int> > > > get_map_pred_obsch_histos(){return map_pred_obsch_histos;};

    int get_obsch_name(TString name);

    void fill_data_histograms(int run, std::map<int, std::vector<TH1F*> >& map_obsch_histos, std::map<TString, std::pair<TH1F*, double> >& map_name_histogram);
    void fill_pred_histograms(int run, std::map<int, std::vector<TH1F*> >& map_obsch_histos, std::map<int, std::vector< std::vector< std::pair<double, double> > > >& map_obsch_bayes, std::map<TString, std::pair<TH1F*, double> >& map_name_histogram, float lee_strength, std::map<int, double> map_data_period_pot);
    
  private:
    TMatrixD* mat_collapse;
    
    // basic information about the channels
    // name, var_name, bin, llmit, hlimit, weight, obs_no, lee
    std::map<int, std::tuple<TString, TString, int, float, float, TString, int, int> > map_ch_hist;
    std::map<TString, int> map_name_ch;

    // information regarding ch and their filetype
    std::map<int, int> map_ch_filetype;
    std::map<int, std::vector<int> > map_filetype_chs;

    // information regarding systematics (xs_flux, det, additional, mc stat...)
    // xs_flux, det, additional relative uncertainties, mc_stat
    std::map<int, std::tuple<int, int, float, int> > map_ch_systematics;
    std::set<int> xfs_filetypes;
    std::set<int> det_filetypes;
    std::set<int> add_filetypes;

    // MC statistics ... same channels  (diagnal to start with ...)
    std::map<int, std::vector<int> > map_mcstat_same_chs;
    std::map<int, std::set<int> > map_filetype_mcstats;
    
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
    std::map<TString, int> map_inputfile_filetype;
    
    // filetype, period, outfile_name, external pot if any, file_no
    std::map<TString, std::tuple<int, int, TString, float, int> > map_inputfile_info;
    std::map<int, int> map_fileno_period;
    std::map<TString, std::vector<TString> > map_inputfile_cuts;

    // histogram infos (name, nbin, lowlimit, highlimit, variable, channel cut, additional cut, weight
    std::map<TString, std::vector< std::tuple<TString, int, float, float, TString, TString, TString, TString > > > map_inputfile_histograms;
    std::map<TString, std::vector< std::tuple<TString, int, float, float, TString, TString, TString, TString > > > map_inputfile_histograms_err2;
    std::map<TString, std::vector< std::tuple<TString, int, float, float, TString, TString, TString, TString > > > map_inputfile_histograms_cros;
    std::map<TString, TString> map_histogram_inputfile;

    // structure of summing histograms together for prediction ...
    std::map<int, std::set<int> > map_pred_obsch_covch; // OK
    std::map<int, std::set<int> > map_pred_covch_ch; // OK
    
    std::map<int, std::set<std::pair<TString, TString> > > map_pred_ch_subch; // OK
    std::map<std::pair<TString, TString> , std::set<std::pair<TString, int> > > map_pred_subch_histos; //OK
    
    
    
    std::map<TString, std::pair<TString, int> > map_pred_histo_histo_err2_lee; //OK

    // Now the cross uncertainty term
    std::map<std::pair<TString, TString>, std::pair<TString, int> > map_pair_histo_histos_cros; // OK

    std::map<int, std::set<std::set<std::pair<TString, int> > > > map_pred_obsch_histos; // total ...
    
   
    
  };
}

#endif
