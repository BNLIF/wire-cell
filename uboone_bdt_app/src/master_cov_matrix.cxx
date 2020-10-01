#include "WCPLEEANA/master_cov_matrix.h"

#include <iostream>
#include <fstream>

LEEana::CovMatrix::CovMatrix(TString filename){
  std::ifstream infile(filename);
  TString name, var_name;
  Int_t bin_num;
  Float_t low_limit, high_limit;
  Int_t obs_no;
  Int_t flag_xs_flux;
  Int_t flag_det;
  Float_t flag_add;
  Int_t flag_same_mc_stat;
  Int_t cov_sec_no;
  Int_t file_no;

  int ch_no = 0;
  
  int covbin = 0;
  int obscovbin = 0;
  
  while(!infile.eof()){
    infile >> name >> var_name >> bin_num >> low_limit >> high_limit >> obs_no >> flag_xs_flux >> flag_det >> flag_add >> flag_same_mc_stat >> cov_sec_no >> file_no;
    //    std::cout << name << " " << var_name << " " << low_limit << " " << bin_num << " " << file_no << std::endl;
    if (bin_num == -1) break;
    
    map_ch_hist[ch_no] = std::make_tuple(name, var_name, bin_num, low_limit, high_limit);
    
    map_ch_filetype[ch_no] = file_no;
    map_filetype_chs[file_no].push_back(ch_no);

    map_ch_systematics[ch_no] = std::make_tuple(flag_xs_flux, flag_det, flag_add, flag_same_mc_stat);
    if (flag_xs_flux == 1) xfs_filetypes.insert(file_no);
    if (flag_det == 1) det_filetypes.insert(file_no);
    if (flag_add != 0) add_filetypes.insert(file_no);

    if (flag_same_mc_stat !=0) map_mcstat_same_covchs[flag_same_mc_stat].push_back(ch_no);
    
    // prepare for the matrix
    map_ch_obsch[ch_no] = obs_no;
    map_ch_covch[ch_no] = cov_sec_no;

    map_obsch_nbin[obs_no]     = bin_num + 1; // add the overflow bin
    map_covch_nbin[cov_sec_no] = bin_num + 1; // add the overflow bin

    map_covch_obsch[cov_sec_no] = obs_no; // ch map

    ch_no ++;
  }

  int start_bin = 0;
  for (auto it = map_obsch_nbin.begin(); it!= map_obsch_nbin.end(); it++){
    map_obsch_startbin[it->first] = start_bin;
    start_bin += it->second;
    //std::cout << it->first << " " << it->second << std::endl;
  }

  start_bin = 0;
  for (auto it = map_covch_nbin.begin(); it != map_covch_nbin.end(); it++){
    map_covch_startbin[it->first] = start_bin;
    start_bin += it->second;
  }
  
  
  // form the large covariance matrix, and the start bin ...
  for (auto it = map_covch_obsch.begin(); it!= map_covch_obsch.end(); it++){
    //std::cout << it->first << " " << it->second << std::endl;
    int nbin = map_covch_nbin[it->first];
    int start_bin_cov = map_covch_startbin[it->first];
    int start_bin_obs = map_obsch_startbin[it->second];
    for (int i=0; i!=nbin; i++){
      map_covchbin_obschbin[start_bin_cov + i] =  start_bin_obs + i;
      //std::cout << start_bin_cov + i << " " << start_bin_obs + i << std::endl;
    }
  }


  
}


LEEana::CovMatrix::~CovMatrix(){
  
}

void LEEana::CovMatrix::print_matrix(){

  for (auto it = map_ch_obsch.begin(); it != map_ch_obsch.end(); it++){
    std::cout << it->first << " " << it->second << " " << map_ch_covch[it->first] << std::endl;
  }
  
  std::cout << "covch: " << std::endl;
  for (auto it = map_covch_nbin.begin(); it != map_covch_nbin.end(); it++){
    std::cout << it->first << " " << it->second << " " << map_covch_startbin[it->first] << std::endl;
    //map_covch_nbin[it->first] = start_bin;
  }

  std::cout << "obsch: " << std::endl;
  for (auto it = map_obsch_nbin.begin(); it != map_obsch_nbin.end(); it++){
    std::cout << it->first << " " << it->second << " " << map_obsch_startbin[it->first] << std::endl;
    //map_covch_nbin[it->first] = start_bin;
  }
  
  std::cout << "matrix: " << std::endl;
  for (auto it = map_covchbin_obschbin.begin(); it != map_covchbin_obschbin.end(); it++){
    std::cout << it->first << " " << it->second << std::endl;
  }
}

void LEEana::CovMatrix::print_systematics(){
  for (auto it = map_ch_systematics.begin(); it!= map_ch_systematics.end(); it++){
    std::cout << it->first << " " << std::get<0>(it->second) << " " << std::get<1>(it->second) << " " << std::get<2>(it->second) << " " << std::get<3>(it->second) << std::endl;
  }

  std::cout << "Xs+Flux: ";
  for (auto it = xfs_filetypes.begin(); it != xfs_filetypes.end(); it++){
    std::cout << (*it) << " ";
  }
  std::cout << std::endl;

  std::cout << "Det    : ";
  for (auto it = det_filetypes.begin(); it != det_filetypes.end(); it++){
    std::cout << (*it) << " ";
  }
  std::cout << std::endl;

  
  std::cout << "Add    : ";
  for (auto it = add_filetypes.begin(); it != add_filetypes.end(); it++){
    std::cout << (*it) << " ";
  }
  std::cout << std::endl;

  std::cout << "same MC stat: " << std::endl;
  for (auto it = map_mcstat_same_covchs.begin(); it != map_mcstat_same_covchs.end(); it++){
    std::cout << it->first << ": ";
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      std::cout <<get_ch_name(*it1) << " ";
    }
    std::cout << std::endl;
  }
  
}

void LEEana::CovMatrix::print_filetype_info(){
  for (auto it = map_filetype_chs.begin(); it != map_filetype_chs.end(); it++){
    std::cout << it->first << ": ";
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      std::cout << get_ch_name(*it1) << " " ;
    }
    std::cout << std::endl;
  }
}



TString LEEana::CovMatrix::get_ch_name(int ch){
  TString result = "Null";
  if (map_ch_hist.find(ch) != map_ch_hist.end()){
    return std::get<0>(map_ch_hist[ch]);
  }
  return result;
}

TString LEEana::CovMatrix::get_ch_var(int ch){
  TString result = "Null";
  if (map_ch_hist.find(ch) != map_ch_hist.end()){
    return std::get<1>(map_ch_hist[ch]);
  }
  return result;
}




void LEEana::CovMatrix::print_ch_info(){
  for (auto it = map_ch_hist.begin(); it!= map_ch_hist.end(); it++){
    std::cout << it->first << " " << std::get<0>(it->second) << " " << std::get<1>(it->second) << " " << std::get<2>(it->second) << " " << std::get<3>(it->second) << " " << std::get<4>(it->second) << std::endl;
  }
}
