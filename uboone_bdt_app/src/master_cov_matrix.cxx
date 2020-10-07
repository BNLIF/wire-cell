#include "WCPLEEANA/master_cov_matrix.h"

#include <iostream>
#include <fstream>

LEEana::CovMatrix::CovMatrix(TString cov_filename, TString cv_filename, TString file_filename){
  std::ifstream infile(cov_filename);
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
  TString weight;
  int lee_strength;
  
  while(!infile.eof()){
    infile >> name >> var_name >> bin_num >> low_limit >> high_limit >> obs_no >> flag_xs_flux >> flag_det >> flag_add >> flag_same_mc_stat >> cov_sec_no >> file_no >> weight >> lee_strength;
    //    std::cout << name << " " << var_name << " " << low_limit << " " << bin_num << " " << file_no << std::endl;
    if (bin_num == -1) break;
    
    map_ch_hist[ch_no] = std::make_tuple(name, var_name, bin_num, low_limit, high_limit, weight, obs_no, lee_strength);
    map_name_ch[name] = ch_no;
    
    map_ch_filetype[ch_no] = file_no;
    map_filetype_chs[file_no].push_back(ch_no);

    if (file_no !=5) { // data
      map_ch_systematics[ch_no] = std::make_tuple(flag_xs_flux, flag_det, flag_add, flag_same_mc_stat);
      if (flag_xs_flux == 1) xfs_filetypes.insert(file_no);
      if (flag_det == 1) det_filetypes.insert(file_no);
      if (flag_add != 0) add_filetypes.insert(file_no);
      
      if (flag_same_mc_stat !=0) {
	map_mcstat_same_chs[flag_same_mc_stat].push_back(ch_no);
	map_filetype_mcstats[file_no].insert(flag_same_mc_stat);
      }
    
      // prepare for the matrix
      map_ch_obsch[ch_no] = obs_no;
      map_ch_covch[ch_no] = cov_sec_no;
      
      map_obsch_nbin[obs_no]     = bin_num + 1; // add the overflow bin
      map_covch_nbin[cov_sec_no] = bin_num + 1; // add the overflow bin
      
      map_covch_obsch[cov_sec_no] = obs_no; // ch map

      // relation among obs_no, cov_no, ch_no ...
      map_pred_obsch_covch[obs_no].insert(cov_sec_no);
      map_pred_covch_ch[cov_sec_no].insert(ch_no);
    }
    
    ch_no ++;
  }

  int total_obs_bin = 0;
  int start_bin = 0;
  for (auto it = map_obsch_nbin.begin(); it!= map_obsch_nbin.end(); it++){
    map_obsch_startbin[it->first] = start_bin;
    start_bin += it->second;
    total_obs_bin += it->second;
    //std::cout << it->first << " " << it->second << std::endl;
  }

  int total_cov_bin = 0;
  start_bin = 0;
  for (auto it = map_covch_nbin.begin(); it != map_covch_nbin.end(); it++){
    map_covch_startbin[it->first] = start_bin;
    start_bin += it->second;
    total_cov_bin += it->second;
  }
  
  mat_collapse = new TMatrixD(total_cov_bin,total_obs_bin);
  // form the large covariance matrix, and the start bin ...
  for (auto it = map_covch_obsch.begin(); it!= map_covch_obsch.end(); it++){
    //std::cout << it->first << " " << it->second << std::endl;
    int nbin = map_covch_nbin[it->first];
    int start_bin_cov = map_covch_startbin[it->first];
    int start_bin_obs = map_obsch_startbin[it->second];
    for (int i=0; i!=nbin; i++){
      map_covchbin_obschbin[start_bin_cov + i] =  start_bin_obs + i;
      (*mat_collapse)(start_bin_cov + i, start_bin_obs + i) = 1;
      //std::cout << start_bin_cov + i << " " << start_bin_obs + i << std::endl;
    }
  }

  int filetype;
  //  TString name;
  int period;
  TString input_filename;
  TString out_filename;
  float ext_pot;
  //std::cout << cv_filename << std::endl;
  std::ifstream infile1(cv_filename);
  while(!infile1.eof()){
    infile1 >> filetype >> name >> period >> input_filename >> out_filename >> ext_pot >> file_no;
    // std::cout << filetype << " " << out_filename << std::endl;
    
    if (filetype == -1) break;
    
    map_filetype_name[filetype] = name;
    map_filetype_inputfiles[filetype].push_back(input_filename);
    map_inputfile_filetype[input_filename] = filetype;
    map_inputfile_info[input_filename] = std::make_tuple(filetype, period, out_filename, ext_pot, file_no);
    map_fileno_period[file_no] = period;
  }

  std::ifstream infile2(file_filename);
  TString cut_name;
  while (!infile2.eof()){
    infile2 >> input_filename >> cut_name;
    if (input_filename == "end") break;
    map_inputfile_cuts[input_filename].push_back(cut_name);

    // cut_name, input_filename --> filetype --> chs ...
    int filetype = map_inputfile_filetype[input_filename];
    std::vector<int> chs = map_filetype_chs[filetype];
    for (auto it = chs.begin(); it!= chs.end(); it++){
      map_pred_ch_subch[*it].insert(std::make_pair(std::get<0>(map_ch_hist[*it]),cut_name));
    }
  }


  // sort out the histograms ...
  for (auto it = map_inputfile_info.begin(); it!= map_inputfile_info.end(); it++){
    TString filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    int file_no = std::get<4>(it->second);
    for (auto it1 = map_inputfile_cuts[filename].begin(); it1 != map_inputfile_cuts[filename].end(); it1++){
      TString add_cut = *it1;
      
      for (auto it2 = map_filetype_chs[filetype].begin(); it2 != map_filetype_chs[filetype].end(); it2++){
	int ch = *it2;
	auto it3 = map_ch_hist.find(ch);

	TString name = std::get<0>(it3->second);
	TString var_name = std::get<1>(it3->second);
	int nbin = std::get<2>(it3->second);
	int llimit = std::get<3>(it3->second);
	int hlimit = std::get<4>(it3->second);
	TString weight = std::get<5>(it3->second);
	int lee_strength = std::get<7>(it3->second);

	//	std::cout << name << " " << lee_strength << std::endl;
	
	TString weight2 = weight + "_" + weight;
	TString histo_name = name + Form("_%d_",file_no) + var_name + "_" + add_cut;
	TString histo_name1 = histo_name + "_err2";

	map_histogram_inputfile[histo_name] = filename;
	map_histogram_inputfile[histo_name1] = filename;
	
	map_inputfile_histograms[filename].push_back(std::make_tuple(histo_name, nbin, llimit, hlimit, var_name, name, add_cut, weight));
	map_inputfile_histograms_err2[filename].push_back(std::make_tuple(histo_name1, nbin, llimit, hlimit, var_name, name, add_cut, weight2));

	map_pred_subch_histos[std::make_pair(name,add_cut)].insert(std::make_pair(histo_name, period));
	map_pred_histo_histo_err2_lee[histo_name] = std::make_pair(histo_name1,lee_strength);
	//std::cout << histo_name << " " << " " << histo_name1 << " " << nbin << " " << llimit << " " << hlimit << " " << var_name << " " << name << " " << add_cut << std::endl;
	
	//	std::cout << filename << " " << add_cut << std::endl;
      }
      std::set<int> mc_stat = map_filetype_mcstats[filetype];
      for (auto it5 = mc_stat.begin(); it5 != mc_stat.end(); it5++){
	std::vector<int> correlated_chs = map_mcstat_same_chs[*it5];
	//std::cout << correlated_chs.size() << " " << filetype << " " << filename << std::endl;
	for (size_t i=0;i!=correlated_chs.size();i++){
	  int ch1 = correlated_chs.at(i);
	  int obsch1 = map_ch_obsch[ch1];
	  auto it3 = map_ch_hist.find(ch1);
	  TString name1 = std::get<0>(it3->second);
	  TString var_name1 = std::get<1>(it3->second);
	  int nbin1 = std::get<2>(it3->second);
	  int llimit1 = std::get<3>(it3->second);
	  int hlimit1 = std::get<4>(it3->second);
	  TString weight1 = std::get<5>(it3->second);
	  TString histo_name1 = name1 + Form("_%d_",file_no) + var_name1 + "_" + add_cut;
	  
	  for (size_t j=i+1;j<correlated_chs.size();j++){
	    int ch2 = correlated_chs.at(j);
	    auto it4 = map_ch_hist.find(ch2);
	    TString name2 = std::get<0>(it4->second);
	    TString var_name2 = std::get<1>(it4->second);
	    int nbin2 = std::get<2>(it4->second);
	    int llimit2 = std::get<3>(it4->second);
	    int hlimit2 = std::get<4>(it4->second);
	    TString weight2 = std::get<5>(it4->second);
	    TString histo_name2 = name2 + Form("_%d_",file_no) + var_name2 + "_" + add_cut;
	    TString histo_name = name1 + "_" + name2 + "_" + add_cut + Form("_%d",file_no);
	    TString weight = weight1 +"_" + weight2;

	    map_histogram_inputfile[histo_name] = filename;
	    
	    map_inputfile_histograms_cros[filename].push_back(std::make_tuple(histo_name, nbin1, llimit1, hlimit1, var_name1, name1, add_cut, weight));
	    
	    map_pair_histo_histos_cros[std::make_pair(histo_name1, histo_name2)] = std::make_pair(histo_name, obsch1);
	    
	  }
	}
      }
      
    }
  }


  // now form the final prediction map ...
  for (auto it = map_pred_obsch_covch.begin(); it!= map_pred_obsch_covch.end(); it++){
    int obsch = it->first;
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      int covch = *it1;
      std::set<int> chs = map_pred_covch_ch[covch];
      for (auto it2 = chs.begin(); it2 != chs.end(); it2++){
	int ch = *it2;
	std::set<std::pair<TString, TString> > subchs = map_pred_ch_subch[ch];
	for (auto it3 = subchs.begin(); it3 != subchs.end(); it3++){
	  std::pair<TString, TString> subch = *it3;
	  std::set<std::pair<TString, int> > histos = map_pred_subch_histos[subch];

	  map_pred_obsch_histos[obsch].insert(histos);
	  //for (auto it4 = histos.begin(); it4 != histos.end(); it4++){
	  //  TString histo = *it4;
	  //map_pred_obsch_histos[obsch].insert(histo);
	  //}
	  
	}
      }
    }
  }

  // for (auto it = map_pred_obsch_histos.begin(); it!=map_pred_obsch_histos.end();it++){
  //   std::cout << it->first << std::endl;
  //   for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
  //     std::cout << "sub: " << (*it1).size() << std::endl;
  //     for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
  // 	std::cout << it->first << " " << (*it2).first << " " << (*it2).second << std::endl;
  //     }
  //   }
  // }

  llimit = 0;
  hlimit = 100;
  
  Double_t x[101]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100};

  
  Double_t yl[101]={0, 0, 0, 0.856632, 1.70317, 2.51005, 3.32075, 4.14046, 4.9693, 5.80646, 6.65117, 7.5025, 8.35978, 9.22237, 10.0898, 10.9615, 11.8372, 12.7165, 13.5992, 14.4849, 15.3734, 16.2646, 17.1583, 18.0543, 18.9524, 19.8526, 20.7547, 21.6586, 22.5642, 23.4715, 24.3803, 25.2906, 26.2023, 27.1153, 28.0297, 28.9452, 29.8619, 30.7797, 31.6987, 32.6187, 33.5396, 34.4616, 35.3845, 36.3083, 37.2329, 38.1584, 39.0847, 40.0118, 40.9396, 41.8682, 42.7975, 43.7275, 44.6581, 45.5895, 46.5215, 47.454, 48.3873, 49.321, 50.2554, 51.1903, 52.1257, 53.0617, 53.9982, 54.9352, 55.8727, 56.8107, 57.7491, 58.6881, 59.6274, 60.5673, 61.5075, 62.4482, 63.3892, 64.3307, 65.2725, 66.2148, 67.1575, 68.1005, 69.0438, 69.9876, 70.9317, 71.8761, 72.8209, 73.766, 74.7114, 75.6572, 76.6033, 77.5497, 78.4964, 79.4434, 80.3907, 81.3383, 82.2862, 83.2342, 84.1827, 85.1314, 86.0804, 87.0296, 87.9791, 88.9288, 89.8788}; 
  Double_t yh[101]={1.1478, 2.35971, 3.51917, 4.72422, 5.98186, 7.21064, 8.41858, 9.61053, 10.7896, 11.9582, 13.1179, 14.27, 15.4155, 16.5552, 17.6898, 18.8197, 19.9454, 21.0673, 22.1858, 23.3011, 24.4133, 25.5229, 26.6299, 27.7346, 28.837, 29.9374, 31.0358, 32.1322, 33.2271, 34.3201, 35.4117, 36.5017, 37.5904, 38.6776, 39.7635, 40.8483, 41.9318, 43.0141, 44.0955, 45.1757, 46.2549, 47.3331, 48.4104, 49.4868, 50.5623, 51.637, 52.7108, 53.7839, 54.8561, 55.9277, 56.9985, 58.0686, 59.1381, 60.2068, 61.275, 62.3425, 63.4094, 64.4757, 65.5415, 66.6066, 67.6713, 68.7354, 69.7989, 70.862, 71.9246, 72.9866, 74.0483, 75.1094, 76.1701, 77.2304, 78.2902, 79.3496, 80.4085, 81.4672, 82.5253, 83.5831, 84.6406, 85.6976, 86.7542, 87.8105, 88.8665, 89.9221, 90.9774, 92.0323, 93.0869, 94.1411, 95.1951, 96.2488, 97.3021, 98.3552, 99.4079, 100.46, 101.513, 102.564, 103.616, 104.667, 105.718, 106.769, 107.82, 108.87, 109.92};

  for (Int_t i=0;i!=101;i++){
    yl[i] = fabs(yl[i] - x[i]);
    yh[i] = fabs(yh[i] - x[i]);
  }
  
  gl = new TGraph(101,x,yl);
  gh = new TGraph(101,x,yh);
  
}


LEEana::CovMatrix::~CovMatrix(){

  delete mat_collapse;
  delete gl;
  delete gh;
}

std::pair<double, double> LEEana::CovMatrix::get_bayes_errors(double num){
  if (num > hlimit) return std::make_pair(sqrt(num), sqrt(num));
  else if (num < llimit) num = llimit;

  //  std::cout << num << " " << gl->Eval(num) << " " << gh->Eval(num) << std::endl;
  
  return std::make_pair(gl->Eval(num),gh->Eval(num));
}


std::vector< std::tuple<TString, int, float, float, TString, TString, TString, TString > > LEEana::CovMatrix::get_histograms(TString filename, int flag){
  if (flag == 0){
    auto it = map_inputfile_histograms.find(filename);
    if (it != map_inputfile_histograms.end()){
      return it->second;
    }else{
      std::vector< std::tuple<TString, int, float, float, TString, TString, TString, TString > > temp;
      return temp;
    }
  }else if (flag==1){
    auto it = map_inputfile_histograms_err2.find(filename);
    if (it != map_inputfile_histograms_err2.end()){
      return it->second;
    }else{
      std::vector< std::tuple<TString, int, float, float, TString, TString, TString, TString > > temp;
      return temp;
    }
  }else if (flag==2){
    auto it = map_inputfile_histograms_cros.find(filename);
    if (it != map_inputfile_histograms_cros.end()){
      return it->second;
    }else{
      std::vector< std::tuple<TString, int, float, float, TString, TString, TString, TString > > temp;
      return temp;
    }
  }
}

void LEEana::CovMatrix::fill_data_histograms(int run, std::map<int, std::vector<TH1F*> >& map_obsch_histos, std::map<TString, std::pair<TH1F*, double> >& map_name_histogram){
  for (auto it = map_inputfile_info.begin(); it != map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    TString out_filename = std::get<2>(it->second);
    int file_no = std::get<4>(it->second);
    if (period == run || run == 0){
      if (filetype == 5){
	std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = get_histograms(input_filename,0);
	for (auto it1 = histo_infos.begin(); it1 != histo_infos.end(); it1++){
	  int obsch = get_obsch_name(std::get<5>(*it1));	
	  TH1F *htemp = map_name_histogram[std::get<0>(*it1)].first;
	  map_obsch_histos[obsch].at(0)->Add(htemp);
	}
      }
    }
  }
  
}


void LEEana::CovMatrix::fill_pred_histograms(int run, std::map<int, std::vector<TH1F*> >& map_obsch_histos, std::map<int, std::vector< std::vector< std::tuple<double, double, double, int> > > >& map_obsch_bayes, std::map<TString, std::pair<TH1F*, double> >& map_name_histogram, float lee_strength, std::map<int, double> map_data_period_pot){
  
  for (auto it = map_pred_obsch_histos.begin(); it!=map_pred_obsch_histos.end();it++){
    //std::cout << it->first << std::endl;
    int obsch = it->first;
    TH1F *hpred = map_obsch_histos[obsch].at(1);
    TH1F *hpred_err2 = map_obsch_histos[obsch].at(2);

    std::map<TString, std::pair<double, double> > temp_map_histo_ratios;

    std::map<TString, std::vector< std::tuple<double, double, double, int> > > map_histoname_values;
    
    // group
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      //std::cout << "sub: " << (*it1).size() << std::endl;

      TH1F *htemp = (TH1F*)hpred->Clone("htemp");
      htemp->Reset();
      TH1F *htemp_err2 = (TH1F*)hpred_err2->Clone("htemp_err2");
      htemp_err2->Reset();
      
      std::map<int, double> temp_map_data_pot;
      std::map<int, double> temp_map_mc_acc_pot;

      // get accumulated POT ...
      for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
	TString histoname = (*it2).first;
	
	TString input_filename = map_histogram_inputfile[histoname];
	auto it3 = map_inputfile_info.find(input_filename);
	int period = std::get<1>(it3->second);  if (period != run && run !=0) continue; // skip ...
	auto it_data = map_data_period_pot.find(period); if (it_data == map_data_period_pot.end()) continue; // no corresponding data ...
	
	double data_pot = it_data->second;
	temp_map_data_pot[period] = data_pot;
	
	std::pair<TString, int> err2_lee = map_pred_histo_histo_err2_lee[histoname];
	TString histoname_err2 = err2_lee.first;
	int flag_lee = err2_lee.second;
	if (flag_lee == 1 && lee_strength == 0) continue; // no need to add ...
	
	TH1F *hmc = map_name_histogram[histoname].first;
	TH1F *hmc_err2 = map_name_histogram[histoname_err2].first;
	double mc_pot = map_name_histogram[histoname].second;
	if (temp_map_mc_acc_pot.find(period) == temp_map_mc_acc_pot.end()){
	  temp_map_mc_acc_pot[period] = mc_pot;
	}else{
	  temp_map_mc_acc_pot[period] += mc_pot;
	}
      }

      for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
	TString histoname = (*it2).first;
	TString input_filename = map_histogram_inputfile[histoname];
	auto it3 = map_inputfile_info.find(input_filename);
	int period = std::get<1>(it3->second);  if (period != run && run !=0) continue; // skip ...
	auto it_data = map_data_period_pot.find(period); if (it_data == map_data_period_pot.end()) continue; // no corresponding data ...
	std::pair<TString, int> err2_lee = map_pred_histo_histo_err2_lee[histoname];
	TString histoname_err2 = err2_lee.first;
	int flag_lee = err2_lee.second;
	if (flag_lee == 1 && lee_strength == 0) continue; // no need to add ...
	TH1F *hmc = map_name_histogram[histoname].first;
	TH1F *hmc_err2 = map_name_histogram[histoname_err2].first;
	
	double ratio = temp_map_data_pot[period]/temp_map_mc_acc_pot[period];
	
	temp_map_histo_ratios[histoname] = std::make_pair(ratio, lee_strength);
	
	if (flag_lee == 1) ratio *= lee_strength;
	
	htemp->Add(hmc, ratio);
	htemp_err2->Add(hmc_err2, ratio*ratio);
	
	std::vector< std::tuple<double, double, double, int> > values;
	for (int i=0;i!=hmc->GetNbinsX()+1;i++){
	  values.push_back(std::make_tuple(hmc->GetBinContent(i+1)*ratio , hmc_err2->GetBinContent(i+1)*ratio*ratio, temp_map_data_pot[period]/temp_map_mc_acc_pot[period], 0));
	}
	map_histoname_values[histoname] = values;

	//	if (obsch==1) std::cout << values.at(2).first << " sep " << values.at(2).second << " " << histoname << std::endl;
	//	std::cout << histoname << " " << period << " " << temp_map_data_pot[period] << " " << temp_map_mc_acc_pot[period] << " " << htemp->GetSum() << " " << ratio << std::endl;
      }

      

      
      hpred->Add(htemp);
      hpred_err2->Add(htemp_err2);
      // delete the histograms ...
      delete htemp;
      delete htemp_err2;
    }

    //    std::cout << map_histoname_values.size() << std::endl;
    // treat cross term ???
    
    for (auto it2 = map_pair_histo_histos_cros.begin(); it2 != map_pair_histo_histos_cros.end(); it2++){
      TString hist1 = it2->first.first;
      TString hist2 = it2->first.second;
      TString histogram = it2->second.first;
      int obsch1 = it2->second.second;
      if (obsch1 != obsch) continue;
      auto it3 = temp_map_histo_ratios.find(hist1);
      auto it4 = temp_map_histo_ratios.find(hist2);
      if (it3 == temp_map_histo_ratios.end() || it4 == temp_map_histo_ratios.end()) continue;
      TH1F *hmc_cros = map_name_histogram[histogram].first;
      // cross term ...
      hpred_err2->Add(hmc_cros, it3->second.first * it4->second.first * 2 * it3->second.second);
      
      std::vector< std::tuple<double, double, double, int> >& values1 = map_histoname_values[hist1];
      std::vector< std::tuple<double, double, double, int> >& values2 = map_histoname_values[hist2];
      for (size_t j=0;j!=values1.size();j++){
	std::get<0>(values1.at(j)) += std::get<0>(values2.at(j));
	std::get<1>(values1.at(j)) += std::get<1>(values2.at(j)) + hmc_cros->GetBinContent(j+1) * it3->second.first * it4->second.first * 2 * it3->second.second; 
      }
      map_histoname_values.erase(hist2);
      //      std::cout << values1.size() << " " << values2.size() << std::endl;
      //      std::cout << hist1 << " " << hist2 << " " << histogram << " " << it3->second.first << " " << it3->second.second << " " << it4->second.first << " " << it4->second.second << std::endl;
    }
    //   std::cout <<  map_histoname_values.size() << std::endl;

    for (auto it2 = map_histoname_values.begin(); it2 != map_histoname_values.end(); it2++){
      //      if (obsch==1) std::cout << it2->second.at(2).first << " " << it2->second.at(2).second << std::endl;
      map_obsch_bayes[obsch].push_back(it2->second);
    }
    
  } // work on obs channel ...
}



int LEEana::CovMatrix::get_obsch_name(TString name){
  int ch = map_name_ch[name];
  // std::cout << ch << std::endl;
  return std::get<6>(map_ch_hist[ch]);
}

float LEEana::CovMatrix::get_ext_pot(TString filename){
  auto it = map_inputfile_info.find(filename);
  if (it != map_inputfile_info.end()){
    return std::get<3>(it->second);
  }else{
    return 0;
  }
}

void LEEana::CovMatrix::print_cvfile_info(){
  
  for (auto it = map_filetype_inputfiles.begin(); it!= map_filetype_inputfiles.end(); it++){
    std::cout << it->first << " " << map_filetype_name[it->first] << std::endl;
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      std::cout << *it1 << " " << std::get<0>(map_inputfile_info[*it1]) << " " << std::get<1>(map_inputfile_info[*it1]) << " " << std::get<2>(map_inputfile_info[*it1]) << " " << std::get<3>(map_inputfile_info[*it1]) << " " << map_inputfile_cuts[*it1].size() << std::endl;
    }
  }
}


int LEEana::CovMatrix::get_ch(TString name){
  auto it = map_name_ch.find(name);
  if (it != map_name_ch.end()){
    return it->second;
  }else{
    return -1;
  }
}


int LEEana::CovMatrix::get_obsch(int ch){
  auto it = map_ch_obsch.find(ch);
  if (it != map_ch_obsch.end()){
    return it->second;
  }else{
    return -1;
  }
}

int LEEana::CovMatrix::get_covch(int ch){
  auto it = map_ch_covch.find(ch);
  if (it != map_ch_covch.end()){
    return it->second;
  }else{
    return -1;
  }
}

int LEEana::CovMatrix::get_obsch_fcov(int covch){
  auto it = map_covch_obsch.find(covch);
  if (it != map_covch_obsch.end()){
    return it->second;
  }else{
    return -1;
  }
}

std::pair<int, int> LEEana::CovMatrix::get_obsch_info(int obsch){
  auto it = map_obsch_startbin.find(obsch);
  if (it != map_obsch_startbin.end()){
    return std::make_pair(it->second, map_obsch_nbin[it->first]);
  }else{
    return std::make_pair(0,0);
  }
}
std::pair<int, int> LEEana::CovMatrix::get_covch_info(int covch){
  auto it = map_covch_startbin.find(covch);
  if (it != map_covch_startbin.end()){
    return std::make_pair(it->second, map_covch_nbin[it->first]);
  }else{
    return std::make_pair(0,0);
  }
}


bool LEEana::CovMatrix::get_sys_xs_flux(int ch){
  auto it = map_ch_systematics.find(ch);
  if (it != map_ch_systematics.end()){
    if (std::get<0>(it->second)==0){
      return false;
    }else{
      return true;
    }
  }else{
    return false;
  }
}

bool LEEana::CovMatrix::get_sys_det(int ch){
  auto it = map_ch_systematics.find(ch);
  if (it != map_ch_systematics.end()){
    if (std::get<1>(it->second)==0){
      return false;
    }else{
      return true;
    }
  }else{
    return false;
  }
}
std::pair<bool, float> LEEana::CovMatrix::get_sys_add(int ch){
  auto it = map_ch_systematics.find(ch);
  if (it != map_ch_systematics.end()){
    if (std::get<2>(it->second)==0){
      return std::make_pair(false,0);
    }else{
      return std::make_pair(true,std::get<2>(it->second) );
    }
  }else{
    return std::make_pair(false,0);
  }
}

int LEEana::CovMatrix::get_sys_mc_same(int ch){
  auto it = map_ch_systematics.find(ch);
  if (it != map_ch_systematics.end()){
    return std::get<3>(it->second);
  }else{
    return 0;
  }
}

std::vector<int> LEEana::CovMatrix::get_filetype_chs(int filetype){
  auto it = map_filetype_chs.find(filetype);
  if (it != map_filetype_chs.end()){
    return it->second;
  }else{
    std::vector<int> tmp;
    return tmp;
  }
}

int LEEana::CovMatrix::get_ch_filetype(int ch){
  auto it = map_ch_filetype.find(ch);

  if (it != map_ch_filetype.end()){
    return it->second;
  }else{
    return -1;
  }
}

std::tuple<int, float, float> LEEana::CovMatrix::get_ch_hist(int ch){
  auto it = map_ch_hist.find(ch);
  if (it != map_ch_hist.end()){
    return std::make_tuple(std::get<2>(it->second), std::get<3>(it->second), std::get<4>(it->second));
  }else{
    return std::make_tuple(-1,0,0);
  }
}


TString LEEana::CovMatrix::get_ch_name(int ch){
  TString result = "Null";
  auto it = map_ch_hist.find(ch);
  if (it != map_ch_hist.end()){
    return std::get<0>(it->second);
  }
  return result;
}

TString LEEana::CovMatrix::get_ch_var(int ch){
  TString result = "Null";
  auto it = map_ch_hist.find(ch);
  if (it != map_ch_hist.end()){
    return std::get<1>(it->second);
  }
  return result;
}




void LEEana::CovMatrix::print_ch_info(){
  for (auto it = map_ch_hist.begin(); it!= map_ch_hist.end(); it++){
    std::cout << it->first << " " << std::get<0>(it->second) << " " << std::get<1>(it->second) << " " << std::get<2>(it->second) << " " << std::get<3>(it->second) << " " << std::get<4>(it->second) << std::endl;
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
  for (auto it = map_mcstat_same_chs.begin(); it != map_mcstat_same_chs.end(); it++){
    std::cout << it->first << ": ";
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      std::cout <<get_ch_name(*it1) << " ";
    }
    std::cout << std::endl;
  }
  
}
