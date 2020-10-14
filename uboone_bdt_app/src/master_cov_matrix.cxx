#include "WCPLEEANA/master_cov_matrix.h"

#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TPrincipal.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TRandom3.h"

#include "WCPLEEANA/cuts.h"
#include "WCPLEEANA/pot.h"
#include "WCPLEEANA/pfeval.h"
#include "WCPLEEANA/weights.h"

//#include "WCPLEEANA/eval.h"
//#include "WCPLEEANA/kine.h"
//#include "WCPLEEANA/tagger.h"





using namespace LEEana;


float leeweight(float Enu)
 {
         if(Enu<200 || Enu>=800) return 0.0;
         else if(Enu>=200 && Enu<250) return 6.37441;
         else if(Enu>=250 && Enu<300) return 5.64554;
         else if(Enu>=300 && Enu<350) return 3.73055;
         else if(Enu>=350 && Enu<400) return 1.50914;
         else if(Enu>=400 && Enu<450) return 1.07428;
         else if(Enu>=450 && Enu<500) return 0.754093;
         else if(Enu>=500 && Enu<600) return 0.476307;
         else if(Enu>=600 && Enu<800) return 0.152327;
         else return 0.0;
 }

#include "mcm_1.h"

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
  mat_add_cov = new TMatrixD(total_cov_bin, total_cov_bin);
  
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
  double norm_pot;
  //std::cout << cv_filename << std::endl;
  std::ifstream infile1(cv_filename);
  while(!infile1.eof()){
    infile1 >> filetype >> name >> period >> input_filename >> out_filename >> ext_pot >> file_no >> norm_pot;
    //std::cout << filetype << " " << out_filename << std::endl;
    
    if (filetype == -1) break;
    
    map_filetype_name[filetype] = name;
    map_filetype_inputfiles[filetype].push_back(input_filename);
    map_inputfile_filetype[input_filename] = filetype;
    map_inputfile_info[input_filename] = std::make_tuple(filetype, period, out_filename, ext_pot, file_no, norm_pot);
    map_fileno_period[file_no] = period;
  }

  std::ifstream infile2(file_filename);
  TString cut_name;
  while (!infile2.eof()){
    infile2 >> input_filename >> cut_name;
    if (input_filename == "end") break;
    map_inputfile_cuts[input_filename].push_back(cut_name);

    //    std::cout << input_filename << " " << cut_name << std::endl;
    
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

      //      std::cout << filename << " " << add_cut << std::endl;
      
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
	
	map_histogram_covch_add[histo_name] = std::make_pair(map_ch_covch[map_name_ch[name]], std::get<2>(map_ch_systematics[map_name_ch[name]]));
	map_histogram_covch_add[histo_name1] = std::make_pair(map_ch_covch[map_name_ch[name]], std::get<2>(map_ch_systematics[map_name_ch[name]]));
	
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
	    map_histogram_covch_add[histo_name] = std::make_pair(-1,0);
	    
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
	  map_pred_covch_histos[covch].insert(histos);
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
  
  Double_t x1[101]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
		    11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
		    21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
		    31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
		    41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
		    51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
		    61, 62, 63, 64, 65, 66, 67, 68, 69, 70,
		    71, 72, 73, 74, 75, 76, 77, 78, 79, 80,
		    81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
		    91, 92, 93, 94, 95, 96, 97, 98, 99, 100};

  
  Double_t yl[101]={0, 0, 0, 0.856632, 1.70317, 2.51005, 3.32075, 4.14046, 4.9693, 5.80646, 6.65117,
		    7.5025, 8.35978, 9.22237, 10.0898, 10.9615, 11.8372, 12.7165, 13.5992, 14.4849, 15.3734,
		    16.2646, 17.1583, 18.0543, 18.9524, 19.8526, 20.7547, 21.6586, 22.5642, 23.4715, 24.3803,
		    25.2906, 26.2023, 27.1153, 28.0297, 28.9452, 29.8619, 30.7797, 31.6987, 32.6187, 33.5396,
		    34.4616, 35.3845, 36.3083, 37.2329, 38.1584, 39.0847, 40.0118, 40.9396, 41.8682, 42.7975,
		    43.7275, 44.6581, 45.5895, 46.5215, 47.454, 48.3873, 49.321, 50.2554, 51.1903, 52.1257,
		    53.0617, 53.9982, 54.9352, 55.8727, 56.8107, 57.7491, 58.6881, 59.6274, 60.5673, 61.5075,
		    62.4482, 63.3892, 64.3307, 65.2725, 66.2148, 67.1575, 68.1005, 69.0438, 69.9876, 70.9317,
		    71.8761, 72.8209, 73.766, 74.7114, 75.6572, 76.6033, 77.5497, 78.4964, 79.4434, 80.3907,
		    81.3383, 82.2862, 83.2342, 84.1827, 85.1314, 86.0804, 87.0296, 87.9791, 88.9288, 89.8788}; 
  Double_t yh[101]={1.1478, 2.35971, 3.51917, 4.72422, 5.98186, 7.21064, 8.41858, 9.61053, 10.7896, 11.9582, 13.1179,
		    14.27, 15.4155, 16.5552, 17.6898, 18.8197, 19.9454, 21.0673, 22.1858, 23.3011, 24.4133,
		    25.5229, 26.6299, 27.7346, 28.837, 29.9374, 31.0358, 32.1322, 33.2271, 34.3201, 35.4117,
		    36.5017, 37.5904, 38.6776, 39.7635, 40.8483, 41.9318, 43.0141, 44.0955, 45.1757, 46.2549,
		    47.3331, 48.4104, 49.4868, 50.5623, 51.637, 52.7108, 53.7839, 54.8561, 55.9277, 56.9985,
		    58.0686, 59.1381, 60.2068, 61.275, 62.3425, 63.4094, 64.4757, 65.5415, 66.6066, 67.6713,
		    68.7354, 69.7989, 70.862, 71.9246, 72.9866, 74.0483, 75.1094, 76.1701, 77.2304, 78.2902,
		    79.3496, 80.4085, 81.4672, 82.5253, 83.5831, 84.6406, 85.6976, 86.7542, 87.8105, 88.8665,
		    89.9221, 90.9774, 92.0323, 93.0869, 94.1411, 95.1951, 96.2488, 97.3021, 98.3552, 99.4079,
		    100.46, 101.513, 102.564, 103.616, 104.667, 105.718, 106.769, 107.82, 108.87, 109.92};

  for (Int_t i=0;i!=101;i++){
    yl[i] = fabs(yl[i] - x1[i]);
    yh[i] = fabs(yh[i] - x1[i]);
  }
  
  gl = new TGraph(101,x1,yl);
  gh = new TGraph(101,x1,yh);
  
}


LEEana::CovMatrix::~CovMatrix(){

  //std::cout << "haha " << std::endl;
  delete mat_collapse;
  delete mat_add_cov;
  delete gl;
  delete gh;

  //  std::cout << "hehe " << std::endl;
}

void LEEana::CovMatrix::gen_xf_cov_matrix(int run, std::map<int, TH1F*>& map_covch_hist, std::map<TString, TH1F*>& map_histoname_hist, TVectorD* vec_mean,  TMatrixD* cov_xf_mat){
  // prepare the maps ... name --> no,  covch, lee
  std::map<TString, std::tuple<int, int, int, TString>> map_histoname_infos ; 
  std::map<int, TString> map_no_histoname; 
  
  int ncount = 0;
  for (auto it = map_inputfile_info.begin(); it != map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    
    if (period != run) continue;
    TString out_filename = std::get<2>(it->second);
    int file_no = std::get<4>(it->second);
    std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = get_histograms(input_filename,0);
    
    for (auto it1 = histo_infos.begin(); it1 != histo_infos.end(); it1++){
      int ch = map_name_ch[std::get<5>(*it1)];
      int obsch = get_obsch_name(std::get<5>(*it1));
      int covch = get_covch_name(std::get<5>(*it1));
      int flag_lee = std::get<7>(map_ch_hist[ch]);
      TString histoname = std::get<0>(*it1);
      TH1F *htemp = map_histoname_hist[histoname];
      //
      map_histoname_infos[histoname] = std::make_tuple(ncount, covch, flag_lee, input_filename);
      map_no_histoname[ncount] = histoname;
      ncount ++;

      // std::cout << histoname << " " << obsch << " " << covch << " " << flag_lee << std::endl;
    }
  }

  // now prepare the output ...
  // filename ...   # events               #weight  #leeweight #difference   #different types
  std::map<TString, std::set<std::tuple<float, float, std::vector<float>, std::vector<int>, std::set<std::pair<int, float> > > > > map_passed_events; // last one is variable name ...
  std::map<TString, double> map_filename_pot;
  for (auto it = map_inputfile_info.begin(); it != map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    //int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    if (period != run) continue;

    //map_all_events[input_filename];
    get_events_weights(input_filename, map_passed_events, map_filename_pot, map_histoname_infos);      
  }
  
  
  
  
}


void LEEana::CovMatrix::get_events_weights(TString input_filename, std::map<TString, std::set<std::tuple<float, float, std::vector<float>, std::vector<int>, std::set<std::pair<int, float> > > > >& map_passed_events, std::map<TString, double>& map_filename_pot, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos){
  TFile *file = new TFile(input_filename);

  TTree *T_BDTvars = (TTree*)file->Get("wcpselection/T_BDTvars");
  TTree *T_eval = (TTree*)file->Get("wcpselection/T_eval");
  TTree *T_pot = (TTree*)file->Get("wcpselection/T_pot");
  TTree *T_PFeval = (TTree*)file->Get("wcpselection/T_PFeval");
  TTree *T_KINEvars = (TTree*)file->Get("wcpselection/T_KINEvars");

  EvalInfo eval;
  POTInfo pot;
  TaggerInfo tagger;
  PFevalInfo pfeval;
  KineInfo kine;

  kine.kine_energy_particle = new std::vector<float>;
  kine.kine_energy_info = new std::vector<int>;
  kine.kine_particle_type = new std::vector<int>;
  kine.kine_energy_included = new std::vector<int>;
  
  
  tagger.pio_2_v_dis2 = new std::vector<float>;
  tagger.pio_2_v_angle2 = new std::vector<float>;
  tagger.pio_2_v_acc_length = new std::vector<float>;
  tagger.pio_2_v_flag = new std::vector<float>;
  tagger.sig_1_v_angle = new std::vector<float>;
  tagger.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger.sig_1_v_energy = new std::vector<float>;
  tagger.sig_1_v_energy_1 = new std::vector<float>;
  tagger.sig_1_v_flag = new std::vector<float>;
  tagger.sig_2_v_energy = new std::vector<float>;
  tagger.sig_2_v_shower_angle = new std::vector<float>;
  tagger.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_flag = new std::vector<float>;
  tagger.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_energy = new std::vector<float>;
  tagger.stw_2_v_angle = new std::vector<float>;
  tagger.stw_2_v_dir_length = new std::vector<float>;
  tagger.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_flag = new std::vector<float>;
  tagger.stw_3_v_angle = new std::vector<float>;
  tagger.stw_3_v_dir_length = new std::vector<float>;
  tagger.stw_3_v_energy = new std::vector<float>;
  tagger.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_3_v_flag = new std::vector<float>;
  tagger.stw_4_v_angle = new std::vector<float>;
  tagger.stw_4_v_dis = new std::vector<float>;
  tagger.stw_4_v_energy = new std::vector<float>;
  tagger.stw_4_v_flag = new std::vector<float>;
  tagger.br3_3_v_energy = new std::vector<float>;
  tagger.br3_3_v_angle = new std::vector<float>;
  tagger.br3_3_v_dir_length = new std::vector<float>;
  tagger.br3_3_v_length = new std::vector<float>;
  tagger.br3_3_v_flag = new std::vector<float>;
  tagger.br3_5_v_dir_length = new std::vector<float>;
  tagger.br3_5_v_total_length = new std::vector<float>;
  tagger.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger.br3_5_v_n_seg = new std::vector<float>;
  tagger.br3_5_v_angle = new std::vector<float>;
  tagger.br3_5_v_sg_length = new std::vector<float>;
  tagger.br3_5_v_energy = new std::vector<float>;
  tagger.br3_5_v_n_main_segs = new std::vector<float>;
  tagger.br3_5_v_n_segs = new std::vector<float>;
  tagger.br3_5_v_shower_main_length = new std::vector<float>;
  tagger.br3_5_v_shower_total_length = new std::vector<float>;
  tagger.br3_5_v_flag = new std::vector<float>;
  tagger.br3_6_v_angle = new std::vector<float>;
  tagger.br3_6_v_angle1 = new std::vector<float>;
  tagger.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger.br3_6_v_direct_length = new std::vector<float>;
  tagger.br3_6_v_length = new std::vector<float>;
  tagger.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger.br3_6_v_energy = new std::vector<float>;
  tagger.br3_6_v_flag = new std::vector<float>;
  tagger.tro_1_v_particle_type = new std::vector<float>;
  tagger.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger.tro_1_v_min_dis = new std::vector<float>;
  tagger.tro_1_v_sg1_length = new std::vector<float>;
  tagger.tro_1_v_shower_main_length = new std::vector<float>;
  tagger.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger.tro_1_v_tmp_length = new std::vector<float>;
  tagger.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger.tro_1_v_flag = new std::vector<float>;
  tagger.tro_2_v_energy = new std::vector<float>;
  tagger.tro_2_v_stem_length = new std::vector<float>;
  tagger.tro_2_v_iso_angle = new std::vector<float>;
  tagger.tro_2_v_max_length = new std::vector<float>;
  tagger.tro_2_v_angle = new std::vector<float>;
  tagger.tro_2_v_flag = new std::vector<float>;
  tagger.tro_4_v_dir2_mag = new std::vector<float>;
  tagger.tro_4_v_angle = new std::vector<float>;
  tagger.tro_4_v_angle1 = new std::vector<float>;
  tagger.tro_4_v_angle2 = new std::vector<float>;
  tagger.tro_4_v_length = new std::vector<float>;
  tagger.tro_4_v_length1 = new std::vector<float>;
  tagger.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_energy = new std::vector<float>;
  tagger.tro_4_v_shower_main_length = new std::vector<float>;
  tagger.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger.tro_4_v_flag = new std::vector<float>;
  tagger.tro_5_v_max_angle = new std::vector<float>;
  tagger.tro_5_v_min_angle = new std::vector<float>;
  tagger.tro_5_v_max_length = new std::vector<float>;
  tagger.tro_5_v_iso_angle = new std::vector<float>;
  tagger.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger.tro_5_v_min_count = new std::vector<float>;
  tagger.tro_5_v_max_count = new std::vector<float>;
  tagger.tro_5_v_energy = new std::vector<float>;
  tagger.tro_5_v_flag = new std::vector<float>;
  tagger.lol_1_v_energy = new std::vector<float>;
  tagger.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_1_v_nseg = new std::vector<float>;
  tagger.lol_1_v_angle = new std::vector<float>;
  tagger.lol_1_v_flag = new std::vector<float>;
  tagger.lol_2_v_length = new std::vector<float>;
  tagger.lol_2_v_angle = new std::vector<float>;
  tagger.lol_2_v_type = new std::vector<float>;
  tagger.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_2_v_energy = new std::vector<float>;
  tagger.lol_2_v_shower_main_length = new std::vector<float>;
  tagger.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger.lol_2_v_flag = new std::vector<float>;
  tagger.cosmict_flag_10 = new std::vector<float>;
  tagger.cosmict_10_flag_inside = new std::vector<float>;
  tagger.cosmict_10_vtx_z = new std::vector<float>;
  tagger.cosmict_10_flag_shower = new std::vector<float>;
  tagger.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger.cosmict_10_angle_beam = new std::vector<float>;
  tagger.cosmict_10_length = new std::vector<float>;
  tagger.numu_cc_flag_1 = new std::vector<float>;
  tagger.numu_cc_1_particle_type = new std::vector<float>;
  tagger.numu_cc_1_length = new std::vector<float>;
  tagger.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger.numu_cc_1_direct_length = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger.numu_cc_flag_2 = new std::vector<float>;
  tagger.numu_cc_2_length = new std::vector<float>;
  tagger.numu_cc_2_total_length = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_all = new std::vector<float>;
  tagger.pio_2_v_dis2 = new std::vector<float>;
  tagger.pio_2_v_angle2 = new std::vector<float>;
  tagger.pio_2_v_acc_length = new std::vector<float>;
  tagger.pio_2_v_flag = new std::vector<float>;
  tagger.sig_1_v_angle = new std::vector<float>;
  tagger.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger.sig_1_v_energy = new std::vector<float>;
  tagger.sig_1_v_energy_1 = new std::vector<float>;
  tagger.sig_1_v_flag = new std::vector<float>;
  tagger.sig_2_v_energy = new std::vector<float>;
  tagger.sig_2_v_shower_angle = new std::vector<float>;
  tagger.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_flag = new std::vector<float>;
  tagger.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_energy = new std::vector<float>;
  tagger.stw_2_v_angle = new std::vector<float>;
  tagger.stw_2_v_dir_length = new std::vector<float>;
  tagger.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_flag = new std::vector<float>;
  tagger.stw_3_v_angle = new std::vector<float>;
  tagger.stw_3_v_dir_length = new std::vector<float>;
  tagger.stw_3_v_energy = new std::vector<float>;
  tagger.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_3_v_flag = new std::vector<float>;
  tagger.stw_4_v_angle = new std::vector<float>;
  tagger.stw_4_v_dis = new std::vector<float>;
  tagger.stw_4_v_energy = new std::vector<float>;
  tagger.stw_4_v_flag = new std::vector<float>;
  tagger.br3_3_v_energy = new std::vector<float>;
  tagger.br3_3_v_angle = new std::vector<float>;
  tagger.br3_3_v_dir_length = new std::vector<float>;
  tagger.br3_3_v_length = new std::vector<float>;
  tagger.br3_3_v_flag = new std::vector<float>;
  tagger.br3_5_v_dir_length = new std::vector<float>;
  tagger.br3_5_v_total_length = new std::vector<float>;
  tagger.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger.br3_5_v_n_seg = new std::vector<float>;
  tagger.br3_5_v_angle = new std::vector<float>;
  tagger.br3_5_v_sg_length = new std::vector<float>;
  tagger.br3_5_v_energy = new std::vector<float>;
  tagger.br3_5_v_n_main_segs = new std::vector<float>;
  tagger.br3_5_v_n_segs = new std::vector<float>;
  tagger.br3_5_v_shower_main_length = new std::vector<float>;
  tagger.br3_5_v_shower_total_length = new std::vector<float>;
  tagger.br3_5_v_flag = new std::vector<float>;
  tagger.br3_6_v_angle = new std::vector<float>;
  tagger.br3_6_v_angle1 = new std::vector<float>;
  tagger.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger.br3_6_v_direct_length = new std::vector<float>;
  tagger.br3_6_v_length = new std::vector<float>;
  tagger.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger.br3_6_v_energy = new std::vector<float>;
  tagger.br3_6_v_flag = new std::vector<float>;
  tagger.tro_1_v_particle_type = new std::vector<float>;
  tagger.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger.tro_1_v_min_dis = new std::vector<float>;
  tagger.tro_1_v_sg1_length = new std::vector<float>;
  tagger.tro_1_v_shower_main_length = new std::vector<float>;
  tagger.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger.tro_1_v_tmp_length = new std::vector<float>;
  tagger.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger.tro_1_v_flag = new std::vector<float>;
  tagger.tro_2_v_energy = new std::vector<float>;
  tagger.tro_2_v_stem_length = new std::vector<float>;
  tagger.tro_2_v_iso_angle = new std::vector<float>;
  tagger.tro_2_v_max_length = new std::vector<float>;
  tagger.tro_2_v_angle = new std::vector<float>;
  tagger.tro_2_v_flag = new std::vector<float>;
  tagger.tro_4_v_dir2_mag = new std::vector<float>;
  tagger.tro_4_v_angle = new std::vector<float>;
  tagger.tro_4_v_angle1 = new std::vector<float>;
  tagger.tro_4_v_angle2 = new std::vector<float>;
  tagger.tro_4_v_length = new std::vector<float>;
  tagger.tro_4_v_length1 = new std::vector<float>;
  tagger.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_energy = new std::vector<float>;
  tagger.tro_4_v_shower_main_length = new std::vector<float>;
  tagger.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger.tro_4_v_flag = new std::vector<float>;
  tagger.tro_5_v_max_angle = new std::vector<float>;
  tagger.tro_5_v_min_angle = new std::vector<float>;
  tagger.tro_5_v_max_length = new std::vector<float>;
  tagger.tro_5_v_iso_angle = new std::vector<float>;
  tagger.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger.tro_5_v_min_count = new std::vector<float>;
  tagger.tro_5_v_max_count = new std::vector<float>;
  tagger.tro_5_v_energy = new std::vector<float>;
  tagger.tro_5_v_flag = new std::vector<float>;
  tagger.lol_1_v_energy = new std::vector<float>;
  tagger.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_1_v_nseg = new std::vector<float>;
  tagger.lol_1_v_angle = new std::vector<float>;
  tagger.lol_1_v_flag = new std::vector<float>;
  tagger.lol_2_v_length = new std::vector<float>;
  tagger.lol_2_v_angle = new std::vector<float>;
  tagger.lol_2_v_type = new std::vector<float>;
  tagger.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_2_v_energy = new std::vector<float>;
  tagger.lol_2_v_shower_main_length = new std::vector<float>;
  tagger.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger.lol_2_v_flag = new std::vector<float>;
  tagger.cosmict_flag_10 = new std::vector<float>;
  tagger.cosmict_10_flag_inside = new std::vector<float>;
  tagger.cosmict_10_vtx_z = new std::vector<float>;
  tagger.cosmict_10_flag_shower = new std::vector<float>;
  tagger.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger.cosmict_10_angle_beam = new std::vector<float>;
  tagger.cosmict_10_length = new std::vector<float>;
  tagger.numu_cc_flag_1 = new std::vector<float>;
  tagger.numu_cc_1_particle_type = new std::vector<float>;
  tagger.numu_cc_1_length = new std::vector<float>;
  tagger.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger.numu_cc_1_direct_length = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger.numu_cc_flag_2 = new std::vector<float>;
  tagger.numu_cc_2_length = new std::vector<float>;
  tagger.numu_cc_2_total_length = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_all = new std::vector<float>;

  set_tree_address(T_BDTvars, tagger,2 );
  set_tree_address(T_eval, eval);
  set_tree_address(T_PFeval, pfeval);
  set_tree_address(T_pot, pot);
  set_tree_address(T_KINEvars, kine);

  double total_pot = 0;
  for (Int_t i=0;i!=T_pot->GetEntries();i++){
    T_pot->GetEntry(i);
    total_pot += pot.pot_tor875;
  }
  // total POT calculations ...
  map_filename_pot[input_filename] = total_pot;

  // fill histogram ...
  T_BDTvars->SetBranchStatus("*",0);
  T_BDTvars->SetBranchStatus("numu_cc_flag",1);
  T_BDTvars->SetBranchStatus("numu_score",1);
  T_BDTvars->SetBranchStatus("nue_score",1);
  T_BDTvars->SetBranchStatus("cosmict_flag",1);

  T_eval->SetBranchStatus("*",0);
  T_eval->SetBranchStatus("run",1);
  T_eval->SetBranchStatus("subrun",1);
  T_eval->SetBranchStatus("event",1);
  
  T_eval->SetBranchStatus("match_isFC",1);
  T_eval->SetBranchStatus("match_found",1);
  if (T_eval->GetBranch("match_found_asInt")) T_eval->SetBranchStatus("match_found_asInt",1); 
  T_eval->SetBranchStatus("stm_eventtype",1);
  T_eval->SetBranchStatus("stm_lowenergy",1);
  T_eval->SetBranchStatus("stm_LM",1);
  T_eval->SetBranchStatus("stm_TGM",1);
  T_eval->SetBranchStatus("stm_STM",1);
  T_eval->SetBranchStatus("stm_FullDead",1);
  T_eval->SetBranchStatus("stm_clusterlength",1);
  
  
  T_eval->SetBranchStatus("weight_spline",1);
  T_eval->SetBranchStatus("weight_cv",1);
  T_eval->SetBranchStatus("weight_lee",1);
  T_eval->SetBranchStatus("weight_change",1);
  // MC enable truth information ...
  T_eval->SetBranchStatus("truth_isCC",1);
  T_eval->SetBranchStatus("truth_nuPdg",1);
  T_eval->SetBranchStatus("truth_vtxInside",1);
  T_eval->SetBranchStatus("truth_nuEnergy",1);
  T_eval->SetBranchStatus("truth_vtxX",1);
  T_eval->SetBranchStatus("truth_vtxY",1);
  T_eval->SetBranchStatus("truth_vtxZ",1);
 
  
  T_KINEvars->SetBranchStatus("*",0);
  T_KINEvars->SetBranchStatus("kine_reco_Enu",1);
  T_KINEvars->SetBranchStatus("kine_pio_mass",1);
  T_KINEvars->SetBranchStatus("kine_pio_flag",1);
  T_KINEvars->SetBranchStatus("kine_pio_vtx_dis",1);
  T_KINEvars->SetBranchStatus("kine_pio_energy_1",1);
  T_KINEvars->SetBranchStatus("kine_pio_theta_1",1);
  T_KINEvars->SetBranchStatus("kine_pio_phi_1",1);
  T_KINEvars->SetBranchStatus("kine_pio_dis_1",1);
  T_KINEvars->SetBranchStatus("kine_pio_energy_2",1);
  T_KINEvars->SetBranchStatus("kine_pio_theta_2",1);
  T_KINEvars->SetBranchStatus("kine_pio_phi_2",1);
  T_KINEvars->SetBranchStatus("kine_pio_dis_2",1);
  T_KINEvars->SetBranchStatus("kine_pio_angle",1);

  T_PFeval->SetBranchStatus("*",0);


  WeightInfo weight;
  TTree *T_weight = (TTree*)file->Get("wcpselection/T_weight");
  weight.file_type = new std::string();
  weight.expskin_FluxUnisim= new std::vector<float>;
  weight.horncurrent_FluxUnisim= new std::vector<float>;
  weight.kminus_PrimaryHadronNormalization= new std::vector<float>;
  weight.kplus_PrimaryHadronFeynmanScaling= new std::vector<float>;
  weight.kzero_PrimaryHadronSanfordWang= new std::vector<float>;
  weight.nucleoninexsec_FluxUnisim= new std::vector<float>;
  weight.nucleonqexsec_FluxUnisim= new std::vector<float>;
  weight.nucleontotxsec_FluxUnisim= new std::vector<float>;
  weight.piminus_PrimaryHadronSWCentralSplineVariation= new std::vector<float>;
  weight.pioninexsec_FluxUnisim= new std::vector<float>;
  weight.pionqexsec_FluxUnisim= new std::vector<float>;
  weight.piontotxsec_FluxUnisim= new std::vector<float>;
  weight.piplus_PrimaryHadronSWCentralSplineVariation= new std::vector<float>;
  
  weight.All_UBGenie= new std::vector<float>;
  weight.AxFFCCQEshape_UBGenie= new std::vector<float>;
  weight.DecayAngMEC_UBGenie= new std::vector<float>;
  weight.NormCCCOH_UBGenie= new std::vector<float>;
  weight.NormNCCOH_UBGenie= new std::vector<float>;
  weight.RPA_CCQE_Reduced_UBGenie= new std::vector<float>;
  weight.RPA_CCQE_UBGenie= new std::vector<float>;
  weight.RootinoFix_UBGenie= new std::vector<float>;
  weight.ThetaDelta2NRad_UBGenie= new std::vector<float>;
  weight.Theta_Delta2Npi_UBGenie= new std::vector<float>;
  weight.TunedCentralValue_UBGenie= new std::vector<float>;
  weight.VecFFCCQEshape_UBGenie= new std::vector<float>;
  weight.XSecShape_CCMEC_UBGenie= new std::vector<float>;
  weight.splines_general_Spline= new std::vector<float>;
  weight.xsr_scc_Fa3_SCC= new std::vector<float>;
  weight.xsr_scc_Fv3_SCC= new std::vector<float>;

  TString option;
  if (T_weight->GetBranch("expskin_FluxUnisim")){
    option = "expskin_FluxUnisim";
  }else if (T_weight->GetBranch("horncurrent_FluxUnisim")){
    option = "horncurrent_FluxUnisim";
  }else if (T_weight->GetBranch("kminus_PrimaryHadronNormalization")){
    option = "kminus_PrimaryHadronNormalization";
  }else if (T_weight->GetBranch("kplus_PrimaryHadronFeynmanScaling")){
    option = "kplus_PrimaryHadronFeynmanScaling";
  }else if (T_weight->GetBranch("kzero_PrimaryHadronSanfordWang")){
    option = "kzero_PrimaryHadronSanfordWang";
  }else if (T_weight->GetBranch("nucleoninexsec_FluxUnisim")){
    option = "nucleoninexsec_FluxUnisim";
  }else if (T_weight->GetBranch("nucleonqexsec_FluxUnisim")){
    option = "nucleonqexsec_FluxUnisim";
  }else if (T_weight->GetBranch("nucleontotxsec_FluxUnisim")){
    option = "nucleontotxsec_FluxUnisim";
  }else if (T_weight->GetBranch("piminus_PrimaryHadronSWCentralSplineVariation")){
    option = "piminus_PrimaryHadronSWCentralSplineVariation";
  }else if (T_weight->GetBranch("pioninexsec_FluxUnisim")){
    option = "pioninexsec_FluxUnisim";
  }else if (T_weight->GetBranch("pionqexsec_FluxUnisim")){
    option = "pionqexsec_FluxUnisim";
  }else if (T_weight->GetBranch("piontotxsec_FluxUnisim")){
    option = "piontotxsec_FluxUnisim";
  }else if (T_weight->GetBranch("piplus_PrimaryHadronSWCentralSplineVariation")){
    option = "piplus_PrimaryHadronSWCentralSplineVariation";
  }else if (T_weight->GetBranch("All_UBGenie")){
    option = "UBGenieFluxSmallUni";
  }
  
  set_tree_address(T_weight, weight, option);
  //std::cout << T_eval->GetEntries() << " " << T_weight->GetEntries() << " " << option << std::endl;

  std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = get_histograms(input_filename,0);

  std::set<std::tuple<float, float, std::vector<float>, std::vector<int>, std::set<std::pair<int, float> > > >& set_events = map_passed_events[input_filename];


  for (size_t i=0;i!=T_eval->GetEntries();i++){
    T_BDTvars->GetEntry(i);
    T_eval->GetEntry(i);
    T_KINEvars->GetEntry(i);
    T_PFeval->GetEntry(i);
    T_weight->GetEntry(i);

    std::tuple<float, float, std::vector<float>, std::vector<int>, std::set<std::pair<int, float> > > event_info;
    std::get<0>(event_info) = eval.weight_cv * eval.weight_spline;
    std::get<1>(event_info) = leeweight(eval.truth_nuEnergy);
    
     for (auto it = histo_infos.begin(); it != histo_infos.end(); it++){
      TString histoname = std::get<0>(*it);

      auto it2 = map_histoname_infos.find(histoname);
      int no = std::get<0>(it2->second);
      
      TString var_name = std::get<4>(*it);
      TString ch_name = std::get<5>(*it);
      TString add_cut = std::get<6>(*it);

      float val = get_kine_var(kine, pfeval, var_name);
      bool flag_pass = get_cut_pass(ch_name, add_cut, false, eval, tagger, kine);

      if (flag_pass) std::get<4>(event_info).insert(std::make_pair(no, val));
     }
    
    if (std::get<4>(event_info).size()>0){
      if (option == "expskin_FluxUnisim"){
	std::get<2>(event_info).resize(weight.expskin_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.expskin_FluxUnisim->size());
	for (size_t j=0;j!=weight.expskin_FluxUnisim->size();j++){
	  std::get<2>(event_info).at(j) = weight.expskin_FluxUnisim->at(j) - 1.0; // relative ...
	}
      }else if (option == "horncurrent_FluxUnisim"){
	std::get<2>(event_info).resize(weight.horncurrent_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.horncurrent_FluxUnisim->size());
	for (size_t j=0;j!= weight.horncurrent_FluxUnisim->size(); j++){
	  std::get<2>(event_info).at(j) = weight.horncurrent_FluxUnisim->at(j) - 1.0; // relative ...
	}
      }else if (option == "kminus_PrimaryHadronNormalization"){
	std::get<2>(event_info).resize(weight.kminus_PrimaryHadronNormalization->size());
	std::get<3>(event_info).push_back(weight.kminus_PrimaryHadronNormalization->size());
	for (size_t j=0;j!= weight.kminus_PrimaryHadronNormalization->size(); j++){
	  std::get<2>(event_info).at(j) = weight.kminus_PrimaryHadronNormalization->at(j) - 1.0; 
	}
      }else if (option == "kplus_PrimaryHadronFeynmanScaling"){
	std::get<2>(event_info).resize(weight.kplus_PrimaryHadronFeynmanScaling->size());
	std::get<3>(event_info).push_back(weight.kplus_PrimaryHadronFeynmanScaling->size());
	for (size_t j=0;j!=weight.kplus_PrimaryHadronFeynmanScaling->size();j++){
	  std::get<2>(event_info).at(j) = weight.kplus_PrimaryHadronFeynmanScaling->at(j) - 1.0;
	}
      }else if (option == "kzero_PrimaryHadronSanfordWang"){
	std::get<2>(event_info).resize(weight.kzero_PrimaryHadronSanfordWang->size());
	std::get<3>(event_info).push_back(weight.kzero_PrimaryHadronSanfordWang->size());
	for (size_t j=0;j!=weight.kzero_PrimaryHadronSanfordWang->size();j++){
	  std::get<2>(event_info).at(j) = weight.kzero_PrimaryHadronSanfordWang->at(j) - 1.0;
	}
      }else if (option == "nucleoninexsec_FluxUnisim"){
	std::get<2>(event_info).resize(weight.nucleoninexsec_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.nucleoninexsec_FluxUnisim->size());
	for (size_t j=0;j!=weight.nucleoninexsec_FluxUnisim->size();j++){
	  std::get<2>(event_info).at(j) = weight.nucleoninexsec_FluxUnisim->at(j) - 1.0;
	}
      }else if (option == "nucleonqexsec_FluxUnisim"){
	std::get<2>(event_info).resize(weight.nucleonqexsec_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.nucleonqexsec_FluxUnisim->size());
	for (size_t j=0; j!= weight.nucleonqexsec_FluxUnisim->size(); j++){
	  std::get<2>(event_info).at(j) = weight.nucleonqexsec_FluxUnisim->at(j) - 1.0;
	}
      }else if (option == "nucleontotxsec_FluxUnisim"){
	std::get<2>(event_info).resize(weight.nucleontotxsec_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.nucleontotxsec_FluxUnisim->size());
	for (size_t j=0; j!= weight.nucleontotxsec_FluxUnisim->size(); j++){
	  std::get<2>(event_info).at(j) = weight.nucleontotxsec_FluxUnisim->at(j) - 1.0;
	}
      }else if (option == "piminus_PrimaryHadronSWCentralSplineVariation"){
	std::get<2>(event_info).resize(weight.piminus_PrimaryHadronSWCentralSplineVariation->size());
	std::get<3>(event_info).push_back(weight.piminus_PrimaryHadronSWCentralSplineVariation->size());
	for (size_t j=0; j!= weight.piminus_PrimaryHadronSWCentralSplineVariation->size(); j++){
	  std::get<2>(event_info).at(j) = weight.piminus_PrimaryHadronSWCentralSplineVariation->at(j) - 1.0;
	}
      }else if (option == "pioninexsec_FluxUnisim"){
	std::get<2>(event_info).resize(weight.pioninexsec_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.pioninexsec_FluxUnisim->size());
	for (size_t j=0;j!=weight.pioninexsec_FluxUnisim->size();j++){
	  std::get<2>(event_info).at(j) = weight.pioninexsec_FluxUnisim->at(j) - 1.0;
	}
      }else if (option == "pionqexsec_FluxUnisim"){
	std::get<2>(event_info).resize(weight.pionqexsec_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.pionqexsec_FluxUnisim->size());
	for (size_t j=0;j!=weight.pionqexsec_FluxUnisim->size();j++){
	  std::get<2>(event_info).at(j) = weight.pionqexsec_FluxUnisim->at(j) - 1.0;
	}
      }else if (option == "piontotxsec_FluxUnisim"){
	std::get<2>(event_info).resize(weight.piontotxsec_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.piontotxsec_FluxUnisim->size());
	for (size_t j=0;j!=weight.piontotxsec_FluxUnisim->size();j++){
	  std::get<2>(event_info).at(j) = weight.piontotxsec_FluxUnisim->at(j) - 1.0;
	}
      }else if (option == "piplus_PrimaryHadronSWCentralSplineVariation"){
	std::get<2>(event_info).resize(weight.piplus_PrimaryHadronSWCentralSplineVariation->size());
	std::get<3>(event_info).push_back(weight.piplus_PrimaryHadronSWCentralSplineVariation->size());
	for (size_t j=0; j!= weight.piplus_PrimaryHadronSWCentralSplineVariation->size(); j++){
	  std::get<2>(event_info).at(j) = weight.piplus_PrimaryHadronSWCentralSplineVariation->at(j) - 1.0;
	}
      }else if (option == "UBGenieFluxSmallUni"){
	
      }
      set_events.insert(event_info);
    }
    
  }
  
  delete file;
}
