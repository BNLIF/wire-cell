

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


void LEEana::CovMatrix::fill_pred_histograms(int run, std::map<int, std::vector<TH1F*> >& map_obsch_histos, std::map<int, std::vector< std::vector< std::tuple<double, double, double, int, double> > > >& map_obsch_bayes, std::map<TString, std::pair<TH1F*, double> >& map_name_histogram, float lee_strength, std::map<int, double> map_data_period_pot){
  
  for (auto it = map_pred_obsch_histos.begin(); it!=map_pred_obsch_histos.end();it++){
    //std::cout << it->first << std::endl;
    int obsch = it->first;
    TH1F *hpred = map_obsch_histos[obsch].at(1);
    TH1F *hpred_err2 = map_obsch_histos[obsch].at(2);

    std::map<TString, std::pair<double, double> > temp_map_histo_ratios;

    std::map<TString, std::vector< std::tuple<double, double, double, int, double> > > map_histoname_values;
    
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

	// mean, sigma2, pot_ratio, covch, add_sys2
	std::vector< std::tuple<double, double, double, int, double> > values;
	for (int i=0;i!=hmc->GetNbinsX()+1;i++){
	  values.push_back(std::make_tuple(hmc->GetBinContent(i+1)*ratio , hmc_err2->GetBinContent(i+1)*ratio*ratio, temp_map_data_pot[period]/temp_map_mc_acc_pot[period], map_histogram_covch_add[histoname].first, pow(hmc->GetBinContent(i+1) *ratio * map_histogram_covch_add[histoname].second,2) ));
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
      
      std::vector< std::tuple<double, double, double, int, double> >& values1 = map_histoname_values[hist1];
      std::vector< std::tuple<double, double, double, int, double> >& values2 = map_histoname_values[hist2];
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
