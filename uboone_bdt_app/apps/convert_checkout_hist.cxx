#include <iostream>

#include "WCPLEEANA/master_cov_matrix.h"

#include "TROOT.h"
#include "TMath.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"

#include "WCPLEEANA/cuts.h"
#include "WCPLEEANA/pot.h"
#include "WCPLEEANA/pfeval.h"

using namespace std;
using namespace LEEana;

int main( int argc, char** argv )
{
  CovMatrix cov;

  TString input_filename = argv[1];
  TString out_filename = argv[2];


  bool flag_data = true;
  TFile *file = new TFile(input_filename);

  TTree *T_BDTvars = (TTree*)file->Get("wcpselection/T_BDTvars");
  TTree *T_eval = (TTree*)file->Get("wcpselection/T_eval");
  TTree *T_pot = (TTree*)file->Get("wcpselection/T_pot");
  TTree *T_PFeval = (TTree*)file->Get("wcpselection/T_PFeval");
  TTree *T_KINEvars = (TTree*)file->Get("wcpselection/T_KINEvars");

  if (T_eval->GetBranch("weight_cv")) flag_data = false;

  EvalInfo eval;
  POTInfo pot;
  TaggerInfo tagger;
  PFevalInfo pfeval;
  KineInfo kine;

#include "init.txt"

  set_tree_address(T_BDTvars, tagger,2 );
  if (flag_data){
    set_tree_address(T_eval, eval,2);
    set_tree_address(T_PFeval, pfeval,2);
  }else{
    set_tree_address(T_eval, eval);
    set_tree_address(T_PFeval, pfeval);
  }
  set_tree_address(T_pot, pot);
  set_tree_address(T_KINEvars, kine);

  double total_pot = 0;
  for (Int_t i=0;i!=T_pot->GetEntries();i++){
    T_pot->GetEntry(i);
    total_pot += pot.pot_tor875;
  }
  double ext_pot = cov.get_ext_pot(input_filename);
  if (ext_pot != 0) total_pot = ext_pot;
  
  std::cout << "Total POT: " << total_pot << " external POT: " << ext_pot << std::endl;

  // prepare histograms ...
  // declare histograms ...
  TH1F *htemp;
  std::map<TString, TH1F*> map_histoname_hist;
  std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > all_histo_infos;
  
  std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = cov.get_histograms(input_filename,0);
  std::copy(histo_infos.begin(), histo_infos.end(), std::back_inserter(all_histo_infos));
  //  std::cout << "CV:" << std::endl;
  for (auto it = histo_infos.begin(); it != histo_infos.end(); it++){
    TString histoname = std::get<0>(*it);
    Int_t nbin = std::get<1>(*it);
    float llimit = std::get<2>(*it);
    float hlimit = std::get<3>(*it);
    TString var_name = std::get<4>(*it);
    TString ch_name = std::get<5>(*it);
    TString add_cut = std::get<6>(*it);
    TString weight = std::get<7>(*it);
    
    //    std::cout << std::get<0>( *it)  << " " << std::get<1>(*it) << " " << std::get<4>(*it) << " " << std::get<5>(*it) << " " << std::get<6>(*it) << " " << std::get<7>(*it) << std::endl;
    htemp = new TH1F(histoname, histoname, nbin, llimit, hlimit);
    map_histoname_hist[histoname] = htemp;
  }
  //  std::cout << std::endl;

  std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos_err2 = cov.get_histograms(input_filename,1);
  std::copy(histo_infos_err2.begin(), histo_infos_err2.end(), std::back_inserter(all_histo_infos));
  //  std::cout << "Error2: " << std::endl;
  for (auto it = histo_infos_err2.begin(); it != histo_infos_err2.end(); it++){
    TString histoname = std::get<0>(*it);
    Int_t nbin = std::get<1>(*it);
    float llimit = std::get<2>(*it);
    float hlimit = std::get<3>(*it);
    TString var_name = std::get<4>(*it);
    TString ch_name = std::get<5>(*it);
    TString add_cut = std::get<6>(*it);
    TString weight = std::get<7>(*it);
    //std::cout << std::get<0>( *it) << " " << std::get<1>(*it) << " " << std::get<4>(*it) << " " << std::get<5>(*it) << " " << std::get<6>(*it) << " " << std::get<7>(*it) << std::endl;
    htemp = new TH1F(histoname, histoname, nbin, llimit, hlimit);
    map_histoname_hist[histoname] = htemp;
  }
  //  std::cout << std::endl;

  std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos_cros = cov.get_histograms(input_filename,2);
  std::copy(histo_infos_cros.begin(), histo_infos_cros.end(), std::back_inserter(all_histo_infos));
  //std::cout << "Cross: " << std::endl;
  for (auto it = histo_infos_cros.begin(); it != histo_infos_cros.end(); it++){
    TString histoname = std::get<0>(*it);
    Int_t nbin = std::get<1>(*it);
    float llimit = std::get<2>(*it);
    float hlimit = std::get<3>(*it);
    TString var_name = std::get<4>(*it);
    TString ch_name = std::get<5>(*it);
    TString add_cut = std::get<6>(*it);
    TString weight = std::get<7>(*it);
    //std::cout << std::get<0>( *it) << " " << std::get<1>(*it) << " " << std::get<4>(*it) << " " << std::get<5>(*it) << " " << std::get<6>(*it) << " " << std::get<7>(*it) << std::endl;
    htemp = new TH1F(histoname, histoname, nbin, llimit, hlimit);
    map_histoname_hist[histoname] = htemp;
  }
  //  std::cout << std::endl;
  
  // fill histogram ...
  T_BDTvars->SetBranchStatus("*",0);
  T_BDTvars->SetBranchStatus("numu_cc_flag",1);
  T_BDTvars->SetBranchStatus("numu_score",1);
  T_BDTvars->SetBranchStatus("nue_score",1);
  T_BDTvars->SetBranchStatus("cosmict_flag",1);

  T_eval->SetBranchStatus("*",0);
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
  
  if (!flag_data){
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
  }

  
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
  
  std::cout << "Total entries: " << T_eval->GetEntries() << std::endl;
  for (Int_t i=0;i!=T_eval->GetEntries();i++){
    T_BDTvars->GetEntry(i);
    T_eval->GetEntry(i);
    T_KINEvars->GetEntry(i);
    T_PFeval->GetEntry(i);
    
    if (!is_preselection(eval)) continue;

    for (auto it = all_histo_infos.begin(); it != all_histo_infos.end(); it++){
      TString histoname = std::get<0>(*it);
      Int_t nbin = std::get<1>(*it);
      float llimit = std::get<2>(*it);
      float hlimit = std::get<3>(*it);
      TString var_name = std::get<4>(*it);
      TString ch_name = std::get<5>(*it);
      TString add_cut = std::get<6>(*it);
      TString weight = std::get<7>(*it);
 
      htemp = map_histoname_hist[histoname];
      // get kinematics variable ...
      double val = get_kine_var(kine, pfeval, var_name);
      // get pass or not
      bool flag_pass = get_cut_pass(ch_name, add_cut, flag_data, eval, tagger, kine);

      // std::cout << weight << std::endl;
      // get weight ...
      double weight_val = get_weight(weight, eval);

      if (flag_pass)
	htemp->Fill(val,weight_val);
    }
  }
  

  // save histograms ...
  TFile *file1 = new TFile(out_filename,"RECREATE");
  file1->cd();
  TTree *T = new TTree("T","T");
  T->SetDirectory(file1);
  T->Branch("pot",&total_pot,"pot/D");
  T->Fill();
  
  for (auto it = map_histoname_hist.begin(); it!= map_histoname_hist.end(); it++){
    it->second->SetDirectory(file1);
  }
  
  file1->Write();
  file1->Close();
  
  return 0;
}
