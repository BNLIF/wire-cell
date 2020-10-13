// cz: code modified from tutorials/tmva/TMVAClassification.C

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"

#include "WCPLEEANA/tagger.h"

#include "WCPLEEANA/eval.h"

using namespace std;
using namespace LEEana;

//#include "WCPLEEANA/bdt.h"
#include "WCPLEEANA/pot.h"
#include "WCPLEEANA/pfeval.h"
#include "WCPLEEANA/kine.h"
#include "WCPLEEANA/weights.h"

int main( int argc, char** argv )
{
  if (argc < 5) {
    std::cout << "merge_xf #input_file_cv #input_file_xf #output_file #option" << std::endl;
    return -1;
  }

  TString input_file_cv = argv[1];
  TString input_file_xf = argv[2];
  TString out_file = argv[3];
  TString option = argv[4];

  TFile *file1 = new TFile(input_file_cv);
  TTree *T_BDTvars = (TTree*)file1->Get("wcpselection/T_BDTvars");
  TTree *T_eval = (TTree*)file1->Get("wcpselection/T_eval");
  TTree *T_pot = (TTree*)file1->Get("wcpselection/T_pot");
  TTree *T_PFeval = (TTree*)file1->Get("wcpselection/T_PFeval");
  TTree *T_KINEvars = (TTree*)file1->Get("wcpselection/T_KINEvars");

  TFile *file2 = new TFile(input_file_xf);
  TTree *T;
  WeightInfo weight;
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
  
 

  
  
  
  if (option == "expskin_FluxUnisim"){
    T = (TTree*)file2->Get("expskin_FluxUnisim");
  }else if (option == "horncurrent_FluxUnisim"){
    T = (TTree*)file2->Get("horncurrent_FluxUnisim");
  }else if (option == "kminus_PrimaryHadronNormalization"){
    T = (TTree*)file2->Get("kminus_PrimaryHadronNormalization");
  }else if (option == "kplus_PrimaryHadronFeynmanScaling"){
    T = (TTree*)file2->Get("kplus_PrimaryHadronFeynmanScaling");
  }else if (option == "kzero_PrimaryHadronSanfordWang"){
    T = (TTree*)file2->Get("kzero_PrimaryHadronSanfordWang");
  }else if (option == "nucleoninexsec_FluxUnisim"){
    T = (TTree*)file2->Get("nucleoninexsec_FluxUnisim");
  }else if (option == "nucleonqexsec_FluxUnisim"){
    T = (TTree*)file2->Get("nucleonqexsec_FluxUnisim");
  }else if (option == "nucleontotxsec_FluxUnisim"){
    T = (TTree*)file2->Get("nucleontotxsec_FluxUnisim");
  }else if (option == "piminus_PrimaryHadronSWCentralSplineVariation"){
    T = (TTree*)file2->Get("piminus_PrimaryHadronSWCentralSplineVariation");
  }else if (option == "pioninexsec_FluxUnisim"){
    T = (TTree*)file2->Get("pioninexsec_FluxUnisim");
  }else if (option == "pionqexsec_FluxUnisim"){
    T = (TTree*)file2->Get("pionqexsec_FluxUnisim");
  }else if (option == "piontotxsec_FluxUnisim"){
    T = (TTree*)file2->Get("piontotxsec_FluxUnisim");
  }else if (option == "piplus_PrimaryHadronSWCentralSplineVariation"){
    T = (TTree*)file2->Get("piplus_PrimaryHadronSWCentralSplineVariation");
  }else if (option == "UBGenieFluxSmallUni"){
    T = (TTree*)file2->Get("UBGenieFluxSmallUni");
  }
  
  TFile *file3 = new TFile(out_file,"RECREATE");
  file3->mkdir("wcpselection");
  file3->cd("wcpselection");
  TTree *t1 = new TTree("T_eval","T_eval");
  TTree *t2 = new TTree("T_pot","T_pot");
  TTree *t3 = new TTree("T_PFeval", "T_PFeval");
  TTree *t5 = new TTree("T_KINEvars", "T_KINEvars");
  TTree *t4 = new TTree("T_BDTvars","T_BDTvars");

  TTree *t6 = new TTree("T_weight","T_weight");
  
  EvalInfo eval;
  eval.file_type = new std::string();
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
  set_tree_address(T, weight, option);
  
  put_tree_address(t4, tagger,2);
  put_tree_address(t1, eval);
  put_tree_address(t3, pfeval);
  put_tree_address(t2, pot);
  put_tree_address(t5, kine);
  put_tree_address(t6, weight, option);


  //std::cout << T->GetEntries() << std::endl;


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
  T_eval->SetBranchStatus("run",1);
  T_eval->SetBranchStatus("subrun",1);
  T_eval->SetBranchStatus("event",1);

  
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

  std::map<std::pair<int, int>, int> map_re_entry;
  std::map<std::pair<int, int>, std::set<std::pair<int, int> > > map_rs_re;

  bool flag_presel = false;

  int num_check = 0;

  for (int i=0;i!=T_eval->GetEntries();i++){
    T_eval->GetEntry(i);
    T_BDTvars->GetEntry(i);
    
    map_rs_re[std::make_pair(eval.run, eval.subrun)].insert(std::make_pair(eval.run, eval.event));
    
    int tmp_match_found = eval.match_found;
    if (eval.is_match_found_int) tmp_match_found = eval.match_found_asInt;
    
    flag_presel = false;
    if (tmp_match_found != 0 && eval.stm_eventtype != 0 && eval.stm_lowenergy ==0 && eval.stm_LM ==0 && eval.stm_TGM ==0 && eval.stm_STM==0 && eval.stm_FullDead == 0 && eval.stm_clusterlength >0) {
      flag_presel = true; // preselection ...
    }
    
    if (tmp_match_found == -1  || (tmp_match_found == 1 && eval.stm_lowenergy == -1) || (flag_presel && tagger.numu_cc_flag == -1)) {
      num_check ++;
      continue;
    }

    //if (eval.run == 7486 && eval.event==3964) std::cout << flag_presel << " " << tagger.numu_cc_flag << std::endl;
    
    map_re_entry[std::make_pair(eval.run, eval.event)] = i;
  }

  std::map<std::pair<int, int>, std::pair<int, double> > map_rs_entry_pot;
  for (Int_t i=0;i!=T_pot->GetEntries();i++){
    T_pot->GetEntry(i);
    map_rs_entry_pot[std::make_pair(pot.runNo,pot.subRunNo)] = std::make_pair(i, pot.pot_tor875);
  }

  std::map<std::pair<int,int>, int> map_rs_failed;
  for (auto it =  map_rs_re.begin(); it !=  map_rs_re.end(); it++){
    int failed_num = 0;
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      if (map_re_entry.find(*it1) == map_re_entry.end()) failed_num++;
    }
    map_rs_failed[it->first] = failed_num;
  }

  std::map<std::pair<int, int>, int> map_re_weight_entry;
  for (size_t i=0;i!=T->GetEntries();i++){
    T->GetEntry(i);
    if (get_size(weight, option)>0){
      map_re_weight_entry[std::make_pair(weight.run, weight.event)] = i;
    }
  }

  for (auto it = map_re_entry.begin(); it != map_re_entry.end(); it++){
    auto it1 = map_re_weight_entry.find(it->first);
    if (it1 != map_re_weight_entry.end()){
      T_BDTvars->GetEntry(it->second);
      T_eval->GetEntry(it->second);
      T_KINEvars->GetEntry(it->second);
      T_PFeval->GetEntry(it->second);
      T->GetEntry(it1->second);
    

      t1->Fill();
      t3->Fill();
      t4->Fill();
      t5->Fill();
      t6->Fill();
    }
  }

  double cv_pot = 0;
  double cv1_pot = 0;
  float pass_ratio;
  
  for (auto it = map_rs_entry_pot.begin(); it != map_rs_entry_pot.end(); it++){
    T_pot->GetEntry(it->second.first);
    cv_pot += it->second.second;

    if (map_rs_re[it->first].size()==0) continue;
    
    pass_ratio = 1-map_rs_failed[it->first] * 1.0 / map_rs_re[it->first].size();
    
    cv1_pot += it->second.second * pass_ratio;

    pot.pot_tor875 *= pass_ratio;
    pot.pot_tor875good *= pass_ratio;
    
    t2->Fill();
  }
  std::cout << out_file << std::endl;
  std::cout << "Events: " << t1->GetEntries()<<"/"<<T_eval->GetEntries() << std::endl;
  std::cout << "POT:    " << cv1_pot << " " << cv_pot << std::endl;
  
  file3->Write();
  file3->Close();
  
  return 0;
}
