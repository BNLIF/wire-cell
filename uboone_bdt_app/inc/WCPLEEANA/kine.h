#ifndef UBOONE_LEE_KINE
#define UBOONE_LEE_KINE

namespace LEEana{
struct KineInfo{
  Float_t kine_reco_Enu;
  Float_t kine_reco_add_energy;
  std::vector<float> *kine_energy_particle;
  std::vector<int> *kine_energy_info; 
  std::vector<int> *kine_particle_type;
  std::vector<int> *kine_energy_included;
  Float_t kine_pio_mass;
  Int_t kine_pio_flag;
  Float_t kine_pio_vtx_dis;
  Float_t kine_pio_energy_1;
  Float_t kine_pio_theta_1;
  Float_t kine_pio_phi_1;
  Float_t kine_pio_dis_1;
  Float_t kine_pio_energy_2;
  Float_t kine_pio_theta_2;
  Float_t kine_pio_phi_2;
  Float_t kine_pio_dis_2;
  Float_t kine_pio_angle;
};

void set_tree_address(TTree *tree0, KineInfo& tagger_info);
void put_tree_address(TTree *Tsig, KineInfo& tagger_info);
}

void LEEana::set_tree_address(TTree *tree0, KineInfo& tagger_info){
  tree0->SetBranchAddress("kine_reco_Enu", &tagger_info.kine_reco_Enu);
  tree0->SetBranchAddress("kine_reco_add_energy", &tagger_info.kine_reco_add_energy);
  tree0->SetBranchAddress("kine_energy_particle", &tagger_info.kine_energy_particle);
  tree0->SetBranchAddress("kine_energy_info", &tagger_info.kine_energy_info); 
  tree0->SetBranchAddress("kine_particle_type", &tagger_info.kine_particle_type);
  tree0->SetBranchAddress("kine_energy_included", &tagger_info.kine_energy_included);
  tree0->SetBranchAddress("kine_pio_mass", &tagger_info.kine_pio_mass);
  tree0->SetBranchAddress("kine_pio_flag", &tagger_info.kine_pio_flag);
  tree0->SetBranchAddress("kine_pio_vtx_dis", &tagger_info.kine_pio_vtx_dis);
  tree0->SetBranchAddress("kine_pio_energy_1", &tagger_info.kine_pio_energy_1);
  tree0->SetBranchAddress("kine_pio_theta_1", &tagger_info.kine_pio_theta_1);
  tree0->SetBranchAddress("kine_pio_phi_1", &tagger_info.kine_pio_phi_1);
  tree0->SetBranchAddress("kine_pio_dis_1", &tagger_info.kine_pio_dis_1);
  tree0->SetBranchAddress("kine_pio_energy_2", &tagger_info.kine_pio_energy_2);
  tree0->SetBranchAddress("kine_pio_theta_2", &tagger_info.kine_pio_theta_2);
  tree0->SetBranchAddress("kine_pio_phi_2", &tagger_info.kine_pio_phi_2);
  tree0->SetBranchAddress("kine_pio_dis_2", &tagger_info.kine_pio_dis_2);
  tree0->SetBranchAddress("kine_pio_angle", &tagger_info.kine_pio_angle);
};

void LEEana::put_tree_address(TTree *tree0, KineInfo& tagger_info){
  tree0->Branch("kine_reco_Enu", &tagger_info.kine_reco_Enu,"data/F");
  tree0->Branch("kine_reco_add_energy", &tagger_info.kine_reco_add_energy,"data/F");
  tree0->Branch("kine_energy_particle", &tagger_info.kine_energy_particle);
  tree0->Branch("kine_energy_info", &tagger_info.kine_energy_info); 
  tree0->Branch("kine_particle_type", &tagger_info.kine_particle_type);
  tree0->Branch("kine_energy_included", &tagger_info.kine_energy_included);
  tree0->Branch("kine_pio_mass", &tagger_info.kine_pio_mass,"data/F");
  tree0->Branch("kine_pio_flag", &tagger_info.kine_pio_flag,"data/I");
  tree0->Branch("kine_pio_vtx_dis", &tagger_info.kine_pio_vtx_dis,"data/F");
  tree0->Branch("kine_pio_energy_1", &tagger_info.kine_pio_energy_1,"data/F");
  tree0->Branch("kine_pio_theta_1", &tagger_info.kine_pio_theta_1,"data/F");
  tree0->Branch("kine_pio_phi_1", &tagger_info.kine_pio_phi_1,"data/F");
  tree0->Branch("kine_pio_dis_1", &tagger_info.kine_pio_dis_1,"data/F");
  tree0->Branch("kine_pio_energy_2", &tagger_info.kine_pio_energy_2,"data/F");
  tree0->Branch("kine_pio_theta_2", &tagger_info.kine_pio_theta_2,"data/F");
  tree0->Branch("kine_pio_phi_2", &tagger_info.kine_pio_phi_2,"data/F");
  tree0->Branch("kine_pio_dis_2", &tagger_info.kine_pio_dis_2,"data/F");
  tree0->Branch("kine_pio_angle", &tagger_info.kine_pio_angle,"data/F");
};


#endif
