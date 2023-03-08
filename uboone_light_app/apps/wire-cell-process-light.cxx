#include "WCPSst/GeomDataSource.h"
#include "WCP2dToy/ToyLightReco.h"
#include "TH1F.h"
#include "TH2F.h"

#include <iostream>

using namespace WCP;
using namespace std;

int main(int argc, char* argv[])
{
  if (argc < 3) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root " << endl;
    return 1;
  }
  TH1::AddDirectory(kFALSE);

  const char* root_file = argv[2];

  WCP2dToy::ToyLightReco uboone_flash(root_file);

  TFile *file1 = new TFile(root_file);
  TTree *T = (TTree*)file1->Get("/Event/Sim");
  
  TClonesArray* cosmic_hg_wf = new TClonesArray;
  TClonesArray* cosmic_lg_wf = new TClonesArray;
  TClonesArray* beam_hg_wf = new TClonesArray;
  TClonesArray* beam_lg_wf = new TClonesArray;
  vector<short> *cosmic_hg_opch = new vector<short>;
  vector<short> *cosmic_lg_opch = new vector<short>;
  vector<short> *beam_hg_opch = new vector<short>;
  vector<short> *beam_lg_opch = new vector<short>;
  vector<double> *cosmic_hg_timestamp = new vector<double>;
  vector<double> *cosmic_lg_timestamp = new vector<double>;
  vector<double> *beam_hg_timestamp = new vector<double>;
  vector<double> *beam_lg_timestamp = new vector<double>;
  std::vector<float> *op_gain = new std::vector<float>;
  std::vector<float> *op_gainerror = new std::vector<float>; 
  double triggerTime;
  
  TClonesArray* op_wf = new TClonesArray("TH1S");
  std::vector<short> *op_femch = new std::vector<short>;
  std::vector<double> *op_timestamp = new std::vector<double>;
  int event_no, run_no, subrun_no;
  
  T->SetBranchAddress("cosmic_hg_wf",&cosmic_hg_wf);
  T->SetBranchAddress("cosmic_lg_wf",&cosmic_lg_wf);
  T->SetBranchAddress("beam_hg_wf",&beam_hg_wf);
  T->SetBranchAddress("beam_lg_wf",&beam_lg_wf);
  T->SetBranchAddress("cosmic_hg_opch",&cosmic_hg_opch);
  T->SetBranchAddress("cosmic_lg_opch",&cosmic_lg_opch);
  T->SetBranchAddress("beam_hg_opch",&beam_hg_opch);
  T->SetBranchAddress("beam_lg_opch",&beam_lg_opch);
  T->SetBranchAddress("cosmic_hg_timestamp",&cosmic_hg_timestamp);
  T->SetBranchAddress("cosmic_lg_timestamp",&cosmic_lg_timestamp);
  T->SetBranchAddress("beam_hg_timestamp",&beam_hg_timestamp);
  T->SetBranchAddress("beam_lg_timestamp",&beam_lg_timestamp);
  T->SetBranchAddress("op_gain",&op_gain);
  T->SetBranchAddress("op_gainerror",&op_gainerror);
  T->SetBranchAddress("triggerTime",&triggerTime);
  T->SetBranchStatus("eventNo",1);
  T->SetBranchAddress("eventNo" , &event_no);
  T->SetBranchStatus("runNo",1);
  T->SetBranchAddress("runNo"   , &run_no);
  T->SetBranchStatus("subRunNo",1);
  T->SetBranchAddress("subRunNo", &subrun_no);

  Float_t elifetime;
  bool flag_elifetime = false;
  if (T->GetBranch("elifetime")){
    flag_elifetime = true;
    T->SetBranchAddress("elifetime",&elifetime);
  }

  
  int nEntries = T->GetEntries();  

  TFile** file = new TFile*[nEntries];    

  for(int i=0; i<nEntries; i++){
    T->GetEntry(i);
    std::cout<<"Run "<<run_no<<"   Subrun "<<subrun_no<<"   Event "<<i<<std::endl;
    file[i] = new TFile(Form("processedFlash_%d_%d_%d.root",run_no, subrun_no,i),"RECREATE");
    TTree *T_flash = new TTree("T_flash","T_flash");
    T_flash->SetDirectory(file[i]);
    int type;
    double low_time, high_time;
    double time;
    double total_PE;
    double PE[32],PE_err[32];
    std::vector<int> fired_channels;
    std::vector<double> l1_fired_time;
    std::vector<double> l1_fired_pe;
    T_flash->Branch("runNo",&run_no);
    T_flash->Branch("subRunNo",&subrun_no);
    T_flash->Branch("eventNo",&event_no);
    T_flash->Branch("triggerTime",&triggerTime);
    T_flash->Branch("type",&type);
    T_flash->Branch("low_time",&low_time);
    T_flash->Branch("high_time",&high_time);
    T_flash->Branch("time",&time);
    T_flash->Branch("total_PE",&total_PE);
    T_flash->Branch("PE",PE,"PE[32]/D");
    T_flash->Branch("PE_err",PE_err,"PE_err[32]/D");
    T_flash->Branch("fired_channels",&fired_channels);
    T_flash->Branch("l1_fired_time",&l1_fired_time);
    T_flash->Branch("l1_fired_pe",&l1_fired_pe);

    if (flag_elifetime)
      T_flash->Branch("elifetime",&elifetime);
    
    uboone_flash.load_event_raw(i);
    WCP::OpflashSelection& flashes = uboone_flash.get_flashes();
    for (auto it = flashes.begin(); it!=flashes.end(); it++){
      fired_channels.clear();
      Opflash *flash = (*it);
      type = flash->get_type();
      low_time = flash->get_low_time();
      high_time = flash->get_high_time();
      time = flash->get_time();
      total_PE = flash->get_total_PE();
      for (int i=0;i!=32;i++){
	PE[i] = flash->get_PE(i);
	PE_err[i] = flash->get_PE_err(i);
	if (flash->get_fired(i))
	  fired_channels.push_back(i);
      }
      l1_fired_time = flash->get_l1_fired_time();
      l1_fired_pe = flash->get_l1_fired_pe();
      T_flash->Fill();

     
    }
    uboone_flash.clear_flashes();
    
    file[i]->Write();
    file[i]->Close();
  }
  return 1;
}
