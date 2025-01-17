#include "WCPSst/GeomDataSource.h"
#include "WCP2dToy/ToyLightReco.h"
#include "TH1F.h"
#include "TH2F.h"

#include <TLeaf.h>
#include <iostream>

using namespace WCP;
using namespace std;

int main(int argc, char* argv[])
{
  if (argc < 4) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root eve_num -i[0,1] -d[0,1,2]" << endl;
    return 1;
  }
  TH1::AddDirectory(kFALSE);

  int imagingoutput = 0;
  int datatier = 0; // data=0, overlay=1, full mc=2
  for(Int_t i=4; i!=argc; i++){
    switch(argv[i][1]){
    case 'i':
      imagingoutput = atoi(&argv[i][2]);
      break;
    case 'd':
      datatier = atoi(&argv[i][2]);
      break;
    }
  }

  double lowerwindow = 3;
  double upperwindow = 6;

  WCPSst::GeomDataSource gds(argv[1]);
  std::vector<double> ex = gds.extent();

  const char* root_file = argv[2];
  int eve_num = atoi(argv[3]);

  bool use_imagingoutput = false;
  if(imagingoutput==1) use_imagingoutput = true;

  WCP2dToy::ToyLightReco uboone_flash(root_file,use_imagingoutput,datatier); 

  uboone_flash.load_event_raw(eve_num);
  TFile *file1 = new TFile(root_file);
  TTree *T;
  if(!use_imagingoutput){ T = (TTree*)file1->Get("/Event/Sim"); } // celltree input
  else{ T = (TTree*)file1->Get("Trun"); }
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
  //Added for mc
  TClonesArray* op_wf = new TClonesArray("TH1S");
  std::vector<short> *op_femch = new std::vector<short>;
  std::vector<double> *op_timestamp = new std::vector<double>;
  int event_no, run_no, subrun_no;
  
  T->SetBranchAddress("cosmic_hg_wf",&cosmic_hg_wf);
  T->SetBranchAddress("cosmic_lg_wf",&cosmic_lg_wf);

  if (datatier == 1){
    T->SetBranchAddress("mixer_beam_hg_wf",&beam_hg_wf);
    T->SetBranchAddress("mixer_beam_lg_wf",&beam_lg_wf);
    T->SetBranchAddress("mixer_beam_hg_opch",&beam_hg_opch);
    T->SetBranchAddress("mixer_beam_lg_opch",&beam_lg_opch);
    T->SetBranchAddress("mixer_beam_hg_timestamp",&beam_hg_timestamp);
    T->SetBranchAddress("mixer_beam_lg_timestamp",&beam_lg_timestamp);
  }else{
    T->SetBranchAddress("beam_hg_wf",&beam_hg_wf);
    T->SetBranchAddress("beam_lg_wf",&beam_lg_wf);
    T->SetBranchAddress("beam_hg_opch",&beam_hg_opch);
    T->SetBranchAddress("beam_lg_opch",&beam_lg_opch);
    T->SetBranchAddress("beam_hg_timestamp",&beam_hg_timestamp);
    T->SetBranchAddress("beam_lg_timestamp",&beam_lg_timestamp);
  }
  T->SetBranchAddress("cosmic_hg_opch",&cosmic_hg_opch);
  T->SetBranchAddress("cosmic_lg_opch",&cosmic_lg_opch);
  
  T->SetBranchAddress("cosmic_hg_timestamp",&cosmic_hg_timestamp);
  T->SetBranchAddress("cosmic_lg_timestamp",&cosmic_lg_timestamp);

  T->SetBranchAddress("op_gain",&op_gain);
  T->SetBranchAddress("op_gainerror",&op_gainerror);
  T->SetBranchAddress("triggerTime",&triggerTime);
  T->SetBranchStatus("eventNo",1);
  T->SetBranchAddress("eventNo" , &event_no);
  T->SetBranchStatus("runNo",1);
  T->SetBranchAddress("runNo"   , &run_no);
  T->SetBranchStatus("subRunNo",1);
  T->SetBranchAddress("subRunNo", &subrun_no);
  T->GetEntry(eve_num);
  std::cout<<"run "<<run_no<<" subrun "<<" event "<<std::endl;
  std::cout<<run_no<<"_"<<subrun_no<<"_"<<event_no<<std::endl;
 
  //TFile *acpttrig_file = new TFile(acpttrig_filename);
  TTree *T_acpttrig;
  int eve_num_acpttrig=-1;
  float _len;
  float _trk_beg_x, _trk_beg_y, _trk_beg_z;
  float _trk_end_x, _trk_end_y, _trk_end_z;
  float _trk_beg_x_off, _trk_beg_y_off, _trk_beg_z_off;
  float _trk_end_x_off, _trk_end_y_off, _trk_end_z_off;
  float _trk_len; //corrected by SCE
  int    _nflash;
  float _flash_t, _flash_pe, _flash_yw, _flash_yc, _flash_zw, _flash_zc;
  unsigned int _timehigh, _timelow;
  vector<double>* _gain_v = new vector<double>;
  vector<double>* _flash_pe_v = new vector<double>;
  vector<double>* _flash_pe_corr_v = new vector<double>;
  if(file1){
    int event_acpttrig, run_acpttrig, subrun_acpttrig;
    TTree *T_acpttrig = (TTree*)file1->Get("/acpttrig/_tree");
    if(T_acpttrig->GetBranch("_event")){
      T_acpttrig->SetBranchAddress("_event" , &event_acpttrig);
      T_acpttrig->SetBranchAddress("_run"   , &run_acpttrig);
      T_acpttrig->SetBranchAddress("_subrun", &subrun_acpttrig);
      T_acpttrig->SetBranchAddress("_nflash",&_nflash);
      T_acpttrig->SetBranchAddress("_flash_t",&_flash_t);
      T_acpttrig->SetBranchAddress("_flash_pe",&_flash_pe);
      T_acpttrig->SetBranchAddress("_gain_v",&_gain_v);
      T_acpttrig->SetBranchAddress("_flash_pe_v",&_flash_pe_v);
      T_acpttrig->SetBranchAddress("_flash_pe_corr_v",&_flash_pe_corr_v);
      T_acpttrig->SetBranchAddress("_timelow" ,&_timelow);
      T_acpttrig->SetBranchAddress("_timehigh",&_timehigh);
      T_acpttrig->SetBranchAddress("_trk_beg_x",&_trk_beg_x);
      T_acpttrig->SetBranchAddress("_trk_beg_y",&_trk_beg_y);
      T_acpttrig->SetBranchAddress("_trk_beg_z",&_trk_beg_z);
      T_acpttrig->SetBranchAddress("_trk_end_x",&_trk_end_x);
      T_acpttrig->SetBranchAddress("_trk_end_y",&_trk_end_y);
      T_acpttrig->SetBranchAddress("_trk_end_z",&_trk_end_z);
      T_acpttrig->SetBranchAddress("_len",&_len);
      T_acpttrig->SetBranchAddress("_trk_beg_x_off",&_trk_beg_x_off);
      T_acpttrig->SetBranchAddress("_trk_beg_y_off",&_trk_beg_y_off);
      T_acpttrig->SetBranchAddress("_trk_beg_z_off",&_trk_beg_z_off);
      T_acpttrig->SetBranchAddress("_trk_end_x_off",&_trk_end_x_off);
      T_acpttrig->SetBranchAddress("_trk_end_y_off",&_trk_end_y_off);
      T_acpttrig->SetBranchAddress("_trk_end_z_off",&_trk_end_z_off);
      T_acpttrig->SetBranchAddress("_trk_len",&_trk_len);
      T_acpttrig->SetBranchAddress("_flash_yw",&_flash_yw);
      T_acpttrig->SetBranchAddress("_flash_zw",&_flash_zw);
      T_acpttrig->SetBranchAddress("_flash_yc",&_flash_yc);
      T_acpttrig->SetBranchAddress("_flash_zc",&_flash_zc);
      for(int i=0; i<T_acpttrig->GetEntries();  i++){
        T_acpttrig->GetEntry(i);
	if(event_acpttrig==event_no && run_acpttrig==run_no && subrun_acpttrig==subrun_no){eve_num_acpttrig=i; break;}
      }
    }else{std::cout<<"WARNING: No acpttrig tree"<<std::endl; return 0;}
  }else{std::cout<<"WARNING: No acpttrig file"<<std::endl; return 0;}
  if(eve_num_acpttrig<0){std::cout<<"No matching acpttrig event"<<std::endl; return 0;}

  TFile *file = new TFile(Form("flash_%d_%d_%d.root",run_no, subrun_no, event_no),"RECREATE");
  
  int type;
  double low_time, high_time;
  double time;
  double total_PE;
  double PE[32],PE_err[32];
  std::vector<int> fired_channels;
  std::vector<double> l1_fired_time;
  std::vector<double> l1_fired_pe;
  int wc_nflash=0;
  
  WCP::OpflashSelection& flashes = uboone_flash.get_flashes();

  for (auto it = flashes.begin(); it!=flashes.end(); it++){

    fired_channels.clear();
    Opflash *flash = (*it);
    double temp_time = flash->get_time();
    if (temp_time <= lowerwindow || temp_time >= upperwindow) continue;
    wc_nflash+=1;

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
    
  }

  if(wc_nflash==0){std::cout<<"WARNING: no WC flashes found"<<std::endl;}
  if(wc_nflash>1){ std::cout<<"WARNING: multiple WC flashes found"<<std::endl;}


  TTree *T_acpttrig_clone = new TTree("T_acpttrig","T_acpttrig");
  T_acpttrig_clone->SetDirectory(file);

  T_acpttrig_clone->Branch("event" , &event_no);
  T_acpttrig_clone->Branch("run"   , &run_no);
  T_acpttrig_clone->Branch("subrun", &subrun_no);

  T_acpttrig_clone->Branch("_timelow" ,&_timelow ,"_timelow/i" );
  T_acpttrig_clone->Branch("_timehigh",&_timehigh,"_timehigh/i");
  T_acpttrig_clone->Branch("_trk_beg_x",&_trk_beg_x,"_trk_beg_x/F");
  T_acpttrig_clone->Branch("_trk_beg_y",&_trk_beg_y,"_trk_beg_y/F");
  T_acpttrig_clone->Branch("_trk_beg_z",&_trk_beg_z,"_trk_beg_z/F");
  T_acpttrig_clone->Branch("_trk_end_x",&_trk_end_x,"_trk_end_x/F");
  T_acpttrig_clone->Branch("_trk_end_y",&_trk_end_y,"_trk_end_y/F");
  T_acpttrig_clone->Branch("_trk_end_z",&_trk_end_z,"_trk_end_z/F");
  T_acpttrig_clone->Branch("_len",&_len,"_len/F");

  T_acpttrig_clone->Branch("_trk_beg_x_off",&_trk_beg_x_off,"_trk_beg_x_off/F");
  T_acpttrig_clone->Branch("_trk_beg_y_off",&_trk_beg_y_off,"_trk_beg_y_off/F");
  T_acpttrig_clone->Branch("_trk_beg_z_off",&_trk_beg_z_off,"_trk_beg_z_off/F");
  T_acpttrig_clone->Branch("_trk_end_x_off",&_trk_end_x_off,"_trk_end_x_off/F");
  T_acpttrig_clone->Branch("_trk_end_y_off",&_trk_end_y_off,"_trk_end_y_off/F");
  T_acpttrig_clone->Branch("_trk_end_z_off",&_trk_end_z_off,"_trk_end_z_off/F");
  T_acpttrig_clone->Branch("_trk_len",&_trk_len,"_trk_len/F");

  T_acpttrig_clone->Branch("_nflash",&_nflash,"_nflash/I");
  T_acpttrig_clone->Branch("_flash_t",&_flash_t,"_flash_t/F");
  T_acpttrig_clone->Branch("_flash_pe",&_flash_pe,"_flash_pe/F");
  T_acpttrig_clone->Branch("_gain_v",&_gain_v);
  T_acpttrig_clone->Branch("_flash_pe_v",&_flash_pe_v);
  T_acpttrig_clone->Branch("_flash_pe_corr_v",&_flash_pe_corr_v);

  T_acpttrig_clone->Branch("_flash_yw",&_flash_yw,"_flash_yw/F");
  T_acpttrig_clone->Branch("_flash_zw",&_flash_zw,"_flash_zw/F");
  T_acpttrig_clone->Branch("_flash_yc",&_flash_yc,"_flash_yc/F");
  T_acpttrig_clone->Branch("_flash_zc",&_flash_zc,"_flash_zc/F");

  T_acpttrig_clone->Branch("wc_type",&type);
  T_acpttrig_clone->Branch("wc_low_time",&low_time);
  T_acpttrig_clone->Branch("wc_high_time",&high_time);
  T_acpttrig_clone->Branch("wc_time",&time);
  T_acpttrig_clone->Branch("wc_total_PE",&total_PE);
  T_acpttrig_clone->Branch("wc_PE",PE,"wc_PE[32]/D");
  T_acpttrig_clone->Branch("wc_PE_err",PE_err,"wc_PE_err[32]/D");
  T_acpttrig_clone->Branch("wc_fired_channels",&fired_channels);
  T_acpttrig_clone->Branch("wc_l1_fired_time",&l1_fired_time);
  T_acpttrig_clone->Branch("wc_l1_fired_pe",&l1_fired_pe);
  T_acpttrig_clone->Branch("wc_nflash",&wc_nflash,"wc_nflash/I");
  T_acpttrig_clone->Fill();

  file->Write();
  file->Close();
  return 1;
}
