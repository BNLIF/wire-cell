#include "WireCellSst/GeomDataSource.h"
#include "WireCellSst/DatauBooNEFrameDataSource.h"
#include "WireCellSst/ToyuBooNESliceDataSource.h"
#include "WireCell2dToy/ToyEventDisplay.h"
#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/BadTiling.h"

#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/TruthToyTiling.h"
#include "WireCell2dToy/SimpleBlobToyTiling.h"

#include "WireCell2dToy/ToyMatrix.h"
#include "WireCell2dToy/ToyMatrixExclusive.h"
#include "WireCell2dToy/ToyMatrixKalman.h"
#include "WireCell2dToy/ToyMatrixIterate.h"
#include "WireCell2dToy/ToyMatrixIterate_SingleWire.h"
#include "WireCell2dToy/ToyMatrixIterate_Only.h"


#include "WireCell2dToy/ToyMatrixMarkov.h"
#include "WireCell2dToy/ToyMetric.h"
#include "WireCell2dToy/BlobMetric.h"
#include "WireCellData/TPCParams.h"
#include "WireCellData/Singleton.h"

#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/MergeGeomWire.h"

#include "WireCellData/GeomCluster.h"
//#include "WireCellNav/SliceDataSource.h"


#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/SimDataSource.h"
#include "WireCellNav/SliceDataSource.h"
#include "WireCellSst/Util.h"
#include "WireCellData/SimTruth.h"
#include "WireCell2dToy/ToyDepositor.h"
#include "WireCellNav/GenerativeFDS.h"
#include "WireCell2dToy/ToySignalSimu.h"
#include "WireCell2dToy/ToySignalSimuTrue.h"
#include "WireCell2dToy/DataSignalGaus_ROI.h"
#include "WireCell2dToy/DataSignalWien_ROI.h"

#include "WireCell2dToy/uBooNE_Data_2D_Deconvolution.h"
#include "WireCell2dToy/uBooNE_Data_ROI.h"
#include "WireCell2dToy/uBooNE_Data_After_ROI.h"
#include "WireCell2dToy/uBooNE_Data_After_ROI_gaus.h"
#include "WireCell2dToy/ExecMon.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TColor.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include <iostream>

using namespace WireCell;
using namespace std;


int main(int argc, char* argv[])
{
  if (argc < 4) {
    cout << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root eve_num -t[0,1] -s[0,1,2]" << endl;
    return 1;
  }

  TH1::AddDirectory(kFALSE);

  int two_plane = 0;
  int save_file = 0;
  int nt_off1 = 0;
  int nt_off2 = 0;
  int dead_region = 0; //enable manully added dead region, -d1
  for(Int_t i = 1; i != argc; i++){
     switch(argv[i][1]){
     case 't':
       two_plane = atoi(&argv[i][2]); 
       break;
     case 's':
       save_file = atoi(&argv[i][2]); 
       break;
     case 'd':
       dead_region = atoi(&argv[i][2]); 
       break;
     case 'a':
       nt_off1 = atoi(&argv[i][2]);
       break;
     case 'b':
       nt_off2 = atoi(&argv[i][2]);
       break;
     }
  }
  
  if (two_plane)
    cout << "Enable Two Plane Reconstruction " << endl; 

  ExecMon em("starting");

  
  WireCellSst::GeomDataSource gds(argv[1]);
  std::vector<double> ex = gds.extent();
  cout << "Extent: "
       << " x:" << ex[0]/units::mm << " mm"
       << " y:" << ex[1]/units::m << " m"
       << " z:" << ex[2]/units::m << " m"
       << endl;

  cout << "Pitch: " << gds.pitch(WirePlaneType_t(0)) 
       << " " << gds.pitch(WirePlaneType_t(1)) 
       << " " << gds.pitch(WirePlaneType_t(2))
       << endl;
  cout << "Angle: " << gds.angle(WirePlaneType_t(0)) 
       << " " << gds.angle(WirePlaneType_t(1)) 
       << " " << gds.angle(WirePlaneType_t(2))
       << endl;
   
  //float unit_dis = 1.01483;  // 58KV @ 226.5 V/cm
  //float unit_dis = 1.14753;  // 70 KV @ 226.5 V/cm
  //float unit_dis = 1.119;  // 70 KV @ 273 V/cm
  //float unit_dis = 1.105;  // doc-db 6683 matched with 256 cm
  //Note: the above one is still on the high end, 
  float unit_dis = 1.114; // match 256 cm

  
  float toffset_1=0.0; //(nt_off1 * 0.2 - 1.0 );  // time offset between u/v 
  float toffset_2=0.0; //(nt_off2 * 0.2 - 1.0); // time offset between u/w
  float toffset_3=0.0;
  
  int save_image_outline_flag = 0; // prescale flag 
  
  int total_time_bin = 9592;
  int recon_threshold = 2000;
  int frame_length = 3200;
  int max_events = 100;
  int eve_num  = atoi(argv[3]);
  int nrebin = 4;

  TPCParams& mp = Singleton<TPCParams>::Instance();
  
  double pitch_u = gds.pitch(WirePlaneType_t(0));
  double pitch_v = gds.pitch(WirePlaneType_t(1));
  double pitch_w = gds.pitch(WirePlaneType_t(2));
  double time_slice_width = nrebin * unit_dis * 0.5 * units::mm;

  mp.set_pitch_u(pitch_u);
  mp.set_pitch_v(pitch_v);
  mp.set_pitch_w(pitch_w);
  mp.set_ts_width(time_slice_width);
  
  std::cout << "Singleton: " << mp.get_pitch_u() << " " << mp.get_pitch_v() << " " << mp.get_pitch_w() << " " << mp.get_ts_width() << std::endl;
  
  float threshold_u = 5.87819e+02 * 4.0;
  float threshold_v = 8.36644e+02 * 4.0;
  float threshold_w = 5.67974e+02 * 4.0;

  float threshold_ug = 755.96;
  float threshold_vg = 822.81;
  float threshold_wg = 510.84;
  
  int time_offset = -90; // Now the time offset is taken care int he signal processing, so we just need the overall offset ... 
  


  const char* root_file = argv[2];
 
  
  int run_no, subrun_no, event_no;
  // sst->SetBranchAddress("eventNo",&event_no);
  // sst->SetBranchAddress("runNo",&run_no);
  // sst->SetBranchAddress("subRunNo",&subrun_no);
  // sst->GetEntry(eve_num);
  
  
  WireCellSst::DatauBooNEFrameDataSource *data_fds = new WireCellSst::DatauBooNEFrameDataSource(root_file,gds,total_time_bin, 1); // 1 means inactive noise filtering. Default: 0
  if (save_file != 2){
    data_fds->jump(eve_num);
    if (save_file == 1)
      data_fds->Save();
  }

  cout << em("noise filtering") << endl;
  
  run_no = data_fds->get_run_no();
  subrun_no = data_fds->get_subrun_no();
  event_no = data_fds->get_event_no();
  
  cout << "Run No: " << run_no << " " << subrun_no << " " << event_no << endl;

  ChirpMap uplane_map = data_fds->get_u_cmap();
  ChirpMap vplane_map = data_fds->get_v_cmap();
  ChirpMap wplane_map = data_fds->get_w_cmap();
  
  std::set<int> lf_noisy_channels = data_fds->get_lf_noisy_channels();
  

  cout << "Deconvolution with Wiener filter" << endl; 
  WireCell2dToy::uBooNEData2DDeconvolutionFDS *wien_fds = new WireCell2dToy::uBooNEData2DDeconvolutionFDS(*data_fds,gds,uplane_map, vplane_map, wplane_map,100,toffset_1,toffset_2,toffset_3, 1); //default: 0 -- for data channel-by-channel correction; 1: for simulation without any channel-by-channel correction
  wien_fds->jump(eve_num);

  cout << em("2D deconvolution") << endl;
  
  WireCell2dToy::uBooNEDataROI *uboone_rois = new WireCell2dToy::uBooNEDataROI(*data_fds,*wien_fds,gds,uplane_map,vplane_map,wplane_map,lf_noisy_channels);
  WireCell2dToy::uBooNEDataAfterROI roi_fds(*wien_fds,gds,*uboone_rois,nrebin);
  roi_fds.jump(eve_num);
  WireCell2dToy::uBooNEDataAfterROI_Gaus roi_gaus_fds(wien_fds, &roi_fds,gds);
  roi_gaus_fds.jump(eve_num);

  cout << em("ROI finder") << endl;
  
  std::vector<float> uplane_rms = uboone_rois->get_uplane_rms();
  std::vector<float> vplane_rms = uboone_rois->get_vplane_rms();
  std::vector<float> wplane_rms = uboone_rois->get_wplane_rms();

  
  delete wien_fds;
  delete uboone_rois;
  roi_fds.Clear();
  
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  int nwire_u = wires_u.size();
  int nwire_v = wires_v.size();
  int nwire_w = wires_w.size();
  
   
  TFile *f = new TFile(root_file);
  TTree *t = (TTree*)f->Get("Event/Sim");
  TBranch *br = t->GetBranch("oh_nHits");

  //save the image before deghosting ... 
  TFile *file = new TFile(Form("data-sp_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");
  
  //save light information in here ... 

  // OP HIT
  Int_t oh_nHits;
  std::vector<Int_t> * oh_channel = 0;
  std::vector<Double_t> * oh_bgtime = 0;
  std::vector<Double_t> * oh_trigtime = 0;
  std::vector<Double_t> * oh_pe = 0;

  // OP FLASH
  Int_t of_nFlash;
  std::vector<Float_t> * of_t = 0;
  std::vector<Float_t> * of_peTotal = 0;
  std::vector<Int_t> * of_multiplicity = 0;
  TClonesArray* pe_opdet = 0;
  
  if (br!=0){
    t->SetBranchAddress("oh_nHits", &oh_nHits);
    t->SetBranchAddress("oh_channel",&oh_channel);
    t->SetBranchAddress("oh_bgtime", &oh_bgtime);
    t->SetBranchAddress("oh_trigtime", &oh_trigtime);
    t->SetBranchAddress("oh_pe", &oh_pe);
    t->SetBranchAddress("of_nFlash", &of_nFlash);
    t->SetBranchAddress("of_t", &of_t);
    t->SetBranchAddress("of_peTotal", &of_peTotal);
    t->SetBranchAddress("of_multiplicity", &of_multiplicity);
    t->SetBranchAddress("pe_opdet", &pe_opdet);
    
    t->GetEntry(eve_num);
    
    TTree *T_op = new TTree("T_op","T_op");
    T_op->SetDirectory(file);
    
    T_op->Branch("oh_nHits",&oh_nHits);  
    T_op->Branch("oh_channel",&oh_channel);
    T_op->Branch("oh_bgtime", &oh_bgtime);
    T_op->Branch("oh_trigtime", &oh_trigtime);
    T_op->Branch("oh_pe", &oh_pe);
    T_op->Branch("of_nFlash", &of_nFlash);
    T_op->Branch("of_t", &of_t);
    T_op->Branch("of_peTotal", &of_peTotal);
    T_op->Branch("of_multiplicity", &of_multiplicity);
    T_op->Branch("pe_opdet", &pe_opdet);
    
    T_op->Fill();
  }
  
  
  TTree *Trun = new TTree("Trun","Trun");
  Trun->SetDirectory(file);

  int detector = 0; // MicroBooNE
  Trun->Branch("detector",&detector,"detector/I");

  Trun->Branch("eventNo",&event_no,"eventNo/I");
  Trun->Branch("runNo",&run_no,"runNo/I");
  Trun->Branch("subRunNo",&subrun_no,"runRunNo/I");
  
  Trun->Branch("unit_dis",&unit_dis,"unit_dis/F");
  Trun->Branch("toffset_uv",&toffset_1,"toffset_uv/F");
  Trun->Branch("toffset_uw",&toffset_2,"toffset_uw/F");
  Trun->Branch("toffset_u",&toffset_3,"toffset_u/F");
  Trun->Branch("total_time_bin",&total_time_bin,"total_time_bin/I");
  Trun->Branch("recon_threshold",&recon_threshold,"recon_threshold/I");
  Trun->Branch("frame_length",&frame_length,"frame_length/I");
  Trun->Branch("max_events",&max_events,"max_events/I");
  Trun->Branch("eve_num",&eve_num,"eve_num/I");
  Trun->Branch("nrebin",&nrebin,"nrebin/I");
  Trun->Branch("threshold_u",&threshold_u,"threshold_u/F");
  Trun->Branch("threshold_v",&threshold_v,"threshold_v/F");
  Trun->Branch("threshold_w",&threshold_w,"threshold_w/F");
  Trun->Branch("threshold_ug",&threshold_ug,"threshold_ug/F");
  Trun->Branch("threshold_vg",&threshold_vg,"threshold_vg/F");
  Trun->Branch("threshold_wg",&threshold_wg,"threshold_wg/F");
  Trun->Branch("time_offset",&time_offset,"time_offset/I");

  pitch_u = pitch_u/units::cm;
  pitch_v = pitch_v/units::cm;
  pitch_w = pitch_w/units::cm;
  Trun->Branch("pitch_u",&pitch_u,"pitch_u/D");
  Trun->Branch("pitch_v",&pitch_v,"pitch_v/D");
  Trun->Branch("pitch_w",&pitch_w,"pitch_w/D");



  Trun->Fill();

  TTree *T_chirp = new TTree("T_bad","T_bad");
  Int_t chid, plane;
  Int_t start_time,end_time;
  T_chirp->Branch("chid",&chid,"chid/I");
  T_chirp->Branch("plane",&plane,"plane/I");
  T_chirp->Branch("start_time",&start_time,"start_time/I");
  T_chirp->Branch("end_time",&end_time,"end_time/I");
  T_chirp->SetDirectory(file);
  for (auto it = uplane_map.begin(); it!=uplane_map.end();it++){
    chid = it->first;
    plane = 0;
    start_time = it->second.first;
    end_time = it->second.second;
    T_chirp->Fill();
  }
  for (auto it = vplane_map.begin(); it!=vplane_map.end();it++){
    chid = it->first + nwire_u;
    plane = 1;
    start_time = it->second.first;
    end_time = it->second.second;
    T_chirp->Fill();
  }
  for (auto it = wplane_map.begin(); it!=wplane_map.end();it++){
    chid = it->first + nwire_u + nwire_v;
    plane = 2;
    start_time = it->second.first;
    end_time = it->second.second;
    T_chirp->Fill();
  }
 
  // If dead region mode, then manually add dead channels.
  if(dead_region)
  {
      for(int ch = 0; ch< 8256; ch++)
      {
            chid = ch;
            start_time = 0;
            end_time = total_time_bin-1;
          // if( (ch>=296 && ch<=671) || (ch>=864 && ch<=911) ) {
          //     plane = 0;
          //     T_chirp->Fill();
          // }
          // if( ch>=3936 && ch<=3983 ) {
          //     plane = 1;
          //     T_chirp->Fill();
          // }
          // if( ch>=7136 && ch<=7263 ) {
          //     plane = 2;
          //     T_chirp->Fill();
          // }

	  if( (ch>=864 && ch<=911) ||
	      (ch>=296 && ch<=296) ||
	      (ch>=298 && ch<=315) ||
	      (ch>=317 && ch<=317) ||
	      (ch>=319 && ch<=327) ||
	      (ch>=336 && ch<=337) ||
	      (ch>=343 && ch<=345) ||
	      (ch>=348 && ch<=351) ||
	      (ch>=376 && ch<=400) ||
	      (ch>=410 && ch<=445) ||
	      (ch>=447 && ch<=484) ||
	      (ch>=501 && ch<=503) ||
	      (ch>=505 && ch<=520) ||
	      (ch>=522 && ch<=524) ||
	      (ch>=536 && ch<=559) ||
	      (ch>=561 && ch<=592) ||
	      (ch>=595 && ch<=598) ||
	      (ch>=600 && ch<=632) ||
	      (ch>=634 && ch<=652) ||
	      (ch>=654 && ch<=654) ||
	      (ch>=656 && ch<=671)
	      ) {
	    plane = 0;
	    T_chirp->Fill();
	  }
	  if( ch>=3936 && ch<=3983 ) {
	    plane = 1;
	    T_chirp->Fill();
          }
	  if( (ch>=4800+2336 && ch<=4800+2399) ||
	      (ch>=4800+2401 && ch<=4800+2414) ||
	      (ch>=4800+2416 && ch<=4800+2463)
	      ) {
	    plane = 2;
	    T_chirp->Fill();
	  }
	  

	  
      }
  }

  TH1F *hu_threshold = new TH1F("hu_threshold","hu_threshold",nwire_u,-0.5,-0.5+nwire_u);
  TH1F *hv_threshold = new TH1F("hv_threshold","hv_threshold",nwire_v,-0.5+nwire_u,-0.5+nwire_u+nwire_v);
  TH1F *hw_threshold = new TH1F("hw_threshold","hw_threshold",nwire_w,-0.5+nwire_u+nwire_v,-0.5+nwire_u+nwire_v+nwire_w);
  hu_threshold->SetDirectory(file);
  hv_threshold->SetDirectory(file);
  hw_threshold->SetDirectory(file);
  
 for (Int_t i=0;i!=uplane_rms.size();i++){
   hu_threshold->SetBinContent(i+1,uplane_rms.at(i));
  }
  for (Int_t i=0;i!=vplane_rms.size();i++){
    hv_threshold->SetBinContent(i+1,vplane_rms.at(i));
  }
  for (Int_t i=0;i!=wplane_rms.size();i++){
    hw_threshold->SetBinContent(i+1,wplane_rms.at(i));
  }


  TH2F *hu_raw = new TH2F("hu_raw","hu_raw",nwire_u,-0.5,nwire_u-0.5,total_time_bin,0,total_time_bin);
  TH2F *hv_raw = new TH2F("hv_raw","hv_raw",nwire_v,-0.5+nwire_u,nwire_v-0.5+nwire_u,total_time_bin,0,total_time_bin);
  TH2F *hw_raw = new TH2F("hw_raw","hw_raw",nwire_w,-0.5+nwire_u+nwire_v,nwire_w-0.5+nwire_u+nwire_v,total_time_bin,0,total_time_bin);
  hu_raw->SetDirectory(file);
  hv_raw->SetDirectory(file);
  hw_raw->SetDirectory(file);
  TH2F *htemp;
  
  const Frame& frame = data_fds->get();
  size_t ntraces = frame.traces.size();
  for (size_t ind=0; ind<ntraces; ++ind) {
    const Trace& trace = frame.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();
    WirePlaneType_t plane = gds.by_channel(chid).at(0)->plane();
    if (plane == WirePlaneType_t(0)){
       htemp = hu_raw;
      //continue;
    }else if (plane == WirePlaneType_t(1)){
      htemp = hv_raw;
      chid -= nwire_u;
      //if (chid < 1166 || chid > 1905) continue;
    }else if (plane == WirePlaneType_t(2)){
      htemp = hw_raw;
      chid -= nwire_u + nwire_v;
      //continue;
    }
    for (int i = tbin;i!=tbin+nbins;i++){
      int tt = i+1;
      htemp->SetBinContent(chid+1,tt,trace.charge.at(i));
    }
  }
  delete data_fds;

  TH2F *hu_decon = new TH2F("hu_decon","hu_decon",nwire_u,-0.5,nwire_u-0.5,total_time_bin/nrebin,0,total_time_bin);
  TH2F *hv_decon = new TH2F("hv_decon","hv_decon",nwire_v,-0.5+nwire_u,nwire_v-0.5+nwire_u,total_time_bin/nrebin,0,total_time_bin);
  TH2F *hw_decon = new TH2F("hw_decon","hw_decon",nwire_w,-0.5+nwire_u+nwire_v,nwire_w-0.5+nwire_u+nwire_v,total_time_bin/nrebin,0,total_time_bin);
  hu_decon->SetDirectory(file);
  hv_decon->SetDirectory(file);
  hw_decon->SetDirectory(file);
  TH2F *htemp1;
  const Frame& frame1 = roi_fds.get();
  ntraces = frame1.traces.size();
  for (size_t ind=0; ind<ntraces; ++ind) {
    const Trace& trace = frame1.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();
    WirePlaneType_t plane = gds.by_channel(chid).at(0)->plane();
    if (plane == WirePlaneType_t(0)){
      htemp1 = hu_decon;
    }else if (plane == WirePlaneType_t(1)){
      htemp1 = hv_decon;
      chid -= nwire_u;
    }else if (plane == WirePlaneType_t(2)){
      htemp1 = hw_decon;
      chid -= nwire_u + nwire_v;
    }
    for (int i = tbin;i!=tbin+nbins;i++){
      int tt = i+1;
      htemp1->SetBinContent(chid+1,tt,int(trace.charge.at(i)));
    }
  }

  TH2F *hu_decon_g = new TH2F("hu_decon_g","hu_decon_g",nwire_u,-0.5,nwire_u-0.5,total_time_bin/nrebin,0,total_time_bin);
  TH2F *hv_decon_g = new TH2F("hv_decon_g","hv_decon_g",nwire_v,-0.5+nwire_u,nwire_v-0.5+nwire_u,total_time_bin/nrebin,0,total_time_bin);
  TH2F *hw_decon_g = new TH2F("hw_decon_g","hw_decon_g",nwire_w,-0.5+nwire_u+nwire_v,nwire_w-0.5+nwire_u+nwire_v,total_time_bin/nrebin,0,total_time_bin);
  hu_decon_g->SetDirectory(file);
  hv_decon_g->SetDirectory(file);
  hw_decon_g->SetDirectory(file);
  const Frame& frame3 = roi_gaus_fds.get();
  ntraces = frame3.traces.size();
  for (size_t ind=0; ind<ntraces; ++ind) {
    const Trace& trace = frame3.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();
    WirePlaneType_t plane = gds.by_channel(chid).at(0)->plane();
    if (plane == WirePlaneType_t(0)){
      htemp1 = hu_decon_g;
    }else if (plane == WirePlaneType_t(1)){
      htemp1 = hv_decon_g;
      chid -= nwire_u;
    }else if (plane == WirePlaneType_t(2)){
      htemp1 = hw_decon_g;
      chid -= nwire_u + nwire_v;
    }
    for (int i = tbin;i!=tbin+nbins;i++){
      int tt = i+1;
      htemp1->SetBinContent(chid+1,tt,int(trace.charge.at(i)));
    }
  }
  

  cout << em("save file") << endl;

  file->Write();
  file->Close();  

  return 0;
  
} // main()
