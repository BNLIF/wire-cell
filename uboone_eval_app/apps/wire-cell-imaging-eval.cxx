#include "WCPSst/GeomDataSource.h"
#include "WCPSst/DatauBooNEFrameDataSource.h"
#include "WCPSst/ToyuBooNESliceDataSource.h"
#include "WCPSst/uBooNESliceDataSource.h"

#include "WCP2dToy/ToyEventDisplay.h"
#include "WCP2dToy/ToyTiling.h"
#include "WCP2dToy/BadTiling.h"
#include "WCP2dToy/LowmemTiling.h"
#include "WCP2dToy/uBooNE_L1SP.h"
#include "WCP2dToy/WCPHolder.h"

#include "WCP2dToy/MergeToyTiling.h"
#include "WCP2dToy/TruthToyTiling.h"
#include "WCP2dToy/SimpleBlobToyTiling.h"

#include "WCP2dToy/ChargeSolving.h"
#include "WCP2dToy/ToyMatrix.h"
#include "WCP2dToy/ToyMatrixExclusive.h"
#include "WCP2dToy/ToyMatrixKalman.h"
#include "WCP2dToy/ToyMatrixIterate.h"
#include "WCP2dToy/ToyMatrixIterate_SingleWire.h"
#include "WCP2dToy/ToyMatrixIterate_Only.h"

#include "WCP2dToy/ToyMatrixMarkov.h"
#include "WCP2dToy/ToyMetric.h"
#include "WCP2dToy/BlobMetric.h"
#include "WCPData/TPCParams.h"
#include "WCPData/Singleton.h"

#include "WCPData/MergeGeomCell.h"
#include "WCPData/MergeGeomWire.h"

#include "WCPData/Slim3DCluster.h"
#include "WCPData/Slim3DDeadCluster.h"
//#include "WCPNav/SliceDataSource.h"


#include "WCPNav/FrameDataSource.h"
#include "WCPNav/SimDataSource.h"
#include "WCPNav/SliceDataSource.h"
#include "WCPSst/Util.h"
#include "WCPData/SimTruth.h"
#include "WCP2dToy/ToyDepositor.h"
#include "WCPNav/GenerativeFDS.h"
#include "WCP2dToy/ToySignalSimu.h"
#include "WCP2dToy/ToySignalSimuTrue.h"
#include "WCP2dToy/DataSignalGaus_ROI.h"
#include "WCP2dToy/DataSignalWien_ROI.h"

#include "WCP2dToy/uBooNE_Data_2D_Deconvolution.h"
#include "WCP2dToy/uBooNE_Data_ROI.h"
#include "WCP2dToy/uBooNE_Data_After_ROI.h"
#include "WCP2dToy/uBooNE_Data_After_ROI_gaus.h"
#include "WCP2dToy/pd_Data_FDS.h"
#include "WCP2dToy/uBooNE_Data_Error.h"
#include "WCP2dToy/ExecMon.h"
#include "WCP2dToy/ToyDataQuality.h"


#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TColor.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TRandom.h"
#include "TF1.h"
#include <iostream>

using namespace WCP;
using namespace std;

bool GeomWireSelectionCompare(GeomWireSelection a, GeomWireSelection b) {
  if (a.size() > b.size()){
    return true;
  }else if (a.size() < b.size()){
    return false;
  }else{
    for (int i=0;i!=a.size();i++){
      if (a.at(i)->ident() > b.at(i)->ident()){
        return true;
      }else if (a.at(i)->ident() < b.at(i)->ident()){
        return false;
      }
    }
    return false;
  }
  return false;
}


int main(int argc, char* argv[])
{
  if (argc < 3) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root /path/to/true.root file_tag -t[0,1] -s[0,1,2]" << endl;
    return 1;
  }

  TH1::AddDirectory(kFALSE);

  int two_plane = 0; //not used
  
  int nt_off1 = 0; // not used
  int nt_off2 = 0; // not used
  int dead_region = 0; // manually add dead channels

  int save_file = 0; //
  // 1 for debug mode for bee ...

  int flag_l1 = 0; // do not run l1sp code 
  
  int landau_fluc = 0; // landau fluctuation for each charge hit
  gRandom->SetSeed(0);
  TF1 *ld = new TF1("ld","[0]*TMath::Landau(x, 1.8, 0.15, 1)",0,10); // MPV: 1.8, psi: 0.15, 4*psi = FMHW
  ld->SetParameter(0, 1);
  
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
     case 'l':
       flag_l1 = atoi(&argv[i][2]);
       break;
     case 'f':
       landau_fluc = atoi(&argv[i][2]);
       break;
     }
  }

  if(save_file==1){
    std::cout << "Save file for bee. " << std::endl;
  }else if (save_file==0){
    std::cout << "Save file for pattern recognition. " << std::endl;
  }else if (save_file==2){
    std::cout << "Save file for clustering and imaging evaluation. " << std::endl;
  }
  
  // if (two_plane)
  //   cout << "Enable Two Plane Reconstruction " << endl; 
  
  ExecMon em("starting");
  cout << em("load geometry") << endl;

  WCPSst::GeomDataSource gds(argv[1]);
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
      
  //Note: the above one is still on the high end, 
  float unit_dis = 1.114; // match 256 cm

  //**** time offset for 58kV ****// 
  float toffset_1=0.0; //(nt_off1 * 0.2 - 1.0 );  // time offset between u/v 
  float toffset_2=0.0; //(nt_off2 * 0.2 - 1.0); // time offset between u/w
  float toffset_3=0.0;
  
  int save_image_outline_flag = 0; // prescale flag 
  
  int total_time_bin = 9592;
  int recon_threshold = 2000;
  int frame_length = 3200;
  int max_events = 100;
  int eve_num = 0;
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
  

  int time_offset = 0; // Now the time offset is taken care int he signal processing, so we just need the overall offset ... 
  
  //root file saving truth info
  const char* true_file;
  if(argv[3][0] != '-'){
    true_file = argv[3];  
  }
  TFile* file0 = new TFile(true_file);
  TTree* T_true = (TTree*)file0->Get("T_true");
  TTree* T_track = (TTree*)file0->Get("T_track");

  const char* root_file = argv[2];
  int run_no, subrun_no, event_no;
  
  cout << em("load data") << endl;

  // load Trun
  TFile *file1 = new TFile(root_file);
  TTree *Trun0 = (TTree*)file1->Get("Trun");
  
  Trun0->SetBranchAddress("eventNo",&event_no);
  Trun0->SetBranchAddress("runNo",&run_no);
  Trun0->SetBranchAddress("subRunNo",&subrun_no);
  Trun0->SetBranchAddress("unit_dis",&unit_dis);
  Trun0->SetBranchAddress("toffset_uv",&toffset_1);
  Trun0->SetBranchAddress("toffset_uw",&toffset_2);
  Trun0->SetBranchAddress("toffset_u",&toffset_3);
  Trun0->SetBranchAddress("total_time_bin",&total_time_bin);
  Trun0->SetBranchAddress("recon_threshold",&recon_threshold);
  Trun0->SetBranchAddress("frame_length",&frame_length);
  Trun0->SetBranchAddress("max_events",&max_events);
  Trun0->SetBranchAddress("eve_num",&eve_num);
  Trun0->SetBranchAddress("nrebin",&nrebin);
  Trun0->SetBranchAddress("threshold_u",&threshold_u);
  Trun0->SetBranchAddress("threshold_v",&threshold_v);
  Trun0->SetBranchAddress("threshold_w",&threshold_w);
  Trun0->SetBranchAddress("threshold_ug",&threshold_ug);
  Trun0->SetBranchAddress("threshold_vg",&threshold_vg);
  Trun0->SetBranchAddress("threshold_wg",&threshold_wg);
  Trun0->SetBranchAddress("time_offset",&time_offset);
  
  Trun0->GetEntry(0);

  TTree *T_op = (TTree*)file1->Get("T_op");
  
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  int nwire_u = wires_u.size();
  int nwire_v = wires_v.size();
  int nwire_w = wires_w.size();


  cout << "Run No: " << run_no << " " << subrun_no << " " << event_no << endl;
  ChirpMap uplane_map;
  ChirpMap vplane_map;
  ChirpMap wplane_map;
  TTree *T_chirp = (TTree*)file1->Get("T_bad");
  Int_t chid=0, plane=0;
  Int_t start_time=0, end_time=0;
  T_chirp->SetBranchAddress("chid",&chid);
  T_chirp->SetBranchAddress("plane",&plane);
  T_chirp->SetBranchAddress("start_time",&start_time);
  T_chirp->SetBranchAddress("end_time",&end_time);
  for (Int_t i=0;i!=T_chirp->GetEntries();i++){
    T_chirp->GetEntry(i);
    std::pair<int,int> abc(start_time,end_time);
    if (plane == 0){
      uplane_map[chid] = abc;
    }else if (plane == 1){
      vplane_map[chid-nwire_u] = abc;
    }else if (plane == 2){
      wplane_map[chid-nwire_u-nwire_v] = abc;
    }
  }
  
  // If dead region mode, then manually add dead channels.
  if(dead_region)
  {
      std::pair<int,int> abc(0, total_time_bin-1);
      for(int ch = 0; ch< 8256; ch++)
      {
          // if( (ch>=296 && ch<=671) || (ch>=864 && ch<=911) ) {
          //   uplane_map[ch] = abc;
          // }
          // if( ch>=3936 && ch<=3983 ) {
          //   vplane_map[ch - nwire_u] = abc;
          // }
          // if( ch>=7136 && ch<=7263 ) {
          //   wplane_map[ch - nwire_u - nwire_v] = abc;
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
	    uplane_map[ch] = abc;
	  }
	  if( ch>=3936 && ch<=3983 ) {
	    vplane_map[ch - nwire_u] = abc;
          }
	  if( (ch>=4800+2336 && ch<=4800+2399) ||
	      (ch>=4800+2401 && ch<=4800+2414) ||
	      (ch>=4800+2416 && ch<=4800+2463)
	      ) {
	    wplane_map[ch - nwire_u - nwire_v] = abc;
	  }	  
	
      }
  }
   
  
  std::vector<float> uplane_rms;
  std::vector<float> vplane_rms;
  std::vector<float> wplane_rms;
  // Note, there is a mismatch here
  // These RMS values are for single ticks
  // Later they are applied to the rebinned data
  // Probably OK as the TPC signal processing already took care the fake hits ... 
  TH1F *hu_threshold = (TH1F*)file1->Get("hu_threshold");
  TH1F *hv_threshold = (TH1F*)file1->Get("hv_threshold");
  TH1F *hw_threshold = (TH1F*)file1->Get("hw_threshold");
  for (int i=0;i!=hu_threshold->GetNbinsX();i++){
    uplane_rms.push_back(hu_threshold->GetBinContent(i+1));
  }
  for (int i=0;i!=hv_threshold->GetNbinsX();i++){
    vplane_rms.push_back(hv_threshold->GetBinContent(i+1));
  }
  for (int i=0;i!=hw_threshold->GetNbinsX();i++){
    wplane_rms.push_back(hw_threshold->GetBinContent(i+1));
  }
  
  TH2F *hu_decon = (TH2F*)file1->Get("hu_decon");
  TH2F *hv_decon = (TH2F*)file1->Get("hv_decon");
  TH2F *hw_decon = (TH2F*)file1->Get("hw_decon");
  
  TH2F *hu_decon_g = (TH2F*)file1->Get("hu_decon_g");
  TH2F *hv_decon_g = (TH2F*)file1->Get("hv_decon_g");
  TH2F *hw_decon_g = (TH2F*)file1->Get("hw_decon_g");


  if(landau_fluc){
    Double_t random_factor = 1.0;
    for(int i=1; i<=hu_decon->GetNbinsX(); i++){
        for(int j=1; j<=hu_decon->GetNbinsY(); j++){
            random_factor = ld->GetRandom()/1.8;
            Double_t content = hu_decon->GetBinContent(i, j);
            Double_t content_g = hu_decon_g->GetBinContent(i, j);
            hu_decon->SetBinContent(i, j, content*random_factor);
            hu_decon_g->SetBinContent(i, j, content_g*random_factor);
        }
    }
    for(int i=1; i<=hv_decon->GetNbinsX(); i++){
        for(int j=1; j<=hv_decon->GetNbinsY(); j++){
            random_factor = ld->GetRandom()/1.8;
            Double_t content = hv_decon->GetBinContent(i, j);
            Double_t content_g = hv_decon_g->GetBinContent(i, j);
            hv_decon->SetBinContent(i, j, content*random_factor);
            hv_decon_g->SetBinContent(i, j, content_g*random_factor);
        }
    }
    for(int i=1; i<=hw_decon->GetNbinsX(); i++){
        for(int j=1; j<=hw_decon->GetNbinsY(); j++){
            random_factor = ld->GetRandom()/1.8;
            Double_t content = hw_decon->GetBinContent(i, j);
            Double_t content_g = hw_decon_g->GetBinContent(i, j);
            hw_decon->SetBinContent(i, j, content*random_factor);
            hw_decon_g->SetBinContent(i, j, content_g*random_factor);
        }
    }
  }

  // add a special treatment here ...
  /* for (int i=296; i!=671;i++){ */
  /*   if (uplane_map.find(i)==uplane_map.end()){ */
  /*     if ((uplane_map.find(i-1)!=uplane_map.end() && */
	   /* uplane_map.find(i+2)!=uplane_map.end()) || */
	  /* (uplane_map.find(i-2)!=uplane_map.end() && */
	   /* uplane_map.find(i+1)!=uplane_map.end()) */
	  /* ){ */
	/* std::cout << "U plane (shorted): " << i << " added to bad chanel list" << std::endl; */
	/* uplane_map[i] = uplane_map[i-1]; */
	/* for (int j=0;j!=hu_decon->GetNbinsY();j++){ */
	  /* hu_decon->SetBinContent(i+1,j+1,0); */
	  /* hu_decon_g->SetBinContent(i+1,j+1,0); */
	/* } */
  /*     } */
  /*   } */
  /* } */
  /* for (int i=2336;i!=2463;i++){ */
  /*   if (wplane_map.find(i)==wplane_map.end()){ */
  /*     if ((wplane_map.find(i-1)!=wplane_map.end() && */
	   /* wplane_map.find(i+2)!=wplane_map.end()) || */
	  /* (wplane_map.find(i-2)!=wplane_map.end() && */
	   /* wplane_map.find(i+1)!=wplane_map.end()) */
	  /* ){ */
	/* wplane_map[i] = wplane_map[i-1]; */
	/* std::cout << "W plane (shorted): " << i << " added to bad channel list" << std::endl; */
	/* for (int j=0;j!=hw_decon->GetNbinsY();j++){ */
	  /* hw_decon->SetBinContent(i+1,j+1,0); */
	  /* hw_decon_g->SetBinContent(i+1,j+1,0); */
	/* } */
  /*     } */
  /*   } */
  /* } */
  //
  
  TH2F *hv_raw = (TH2F*)file1->Get("hv_raw");
  
  // V wire noisy channels 10 vetoed ...
  /* for (int i=3684;i!=3699;i++){ */
  /*   if (vplane_map.find(i-2400)==vplane_map.end()){ */
  /*     vplane_map[i-2400] = std::make_pair(0,hv_raw->GetNbinsY()-1); */
  /*     std::cout << "V plane (noisy): " << i -2400 << " added to bad channel list" << std::endl; */
  /*   }else{ */
  /*     vplane_map[i-2400] = std::make_pair(0,hv_raw->GetNbinsY()-1); */
  /*   } */
  /*   for (int j=0;j!=hv_decon->GetNbinsY();j++){ */
  /*     hv_decon->SetBinContent(i+1-2400,j+1,0); */
  /*     hv_decon_g->SetBinContent(i+1-2400,j+1,0); */
  /*   } */
  /*   for (int j=0;j!=hv_raw->GetNbinsY();j++){ */
  /*     hv_raw->SetBinContent(i+1-2400,j+1,0); */
  /*   } */
  /* } */
  // U wire plane bad channels
  /* for (int i=2160;i!=2176;i++){ */
  /*   if (uplane_map.find(i)==uplane_map.end()){ */
  /*     uplane_map[i] = std::make_pair(0,hv_raw->GetNbinsY()-1); */
  /*     std::cout << "U plane (bad): " << i << " added to bad channel list" << std::endl; */
  /*   }else{ */
  /*     uplane_map[i] = std::make_pair(0,hv_raw->GetNbinsY()-1); */
  /*   } */
  /*   for (int j=0;j!=hu_decon->GetNbinsY();j++){ */
  /*     hu_decon->SetBinContent(i+1,j+1,0); */
  /*     hu_decon_g->SetBinContent(i+1,j+1,0); */
  /*   } */
  /* } */

  // loop through U/V/W plane to disable the bad channels completely
  for (auto it = uplane_map.begin(); it!=uplane_map.end(); it++){
    int ch = it->first;
    int start = it->second.first/nrebin;
    int end = it->second.second/nrebin;
    for (int j=start; j!=end+1;j++){
      hu_decon->SetBinContent(ch+1,j+1,0);
      hu_decon_g->SetBinContent(ch+1,j+1,0);
    }
  }
  for (auto it = vplane_map.begin(); it!=vplane_map.end(); it++){
    int ch = it->first;
    int start = it->second.first/nrebin;
    int end = it->second.second/nrebin;
    for (int j=start; j!=end+1;j++){
      hv_decon->SetBinContent(ch+1,j+1,0);
      hv_decon_g->SetBinContent(ch+1,j+1,0);
    }
    for (int j=it->second.first;j!=it->second.second+1;j++){
      hv_raw->SetBinContent(ch+1,j+1,0);
    }
  }
  for (auto it = wplane_map.begin(); it!=wplane_map.end(); it++){
    int ch = it->first;
    int start = it->second.first/nrebin;
    int end = it->second.second/nrebin;
    for (int j=start; j!=end+1;j++){
      hw_decon->SetBinContent(ch+1,j+1,0);
      hw_decon_g->SetBinContent(ch+1,j+1,0);
    }
  }

  



  

  int tpc_status = WCP2dToy::Noisy_Event_ID(hu_decon, hv_decon, hw_decon, uplane_rms, vplane_rms, wplane_rms, uplane_map, vplane_map, wplane_map, hu_decon_g, hv_decon_g, hw_decon_g, nrebin, hv_raw, true);
  WCP2dToy::Organize_Dead_Channels(uplane_map, vplane_map, wplane_map, hv_raw->GetNbinsY()-1,nrebin);

    // loop through U/V/W plane to disable the bad channels completely
  for (auto it = uplane_map.begin(); it!=uplane_map.end(); it++){
    int ch = it->first;
    int start = it->second.first/nrebin;
    int end = it->second.second/nrebin;
    for (int j=start; j!=end+1;j++){
      hu_decon->SetBinContent(ch+1,j+1,0);
      hu_decon_g->SetBinContent(ch+1,j+1,0);
    }
  }
  for (auto it = vplane_map.begin(); it!=vplane_map.end(); it++){
    int ch = it->first;
    int start = it->second.first/nrebin;
    int end = it->second.second/nrebin;
    for (int j=start; j!=end+1;j++){
      hv_decon->SetBinContent(ch+1,j+1,0);
      hv_decon_g->SetBinContent(ch+1,j+1,0);
    }
    for (int j=it->second.first;j!=it->second.second+1;j++){
      hv_raw->SetBinContent(ch+1,j+1,0);
    }
  }
  for (auto it = wplane_map.begin(); it!=wplane_map.end(); it++){
    int ch = it->first;
    int start = it->second.first/nrebin;
    int end = it->second.second/nrebin;
    for (int j=start; j!=end+1;j++){
      hw_decon->SetBinContent(ch+1,j+1,0);
      hw_decon_g->SetBinContent(ch+1,j+1,0);
    }
  }


  
  
  WCP2dToy::pdDataFDS roi_fds(gds,hu_decon,hv_decon,hw_decon,eve_num);
  roi_fds.jump(eve_num);
  
  
  
  WCP2dToy::uBooNE_L1SP l1sp(hv_raw,hv_decon,hv_decon_g,nrebin, 0.0);
  
  WCP2dToy::pdDataFDS roi_gaus_fds(gds,hu_decon_g,hv_decon_g,hw_decon_g,eve_num);
  roi_gaus_fds.jump(eve_num);

  WCP2dToy::uBooNEDataError error_fds(gds,hu_decon_g, hv_decon_g, hw_decon_g, eve_num, nrebin);
  error_fds.jump(eve_num);
  
  
  // WCPSst::ToyuBooNESliceDataSource sds(roi_fds,roi_gaus_fds,threshold_u, 
  // 					    threshold_v, threshold_w, 
  // 					    threshold_ug, 
  // 					    threshold_vg, threshold_wg, 
  // 					    nwire_u, 
  // 					    nwire_v, nwire_w,
  // 					    &uplane_rms, &vplane_rms, &wplane_rms); 

  WCPSst::uBooNESliceDataSource sds(roi_fds,roi_gaus_fds,error_fds,
					 threshold_u, threshold_v, threshold_w,
					 nwire_u, nwire_v, nwire_w,
					 &uplane_rms, &vplane_rms, &wplane_rms); 
  
  // sds.jump(100);
  // full_sds.jump(100);
  // WCP::Slice slice = sds.get();
  // WCP::Slice slice1 = full_sds.get();
  // WCP::Slice slice2 = full_sds.get_error();

  // WCP::Channel::Group group = slice.group();
  // WCP::Channel::Group group1 = slice1.group();
  // WCP::Channel::Group group2 = slice2.group();
  // for (int i=0;i!=group.size();i++){
  //   std::cout << group.at(i).second << " " << group1.at(i).second << " " << group2.at(i).second << std::endl;
  // }
  
  cout << em("begin tiling") << endl;

  
  int ncount = 0;
  int ncount1 = 0;
  int ncount2 = 0;

  int ncount_t = 0;
  

  // WCP2dToy::ToyTiling **toytiling = new WCP2dToy::ToyTiling*[2400];
  // WCP2dToy::BadTiling **badtiling = new WCP2dToy::BadTiling*[2400];
  // WCP2dToy::MergeToyTiling **mergetiling = new WCP2dToy::MergeToyTiling*[2400];
  // WCP2dToy::ToyMatrix **toymatrix = new WCP2dToy::ToyMatrix*[2400];
  WCP2dToy::LowmemTiling **lowmemtiling = new WCP2dToy::LowmemTiling*[2400];
  WCP2dToy::ChargeSolving **chargesolver = new WCP2dToy::ChargeSolving*[2400];
  
  WCP2dToy::WCPHolder WCholder;

  //add in cluster
  Slim3DClusterSet cluster_set, cluster_delset, cluster_set_save;
  Slim3DDeadClusterSet dead_cluster_set, dead_cluster_delset, dead_cluster_set_save;
  
  int ncount_mcell = 0;

 
  int start_num = 0 ;
  int end_num = sds.size()-1;

  
  // start_num = 1925;
  // end_num = 1925;

  // start_num = 1665;
  // end_num = 1665;

  // start_num = 50;
  // end_num = 100;
 





  TFile *file = new TFile(Form("result_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");

  if(T_true!=0 && save_file!=0){
    T_true->CloneTree();
  }
  if(T_track!=0 && save_file!=0){
    T_track->CloneTree();
  }
  file0->Close();
  file0->Delete();

  if (T_op!=0){
    T_op->CloneTree();
  }
  
  Int_t n_cells;
  Int_t n_good_wires;
  Int_t n_bad_wires;
  Int_t time_slice;
  Int_t n_single_cells;
  Int_t ndirect_solved;
  Int_t nL1_solved;
  Int_t total_ndf;
  Double_t total_chi2;
  
  Int_t L1_ndf[1000];
  Int_t direct_ndf[1000];
  Double_t L1_chi2_base[1000];
  Double_t L1_chi2_penalty[1000];
  Double_t direct_chi2[1000];
  TTree *Twc;
  if (save_file==1){
    Twc = new TTree("Twc","Twc");
    Twc->SetDirectory(file);

    Twc->Branch("time_slice",&time_slice,"time_slice/I");
    Twc->Branch("n_cells",&n_cells,"n_cells/I");
    Twc->Branch("n_good_wires",&n_good_wires,"n_good_wires/I");  
    Twc->Branch("n_bad_wires",&n_bad_wires,"n_bad_wires/I");
    Twc->Branch("n_single_cells",&n_single_cells,"n_single_cells/I");
    
    Twc->Branch("ndirect_solved",&ndirect_solved,"ndirect_solved/I");
    Twc->Branch("nL1_solved",&nL1_solved,"nL1_solved/I");
    Twc->Branch("total_ndf",&total_ndf,"total_ndf/I");
    Twc->Branch("total_chi2",&total_chi2,"total_chi2/D");
    
    Twc->Branch("L1_ndf",L1_ndf,"L1_ndf[nL1_solved]/I");
    Twc->Branch("L1_chi2_base",L1_chi2_base,"L1_chi2_base[nL1_solved]/D");
    Twc->Branch("L1_chi2_penalty",L1_chi2_penalty,"L1_chi2_penalty[nL1_solved]/D");
    Twc->Branch("direct_ndf",direct_ndf,"direct_ndf[ndirect_solved]/I");
    Twc->Branch("direct_chi2",direct_chi2,"direct_chi2[ndirect_solved]/D");
  }
  
  //test 
  // uplane_map.begin()->second.second=5000;

  for (int i=start_num;i!=end_num+1;i++){
    if (i%200==0)
      std::cout << "Tiling: " << i << std::endl;

    sds.jump(i);
    WCP::Slice& slice = sds.get();
    WCP::Slice& slice_err = sds.get_error();
    
    lowmemtiling[i] = new WCP2dToy::LowmemTiling(i,nrebin,gds,WCholder);
    if (i==start_num){
      lowmemtiling[i]->init_bad_cells(uplane_map,vplane_map,wplane_map);
    }else{
      lowmemtiling[i]->check_bad_cells(lowmemtiling[i-1],uplane_map,vplane_map,wplane_map);
    }
    lowmemtiling[i]->init_good_cells(slice,slice_err,uplane_rms,vplane_rms,wplane_rms);

    //  std::cout << lowmemtiling[i]->get_two_bad_wire_cells().size() << std::endl;
    if (flag_l1){
      GeomWireSelection wires = lowmemtiling[i]->find_L1SP_wires();
      l1sp.AddWires(i,wires);
    }
    
  }

  // for (int i=start_num;i!=end_num+1;i++){
  //   delete lowmemtiling[i] ;
  // }

  cout << em("finish initial tiling") << endl;

  if (flag_l1){
    l1sp.AddWireTime_Raw();
    l1sp.Form_rois(6);
    
    roi_fds.refresh(hu_decon,hv_decon,hw_decon,eve_num);
    roi_gaus_fds.refresh(hu_decon_g,hv_decon_g,hw_decon_g,eve_num);
    error_fds.refresh(hu_decon_g, hv_decon_g, hw_decon_g, eve_num);
    
    //    cout << em("finish L1SP ") << endl;
    
    std::set<int>& time_slice_set = l1sp.get_time_slice_set();
    for (auto it = time_slice_set.begin(); it!= time_slice_set.end(); it++){
      int time_slice = *it;
      //std::cout << time_slice << std::endl;
      if (time_slice >= start_num && time_slice <=end_num){
	
	sds.jump(time_slice);
	WCP::Slice& slice = sds.get();
	WCP::Slice& slice_err = sds.get_error();
	
	// std::cout << lowmemtiling[time_slice]->get_wire_charge_error_map().size() << std::endl;
	// WCP::WireChargeMap& wire_charge_err_map = lowmemtiling[time_slice]->get_wire_charge_error_map();
	// for (auto it1= wire_charge_err_map.begin(); it1 != wire_charge_err_map.end(); it1++){
	// 	if ((*it1).second==0)
	// 	  std::cout << "A: " << (*it1).second << std::endl;
	// }
	//	lowmemtiling[time_slice]->reset_cells();

	delete lowmemtiling[time_slice];
	lowmemtiling[time_slice] = new WCP2dToy::LowmemTiling(time_slice,nrebin,gds,WCholder);
	if (time_slice==start_num){
	  lowmemtiling[time_slice]->init_bad_cells(uplane_map,vplane_map,wplane_map);
	}else{
	  lowmemtiling[time_slice]->check_bad_cells(lowmemtiling[time_slice-1],uplane_map,vplane_map,wplane_map);
	}
	lowmemtiling[time_slice]->init_good_cells(slice,slice_err,uplane_rms,vplane_rms,wplane_rms);
	
	//std::cout << lowmemtiling[time_slice]->get_two_bad_wire_cells().size() << std::endl;
	// std::cout << lowmemtiling[time_slice]->get_wire_charge_error_map().size() << std::endl;
	// {
	// 	WCP::WireChargeMap& wire_charge_err_map = lowmemtiling[time_slice]->get_wire_charge_error_map();
	// 	for (auto it1= wire_charge_err_map.begin(); it1 != wire_charge_err_map.end(); it1++){
	// 	  if ((*it1).second==0)
	// 	    std::cout << "B: " << (*it1).second << std::endl;
	// 	} 
	// }
      }
    }
    cout << em("finish L1SP and retiling") << endl;
  }
  
  delete hu_decon;
  delete hv_decon;
  delete hw_decon;
  delete hu_decon_g;
  delete hv_decon_g;
  delete hw_decon_g;
  delete hv_raw;
  delete hu_threshold;
  delete hv_threshold;
  delete hw_threshold;
  

  
  // to save original charge info (after L1SP) into output Trun
  std::vector<int> *timesliceId = new std::vector<int>;
  std::vector<std::vector<int>> *timesliceChannel = new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *raw_charge = new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *raw_charge_err = new std::vector<std::vector<int>>;
  //TClonesArray *raw_charge = new TClonesArray("TH1I");
  //TH1::AddDirectory(kFALSE);
  int raw_charge_ind = 0;
  for (int i=start_num; i<=end_num; i++){
      sds.jump(i);
      //WCP::Slice& slice = sds.get();
      //WCP::Slice& slice_err = sds.get_error();
      WCP::Channel::Group group = sds.get().group();
      //cout<<"slice charge group size: "<<group.size()<<endl;
      if(group.size()!=0){
        WCP::Channel::Group group_err = sds.get_error().group();
        //TH1I *htemp = new ( (*raw_charge)[raw_charge_ind] ) TH1I("","", 8256, 0, 8256);
        timesliceId->push_back(i);
        std::vector<int> channel;
        std::vector<int> charge;
        std::vector<int> charge_err;
        for(int j=0; j<group.size(); j++){
            //int channel = group.at(j).first;
            //float charge = group.at(j).second;
            //float charge_err = group_err.at(j).second;
                //cout <<"slice: "<<i<<" channel: "<<channel<<" charge: "<<charge<<" error: "<<charge_err<<endl;
            if (group.at(j).first >= 8256) { cout << "ERROR: out of bound of MicroBooNE channels." << endl;}
            //htemp->SetBinContent(channel+1, charge);    
            //htemp->SetBinError(channel+1, charge_err);
            channel.push_back(group.at(j).first);
            charge.push_back(group.at(j).second);
            charge_err.push_back(group_err.at(j).second);
        }
        //raw_charge_ind ++;
        timesliceChannel->push_back(channel);
        raw_charge->push_back(charge);
        raw_charge_err->push_back(charge_err);
      }
  }

  TTree *Trun = Trun0->CloneTree(0);
  Trun->SetTitle("Trun");
  Trun->SetName("Trun");
  Trun->SetDirectory(file);
  Trun0->GetEntry(eve_num);

  /* int detector = 0; // MicroBooNE */ 
  /* Trun->Branch("triggerBits",&triggerbits,"triggerBits/i"); */
  /* Trun->Branch("triggerTime",&triggerTime, "triggerTime/D"); */
  /* Trun->Branch("detector",&detector,"detector/I"); */

  if(argv[4][0]!='-') event_no = atoi(argv[4]);

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
  Trun->Branch("tpc_status",&tpc_status,"tpc_status/I");
  
  Trun->Branch("timesliceId",&timesliceId);
  //Trun->Branch("raw_charge",&raw_charge, 256000, 0);
  Trun->Branch("timesliceChannel",&timesliceChannel);
  Trun->Branch("raw_charge",&raw_charge);
  Trun->Branch("raw_charge_err",&raw_charge_err);

  Trun->Fill();
  delete raw_charge;
  delete raw_charge_err;
  delete timesliceId;
  delete timesliceChannel;

  file1->Close();
  delete file1;
 

  for (int i=start_num;i!=end_num+1;i++){
    
    // form clusters
    WCP::GeomCellMap& cell_wires_map = lowmemtiling[i]->get_cell_wires_map();
    GeomCellSelection allmcell;
    for (auto it=cell_wires_map.begin(); it!= cell_wires_map.end(); it++){
      SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)it->first;
      allmcell.push_back(mcell);
    }
    
    if (cluster_set.empty()){
      // if cluster is empty, just insert all the mcell, each as a cluster
      for (int j=0;j!=allmcell.size();j++){
	Slim3DCluster *cluster = new Slim3DCluster(*((SlimMergeGeomCell*)allmcell[j]));
	cluster_set.insert(cluster);
      }
    }else{
      for (int j=0;j!=allmcell.size();j++){
	int flag = 0;
	int flag_save = 0;
	Slim3DCluster *cluster_save = 0;
	cluster_delset.clear();

	// loop through merged cell
	for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
	  //loop through clusters
	  flag += (*it)->AddCell(*((SlimMergeGeomCell*)allmcell[j]),2);
	  if (flag==1 && flag != flag_save){
	    cluster_save = *it;
	  }else if (flag>1 && flag != flag_save){
	    cluster_save->MergeCluster(*(*it));
	    cluster_delset.insert(*it);
	  }
	  flag_save = flag;
	}
	
	for (auto it = cluster_delset.begin();it!=cluster_delset.end();it++){
	  cluster_set.erase(*it);
	  delete (*it);
	}
	
	if (flag==0){
	  Slim3DCluster *cluster = new Slim3DCluster(*((SlimMergeGeomCell*)allmcell[j]));
	  cluster_set.insert(cluster);
	}
      }

      // int ncount_mcell_cluster = 0;
      // for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
      // 	ncount_mcell_cluster += (*it)->get_allcell().size();
      // }
      // ncount_mcell += allmcell.size();
      // cout << i << " " << allmcell.size()  << " " << cluster_set.size()  << endl;
    }
    
    for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
      std::set<SlimMergeGeomCell*>&  abc = ((*it)->get_ordercell()).back();
      if (i - (*abc.begin())->GetTimeSlice()>1){
	cluster_set_save.insert(*it);
	cluster_set.erase(*it);
      }
    }
  }
  for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
    cluster_set_save.insert(*it);
  }
  cluster_set = cluster_set_save;
  cluster_set_save.clear();

 
  
  cout << em("finish initial clustering") << endl;

  //return 0;




  // first round only deal with the absolute match ... 
  std::map<Projected2DCluster*, std::vector<Slim3DCluster*>> u_2D_3D_clus_map;
  std::map<Projected2DCluster*, std::vector<Slim3DCluster*>> v_2D_3D_clus_map;
  std::map<Projected2DCluster*, std::vector<Slim3DCluster*>> w_2D_3D_clus_map;
  
  for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
    (*it)->Calc_Projection();
    // get U
    Projected2DCluster *u_2Dclus = (*it)->get_projection(WirePlaneType_t(0));
    if (u_2Dclus->get_number_time_slices() >0){
      bool flag_save = true;
      std::vector<Projected2DCluster*> to_be_removed;
      for (auto it1 = u_2D_3D_clus_map.begin(); it1!= u_2D_3D_clus_map.end(); it1++){
	Projected2DCluster *comp_2Dclus = it1->first;
	std::vector<Slim3DCluster*>& vec_3Dclus = it1->second;
	int comp_score = comp_2Dclus->judge_coverage(u_2Dclus);

	if (comp_score == 1){
	  // u_2Dclus is part of comp_2Dclus
	  flag_save = false;
	  break;
	}else if (comp_score == 2){
	  // u_2D_clus is the same as comp_2Dclus
	  flag_save = false;
	  vec_3Dclus.push_back((*it));
	  break;
	}else if (comp_score == -1){
	  // comp_2Dclus is part of u_2Dclus
	  to_be_removed.push_back(comp_2Dclus);
	}else if (comp_score == 0){
	  //they do not match ...
	  // do nothing ... 
	}
      }

      // remove the small stuff ...
      for (auto it1 = to_be_removed.begin(); it1!=to_be_removed.end(); it1++){
	u_2D_3D_clus_map.erase((*it1));
      }
      // save it 
      if (flag_save){
      	std::vector<Slim3DCluster*> vec_3Dclus;
      	vec_3Dclus.push_back((*it));
      	u_2D_3D_clus_map[u_2Dclus] = vec_3Dclus;
      }
      //std::cout << u_2D_3D_clus_map.size() << std::endl;
    }
  
    // get V
    Projected2DCluster *v_2Dclus = (*it)->get_projection(WirePlaneType_t(1));
    if (v_2Dclus->get_number_time_slices() >0){
      bool flag_save = true;
      std::vector<Projected2DCluster*> to_be_removed;
      for (auto it1 = v_2D_3D_clus_map.begin(); it1!= v_2D_3D_clus_map.end(); it1++){
  	Projected2DCluster *comp_2Dclus = it1->first;
  	std::vector<Slim3DCluster*>& vec_3Dclus = it1->second;
  	int comp_score = comp_2Dclus->judge_coverage(v_2Dclus);
  	if (comp_score == 1){
  	  // v_2Dclus is part of comp_2Dclus
  	  flag_save = false;
  	  break;
  	}else if (comp_score == 2){
  	  // v_2D_clus is the same as comp_2Dclus
  	  flag_save = false;
  	  vec_3Dclus.push_back((*it));
  	  break;
  	}else if (comp_score == -1){
  	  // comp_2Dclus is part of v_2Dclus
  	  to_be_removed.push_back(comp_2Dclus);
  	}else if (comp_score == 0){
  	  //they do not match ...
  	  // do nothing ... 
  	}
      }
      // remove the small stuff ...
      for (auto it1 = to_be_removed.begin(); it1!=to_be_removed.end(); it1++){
  	v_2D_3D_clus_map.erase((*it1));
      }
      // save it 
      if (flag_save){
      	std::vector<Slim3DCluster*> vec_3Dclus;
      	vec_3Dclus.push_back((*it));
      	v_2D_3D_clus_map[v_2Dclus] = vec_3Dclus;
      }
      //std::cout << v_2D_3D_clus_map.size() << std::endl;
    }
    
    // get W
    Projected2DCluster *w_2Dclus = (*it)->get_projection(WirePlaneType_t(2));
    if (w_2Dclus->get_number_time_slices() >0){
      bool flag_save = true;
      std::vector<Projected2DCluster*> to_be_removed;
      for (auto it1 = w_2D_3D_clus_map.begin(); it1!= w_2D_3D_clus_map.end(); it1++){
  	Projected2DCluster *comp_2Dclus = it1->first;
  	std::vector<Slim3DCluster*>& vec_3Dclus = it1->second;
	
  	int comp_score = comp_2Dclus->judge_coverage(w_2Dclus);
	if (comp_score == 1){
      	  // w_2Dclus is part of comp_2Dclus
      	  flag_save = false;
      	  break;
      	}else if (comp_score == 2){
      	  // w_2D_clus is the same as comp_2Dclus
      	  flag_save = false;
      	  vec_3Dclus.push_back((*it));
      	  break;
      	}else if (comp_score == -1){
      	  // comp_2Dclus is part of w_2Dclus
      	  to_be_removed.push_back(comp_2Dclus);
      	}else if (comp_score == 0){
      	  //they do not match ...
      	  // do nothing ... 
      	}
      }
      // remove the small stuff ...
      for (auto it1 = to_be_removed.begin(); it1!=to_be_removed.end(); it1++){
  	w_2D_3D_clus_map.erase((*it1));
      }
      // save it 
      if (flag_save){
      	std::vector<Slim3DCluster*> vec_3Dclus;
      	vec_3Dclus.push_back((*it));
      	w_2D_3D_clus_map[w_2Dclus] = vec_3Dclus;
      }
      //std::cout << w_2D_3D_clus_map.size() << std::endl;
    }
  }


  // TTree *T_2Dcluster = new TTree("T_2Dcluster","T_2Dcluster");
  // T_2Dcluster->SetDirectory(file);
  Int_t cluster1_ID, cluster2_ID;
  Int_t plane_no;
  Int_t cluster1_wire, cluster2_wire;
  Int_t cluster1_dead_wire, cluster2_dead_wire;
  Int_t common_wire;
  Float_t cluster1_charge, cluster2_charge;
  Float_t cluster1_charge_estimated, cluster2_charge_estimated;
  Float_t common_charge;
  Int_t value;
  // T_2Dcluster->Branch("cluster1_ID",&cluster1_ID,"cluster1_ID/I");
  // T_2Dcluster->Branch("cluster2_ID",&cluster2_ID,"cluster2_ID/I");
  // T_2Dcluster->Branch("plane_no",&plane_no,"plane_no/I");
  // T_2Dcluster->Branch("cluster1_wire",&cluster1_wire,"cluster1_wire/I");
  // T_2Dcluster->Branch("cluster2_wire",&cluster2_wire,"cluster2_wire/I");
  // T_2Dcluster->Branch("cluster1_dead_wire",&cluster1_dead_wire,"cluster1_dead_wire/I");
  // T_2Dcluster->Branch("cluster2_dead_wire",&cluster2_dead_wire,"cluster2_dead_wire/I");
  // T_2Dcluster->Branch("common_wire",&common_wire,"common_wire/I");
  // T_2Dcluster->Branch("cluster1_charge",&cluster1_charge,"cluster1_charge/F");
  // T_2Dcluster->Branch("cluster2_charge",&cluster2_charge,"cluster2_charge/F");
  // T_2Dcluster->Branch("cluster1_charge_estimated",&cluster1_charge_estimated,"cluster1_charge_estimated/F");
  // T_2Dcluster->Branch("cluster2_charge_estimated",&cluster2_charge_estimated,"cluster2_charge_estimated/F");
  // T_2Dcluster->Branch("common_charge",&common_charge,"common_charge/F");
  // T_2Dcluster->Branch("value",&value,"value/I");

  
  std::cout << cluster_set.size() << " " << u_2D_3D_clus_map.size() << " " << v_2D_3D_clus_map.size() << " " << w_2D_3D_clus_map.size() << std::endl;
  {
    std::vector<Projected2DCluster*> to_be_removed;
    for (auto it = u_2D_3D_clus_map.begin(); it!= u_2D_3D_clus_map.end(); it++){
      Projected2DCluster *u_2Dclus = it->first;
      if (find(to_be_removed.begin(), to_be_removed.end(), u_2Dclus) == to_be_removed.end()){
  	cluster2_ID = u_2Dclus->get_parent_cluster_id();
  	plane_no = 0;
  	auto it1 = it; it1++;
  	for (auto it2 = it1; it2!= u_2D_3D_clus_map.end(); it2++){
  	  Projected2DCluster *comp_2Dclus = it2->first;
	  
  	  std::vector<int> comp_results = comp_2Dclus->calc_coverage(u_2Dclus);
  	  cluster1_ID = comp_2Dclus->get_parent_cluster_id();
  	  cluster1_wire = comp_results.at(0);
  	  cluster2_wire = comp_results.at(1);
  	  cluster1_dead_wire = comp_results.at(2);
  	  cluster2_dead_wire = comp_results.at(3);
  	  common_wire = comp_results.at(4);
  	  cluster1_charge = comp_results.at(5);
  	  cluster2_charge = comp_results.at(6);
  	  cluster1_charge_estimated = comp_results.at(7);
  	  cluster2_charge_estimated = comp_results.at(8);
  	  common_charge = comp_results.at(9);
  	  value = comp_2Dclus->judge_coverage_alt(u_2Dclus);
  	  //T_2Dcluster->Fill();

  	  if (value==1){
  	    to_be_removed.push_back(u_2Dclus);
  	  }else if (value==-1){
  	    to_be_removed.push_back(comp_2Dclus);
  	  }
  	}
      }
    }
    for (auto it = to_be_removed.begin(); it!= to_be_removed.end(); it++){
      u_2D_3D_clus_map.erase((*it));
    }
  }

  {
    std::vector<Projected2DCluster*> to_be_removed;
    for (auto it = v_2D_3D_clus_map.begin(); it!= v_2D_3D_clus_map.end(); it++){
      Projected2DCluster *v_2Dclus = it->first;
      if (find(to_be_removed.begin(), to_be_removed.end(), v_2Dclus) == to_be_removed.end()){
  	cluster2_ID = v_2Dclus->get_parent_cluster_id();
  	plane_no = 1;
  	auto it1 = it; it1++;
  	for (auto it2 = it1; it2!= v_2D_3D_clus_map.end(); it2++){
  	  Projected2DCluster *comp_2Dclus = it2->first;
  	  std::vector<int> comp_results = comp_2Dclus->calc_coverage(v_2Dclus);
	  
  	  cluster1_ID = comp_2Dclus->get_parent_cluster_id();
  	  cluster1_wire = comp_results.at(0);
  	  cluster2_wire = comp_results.at(1);
  	  cluster1_dead_wire = comp_results.at(2);
  	  cluster2_dead_wire = comp_results.at(3);
  	  common_wire = comp_results.at(4);
  	  cluster1_charge = comp_results.at(5);
  	  cluster2_charge = comp_results.at(6);
  	  cluster1_charge_estimated = comp_results.at(7);
  	  cluster2_charge_estimated = comp_results.at(8);
  	  common_charge = comp_results.at(9);
  	  value = comp_2Dclus->judge_coverage_alt(v_2Dclus);
  	  //T_2Dcluster->Fill();

  	  if (value==1){
  	    to_be_removed.push_back(v_2Dclus);
  	  }else if (value==-1){
  	    to_be_removed.push_back(comp_2Dclus);
  	  }
  	}
      }
    }
    for (auto it = to_be_removed.begin(); it!= to_be_removed.end(); it++){
      v_2D_3D_clus_map.erase((*it));
    }
  }

  
  {
    std::vector<Projected2DCluster*> to_be_removed;
    for (auto it = w_2D_3D_clus_map.begin(); it!= w_2D_3D_clus_map.end(); it++){
      Projected2DCluster *w_2Dclus = it->first;
      if (find(to_be_removed.begin(), to_be_removed.end(), w_2Dclus) == to_be_removed.end()){
  	cluster2_ID = w_2Dclus->get_parent_cluster_id();
  	plane_no = 2;
  	auto it1 = it; it1++;
  	for (auto it2 = it1; it2!= w_2D_3D_clus_map.end(); it2++){
  	  Projected2DCluster *comp_2Dclus = it2->first;
  	  std::vector<int> comp_results = comp_2Dclus->calc_coverage(w_2Dclus);
  	  cluster1_ID = comp_2Dclus->get_parent_cluster_id();
  	  cluster1_wire = comp_results.at(0);
  	  cluster2_wire = comp_results.at(1);
  	  cluster1_dead_wire = comp_results.at(2);
  	  cluster2_dead_wire = comp_results.at(3);
  	  common_wire = comp_results.at(4);
  	  cluster1_charge = comp_results.at(5);
  	  cluster2_charge = comp_results.at(6);
  	  cluster1_charge_estimated = comp_results.at(7);
  	  cluster2_charge_estimated = comp_results.at(8);
  	  common_charge = comp_results.at(9);
  	  value = comp_2Dclus->judge_coverage_alt(w_2Dclus);
  	  //T_2Dcluster->Fill();

  	  if (value==1){
  	    to_be_removed.push_back(w_2Dclus);
  	  }else if (value==-1){
  	    to_be_removed.push_back(comp_2Dclus);
  	  }
  	}
      }
    }
    for (auto it = to_be_removed.begin(); it!= to_be_removed.end(); it++){
      w_2D_3D_clus_map.erase((*it));
    }
  }

  std::cout << cluster_set.size() << " " << u_2D_3D_clus_map.size() << " " << v_2D_3D_clus_map.size() << " " << w_2D_3D_clus_map.size() << std::endl;
  

  //label the cluster ...
  for (auto it = u_2D_3D_clus_map.begin(); it!= u_2D_3D_clus_map.end(); it++){
    // if (it->second.size()>1)
    //   std::cout << "U: " << it->second.size() << std::endl;
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      (*it1)->set_flag_saved((*it1)->get_flag_saved()+1);
      // (*it1)->set_flag_saved(1);
    }
  }
  for (auto it = v_2D_3D_clus_map.begin(); it!= v_2D_3D_clus_map.end(); it++){
    // if (it->second.size()>1)
    //   std::cout << "V: " << it->second.size() << std::endl;
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      (*it1)->set_flag_saved((*it1)->get_flag_saved()+1);
      //(*it1)->set_flag_saved(1);
    }
  }
  for (auto it = w_2D_3D_clus_map.begin(); it!= w_2D_3D_clus_map.end(); it++){
    // if (it->second.size()>1)
    //   std::cout << "W: " << it->second.size() << std::endl;
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      (*it1)->set_flag_saved((*it1)->get_flag_saved()+1);
      //(*it1)->set_flag_saved(1);
    }
  }
  
  // rescan it
  for (auto it = u_2D_3D_clus_map.begin(); it!= u_2D_3D_clus_map.end(); it++){
    int max_flag_saved = -1;
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      if ((*it1)->get_flag_saved()>max_flag_saved)
  	max_flag_saved = (*it1)->get_flag_saved();
    }
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      if ((*it1)->get_flag_saved()!=max_flag_saved)
  	(*it1)->set_flag_saved_1((*it1)->get_flag_saved_1()+1);
    }
  }
  for (auto it = v_2D_3D_clus_map.begin(); it!= v_2D_3D_clus_map.end(); it++){
    int max_flag_saved = -1;
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      if ((*it1)->get_flag_saved()>max_flag_saved)
  	max_flag_saved = (*it1)->get_flag_saved();
    }
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      if ((*it1)->get_flag_saved()!=max_flag_saved)
  	(*it1)->set_flag_saved_1((*it1)->get_flag_saved_1()+1);
    }
  }
  for (auto it = w_2D_3D_clus_map.begin(); it!= w_2D_3D_clus_map.end(); it++){
    int max_flag_saved = -1;
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      if ((*it1)->get_flag_saved()>max_flag_saved)
  	max_flag_saved = (*it1)->get_flag_saved();
    }
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      if ((*it1)->get_flag_saved()!=max_flag_saved)
  	(*it1)->set_flag_saved_1((*it1)->get_flag_saved_1()+1);
    }
  }

  TTree *T_3Dcluster;
  Int_t cluster_id;
  Int_t saved;
  Float_t total_charge;
  Float_t min_charge;
  Int_t flag_saved;
  Int_t flag_saved_1;
  Int_t n_mcells;
  Int_t n_timeslices;
  
  if (save_file ==1){
    T_3Dcluster = new TTree("T_3Dcluster","T_3Dcluster");
    T_3Dcluster->Branch("cluster_id",&cluster_id,"cluster_id/I");
    T_3Dcluster->Branch("saved",&saved,"saved/I");
    T_3Dcluster->Branch("total_charge",&total_charge,"total_charge/F");
    T_3Dcluster->Branch("min_charge",&min_charge,"min_charge/F");
    T_3Dcluster->Branch("flag_saved",&flag_saved,"flag_saved/I");
    T_3Dcluster->Branch("flag_saved_1",&flag_saved_1,"flag_saved_1/I");
    T_3Dcluster->Branch("n_mcells",&n_mcells,"n_mcells/I");
    T_3Dcluster->Branch("n_timeslices",&n_timeslices,"n_timeslices/I");
    T_3Dcluster->SetDirectory(file);
  }
  
  // remove the mcell from tiling ...
  int ncluster_saved = 0;
  int ncluster_deleted = 0;
  int nmcell_saved = 0;
  int nmcell_deleted = 0;
  int nmcell_before = 0, nmcell_after = 0;
  for (int i=start_num;i!=end_num+1;i++){
    nmcell_before += lowmemtiling[i]->get_cell_wires_map().size();
  }
  
  for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
    GeomCellSelection mcells =(*it)->get_allcell();
    int num = 0;
    if ((*it)->get_projection(WirePlaneType_t(0))->get_number_time_slices()!=0) num++;
    if ((*it)->get_projection(WirePlaneType_t(1))->get_number_time_slices()!=0) num++;
    if ((*it)->get_projection(WirePlaneType_t(2))->get_number_time_slices()!=0) num++;
    
    cluster_id = (*it)->get_id();
    total_charge = (*it)->get_total_charge();
    min_charge = (*it)->get_min_total_charge();
    flag_saved = (*it)->get_flag_saved();
    flag_saved_1 = (*it)->get_flag_saved_1();
    n_mcells = (*it)->get_allcell().size();
    n_timeslices = (*it)->get_ordercell().size();
    
    if ((*it)->get_flag_saved()-(*it)->get_flag_saved_1() ==3){
      // look at each cell level ...
      if ( sqrt(pow(n_timeslices/3.,2) + pow(min_charge/n_mcells/3000.,2))<1 || min_charge/n_mcells/2000.<1.){
	saved = 0;
      }else{
	saved = 1;
      }
    }else if ((*it)->get_flag_saved()-(*it)->get_flag_saved_1()  ==2){
      if ( sqrt(pow(n_timeslices/8.,2) + pow(min_charge/n_mcells/8000.,2))<1 ||  min_charge/n_mcells/4000.<1.){
	saved = 0;
      }else{
	saved = 1;
      }
    }else if ((*it)->get_flag_saved()-(*it)->get_flag_saved_1()  ==1){
      if ( sqrt(pow(n_timeslices/8.,2) + pow(min_charge/n_mcells/8000.,2))<1 || min_charge/n_mcells/6000.<1.){
	saved = 0;
      }else{
	saved = 1;
      }
    }else{
      saved = 0;
    }

    // test
    //    saved = 1;
    
    // if (min_charge/n_mcells < 5000) saved = 0;

    if (saved==1){
      ncluster_saved ++;
      nmcell_saved += mcells.size();
    }else{
      ncluster_deleted ++;
      nmcell_deleted += mcells.size();
      // remove them ...
      for (auto it1 = mcells.begin(); it1!=mcells.end(); it1++){
	SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)(*it1);
	lowmemtiling[mcell->GetTimeSlice()]->Erase_Cell(mcell);
      }
    }
    
    //T_3Dcluster->Fill();
  }

  for (int i=start_num;i!=end_num+1;i++){
    nmcell_after += lowmemtiling[i]->get_cell_wires_map().size();
  }
  
  std::cout << ncluster_saved << " " << nmcell_saved << " "
	    << ncluster_deleted << " " << nmcell_deleted << " "
	    << nmcell_before << " " << nmcell_after << " "
	    << std::endl;
  
  cout << em("finish 1st round of deghosting") << endl;
  
  
  for (int i=start_num;i!=end_num+1;i++){
    if (i%400==0)
      std::cout << "1st Solving: " << i << std::endl;

    // tiling after the firs round of deghosting ... 
    lowmemtiling[i]->MergeWires();
       
    // create individual cells ...
    //GeomCellSelection single_cells = lowmemtiling[i]->create_single_cells();
    time_slice = i;
    n_cells = lowmemtiling[i]->get_cell_wires_map().size();//get_all_cell_centers().size();
    n_good_wires = lowmemtiling[i]->get_all_good_wires().size();
    n_bad_wires = lowmemtiling[i]->get_all_bad_wires().size();
    //n_single_cells = single_cells.size();
    // L1 solving
    chargesolver[i] = new WCP2dToy::ChargeSolving(gds, *lowmemtiling[i]);
    ndirect_solved = chargesolver[i]->get_ndirect_solved();
    nL1_solved = chargesolver[i]->get_nL1_solved();
    chargesolver[i]->Update_ndf_chi2();
    total_chi2 = chargesolver[i]->get_ndf();
    total_ndf = chargesolver[i]->get_chi2();
    for (Int_t k=0;k!=nL1_solved;k++){
      L1_ndf[k] = chargesolver[i]->get_L1_ndf(k);
      L1_chi2_base[k] = chargesolver[i]->get_L1_chi2_base(k);
      L1_chi2_penalty[k] = chargesolver[i]->get_L1_chi2_penalty(k);
    }
    for (Int_t k=0;k!=ndirect_solved;k++){
      direct_ndf[k] = chargesolver[i]->get_direct_ndf(k);
      direct_chi2[k] = chargesolver[i]->get_direct_chi2(k);
    }
    //Twc->Fill();
  }
  // cout << em("finish 1st round of solving") << endl;

  // Int_t nc_mcells = 0;
  // for (int i=start_num; i!=end_num+1;i++){
  //   GeomCellMap cell_wires_map = lowmemtiling[i]->get_cell_wires_map();
  //   for (auto it = cell_wires_map.begin(); it!= cell_wires_map.end(); it++){
  //     SlimMergeGeomCell *mcell = (SlimMergeGeomCell*) it->first;
  //      if (chargesolver[i]->get_mcell_charge(mcell)>300){
  // 	 nc_mcells ++;
  //      }
  //   }
  // }
  // std::cout << nc_mcells << std::endl;
  // label merge cells according to connectivities, going through clusters ..
  std::map<const GeomCell*, GeomCellSelection> front_cell_map;
  std::map<const GeomCell*, GeomCellSelection> back_cell_map;
  for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
    (*it)->Form_maps(2, front_cell_map,back_cell_map);
  }
  // std::cout << front_cell_map.size() << " " << back_cell_map.size() << std::endl;

  // repeat solving by changing the weight // getting weight ...

  for (auto it = front_cell_map.begin(); it!=front_cell_map.end();it++){
    SlimMergeGeomCell *mcell1 = (SlimMergeGeomCell*) it->first;
    int time_slice1 = mcell1->GetTimeSlice();
    bool flag1 = chargesolver[time_slice1]->get_mcell_charge(mcell1)>300;
    //std::cout << flag1 << std::endl;
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      SlimMergeGeomCell *mcell2 = (SlimMergeGeomCell*)(*it1);
      int time_slice2 = mcell2->GetTimeSlice();
      bool flag2 = chargesolver[time_slice2]->get_mcell_charge(mcell2)>300;
      if (flag2)
	chargesolver[time_slice1]->add_front_connectivity(mcell1);
      if(flag1)
	chargesolver[time_slice2]->add_back_connectivity(mcell2);
    }
  }

  // cout << em("finish 1st solving") << endl;
  
  for (int i=start_num; i!=end_num+1;i++){
    if (i%400==0)
      std::cout << "1st Solving with connectivity: " << i << std::endl;
    chargesolver[i]->L1_resolve(9,3);
  }
  
  // cout << em("finish 1st solving with connectivities") << endl;
   
  std::set<SlimMergeGeomCell*> potential_good_mcells, good_mcells;
  for (auto it = front_cell_map.begin(); it!=front_cell_map.end();it++){
    SlimMergeGeomCell *mcell1 = (SlimMergeGeomCell*) it->first;
    int time_slice1 = mcell1->GetTimeSlice();
    bool flag1 = chargesolver[time_slice1]->get_mcell_charge(mcell1)>300;
    if (flag1)
      good_mcells.insert(mcell1);
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      SlimMergeGeomCell *mcell2 = (SlimMergeGeomCell*)(*it1);
      int time_slice2 = mcell2->GetTimeSlice();
      bool flag2 = chargesolver[time_slice2]->get_mcell_charge(mcell2)>300;
      if (flag1)
  	potential_good_mcells.insert(mcell2);
      if (flag2){
	potential_good_mcells.insert(mcell1);
	good_mcells.insert(mcell2);
      }
    }
  }
  for (int i=start_num;i!=end_num+1;i++){
    WCP::GeomCellMap& cell_map = lowmemtiling[i]->get_cell_wires_map();
    for (auto it = cell_map.begin(); it!=cell_map.end(); it++){
      SlimMergeGeomCell *mcell = (SlimMergeGeomCell*) it->first;
      bool flag1 = chargesolver[i]->get_mcell_charge(mcell)>300;
      if (flag1){
	good_mcells.insert(mcell);
	potential_good_mcells.insert(mcell);
      }
    }
  }
  
  nmcell_before = 0;
  for (int i=start_num;i!=end_num+1;i++){
    nmcell_before += lowmemtiling[i]->get_cell_wires_map().size();
  }
  
  // removel absolute can be removed ...
  // completely overlapped with the good three-wire-cells ... 
  std::set<SlimMergeGeomCell*> potential_bad_mcells;
  for (int i=start_num; i!=end_num+1;i++){
    GeomCellSelection mcells = lowmemtiling[i]->local_deghosting(potential_good_mcells,good_mcells,false);
    for (auto it = mcells.begin(); it!=mcells.end(); it++){
      potential_bad_mcells.insert((SlimMergeGeomCell*)(*it));
    }
    //  if (i==1681){
    //   //draw ...
    //   sds.jump(i);
    //   WCP::Slice slice = sds.get();
    //   TApplication theApp("theApp",&argc,argv);
    //   theApp.SetReturnFromRun(true);
      
    //   TCanvas c1("ToyMC","ToyMC",800,600);
    //   c1.Draw();
      
    //   WCP2dToy::ToyEventDisplay display(c1, gds);
    //   display.charge_min = 0;
    //   display.charge_max = 5e4;
      
      
    //   gStyle->SetOptStat(0);
      
    //   const Int_t NRGBs = 5;
    //   const Int_t NCont = 255;
    //   Int_t MyPalette[NCont];
    //   Double_t stops[NRGBs] = {0.0, 0.34, 0.61, 0.84, 1.0};
    //   Double_t red[NRGBs] = {0.0, 0.0, 0.87 ,1.0, 0.51};
    //   Double_t green[NRGBs] = {0.0, 0.81, 1.0, 0.2 ,0.0};
    //   Double_t blue[NRGBs] = {0.51, 1.0, 0.12, 0.0, 0.0};
    //   Int_t FI = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    //   gStyle->SetNumberContours(NCont);
    //   for (int kk=0;kk!=NCont;kk++) MyPalette[kk] = FI+kk;
    //   gStyle->SetPalette(NCont,MyPalette);
      
    //   GeomCellSelection single_cells = lowmemtiling[i]->create_single_cells();
      
    //   display.init(0,10.3698,-2.33/2.,2.33/2.);
    //   display.draw_mc(1,WCP::PointValueVector(),"colz");
    //   display.draw_slice(slice,""); // draw wire 
    //   // display.draw_wires(vec1_wires.at(64),"same"); // draw wire 
    //   // // display.draw_bad_region(uplane_map,i,nrebin,0,"same");
    //   // // display.draw_bad_region(vplane_map,i,nrebin,1,"same");
    //   // // display.draw_bad_region(wplane_map,i,nrebin,2,"same");
    //   // display.draw_bad_cell(badtiling[i]->get_cell_all());
    //   display.draw_cells(single_cells,"*same");
    //   //display.draw_points(lowmemtiling[i]->get_all_cell_centers(),"*");
    //   //display.draw_merged_wires(lowmemtiling[i]->get_all_good_wires(),"same",2);
    //   //display.draw_merged_wires(lowmemtiling[i]->get_all_bad_wires(),"same",1);
      
    //   //display.draw_mergecells(mergetiling[i]->get_allcell(),"*same",0); //0 is normal, 1 is only draw the ones containt the truth cell
      
    //   // display.draw_wires_charge(toytiling[i]->wcmap(),"Fsame",FI);
    //   // display.draw_cells_charge(toytiling[i]->get_allcell(),"Fsame");
    //   theApp.Run();
    // }
  }
  nmcell_after = 0;
  for (int i=start_num;i!=end_num+1;i++){
    nmcell_after += lowmemtiling[i]->get_cell_wires_map().size();
  }
  std::cout << nmcell_before << " " << nmcell_after << " " << potential_bad_mcells.size() << std::endl;
  
  
  //std::cout << good_mcells.size() << std::endl;
  // // see the difference
  // nc_mcells = 0;
  // for (int i=start_num; i!=end_num+1;i++){
  //   GeomCellMap cell_wires_map = lowmemtiling[i]->get_cell_wires_map();
  //   for (auto it = cell_wires_map.begin(); it!= cell_wires_map.end(); it++){
  //     SlimMergeGeomCell *mcell = (SlimMergeGeomCell*) it->first;
  //      if (chargesolver[i]->get_mcell_charge(mcell)>300){
  // 	 nc_mcells ++;
  //      }
  //   }
  // }
  // std::cout << nc_mcells << std::endl;
  cout << em("finish 2nd round of solving with connectivities") << endl;

  
  // delete clusters here ... 
  for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
    delete *it;
  }
  cluster_set.clear();
  cluster_delset.clear();
  u_2D_3D_clus_map.clear();
  v_2D_3D_clus_map.clear();
  w_2D_3D_clus_map.clear();

  {
    Slim3DClusterSet temp_cluster_set, temp_cluster_delset, temp_cluster_set_save;
    // 2nd round of clustering
    for (int i=start_num; i!=end_num+1;i++){
      // if (i%400==0)
      //   std::cout << "2nd Clustering: " << i << std::endl;
      // form clusters
      WCP::GeomCellMap& cell_wires_map = lowmemtiling[i]->get_cell_wires_map();
      GeomCellSelection allmcell;
      for (auto it=cell_wires_map.begin(); it!= cell_wires_map.end(); it++){
	SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)it->first;
	if (potential_good_mcells.find(mcell)!=potential_good_mcells.end() && potential_bad_mcells.find(mcell)==potential_bad_mcells.end())
	  allmcell.push_back(mcell);
      }
      if (temp_cluster_set.empty()){
	// if cluster is empty, just insert all the mcell, each as a cluster
	for (int j=0;j!=allmcell.size();j++){
	  Slim3DCluster *cluster = new Slim3DCluster(*((SlimMergeGeomCell*)allmcell[j]));
	  temp_cluster_set.insert(cluster);
	}
      }else{
	for (int j=0;j!=allmcell.size();j++){
	  int flag = 0;
	  int flag_save = 0;
	  Slim3DCluster *cluster_save = 0;
	  temp_cluster_delset.clear();
	  
	  // loop through merged cell
	  for (auto it = temp_cluster_set.begin();it!=temp_cluster_set.end();it++){
	    //loop through clusters
	    
	    flag += (*it)->AddCell(*((SlimMergeGeomCell*)allmcell[j]),2);
	    if (flag==1 && flag != flag_save){
	      cluster_save = *it;
	    }else if (flag>1 && flag != flag_save){
	      cluster_save->MergeCluster(*(*it));
	      temp_cluster_delset.insert(*it);
	    }
	    flag_save = flag;
	    
	  }
	  
	  for (auto it = temp_cluster_delset.begin();it!=temp_cluster_delset.end();it++){
	    temp_cluster_set.erase(*it);
	    delete (*it);
	  }
	  
	  if (flag==0){
	    Slim3DCluster *cluster = new Slim3DCluster(*((SlimMergeGeomCell*)allmcell[j]));
	    temp_cluster_set.insert(cluster);
	  }
	}
	
	// int ncount_mcell_cluster = 0;
	// for (auto it = temp_cluster_set.begin();it!=temp_cluster_set.end();it++){
	// 	ncount_mcell_cluster += (*it)->get_allcell().size();
	// }
	// ncount_mcell += allmcell.size();
	// cout << i << " " << allmcell.size()  << " " << temp_cluster_set.size()  << endl;
      }
      
      for (auto it = temp_cluster_set.begin();it!=temp_cluster_set.end();it++){
	std::set<SlimMergeGeomCell*>&  abc = ((*it)->get_ordercell()).back();
	if (i - (*abc.begin())->GetTimeSlice()>1){
	  temp_cluster_set_save.insert(*it);
	  temp_cluster_set.erase(*it);
	}
      }
    }
    for (auto it = temp_cluster_set.begin();it!=temp_cluster_set.end();it++){
      temp_cluster_set_save.insert(*it);
    }
    temp_cluster_set = temp_cluster_set_save;
    temp_cluster_set_save.clear();
    
    for (auto it = temp_cluster_set.begin(); it!= temp_cluster_set.end(); it++){
      cluster_set.insert(*it);
    }
  }

  {
    Slim3DClusterSet temp_cluster_set, temp_cluster_delset, temp_cluster_set_save;
    for (int i=start_num; i!=end_num+1;i++){
      //    if (i%400==0)
      // std::cout << "2nd Clustering (bad cells): " << i << std::endl;
      
      // form clusters
      WCP::GeomCellMap& cell_wires_map = lowmemtiling[i]->get_cell_wires_map();
      GeomCellSelection allmcell;
      for (auto it=cell_wires_map.begin(); it!= cell_wires_map.end(); it++){
	SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)it->first;
	if (potential_good_mcells.find(mcell)!=potential_good_mcells.end() && potential_bad_mcells.find(mcell)!=potential_bad_mcells.end())
	  allmcell.push_back(mcell);
      }
      if (temp_cluster_set.empty()){
	// if cluster is empty, just insert all the mcell, each as a cluster
	for (int j=0;j!=allmcell.size();j++){
	  Slim3DCluster *cluster = new Slim3DCluster(*((SlimMergeGeomCell*)allmcell[j]));
	  temp_cluster_set.insert(cluster);
	}
      }else{
	for (int j=0;j!=allmcell.size();j++){
	  int flag = 0;
	  int flag_save = 0;
	  Slim3DCluster *cluster_save = 0;
	  temp_cluster_delset.clear();
	  
	  // loop through merged cell
	  for (auto it = temp_cluster_set.begin();it!=temp_cluster_set.end();it++){
	    //loop through clusters
	    
	    flag += (*it)->AddCell(*((SlimMergeGeomCell*)allmcell[j]),2);
	    if (flag==1 && flag != flag_save){
	      cluster_save = *it;
	    }else if (flag>1 && flag != flag_save){
	      cluster_save->MergeCluster(*(*it));
	      temp_cluster_delset.insert(*it);
	    }
	    flag_save = flag;
	    
	  }
	  
	  for (auto it = temp_cluster_delset.begin();it!=temp_cluster_delset.end();it++){
	    temp_cluster_set.erase(*it);
	    delete (*it);
	  }
	  
	  if (flag==0){
	    Slim3DCluster *cluster = new Slim3DCluster(*((SlimMergeGeomCell*)allmcell[j]));
	    temp_cluster_set.insert(cluster);
	  }
	}
	
	// int ncount_mcell_cluster = 0;
	// for (auto it = temp_cluster_set.begin();it!=temp_cluster_set.end();it++){
	// 	ncount_mcell_cluster += (*it)->get_allcell().size();
	// }
	// ncount_mcell += allmcell.size();
	// cout << i << " " << allmcell.size()  << " " << temp_cluster_set.size()  << endl;
      }

      for (auto it = temp_cluster_set.begin();it!=temp_cluster_set.end();it++){
	std::set<SlimMergeGeomCell*>&  abc = ((*it)->get_ordercell()).back();
	if (i - (*abc.begin())->GetTimeSlice()>1){
	  temp_cluster_set_save.insert(*it);
	  temp_cluster_set.erase(*it);
	}
      }
    }

    for (auto it = temp_cluster_set.begin();it!=temp_cluster_set.end();it++){
      temp_cluster_set_save.insert(*it);
    }
    temp_cluster_set = temp_cluster_set_save;
    temp_cluster_set_save.clear();

    for (auto it = temp_cluster_set.begin(); it!= temp_cluster_set.end(); it++){
      cluster_set.insert(*it);
    }
  }

  {
    Slim3DClusterSet temp_cluster_set, temp_cluster_delset, temp_cluster_set_save;
    for (int i=start_num; i!=end_num+1;i++){
      //    if (i%400==0)
      // std::cout << "2nd Clustering (bad cells): " << i << std::endl;
      
      // form clusters
      WCP::GeomCellMap& cell_wires_map = lowmemtiling[i]->get_cell_wires_map();
      GeomCellSelection allmcell;
      for (auto it=cell_wires_map.begin(); it!= cell_wires_map.end(); it++){
	SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)it->first;
	if (potential_good_mcells.find(mcell)==potential_good_mcells.end() && potential_bad_mcells.find(mcell)==potential_bad_mcells.end())
	  allmcell.push_back(mcell);
      }
      if (temp_cluster_set.empty()){
	// if cluster is empty, just insert all the mcell, each as a cluster
	for (int j=0;j!=allmcell.size();j++){
	  Slim3DCluster *cluster = new Slim3DCluster(*((SlimMergeGeomCell*)allmcell[j]));
	  temp_cluster_set.insert(cluster);
	}
      }else{
	for (int j=0;j!=allmcell.size();j++){
	  int flag = 0;
	  int flag_save = 0;
	  Slim3DCluster *cluster_save = 0;
	  temp_cluster_delset.clear();
	  
	  // loop through merged cell
	  for (auto it = temp_cluster_set.begin();it!=temp_cluster_set.end();it++){
	    //loop through clusters
	    
	    flag += (*it)->AddCell(*((SlimMergeGeomCell*)allmcell[j]),2);
	    if (flag==1 && flag != flag_save){
	      cluster_save = *it;
	    }else if (flag>1 && flag != flag_save){
	      cluster_save->MergeCluster(*(*it));
	      temp_cluster_delset.insert(*it);
	    }
	    flag_save = flag;
	    
	  }
	  
	  for (auto it = temp_cluster_delset.begin();it!=temp_cluster_delset.end();it++){
	    temp_cluster_set.erase(*it);
	    delete (*it);
	  }
	  
	  if (flag==0){
	    Slim3DCluster *cluster = new Slim3DCluster(*((SlimMergeGeomCell*)allmcell[j]));
	    temp_cluster_set.insert(cluster);
	  }
	}

	// int ncount_mcell_cluster = 0;
	// for (auto it = temp_cluster_set.begin();it!=temp_cluster_set.end();it++){
	// 	ncount_mcell_cluster += (*it)->get_allcell().size();
	// }
	// ncount_mcell += allmcell.size();
	// cout << i << " " << allmcell.size()  << " " << temp_cluster_set.size()  << endl;
      }

      for (auto it = temp_cluster_set.begin();it!=temp_cluster_set.end();it++){
	std::set<SlimMergeGeomCell*>&  abc = ((*it)->get_ordercell()).back();
	if (i - (*abc.begin())->GetTimeSlice()>1){
	  temp_cluster_set_save.insert(*it);
	  temp_cluster_set.erase(*it);
	}
      }
    }

    for (auto it = temp_cluster_set.begin();it!=temp_cluster_set.end();it++){
      temp_cluster_set_save.insert(*it);
    }
    temp_cluster_set = temp_cluster_set_save;
    temp_cluster_set_save.clear();
    
    for (auto it = temp_cluster_set.begin(); it!= temp_cluster_set.end(); it++){
      cluster_set.insert(*it);
    }
  }

  {
    Slim3DClusterSet temp_cluster_set, temp_cluster_delset, temp_cluster_set_save;
    
    for (int i=start_num; i!=end_num+1;i++){
      //    if (i%400==0)
      // std::cout << "2nd Clustering (bad cells): " << i << std::endl;
      
      // form clusters
      WCP::GeomCellMap& cell_wires_map = lowmemtiling[i]->get_cell_wires_map();
      GeomCellSelection allmcell;
      for (auto it=cell_wires_map.begin(); it!= cell_wires_map.end(); it++){
	SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)it->first;
	if (potential_good_mcells.find(mcell)==potential_good_mcells.end() && potential_bad_mcells.find(mcell)!=potential_bad_mcells.end())
	  allmcell.push_back(mcell);
      }
      if (temp_cluster_set.empty()){
	// if cluster is empty, just insert all the mcell, each as a cluster
	for (int j=0;j!=allmcell.size();j++){
	  Slim3DCluster *cluster = new Slim3DCluster(*((SlimMergeGeomCell*)allmcell[j]));
	  temp_cluster_set.insert(cluster);
	}
      }else{
	for (int j=0;j!=allmcell.size();j++){
	  int flag = 0;
	  int flag_save = 0;
	  Slim3DCluster *cluster_save = 0;
	  temp_cluster_delset.clear();
	  
	  // loop through merged cell
	  for (auto it = temp_cluster_set.begin();it!=temp_cluster_set.end();it++){
	    //loop through clusters
	    
	    flag += (*it)->AddCell(*((SlimMergeGeomCell*)allmcell[j]),2);
	    if (flag==1 && flag != flag_save){
	      cluster_save = *it;
	    }else if (flag>1 && flag != flag_save){
	      cluster_save->MergeCluster(*(*it));
	      temp_cluster_delset.insert(*it);
	    }
	    flag_save = flag;
	    
	  }
	  
	  for (auto it = temp_cluster_delset.begin();it!=temp_cluster_delset.end();it++){
	    temp_cluster_set.erase(*it);
	    delete (*it);
	  }
	  
	  if (flag==0){
	    Slim3DCluster *cluster = new Slim3DCluster(*((SlimMergeGeomCell*)allmcell[j]));
	    temp_cluster_set.insert(cluster);
	  }
	}
	
	// int ncount_mcell_cluster = 0;
	// for (auto it = temp_cluster_set.begin();it!=temp_cluster_set.end();it++){
	// 	ncount_mcell_cluster += (*it)->get_allcell().size();
	// }
	// ncount_mcell += allmcell.size();
	// cout << i << " " << allmcell.size()  << " " << temp_cluster_set.size()  << endl;
      }


      for (auto it = temp_cluster_set.begin();it!=temp_cluster_set.end();it++){
	std::set<SlimMergeGeomCell*>&  abc = ((*it)->get_ordercell()).back();
	if (i - (*abc.begin())->GetTimeSlice()>1){
	  temp_cluster_set_save.insert(*it);
	  temp_cluster_set.erase(*it);
	}
      }
    }

    for (auto it = temp_cluster_set.begin();it!=temp_cluster_set.end();it++){
      temp_cluster_set_save.insert(*it);
    }
    temp_cluster_set = temp_cluster_set_save;
    temp_cluster_set_save.clear();

    
    for (auto it = temp_cluster_set.begin(); it!= temp_cluster_set.end(); it++){
      cluster_set.insert(*it);
    }
  }

  
  
  
 

  
  for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
    (*it)->Calc_Projection();
    // get U
    Projected2DCluster *u_2Dclus = (*it)->get_projection(WirePlaneType_t(0));
    if (u_2Dclus->get_number_time_slices() >0){
      bool flag_save = true;
      std::vector<Projected2DCluster*> to_be_removed;
      for (auto it1 = u_2D_3D_clus_map.begin(); it1!= u_2D_3D_clus_map.end(); it1++){
  	Projected2DCluster *comp_2Dclus = it1->first;
  	std::vector<Slim3DCluster*>& vec_3Dclus = it1->second;
  	int comp_score = comp_2Dclus->judge_coverage(u_2Dclus);

  	if (comp_score == 1){
  	  // u_2Dclus is part of comp_2Dclus
  	  flag_save = false;
  	  break;
  	}else if (comp_score == 2){
  	  // u_2D_clus is the same as comp_2Dclus
  	  flag_save = false;
  	  vec_3Dclus.push_back((*it));
  	  break;
  	}else if (comp_score == -1){
  	  // comp_2Dclus is part of u_2Dclus
  	  to_be_removed.push_back(comp_2Dclus);
  	}else if (comp_score == 0){
  	  //they do not match ...
  	  // do nothing ... 
  	}
      }

      // remove the small stuff ...
      for (auto it1 = to_be_removed.begin(); it1!=to_be_removed.end(); it1++){
  	u_2D_3D_clus_map.erase((*it1));
      }
      // save it 
      if (flag_save){
      	std::vector<Slim3DCluster*> vec_3Dclus;
      	vec_3Dclus.push_back((*it));
      	u_2D_3D_clus_map[u_2Dclus] = vec_3Dclus;
      }
      //std::cout << u_2D_3D_clus_map.size() << std::endl;
    }
  
    // get V
    Projected2DCluster *v_2Dclus = (*it)->get_projection(WirePlaneType_t(1));
    if (v_2Dclus->get_number_time_slices() >0){
      bool flag_save = true;
      std::vector<Projected2DCluster*> to_be_removed;
      for (auto it1 = v_2D_3D_clus_map.begin(); it1!= v_2D_3D_clus_map.end(); it1++){
  	Projected2DCluster *comp_2Dclus = it1->first;
  	std::vector<Slim3DCluster*>& vec_3Dclus = it1->second;
  	int comp_score = comp_2Dclus->judge_coverage(v_2Dclus);
  	if (comp_score == 1){
  	  // v_2Dclus is part of comp_2Dclus
  	  flag_save = false;
  	  break;
  	}else if (comp_score == 2){
  	  // v_2D_clus is the same as comp_2Dclus
  	  flag_save = false;
  	  vec_3Dclus.push_back((*it));
  	  break;
  	}else if (comp_score == -1){
  	  // comp_2Dclus is part of v_2Dclus
  	  to_be_removed.push_back(comp_2Dclus);
  	}else if (comp_score == 0){
  	  //they do not match ...
  	  // do nothing ... 
  	}
      }
      // remove the small stuff ...
      for (auto it1 = to_be_removed.begin(); it1!=to_be_removed.end(); it1++){
  	v_2D_3D_clus_map.erase((*it1));
      }
      // save it 
      if (flag_save){
      	std::vector<Slim3DCluster*> vec_3Dclus;
      	vec_3Dclus.push_back((*it));
      	v_2D_3D_clus_map[v_2Dclus] = vec_3Dclus;
      }
      //std::cout << v_2D_3D_clus_map.size() << std::endl;
    }
    
    // get W
    Projected2DCluster *w_2Dclus = (*it)->get_projection(WirePlaneType_t(2));
    if (w_2Dclus->get_number_time_slices() >0){
      bool flag_save = true;
      std::vector<Projected2DCluster*> to_be_removed;
      for (auto it1 = w_2D_3D_clus_map.begin(); it1!= w_2D_3D_clus_map.end(); it1++){
  	Projected2DCluster *comp_2Dclus = it1->first;
  	std::vector<Slim3DCluster*>& vec_3Dclus = it1->second;
	
  	int comp_score = comp_2Dclus->judge_coverage(w_2Dclus);
  	if (comp_score == 1){
      	  // w_2Dclus is part of comp_2Dclus
      	  flag_save = false;
      	  break;
      	}else if (comp_score == 2){
      	  // w_2D_clus is the same as comp_2Dclus
      	  flag_save = false;
      	  vec_3Dclus.push_back((*it));
      	  break;
      	}else if (comp_score == -1){
      	  // comp_2Dclus is part of w_2Dclus
      	  to_be_removed.push_back(comp_2Dclus);
      	}else if (comp_score == 0){
      	  //they do not match ...
      	  // do nothing ... 
      	}
      }
      // remove the small stuff ...
      for (auto it1 = to_be_removed.begin(); it1!=to_be_removed.end(); it1++){
  	w_2D_3D_clus_map.erase((*it1));
      }
      // save it 
      if (flag_save){
      	std::vector<Slim3DCluster*> vec_3Dclus;
      	vec_3Dclus.push_back((*it));
      	w_2D_3D_clus_map[w_2Dclus] = vec_3Dclus;
      }
      //std::cout << w_2D_3D_clus_map.size() << std::endl;
    }
  }


   {
    std::vector<Projected2DCluster*> to_be_removed;
    for (auto it = u_2D_3D_clus_map.begin(); it!= u_2D_3D_clus_map.end(); it++){
      Projected2DCluster *u_2Dclus = it->first;
      if (find(to_be_removed.begin(), to_be_removed.end(), u_2Dclus) == to_be_removed.end()){
  	cluster2_ID = u_2Dclus->get_parent_cluster_id();
  	plane_no = 0;
  	auto it1 = it; it1++;
  	for (auto it2 = it1; it2!= u_2D_3D_clus_map.end(); it2++){
  	  Projected2DCluster *comp_2Dclus = it2->first;
	  
  	  std::vector<int> comp_results = comp_2Dclus->calc_coverage(u_2Dclus);
  	  cluster1_ID = comp_2Dclus->get_parent_cluster_id();
  	  cluster1_wire = comp_results.at(0);
  	  cluster2_wire = comp_results.at(1);
  	  cluster1_dead_wire = comp_results.at(2);
  	  cluster2_dead_wire = comp_results.at(3);
  	  common_wire = comp_results.at(4);
  	  cluster1_charge = comp_results.at(5);
  	  cluster2_charge = comp_results.at(6);
  	  cluster1_charge_estimated = comp_results.at(7);
  	  cluster2_charge_estimated = comp_results.at(8);
  	  common_charge = comp_results.at(9);
  	  value = comp_2Dclus->judge_coverage_alt(u_2Dclus);
  	  //T_2Dcluster->Fill();

  	  if (value==1){
  	    to_be_removed.push_back(u_2Dclus);
  	  }else if (value==-1){
  	    to_be_removed.push_back(comp_2Dclus);
  	  }
  	}
      }
    }
    for (auto it = to_be_removed.begin(); it!= to_be_removed.end(); it++){
      u_2D_3D_clus_map.erase((*it));
    }
  }

  {
    std::vector<Projected2DCluster*> to_be_removed;
    for (auto it = v_2D_3D_clus_map.begin(); it!= v_2D_3D_clus_map.end(); it++){
      Projected2DCluster *v_2Dclus = it->first;
      if (find(to_be_removed.begin(), to_be_removed.end(), v_2Dclus) == to_be_removed.end()){
  	cluster2_ID = v_2Dclus->get_parent_cluster_id();
  	plane_no = 1;
  	auto it1 = it; it1++;
  	for (auto it2 = it1; it2!= v_2D_3D_clus_map.end(); it2++){
  	  Projected2DCluster *comp_2Dclus = it2->first;
  	  std::vector<int> comp_results = comp_2Dclus->calc_coverage(v_2Dclus);
	  
  	  cluster1_ID = comp_2Dclus->get_parent_cluster_id();
  	  cluster1_wire = comp_results.at(0);
  	  cluster2_wire = comp_results.at(1);
  	  cluster1_dead_wire = comp_results.at(2);
  	  cluster2_dead_wire = comp_results.at(3);
  	  common_wire = comp_results.at(4);
  	  cluster1_charge = comp_results.at(5);
  	  cluster2_charge = comp_results.at(6);
  	  cluster1_charge_estimated = comp_results.at(7);
  	  cluster2_charge_estimated = comp_results.at(8);
  	  common_charge = comp_results.at(9);
  	  value = comp_2Dclus->judge_coverage_alt(v_2Dclus);
  	  //T_2Dcluster->Fill();

  	  if (value==1){
  	    to_be_removed.push_back(v_2Dclus);
  	  }else if (value==-1){
  	    to_be_removed.push_back(comp_2Dclus);
  	  }
  	}
      }
    }
    for (auto it = to_be_removed.begin(); it!= to_be_removed.end(); it++){
      v_2D_3D_clus_map.erase((*it));
    }
  }

  
  {
    std::vector<Projected2DCluster*> to_be_removed;
    for (auto it = w_2D_3D_clus_map.begin(); it!= w_2D_3D_clus_map.end(); it++){
      Projected2DCluster *w_2Dclus = it->first;
      if (find(to_be_removed.begin(), to_be_removed.end(), w_2Dclus) == to_be_removed.end()){
  	cluster2_ID = w_2Dclus->get_parent_cluster_id();
  	plane_no = 2;
  	auto it1 = it; it1++;
  	for (auto it2 = it1; it2!= w_2D_3D_clus_map.end(); it2++){
  	  Projected2DCluster *comp_2Dclus = it2->first;
  	  std::vector<int> comp_results = comp_2Dclus->calc_coverage(w_2Dclus);
  	  cluster1_ID = comp_2Dclus->get_parent_cluster_id();
  	  cluster1_wire = comp_results.at(0);
  	  cluster2_wire = comp_results.at(1);
  	  cluster1_dead_wire = comp_results.at(2);
  	  cluster2_dead_wire = comp_results.at(3);
  	  common_wire = comp_results.at(4);
  	  cluster1_charge = comp_results.at(5);
  	  cluster2_charge = comp_results.at(6);
  	  cluster1_charge_estimated = comp_results.at(7);
  	  cluster2_charge_estimated = comp_results.at(8);
  	  common_charge = comp_results.at(9);
  	  value = comp_2Dclus->judge_coverage_alt(w_2Dclus);
  	  //T_2Dcluster->Fill();

  	  if (value==1){
  	    to_be_removed.push_back(w_2Dclus);
  	  }else if (value==-1){
  	    to_be_removed.push_back(comp_2Dclus);
  	  }
  	}
      }
    }
    for (auto it = to_be_removed.begin(); it!= to_be_removed.end(); it++){
      w_2D_3D_clus_map.erase((*it));
    }
  }

  
  std::cout << cluster_set.size() << " " << u_2D_3D_clus_map.size() << " " << v_2D_3D_clus_map.size() << " " << w_2D_3D_clus_map.size() << std::endl;

  //label the cluster ...
  for (auto it = u_2D_3D_clus_map.begin(); it!= u_2D_3D_clus_map.end(); it++){
    // if (it->second.size()>1)
    //   std::cout << "U: " << it->second.size() << std::endl;
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      (*it1)->set_flag_saved((*it1)->get_flag_saved()+1);
      // (*it1)->set_flag_saved(1);
    }
  }
  for (auto it = v_2D_3D_clus_map.begin(); it!= v_2D_3D_clus_map.end(); it++){
    // if (it->second.size()>1)
    //   std::cout << "V: " << it->second.size() << std::endl;
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      (*it1)->set_flag_saved((*it1)->get_flag_saved()+1);
      //(*it1)->set_flag_saved(1);
    }
  }
  for (auto it = w_2D_3D_clus_map.begin(); it!= w_2D_3D_clus_map.end(); it++){
    // if (it->second.size()>1)
    //   std::cout << "W: " << it->second.size() << std::endl;
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      (*it1)->set_flag_saved((*it1)->get_flag_saved()+1);
      //(*it1)->set_flag_saved(1);
    }
  }
  
  // rescan it
  for (auto it = u_2D_3D_clus_map.begin(); it!= u_2D_3D_clus_map.end(); it++){
    int max_flag_saved = -1;
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      if ((*it1)->get_flag_saved()>max_flag_saved)
  	max_flag_saved = (*it1)->get_flag_saved();
    }
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      if ((*it1)->get_flag_saved()!=max_flag_saved)
  	(*it1)->set_flag_saved_1((*it1)->get_flag_saved_1()+1);
    }
  }
  for (auto it = v_2D_3D_clus_map.begin(); it!= v_2D_3D_clus_map.end(); it++){
    int max_flag_saved = -1;
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      if ((*it1)->get_flag_saved()>max_flag_saved)
  	max_flag_saved = (*it1)->get_flag_saved();
    }
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      if ((*it1)->get_flag_saved()!=max_flag_saved)
  	(*it1)->set_flag_saved_1((*it1)->get_flag_saved_1()+1);
    }
  }
  for (auto it = w_2D_3D_clus_map.begin(); it!= w_2D_3D_clus_map.end(); it++){
    int max_flag_saved = -1;
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      if ((*it1)->get_flag_saved()>max_flag_saved)
  	max_flag_saved = (*it1)->get_flag_saved();
    }
    for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
      if ((*it1)->get_flag_saved()!=max_flag_saved)
  	(*it1)->set_flag_saved_1((*it1)->get_flag_saved_1()+1);
    }
  }


  ncluster_saved=0;
  nmcell_saved=0;
  ncluster_deleted=0;
  nmcell_deleted=0;
  nmcell_before=0;
  nmcell_after=0;
  
  for (int i=start_num;i!=end_num+1;i++){
    // re-establish map ...
    lowmemtiling[i]->re_establish_maps();
    nmcell_before += lowmemtiling[i]->get_cell_wires_map().size();
  }
  
  for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
    GeomCellSelection& mcells =(*it)->get_allcell();
    int num = 0;
    if ((*it)->get_projection(WirePlaneType_t(0))->get_number_time_slices()!=0) num++;
    if ((*it)->get_projection(WirePlaneType_t(1))->get_number_time_slices()!=0) num++;
    if ((*it)->get_projection(WirePlaneType_t(2))->get_number_time_slices()!=0) num++;
    
    cluster_id = (*it)->get_id();
    total_charge = (*it)->get_total_charge();
    min_charge = (*it)->get_min_total_charge();
    flag_saved = (*it)->get_flag_saved();
    flag_saved_1 = (*it)->get_flag_saved_1();
    n_mcells = (*it)->get_allcell().size();
    n_timeslices = (*it)->get_ordercell().size();
    
    if ((*it)->get_flag_saved()-(*it)->get_flag_saved_1() ==3){
      // look at each cell level ...
      if ( sqrt(pow(n_timeslices/3.,2) + pow(min_charge/n_mcells/3000.,2))<1 || min_charge/n_mcells/2000.<1.){
  	saved = 0;
      }else{
  	saved = 1;
      }
    }else if ((*it)->get_flag_saved()-(*it)->get_flag_saved_1()  ==2){
      if ( sqrt(pow(n_timeslices/8.,2) + pow(min_charge/n_mcells/8000.,2))<1 ||  min_charge/n_mcells/4000.<1.){
  	saved = 0;
      }else{
  	saved = 1;
      }
    }else if ((*it)->get_flag_saved()-(*it)->get_flag_saved_1()  ==1){
      if (sqrt(pow(n_timeslices/8.,2) + pow(min_charge/n_mcells/8000.,2))<1 || min_charge/n_mcells/6000.<1.){
  	saved = 0;
      }else{
  	saved = 1;
      }
    }else{
      saved = 0;
    }

    // test
    // saved = 1;
    
    // if (min_charge/n_mcells < 5000) saved = 0;
    if (saved==1){
      ncluster_saved ++;
      nmcell_saved += mcells.size();
    }else{
      ncluster_deleted ++;
      nmcell_deleted += mcells.size();
      // remove them ...
      for (auto it1 = mcells.begin(); it1!=mcells.end(); it1++){
      	SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)(*it1);
      	lowmemtiling[mcell->GetTimeSlice()]->Erase_Cell(mcell);
      }
    }
    
    //    T_3Dcluster->Fill();
  }

  for (int i=start_num;i!=end_num+1;i++){
    nmcell_after += lowmemtiling[i]->get_cell_wires_map().size();
  }
  
  std::cout << ncluster_saved << " " << nmcell_saved << " "
  	    << ncluster_deleted << " " << nmcell_deleted << " "
  	    << nmcell_before << " " << nmcell_after << " "
  	    << std::endl;
  
  cout << em("finish 2nd round of clustering and deghosting") << std::endl;


   for (int i=start_num;i!=end_num+1;i++){
    if (i%400==0)
      std::cout << "2nd Solving: " << i << std::endl;
    // tiling after the firs round of deghosting ... 
    lowmemtiling[i]->MergeWires();
    // create individual cells ...
    //GeomCellSelection single_cells = lowmemtiling[i]->create_single_cells();
    time_slice = i;
    n_cells = lowmemtiling[i]->get_cell_wires_map().size();//get_all_cell_centers().size();
    n_good_wires = lowmemtiling[i]->get_all_good_wires().size();
    n_bad_wires = lowmemtiling[i]->get_all_bad_wires().size();
    //n_single_cells = single_cells.size();
    // L1 solving
    chargesolver[i]->clear_connectivity();
    chargesolver[i]->L1_resolve(9,1);
    
    ndirect_solved = chargesolver[i]->get_ndirect_solved();
    nL1_solved = chargesolver[i]->get_nL1_solved();
    chargesolver[i]->Update_ndf_chi2();
    total_chi2 = chargesolver[i]->get_ndf();
    total_ndf = chargesolver[i]->get_chi2();
    for (Int_t k=0;k!=nL1_solved;k++){
      L1_ndf[k] = chargesolver[i]->get_L1_ndf(k);
      L1_chi2_base[k] = chargesolver[i]->get_L1_chi2_base(k);
      L1_chi2_penalty[k] = chargesolver[i]->get_L1_chi2_penalty(k);
    }
    for (Int_t k=0;k!=ndirect_solved;k++){
      direct_ndf[k] = chargesolver[i]->get_direct_ndf(k);
      direct_chi2[k] = chargesolver[i]->get_direct_chi2(k);
    }
    //Twc->Fill();
  }

   
   // label merge cells according to connectivities, going through clusters ..
   front_cell_map.clear();
   back_cell_map.clear();
   for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
     (*it)->Form_maps(2, front_cell_map,back_cell_map);
   }
   // std::cout << front_cell_map.size() << " " << back_cell_map.size() << std::endl;
   
   // repeat solving by changing the weight // getting weight ...
   
   for (auto it = front_cell_map.begin(); it!=front_cell_map.end();it++){
     SlimMergeGeomCell *mcell1 = (SlimMergeGeomCell*) it->first;
     int time_slice1 = mcell1->GetTimeSlice();
     bool flag1 = chargesolver[time_slice1]->get_mcell_charge(mcell1)>300;
     //std::cout << flag1 << std::endl;
     for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
       SlimMergeGeomCell *mcell2 = (SlimMergeGeomCell*)(*it1);
       int time_slice2 = mcell2->GetTimeSlice();
       bool flag2 = chargesolver[time_slice2]->get_mcell_charge(mcell2)>300;
       if (flag2)
	 chargesolver[time_slice1]->add_front_connectivity(mcell1);
       if(flag1)
	 chargesolver[time_slice2]->add_back_connectivity(mcell2);
     }
   }
  
   for (int i=start_num; i!=end_num+1;i++){
     if (i%400==0)
       std::cout << "2nd Solving with connectivity: " << i << std::endl;
     chargesolver[i]->L1_resolve(9,3);
   }
   
   potential_good_mcells.clear();
   good_mcells.clear();
   for (auto it = front_cell_map.begin(); it!=front_cell_map.end();it++){
     SlimMergeGeomCell *mcell1 = (SlimMergeGeomCell*) it->first;
     int time_slice1 = mcell1->GetTimeSlice();
     bool flag1 = chargesolver[time_slice1]->get_mcell_charge(mcell1)>300;
     if (flag1)
       good_mcells.insert(mcell1);
     for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
       SlimMergeGeomCell *mcell2 = (SlimMergeGeomCell*)(*it1);
       int time_slice2 = mcell2->GetTimeSlice();
       bool flag2 = chargesolver[time_slice2]->get_mcell_charge(mcell2)>300;
       if (flag1)
	 potential_good_mcells.insert(mcell2);
       if (flag2){
	 potential_good_mcells.insert(mcell1);
	 good_mcells.insert(mcell2);
       }
     }
   }
   for (int i=start_num;i!=end_num+1;i++){
     WCP::GeomCellMap& cell_map = lowmemtiling[i]->get_cell_wires_map();
     for (auto it = cell_map.begin(); it!=cell_map.end(); it++){
       SlimMergeGeomCell *mcell = (SlimMergeGeomCell*) it->first;
       bool flag1 = chargesolver[i]->get_mcell_charge(mcell)>300;
       if (flag1){
	 good_mcells.insert(mcell);
	 potential_good_mcells.insert(mcell);
       }
     }
   }
   
   cout << em("finish the 2nd round of solving") << std::endl;
   
   
   nmcell_before = 0;
   for (int i=start_num;i!=end_num+1;i++){
     nmcell_before += lowmemtiling[i]->get_cell_wires_map().size();
   }

   std::map<SlimMergeGeomCell*, double> map_mcell_charge;
   for (auto it = potential_good_mcells.begin(); it!=potential_good_mcells.end();it++){
     map_mcell_charge[*it] = chargesolver[(*it)->GetTimeSlice()]->get_mcell_charge(*it);
   }
   
   for (int i=start_num; i!=end_num+1;i++){
     lowmemtiling[i]->local_deghosting1(potential_good_mcells, map_mcell_charge);//(potential_good_mcells,false);
   }
   
   nmcell_after = 0;
   for (int i=start_num;i!=end_num+1;i++){
     lowmemtiling[i]->re_establish_maps();
     nmcell_after += lowmemtiling[i]->get_cell_wires_map().size();
   }
   std::cout << nmcell_before << " " << nmcell_after << std::endl;

   //   cout << em("finish the local deg ... ") << std::endl;
   
   
   //solve again ... 
   for (int i=start_num;i!=end_num+1;i++){
     if (i%400==0)
       std::cout << "3rd Solving: " << i << std::endl;
     // tiling after the firs round of deghosting ... 
     lowmemtiling[i]->MergeWires();
     // create individual cells ...
     //GeomCellSelection single_cells = lowmemtiling[i]->create_single_cells();
     time_slice = i;
     n_cells = lowmemtiling[i]->get_cell_wires_map().size();//get_all_cell_centers().size();
     n_good_wires = lowmemtiling[i]->get_all_good_wires().size();
     n_bad_wires = lowmemtiling[i]->get_all_bad_wires().size();
     //n_single_cells = single_cells.size();
     // L1 solving
     chargesolver[i]->clear_connectivity();
     chargesolver[i]->L1_resolve(9,1);
     
     ndirect_solved = chargesolver[i]->get_ndirect_solved();
     nL1_solved = chargesolver[i]->get_nL1_solved();
     chargesolver[i]->Update_ndf_chi2();
     total_chi2 = chargesolver[i]->get_ndf();
     total_ndf = chargesolver[i]->get_chi2();
     for (Int_t k=0;k!=nL1_solved;k++){
       L1_ndf[k] = chargesolver[i]->get_L1_ndf(k);
       L1_chi2_base[k] = chargesolver[i]->get_L1_chi2_base(k);
       L1_chi2_penalty[k] = chargesolver[i]->get_L1_chi2_penalty(k);
     }
     for (Int_t k=0;k!=ndirect_solved;k++){
       direct_ndf[k] = chargesolver[i]->get_direct_ndf(k);
       direct_chi2[k] = chargesolver[i]->get_direct_chi2(k);
     }
     //Twc->Fill();
   }

   
   
   
   // label merge cells according to connectivities, going through clusters ..
   front_cell_map.clear();
   back_cell_map.clear();
   for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
     (*it)->Form_maps(2, front_cell_map,back_cell_map);
   }
   // std::cout << front_cell_map.size() << " " << back_cell_map.size() << std::endl;
   
   // repeat solving by changing the weight // getting weight ...
   
   for (auto it = front_cell_map.begin(); it!=front_cell_map.end();it++){
     SlimMergeGeomCell *mcell1 = (SlimMergeGeomCell*) it->first;
     int time_slice1 = mcell1->GetTimeSlice();
     bool flag1 = chargesolver[time_slice1]->get_mcell_charge(mcell1)>300;
     //std::cout << flag1 << std::endl;
     for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
       SlimMergeGeomCell *mcell2 = (SlimMergeGeomCell*)(*it1);
       int time_slice2 = mcell2->GetTimeSlice();
       bool flag2 = chargesolver[time_slice2]->get_mcell_charge(mcell2)>300;
       if (flag2)
	 chargesolver[time_slice1]->add_front_connectivity(mcell1);
       if(flag1)
	 chargesolver[time_slice2]->add_back_connectivity(mcell2);
     }
   }
   
   for (int i=start_num; i!=end_num+1;i++){
     if (i%400==0)
       std::cout << "3rd Solving with connectivity: " << i << std::endl;
     chargesolver[i]->L1_resolve(9,3);
   }

   //   cout << em("finish 3rd solving ... ") << std::endl;
   
   potential_good_mcells.clear();
   good_mcells.clear();
   for (auto it = front_cell_map.begin(); it!=front_cell_map.end();it++){
     SlimMergeGeomCell *mcell1 = (SlimMergeGeomCell*) it->first;
     int time_slice1 = mcell1->GetTimeSlice();
     bool flag1 = chargesolver[time_slice1]->get_mcell_charge(mcell1)>300;
     if (flag1)
       good_mcells.insert(mcell1);
     for (auto it1 = it->second.begin(); it1!=it->second.end(); it1++){
       SlimMergeGeomCell *mcell2 = (SlimMergeGeomCell*)(*it1);
       int time_slice2 = mcell2->GetTimeSlice();
       bool flag2 = chargesolver[time_slice2]->get_mcell_charge(mcell2)>300;
       if (flag1)
	 potential_good_mcells.insert(mcell2);
       if (flag2){
	 potential_good_mcells.insert(mcell1);
	 good_mcells.insert(mcell2);
       }
     }
   }
   for (int i=start_num;i!=end_num+1;i++){
     WCP::GeomCellMap& cell_map = lowmemtiling[i]->get_cell_wires_map();
     for (auto it = cell_map.begin(); it!=cell_map.end(); it++){
       SlimMergeGeomCell *mcell = (SlimMergeGeomCell*) it->first;
       bool flag1 = chargesolver[i]->get_mcell_charge(mcell)>300;
       if (flag1){
	 good_mcells.insert(mcell);
	 potential_good_mcells.insert(mcell);
       }
     }
   }
   std::cout << good_mcells.size() << std::endl;
   
   
   
   // cluster again
   // delete clusters here ... 
   for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
     delete *it;
   }
   cluster_set.clear();
   cluster_delset.clear();
   
   
   {
     Slim3DClusterSet temp_cluster_set, temp_cluster_delset;
     // 2nd round of clustering
     for (int i=start_num; i!=end_num+1;i++){
       WCP::GeomCellMap& cell_wires_map = lowmemtiling[i]->get_cell_wires_map();
       GeomCellSelection allmcell;
       for (auto it=cell_wires_map.begin(); it!= cell_wires_map.end(); it++){
	 SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)it->first;
	 if (potential_good_mcells.find(mcell)!=potential_good_mcells.end())
	   allmcell.push_back(mcell);
       }
       if (temp_cluster_set.empty()){
	 // if cluster is empty, just insert all the mcell, each as a cluster
	 for (int j=0;j!=allmcell.size();j++){
	   Slim3DCluster *cluster = new Slim3DCluster(*((SlimMergeGeomCell*)allmcell[j]));
	   temp_cluster_set.insert(cluster);
	 }
       }else{
	 for (int j=0;j!=allmcell.size();j++){
	   int flag = 0;
	   int flag_save = 0;
	   Slim3DCluster *cluster_save = 0;
	   temp_cluster_delset.clear();
	   
	   // loop through merged cell
	   for (auto it = temp_cluster_set.begin();it!=temp_cluster_set.end();it++){
	     //loop through clusters
	     
	     flag += (*it)->AddCell(*((SlimMergeGeomCell*)allmcell[j]),2);
	     if (flag==1 && flag != flag_save){
	       cluster_save = *it;
	     }else if (flag>1 && flag != flag_save){
	       cluster_save->MergeCluster(*(*it));
	       temp_cluster_delset.insert(*it);
	     }
	     flag_save = flag;
	   }
	   
	   for (auto it = temp_cluster_delset.begin();it!=temp_cluster_delset.end();it++){
	     temp_cluster_set.erase(*it);
	     delete (*it);
	   }
	   
	   if (flag==0){
	     Slim3DCluster *cluster = new Slim3DCluster(*((SlimMergeGeomCell*)allmcell[j]));
	     temp_cluster_set.insert(cluster);
	   }
	 }
	 
	 // int ncount_mcell_cluster = 0;
	 // for (auto it = temp_cluster_set.begin();it!=temp_cluster_set.end();it++){
	 // 	ncount_mcell_cluster += (*it)->get_allcell().size();
	// }
	// ncount_mcell += allmcell.size();
	// cout << i << " " << allmcell.size()  << " " << temp_cluster_set.size()  << endl;
      }
    }

     for (auto it = temp_cluster_set.begin(); it!= temp_cluster_set.end(); it++){
       cluster_set.insert(*it);
     }
   }

   //cout << em("finish local clustering ... ") << std::endl;
   
   
   // std::cout << cluster_set.size() << std::endl;
   // save cluster into the output file ... 
   Double_t x,y,z,q,nq;
   //save cluster
   Int_t ncluster = 0;
   TTree *T_cluster ;
   if (save_file!=0){
     T_cluster = new TTree("T_cluster","T_cluster");
     T_cluster->Branch("cluster_id",&ncluster,"cluster_id/I");
     T_cluster->Branch("x",&x,"x/D");
     T_cluster->Branch("y",&y,"y/D");
     T_cluster->Branch("z",&z,"z/D");
     T_cluster->Branch("q",&q,"q/D");
     T_cluster->Branch("nq",&nq,"nq/D");
     
     T_cluster->SetDirectory(file);
   }
   
   good_mcells.clear();
   for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
     (*it)->set_id(ncluster);
     //   TGraph2D *g1 = new TGraph2D();
     for (int i=0; i!=(*it)->get_allcell().size();i++){
       SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)((*it)->get_allcell().at(i));
       good_mcells.insert(mcell);
       if (save_file!=0){
	 int time_slice_save = mcell->GetTimeSlice();
	 GeomCellSelection temp_cells = lowmemtiling[time_slice_save]->create_single_cells((SlimMergeGeomCell*)mcell);

	 nq = temp_cells.size();
	 if (nq>0){
	   q = chargesolver[time_slice_save]->get_mcell_charge(mcell) / (temp_cells.size() * 1.0) ;
	 }else{
	   q = 0;
	 }
	 
	 
	 for (int j=0; j!=temp_cells.size();j++){
	   Point p = temp_cells.at(j)->center();
	   x = mcell->GetTimeSlice()*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
	   y = p.y/units::cm;
	   z = p.z/units::cm;
	   T_cluster->Fill();
	   // 	g1->SetPoint(ncount,x,y,z);
	   // 	ncount ++;
	 }
       }
     }
     // cout << ncount << endl;
     //   g1->Write(Form("cluster_%d",ncluster));
     ncluster ++;
   }
     
  cout << em("finish the local deghosting ... ") << std::endl;

  // try to group the dead cells ...
  for (int i=start_num;i!=end_num+1;i++){
    GeomCellSelection& allmcell = lowmemtiling[i]->get_two_bad_wire_cells();
    if (lowmemtiling[i]->get_regen_two_bad_wire_cells() || dead_cluster_set.empty()){
      if (dead_region && dead_cluster_set.empty()){
	Slim3DDeadCluster *cluster = new Slim3DDeadCluster(*((SlimMergeGeomCell*)allmcell[0]),i);
	dead_cluster_set.insert(cluster);
      }
      for (int j=0;j<allmcell.size();j++){
	int flag = 0;
	int flag_save = 0;
	Slim3DDeadCluster *cluster_save = 0;
	dead_cluster_delset.clear();

	// loop through merged cell
	for (auto it = dead_cluster_set.begin();it!=dead_cluster_set.end();it++){
	  //loop through clusters
	  flag += (*it)->AddCell(*((SlimMergeGeomCell*)allmcell[j]),i);
	  if (flag==1 && flag != flag_save){
	    cluster_save = *it;
	  }else if (flag>1 && flag != flag_save){
	    cluster_save->MergeCluster(*(*it));
	    dead_cluster_delset.insert(*it);
	  }
	  flag_save = flag;
	}
	
	//std::cout << i << " " << j << " " << flag << " " << flag_save << std::endl;
	for (auto it = dead_cluster_delset.begin();it!=dead_cluster_delset.end();it++){
	  dead_cluster_set.erase(*it);
	  delete (*it);
	}
	
	if (flag==0){
	  Slim3DDeadCluster *cluster = new Slim3DDeadCluster(*((SlimMergeGeomCell*)allmcell[j]),i);
	  dead_cluster_set.insert(cluster);
	}
      }
    }else{
      for (auto it = dead_cluster_set.begin();it!=dead_cluster_set.end();it++){
	(*it)->Extend(i);
      }
    }
    //std::cout << i << " " << dead_cluster_set.size() << " " << allmcell.size() << std::endl;

    for (auto it = dead_cluster_set.begin();it!=dead_cluster_set.end();it++){
      if (i - ((*it)->get_cluster()).rbegin()->first >0){
	dead_cluster_set_save.insert(*it);
	dead_cluster_set.erase(*it);
      }
    }
  }

  for (auto it = dead_cluster_set.begin();it!=dead_cluster_set.end();it++){
    dead_cluster_set_save.insert(*it);
  }
  dead_cluster_set = dead_cluster_set_save;
  dead_cluster_set_save.clear();

  int cluster_id1 = 0;
  for (auto it = dead_cluster_set.begin();it!=dead_cluster_set.end();it++){
    (*it)->set_id(cluster_id1);
    cluster_id1++;
  }
    
  // for (int i=start_num;i!=end_num+1;i++){
  //   GeomCellSelection allmcell = lowmemtiling[i]->get_two_bad_wire_cells();
  //   for (int j=0;j<allmcell.size();j++){
  //     bool is_contain = false;
  //     for (auto it = dead_cluster_set.begin();it!=dead_cluster_set.end();it++){
  // 	is_contain = (*it)->IsContain(*((SlimMergeGeomCell*)allmcell[j]),i);
  // 	if (is_contain) break;
  //     }
  //     if(!is_contain) std::cout << i << " " << j << " " << "Wrong!" << std::endl;
  //   }
  // }
  // std::cout << "Done ... " << std::endl;
  // end to group the dead cells ... 
  cout << em("finish cluster dead region ... ") << std::endl;
  // for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
  //   for (auto it1 = dead_cluster_set.begin();it1!=dead_cluster_set.end();it1++){
  //     GeomCellSelection mcells = (*it)->Is_Connected(*it1,2);
  //     if (mcells.size()>0)
  // 	std::cout << (*it)->get_id() << " " << (*it)->get_allcell().size() << " " << (*it1)->get_id() << " " << (*it1)->get_mcells().size() << " " << mcells.size() << std::endl;
  //   }
  // }
  // cout << em("form map between dead and live clusters ... ") << std::endl;
  
  if (save_file==1){
    Double_t x_save, y_save, z_save;
    Double_t charge_save;
    Double_t ncharge_save;
    Double_t chi2_save;
    Double_t ndf_save;
    
    Double_t uq, vq, wq, udq, vdq, wdq;
    Int_t time_slice_save;
    
    TGraph2D *g = new TGraph2D();
    TTree *t_rec_simple = new TTree("T_rec","T_rec");
    t_rec_simple->SetDirectory(file);
    t_rec_simple->Branch("x",&x_save,"x/D");
    t_rec_simple->Branch("y",&y_save,"y/D");
    t_rec_simple->Branch("z",&z_save,"z/D");
    
    TGraph2D *g_rec = new TGraph2D();
    TTree *t_rec_charge = new TTree("T_rec_charge","T_rec_charge");
    t_rec_charge->SetDirectory(file);
    t_rec_charge->Branch("x",&x_save,"x/D");
    t_rec_charge->Branch("y",&y_save,"y/D");
    t_rec_charge->Branch("z",&z_save,"z/D");
    t_rec_charge->Branch("q",&charge_save,"q/D");
    t_rec_charge->Branch("nq",&ncharge_save,"nq/D");
    t_rec_charge->Branch("chi2",&chi2_save,"chi2/D");
    t_rec_charge->Branch("ndf",&ndf_save,"ndf/D");

    
    TGraph2D *g_deblob = new TGraph2D();
    TTree *t_rec_deblob = new TTree("T_rec_charge_blob","T_rec_charge_blob");
    t_rec_deblob->SetDirectory(file);
    t_rec_deblob->Branch("x",&x_save,"x/D");
    t_rec_deblob->Branch("y",&y_save,"y/D");
    t_rec_deblob->Branch("z",&z_save,"z/D");
    t_rec_deblob->Branch("q",&charge_save,"q/D");
    t_rec_deblob->Branch("nq",&ncharge_save,"nq/D");
    t_rec_deblob->Branch("chi2",&chi2_save,"chi2/D");
    t_rec_deblob->Branch("ndf",&ndf_save,"ndf/D");

    TTree *t_mcell = new TTree("T_mcell","T_mcell");
    t_mcell->SetDirectory(file);
    t_mcell->Branch("t",&time_slice_save,"t/I");
    t_mcell->Branch("x",&x_save,"x/D");
    t_mcell->Branch("y",&y_save,"y/D");
    t_mcell->Branch("z",&z_save,"z/D");
    t_mcell->Branch("q",&charge_save,"q/D");
    t_mcell->Branch("uq",&uq,"uq/D");
    t_mcell->Branch("vq",&vq,"vq/D");
    t_mcell->Branch("wq",&wq,"wq/D");
    t_mcell->Branch("udq",&udq,"udq/D");
    t_mcell->Branch("vdq",&vdq,"vdq/D");
    t_mcell->Branch("wdq",&wdq,"wdq/D");
    
    for (int i=start_num; i!=end_num+1;i++){
      GeomCellMap cell_wires_map = lowmemtiling[i]->get_cell_wires_map();
      GeomWireMap wire_cells_map = lowmemtiling[i]->get_wire_cells_map();
      GeomCellSelection three_good_wire_cells = lowmemtiling[i]->get_three_good_wire_cells();
      WCP::WireChargeMap& wire_charge = lowmemtiling[i]->get_wire_charge_map();
      WCP::WireChargeMap& wire_charge_error = lowmemtiling[i]->get_wire_charge_error_map();
      
      chi2_save = chargesolver[i]->get_chi2();
      ndf_save = chargesolver[i]->get_ndf();
      
      for (auto it = cell_wires_map.begin(); it!= cell_wires_map.end(); it++){
	SlimMergeGeomCell *mcell = (SlimMergeGeomCell*) it->first;
	GeomCellSelection temp_cells = lowmemtiling[i]->create_single_cells((SlimMergeGeomCell*)it->first);
	for (auto it1 = temp_cells.begin(); it1!=temp_cells.end(); it1++){
	  Point p = (*it1)->center();
	  x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
	  y_save = p.y/units::cm;
	  z_save = p.z/units::cm;
	  g->SetPoint(ncount,x_save, y_save, z_save);
	  t_rec_simple->Fill();
	  ncount ++;
	}
	
	// need the cell to be very good,
	bool flag_save = true;
	// if (cell_wires_map[mcell].size()==3 && find(three_good_wire_cells.begin(), three_good_wire_cells.end(), mcell)!= three_good_wire_cells.end()){
	// 	for (auto it1 = cell_wires_map[mcell].begin(); it1!= cell_wires_map[mcell].end(); it1++){
	// 	  MergeGeomWire *mwire = (MergeGeomWire*)(*it1);
	// 	  if (wire_cells_map[mwire].size()!=1){
	// 	    flag_save = false;
	// 	    break;
	// 	  }
	// 	}
	// }else{
	// 	flag_save = false;
	// }
	// //  std::cout << flag_save << std::endl;
	
	// if (chargesolver[i]->get_mcell_charge(mcell)>300 && flag_save)
	{
	  charge_save = chargesolver[i]->get_mcell_charge(mcell);
	  time_slice_save = mcell->GetTimeSlice();
	  x_save = 0;
	  y_save = 0;
	  z_save = 0;
	  
	  for (auto it1 = temp_cells.begin(); it1!=temp_cells.end(); it1++){
	    Point p = (*it1)->center();
	    x_save += i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
	    y_save += p.y/units::cm;
	    z_save += p.z/units::cm;
	  }
	  x_save /= temp_cells.size();
	  y_save /= temp_cells.size();
	  z_save /= temp_cells.size();
	  
	  for (auto it1 = cell_wires_map[mcell].begin(); it1!= cell_wires_map[mcell].end(); it1++){
	    MergeGeomWire *mwire = (MergeGeomWire*)(*it1);
	    if (mwire->get_allwire().front()->iplane()==0){
	      uq = wire_charge[mwire];
	      udq = wire_charge_error[mwire];
	    }else if(mwire->get_allwire().front()->iplane()==1){
	      vq = wire_charge[mwire];
	      vdq = wire_charge_error[mwire];
	    }else if(mwire->get_allwire().front()->iplane()==2){
	      wq = wire_charge[mwire];
	      wdq = wire_charge_error[mwire];
	    }
	  }
	  t_mcell->Fill();
	}
	
	// fill the charge ... 
	if (chargesolver[i]->get_mcell_charge(mcell)>300){
	  charge_save = chargesolver[i]->get_mcell_charge(mcell) / (temp_cells.size() * 1.0) ;
	  ncharge_save = temp_cells.size();
	  
	  for (auto it1 = temp_cells.begin(); it1!=temp_cells.end(); it1++){
	    Point p = (*it1)->center();
	    x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
	    y_save = p.y/units::cm;
	    z_save = p.z/units::cm;
	    g_rec->SetPoint(ncount1,x_save, y_save, z_save);
	    t_rec_charge->Fill();
	    
	    ncount1 ++;
	  }
	}
	
	if (good_mcells.find(mcell)!=good_mcells.end()){
	  if (chargesolver[i]->get_mcell_charge(mcell)>300){
	    charge_save = chargesolver[i]->get_mcell_charge(mcell) / (temp_cells.size() * 1.0) ;
	    ncharge_save = temp_cells.size();
	  }else{
	    charge_save = 300 / (temp_cells.size() * 1.0) ;
	    ncharge_save = temp_cells.size();
	  }
	  
	  for (auto it1 = temp_cells.begin(); it1!=temp_cells.end(); it1++){
	    Point p = (*it1)->center();
	    x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
	    y_save = p.y/units::cm;
	    z_save = p.z/units::cm;
	    g_deblob->SetPoint(ncount1,x_save, y_save, z_save);
	    t_rec_deblob->Fill();
	    
	    ncount1 ++;
	  }
	}
	
      }
      // 
    }
    g->Write("g");
    g_rec->Write("g_rec");
    g_deblob->Write("g_blob");
    
    TTree *t_bad = new TTree("T_bad","T_bad");
    t_bad->SetDirectory(file);
    Int_t bad_npoints;
    Double_t bad_y[100],bad_z[100];
    t_bad->Branch("bad_npoints",&bad_npoints,"bad_npoints/I");
    t_bad->Branch("bad_y",bad_y,"bad_y[bad_npoints]/D");
    t_bad->Branch("bad_z",bad_z,"bad_z[bad_npoints]/D");

    for (int i=0; i!=lowmemtiling[start_num]->get_two_bad_wire_cells().size();i++){
      const SlimMergeGeomCell *cell = (SlimMergeGeomCell*)lowmemtiling[start_num]->get_two_bad_wire_cells().at(i);
      PointVector ps = cell->boundary();
      bad_npoints = ps.size();
      for (int j=0;j!=bad_npoints;j++){
	bad_y[j] = ps.at(j).y/units::cm;
	bad_z[j] = ps.at(j).z/units::cm;
      }
      t_bad->Fill();
    }
    
  }

  // save
  if (save_file==0){
    TTree *TC = new TTree("TC","TC");
    TC->SetDirectory(file);
    
    Int_t cluster_id;
    Int_t time_slice;
    Double_t x_save, y_save, z_save;
    Double_t q, uq, vq, wq, udq, vdq, wdq;

    TC->Branch("cluster_id",&cluster_id,"cluster_id/I");
    TC->Branch("time_slice",&time_slice,"time_slice/I");

    // TC->Branch("x",&x_save,"x/D");
    // TC->Branch("y",&y_save,"y/D");
    // TC->Branch("z",&z_save,"z/D");
    TC->Branch("q",&q,"q/D");
    TC->Branch("uq",&uq,"uq/D");
    TC->Branch("vq",&vq,"vq/D");
    TC->Branch("wq",&wq,"wq/D");
    TC->Branch("udq",&udq,"udq/D");
    TC->Branch("vdq",&vdq,"vdq/D");
    TC->Branch("wdq",&wdq,"wdq/D");
    
    Int_t nwire_u=0, flag_u; //number of wires, dead?
    Int_t nwire_v=0, flag_v;
    Int_t nwire_w=0, flag_w;
    Int_t wire_index_u[2400];
    Int_t wire_index_v[2400];
    Int_t wire_index_w[3256];
    Double_t wire_charge_u[2400];
    Double_t wire_charge_v[2400];
    Double_t wire_charge_w[2400];
    Double_t wire_charge_err_u[2400];
    Double_t wire_charge_err_v[2400];
    Double_t wire_charge_err_w[2400];
    TC->Branch("nwire_u",&nwire_u,"nwire_u/I");
    TC->Branch("nwire_v",&nwire_v,"nwire_v/I");
    TC->Branch("nwire_w",&nwire_w,"nwire_w/I");
    TC->Branch("flag_u",&flag_u,"flag_u/I");
    TC->Branch("flag_v",&flag_v,"flag_v/I");
    TC->Branch("flag_w",&flag_w,"flag_w/I");
    TC->Branch("wire_index_u",wire_index_u,"wire_index_u[nwire_u]/I");
    TC->Branch("wire_index_v",wire_index_v,"wire_index_v[nwire_v]/I");
    TC->Branch("wire_index_w",wire_index_w,"wire_index_w[nwire_w]/I");
    TC->Branch("wire_charge_u",wire_charge_u,"wire_charge_u[nwire_u]/D");
    TC->Branch("wire_charge_v",wire_charge_v,"wire_charge_v[nwire_v]/D");
    TC->Branch("wire_charge_w",wire_charge_w,"wire_charge_w[nwire_w]/D");
    TC->Branch("wire_charge_err_u",wire_charge_err_u,"wire_charge_err_u[nwire_u]/D");
    TC->Branch("wire_charge_err_v",wire_charge_err_v,"wire_charge_err_v[nwire_v]/D");
    TC->Branch("wire_charge_err_w",wire_charge_err_w,"wire_charge_err_w[nwire_w]/D");
    
    for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
      cluster_id = (*it)->get_id();
    
      for (int i=0; i!=(*it)->get_allcell().size();i++){
	SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)((*it)->get_allcell().at(i));
	time_slice = mcell->GetTimeSlice();
	q = chargesolver[time_slice]->get_mcell_charge(mcell);
	
	GeomCellMap& cell_wires_map = lowmemtiling[time_slice]->get_cell_wires_map();
	WCP::WireChargeMap& wire_charge = lowmemtiling[time_slice]->get_wire_charge_map();
	WCP::WireChargeMap& wire_charge_error = lowmemtiling[time_slice]->get_wire_charge_error_map();
	for (auto it1 = cell_wires_map[mcell].begin(); it1!= cell_wires_map[mcell].end(); it1++){
	  MergeGeomWire *mwire = (MergeGeomWire*)(*it1);
	  if (mwire->get_allwire().front()->iplane()==0){
	    uq = wire_charge[mwire];
	    udq = wire_charge_error[mwire];
	  }else if(mwire->get_allwire().front()->iplane()==1){
	    vq = wire_charge[mwire];
	    vdq = wire_charge_error[mwire];
	  }else if(mwire->get_allwire().front()->iplane()==2){
	    wq = wire_charge[mwire];
	    wdq = wire_charge_error[mwire];
	  }
	}
	
      
	GeomWireSelection& uwires = mcell->get_uwires();
	GeomWireSelection& vwires = mcell->get_vwires();
	GeomWireSelection& wwires = mcell->get_wwires();
	nwire_u = uwires.size();
	nwire_v = vwires.size();
	nwire_w = wwires.size();

	std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();
	flag_u = 1;
	flag_v = 1;
	flag_w = 1;
	for (size_t j= 0 ; j!=bad_planes.size(); j++){
	  if (bad_planes.at(j)==WirePlaneType_t(0)){
	    flag_u = 0;
	  }else if (bad_planes.at(j)==WirePlaneType_t(1)){
	    flag_v = 0;
	  }else if (bad_planes.at(j)==WirePlaneType_t(2)){
	    flag_w = 0;
	  }
	} 
	for (int j=0;j!=nwire_u;j++){
	  const GeomWire *wire = uwires.at(j);
	  wire_index_u[j] = wire->index();
	  wire_charge_u[j] = wire_charge[wire];
	  wire_charge_err_u[j] = wire_charge_error[wire];
	}
	for (int j=0;j!=nwire_v;j++){
	  const GeomWire *wire = vwires.at(j);
	  wire_index_v[j] = wire->index();
	  wire_charge_v[j] = wire_charge[wire];
	  wire_charge_err_v[j] = wire_charge_error[wire];
	}
	for (int j=0;j!=nwire_w;j++){
	  const GeomWire *wire = wwires.at(j);
	  wire_index_w[j] = wire->index();
	  wire_charge_w[j] = wire_charge[wire];
	  wire_charge_err_w[j] = wire_charge_error[wire];
	}
	
	
	
	TC->Fill();
      }
    } // end saving live stuff ...

    TTree *TDC = new TTree("TDC","TDC");
    TDC->SetDirectory(file);
    TDC->Branch("cluster_id",&cluster_id,"cluster_id/I");
    int ntime_slice = 0,time_slices[2400];
    TDC->Branch("ntime_slice",&ntime_slice,"ntime_slice/I");
    TDC->Branch("time_slice",time_slices,"time_slice[ntime_slice]/I");

    TDC->Branch("nwire_u",&nwire_u,"nwire_u/I");
    TDC->Branch("nwire_v",&nwire_v,"nwire_v/I");
    TDC->Branch("nwire_w",&nwire_w,"nwire_w/I");
    TDC->Branch("flag_u",&flag_u,"flag_u/I");
    TDC->Branch("flag_v",&flag_v,"flag_v/I");
    TDC->Branch("flag_w",&flag_w,"flag_w/I");
    TDC->Branch("wire_index_u",wire_index_u,"wire_index_u[nwire_u]/I");
    TDC->Branch("wire_index_v",wire_index_v,"wire_index_v[nwire_v]/I");
    TDC->Branch("wire_index_w",wire_index_w,"wire_index_w[nwire_w]/I");

    
    for (auto it = dead_cluster_set.begin();it!=dead_cluster_set.end();it++){
      cluster_id = (*it)->get_id();
      std::map<SlimMergeGeomCell*,std::set<int>>& results = (*it)->get_mcell_time_map();
      for (auto it1 = results.begin(); it1!=results.end(); it1++){
	SlimMergeGeomCell* mcell = it1->first;
	std::set<int>& times = it1->second;
	ntime_slice = times.size();
	int temp_num = 0;
	for (auto it2 = times.begin(); it2!=times.end();it2++){
	  time_slices[temp_num] = *it2;
	  temp_num ++;
	}

	GeomWireSelection& uwires = mcell->get_uwires();
    	GeomWireSelection& vwires = mcell->get_vwires();
    	GeomWireSelection& wwires = mcell->get_wwires();
    	nwire_u = uwires.size();
    	nwire_v = vwires.size();
    	nwire_w = wwires.size();

    	std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();
    	flag_u = 1;
    	flag_v = 1;
    	flag_w = 1;
    	for (size_t j= 0 ; j!=bad_planes.size(); j++){
    	  if (bad_planes.at(j)==WirePlaneType_t(0)){
    	    flag_u = 0;
    	  }else if (bad_planes.at(j)==WirePlaneType_t(1)){
    	    flag_v = 0;
    	  }else if (bad_planes.at(j)==WirePlaneType_t(2)){
    	    flag_w = 0;
    	  }
    	}
	
    	for (int j=0;j!=nwire_u;j++){
    	  const GeomWire *wire = uwires.at(j);
    	  wire_index_u[j] = wire->index();
    	}
    	for (int j=0;j!=nwire_v;j++){
    	  const GeomWire *wire = vwires.at(j);
    	  wire_index_v[j] = wire->index();
    	}
    	for (int j=0;j!=nwire_w;j++){
    	  const GeomWire *wire = wwires.at(j);
    	  wire_index_w[j] = wire->index();
    	}
	TDC->Fill();
	
      }
    }
    
    
  }

  // save one event TC, TDC into one entry
  if (save_file==2){
    TTree *TC = new TTree("TC","TC");
    TC->SetDirectory(file);
    
    std::vector<int> *cluster_id = new std::vector<int>;
    std::vector<int> *time_slice = new std::vector<int>;
    //std::vector<double> *x_save = new std::vector<double>;    
    //std::vector<double> *y_save = new std::vector<double>;
    //std::vector<double> *z_save = new std::vector<double>;
    std::vector<double> *q = new std::vector<double>;    
    std::vector<double> *uq = new std::vector<double>;
    std::vector<double> *vq = new std::vector<double>;
    std::vector<double> *wq = new std::vector<double>;    
    std::vector<double> *udq = new std::vector<double>;
    std::vector<double> *vdq = new std::vector<double>;   
    std::vector<double> *wdq = new std::vector<double>;    

    TC->Branch("cluster_id",&cluster_id);

    TC->Branch("time_slice",&time_slice);

    // TC->Branch("x",&x_save,"x/D");
    // TC->Branch("y",&y_save,"y/D");
    // TC->Branch("z",&z_save,"z/D");
    TC->Branch("q",&q);
    TC->Branch("uq",&uq);
    TC->Branch("vq",&vq);
    TC->Branch("wq",&wq);
    TC->Branch("udq",&udq);
    TC->Branch("vdq",&vdq);
    TC->Branch("wdq",&wdq);
    
    std::vector<int> *nwire_u = new  std::vector<int>;
    std::vector<int> *nwire_v = new  std::vector<int>;
    std::vector<int> *nwire_w = new  std::vector<int>;
    std::vector<int> *flag_u = new  std::vector<int>;
    std::vector<int> *flag_v = new  std::vector<int>;
    std::vector<int> *flag_w = new  std::vector<int>;
    int flag_u_i;
    int flag_v_i;
    int flag_w_i;

    std::vector<std::vector<int>> *wire_index_u = new std::vector<std::vector<int>>;
    std::vector<std::vector<int>> *wire_index_v = new std::vector<std::vector<int>>;
    std::vector<std::vector<int>> *wire_index_w = new std::vector<std::vector<int>>;
    std::vector<std::vector<double>> *wire_charge_u = new std::vector<std::vector<double>>;
    std::vector<std::vector<double>> *wire_charge_v = new std::vector<std::vector<double>>;
    std::vector<std::vector<double>> *wire_charge_w = new std::vector<std::vector<double>>;
    std::vector<std::vector<double>> *wire_charge_err_u = new std::vector<std::vector<double>>;
    std::vector<std::vector<double>> *wire_charge_err_v = new std::vector<std::vector<double>>;
    std::vector<std::vector<double>> *wire_charge_err_w = new std::vector<std::vector<double>>;

    std::vector<int> wire_index_ui;
    std::vector<int> wire_index_vi;
    std::vector<int> wire_index_wi;
    std::vector<double> wire_charge_ui;
    std::vector<double> wire_charge_vi;
    std::vector<double> wire_charge_wi;
    std::vector<double> wire_charge_err_ui;
    std::vector<double> wire_charge_err_vi;
    std::vector<double> wire_charge_err_wi;

    TC->Branch("nwire_u",&nwire_u);
    TC->Branch("nwire_v",&nwire_v);
    TC->Branch("nwire_w",&nwire_w);
    TC->Branch("flag_u",&flag_u);
    TC->Branch("flag_v",&flag_v);
    TC->Branch("flag_w",&flag_w);
    TC->Branch("wire_index_u",&wire_index_u);
    TC->Branch("wire_index_v",&wire_index_v);
    TC->Branch("wire_index_w",&wire_index_w);
    TC->Branch("wire_charge_u",&wire_charge_u);
    TC->Branch("wire_charge_v",&wire_charge_v);
    TC->Branch("wire_charge_w",&wire_charge_w);
    TC->Branch("wire_charge_err_u",&wire_charge_err_u);
    TC->Branch("wire_charge_err_v",&wire_charge_err_v);
    TC->Branch("wire_charge_err_w",&wire_charge_err_w);
    
    for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
    //  cluster_id->push_back((*it)->get_id());
    
      for (int i=0; i!=(*it)->get_allcell().size();i++){
    cluster_id->push_back((*it)->get_id());
	SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)((*it)->get_allcell().at(i));
	int time_slice_temp = mcell->GetTimeSlice();
    time_slice->push_back(time_slice_temp);
	q->push_back(chargesolver[time_slice_temp]->get_mcell_charge(mcell));
	
	GeomCellMap cell_wires_map = lowmemtiling[time_slice_temp]->get_cell_wires_map();
	WCP::WireChargeMap& wire_charge = lowmemtiling[time_slice_temp]->get_wire_charge_map();
	WCP::WireChargeMap& wire_charge_error = lowmemtiling[time_slice_temp]->get_wire_charge_error_map();
	for (auto it1 = cell_wires_map[mcell].begin(); it1!= cell_wires_map[mcell].end(); it1++){
	  MergeGeomWire *mwire = (MergeGeomWire*)(*it1);
	  if (mwire->get_allwire().front()->iplane()==0){
	    uq->push_back(wire_charge[mwire]);
	    udq->push_back(wire_charge_error[mwire]);
	  }else if(mwire->get_allwire().front()->iplane()==1){
	    vq->push_back(wire_charge[mwire]);
	    vdq->push_back(wire_charge_error[mwire]);
	  }else if(mwire->get_allwire().front()->iplane()==2){
	    wq->push_back(wire_charge[mwire]);
	    wdq->push_back(wire_charge_error[mwire]);
	  }
	}
	
      
	GeomWireSelection& uwires = mcell->get_uwires();
	GeomWireSelection& vwires = mcell->get_vwires();
	GeomWireSelection& wwires = mcell->get_wwires();
	nwire_u->push_back(uwires.size());
	nwire_v->push_back(vwires.size());
	nwire_w->push_back(wwires.size());

	std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();
	flag_u_i = 1;
      flag_v_i = 1;
      flag_w_i = 1;
	for (size_t j= 0 ; j!=bad_planes.size(); j++){
	  if (bad_planes.at(j)==WirePlaneType_t(0)){
	    flag_u_i = 0;
	  }else if (bad_planes.at(j)==WirePlaneType_t(1)){
	    flag_v_i = 0;
	  }else if (bad_planes.at(j)==WirePlaneType_t(2)){
	    flag_w_i = 0;
	  }
	} 
      flag_u->push_back(flag_u_i);
 	flag_v->push_back(flag_v_i);
 	flag_w->push_back(flag_w_i);

	for (int j=0;j!=uwires.size();j++){
	  const GeomWire *wire = uwires.at(j);
	  wire_index_ui.push_back(wire->index());
	  wire_charge_ui.push_back(wire_charge[wire]);
	  wire_charge_err_ui.push_back(wire_charge_error[wire]);
	}
	for (int j=0;j!=vwires.size();j++){
	  const GeomWire *wire = vwires.at(j);
	  wire_index_vi.push_back(wire->index());
	  wire_charge_vi.push_back(wire_charge[wire]);
	  wire_charge_err_vi.push_back(wire_charge_error[wire]);
	}
	for (int j=0;j!=wwires.size();j++){
	  const GeomWire *wire = wwires.at(j);
	  wire_index_wi.push_back(wire->index());
	  wire_charge_wi.push_back(wire_charge[wire]);
	  wire_charge_err_wi.push_back(wire_charge_error[wire]);
	}
		
    wire_index_u->push_back(wire_index_ui);	
    wire_index_v->push_back(wire_index_vi);	
    wire_index_w->push_back(wire_index_wi);
    wire_charge_u->push_back(wire_charge_ui);
    wire_charge_v->push_back(wire_charge_vi);
    wire_charge_w->push_back(wire_charge_wi);
    wire_charge_err_u->push_back(wire_charge_err_ui);
    wire_charge_err_v->push_back(wire_charge_err_vi);
    wire_charge_err_w->push_back(wire_charge_err_wi);
    wire_index_ui.clear();
    wire_index_vi.clear();
    wire_index_wi.clear();
    wire_charge_ui.clear();
    wire_charge_vi.clear();
    wire_charge_wi.clear();
    wire_charge_err_ui.clear();
    wire_charge_err_vi.clear();
    wire_charge_err_wi.clear();
      } // each cell
    } // end saving live stuff ...
    
    TC->Fill();
    cluster_id->clear();
    nwire_u->clear();
    nwire_v->clear();
    nwire_w->clear();
    flag_u->clear();
    flag_v->clear();
    flag_w->clear();
    wire_index_u->clear();
    wire_index_v->clear();
    wire_index_w->clear();

    TTree *TDC = new TTree("TDC","TDC");
    TDC->SetDirectory(file);

    std::vector<int> *ntime_slice = new std::vector<int>;
    std::vector<std::vector<int>> *time_slices = new std::vector<std::vector<int>>;
    std::vector<int> time_slices_i;

    TDC->Branch("cluster_id",&cluster_id);
    TDC->Branch("ntime_slice",&ntime_slice);
    TDC->Branch("time_slice",&time_slices);
    
    TDC->Branch("nwire_u",&nwire_u);
    TDC->Branch("nwire_v",&nwire_v);
    TDC->Branch("nwire_w",&nwire_w);
    TDC->Branch("flag_u",&flag_u);
    TDC->Branch("flag_v",&flag_v);
    TDC->Branch("flag_w",&flag_w);
    TDC->Branch("wire_index_u",&wire_index_u);
    TDC->Branch("wire_index_v",&wire_index_v);
    TDC->Branch("wire_index_w",&wire_index_w);

    
    for (auto it = dead_cluster_set.begin();it!=dead_cluster_set.end();it++){
      //cluster_id = (*it)->get_id();
      std::map<SlimMergeGeomCell*,std::set<int>>& results = (*it)->get_mcell_time_map();
      for (auto it1 = results.begin(); it1!=results.end(); it1++){
    cluster_id->push_back((*it)->get_id());
	SlimMergeGeomCell* mcell = it1->first;
	std::set<int>& times = it1->second;
	ntime_slice->push_back(times.size());
	int temp_num = 0;
	for (auto it2 = times.begin(); it2!=times.end();it2++){
	  time_slices_i.push_back(*it2);
	  temp_num ++;
	}
    time_slices->push_back(time_slices_i);
    time_slices_i.clear();

	GeomWireSelection& uwires = mcell->get_uwires();
    	GeomWireSelection& vwires = mcell->get_vwires();
    	GeomWireSelection& wwires = mcell->get_wwires();
    	nwire_u->push_back(uwires.size());
    	nwire_v->push_back(vwires.size());
    	nwire_w->push_back(wwires.size());

    	std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();

	flag_u_i = 1;
      flag_v_i = 1;
      flag_w_i = 1; 
    	for (size_t j= 0 ; j!=bad_planes.size(); j++){
    	  if (bad_planes.at(j)==WirePlaneType_t(0)){
    	    flag_u_i = 0;
    	  }else if (bad_planes.at(j)==WirePlaneType_t(1)){
    	    flag_v_i = 0;
    	  }else if (bad_planes.at(j)==WirePlaneType_t(2)){
    	    flag_w_i = 0;
    	  }
    	}
      flag_u->push_back(flag_u_i);
 	flag_v->push_back(flag_v_i);
 	flag_w->push_back(flag_w_i);
	
    	for (int j=0;j!=uwires.size();j++){
    	  const GeomWire *wire = uwires.at(j);
    	  wire_index_ui.push_back(wire->index());
    	}
    	for (int j=0;j!=vwires.size();j++){
    	  const GeomWire *wire = vwires.at(j);
    	  wire_index_vi.push_back(wire->index());
    	}
    	for (int j=0;j!=wwires.size();j++){
    	  const GeomWire *wire = wwires.at(j);
    	  wire_index_wi.push_back(wire->index());
    	}

        wire_index_u->push_back(wire_index_ui);
        wire_index_v->push_back(wire_index_vi);
        wire_index_w->push_back(wire_index_wi);
        wire_index_ui.clear();
        wire_index_vi.clear();
        wire_index_wi.clear();
      }
    }
    
  TDC->Fill();
  }

   
  file->Write();
  file->Close();
  
  return 0;
  
} // main()
