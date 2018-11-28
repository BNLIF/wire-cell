#include "WireCellSst/GeomDataSource.h"
#include "WireCellSst/DatauBooNEFrameDataSource.h"
#include "WireCellSst/ToyuBooNESliceDataSource.h"
#include "WireCellSst/uBooNESliceDataSource.h"

#include "WireCell2dToy/ToyEventDisplay.h"
#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/BadTiling.h"
#include "WireCell2dToy/LowmemTiling.h"
#include "WireCell2dToy/uBooNE_L1SP.h"
#include "WireCell2dToy/WireCellHolder.h"

#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/TruthToyTiling.h"
#include "WireCell2dToy/SimpleBlobToyTiling.h"

#include "WireCell2dToy/ChargeSolving.h"
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

#include "WireCellData/Slim3DCluster.h"
#include "WireCellData/Slim3DDeadCluster.h"
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
#include "WireCell2dToy/pd_Data_FDS.h"
#include "WireCell2dToy/uBooNE_Data_Error.h"
#include "WireCell2dToy/ToyDataQuality.h"

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
  if (argc < 3) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root -t[0,1] -s[0,1,2]" << endl;
    return 1;
  }

  TH1::AddDirectory(kFALSE);

  int two_plane = 0; //not used
  
  int nt_off1 = 0; // not used
  int nt_off2 = 0; // not used
  int solve_charge = 1; // not used

  int save_file = 0; //
  // 1 for debug mode for bee ...

  for(Int_t i = 1; i != argc; i++){
     switch(argv[i][1]){
     case 't':
       two_plane = atoi(&argv[i][2]); 
       break;
     case 's':
       save_file = atoi(&argv[i][2]); 
       break;
     case 'd':
       solve_charge = atoi(&argv[i][2]); 
       break;
     case 'a':
       nt_off1 = atoi(&argv[i][2]);
       break;
     case 'b':
       nt_off2 = atoi(&argv[i][2]);
       break;
     }
  }

  // if(save_file==1){
  //   std::cout << "Save file for bee. " << std::endl;
  // }else if (save_file==0){
  //   std::cout << "Save file for pattern recognition. " << std::endl;
  // }
  
  // if (two_plane)
  //   cout << "Enable Two Plane Reconstruction " << endl; 
  

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
      
  //Note: the above one is still on the high end, 
  float unit_dis = 1.101; // match 256 cm

  //**** time offset for 58kV ****// 
  float toffset_1=0.0; //(nt_off1 * 0.2 - 1.0 );  // time offset between u/v 
  float toffset_2=0.0; //(nt_off2 * 0.2 - 1.0); // time offset between u/w
  float toffset_3=0.0;
  
  int save_image_outline_flag = 0; // prescale flag 
  
  int total_time_bin = 9592;
  int recon_threshold = 2000;
  int frame_length = 3200;
  int max_events = 100;
  int eve_num;
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
  

  int time_offset = 4; // Now the time offset is taken care int he signal processing, so we just need the overall offset ... 
  

  const char* root_file = argv[2];  
  int run_no, subrun_no, event_no;
  

  // load Trun
  TFile *file1 = new TFile(root_file);
  TTree *Trun = (TTree*)file1->Get("Trun");
  
  Trun->SetBranchAddress("eventNo",&event_no);
  Trun->SetBranchAddress("runNo",&run_no);
  Trun->SetBranchAddress("subRunNo",&subrun_no);
  Trun->SetBranchAddress("unit_dis",&unit_dis);
  Trun->SetBranchAddress("toffset_uv",&toffset_1);
  Trun->SetBranchAddress("toffset_uw",&toffset_2);
  Trun->SetBranchAddress("toffset_u",&toffset_3);
  Trun->SetBranchAddress("total_time_bin",&total_time_bin);
  Trun->SetBranchAddress("recon_threshold",&recon_threshold);
  Trun->SetBranchAddress("frame_length",&frame_length);
  Trun->SetBranchAddress("max_events",&max_events);
  Trun->SetBranchAddress("eve_num",&eve_num);
  Trun->SetBranchAddress("nrebin",&nrebin);
  Trun->SetBranchAddress("threshold_u",&threshold_u);
  Trun->SetBranchAddress("threshold_v",&threshold_v);
  Trun->SetBranchAddress("threshold_w",&threshold_w);
  Trun->SetBranchAddress("threshold_ug",&threshold_ug);
  Trun->SetBranchAddress("threshold_vg",&threshold_vg);
  Trun->SetBranchAddress("threshold_wg",&threshold_wg);
  Trun->SetBranchAddress("time_offset",&time_offset);
  
  Trun->GetEntry(0);

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
  TTree *T_chirp = (TTree*)file1->Get("T_chirp");
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

   // add a special treatment here ...
  for (int i=296; i!=671;i++){
    if (uplane_map.find(i)==uplane_map.end()){
      if ((uplane_map.find(i-1)!=uplane_map.end() &&
	   uplane_map.find(i+2)!=uplane_map.end()) ||
	  (uplane_map.find(i-2)!=uplane_map.end() &&
	   uplane_map.find(i+1)!=uplane_map.end())
	  ){
	std::cout << "U plane (shorted): " << i << " added to bad chanel list" << std::endl;
	uplane_map[i] = uplane_map[i-1];
	for (int j=0;j!=hu_decon->GetNbinsY();j++){
	  hu_decon->SetBinContent(i+1,j+1,0);
	  hu_decon_g->SetBinContent(i+1,j+1,0);
	}
      }
    }
  }
  for (int i=2336;i!=2463;i++){
    if (wplane_map.find(i)==wplane_map.end()){
      if ((wplane_map.find(i-1)!=wplane_map.end() &&
	   wplane_map.find(i+2)!=wplane_map.end()) ||
	  (wplane_map.find(i-2)!=wplane_map.end() &&
	   wplane_map.find(i+1)!=wplane_map.end())
	  ){
	wplane_map[i] = wplane_map[i-1];
	std::cout << "W plane (shorted): " << i << " added to bad channel list" << std::endl;
	for (int j=0;j!=hw_decon->GetNbinsY();j++){
	  hw_decon->SetBinContent(i+1,j+1,0);
	  hw_decon_g->SetBinContent(i+1,j+1,0);
	}
      }
    }
  }
  // V wire noisy channels 10 vetoed ...
  for (int i=3686;i!=3697;i++){
    if (vplane_map.find(i-2400)==vplane_map.end()){
      vplane_map[i-2400] = std::make_pair(0,hv_decon->GetNbinsY()-1);
      std::cout << "V plane (noisy): " << i -2400 << " added to bad channel list" << std::endl;
    }else{
      vplane_map[i-2400] = std::make_pair(0,hv_decon->GetNbinsY()-1);
    }
    for (int j=0;j!=hv_decon->GetNbinsY();j++){
      hv_decon->SetBinContent(i+1-2400,j+1,0);
      hv_decon_g->SetBinContent(i+1-2400,j+1,0);
      
    }
  }
  //

  //  std::cout << uplane_map.size() << " " << vplane_map.size() << " " << wplane_map.size() << std::endl;

  TH2F *hv_raw = (TH2F*)file1->Get("hv_raw");
  
  WireCell2dToy::Noisy_Event_ID(hu_decon, hv_decon, hw_decon, uplane_rms, vplane_rms, wplane_rms, uplane_map, vplane_map, wplane_map, hu_decon_g, hv_decon_g, hw_decon_g, nrebin, hv_raw, true);

  
  if (save_file==1){
    TFile *file1 = new TFile("temp.root","RECREATE");
    hu_decon->SetDirectory(file1);
    hv_decon->SetDirectory(file1);
    hw_decon->SetDirectory(file1);

    hu_decon_g->SetDirectory(file1);
    hv_decon_g->SetDirectory(file1);
    hw_decon_g->SetDirectory(file1);

    hv_raw->SetDirectory(file1);

    hu_threshold->SetDirectory(file1);
    hv_threshold->SetDirectory(file1);
    hw_threshold->SetDirectory(file1);

    Trun->CloneTree()->Write();
    T_chirp->CloneTree()->Write();
    file1->Write();
    file1->Close();
  }
  
  //std::cout << uplane_map.size() << " " << vplane_map.size() << " " << wplane_map.size() << std::endl;
}
