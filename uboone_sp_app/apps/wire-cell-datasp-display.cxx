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
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root -s[0,1,2]" << endl;
    return 1;
  }

  TH1::AddDirectory(kFALSE);
  
  int save_file = 0; //
  
  for(Int_t i = 1; i != argc; i++){
     switch(argv[i][1]){
     case 's':
       save_file = atoi(&argv[i][2]); 
       break;
     }
  }

  
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
  
  cout << em("load data") << endl;

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
  TH2F *hv_raw = (TH2F*)file1->Get("hv_raw");
  
  if (save_file==0){
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
    //
   
    // V wire noisy channels 10 vetoed ...
    for (int i=3684;i!=3699;i++){
      if (vplane_map.find(i-2400)==vplane_map.end()){
	vplane_map[i-2400] = std::make_pair(0,hv_raw->GetNbinsY()-1);
	std::cout << "V plane (noisy): " << i -2400 << " added to bad channel list" << std::endl;
      }else{
	vplane_map[i-2400] = std::make_pair(0,hv_raw->GetNbinsY()-1);
      }
      for (int j=0;j!=hv_decon->GetNbinsY();j++){
	hv_decon->SetBinContent(i+1-2400,j+1,0);
	hv_decon_g->SetBinContent(i+1-2400,j+1,0);
      }
      for (int j=0;j!=hv_raw->GetNbinsY();j++){
	hv_raw->SetBinContent(i+1-2400,j+1,0);
      }
    }
    // U wire plane bad channels
    for (int i=2160;i!=2176;i++){
      if (uplane_map.find(i)==uplane_map.end()){
	uplane_map[i] = std::make_pair(0,hv_raw->GetNbinsY()-1);
	std::cout << "U plane (bad): " << i << " added to bad channel list" << std::endl;
      }else{
	uplane_map[i] = std::make_pair(0,hv_raw->GetNbinsY()-1);
      }
      for (int j=0;j!=hu_decon->GetNbinsY();j++){
	hu_decon->SetBinContent(i+1,j+1,0);
	hu_decon_g->SetBinContent(i+1,j+1,0);
      }
    }
    
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
    
    
    WCP2dToy::Noisy_Event_ID(hu_decon, hv_decon, hw_decon, uplane_rms, vplane_rms, wplane_rms, uplane_map, vplane_map, wplane_map, hu_decon_g, hv_decon_g, hw_decon_g, nrebin, hv_raw, true);
    WCP2dToy::Organize_Dead_Channels(uplane_map, vplane_map, wplane_map, hv_raw->GetNbinsY()-1,nrebin);
  }
  
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


   TFile *file = new TFile(Form("nsp_2D_display_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");

   Trun->CloneTree()->Write();

   TTree *T_bad = new TTree("T_bad","T_bad");
   // Int_t chid, plane;
   // Int_t start_time,end_time;
   T_bad->Branch("chid",&chid,"chid/I");
   T_bad->Branch("plane",&plane,"plane/I");
   T_bad->Branch("start_time",&start_time,"start_time/I");
   T_bad->Branch("end_time",&end_time,"end_time/I");
   T_bad->SetDirectory(file);
   
   for (auto it = uplane_map.begin(); it!=uplane_map.end(); it++){
     chid = it->first;
     plane = 0;
     start_time = it->second.first;
     end_time = it->second.second;
     T_bad->Fill();
   }
   for (auto it = vplane_map.begin(); it!=vplane_map.end(); it++){
     chid = it->first + 2400;
     plane = 1;
     start_time = it->second.first;
     end_time = it->second.second;
     T_bad->Fill();
   }
   for (auto it = wplane_map.begin(); it!=wplane_map.end(); it++){
     chid = it->first + 4800;
     plane = 2;
     start_time = it->second.first;
     end_time = it->second.second;
     T_bad->Fill();
   }
   
   
   
   hu_threshold->SetDirectory(file);
   hv_threshold->SetDirectory(file);
   hw_threshold->SetDirectory(file);

   hu_decon->SetDirectory(file);
   hv_decon->SetDirectory(file);
   hw_decon->SetDirectory(file);

   hu_decon_g->SetDirectory(file);
   hv_decon_g->SetDirectory(file);
   hw_decon_g->SetDirectory(file);

   hv_raw->SetDirectory(file);
   
  file->Write();
  file->Close();
  
  return 0;
  
} 
