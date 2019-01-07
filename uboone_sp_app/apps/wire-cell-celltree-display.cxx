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
#include "WireCell2dToy/ToyDataQuality.h"

#include "WireCell2dToy/uBooNE_Data_2D_Deconvolution.h"
#include "WireCell2dToy/uBooNE_Data_ROI.h"
#include "WireCell2dToy/uBooNE_Data_After_ROI.h"
#include "WireCell2dToy/uBooNE_Data_After_ROI_gaus.h"

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
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root eve_num " << endl;
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
  int eve_num = atoi(argv[3]);
  int nrebin = 4;



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
  TTree *T = (TTree*)file1->Get("/Event/Sim");
  T->SetBranchAddress("eventNo",&event_no);
  T->SetBranchAddress("runNo",&run_no);
  T->SetBranchAddress("subRunNo",&subrun_no);

  // std::vector<int> *badChannelList = new std::vector<int>;
  // T->SetBranchAddress("badChannelList",&badChannelList);

  std::vector<int> *badChannel = new std::vector<int>;
  std::vector<int> *badBegin = new std::vector<int>;
  std::vector<int> *badEnd = new std::vector<int>;
  T->SetBranchAddress("badChannel",&badChannel);
  T->SetBranchAddress("badBegin",&badBegin);
  T->SetBranchAddress("badEnd",&badEnd);

  Short_t nf_shift, nf_scale;
  T->SetBranchAddress("nf_shift",&nf_shift);
  T->SetBranchAddress("nf_scale",&nf_scale);
  
  std::vector<double> *channelThreshold = new std::vector<double>;
  T->SetBranchAddress("channelThreshold",&channelThreshold);
  std::vector<int> *calibGaussian_channelId = new std::vector<int>;
  std::vector<int> *calibWiener_channelId = new std::vector<int>;
  std::vector<int> *nf_channelId = new std::vector<int>;
  T->SetBranchAddress("calibGaussian_channelId",&calibGaussian_channelId);
  T->SetBranchAddress("calibWiener_channelId",&calibWiener_channelId);
  T->SetBranchAddress("nf_channelId",&nf_channelId);
  TClonesArray* calibWiener_wf = new TClonesArray;
  TClonesArray* calibGaussian_wf = new TClonesArray;
  TClonesArray* nf_wf = new TClonesArray;
  // TH1F  ... 
  T->SetBranchAddress("calibWiener_wf",&calibWiener_wf);
  T->SetBranchAddress("calibGaussian_wf",&calibGaussian_wf);
  T->SetBranchAddress("nf_wf",&nf_wf);
  
  T->GetEntry(eve_num);
  cout << "Run No: " << run_no << " " << subrun_no << " " << event_no << endl;


  TTree *T_op = 0;//(TTree*)file1->Get("T_op");
  
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  int nwire_u = wires_u.size();
  int nwire_v = wires_v.size();
  int nwire_w = wires_w.size();
  
  ChirpMap uplane_map;
  ChirpMap vplane_map;
  ChirpMap wplane_map;

  Int_t chid=0, plane=0;
  Int_t start_time=0, end_time=9594;

  

  
  for (size_t i=0;i!=badChannel->size();i++){
    std::pair<int,int> abc(badBegin->at(i),badEnd->at(i));
    chid = badChannel->at(i);
    if (chid < nwire_u){
      uplane_map[chid] = abc;
    }else if (chid < nwire_u + nwire_v){
      vplane_map[chid-nwire_u] = abc;
    }else if (chid < nwire_u + nwire_v + nwire_w){
      wplane_map[chid-nwire_u-nwire_v] = abc;
    }
  }

  // std::cout << uplane_map.size() << " " << vplane_map.size() << " " << wplane_map.size() << std::endl;
  
  std::vector<float> uplane_rms;
  std::vector<float> vplane_rms;
  std::vector<float> wplane_rms;
  // Note, there is a mismatch here
  // These RMS values are for single ticks
  // Later they are applied to the rebinned data
  // Probably OK as the TPC signal processing already took care the fake hits ...
  for (size_t i=0; i!= channelThreshold->size(); i++){
    if (i<nwire_u){
      uplane_rms.push_back(channelThreshold->at(i));
    }else if (i<nwire_u + nwire_v){
      vplane_rms.push_back(channelThreshold->at(i));
    }else if (i<nwire_u + nwire_v + nwire_w){
      wplane_rms.push_back(channelThreshold->at(i));
    }
  }

  //  std::cout <<uplane_rms.size() << " " << vplane_rms.size() << " " << wplane_rms.size() << std::endl;
  
  const int nwindow_size = ((TH1F*)calibWiener_wf->At(0))->GetNbinsX();

  
  
  TH2F *hu_decon = new TH2F("hu_decon","hu_decon",nwire_u,-0.5,nwire_u-0.5,nwindow_size,0,nwindow_size*nrebin);
  TH2F *hv_decon = new TH2F("hv_decon","hv_decon",nwire_v,-0.5+nwire_u,-0.5+nwire_u+nwire_v,nwindow_size,0,nwindow_size*nrebin);
  TH2F *hw_decon = new TH2F("hw_decon","hw_decon",nwire_w,-0.5+nwire_u+nwire_v,-0.5+nwire_u+nwire_v+nwire_w,nwindow_size,0,nwindow_size*nrebin);
  
  TH2F *hu_decon_g = new TH2F("hu_decon_g","hu_decon_g",nwire_u,-0.5,nwire_u-0.5,nwindow_size,0,nwindow_size*nrebin);
  TH2F *hv_decon_g = new TH2F("hv_decon_g","hv_decon_g",nwire_v,-0.5+nwire_u,-0.5+nwire_u+nwire_v,nwindow_size,0,nwindow_size*nrebin);
  TH2F *hw_decon_g = new TH2F("hw_decon_g","hw_decon_g",nwire_w,-0.5+nwire_u+nwire_v,-0.5+nwire_u+nwire_v+nwire_w,nwindow_size,0,nwindow_size*nrebin);
  
  for (size_t i=0;i!=calibWiener_channelId->size();i++){
    int chid = calibWiener_channelId->at(i);
    TH1F *htemp = (TH1F*)calibWiener_wf->At(i);
    if (chid < nwire_u){
      for (size_t j=0;j!=nwindow_size;j++){
	hu_decon->SetBinContent(chid+1,j+1,htemp->GetBinContent(j+1));
      }
    }else if (chid < nwire_v+nwire_u){
      for (size_t j=0;j!=nwindow_size;j++){
	hv_decon->SetBinContent(chid-nwire_u+1,j+1,htemp->GetBinContent(j+1));
      }
    }else if (chid < nwire_w+nwire_v+nwire_u){
      for (size_t j=0;j!=nwindow_size;j++){
	hw_decon->SetBinContent(chid-nwire_u-nwire_v+1,j+1,htemp->GetBinContent(j+1));
      }
    }

    chid = calibGaussian_channelId->at(i);
    htemp = (TH1F*)calibGaussian_wf->At(i);
    if (chid < nwire_u){
      for (size_t j=0;j!=nwindow_size;j++){
	hu_decon_g->SetBinContent(chid+1,j+1,htemp->GetBinContent(j+1));
      }
    }else if (chid < nwire_v+nwire_u){
      for (size_t j=0;j!=nwindow_size;j++){
	hv_decon_g->SetBinContent(chid-nwire_u+1,j+1,htemp->GetBinContent(j+1));
      }
    }else if (chid < nwire_w+nwire_v+nwire_u){
      for (size_t j=0;j!=nwindow_size;j++){
	hw_decon_g->SetBinContent(chid-nwire_u-nwire_v+1,j+1,htemp->GetBinContent(j+1));
      }
    }
    
  }

  

  TFile *file = new TFile(Form("celltree_2D_display_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");
  
  
  TTree *Trun = new TTree("Trun","Trun");
  Trun->SetDirectory(file);

  int detector = 0; // MicroBooNE
  Trun->Branch("detector",&detector,"detector/I");
  Trun->Branch("eventNo",&event_no,"eventNo/I");
  Trun->Branch("runNo",&run_no,"runNo/I");
  Trun->Branch("subRunNo",&subrun_no,"subRunNo/I");
  Trun->Branch("unit_dis",&unit_dis,"unit_dis/F");
  Trun->Branch("toffset_uv",&toffset_1,"toffset_uv/F");
  Trun->Branch("toffset_uw",&toffset_2,"toffset_uw/F");
  Trun->Branch("toffset_u",&toffset_3,"toffset_u/F");
  Trun->Branch("total_time_bin",&total_time_bin,"total_time_bin/I");
  //Trun->Branch("recon_threshold",&recon_threshold,"recon_threshold/I");
  Trun->Branch("frame_length",&frame_length,"frame_length/I");
  //Trun->Branch("max_events",&max_events,"max_events/I");
  Trun->Branch("eve_num",&eve_num,"eve_num/I");
  Trun->Branch("nrebin",&nrebin,"nrebin/I");
  //Trun->Branch("threshold_u",&threshold_u,"threshold_u/F");
  //Trun->Branch("threshold_v",&threshold_v,"threshold_v/F");
  //Trun->Branch("threshold_w",&threshold_w,"threshold_w/F");
  Trun->Branch("time_offset",&time_offset,"time_offset/I");
  Trun->Fill();

  TTree *T_lf = new TTree("T_lf","T_lf");
  Int_t channel;
  T_lf->SetDirectory(file);
  T_lf->Branch("channel",&channel,"channel/I");


   
   //
    TH2F *hv_raw = new TH2F("hv_raw","hv_raw",nwire_v,-0.5+nwire_u,-0.5+nwire_u+nwire_v,nwindow_size*nrebin,0,nwindow_size*nrebin);
    for (size_t i=0;i!=nf_channelId->size();i++){
      int chid = nf_channelId->at(i);
      TH1S *htemp = (TH1S*)nf_wf->At(i);
      if (chid < nwire_u){
      }else if (chid < nwire_v+nwire_u){
	for (size_t j=0;j!=nwindow_size*nrebin;j++){
	  hv_raw->SetBinContent(chid-nwire_u+1,j+1,(htemp->GetBinContent(j+1)+nf_shift)*1.0/nf_scale);
	}
      }else if (chid < nwire_w+nwire_v+nwire_u){
      }
    }


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
	  std::cout << "W plane (shorted): " << i << " added to bad chanel list" << std::endl;
	  for (int j=0;j!=hw_decon->GetNbinsY();j++){
	    hw_decon->SetBinContent(i+1,j+1,0);
	    hw_decon_g->SetBinContent(i+1,j+1,0);
	  }
	}
      }
    }
   
    
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

  //  std::cout << "Xin1 " << " " << uplane_rms.size() << " " << std::endl;  
 

  int tpc_status = WireCell2dToy::Noisy_Event_ID(hu_decon, hv_decon, hw_decon, uplane_rms, vplane_rms, wplane_rms, uplane_map, vplane_map, wplane_map, hu_decon_g, hv_decon_g, hw_decon_g, nrebin, hv_raw, true);
  WireCell2dToy::Organize_Dead_Channels(uplane_map, vplane_map, wplane_map, hv_raw->GetNbinsY()-1,nrebin);
  
  }


    // loop through U/V/W plane to disable the bad channels completely
  for (auto it = uplane_map.begin(); it!=uplane_map.end(); it++){
    int ch = it->first;
    int start = it->second.first/nrebin;
    int end = it->second.second/nrebin;
    for (int j=start; j!=end+1;j++){
      hu_decon->SetBinContent(ch+1,j+1,0);
    }
  }
  for (auto it = vplane_map.begin(); it!=vplane_map.end(); it++){
    int ch = it->first;
    int start = it->second.first/nrebin;
    int end = it->second.second/nrebin;
    for (int j=start; j!=end+1;j++){
      hv_decon->SetBinContent(ch+1,j+1,0);
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
    }
  }

  

  TTree *T_bad = new TTree("T_bad","T_bad");
  // Int_t chid, plane;
  // Int_t start_time,end_time;
  T_bad->Branch("chid",&chid,"chid/I");
  T_bad->Branch("plane",&plane,"plane/I");
  T_bad->Branch("start_time",&start_time,"start_time/I");
  T_bad->Branch("end_time",&end_time,"end_time/I");
  T_bad->SetDirectory(file);

  // for (size_t i=0; i!=badChannel->size();i++){
  //   chid = badChannel->at(i);
  //   if (chid < nwire_u){
  //     plane = 0;
  //   }else if (chid < nwire_u + nwire_v){
  //     plane = 1;
  //   }else{
  //     plane = 2;
  //   }
  //   start_time = badBegin->at(i);
  //   end_time = badEnd->at(i);
  //   T_bad->Fill();
  // }
  
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

   hu_decon->SetDirectory(file);
   hv_decon->SetDirectory(file);
   hw_decon->SetDirectory(file);

   hu_decon_g->SetDirectory(file);
   hv_decon_g->SetDirectory(file);
   hw_decon_g->SetDirectory(file);

   hv_raw->SetDirectory(file);

    TH1I *hu_threshold = new TH1I("hu_threshold","hu_threshold",nwire_u,-0.5,-0.5+nwire_u);
  TH1I *hv_threshold = new TH1I("hv_threshold","hv_threshold",nwire_v,-0.5+nwire_u,-0.5+nwire_u+nwire_v);
  TH1I *hw_threshold = new TH1I("hw_threshold","hw_threshold",nwire_w,-0.5+nwire_u+nwire_v,-0.5+nwire_u+nwire_v+nwire_w);

  for (Int_t i=0;i!=uplane_rms.size();i++){
    hu_threshold->SetBinContent(i+1,uplane_rms.at(i)*3.0 * nrebin);
  }
  for (Int_t i=0;i!=vplane_rms.size();i++){
    hv_threshold->SetBinContent(i+1,vplane_rms.at(i)*3.0 * nrebin);
  }
  for (Int_t i=0;i!=wplane_rms.size();i++){
    hw_threshold->SetBinContent(i+1,wplane_rms.at(i)*3.0 * nrebin);
  }
  hu_threshold->SetDirectory(file);
  hv_threshold->SetDirectory(file);
  hw_threshold->SetDirectory(file);
  

   
  file->Write();
  file->Close();
  
  
}
