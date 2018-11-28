#include "WireCellNav/DetectorGDS.h"
#include "WireCellNav/GeomDataSource.h"
#include "WireCellNav/DetGenerativeFDS.h"
#include "WireCellNav/FrameDataSource.h"
#include "WireCell2dToy/ToyDepositor.h"
#include "WireCellSst/Util.h"
#include "WireCell2dToy/ToySignalSimu.h"
#include "WireCell2dToy/ToySignalSimuTrue.h"
#include "WireCell2dToy/ToySignalGaus.h"
#include "WireCell2dToy/ToySignalWien.h"
#include "WireCellSst/ToyuBooNESliceDataSource.h"
#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/TotalTiling.h"
#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/TruthToyTiling.h"
#include "WireCellData/MergeGeomCell.h"

#include "WireCell2dToy/ToyEventDisplay.h"

#include "WireCell2dToy/ToyMatrix.h"
#include "WireCell2dToy/ToyMatrixExclusive.h"
#include "WireCell2dToy/ToyMatrixKalman.h"
#include "WireCell2dToy/ToyMatrixIterate.h"
#include "WireCell2dToy/ToyMatrixMarkov.h"
#include "WireCell2dToy/ToyMetric.h"
#include "WireCellData/TPCParams.h"
#include "WireCellData/Singleton.h"
#include "WireCellData/GeomCluster.h"
#include "WireCellSst/MCTruth.h"

#include "WireCellSignal/ElectronicsConfig.h"
#include "WireCellSignal/ConvolutedResponse.h"
#include "WireCellSignal/DetGenerativeFDS.h"
//#include "WireCellSignal/GenerativeFDS.h"


#include "WireCellNav/FrameDataSource.h"
#include "WireCellNav/GenerativeFDS.h"
#include "WireCellNav/PepperDepositor.h"
#include "WireCell2dToy/ToySignalSimu.h"
#include "WireCell2dToy/DataSignalGaus.h"
//#include "WireCell2dToy/SimuSignalWien_ROI.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TBenchmark.h"
#include "TH1F.h"
#include "TLine.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph2D.h"
#include "TColor.h"
#include "TStyle.h"
#include "TMath.h"
#include "TRandom3.h"

#include <iostream>
#include <string>
#include <vector>

//#define MAX_TRACKS 10
using namespace WireCellSignal;
using namespace WireCell;
using namespace std;

int main(int argc, char* argv[])
{

  if (argc < 3){
    cerr << "usage: wire-cell-dune_2Dsignal /path/to/ChannelWireGeometry.txt /path/to/celltree.root eve_num" << endl;
    return 1;
  }

  int two_plane = 0;
  int save_file = 0;
  int nt_off1 = 0;
  int nt_off2 = 0;

  /*
  WireCellSst::GeomDataSource gds(argv[1]);
  std::vector<double> ex = gds.extent();
  cerr << "Extent: "
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
  */
  
  double pitchU, pitchV, pitchW;
  pitchU = 4.667*units::mm;
  pitchV = 4.667*units::mm;
  pitchW = 4.79*units::mm;
  
  //build GDS ... 
  DetectorGDS gds;
  gds.set_ncryos(1);
  gds.set_napas(0,1);
  Vector center0(-0*units::cm, -300.05*units::cm, 115.318875*units::cm);
  Vector halves0(3.995355*units::cm, 299.5*units::cm, 115.318875*units::cm);
  gds.set_apa(0, 0, 35.71*units::deg, 35.71*units::deg, pitchU, pitchV, pitchW, center0, halves0);
  gds.buildGDS();

  
  float unit_dis = 1.119;  
  float toffset_1=1.647;
  float toffset_2=1.539+1.647;
  float toffset_3=0;

  const char* root_file = argv[2];
  const char* tpath = "/Event/Sim";
  int eve_num = atoi(argv[3]);

  int total_time_bin=9600;
  int recon_threshold = 2000;
  int frame_length = 3200;
  int max_events = 100;
  int nrebin = 6;

  /*
  TPCParams& mp = Singleton<TPCParams>::Instance();
  
  double pitch_u = gds.get_pitch(0,WirePlaneType_t(0));
  double pitch_v = gds.get_pitch(0,WirePlaneType_t(1));
  double pitch_w = gds.get_pitch(0,WirePlaneType_t(2));
  double time_slice_width = nrebin * unit_dis * 0.5 * units::mm;

  mp.set_pitch_u(pitch_u);
  mp.set_pitch_v(pitch_v);
  mp.set_pitch_w(pitch_w);
  mp.set_ts_width(time_slice_width);
  */
  //std::cout << "Singleton: " << mp.get_pitch_u() << " " << mp.get_pitch_v() << " " << mp.get_pitch_w() << " " << mp.get_ts_width() << std::endl;
  
  float threshold_u = 5.87819e+02 * 4.0;
  float threshold_v = 8.36644e+02 * 4.0;
  float threshold_w = 5.67974e+02 * 4.0;

  float threshold_ug = 755.96;
  float threshold_vg = 822.81;
  float threshold_wg = 510.84;
  
  int time_offset = 0;

  
  TFile *tfile = TFile::Open(root_file);
  TTree* sst = dynamic_cast<TTree*>(tfile->Get(tpath));

  int run_no, subrun_no, event_no;
  sst->SetBranchAddress("eventNo",&event_no);
  sst->SetBranchAddress("runNo",&run_no);
  sst->SetBranchAddress("subRunNo",&subrun_no);
  sst->GetEntry(eve_num);
  
   //save MC truth ...
  WireCellSst::MCTruth *mctruth = new WireCellSst::MCTruth(root_file);
  mctruth->GetEntry(eve_num);
  
  std::cout << "Energy: " << mctruth->mc_nu_mom[3] << " CC/NC: " << mctruth->mc_nu_ccnc << " Neutrino: " << mctruth->mc_nu_pdg << std::endl;
  cout << "Run No: " << mctruth->runNo << " " << mctruth->subrunNo << " " << mctruth->eventNo << endl;
  
  WireCell::FrameDataSource* fds = 0;
  fds = WireCellSst::make_fds(*tfile);
  if (!fds) {
    cerr << "ERROR: failed to get FDS from " << root_file << endl;
    return 1;
  }
  
  TH1::AddDirectory(kFALSE);
  
  WireCell::ToyDepositor *toydep;
  //WireCell::PepperDepositor *toydep;

  Double_t x_center = mctruth->mc_startXYZT[0][0];
  Double_t y_center = mctruth->mc_startXYZT[0][1];
  Double_t z_center = mctruth->mc_startXYZT[0][2];
  Double_t x_shift = -mctruth->mc_startXYZT[0][0] + 10;
  Double_t y_shift = -mctruth->mc_startXYZT[0][1] - 300;
  Double_t z_shift = -mctruth->mc_startXYZT[0][2] + 115;

  std::pair<float,float> xmm, ymm, zmm, qmm;
  xmm.first = 50*units::cm;
  xmm.second = 200*units::cm;
  ymm.first = -500*units::cm;
  ymm.second = -100*units::cm;
  zmm.first = 50*units::cm;
  zmm.second = 180*units::cm;
  qmm.first = 2000;
  qmm.second = 2001;
  //toydep = new WireCell::PepperDepositor(xmm, ymm, zmm, qmm, 1);

  toydep = new WireCell::ToyDepositor(fds,0,unit_dis,frame_length);
  cout << "Primary vertex is (" << mctruth->mc_startXYZT[0][0] << "," << mctruth->mc_startXYZT[0][1] << "," << mctruth->mc_startXYZT[0][2] <<")" << endl;
  toydep->translation(x_shift,y_shift,z_shift);
  mctruth->Rotate_Shift(0,0,0, 0,
  			x_shift,y_shift,z_shift);

  const PointValueVector& pvv = toydep->depositions(eve_num);

  std::cout << "Points deposited: " << pvv.size() << std::endl;

  WireCellSignal::ElectronicsConfig *conf = new WireCellSignal::ElectronicsConfig();
  WireCellSignal::ConvolutedResponse *fRsp = new WireCellSignal::ConvolutedResponse(conf, "./signal/dune.root", time_offset, toffset_1, toffset_2);
  //WireCellSignal::GenerativeFDS gfds(*toydep, gds, conf->NTDC(), max_events, 0.5*1.6*units::millimeter);
  WireCellSignal::DetGenerativeFDS gfds(*toydep, gds, conf->NTDC(), max_events, 0.5*1.6*units::millimeter);// 3m maximum drift distance
  gfds.SetResponseFunctions(fRsp);
  gfds.jump(eve_num);
  tfile->Close("R");
  delete tfile;
  //return 0;
  cout << "Simulate Raw WaveForm " << endl; 
  WireCell2dToy::ToySignalSimuFDS simu_fds(gfds,gds,*conf,max_events,toffset_1,toffset_2,1); // time offset among different planes for the time electrons travel among different planes
  simu_fds.jump(eve_num);

  //return 0;
  TFile *file = new TFile(Form("2D_display_%d_%d_%d.root",run_no,subrun_no,eve_num),"RECREATE");
  TH2F *hsig[4][100][3][2];
  for (int i=0; i<gds.ncryos(); ++i) {
    for (int j=0; j<gds.napa(i); ++j) {
      const WrappedGDS *apa_gds = gds.get_apaGDS(i,j);
      for (int p=0; p<3; ++p) {
	WirePlaneType_t plane = static_cast<WirePlaneType_t>(p);
	for (int f = 0; f < 2; ++f) {
	  int nwires = apa_gds->wires_in_plane(f,plane).size();
	  hsig[i][j][p][f] = new TH2F(Form("hsig_%d_%d_%d_%d", i, j, p, f),
				      Form("hsig_%d_%d_%d_%d", i, j, p, f), 
				      nwires, 0, nwires,
				      conf->NTDC(), 0, conf->NTDC());
	}
      }
    }
  }    
  
  const Frame& frame = simu_fds.get();
  size_t ntraces = frame.traces.size();    
  for (size_t ind=0; ind<ntraces; ++ind) {
    //std::cout<<ind<<std::endl;
    const Trace& trace = frame.traces[ind];
    int chid = trace.chid;
    const GeomWireSelection& wire_in_channel = gds.by_channel(chid);     
    int cryo = wire_in_channel.at(0)->cryo();
    int apa = wire_in_channel.at(0)->apa();
    WirePlaneType_t plane = wire_in_channel.at(0)->plane();
    for (int w=0; w<wire_in_channel.size(); ++w) {      
      int face = wire_in_channel.at(w)->face();
      const GeomWire *wire = gds.by_channel(chid).at(w);
      const GeomWireSelection& wires_in_plane = gds.get_apaGDS(cryo, apa)->wires_in_plane(face, plane);
      int wirebin = (int)(find(wires_in_plane.begin(), wires_in_plane.end(), wire)-wires_in_plane.begin());
      for (int t=0; t<conf->NTDC(); ++t) {
	hsig[cryo][apa][(int)plane][face]->SetBinContent(wirebin+1, t+1, trace.charge.at(t));
      }
    }
  }

  file->cd();
  for (int i=0; i<gds.ncryos(); ++i) {
    for (int j=0; j<gds.napa(i); ++j) {
      for (int p=0; p<3; ++p) {
	for (int f = 0; f < 2; ++f) {
	  if (hsig[i][j][p][f]!=NULL)hsig[i][j][p][f]->Write();
	  //delete hsig[i][j][p][f];
	}
      }
    }
  }
  
  file->Close();
}

