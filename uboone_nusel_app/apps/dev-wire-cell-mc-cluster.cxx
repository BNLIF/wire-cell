#include "WCPSst/GeomDataSource.h"
#include "WCPData/PR3DCluster.h"
#include "WCPData/SlimMergeGeomCell.h"
#include "WCPData/TPCParams.h"
#include "WCPData/Singleton.h"
#include "WCPData/ToyCTPointCloud.h"

#include "WCP2dToy/ExecMon.h"
#include "WCP2dToy/CalcPoints.h"
#include "WCP2dToy/ToyClustering.h"
//#include "WCP2dToy/uBooNE_light_reco.h"
#include "WCP2dToy/ToyLightReco.h"

#include "WCP2dToy/ToyMatching.h"
#include "WCP2dToy/ToyFiducial.h"
#include "WCP2dToy/ImprovePR3DCluster.h"
#include "WCP2dToy/ExamineBundles.h"

#include "TH1F.h"
#include "TH2F.h"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"

using namespace WCP;
using namespace std;

int main(int argc, char* argv[])
{
  if (argc < 3) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/mc_imaging.root " << endl;
    return 1;
  }
  TH1::AddDirectory(kFALSE);
  
  int flag_pos_corr = 0; // correct X position after matching ...
  int bee_debug = 0;
  
  for(Int_t i = 1; i != argc; i++){
    switch(argv[i][1]){
    case 'c':
      flag_pos_corr = atoi(&argv[i][2]); 
      break;
    case 'b':
      bee_debug = atoi(&argv[i][2]);
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
  
  
  
  // test geometry ...
  const GeomWire *uwire = gds.by_planeindex(WirePlaneType_t(0),0);
  const GeomWire *vwire = gds.by_planeindex(WirePlaneType_t(1),0);
  const GeomWire *wwire = gds.by_planeindex(WirePlaneType_t(2),0);
  double first_u_dis = gds.wire_dist(*uwire) ; // first U wire center ...
  double first_v_dis = gds.wire_dist(*vwire) ; // first V wire center ...
  double first_w_dis = gds.wire_dist(*wwire) ; // first W wire center ... 
  
  
  TString filename = argv[2];
  TFile *file = new TFile(filename);
  TTree *Trun = (TTree*)file->Get("Trun");
  
  int run_no, subrun_no, event_no;
  int time_offset;
  int nrebin;
  int frame_length;
  int eve_num;
  float unit_dis;
  
  std::vector<int> *timesliceId = new std::vector<int>;
  std::vector<std::vector<int>> *timesliceChannel = new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *raw_charge = new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *raw_charge_err = new std::vector<std::vector<int>>;
  
  Trun->SetBranchAddress("eventNo",&event_no);
  Trun->SetBranchAddress("runNo",&run_no);
  Trun->SetBranchAddress("subRunNo",&subrun_no);
  Trun->SetBranchAddress("unit_dis",&unit_dis);
  Trun->SetBranchAddress("frame_length",&frame_length);
  Trun->SetBranchAddress("eve_num",&eve_num);
  Trun->SetBranchAddress("nrebin",&nrebin);
  Trun->SetBranchAddress("time_offset",&time_offset);
  
  Trun->SetBranchAddress("timesliceId",&timesliceId);
  Trun->SetBranchAddress("timesliceChannel",&timesliceChannel);
  Trun->SetBranchAddress("raw_charge",&raw_charge);
  Trun->SetBranchAddress("raw_charge_err",&raw_charge_err);
  
  
  double triggerTime;
  Trun->SetBranchAddress("triggerTime",&triggerTime); 
  Trun->GetEntry(0);
  
    
  // define singleton ... 
  TPCParams& mp = Singleton<TPCParams>::Instance();
  
  double pitch_u = gds.pitch(WirePlaneType_t(0));
  double pitch_v = gds.pitch(WirePlaneType_t(1));
  double pitch_w = gds.pitch(WirePlaneType_t(2));
  double time_slice_width = nrebin * unit_dis * 0.5 * units::mm;
  
  double angle_u = gds.angle(WirePlaneType_t(0));
  double angle_v = gds.angle(WirePlaneType_t(1));
  double angle_w = gds.angle(WirePlaneType_t(2));
  
  //std::cout << angle_u << " " << angle_v << " " << angle_w << std::endl;
  
  mp.set_pitch_u(pitch_u);
  mp.set_pitch_v(pitch_v);
  mp.set_pitch_w(pitch_w);
  mp.set_angle_u(angle_u);
  mp.set_angle_v(angle_v);
  mp.set_angle_w(angle_w);
  mp.set_ts_width(time_slice_width);
  mp.set_first_u_dis(first_u_dis);
  mp.set_first_v_dis(first_v_dis);
  mp.set_first_w_dis(first_w_dis);
  
  
  
  
  
  std::map<int,std::pair<double,double>> dead_u_index;
  std::map<int,std::pair<double,double>> dead_v_index;
  std::map<int,std::pair<double,double>> dead_w_index;
  // load mcell
  
  TTree *TC = (TTree*)file->Get("TC");
  std::vector<int> *cluster_id_vec = new std::vector<int>;
  std::vector<int> *time_slice_vec = new std::vector<int>;
  std::vector<double> *q_vec = new std::vector<double>;
  std::vector<double> *uq_vec = new std::vector<double>;
  std::vector<double> *vq_vec = new std::vector<double>;
  std::vector<double> *wq_vec = new std::vector<double>;
  std::vector<double> *udq_vec = new std::vector<double>;
  std::vector<double> *vdq_vec = new std::vector<double>;
  std::vector<double> *wdq_vec = new std::vector<double>;
  
  std::vector<int> *nwire_u_vec = new  std::vector<int>;
  std::vector<int> *nwire_v_vec = new  std::vector<int>;
  std::vector<int> *nwire_w_vec = new  std::vector<int>;
  std::vector<int> *flag_u_vec = new  std::vector<int>;
  std::vector<int> *flag_v_vec = new  std::vector<int>;
  std::vector<int> *flag_w_vec = new  std::vector<int>;
  
  std::vector<std::vector<int>> *wire_index_u_vec = new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *wire_index_v_vec = new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *wire_index_w_vec = new std::vector<std::vector<int>>;
  std::vector<std::vector<double>> *wire_charge_u_vec = new std::vector<std::vector<double>>;
  std::vector<std::vector<double>> *wire_charge_v_vec = new std::vector<std::vector<double>>;
  std::vector<std::vector<double>> *wire_charge_w_vec = new std::vector<std::vector<double>>;
  std::vector<std::vector<double>> *wire_charge_err_u_vec = new std::vector<std::vector<double>>;
  std::vector<std::vector<double>> *wire_charge_err_v_vec = new std::vector<std::vector<double>>;
  std::vector<std::vector<double>> *wire_charge_err_w_vec = new std::vector<std::vector<double>>;
  

  TC->SetBranchAddress("cluster_id",&cluster_id_vec);
  TC->SetBranchAddress("time_slice",&time_slice_vec);
  TC->SetBranchAddress("q",&q_vec);
  TC->SetBranchAddress("uq",&uq_vec);
  TC->SetBranchAddress("vq",&vq_vec);
  TC->SetBranchAddress("wq",&wq_vec);
  TC->SetBranchAddress("udq",&udq_vec);
  TC->SetBranchAddress("vdq",&vdq_vec);
  TC->SetBranchAddress("wdq",&wdq_vec);
  TC->SetBranchAddress("nwire_u",&nwire_u_vec);
  TC->SetBranchAddress("nwire_v",&nwire_v_vec);
  TC->SetBranchAddress("nwire_w",&nwire_w_vec);
  TC->SetBranchAddress("flag_u",&flag_u_vec);
  TC->SetBranchAddress("flag_v",&flag_v_vec);
  TC->SetBranchAddress("flag_w",&flag_w_vec);
  TC->SetBranchAddress("wire_index_u",&wire_index_u_vec);
  TC->SetBranchAddress("wire_index_v",&wire_index_v_vec);
  TC->SetBranchAddress("wire_index_w",&wire_index_w_vec);
  TC->SetBranchAddress("wire_charge_u",&wire_charge_u_vec);
  TC->SetBranchAddress("wire_charge_v",&wire_charge_v_vec);
  TC->SetBranchAddress("wire_charge_w",&wire_charge_w_vec);
  TC->SetBranchAddress("wire_charge_err_u",&wire_charge_err_u_vec);
  TC->SetBranchAddress("wire_charge_err_v",&wire_charge_err_v_vec);
  TC->SetBranchAddress("wire_charge_err_w",&wire_charge_err_w_vec);
  
  
  //load mcell
  
  TTree *TDC = (TTree*)file->Get("TDC");
  std::vector<int> *ntime_slice_vec = new std::vector<int>;
  std::vector<std::vector<int>> *time_slices_vec = new std::vector<std::vector<int>>;
  
  TDC->SetBranchAddress("cluster_id",&cluster_id_vec);
  TDC->SetBranchAddress("ntime_slice",&ntime_slice_vec);
  TDC->SetBranchAddress("time_slice",&time_slices_vec);
  
  TDC->SetBranchAddress("nwire_u",&nwire_u_vec);
  TDC->SetBranchAddress("nwire_v",&nwire_v_vec);
  TDC->SetBranchAddress("nwire_w",&nwire_w_vec);
  TDC->SetBranchAddress("flag_u",&flag_u_vec);
  TDC->SetBranchAddress("flag_v",&flag_v_vec);
  TDC->SetBranchAddress("flag_w",&flag_w_vec);
  TDC->SetBranchAddress("wire_index_u",&wire_index_u_vec);
  TDC->SetBranchAddress("wire_index_v",&wire_index_v_vec);
  TDC->SetBranchAddress("wire_index_w",&wire_index_w_vec);
  
  
  WCP2dToy::ToyFiducial *fid = new WCP2dToy::ToyFiducial(3,800,-first_u_dis/pitch_u, -first_v_dis/pitch_v, -first_w_dis/pitch_w,
								   1./time_slice_width, 1./pitch_u, 1./pitch_v, 1./pitch_w, // slope
								   angle_u,angle_v,angle_w,// angle
								   3*units::cm, 117*units::cm, -116*units::cm, 0*units::cm, 1037*units::cm, 0*units::cm, 256*units::cm);
  
  // load cells ... 
  GeomCellSelection mcells;
  PR3DClusterSelection live_clusters;
  PR3DClusterSelection dead_clusters;
  PR3DCluster *cluster;
  int prev_cluster_id=-1;
  int ident = 0;
  TC->GetEntry(0);
  for (int i=0;i!=cluster_id_vec->size();i++){
    int cluster_id = cluster_id_vec->at(i);
    SlimMergeGeomCell *mcell = new SlimMergeGeomCell(ident);
    int time_slice = time_slice_vec->at(i);
    std::vector<int> wire_index_u = wire_index_u_vec->at(i);
    std::vector<int> wire_index_v = wire_index_v_vec->at(i);
    std::vector<int> wire_index_w = wire_index_w_vec->at(i);
    std::vector<double> wire_charge_u = wire_charge_u_vec->at(i);
    std::vector<double> wire_charge_v = wire_charge_v_vec->at(i);
    std::vector<double> wire_charge_w = wire_charge_w_vec->at(i);
    std::vector<double> wire_charge_err_u = wire_charge_err_u_vec->at(i);
    std::vector<double> wire_charge_err_v = wire_charge_err_v_vec->at(i);
    std::vector<double> wire_charge_err_w = wire_charge_err_w_vec->at(i);
    
    mcell->SetTimeSlice(time_slice);
    
    mcell->set_uq(uq_vec->at(i));
    mcell->set_vq(vq_vec->at(i));
    mcell->set_wq(wq_vec->at(i));
    
    mcell->set_udq(udq_vec->at(i));
    mcell->set_vdq(vdq_vec->at(i));
    mcell->set_wdq(wdq_vec->at(i));
    
    mcell->set_q(q_vec->at(i));
    
    double temp_x = (time_slice*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) * units::cm;
    
    if (flag_u_vec->at(i)==0){
      mcell->add_bad_planes(WirePlaneType_t(0));
      for (int j=0;j!=nwire_u_vec->at(i);j++){
	if (dead_u_index.find(wire_index_u.at(j))==dead_u_index.end()){
	  dead_u_index[wire_index_u.at(j)] = std::make_pair(temp_x-0.1*units::cm,temp_x+0.1*units::cm);
	}else{
	  if (temp_x-0.1*units::cm < dead_u_index[wire_index_u.at(j)].first){
	    dead_u_index[wire_index_u.at(j)].first = temp_x - 0.1*units::cm;
	  }else if (temp_x+0.1*units::cm > dead_u_index[wire_index_u.at(j)].second){
	    dead_u_index[wire_index_u.at(j)].second = temp_x + 0.1*units::cm;
	  }
	}
	
      }
    }
    if (flag_v_vec->at(i)==0){
      mcell->add_bad_planes(WirePlaneType_t(1));
      for (int j=0;j!=nwire_v_vec->at(i);j++){
	if (dead_v_index.find(wire_index_v.at(j))==dead_v_index.end()){
	  dead_v_index[wire_index_v.at(j)] = std::make_pair(temp_x-0.1*units::cm,temp_x+0.1*units::cm);
	}else{
	  if (temp_x-0.1*units::cm < dead_v_index[wire_index_v.at(j)].first){
	    dead_v_index[wire_index_v.at(j)].first = temp_x-0.1*units::cm;
	  }else if (temp_x+0.1*units::cm > dead_v_index[wire_index_v.at(j)].second){
	    dead_v_index[wire_index_v.at(j)].second = temp_x + 0.1*units::cm;
	  }
	}
      }
    }
    if (flag_w_vec->at(i)==0){
      mcell->add_bad_planes(WirePlaneType_t(2));
      for (int j=0;j!=nwire_w_vec->at(i);j++){
	if (dead_w_index.find(wire_index_w.at(j))==dead_w_index.end()){
	  dead_w_index[wire_index_w.at(j)] = std::make_pair(temp_x-0.1*units::cm,temp_x+0.1*units::cm);
	}else{
	  if (temp_x-0.1*units::cm < dead_w_index[wire_index_w.at(j)].first){
	    dead_w_index[wire_index_w.at(j)].first = temp_x-0.1*units::cm;
	  }else if (temp_x+0.1*units::cm > dead_w_index[wire_index_w.at(j)].second){
	    dead_w_index[wire_index_w.at(j)].second = temp_x + 0.1*units::cm;
	  }
	}
      }
    }
    for (int j=0;j!=nwire_u_vec->at(i);j++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(0),wire_index_u.at(j));
      mcell->AddWire(wire,WirePlaneType_t(0),wire_charge_u.at(j),wire_charge_err_u.at(j));
    }
    for (int j=0;j!=nwire_v_vec->at(i);j++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(1),wire_index_v.at(j));
      mcell->AddWire(wire,WirePlaneType_t(1),wire_charge_v.at(j),wire_charge_err_v.at(j));
    }
    for (int j=0;j!=nwire_w_vec->at(i);j++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(2),wire_index_w.at(j));
      mcell->AddWire(wire,WirePlaneType_t(2),wire_charge_w.at(j),wire_charge_err_w.at(j));
    }
    mcells.push_back(mcell);
    
    if (cluster_id != prev_cluster_id){
      cluster = new PR3DCluster(cluster_id);
      live_clusters.push_back(cluster);
    }
    cluster->AddCell(mcell,time_slice);
    
    prev_cluster_id = cluster_id;
    ident++;
  }
  //  std::cout << live_clusters.size() << std::endl;
  
  prev_cluster_id = -1;
  // TDC
  cluster_id_vec->clear();
  wire_index_u_vec->clear();
  wire_index_v_vec->clear();
  wire_index_w_vec->clear();
  nwire_u_vec->clear();
  nwire_v_vec->clear();
  nwire_w_vec->clear();
  flag_u_vec->clear();
  flag_v_vec->clear();
  flag_w_vec->clear();
  
  TDC->GetEntry(0);
  for (int i=0;i!=cluster_id_vec->size();i++){
    int cluster_id = cluster_id_vec->at(i);
    SlimMergeGeomCell *mcell = new SlimMergeGeomCell(ident);
    std::vector<int> time_slices = time_slices_vec->at(i);
    std::vector<int> wire_index_u = wire_index_u_vec->at(i);
    std::vector<int> wire_index_v = wire_index_v_vec->at(i);
    std::vector<int> wire_index_w = wire_index_w_vec->at(i);
    
    mcell->SetTimeSlice(time_slices.at(0));
    
    double temp_x1 = (time_slices.front()*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) * units::cm;
    double temp_x2 = (time_slices.back()*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) * units::cm;
    // std::cout << temp_x1/units::cm << " " << temp_x2/units::cm << std::endl;
    
    if (flag_u_vec->at(i)==0){
      mcell->add_bad_planes(WirePlaneType_t(0));
      for (int j=0;j!=nwire_u_vec->at(i);j++){
	if (dead_u_index.find(wire_index_u.at(j))==dead_u_index.end()){
	  dead_u_index[wire_index_u.at(j)] = std::make_pair(temp_x1-0.1*units::cm,temp_x2+0.1*units::cm);
	}else{
	  if (temp_x1-0.1*units::cm < dead_u_index[wire_index_u.at(j)].first){
	    dead_u_index[wire_index_u.at(j)].first = temp_x1 - 0.1*units::cm;
	  }else if (temp_x2+0.1*units::cm > dead_u_index[wire_index_u.at(j)].second){
	    dead_u_index[wire_index_u.at(j)].second = temp_x2 + 0.1*units::cm;
	  }
	}
      }
    }
    if (flag_v_vec->at(i)==0){
      mcell->add_bad_planes(WirePlaneType_t(1));
      for (int j=0;j!=nwire_v_vec->at(i);j++){
	if (dead_v_index.find(wire_index_v.at(j))==dead_v_index.end()){
	  dead_v_index[wire_index_v.at(j)] = std::make_pair(temp_x1-0.1*units::cm,temp_x2+0.1*units::cm);
	}else{
	  if (temp_x1-0.1*units::cm < dead_v_index[wire_index_v.at(j)].first){
	    dead_v_index[wire_index_v.at(j)].first = temp_x1 - 0.1*units::cm;
	  }else if (temp_x2+0.1*units::cm > dead_v_index[wire_index_v.at(j)].second){
	    dead_v_index[wire_index_v.at(j)].second = temp_x2 + 0.1*units::cm;
	  }
	}
      }
    }
    if (flag_w_vec->at(i)==0){
      mcell->add_bad_planes(WirePlaneType_t(2));
      for (int j=0;j!=nwire_w_vec->at(i);j++){
	if (dead_w_index.find(wire_index_w.at(j))==dead_w_index.end()){
	  dead_w_index[wire_index_w.at(j)] = std::make_pair(temp_x1-0.1*units::cm,temp_x2+0.1*units::cm);
	}else{
	  if (temp_x1-0.1*units::cm < dead_w_index[wire_index_w.at(j)].first){
	    dead_w_index[wire_index_w.at(j)].first = temp_x1 - 0.1*units::cm;
	  }else if (temp_x2+0.1*units::cm > dead_w_index[wire_index_w.at(j)].second){
	    dead_w_index[wire_index_w.at(j)].second = temp_x2 + 0.1*units::cm;
	  }
	}
      }
    }
    for (int j=0;j!=nwire_u_vec->at(i);j++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(0),wire_index_u.at(j));
      mcell->AddWire(wire,WirePlaneType_t(0));
    }
    for (int j=0;j!=nwire_v_vec->at(i);j++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(1),wire_index_v.at(j));
      mcell->AddWire(wire,WirePlaneType_t(1));
    }
    for (int j=0;j!=nwire_w_vec->at(i);j++){
      const GeomWire *wire = gds.by_planeindex(WirePlaneType_t(2),wire_index_w.at(j));
      mcell->AddWire(wire,WirePlaneType_t(2));
    }
    mcells.push_back(mcell);
    
    if (cluster_id!=prev_cluster_id){
      cluster = new PR3DCluster(cluster_id);
      dead_clusters.push_back(cluster);
    }
    for (int j=0;j!=ntime_slice_vec->at(i);j++){
      cluster->AddCell(mcell,time_slices.at(j));
    }
    fid->AddDeadRegion(mcell,time_slices);
    
    prev_cluster_id=cluster_id;
    ident++;
  }
  
  for (size_t i=0;i!=dead_clusters.size();i++){
    dead_clusters.at(i)->Remove_duplicated_mcells();
    // std::cout << dead_clusters.at(i)->get_cluster_id() << " " 
    // 	       << dead_clusters.at(i)->get_num_mcells() << " "
    // 	       << dead_clusters.at(i)->get_num_time_slices() << std::endl;
  }
  
  if (timesliceId->size()==0){
    std::cout << "No points! Quit! " << std::endl;
    return 0;
  }
  
  // veto 16 channels in U ...
  // for (int i=2080; i!=2096;i++){
  //   if (dead_u_index.find(i)==dead_u_index.end()){
  //     dead_u_index[i] = std::make_pair((0*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) * units::cm-0.1*units::cm, (2400*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) * units::cm+0.1*units::cm);
  //   }
  // }
  
  // Load T_ch_bad tree ...
  TTree *T_bad_ch = (TTree*)file->Get("T_bad_ch");
  if (T_bad_ch!=0){
    Int_t chid, plane;
    Int_t start_time,end_time;
    T_bad_ch->Branch("chid",&chid,"chid/I");
    T_bad_ch->Branch("plane",&plane,"plane/I");
    T_bad_ch->Branch("start_time",&start_time,"start_time/I");
    T_bad_ch->Branch("end_time",&end_time,"end_time/I");
    
    for (int i=0;i!=T_bad_ch->GetEntries();i++){
      T_bad_ch->GetEntry(i);
      double temp_x1 = (start_time/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) * units::cm;
      double temp_x2 = (end_time/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) * units::cm;
      if (plane==1){
	chid -= 2400;
      }else if (plane==2){
	chid -=4800;
      }
      if (plane==0){
	if (dead_u_index.find(chid)==dead_u_index.end()){
	  dead_u_index[chid] = std::make_pair(temp_x1-0.1*units::cm,temp_x2+0.1*units::cm);
	}else{
	  if (temp_x1-0.1*units::cm < dead_u_index[chid].first){
	    dead_u_index[chid].first = temp_x1 - 0.1*units::cm;
	  }else if (temp_x2+0.1*units::cm > dead_u_index[chid].second){
	    dead_u_index[chid].second = temp_x2 + 0.1*units::cm;
	  }
	}
      }else if (plane==1){
	if (dead_v_index.find(chid)==dead_v_index.end()){
	  dead_v_index[chid] = std::make_pair(temp_x1-0.1*units::cm,temp_x2+0.1*units::cm);
	  //std::cout << plane << " " << chid << std::endl;
	}else{
	  if (temp_x1-0.1*units::cm < dead_v_index[chid].first){
	    dead_v_index[chid].first = temp_x1 - 0.1*units::cm;
	    //std::cout << plane << " a " << chid << std::endl;
	  }else if (temp_x2+0.1*units::cm > dead_v_index[chid].second){
	    dead_v_index[chid].second = temp_x2 + 0.1*units::cm;
	    //std::cout << plane << " b " << chid << std::endl;
	  }
	}
      }else if (plane==2){
	if (dead_w_index.find(chid)==dead_w_index.end()){
	  dead_w_index[chid] = std::make_pair(temp_x1-0.1*units::cm,temp_x2+0.1*units::cm);
	  //std::cout << plane << " " << chid << std::endl;
	}else{
	  if (temp_x1-0.1*units::cm < dead_w_index[chid].first){
	    dead_w_index[chid].first = temp_x1 - 0.1*units::cm;
	    //std::cout << plane << " a " << chid << std::endl;
	  }else if (temp_x2+0.1*units::cm > dead_w_index[chid].second){
	    dead_w_index[chid].second = temp_x2 + 0.1*units::cm;
	    //std::cout << plane << " b " << chid << std::endl;
	  }
	}
      }
    }
    
  }
  
  
  cout << em("load clusters from file") << endl;
  
  // Start to add X, Y, Z points
  // form boundaries for the bad cells ... 
  for (size_t j = 0; j!= dead_clusters.size(); j++){
    WCP2dToy::calc_boundary_points_dead(gds,dead_clusters.at(j));
  }
  // form sampling points for the normal cells ...
  DynamicToyPointCloud global_point_cloud(angle_u,angle_v,angle_w);
  for (size_t i=0; i!=live_clusters.size();i++){
    WCP2dToy::calc_sampling_points(gds,live_clusters.at(i),nrebin, frame_length, unit_dis);
    live_clusters.at(i)->Create_point_cloud();
    global_point_cloud.AddPoints(live_clusters.at(i),0);
    //live_clusters.at(i)->Calc_PCA();
  }
  cout << em("Add X, Y, Z points") << std::endl;
  
  
  // create global CT point cloud ...
  double first_t_dis = live_clusters.at(0)->get_mcells().front()->GetTimeSlice()*time_slice_width - live_clusters.at(0)->get_mcells().front()->get_sampling_points().front().x;
  double offset_t = first_t_dis/time_slice_width;
  
  // test the fiducial volume cut 
  fid->set_offset_t(offset_t);
  
  ToyCTPointCloud ct_point_cloud(0,2399,2400,4799,4800,8255, // channel range
				 offset_t, -first_u_dis/pitch_u, -first_v_dis/pitch_v, -first_w_dis/pitch_w, // offset
				 1./time_slice_width, 1./pitch_u, 1./pitch_v, 1./pitch_w, // slope
				 angle_u,angle_v,angle_w// angle
				 );
  ct_point_cloud.AddPoints(timesliceId,timesliceChannel,raw_charge,raw_charge_err);
  ct_point_cloud.AddDeadChs(dead_u_index, dead_v_index, dead_w_index);
  ct_point_cloud.build_kdtree_index();
   
  // finish creating global CT point cloud 
  
  
  
  // WCP2dToy::Clustering_live_dead(live_clusters, dead_clusters);
  // cerr << em("Clustering live and dead clusters") << std::endl;
  
  std::map<PR3DCluster*,std::vector<std::pair<PR3DCluster*,double>>> group_clusters = WCP2dToy::Clustering_jump_gap_cosmics(live_clusters, dead_clusters,dead_u_index, dead_v_index, dead_w_index, global_point_cloud, ct_point_cloud);
  cout << em("Clustering to jump gap in cosmics") << std::endl;

  // form a global map with the current map information
  std::map<int,std::map<const GeomWire*, SMGCSelection > > global_wc_map;
  for (size_t i=0; i!=live_clusters.size();i++){
    PR3DCluster *cluster = live_clusters.at(i);
    SMGCSelection& mcells = cluster->get_mcells();
    for (auto it = mcells.begin(); it!= mcells.end(); it++){
      SlimMergeGeomCell *mcell = *it;
      int time_slice = mcell->GetTimeSlice();
      if (global_wc_map.find(time_slice)==global_wc_map.end()){
	std::map<const GeomWire*, SMGCSelection> temp_wc_map;
	global_wc_map[time_slice] = temp_wc_map;
      }
      std::map<const GeomWire*, SMGCSelection>& timeslice_wc_map = global_wc_map[time_slice];
      
      GeomWireSelection& uwires = mcell->get_uwires();
      GeomWireSelection& vwires = mcell->get_vwires();
      GeomWireSelection& wwires = mcell->get_wwires();
      std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();
      if (find(bad_planes.begin(),bad_planes.end(),WirePlaneType_t(0))==bad_planes.end()){
	for (int j=0;j!=uwires.size();j++){
	  const GeomWire *wire = uwires.at(j);
	  if (timeslice_wc_map.find(wire)==timeslice_wc_map.end()){
	    SMGCSelection temp_mcells;
	    temp_mcells.push_back(mcell);
	    timeslice_wc_map[wire] = temp_mcells;
	  }else{
	    timeslice_wc_map[wire].push_back(mcell);
	  }
	}
      }
      if (find(bad_planes.begin(),bad_planes.end(),WirePlaneType_t(1))==bad_planes.end()){
	for (int j=0;j!=vwires.size();j++){
	  const GeomWire *wire = vwires.at(j);
	  if (timeslice_wc_map.find(wire)==timeslice_wc_map.end()){
	    SMGCSelection temp_mcells;
	    temp_mcells.push_back(mcell);
	    timeslice_wc_map[wire] = temp_mcells;
	  }else{
	    timeslice_wc_map[wire].push_back(mcell);
	  }
	}
      }
      if (find(bad_planes.begin(),bad_planes.end(),WirePlaneType_t(2))==bad_planes.end()){
	for (int j=0;j!=wwires.size();j++){
	  const GeomWire *wire = wwires.at(j);
	  if (timeslice_wc_map.find(wire)==timeslice_wc_map.end()){
	    SMGCSelection temp_mcells;
	    temp_mcells.push_back(mcell);
	    timeslice_wc_map[wire] = temp_mcells;
	  }else{
	    timeslice_wc_map[wire].push_back(mcell);
	  }
	}
      }
    }
  }
  
  
  //std::cout << group_clusters.size() << std::endl;
 
  live_clusters.clear();
  for (auto it = group_clusters.begin(); it!=group_clusters.end(); it++){
    live_clusters.push_back(it->first);
    for (auto it1 = it->second.begin(); it1!=it->second.end();it1++){
      live_clusters.push_back((*it1).first);
    }
  }
  
  for (size_t i=0;i!=live_clusters.size();i++){
    live_clusters.at(i)->Create_graph(ct_point_cloud);
    std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = live_clusters.at(i)->get_highest_lowest_wcps();
    live_clusters.at(i)->dijkstra_shortest_paths(wcps.first);
    live_clusters.at(i)->cal_shortest_path(wcps.second);

    live_clusters.at(i)->collect_charge_trajectory(ct_point_cloud);
  }


   
  cerr << em("Create Graph in all clusters") << std::endl;
   
  TFile *file1 = new TFile(Form("match_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");
  TTree *T_match = new TTree("T_match","T_match");
  T_match->SetDirectory(file1);
  Int_t ncluster=live_clusters.at(0)->get_cluster_id();
  Int_t flash_id=0;
  Double_t strength=1;
  Double_t pe_pred[32]={1,1,1,1,1,1,1,1,1,1,
			1,1,1,1,1,1,1,1,1,1,
			1,1,1,1,1,1,1,1,1,1,
			1,1};
  Double_t pe_meas[32]={1,1,1,1,1,1,1,1,1,1,
			1,1,1,1,1,1,1,1,1,1,
			1,1,1,1,1,1,1,1,1,1,
			1,1};
  Double_t pe_meas_err[32]={1,1,1,1,1,1,1,1,1,1,
			1,1,1,1,1,1,1,1,1,1,
			1,1,1,1,1,1,1,1,1,1,
			1,1};
  Int_t event_type=0;
  
  // 1: neutrino candidate
  // 2: light mismatch
  // 3: through-going muon
  // 4: stopped muon
  
  T_match->Branch("tpc_cluster_id",&ncluster,"tpc_cluster_id/I");
  T_match->Branch("flash_id",&flash_id,"flash_id/I");
  T_match->Branch("strength",&strength,"strength/D");
  T_match->Branch("pe_pred",pe_pred,"pe_pred[32]/D");
  T_match->Branch("pe_meas",pe_meas,"pe_meas[32]/D");
  T_match->Branch("pe_meas_err",pe_meas_err,"pe_meas_err[32]/D");
  T_match->Branch("event_type",&event_type,"event_type/I");
  
   bool flag_close_to_PMT=false;
   bool flag_at_x_boundary=false;
   double ks_dis=1;
   double chi2=1;
   int ndf=1;
   double cluster_length=1;

   T_match->Branch("flag_close_to_PMT",&flag_close_to_PMT,"flag_close_to_PMT/B");
   T_match->Branch("flag_at_x_boundary",&flag_at_x_boundary,"flag_at_x_boundary/B");
   T_match->Branch("ks_dis",&ks_dis,"ks_dis/D");
   T_match->Branch("chi2",&chi2,"chi2/D");
   T_match->Branch("ndf",&ndf,"ndf/I");
   T_match->Branch("cluster_length",&cluster_length,"cluster_length/D");

   T_match->Fill();
   
   
   if (T_bad_ch!=0){
     T_bad_ch->CloneTree(-1,"fast");
   }

   if (bee_debug!=0){
     TTree *t_bad = new TTree("T_bad","T_bad");
     t_bad->SetDirectory(file1);
     Int_t bad_npoints;
     
     Double_t bad_y[100],bad_z[100];
     t_bad->Branch("cluster_id",&ncluster,"cluster_id/I");
     t_bad->Branch("bad_npoints",&bad_npoints,"bad_npoints/I");
     t_bad->Branch("bad_y",bad_y,"bad_y[bad_npoints]/D");
     t_bad->Branch("bad_z",bad_z,"bad_z[bad_npoints]/D");
     
     for (size_t j = 0; j!= dead_clusters.size(); j++){
       SMGCSelection& mcells = dead_clusters.at(j)->get_mcells();
       ncluster = dead_clusters.at(j)->get_cluster_id();
       for (size_t i=0;i!=mcells.size();i++){
  	 PointVector ps = mcells.at(i)->boundary();
  	 bad_npoints = ps.size();
  	 for (int k=0;k!=bad_npoints;k++){
  	   bad_y[k] = ps.at(k).y/units::cm;
  	   bad_z[k] = ps.at(k).z/units::cm;
  	 }
  	 t_bad->Fill();
       }
     }
     
     TTree *T_cluster ;
     Double_t x,y,z,q,nq;
     Int_t temp_time_slice, ch_u, ch_v, ch_w;
     T_cluster = new TTree("T_cluster","T_cluster");
     T_cluster->Branch("cluster_id",&ncluster,"cluster_id/I");
     T_cluster->Branch("x",&x,"x/D");
     T_cluster->Branch("y",&y,"y/D");
     T_cluster->Branch("z",&z,"z/D");
     T_cluster->Branch("q",&q,"q/D");
     T_cluster->Branch("nq",&nq,"nq/D");
     T_cluster->Branch("time_slice",&temp_time_slice,"time_slice/I");
     T_cluster->Branch("ch_u",&ch_u,"ch_u/I");
     T_cluster->Branch("ch_v",&ch_v,"ch_v/I");
     T_cluster->Branch("ch_w",&ch_w,"ch_w/I");
     T_cluster->SetDirectory(file1);
     
     Double_t pu, pv, pw, pt;
     Double_t charge_save=1, ncharge_save=1, chi2_save=1, ndf_save=1;
     TTree *T_rec = new TTree("T_rec","T_rec");
     T_rec->Branch("x",&x,"x/D");
     T_rec->Branch("y",&y,"y/D");
     T_rec->Branch("z",&z,"z/D");
     T_rec->Branch("q",&charge_save,"q/D");
     T_rec->Branch("nq",&ncharge_save,"nq/D");
     T_rec->Branch("chi2",&chi2_save,"chi2/D");
     T_rec->Branch("ndf",&ndf_save,"ndf/D");
     T_rec->Branch("pu",&pu,"pu/D");
     T_rec->Branch("pv",&pv,"pv/D");
     T_rec->Branch("pw",&pw,"pw/D");
     T_rec->Branch("pt",&pt,"pt/D");
     T_rec->SetDirectory(file1);
     
     
     
     TTree *t_rec_charge = new TTree("T_rec_charge","T_rec_charge");
     t_rec_charge->SetDirectory(file1);
     t_rec_charge->Branch("x",&x,"x/D");
     t_rec_charge->Branch("y",&y,"y/D");
     t_rec_charge->Branch("z",&z,"z/D");
     t_rec_charge->Branch("q",&charge_save,"q/D");
     t_rec_charge->Branch("nq",&ncharge_save,"nq/D");
     t_rec_charge->Branch("chi2",&chi2_save,"chi2/D");
     t_rec_charge->Branch("ndf",&ndf_save,"ndf/D");
     t_rec_charge->Branch("pu",&pu,"pu/D");
     t_rec_charge->Branch("pv",&pv,"pv/D");
     t_rec_charge->Branch("pw",&pw,"pw/D");
     t_rec_charge->Branch("pt",&pt,"pt/D");

   

     ncluster = live_clusters.at(0)->get_cluster_id();
     
       
 
     double offset_x =0;
         
     PR3DClusterSelection temp_clusters = live_clusters;
     for (size_t j = 0; j!= temp_clusters.size(); j++){
       // show individual clusters ... 
       //	 ncluster = temp_clusters.at(j)->get_cluster_id();
	 
       SMGCSelection& mcells = temp_clusters.at(j)->get_mcells();
       for (size_t i=0;i!=mcells.size();i++){
	 PointVector ps = mcells.at(i)->get_sampling_points();
	 int time_slice = mcells.at(i)->GetTimeSlice();
	 
	 if (ps.size()==0){
	   std::cout << "zero sampling points!" << std::endl;
	 }else{
	   q = mcells.at(i)->get_q() / ps.size();
	   nq = ps.size();
	   for (int k=0;k!=ps.size();k++){
	     x = (ps.at(k).x- offset_x)/units::cm ;
	     y = ps.at(k).y/units::cm;
	     z = ps.at(k).z/units::cm;
	     std::vector<int> time_chs = ct_point_cloud.convert_3Dpoint_time_ch(ps.at(k));
	     temp_time_slice = time_chs.at(0);
	     ch_u = time_chs.at(1);
	     ch_v = time_chs.at(2);
	     ch_w = time_chs.at(3);
	     
	     T_cluster->Fill();
	   }
	 }
       }
     }

     
     for (size_t j = 0; j!= live_clusters.size(); j++){
       
       // save wcps
       PR3DCluster *new_cluster = live_clusters.at(j);//old_new_cluster_map[live_clusters.at(j)];
       std::list<WCPointCloud<double>::WCPoint>& wcps_list = new_cluster->get_path_wcps();
       ndf_save = live_clusters.at(j)->get_cluster_id();
       charge_save = 0;
       ncharge_save = 0;
       chi2_save = 0;
       for (auto it = wcps_list.begin(); it!=wcps_list.end(); it++){
  	 x = (*it).x/units::cm;
  	 y = (*it).y/units::cm;
  	 z = (*it).z/units::cm;
	 
  	 Point temp_p((*it).x,(*it).y,(*it).z);
  	 std::vector<int> time_chs = ct_point_cloud.convert_3Dpoint_time_ch(temp_p);
  	 pt = time_chs.at(0);
  	 pu = time_chs.at(1);
  	 pv = time_chs.at(2);
  	 pw = time_chs.at(3);
	 
  	 T_rec->Fill();
       }
       
     }
   }

   
   Trun->CloneTree(-1,"fast");
   
   TTree *T_flash = new TTree("T_flash","T_flash");
   T_flash->SetDirectory(file1);
   int type=2;
   double low_time=3, high_time=5;
   double time=4;
   double total_PE=100;
   double PE[32]={1,1,1,1,1,1,1,1,1,1,
		  1,1,1,1,1,1,1,1,1,1,
		  1,1,1,1,1,1,1,1,1,1,
		  1,1};
   double PE_err[32]={1,1,1,1,1,1,1,1,1,1,
		  1,1,1,1,1,1,1,1,1,1,
		  1,1,1,1,1,1,1,1,1,1,
		  1,1};
   std::vector<int> fired_channels;
   std::vector<double> l1_fired_time;
   std::vector<double> l1_fired_pe;
   flash_id = 0;
   T_flash->Branch("type",&type);
   T_flash->Branch("flash_id",&flash_id);
   T_flash->Branch("low_time",&low_time);
   T_flash->Branch("high_time",&high_time);
   T_flash->Branch("time",&time);
   T_flash->Branch("total_PE",&total_PE);
   T_flash->Branch("PE",PE,"PE[32]/D");
   T_flash->Branch("PE_err",PE_err,"PE_err[32]/D");
   T_flash->Branch("fired_channels",&fired_channels);
   T_flash->Branch("l1_fired_time",&l1_fired_time);
   T_flash->Branch("l1_fired_pe",&l1_fired_pe);
   
   T_flash->Fill();
    
   TDC->CloneTree(-1,"fast");
   {
    TTree *TC_n = new TTree("TC","TC");
    TC_n->SetDirectory(file1);
    std::vector<int> *cluster_id = new std::vector<int>;
    //
    std::vector<int> *parent_cluster_id = new std::vector<int>;
    //
    std::vector<int> *time_slice = new std::vector<int>;
    std::vector<double> *q = new std::vector<double>;    
    std::vector<double> *uq = new std::vector<double>;
    std::vector<double> *vq = new std::vector<double>;
    std::vector<double> *wq = new std::vector<double>;    
    std::vector<double> *udq = new std::vector<double>;
    std::vector<double> *vdq = new std::vector<double>;   
    std::vector<double> *wdq = new std::vector<double>; 
    
    TC_n->Branch("cluster_id",&cluster_id);
    TC_n->Branch("parent_cluster_id",&parent_cluster_id);
    TC_n->Branch("time_slice",&time_slice);
    TC_n->Branch("q",&q);
    TC_n->Branch("uq",&uq);
    TC_n->Branch("vq",&vq);
    TC_n->Branch("wq",&wq);
    TC_n->Branch("udq",&udq);
    TC_n->Branch("vdq",&vdq);
    TC_n->Branch("wdq",&wdq);

    std::vector<int> *nwire_u = new  std::vector<int>;
    std::vector<int> *nwire_v = new  std::vector<int>;
    std::vector<int> *nwire_w = new  std::vector<int>;
    std::vector<int> *flag_u = new  std::vector<int>;
    std::vector<int> *flag_v = new  std::vector<int>;
    std::vector<int> *flag_w = new  std::vector<int>;
    

    std::vector<std::vector<int>> *wire_index_u = new std::vector<std::vector<int>>;
    std::vector<std::vector<int>> *wire_index_v = new std::vector<std::vector<int>>;
    std::vector<std::vector<int>> *wire_index_w = new std::vector<std::vector<int>>;
    std::vector<std::vector<double>> *wire_charge_u = new std::vector<std::vector<double>>;
    std::vector<std::vector<double>> *wire_charge_v = new std::vector<std::vector<double>>;
    std::vector<std::vector<double>> *wire_charge_w = new std::vector<std::vector<double>>;
    std::vector<std::vector<double>> *wire_charge_err_u = new std::vector<std::vector<double>>;
    std::vector<std::vector<double>> *wire_charge_err_v = new std::vector<std::vector<double>>;
    std::vector<std::vector<double>> *wire_charge_err_w = new std::vector<std::vector<double>>;

    

    TC_n->Branch("nwire_u",&nwire_u);
    TC_n->Branch("nwire_v",&nwire_v);
    TC_n->Branch("nwire_w",&nwire_w);
    TC_n->Branch("flag_u",&flag_u);
    TC_n->Branch("flag_v",&flag_v);
    TC_n->Branch("flag_w",&flag_w);
    TC_n->Branch("wire_index_u",&wire_index_u);
    TC_n->Branch("wire_index_v",&wire_index_v);
    TC_n->Branch("wire_index_w",&wire_index_w);
    TC_n->Branch("wire_charge_u",&wire_charge_u);
    TC_n->Branch("wire_charge_v",&wire_charge_v);
    TC_n->Branch("wire_charge_w",&wire_charge_w);
    TC_n->Branch("wire_charge_err_u",&wire_charge_err_u);
    TC_n->Branch("wire_charge_err_v",&wire_charge_err_v);
    TC_n->Branch("wire_charge_err_w",&wire_charge_err_w);

    

    int temp_parent_cluster_id = live_clusters.at(0)->get_cluster_id();
    PR3DClusterSelection temp_clusters = live_clusters;
    for (auto it1 = temp_clusters.begin(); it1!=temp_clusters.end(); it1++){
      SMGCSelection& mcells = (*it1)->get_mcells();
      
      for (int i=0; i!=mcells.size();i++){
	cluster_id->push_back((*it1)->get_cluster_id());
	parent_cluster_id->push_back(temp_parent_cluster_id);
	
	SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)mcells.at(i);
	int time_slice_temp = mcell->GetTimeSlice();
	time_slice->push_back(time_slice_temp);
	q->push_back(mcell->get_q());
	uq->push_back(mcell->get_uq());
	vq->push_back(mcell->get_vq());
	wq->push_back(mcell->get_wq());
	udq->push_back(mcell->get_udq());
	vdq->push_back(mcell->get_vdq());
	wdq->push_back(mcell->get_wdq());
	
	std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();
	flag_u->push_back(1);
	flag_v->push_back(1);
	flag_w->push_back(1);
	for (size_t i=0;i!=bad_planes.size();i++){
	  if (bad_planes.at(i)==(WirePlaneType_t)(0)){
	    flag_u->back() = 0;
	  }else if (bad_planes.at(i)==(WirePlaneType_t)(1)){
	    flag_v->back() = 0;
	  }else if (bad_planes.at(i)==(WirePlaneType_t)(2)){
	    flag_w->back() = 0;
	  }else{
	  }
	}
	GeomWireSelection& uwires = mcell->get_uwires();
	GeomWireSelection& vwires = mcell->get_vwires();
	GeomWireSelection& wwires = mcell->get_wwires();
	nwire_u->push_back(uwires.size());
	nwire_v->push_back(vwires.size());
	nwire_w->push_back(wwires.size());
	
	std::vector<int> wire_index_ui;
	std::vector<int> wire_index_vi;
	std::vector<int> wire_index_wi;
	std::vector<double> wire_charge_ui;
	std::vector<double> wire_charge_vi;
	std::vector<double> wire_charge_wi;
	std::vector<double> wire_charge_err_ui;
	std::vector<double> wire_charge_err_vi;
	std::vector<double> wire_charge_err_wi;
	
	for (auto it2 = uwires.begin(); it2!=uwires.end(); it2++){
	  const GeomWire *wire = (*it2);
	  wire_index_ui.push_back(wire->index());
	  wire_charge_ui.push_back(mcell->Get_Wire_Charge(wire));
	  wire_charge_err_ui.push_back(mcell->Get_Wire_Charge_Err(wire));
	}
	for (auto it2 = vwires.begin(); it2!=vwires.end(); it2++){
	  const GeomWire *wire = (*it2);
	  wire_index_vi.push_back(wire->index());
	  wire_charge_vi.push_back(mcell->Get_Wire_Charge(wire));
	  wire_charge_err_vi.push_back(mcell->Get_Wire_Charge_Err(wire));
	}
	for (auto it2 = wwires.begin(); it2!=wwires.end(); it2++){
	  const GeomWire *wire = (*it2);
	  wire_index_wi.push_back(wire->index());
	  wire_charge_wi.push_back(mcell->Get_Wire_Charge(wire));
	  wire_charge_err_wi.push_back(mcell->Get_Wire_Charge_Err(wire));
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
      }
    }
    
    
    TC_n->Fill();
    cluster_id->clear();
    parent_cluster_id->clear();
    nwire_u->clear();
    nwire_v->clear();
    nwire_w->clear();
    flag_u->clear();
    flag_v->clear();
    flag_w->clear();
    wire_index_u->clear();
    wire_index_v->clear();
    wire_index_w->clear();
  }

     
   file1->Write();
   file1->Close();
}
