#include "WireCellSst/GeomDataSource.h"
#include "WireCellData/PR3DCluster.h"
#include "WireCellData/SlimMergeGeomCell.h"
#include "WireCellData/TPCParams.h"
#include "WireCellData/Singleton.h"
#include "WireCellData/ToyCTPointCloud.h"

#include "WireCell2dToy/ExecMon.h"
#include "WireCell2dToy/CalcPoints.h"
#include "WireCell2dToy/ToyClustering.h"
#include "WireCell2dToy/ToyLightReco.h"


#include "WireCell2dToy/ToyMatching.h"
#include "WireCell2dToy/ToyFiducial.h"

#include "WireCell2dToy/ImprovePR3DCluster.h"


#include "TH1F.h"
#include "TH2F.h"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"

using namespace WireCell;
using namespace std;

int main(int argc, char* argv[])
{
  if (argc < 3) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/imaging.root" << endl;
    return 1;
  }
  TH1::AddDirectory(kFALSE);
  
  int flag_pos_corr = 0; // correct X position after matching ... 
  for(Int_t i = 1; i != argc; i++){
     switch(argv[i][1]){
     case 'c':
       flag_pos_corr = atoi(&argv[i][2]); 
       break;
     }
  }
  
  
  ExecMon em("starting");
  cout << em("load geometry") << endl;

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
  
  
  std::vector<float> *op_gain = new std::vector<float>;
  std::vector<float> *op_gainerror = new std::vector<float>;
  double triggerTime;
  
  Trun->SetBranchAddress("op_gain",&op_gain); 
  Trun->SetBranchAddress("op_gainerror",&op_gainerror); 
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


  WireCell2dToy::ToyFiducial *fid = new WireCell2dToy::ToyFiducial(3,800,-first_u_dis/pitch_u, -first_v_dis/pitch_v, -first_w_dis/pitch_w,
								   1./time_slice_width, 1./pitch_u, 1./pitch_v, 1./pitch_w, // slope
								   angle_u,angle_v,angle_w,// angle
								   3*units::cm);
  
  
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
  }
  
  if (timesliceId->size()==0){
    std::cout << "No points! Quit! " << std::endl;
    return 0;
  }
  
  cout << em("load clusters from file") << endl;
    
  
  // Start to add X, Y, Z points
  // form boundaries for the bad cells ... 
  for (size_t j = 0; j!= dead_clusters.size(); j++){
    WireCell2dToy::calc_boundary_points_dead(gds,dead_clusters.at(j));
  }
  // form sampling points for the normal cells ...
  DynamicToyPointCloud global_point_cloud(angle_u,angle_v,angle_w);
  for (size_t i=0; i!=live_clusters.size();i++){
    WireCell2dToy::calc_sampling_points(gds,live_clusters.at(i),nrebin, frame_length, unit_dis);
    live_clusters.at(i)->Create_point_cloud();
    global_point_cloud.AddPoints(live_clusters.at(i),0);
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
  
  std::map<PR3DCluster*,std::vector<std::pair<PR3DCluster*,double>>> group_clusters = WireCell2dToy::Clustering_jump_gap_cosmics(live_clusters, dead_clusters,dead_u_index, dead_v_index, dead_w_index, global_point_cloud, ct_point_cloud);
  cout << em("Clustering to jump gap in cosmics") << std::endl;
    
  // process light information
  WireCell2dToy::ToyLightReco uboone_flash(filename, 1); // 1: imagingoutput, "Trun"; default/not specfified: path "Event/Sim"
  uboone_flash.load_event_raw(0);
  cout << em("flash reconstruction") << std::endl;
  
  // prepare light matching ....
  WireCell::OpflashSelection& flashes = uboone_flash.get_flashes();
  for (size_t i=0;i!=flashes.size(); i++){
    flashes.at(i)->set_flash_id(i);
  }
  
  // // form a global map with the current map information
  // std::map<int,std::map<const GeomWire*, SMGCSelection > > global_wc_map;
  // for (size_t i=0; i!=live_clusters.size();i++){
  //   PR3DCluster *cluster = live_clusters.at(i);
  //   SMGCSelection& mcells = cluster->get_mcells();
  //   for (auto it = mcells.begin(); it!= mcells.end(); it++){
  //     SlimMergeGeomCell *mcell = *it;
  //     int time_slice = mcell->GetTimeSlice();
  //     if (global_wc_map.find(time_slice)==global_wc_map.end()){
  // 	std::map<const GeomWire*, SMGCSelection> temp_wc_map;
  // 	global_wc_map[time_slice] = temp_wc_map;
  //     }
  //     std::map<const GeomWire*, SMGCSelection>& timeslice_wc_map = global_wc_map[time_slice];
      
  //     GeomWireSelection& uwires = mcell->get_uwires();
  //     GeomWireSelection& vwires = mcell->get_vwires();
  //     GeomWireSelection& wwires = mcell->get_wwires();
  //     std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();
  //     if (find(bad_planes.begin(),bad_planes.end(),WirePlaneType_t(0))==bad_planes.end()){
  // 	for (int j=0;j!=uwires.size();j++){
  // 	  const GeomWire *wire = uwires.at(j);
  // 	  if (timeslice_wc_map.find(wire)==timeslice_wc_map.end()){
  // 	    SMGCSelection temp_mcells;
  // 	    temp_mcells.push_back(mcell);
  // 	    timeslice_wc_map[wire] = temp_mcells;
  // 	  }else{
  // 	    timeslice_wc_map[wire].push_back(mcell);
  // 	  }
  // 	}
  //     }
  //     if (find(bad_planes.begin(),bad_planes.end(),WirePlaneType_t(1))==bad_planes.end()){
  // 	for (int j=0;j!=vwires.size();j++){
  // 	  const GeomWire *wire = vwires.at(j);
  // 	  if (timeslice_wc_map.find(wire)==timeslice_wc_map.end()){
  // 	    SMGCSelection temp_mcells;
  // 	    temp_mcells.push_back(mcell);
  // 	    timeslice_wc_map[wire] = temp_mcells;
  // 	  }else{
  // 	    timeslice_wc_map[wire].push_back(mcell);
  // 	  }
  // 	}
  //     }
  //     if (find(bad_planes.begin(),bad_planes.end(),WirePlaneType_t(2))==bad_planes.end()){
  // 	for (int j=0;j!=wwires.size();j++){
  // 	  const GeomWire *wire = wwires.at(j);
  // 	  if (timeslice_wc_map.find(wire)==timeslice_wc_map.end()){
  // 	    SMGCSelection temp_mcells;
  // 	    temp_mcells.push_back(mcell);
  // 	    timeslice_wc_map[wire] = temp_mcells;
  // 	  }else{
  // 	    timeslice_wc_map[wire].push_back(mcell);
  // 	  }
  // 	}
  //     }
  //   }
  // }
  // //
   
  FlashTPCBundleSelection matched_bundles = WireCell2dToy::tpc_light_match(time_offset,nrebin,group_clusters,flashes);
  cout << em("TPC Light Matching") << std::endl;
  
  // create the live clusters ...
  live_clusters.clear();
  
  for (auto it = matched_bundles.begin(); it!= matched_bundles.end(); it++){
    FlashTPCBundle *bundle = *it;
    PR3DCluster *main_cluster = bundle->get_main_cluster();
    live_clusters.push_back(main_cluster);
    for (auto it1 = bundle->get_other_clusters().begin(); it1!=bundle->get_other_clusters().end();it1++){
      live_clusters.push_back(*it1);
    }
  }
  
  std::map<PR3DCluster*, PR3DCluster*> old_new_cluster_map;
  for (size_t i=0;i!=live_clusters.size();i++){
    live_clusters.at(i)->Create_graph(ct_point_cloud);
    
    std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = live_clusters.at(i)->get_highest_lowest_wcps();
    
    live_clusters.at(i)->dijkstra_shortest_paths(wcps.first);
    live_clusters.at(i)->cal_shortest_path(wcps.second);
    {
      // add dead channels in??? 
      PR3DCluster *new_cluster = WireCell2dToy::Improve_PR3DCluster(live_clusters.at(i),ct_point_cloud, gds);
      WireCell2dToy::calc_sampling_points(gds,new_cluster,nrebin, frame_length, unit_dis);
      new_cluster->Create_point_cloud();
      old_new_cluster_map[live_clusters.at(i)] = new_cluster;
      
      //      new_cluster->Create_graph(ct_point_cloud);
      //std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> new_wcps = new_cluster->get_highest_lowest_wcps();
      //new_cluster->dijkstra_shortest_paths(new_wcps.first);
      //new_cluster->cal_shortest_path(new_wcps.second);
    }
    
    //    live_clusters.at(i)->fine_tracking(global_wc_map);
    //    live_clusters.at(i)->collect_charge_trajectory(ct_point_cloud);
  }


  // do the dQ/dx fitting ... 
  for (auto it = matched_bundles.begin(); it!= matched_bundles.end(); it++){
    FlashTPCBundle *bundle = *it;
    Opflash *flash = bundle->get_flash();
    PR3DCluster *main_cluster = bundle->get_main_cluster();
    
    if (flash!=0){
      
      // if (flash->get_time() > 2 && flash->get_time() < 6){
      // std::cout << flash->get_time() << std::endl;
      // main_cluster->dQ_dx_fit(global_wc_map, flash->get_time()*units::microsecond);
      // }
    }
    
  }
  
  cerr << em("Create Graph in all clusters") << std::endl;
  
  
  TFile *file1 = new TFile(Form("matching_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");
  TTree *T_match = new TTree("T_match","T_match");
  T_match->SetDirectory(file1);
  Int_t ncluster=0;
  Int_t flash_id;
  Double_t strength;
  Double_t pe_pred[32];
  Double_t pe_meas[32];
  Double_t pe_meas_err[32];
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
  
  bool flag_close_to_PMT;
  bool flag_at_x_boundary;
  double ks_dis;
  double chi2;
  int ndf;
  double cluster_length;
  
  T_match->Branch("flag_close_to_PMT",&flag_close_to_PMT,"flag_close_to_PMT/B");
  T_match->Branch("flag_at_x_boundary",&flag_at_x_boundary,"flag_at_x_boundary/B");
  T_match->Branch("ks_dis",&ks_dis,"ks_dis/D");
  T_match->Branch("chi2",&chi2,"chi2/D");
  T_match->Branch("ndf",&ndf,"ndf/I");
  T_match->Branch("cluster_length",&cluster_length,"cluster_length/D");
  
  for (auto it = matched_bundles.begin(); it!=matched_bundles.end(); it++){
    FlashTPCBundle *bundle = *it;
    
    Opflash *flash = bundle->get_flash();
    PR3DCluster *main_cluster = bundle->get_main_cluster();
    cluster_length = -1;
    
    if (flash!=0){
      auto it1 = find(flashes.begin(),flashes.end(),flash);
      flash_id = flash->get_flash_id();
      strength = bundle->get_strength();
      std::vector<double> temp = bundle->get_pred_pmt_light();
      for (int i=0;i!=32;i++){
	pe_pred[i] = temp.at(i);
	pe_meas[i] = flash->get_PE(i);
	pe_meas_err[i] = flash->get_PE_err(i);
      }
      flag_close_to_PMT = bundle->get_flag_close_to_PMT();
      flag_at_x_boundary = bundle->get_flag_at_x_boundary();
      ks_dis = bundle->get_ks_dis();
      chi2 = bundle->get_chi2();
      ndf = bundle->get_ndf();
      
      
    }else{
      flash_id = -1;
      strength = 0;
      for (int i=0;i!=32;i++){
	pe_pred[i] = 0;
	pe_meas[i] = 0;
	pe_meas_err[i] = 0.;
      }
      flag_close_to_PMT = false;
      flag_at_x_boundary = false;
      ks_dis = 1;
      chi2 = 1e9;
      ndf = 32;
    }
    ncluster = main_cluster->get_cluster_id();
    
    // check if this is through going muon ...
    event_type = 0;
    if (flash!=0){
      //std::cout << "Flash: " << flash->get_flash_id() << " " << flash->get_time() << std::endl;
      
      // Note, for BNB, flash time is aroudn 4us
      // time offset happened to be 4us ...
      double offset_x = (flash->get_time() - time_offset)*2./nrebin*time_slice_width;
      if (fid->check_tgm(bundle,offset_x, ct_point_cloud,old_new_cluster_map))
	event_type |= 1UL << 3; // 3rd bit for TGM
      
      int temp_flag = fid->check_LM(bundle,cluster_length);
      if (temp_flag==1){
	event_type |= 1UL << 4; // 4th bit for low energy ...
      }else if (temp_flag==2){
	event_type |= 1UL << 1; // 1st bit for light mismatch ...
      }
    }
    
    T_match->Fill();
  }
  
  cerr << em("test TGM") << std::endl;
  
  // clone TDC tree ... 
  TDC->CloneTree()->Write();
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

    
    for (auto it = matched_bundles.begin(); it!= matched_bundles.end(); it++){
      FlashTPCBundle *bundle = *it;
      PR3DCluster *main_cluster = bundle->get_main_cluster();
      Opflash *flash = bundle->get_flash();
      // now prepare saving it
      int temp_parent_cluster_id = main_cluster->get_cluster_id();
      PR3DClusterSelection temp_clusters;
      temp_clusters.push_back(main_cluster);
      for (auto it1 = bundle->get_other_clusters().begin(); it1!=bundle->get_other_clusters().end();it1++){
	temp_clusters.push_back(*it1);
      }

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
  //

  
  TTree *T_bad_ch = (TTree*)file->Get("T_bad_ch");
  if (T_bad_ch!=0){
    T_bad_ch->CloneTree()->Write();
  }
  
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
  // TTree *T_rec = new TTree("T_rec","T_rec");
  // T_rec->Branch("x",&x,"x/D");
  // T_rec->Branch("y",&y,"y/D");
  // T_rec->Branch("z",&z,"z/D");
  // T_rec->Branch("q",&charge_save,"q/D");
  // T_rec->Branch("nq",&nchareg_save,"nq/D");
  // T_rec->Branch("chi2",&chi2_save,"chi2/D");
  // T_rec->Branch("ndf",&ndf_save,"ndf/D");
  // T_rec->Branch("pu",&pu,"pu/D");
  // T_rec->Branch("pv",&pv,"pv/D");
  // T_rec->Branch("pw",&pw,"pw/D");
  // T_rec->Branch("pt",&pt,"pt/D");
  // T_rec->SetDirectory(file1);
  
  // TTree *t_rec_charge = new TTree("T_rec_charge","T_rec_charge");
  // t_rec_charge->SetDirectory(file1);
  // t_rec_charge->Branch("x",&x,"x/D");
  // t_rec_charge->Branch("y",&y,"y/D");
  // t_rec_charge->Branch("z",&z,"z/D");
  // t_rec_charge->Branch("q",&charge_save,"q/D");
  // t_rec_charge->Branch("nq",&ncharge_save,"nq/D");
  // t_rec_charge->Branch("chi2",&chi2_save,"chi2/D");
  // t_rec_charge->Branch("ndf",&ndf_save,"ndf/D");
  // t_rec_charge->Branch("pu",&pu,"pu/D");
  // t_rec_charge->Branch("pv",&pv,"pv/D");
  // t_rec_charge->Branch("pw",&pw,"pw/D");
  // t_rec_charge->Branch("pt",&pt,"pt/D");
  
  // TTree *T_proj_data = new TTree("T_proj_data","T_proj_data");
  // std::vector<int> *proj_data_cluster_id = new std::vector<int>;
  // std::vector<std::vector<int>> *proj_data_cluster_channel = new std::vector<std::vector<int>>;
  // std::vector<std::vector<int>> *proj_data_cluster_timeslice= new std::vector<std::vector<int>>;
  // std::vector<std::vector<int>> *proj_data_cluster_charge= new std::vector<std::vector<int>>;
  // std::vector<std::vector<int>> *proj_data_cluster_charge_err= new std::vector<std::vector<int>>;
  // std::vector<std::vector<int>> *proj_data_cluster_charge_pred= new std::vector<std::vector<int>>;
  
  // T_proj_data->Branch("cluster_id",&proj_data_cluster_id);
  // T_proj_data->Branch("channel",&proj_data_cluster_channel);
  // T_proj_data->Branch("time_slice",&proj_data_cluster_timeslice);
  // T_proj_data->Branch("charge",&proj_data_cluster_charge);
  // T_proj_data->Branch("charge_err",&proj_data_cluster_charge_err);
  // T_proj_data->Branch("charge_pred",&proj_data_cluster_charge_pred);
  // T_proj_data->SetDirectory(file1);
  
  // note did not save the unmatched cluster ... 
  ncluster = 0;
  //   for (auto it = matched_results.begin(); it!= matched_results.end(); it++){
  for (auto it = matched_bundles.begin(); it!= matched_bundles.end(); it++){
    FlashTPCBundle *bundle = *it;
    PR3DCluster *main_cluster = bundle->get_main_cluster();//std::get<0>(*it);
    Opflash *flash = bundle->get_flash();//std::get<1>(*it);
    double offset_x ;
    if (flash!=0){
      offset_x = (flash->get_time() - time_offset)*2./nrebin*time_slice_width;
    }else{
      offset_x = 0;
    }
    
    if (flag_pos_corr==0)
      offset_x = 0;
    
    ncluster = main_cluster->get_cluster_id();
    
    PR3DClusterSelection temp_clusters;
    temp_clusters.push_back(main_cluster);
    for (auto it1 = bundle->get_other_clusters().begin(); it1!=bundle->get_other_clusters().end();it1++){
      temp_clusters.push_back(*it1);
    }
    for (size_t j = 0; j!= temp_clusters.size(); j++){
      
      // show individual clusters ... 
      SMGCSelection& mcells = temp_clusters.at(j)->get_mcells();
      //ncluster = temp_clusters.at(0)->get_cluster_id();
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
    //     ncluster ++;
  }
  
  for (size_t j = 0; j!= live_clusters.size(); j++){
    
    // save wcps
    PR3DCluster *new_cluster = old_new_cluster_map[live_clusters.at(j)];
    // std::list<WCPointCloud<double>::WCPoint>& wcps_list = new_cluster->get_path_wcps();
  //   //ncluster = -1 * ncluster-100;
  //   ndf_save = live_clusters.at(j)->get_cluster_id();
  //   charge_save = 0;
  //   ncharge_save = 0;
  //   chi2_save = 0;
  //   for (auto it = wcps_list.begin(); it!=wcps_list.end(); it++){
  //     x = (*it).x/units::cm;
  //     y = (*it).y/units::cm;
  //     z = (*it).z/units::cm;
      
  //     Point temp_p((*it).x,(*it).y,(*it).z);
  //     std::vector<int> time_chs = ct_point_cloud.convert_3Dpoint_time_ch(temp_p);
  //     pt = time_chs.at(0);
  //     pu = time_chs.at(1);
  //     pv = time_chs.at(2);
  //     pw = time_chs.at(3);
      
  //     T_rec->Fill();
  //   }
    
    PointVector& pts = live_clusters.at(j)->get_fine_tracking_path();
    std::vector<double>& dQ = live_clusters.at(j)->get_dQ();
    std::vector<double>& dx = live_clusters.at(j)->get_dx();
    std::vector<double>& tpu = live_clusters.at(j)->get_pu();
    std::vector<double>& tpv = live_clusters.at(j)->get_pv();
    std::vector<double>& tpw = live_clusters.at(j)->get_pw();
    std::vector<double>& tpt = live_clusters.at(j)->get_pt();
    
    std::map<std::pair<int,int>, std::tuple<double,double,double> > & proj_data_u_map = live_clusters.at(j)->get_proj_data_u_map();
    std::map<std::pair<int,int>, std::tuple<double,double,double> > & proj_data_v_map = live_clusters.at(j)->get_proj_data_v_map();
    std::map<std::pair<int,int>, std::tuple<double,double,double> > & proj_data_w_map = live_clusters.at(j)->get_proj_data_w_map();
    
    // ndf_save = live_clusters.at(j)->get_cluster_id();
    // if (pts.size()==dQ.size()){
    //   proj_data_cluster_id->push_back(ndf_save);
    //   std::vector<int> temp_channel;
    //   std::vector<int> temp_timeslice;
    //   std::vector<int> temp_charge;
    //   std::vector<int> temp_charge_err;
    //   std::vector<int> temp_charge_pred;
    //   for (auto it = proj_data_u_map.begin(); it!=proj_data_u_map.end(); it++){
    // 	temp_channel.push_back(it->first.first);
    // 	temp_timeslice.push_back(it->first.second);
    // 	temp_charge.push_back(std::get<0>(it->second));
    // 	temp_charge_err.push_back(std::get<1>(it->second));
    // 	temp_charge_pred.push_back(std::get<2>(it->second));
    //   }
    //   for (auto it = proj_data_v_map.begin(); it!=proj_data_v_map.end(); it++){
    // 	temp_channel.push_back(it->first.first);
    // 	temp_timeslice.push_back(it->first.second);
    // 	temp_charge.push_back(std::get<0>(it->second));
    // 	temp_charge_err.push_back(std::get<1>(it->second));
    // 	temp_charge_pred.push_back(std::get<2>(it->second));
    //   }
    //   for (auto it = proj_data_w_map.begin(); it!=proj_data_w_map.end(); it++){
    // 	temp_channel.push_back(it->first.first);
    // 	temp_timeslice.push_back(it->first.second);
    // 	temp_charge.push_back(std::get<0>(it->second));
    // 	temp_charge_err.push_back(std::get<1>(it->second));
    // 	temp_charge_pred.push_back(std::get<2>(it->second));
    //   }
    //   proj_data_cluster_channel->push_back(temp_channel);
    //   proj_data_cluster_timeslice->push_back(temp_timeslice);
    //   proj_data_cluster_charge->push_back(temp_charge);
    //   proj_data_cluster_charge_err->push_back(temp_charge_err);
    //   proj_data_cluster_charge_pred->push_back(temp_charge_pred);
    // }
    
    
    // for (size_t i=0; i!=pts.size(); i++){
    //   x = pts.at(i).x/units::cm;
    //   y = pts.at(i).y/units::cm;
    //   z = pts.at(i).z/units::cm;
      
    //   if (pts.size()==dQ.size()){
    // 	charge_save = dQ.at(i);
    // 	ncharge_save = dx.at(i)/units::cm;
    // 	pu = tpu.at(i);
    // 	pv = tpv.at(i);
    // 	pw = tpw.at(i);
    // 	pt = tpt.at(i);
	
	
    // 	t_rec_charge->Fill();
	
    //   }else{
    // 	charge_save = 0;
    // 	ncharge_save = -1;
    // 	pu = -1;
    // 	pv = -1;
    // 	pw = -1;
    // 	pt = -1;
    //   }
    // }
  }
  
  
  Trun->CloneTree()->Write();
  
  
  // TTree *t1 = new TTree("T_data","T_data");
  // t1->SetDirectory(file1);
  
  // TClonesArray* op_wf = new TClonesArray("TH1S");
  // std::vector<short> *op_femch = new std::vector<short>;
  // std::vector<double> *op_timestamp = new std::vector<double>;
  
  // t1->Branch("op_gain",&op_gain);
  // t1->Branch("op_gainerror",&op_gainerror);
  // t1->Branch("op_femch",&op_femch);
  // t1->Branch("op_timestamp",&op_timestamp);
  // t1->Branch("op_wf",&op_wf,256000,0);
  // t1->Branch("triggerTime",&triggerTime);
  
  // t1->Branch("runNo",&run_no);
  // t1->Branch("subRunNo",&subrun_no);
  // t1->Branch("eventNo",&event_no);
  
  // op_wf = uboone_flash.get_rawWfm();
  // op_femch = uboone_flash.get_rawChan();
  // op_timestamp = uboone_flash.get_rawTimestamp();
  
  // t1->Fill();
  
  // TH2F *h1 = new TH2F("hraw","hraw",1500,0,1500,32,0,32);
  // TH2F *h2 = new TH2F("hdecon","hdecon",250,0,250,32,0,32);
  // TH2F *h3 = new TH2F("hl1","hl1",250,0,250,32,0,32);
  // h1->SetDirectory(file1);
  // h2->SetDirectory(file1);
  // h3->SetDirectory(file1);
  // for (int i=0;i!=32;i++){
  //   TH1F *h10 = uboone_flash.get_raw_hist(i);
  //   TH1F *h20 = uboone_flash.get_decon_hist(i);
  //   TH1F *h30 = uboone_flash.get_l1_hist(i);
  //   for (int j=0;j!=1500;j++){
  //     h1->SetBinContent(j+1,i+1,h10->GetBinContent(j+1));
  //   }
  //   for (int j=0;j!=250;j++){
  //     h2->SetBinContent(j+1,i+1,h20->GetBinContent(j+1));
  //     h3->SetBinContent(j+1,i+1,h30->GetBinContent(j+1));
  //   }
  // }
  
  // TH1F *h_totPE = (TH1F*)uboone_flash.get_totPE()->Clone("totPE");
  // TH1F *h_mult = (TH1F*)uboone_flash.get_mult()->Clone("mult");
  // TH1F *h_l1_mult = (TH1F*)uboone_flash.get_l1_mult()->Clone("l1_mult");
  // TH1F *h_l1_totPE = (TH1F*)uboone_flash.get_l1_totPE()->Clone("l1_totPE");
  
  // h_totPE->SetDirectory(file1);
  // h_mult->SetDirectory(file1);
  // h_l1_mult->SetDirectory(file1);
  // h_l1_totPE->SetDirectory(file1);
  
  TTree *T_flash = new TTree("T_flash","T_flash");
  T_flash->SetDirectory(file1);
  int type;
  double low_time, high_time;
  double time;
  double total_PE;
  double PE[32],PE_err[32];
  // std::vector<int> matched_tpc_ids;
  std::vector<int> fired_channels;
  std::vector<double> l1_fired_time;
  std::vector<double> l1_fired_pe;
  
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
  
  
  for (auto it = flashes.begin(); it!=flashes.end(); it++){
    fired_channels.clear();
    Opflash *flash = (*it);
    flash_id = flash->get_flash_id();
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
  
  // // now save the projected charge information ... 
  // TTree *T_proj = new TTree("T_proj","T_proj");
  // std::vector<int> *proj_cluster_id = new std::vector<int>;
  // std::vector<std::vector<int>> *proj_cluster_channel = new std::vector<std::vector<int>>;
  // std::vector<std::vector<int>> *proj_cluster_timeslice= new std::vector<std::vector<int>>;
  // std::vector<std::vector<int>> *proj_cluster_charge= new std::vector<std::vector<int>>;
  // std::vector<std::vector<int>> *proj_cluster_charge_err= new std::vector<std::vector<int>>;
  // T_proj->Branch("cluster_id",&proj_cluster_id);
  // T_proj->Branch("channel",&proj_cluster_channel);
  // T_proj->Branch("time_slice",&proj_cluster_timeslice);
  // T_proj->Branch("charge",&proj_cluster_charge);
  // T_proj->Branch("charge_err",&proj_cluster_charge_err);
  
  // T_proj->SetDirectory(file1);
  
  // for (auto it = matched_bundles.begin(); it!= matched_bundles.end(); it++){
  //   FlashTPCBundle *bundle = *it;
  //   PR3DCluster *main_cluster = bundle->get_main_cluster();
  //   Opflash *flash = bundle->get_flash();
  //   //  if (flash!=0){
    
  //   // now prepare saving it
  //   int cluster_id = main_cluster->get_cluster_id();
    
  //   PR3DClusterSelection temp_clusters;
  //   temp_clusters.push_back(main_cluster);
  //   for (auto it1 = bundle->get_other_clusters().begin(); it1!=bundle->get_other_clusters().end();it1++){
  //     temp_clusters.push_back(*it1);
  //   }
    
  //   std::vector<int> proj_channel;
  //   std::vector<int> proj_timeslice;
  //   std::vector<int> proj_charge;
  //   std::vector<int> proj_charge_err;
  //   std::vector<int> proj_flag;
    
  //   for (size_t j = 0; j!= temp_clusters.size(); j++){
  //     PR3DCluster *cluster = temp_clusters.at(j);
  //     cluster->get_projection(proj_channel,proj_timeslice,proj_charge, proj_charge_err, proj_flag, global_wc_map);
  //   }
  //   proj_cluster_id->push_back(cluster_id);
  //   proj_cluster_channel->push_back(proj_channel);
  //   proj_cluster_timeslice->push_back(proj_timeslice);
  //   proj_cluster_charge->push_back(proj_charge);
  //   proj_cluster_charge_err->push_back(proj_charge_err);
    
  // }
  // T_proj->Fill();
  //  T_proj_data->Fill();
  
  file1->Write();
  file1->Close();
}
