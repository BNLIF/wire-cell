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
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/imaging.root -d[0,1,2]" << endl;
    return 1;
  }
  TH1::AddDirectory(kFALSE);
  
  int flag_pos_corr = 0; // correct X position after matching ...
  int datatier = 0; // data=0, overlay=1, full mc=2
  int save_light = 0;
  int bee_debug = 0;
  int save_proj = 0;
  
  for(Int_t i = 1; i != argc; i++){
    switch(argv[i][1]){
    case 'c':
      flag_pos_corr = atoi(&argv[i][2]); 
      break;
    case 'd':
      datatier = atoi(&argv[i][2]);
      break;
    case 'l':
      save_light = atoi(&argv[i][2]);
      break;
    case 'b':
      bee_debug = atoi(&argv[i][2]);
      break;
    case 'p':
      save_proj = atoi(&argv[i][2]);
    }
  }
  bool flag_timestamp = false;
  
  int flag_data = 1; // data
  if (datatier==1 || datatier==2) flag_data=0; // overlay, full mc
  bool flag_match_data = true;
  if (datatier == 2) flag_match_data = false; // if MC we do not take into account the dead PMT

  
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
  double eventTime;
  int run_no, subrun_no, event_no;
  int time_offset;
  int nrebin;
  int frame_length;
  int eve_num;
  float unit_dis;

  // get electron lifetime
  
  Float_t elifetime = 1000; // large number 
  if (Trun->GetBranch("elifetime")){
    Trun->SetBranchAddress("elifetime",&elifetime);
  }

  
  std::vector<int> *timesliceId = new std::vector<int>;
  std::vector<std::vector<int>> *timesliceChannel = new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *raw_charge = new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *raw_charge_err = new std::vector<std::vector<int>>;
  
  Trun->SetBranchAddress("eventNo",&event_no);
  Trun->SetBranchAddress("runNo",&run_no);
  Trun->SetBranchAddress("subRunNo",&subrun_no);
  Trun->SetBranchAddress("eventTime",&eventTime);
  Trun->SetBranchAddress("unit_dis",&unit_dis);
  Trun->SetBranchAddress("frame_length",&frame_length);
  Trun->SetBranchAddress("eve_num",&eve_num);
  Trun->SetBranchAddress("nrebin",&nrebin);
  Trun->SetBranchAddress("time_offset",&time_offset);
  //triggerbits for Beam window selection
  unsigned int triggerbits;
  Trun->SetBranchAddress("triggerBits",&triggerbits); 
  
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
  
  
  
  
  //std::cout << nrebin << " " << time_offset << std::endl;
  
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

  if (elifetime < 1000){
    // read the variable from the Trun tree ...
    mp.set_electron_lifetime(elifetime);

    std::cout << "Electron Lifetime Read in: " << elifetime << " ms" << std::endl;
  }
  
  
  
  
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
								   3*units::cm, 117*units::cm, -116*units::cm, 0*units::cm, 1037*units::cm, 0*units::cm, 256*units::cm, flag_data);
  
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
    T_bad_ch->SetBranchAddress("chid",&chid);
    T_bad_ch->SetBranchAddress("plane",&plane);
    T_bad_ch->SetBranchAddress("start_time",&start_time);
    T_bad_ch->SetBranchAddress("end_time",&end_time);
    
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
  

  // for (size_t i = 0 ;i!=live_clusters.size();i++){
  //   ToyPointCloud *pcloud = live_clusters.at(i)->get_point_cloud();
  //   Point test_p(236.3*units::cm,13.4*units::cm,392.6*units::cm);
  //   if (pcloud->get_closest_dis(test_p) < 1*units::cm)
  //     std::cout << i << " " << live_clusters.at(i)->get_cluster_id() << " " << pcloud->get_closest_dis(test_p)/units::cm << std::endl;
  // }
  
  
  //WCP2dToy::Clustering_live_dead(live_clusters, dead_clusters);
  // cerr << em("Clustering live and dead clusters") << std::endl;
  map_cluster_cluster_vec group_clusters = WCP2dToy::Clustering_jump_gap_cosmics(live_clusters, dead_clusters,dead_u_index, dead_v_index, dead_w_index, global_point_cloud, ct_point_cloud);
  cout << em("Clustering to jump gap in cosmics") << std::endl;

  
  // // reset live_clusters;
  // live_clusters.clear();
  // for (auto it = group_clusters.begin(); it!= group_clusters.end(); it++){
  //   live_clusters.push_back(it->first);
  // }
  
  
  
   
   TFile *file1 = new TFile(Form("match_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");
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

    //  for (size_t j = 0; j!= live_clusters.size(); j++){
    //    SMGCSelection& mcells = live_clusters.at(j)->get_mcells();
    //    ncluster = live_clusters.at(j)->get_cluster_id();
    //    for (size_t i=0;i!=mcells.size();i++){
    //     PointVector ps = mcells.at(i)->get_sampling_points();
    //     int time_slice = mcells.at(i)->GetTimeSlice();
    //     if (ps.size()==0){
    //       std::cout << "zero sampling points!" << std::endl;
    //     }else{
    //       q = mcells.at(i)->get_q() / ps.size();
    //       nq = ps.size();
    //       for (int k=0;k!=ps.size();k++){
    //         x = ps.at(k).x/units::cm ;
    //         y = ps.at(k).y/units::cm;
    //         z = ps.at(k).z/units::cm;
    //         T_cluster->Fill();
    //       }
    //     }
    //    }
    //  }

    for (auto it = group_clusters.begin(); it != group_clusters.end(); ++it) {
       PR3DCluster* main_cluster = it->first;
       SMGCSelection& mcells = main_cluster->get_mcells();
       ncluster = main_cluster->get_cluster_id();
       for (size_t i = 0; i != mcells.size(); i++) {
        PointVector ps = mcells.at(i)->get_sampling_points();
        int time_slice = mcells.at(i)->GetTimeSlice();
        if (ps.size() == 0) {
          std::cout << "zero sampling points!" << std::endl;
        } else {
          q = mcells.at(i)->get_q() / ps.size();
          nq = ps.size();
          for (int k = 0; k != ps.size(); k++) {
          x = ps.at(k).x / units::cm;
          y = ps.at(k).y / units::cm;
          z = ps.at(k).z / units::cm;
          T_cluster->Fill();
          }
        }
       }
       for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++) {
        PR3DCluster* temp_cluster = (*it1).first;
        if (temp_cluster == main_cluster) continue;
        SMGCSelection& mcells = temp_cluster->get_mcells();
        for (size_t i = 0; i != mcells.size(); i++) {
          PointVector ps = mcells.at(i)->get_sampling_points();
          int time_slice = mcells.at(i)->GetTimeSlice();
          if (ps.size() == 0) {
            std::cout << "zero sampling points!" << std::endl;
          } else {
            q = mcells.at(i)->get_q() / ps.size();
            nq = ps.size();
            for (int k = 0; k != ps.size(); k++) {
              x = ps.at(k).x / units::cm;
              y = ps.at(k).y / units::cm;
              z = ps.at(k).z / units::cm;
              T_cluster->Fill();
            }
          }
       }
    }
   }
   
   }

   Trun->CloneTree(-1,"fast");


  
	


  
  // T_proj_data->Fill();
   
  file1->Write();
  file1->Close();

}

