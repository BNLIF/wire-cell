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

#include "TH1F.h"
#include "TH2F.h"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"

using namespace WCP;
using namespace std;

/* 
Mix matchng results within T_match tree 
ensuring that correct matches are not paired
--> provides a sample to test for LM event identification 
*/


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
								   3*units::cm);
  

  // {
  //   Point test_p(166.691*units::cm,-104.599*units::cm, 248.05*units::cm);
  //   // Point test_p(166.691*units::cm,-104.599*units::cm, 248.05*units::cm);
  //   //Point test_p(50.691*units::cm, 50.599*units::cm, 248.05*units::cm);
  //   double offset_x = 0.0636516*units::cm;
  //   fid->inside_fiducial_volume(test_p,offset_x);
  // }
  
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
  
  // std::cout << live_clusters.size() << std::endl;
  // for (size_t i=0;i!=live_clusters.size();i++){
  //   std::cout << live_clusters.at(i)->get_cluster_id() << " " 
  // 	       << live_clusters.at(i)->get_num_mcells() << " "
  // 	       << live_clusters.at(i)->get_num_time_slices() << std::endl;
  // }
  // std::cout << dead_clusters.size() << std::endl;
  for (size_t i=0;i!=dead_clusters.size();i++){
    dead_clusters.at(i)->Remove_duplicated_mcells();
    // std::cout << dead_clusters.at(i)->get_cluster_id() << " " 
    // 	       << dead_clusters.at(i)->get_num_mcells() << " "
    // 	       << dead_clusters.at(i)->get_num_time_slices() << std::endl;
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
   // {
   //   WCP::Point p(30.0*units::cm,30*units::cm,30*units::cm);
   //   WCP::Point p1(110.0*units::cm,0*units::cm,0*units::cm);
     
   //   std::cout << fid->inside_fiducial_volume(p) << " " << fid->inside_fiducial_volume(p1) << std::endl;
   //   WCP::Point p2(-100*units::cm, 42.5*units::cm+3*units::cm, 738.4*units::cm);
   //   std::cout << fid->inside_dead_region(p2) << std::endl;
     
   //   // for (int i=0;i!=1000;i++){
   //   //   for (int j=0;j!=1000;j++){
   //   // 	 //	 
   //   // 	 WCP::Point p2(302.8*units::cm, -116*units::cm + 233*units::cm/1000.*j , 1037*units::cm/1000.*i);
   //   // 	 if (fid->inside_dead_region(p2))
   //   // 	   std::cout << "Xin: " << p2.y/units::cm << " " << p2.z/units::cm << std::endl;
	       	 
   //   //   }
   //   // }
   // }
			
   
   
   ToyCTPointCloud ct_point_cloud(0,2399,2400,4799,4800,8255, // channel range
				  offset_t, -first_u_dis/pitch_u, -first_v_dis/pitch_v, -first_w_dis/pitch_w, // offset
				  1./time_slice_width, 1./pitch_u, 1./pitch_v, 1./pitch_w, // slope
				  angle_u,angle_v,angle_w// angle
				  );
   ct_point_cloud.AddPoints(timesliceId,timesliceChannel,raw_charge,raw_charge_err);
   ct_point_cloud.AddDeadChs(dead_u_index, dead_v_index, dead_w_index);
   ct_point_cloud.build_kdtree_index();


   
   
   // // test the usage of this CT point cloud
   // {
   //   std::cout << live_clusters.at(0)->get_mcells().front()->get_sampling_points().front().x/units::cm << " " << live_clusters.at(0)->get_mcells().front()->get_sampling_points().front().y/units::cm << " " << live_clusters.at(0)->get_mcells().front()->get_sampling_points().front().z/units::cm << std::endl;
   //   std::cout << live_clusters.at(0)->get_mcells().front()->GetTimeSlice() << std::endl;
   //   for (auto it = live_clusters.at(0)->get_mcells().front()->get_uwires().begin(); it!= live_clusters.at(0)->get_mcells().front()->get_uwires().end(); it++){
   //     std::cout << "U: " << (*it)->index() << " " << live_clusters.at(0)->get_mcells().front()->Get_Wire_Charge(*it) <<std::endl;
   //   }
   //   for (auto it = live_clusters.at(0)->get_mcells().front()->get_vwires().begin(); it!= live_clusters.at(0)->get_mcells().front()->get_vwires().end(); it++){
   //     std::cout << "V: " << 2400+(*it)->index() << " " << live_clusters.at(0)->get_mcells().front()->Get_Wire_Charge(*it) <<std::endl;
   //   }
   //   for (auto it = live_clusters.at(0)->get_mcells().front()->get_wwires().begin(); it!= live_clusters.at(0)->get_mcells().front()->get_wwires().end(); it++){
   //     std::cout << "W: " << 4800+(*it)->index() << " " << live_clusters.at(0)->get_mcells().front()->Get_Wire_Charge(*it) << std::endl;
   //   }

   //   ct_point_cloud.Print(live_clusters.at(0)->get_mcells().front()->get_sampling_points().front());
   //   std::cout << ct_point_cloud.get_num_points(0) << " " << ct_point_cloud.get_num_points(1) << " " << ct_point_cloud.get_num_points(2) << std::endl;
     
   //   WCP::CTPointCloud<double> nearby_points = ct_point_cloud.get_closest_points(live_clusters.at(0)->get_mcells().front()->get_sampling_points().front(),1*units::cm,0);
   //   for (size_t i=0;i!=nearby_points.pts.size();i++){
   //     std::cout << "U1: " << nearby_points.pts.at(i).channel << " " << nearby_points.pts.at(i).time_slice << " " << nearby_points.pts.at(i).charge << std::endl;
   //   }
   //   nearby_points = ct_point_cloud.get_closest_points(live_clusters.at(0)->get_mcells().front()->get_sampling_points().front(),1*units::cm,1);
   //   for (size_t i=0;i!=nearby_points.pts.size();i++){
   //     std::cout << "V1: " << nearby_points.pts.at(i).channel << " " << nearby_points.pts.at(i).time_slice << " " << nearby_points.pts.at(i).charge << std::endl;
   //   }
   //   nearby_points = ct_point_cloud.get_closest_points(live_clusters.at(0)->get_mcells().front()->get_sampling_points().front(),1*units::cm,2);
   //   for (size_t i=0;i!=nearby_points.pts.size();i++){
   //     std::cout << "W1: " << nearby_points.pts.at(i).channel << " " << nearby_points.pts.at(i).time_slice << " " << nearby_points.pts.at(i).charge << std::endl;
   //   }
   // }
   
   
   // finish creating global CT point cloud 

   
   
   // WCP2dToy::Clustering_live_dead(live_clusters, dead_clusters);
   // cerr << em("Clustering live and dead clusters") << std::endl;

   std::map<PR3DCluster*,std::vector<std::pair<PR3DCluster*,double>>> group_clusters = WCP2dToy::Clustering_jump_gap_cosmics(live_clusters, dead_clusters,dead_u_index, dead_v_index, dead_w_index, global_point_cloud, ct_point_cloud);
   cout << em("Clustering to jump gap in cosmics") << std::endl;
   
   // // need to further cluster things ...
   // std::map<PR3DCluster*,std::vector<std::pair<PR3DCluster*,double>>> group_clusters =  WCP2dToy::Clustering_isolated(live_clusters);
   // cerr << em("Clustering isolated") << std::endl;
   
   
   
  
   
   // processing light information
   //const char* root_file = argv[3];
   //  WCP2dToy::uBooNE_light_reco uboone_flash(root_file);
   WCP2dToy::ToyLightReco uboone_flash(filename, 1); // 1: imagingoutput, "Trun"; default/not specfified: path "Event/Sim"
   uboone_flash.load_event_raw(0);
   cout << em("flash reconstruction") << std::endl;

   // prepare light matching ....
   WCP::OpflashSelection& flashes = uboone_flash.get_flashes();
   for (size_t i=0;i!=flashes.size(); i++){
     flashes.at(i)->set_flash_id(i);
   }



   

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
  
  //

   
   //   cout<<"BUGGGG"<<endl;
   

   //   std::vector<std::tuple<WCP::PR3DCluster*, WCP::Opflash*, double, std::vector<double>>> matched_results = WCP2dToy::tpc_light_match(time_offset,nrebin,group_clusters,flashes);
  //FlashTPCBundleSelection matched_bundles = WCP2dToy::tpc_light_match(time_offset,nrebin,group_clusters,flashes, run_no, true); // data
  WCP::Photon_Library pl(run_no,true);
  FlashTPCBundleSelection matched_bundles = WCP2dToy::tpc_light_match(time_offset,nrebin,&pl,group_clusters,flashes, run_no, true); // data
   cout << em("TPC Light Matching") << std::endl;

   // create the live clusters ...
   //std::cout << live_clusters.size() << std::endl;
   live_clusters.clear();

   for (auto it = matched_bundles.begin(); it!= matched_bundles.end(); it++){
     FlashTPCBundle *bundle = *it;
     PR3DCluster *main_cluster = bundle->get_main_cluster();
     live_clusters.push_back(main_cluster);
     for (auto it1 = bundle->get_other_clusters().begin(); it1!=bundle->get_other_clusters().end();it1++){
       live_clusters.push_back(*it1);
     }
   }
   //std::cout << live_clusters.size() << std::endl;

   std::map<PR3DCluster*, PR3DCluster*> old_new_cluster_map;
   
   for (size_t i=0;i!=live_clusters.size();i++){

     //if (live_clusters.at(i)->get_cluster_id()!=3) continue;

     //std::cout << i << " " << live_clusters.at(i)->get_cluster_id() << " " << live_clusters.at(i)->get_mcells().size() << " " << live_clusters.at(i)->get_num_time_slices() << std::endl;
     
     live_clusters.at(i)->Create_graph(ct_point_cloud);

     std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = live_clusters.at(i)->get_highest_lowest_wcps();
     // std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = live_clusters.at(i)->get_extreme_wcps();

     // std::cout << wcps.first.x/units::cm << " " << wcps.first.y/units::cm << " " << wcps.first.z/units::cm << " " << wcps.second.x/units::cm << " " << wcps.second.y/units::cm << " " << wcps.second.z/units::cm << std::endl;
     
     live_clusters.at(i)->dijkstra_shortest_paths(wcps.first);

     //std::cout << "shortest path start point" << std::endl;
     
     live_clusters.at(i)->cal_shortest_path(wcps.second);

     {
       // temp ... 
       PR3DCluster *new_cluster = WCP2dToy::Improve_PR3DCluster(live_clusters.at(i),ct_point_cloud, gds);
       WCP2dToy::calc_sampling_points(gds,new_cluster,nrebin, frame_length, unit_dis);
       new_cluster->Create_point_cloud();
       old_new_cluster_map[live_clusters.at(i)] = new_cluster;

       new_cluster->Create_graph(ct_point_cloud);
       std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> new_wcps = new_cluster->get_highest_lowest_wcps();
       new_cluster->dijkstra_shortest_paths(new_wcps.first);
       new_cluster->cal_shortest_path(new_wcps.second);
     }
     
     //std::cout << "shortest path end point" << std::endl;
     
     live_clusters.at(i)->fine_tracking(global_wc_map);

     //std::cout << "fine tracking" << std::endl;
     
     live_clusters.at(i)->collect_charge_trajectory(ct_point_cloud);

     //std::cout << "Collect points" << std::endl;
   }
   
   cerr << em("Create Graph in all clusters") << std::endl;

   
   
   TFile *file1 = new TFile(Form("lmTypeA_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");
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
     //PR3DCluster *main_cluster = bundle->get_main_cluster();
     FlashTPCBundle *bundleAlt = *(it+1);
     PR3DCluster *main_cluster = bundleAlt->get_main_cluster();
     cluster_length = -1;
     
     if (flash!=0){
       auto it1 = find(flashes.begin(),flashes.end(),flash);
       //std::cout<<"flash ID:   "<<flash->get_flash_id()<<std::endl;
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
       double offset_x = (flash->get_time() - time_offset)*2./nrebin*time_slice_width;
       if (fid->check_tgm(bundle,offset_x, ct_point_cloud,old_new_cluster_map))
	 event_type |= 1UL << 3; // 3rd bit for TGM

       /* X.Q. cuts-based LM ID */
       int temp_flag = fid->check_LM(bundle,cluster_length);
       if (temp_flag==1){
	 event_type |= 1UL << 4; // 4th bit for low energy ...
       }else if (temp_flag==2){
	 event_type |= 1UL << 1; // 1st bit for light mismatch ...
       }

       /* B.E.R.M. cuts-based LM ID */
       int temp_flag_cut = fid->check_LM_cuts(bundle,cluster_length);
       if(temp_flag!=temp_flag_cut){
	 std::cout<<"temp_flag: XQ|BERM  "<<temp_flag<<"|"<<temp_flag_cut<<std::endl;
	 if (temp_flag_cut==1){
	   event_type |= 1UL << 5; // 5th bit for low energy ...
	 }else if (temp_flag_cut==2){
	   event_type |= 1UL << 6; // 6th bit for light mismatch ...
	   //event_type |= 1UL; // 5
	 }
       }
       
       /* B.E.R.M. BDT-based LM ID */
       int temp_flag_bdt = fid->check_LM_bdt(bundle,cluster_length);
       if(temp_flag!=temp_flag_bdt || temp_flag_cut!=temp_flag_bdt){
	 std::cout<<"temp_flag: XQ|BERM|BERM "<<temp_flag<<"|"<<temp_flag_cut<<"|"<<temp_flag_bdt<<std::endl;
     	 event_type = 0;
	 if (temp_flag_bdt==1){
	   event_type |= 1UL << 7; // 7th bit for low energy ...
	 }else if (temp_flag_bdt==2){
	   event_type |= 1UL << 8; // 8th bit for light mismatch ...
	   //int temp = 0;
	   //temp |= 1UL << 1;
	   //event_type |= temp; // 6
	 }
       }
       std::cout<<"event designation: "<<temp_flag<<" "<<temp_flag_cut<<" "<<temp_flag_bdt<<std::endl;
       
     }
     
     T_match->Fill();
    }

   cerr << em("test TGM and LM") << std::endl;

   
   // for (auto it = matched_results.begin(); it!=matched_results.end(); it++){
   //   Opflash *flash = std::get<1>(*it);
   //   PR3DCluster *main_cluster = std::get<0>(*it);
   //   if (flash!=0){
   //     auto it1 = find(flashes.begin(),flashes.end(),flash);
   //     // flash_id = it1 - flashes.begin();
   //     flash_id = flash->get_flash_id();
   //     strength = std::get<2>(*it);
   //     std::vector<double> temp = std::get<3>(*it);
   //     for (int i=0;i!=32;i++){
   //   	 pe_pred[i] = temp.at(i);
   //   	 pe_meas[i] = flash->get_PE(i);
   //   	 pe_meas_err[i] = flash->get_PE_err(i);
   //     }
   //   }else{
   //     flash_id = -1;
   //     strength = 0;
   //     for (int i=0;i!=32;i++){
   //   	 pe_pred[i] = 0;
   //   	 pe_meas[i] = 0;
   //   	 pe_meas_err[i] = 0.;
   //     }
   //   }
   //   ncluster = main_cluster->get_cluster_id();
   //   T_match->Fill();
   //   //ncluster++;
   // }
   
   
   
   
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
   
   T_cluster = new TTree("T_cluster","T_cluster");
   T_cluster->Branch("cluster_id",&ncluster,"cluster_id/I");
   T_cluster->Branch("x",&x,"x/D");
   T_cluster->Branch("y",&y,"y/D");
   T_cluster->Branch("z",&z,"z/D");
   T_cluster->Branch("q",&q,"q/D");
   T_cluster->Branch("nq",&nq,"nq/D");
   T_cluster->SetDirectory(file1);

   TTree *T_rec = new TTree("T_rec","T_rec");
   T_rec->Branch("x",&x,"x/D");
   T_rec->Branch("y",&y,"y/D");
   T_rec->Branch("z",&z,"z/D");
   T_rec->SetDirectory(file1);

   Double_t charge_save=1, ncharge_save=1, chi2_save=1, ndf_save=1;
   TTree *t_rec_charge = new TTree("T_rec_charge","T_rec_charge");
   t_rec_charge->SetDirectory(file1);
   t_rec_charge->Branch("x",&x,"x/D");
   t_rec_charge->Branch("y",&y,"y/D");
   t_rec_charge->Branch("z",&z,"z/D");
   t_rec_charge->Branch("q",&charge_save,"q/D");
   t_rec_charge->Branch("nq",&ncharge_save,"nq/D");
   t_rec_charge->Branch("chi2",&chi2_save,"chi2/D");
   t_rec_charge->Branch("ndf",&ndf_save,"ndf/D");
   
   
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
       //     for (auto it1 = group_clusters[main_cluster].begin(); it1!=group_clusters[main_cluster].end(); it1++){
       //temp_clusters.push_back((*it1).first);
       //std::cout << (*it1).second/units::cm << std::endl;
     }
     for (size_t j = 0; j!= temp_clusters.size(); j++){

       // show individual clusters ... 
       //ncluster = temp_clusters.at(j)->get_cluster_id();
       
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
     std::list<WCPointCloud<double>::WCPoint>& wcps_list = new_cluster->get_path_wcps();
     //ncluster = -1 * ncluster-100;
     for (auto it = wcps_list.begin(); it!=wcps_list.end(); it++){
       x = (*it).x/units::cm;
       y = (*it).y/units::cm;
       z = (*it).z/units::cm;
       T_rec->Fill();
     }

     // PR3DCluster *new_cluster = old_new_cluster_map[live_clusters.at(j)];
     // SMGCSelection& mcells = new_cluster->get_mcells();
     // //ncluster = temp_clusters.at(0)->get_cluster_id();
     // for (size_t i=0;i!=mcells.size();i++){
     //   PointVector ps = mcells.at(i)->get_sampling_points();
     //   int time_slice = mcells.at(i)->GetTimeSlice();
     //   for (int k=0;k!=ps.size();k++){
     // 	 x = ps.at(k).x/units::cm ;
     // 	 y = ps.at(k).y/units::cm;
     // 	 z = ps.at(k).z/units::cm;
     // 	 T_rec->Fill();
     //   }
     // }
     
     PointVector& pts = live_clusters.at(j)->get_fine_tracking_path();
     for (size_t i=0; i!=pts.size(); i++){
       x = pts.at(i).x/units::cm;
       y = pts.at(i).y/units::cm;
       z = pts.at(i).z/units::cm;
       t_rec_charge->Fill();
     }

     

     // // save mcells
     // std::list<SlimMergeGeomCell*>& mcells_list = live_clusters.at(j)->get_path_mcells();
     // ncluster = -1 * ncluster-100;
     // for (auto it = mcells_list.begin(); it!=mcells_list.end(); it++){
     //   Point p = (*it)->center();
     //   x = p.x/units::cm;
     //   y = p.y/units::cm;
     //   z = p.z/units::cm;
     //   T_cluster->Fill();
     // }
     
     
     // if (live_clusters.at(j)->get_num_mcells()>30){
     //   // add PCA axis point
     //   Vector center = live_clusters.at(j)->get_center();
     //   Vector dir = live_clusters.at(j)->get_PCA_axis(0);
     //   for (int i=-200;i!=200;i++){
     // 	 x = (center.x + dir.x *(i*units::cm) )/units::cm;
     // 	 y = (center.y + dir.y *(i*units::cm) )/units::cm;
     // 	 z = (center.z + dir.z *(i*units::cm) )/units::cm;
     // 	 T_cluster->Fill();
     //   }
     // }
    }

   // ncluster = 0;
   // for (auto it = dead_live_cluster_mapping.begin(); it!= dead_live_cluster_mapping.end(); it++){
   //   std::vector<PR3DCluster*> clusters = (*it).second;
   //   if (clusters.size()>1){
   //     //std::cout << clusters.size() << std::endl;
   //     for (auto it1 = clusters.begin(); it1!=clusters.end(); it1++){
   // 	 PR3DCluster* cluster = (*it1);
   // 	 ncluster = cluster->get_cluster_id();
   // 	 SMGCSelection& mcells = cluster->get_mcells();
   // 	 for (size_t i=0;i!=mcells.size();i++){
   // 	   PointVector ps = mcells.at(i)->get_sampling_points();
   // 	   for (int k=0;k!=ps.size();k++){
   // 	     x = ps.at(k).x/units::cm;//time_slice*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
   // 	     y = ps.at(k).y/units::cm;
   // 	     z = ps.at(k).z/units::cm;
   // 	     T_cluster->Fill();
   // 	   }
   // 	 }
   //     }
   //   }
   //   // ncluster++;
   // }
   // for (auto it = dead_live_mcells_mapping.begin(); it!= dead_live_mcells_mapping.end(); it++){
   //   std::vector<std::vector<SlimMergeGeomCell*>> mcellss = (*it).second;
   //   // std::cout << mcellss.size() << std::endl;
   //   if (mcellss.size()>1){
   //     for (auto it1 = mcellss.begin(); it1!=mcellss.end(); it1++){
   // 	 std::vector<SlimMergeGeomCell*> mcells = (*it1);
   // 	 for (size_t i=0;i!=mcells.size();i++){
   // 	   PointVector ps = mcells.at(i)->get_sampling_points();
   // 	   for (int k=0;k!=ps.size();k++){
   // 	     x = ps.at(k).x/units::cm;//time_slice*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
   // 	     y = ps.at(k).y/units::cm;
   // 	     z = ps.at(k).z/units::cm;
   // 	     T_cluster->Fill();
   // 	   }
   // 	 }
   //     }
   //   }
   //   ncluster++;
   // }
   
   Trun->CloneTree()->Write();
   //Trun->CloneTree();
   //Trun->SetDirectory(file1);

   
   TTree *t1 = new TTree("T_data","T_data");
   t1->SetDirectory(file1);

   TClonesArray* op_wf = new TClonesArray("TH1S");
   std::vector<short> *op_femch = new std::vector<short>;
   std::vector<double> *op_timestamp = new std::vector<double>;

  t1->Branch("op_gain",&op_gain);
  t1->Branch("op_gainerror",&op_gainerror);
  t1->Branch("op_femch",&op_femch);
  t1->Branch("op_timestamp",&op_timestamp);
  t1->Branch("op_wf",&op_wf,256000,0);
  t1->Branch("triggerTime",&triggerTime);
  
  t1->Branch("runNo",&run_no);
  t1->Branch("subRunNo",&subrun_no);
  t1->Branch("eventNo",&event_no);

  op_wf = uboone_flash.get_rawWfm();
  op_femch = uboone_flash.get_rawChan();
  op_timestamp = uboone_flash.get_rawTimestamp();

  t1->Fill();

  TH2F *h1 = new TH2F("hraw","hraw",1500,0,1500,32,0,32);
  TH2F *h2 = new TH2F("hdecon","hdecon",250,0,250,32,0,32);
  TH2F *h3 = new TH2F("hl1","hl1",250,0,250,32,0,32);
  h1->SetDirectory(file1);
  h2->SetDirectory(file1);
  h3->SetDirectory(file1);
  for (int i=0;i!=32;i++){
    TH1F *h10 = uboone_flash.get_raw_hist(i);
    TH1F *h20 = uboone_flash.get_decon_hist(i);
    TH1F *h30 = uboone_flash.get_l1_hist(i);
    for (int j=0;j!=1500;j++){
      h1->SetBinContent(j+1,i+1,h10->GetBinContent(j+1));
    }
    for (int j=0;j!=250;j++){
      h2->SetBinContent(j+1,i+1,h20->GetBinContent(j+1));
      h3->SetBinContent(j+1,i+1,h30->GetBinContent(j+1));
    }
  }

  TH1F *h_totPE = (TH1F*)uboone_flash.get_totPE()->Clone("totPE");
  TH1F *h_mult = (TH1F*)uboone_flash.get_mult()->Clone("mult");
  TH1F *h_l1_mult = (TH1F*)uboone_flash.get_l1_mult()->Clone("l1_mult");
  TH1F *h_l1_totPE = (TH1F*)uboone_flash.get_l1_totPE()->Clone("l1_totPE");

  h_totPE->SetDirectory(file1);
  h_mult->SetDirectory(file1);
  h_l1_mult->SetDirectory(file1);
  h_l1_totPE->SetDirectory(file1);

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
  //  T_flash->Branch("matched_tpc_ids",&matched_tpc_ids);

  
  
  
  for (auto it = flashes.begin(); it!=flashes.end(); it++){
    fired_channels.clear();
    //it - flashes.begin();
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

 
 
  
  
  // now save the projected charge information ... 
  TTree *T_proj = new TTree("T_proj","T_proj");
  std::vector<int> *proj_cluster_id = new std::vector<int>;
  std::vector<std::vector<int>> *proj_cluster_channel = new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *proj_cluster_timeslice= new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *proj_cluster_charge= new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *proj_cluster_charge_err= new std::vector<std::vector<int>>;
  T_proj->Branch("cluster_id",&proj_cluster_id);
  T_proj->Branch("channel",&proj_cluster_channel);
  T_proj->Branch("time_slice",&proj_cluster_timeslice);
  T_proj->Branch("charge",&proj_cluster_charge);
  T_proj->Branch("charge_err",&proj_cluster_charge_err);
  
  T_proj->SetDirectory(file1);
  
  for (auto it = matched_bundles.begin(); it!= matched_bundles.end(); it++){
    FlashTPCBundle *bundle = *it;
    PR3DCluster *main_cluster = bundle->get_main_cluster();
    Opflash *flash = bundle->get_flash();
    //  if (flash!=0){
      
      // now prepare saving it
    int cluster_id = main_cluster->get_cluster_id();
    
    PR3DClusterSelection temp_clusters;
    temp_clusters.push_back(main_cluster);
    for (auto it1 = bundle->get_other_clusters().begin(); it1!=bundle->get_other_clusters().end();it1++){
      temp_clusters.push_back(*it1);
       //for (auto it1 = group_clusters[main_cluster].begin(); it1!=group_clusters[main_cluster].end(); it1++){
       //  temp_clusters.push_back((*it1).first);
    }
    
    std::vector<int> proj_channel;
    std::vector<int> proj_timeslice;
    std::vector<int> proj_charge;
    std::vector<int> proj_charge_err;
    std::vector<int> proj_flag;
    
    for (size_t j = 0; j!= temp_clusters.size(); j++){
      PR3DCluster *cluster = temp_clusters.at(j);
      cluster->get_projection(proj_channel,proj_timeslice,proj_charge, proj_charge_err, proj_flag, global_wc_map);
    }
    proj_cluster_id->push_back(cluster_id);
    proj_cluster_channel->push_back(proj_channel);
    proj_cluster_timeslice->push_back(proj_timeslice);
    proj_cluster_charge->push_back(proj_charge);
    proj_cluster_charge_err->push_back(proj_charge_err);
    
      // }
  }
  T_proj->Fill();

   
  file1->Write();
  file1->Close();
}
