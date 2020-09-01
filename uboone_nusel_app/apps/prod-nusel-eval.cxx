#include "WCPSst/GeomDataSource.h"
#include "WCPData/PR3DCluster.h"
#include "WCPData/SlimMergeGeomCell.h"
#include "WCPData/TPCParams.h"
#include "WCPData/Singleton.h"
#include "WCPData/ToyCTPointCloud.h"
#include "WCP2dToy/ExecMon.h"
#include "WCP2dToy/CalcPoints.h"
#include "WCP2dToy/ToyClustering.h"
#include "WCP2dToy/ToyLightReco.h"
#include "WCP2dToy/ToyMatching.h"
#include "WCP2dToy/ToyFiducial.h"
#include "WCP2dToy/ImprovePR3DCluster.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"

using namespace WCP;
using namespace std;

int main(int argc, char* argv[])
{
  if (argc < 3) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/imaging.root -d[0,1,2] -p[0,1]" << endl;
    return 1;
  }
  TH1::AddDirectory(kFALSE);
  
  int flag_pos_corr = 0; // correct X position after matching ...
  int port = -1; // saved if =1
  int datatier = 0; // data=0, overlay=1, full mc=2
  for(Int_t i = 1; i != argc; i++){
    switch(argv[i][1]){
    case 'c':
      flag_pos_corr = atoi(&argv[i][2]); 
      break;
    case 'd':
      datatier = atoi(&argv[i][2]);
      break;
    case 'p':
      port = atoi(&argv[i][2]);
      break;
    }
  }
  bool flag_match_data = true;
  bool flag_timestamp = false;
  
  if(datatier==2)
    flag_match_data = false;
  ExecMon em("starting");
  cout << em("load geometry") << endl;

  WCPSst::GeomDataSource gds(argv[1]);
  std::vector<double> ex = gds.extent();
  const GeomWire *uwire = gds.by_planeindex(WirePlaneType_t(0),0);
  const GeomWire *vwire = gds.by_planeindex(WirePlaneType_t(1),0);
  const GeomWire *wwire = gds.by_planeindex(WirePlaneType_t(2),0);
  double first_u_dis = gds.wire_dist(*uwire) ; // first U wire center ...
  double first_v_dis = gds.wire_dist(*vwire) ; // first V wire center ...
  double first_w_dis = gds.wire_dist(*wwire) ; // first W wire center ... 
  
  TString filename = argv[2];
  //TFile *file = new TFile(filename);
  TFile *file = TFile::Open(filename); //enable xrootd url
  TTree *Trun = (TTree*)file->Get("Trun");
  double eventTime;
  int run_no, subrun_no, event_no;
  int time_offset;
  int nrebin;
  int frame_length;
  int eve_num;
  float unit_dis;
  unsigned int triggerbits;
  std::vector<int> *timesliceId = new std::vector<int>;
  std::vector<std::vector<int>> *timesliceChannel = new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *raw_charge = new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *raw_charge_err = new std::vector<std::vector<int>>;
  std::vector<float> *op_gain = new std::vector<float>;
  std::vector<float> *op_gainerror = new std::vector<float>;
  double triggerTime;  

  int nmax = 10000;
  int Ntrack;
  float nu_mom[4], nu_pos[4];
  int nu_pdg, nu_ccnc;
  vector<double> *i_nelectrons = new vector<double>;
  vector<double> *i_nphotons = new vector<double>;
  vector<double> *i_time_start = new vector<double>;
  vector<double> *i_time_end = new vector<double>;
  vector<double> *i_x_start = new vector<double>;
  vector<double> *i_x_end = new vector<double>;
  vector<double> *i_y_start = new vector<double>;
  vector<double> *i_y_end = new vector<double>;
  vector<double> *i_z_start = new vector<double>;
  vector<double> *i_z_end = new vector<double>;
  vector<int> *i_pdg = new vector<int>;
  vector<int> *i_trackId = new vector<int>;
  vector<double> *i_energy = new vector<double>;

  Trun->SetBranchAddress("eventNo",&event_no);
  Trun->SetBranchAddress("runNo",&run_no);
  Trun->SetBranchAddress("subRunNo",&subrun_no);
  Trun->SetBranchAddress("eventTime",&eventTime);
  Trun->SetBranchAddress("unit_dis",&unit_dis);
  Trun->SetBranchAddress("frame_length",&frame_length);
  Trun->SetBranchAddress("eve_num",&eve_num);
  Trun->SetBranchAddress("nrebin",&nrebin);
  Trun->SetBranchAddress("time_offset",&time_offset);
  Trun->SetBranchAddress("triggerBits",&triggerbits);
  Trun->SetBranchAddress("timesliceId",&timesliceId);
  Trun->SetBranchAddress("timesliceChannel",&timesliceChannel);
  Trun->SetBranchAddress("raw_charge",&raw_charge);
  Trun->SetBranchAddress("raw_charge_err",&raw_charge_err);
  Trun->SetBranchAddress("op_gain",&op_gain); 
  Trun->SetBranchAddress("op_gainerror",&op_gainerror); 
  Trun->SetBranchAddress("triggerTime",&triggerTime);
  if(datatier==1 || datatier==2){
    Trun->SetBranchAddress("mc_Ntrack",&Ntrack);
    Trun->SetBranchAddress("mc_nu_mom",nu_mom);
    Trun->SetBranchAddress("mc_nu_pos",nu_pos);
    Trun->SetBranchAddress("mc_nu_pdg",&nu_pdg);
    Trun->SetBranchAddress("mc_nu_ccnc",&nu_ccnc);
    Trun->SetBranchAddress("sedi_nelectrons",&i_nelectrons);
    Trun->SetBranchAddress("sedi_nphotons",&i_nphotons);
    Trun->SetBranchAddress("sedi_time_start",&i_time_start);
    Trun->SetBranchAddress("sedi_time_end",&i_time_end);
    Trun->SetBranchAddress("sedi_x_start",&i_x_start);
    Trun->SetBranchAddress("sedi_x_end",&i_x_end);
    Trun->SetBranchAddress("sedi_y_start",&i_y_start);
    Trun->SetBranchAddress("sedi_y_end",&i_y_end);
    Trun->SetBranchAddress("sedi_z_start",&i_z_start);
    Trun->SetBranchAddress("sedi_z_end",&i_z_end);
    Trun->SetBranchAddress("sedi_pdg",&i_pdg);
    Trun->SetBranchAddress("sedi_trackId",&i_trackId);
    Trun->SetBranchAddress("sedi_energy",&i_energy);
  }
  Trun->GetEntry(0);

  
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

  int flag_data=1;
  if(datatier==1 || datatier==2){ flag_data=0; } // overlay and full MC
  WCP2dToy::ToyFiducial *fid = new WCP2dToy::ToyFiducial(3,800,-first_u_dis/pitch_u, -first_v_dis/pitch_v, -first_w_dis/pitch_w,
								   1./time_slice_width, 1./pitch_u, 1./pitch_v, 1./pitch_w, // slope
								   angle_u,angle_v,angle_w,// angle
                                   				   3*units::cm, 117*units::cm, -116*units::cm, 0*units::cm, 1037*units::cm, 0*units::cm, 256*units::cm, flag_data);

  std::vector<int> nuFinalState_pdg;
  std::vector<double> nuFinalState_energy, nuFinalState_x, nuFinalState_y, nuFinalState_z, nuFinalState_nelectrons;
  nuFinalState_pdg.clear();
  nuFinalState_energy.clear();
  nuFinalState_z.clear();
  nuFinalState_y.clear();
  nuFinalState_z.clear();  
  nuFinalState_nelectrons.clear();  
  
  TFile *file1 = new TFile(Form("nuselEval_old_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");
  TTree *T_true;
  int nufs_pdg;
  double nufs_x, nufs_y, nufs_z, nufs_q, nufs_e;
  bool truth_isFC = true;
  bool check_status = true;
  if(datatier==1 || datatier==2){
    T_true = new TTree("T_true","T_true");
    T_true->SetDirectory(file1);
    T_true->Branch("pdg", &nufs_pdg);  
    T_true->Branch("x", &nufs_x);  
    T_true->Branch("y", &nufs_y);  
    T_true->Branch("z", &nufs_z);  
    T_true->Branch("q", &nufs_q);  
    T_true->Branch("e", &nufs_e);  
  
    for(size_t i=0; i<i_time_start->size(); i++){
      if(i_time_start->at(i)>nu_pos[3] && i_time_end->at(i)<nu_pos[3]+1000.){ // 1us window, negligible cosmic muon (5kHz)
	if(i_pdg->at(i)<10000){ // Argon excitation would deposi energy but not invisible to TPC
	  nuFinalState_pdg.push_back(i_pdg->at(i));
	  nuFinalState_energy.push_back(i_energy->at(i));
	  nuFinalState_nelectrons.push_back(i_nelectrons->at(i));
	  double xx = (i_x_start->at(i)+i_x_end->at(i))/2.;
	  double yy = (i_y_start->at(i)+i_y_end->at(i))/2.;
	  double zz = (i_z_start->at(i)+i_z_end->at(i))/2.;
	  nuFinalState_x.push_back( xx );
	  nuFinalState_y.push_back( yy );
	  nuFinalState_z.push_back( zz );
	  
	  nufs_x = xx; //cm
	  nufs_y = yy;
	  nufs_z = zz;
	  nufs_q = i_nelectrons->at(i);
	  nufs_e = i_energy->at(i);
	  nufs_pdg = i_pdg->at(i);
	  T_true->Fill();

	  // is fully contained?
	  WCP::Point pp(xx*units::cm, yy*units::cm, zz*units::cm);
	  if(check_status && !fid->inside_fiducial_volume(pp, 0) && i_energy->at(i)>0.001 /* MeV, dot-like depo */){
	    truth_isFC = false;
	    check_status = false;
	  }
	}
      }
    }
  } // datatier
   
  std::map<int,std::pair<double,double>> dead_u_index;
  std::map<int,std::pair<double,double>> dead_v_index;
  std::map<int,std::pair<double,double>> dead_w_index;
  
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

  // veto 16 channels in U ...
  // for (int i=2080; i!=2096;i++){
  //   if (dead_u_index.find(i)==dead_u_index.end()){
  //     dead_u_index[i] = std::make_pair((0*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) * units::cm-0.1*units::cm, (2400*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) * units::cm+0.1*units::cm);
  //   }
  // }
  //
  
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
   }
   cout << em("Add X, Y, Z points") << std::endl;

   // create global CT point cloud ...
   double first_t_dis = live_clusters.at(0)->get_mcells().front()->GetTimeSlice()*time_slice_width - live_clusters.at(0)->get_mcells().front()->get_sampling_points().front().x;
   double offset_t = first_t_dis/time_slice_width;

   // test the fiducial volume cut 
   // fid->set_offset_t(offset_t);
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

   map_cluster_cluster_vec group_clusters = WCP2dToy::Clustering_jump_gap_cosmics(live_clusters, dead_clusters,dead_u_index, dead_v_index, dead_w_index, global_point_cloud, ct_point_cloud);
   cout << em("Clustering to jump gap in cosmics") << std::endl;
   
   double lowerwindow = 0., upperwindow = 0.;
   if(datatier==0 && ((triggerbits>>11) & 1U)){ lowerwindow=3.1875; upperwindow=4.96876; } // BNB
   if((datatier==0 && ((triggerbits>>9) & 1U)) || datatier==1){ lowerwindow=3.5625; upperwindow=5.34376; } // BNB EXT & overlay
   if(datatier==2){ lowerwindow=3.1718; upperwindow=4.95306; } // full mc
   
   WCP2dToy::ToyLightReco uboone_flash(filename, 1, datatier);
   uboone_flash.load_event_raw(0, lowerwindow, upperwindow);
   cout << em("flash reconstruction") << std::endl;

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
  
   //FlashTPCBundleSelection matched_bundles = WCP2dToy::tpc_light_match(time_offset,nrebin,group_clusters,flashes, run_no, flag_match_data);
   WCP::Photon_Library pl(eventTime,run_no,flag_match_data);
   FlashTPCBundleSelection matched_bundles = WCP2dToy::tpc_light_match(eventTime,time_offset,nrebin,&pl,group_clusters,flashes, run_no, flag_match_data);
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
       PR3DCluster *new_cluster = WCP2dToy::Improve_PR3DCluster(live_clusters.at(i),ct_point_cloud, gds);
       WCP2dToy::calc_sampling_points(gds,new_cluster,nrebin, frame_length, unit_dis);
       new_cluster->Create_point_cloud();
       old_new_cluster_map[live_clusters.at(i)] = new_cluster;

       new_cluster->Create_graph(ct_point_cloud);
       std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> new_wcps = new_cluster->get_highest_lowest_wcps();
       new_cluster->dijkstra_shortest_paths(new_wcps.first);
       new_cluster->cal_shortest_path(new_wcps.second);
     }
     live_clusters.at(i)->fine_tracking(global_wc_map);
     live_clusters.at(i)->collect_charge_trajectory(ct_point_cloud);
   }

   // do the dQ/dx fitting ... 
   //for (auto it = matched_bundles.begin(); it!= matched_bundles.end(); it++){
   //FlashTPCBundle *bundle = *it;
   //Opflash *flash = bundle->get_flash();
   //PR3DCluster *main_cluster = bundle->get_main_cluster();
   //if (flash!=0){
       //
       //       if (flash->get_time() > 2 && flash->get_time() < 6){
       //std::cout << flash->get_time() << std::endl;
       //main_cluster->dQ_dx_fit(global_wc_map, flash->get_time()*units::microsecond);
	 // }
   //}
   //}
  
   cerr << em("Create Graph in all clusters") << std::endl;
   
   TTree *T_eval = new TTree("T_eval","T_eval");
   T_eval->SetDirectory(file1);
   bool truth_isCC=false, truth_isEligible=false, truth_vtxInside=false;
   int truth_nuPdg=-1;
   float truth_nuTime=-1, truth_nuEnergy=-1., truth_energyInside=0., truth_electronInside=0.;
   float truth_vtxX=-1., truth_vtxY=-1., truth_vtxZ=-1.;
   float flash_time=-1, flash_measPe=-1., flash_predPe=-1.;
   bool flash_found=false, match_found=false;
   float match_charge=0., match_energy=0.;
   // 3D matching
   float match_completeness=0./* matched/true */, match_completeness_energy=0., match_purity=0./* matched/reco */;
   float match_purity_xz=0, match_purity_xy=0;
   Int_t match_type=0;
   bool match_notFC_FV=false, match_notFC_SP=false, match_notFC_DC=false, match_isTgm=false, match_isFC=false;

     T_eval->Branch("run", &run_no);
     T_eval->Branch("subrun", &subrun_no);
     T_eval->Branch("event", &event_no);
     T_eval->Branch("flash_found", &flash_found);
     T_eval->Branch("flash_time", &flash_time);
     T_eval->Branch("flash_measPe", &flash_measPe);
     T_eval->Branch("flash_predPe", &flash_predPe);
     T_eval->Branch("match_found", &match_found);
     T_eval->Branch("match_type", &match_type);
     T_eval->Branch("match_isFC", &match_isFC);
     T_eval->Branch("match_isTgm", &match_isTgm);
     T_eval->Branch("match_notFC_FV", &match_notFC_FV);
     T_eval->Branch("match_notFC_SP", &match_notFC_SP);
     T_eval->Branch("match_notFC_DC", &match_notFC_DC);

   if(datatier==1 || datatier==2){
     T_eval->Branch("truth_nuEnergy", &truth_nuEnergy);
     T_eval->Branch("truth_energyInside", &truth_energyInside);
     T_eval->Branch("truth_electronInside", &truth_electronInside);
     T_eval->Branch("truth_nuPdg", &truth_nuPdg);
     T_eval->Branch("truth_isCC", &truth_isCC);
     T_eval->Branch("truth_isEligible", &truth_isEligible);
     T_eval->Branch("truth_isFC", &truth_isFC);
     T_eval->Branch("truth_vtxInside", &truth_vtxInside);
     T_eval->Branch("truth_vtxX", &truth_vtxX);
     T_eval->Branch("truth_vtxY", &truth_vtxY);
     T_eval->Branch("truth_vtxZ", &truth_vtxZ);
     T_eval->Branch("truth_nuTime", &truth_nuTime);
     T_eval->Branch("match_completeness", &match_completeness);
     T_eval->Branch("match_completeness_energy", &match_completeness_energy);
     T_eval->Branch("match_purity", &match_purity);
     T_eval->Branch("match_purity_xz", &match_purity_xz);
     T_eval->Branch("match_purity_xy", &match_purity_xy);
     T_eval->Branch("match_charge", &match_charge);
     T_eval->Branch("match_energy", &match_energy);

     truth_vtxX = nu_pos[0];
     truth_vtxY = nu_pos[1];
     truth_vtxZ = nu_pos[2];
     truth_nuTime = nu_pos[3]/1000.; // us
     truth_nuEnergy = nu_mom[3]*1000; //[MeV]
     if(nu_ccnc==0){ truth_isCC=true; }
     truth_nuPdg = nu_pdg;
     for(int a=0; a<nuFinalState_energy.size(); a++){
       {
	 truth_energyInside += nuFinalState_energy.at(a);	   
	 truth_electronInside += nuFinalState_nelectrons.at(a);	   
       }
     }
     if(truth_vtxX > 5.0 && truth_vtxX < 251.0 &&
	truth_vtxY > -101.0 && truth_vtxY < 101.0 &&
	truth_vtxZ > 10.0 && truth_vtxZ < 1026.0){
       truth_vtxInside=true;
       if(truth_isCC==true){ truth_isEligible=true; }
       else{
	 if(truth_energyInside > 50.){ truth_isEligible=true; }
       }
     }
   } // datatier==1 || datatier==2
   
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

   
   bool is_inTime=false;
   for (auto it = matched_bundles.begin(); it!=matched_bundles.end(); it++){
     FlashTPCBundle *bundle = *it;
     
     Opflash *flash = bundle->get_flash();
     PR3DCluster *main_cluster = bundle->get_main_cluster();
     cluster_length = -1;
     
     if (flash!=0){
       auto it1 = find(flashes.begin(),flashes.end(),flash);
       flash_id = flash->get_flash_id();
       strength = bundle->get_strength();
       if(flash->get_time()>lowerwindow && flash->get_time()<upperwindow){
	 is_inTime=true;
       }
       std::vector<double> temp = bundle->get_pred_pmt_light();
       for (int i=0;i!=32;i++){
     	 pe_pred[i] = temp.at(i);
     	 pe_meas[i] = flash->get_PE(i);
     	 pe_meas_err[i] = flash->get_PE_err(i);
	 if(is_inTime==true){
	   flash_measPe += flash->get_PE(i);
	   flash_predPe += temp.at(i);
	 }
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
     if(is_inTime==true){
       double offset_x = (flash->get_time() - time_offset)*2./nrebin*time_slice_width;
       if (fid->check_tgm(bundle,offset_x, ct_point_cloud,old_new_cluster_map)){ match_isTgm = true; }

       unsigned int fc_breakdown=0;// should be initilized to ZERO
       if (fid->check_fully_contained(bundle,offset_x, ct_point_cloud,old_new_cluster_map, &fc_breakdown)){  
           match_isFC = true;
           match_notFC_FV = false;
           match_notFC_SP = false;
           match_notFC_DC = false;
       }
       else {
           match_isFC = false;
           match_notFC_FV = (fc_breakdown>>2) & 1U; // outside fiducial volume
           match_notFC_SP = (fc_breakdown>>1) & 1U; // SP gaps
           match_notFC_DC = (fc_breakdown) & 1U; // dead regions
       }
       is_inTime=false;
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
   match_type = event_type;
   
   cerr << em("test TGM") << std::endl;
   /*
   TTree *T_bad_ch = (TTree*)file->Get("T_bad_ch");
   //if (T_bad_ch!=0){
   //T_bad_ch->CloneTree()->Write();
   //}
   
   TTree *t_bad = new TTree("T_bad","T_bad");
   //t_bad->SetDirectory(file1);
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
       //t_bad->Fill();
     }
   }
*/   
   TTree *T_cluster ;
   Double_t x,y,z,q,nq;
   Int_t temp_time_slice, ch_u, ch_v, ch_w;
   int flag_main;
   T_cluster = new TTree("T_cluster","T_cluster");
   T_cluster->Branch("cluster_id",&ncluster,"cluster_id/I");
   T_cluster->Branch("flag_main", &flag_main, "flag_main/I");
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
   //T_rec->SetDirectory(file1);
   
   TTree *t_rec_charge = new TTree("T_rec_charge","T_rec_charge");
   //t_rec_charge->SetDirectory(file1);
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
   //t_rec_charge->Branch("pt",&pt,"pt/D");

   TTree *T_proj_data = new TTree("T_proj_data","T_proj_data");
   std::vector<int> *proj_data_cluster_id = new std::vector<int>;
   std::vector<std::vector<int>> *proj_data_cluster_channel = new std::vector<std::vector<int>>;
   std::vector<std::vector<int>> *proj_data_cluster_timeslice= new std::vector<std::vector<int>>;
   std::vector<std::vector<int>> *proj_data_cluster_charge= new std::vector<std::vector<int>>;
   std::vector<std::vector<int>> *proj_data_cluster_charge_err= new std::vector<std::vector<int>>;
   std::vector<std::vector<int>> *proj_data_cluster_charge_pred= new std::vector<std::vector<int>>;
   
   T_proj_data->Branch("cluster_id",&proj_data_cluster_id);
   T_proj_data->Branch("channel",&proj_data_cluster_channel);
   T_proj_data->Branch("time_slice",&proj_data_cluster_timeslice);
   T_proj_data->Branch("charge",&proj_data_cluster_charge);
   T_proj_data->Branch("charge_err",&proj_data_cluster_charge_err);
   T_proj_data->Branch("charge_pred",&proj_data_cluster_charge_pred);
   //T_proj_data->SetDirectory(file1);
   
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
       if(j==0) flag_main = 1; // main cluster flag 
       else flag_main = 0;

       for (size_t i=0;i!=mcells.size();i++){
	 PointVector ps = mcells.at(i)->get_sampling_points();
	 int time_slice = mcells.at(i)->GetTimeSlice();

	 //	 temp_time_slice = time_slice;
	 
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
     std::list<WCPointCloud<double>::WCPoint>& wcps_list = new_cluster->get_path_wcps();
     //ncluster = -1 * ncluster-100;
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
       
       //T_rec->Fill();
     }

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
     
     ndf_save = live_clusters.at(j)->get_cluster_id();

     if (pts.size()==dQ.size()){
       proj_data_cluster_id->push_back(ndf_save);
       std::vector<int> temp_channel;
       std::vector<int> temp_timeslice;
       std::vector<int> temp_charge;
       std::vector<int> temp_charge_err;
       std::vector<int> temp_charge_pred;
       for (auto it = proj_data_u_map.begin(); it!=proj_data_u_map.end(); it++){
	 temp_channel.push_back(it->first.first);
	 temp_timeslice.push_back(it->first.second);
	 temp_charge.push_back(std::get<0>(it->second));
	 temp_charge_err.push_back(std::get<1>(it->second));
	 temp_charge_pred.push_back(std::get<2>(it->second));
       }
       for (auto it = proj_data_v_map.begin(); it!=proj_data_v_map.end(); it++){
	 temp_channel.push_back(it->first.first);
	 temp_timeslice.push_back(it->first.second);
	 temp_charge.push_back(std::get<0>(it->second));
	 temp_charge_err.push_back(std::get<1>(it->second));
	 temp_charge_pred.push_back(std::get<2>(it->second));
       }
       for (auto it = proj_data_w_map.begin(); it!=proj_data_w_map.end(); it++){
	 temp_channel.push_back(it->first.first);
	 temp_timeslice.push_back(it->first.second);
	 temp_charge.push_back(std::get<0>(it->second));
	 temp_charge_err.push_back(std::get<1>(it->second));
	 temp_charge_pred.push_back(std::get<2>(it->second));
       }
       proj_data_cluster_channel->push_back(temp_channel);
       proj_data_cluster_timeslice->push_back(temp_timeslice);
       proj_data_cluster_charge->push_back(temp_charge);
       proj_data_cluster_charge_err->push_back(temp_charge_err);
       proj_data_cluster_charge_pred->push_back(temp_charge_pred);
     }

     for (size_t i=0; i!=pts.size(); i++){
       x = pts.at(i).x/units::cm;
       y = pts.at(i).y/units::cm;
       z = pts.at(i).z/units::cm;
       
       if (pts.size()==dQ.size()){
	 charge_save = dQ.at(i);
	 ncharge_save = dx.at(i)/units::cm;
	 pu = tpu.at(i);
	 pv = tpv.at(i);
	 pw = tpw.at(i);
	 pt = tpt.at(i);

	 
	 //t_rec_charge->Fill();
	 
       }else{
	 charge_save = 0;
	 ncharge_save = -1;
	 pu = -1;
	 pv = -1;
	 pw = -1;
	 pt = -1;
       }
       
     }

    }

   //Trun->CloneTree()->Write();
   
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

  // now save the projected charge information ... 
  TTree *T_proj = new TTree("T_proj","T_proj");
  std::vector<int> *proj_cluster_id = new std::vector<int>;
  std::vector<std::vector<int>> *proj_cluster_channel = new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *proj_cluster_timeslice= new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *proj_cluster_charge= new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *proj_cluster_charge_err= new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *proj_cluster_main_flag = new std::vector<std::vector<int>>;
  T_proj->Branch("cluster_id",&proj_cluster_id);
  T_proj->Branch("channel",&proj_cluster_channel);
  T_proj->Branch("time_slice",&proj_cluster_timeslice);
  T_proj->Branch("charge",&proj_cluster_charge);
  T_proj->Branch("charge_err",&proj_cluster_charge_err);
  T_proj->Branch("flag_main",&proj_cluster_main_flag); 
  T_proj->SetDirectory(file1);
  
  for (auto it = matched_bundles.begin(); it!= matched_bundles.end(); it++){
    FlashTPCBundle *bundle = *it;
    PR3DCluster *main_cluster = bundle->get_main_cluster();
    Opflash *flash = bundle->get_flash();
      
      // now prepare saving it
    int cluster_id = main_cluster->get_cluster_id();
    
    PR3DClusterSelection temp_clusters;
    temp_clusters.push_back(main_cluster);
    for (auto it1 = bundle->get_other_clusters().begin(); it1!=bundle->get_other_clusters().end();it1++){
      temp_clusters.push_back(*it1);
    }
    
    std::vector<int> proj_channel;
    std::vector<int> proj_timeslice;
    std::vector<int> proj_charge;
    std::vector<int> proj_charge_err;
    std::vector<int> proj_flag;
    std::vector<int> proj_flag_main;
    int size_main = 0;
    for (size_t j = 0; j!= temp_clusters.size(); j++){
      PR3DCluster *cluster = temp_clusters.at(j);
      cluster->get_projection(proj_channel,proj_timeslice,proj_charge, proj_charge_err, proj_flag, global_wc_map);
      if (j==0)
          size_main = proj_channel.size();
    }
    proj_flag_main.resize(proj_channel.size(),0);
    for (size_t i=0;i!=size_main;i++){
        proj_flag_main.at(i)=1; 
    }
    proj_cluster_id->push_back(cluster_id);
    proj_cluster_channel->push_back(proj_channel);
    proj_cluster_timeslice->push_back(proj_timeslice);
    proj_cluster_charge->push_back(proj_charge);
    proj_cluster_charge_err->push_back(proj_charge_err);
    proj_cluster_main_flag->push_back(proj_flag_main);
  }
  T_proj->Fill();
  T_proj_data->Fill();


  TFile *file2;
  if(port==1){
    file2 = new TFile(Form("port_old_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");

    TTree *T_port_flash = new TTree("T_port_flash","T_port_flash");
    int f_type;
    double f_time, f_low_time, f_high_time;
    double f_totalPE;
    std::vector<double> f_pe, f_pe_err;
    unsigned int f_run, f_subrun, f_event;
    T_port_flash->Branch("run",&f_run,"run/I");
    T_port_flash->Branch("subrun",&f_subrun,"subrun/I");
    T_port_flash->Branch("event",&f_event,"event/I");
    T_port_flash->Branch("totalPE",&f_totalPE,"totalPE/D");
    T_port_flash->Branch("pe",&f_pe);
    T_port_flash->Branch("pe_err",&f_pe_err);
    T_port_flash->Branch("time",&f_time,"time/D");
    T_port_flash->Branch("low_time",&f_low_time,"low_time/D");
    T_port_flash->Branch("high_time",&f_high_time,"high_time/D");
    T_port_flash->Branch("type",&f_type,"type/I");
    T_port_flash->SetDirectory(file2);

    T_flash->SetBranchAddress("type",&type);    
    T_flash->SetBranchAddress("flash_id",&flash_id);
    T_flash->SetBranchAddress("time",&time);
    T_flash->SetBranchAddress("low_time",&low_time);
    T_flash->SetBranchAddress("high_time",&high_time);
    T_flash->SetBranchAddress("total_PE",&total_PE);
    T_flash->SetBranchAddress("PE",PE);
    T_flash->SetBranchAddress("PE_err",PE_err);
    std::set<int> intime_light_id;
    //if(datatier==0 && triggerbits==2048){ lowerwindow=3.1875; upperwindow=4.96876; } // BNB
    //if(datatier==0 && triggerbits==512){ lowerwindow=3.5625; upperwindow=5.34376; } // BNB EXT
    //if(datatier==2 && triggerbits==2048){ lowerwindow=3.1718; upperwindow=4.95306; } // full mc
    for(int i=0; i<T_flash->GetEntries(); i++){
      T_flash->GetEntry(i);
      if(time>lowerwindow && time<upperwindow){
	intime_light_id.insert(flash_id);
	flash_found=true;
	flash_time=time;

	f_type = type;
	f_time = time;
	f_low_time = low_time;
	f_high_time = high_time;
	f_totalPE = total_PE;
	for(int j=0; j<32; j++){
	  f_pe.push_back(PE[j]);
	  f_pe_err.push_back(PE_err[j]);
	}
	f_run=run_no;
	f_subrun=subrun_no;
	f_event=event_no;
	T_port_flash->Fill();
      }
    }
    
    int tpc_cluster_id;
    T_match->SetBranchAddress("flash_id",&flash_id);
    T_match->SetBranchAddress("tpc_cluster_id",&tpc_cluster_id);
    std::set<int> assoc_cluster;
    int selected_cluster;
    for(int i=0; i<T_match->GetEntries(); i++){
      T_match->GetEntry(i);
      for(std::set<int>::iterator it=intime_light_id.begin(); it!=intime_light_id.end(); it++){
	if(flash_id != *it) continue;
	assoc_cluster.insert(tpc_cluster_id);
	selected_cluster = tpc_cluster_id;
      }
    }

    TTree *T_port_3d = new TTree("T_port_3d","T_port_3d");
    int p_cluster_id, p_main_flag, p_time_slice, p_ch_u, p_ch_v, p_ch_w;
    double p_x, p_y, p_z, p_q, p_nq;
    unsigned int p_run, p_subrun, p_event;
    T_port_3d->Branch("run",&p_run,"run/I");
    T_port_3d->Branch("subrun",&p_subrun,"subrun/I");
    T_port_3d->Branch("event",&p_event,"event/I");
    T_port_3d->Branch("cluster_id",&p_cluster_id,"cluster_id/I");
    T_port_3d->Branch("main_flag",&p_main_flag,"main_flag/I");
    T_port_3d->Branch("time_slice",&p_time_slice,"time_slice/I");
    T_port_3d->Branch("ch_u",&p_ch_u,"ch_u/I");
    T_port_3d->Branch("ch_v",&p_ch_v,"ch_v/I");
    T_port_3d->Branch("ch_w",&p_ch_w,"ch_w/I");
    T_port_3d->Branch("x",&p_x,"x/D");
    T_port_3d->Branch("y",&p_y,"y/D");
    T_port_3d->Branch("z",&p_z,"z/D");
    T_port_3d->Branch("q",&p_q,"q/D");
    T_port_3d->Branch("nq",&p_nq,"nq/D");
    T_port_3d->SetDirectory(file2);
    
    int c_id,c_ch_u,c_ch_v,c_ch_w;
    double c_x,c_y,c_z,c_q,c_nq;
    int c_flag_main,c_time_slice;
    std::vector<double> intime_charge_x, intime_charge_y, intime_charge_z;
    intime_charge_x.clear();
    intime_charge_y.clear();
    intime_charge_z.clear();
    T_cluster->SetBranchAddress("cluster_id",&c_id);
    T_cluster->SetBranchAddress("time_slice",&c_time_slice);
    T_cluster->SetBranchAddress("ch_u",&c_ch_u);
    T_cluster->SetBranchAddress("ch_v",&c_ch_v);
    T_cluster->SetBranchAddress("ch_w",&c_ch_w);
    T_cluster->SetBranchAddress("q",&c_q);
    T_cluster->SetBranchAddress("nq",&c_nq);
    T_cluster->SetBranchAddress("x",&c_x);
    T_cluster->SetBranchAddress("y",&c_y);
    T_cluster->SetBranchAddress("z",&c_z);
    T_cluster->SetBranchAddress("flag_main",&c_flag_main);
    for(int i=0; i<T_cluster->GetEntries(); i++){
      T_cluster->GetEntry(i);
      for(std::set<int>::iterator it=assoc_cluster.begin(); it!=assoc_cluster.end(); it++){
	if(*it == c_id && c_flag_main==1) //main cluster
        {
	  intime_charge_x.push_back(c_x);
	  intime_charge_y.push_back(c_y);
	  intime_charge_z.push_back(c_z);
	}
	if(*it == c_id)
	{
	p_run=run_no;
	p_subrun=subrun_no;
	p_event=event_no;
	p_cluster_id = c_id;
	p_main_flag = c_flag_main;
	p_time_slice = c_time_slice;
	p_ch_u = c_ch_u;
	p_ch_v = c_ch_v;
	p_ch_w = c_ch_w;
	p_x = c_x;
	p_y = c_y;
	p_z = c_z;
	p_q = c_q;
	p_nq = c_nq;
	T_port_3d->Fill();
	}
      }
    }
    if(intime_charge_x.size()!=0 && nuFinalState_x.size()!=0){
      double minX_element = *std::min_element(intime_charge_x.begin(),intime_charge_x.end());
      double minY_element = *std::min_element(intime_charge_y.begin(),intime_charge_y.end());
      double minZ_element = *std::min_element(intime_charge_z.begin(),intime_charge_z.end());
      double maxX_element = *std::max_element(intime_charge_x.begin(),intime_charge_x.end());
      double maxY_element = *std::max_element(intime_charge_y.begin(),intime_charge_y.end());
      double maxZ_element = *std::max_element(intime_charge_z.begin(),intime_charge_z.end());

	double true_minX_element = *std::min_element(nuFinalState_x.begin(), nuFinalState_x.end());
	double true_minY_element = *std::min_element(nuFinalState_y.begin(), nuFinalState_y.end());
	double true_minZ_element = *std::min_element(nuFinalState_z.begin(), nuFinalState_z.end());
	double true_maxX_element = *std::max_element(nuFinalState_x.begin(), nuFinalState_x.end());
	double true_maxY_element = *std::max_element(nuFinalState_y.begin(), nuFinalState_y.end());
	double true_maxZ_element = *std::max_element(nuFinalState_z.begin(), nuFinalState_z.end());

    /* 3D histogram voxelization */
    // Reco and True --> match completeness
    // Reco and BlurTrue --> match purity 
    double minX = true_minX_element<minX_element ? true_minX_element:minX_element; 
    double minY = true_minY_element<minY_element ? true_minY_element:minY_element; 
    double minZ = true_minZ_element<minZ_element ? true_minZ_element:minZ_element; 
    double maxX = true_maxX_element>maxX_element ? true_maxX_element:maxX_element; 
    double maxY = true_maxY_element>maxY_element ? true_maxY_element:maxY_element; 
    double maxZ = true_maxZ_element>maxZ_element ? true_maxZ_element:maxZ_element; 
    
    TH3D *h3RecoCluster = new TH3D("h3RecoCluster", "In-time matched recostructed cluster", 
				    ((int)maxX)-((int)minX)+6, ((int)minX)-3, ((int)maxX)+3,
				    ((int)maxY)-((int)minY)+6, ((int)minY)-3, ((int)maxY)+3,
				    ((int)maxZ)-((int)minZ)+6, ((int)minZ)-3, ((int)maxZ)+3);
    TH3D *h3TrueCluster = new TH3D("h3TrueCluster", "In-time true neutrino cluster", 
				    ((int)maxX)-((int)minX)+6, ((int)minX)-3, ((int)maxX)+3,
				    ((int)maxY)-((int)minY)+6, ((int)minY)-3, ((int)maxY)+3,
				    ((int)maxZ)-((int)minZ)+6, ((int)minZ)-3, ((int)maxZ)+3);
    TH3D *h3TrueBlurCluster = new TH3D("h3TrueBlurCluster", "In-time true neutrino cluster blurred", 
				    ((int)maxX)-((int)minX)+6, ((int)minX)-3, ((int)maxX)+3,
				    ((int)maxY)-((int)minY)+6, ((int)minY)-3, ((int)maxY)+3,
				    ((int)maxZ)-((int)minZ)+6, ((int)minZ)-3, ((int)maxZ)+3);
    h3RecoCluster->SetDirectory(file1);
    h3TrueCluster->SetDirectory(file1);
    h3TrueBlurCluster->SetDirectory(file1);
    

    for(size_t i=0; i<intime_charge_x.size(); i++){
	    h3RecoCluster->Fill( intime_charge_x.at(i), intime_charge_y.at(i), intime_charge_z.at(i), 1.0);
    }  	
	for(size_t i=0; i<nuFinalState_x.size(); i++){
        // truth x position considering offset in imaging
        float corrected_x =  ((nuFinalState_x.at(i) + 0.6 /*6 mm offset from U plane to Y plane*/)*unit_dis/1.098 /*drift velocity in sim*/
                                              + (truth_nuTime+0.0+0.0)*unit_dis*0.1 /*unit: cm, t0 + decon residual + rebin*/);
	    h3TrueCluster->Fill( corrected_x, nuFinalState_y.at(i), nuFinalState_z.at(i), nuFinalState_energy.at(i) ); // also for calculating true energy match fraction
        for(int indx = -1; indx<=1; indx++)
        {
            for(int indy = -1; indy<=1; indy++)
            {
                for(int indz = -1; indz<=1; indz++)
                {
	                h3TrueBlurCluster->Fill( corrected_x+indx, nuFinalState_y.at(i)+indy, nuFinalState_z.at(i)+indz, 1./27);
                }
            }
        }
	}
    // calculate match completeness, completeness energy, purity 
    float match_reco_elements = 0;
    float match_true_elements = 0;
    for(int xind = 0; xind < h3RecoCluster->GetNbinsX(); xind++)
    {
        for(int yind = 0; yind < h3RecoCluster->GetNbinsY(); yind++)
        {
            for(int zind = 0; zind < h3RecoCluster->GetNbinsZ(); zind++)
            {
                float recoc = h3RecoCluster->GetBinContent(xind, yind, zind);          
                float truec = h3TrueCluster->GetBinContent(xind, yind, zind);          
                float trueblurc = h3TrueBlurCluster->GetBinContent(xind, yind, zind);
                if(recoc!=0) match_reco_elements += 1.0;
                if(truec!=0) match_true_elements += 1.0;
                if(recoc!=0 && truec!=0) { match_completeness += 1.0; match_completeness_energy += truec; }
                if(recoc!=0 && trueblurc!=0) { match_purity += 1.0; }
            }
        }
    }
    match_completeness = match_completeness/match_true_elements; 
    match_purity = match_purity/match_reco_elements; 
    
    // 2D purity 
    TH2D* hrecoxz = (TH2D*)h3RecoCluster->Project3D("xz");
    TH2D* hrecoxy = (TH2D*)h3RecoCluster->Project3D("xy");
    TH2D* htruexz = (TH2D*)h3TrueBlurCluster->Project3D("xz");
    TH2D* htruexy = (TH2D*)h3TrueBlurCluster->Project3D("xy");

    match_reco_elements=0;
    for(int xind = 0; xind < hrecoxz->GetNbinsX(); xind++)
    {
	    for(int yind = 0; yind < hrecoxz->GetNbinsY(); yind++)
	    {
		    float recoc = hrecoxz->GetBinContent(xind, yind);
		    float truec = htruexz->GetBinContent(xind, yind);
		    if(recoc!=0) match_reco_elements += 1.0;
		    if(recoc!=0 && truec!=0) { match_purity_xz += 1.0; }
	    }
    }
    match_purity_xz = match_purity_xz/match_reco_elements;

    match_reco_elements = 0;
    for(int xind = 0; xind < hrecoxy->GetNbinsX(); xind++)
    {
	    for(int yind = 0; yind < hrecoxy->GetNbinsY(); yind++)
	    {
		    float recoc = hrecoxy->GetBinContent(xind, yind);
		    float truec = htruexy->GetBinContent(xind, yind);
		    if(recoc!=0) match_reco_elements += 1.0;
		    if(recoc!=0 && truec!=0) { match_purity_xy += 1.0; }
	    }
    }
    match_purity_xy = match_purity_xy/match_reco_elements;


    } // intime_charge.size != 0

    TTree *T_port_2d = new TTree("T_port_2d","T_port_2d");
    unsigned int run, subrun, event;
    int channel, start_tick, main_flag;
    float charge, charge_error;
    T_port_2d->Branch("run",&run,"run/I");
    T_port_2d->Branch("subrun",&subrun,"subrun/I");
    T_port_2d->Branch("event",&event,"event/I");
    T_port_2d->Branch("channel",&channel,"channel/I");
    T_port_2d->Branch("start_tick",&start_tick,"start_tick/I");
    T_port_2d->Branch("main_flag",&main_flag,"main_flag/I");
    T_port_2d->Branch("charge",&charge,"charge/F");
    T_port_2d->Branch("charge_error",&charge_error,"charge_error/F");
    T_port_2d->SetDirectory(file2);

    std::set<int> intime_charge_id;
    for(size_t i=0; i!=proj_cluster_id->size(); i++){
      for(std::set<int>::iterator it=assoc_cluster.begin(); it!=assoc_cluster.end(); it++){
	if(proj_cluster_id->at(i)==*it){ intime_charge_id.insert(*it); }
      }
    } 
    
    for(std::set<int>::iterator it=intime_charge_id.begin(); it!=intime_charge_id.end(); it++){
      for(size_t a=0; a!=proj_cluster_id->size(); a++){
	if(proj_cluster_id->at(a)!=*it) continue;
	for(size_t i=0; i!=proj_cluster_channel->at(a).size(); i++){
      // select "main cluster"
	  if(proj_cluster_charge->at(a).at(i)*0.005<1.0) continue;
	  run=run_no;
	  subrun=subrun_no;
	  event=event_no;
	  channel = proj_cluster_channel->at(a).at(i);
	  start_tick = (int)4*proj_cluster_timeslice->at(a).at(i);
	  main_flag = proj_cluster_main_flag->at(a).at(i);
	  charge = proj_cluster_charge->at(a).at(i);
	  charge_error = proj_cluster_charge_err->at(a).at(i);
	  T_port_2d->Fill();
      // only collection plane
	  if(channel>=4800 && main_flag==1) match_charge += charge;
	} // i
      } // a
    } // it
    //file2->Write();
    //file2->Close();
  } 
  cout << em("matching evaluation finished") << endl;
  
  //match_charge /= 3.0; all planes
  match_energy = match_charge / 0.826 * 23.6 * 1e-6 / 0.55; // [e-]*[23.6 eV / e-]*[10^-6 MeV / eV]*[1/recombination]
  cout<< match_charge << endl;
  
  T_eval->Fill();

  TTree *T_port_eval = T_eval->CloneTree();
  T_port_eval->SetDirectory(file2);
  file2->Write();
  file2->Close();
    
  file1->Write();
  file1->Close();
}
