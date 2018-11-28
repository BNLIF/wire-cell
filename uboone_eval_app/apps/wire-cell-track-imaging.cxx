#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPolyLine3D.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TMath.h"
using namespace std;

#include <iostream>
#include <map>

#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVector3.h"

#include "TStyle.h"
#include "TColor.h"
#include "TROOT.h"

#include "TString.h"
#include "TGraph2D.h"
#include "TPolyLine3D.h"
#include "TGraph.h"


/*
  short cluster: satisfy "distance"
  Flag: multiple cluters to one track
*/


void set_plot_style()
{
  TStyle* myStyle = new TStyle("myStyle","My ROOT plot style");
  // plot style
  //myStyle->SetPalette(kInvertedDarkBodyRadiator); 
  //myStyle->SetPalette(kRainBow);
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  myStyle->SetNumberContours(NCont);
  myStyle->SetOptStat(0);
    
  myStyle->SetLabelFont(62,"xyz");
  myStyle->SetLabelSize(0.05,"xyz");
  myStyle->SetTitleFont(62, "xyz");
  myStyle->SetTitleSize(0.06,"xyz");
  myStyle->SetTitleOffset(1.3, "y");
  myStyle->SetTitleOffset(0.8, "x");
  myStyle->SetTitleOffset(1.5, "z");
    
  // only 5 in x to avoid label overlaps
  myStyle->SetNdivisions(505, "x");

  // set the margin sizes
  myStyle->SetPadTopMargin(0.05);
  myStyle->SetPadRightMargin(0.2); // increase for colz plots
  myStyle->SetPadBottomMargin(0.1);
  myStyle->SetPadLeftMargin(0.15);
  myStyle->SetFrameBorderMode(0);
  myStyle->SetCanvasBorderMode(0);
  myStyle->SetPadBorderMode(0);
  myStyle->SetCanvasColor(0);
  myStyle->SetPadColor(0);

  gROOT->SetStyle("myStyle");
  gROOT->ForceStyle();
}

int main(int argc, char* argv[])
{
  if(argc < 3){
    cout<<"Usage: wire-cell-track-eval input.root outputname.root [option: cluster_id]"<<endl;
    return 1;
  }
  const char* inputroot = argv[1];
  const char* outputroot = argv[2];
  int cluster_check=0; // track classification to check + clusters
  //if(argc==4) cluster_check=atoi(argv[3]);
    
  for(int i=1; i!=argc; i++){
    switch(argv[i][1]){
    case 'c': {
      cluster_check = atoi(&argv[i][2]); // 1
      break;
    }
    default: break; // 0
    }
  }

  //cout<<" ---> test track-imaging A"<<endl;
  //cout<<" ---> test track-imaging B"<<endl;

  TString roostr = "";

  //////////////////////////////////////////////////////////////////////////////////////////////// variable histgram start

  roostr = "graph_ghost_yz";
  TGraph *graph_ghost_yz = new TGraph();
  graph_ghost_yz->SetName(roostr);

  std::map< int, std::vector<int> > map_non_ghost_track_clusters;
  std::map<int,int> map_flag_cluster_used;

  TGraph2D *graph_track_each[21];
  for(int idx=0; idx<21; idx++) {
    roostr = TString::Format("graph_track_each_%2d", idx);
    graph_track_each[idx] = new TGraph2D();
    graph_track_each[idx]->SetName(roostr);
  }
  int count_graph_track_each[21] = {0};

  TGraph2D *graph_cluster_all = new TGraph2D();
  graph_cluster_all->SetName("graph_cluster_all");

  std::map< int, TVector3 > map_cluster_center;

  //////////////////////////////////////////////////////////////////////////////////////////////// variable histgram end
  
  // rotation matrix to convert collection plane (default detector) coordinate to inductions plane
  TMatrixD Muplane(3,3);
  TMatrixD Mvplane(3,3);
  double angle = 60./180*TMath::Pi();
  double cosine = TMath::Cos(angle);
  double sine = TMath::Sin(angle);
  Muplane(0,0) = 1.;
  Muplane(0,1) = 0;
  Muplane(0,2) = 0;
  Muplane(1,0) = 0;
  Muplane(1,1) = cosine;
  Muplane(1,2) = sine;
  Muplane(2,0) = 0;
  Muplane(2,1) = -sine;
  Muplane(2,2) = cosine;

  Mvplane(0,0) = 1.;
  Mvplane(0,1) = 0;
  Mvplane(0,2) = 0;
  Mvplane(1,0) = 0;
  Mvplane(1,1) = cosine;
  Mvplane(1,2) = -sine;
  Mvplane(2,0) = 0;
  Mvplane(2,1) = sine;
  Mvplane(2,2) = cosine;

  // read in the point cloud for each cluster_id
  TFile* f = new TFile(inputroot, "READ");

  ///////////////////////////////////////////////////////////

  TTree* track_true = (TTree*)f->Get("T_track");
  double x0,y0,z0;
  double x1,y1,z1;
  int track_id;
  track_true->SetBranchAddress("x0",&x0);
  track_true->SetBranchAddress("y0",&y0);
  track_true->SetBranchAddress("z0",&z0);
  track_true->SetBranchAddress("x1",&x1);
  track_true->SetBranchAddress("y1",&y1);
  track_true->SetBranchAddress("z1",&z1);  
  track_true->SetBranchAddress("track_id",&track_id);
  std::map<int, TVector3> vc_track_dir;
  std::map<int, TVector3> vc_track_start;
  std::map<int, TVector3> vc_track_end;
  std::map<int, double> vc_track_length;

  vector<int> vc_track_id;
  
  int num_track_clusters = track_true->GetEntries();
  for(int ientry=0; ientry<num_track_clusters; ientry++) {
    track_true->GetEntry( ientry );
    
    TVector3 v3_start(x0, y0, z0);
    TVector3 v3_end(x1, y1, z1);

    TVector3 dir = v3_end - v3_start;
    vc_track_length[track_id] = dir.Mag();
    dir = 1./dir.Mag() * dir;
    if( dir.Y()<0 ) dir *= -1.;
    vc_track_dir[track_id] = dir;
    
    if( v3_start.Y() < v3_end.Y() ) {
      vc_track_start[track_id] = v3_start;
      vc_track_end[track_id] = v3_end;
    }
    else {
      vc_track_start[track_id] = v3_end;
      vc_track_end[track_id] = v3_start;
    }

    vc_track_id.push_back( track_id );
  }

  
  for(int i=0; i<vc_track_id.size(); i++) {
    for(int j=i+1; j<vc_track_id.size(); j++) {

      int this_id = vc_track_id.at( i );
      int next_id = vc_track_id.at( j );

      double step = 0.1;// cm
      int num_cycle =  vc_track_length[this_id]/step;
      
      for(int idx=0; idx<=num_cycle; idx++) {
	TVector3 this_point = vc_track_start[this_id] + step*idx*vc_track_dir[this_id];
	TVector3 this_point_2_next_start = vc_track_start[next_id] - this_point;
        double distance_vertical = ( this_point_2_next_start.Cross( vc_track_dir[next_id] ) ).Mag();
	double distance_2start = ( this_point-vc_track_start[next_id] ).Mag();
	double distance_2end = ( this_point-vc_track_end[next_id] ).Mag();
	double distance = 10000;
	
	/*         

		   - test_A
		   -
		   -
		   <---------
		   start   end



		   - test_B
		   -
		   -
		   ---------->
		   start     end

	*/

	double angle_end = (this_point-vc_track_end[next_id]).Dot( -1.*vc_track_dir[next_id] );
	double angle_start = (this_point-vc_track_start[next_id]).Dot( vc_track_dir[next_id] );
	if( angle_end * angle_start<0 ) {
	  distance = (distance_2start<distance_2end) ? distance_2start : distance_2end;
	}
	else {
	  distance = distance_vertical;
	}
	
	if( distance < 2 ) {// cm
	  cout<<endl<<" WARNING : track crossing: distance between two tracks is very close ( < 2 cm )"<<endl<<endl;
	  cout<<this_id<<"\t"<<next_id<<"\t"<<distance<<endl;
	  exit(0);
	}
	
      }
    }
  }

  //exit(0);
  
  //////////////////////////////////////////////////////////
  
  TTree* clusters = (TTree*)f->Get("T_cluster");
  Double_t x=0;
  Double_t y=0;
  Double_t z=0;
  Double_t q=0;
  Int_t cluster_id=0;
  clusters->SetBranchAddress("x",&x);
  clusters->SetBranchAddress("y",&y);
  clusters->SetBranchAddress("z",&z);
  clusters->SetBranchAddress("q",&q);
  clusters->SetBranchAddress("cluster_id",&cluster_id);

  std::map<int, std::vector<double>> xpt;
  std::map<int, std::vector<double>> ypt;
  std::map<int, std::vector<double>> zpt;
  std::map<int, std::vector<double>> qpt;

  ///////////////////////////////////////////////////////////
  
  int track_id_debug = 0;
  TTree *truth = (TTree*)f->Get("T_true");
  truth->SetBranchAddress("x",&x);
  truth->SetBranchAddress("y",&y);
  truth->SetBranchAddress("z",&z);
  truth->SetBranchAddress("track_id",&track_id_debug);
  
  TGraph2D *gh_truth = new TGraph2D();
  gh_truth->SetName("gh_truth");
  
  for(long ientry=0; ientry<truth->GetEntries(); ientry++) {
    truth->GetEntry(ientry);
    gh_truth->SetPoint( ientry, x,y,z );

    //////////
    count_graph_track_each[track_id_debug]++;
    graph_track_each[track_id_debug]->SetPoint( count_graph_track_each[track_id_debug]-1, x,y,z );    
  }
  
  
  ///////
  std::map<int, TVector3> vc_cluster_dir;
  std::map<int, TVector3> vc_cluster_start;
  std::map<int, TVector3> vc_cluster_end;
  std::map<int, double> vc_cluster_length;
  
  TGraph2D *gh_ghost_all = new TGraph2D();
  gh_ghost_all->SetName("gh_ghost_all");
  int count_gh_ghost_all = 0;
  
  cout<<endl;
  cout<<" -------> Points in cluster, entries: "<<clusters->GetEntries()<<endl;
  cout<<endl;
  
  for(long i=0; i<clusters->GetEntries();i++) {
    clusters->GetEntry(i);
    auto it = xpt.find(cluster_id);
    if( it == xpt.end() ){
      std::vector<double> xvec;
      xvec.push_back(x);
      std::vector<double> yvec;
      yvec.push_back(y);
      std::vector<double> zvec;
      zvec.push_back(z);
      std::vector<double> qvec;
      qvec.push_back(q);
	
      xpt[cluster_id] = xvec;
      ypt[cluster_id] = yvec;
      zpt[cluster_id] = zvec;
      qpt[cluster_id] = qvec;
    }
    else{
      xpt[cluster_id].push_back(x);
      ypt[cluster_id].push_back(y);
      zpt[cluster_id].push_back(z);
      qpt[cluster_id].push_back(q);
    }

    graph_cluster_all->SetPoint(i, x,y,z);
  }

  // PCA (Priciple Component Analysis): Matrix Decomposition --> Coordinate Rotation
  // n-dimentional points in nxn position covariance matrix;
  // diagonalization (rotation);
  // the biggest eigen value corresponds to the main axis --> primary direction of the point cloud
    
  // Step 1: calculate average position
  // Step 2: construct covariance matrix
  // merge the two steps into a single loop of the point cloud
  
  for( auto it=xpt.begin(); it!=xpt.end(); it++ ) { // each cluster
    int i = it->first;
  
    int npts = xpt[i].size();
    double xx=0;
    double yy=0;
    double zz=0;
    double xy=0;
    double yz=0;
    double zx=0;
    double mx=0;
    double my=0;
    double mz=0;
    for(int j=0; j<npts; j++) { // each point      
      double x = xpt[i].at(j);
      double y = ypt[i].at(j);
      double z = zpt[i].at(j);

      xx += x*x/npts;
      yy += y*y/npts;
      zz += z*z/npts;
      xy += x*y/npts;
      yz += y*z/npts;
      zx += z*x/npts;
      mx += x/npts;
      my += y/npts;
      mz += z/npts;

      //polyline[i]->SetPoint(j, x,y,z);
    } //each point   
    
    TMatrixDSym m(3);
    m(0,0) = xx - mx*mx; 
    m(0,1) = xy - mx*my;
    m(0,2) = zx - mz*mx;
    m(1,0) = xy - mx*my;
    m(1,1) = yy - my*my;
    m(1,2) = yz - my*mz;
    m(2,0) = zx - mz*mx;
    m(2,1) = yz - my*mz;
    m(2,2) = zz - mz*mz;

    TMatrixDSymEigen me(m);
    auto& eigenval = me.GetEigenValues();
    auto& eigenvec = me.GetEigenVectors();

    double maxeval =  eigenval(0);
    int maxevalaxis = 0;

    for(int k=1; k<3; k++) {
      if(eigenval(k)>maxeval)
	{
	  maxevalaxis = k;
	  maxeval = eigenval(k);
	}
    }
        
    // PCA main direction
    double pcadir_x = eigenvec(0, maxevalaxis);
    double pcadir_y = eigenvec(1, maxevalaxis);
    double pcadir_z = eigenvec(2, maxevalaxis);

    // track direction
    TVector3 pca_dir(pcadir_x, pcadir_y, pcadir_z);
    TVector3 track_dir;
    track_dir = (1./pca_dir.Mag())*pca_dir;
      
    // theta_y: 0-90 degree, upgoing
    if(track_dir.Y()<0) track_dir *= -1.0;
      
    // Loop over point cloud to find the edges (start, end) along PCA main direction        
    TVector3 start(xpt[i].at(0), ypt[i].at(0), zpt[i].at(0));
    TVector3 end(start);
    double start_proj = track_dir.Dot(start);
    double end_proj = track_dir.Dot(end);

    for(int j=1; j<npts; j++) {//each point
      TVector3 point(xpt[i].at(j), ypt[i].at(j), zpt[i].at(j));
      double point_proj = track_dir.Dot(point);
      if(point_proj < start_proj)
	{
	  start=point;
	  start_proj = point_proj;
	}
      if(point_proj > end_proj)
	{
	  end=point;
	  end_proj = point_proj;
	} 
    }//each point

    // the start, end points are more precise than the PCA direction in case of a big blob (isochronous track)

    ///////
    TVector3 cluster_dir = end - start;// direction is alway "Up" in Y-axis
    vc_cluster_dir[i] = cluster_dir;
    vc_cluster_dir[i] = 1./vc_cluster_dir[i].Mag() * vc_cluster_dir[i];
    if( vc_cluster_dir[i].Y()<0 ) vc_cluster_dir[i] *= -1.;
    
    if( start.Y()<=end.Y() ) {
      vc_cluster_start[i] = start;
      vc_cluster_end[i] = end;
    }
    else {
      vc_cluster_start[i] = end;
      vc_cluster_end[i] = start;
    }
    vc_cluster_length[i] = cluster_dir.Mag();

    ///////
    TVector3 cluster_dir_y = cluster_dir;
    if(cluster_dir_y.Y()<0) cluster_dir_y *= -1.0;
    
    TVector3 cluster_dir_u = Muplane*cluster_dir;
    if(cluster_dir_u.Y()<0) cluster_dir_u *= -1.0;

    TVector3 cluster_dir_v = Mvplane*cluster_dir;
    if(cluster_dir_v.Y()<0) cluster_dir_v *= -1.0;


    ////////////////////////////////////////////

    TVector3 cluster_center( mx, my, mz );
    
    map_cluster_center[ i ] = cluster_center;
   
  }

  ///////////////////////////////////////////////////////////////
  ///////////// clusters in one track, and ghost ////////////////
  ///////////////////////////////////////////////////////////////
  
  cout<<endl<<" -------> processing clusters in one track, and ghost"<<endl<<endl;

  int num_clusters = vc_cluster_dir.size();
  cout<<" -------> Cluster numbers: "<<num_clusters<<endl;
  cout<<" -------> Track numbers:   "<<num_track_clusters<<endl;
  cout<<endl;

  /*
    cell_unit = 3 mm
    ---> +/- uncertainty = 6 mm
    
    arc_length = radius * angle(in PI unit)

    delta_arc_length = sqrt( 6*6 + 6*6 ) = 8.5 mm ---> 0.85 cm

    tolerance = delta_arc_length/radius

    sin(0.05) = 0.04998
    sin(0.10) = 0.09983
    sin(0.20) = 0.19867
    sin(0.50) = 0.47943
  */

  const double track_length_setting = 1;// cm

  const double delta_arc_length = 0.85;
  double tolerance_parallel = delta_arc_length/vc_track_length[track_id];      
  tolerance_parallel = 0.17365;// ---> 10 degree

  double angle_15 = 0.25882;// sin(15 degree) = 0.2588;
  double angle_10 = 0.17365;
  double tolerance_distance2track = 2; // cm
  double tolerance_contain = 2;// cm
  double short_cluster = 5;
  
  std::map<int, int> vc_ghost_id;
  std::map<int, int> vc_ghost_angle;
  std::map<int, int> vc_ghost_distance;
  std::map<int, int> vc_ghost_contain;

  for( auto it_track=vc_track_dir.begin(); it_track!=vc_track_dir.end(); it_track++ ) {
    map_non_ghost_track_clusters[it_track->first].push_back(-1);//
  }
  
  for( auto it_cluster=vc_cluster_dir.begin(); it_cluster!=vc_cluster_dir.end(); it_cluster++ ) {// By default, all the cluster are ghosts
    vc_ghost_id[it_cluster->first] = it_cluster->first;
  }
  

  for( auto it_cluster=vc_cluster_dir.begin(); it_cluster!=vc_cluster_dir.end(); it_cluster++ ) {// each cluster cluster
    for( auto it_track=vc_track_dir.begin(); it_track!=vc_track_dir.end(); it_track++ ) {// each track cluster
      
      int track_id = it_track->first;
      int cluster_id = it_cluster->first;
      
      // if( it_cluster->first!=9 && it_cluster->first!=26 ) continue;
      // int id = cluster_id;
      // if( id!=20 ) continue;
      // if( id!=15 && id!=19 && id!=22 && id!=26 ) continue;
      
      /// angle
      TVector3 v3_trackXcluster = vc_track_dir[track_id].Cross( vc_cluster_dir[cluster_id] );

      /// distance
      TVector3 v3_cluster_start_2_track_end = vc_track_end[track_id] - vc_cluster_start[cluster_id];
      TVector3 v3_startXend = v3_cluster_start_2_track_end.Cross( vc_track_dir[track_id] );
      double distance_startXend = v3_startXend.Mag();
      TVector3 v3_cluster_end_2_track_start = vc_track_start[track_id] - vc_cluster_end[cluster_id];
      TVector3 v3_endXstart = v3_cluster_end_2_track_start.Cross( vc_track_dir[track_id] );
      double distance_endXstart = v3_endXstart.Mag();
      //double distance_center = 0.5*(distance_startXend + distance_endXstart);      
      TVector3 v3_cluster_center_2_track_start = vc_track_start[track_id] - map_cluster_center[cluster_id];
      TVector3 v3_centerXstart = v3_cluster_center_2_track_start.Cross( vc_track_dir[track_id] );
      double distance_center = v3_centerXstart.Mag();
      
      /// contain
      double value_track_start = vc_track_start[track_id].Dot( vc_track_dir[track_id] );//A
      double value_track_end = vc_track_end[track_id].Dot( vc_track_dir[track_id] );//B
      double value_cluster_start = vc_cluster_start[cluster_id].Dot( vc_track_dir[track_id] );//C
      double value_cluster_end = vc_cluster_end[cluster_id].Dot( vc_track_dir[track_id] );//D
      double min = (value_track_start>value_track_end)?(value_track_end-tolerance_contain):(value_track_start-tolerance_contain);
      double max = (value_track_start>value_track_end)?(value_track_start+tolerance_contain):(value_track_end+tolerance_contain);

      /* cout<<TString::Format(" cluster %2d , track %2d, angleT %5.3f, angle %5.3f, distance center %8.2f, contain min/max %8.2f %8.2f, start/end, %8.2f %8.2f", */
      /* 			    cluster_id, track_id, */
      /* 			    tolerance_parallel, v3_trackXcluster.Mag(), */
      /* 			    distance_center, */
      /* 			    min, max, value_cluster_start, value_cluster_end )<<endl; */


      bool flag_contain = false;
      if( value_cluster_start>min && value_cluster_start<max && value_cluster_end>min && value_cluster_end<max ) flag_contain = true;

      bool flag_length = false;
      if( vc_cluster_length[cluster_id]<short_cluster ) flag_length = true;

      bool flag_angle_LE_10 = false;
      if( v3_trackXcluster.Mag()<angle_10 ) flag_angle_LE_10 = true;

      bool flag_angle_LE_15 = false;
      if( v3_trackXcluster.Mag()<angle_15 ) flag_angle_LE_15 = true;
      
      bool flag_distance_2_3cm = false;
      if( distance_center< 3 && distance_center>2 ) flag_distance_2_3cm = true;

      bool flag_distance_0_2cm = false;
      if( distance_center<= 2 ) flag_distance_0_2cm = true;

      bool case_A = flag_contain && flag_length && ( distance_center<2 );
      bool case_B = flag_contain && flag_angle_LE_10 && flag_distance_2_3cm;
      bool case_C = flag_contain && flag_angle_LE_15 && flag_distance_0_2cm;
      
      if( case_A || case_B || case_C ) {
      
	// std::map<int, TVector3> vc_cluster_dir;
	// std::map<int, TVector3> vc_cluster_start;
	// std::map<int, TVector3> vc_cluster_end;
	// std::map<int, double> vc_cluster_length;
 
	// std::map<int, TVector3> vc_track_dir;
	// std::map<int, TVector3> vc_track_start;
	// std::map<int, TVector3> vc_track_end;
	// std::map<int, double> vc_track_length;
	
	// not ghost
	vc_ghost_id.erase( cluster_id );
	  
	if( map_flag_cluster_used.find( cluster_id )!=map_flag_cluster_used.end() ) {
	  map_flag_cluster_used[cluster_id] += 1;
	  continue;
	}
	    
	map_flag_cluster_used[cluster_id] = 1;
	map_non_ghost_track_clusters[track_id].push_back( cluster_id );
	
      }// if( case_A || case_B || case_C )
    }// each cluster cluster
  }// each track cluster


  map<int, int>flag_track2recon;// 0: inefficiency, 1: broken, 2: good
  for( auto it=map_non_ghost_track_clusters.begin(); it!=map_non_ghost_track_clusters.end(); it++ ) {
    int user_size = it->second.size() - 1;

    //cout<<it->first<<"\t"<<user_size<<endl;
      
    if( user_size==0 ) flag_track2recon[it->first] = 0;
    if( user_size>1  ) {

      flag_track2recon[it->first] = 1;
      
      for(int idx=1; idx<=user_size; idx++) {
	int cluster_id = it->second.at(idx);
	double cluster_length = vc_cluster_length[cluster_id];
	if( cluster_length>track_length_setting*0.95 ) {
	  flag_track2recon[it->first] = 2;
	  break;
	}
      }// for(int idx=1; idx<=user_size; idx++)
    }// if( user_size>1  )
    
    if( user_size==1 ) {
      flag_track2recon[it->first] = 2;

      /// get cluster_id
      int cluster_id = it->second.at(1);
      double cluster_length = vc_cluster_length[cluster_id];

      // if( cluster_length<5 ) {// ghost ---> inefficinecy
      // 	//flag_track2recon[it->first] = 0;
      // 	//vc_ghost_id[cluster_id] = cluster_id;// pay attention
      // 	bool flag = false;
      // }
      
      if( cluster_length<track_length_setting*0.95 ) {// broken
	flag_track2recon[it->first] = 1;
      }
      
    }
    
  }
  
  for( auto it=map_flag_cluster_used.begin(); it!=map_flag_cluster_used.end(); it++ ) {
    if(it->second > 1) {
      cout<<" -----------> WARNING : cluster "<<it->first<<" used "<<it->second<< " times." <<endl;
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////

  double x_sigma_cut = 1;// cm
  
  for( auto it=vc_cluster_dir.begin(); it!=vc_cluster_dir.end(); it++ ) {
    
    if( vc_ghost_id.find( it->first )==vc_ghost_id.end() ) continue;
    
    int id = it->first;
    int cluster_id = id;
    
    if( vc_cluster_length[id] > 10 ) {// cm      
      int npts = xpt[id].size();
      double x2_bar = 0;
      double x_bar = 0;
      for(int i=0; i<npts; i++) {
	double x =  xpt[id].at(i);
	x2_bar += x*x/npts;
	x_bar += x/npts;
      }
      
      double var_x =  x2_bar - x_bar*x_bar;
      if(var_x<=0) var_x = 0;
      double sigma_x = sqrt( var_x );// standared deviation
      //cout<<"DEBUG: x2_bar "<<x2_bar<<" x_bar*x_bar "<<x_bar*x_bar<<" difference "<<x2_bar - x_bar*x_bar<<endl;
      cout<<" WARNING : big blob "<<id<<"\t"<<vc_cluster_length[id]<<"\t"<<npts<<"\t"<<sigma_x<<endl;

      //////////////////////////////////////
      
      if( sigma_x<x_sigma_cut ) {// cm, big blob, isochronous track
	
	for( auto it_track=vc_track_dir.begin(); it_track!=vc_track_dir.end(); it_track++ ) {// each track cluster
	  int track_id = it_track->first;

	  if( flag_track2recon[track_id]!=0 ) continue;
	  
	  double x_diff_start_end = fabs( vc_track_end[track_id].X() - vc_track_start[track_id].X() );
	  if( x_diff_start_end<2*x_sigma_cut ) {

	    double distance_center2center = fabs( x_bar - 0.5*(vc_track_end[track_id].X() + vc_track_start[track_id].X()) );
	    
	    if(distance_center2center < 2) {// cm
	      cout<<" WARNING : distance between big blob and track < 2 cm"<<endl;


	      bool flag_start_contain = false;
	      bool flag_end_contain = false;
	      
	      for(int i=0; i<npts; i++) {
		double blob_y =  ypt[id].at(i);
		double blob_z =  zpt[id].at(i);
		double track_start_y = vc_track_start[track_id].Y();
		double track_start_z = vc_track_start[track_id].z();
		double track_end_y = vc_track_end[track_id].Y();
		double track_end_z = vc_track_end[track_id].z();

		double distance_start = sqrt( pow(track_start_y-blob_y,2) + pow(track_start_z-blob_z,2) );
		double distance_end = sqrt( pow(track_end_y-blob_y,2) + pow(track_end_z-blob_z,2) );

		if( distance_start<1 ) flag_start_contain = true; // cm
		if( distance_end<1 ) flag_end_contain = true; // cm
		  
		if( flag_start_contain && flag_end_contain ) break;
	      }

	      if( flag_start_contain && flag_end_contain ) {// good track
		flag_track2recon[track_id] = 2;		
		vc_ghost_id.erase( cluster_id );
		map_non_ghost_track_clusters[track_id].push_back( cluster_id );

		///
		int i = cluster_id;
		int npts = xpt[i].size();
		double start_proj = 1e6;
		double end_proj  = -1e6;
		
		for(int j=0; j<npts; j++) {//each point
		  TVector3 point(xpt[i].at(j), ypt[i].at(j), zpt[i].at(j));
		  double point_proj = vc_track_dir[track_id].Dot(point);
		  if(point_proj < start_proj)
		    {
		      start_proj = point_proj;
		    }
		  if(point_proj > end_proj)
		    {
		      end_proj = point_proj;
		    } 
		}//each point
		vc_cluster_length[cluster_id] = (end_proj-start_proj);
		
	      }
	      else if(flag_start_contain || flag_end_contain) {// broken
		flag_track2recon[track_id] = 1;		
		vc_ghost_id.erase( cluster_id );
		map_non_ghost_track_clusters[track_id].push_back( cluster_id );	
	      }
	      
	    }// if(distance_center2center < 2) 
	  }// if( x_diff_start_end<2*x_sigma_cut )
	}// for( auto it_track=vc_track_dir.begin(); it_track!=vc_track_dir.end(); it_track++ )
      }// if( sgima_x<1 )
    }// if( vc_cluster_length[id] > 10 ) 
  }// 
  
  cout<<endl;

  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////

 
  for( auto it=vc_cluster_dir.begin(); it!=vc_cluster_dir.end(); it++ ) {
    
    int id = it->first;
    int cluster_id = id;
    if( vc_ghost_id.find(id)==vc_ghost_id.end() ) continue;

    //cout<<" ---> debug "<<cluster_id<<"\t"<<vc_ghost_id.size()<<endl;
    
    for( auto it_track=vc_track_dir.begin(); it_track!=vc_track_dir.end(); it_track++ ) {// each track cluster
      int track_id = it_track->first;
      //if( flag_track2recon[track_id]==2 ) continue;

      bool flag_start_near = false;
      bool flag_end_near = false;
      
      int npts = xpt[id].size();
      int count_npts = 0;
      
      for(int idx=0; idx<npts; idx++) {
	TVector3 ghost_point( xpt[id].at(idx), ypt[id].at(idx), zpt[id].at(idx) );
	double distance_2start = ( vc_track_start[track_id]-ghost_point ).Mag();
	double distance_2end = ( vc_track_end[track_id]-ghost_point ).Mag();
	if( distance_2start<4 ) flag_start_near = true;// cm
	if( distance_2end<4 ) flag_end_near = true;// cm
	
	/*         

		   - test_A
		   -
		   -
		   <---------
		   start   end



		   - test_B
		   -
		   -
		   ---------->
		   start     end

	*/

	
	TVector3 this_point = ghost_point;
	TVector3 this_point_2_next_start = vc_track_start[track_id] - this_point;
        double distance_vertical = ( this_point_2_next_start.Cross( vc_track_dir[track_id] ) ).Mag();
	double angle_end = (this_point-vc_track_end[track_id]).Dot( -1.*vc_track_dir[track_id] );
	double angle_start = (this_point-vc_track_start[track_id]).Dot( vc_track_dir[track_id] );
	double distance = 1e6;
	if( angle_end * angle_start<0 ) distance = (distance_2start<distance_2end) ? distance_2start : distance_2end;
	else distance = distance_vertical;

	if( distance < 4 ) {// cm
	  count_npts++;
	}
	
      }// for(int idx=0; idx<npts; idx++)

      if( count_npts==0 ) continue;
      cout<<" ------------> percentage of points of cluster near to track : "<<count_npts*1./npts<<endl;

      
      if( (count_npts*1./npts < 0.1 || vc_cluster_length[cluster_id]>5) && count_npts*1./npts < 0.15 )  continue; // equivalent to an angle cut


      if( flag_start_near && flag_end_near ) {// good track
	flag_track2recon[track_id] = 2;
	vc_ghost_id.erase( cluster_id );
	map_non_ghost_track_clusters[track_id].push_back( cluster_id );
	//cout<<" WARNING : big ghost in good track "<<cluster_id<<"\t"<<vc_cluster_length[cluster_id]<<endl;


	///
	int i = cluster_id;
	int npts = xpt[i].size();
	double start_proj = 1e6;
	double end_proj  = -1e6;
		
	for(int j=0; j<npts; j++) {//each point
	  TVector3 point(xpt[i].at(j), ypt[i].at(j), zpt[i].at(j));
	  double point_proj = vc_track_dir[track_id].Dot(point);
	  if(point_proj < start_proj)
	    {
	      start_proj = point_proj;
	    }
	  if(point_proj > end_proj)
	    {
	      end_proj = point_proj;
	    } 
	}//each point
	vc_cluster_length[cluster_id] = (end_proj-start_proj);
	cout<<" WARNING : big ghost "<<cluster_id<<" in good track "<<track_id<<"\t"<<vc_cluster_length[cluster_id]<<endl;
	
      }
      else {// broken
	if( flag_track2recon[track_id]!=2 ) flag_track2recon[track_id] = 1;
	vc_ghost_id.erase( cluster_id );
	map_non_ghost_track_clusters[track_id].push_back( cluster_id );
	cout<<" WARNING : big ghost "<<cluster_id<<" in broken track "<<track_id<<"\t"<<vc_cluster_length[cluster_id]<<endl;
      }
      
    }
    
  }

 
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
 
  for( auto it=map_non_ghost_track_clusters.begin(); it!=map_non_ghost_track_clusters.end(); it++ ) {
    int user_size = it->second.size() - 1;
    cout<<it->first<<"\t"<<user_size<<endl;
  }
  
  /////////////////////////////////////////////////////////////////////////////////// Fill
  ///////////////////////////////////////////////////////////////////////////////////
     
  // std::map<int, TVector3> vc_cluster_dir;
  // std::map<int, TVector3> vc_cluster_start;
  // std::map<int, TVector3> vc_cluster_end;
  // std::map<int, double> vc_cluster_length;
 
  // std::map<int, TVector3> vc_track_dir;
  // std::map<int, TVector3> vc_track_start;
  // std::map<int, TVector3> vc_track_end;
  // std::map<int, double> vc_track_length;
	
  roostr = "h1_flag_track_ghost";
  TH1I *h1_flag_track_ghost = new TH1I( roostr, roostr, 5, 0.5, 5.5);// ghost, total_track, good, broken, inefficiency
  
  roostr = "h1_good_track_length";
  TH1D *h1_good_track_length = new TH1D(roostr, roostr, 150, 0, 150);
 
  roostr = "h1_broken_track_length";
  TH1D *h1_broken_track_length = new TH1D(roostr, roostr, 150, 0, 150);

  roostr = "h2_good_track_phi_costheta";
  TH2D *h2_good_track_phi_costheta = new TH2D(roostr, roostr, 10, 0, TMath::Pi(), 10, 0, 1);
  roostr = "h2_broken_track_phi_costheta";
  TH2D *h2_broken_track_phi_costheta = new TH2D(roostr, roostr, 10, 0, TMath::Pi(), 10, 0, 1);
  roostr = "h2_ineff_track_phi_costheta";
  TH2D *h2_ineff_track_phi_costheta = new TH2D(roostr, roostr, 10, 0, TMath::Pi(), 10, 0, 1);
  
  roostr = "h2_good_track_u_phi_costheta";
  TH2D *h2_good_track_u_phi_costheta = new TH2D(roostr, roostr, 10, 0, TMath::Pi(), 10, 0, 1);
  roostr = "h2_broken_track_u_phi_costheta";
  TH2D *h2_broken_track_u_phi_costheta = new TH2D(roostr, roostr, 10, 0, TMath::Pi(), 10, 0, 1);
  roostr = "h2_ineff_track_u_phi_costheta";
  TH2D *h2_ineff_track_u_phi_costheta = new TH2D(roostr, roostr, 10, 0, TMath::Pi(), 10, 0, 1);
  
  roostr = "h2_good_track_v_phi_costheta";
  TH2D *h2_good_track_v_phi_costheta = new TH2D(roostr, roostr, 10, 0, TMath::Pi(), 10, 0, 1);
  roostr = "h2_broken_track_v_phi_costheta";
  TH2D *h2_broken_track_v_phi_costheta = new TH2D(roostr, roostr, 10, 0, TMath::Pi(), 10, 0, 1);
  roostr = "h2_ineff_track_v_phi_costheta";
  TH2D *h2_ineff_track_v_phi_costheta = new TH2D(roostr, roostr, 10, 0, TMath::Pi(), 10, 0, 1);
  
  int num_ghost = vc_ghost_id.size();
  int num_tracks = flag_track2recon.size();
  int num_flag[3] = {0};
  for(auto it=flag_track2recon.begin(); it!=flag_track2recon.end(); it++ ) {
    int flag = it->second;
    num_flag[flag]++;

    int track_id = it->first;

    TVector3 dir = vc_track_dir[track_id];
    double cosy = dir.Y()/dir.Mag();
    double phiz = TMath::ACos(dir.Z()/TMath::Sqrt(dir.X()*dir.X()+dir.Z()*dir.Z()));
 
    /////// Muplane
    TVector3 u_dir = Muplane * dir;
    if( u_dir.Y()<0 ) u_dir *= -1.;
    double u_cosy = u_dir.Y()/u_dir.Mag();
    double u_phiz = TMath::ACos(u_dir.Z()/TMath::Sqrt(u_dir.X()*u_dir.X()+u_dir.Z()*u_dir.Z()));
    
    /////// Mvplane
    TVector3 v_dir = Mvplane * dir;
    if( v_dir.Y()<0 ) v_dir *= -1.;
    double v_cosy = v_dir.Y()/v_dir.Mag();
    double v_phiz = TMath::ACos(v_dir.Z()/TMath::Sqrt(v_dir.X()*v_dir.X()+v_dir.Z()*v_dir.Z()));
   
    ///////
    if( flag==2 ) {
      // int cluster_id = map_non_ghost_track_clusters[it->first].at(1);
      // double length = vc_cluster_length[cluster_id] * fabs( vc_cluster_dir[cluster_id].Dot(vc_track_dir[track_id]) );      
      // h1_good_track_length->Fill( length );

      int user_size = map_non_ghost_track_clusters[it->first].size() - 1;
      
      if( user_size==1 ) {
	int cluster_id = map_non_ghost_track_clusters[it->first].at(user_size);
	//double length = vc_cluster_length[cluster_id] * fabs( vc_cluster_dir[cluster_id].Dot(vc_track_dir[track_id]) );
	h1_good_track_length->Fill( vc_cluster_length[cluster_id] );
      }
      else {
	double length_max = 0;
	for(int idx=1; idx<=user_size; idx++) {
	  int cluster_id = map_non_ghost_track_clusters[it->first].at(idx);
	  //double length = vc_cluster_length[cluster_id] * fabs( vc_cluster_dir[cluster_id].Dot(vc_track_dir[track_id]) );
	  double length = vc_cluster_length[cluster_id];
	  if( length>length_max ) length_max = length;
	}

	h1_good_track_length->Fill( length_max );
      }
	
    }

    if( flag==1 ) {
      //store each cluster's min/max projective position along track direction
      //then sorted by the min
      std::vector< std::pair<double, double> > vec_min_max_proj;  
      for(int i=1; i<map_non_ghost_track_clusters[it->first].size(); i++) {
	int cluster_id = map_non_ghost_track_clusters[it->first].at(i);
	//length += vc_cluster_length[cluster_id] * fabs( vc_cluster_dir[cluster_id].Dot(vc_track_dir[track_id]) );
	//cout<<" debug : "<< it->first<<"\t"<<vc_cluster_length[ map_non_ghost_track_clusters[it->first].at(i) ]<<" "<<map_non_ghost_track_clusters[it->first].at(i)<<endl;
        double start_proj = vc_cluster_start[cluster_id].Dot(vc_track_dir[track_id]);
        double end_proj = vc_cluster_end[cluster_id].Dot(vc_track_dir[track_id]); 
        double min = 1e6;
        double max = -1e6;
        min = min<start_proj?min:start_proj;
        min = min<end_proj?min:end_proj;
        max = max>start_proj?max:start_proj;
        max = max>end_proj?max:end_proj;
       
        //cout<<"Debug: cluster_id: "<<cluster_id<<endl;
        //cout<<"Debug: start: "<<vc_cluster_start[cluster_id].X()<<" "<<vc_cluster_start[cluster_id].Y()<<" "<<vc_cluster_start[cluster_id].Z()<<endl;
        //cout<<"Debug: end: "<<vc_cluster_end[cluster_id].X()<<" "<<vc_cluster_end[cluster_id].Y()<<" "<<vc_cluster_end[cluster_id].Z()<<endl;
        //cout<<"Debug: min: "<<min<<" max: "<<max<<endl;
        std::pair<double, double> p(min, max);
        vec_min_max_proj.push_back(p);
      }

      //sorted by ascending order of fist element of pair, i.e. min
      std::sort(vec_min_max_proj.begin(), vec_min_max_proj.end());
      double previous_min=-1e6;
      double previous_max=-1e6;
      double length = 0;
      cout<<"Broken: Debug: track: "<<track_id<<endl;
      for(int i=0; i<vec_min_max_proj.size(); i++)
	{
	  double cluster_min = vec_min_max_proj.at(i).first;
	  double cluster_max = vec_min_max_proj.at(i).second;
	  cout<<"Debug: min: "<<cluster_min<<" max: "<<cluster_max<<endl;
	  if(previous_max>=cluster_min){
            if(previous_max<cluster_max) previous_max = cluster_max;
	  }
	  else{
            length += previous_max - previous_min;
            previous_max = cluster_max;
            previous_min = cluster_min;
	  }
	}
      length += previous_max - previous_min;

      cout<<"Debug: length: "<<length<<endl;
      h1_broken_track_length->Fill( length );
      if(length>110) cout<<" WARNING : "<<"broken track length too large!"<<endl;
    }

    double val_u = u_cosy * fabs( TMath::Pi()/2 - u_phiz );
    double val_v = v_cosy * fabs( TMath::Pi()/2 - v_phiz );
    bool flag_u = false;
    if( val_u < val_v ) flag_u = true;
    
    ///////
    if( flag==2 ) {
      h2_good_track_phi_costheta->Fill( phiz, cosy );
      if( flag_u ) h2_good_track_u_phi_costheta->Fill( u_phiz, u_cosy );
      else h2_good_track_v_phi_costheta->Fill( v_phiz, v_cosy );
    }
    if( flag==1 ) {
      h2_broken_track_phi_costheta->Fill( phiz, cosy );
      if( flag_u ) h2_broken_track_u_phi_costheta->Fill( u_phiz, u_cosy );
      else h2_broken_track_v_phi_costheta->Fill( v_phiz, v_cosy );
    }
    if( flag==0 ) {
      h2_ineff_track_phi_costheta->Fill( phiz, cosy );
      if( flag_u ) h2_ineff_track_u_phi_costheta->Fill( u_phiz, u_cosy );
      else h2_ineff_track_v_phi_costheta->Fill( v_phiz, v_cosy );      
    }
      
  }
  int num_good = num_flag[2];
  int num_broken = num_flag[1];
  int num_ineff = num_flag[0];

  h1_flag_track_ghost->SetBinContent( 1, num_ghost );
  h1_flag_track_ghost->SetBinContent( 2, num_tracks );
  h1_flag_track_ghost->SetBinContent( 3, num_good );
  h1_flag_track_ghost->SetBinContent( 4, num_broken );
  h1_flag_track_ghost->SetBinContent( 5, num_ineff );

  //////
  //////

  
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  
  // TCanvas *canv = new TCanvas("canv", "canv", 1000, 800);
  // TH3D *h3 = new TH3D("h3", "h3", 100, 0, 300, 100, -150, 150, 100, 0, 1500);
  // h3->SetStats(0);
  // h3->Draw();
  // h3->SetXTitle("X axis");
  // h3->SetYTitle("Y axis");
  // h3->SetZTitle("Z axis");

  roostr = "h1_ghost_track_length";
  TH1D *h1_ghost_track_length = new TH1D(roostr, roostr, 150, 0, 150);

  roostr = "h2_long_ghost_y_phi_costheta";
  TH2D *h2_long_ghost_y_phi_costheta = new TH2D(roostr, roostr, 10, 0, TMath::Pi(), 10, 0, 1);
  roostr = "h2_long_ghost_u_phi_costheta";
  TH2D *h2_long_ghost_u_phi_costheta = new TH2D(roostr, roostr, 10, 0, TMath::Pi(), 10, 0, 1);
  roostr = "h2_long_ghost_v_phi_costheta";
  TH2D *h2_long_ghost_v_phi_costheta = new TH2D(roostr, roostr, 10, 0, TMath::Pi(), 10, 0, 1);
  
  roostr = "h2_ghost_y_phi_costheta";
  TH2D *h2_ghost_y_phi_costheta = new TH2D(roostr, roostr, 10, 0, TMath::Pi(), 10, 0, 1);
  roostr = "h2_ghost_u_phi_costheta";
  TH2D *h2_ghost_u_phi_costheta = new TH2D(roostr, roostr, 10, 0, TMath::Pi(), 10, 0, 1);
  roostr = "h2_ghost_v_phi_costheta";
  TH2D *h2_ghost_v_phi_costheta = new TH2D(roostr, roostr, 10, 0, TMath::Pi(), 10, 0, 1);
  
  int size_ghost = vc_ghost_id.size();
  cout<<" ---> number of ghost "<<size_ghost<<endl;

  long count_graph_ghost_yz = 0;
  for( auto it=vc_ghost_id.begin(); it!=vc_ghost_id.end(); it++ ) {
    int id = it->first;
    cout<<TString::Format( " ---> ghost %3d, length %7.4f", id, vc_cluster_length[id] )<<endl;

    int npts = xpt[id].size();
    for(int j=0; j<npts; j++) { // each point      
      double x = xpt[id].at(j);
      double y = ypt[id].at(j);
      double z = zpt[id].at(j);
      count_gh_ghost_all++;
      gh_ghost_all->SetPoint( count_gh_ghost_all-1, x,y,z );

      count_graph_ghost_yz++;
      graph_ghost_yz->SetPoint(count_graph_ghost_yz-1, z,y);      
    }


    // if( id!=9 && id!=26 ) continue;
    // if( id!=13 && id!=16 && id!=18 && id!=19 ) continue;

    /////// ghost histgram
    /////// ghost histgram
    h1_ghost_track_length->Fill( vc_cluster_length[id] );

    cout<<"length_ghost "<<vc_cluster_length[id]<<"\t"<<inputroot<<endl;
    
    /////// ghost-angle
    /////// ghost-angle
    TVector3 dir = vc_cluster_dir[id];
    double cosy = dir.Y()/dir.Mag();
    double phiz = TMath::ACos(dir.Z()/TMath::Sqrt(dir.X()*dir.X()+dir.Z()*dir.Z()));
 
    /////// Muplane
    TVector3 u_dir = Muplane * dir;
    if( u_dir.Y()<0 ) u_dir *= -1.;
    double u_cosy = u_dir.Y()/u_dir.Mag();
    double u_phiz = TMath::ACos(u_dir.Z()/TMath::Sqrt(u_dir.X()*u_dir.X()+u_dir.Z()*u_dir.Z()));
    
    /////// Mvplane
    TVector3 v_dir = Mvplane * dir;
    if( v_dir.Y()<0 ) v_dir *= -1.;
    double v_cosy = v_dir.Y()/v_dir.Mag();
    double v_phiz = TMath::ACos(v_dir.Z()/TMath::Sqrt(v_dir.X()*v_dir.X()+v_dir.Z()*v_dir.Z()));

    ///////
    h2_ghost_y_phi_costheta->Fill( phiz, cosy );
    h2_ghost_u_phi_costheta->Fill( u_phiz, u_cosy );
    h2_ghost_v_phi_costheta->Fill( v_phiz, v_cosy );
    
    /////// long ghost: length > 20 cm
    if( vc_cluster_length[id]>20 ) {
      h2_long_ghost_y_phi_costheta->Fill( phiz, cosy );
      h2_long_ghost_u_phi_costheta->Fill( u_phiz, u_cosy );
      h2_long_ghost_v_phi_costheta->Fill( v_phiz, v_cosy );
    }
    

  }

  // if( gh_ghost_all->GetN()!=0 ) gh_ghost_all->Draw("same p");
    
  // gh_truth->SetMarkerColor(kRed);
  // gh_truth->Draw("same p");
    
  // canv->SaveAs("canv.root");



  /////////////////////////////////////////////////////////////////////////////////// drawgood
  ///////////////////////////////////////////////////////////////////////////////////

  roostr = "graph_good";
  TH2D *graph_good = new TH2D(roostr, roostr, 550, 0, 1100, 120, -120, 120);
  // TGraph *graph_good = new TGraph();
  // graph_good->SetName(roostr);
  
  roostr = "graph_broken";
  TH2D *graph_broken = new TH2D(roostr, roostr, 550, 0, 1100, 120, -120, 120);
  // TGraph *graph_broken = new TGraph();
  // graph_broken->SetName(roostr);
  
  // roostr = "graph_noGhost";
  // TGraph *graph_noGhost = new TGraph();
  // graph_noGhost->SetName(roostr);

  long count_good = 0;
  long count_broken = 0;
  long count_noGhost = 0;

  // std::map< int, std::vector<int> > map_non_ghost_track_clusters;
  // map<int, int>flag_track2recon;// 0: inefficiency, 1: broken, 2: good
  // std::map<int, std::vector<double>> xpt;

  /// chargeTTT
  map<int, double>map_charge_TotMean_goodtrk;// [cluster_id]
  map<int, int>map_charge_TotMean_goodtrk_from_clusterid;
  

  for(auto it_tk=flag_track2recon.begin(); it_tk!=flag_track2recon.end(); it_tk++ ) {
    int track_id = it_tk->first;
    int flag = it_tk->second;
    if( flag==0 ) continue;

    int user_size = map_non_ghost_track_clusters[track_id].size() - 1;
    //cout<<TString::Format(" ---> %2d %2d %3d", track_id, flag, user_size)<<endl;

    for(int ic=1; ic<=user_size; ic++) {
      int cluster_id = map_non_ghost_track_clusters[track_id].at(ic);

      int cluster_size = xpt[cluster_id].size();
      for(int idx=0; idx<cluster_size; idx++) {
	
	if( flag==1 ) {
	  count_broken++;
	  //graph_broken->SetPoint(count_broken-1, zpt[cluster_id].at(idx), ypt[cluster_id].at(idx) );
	  graph_broken->Fill( zpt[cluster_id].at(idx), ypt[cluster_id].at(idx) );
	}
	
	if( flag==2 ) {
	  count_good++;
	  //graph_good->SetPoint(count_good-1, zpt[cluster_id].at(idx), ypt[cluster_id].at(idx) );
	  graph_good->Fill( zpt[cluster_id].at(idx), ypt[cluster_id].at(idx) );

	  /// chargeTTT
	  //map_charge_TotMean_goodtrk[cluster_id] += qpt[cluster_id].at(idx)/vc_cluster_length[cluster_id]/10; // 1m = 1000 mm; cm --> mm
	  map_charge_TotMean_goodtrk[cluster_id] += qpt[cluster_id].at(idx);
	  map_charge_TotMean_goodtrk_from_clusterid[cluster_id] = track_id;
	}
	
	count_noGhost++;
	//graph_noGhost->SetPoint(count_noGhost-1, zpt[cluster_id].at(idx), ypt[cluster_id].at(idx) );
      }
    }

    
  }// for(auto it_tk=flag_track2recon.begin(); it_tk!=flag_track2recon.end(); it_tk++ )


  
  roostr = "graph_truth";
  TH2D *graph_truth = new TH2D(roostr, roostr, 550, 0, 1100, 120, -120, 120);

  roostr = "graph_truth_good";
  TH2D *graph_truth_good = new TH2D(roostr, roostr, 550, 0, 1100, 120, -120, 120);
  
  roostr = "graph_truth_broken";
  TH2D *graph_truth_broken = new TH2D(roostr, roostr, 550, 0, 1100, 120, -120, 120);
  
  
  for(long ientry=0; ientry<truth->GetEntries(); ientry++) {
    truth->GetEntry(ientry);
    
    graph_truth->Fill( z,y );

    if( flag_track2recon[track_id]==2 ) graph_truth_good->Fill( z,y );
    if( flag_track2recon[track_id]==1 ) graph_truth_broken->Fill( z,y );
  }

  /// chargeTTT
  roostr = "h1_charge_TotMean_goodtrack";
  TH1D *h1_charge_TotMean_goodtrack = new TH1D(roostr, roostr, 600, 0, 60000);

  roostr = "h1_charge_TotMean_goodtrack_y_phi_costheta";
  TH2D *h1_charge_TotMean_goodtrack_y_phi_costheta = new TH2D(roostr, roostr,  10, 0, TMath::Pi(), 10, 0, 1);
  
  for(auto it_tk=map_charge_TotMean_goodtrk.begin(); it_tk!=map_charge_TotMean_goodtrk.end(); it_tk++ ) {
    // cout<<" ---> "<<it_tk->first<<"\t"<<it_tk->second<<endl;
    int cluster_id = it_tk->first;
    double TotMean = it_tk->second;

    h1_charge_TotMean_goodtrack->Fill( TotMean );

    cout<<" ------> chargeTTT, length and meanCharge :"<<vc_cluster_length[cluster_id]/10<<"\t"<<TotMean<<endl;
    
    //////////////////////////////////////
    
    int track_id = map_charge_TotMean_goodtrk_from_clusterid[cluster_id];
    
    TVector3 dir = vc_track_dir[track_id];
    double cosy = dir.Y()/dir.Mag();
    double phiz = TMath::ACos(dir.Z()/TMath::Sqrt(dir.X()*dir.X()+dir.Z()*dir.Z()));
 
    /////// Muplane
    TVector3 u_dir = Muplane * dir;
    if( u_dir.Y()<0 ) u_dir *= -1.;
    double u_cosy = u_dir.Y()/u_dir.Mag();
    double u_phiz = TMath::ACos(u_dir.Z()/TMath::Sqrt(u_dir.X()*u_dir.X()+u_dir.Z()*u_dir.Z()));
    
    /////// Mvplane
    TVector3 v_dir = Mvplane * dir;
    if( v_dir.Y()<0 ) v_dir *= -1.;
    double v_cosy = v_dir.Y()/v_dir.Mag();
    double v_phiz = TMath::ACos(v_dir.Z()/TMath::Sqrt(v_dir.X()*v_dir.X()+v_dir.Z()*v_dir.Z()));

    ///////
    h1_charge_TotMean_goodtrack_y_phi_costheta->Fill( phiz, cosy, TotMean );
    
  }
  

  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////

  if(cluster_check){
    roostr = "canv_visuliaztion_tracks_cluster";
    TCanvas *canv_visuliaztion_tracks_cluster = new TCanvas(roostr, roostr, 1000, 800);

    roostr = "h3_visuliaztion_tracks_cluster";
    TH3D *h3_visuliaztion_tracks_cluster = new TH3D(roostr, roostr, 100, 0, 300, 100, -150, 150, 100, 0, 1500);
    h3_visuliaztion_tracks_cluster->SetStats(0);
    h3_visuliaztion_tracks_cluster->Draw();
    h3_visuliaztion_tracks_cluster->SetXTitle("X axis");
    h3_visuliaztion_tracks_cluster->SetYTitle("Y axis");
    h3_visuliaztion_tracks_cluster->SetZTitle("Z axis");

    int colors_vis[3] = {kRed, kOrange-3, kGreen};// inefficiency, broken, good
  
    for( auto it=vc_track_dir.begin(); it!=vc_track_dir.end(); it++ ) {
      int id = it->first;
    
      graph_track_each[id]->Draw("same p");
      graph_track_each[id]->SetMarkerColor( colors_vis[flag_track2recon[id]] );
      graph_track_each[id]->SetMarkerStyle(4);
      graph_track_each[id]->SetMarkerSize(0.5);
    }

    graph_cluster_all->Draw("same p");
    graph_cluster_all->SetMarkerColor(kGray+2);

    if( gh_ghost_all->GetN()!=0 ) {
      gh_ghost_all->Draw("same p");
      gh_ghost_all->SetMarkerColor(kBlue);
      gh_ghost_all->SetMarkerStyle(4);
      gh_ghost_all->SetMarkerSize(0.3);
    }

    canv_visuliaztion_tracks_cluster->SaveAs("canv_visuliaztion_tracks_cluster.root");
  }
  /////////////////////////////////////////////////////////////////////////////////// WRITE FILE
  /////////////////////////////////////////////////////////////////////////////////// WRITE FILE
  

  cout<<" --------> Number of tracks "<< vc_track_dir.size() <<endl;
  cout<<" --------> Number of ghosts "<<size_ghost<<endl;

  cout<<" file "<<inputroot<<endl;
  
  TFile* output = new TFile(outputroot, "RECREATE");

  ///////
    
  graph_good->SetMarkerColor(kGreen);
  graph_broken->SetMarkerColor(kBlue);
  
  graph_truth->Write();
  graph_truth_good->Write();
  graph_truth_broken->Write();
  
  graph_good->Write();
  graph_broken->Write();
  //graph_noGhost->Write();
  graph_ghost_yz->Write();
  h1_ghost_track_length->Write();
  
  ///////
  h1_flag_track_ghost->Write();
  h1_good_track_length->Write();
  h1_broken_track_length->Write();

  h2_good_track_phi_costheta->Write();
  h2_broken_track_phi_costheta->Write();
  h2_ineff_track_phi_costheta->Write();
    
  h2_good_track_u_phi_costheta->Write();
  h2_broken_track_u_phi_costheta->Write();
  h2_ineff_track_u_phi_costheta->Write();
    
  h2_good_track_v_phi_costheta->Write();
  h2_broken_track_v_phi_costheta->Write();
  h2_ineff_track_v_phi_costheta->Write();

  h2_long_ghost_y_phi_costheta->Write();
  h2_long_ghost_u_phi_costheta->Write();
  h2_long_ghost_v_phi_costheta->Write();
  
  h2_ghost_y_phi_costheta->Write();
  h2_ghost_u_phi_costheta->Write();
  h2_ghost_v_phi_costheta->Write();

  h1_charge_TotMean_goodtrack->Write();
  h1_charge_TotMean_goodtrack_y_phi_costheta->Write();
    
  output->Close();

  return 0;

 
  // numbers of cycle = C_n^2 =  num_clusters*(num_clusters-1)/2 

  // int count_it_a = 0;
  // for( auto it_a=vc_cluster_dir.begin(); it_a!=std::prev( vc_cluster_dir.end(), 1); it_a++ ) {
  //   int cluster_id_a = it_a->first;
  //   cout<<TString::Format(" ---> cluster and length: %3d, %6.2f", cluster_id_a, vc_cluster_length[cluster_id_a] )<<endl;

  //   ///////
  //   count_it_a++;
  //   for( auto it_b=std::prev( vc_cluster_dir.end(), vc_cluster_dir.size()-count_it_a ); it_b!=vc_cluster_dir.end(); it_b++ ) {
  //     int cluster_id_b = it_b->first;
  //     cout<<TString::Format("   -------> %3d, %6.2f", cluster_id_b, vc_cluster_length[cluster_id_b])<<endl;
  //   }
  // }

  // int count_it_a = 0;
  // for( auto it_a=vc_cluster_dir.begin(); it_a!=vc_cluster_dir.end(); it_a++ ) {
  //   int cluster_id_a = it_a->first;
  //   //cout<<TString::Format(" ---> cluster and length: %3d, %6.2f", cluster_id_a, vc_cluster_length[cluster_id_a] )<<endl;

  //   TVector3 vector_unit = 1./vc_cluster_dir[cluster_id_a].Mag() * vc_cluster_dir[cluster_id_a];
  //   vector_unit.Print();
    
  //   for( auto it_b=vc_cluster_dir.begin(); it_b!=vc_cluster_dir.end(); it_b++ ) {
  //     int cluster_id_b = it_b->first;
  //     if( cluster_id_a==cluster_id_b ) continue;
  //     //cout<<TString::Format("   -------> %3d, %6.2f", cluster_id_b, vc_cluster_length[cluster_id_b])<<endl;
      
  //     TVector3 vector_a_unit = 1./vc_cluster_dir[cluster_id_a].Mag() * vc_cluster_dir[cluster_id_a];
  //     TVector3 vector_b_unit = 1./vc_cluster_dir[cluster_id_b].Mag() * vc_cluster_dir[cluster_id_b];

  //     TVector3 vector_axb = vector_a_unit.Cross( vector_b_unit );
  //     if( vector_axb.Mag()==0 ) cout<<TString::Format(" ---> Parallel of cluster %3d and %3d", cluster_id_a, cluster_id_b )<<endl;
      
  //   }
  // }
  
}
