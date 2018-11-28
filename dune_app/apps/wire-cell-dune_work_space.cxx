#include "WireCellNav/DetectorGDS.h"
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

using namespace WireCell;
using namespace std;

int main(int argc, char* argv[])
{

  if (argc < 5){
    cerr << "usage: wire-cell-dune_work_space /path/to/celltree.root eve_num -o[0,1](do_rotation) -p[0,1](is_3mm) -r[0,1](randomly place vertices) -s[int](seed) -v[0,1](original vertex?)" << endl;
    return 1;
  }

  bool rotate_90deg=false;
  bool is_3mm=false;
  bool random_vertices=false;
  bool orig_vertex = false;
  bool cut_events = true;

  int seed=0;
  for (Int_t i = 1; i != argc; i++){
     switch(argv[i][1]){
     case 'o':
       rotate_90deg = atoi(&argv[i][2]); 
       break;
     case 'p':
       is_3mm = atoi(&argv[i][2]); 
       break;
     case 'r':
       random_vertices = atoi(&argv[i][2]);
       break;
     case 's':
       seed = atoi(&argv[i][2]);
       break;
     case 'v':
       orig_vertex = atoi(&argv[i][2]);
       break;
     case 'c':
       cut_events = atoi(&argv[i][2]);
       break;
     }
  }

  if (rotate_90deg) cout<<"Beam is perpendicular to wire planes. ";
  else cout<<"Beam is parallel to wire planes (default). ";
  if (is_3mm) cout<<"Wire pitch is 3 mm."<<endl;
  else cout<<"Wire pitch is 5 mm (default)."<<endl;
  if (orig_vertex) cout << "Original Vertex " << endl;

  double pitchU, pitchV, pitchW;
  if (!is_3mm) {
    pitchU = 0.4667*units::cm;
    pitchV = 0.4667*units::cm;
    pitchW = 0.479*units::cm;
  } else {
    pitchU = 0.3*units::cm;
    pitchV = 0.3*units::cm;
    pitchW = 0.3*units::cm;
  }
  
  //build GDS ... 
  DetectorGDS gds;
  gds.set_ncryos(1);
  gds.set_napas(0,4);
  Vector center0(-0*units::cm, -300.05*units::cm, 115.318875*units::cm);
  Vector halves0(3.995355*units::cm, 299.5*units::cm, 115.318875*units::cm);
  gds.set_apa(0, 0, 35.71*units::deg, 35.71*units::deg, pitchU, pitchV, pitchW, center0, halves0);
  Vector center1(-0*units::cm, 300.05*units::cm, 115.318875*units::cm);
  Vector halves1(3.995355*units::cm, 299.5*units::cm, 115.318875*units::cm);
  gds.set_apa(0, 1, 35.71*units::deg, 35.71*units::deg, pitchU, pitchV, pitchW, center1, halves1);
  Vector center2(-0*units::cm, -300.05*units::cm, 347.709*units::cm);
  Vector halves2(3.995355*units::cm, 299.5*units::cm, 115.318875*units::cm);
  gds.set_apa(0, 2, 35.71*units::deg, 35.71*units::deg, pitchU, pitchV, pitchW, center2, halves2);
  Vector center3(-0*units::cm, 300.05*units::cm, 347.709*units::cm);
  Vector halves3(3.995355*units::cm, 299.5*units::cm, 115.318875*units::cm);
  gds.set_apa(0, 3, 35.71*units::deg, 35.71*units::deg, pitchU, pitchV, pitchW, center3, halves3);
  gds.buildGDS();

  Vector vertex(200,300,100);
  if (random_vertices) {
    TRandom3 rnd(seed);
    vertex.x = rnd.Uniform(-300,300);
    vertex.y = rnd.Uniform(-550,550);
    vertex.z = rnd.Uniform(50,410);
    std::cout<<"new vertex is ("<<vertex.x<<","<<vertex.y<<","<<vertex.z<<")"<<std::endl;
  }
  
  const char* root_file = argv[1];
  const char* tpath = "/Event/Sim";

  int recon_threshold = 2000;
  int max_events = 100;
  int eve_num = atoi(argv[2]);
  float unit_dis = 1.6;  
  float toffset_1=1.647;
  float toffset_2=1.539+1.647;
  float toffset_3=0;

  int save_image_outline_flag = 1; // prescale flag 
  
  int total_time_bin=9600;
  int frame_length = 3200;
  //int frame_length = 800;  // hack for now
  int nrebin = 4;

  TPCParams& mp = Singleton<TPCParams>::Instance();
  
  double pitch_u = gds.get_pitch(0,WirePlaneType_t(0));
  double pitch_v = gds.get_pitch(0,WirePlaneType_t(1));
  double pitch_w = gds.get_pitch(0,WirePlaneType_t(2));
  double time_slice_width = nrebin * unit_dis * 0.5 * units::mm;

  mp.set_pitch_u(pitch_u);
  mp.set_pitch_v(pitch_v);
  mp.set_pitch_w(pitch_w);
  mp.set_ts_width(time_slice_width);
  
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

  if (cut_events){
    //Now judge which event to replay 
    int flag_replay = 0;
    // replay all nue cc for energy < 10 GeV 
    if (mctruth->mc_nu_mom[3]<10) flag_replay = 1;
    // replay all nc events
    if (mctruth->mc_nu_ccnc==1) flag_replay = 1;
    // replay all nu_tau events 
    if (fabs(mctruth->mc_nu_pdg)==16) flag_replay = 1;
    
    if (flag_replay==0)
      return 0;
  }
  /*
  int mc_Ntrack;  // number of tracks in MC
  int mc_id[MAX_TRACKS];  // track id; size == mc_Ntrack
  int mc_pdg[MAX_TRACKS];  // track particle pdg; size == mc_Ntrack
  int mc_process[MAX_TRACKS];  // track generation process code; size == mc_Ntrack
  int mc_mother[MAX_TRACKS];  // mother id of this track; size == mc_Ntrack
  float mc_startXYZT[MAX_TRACKS][4];  // start position of this track; size == mc_Ntrack
  float mc_endXYZT[MAX_TRACKS][4];  // end position of this track; size == mc_Ntrack
  float mc_startMomentum[MAX_TRACKS][4];  // start momentum of this track; size == mc_Ntrack
  float mc_endMomentum[MAX_TRACKS][4];  // end momentum of this track; size == mc_Ntrack
  std::vector<std::vector<int> > *mc_daughters= new std::vector<std::vector<int> >;  // daughters id of this track; vector
  TObjArray* mc_trackPosition = new TObjArray();

  sst->SetBranchAddress("mc_Ntrack", &mc_Ntrack);  // number of tracks in MC
  sst->SetBranchAddress("mc_id", &mc_id);  // track id; size == mc_Ntrack
  sst->SetBranchAddress("mc_pdg", &mc_pdg);  // track particle pdg; size == mc_Ntrack
  sst->SetBranchAddress("mc_process", &mc_process);  // track generation process code; size == mc_Ntrack
  sst->SetBranchAddress("mc_mother", &mc_mother);  // mother id of this track; size == mc_Ntrack
  sst->SetBranchAddress("mc_daughters", &mc_daughters);  // daughters id of this track; vector
  sst->SetBranchAddress("mc_startXYZT", &mc_startXYZT);  // start position of this track; size == mc_Ntrack
  sst->SetBranchAddress("mc_endXYZT", &mc_endXYZT);  // start position of this track; size == mc_Ntrack
  sst->SetBranchAddress("mc_startMomentum", &mc_startMomentum);  // start momentum of this track; size == mc_Ntrack
  sst->SetBranchAddress("mc_endMomentum", &mc_endMomentum);  // start momentum of this track; size == mc_Ntrack
  sst->SetBranchAddress("mc_trackPosition",&mc_trackPosition);
  */
  //sst->GetEntry(eve_num);
  
  cout << "Run No: " << mctruth->runNo << " " << mctruth->subrunNo << " " << mctruth->eventNo << endl;
  
  WireCell::FrameDataSource* fds = 0;
  fds = WireCellSst::make_fds(*tfile);
  if (!fds) {
    cerr << "ERROR: failed to get FDS from " << root_file << endl;
    return 1;
  }
  
  TH1::AddDirectory(kFALSE);
  
  WireCell::ToyDepositor *toydep;

  Double_t x_center = mctruth->mc_startXYZT[0][0];
  Double_t y_center = mctruth->mc_startXYZT[0][1];
  Double_t z_center = mctruth->mc_startXYZT[0][2];
  Double_t x_shift = -mctruth->mc_startXYZT[0][0]+vertex.x;
  Double_t y_shift = -mctruth->mc_startXYZT[0][1]+vertex.y;
  Double_t z_shift = -mctruth->mc_startXYZT[0][2]+vertex.z;
  
  if (orig_vertex){
    x_shift = 0;
    y_shift = 0;
    z_shift = 0;
  }

  if (rotate_90deg) {
    toydep = new WireCell::ToyDepositor(fds,0,unit_dis,frame_length, x_center, y_center, z_center, TMath::Pi()/2.);
  } else {
    toydep = new WireCell::ToyDepositor(fds,0,unit_dis,frame_length);
  }
  cout << "Primary vertex is (" << mctruth->mc_startXYZT[0][0] << "," << mctruth->mc_startXYZT[0][1] << "," << mctruth->mc_startXYZT[0][2] <<")" << endl;
  // if (random_vertices) {
  toydep->translation(x_shift,y_shift,z_shift);
  //}

  if (rotate_90deg){
    mctruth->Rotate_Shift(x_center, y_center, z_center, TMath::Pi()/2.,
			  x_shift,y_shift,z_shift);
  }else{
    mctruth->Rotate_Shift(0,0,0, 0,
			  x_shift,y_shift,z_shift);
  }
  
  const PointValueVector& pvv = toydep->depositions(eve_num);

  std::cout << "Points deposited: " << pvv.size() << std::endl;
  //return 0;
  // DetGenerativeFDS gfds(toydep,gds, 2400,max_events,2.0*1.6*units::millimeter);
  DetGenerativeFDS gfds(*toydep, gds, total_time_bin,max_events,0.5*1.6*units::millimeter);
  gfds.jump(eve_num);

  tfile->Close("R");
  delete tfile;

  //WireCellSst::ToyuBooNESliceDataSource sds(gfds,1500); //set threshold at 2000 electrons

  
  
  WireCell2dToy::ToySignalSimuTrueFDS st_fds(gfds,gds,total_time_bin/nrebin,max_events,0); //truth  
  cout << "Put in Truth " << endl; 
  st_fds.jump(eve_num);
  
  cout << "Simulate Raw WaveForm " << endl; 
  WireCell2dToy::ToySignalSimuFDS simu_fds(gfds,gds,total_time_bin,max_events,toffset_1,toffset_2,1); // time offset among different planes for the time electrons travel among different planes
  simu_fds.jump(eve_num);
  //simu_fds.Save();
  
  cout << "Deconvolution with Gaussian filter" << endl;
  WireCell2dToy::ToySignalGausFDS gaus_fds(simu_fds,gds,total_time_bin/nrebin,max_events,toffset_1,toffset_2); // gaussian smearing for charge estimation
  gaus_fds.jump(eve_num);
  //gaus_fds.Save();

  cout << "Deconvolution with Wiener filter" << endl;
  WireCell2dToy::ToySignalWienFDS wien_fds(simu_fds,gds,total_time_bin/nrebin,max_events,toffset_1,toffset_2); // weiner smearing for hit identification
  wien_fds.jump(eve_num);
  //wien_fds.Save();
  
  
  

  int nwire_u = gds.get_total_nwires(WirePlaneType_t(0));
  int nwire_v = gds.get_total_nwires(WirePlaneType_t(1));
  int nwire_w = gds.get_total_nwires(WirePlaneType_t(2));
  
 
  
  // cin >> abc;

  WireCellSst::ToyuBooNESliceDataSource sds(gds,wien_fds,gaus_fds,threshold_u, 
  					    threshold_v, threshold_w, 
  					    threshold_ug, 
  					    threshold_vg, threshold_wg, 
  					    nwire_u, 
  					    nwire_v, nwire_w); 

  WireCellSst::ToyuBooNESliceDataSource sds_th(gds,st_fds,st_fds,500, 
  					    500, 500, 
  					    threshold_ug, 
  					    threshold_vg, threshold_wg, 
  					    nwire_u, 
  					    nwire_v, nwire_w); 

  sds.set_flag(2);
  sds_th.set_flag(2);
  

   int ncount = 0;
  int ncount1 = 0;  
  int ncount2 = 0;

  int ncount_t = 0;


   WireCell2dToy::ToyTiling **toytiling = new WireCell2dToy::ToyTiling*[2400];
   WireCell2dToy::MergeToyTiling **mergetiling = new WireCell2dToy::MergeToyTiling*[2400];
   WireCell2dToy::TruthToyTiling **truthtiling = new WireCell2dToy::TruthToyTiling*[2400];
   WireCell2dToy::ToyMatrix **toymatrix = new WireCell2dToy::ToyMatrix*[2400];
   
    //save truth ...
   // WireCell2dToy::ToyTiling **toytiling_th = new WireCell2dToy::ToyTiling*[2400];
   // WireCell2dToy::TruthToyTiling **truthtiling_th = new WireCell2dToy::TruthToyTiling*[2400];

   WireCell2dToy::ToyMetric toymetric;
   

   //add in cluster
   GeomClusterSet cluster_set, cluster_delset;
   
   int ncount_mcell = 0;
   
   delete fds;

   int start_num = 0 ;
   int end_num = sds.size()-1;
   
   // start_num = 841;
   // end_num = 842;

  //int i = 317+800;{
   //int i = 292+800;
   for (int i=start_num;i!=end_num+1;i++){
     sds.jump(i);
     sds_th.jump(i);
     
     WireCell::Slice slice = sds.get();
     WireCell::Slice slice_th = sds_th.get();
     cout << i << " " << slice.group().size() << " " << slice_th.group().size() << endl;

     // if ( slice.group().size() >0){
     // cout << i << " " << slice.group().size() << endl;
     toytiling[i] = new WireCell2dToy::ToyTiling(slice,gds,0,0,0,threshold_ug,threshold_vg, threshold_wg);
     //allcell = toytiling[i]->get_allcell();
     GeomCellSelection allcell = toytiling[i]->get_allcell();
     GeomWireSelection allwire = toytiling[i]->get_allwire();
     cout << "Single Cell: " << i << " "  << allcell.size() << " " << allwire.size() << endl;

     mergetiling[i] = new WireCell2dToy::MergeToyTiling(gds,*toytiling[i],i); 
     GeomCellSelection allmcell = mergetiling[i]->get_allcell();
     GeomWireSelection allmwire = mergetiling[i]->get_allwire();
     
     cout <<"Blob: " << i << " " << allmcell.size() << " " << allmwire.size() << endl;
     
     
     truthtiling[i] = new WireCell2dToy::TruthToyTiling(*toytiling[i],pvv,i,gds,frame_length/nrebin,unit_dis);
     
     
     toymatrix[i] = new WireCell2dToy::ToyMatrix(gds,*toytiling[i],*mergetiling[i]);
   
     // origininal iteration method, now want to add penalty ... 
     // if (toymatrix[i]->Get_Solve_Flag()==0){
     //   WireCell2dToy::ToyMatrixIterate toymatrix_it(*toymatrix[i],2000,1e5);
     // }
     
     cout << i << " chi2: " << toymatrix[i]->Get_Chi2() <<
       " NDF: " << toymatrix[i]->Get_ndf() << endl;
     
     // toytiling_th[i] = new WireCell2dToy::ToyTiling(slice_th,gds,0,0,0,threshold_ug,threshold_vg, threshold_wg);
     // truthtiling_th[i] = new WireCell2dToy::TruthToyTiling(*toytiling_th[i],pvv,i,gds,frame_length/nrebin,unit_dis);
    
     CellChargeMap ccmap = truthtiling[i]->ccmap();
     if (toymatrix[i]->Get_Solve_Flag()!=0)
       toymetric.Add(allmcell,*toymatrix[i],ccmap);
     
     toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
     
     //}
     
     //  TApplication theApp("theApp",&argc,argv);
     // theApp.SetReturnFromRun(true);
     
     // TCanvas c1("ToyMC","ToyMC",1200,600);
     // c1.Divide(2,1);
    // c1.Draw();
    
    // float charge_min = 0;
    // float charge_max = 1e5;


    // WireCell2dToy::ToyEventDisplay display(c1, gds);
    // display.charge_min = charge_min;
    // display.charge_max = charge_max;


    // gStyle->SetOptStat(0);
    
    // const Int_t NRGBs = 5;
    // const Int_t NCont = 255;
    // Int_t MyPalette[NCont];
    // Double_t stops[NRGBs] = {0.0, 0.34, 0.61, 0.84, 1.0};
    // Double_t red[NRGBs] = {0.0, 0.0, 0.87 ,1.0, 0.51};
    // Double_t green[NRGBs] = {0.0, 0.81, 1.0, 0.2 ,0.0};
    // Double_t blue[NRGBs] = {0.51, 1.0, 0.12, 0.0, 0.0};
    // Int_t FI = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    // gStyle->SetNumberContours(NCont);
    // for (int kk=0;kk!=NCont;kk++) MyPalette[kk] = FI+kk;
    // gStyle->SetPalette(NCont,MyPalette);
    
    // display.init(-0.03,1.568,-0.845,1.151);
    // display.draw_mc(1,WireCell::PointValueVector(),"colz");
    // display.draw_slice(slice,"");
    // display.draw_cells(toytiling[i]->get_allcell(),"*same");

    
    // //std::cout << ccmap.size() << std::endl;
    // display.draw_truthcells(ccmap,"*same");
    // display.draw_mergecells(mergetiling[i]->get_allcell(),"*same",0); //0 is normal, 1 is only draw the ones containt the truth cell

    
    // theApp.Run();
   }
   
   // WireCell2dToy::TotalTiling& totaltiling = Singleton<WireCell2dToy::TotalTiling>::Instance();
   // totaltiling.Clear();
   // int qx_abc;
   // std::cout << "Stop, check memory!" << std::endl;
   // std::cin >> qx_abc;




   double penalty = 6;
   double penalty_ncpt = 0.1;

   
   // figure out how to add chi^2 in the solving ... 
   std::cout << "Starting to use connectivitiy" << std::endl;
   
   std::list<int> solve_list;
   
   if (start_num != end_num){
     int first_solve=-1;
     for (int i=start_num; i!=end_num+1;i++){
       if (toymatrix[i]->Get_Solve_Flag()!=0){
   	 first_solve = i;
   	 break;
       }
     }
     if (first_solve == -1) first_solve = start_num;

    
     for (int i=first_solve+1;i<=end_num-1;i++){
       if (toymatrix[i]->Get_Solve_Flag()==0 && toymatrix[i]->Get_mcindex()>0){ 
	 if ((toymatrix[i-1]->Get_Solve_Flag() + toymatrix[i+1]->Get_Solve_Flag()) > 0){
	   WireCell2dToy::ToyMatrixIterate toymatrix_it(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],2000,1e5,penalty,penalty_ncpt);
	   
	   GeomCellSelection allmcell = mergetiling[i]->get_allcell();
	   CellChargeMap ccmap = truthtiling[i]->ccmap();
	   if (toymatrix[i]->Get_Solve_Flag()!=0)
	     toymetric.Add(allmcell,*toymatrix[i],ccmap);
	   toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
	   
	   cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
	   cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
	 }else{
	   solve_list.push_back(i); 
	 }
       }
     }
     
     for (int i=first_solve-1;i>=start_num+1;i--){
       if (toymatrix[i]->Get_Solve_Flag()==0&& toymatrix[i]->Get_mcindex()>0){
	 if ((toymatrix[i-1]->Get_Solve_Flag() + toymatrix[i+1]->Get_Solve_Flag()) > 0){
	   WireCell2dToy::ToyMatrixIterate toymatrix_it(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],2000,1e5,penalty,penalty_ncpt);
	 
	   GeomCellSelection allmcell = mergetiling[i]->get_allcell();
	   CellChargeMap ccmap = truthtiling[i]->ccmap();
	   if (toymatrix[i]->Get_Solve_Flag()!=0)
	     toymetric.Add(allmcell,*toymatrix[i],ccmap);
	   toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
	 
	   cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
	   cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
	 }else{
	   solve_list.push_back(i);
	 }
       }
     }
   }
   
   // start second round ...
   // std::cout << "Connectivitiy rest " << solve_list.size() << std::endl;
   int prev_count = 0;
   while (solve_list.size() >0){
     int curr_count = solve_list.size();
     std::cout << "Connectivitiy rest " << solve_list.size() << std::endl;

     if (curr_count == prev_count){
       int i = solve_list.front(); // pick the first element ... 
       if (toymatrix[i]->Get_Solve_Flag()==0&& toymatrix[i]->Get_mcindex()>0){
	 WireCell2dToy::ToyMatrixIterate toymatrix_it(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],2000,1e5,penalty,penalty_ncpt);
	 GeomCellSelection allmcell = mergetiling[i]->get_allcell();
	 CellChargeMap ccmap = truthtiling[i]->ccmap();
	 if (toymatrix[i]->Get_Solve_Flag()!=0)
	   toymetric.Add(allmcell,*toymatrix[i],ccmap);
	 toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
	 cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
	 cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
	 solve_list.erase(solve_list.begin());
       }
     }else{
       for (auto it = solve_list.begin(); it!= solve_list.end(); it++){
	 int i = *it;
	 if (toymatrix[i]->Get_Solve_Flag()==0&& toymatrix[i]->Get_mcindex()>0){
	   if ((toymatrix[i-1]->Get_Solve_Flag() + toymatrix[i+1]->Get_Solve_Flag()) > 0){
	     WireCell2dToy::ToyMatrixIterate toymatrix_it(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],2000,1e5,penalty,penalty_ncpt);
	     GeomCellSelection allmcell = mergetiling[i]->get_allcell();
	     CellChargeMap ccmap = truthtiling[i]->ccmap();
	     if (toymatrix[i]->Get_Solve_Flag()!=0)
	       toymetric.Add(allmcell,*toymatrix[i],ccmap);
	     toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
	     cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
	     cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
	     it = solve_list.erase(it);
	   }
	 }
       }
     }
     
     prev_count = curr_count;
   }


   
   // by the end do the final two
   if (toymatrix[end_num]->Get_Solve_Flag()==0&& toymatrix[end_num]->Get_mcindex()>0){
     WireCell2dToy::ToyMatrixIterate toymatrix_it(toymatrix[end_num-1],toymatrix[end_num],0,mergetiling[end_num-1],mergetiling[end_num],0,2000,1e5,penalty,penalty_ncpt);
     GeomCellSelection allmcell = mergetiling[end_num]->get_allcell();
     CellChargeMap ccmap = truthtiling[end_num]->ccmap();
     if (toymatrix[end_num]->Get_Solve_Flag()!=0)
       toymetric.Add(allmcell,*toymatrix[end_num],ccmap);
     toymetric.AddSolve(toymatrix[end_num]->Get_Solve_Flag());
     
     cout << "chi2: " << end_num << " " << toymatrix[end_num]->Get_Chi2() << endl;
     cout << "NDF: " << toymatrix[end_num]->Get_ndf() << endl;
   }
   
   if (toymatrix[start_num]->Get_Solve_Flag()==0){
     WireCell2dToy::ToyMatrixIterate toymatrix_it(0,toymatrix[start_num],toymatrix[start_num+1],0,mergetiling[start_num],mergetiling[start_num+1],2000,1e5,penalty,penalty_ncpt);
     
     GeomCellSelection allmcell = mergetiling[start_num]->get_allcell();
     CellChargeMap ccmap = truthtiling[start_num]->ccmap();
     if (toymatrix[start_num]->Get_Solve_Flag()!=0&& toymatrix[start_num]->Get_mcindex()>0)
       toymetric.Add(allmcell,*toymatrix[start_num],ccmap);
     toymetric.AddSolve(toymatrix[start_num]->Get_Solve_Flag());
     
     cout << "chi2: " << start_num << " " << toymatrix[start_num]->Get_Chi2() << endl;
     cout << "NDF: " << toymatrix[start_num]->Get_ndf() << endl;
   }

   
   toymetric.Print();
   std::cout << "Starting MCMC" << std::endl;
   solve_list.clear();
   
   // //without  time information
   // for (int i=start_num;i!=end_num+1;i++){
   //   if (toymatrix[i]->Get_Solve_Flag()==0){
   //     GeomCellSelection allmcell = mergetiling[i]->get_allcell();
   //     //WireCell2dToy::ToyMatrixMarkov toymatrix_markov(*toymatrix[i],*mergetiling[i],&allmcell,Good_MCells.at(i-start_num));
   //     WireCell2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i],&allmcell);
   //     CellChargeMap ccmap = truthtiling[i]->ccmap();
   //     if (toymatrix[i]->Get_Solve_Flag()!=0)
   // 	toymetric.Add(allmcell,*toymatrix[i],ccmap);
   //     toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
   //     cout << " chi2: " << i << " " << toymatrix[i]->Get_Chi2() 
   // 	   << " NDF: " << toymatrix[i]->Get_ndf() << endl;
   //   }
   // }
   
   
   //with time information
   if (start_num != end_num){
     int first_solve=-1;
     for (int i=start_num; i!=end_num+1;i++){
       if (toymatrix[i]->Get_Solve_Flag()!=0){
   	 first_solve = i;
   	 break;
       }
     }
     if (first_solve == -1) first_solve = start_num;
     
     
     for (int i=first_solve+1;i<=end_num-1;i++){
       if (toymatrix[i]->Get_Solve_Flag()==0&& toymatrix[i]->Get_mcindex()>0){
   	 if ((toymatrix[i-1]->Get_Solve_Flag() + toymatrix[i+1]->Get_Solve_Flag()) > 0){
   	   GeomCellSelection allmcell = mergetiling[i]->get_allcell();
   	   WireCell2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],&allmcell,1500,2000,penalty,penalty_ncpt);
   	   CellChargeMap ccmap = truthtiling[i]->ccmap();
   	   if (toymatrix[i]->Get_Solve_Flag()!=0)
   	     toymetric.Add(allmcell,*toymatrix[i],ccmap);
   	   toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
	   
   	   cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
   	   cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
   	 }else{
   	   solve_list.push_back(i);
   	 }
   	 //toymetric.Print();
       }
     }
     
     // go to early ones 
     for (int i=first_solve-1;i>=start_num+1;i--){
       if (toymatrix[i]->Get_Solve_Flag()==0&& toymatrix[i]->Get_mcindex()>0){
   	 if ((toymatrix[i-1]->Get_Solve_Flag() + toymatrix[i+1]->Get_Solve_Flag()) > 0){
   	   GeomCellSelection allmcell = mergetiling[i]->get_allcell();
   	   WireCell2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],&allmcell,1500,2000,penalty,penalty_ncpt);
	   
   	   CellChargeMap ccmap = truthtiling[i]->ccmap();
   	   if (toymatrix[i]->Get_Solve_Flag()!=0)
   	     toymetric.Add(allmcell,*toymatrix[i],ccmap);
   	   toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
	   
   	   cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
   	   cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
   	 }else{
   	   solve_list.push_back(i);
   	 }
       }
     }
      
     // do the while ... 
     int prev_count = 0;
   while (solve_list.size() >0){
     int curr_count = solve_list.size();
     std::cout << "Connectivitiy rest " << solve_list.size() << std::endl;

     if (curr_count == prev_count){
       int i = solve_list.front(); // pick the first element ... 
       if (toymatrix[i]->Get_Solve_Flag()==0&& toymatrix[i]->Get_mcindex()>0){
   	 GeomCellSelection allmcell = mergetiling[i]->get_allcell();
   	 WireCell2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],&allmcell,1500,2000,penalty,penalty_ncpt);
	 
   	 CellChargeMap ccmap = truthtiling[i]->ccmap();
   	 if (toymatrix[i]->Get_Solve_Flag()!=0)
   	   toymetric.Add(allmcell,*toymatrix[i],ccmap);
   	 toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
   	 cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
   	 cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
   	 solve_list.erase(solve_list.begin());
       }
     }else{
       for (auto it = solve_list.begin(); it!= solve_list.end(); it++){
   	 int i = *it;
   	 if (toymatrix[i]->Get_Solve_Flag()==0&& toymatrix[i]->Get_mcindex()>0){
   	   if ((toymatrix[i-1]->Get_Solve_Flag() + toymatrix[i+1]->Get_Solve_Flag()) > 0){
   	     GeomCellSelection allmcell = mergetiling[i]->get_allcell();
   	     WireCell2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[i-1],toymatrix[i],toymatrix[i+1],mergetiling[i-1],mergetiling[i],mergetiling[i+1],&allmcell,1500,2000,penalty,penalty_ncpt);
	     
   	     CellChargeMap ccmap = truthtiling[i]->ccmap();
   	     if (toymatrix[i]->Get_Solve_Flag()!=0)
   	       toymetric.Add(allmcell,*toymatrix[i],ccmap);
   	     toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());
   	     cout << "chi2: " << i << " " << toymatrix[i]->Get_Chi2() << endl;
   	     cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
   	     it = solve_list.erase(it);
   	   }
   	 }
       }
     }
     
     prev_count = curr_count;
   }


     //deal with the start/end ones ... 
     if (toymatrix[end_num]->Get_Solve_Flag()==0&& toymatrix[end_num]->Get_mcindex()>0){
       GeomCellSelection allmcell = mergetiling[end_num]->get_allcell();
       WireCell2dToy::ToyMatrixMarkov toymatrix_markov(toymatrix[end_num-1],toymatrix[end_num],0,mergetiling[end_num-1],mergetiling[end_num],0,&allmcell,1500,2000,penalty,penalty_ncpt);
       
       
       CellChargeMap ccmap = truthtiling[end_num]->ccmap();
       if (toymatrix[end_num]->Get_Solve_Flag()!=0)
   	 toymetric.Add(allmcell,*toymatrix[end_num],ccmap);
       toymetric.AddSolve(toymatrix[end_num]->Get_Solve_Flag());
       
       cout << "chi2: " << end_num << " " << toymatrix[end_num]->Get_Chi2() << endl;
       cout << "NDF: " << toymatrix[end_num]->Get_ndf() << endl;
     }
     
     if (toymatrix[start_num]->Get_Solve_Flag()==0&& toymatrix[start_num]->Get_mcindex()>0){
       GeomCellSelection allmcell = mergetiling[start_num]->get_allcell();
       WireCell2dToy::ToyMatrixMarkov toymatrix_markov(0,toymatrix[start_num],toymatrix[start_num+1],0,mergetiling[start_num],mergetiling[start_num+1],&allmcell,1500,2000,penalty,penalty_ncpt);
       
       
       CellChargeMap ccmap = truthtiling[start_num]->ccmap();
       if (toymatrix[start_num]->Get_Solve_Flag()!=0)
   	 toymetric.Add(allmcell,*toymatrix[start_num],ccmap);
       toymetric.AddSolve(toymatrix[start_num]->Get_Solve_Flag());
       
       cout << "chi2: " << start_num << " " << toymatrix[start_num]->Get_Chi2() << endl;
       cout << "NDF: " << toymatrix[start_num]->Get_ndf() << endl;
     }
   }


   // form a map to illustrate connectivities 
   GeomCellCellsMap cell_prev_map;
   GeomCellCellsMap cell_next_map;
   
   for (int i=start_num;i!=end_num;i++){
     GeomCellSelection curr_mcell = mergetiling[i]->get_allcell();
     GeomCellSelection next_mcell = mergetiling[i+1]->get_allcell();
     for (int j=0;j!=curr_mcell.size();j++){
       const MergeGeomCell *curr_cell = (MergeGeomCell*)curr_mcell.at(j);
       for (int k=0;k!=next_mcell.size();k++){
	 const MergeGeomCell *next_cell = (MergeGeomCell*)next_mcell.at(k);
	 if (curr_cell->Overlap(*next_cell) && curr_cell->get_allcell().at(0)->get_face() == next_cell->get_allcell().at(0)->get_face()){
	   if (cell_next_map.find(curr_cell)==cell_next_map.end()){
	     GeomCellSelection cells;
	     cell_next_map[curr_cell] = cells;
	   }
	   cell_next_map[curr_cell].push_back(next_cell);
	   
	   if (cell_prev_map.find(next_cell)==cell_prev_map.end()){
	     GeomCellSelection cells;
	     cell_prev_map[next_cell] = cells;
	   }
	   cell_prev_map[next_cell].push_back(curr_cell);
	   
	 }
       }
     }
   }
   // 

   // save good cluster_cells;
   GeomCellSelection good_cluster_cells;
   for (int i=start_num;i!=end_num+1;i++){
     GeomCellSelection curr_mcell = mergetiling[i]->get_allcell();
     for (int j=0;j!=curr_mcell.size();j++){
       const MergeGeomCell *curr_cell = (MergeGeomCell*)curr_mcell.at(j);
       int flag_good_cluster_cell = 0;
     
       if (cell_next_map.find(curr_cell) != cell_next_map.end()){
	 for (int k=0;k!=cell_next_map[curr_cell].size();k++){
	   if (toymatrix[i+1]->Get_Cell_Charge(cell_next_map[curr_cell].at(k))> recon_threshold || toymatrix[i+1]->Get_Solve_Flag()==0){
	     flag_good_cluster_cell = 1;
	     break;
	   }
	 }
       }
       if (flag_good_cluster_cell == 0){
	 if (cell_prev_map.find(curr_cell) != cell_prev_map.end()){
	   for (int k=0;k!=cell_prev_map[curr_cell].size();k++){
	     if (toymatrix[i-1]->Get_Cell_Charge(cell_prev_map[curr_cell].at(k))> recon_threshold || toymatrix[i-1]->Get_Solve_Flag()==0 ){
	       flag_good_cluster_cell = 1;
	       break;
	     }
	   }
	 }
       }

       if (flag_good_cluster_cell == 1)
	 good_cluster_cells.push_back(curr_cell);
     }
   }
   


   //do clustering ... 
   for (int i=start_num;i!=end_num+1;i++){
     GeomCellSelection pallmcell = mergetiling[i]->get_allcell();
     GeomCellSelection allmcell;
     for (int j=0;j!=pallmcell.size();j++){
       const GeomCell* mcell = pallmcell[j];
       
       // if (toymatrix[i]->Get_Cell_Charge(mcell)> recon_threshold || toymatrix[i]->Get_Solve_Flag()==0){
       // 	 allmcell.push_back(mcell);
       // }
       
        int flag_save_cell = 0;

      if (toymatrix[i]->Get_Solve_Flag()==0){
      	flag_save_cell = 1;
      }else{
      	if (toymatrix[i]->Get_Cell_Charge(mcell)> recon_threshold){
      	  flag_save_cell = 1;
      	}else{
      	  if (i == start_num || i == end_num + 1) continue;
	  
      	  // there are good cells from the prev and next
      	  flag_save_cell = 0;
      	  //	  std::cout << "Xin: " << cell_next_map[mcell].size() <<  " " << cell_prev_map[mcell].size() << std::endl;

      	  for (int k=0;k!=cell_next_map[mcell].size();k++){
	    if (find(good_cluster_cells.begin(),good_cluster_cells.end(),cell_next_map[mcell].at(k))!=good_cluster_cells.end()){
	      //if (toymatrix[i+1]->Get_Cell_Charge(cell_next_map[mcell].at(k))> recon_threshold){
      	      flag_save_cell = 1;
      	      break;
      	    }
      	  }
      	  if (flag_save_cell==1){
      	    flag_save_cell = 0;
      	    for (int k=0;k!=cell_prev_map[mcell].size();k++){
	      if (find(good_cluster_cells.begin(),good_cluster_cells.end(),cell_prev_map[mcell].at(k))!=good_cluster_cells.end()){
		//if (toymatrix[i-1]->Get_Cell_Charge(cell_prev_map[mcell].at(k))> recon_threshold){
      		flag_save_cell = 1;
      		break;
      	      }
      	    }
      	  }
      	}
      }
      
      
      if (flag_save_cell == 1)
      	allmcell.push_back(mcell);
     }
     
     
     
     if (cluster_set.empty()){
       // if cluster is empty, just insert all the mcell, each as a cluster
       for (int j=0;j!=allmcell.size();j++){
  	 GeomCluster *cluster = new GeomCluster(*((MergeGeomCell*)allmcell[j]));
  	 cluster_set.insert(cluster);
       }
     }else{
       for (int j=0;j!=allmcell.size();j++){
  	 int flag = 0;
  	 int flag_save = 0;
  	 GeomCluster *cluster_save = 0;
	 
  	 cluster_delset.clear();
	 
  	 // loop through merged cell
  	 for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
  	   //loop through clusters
	   
  	   flag += (*it)->AddCell(*((MergeGeomCell*)allmcell[j]));
  	   if (flag==1 && flag != flag_save){
  	     cluster_save = *it;
  	   }else if (flag>1 && flag != flag_save){
  	     cluster_save->MergeCluster(*(*it));
  	     cluster_delset.insert(*it);
  	   }
  	   flag_save = flag;
  	   
  	 }
	 
  	 for (auto it = cluster_delset.begin();it!=cluster_delset.end();it++){
  	   cluster_set.erase(*it);
  	   delete (*it);
  	 }
	 
  	 if (flag==0){
  	   GeomCluster *cluster = new GeomCluster(*((MergeGeomCell*)allmcell[j]));
  	   cluster_set.insert(cluster);
  	 }
	 
       }
     }

     int ncount_mcell_cluster = 0;
      for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
  	ncount_mcell_cluster += (*it)->get_allcell().size();
      }
      ncount_mcell += allmcell.size();
      cout << i << " " << allmcell.size()  << " " << cluster_set.size()  << endl;
   }

   int ncount_mcell_cluster = 0;
   for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
     ncount_mcell_cluster += (*it)->get_allcell().size();
   }
   cout << "Summary: " << ncount << " " << ncount_mcell << " " << ncount_mcell_cluster << endl;
   

 
  TFile *file = new TFile(Form("dune_work_space_imaging_%d.root",eve_num),"RECREATE");
  TTree *t_true = new TTree("T_true","T_true");
  TTree *t_rec = new TTree("T_rec","T_rec");
  TTree *t_rec_charge = new TTree("T_rec_charge","T_rec_charge");
  TTree *t_rec_charge_blob = new TTree("T_rec_charge_blob","T_rec_charge_blob");

  Double_t x_save, y_save, z_save;
  Double_t charge_save;
  Double_t ncharge_save;
  Double_t chi2_save;
  Double_t ndf_save;
  Double_t type_save;

  t_true->SetDirectory(file);
  t_true->Branch("x",&x_save,"x/D");
  t_true->Branch("y",&y_save,"y/D");
  t_true->Branch("z",&z_save,"z/D");
  t_true->Branch("q",&charge_save,"q/D");
  
  t_rec->SetDirectory(file);
  t_rec->Branch("x",&x_save,"x/D");
  t_rec->Branch("y",&y_save,"y/D");
  t_rec->Branch("z",&z_save,"z/D");
  
  t_rec_charge->SetDirectory(file);
  t_rec_charge->Branch("x",&x_save,"x/D");
  t_rec_charge->Branch("y",&y_save,"y/D");
  t_rec_charge->Branch("z",&z_save,"z/D");
  t_rec_charge->Branch("q",&charge_save,"q/D");
  t_rec_charge->Branch("nq",&ncharge_save,"nq/D");
  t_rec_charge->Branch("chi2",&chi2_save,"chi2/D");
  t_rec_charge->Branch("ndf",&ndf_save,"ndf/D");

  //blob stuff
  t_rec_charge_blob->SetDirectory(file);
  t_rec_charge_blob->Branch("x",&x_save,"x/D");
  t_rec_charge_blob->Branch("y",&y_save,"y/D");
  t_rec_charge_blob->Branch("z",&z_save,"z/D");
  t_rec_charge_blob->Branch("q",&charge_save,"q/D");
  t_rec_charge_blob->Branch("nq",&ncharge_save,"nq/D");
  
  if (save_image_outline_flag==1){
    t_rec->Branch("type",&type_save,"type/D");
    t_rec_charge->Branch("type",&type_save,"type/D");
    t_rec_charge_blob->Branch("type",&type_save,"type/D");
  }
  
  TGraph2D *g = new TGraph2D();
  TGraph2D *gt = new TGraph2D();
  TGraph2D *g_rec = new TGraph2D();
  TGraph2D *g_rec_blob = new TGraph2D();

  
  
  


   //save results 
  for (int i=start_num;i!=end_num+1;i++){
    //truth
    CellChargeMap ccmap = truthtiling[i]->ccmap();
    for (auto it = ccmap.begin();it!=ccmap.end(); it++){
      Point p = it->first->center();
      
      int cryo = it->first->get_cryo();
      int apa = it->first->get_apa();
      int face = it->first->get_face();
      const WrappedGDS *apa_gds = gds.get_apaGDS(cryo,apa);
      std::pair<double, double> xmm = apa_gds->minmax(0); 
      
      if (face == 1){
	x_save = (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) + xmm.second/units::cm; // *4 is temporary
      }else if (face == 0){
	x_save = xmm.first/units::cm - (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.);
      }
      y_save = p.y/units::cm;
      z_save = p.z/units::cm;
      charge_save = it->second;
      
      gt->SetPoint(ncount_t,x_save,y_save,z_save);
      t_true->Fill();
            
      ncount_t ++;
    }

     //recon 1
    GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    for (int j=0;j!=allmcell.size();j++){
      MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];

      if (save_image_outline_flag==1){
	Point p = mcell->get_allcell().at(0)->center();
	int cryo = mcell->get_allcell().at(0)->get_cryo();
	int apa = mcell->get_allcell().at(0)->get_apa();
	int face = mcell->get_allcell().at(0)->get_face();
	const WrappedGDS *apa_gds = gds.get_apaGDS(cryo,apa);
	std::pair<double, double> xmm = apa_gds->minmax(0); 
	
	if (face == 1){
	  x_save = (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) + xmm.second/units::cm; // *4 is temporary
	}else if (face == 0){
	  x_save = xmm.first/units::cm - (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.);
	}
	y_save = mcell->center().y/units::cm;
	z_save = mcell->center().z/units::cm;
	
	g->SetPoint(ncount,x_save,y_save,z_save);
	ncount ++;
	type_save = 1; //center
	t_rec->Fill();
	
	for (int k=0;k!=mcell->get_edgecells().size();k++){
	  Point p = mcell->get_edgecells().at(k)->center();
	  y_save = p.y/units::cm;
    	  z_save = p.z/units::cm;
	  g->SetPoint(ncount,x_save,y_save,z_save);
	  ncount ++;
	  type_save = 2; //boundary ...
	  t_rec->Fill();
	}
	
	for (int k=0;k!=mcell->get_allcell().size();k++){
    	  Point p = mcell->get_allcell().at(k)->center();
	  
	  
	  int cryo = mcell->get_allcell().at(k)->get_cryo();
	  int apa = mcell->get_allcell().at(k)->get_apa();
	  int face = mcell->get_allcell().at(k)->get_face();
	  const WrappedGDS *apa_gds = gds.get_apaGDS(cryo,apa);
	  std::pair<double, double> xmm = apa_gds->minmax(0); 
	  
	  if (face == 1){
	    x_save = (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) + xmm.second/units::cm; // *4 is temporary
	  }else if (face == 0){
	    x_save = xmm.first/units::cm - (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.);
	  }
	  
    	  y_save = p.y/units::cm;
    	  z_save = p.z/units::cm;
	  
	  type_save = 0;
	  t_rec->Fill();
	}
      }else{
	for (int k=0;k!=mcell->get_allcell().size();k++){
    	  Point p = mcell->get_allcell().at(k)->center();
	  
	  
	  int cryo = mcell->get_allcell().at(k)->get_cryo();
	  int apa = mcell->get_allcell().at(k)->get_apa();
	  int face = mcell->get_allcell().at(k)->get_face();
	  const WrappedGDS *apa_gds = gds.get_apaGDS(cryo,apa);
	  std::pair<double, double> xmm = apa_gds->minmax(0); 
	  
	  if (face == 1){
	    x_save = (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) + xmm.second/units::cm; // *4 is temporary
	  }else if (face == 0){
	    x_save = xmm.first/units::cm - (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.);
	  }
	  
    	  y_save = p.y/units::cm;
    	  z_save = p.z/units::cm;
	  
	  // int fill_flag = 1;
	  // if (pre_scale_flag ==1 ){
	  //   if (mcell->get_allcell().size()>pre_scale_ncells){
	  //     fill_flag = 0;
	  //   }
	  // }

	  // if (fill_flag==1){
	  g->SetPoint(ncount,x_save,y_save,z_save);
	  ncount ++;
	  t_rec->Fill();
	  // }
	}
      }
    }
    // GeomCellSelection allcell = toytiling[i]->get_allcell();
    // for (int j=0;j!=allcell.size();j++){
    //   Point p = allcell[j]->center();

    //   int cryo = allcell[j]->get_cryo();
    //   int apa = allcell[j]->get_apa();
    //   int face = allcell[j]->get_face();
    //   const WrappedGDS *apa_gds = gds.get_apaGDS(cryo,apa);
    //   std::pair<double, double> xmm = apa_gds->minmax(0); 
      
    //   if (face == 1){
    // 	x_save = (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) + xmm.second/units::cm; // *4 is temporary
    //   }else if (face == 0){
    // 	x_save = xmm.first/units::cm - (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.);
    //   }

    //   //      x_save = i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.*4;
    //   y_save = p.y/units::cm;
    //   z_save = p.z/units::cm;
      

    //   g->SetPoint(ncount,x_save,y_save,z_save);
    //   t_rec->Fill();

    //   ncount ++;
    // }

      //recon 2 with charge
    // GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    for (int j=0;j!=allmcell.size();j++){
      MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
      double charge = toymatrix[i]->Get_Cell_Charge(mcell,1);
      //if (charge> recon_threshold || toymatrix[i]->Get_Solve_Flag()==0){
      if (charge> recon_threshold ){

    	if (toymatrix[i]->Get_Solve_Flag()==0)
	  charge = toytiling[i]->get_ave_charge();

	if (save_image_outline_flag==1){
	  
	  Point p = mcell->get_allcell().at(0)->center();
	  int cryo = mcell->get_allcell().at(0)->get_cryo();
	  int apa = mcell->get_allcell().at(0)->get_apa();
	  int face = mcell->get_allcell().at(0)->get_face();
	  const WrappedGDS *apa_gds = gds.get_apaGDS(cryo,apa);
	  std::pair<double, double> xmm = apa_gds->minmax(0); 
	  if (face == 1){
	    x_save = (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) + xmm.second/units::cm; // *4 is temporary
	  }else if (face == 0){
	    x_save = xmm.first/units::cm - (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.);
	  }
	  y_save = mcell->center().y/units::cm;
	  z_save = mcell->center().z/units::cm;
	  charge_save = charge;
	  ncharge_save = mcell->get_allcell().size();
	  chi2_save = toymatrix[i]->Get_Chi2();
	  ndf_save = toymatrix[i]->Get_ndf();
	  g_rec->SetPoint(ncount1,x_save,y_save,z_save);
	  ncount1 ++;
	  type_save = 1; //center
	  t_rec_charge->Fill();

	  charge_save = 0;
	  ncharge_save = 0;
	  chi2_save= 0;
	  ndf_save = 0;
	  
	  for (int k=0;k!=mcell->get_edgecells().size();k++){
	    Point p = mcell->get_edgecells().at(k)->center();
	    y_save = p.y/units::cm;
	    z_save = p.z/units::cm;
	    g_rec->SetPoint(ncount1,x_save,y_save,z_save);
	    ncount1 ++;
	    type_save = 2; //boundary ...
	    t_rec_charge->Fill();
	  }

	   for (int k=0;k!=mcell->get_allcell().size();k++){
	    Point p = mcell->get_allcell().at(k)->center();
	    
	    
	    int cryo = mcell->get_allcell().at(k)->get_cryo();
	    int apa = mcell->get_allcell().at(k)->get_apa();
	    int face = mcell->get_allcell().at(k)->get_face();
	    const WrappedGDS *apa_gds = gds.get_apaGDS(cryo,apa);
	    std::pair<double, double> xmm = apa_gds->minmax(0); 
	    
	    if (face == 1){
	      x_save = (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) + xmm.second/units::cm; // *4 is temporary
	    }else if (face == 0){
	      x_save = xmm.first/units::cm - (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.);
	    }
	    
	    y_save = p.y/units::cm;
	    z_save = p.z/units::cm;
	    charge_save = charge/mcell->get_allcell().size();
	    ncharge_save = mcell->get_allcell().size();
	    chi2_save = toymatrix[i]->Get_Chi2();
	    ndf_save = toymatrix[i]->Get_ndf();
	    type_save = 0;
	    t_rec_charge->Fill();
	    
	   }

	}else{
	  for (int k=0;k!=mcell->get_allcell().size();k++){
	    Point p = mcell->get_allcell().at(k)->center();
	    
	    
	    int cryo = mcell->get_allcell().at(k)->get_cryo();
	    int apa = mcell->get_allcell().at(k)->get_apa();
	    int face = mcell->get_allcell().at(k)->get_face();
	    const WrappedGDS *apa_gds = gds.get_apaGDS(cryo,apa);
	    std::pair<double, double> xmm = apa_gds->minmax(0); 
	    
	    if (face == 1){
	      x_save = (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) + xmm.second/units::cm; // *4 is temporary
	    }else if (face == 0){
	      x_save = xmm.first/units::cm - (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.);
	    }
	    
	    y_save = p.y/units::cm;
	    z_save = p.z/units::cm;
	    charge_save = charge/mcell->get_allcell().size();
	    ncharge_save = mcell->get_allcell().size();
	    chi2_save = toymatrix[i]->Get_Chi2();
	    ndf_save = toymatrix[i]->Get_ndf();
	    
	    g_rec->SetPoint(ncount1,x_save,y_save,z_save);
	    t_rec_charge->Fill();
	    
	    ncount1 ++;
	  }
	}
      }
    }
    
    // for (int j=0;j!=allmcell.size();j++){
    //   MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
    //   double charge = toymatrix[i]->Get_Cell_Charge(mcell,1);
    //   if (charge> recon_threshold || toymatrix[i]->Get_Solve_Flag()==0){
    // 	if (toymatrix[i]->Get_Solve_Flag()==0)
    // 	  charge = toytiling[i]->get_ave_charge();

    // 	for (int k=0;k!=mcell->get_allcell().size();k++){
    // 	  Point p = mcell->get_allcell().at(k)->center();

    // 	   int cryo = mcell->get_allcell().at(k)->get_cryo();
    // 	  int apa = mcell->get_allcell().at(k)->get_apa();
    // 	  int face = mcell->get_allcell().at(k)->get_face();
    // 	  const WrappedGDS *apa_gds = gds.get_apaGDS(cryo,apa);
    // 	  std::pair<double, double> xmm = apa_gds->minmax(0); 
	  
    // 	  if (face == 1){
    // 	    x_save = (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) + xmm.second/units::cm; // *4 is temporary
    // 	  }else if (face == 0){
    // 	    x_save = xmm.first/units::cm - (i*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.);
    // 	  }

    	  
    // 	  y_save = p.y/units::cm;
    // 	  z_save = p.z/units::cm;
    // 	  charge_save = charge/mcell->get_allcell().size();
    // 	  ncharge_save = mcell->get_allcell().size();
	  
    // 	  g_rec_blob->SetPoint(ncount2,x_save,y_save,z_save);
    // 	  t_rec_charge_blob->Fill();
	  
    // 	  ncount2 ++;
    // 	}
    //   }
    // }

  }

 

  
   TTree *ttree1 = new TTree("TC","TC");
  // To save cluster, we need to save
  // 1. time slice
  // 2. single cell
  // 3. charge
  // 4. cluster number
  const GeomCell* cell_save = 0;
  int cluster_num = -1;
  int mcell_id = -1;
  int time_slice;
  
  ttree1->Branch("time_slice",&time_slice,"time_slice/I"); // done
  ttree1->Branch("cell",&cell_save);
  ttree1->Branch("ncluster",&cluster_num,"cluster_num/I"); //done
  ttree1->Branch("mcell_id",&mcell_id,"mcell_id/I");
  ttree1->Branch("charge",&charge_save,"charge/D"); 
  
  double xx,yy,zz;

  ttree1->Branch("xx",&xx,"xx/D");    //don
  ttree1->Branch("yy",&yy,"yy/D");    //don
  ttree1->Branch("zz",&zz,"zz/D");    //don
  // ttree1->Branch("x",&x,"x/D");    //done
  // ttree1->Branch("y",&y,"y/D");
  // ttree1->Branch("z",&z,"z/D");

  // save information to reconstruct the toytiling
  int u_index, v_index, w_index;
  double u_charge, v_charge, w_charge;
  double u_charge_err, v_charge_err, w_charge_err;
  
  int apa_no=0, cryostat_no=0;
  int face = 0;
  ttree1->Branch("face",&face,"face/I");
  ttree1->Branch("apa_no",&apa_no,"apa_no/I");
  ttree1->Branch("cryostat_no",&cryostat_no,"cryostat_no/I");

  ttree1->Branch("u_index",&u_index,"u_index/I");
  ttree1->Branch("v_index",&v_index,"v_index/I");
  ttree1->Branch("w_index",&w_index,"w_index/I");
  
  ttree1->Branch("u_charge",&u_charge,"u_charge/D");
  ttree1->Branch("v_charge",&v_charge,"v_charge/D");
  ttree1->Branch("w_charge",&w_charge,"w_charge/D");

  ttree1->Branch("u_charge_err",&u_charge_err,"u_charge_err/D");
  ttree1->Branch("v_charge_err",&v_charge_err,"v_charge_err/D");
  ttree1->Branch("w_charge_err",&w_charge_err,"w_charge_err/D");

  //end save 

  ttree1->SetDirectory(file);
  
  for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
    cluster_num ++;
    //loop merged cell
    for (int i=0; i!=(*it)->get_allcell().size();i++){
      const MergeGeomCell *mcell = (const MergeGeomCell*)((*it)->get_allcell().at(i));
      mcell_id ++;
      time_slice = mcell->GetTimeSlice();

      //      
      //loop single cell
      for (int j=0; j!=mcell->get_allcell().size();j++){
	cell_save = mcell->get_allcell().at(j);

	//fill the information needed for toytiling
	GeomWireSelection wires = toytiling[time_slice]->wires(*cell_save);
	
	//	if (i==0 && j==0) cout << "abc: " << time_slice << " " << toytiling[time_slice]->get_allcell().size() << " " << wires.size() << endl;

	for (int k=0;k!=wires.size();k++){
	  const GeomWire *wire = wires.at(k);
	  WirePlaneType_t plane = wire->plane();
	  if (plane==0){
	    u_index = wire->index();
	    u_charge = toytiling[time_slice]->wcmap()[wire];
	    u_charge_err = toytiling[time_slice]->wcemap()[wire];
	  }else if (plane==1){
	    v_index = wires.at(k)->index();
	    v_charge = toytiling[time_slice]->wcmap()[wire];
	    v_charge_err = toytiling[time_slice]->wcemap()[wire];
	  }else if (plane==2){
	    w_index = wire->index();
	    w_charge = toytiling[time_slice]->wcmap()[wire];
	    w_charge_err = toytiling[time_slice]->wcemap()[wire];
	  }
	}

	//end fill

	Point p = mcell->get_allcell().at(j)->center();

	int cryo = mcell->get_allcell().at(j)->get_cryo();
	apa_no = mcell->get_allcell().at(j)->get_apa();
	face = mcell->get_allcell().at(j)->get_face();
	const WrappedGDS *apa_gds = gds.get_apaGDS(cryo,apa_no);
	std::pair<double, double> xmm = apa_gds->minmax(0); 
	
	if (face == 1){
	  xx = (time_slice*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) + xmm.second/units::cm; // *4 is temporary
	}else if (face == 0){
	  xx = xmm.first/units::cm - (time_slice*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.);
	}
	//xx = time_slice*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.;
	
	charge_save = toymatrix[time_slice]->Get_Cell_Charge(mcell,1)/mcell->cross_section() * cell_save->cross_section();
	yy = p.y/units::cm;
  	zz = p.z/units::cm;
	ttree1->Fill();
	
	if (save_image_outline_flag==0){
	  //save the g_rec_blob tree ... 
	  x_save = xx;
	  y_save = yy;
	  z_save = zz;
	  ncharge_save = mcell->get_allcell().size();
	  g_rec_blob->SetPoint(ncount2,x_save,y_save,z_save);
	  ncount2++;
	  t_rec_charge_blob->Fill();
	}else if (save_image_outline_flag==1){
	  //save the g_rec_blob tree ... 
	  x_save = xx;
	  y_save = yy;
	  z_save = zz;
	  ncharge_save = mcell->get_allcell().size();
	  type_save = 0;
	  t_rec_charge_blob->Fill();
	}
      }
      
      if (save_image_outline_flag==1){
	Point p = mcell->get_allcell().at(0)->center();

	int cryo = mcell->get_allcell().at(0)->get_cryo();
	apa_no = mcell->get_allcell().at(0)->get_apa();
	face = mcell->get_allcell().at(0)->get_face();
	const WrappedGDS *apa_gds = gds.get_apaGDS(cryo,apa_no);
	std::pair<double, double> xmm = apa_gds->minmax(0); 
	if (face == 1){
	  xx = (time_slice*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.) + xmm.second/units::cm; // *4 is temporary
	}else if (face == 0){
	  xx = xmm.first/units::cm - (time_slice*nrebin/2.*unit_dis/10. - frame_length/2.*unit_dis/10.);
	}
	
	x_save = xx;
	y_save = mcell->center().y/units::cm;
	z_save = mcell->center().z/units::cm;
	charge_save = toymatrix[time_slice]->Get_Cell_Charge(mcell,1);
	ncharge_save = mcell->get_allcell().size();
	g_rec_blob->SetPoint(ncount2,x_save,y_save,z_save);
	ncount2++;
	type_save = 1;
	t_rec_charge_blob->Fill();

	charge_save = 0;
	ncharge_save = 0;
	
	for (int k=0;k!=mcell->get_edgecells().size();k++){
	  Point p = mcell->get_edgecells().at(k)->center();
	  y_save = p.y/units::cm;
	  z_save = p.z/units::cm;
	  g_rec_blob->SetPoint(ncount2,x_save,y_save,z_save);
	  ncount2 ++;
	  type_save = 2; //boundary ...
	  t_rec_charge_blob->Fill();
	}

      }

    }
  }
  ttree1->Write();
  g->Write("shower3D");
  gt->Write("shower3D_truth");
  g_rec->Write("shower3D_charge");
  g_rec_blob->Write("shower3D_charge_blob");

  TTree *Trun = new TTree("Trun","Trun");
  Trun->SetDirectory(file);

  int detector = 3; // far detector dune work space
  Trun->Branch("detector",&detector,"detector/I");

  Trun->Branch("eventNo",&event_no,"eventNo/I");
  Trun->Branch("runNo",&run_no,"runNo/I");
  Trun->Branch("subRunNo",&subrun_no,"runRunNo/I");

  Trun->Branch("unit_dis",&unit_dis,"unit_dis/F");
  Trun->Branch("toffset_uv",&toffset_1,"toffset_uv/F");
  Trun->Branch("toffset_uw",&toffset_2,"toffset_uw/F");
  Trun->Branch("toffset_u",&toffset_3,"toffset_u/F");
  Trun->Branch("total_time_bin",&total_time_bin,"total_time_bin/I");
  Trun->Branch("recon_threshold",&recon_threshold,"recon_threshold/I");
  Trun->Branch("frame_length",&frame_length,"frame_length/I");
  Trun->Branch("max_events",&max_events,"max_events/I");
  Trun->Branch("eve_num",&eve_num,"eve_num/I");
  Trun->Branch("nrebin",&nrebin,"nrebin/I");
  Trun->Branch("threshold_u",&threshold_u,"threshold_u/F");
  Trun->Branch("threshold_v",&threshold_v,"threshold_v/F");
  Trun->Branch("threshold_w",&threshold_w,"threshold_w/F");
  Trun->Branch("time_offset",&time_offset,"time_offset/I");

 pitch_u = pitch_u/units::cm;
  pitch_v = pitch_v/units::cm;
  pitch_w = pitch_w/units::cm;
  Trun->Branch("pitch_u",&pitch_u,"pitch_u/D");
  Trun->Branch("pitch_v",&pitch_v,"pitch_v/D");
  Trun->Branch("pitch_w",&pitch_w,"pitch_w/D");


  Trun->Fill();

  TTree *TMC = new TTree("TMC","TMC");
  TMC->SetDirectory(file);
  Float_t mc_oldVertex[3];
  mc_oldVertex[0] = vertex.x - x_shift;
  mc_oldVertex[1] = vertex.y - y_shift;
  mc_oldVertex[2] = vertex.z - z_shift;
  
  TMC->Branch("mc_Ntrack", &mctruth->mc_Ntrack);  // number of tracks in MC
  TMC->Branch("mc_id", &mctruth->mc_id, "mc_id[mc_Ntrack]/I");  // track id; size == mc_Ntrack
  TMC->Branch("mc_pdg", &mctruth->mc_pdg, "mc_id[mc_Ntrack]/I");  // track particle pdg; size == mc_Ntrack
  TMC->Branch("mc_process", &mctruth->mc_process, "mc_process[mc_Ntrack]/I");  // track generation process code; size == mc_Ntrack
  TMC->Branch("mc_mother", &mctruth->mc_mother, "mc_mother[mc_Ntrack]/I");  // mother id of this track; size == mc_Ntrack
  TMC->Branch("mc_daughters", mctruth->mc_daughters);  // daughters id of this track; vector
  TMC->Branch("mc_startXYZT", &mctruth->mc_startXYZT, "mc_startXYZT[mc_Ntrack][4]/F");  // start position of this track; size == mc_Ntrack
  TMC->Branch("mc_endXYZT", &mctruth->mc_endXYZT, "mc_endXYZT[mc_Ntrack][4]/F");  // start position of this track; size == mc_Ntrack
  TMC->Branch("mc_startMomentum", &mctruth->mc_startMomentum, "mc_startMomentum[mc_Ntrack][4]/F");  // start momentum of this track; size == mc_Ntrack
  TMC->Branch("mc_endMomentum", &mctruth->mc_endMomentum, "mc_endMomentum[mc_Ntrack][4]/F");  // start momentum of this track; size == mc_Ntrack
  if (mctruth->mc_trackPosition != NULL) TMC->Branch("mc_trackPosition",mctruth->mc_trackPosition);
  
  
  TMC->Branch("mc_isnu", &mctruth->mc_isnu);
  TMC->Branch("mc_nGeniePrimaries", &mctruth->mc_nGeniePrimaries);
  TMC->Branch("mc_nu_pdg", &mctruth->mc_nu_pdg);
  TMC->Branch("mc_nu_ccnc", &mctruth->mc_nu_ccnc);
  TMC->Branch("mc_nu_mode", &mctruth->mc_nu_mode);
  TMC->Branch("mc_nu_intType", &mctruth->mc_nu_intType);
  TMC->Branch("mc_nu_target", &mctruth->mc_nu_target);
  TMC->Branch("mc_hitnuc", &mctruth->mc_hitnuc);
  TMC->Branch("mc_hitquark", &mctruth->mc_hitquark);
  TMC->Branch("mc_nu_Q2", &mctruth->mc_nu_Q2);
  TMC->Branch("mc_nu_W", &mctruth->mc_nu_W);
  TMC->Branch("mc_nu_X", &mctruth->mc_nu_X);
  TMC->Branch("mc_nu_Y", &mctruth->mc_nu_Y);
  TMC->Branch("mc_nu_Pt", &mctruth->mc_nu_Pt);
  TMC->Branch("mc_nu_Theta", &mctruth->mc_nu_Theta);
  TMC->Branch("mc_nu_pos", &mctruth->mc_nu_pos, "mc_nu_pos[4]/F");
  TMC->Branch("mc_nu_mom", &mctruth->mc_nu_mom, "mc_nu_mom[4]/F");



  TMC->Branch("mc_oldVertex", &mc_oldVertex, "mc_oldVertex[3]/F");
  TMC->Fill();


  file->Write();
  file->Close();

  return 0;

  // TCanvas *c = new TCanvas();
  // c->Range(-5*units::cm, -90*units::cm, 160*units::cm, 120*units::cm);    
  
  // TLine *l = new TLine();
  // l->SetLineWidth(0);
  // c->cd();
  
  
  
  // int colors[] = {2,4,1};
  // for (short cryo = 0; cryo < gds.ncryos(); cryo++) {
  //   for (short apa = 0; apa < gds.napa(cryo); apa++) {
  // 	std::cout << cryo << " " << apa << std::endl;
  
  // 	    const WrappedGDS *apa_gds = gds.get_apaGDS(cryo, apa);
  // 	    for (int iplane=0; iplane<3; ++iplane) {
  // 	        WirePlaneType_t plane = (WirePlaneType_t)iplane;
  // 		GeomWireSelection wip = apa_gds->wires_in_plane(plane);
  // 	        // std::cout<<"\n[CRYO] "<<cryo<<" [APA] "<<apa<<" [PLANE] "<<iplane
  // 		// 	 <<" has "<< wip.size()<<" wires, wire angle is "<<apa_gds->angle(plane)*180/TMath::Pi()<<std::endl;
  // 		//for (auto wit = wip.begin(); wit != wip.end(); ++wit) {
  // 		//  const GeomWire& wire = **wit;
  // 		for (int index=0; index<(int)wip.size(); ++index) {
  // 		    const GeomWire* wire = apa_gds->by_planeindex(plane, index);
  // 		    //if (wire.face() == 0) continue;
  // 		    const Vector& p1 = wire->point1();
  // 		    const Vector& p2 = wire->point2();
  // 		    //  std::cout<<*wire<<" ("<<p1.x<<","<<p1.y<<","<<p1.z<<") ("<<p2.x<<","<<p2.y<<","<<p2.z<<")\n";
  // 		    l->SetLineColor(colors[iplane]);
  // 		    l->DrawLine(p1.z, p1.y, p2.z, p2.y);
  // 		}
  // 	    }	    
  // 	}
  // }
  
  // c->SaveAs("./test_detectorgds_35t.pdf");
  }
