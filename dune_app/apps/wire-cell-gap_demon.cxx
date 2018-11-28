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
#include "WireCell2dToy/DataSignalGaus.h"
#include "WireCell2dToy/DataSignalWien.h"

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
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root eve_num -t[0,1] -s[0,1,2] -n[time_slice]" << endl;
    return 1;
  }

  int two_plane = 0;
  int save_file = 0;
  int time_slice = 0;
  
  for(Int_t i = 1; i != argc; i++){
     switch(argv[i][1]){
     case 't':
       two_plane = atoi(&argv[i][2]); 
       break;
     case 's':
       save_file = atoi(&argv[i][2]); 
       break;
     case 'n':
       time_slice = atoi(&argv[i][2]);
       break;
     }
  }
  
  if (two_plane)
    cout << "Enable Two Plane Reconstruction " << endl; 
  
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
  
  
  
  
  // WireCell::FrameDataSource* fds = 0;
  // fds = WireCellSst::make_fds(root_file);
  // if (!fds) {
  //   cerr << "ERROR: failed to get FDS from " << root_file << endl;
  //   return 1;
  // }
  
  //float unit_dis = 1.01483;  // 58KV @ 226.5 V/cm
  float unit_dis = 1.14753;  // 70 KV @ 226.5 V/cm

  //**** time offset for 58kV ****// 
  //float toffset_1=-(1.834-1.647) -0.1;//+ 0.1 + (nt_off1 * 0.1 - 0.5 );  // time offset between u/v 
  //float toffset_2=-(1.834+1.555-1.539-1.647) +0.1;//+ 0.3 + (nt_off2*0.2 - 1); // time offset between u/w
  //float toffset_3=-0.5; //overall time shift

  //final offset after time scan (70kV)
  float toffset_1=-(1.834-1.647) -0.1 - 0.5;
  float toffset_2=-(1.834+1.555-1.539-1.647) +0.1 - 0.5;
  float toffset_3=0.0;
  

  int total_time_bin=9594;
  int recon_threshold = 2000;
  int frame_length = 3200;
  int max_events = 100;
  int eve_num  = atoi(argv[3]);
  int nrebin = 6;
  float threshold_u = 5.87819e+02 * 4.0;
  float threshold_v = 8.36644e+02 * 4.0;
  float threshold_w = 5.67974e+02 * 4.0;

  float threshold_ug = 755.96;
  float threshold_vg = 822.81;
  float threshold_wg = 510.84;
  
  int time_offset = -92.;
  


  const char* root_file = argv[2];
 
  
  int run_no, subrun_no, event_no;
  // sst->SetBranchAddress("eventNo",&event_no);
  // sst->SetBranchAddress("runNo",&run_no);
  // sst->SetBranchAddress("subRunNo",&subrun_no);
  // sst->GetEntry(eve_num);
  
  
  WireCellSst::DatauBooNEFrameDataSource data_fds(root_file,gds,total_time_bin);
  if (save_file != 2){
    data_fds.jump(eve_num);
    if (save_file == 1)
      data_fds.Save();
  }
  
  run_no = data_fds.get_run_no();
  subrun_no = data_fds.get_subrun_no();
  event_no = data_fds.get_event_no();
  
  cout << "Run No: " << run_no << " " << subrun_no << " " << eve_num << endl;


  // WireMap& uplane_map = data_fds.get_u_map();
  // WireMap& vplane_map = data_fds.get_v_map();
  // WireMap& wplane_map = data_fds.get_w_map();

  ChirpMap& uplane_map = data_fds.get_u_cmap();
  ChirpMap& vplane_map = data_fds.get_v_cmap();
  ChirpMap& wplane_map = data_fds.get_w_cmap();

  // std::cout << uplane_map.size() << " " << vplane_map.size() << " " << wplane_map.size() << std::endl;
    
  cout << "Deconvolution with Gaussian filter" << endl;
  WireCell2dToy::DataSignalGausFDS gaus_fds(data_fds,gds,total_time_bin/nrebin,max_events,toffset_1,toffset_2,toffset_3); // gaussian smearing for charge estimation
  if (save_file != 2){
    gaus_fds.jump(eve_num);
    if (save_file == 1)
      gaus_fds.Save();
  }else{
  
  }

  cout << "Deconvolution with Wiener filter" << endl; 
  WireCell2dToy::DataSignalWienFDS wien_fds(data_fds,gds,uplane_map, vplane_map, wplane_map, total_time_bin/nrebin,max_events,toffset_1,toffset_2,toffset_3); // weiner smearing for hit identification
  if (save_file !=2 ){
    wien_fds.jump(eve_num);
    if (save_file == 1)
      wien_fds.Save();
  }else{
    
  }


  data_fds.Clear();

  std::vector<float>& uplane_rms = wien_fds.get_uplane_rms();
  std::vector<float>& vplane_rms = wien_fds.get_vplane_rms();
  std::vector<float>& wplane_rms = wien_fds.get_wplane_rms();

  // // hack for now ...  remove the very busy wires ... 
  // for (int i=0;i!=uplane_rms.size();i++){
  //   //cout << "U " << i << " " << uplane_rms.at(i) << endl;
  //   if (uplane_rms.at(i) > 1500) {
  //     uplane_rms.at(i) *=2;
  //     uplane_map.erase(i);
  //   }
  // }
  // for (int i=0;i!=vplane_rms.size();i++){
  //   //cout << "V " << i << " " << vplane_rms.at(i) << endl;
  //   if (vplane_rms.at(i) > 2000 && vplane_rms.at(i)<3000){
  //     vplane_rms.at(i) *=2;
  //     vplane_map.erase(i);
  //   }else if (vplane_rms.at(i)>=3000){
  //     vplane_rms.at(i) *=10;
  //     vplane_map.erase(i);
  //   }
  // }
  // for (int i=0;i!=wplane_rms.size();i++){
  //   //cout << "W " << i << " " << wplane_rms.at(i) << endl;
  //   if (wplane_rms.at(i) > 1000) {
  //     wplane_rms.at(i) *=2;
  //     wplane_map.erase(i);
  //   }
  // }
  

  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  int nwire_u = wires_u.size();
  int nwire_v = wires_v.size();
  int nwire_w = wires_w.size();
  
  
 
  

 

  // WireCellSst::ToyuBooNESliceDataSource sds(wien_fds,gaus_fds,threshold_u, 
  // 					    threshold_v, threshold_w, 
  // 					    threshold_ug, 
  // 					    threshold_vg, threshold_wg, 
  // 					    nwire_u, 
  // 					    nwire_v, nwire_w); 

  WireCellSst::ToyuBooNESliceDataSource sds(wien_fds,gaus_fds,threshold_u, 
  					    threshold_v, threshold_w, 
  					    threshold_ug, 
  					    threshold_vg, threshold_wg, 
  					    nwire_u, 
  					    nwire_v, nwire_w,
  					    &uplane_rms, &vplane_rms, &wplane_rms); 
    
  
  
  int ncount = 0;
  int ncount1 = 0;
  int ncount2 = 0;

  int ncount_t = 0;
  

  WireCell2dToy::ToyTiling **toytiling = new WireCell2dToy::ToyTiling*[2400];
  WireCell2dToy::BadTiling **badtiling = new WireCell2dToy::BadTiling*[2400];
  WireCell2dToy::MergeToyTiling **mergetiling = new WireCell2dToy::MergeToyTiling*[2400];
  WireCell2dToy::ToyMatrix **toymatrix = new WireCell2dToy::ToyMatrix*[2400];
    
  //add in cluster
  GeomClusterSet cluster_set, cluster_delset;
  
  int ncount_mcell = 0;

  //  delete fds;

  // int start_num = 0 ;
  // int end_num = sds.size()-1;

  int start_num = time_slice;
  int end_num = time_slice;

  for (int i=start_num;i!=end_num+1;i++){
 
    sds.jump(i);
    WireCell::Slice slice = sds.get();

    //cout << i << " " << slice.group().size() << std::endl;
    toytiling[i] = new WireCell2dToy::ToyTiling(slice,gds,0.15,0.2,0.1,threshold_ug,threshold_vg, threshold_wg, &uplane_rms, &vplane_rms, &wplane_rms);

    if (two_plane)
      toytiling[i]->twoplane_tiling(i,nrebin,gds,uplane_rms,vplane_rms,wplane_rms, uplane_map, vplane_map, wplane_map);


    GeomCellSelection allcell = toytiling[i]->get_allcell();
    GeomWireSelection allwire = toytiling[i]->get_allwire();

    cout << i << " " << allcell.size() << " " << allwire.size() << endl;

    mergetiling[i] = new WireCell2dToy::MergeToyTiling(*toytiling[i],i,3);
    
    //    if (two_plane) 
    // mergetiling[i]->deghost();
   
    GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    GeomWireSelection allmwire = mergetiling[i]->get_allwire();
    
    cout << i << " " << allmcell.size() << " " << allmwire.size() << endl;
    

    toymatrix[i] = new WireCell2dToy::ToyMatrix(*toytiling[i],*mergetiling[i]);
    
    GeomCellSelection& two_wires_cells = mergetiling[i]->get_two_wires_cells();
    GeomCellSelection& three_wires_cells = mergetiling[i]->get_three_wires_cells();

 
    
    //comebine them 
    
	
    
    
    cout << "chi2: " << toymatrix[i]->Get_Chi2() << endl;
    cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
    
    badtiling[i] = new WireCell2dToy::BadTiling(i,nrebin,uplane_map,vplane_map,wplane_map,gds);

    // if (toymatrix[i]->Get_Solve_Flag()!=0){
    //   toymatrix[i]->Update_pred();
    //   toymatrix[i]->Print();
    // }

    //draw ... 
    TApplication theApp("theApp",&argc,argv);
    theApp.SetReturnFromRun(true);
    
    TCanvas c1("ToyMC","ToyMC",800,600);
    c1.Draw();
    
    WireCell2dToy::ToyEventDisplay display(c1, gds);
    display.charge_min = 0;
    display.charge_max = 5e4;


    gStyle->SetOptStat(0);
    
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Int_t MyPalette[NCont];
    Double_t stops[NRGBs] = {0.0, 0.34, 0.61, 0.84, 1.0};
    Double_t red[NRGBs] = {0.0, 0.0, 0.87 ,1.0, 0.51};
    Double_t green[NRGBs] = {0.0, 0.81, 1.0, 0.2 ,0.0};
    Double_t blue[NRGBs] = {0.51, 1.0, 0.12, 0.0, 0.0};
    Int_t FI = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    for (int kk=0;kk!=NCont;kk++) MyPalette[kk] = FI+kk;
    gStyle->SetPalette(NCont,MyPalette);

    

    display.init(0,10.3698,-2.33/2.,2.33/2.);
    display.draw_mc(1,WireCell::PointValueVector(),"colz");
    
    // display.draw_bad_region(uplane_map,i,nrebin,0,"same");
    // display.draw_bad_region(vplane_map,i,nrebin,1,"same");
    // display.draw_bad_region(wplane_map,i,nrebin,2,"same");
    display.draw_bad_cell(badtiling[i]->get_cell_all());
  
    double bad_area = 0;
    for (int k = 0; k !=  badtiling[i]->get_cell_all().size(); k++){
      bad_area += badtiling[i]->get_cell_all().at(k)->cross_section();
    }
    std::cout << "Bad Area: " << bad_area/units::m/units::m << std::endl;
    
    display.draw_cells(toytiling[i]->get_allcell(),"*same");
    display.draw_mergecells(mergetiling[i]->get_allcell(),"*same",0); //0 is normal, 1 is only draw the ones containt the truth cell
    
    display.draw_wires_charge(toytiling[i]->wcmap(),"Fsame",FI);
    display.draw_cells_charge(toytiling[i]->get_allcell(),"Fsame");
  
    
    
    theApp.Run();
  }
  

  
 

  

  return 0;
  
} // main()
