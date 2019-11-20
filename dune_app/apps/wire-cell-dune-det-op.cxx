#include "WCP2dToy/ToyTiling.h"
#include "WCP2dToy/ClusterDisplay.h"
#include "WCP2dToy/ToyCrawler.h"
#include "WCP2dToy/ToyTracking.h"
#include "WCPData/TPCParams.h"
#include "WCPData/Singleton.h"

#include "WCPData/MergeGeomCell.h"
#include "WCPData/MergeGeomWire.h"
#include "WCPData/SpaceCell.h"
#include "WCPData/MergeSpaceCell.h"

#include "WCPSst/MCTruth.h"

#include "TApplication.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TColor.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TVector3.h"
#include <iostream>

using namespace WCP;
using namespace std;

int main(int argc, char* argv[])
{

  if (argc < 5){
    cerr << "usage: wire-cell-dune-det-op /path/to/celltree.root cluster_num -o[0,1](do_rotation) -p[0,1](is_3mm) " << endl;
    return 1;
  }

  bool rotate_90deg=false;
  bool is_3mm=false;
  bool random_vertices=false;
  int seed=0;
  
  int section_no=0;
  int event_no=0;
  int condition_no=0;

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
     case 'e':
       event_no = atoi(&argv[i][2]);
       break;
     case 's':
       section_no = atoi(&argv[i][2]);
       break;
     case 'c':
       condition_no = atoi(&argv[i][2]);
       break;
       
     }
  }


  if (rotate_90deg) cout<<"Beam is perpendicular to wire planes. ";
  else cout<<"Beam is parallel to wire planes (default). ";
  if (is_3mm) cout<<"Wire pitch is 3 mm."<<endl;
  else cout<<"Wire pitch is 5 mm (default)."<<endl;
  
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

  
  TString filename = argv[1];
  int ncluster = atoi(argv[2]);


  TFile *file = new TFile(filename);
  TTree *T = (TTree*)file->Get("T");
  TTree *TC = (TTree*)file->Get("TC");
  TTree *Trun = (TTree*)file->Get("Trun");
  TTree *TMC = (TTree*)file->Get("TMC");
  
  TGraph2D *shower3D = (TGraph2D*)file->Get("shower3D");

  float unit_dis=0;
  int nrebin=0;
  int total_time_bin=0;
  Trun->SetBranchAddress("nrebin",&nrebin);
  Trun->SetBranchAddress("unit_dis",&unit_dis);
  Trun->SetBranchAddress("total_time_bin",&total_time_bin);
  Trun->GetEntry(0);
  
  TPCParams& mp = Singleton<TPCParams>::Instance();
  
  double pitch_u = pitchU;
  double pitch_v = pitchV;
  double pitch_w = pitchW;
  double time_slice_width = nrebin * unit_dis * 0.5 * units::mm;

  mp.set_pitch_u(pitch_u);
  mp.set_pitch_v(pitch_v);
  mp.set_pitch_w(pitch_w);
  mp.set_ts_width(time_slice_width);
  
  std::cout << "Singleton: " << mp.get_pitch_u() << " " << mp.get_pitch_v() << " " << mp.get_pitch_w() << " " << mp.get_ts_width() << std::endl;

  const int ntime = total_time_bin/nrebin;
  WCP2dToy::ToyTiling **toytiling = new WCP2dToy::ToyTiling*[ntime];
  for (int i=0;i!=ntime;i++){
    toytiling[i] = new WCP2dToy::ToyTiling();
  }

  int time_slice=0;

  const GeomCell *cell = 0;//tt->get_allcell().at(5);
  const GeomWire *wire = 0;//tt->get_allwire().at(3);
  
  double charge=0, x=0,y=0,z=0;
  int cluster_num=0;
  int mcell_id=0;
  TC->SetBranchAddress("time_slice",&time_slice);
  TC->SetBranchAddress("charge",&charge);
  TC->SetBranchAddress("xx",&x);
  TC->SetBranchAddress("yy",&y);
  TC->SetBranchAddress("zz",&z);
  TC->SetBranchAddress("ncluster",&cluster_num);
  TC->SetBranchAddress("mcell_id",&mcell_id);
  TC->SetBranchAddress("cell",&cell);
  
  int u_index=0, v_index=0, w_index=0;
  double u_charge=0, v_charge=0, w_charge=0;
  double u_charge_err=0, v_charge_err=0, w_charge_err=0;

  TC->SetBranchAddress("u_index",&u_index);
  TC->SetBranchAddress("v_index",&v_index);
  TC->SetBranchAddress("w_index",&w_index);
  TC->SetBranchAddress("u_charge",&u_charge);
  TC->SetBranchAddress("v_charge",&v_charge);
  TC->SetBranchAddress("w_charge",&w_charge);
  TC->SetBranchAddress("u_charge_err",&u_charge_err);
  TC->SetBranchAddress("v_charge_err",&v_charge_err);
  TC->SetBranchAddress("w_charge_err",&w_charge_err);

  int apa_no=0, cryostat_no=0;
  int face = 0;
  TC->SetBranchAddress("face",&face);
  TC->SetBranchAddress("apa_no",&apa_no);
  TC->SetBranchAddress("cryostat_no",&cryostat_no);
    
  
  int prev_mcell_id = -1;
  int prev_cluster_num = -1;

  MergeSpaceCellSelection mcells; // save all the cells
  MergeSpaceCellSelection mcells_all; // save all the cells
  int flag = 0;
  MergeSpaceCell *mcell=0;
  SpaceCellSelection cells;
  
  for (int i=0;i!=TC->GetEntries();i++){
    TC->GetEntry(i);
    
    if (cluster_num != prev_cluster_num){
      if (prev_cluster_num!=-1){
	if (mcell->Get_all_spacecell().size()>0){
	  mcells.push_back(mcell);  
	  mcells_all.push_back(mcell);
	}
	// this is a cluster
	mcells.clear();
	cells.clear();
	flag = 0;
      }
    }
    
    if (flag == 0){
      mcell = new MergeSpaceCell();
      mcell->set_id(mcell_id);
      flag = 1;
    }else if (flag==1 && (mcell_id!=prev_mcell_id)){
      mcells.push_back(mcell);
      mcells_all.push_back(mcell);
      mcell = new MergeSpaceCell();
      mcell->set_id(mcell_id);
    }
    GeomCell *cell1 = new GeomCell(cell );
    cell1->set_tpc_no(face + cryostat_no*10000 + apa_no*10);

    toytiling[time_slice]->AddCell(gds,cryostat_no,apa_no,cell1,u_index,v_index,w_index,u_charge,v_charge,w_charge,u_charge_err,v_charge_err,w_charge_err);

    SpaceCell *space_cell = new SpaceCell(cluster_num,*cell1,x*units::cm,charge,0.32*units::cm);
    mcell->AddSpaceCell(space_cell);
    cells.push_back(space_cell);
    
    prev_cluster_num = cluster_num;
    prev_mcell_id = mcell_id;

  }

  // last cluster ... 
  if (mcell!=0){
    if (mcell->Get_all_spacecell().size()>0){
      mcells.push_back(mcell);  
      mcells_all.push_back(mcell);
    }
  }
  
  std::vector<MergeSpaceCellSelection> mcells_time;
  for (int i=0;i!=ntime;i++){
    MergeSpaceCellSelection temp_mcells;
    mcells_time.push_back(temp_mcells);
  }
  
  const WrappedGDS *temp_apa_gds = gds.get_apaGDS(0,0);
  std::pair<double, double> temp_xmm = temp_apa_gds->minmax(0); 
  
  for (int i=0;i!=mcells_all.size();i++){
    MergeSpaceCell *mcell = mcells_all.at(i);
    SpaceCell *cell = mcell->Get_all_spacecell().at(0);
    
    int face = 0;
    float drift_dist=0;
    if (cell->x()>temp_xmm.first && cell->x()<temp_xmm.second) {
    }else{
      if (TMath::Abs(cell->x()-temp_xmm.first) > TMath::Abs(cell->x()-temp_xmm.second)) {
	drift_dist = TMath::Abs(cell->x()-temp_xmm.second);
	face = 1; //"B face +x"
      }else{
	drift_dist = TMath::Abs(cell->x()-temp_xmm.first);
      }
      int time_slice = round(drift_dist / 2.0 / (1.6*units::mm))+800;
      mcells_time.at(time_slice).push_back(mcell);
    }
  }
  
  
   
  
  // for (int i=0;i!=ntime;i++){
  //   GeomCellSelection allcell = toytiling[i]->get_allcell();
  //   GeomWireSelection allwire = toytiling[i]->get_allwire();
  //   cout << "Single Cell: " << i << " "  << allcell.size() << " " << allwire.size() << endl;
  // }
  
  
  TGraph2D *grec = new TGraph2D();
  int num = 0;
  for (int i=0;i!=mcells_all.size();i++){
    // for (int j=0;j!=mcells_all.at(i)->Get_all_spacecell().size();j++){
    //   grec->SetPoint(num,mcells_all.at(i)->Get_all_spacecell().at(j)->x()/units::cm,
    // 		     mcells_all.at(i)->Get_all_spacecell().at(j)->y()/units::cm,
    // 		     mcells_all.at(i)->Get_all_spacecell().at(j)->z()/units::cm
    // 		     );
    //   num ++;
    // }
    grec->SetPoint(i,mcells_all.at(i)->Get_Center().x/units::cm,
    		    mcells_all.at(i)->Get_Center().y/units::cm,
    		    mcells_all.at(i)->Get_Center().z/units::cm);
    
		    
  }

  std::cout << TC->GetEntries() << " " << num << std::endl;

  
  // deal with MC truth ... 
  WCPSst::MCTruth *mctruth = new WCPSst::MCTruth(TMC);
  //mctruth->GetEntry(0);

  // Need the (original) primary vertex to be in the fiducial volume ...
  Point neutrino_vertex = mctruth->find_neutrino_vertex();
  bool contained_flag = false;
  if (fabs(neutrino_vertex[0]) < 360-10 && fabs(neutrino_vertex[1]) < 600-10 && neutrino_vertex[2]>10 && neutrino_vertex[2] < 360-20){
    contained_flag = true;
  }
  std::cout << "neutrino Vertex: " << neutrino_vertex[0] << " " << neutrino_vertex[1] << " " << neutrino_vertex[2] << " " << contained_flag << std::endl;

  // Now need to find the primary electron trajectory ...
  Point primary_vertex = mctruth->find_primary_vertex();
  std::cout << "Primary Vertex: " << primary_vertex[0] << " " << primary_vertex[1] << " " << primary_vertex[2] << std::endl;
  
  // use the vertex and the main point
  // use the direction of the electron
  MCParticle* electron = mctruth->find_primary_electron();
  if (electron != 0){
    std::cout << "Electron's mom: " << electron->startMomentum[0] << " " << electron->startMomentum[1] << " " << electron->startMomentum[2] << std::endl;
  }else{
    std::cout << "No Electrons!" << std::endl;
  }

  // Need to find all the seconary gammas' (energy deposition) starting point and trajectory... 
  // position, the energy of secondary electron etc
  // direction is just gamma's direction?
  MCParticleSelection photons = mctruth->find_primary_photons();
  std::cout << "Photons: " << photons.size() << std::endl;
  
  std::cout << "Neutrino Energies (true, visible): " << mctruth->find_neutrino_true_energy() << " " << mctruth->find_neutrino_visible_energy() << std::endl;
  
  // find gap
  // start from the vertex, and check along the photon direction 
  // the first distance would be the length before a gap
  // the second ditance would be the length after a gap


  // save wire and cell directly? (in addition to the point)
  // gap can be defined as a line between the position and vertex location 

  // Given a position, need to find the corresponding merged cells' location (through merged wire)... 
  // can also cut distance betweem wires ... 
  
  // Given a merged cell, need to find all the wires and then energy ... 
  // dE/dx through projection? 
  

  TString new_file; 
  if (section_no < 10){
    new_file = Form("./out_rootfiles/000%d/%d_%d.root",section_no,event_no,condition_no);
  }else if (section_no < 100){
    new_file = Form("./out_rootfiles/00%d/%d_%d.root",section_no,event_no,condition_no);
  }else if (section_no < 1000){
    new_file = Form("./out_rootfiles/0%d/%d_%d.root",section_no,event_no,condition_no);
  }else{
    new_file = Form("./out_rootfiles/%d/%d_%d.root",section_no,event_no,condition_no);
  }


  if (contained_flag){
    TFile *file1 = new TFile(new_file,"RECREATE");
    TTree *t1 = new TTree("T","T");
    t1->SetDirectory(file1);
    Float_t Enu_true=0, Enu_reco=0;
    Enu_true = mctruth->find_neutrino_true_energy() ;
    Enu_reco = mctruth->find_neutrino_visible_energy();
    
    t1->Branch("Enu_true",&Enu_true,"Enu_true/F");
    t1->Branch("Enu_reco",&Enu_reco,"Enu_reco/F");
    
    t1->Branch("section_no",&section_no,"section_no/I");
    t1->Branch("event_no",&event_no,"event_no/I");
    t1->Branch("condition_no",&condition_no,"condition_no/I");
        
    Int_t nphoton,nele;
    nphoton =  photons.size();
    if (electron == 0){
      nele = 0;
    }else{
      nele = 1;
    }
    t1->Branch("nphoton",&nphoton,"nphoton/I");
    t1->Branch("nele",&nele,"nele/I");
    // save gap related variables for electron (to be added)
    Int_t ele_gap_time_slice[200]={0};
    Float_t ele_gap_dis[200]={0.0}; 
    Float_t ele_gap_charge[200]={0.0};
    Float_t ele_gap_area[200]={0.0};

    Float_t E_ele=0;
    Float_t Theta_ele=0;
    Float_t Phi_ele=0;
    t1->Branch("E_ele",&E_ele,"E_ele/F");
    t1->Branch("Theta_ele",&Theta_ele,"Theta_ele/F");
    t1->Branch("Phi_ele",&Phi_ele,"Phi_ele/F");
   
    t1->Branch("ele_gap_time_slice",ele_gap_time_slice,"ele_gap_time_slice[200]/I");
    t1->Branch("ele_gap_dis",ele_gap_dis,"ele_gap_dis[200]/F");
    // t1->Branch("ele_gap_charge",ele_gap_charge,"ele_gap_charge[200]/F");
    // t1->Branch("ele_gap_area",ele_gap_area,"ele_gap_area[200]/F");

    //std::cout << nele<< std::endl;
    for (int i=0;i!=nele;i++){
      Point p;
      TVector3 dir(electron->startMomentum[0],electron->startMomentum[1],electron->startMomentum[2]);
      
      E_ele = electron->startMomentum[3];
      Theta_ele = dir.Theta();
      Phi_ele = dir.Phi();

      for (int j=0;j!=200;j++){
	p.x = primary_vertex.x*units::cm + j * units::mm * sin(Theta_ele)*cos(Phi_ele);
	p.y = primary_vertex.y*units::cm + j * units::mm * sin(Theta_ele)*sin(Phi_ele);
	p.z = primary_vertex.z*units::cm + j * units::mm * cos(Theta_ele);
	
	//std::cout << p.x/units::cm << " " << p.y/units::cm << " " << p.z/units::cm << std::endl;
	
	int flag = 0;
	short which_cryo = gds.in_which_cryo(p);
	short which_apa = gds.in_which_apa(p);
	const WrappedGDS *apa_gds = gds.get_apaGDS(which_cryo,which_apa);
	if (apa_gds!=NULL) {
	  std::pair<double, double> xmm = apa_gds->minmax(0); 
	  int face = 0;
	  float drift_dist;
	  if (p.x>xmm.first && p.x<xmm.second) {
	  }else{
	    if (TMath::Abs(p.x-xmm.first) > TMath::Abs(p.x-xmm.second)) {
	      drift_dist = TMath::Abs(p.x-xmm.second);
	      face = 1; //"B face +x"
	    }else{
	      drift_dist = TMath::Abs(p.x-xmm.first);
	    }
	    flag = 1;
	    ele_gap_time_slice[j] = round(drift_dist / 2.0 / (1.6*units::mm))+800;
	    
	    const GeomWire* uwire = apa_gds->closest(p,(WirePlaneType_t)0,face);
	    const GeomWire* vwire = apa_gds->closest(p,(WirePlaneType_t)1,face);
	    const GeomWire* wwire = apa_gds->closest(p,(WirePlaneType_t)2,face);
	    
	    int time_slice = ele_gap_time_slice[j];
	    
	    // check with merged cells
	    //std::cout << i << " " << j << " " << mcells_time.at(time_slice).size() << std::endl;
	    
	    ele_gap_dis[j] = 1e9;
	    ele_gap_area[j] = -1;
	    ele_gap_charge[j] = -1;

	    for (int k=0;k!=mcells_time.at(time_slice).size();k++){
	      MergeSpaceCell *mcell = mcells_time.at(time_slice).at(k);

	      // std::cout << mcell->Get_all_spacecell().at(0)->get_cell()->get_face() << " " << face << " " << 
	      // 	mcell->Get_all_spacecell().at(0)->get_cell()->get_cryo() << " " <<  (int)which_cryo << " " <<
	      // 	mcell->Get_all_spacecell().at(0)->get_cell()->get_apa() << " " <<  (int)which_apa << std::endl;
	      
	      if (mcell->Get_all_spacecell().at(0)->get_cell()->get_face()!=face
		  || mcell->Get_all_spacecell().at(0)->get_cell()->get_cryo() != (int)which_cryo
		  || mcell->Get_all_spacecell().at(0)->get_cell()->get_apa() != (int)which_apa
		  ) continue;

	      // check with merged wires
	      GeomWireSelection uwires = mcell->get_uwires();
	      GeomWireSelection vwires = mcell->get_vwires();
	      GeomWireSelection wwires = mcell->get_wwires();
	      
	      double closest_dist_u = 1e9;
	      double closest_dist_v = 1e9;
	      double closest_dist_w = 1e9;

	      //find the closest wire for uplane
	      for (int kk = 0;kk!=uwires.size();kk++){
		if (fabs(apa_gds->wire_dist(*uwire) - apa_gds->wire_dist(*uwires.at(kk))) < closest_dist_u)
		  closest_dist_u = fabs(apa_gds->wire_dist(*uwire) - apa_gds->wire_dist(*uwires.at(kk)));
	      }
	      //find the closest wire for vplane
	      for (int kk = 0;kk!=vwires.size();kk++){
		if (fabs(apa_gds->wire_dist(*vwire) - apa_gds->wire_dist(*vwires.at(kk))) < closest_dist_v)
		  closest_dist_v = fabs(apa_gds->wire_dist(*vwire) - apa_gds->wire_dist(*vwires.at(kk)));
	      }
	      //find the closest wire for wplane
	      for (int kk = 0;kk!=wwires.size();kk++){
		if (fabs(apa_gds->wire_dist(*wwire) - apa_gds->wire_dist(*wwires.at(kk))) < closest_dist_w)
		  closest_dist_w = fabs(apa_gds->wire_dist(*wwire) - apa_gds->wire_dist(*wwires.at(kk)));
	      }
	      if (closest_dist_u+closest_dist_v+closest_dist_w < ele_gap_dis[j])
		ele_gap_dis[j] = (closest_dist_u+closest_dist_v+closest_dist_w)/units::cm;
	    }
	    
	    // std::cout << i << " " << j << " " << time_slice << " " << mcells_time.at(time_slice).size() << " " << ele_gap_dis[j] << std::endl;
	    
	  }
	}
	
	if (flag == 0){
	  ele_gap_time_slice[j] = -1;
	  ele_gap_dis[j] = -1;
	  ele_gap_charge[j] = -1;
	  ele_gap_area[j] = -1;
	}

      }
    }


    //save gap related variables for photons
    // save up to 20 photons
    // 1 mm a point and 200 points = 20 cm
    Int_t gap_time_slice[20][200]={0};
    Float_t gap_dis[20][200]={0.0}; 
    // Float_t gap_charge[200][20];
    // Float_t gap_area[200][20];

    Float_t E_photon[20]={0.0};
    Float_t Theta_photon[20]={0.0};
    Float_t Phi_photon[20]={0.0};

    t1->Branch("E_photon",E_photon,"E_photon[nphoton]/F");
    t1->Branch("Theta_photon",Theta_photon,"Theta_photon[nphoton]/F");
    t1->Branch("Phi_photon",Phi_photon,"Phi_photon[nphoton]/F");
   
    t1->Branch("gap_time_slice",gap_time_slice,"gap_time_slice[nphoton][200]/I");
    t1->Branch("gap_dis",gap_dis,"gap_dis[nphoton][200]/F");
    // t1->Branch("gap_charge",gap_charge,"gap_charge[200][nphoton]/F");
    // t1->Branch("gap_area",gap_area,"gap_area[200][nphoton]/F");
    
    

    for (int i=0;i!=nphoton;i++){
      Point p;
      TVector3 dir(photons.at(i)->startMomentum[0],photons.at(i)->startMomentum[1],photons.at(i)->startMomentum[2]);
      
      E_photon[i] = photons.at(i)->startMomentum[3];
      Theta_photon[i] = dir.Theta();
      Phi_photon[i] = dir.Phi();

      for (int j=0;j!=200;j++){
	p.x = primary_vertex.x*units::cm + j * units::mm * sin(Theta_photon[i])*cos(Phi_photon[i]);
	p.y = primary_vertex.y*units::cm + j * units::mm * sin(Theta_photon[i])*sin(Phi_photon[i]);
	p.z = primary_vertex.z*units::cm + j * units::mm * cos(Theta_photon[i]);
	
	//std::cout << p.x/units::cm << " " << p.y/units::cm << " " << p.z/units::cm << std::endl;
	
	int flag = 0;
	short which_cryo = gds.in_which_cryo(p);
	short which_apa = gds.in_which_apa(p);
	const WrappedGDS *apa_gds = gds.get_apaGDS(which_cryo,which_apa);
	if (apa_gds!=NULL) {
	  std::pair<double, double> xmm = apa_gds->minmax(0); 
	  int face = 0;
	  float drift_dist=0;
	  if (p.x>xmm.first && p.x<xmm.second) {
	  }else{
	    if (TMath::Abs(p.x-xmm.first) > TMath::Abs(p.x-xmm.second)) {
	      drift_dist = TMath::Abs(p.x-xmm.second);
	      face = 1; //"B face +x"
	    }else{
	      drift_dist = TMath::Abs(p.x-xmm.first);
	    }
	    flag = 1;
	    gap_time_slice[i][j] = round(drift_dist / 2.0 / (1.6*units::mm))+800;
	    //std::cout << gap_time_slice[j][i] << std::endl;

	    const GeomWire* uwire = apa_gds->closest(p,(WirePlaneType_t)0,face);
	    const GeomWire* vwire = apa_gds->closest(p,(WirePlaneType_t)1,face);
	    const GeomWire* wwire = apa_gds->closest(p,(WirePlaneType_t)2,face);
	    
	    int time_slice = gap_time_slice[i][j];
	    
	    // check with merged cells
	    //std::cout << i << " " << j << " " << mcells_time.at(time_slice).size() << std::endl;
	    
	    gap_dis[i][j] = 1e9;
	    // gap_area[j][i] = -1;
	    // gap_charge[j][i] = -1;

	    for (int k=0;k!=mcells_time.at(time_slice).size();k++){
	      MergeSpaceCell *mcell = mcells_time.at(time_slice).at(k);

	      // std::cout << mcell->Get_all_spacecell().at(0)->get_cell()->get_face() << " " << face << " " << 
	      // 	mcell->Get_all_spacecell().at(0)->get_cell()->get_cryo() << " " <<  (int)which_cryo << " " <<
	      // 	mcell->Get_all_spacecell().at(0)->get_cell()->get_apa() << " " <<  (int)which_apa << std::endl;
	      
	      if (mcell->Get_all_spacecell().at(0)->get_cell()->get_face()!=face
		  || mcell->Get_all_spacecell().at(0)->get_cell()->get_cryo() != (int)which_cryo
		  || mcell->Get_all_spacecell().at(0)->get_cell()->get_apa() != (int)which_apa
		  ) continue;

	      // check with merged wires
	      GeomWireSelection uwires = mcell->get_uwires();
	      GeomWireSelection vwires = mcell->get_vwires();
	      GeomWireSelection wwires = mcell->get_wwires();
	      
	      double closest_dist_u = 1e9;
	      double closest_dist_v = 1e9;
	      double closest_dist_w = 1e9;

	      //find the closest wire for uplane
	      for (int kk = 0;kk!=uwires.size();kk++){
		if (fabs(apa_gds->wire_dist(*uwire) - apa_gds->wire_dist(*uwires.at(kk))) < closest_dist_u)
		  closest_dist_u = fabs(apa_gds->wire_dist(*uwire) - apa_gds->wire_dist(*uwires.at(kk)));
	      }
	      //find the closest wire for vplane
	      for (int kk = 0;kk!=vwires.size();kk++){
		if (fabs(apa_gds->wire_dist(*vwire) - apa_gds->wire_dist(*vwires.at(kk))) < closest_dist_v)
		  closest_dist_v = fabs(apa_gds->wire_dist(*vwire) - apa_gds->wire_dist(*vwires.at(kk)));
	      }
	      //find the closest wire for wplane
	      for (int kk = 0;kk!=wwires.size();kk++){
		if (fabs(apa_gds->wire_dist(*wwire) - apa_gds->wire_dist(*wwires.at(kk))) < closest_dist_w)
		  closest_dist_w = fabs(apa_gds->wire_dist(*wwire) - apa_gds->wire_dist(*wwires.at(kk)));
	      }
	      if (closest_dist_u+closest_dist_v+closest_dist_w < gap_dis[i][j])
		gap_dis[i][j] = (closest_dist_u+closest_dist_v+closest_dist_w)/units::cm;
	    }
	    
	    std::cout << i << " " << j << " " << gap_time_slice[i][j] << " " << gap_dis[i][j] << std::endl;
	    
	  }
	}
	
	if (flag == 0){
	  gap_time_slice[i][j] = -1;
	  gap_dis[i][j] = -1;
	  // gap_charge[j][i] = -1;
	  // gap_area[j][i] = -1;
	}

      }
    }

    //save dE/dx for photons 
    
    //save dE/dx for electrons

    

    
    t1->Fill();
    file1->Write();
    file1->Close();
  }

  
  // TApplication theApp("theApp",&argc,argv);
  // TCanvas *c1 = new TCanvas("c1","c1",800,600);
  // shower3D->Draw("p");
  // grec->SetMarkerColor(2);
  // grec->Draw("psame");
  // grec->SetMarkerStyle(21);
  // grec->SetMarkerSize(0.5);
  // theApp.Run();


  // //cout << mcells.size() << endl;

  // // do the Toy Crawler
  // std::cout << "Crawling " << std::endl;
  // WCP2dToy::ToyCrawler toycrawler(mcells);
  // //WCP2dToy::ToyCrawler toycrawler(mcells,1,2); //cosmic tune?

  // // test
  // std::cout << "Tracking " << std::endl;
  // WCP2dToy::ToyTracking toytracking(toycrawler);
  // MergeSpaceCellMap& mcells_map = toycrawler.Get_mcells_map();
  // //WCP2dToy::ToyTracking toytracking(toycrawler,1); //cosmic tune?
  // toytracking.IterateMergeTracks(mcells_map);
  
  // //std:cout << mcells.size() << " " << mcells_map.size() << std::endl;
  // // for (int i=0;i!=mcells.size();i++){
  // //   if (mcells_map.find(mcells.at(i)) == mcells_map.end()){
  // //     std::cout << i << std::endl;
  // //   }else{
  // //     std::cout << i << " " << mcells_map[mcells.at(i)].size() << std::endl;
  // //   }
  // // }
  // // for (auto it = mcells_map.begin();it!=mcells_map.end();it++){
  // //   auto it1 = find(mcells.begin(),mcells.end(),it->first);
  // //   if (it1 == mcells.end()){
  // //     std::cout << it->first->Get_all_spacecell().size() << " " 
  // //   		<< it->first->Get_Center().x/units::cm << " " 
  // //   		<< it->first->Get_Center().y/units::cm << " " 
  // //   		<< it->first->Get_Center().z/units::cm << " " 
  // // 		<< it->second.size() << std::endl;
  // //   }
  // // }


  // std::cout << "Good Tracks:     " << toytracking.get_good_tracks().size() <<std::endl;
  // std::cout << "Good Vertices:        " << toytracking.get_good_vertices().size() << std::endl;
  // std::cout << "Bad Tracks:      " << toytracking.get_bad_tracks().size() << std::endl;
  // std::cout << "Parallel Tracks: " << toytracking.get_parallel_tracks().size() << std::endl;
  // std::cout << "Showers:         " << toytracking.get_showers().size() << std::endl;


  // std::cout << "Drawing " << std::endl; 
  // TApplication theApp("theApp",&argc,argv);
  // theApp.SetReturnFromRun(true);
  
  // TCanvas c1("ToyMC","ToyMC",800,600);
  // c1.Draw();
  
  // WCP2dToy::ClusterDisplay display(c1);
  // shower3D_charge->Draw("p");

  // // display.DrawCluster(cells);
  // // display.DrawCluster(mcells);
  // //display.DrawCluster(mcells,toytracking);
  // //display.DrawCrawler(toycrawler,"psame",1);

  // WCVertexSelection& vertices = toytracking.get_good_vertices();
  // //WCVertexSelection& vertices = toytracking.get_vertices();
  // display.DrawVertex(vertices,"psame");
  
  // WCTrackSelection& bad_tracks = toytracking.get_bad_tracks();
  // //display.DrawTracks(bad_tracks,"same",2);

  // WCTrackSelection& short_tracks = toytracking.get_short_tracks();
  // //display.DrawTracks(short_tracks,"psame",4);

  // WCShowerSelection& showers =toytracking.get_showers();
  // if (showers.size() > 0)
  //   display.DrawShower(showers.at(0),"psame",8);
  // // Point p;
  // // p.x = cells.at(0)->x();
  // // p.y = cells.at(0)->y();
  // // p.z = cells.at(0)->z();
  // // display.DrawHough(cells,p,-1,10*units::m);
  

  // theApp.Run();
  // //std::cout << cells.size() << std::endl;
  // //successfully read the TC tree 
  // // TC->GetEntry(0);
  // //cout << x << " " << y << " " << z << " " << charge << " " << time_slice << " " << cluster_num << " " << cell->cross_section() << " " << cell->center().y << endl;
    
    
}
