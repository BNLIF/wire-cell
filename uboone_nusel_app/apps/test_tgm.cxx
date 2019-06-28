#include "WireCellSst/GeomDataSource.h"
#include "WireCellSst/DatauBooNEFrameDataSource.h"
#include "WireCellSst/ToyuBooNESliceDataSource.h"
#include "WireCellSst/uBooNESliceDataSource.h"

#include "WireCell2dToy/ToyEventDisplay.h"
#include "WireCell2dToy/ToyTiling.h"
#include "WireCell2dToy/BadTiling.h"
#include "WireCell2dToy/LowmemTiling.h"
#include "WireCell2dToy/uBooNE_L1SP.h"
#include "WireCell2dToy/WireCellHolder.h"
#include "WireCell2dToy/ToyLightReco.h"

#include "WireCell2dToy/MergeToyTiling.h"
#include "WireCell2dToy/TruthToyTiling.h"
#include "WireCell2dToy/SimpleBlobToyTiling.h"

#include "WireCell2dToy/ChargeSolving.h"
#include "WireCell2dToy/ToyMatrix.h"
#include "WireCell2dToy/ToyMatrixExclusive.h"
#include "WireCell2dToy/ToyMatrixKalman.h"
#include "WireCell2dToy/ToyMatrixIterate.h"
#include "WireCell2dToy/ToyMatrixIterate_SingleWire.h"
#include "WireCell2dToy/ToyMatrixIterate_Only.h"

#include "WireCell2dToy/ToyMatrixMarkov.h"
#include "WireCell2dToy/ToyMetric.h"
#include "WireCell2dToy/BlobMetric.h"
#include "WireCellData/TPCParams.h"
#include "WireCellData/Singleton.h"

#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/MergeGeomWire.h"

#include "WireCellData/Slim3DCluster.h"
#include "WireCellData/Slim3DDeadCluster.h"
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
#include "WireCell2dToy/DataSignalGaus_ROI.h"
#include "WireCell2dToy/DataSignalWien_ROI.h"

#include "WireCell2dToy/uBooNE_Data_2D_Deconvolution.h"
#include "WireCell2dToy/uBooNE_Data_ROI.h"
#include "WireCell2dToy/uBooNE_Data_After_ROI.h"
#include "WireCell2dToy/uBooNE_Data_After_ROI_gaus.h"
#include "WireCell2dToy/pd_Data_FDS.h"
#include "WireCell2dToy/uBooNE_Data_Error.h"
#include "WireCell2dToy/ExecMon.h"
#include "WireCell2dToy/ToyDataQuality.h"

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
  if (argc < 2){
    cerr << "usage: wire-cell-uboone /path/to/match.root" << endl;
    return 1;
  }

  const char* root_file = argv[1];
  TFile *file1 = new TFile(root_file);
  TTree *T = (TTree*)file1->Get("Trun");
  Int_t runNo, subRunNo, eventNo;
  T->SetBranchAddress("runNo",&runNo);
  T->SetBranchAddress("subRunNo",&subRunNo);
  T->SetBranchAddress("eventNo",&eventNo);
  unsigned int triggerbits;
  T->SetBranchAddress("triggerBits",&triggerbits);
  T->GetEntry(0);
  double lowerwindow = 0;
  double upperwindow = 0;

  //  if(triggerbits==2048) { lowerwindow = 3.1875; upperwindow = 4.96875; }// bnb  
  //if(triggerbits==512) { lowerwindow = 3.5625; upperwindow = 5.34375; } // extbnb
  // better timing ... 
  // if(triggerbits==2048) { lowerwindow = 3.1625; upperwindow = 4.96875; }// bnb  
  //if(triggerbits==512) { lowerwindow = 3.5375; upperwindow = 5.34375; } // extbnb
  
  //enlarge window ... 
  if(triggerbits==2048) { lowerwindow = 3.0; upperwindow = 5.0; }// bnb  
  if(triggerbits==512) { lowerwindow = 3.45; upperwindow = 5.45; } // extbnb

  TTree *T_flash = (TTree*)file1->Get("T_flash");
  Double_t time;
  Int_t type;
  Int_t flash_id;
  T_flash->SetBranchAddress("time",&time);
  T_flash->SetBranchAddress("type",&type);
  T_flash->SetBranchAddress("flash_id",&flash_id);
  
  TTree *T_match = (TTree*)file1->Get("T_match");
  Int_t tpc_cluster_id;
  Int_t event_type;
  T_match->SetBranchAddress("tpc_cluster_id",&tpc_cluster_id);
  T_match->SetBranchAddress("flash_id",&flash_id);
  T_match->SetBranchAddress("event_type",&event_type);

  std::set<int> saved_flash_ids;
  std::map<int,double> saved_id_time_map;
  //Int_t saved_flash_id = -1;
  //Double_t saved_time = 1e9;
  
  for (int i=0;i!=T_flash->GetEntries();i++){
    T_flash->GetEntry(i);
    if (type==2){
      if ( time >= lowerwindow && time <=upperwindow){
	saved_flash_ids.insert(flash_id);
	saved_id_time_map[flash_id]=time;
	//saved_times.push_back(time);
	//saved_flash_id = flash_id;
	//saved_time = time;
      }
    }
  }

  if (saved_flash_ids.size()>0){
    for (int i=0;i!=T_match->GetEntries();i++){
      T_match->GetEntry(i);
      
      if (saved_flash_ids.find(flash_id)!=saved_flash_ids.end()){
	int flag_tgm = (event_type >> 3) & 1U;
	int flag_low_energy = (event_type >> 4) & 1U;
	int flag_lm = (event_type >> 1) & 1U;
	int flag_fully_contained = (event_type >> 2) & 1U;

	std::cout << runNo << "_" << subRunNo << "_" << eventNo << " " << flash_id << " " << tpc_cluster_id << " " << saved_id_time_map[flash_id] << " " << event_type << " " << flag_low_energy << " " << flag_lm << " " << flag_tgm << " " << flag_fully_contained << std::endl;

      }
      
    }
  }
  
  
  //T_flash->GetEntries(0);
  //T_match->GetEntries(0);
  
  // std::cout << lowerwindow << " " << upperwindow << std::endl;
}
