#include "WCPSst/GeomDataSource.h"
#include "WCPSst/DatauBooNEFrameDataSource.h"
#include "WCPSst/ToyuBooNESliceDataSource.h"
#include "WCPSst/uBooNESliceDataSource.h"

#include "WCP2dToy/ToyEventDisplay.h"
#include "WCP2dToy/ToyTiling.h"
#include "WCP2dToy/BadTiling.h"
#include "WCP2dToy/LowmemTiling.h"
#include "WCP2dToy/uBooNE_L1SP.h"
#include "WCP2dToy/WCPHolder.h"
#include "WCP2dToy/ToyLightReco.h"

#include "WCP2dToy/MergeToyTiling.h"
#include "WCP2dToy/TruthToyTiling.h"
#include "WCP2dToy/SimpleBlobToyTiling.h"

#include "WCP2dToy/ChargeSolving.h"
#include "WCP2dToy/ToyMatrix.h"
#include "WCP2dToy/ToyMatrixExclusive.h"
#include "WCP2dToy/ToyMatrixKalman.h"
#include "WCP2dToy/ToyMatrixIterate.h"
#include "WCP2dToy/ToyMatrixIterate_SingleWire.h"
#include "WCP2dToy/ToyMatrixIterate_Only.h"

#include "WCP2dToy/ToyMatrixMarkov.h"
#include "WCP2dToy/ToyMetric.h"
#include "WCP2dToy/BlobMetric.h"
#include "WCPData/TPCParams.h"
#include "WCPData/Singleton.h"

#include "WCPData/MergeGeomCell.h"
#include "WCPData/MergeGeomWire.h"

#include "WCPData/Slim3DCluster.h"
#include "WCPData/Slim3DDeadCluster.h"
//#include "WCPNav/SliceDataSource.h"


#include "WCPNav/FrameDataSource.h"
#include "WCPNav/SimDataSource.h"
#include "WCPNav/SliceDataSource.h"
#include "WCPSst/Util.h"
#include "WCPData/SimTruth.h"
#include "WCP2dToy/ToyDepositor.h"
#include "WCPNav/GenerativeFDS.h"
#include "WCP2dToy/ToySignalSimu.h"
#include "WCP2dToy/ToySignalSimuTrue.h"
#include "WCP2dToy/DataSignalGaus_ROI.h"
#include "WCP2dToy/DataSignalWien_ROI.h"

#include "WCP2dToy/uBooNE_Data_2D_Deconvolution.h"
#include "WCP2dToy/uBooNE_Data_ROI.h"
#include "WCP2dToy/uBooNE_Data_After_ROI.h"
#include "WCP2dToy/uBooNE_Data_After_ROI_gaus.h"
#include "WCP2dToy/pd_Data_FDS.h"
#include "WCP2dToy/uBooNE_Data_Error.h"
#include "WCP2dToy/ExecMon.h"
#include "WCP2dToy/ToyDataQuality.h"

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

using namespace WCP;
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
  int time_offset;
  T->SetBranchAddress("time_offset",&time_offset);
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
  if((triggerbits>>11) & 1U) { lowerwindow = 3.0; upperwindow = 5.0; }// bnb
  if ((triggerbits>>12) & 1U) { lowerwindow=4.9295; upperwindow=16.6483;} //NUMI
  if( ((triggerbits>>9) & 1U) && time_offset != 3) { lowerwindow = 3.45; upperwindow = 5.45; } // extbnb
  if (((triggerbits>>9) & 1U) && time_offset == 3) {lowerwindow=5.3045; upperwindow=17.0233;} // EXTNUMI

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
      if ( time > lowerwindow && time <upperwindow){
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
