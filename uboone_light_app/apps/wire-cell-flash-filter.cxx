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


bool flashFilter(const char* file, int eve_num, unsigned int triggerbits)
{
  // light reco and apply a [3, 5] ([3.45, 5.45]) us cut on BNB (extBNB) trigger
  WireCell2dToy::ToyLightReco flash(file, 1);
  flash.load_event_raw(eve_num);
  
  // flash.clear_flashes();
  // flash.load_event_raw(eve_num);

  // flash.clear_flashes();
  // flash.clear_flashes();
  // flash.load_event_raw(eve_num);
  
  WireCell::OpflashSelection& flashes = flash.get_flashes();
  bool beamspill = false;
  for(auto it = flashes.begin(); it!=flashes.end(); it++){
      Opflash *flash = (*it);
      int type = flash->get_type();
      double time = flash->get_time();
      //cout<<"Flash time: "<<time<<" Type: "<<type<<endl;
      double lowerwindow = 0;
      double upperwindow = 0;
      if(triggerbits==2048) { lowerwindow = 3; upperwindow = 5; }// bnb
      if(triggerbits==512) { lowerwindow = 3.45; upperwindow = 5.45; } // extbnb
      if(type == 2 && time > lowerwindow && time < upperwindow)
      {
          beamspill = true;
      }
      if(beamspill){ break; }
  }
  return beamspill;
}


int main(int argc, char* argv[])
{
  if (argc < 3) {
    cerr << "usage: wire-cell-uboone /path/to/celltree.root eve_num" << endl;
    return 1;
  }
   TH1::AddDirectory(kFALSE);
  
  const char* root_file = argv[1];  
  int eve_num = atoi(argv[2]);

  //ExecMon em("Starting");
  // load Trun
  TFile *file1 = new TFile(root_file);
  TTree *T = (TTree*)file1->Get("Trun");
  //int event_no=0, run_no=0, subrun_no=0;
  //T->SetBranchAddress("eventNo",&event_no);
  //T->SetBranchAddress("runNo",&run_no);
  //T->SetBranchAddress("subRunNo",&subrun_no);
  unsigned int triggerbits;
  T->SetBranchAddress("triggerBits",&triggerbits);

  T->GetEntry(eve_num);
  //cout << em("load data") << endl;

  //  std::cout << triggerbits << " " << argv[1] << std::endl;
  
  // flash time filter
  bool beamspill = false;
  beamspill = flashFilter(root_file, eve_num, triggerbits);
  //cout << em("Flash filter") <<endl;
  //cout << "Run No: " << run_no << " " << subrun_no << " " << event_no << endl;
  if(beamspill){
	if(triggerbits==2048){  
	  cout << "BNB Flash time filter: PASS " << argv[1] << endl;
	}
	if(triggerbits==512){  
      	cout << "extBNB Flash time filter: PASS "<< argv[1] << endl;
	}
  }
  else{
	if(triggerbits==2048){  
    	cout << "BNB Flash time filter: FAIL " << argv[1]<< endl;
	}  
	if(triggerbits==512){  
    	cout << "extBNB Flash time filter: FAIL "<< argv[1] << endl;
	}  
  }
  return 0;
  
} // main()
