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


bool flashFilter(const char* file, int eve_num, unsigned int triggerbits)
{
  // light reco and apply a [3, 5] ([3.45, 5.45]) us cut on BNB (extBNB) trigger
  WCP2dToy::ToyLightReco flash(file, 1);
  flash.load_event_raw(eve_num);
  
  // flash.clear_flashes();
  // flash.load_event_raw(eve_num);

  // flash.clear_flashes();
  // flash.clear_flashes();
  // flash.load_event_raw(eve_num);
  
  WCP::OpflashSelection& flashes = flash.get_flashes();
  bool beamspill = false;
  for(auto it = flashes.begin(); it!=flashes.end(); it++){
      Opflash *flash = (*it);
      int type = flash->get_type();
      double time = flash->get_time();
      //cout<<"Flash time: "<<time<<" Type: "<<type<<endl;
      double lowerwindow = 0;
      double upperwindow = 0;
      if((triggerbits>>11) & 1U) { lowerwindow = 3.0; upperwindow = 5.0; }// bnb  
      if((triggerbits>>9) & 1U) { lowerwindow = 3.45; upperwindow = 5.45; } // extbnb
     
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
    if((triggerbits>>11) & 1U){  
	  cout << "BNB Flash time filter: PASS " << argv[1] << endl;
	}
	if((triggerbits>>9) & 1U){  
      	cout << "extBNB Flash time filter: PASS "<< argv[1] << endl;
	}
  }
  else{
	if((triggerbits>>11) & 1U){  
    	cout << "BNB Flash time filter: FAIL " << argv[1]<< endl;
	}  
	if((triggerbits>>9) & 1U){  
    	cout << "extBNB Flash time filter: FAIL "<< argv[1] << endl;
	}  
  }
  return 0;
  
} // main()
