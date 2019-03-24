
/***************************************************************************/
/*             Magnify-to-celltree converter for protoDUNE-SP              */
/***************************************************************************/
/* features:                                                               */
/* - standalone                                                            */

#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH2I.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TDirectory.h"
#include "stdlib.h"

using namespace std;

int main(int argc, char* argv[])
{

  if (argc < 4) {
    cerr << "usage: wire-cell-pdsp-magnify2ctree magnify_file.root celltree_file_name.root frame_tag [run_number (default:0)] [event_number (default:0)]" << endl;
    return 1;
  }
  const char* infile = argv[1];
  const char* outfile = argv[2];
  const char* intag = argv[3];
  std::cerr <<"[wgu] input magnify file: " << infile << std::endl;
  std::cerr <<"[wgu] output celltree file: " << outfile << std::endl;
  std::cerr <<"[wgu] converting frame tag: " << intag << std::endl;
  int fEvent=0, fRun=0, fSubRun=0;
  if(argc>4){
    fRun = atoi(argv[4]);
    fEvent = atoi(argv[5]);
    std::cerr <<"[wgu] run number: " << fRun << " , event number: " << fEvent << std::endl;
  }

  TFile* fInFile = TFile::Open(infile);
  auto hu_orig = (TH2F*)fInFile->Get(TString::Format("hu_%s",intag));
  auto hv_orig = (TH2F*)fInFile->Get(TString::Format("hv_%s",intag));
  auto hw_orig = (TH2F*)fInFile->Get(TString::Format("hw_%s",intag));
  hu_orig->Add(hv_orig);
  hu_orig->Add(hw_orig);

  TFile* fOutFile = new TFile(outfile, "RECREATE");
  TDirectory* subDir = fOutFile->mkdir("Event");
  subDir->cd();

  int fRaw_nChannel;
  std::vector<int> fRaw_channelId;
  TClonesArray *fRaw_wf = new TClonesArray("TH1F");
  auto fEventTree = new TTree("Sim", "Event Tree from Simulation");
  fEventTree->Branch("eventNo", &fEvent);
  fEventTree->Branch("runNo", &fRun);
  fEventTree->Branch("subRunNo", &fSubRun);
  fEventTree->Branch("raw_nChannel", &fRaw_nChannel);  // number of hit channels above threshold
  fEventTree->Branch("raw_channelId" , &fRaw_channelId); // hit channel id; size == raw_nChannel
  fEventTree->Branch("raw_wf", &fRaw_wf, 256000, 0);  // raw waveform adc of each channel

  fRaw_nChannel = hu_orig->GetNbinsX();
  int nticks = hu_orig->GetNbinsY();
  for(int ich=0; ich<fRaw_nChannel; ich++){
    fRaw_channelId.push_back(ich);
    TH1F* h = new((*fRaw_wf)[ich]) TH1F("","", nticks, 0, nticks);
    for(int j=0; j<nticks; j++){
      h->SetBinContent(j+1, hu_orig->GetBinContent(ich+1, j+1));
    }
  }
  fEventTree->Fill();
  
  fOutFile->cd("/Event");
  fEventTree->Write();
  fOutFile->Close();
  fInFile->Close();
  return 0;
  
}
