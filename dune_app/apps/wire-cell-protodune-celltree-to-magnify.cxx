/*
* 09/06/2018, Wenqiang Gu (wgu@bnl.gov)
* A standalone converter for celltree-to-magnify
* currently using standalong classes while not classes from wirecell framework
*/

#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH2I.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TF1.h"
#include "TDirectory.h"

using namespace std;

// 
class ChannelGeo{
  public:
    int chn;
    int tpc;
    int plane;
    int wire;
    double sx; // start of x
    double sy;
    double sz;
    double ex; // end of x
    double ey;
    double ez;
};

class ChannelMap{
  public:
    std::map<int, std::vector<ChannelGeo> > m_map;
    void add(ChannelGeo cg);
    std::vector<ChannelGeo> find(int chn);
    int channel_to_plane(int);
};

void ChannelMap::add(ChannelGeo cg){
    if(0==m_map.count(cg.chn)){
      std::vector<ChannelGeo> vv;
      vv.push_back(cg);
      m_map[cg.chn] =  vv;
    }
    else{
      // some-to-one mapping for wire-channel 
      m_map[cg.chn].push_back(cg);
    }
}

std::vector<ChannelGeo> ChannelMap::find(int chn){
  if(0==m_map.count(chn)) {
    std::cout << "Error. No channel found for " << chn << std::endl;
  }
  return m_map[chn];
}

int ChannelMap::channel_to_plane(int chn){
  std::vector<ChannelGeo> v_cg = find(chn);
  ChannelGeo cg = v_cg[0]; // only one channel-to-tpc mapping
  return cg.plane;
}

void fill_2dw1d(TH2I* h_orig, TH1F* hwf, int chn){
  // cout << "fill chn= " << chn << endl; 
  for(int ind=1; ind<=6000; ind++){
    int rawAdc = hwf->GetBinContent(ind);
    h_orig->SetBinContent(chn+1, ind, rawAdc);
  }
}
void fill_2dw1d(TH2F* h_orig, TH1F* hwf, int chn, double baseline, double scale=1){
  // cout << "fill chn= " << chn << endl; 
  for(int ind=1; ind<=6000; ind++){
    int rawAdc = hwf->GetBinContent(ind) - baseline;
    h_orig->SetBinContent(chn+1, ind, rawAdc * scale);
  }
}



int main(int argc, char* argv[])
{

  if (argc < 4) {
    cerr << "usage: wire-cell-protodune-celltree-to-magnify /path/to/wire-cell/input_data_files/protodune-wires-larsoft-v1.txt /path/to/celltree.root local_entry" << endl;
    return 1;
  }

  const char* map_file = argv[1];
  const char* root_file = argv[2];
  int fEntry = atoi(argv[3]); // local entry in celltree root file
  //const char* dltag = argv[4]; // additional tag

  //ifstream in("protodune-wires-larsoft-v1.txt");
  ifstream in(map_file);
  string dummyLine;
  getline(in, dummyLine);
  ChannelMap cmap;
  ChannelGeo cg;
  while(in >> cg.chn){
    in >> cg.tpc; in >> cg.plane; in >> cg.wire;
    in >> cg.sx; in >> cg.sy; in >> cg.sz;
    in >> cg.ex; in >> cg.ey; in >> cg.ez;
    cmap.add(cg); 
  }
  in.close();


  //TFile* file = TFile::Open("celltree_run003983_0001_dl2.root");
  TFile* file = TFile::Open(root_file);
  TTree* tt = dynamic_cast<TTree*>(file->Get("Event/Sim"));
  tt->SetBranchStatus("*",0);
  tt->SetBranchStatus("eventNo",1);
  tt->SetBranchStatus("runNo",1);
  tt->SetBranchStatus("subRunNo",1);
  tt->SetBranchStatus("raw_nChannel",1);
  tt->SetBranchStatus("raw_channelId",1);
  tt->SetBranchStatus("raw_wf",1);

  int eventNo, runNo, subRunNo;
  tt->SetBranchAddress("eventNo", &eventNo);
  tt->SetBranchAddress("runNo", &runNo);
  tt->SetBranchAddress("subRunNo", &subRunNo);
  int raw_nChannel;
  //std::vector<int>* raw_channelId = new std::vector<int>;
  //TClonesArray* raw_wf = new TClonesArray;
  std::vector<int>* raw_channelId =0;
  TClonesArray* raw_wf =0;
  tt->SetBranchAddress("raw_nChannel", &raw_nChannel);
  tt->SetBranchAddress("raw_channelId", &raw_channelId);
  tt->SetBranchAddress("raw_wf", &raw_wf);

  //TH1::AddDirectory(kFALSE);
  int nchn_u=15360;
  int nchn_v=15360;
  int nchn_w=15360;
  int total_time_bin=6000;
  TH2I* hu_orig = new TH2I("hu_orig","hu_orig",nchn_u,-0.5,nchn_u-0.5,total_time_bin,0,total_time_bin);
  TH2I* hv_orig = new TH2I("hv_orig","hv_orig",nchn_v,-0.5,nchn_v-0.5,total_time_bin,0,total_time_bin);
  TH2I* hw_orig = new TH2I("hw_orig","hw_orig",nchn_w,-0.5,nchn_w-0.5,total_time_bin,0,total_time_bin);
  //
  TH1I* hu_baseline = new TH1I("hu_baseline","hu_baseline",nchn_u,-0.5,nchn_u-0.5);
  TH1I* hv_baseline = new TH1I("hv_baseline","hv_baseline",nchn_v,-0.5,nchn_v-0.5);
  TH1I* hw_baseline = new TH1I("hw_baseline","hw_baseline",nchn_w,-0.5,nchn_w-0.5);
  TH1I* hu_threshold = new TH1I("hu_threshold","hu_threshold",nchn_u,-0.5,nchn_u-0.5);
  TH1I* hv_threshold = new TH1I("hv_threshold","hv_threshold",nchn_v,-0.5,nchn_v-0.5);
  TH1I* hw_threshold = new TH1I("hw_threshold","hw_threshold",nchn_w,-0.5,nchn_w-0.5);
  //
  TH2F* hu_raw = new TH2F("hu_raw","hu_raw",nchn_u,-0.5,nchn_u-0.5,total_time_bin,0,total_time_bin);
  TH2F* hv_raw = new TH2F("hv_raw","hv_raw",nchn_v,-0.5,nchn_v-0.5,total_time_bin,0,total_time_bin);
  TH2F* hw_raw = new TH2F("hw_raw","hw_raw",nchn_w,-0.5,nchn_w-0.5,total_time_bin,0,total_time_bin);
  TH2F* hu_decon = new TH2F("hu_decon","hu_decon",nchn_u,-0.5,nchn_u-0.5,total_time_bin,0,total_time_bin);
  TH2F* hv_decon = new TH2F("hv_decon","hv_decon",nchn_v,-0.5,nchn_v-0.5,total_time_bin,0,total_time_bin);
  TH2F* hw_decon = new TH2F("hw_decon","hw_decon",nchn_w,-0.5,nchn_w-0.5,total_time_bin,0,total_time_bin);


  bool first_u=true;
  bool first_v=true;
  bool first_w=true;
  TF1* fp0 = new TF1("fp0", "pol0", 0, 15360);
  
  tt->GetEntry(fEntry);
  cout << "Run: " << runNo << endl;
  cout << "Event: " << eventNo << endl;
  cout << "subRun: " << subRunNo << endl;
  cout << "raw_nChannel: " << raw_nChannel << endl;

  TString title;
  //title.Form("magnify_%04d_%d_%s.root", runNo, eventNo, dltag);
  title.Form("magnify_%04d_%d.root", runNo, eventNo);
  cout << "Output file: " << title << endl;
  TFile* ofile = new TFile(title.Data(),"recreate");

  for(int ich=0; ich<raw_nChannel; ich++){
    int chId = raw_channelId->at(ich);
    int planeId = cmap.channel_to_plane(chId);
    TH1F* hwf = (TH1F*)raw_wf->At(ich);
    hwf->Fit(fp0,"RQ");
    double baseline = fp0->GetParameter(0);
    if(planeId==0){
      fill_2dw1d(hu_orig, hwf, chId);
      fill_2dw1d(hu_raw, hwf, chId, baseline);
      fill_2dw1d(hu_decon, hwf, chId, baseline, 125);
      hu_baseline->SetBinContent(chId+1, baseline);
      hu_threshold->SetBinContent(chId+1, baseline);
    }
    else if(planeId==1){
      fill_2dw1d(hv_orig, hwf, chId);
      fill_2dw1d(hv_raw, hwf, chId, baseline);
      fill_2dw1d(hv_decon, hwf, chId, baseline, 125);
      hv_baseline->SetBinContent(chId+1, baseline);
      hv_threshold->SetBinContent(chId+1, baseline);
    }
    else if(planeId==2){
      fill_2dw1d(hw_orig, hwf, chId);
      fill_2dw1d(hw_raw, hwf, chId, baseline);
      fill_2dw1d(hw_decon, hwf, chId, baseline, 125);
      hw_baseline->SetBinContent(chId+1, baseline);
      hw_threshold->SetBinContent(chId+1, baseline);
    }
  }

  hu_orig->Write();
  hv_orig->Write();
  hw_orig->Write();


  hu_baseline->Write();
  hv_baseline->Write();
  hw_baseline->Write();
  hu_threshold->Write();
  hv_threshold->Write();
  hw_threshold->Write();

  hu_raw->Write();
  hv_raw->Write();
  hw_raw->Write();


  hu_decon->Rebin2D(1,4);
  hv_decon->Rebin2D(1,4);
  hw_decon->Rebin2D(1,4);
  hu_decon->Write();
  hv_decon->Write();
  hw_decon->Write();

  ofile->Close();

  file->Close();
  return 0;
  
}
