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


std::pair<double,double> cal_mean_rms(TH1F *hist){
  double mean=0, rms=0;
  if(hist->GetEntries()==0){
    return std::make_pair(0,0);
  }

  int nbin = hist->GetNbinsX();
  mean = hist->GetSum()/(float)nbin;
  for(int i=0; i!=nbin; i++)
  {
    rms += pow(hist->GetBinContent(i+1)-mean,2);
  }
  rms = sqrt(rms/(float)nbin);

  TH1F* htemp = new TH1F("pc3erslj","",100,mean-10*rms,mean+10*rms);
  for (int j=0;j!=nbin;j++){
   if (fabs(hist->GetBinContent(j+1)-mean)<6*rms)
     htemp->Fill(hist->GetBinContent(j+1));
  }
  if(htemp->GetEntries()==0){
    delete htemp;
    return std::make_pair(0,0);
  }
  
  double xq = 0.5;
  double par[3];
  htemp->GetQuantiles(1,&par[1],&xq);
  xq = 0.5 - 0.34;
  htemp->GetQuantiles(1,&par[0],&xq);
  xq = 0.5 + 0.34;
  htemp->GetQuantiles(1,&par[2],&xq);
  rms = sqrt((pow(par[0]-par[1],2)+pow(par[2]-par[1],2))/2.);

  delete htemp;
  return std::make_pair(par[1],rms);
}

int planeid(int chid){
  int chid1 = chid % 2560;
  if(chid1<800) return 0;
  else if(chid1<1600) return 1;
  else return 2;
}

void fill_2dw1d(TH2I* h_orig, TH1F* hwf, int chn){
  // cout << "fill chn= " << chn << endl; 
  int nticks = hwf->GetNbinsX();
  for(int ind=1; ind<=nticks; ind++){
    int rawAdc = hwf->GetBinContent(ind);
    h_orig->SetBinContent(chn+1, ind, rawAdc);
  }
}
void fill_2dw1d(TH2F* h_orig, TH1F* hwf, int chn, double baseline, double scale=1){
  // cout << "fill chn= " << chn << endl; 
  int nticks = hwf->GetNbinsX();
  for(int ind=1; ind<=nticks; ind++){
    int rawAdc = hwf->GetBinContent(ind) - baseline;
    h_orig->SetBinContent(chn+1, ind, rawAdc * scale);
  }
}



int main(int argc, char* argv[])
{

  if (argc !=3) {
    cerr << "usage: wire-cell-pdsp-ctree2magnify celltree.root fEntry" << endl;
    return 1;
  }

  const char* root_file = argv[1];
  int fEntry = atoi(argv[2]); // local entry in celltree root file
  //const char* dltag = argv[4]; // additional tag

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


  tt->GetEntry(fEntry);
  cout << "Run: " << runNo << endl;
  cout << "Event: " << eventNo << endl;
  cout << "subRun: " << subRunNo << endl;
  cout << "raw_nChannel: " << raw_nChannel << endl;

  TString title;
  //title.Form("magnify_%04d_%d_%s.root", runNo, eventNo, dltag);
  title.Form("cmagnify_%04d_%d.root", runNo, eventNo);
  cout << "Output file: " << title << endl;
  TFile* ofile = new TFile(title.Data(),"recreate");

  for(int ich=0; ich<raw_nChannel; ich++){
    int chId = raw_channelId->at(ich);
    int planeId = planeid(chId);
    TH1F* hwf = (TH1F*)raw_wf->At(ich);

    auto temp = cal_mean_rms(hwf);
    double baseline = temp.first;

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

  // hu_decon->Rebin2D(1,4);
  // hv_decon->Rebin2D(1,4);
  // hw_decon->Rebin2D(1,4);
  hu_decon->Write();
  hv_decon->Write();
  hw_decon->Write();

  ofile->Close();

  file->Close();
  return 0;
  
}
