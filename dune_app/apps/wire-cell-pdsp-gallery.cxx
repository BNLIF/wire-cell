/**
 * Example for looping over protoDUNE channels and quality check
 */
#include "WireCellNav/SliceDataSource.h"
#include "WireCellSst/FrameDataSource.h"
#include "WireCellTiling/TileMaker.h"

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"

#include <ctime>
#include <fstream>
#include <iostream>

#include "LedgeIdentify.h"
#include "StickyCodeIdent.h"

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

//
int main(int argc, char* argv[])
{

  if (argc < 3) {
  	cerr << "usage: wire-cell-pdune-chan-qual celltree_file_name.root output.root" << endl;
  	exit (1);
  }

  TH1F* hfft[3];
  hfft[0]= new TH1F("hfft_u","", 3000,0,1000);
  hfft[1]= new TH1F("hfft_v","", 3000,0,1000);
  hfft[2]= new TH1F("hfft_w","", 3000,0,1000);
  hfft[0]->GetXaxis()->SetTitle("Frequency [kHz]");
  hfft[0]->SetLineColor(1);
  hfft[1]->SetLineColor(2);
  hfft[2]->SetLineColor(4);  
  int nfft[3] = {0,0,0};
  TH2F* hfft2= new TH2F("hfft2d","",15360,0,15360, 3000,0,1.0);
  hfft2->GetXaxis()->SetTitle("Channel");
  hfft2->GetYaxis()->SetTitle("Frequency [MHz]");
  TH1F* hEnc[3];
  hEnc[0]= new TH1F("hEnc_u","", 150,0,1500);
  hEnc[1]= new TH1F("hEnc_v","", 150,0,1500);
  hEnc[2]= new TH1F("hEnc_w","", 150,0,1500);
  hEnc[0]->GetXaxis()->SetTitle("ENC [e^{-}]");
  hEnc[0]->SetLineColor(1);
  hEnc[1]->SetLineColor(2);
  hEnc[2]->SetLineColor(4);  

  const char* infile = argv[1];
  const char* outfile = argv[2];
  TFile* tfile = TFile::Open(infile);
  TTree* tree = dynamic_cast<TTree*>(tfile->Get("/Event/Sim"));
  WireCellSst::FrameDataSource fds(*tree, "raw");
  int eventNo, runNo;
  tree->SetBranchAddress("eventNo",&eventNo);
  tree->SetBranchAddress("runNo",&runNo);
    
  // Loop over frames (aka "events")
  size_t nframes = fds.size();
  cout << "FDS: " << nframes << " frames" << endl;
  for (size_t iframe = 0; iframe < nframes; ++iframe) {
  	int iframe_got = fds.jump(iframe);
  	if (iframe_got < 0) {
  	    cerr << "Failed to get frame " << iframe << endl;
  	    exit(1); // real code may want to do something less drastic
  	}
  	cout << "FDS: jumping to frame #" << iframe << " eventNo " << eventNo << endl;

    const WireCell::Frame& frame = fds.get();
    size_t ntraces = frame.traces.size();
    cout << "protoDUNE::frame.traces.size(): " << ntraces << endl;
    for (size_t ind=0; ind<ntraces; ++ind) {
      const WireCell::Trace& trace = frame.traces[ind];
      int tbin = trace.tbin;
      int chid = trace.chid;
      int pid = planeid(chid);
      // cout << "pdsp ch: " << chid << endl;
      int nbins = trace.charge.size();
      TH1F* hCharge = dynamic_cast<TH1F*>(trace.hCharge);
      std::pair<double, double> temp = cal_mean_rms(hCharge);
      // exclude signals and fill with mean
      for (int j=0;j!=nbins;j++){
        if(fabs(hCharge->GetBinContent(j+1) - temp.first) > 3.5*temp.second){
          hCharge->SetBinContent(j+1, temp.first);
        }
      }

      hEnc[pid]->Fill(temp.second * 130.); // FIXME: what is the ENC/ADC factor?

      TH1F* hmag = (TH1F*)hCharge->FFT(0, "MAG");
      if(hmag->GetBinContent(2)<5E3){ // 5E3 avoids the extreme case
        nfft[pid]++;
        for(int i=1; i<3000; i++){
        hfft[pid]->AddBinContent(i+1, hmag->GetBinContent(i+1));
        // double prevcont = hfft2->GetBinContent(chid%2560+1,i+1);
        hfft2->SetBinContent(chid+1, i+1, hmag->GetBinContent(i+1));
        }
      }
      hmag->Delete();

    }// traces

  }// frames

  TFile* ofile = new TFile(outfile,"recreate");
  for(int i=0; i<3; i++){
    hfft[i]->Scale(1./nfft[i]);
    hfft[i]->Write();
    hEnc[i]->Write();
  }
  hfft2->Rebin2D(4,4);
  hfft2->SetMaximum(15000);
  hfft2->Write();
  ofile->Close();
  return 0;

} // main()
