/**
 * Example for looping over protoDUNE channels and quality check
 */
#include "WireCellNav/SliceDataSource.h"
#include "WireCellSst/FrameDataSource.h"
#include "WireCellTiling/TileMaker.h"

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TSpectrum.h"

#include <ctime>
#include <fstream>
#include <iostream>


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

int minBin(TH1F* hsig, int bin1, int bin2){
  int minbin = bin1;
  double min = hsig->GetBinContent(bin1);
  for(int i=bin1+1; i<=bin2; i++){
    if( hsig->GetBinContent(i) < min ){
      min = hsig->GetBinContent(i);
      minbin = i;
    }
  }
  return minbin;
}

int crossingBin(TH1F* hsig, int bin1, int bin2, double baseline){
  int mid = (bin1+bin2)/2;
  int crx = 0;
  for(int i=bin1+1; i<bin2; i++){
    if(hsig->GetBinContent(i-1)>baseline
       && hsig->GetBinContent(i) < baseline){
      // if(fabs(crx-mid)<15) 
      crx = i;
    }
  }

  return crx;
}

double area(TH1F* hsig, int bin1, int bin2, double baseline){
  double ret=0;
  for(int i=bin1; i<=bin2; i++){
    ret += hsig->GetBinContent(i) - baseline;
  }
  return fabs(ret);
}


int planeid(int chid){
  int chid1 = chid % 2560;
  if(chid1<800) return 0;
  else if(chid1<1600) return 1;
  else return 2;
}

typedef std::pair<int, int> BinRange;
//
int main(int argc, char* argv[])
{

  if (argc < 3) {
  	cerr << "usage: wire-cell-pdsp-imbalance celltree_file_name.root output.root" << endl;
  	exit (1);
  }

  const char* infile = argv[1];
  const char* outfile = argv[2];
  TFile* tfile = TFile::Open(infile);
  TTree* tree = dynamic_cast<TTree*>(tfile->Get("/Event/Sim"));
  WireCellSst::FrameDataSource fds(*tree, "raw");
  int eventNo, runNo;
  tree->SetBranchAddress("eventNo",&eventNo);
  tree->SetBranchAddress("runNo",&runNo);

  TFile* ofile = new TFile(outfile,"recreate");
  int m_run, m_event, m_chid, m_wpid;
  // double m_sump, m_sumn, m_hp, m_hn;
  int m_roi_t, m_roi_dt;
  double m_roi_q, m_ch_rms;
  TTree* m_tree = new TTree("t","ROOT tree for induction pulse");
  m_tree->Branch("run", &m_run, "run/I");
  m_tree->Branch("event", &m_event, "event/I");
  m_tree->Branch("chid", &m_chid, "chid/I");
  m_tree->Branch("wpid", &m_wpid, "wpid/I");
  // m_tree->Branch("sump", &m_sump, "sump/D");
  // m_tree->Branch("sumn", &m_sumn, "sumn/D");
  // m_tree->Branch("hp", &m_hp, "hp/D");
  // m_tree->Branch("hn", &m_hn, "hn/D");
  m_tree->Branch("roi_t", &m_roi_t, "roi_t/I");
  m_tree->Branch("roi_dt", &m_roi_dt, "roi_dt/I");
  m_tree->Branch("roi_q", &m_roi_q, "roi_q/D");
  m_tree->Branch("ch_rms", &m_ch_rms, "ch_rms/D");


  TSpectrum* ss = new TSpectrum(10);
  // Loop over frames (aka "events")
  size_t nframes = fds.size();
  // nframes = 10;
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
      if (pid==2){
        int nbins = trace.charge.size();
        TH1F* hsignal = dynamic_cast<TH1F*>(trace.hCharge);
        // int nfound = ss->Search(hsignal,2,"",0.2);
        // if(nfound>0) {
        //   ofile->cd();
        //   std::string chname("ch");
        //   chname += std::to_string(chid);
        //   hsignal->SetName(chname.c_str());
        //   hsignal->SetTitle(chname.c_str());
        //   hsignal->Write();
        // }
        // double* xpeaks = ss->GetPositionX();

        m_run = runNo;
        m_event = eventNo;
        m_wpid = pid;
        m_chid = chid;
        std::pair<double, double> temp = cal_mean_rms(hsignal); // mean & rms
        m_ch_rms = temp.second; 
        std::vector<BinRange> RoiList;
        for(int i=0; i<nbins; i++){
          if (trace.charge.at(i) - temp.first > 5*temp.second){
            if (RoiList.empty()) RoiList.push_back(std::make_pair(i,i));
            else{

              if (RoiList.back().second + 1 == i){
                RoiList.back().second = i;
              }
              else{
                RoiList.push_back(std::make_pair(i,i));
              }

            }
          }
        }
        // print all ROIs
        for(auto& rng: RoiList){
          m_roi_t = rng.first;
          m_roi_dt = rng.second - rng.first;
          m_roi_q = 0;
          for(int i=rng.first; i<=rng.second; i++){
              m_roi_q += trace.charge.at(i) - temp.first;
          }
          m_tree->Fill();
          // std::cerr << "[wgu] ch: " << chid << " ROI: " << rng.first << " " << rng.second << " charge: " << m_roi_q << std::endl;
        }


        // int maxbin = hsignal->GetMaximumBin();
        // int minbin = minBin(hsignal, maxbin, maxbin+50);
        // if (minbin!=maxbin
        //     && hsignal->GetBinContent(maxbin) - temp.first > 6*temp.second
        //     && hsignal->GetBinContent(minbin) - temp.first < -6*temp.second){

        //     int crx = crossingBin(hsignal, maxbin, minbin, temp.first);
        //     // std::cerr << "[wgu] ch: " << chid << " maxbin: " << maxbin << " minbin: " << minbin << " crx: " << crx << std::endl;


        //     m_sump = area(hsignal, maxbin-10, crx, temp.first);
        //     m_sumn = area(hsignal, crx, minbin+10, temp.first);
        //     m_hp = hsignal->GetBinContent(maxbin) - temp.first;
        //     m_hn = temp.first - hsignal->GetBinContent(minbin);
        //     m_tree->Fill();
        //     std::string chname("ch"); chname += std::to_string(chid);
        //     hsignal->SetName(chname.c_str());
        //     hsignal->SetTitle(chname.c_str());
        //     hsignal->Write();
        // }
        // for (int j=0;j!=nfound;j++){
        //   double xp = xpeaks[j];
        //   int bin = hsignal->FindBin(xp);
        //   if( hsignal->GetBinContent(bin) - temp.first < 3*temp.second){
        //     std::cerr << "[wgu] ch: " << chid << " peak: " << xp << std::endl;
        //   }

        // }
      }

    }// traces

  }// frames

  m_tree->Write();
  ofile->Close();
  delete ss;
  return 0;

} // main()
