#include <iostream>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TH2F.h"
#include "TTree.h"
#include "TString.h"
#include "TMatrixD.h"
#include "TClonesArray.h"

using namespace std;

void restore_baseline(TH1F *h, int nbin=4096);

int main(int argc, char* argv[])
{
  if (argc < 4) {
    cerr << "usage: icarus-coherent #filename #outfile #event" << endl;
  }else{
    TString filename = argv[1];
    TFile *file = new TFile(filename);

    TString outfile = argv[2];
    int nticks = 4096;

    int neve = atoi(argv[3]);

    const char* tpath = "/Event/Sim";
    TTree* tree = dynamic_cast<TTree*>(file->Get(tpath));
    
    tree->SetBranchStatus("*",0);

    int event_no, run_no, subrun_no;
    tree->SetBranchStatus("eventNo",1);
    tree->SetBranchAddress("eventNo" , &event_no);
    tree->SetBranchStatus("runNo",1);
    tree->SetBranchAddress("runNo"   , &run_no);
    tree->SetBranchStatus("subRunNo",1);
    tree->SetBranchAddress("subRunNo", &subrun_no);
    
    std::vector<int> *channelid = new std::vector<int>;
    TClonesArray* esignal = new TClonesArray;
          
    tree->SetBranchStatus("raw_channelId",1);
    tree->SetBranchAddress("raw_channelId", &channelid);
    tree->SetBranchStatus("raw_wf",1);
    tree->SetBranchAddress("raw_wf", &esignal);

    int n_u1 = 33; int start_u1[2] = {0,13824};      
    int n_u2 = 33; int start_u2[2] = {1152,14976};
    int n_v = 175; int start_v[2] = {2400,16192};
    int n_w = 175; int start_w[2] = {8128,21984};
    
    // save 
    if (neve >= tree->GetEntries()) return 0;

    tree->GetEntry(neve);

    int nchannels = channelid->size();
    std::map<int, TH1F*> map_ch_hist;
    
    for (size_t ind=0; ind < nchannels; ++ind) {
      TH1F* signal = dynamic_cast<TH1F*>(esignal->At(ind));
      restore_baseline(signal);
      map_ch_hist[channelid->at(ind)] = signal;
      //	std::cout << ind << " " << signal->GetBinContent(1000) << std::endl;
    }
    
    std::vector<float> temp_v(32);
    size_t mid = 16;

    
    
    TH2F *hu1 = new TH2F("hu1","hu1",66,0,66,4096,0,4096);
    TH2F *hu2 = new TH2F("hu2","hu2",66,0,66,4096,0,4096);
    TH2F *hv = new TH2F("hv","hv",350,0,350,4096,0,4096);
    TH2F *hw = new TH2F("hw","hw",350,0,350,4096,0,4096);

    TH2F *hu1_m = new TH2F("hu1_m","hu1_m",66,0,66,4096,0,4096);
    TH2F *hu2_m = new TH2F("hu2_m","hu2_m",66,0,66,4096,0,4096);
    TH2F *hv_m = new TH2F("hv_m","hv_m",350,0,350,4096,0,4096);
    TH2F *hw_m = new TH2F("hw_m","hw_m",350,0,350,4096,0,4096);
    
    
    
    // TH2F *hu1_orig = new TH2F("hu1_orig","hu1_orig",2112,0,2112,4096,0,4096);
    // TH2F *hu2_orig = new TH2F("hu2_orig","hu2_orig",2112,0,2112,4096,0,4096);
    // TH2F *hv_orig = new TH2F("hv_orig","hv_orig",11200,0,11200,4096,0,4096);
    // TH2F *hw_orig = new TH2F("hw_orig","hw_orig",11200,0,11200,4096,0,4096);

    // TH2F *hu1_raw = new TH2F("hu1_raw","hu1_raw",2112,0,2112,4096,0,4096);
    // TH2F *hu2_raw = new TH2F("hu2_raw","hu2_raw",2112,0,2112,4096,0,4096);
    // TH2F *hv_raw = new TH2F("hv_raw","hv_raw",11200,0,11200,4096,0,4096);
    // TH2F *hw_raw = new TH2F("hw_raw","hw_raw",11200,0,11200,4096,0,4096);
    
    TH1F *htemp = new TH1F("htemp","htemp", nticks,0,nticks);
    TH1 *h_real = 0;
    TH1 *h_imag = 0;

    TH1 *h_real1 = 0;
    TH1 *h_imag1 = 0;
    
    // save the middle one ...
    
    // U plane ...
    for (Int_t j=0;j!=2;j++){
      //
      for (Int_t k=0;k!=n_u1;k++){
	// look ticks ...
	htemp->Clear();
	for(Int_t k1 = 0; k1 != 4096; k1 ++){
	  // ... chs
	  for (Int_t q = 0; q != 32; q++){
	    temp_v.at(q) = map_ch_hist[start_u1[j]+32*k+q]->GetBinContent(k1+1);
	    // hu1_orig->SetBinContent(32*n_u1*j+32*k+q+1,k1+1,temp_v.at(q));
	  }
	  
	  std::nth_element(temp_v.begin(),temp_v.begin()+mid, temp_v.end());
	  htemp->SetBinContent(k1+1,temp_v.at(mid));
	  //	  hu1->SetBinContent(j*n_u1+k+1,k1+1,temp_v.at(mid));
	  // for (Int_t q = 0; q!=32; q++){
	  //   hu1_raw->SetBinContent(32*n_u1*j+32*k+q+1,k1+1, map_ch_hist[start_u1[j]+32*k+q]->GetBinContent(k1+1) - temp_v.at(mid));
	  // }
	  
	} // loop over k1 ...

	
	h_real = htemp->FFT(0,"RE");
	h_imag = htemp->FFT(0,"IM");

	for (Int_t q1 = 0; q1!=2048;q1++){
	  hu1->SetBinContent(j*n_u1+k+1,2*q1+1,h_real->GetBinContent(q1+1));
	  hu1->SetBinContent(j*n_u1+k+1,2*q1+1+1,h_imag->GetBinContent(q1+1));
	}

	delete h_real;
	delete h_imag;
	
	for (Int_t q = 0; q != 32; q++){
	  h_real1 =  map_ch_hist[start_u1[j]+32*k+q]->FFT(0,"RE");
	  h_imag1 =  map_ch_hist[start_u1[j]+32*k+q]->FFT(0,"IM");

	  //	  std::cout << j*n_u1 + k + 1 << std::endl;
	  
	  for (Int_t q1 =0 ; q1!=2048; q1++){ 
	    hu1_m->SetBinContent(j*n_u1+k+1,2*q1+1, hu1_m->GetBinContent(j*n_u1+k+1,2*q1+1) + pow(h_real1->GetBinContent(q1+1) - hu1->GetBinContent(j*n_u1+k+1,2*q1+1),2));
	    hu1_m->SetBinContent(j*n_u1+k+1,2*q1+1+1, hu1_m->GetBinContent(j*n_u1+k+1,2*q1+1+1) + pow(h_imag1->GetBinContent(q1+1) - hu1->GetBinContent(j*n_u1+k+1,2*q1+1+1),2));
	  }
	  
	  delete h_real1;
	  delete h_imag1;
	}      	
      } // k 
    } //j


    // U2
    for (Int_t j=0;j!=2;j++){
      //
      for (Int_t k=0;k!=n_u2;k++){
	// look ticks ...
	htemp->Clear();
	for(Int_t k1 = 0; k1 != 4096; k1 ++){
	  // ... chs
	  for (Int_t q = 0; q != 32; q++){
	    temp_v.at(q) = map_ch_hist[start_u2[j]+32*k+q]->GetBinContent(k1+1);
	    // hu2_orig->SetBinContent(32*n_u2*j+32*k+q+1,k1+1,temp_v.at(q));
	  }
	  
	  std::nth_element(temp_v.begin(),temp_v.begin()+mid, temp_v.end());
	  htemp->SetBinContent(k1+1,temp_v.at(mid));
	  //	  hu2->SetBinContent(j*n_u2+k+1,k1+1,temp_v.at(mid));
	  // for (Int_t q = 0; q!=32; q++){
	  //   hu2_raw->SetBinContent(32*n_u2*j+32*k+q+1,k1+1, map_ch_hist[start_u2[j]+32*k+q]->GetBinContent(k1+1) - temp_v.at(mid));
	  // }
	  
	} // loop over k1 ...

	
	h_real = htemp->FFT(0,"RE");
	h_imag = htemp->FFT(0,"IM");

	for (Int_t q1 = 0; q1!=2048;q1++){
	  hu2->SetBinContent(j*n_u2+k+1,2*q1+1,h_real->GetBinContent(q1+1));
	  hu2->SetBinContent(j*n_u2+k+1,2*q1+1+1,h_imag->GetBinContent(q1+1));
	}

	delete h_real;
	delete h_imag;
	
	for (Int_t q = 0; q != 32; q++){
	  h_real1 =  map_ch_hist[start_u2[j]+32*k+q]->FFT(0,"RE");
	  h_imag1 =  map_ch_hist[start_u2[j]+32*k+q]->FFT(0,"IM");

	  //	  std::cout << j*n_u2 + k + 1 << std::endl;
	  
	  for (Int_t q1 =0 ; q1!=2048; q1++){ 
	    hu2_m->SetBinContent(j*n_u2+k+1,2*q1+1, hu2_m->GetBinContent(j*n_u2+k+1,2*q1+1) + pow(h_real1->GetBinContent(q1+1) - hu2->GetBinContent(j*n_u2+k+1,2*q1+1),2));
	    hu2_m->SetBinContent(j*n_u2+k+1,2*q1+1+1, hu2_m->GetBinContent(j*n_u2+k+1,2*q1+1+1) + pow(h_imag1->GetBinContent(q1+1) - hu2->GetBinContent(j*n_u2+k+1,2*q1+1+1),2));
	  }
	  
	  delete h_real1;
	  delete h_imag1;
	}      	
      } // k 
    } //j

    // V
    for (Int_t j=0;j!=2;j++){
      //
      for (Int_t k=0;k!=n_v;k++){
	// look ticks ...
	htemp->Clear();
	for(Int_t k1 = 0; k1 != 4096; k1 ++){
	  // ... chs
	  for (Int_t q = 0; q != 32; q++){
	    temp_v.at(q) = map_ch_hist[start_v[j]+32*k+q]->GetBinContent(k1+1);
	    // hv_orig->SetBinContent(32*n_v*j+32*k+q+1,k1+1,temp_v.at(q));
	  }
	  
	  std::nth_element(temp_v.begin(),temp_v.begin()+mid, temp_v.end());
	  htemp->SetBinContent(k1+1,temp_v.at(mid));
	  //	  hv->SetBinContent(j*n_v+k+1,k1+1,temp_v.at(mid));
	  // for (Int_t q = 0; q!=32; q++){
	  //   hv_raw->SetBinContent(32*n_v*j+32*k+q+1,k1+1, map_ch_hist[start_v[j]+32*k+q]->GetBinContent(k1+1) - temp_v.at(mid));
	  // }
	  
	} // loop over k1 ...

	
	h_real = htemp->FFT(0,"RE");
	h_imag = htemp->FFT(0,"IM");

	for (Int_t q1 = 0; q1!=2048;q1++){
	  hv->SetBinContent(j*n_v+k+1,2*q1+1,h_real->GetBinContent(q1+1));
	  hv->SetBinContent(j*n_v+k+1,2*q1+1+1,h_imag->GetBinContent(q1+1));
	}

	delete h_real;
	delete h_imag;
	
	for (Int_t q = 0; q != 32; q++){
	  h_real1 =  map_ch_hist[start_v[j]+32*k+q]->FFT(0,"RE");
	  h_imag1 =  map_ch_hist[start_v[j]+32*k+q]->FFT(0,"IM");

	  //	  std::cout << j*n_v + k + 1 << std::endl;
	  
	  for (Int_t q1 =0 ; q1!=2048; q1++){ 
	    hv_m->SetBinContent(j*n_v+k+1,2*q1+1, hv_m->GetBinContent(j*n_v+k+1,2*q1+1) + pow(h_real1->GetBinContent(q1+1) - hv->GetBinContent(j*n_v+k+1,2*q1+1),2));
	    hv_m->SetBinContent(j*n_v+k+1,2*q1+1+1, hv_m->GetBinContent(j*n_v+k+1,2*q1+1+1) + pow(h_imag1->GetBinContent(q1+1) - hv->GetBinContent(j*n_v+k+1,2*q1+1+1),2));
	  }
	  
	  delete h_real1;
	  delete h_imag1;
	}      	
      } // k 
    } //j

    // W
    for (Int_t j=0;j!=2;j++){
      //
      for (Int_t k=0;k!=n_w;k++){
	// look ticks ...
	htemp->Clear();
	for(Int_t k1 = 0; k1 != 4096; k1 ++){
	  // ... chs
	  for (Int_t q = 0; q != 32; q++){
	    temp_v.at(q) = map_ch_hist[start_w[j]+32*k+q]->GetBinContent(k1+1);
	    // hw_orig->SetBinContent(32*n_w*j+32*k+q+1,k1+1,temp_v.at(q));
	  }
	  
	  std::nth_element(temp_v.begin(),temp_v.begin()+mid, temp_v.end());
	  htemp->SetBinContent(k1+1,temp_v.at(mid));
	  //	  hw->SetBinContent(j*n_w+k+1,k1+1,temp_v.at(mid));
	  // for (Int_t q = 0; q!=32; q++){
	  //   hw_raw->SetBinContent(32*n_w*j+32*k+q+1,k1+1, map_ch_hist[start_w[j]+32*k+q]->GetBinContent(k1+1) - temp_v.at(mid));
	  // }
	  
	} // loop over k1 ...

	
	h_real = htemp->FFT(0,"RE");
	h_imag = htemp->FFT(0,"IM");

	for (Int_t q1 = 0; q1!=2048;q1++){
	  hw->SetBinContent(j*n_w+k+1,2*q1+1,h_real->GetBinContent(q1+1));
	  hw->SetBinContent(j*n_w+k+1,2*q1+1+1,h_imag->GetBinContent(q1+1));
	}

	delete h_real;
	delete h_imag;
	
	for (Int_t q = 0; q != 32; q++){
	  h_real1 =  map_ch_hist[start_w[j]+32*k+q]->FFT(0,"RE");
	  h_imag1 =  map_ch_hist[start_w[j]+32*k+q]->FFT(0,"IM");

	  //	  std::cout << j*n_w + k + 1 << std::endl;
	  
	  for (Int_t q1 =0 ; q1!=2048; q1++){ 
	    hw_m->SetBinContent(j*n_w+k+1,2*q1+1, hw_m->GetBinContent(j*n_w+k+1,2*q1+1) + pow(h_real1->GetBinContent(q1+1) - hw->GetBinContent(j*n_w+k+1,2*q1+1),2));
	    hw_m->SetBinContent(j*n_w+k+1,2*q1+1+1, hw_m->GetBinContent(j*n_w+k+1,2*q1+1+1) + pow(h_imag1->GetBinContent(q1+1) - hw->GetBinContent(j*n_w+k+1,2*q1+1+1),2));
	  }
	  
	  delete h_real1;
	  delete h_imag1;
	}      	
      } // k 
    } //j
   
    
    
    outfile += Form("_%d.root",neve);
    TFile *file2 = new TFile(outfile,"RECREATE");
    hu1->SetDirectory(file2);
    hu2->SetDirectory(file2);
    hv->SetDirectory(file2);
    hw->SetDirectory(file2);

    hu1_m->Scale(1./32.);
    hu2_m->Scale(1./32.);
    hv_m->Scale(1./32.);
    hw_m->Scale(1./32.);
    
    hu1_m->SetDirectory(file2);
    hu2_m->SetDirectory(file2);
    hv_m->SetDirectory(file2);
    hw_m->SetDirectory(file2);
    
    // hu1_orig->SetDirectory(file2);
    // hu2_orig->SetDirectory(file2);
    // hv_orig->SetDirectory(file2);
    // hw_orig->SetDirectory(file2);

    // hu1_raw->SetDirectory(file2);
    // hu2_raw->SetDirectory(file2);
    // hv_raw->SetDirectory(file2);
    // hw_raw->SetDirectory(file2);
    
    
    file2->Write();
    file2->Close();
    
    
  }

  return 0;
}


void restore_baseline(TH1F *htemp, int nbin){
  //correct baseline 
  double max = htemp->GetMaximum();
  double min = htemp->GetMinimum();
  int nbin_b = max - min;
  if (nbin_b <=0) nbin_b = 1;
  TH1F *h1 = new TH1F("h1","h1",nbin_b,min,max);
  for (int j=0;j!=nbin;j++){
    //    if (j%100==0) std::cout << j << " " << nbin << std::endl;
    if (htemp->GetBinContent(j+1)!=0)
      h1->Fill(htemp->GetBinContent(j+1));
  }
  float ped = h1->GetMaximumBin()*(max-min)/(nbin_b*1.) + min;
  float ave=0,ncount = 0;
  for (int j=0;j!=nbin;j++){
    if (fabs(htemp->GetBinContent(j+1)-ped)<400 && htemp->GetBinContent(j+1)!=0){
      ave +=htemp->GetBinContent(j+1);
      ncount ++;
      }
    }
    if (ncount==0) ncount=1;
    ave = ave/ncount;
    
    for (int j=0;j!=nbin;j++){
      double content = htemp->GetBinContent(j+1);
      content -= ave;
      if (htemp->GetBinContent(j+1)!=0)
	htemp->SetBinContent(j+1,content);
    }
    delete h1;
}
