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
  if (argc < 3) {
    cerr << "usage: icarus-coherent #filename #outfile" << endl;
  }else{
    TString filename = argv[1];
    TFile *file = new TFile(filename);
    
    TString outfile = argv[2];
    int nticks = 4096;

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

    //std::cout << tree->GetEntries() << std::endl;

    int n_u1 = 33; int start_u1[2] = {0,13824};      
    int n_u2 = 33; int start_u2[2] = {1152,14976};
    int n_v = 175; int start_v[2] = {2400,16192};
    int n_w = 175; int start_w[2] = {8128,21984};

    TH1F *htemp = new TH1F("htemp","htemp", nticks,0,nticks);
    TH1 *h_real = 0;
    TH1 *h_imag = 0;

    TMatrixD mat_u1(4096,4096); //u1
    int count_u1 = 0;
    TMatrixD mat_u2(4096,4096);
    int count_u2 = 0;
    TMatrixD mat_v(4096,4096);
    int count_v = 0;
    TMatrixD mat_w(4096,4096);
    int count_w = 0;
    
    
    for (Int_t i=0;i!=tree->GetEntries();i++){
      tree->GetEntry(i);

     
      int nchannels = channelid->size();

      //      std::cout << nchannels << std::endl;
      std::map<int, TH1F*> map_ch_hist;
      
      for (size_t ind=0; ind < nchannels; ++ind) {
	TH1F* signal = dynamic_cast<TH1F*>(esignal->At(ind));
	restore_baseline(signal);
	map_ch_hist[channelid->at(ind)] = signal;
	//	std::cout << ind << " " << signal->GetBinContent(1000) << std::endl;
      }
      std::vector<float> temp_v(32);
      size_t mid = 16;

      
      // U plane ...
      for (Int_t j=0;j!=2;j++){
	//
	for (Int_t k=0;k!=n_u1;k++){
	  // look ticks ...
	  //std::vector<float> results;
	  htemp->Clear();
	  for(Int_t k1 = 0; k1 != 4096; k1 ++){
	    // ... chs
	    for (Int_t q = 0; q != 32; q++){
	      temp_v.at(q) = map_ch_hist[start_u1[j]+32*k+q]->GetBinContent(k1+1);
	    }
	    std::nth_element(temp_v.begin(),temp_v.begin()+mid, temp_v.end());
	    //	    results.at(k1) = temp_v.at(mid);
	    htemp->SetBinContent(k1+1,temp_v.at(mid));
	  } // loop over k1 ...
	  h_real = htemp->FFT(0,"RE");
	  h_imag = htemp->FFT(0,"IM");


	  for (Int_t q = 0; q != 2048 ;q++){
	    mat_u1(2*q,2*q) += pow(h_real->GetBinContent(q+1),2);
	    mat_u1(2*q+1,2*q+1) += pow(h_imag->GetBinContent(q+1),2);
	    mat_u1(2*q,2*q+1) += h_real->GetBinContent(q+1) * h_imag->GetBinContent(q+1);
	    
	    for (Int_t p=q+1; p!=2048; p++){
	      mat_u1(2*q,2*p) += h_real->GetBinContent(p+1) * h_real->GetBinContent(q+1);
	      mat_u1(2*q+1,2*p) += h_real->GetBinContent(p+1) * h_imag->GetBinContent(q+1);
	      mat_u1(2*q,2*p+1) += h_imag->GetBinContent(p+1) * h_real->GetBinContent(q+1);
	      mat_u1(2*q+1,2*p+1) += h_imag->GetBinContent(p+1) * h_imag->GetBinContent(q+1);	    
	    }
	  }
	  count_u1 ++; 
	  
	  
	  delete h_real;
	  delete h_imag;
	}
      }


      
      // // U 2 plane ...
      // for (Int_t j=0;j!=2;j++){
      // 	//
      // 	for (Int_t k=0;k!=n_u2;k++){
      // 	  // look ticks ...
      // 	  //std::vector<float> results;
      // 	  htemp->Clear();
      // 	  for(Int_t k1 = 0; k1 != 4096; k1 ++){
      // 	    // ... chs
      // 	    for (Int_t q = 0; q != 32; q++){
      // 	      temp_v.at(q) = map_ch_hist[start_u2[j]+32*k+q]->GetBinContent(k1+1);
      // 	    }
      // 	    std::nth_element(temp_v.begin(),temp_v.begin()+mid, temp_v.end());
      // 	    //	    results.at(k1) = temp_v.at(mid);
      // 	    htemp->SetBinContent(k1+1,temp_v.at(mid));
      // 	  } // loop over k1 ...
      // 	  h_real = htemp->FFT(0,"RE");
      // 	  h_imag = htemp->FFT(0,"IM");


      // 	  for (Int_t q = 0; q != 2048 ;q++){
      // 	    mat_u2(2*q,2*q) += pow(h_real->GetBinContent(q+1),2);
      // 	    mat_u2(2*q+1,2*q+1) += pow(h_imag->GetBinContent(q+1),2);
      // 	    mat_u2(2*q,2*q+1) += h_real->GetBinContent(q+1) * h_imag->GetBinContent(q+1);
	    
      // 	    for (Int_t p=q+1; p!=2048; p++){
      // 	      mat_u2(2*q,2*p) += h_real->GetBinContent(p+1) * h_real->GetBinContent(q+1);
      // 	      mat_u2(2*q+1,2*p) += h_real->GetBinContent(p+1) * h_imag->GetBinContent(q+1);
      // 	      mat_u2(2*q,2*p+1) += h_imag->GetBinContent(p+1) * h_real->GetBinContent(q+1);
      // 	      mat_u2(2*q+1,2*p+1) += h_imag->GetBinContent(p+1) * h_imag->GetBinContent(q+1);	    
      // 	    }
      // 	  }
      // 	  count_u2 ++; 
	  
	  
      // 	  delete h_real;
      // 	  delete h_imag;
      // 	}
      // }


      // // V
      // for (Int_t j=0;j!=2;j++){
      // 	//
      // 	for (Int_t k=0;k!=n_v;k++){
      // 	  // look ticks ...
      // 	  //std::vector<float> results;
      // 	  htemp->Clear();
      // 	  for(Int_t k1 = 0; k1 != 4096; k1 ++){
      // 	    // ... chs
      // 	    for (Int_t q = 0; q != 32; q++){
      // 	      temp_v.at(q) = map_ch_hist[start_v[j]+32*k+q]->GetBinContent(k1+1);
      // 	    }
      // 	    std::nth_element(temp_v.begin(),temp_v.begin()+mid, temp_v.end());
      // 	    //	    results.at(k1) = temp_v.at(mid);
      // 	    htemp->SetBinContent(k1+1,temp_v.at(mid));
      // 	  } // loop over k1 ...
      // 	  h_real = htemp->FFT(0,"RE");
      // 	  h_imag = htemp->FFT(0,"IM");


      // 	  for (Int_t q = 0; q != 2048 ;q++){
      // 	    mat_v(2*q,2*q) += pow(h_real->GetBinContent(q+1),2);
      // 	    mat_v(2*q+1,2*q+1) += pow(h_imag->GetBinContent(q+1),2);
      // 	    mat_v(2*q,2*q+1) += h_real->GetBinContent(q+1) * h_imag->GetBinContent(q+1);
	    
      // 	    for (Int_t p=q+1; p!=2048; p++){
      // 	      mat_v(2*q,2*p) += h_real->GetBinContent(p+1) * h_real->GetBinContent(q+1);
      // 	      mat_v(2*q+1,2*p) += h_real->GetBinContent(p+1) * h_imag->GetBinContent(q+1);
      // 	      mat_v(2*q,2*p+1) += h_imag->GetBinContent(p+1) * h_real->GetBinContent(q+1);
      // 	      mat_v(2*q+1,2*p+1) += h_imag->GetBinContent(p+1) * h_imag->GetBinContent(q+1);	    
      // 	    }
      // 	  }
      // 	  count_v ++; 
	  
	  
      // 	  delete h_real;
      // 	  delete h_imag;
      // 	}
      // }

      // // W
      // for (Int_t j=0;j!=2;j++){
      // 	//
      // 	for (Int_t k=0;k!=n_w;k++){
      // 	  // look ticks ...
      // 	  //std::vector<float> results;
      // 	  htemp->Clear();
      // 	  for(Int_t k1 = 0; k1 != 4096; k1 ++){
      // 	    // ... chs
      // 	    for (Int_t q = 0; q != 32; q++){
      // 	      temp_v.at(q) = map_ch_hist[start_w[j]+32*k+q]->GetBinContent(k1+1);
      // 	    }
      // 	    std::nth_element(temp_v.begin(),temp_v.begin()+mid, temp_v.end());
      // 	    //	    results.at(k1) = temp_v.at(mid);
      // 	    htemp->SetBinContent(k1+1,temp_v.at(mid));
      // 	  } // loop over k1 ...
      // 	  h_real = htemp->FFT(0,"RE");
      // 	  h_imag = htemp->FFT(0,"IM");


      // 	  for (Int_t q = 0; q != 2048 ;q++){
      // 	    mat_w(2*q,2*q) += pow(h_real->GetBinContent(q+1),2);
      // 	    mat_w(2*q+1,2*q+1) += pow(h_imag->GetBinContent(q+1),2);
      // 	    mat_w(2*q,2*q+1) += h_real->GetBinContent(q+1) * h_imag->GetBinContent(q+1);
      // 	    for (Int_t p=q+1; p!=2048; p++){
      // 	      mat_w(2*q,2*p) += h_real->GetBinContent(p+1) * h_real->GetBinContent(q+1);
      // 	      mat_w(2*q+1,2*p) += h_real->GetBinContent(p+1) * h_imag->GetBinContent(q+1);
      // 	      mat_w(2*q,2*p+1) += h_imag->GetBinContent(p+1) * h_real->GetBinContent(q+1);
      // 	      mat_w(2*q+1,2*p+1) += h_imag->GetBinContent(p+1) * h_imag->GetBinContent(q+1);	    
      // 	    }
      // 	  }
      // 	  count_w ++; 
	  
	  
      // 	  delete h_real;
      // 	  delete h_imag;
      // 	}
      // }
      
      
    }


    TFile *file1 = new TFile(outfile,"RECREATE");
    mat_u1.Write("mat_u1");
    mat_u2.Write("mat_u2");
    mat_v.Write("mat_v");
    mat_w.Write("mat_w");
    TTree *T = new TTree("T","T");
    T->SetDirectory(file1);
    T->Branch("count_u1",&count_u1,"count_u1/I");
    T->Branch("count_u2",&count_u2,"count_u2/I");
    T->Branch("count_v",&count_v,"count_v/I");
    T->Branch("count_w",&count_w,"count_w/I");
    T->Fill();
    file1->Write();
    file1->Close();
    
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
