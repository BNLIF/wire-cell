#include <iostream>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TH2F.h"
#include "TString.h"
#include "TMatrixD.h"

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

    TH2F **hu = new TH2F*[4];
    TH2F **hv = new TH2F*[4];
    TH2F **hw = new TH2F*[4];
    
    // 111
    hu[0] = (TH2F*)file->Get("hu_orig111"); // 1056 chs, from 13824 to 14880
    hv[0] = (TH2F*)file->Get("hv_orig111"); // 5600 chs, from 16192 to 21792
    hw[0] = (TH2F*)file->Get("hw_orig111"); // 5600 chs, from 21984 to 27584
    
    // 120
    hu[1] = (TH2F*)file->Get("hu_orig120"); // 1056 chs, from 1152 to 2208
    hv[1] = (TH2F*)file->Get("hv_orig120"); // 5600 chs, from 2400 to 8000
    hw[1] = (TH2F*)file->Get("hw_orig120"); // 5600 chs, from 8128 to 13728
    
    // 110
    hu[2] = (TH2F*)file->Get("hu_orig110"); // 1056 chs, from 0 to 1056
    hv[2] = (TH2F*)file->Get("hv_orig110"); // 5600 chs, from 2400 to 8000
    hw[2] = (TH2F*)file->Get("hw_orig110"); // 5600 chs, from 8128 to 13728
    
    // 121
    hu[3] = (TH2F*)file->Get("hu_orig121"); // 1056 chs, from 14976 to 16032
    hv[3] = (TH2F*)file->Get("hv_orig121"); // 5600 chs, from 16192 to 21792
    hw[3] = (TH2F*)file->Get("hw_orig121"); // 5600 chs, from 21984 to 27584

    int n_u = 33;
    int n_v = 175;
    int n_w = 175;

    TH2F **hu1 = new TH2F*[4];
    TH2F **hv1 = new TH2F*[4];
    TH2F **hw1 = new TH2F*[4];
    for (Int_t i=0;i!=4;i++){
      hu1[i] = new TH2F(Form("hu_%d",i),Form("hu_%d",i),n_u,0-0.5,n_u-0.5,4096,0,4096);
      hv1[i] = new TH2F(Form("hv_%d",i),Form("hv_%d",i),n_v,0-0.5,n_v-0.5,4096,0,4096);
      hw1[i] = new TH2F(Form("hw_%d",i),Form("hw_%d",i),n_w,0-0.5,n_w-0.5,4096,0,4096);
    }
    
    std::vector<float> temp_v(32);
    size_t mid = 16;
    TH1F *ht = new TH1F("ht","ht",4096,0,4096);
    // calculate the baseline ...
    for (Int_t i=0;i!=4;i++){
      // u
      for (Int_t j=0;j!=1056;j++){
	for (Int_t k=0;k!=4096;k++){
	  ht->SetBinContent(k+1,hu[i]->GetBinContent(j+1,k+1));
	}
	restore_baseline(ht);
	for (Int_t k=0;k!=4096;k++){
	  hu[i]->SetBinContent(j+1,k+1,ht->GetBinContent(k+1));
	}
      }

      // calculate medium waveform ...
      for (Int_t j=0;j!=n_u;j++){
	for (Int_t k=0;k!=4096;k++){
	  for (Int_t q = 0; q!=32;q++){
	    temp_v.at(q) = hu[i]->GetBinContent(j*32+q+1,k+1);
	  }
	  // calculate medium ...
	  std::nth_element(temp_v.begin(),temp_v.begin()+mid, temp_v.end());
	  hu1[i]->SetBinContent(j+1,k+1,temp_v.at(mid));
	}
      }
      
      // v
      for (Int_t j=0;j!=5600;j++){
	for (Int_t k=0;k!=4096;k++){
	  ht->SetBinContent(k+1,hv[i]->GetBinContent(j+1,k+1));
	}
	restore_baseline(ht);
	for (Int_t k=0;k!=4096;k++){
	  hv[i]->SetBinContent(j+1,k+1,ht->GetBinContent(k+1));
	}
      }
      for (Int_t j=0;j!=n_v;j++){
	for (Int_t k=0;k!=4096;k++){
	  for (Int_t q = 0; q!=32;q++){
	    temp_v.at(q) = hv[i]->GetBinContent(j*32+q+1,k+1);
	  }
	  // calculate medium ...
	  std::nth_element(temp_v.begin(),temp_v.begin()+mid, temp_v.end());
	  hv1[i]->SetBinContent(j+1,k+1,temp_v.at(mid));
	}
      }
      
      
      // w
      for (Int_t j=0;j!=5600;j++){
	for (Int_t k=0;k!=4096;k++){
	  ht->SetBinContent(k+1,hw[i]->GetBinContent(j+1,k+1));
	}
	restore_baseline(ht);
	for (Int_t k=0;k!=4096;k++){
	  hw[i]->SetBinContent(j+1,k+1,ht->GetBinContent(k+1));
	}
      }

      for (Int_t j=0;j!=n_w;j++){
	for (Int_t k=0;k!=4096;k++){
	  for (Int_t q = 0; q!=32;q++){
	    temp_v.at(q) = hw[i]->GetBinContent(j*32+q+1,k+1);
	  }
	  // calculate medium ...
	  std::nth_element(temp_v.begin(),temp_v.begin()+mid, temp_v.end());
	  hw1[i]->SetBinContent(j+1,k+1,temp_v.at(mid));
	}
      }
    }

    TH1 *h_real = 0;
    TH1 *h_imag = 0;

    // all U plane ... 1056 * 4 ...
    TMatrixD mat(4096,4096);
    Int_t count = 0;
    // Now look at the covariance matrix ... 
    for (Int_t i=0;i!=4;i++){
      for (Int_t j=0;j!=n_u;j++){
	for (Int_t k=0;k!=4096;k++){
	  ht->SetBinContent(k+1,hu1[i]->GetBinContent(j+1,k+1));
	}
	h_real = ht->FFT(0,"RE");
	h_imag = ht->FFT(0,"IM");

	for (Int_t q = 0; q != 2048 ;q++){
	  mat(2*q,2*q) += pow(h_real->GetBinContent(q+1),2);
	  mat(2*q+1,2*q+1) += pow(h_imag->GetBinContent(q+1),2);
	  mat(2*q,2*q+1) += h_real->GetBinContent(q+1) * h_imag->GetBinContent(q+1);
	  
	  for (Int_t p=q+1; p!=2048; p++){
	    mat(2*q,2*p) = h_real->GetBinContent(p+1) * h_real->GetBinContent(q+1);
	    mat(2*q+1,2*p) = h_real->GetBinContent(p+1) * h_imag->GetBinContent(q+1);
	    mat(2*q,2*p+1) = h_imag->GetBinContent(p+1) * h_real->GetBinContent(q+1);
	    mat(2*q+1,2*p+1) = h_imag->GetBinContent(p+1) * h_imag->GetBinContent(q+1);	    
	  }
	}
	count ++; 
	delete h_real;
	delete h_imag;
      }
    }

    // std::cout << count << std::endl;
    mat *= 1./count;
    for (Int_t q = 0; q!=4096;q++){
      for (Int_t p = q+1; p!=4096;p++){
	mat(p,q) = mat(q,p);
      }
    }
    
    
    TFile *file2 = new TFile(outfile,"RECREATE");
    hu[0]->SetDirectory(file2);
    hv[0]->SetDirectory(file2);
    hw[0]->SetDirectory(file2);
    
    hu[1]->SetDirectory(file2);
    hv[1]->SetDirectory(file2);
    hw[1]->SetDirectory(file2);
    
    hu[2]->SetDirectory(file2);
    hv[2]->SetDirectory(file2);
    hw[2]->SetDirectory(file2);
    
    hu[3]->SetDirectory(file2);
    hv[3]->SetDirectory(file2);
    hw[3]->SetDirectory(file2);

    for (Int_t i=0;i!=4;i++){
      hu1[i]->SetDirectory(file2);
      hv1[i]->SetDirectory(file2);
      hw1[i]->SetDirectory(file2);
    }
    mat.Write("mat");
    
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
