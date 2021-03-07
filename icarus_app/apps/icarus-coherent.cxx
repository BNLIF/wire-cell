#include <iostream>

#include "TFile.h"
#include "TH2F.h"
#include "TString.h"

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
    
    //120
    hu[1] = (TH2F*)file->Get("hu_orig120"); // 1056 chs, from 1152 to 2208
    hv[1] = (TH2F*)file->Get("hv_orig120"); // 5600 chs, from 2400 to 8000
    hw[1] = (TH2F*)file->Get("hw_orig120"); // 5600 chs, from 8128 to 13728
    
    // 110
    hu[2] = (TH2F*)file->Get("hu_orig110"); // 1056 chs, from 0 to 1056
    hv[2] = (TH2F*)file->Get("hv_orig110"); // 5600 chs, from 2400 to 8000
    hw[2] = (TH2F*)file->Get("hw_orig110"); // 5600 chs, from 8128 to 13728
    
    //121
    hu[3] = (TH2F*)file->Get("hu_orig121"); // 1056 chs, from 14976 to 16032
    hv[3] = (TH2F*)file->Get("hv_orig121"); // 5600 chs, from 16192 to 21792
    hw[3] = (TH2F*)file->Get("hw_orig121"); // 5600 chs, from 21984 to 27584
    

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
