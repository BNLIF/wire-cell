#include "TH1.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
  TString filename = argv[1];
  TFile *file = new TFile(filename);
  TTree *T = (TTree*)file->Get("Trun");
  int num = T->GetEntries();
  // Float_t elifetime;
  //T->SetBranchAddress("elifetime",&elifetime);
  Int_t runNo, subrunNo, eventNo;
  T->SetBranchAddress("runNo",&runNo);
  T->SetBranchAddress("subRunNo",&subrunNo);
  T->SetBranchAddress("eventNo",&eventNo);
  if (num >0){
    for (Int_t i=0;i!=num; i++){
      T->GetEntry(i);
      std::cout << filename << " " << runNo << " " << subrunNo << " " << eventNo << " "  << i << std::endl;
      //  std::cout << num << " " << runNo << " " << elifetime << std::endl;
    }
  }
}
