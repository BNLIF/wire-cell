#include "WireCell2dToy/ToyLightReco.h"
#include "TH1F.h"
#include "TH2F.h"

#include <iostream>
#include <fstream>
#include <ios>

using namespace WireCell;
using namespace std;

int main(int argc, char* argv[])
{
  if (argc < 2) {
    cerr << "usage: wire-cell-uboone /path/to/celltree.root " << endl;
    return 1;
  }
  TH1::AddDirectory(kFALSE);

  string suffix (".root");
  string prefix ("_");
  vector<size_t> positions;
  
  const char* outName = "test.txt";
  ofstream outputFile;
  outputFile.open(outName, ios::out | ios::app);

  const char* bookName = "book.txt";
  ofstream bookFile;
  bookFile.open(bookName, ios::out | ios::app);

  const char* bermName = "berm.txt";
  ofstream bermFile;
  bermFile.open(bermName, ios::out | ios::app);

  string myArray = argv[1];
  const char* root_file = argv[1];

  cout<<myArray<<endl;
  size_t suffixPos = myArray.find(suffix);
  if(suffixPos == std::string::npos){
    cout <<"incorrect input... BREAKING!!!"<<endl;
    return 1;
  }

  size_t prefixPos = myArray.find(prefix,0);      
  while(prefixPos != string::npos){
    positions.push_back(prefixPos);
    prefixPos = myArray.find(prefix,prefixPos+1);
  }

  string fileID;
  for(int b=positions.back()+1; b<suffixPos; b++){
    fileID += myArray[b];
  }
      
  WireCell2dToy::ToyLightReco uboone_flash(root_file);
      
  TFile *file1 = new TFile(root_file);
  TTree *T = (TTree*)file1->Get("/Event/Sim");
  int berm_entries = T->GetEntries();
  int event_no, run_no, subrun_no;
  unsigned int triggerbits;
  vector<unsigned int> *phmax = new vector<unsigned int>;
  vector<unsigned int> *multiplicity = new vector<unsigned int>;
  T->SetBranchStatus("eventNo",1);
  T->SetBranchAddress("eventNo" , &event_no);
  T->SetBranchStatus("runNo",1);
  T->SetBranchAddress("runNo"   , &run_no);
  T->SetBranchStatus("subRunNo",1);
  T->SetBranchAddress("subRunNo", &subrun_no);
  T->SetBranchStatus("triggerBits",1);
  T->SetBranchAddress("triggerBits", &triggerbits);
  T->SetBranchStatus("PHMAX",1);
  T->SetBranchAddress("PHMAX", &phmax);
  T->SetBranchStatus("multiplicity",1);
  T->SetBranchAddress("multiplicity", &multiplicity);
      
  for(int z=0; z<berm_entries; z++){
    //for(int z=0; z<1; z++){
    T->GetEntry(z);
    // if(run_no != 5124) continue;
    // if(subrun_no != 7) continue;
    // if(event_no != 396) continue;

    uboone_flash.load_event_raw(z);
    bool beamspill = false;
    WireCell::OpflashSelection& flashes = uboone_flash.get_flashes();

    unsigned int ph = 10000, mult = 10000;
    if(phmax->size()>0){
      ph = phmax->at(0);
    }
    if(multiplicity->size()>0){
      mult = multiplicity->at(0);
    }

    for(auto it = flashes.begin(); it!=flashes.end(); it++){
      Opflash *flash = (*it);
      double time = flash->get_time();
      double pe = flash->get_total_PE();
      double low = 0, high = 0;
      int type = flash->get_type();
      if((triggerbits>>11) & 1U) { low = 3.0; high = 5.0; }// bnb  
      if((triggerbits>>9) & 1U) { low = 3.45; high = 5.45; } // extbnb
     
      if(time > low && time < high && type == 2){ beamspill = true; }
      bermFile<<run_no<<" "<<subrun_no<<" "<<event_no<<" "<<ph<<" "<<mult<<" "<<time<<" "<<pe<<endl;
    }

    if(beamspill){
      cout<<"PASS"<<endl;
      outputFile<<run_no<<" "<<subrun_no<<" "<<event_no<<endl;
    }
    else{
      cout<<"FAIL"<<endl;
    }
    
    uboone_flash.clear_flashes();
    bookFile<<run_no<<" "<<subrun_no<<" "<<event_no<<" "<<fileID<<" "<<z<<endl;
  }
  T->Delete();
  file1->Close();
  file1->Delete();
  
  outputFile.close();
  bookFile.close();
  bermFile.close();
  return 1;
}
