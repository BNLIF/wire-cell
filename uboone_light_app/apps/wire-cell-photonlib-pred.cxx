#include "WireCellSst/GeomDataSource.h"
#include "WireCellData/TPCParams.h"
#include "WireCellData/Singleton.h"

#include "WireCell2dToy/ToyLibPe.h"

#include "TH1F.h"
#include "TH2F.h"

#include "TStopwatch.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

using namespace WireCell;
using namespace std;

int main(int argc, char* argv[]){

  if(argc < 3){
    cerr << "usage: wire-cell-photonlib-pred /path/to/ChannelWireGeometry.txt sim.root ";
    return 1;
  }

  // global switch disabling the reference to histo (i.e. not in list of in-memory objects)
  TH1::AddDirectory(kFALSE);

  int setX =  128, minX = 0, maxX = 256, allX=1;
  
  for(int i=1; i!=argc; i++){
    switch(argv[i][1]){
    case 'a':
      setX = atoi(&argv[i][2]);
      break;
    case 'b':
      minX = atoi(&argv[i][2]);
      break;
    case 'c':
      maxX = atoi(&argv[i][2]);
      break;
    case 'd':
      allX = atoi(&argv[i][2]);
      break;
    }
  }

  // sampling points per dimension;
  //int nBinsX=10, nBinsY=20, nBinsZ=100;
  int nBinsX=1, nBinsY=1, nBinsZ=10;  
  std::cout<<"nBinsX: "<<nBinsX<<"  nBinsY: "<<nBinsY<<"  nBinsZ: "<<nBinsZ<<std::endl;
  
  // detector active volume (from flashmatchalg.fcl)
  float xActVolMin = 0.0, xActVolMax = 256.35;
  float yActVolMin = -116.5, yActVolMax = 116.5;
  float zActVolMin = 0.0, zActVolMax = 1036.8;
  
  // setup x sampling points
  float xSample[nBinsX+1] = {0.};
  float ySample[nBinsY+1] = {0.};
  float zSample[nBinsZ+1] = {0.};
  int a = 0; float incrementX = (xActVolMax-xActVolMin)/nBinsX;
  for(float i=xActVolMin; i<=xActVolMax; i+=incrementX){ xSample[a] = i; a++; }
  a = 0; float incrementY = (yActVolMax-yActVolMin)/nBinsY;
  for(float i=yActVolMin; i<=yActVolMax; i+=incrementY){ ySample[a] = i; a++; }
  a = 0; float incrementZ = (zActVolMax-zActVolMin)/nBinsZ;
  for(float i=zActVolMin; i<=zActVolMax; i+=incrementZ){ zSample[a] = i; a++; }
  
  //int eve_num = atoi(argv[3]);
  int eve_num = 50;

  // set up input file
  const char* root_file = argv[2];
  TFile *file = new TFile(root_file);
  TTree *tout = (TTree*)file->Get("tout");
  vector<int> *trackId = new vector<int>;
  vector<float> *energy = new vector<float>;
  vector<float> *numElectrons = new vector<float>;
  vector<float> *x = new vector<float>;
  vector<float> *y = new vector<float>;
  vector<float> *z = new vector<float>;
  tout->SetBranchAddress("trackId",&trackId);
  tout->SetBranchAddress("energy",&energy);
  tout->SetBranchAddress("numElectrons",&numElectrons);
  tout->SetBranchAddress("x",&x);
  tout->SetBranchAddress("y",&y);
  tout->SetBranchAddress("z",&z);
  //tout->GetEntry(eve_num);
  int nEntries = tout->GetEntries();
  /*  
  // set up output file
  TFile *fout = new TFile("libexp.root","RECREATE");
  TTree *T_exp = new TTree("T_exp","T_exp");
  int index = 0;
  Double_t totalEnergy = 0;
  bool fidVol = true;
  bool atBoundary = true;
  std::vector<double> pe_pred_v;
  std::vector<double> shifted_v;
  T_exp->Branch("index",&index,"index/I");
  T_exp->Branch("totalEnergy",&totalEnergy,"totalEnergy/D");
  T_exp->Branch("pe_pred_v",&pe_pred_v);
  T_exp->Branch("shifted_v",&shifted_v);
  T_exp->Branch("fidVol",&fidVol);
  T_exp->Branch("atBoundary",&atBoundary);

  pe_pred_v.reserve(32);
  shifted_v.reserve(4);
  */
  std::cout<<"zeroth   "<<xSample[0]<<std::endl;
  std::cout<<"first    "<<xSample[1]<<std::endl;
  std::cout<<"second   "<<xSample[2]<<std::endl;

  WireCell2dToy::ToyLibPe simCharge(root_file);
  float shiftX = 0, shiftY = 0, shiftZ = 0;
  TStopwatch tTime;

  TFile** fout = new TFile*[nEntries];
  
  for(int i=0; i<nEntries; i++){
    tout->GetEntry(i);
    if(i>=2) break;
    std::cout<<"===== EVENT "<<i<<" ====="<<std::endl;
    TString fileName;
    fileName.Form("nuE_event%02d.root",i);
    TFile *fout = new TFile(fileName,"RECREATE");
    TTree *T_exp = new TTree("T_exp","T_exp");
    int index = 0;
    Double_t totalEnergy = 0;
    bool fidVol = true;
    bool atBoundary = true;
    std::vector<double> pe_pred_v;
    std::vector<double> shifted_v;
    T_exp->Branch("index",&index,"index/I");
    T_exp->Branch("totalEnergy",&totalEnergy,"totalEnergy/D");
    T_exp->Branch("pe_pred_v",&pe_pred_v);
    T_exp->Branch("shifted_v",&shifted_v);
    T_exp->Branch("fidVol",&fidVol);
    T_exp->Branch("atBoundary",&atBoundary);

    pe_pred_v.reserve(32);
    shifted_v.reserve(4);
        
    int nE = (int)energy->size();
    for(int c=0; c<nE; c++){ totalEnergy += energy->at(c); }

    index = i;

    for(int a=0; a<nBinsX; a++){
      for(int b=0; b<nBinsY; b++){
	for(int c=0; c<nBinsZ; c++){

	  tTime.Start(); 
	  pe_pred_v.clear();
	  shifted_v.clear();
	  std::cout<<"----------------------------------------"<<std::endl;
	  std::cout<<"sampling position:\n x: "<<xSample[a]<<"\n y: "<<ySample[b]<<"\n z: "<<zSample[c]<<std::endl;

	  std::vector<double> original_v = {0.0, 0.0, 0.0, 0.0};
	  original_v = simCharge.getXYZQ(x,y,z,numElectrons);
	  std::cout<<"original position:\n x: "<<original_v[0]<<"\n y: "<<original_v[1]<<"\n z: "<<original_v[2]<<std::endl;
	  
	  shiftX = original_v[0]-xSample[a];
	  shiftY = original_v[1]-ySample[b];
	  shiftZ = original_v[2]-zSample[c];
	  std::cout<<"shift distance:\n x: "<<shiftX<<"\n y: "<<shiftY<<"\n z:"<<shiftZ<<std::endl;
	  
	  shifted_v = simCharge.getShiftedXYZQ(x,y,z,numElectrons,shiftX,shiftY,shiftZ);
	  std::cout<<"weighted average position:\n shifted x="<<shifted_v[0]<<"\n shifted y="<<shifted_v[1]<<"\n shifted z="<<shifted_v[2]<<"\n q="<<shifted_v[3]<<std::endl;
	  
	  std::pair<float,float> shiftedMinMax = simCharge.getShiftedXminXmax(x,shiftX);
	  std::cout<<"minimum extent shifted x: "<<shiftedMinMax.first<<"  maximum extent shifted x:"<<shiftedMinMax.second<<std::endl;
	  
	  fidVol = simCharge.xInFidVol(shiftedMinMax);
	  atBoundary = simCharge.xAtBoundary(shiftedMinMax);
	  std::cout<<"Is cluster inside fiducial volume? "<< fidVol 
		   <<"\nIs cluster at boundary? "<< atBoundary<<std::endl;

	  pe_pred_v= simCharge.fromClusterFindPeDist(x,y,z,numElectrons,shiftX,shiftY,shiftZ);

	  T_exp->Fill(); 
	  tTime.Stop();
	  std::cout<<"iteration took "<<tTime.CpuTime()<<" CPU seconds\n"<<std::endl;
	  tTime.Reset();
	} // zSamp
      } // ySamp
    } // xSamp
    fout->Write();
    fout->Close();
  } // i

} // main
