// Created on 10/31/2018
// by Wenqiang Gu (wgu@bnl.gov)
// to compile it standalone: g++ ToyMCS.cxx -o ToyMCS.exe `root-config --cflags --libs`

#include <iostream>
#include <vector>
#include <cmath>
#include "TRandom.h"
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"

namespace MCS{
  double mrad = 1E-3;
  double cm = 1E-2;
  double m = 1.;
  double stepLen = 2*cm;
  double varTheta = 3*mrad;

  struct track{
    int N; // number of vertices 
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;

    void clear(){
      N = 0;
      x.clear(); y.clear(); z.clear();
    }
  };
};

void RotateUz(TVector3 direction, TVector3& v1){
  TVector3 unit = direction.Unit();
  v1.RotateUz(unit); 
}

// c.f. http://pdg.lbl.gov/2018/reviews/rpp2018-rev-passage-particles-matter.pdf
void MultiScattSim(MCS::track& atrack, int N, std::vector<double> initpos, std::vector<double> initdir){
  atrack.clear();
  double x=0;
  double y=0;
  double z=0; // moving in z'-axis
  double thetaX =0;
  double thetaY =0;
  double varThetaPlane = MCS::varTheta / std::sqrt(2);
  TVector3 direction;
  direction.SetMagThetaPhi(1, initdir[0], initdir[1]);
  while(atrack.N < N){
    double z1= gRandom->Gaus(0,1);
    double z2= gRandom->Gaus(0,1);
    double z3= gRandom->Gaus(0,1);
    double z4= gRandom->Gaus(0,1);

    double dx = MCS::stepLen *thetaX + z1*MCS::stepLen*varThetaPlane/std::sqrt(12) + z2* MCS::stepLen* varThetaPlane*0.5;
    double dy = MCS::stepLen *thetaY + z3*MCS::stepLen*varThetaPlane/std::sqrt(12) + z4* MCS::stepLen* varThetaPlane*0.5;
    z += MCS::stepLen;
    x += dx;
    y += dy;
    thetaX += z2 * varThetaPlane;
    thetaY += z4 * varThetaPlane;
    atrack.N ++;
    TVector3 v1(x,y,z);
    RotateUz(direction, v1);
    atrack.x.push_back(v1.X() + initpos[0]);
    atrack.y.push_back(v1.Y() + initpos[1]);
    atrack.z.push_back(v1.Z() + initpos[2]);

  }
}

int main(){
  MCS::track atrack;

  TFile* ofile = new TFile("mcs-tracks.root","RECREATE");
  TTree* T = new TTree("T","tracks and vertices");
  T->Branch("N", &atrack.N, "N/I");
  T->Branch("x", &atrack.x);
  T->Branch("y", &atrack.y);
  T->Branch("z", &atrack.z);

  for(int i=0; i<1000/*ntracks*/; i++){
    double x = gRandom->Uniform(2.56*MCS::m);
    double y = gRandom->Uniform(2.3*MCS::m);
    double z = gRandom->Uniform(10.4*MCS::m);
    double cosTheta = gRandom->Uniform(-1,1);
    double phi = gRandom->Uniform(0,2*3.1415926);
    std::vector<double> initpos = {x,y,z};
    std::vector<double> initdir = {acos(cosTheta), phi};

    MultiScattSim(atrack, 100/*nvertices*/, initpos, initdir);
    T->Fill();
  }
  T->Write();
  ofile->Close();
}
