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
#include "WireCellMCSSim/MCSSim.h"

using namespace WireCell;

// namespace MCS{
//   double mrad = 1E-3;
//   double cm = 1E-2;
//   double m = 1.;
//   double stepLen = 2*cm;
//   double varTheta = 3*mrad;
// };



// c.f. http://pdg.lbl.gov/2018/reviews/rpp2018-rev-passage-particles-matter.pdf


int main(){
  MCStrack atrack;

  TFile* ofile = new TFile("mcs-tracks.root","RECREATE");
  TTree* T = new TTree("T","tracks and vertices");
  T->Branch("N", &atrack.N, "N/I");
  T->Branch("x", &atrack.x);
  T->Branch("y", &atrack.y);
  T->Branch("z", &atrack.z);

  for(int i=0; i<1000/*ntracks*/; i++){
    double x = gRandom->Uniform(2.56*units::m);
    double y = gRandom->Uniform(2.3*units::m);
    double z = gRandom->Uniform(10.4*units::m);
    double cosTheta = gRandom->Uniform(-1,1);
    double phi = gRandom->Uniform(0,2*3.1415926);
    std::vector<double> initpos = {x,y,z};
    std::vector<double> initdir = {acos(cosTheta), phi};

    MultiScattSim(atrack, 100/*nvertices*/, initpos, initdir);

    for (size_t j=0; j!=atrack.x.size(); j++){
      atrack.x.at(j) /= units::m;
      atrack.y.at(j) /= units::m;
      atrack.z.at(j) /= units::m;
    }
    
    T->Fill();
  }
  T->Write();
  ofile->Close();
}
