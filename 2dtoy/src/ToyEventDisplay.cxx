#include "WireCellData/Units.h"
#include "WireCell2dToy/ToyEventDisplay.h"
#include "TGraph.h"
#include <iostream>

using namespace WireCell2dToy;

ToyEventDisplay::ToyEventDisplay(){
  c1 = new TCanvas("ToyMC","ToyMC",800,600);
  c1->Draw();
}

ToyEventDisplay::~ToyEventDisplay(){
  delete c1;
  delete h1;
  delete g1;
}

int ToyEventDisplay::init(){
  h1 = new TH2F("h1","h1",100,4.9,6.1,100,-1.1,1.1);
  h1->SetTitle("Wires and True Hits");
  h1->GetYaxis()->SetNdivisions(506);
  h1->GetXaxis()->SetNdivisions(506);
  h1->SetYTitle("Z (m)");
  h1->SetXTitle("Y (m)");
  return 0;
}

int ToyEventDisplay::draw_mc(int flag, WireCell::PointCVector mctruth, TString option){
  
  if (flag==1){
    h1->Draw(option);
  }else if (flag==2){
    h1->Reset();
    for (int i=0;i!=mctruth.size();i++){
      h1->Fill(mctruth[i].z/units::m,mctruth[i].y/units::m-0.02,int(mctruth[i].charge*10)/10.);
    }
    h1->Draw(option);
  }else if (flag==3){
    g1 = new TGraph();
    for (int i=0;i!=mctruth.size();i++){
      g1->SetPoint(i,mctruth[i].z/units::m,mctruth[i].y/units::m);
    }
    g1->SetMarkerColor(2);
    g1->SetMarkerSize(0.6);
    g1->Draw(option);
    g1->SetMarkerStyle(20);
  }
  
  return 0;
}


