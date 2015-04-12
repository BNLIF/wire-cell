#include "WireCellData/Units.h"
#include "WireCell2dToy/ToyEventDisplay.h"
#include "TGraph.h"
#include "TLine.h"
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
  delete g2;
}

int ToyEventDisplay::init(){
  h1 = new TH2F("h1","h1",1000,4.9,6.1,1000,-1.1,1.1);
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

int ToyEventDisplay::draw_slice(WireCell::Slice slice, WireCellSst::GeomDataSource gds, TString option){
  
  WireCell::Channel::Group group = slice.group();
  //std::cout << group.size() << std::endl;
  
  for (int i=0;i!=group.size();i++){
    //std::cout << group.at(i).first << std::endl;
    const WireCell::GeomWire *wire = gds.by_channel_segment(group.at(i).first,0);
    //std::cout << wire->point1().y << " " << wire->point1().z << std::endl;
    TLine *l1 = new TLine(wire->point1().z/units::m,wire->point1().y/units::m,
			  wire->point2().z/units::m,wire->point2().y/units::m);
    l1->Draw(option);
  }
  return 0;
}


int ToyEventDisplay::draw_cells(WireCell::GeomCellSelection cellall ,TString option){
  g2 = new TGraph();
  for (int i=0;i!=cellall.size();i++){
    g2->SetPoint(i,cellall[i]->center().z/units::m,cellall[i]->center().y/units::m);
  }
  g2->SetMarkerColor(4);
  g2->SetMarkerSize(0.8);
  g2->Draw(option);
  g2->SetMarkerStyle(21);
  
  return 0;
}
