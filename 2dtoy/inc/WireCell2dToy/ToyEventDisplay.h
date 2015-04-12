#ifndef WIRECELL2DTOY_TOYEVENTDISPLAY_H
#define WIRECELL2DTOY_TOYEVENTDISPLAY_H
#include "WireCellData/Point.h"
#include "WireCellData/GeomCell.h"

#include "WireCellNav/SliceDataSource.h"
#include "WireCellSst/GeomDataSource.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TString.h"
#include "TGraph.h"

namespace WireCell2dToy {

  class ToyEventDisplay {
  private:
    
    TCanvas *c1;
    TH2F *h1;
    TGraph *g1;
    TGraph *g2;
    
  public:
    ToyEventDisplay();
    virtual ~ToyEventDisplay();
    
    virtual int init(){ init(4.9,6.1,-1.1,1.1);
      return 0;}

    virtual int init(float x_min, float x_max, float y_min, float y_max);
    
    virtual int draw_mc(int flag, WireCell::PointCVector mctruth, TString option);
    
    virtual int draw_slice(WireCell::Slice slice, WireCellSst::GeomDataSource gds,TString option);

    virtual int draw_cells(WireCell::GeomCellSelection cellall ,TString option);
  };

}

#endif
