#ifndef WIRECELL2DTOY_TOYEVENTDISPLAY_H
#define WIRECELL2DTOY_TOYEVENTDISPLAY_H
#include "WireCellData/Point.h"
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
    
  public:
    ToyEventDisplay();
    virtual ~ToyEventDisplay();
    
    virtual int init();
    
    virtual int draw_mc(int flag, WireCell::PointCVector mctruth, TString option);
    
    
  };

}

#endif
