#ifndef WIRECELLSST_FRAMEDATASOURCE_H
#define  WIRECELLSST_FRAMEDATASOURCE_H

#include "WireCellNav/FrameDataSource.h"
#include "WireCellSst/GeomDataSource.h"

namespace WireCell2dToy {
  
  class FrameDataSource : public WireCell::FrameDataSource {
    int Nevent;
    WireCellSst::GeomDataSource *gds;

  public:
    FrameDataSource(int nevents, WireCellSst::GeomDataSource *gds1);
    virtual ~FrameDataSource(); 
    
    
    virtual int size() const;
    virtual int jump(int frame_number); 
    
  };
  
}

#endif
