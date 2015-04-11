#ifndef WIRECELLSST_FRAMEDATASOURCE_H
#define  WIRECELLSST_FRAMEDATASOURCE_H

#include "WireCellNav/FrameDataSource.h"

namespace WireCell2dToy {

    class FrameDataSource : public WireCell::FrameDataSource {
      int Nevent;
      
      public:
        FrameDataSource(int nevents);
	virtual ~FrameDataSource(); 


	virtual int size() const;
	virtual int jump(int frame_number); 
	
    };

}

#endif
