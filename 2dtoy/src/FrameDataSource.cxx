#include "WireCell2dToy/FrameDataSource.h"
#include "TRandom.h"

WireCell2dToy::FrameDataSource::FrameDataSource(int nevents,
		WireCellSst::GeomDataSource *gds1   )
    : WireCell::FrameDataSource()
{
  Nevent = nevents;
  gds = gds1;
}
WireCell2dToy::FrameDataSource::~FrameDataSource()
{
}

int WireCell2dToy::FrameDataSource::size() const
{
  return Nevent;
}

int WireCell2dToy::FrameDataSource::jump(int frame_number)
{
  if (frame_number >= Nevent) frame_number = Nevent;

  gRandom->SetSeed(frame_number);
  frame.clear();
  
  //main code to construct a frame


  frame.index = frame_number;
  return frame.index;
}
