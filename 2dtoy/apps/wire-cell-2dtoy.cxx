#include "WireCellSst/GeomDataSource.h"
#include "WireCell2dToy/FrameDataSource.h"

#include <iostream>
using namespace WireCell;
using namespace std;


int main(int argc, char* argv[])
{
  WireCellSst::GeomDataSource gds(argv[1]);
  WireCell2dToy::FrameDataSource fds(10,&gds);

  std::vector<float> ex = gds.extent();
  cerr << "Extent: "
       << " x:" << ex[0]/units::mm << " mm"
       << " y:" << ex[1]/units::m << " m"
       << " z:" << ex[2]/units::m << " m"
       << endl;

  cout << fds.size() << endl;
  
  fds.jump(1);
  WireCell::Frame frame = fds.get();
  cout << frame.traces.size() << endl;
  
  cout << fds.mctruth.size() << endl;

  return 0;
  
} // main()
