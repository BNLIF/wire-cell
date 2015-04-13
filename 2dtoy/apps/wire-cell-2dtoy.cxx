#include "WireCellSst/GeomDataSource.h"
#include "WireCell2dToy/FrameDataSource.h"
#include "WireCell2dToy/ToyEventDisplay.h"
#include "WireCell2dToy/ToyTiling.h"

#include "WireCellNav/SliceDataSource.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TH1F.h"
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

  WireCell::SliceDataSource sds(fds);
  sds.jump(0);
  WireCell::Slice slice = sds.get();

  WireCell::ToyTiling toytiling(slice,gds);

  GeomCellSelection allcell = toytiling.get_allcell();
  GeomWireSelection allwire = toytiling.get_allwire();
  //  cout << toytiling.wiremap[allwire.at(0)].size() << endl;
  //  cout << toytiling.cellmap[allcell.at(0)].size() << endl;
  //  cout << sds.size() << endl;

  
  TApplication theApp("theApp",&argc,argv);
  theApp.SetReturnFromRun(true);
  WireCell2dToy::ToyEventDisplay display;
  
  gStyle->SetOptStat(0);

  // // display.init(0,10.3698,-2.33/2.,2.33/2.);
  display.init();
  display.draw_mc(1,fds.mctruth,"");
  display.draw_mc(2,fds.mctruth,"TEXTsame");
  
  
  display.draw_slice(slice,gds,"same");
 
  display.draw_cells(toytiling.get_allcell(),"*same");
  display.draw_mc(3,fds.mctruth,"*same");

  theApp.Run();
  

  return 0;
  
} // main()
