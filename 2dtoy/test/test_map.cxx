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
 

  std::map<int,float> abc;
  abc[1] = 10;
  abc[2] = 20;
  
  cout << abc[3] << endl;
  return 0;
}
