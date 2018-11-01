#ifndef WIRECELL_MCSSIM_H
#define WIRECELL_MCSSIM_H

#include "WireCellData/Units.h"
#include "TVector3.h"
#include "TRandom.h"

#include <vector>

namespace WireCell {

  struct MCStrack{
    int N; // number of vertices 
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    
    void clear(){
      N = 0;
      x.clear(); y.clear(); z.clear();
    }
  };

  void RotateUz(TVector3& direction, TVector3& v1);

  void MultiScattSim(MCStrack& atrack, int N, std::vector<double> initpos, std::vector<double> initdir);
  
};

#endif
