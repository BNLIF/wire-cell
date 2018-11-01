#include "WireCellMCSSim/MCSSim.h"

using namespace WireCell;

void WireCell::RotateUz(TVector3& direction, TVector3& v1){
  TVector3 unit = direction.Unit();
  v1.RotateUz(unit); 
}
   
void WireCell::MultiScattSim(MCStrack& atrack, int N, std::vector<double> initpos, std::vector<double> initdir){
  atrack.clear();
  double x=0;
  double y=0;
  double z=0; // moving in z'-axis
  double thetaX =0;
  double thetaY =0;

  double stepLen = 2*units::cm;
  double varTheta = 3*units::mrad;
  
  double varThetaPlane = varTheta / std::sqrt(2);
  TVector3 direction;
  direction.SetMagThetaPhi(1, initdir[0], initdir[1]);
  while(atrack.N < N){
    double z1= gRandom->Gaus(0,1);
    double z2= gRandom->Gaus(0,1);
    double z3= gRandom->Gaus(0,1);
    double z4= gRandom->Gaus(0,1);

    double dx = stepLen *thetaX + z1*stepLen*varThetaPlane/std::sqrt(12) + z2* stepLen* varThetaPlane*0.5;
    double dy = stepLen *thetaY + z3*stepLen*varThetaPlane/std::sqrt(12) + z4* stepLen* varThetaPlane*0.5;
    z += stepLen;
    x += dx;
    y += dy;
    thetaX += z2 * varThetaPlane;
    thetaY += z4 * varThetaPlane;
    atrack.N ++;
    TVector3 v1(x,y,z);
    RotateUz(direction, v1);
    atrack.x.push_back(v1.X() + initpos[0]);
    atrack.y.push_back(v1.Y() + initpos[1]);
    atrack.z.push_back(v1.Z() + initpos[2]);

  }
}
