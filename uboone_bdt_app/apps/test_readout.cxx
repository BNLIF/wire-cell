#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <string>

#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TH1F.h"


#include "WCPLEEANA/master_cov_matrix.h"
#include "WCPLEEANA/bayes.h"

using namespace std;
using namespace LEEana;

int main( int argc, char** argv )
{
  //  double num = atoi(argv[1]);
  
  CovMatrix cov;

  // cov.print_ch_info();
  // cov.print_filetype_info();
  // cov.print_systematics();
  // cov.print_matrix();
  // int ch = 2;
  // std::cout << cov.get_sys_xs_flux(ch) << " " << cov.get_sys_det(ch) << " " << cov.get_sys_add(ch).first << " " << cov.get_sys_add(ch).second << " " << cov.get_sys_mc_same(ch) << std::endl;
  // TFile *file = new TFile("temp.root","RECREATE");
  // TMatrixD* mat_collapse = cov.get_mat_collapse();
  // mat_collapse->Write("mat_collapse");
  // file->Write();
  // file->Close();
  //cov.print_cvfile_info();

  // double x[101], y[2][101];
  // for (Int_t i=0;i!=101;i++){
  Bayes bayes;
  bayes.add_meas_component(0, 0, 0.184271 );
  bayes.add_meas_component(0 ,0, 0.75); 
  bayes.add_meas_component(0, 0, 0.0317704 );
  bayes.add_meas_component(0.315891, 0.0250575, 0.0759359 );
  bayes.add_meas_component(0, 0, 0.0317704 );
  bayes.add_meas_component(0, 0, 0.000132375 );
  bayes.add_meas_component(1.39223, 0.00129121, 0.000790154 );
  bayes.add_meas_component(0, 0, 0.000132375);
  
  // bayes.add_meas_component(3.64165, 0.746753, 0.184271); 
  // bayes.add_meas_component(10.5, 7.875, 0.75); 
  // bayes.add_meas_component(2.84751, 0.162341, 0.0317704);
  // bayes.add_meas_component(456.453, 46.4274, 0.0759359 );
  // bayes.add_meas_component(4.73775, 0.239406, 0.0317704);

  //  bayes.add_meas_component(3, 3 , 1);
  //bayes.add_meas_component(450,450, 1);
  
  
  bayes.do_convolution();
  // for (int i=0;i!=100;i++){
  //   double x = i/10.;
  //   double par[3]={1,1,1};
  //   //std::cout << bayes.Prop_Poisson_Pdf(&x,par) << std::endl;
  //   //std::cout << x << " " << f->Eval(x) << std::endl;
  // }
  
  std::cout << bayes.get_covariance_mc() << " " << bayes.get_covariance() << std::endl;

  // auto result = bayes.calculate_lower_upper(1);
  // x[i] = i;
  // y[0][i] = result.first;
  // y[1][i] = result.second;
  //    std::cout << i << " " << result.first << " " << result.second << std::endl;
  //  }

  // for (Int_t i=0;i!=101;i++){
  //   std::cout << x[i] << ", ";
  // }
  // std::cout << std::endl;

  // for (Int_t i=0;i!=101;i++){
  //   std::cout << y[0][i] << ", ";
  // }
  // std::cout << std::endl;

  // for (Int_t i=0;i!=101;i++){
  //   std::cout << y[1][i] << ", ";
  // }
  // std::cout << std::endl;
  
  TApplication theApp("theApp",&argc,argv);
  theApp.SetReturnFromRun(true);
  gStyle->SetOptStat(0);

  TCanvas c1("ToyMC","ToyMC",1200,800);
  c1.Draw();
  c1.Divide(3,2);
  for (int i=0;i!=1;i++){
    c1.cd(i+1);
    //    TF1 *f = bayes.get_f(i);
    //    TF1 *f = bayes.get_f_conv();
    TF1 *f = bayes.get_f_conv_num();
    f->Draw();
    //    TGraph *g1 = bayes.get_graph(i);
    //g1->Draw("AL");
    //bayes.get_f_num(i)->Draw();
  }
  c1.cd(2);
  bayes.get_f(0)->Draw();
  c1.cd(3);
  bayes.get_f(1)->Draw();
    
  theApp.Run();
  
  return 0;
}
