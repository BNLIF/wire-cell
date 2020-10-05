#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<map>
#include<vector>
#include<set>
#include<algorithm> // std::sort

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THStack.h"
#include "TF1.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"

#include "TRandom3.h"
#include "TGaxis.h"
#include "TStyle.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TROOT.h"

#include "TF1Convolution.h" // "root-config --features" --> should have "fftw"

///////
#include <chrono>
using namespace std::chrono;

//auto time0 = high_resolution_clock::now();
//auto time1 = high_resolution_clock::now();
//auto duration = duration_cast<microseconds>(stop - start);
//cout << "Time taken by function: "
//<< duration.count() << " microseconds" << endl;

//std::chrono::nanoseconds
//std::chrono::microseconds
//std::chrono::milliseconds
//std::chrono::seconds
//std::chrono::minutes
//std::chrono::hours

///////

void func_canv_margin(TCanvas *canv, double left, double right, double top, double bot)
{
  canv->SetLeftMargin(left);
  canv->SetRightMargin(right);
  canv->SetTopMargin(top);
  canv->SetBottomMargin(bot);
}

/////////////////////////////////////////////////////////////////////////////////// class

//
// Calcualte the credible interval of ratio=meas/pred
//
// Two solutions: analytic and ToyMC
// ROOT package needed for analytic solution: #include "TF1Convolution.h" // "root-config --features" --> should have "fftw"
//
// usage: root -l read_ratio_beta.cc
//
//
// exampleA->Set_analytic_result_precision( 1e-4 );
// exampleA->Set_analytic_gh_digitize( 5000 );// "5000" takes 1sec, "1000" takes 0.4sec
// exampleA->Calculate_ratio_lower_upper(1, true);// nSigma, digitize func_raito. If false --> no digitization
//
//

class TRatio
{
public:
  TRatio() {
    flag_clear = 0;
  }

  void Clear();
  
  static double pdf_func_sterling_poisson(double *x, double *par);
  static double com_func_Y2X(double *x, double *par);
  static double pdf_func_Y2X(double *x, double *par);
  
  void Set_data_range(double low, double hgh);
  void Set_ratio_range(double low, double hgh);
  void Set_analytic_result_precision(double val) { analytic_result_precision = val; }
  double Get_data_range_low() { return data_range_low; }
  double Get_data_range_hgh() { return data_range_hgh; }
  double Get_ratio_range_low() { return ratio_range_low; }
  double Get_ratio_range_hgh() { return ratio_range_hgh; }
  
  void Add_meas_component(double value, double weight);
  void Add_pred_component(double value, double weight); 

  double Get_sum_meas() { return sum_meas; }
  double Get_sum_pred() { return sum_pred; }
  
  TF1 *Get_meas_component(int i) { return map_pdf_func_meas_component[i]; }
  TF1 *Get_pred_component(int i) { return map_pdf_func_pred_component[i]; }
  
  int Get_size_meas() {return size_meas;}
  int Get_size_pred() {return size_pred;}

  TF1 *Get_min_component() {
    int minX_meas = 1;
    double minY_meas = 1e10;
    for(int idx=1; idx<=size_meas; idx++) {
      double Y_meas = map_meas_pars[idx][0]*map_meas_pars[idx][1];
      if( Y_meas < minY_meas ) { minX_meas=idx; minY_meas=Y_meas; }
    }
    
    int minX_pred = 1;
    double minY_pred = 1e10;
    for(int idx=1; idx<=size_pred; idx++) {
      double Y_pred = map_pred_pars[idx][0]*map_pred_pars[idx][1];
      if( Y_pred < minY_pred ) { minX_pred=idx; minY_pred=Y_pred; }
    }
    
    if( minY_meas<minY_pred ) {return map_pdf_func_meas_component[minX_meas];}
    else { return map_pdf_func_pred_component[minX_pred];}  
  }
  
  void Summation_meas_func();
  void Summation_pred_func();

  TF1 *Get_summation_meas_func() { return func_summation_meas; }
  TF1 *Get_summation_pred_func() { return func_summation_pred; }

  void Set_analytic_gh_digitize( int n_digitize ) { gh_digitize = n_digitize; }
  
  void Func_ratio_meas2pred();
  TF1 *Get_func_ratio_meas2pred() {
    if( !flag_ghh ) {return func_ratio_meas2pred;}
    else {return func_gh_func_ratio_meas2pred;}
  }
  
  void Calculate_ratio_lower_upper(int nSigma, bool flag_gh);
  double Get_ratio_lower() { return ratio_lower; }
  double Get_ratio_upper() { return ratio_upper; }
  double Get_ratio() { ratio_value = sum_meas/sum_pred; return ratio_value; }

  void Print_inputs() {
    cout<<endl<<" -------------- Input Information"<<endl<<endl;    
    cout<<" --->  data_range low/hgh: "<<data_range_low<<"\t"<<data_range_hgh<<endl;
    cout<<" ---> ratio_range low/hgh: "<<ratio_range_low<<"\t"<<ratio_range_hgh<<endl;

    cout<<endl;
    cout<<" ---> analytic CI precision: "<<analytic_result_precision<<endl;	
    cout<<" ---> gh_digitize number of pdf_func_ratio: "<<gh_digitize<<endl;
    
    cout<<endl;
    for(int i=1; i<=size_meas; i++) {
      cout<<TString::Format(" ---> %2d measurement (event and weight): %8.2f %8.6f",
                            i, map_meas_pars[i][0], map_meas_pars[i][1])<<endl;
    }
    cout<<endl;
    
    for(int i=1; i<=size_pred; i++) {
      cout<<TString::Format(" ---> %2d prediction  (event and weight): %8.2f %8.6f",
                            i, map_pred_pars[i][0], map_pred_pars[i][1])<<endl;
    }
    cout<<endl<<" -------------- "<<endl<<endl;
  }

  void Toy_ResultsOfRatio(int nSigma, int nToy);
  TH1D *Get_toy_ratio_PDF() { return h1_toy_ratio_PDF; }
  double Get_toy_ratio_lower() { return toy_ratio_lower; }
  double Get_toy_ratio_upper() { return toy_ratio_upper; }

  void Self_Check() {
    
    ////////////////////////////
    
    if( sum_meas==0 ) {  
      cout<<endl<<" -------------> Special case: sum_meas==0"<<endl;
      cout<<" -------------> Special case: sum_meas==0"<<endl;

      flag_meas0 = true;
    }
    
    ////////////////////////////
    
    if( sum_pred<=0 || sum_meas<0 ) {
      cerr<<endl<<" -------------------> ERROR: sum_pred<=0 || sum_meas<0"<<endl<<endl;
      exit(1);
    }
    
    ////////////////////////////

    if( ratio_range_low > Get_ratio() || ratio_range_hgh<Get_ratio() ) {
      cerr<<endl<<" -------------------> WARNING: ratio_range_low > ratio_nominal || ratio_range_hgh<ratioi_nominal"<<endl<<endl;
    }
      
  }
  
private:  
  int size_meas;
  int size_pred;

  map<int, TF1*>map_pdf_func_meas_component;
  map<int, TF1*>map_pdf_func_pred_component;

  double map_meas_pars[50][2];// maximum 50 components
  double map_pred_pars[50][2];

  double data_range_low;
  double data_range_hgh;
  double ratio_range_low;
  double ratio_range_hgh;
  
  double analytic_result_precision;
  
  double sum_meas;
  double sum_pred;

  TF1 *func_summation_meas;
  TF1 *func_summation_pred;
  
  map<int, TF1Convolution*>fconv_meas;
  map<int, TF1*>func_conv_meas;  
  map<int, TF1Convolution*>fconv_pred;
  map<int, TF1*>func_conv_pred;

  TF1 *func_ratio_meas2pred;
  int gh_digitize;
  TGraph *gh_func_ratio_meas2pred;
  TF1 *func_gh_func_ratio_meas2pred;
  bool flag_ghh;
  
  double ratio_value;
  double ratio_lower;
  double ratio_upper;

  double toy_ratio_lower;
  double toy_ratio_upper;

  TH1D *h1_toy_ratio_PDF;

  bool flag_meas0;

  int flag_clear;
  
  static TF1 *func_meas;
  static TF1 *func_pred;
};

TF1* TRatio::func_meas = NULL;
TF1* TRatio::func_pred = NULL;

//////

void TRatio::Toy_ResultsOfRatio(int nSigma, int nToy)
{
  vector<double>vec_ratio;
  
  h1_toy_ratio_PDF = new TH1D("h1_toy_ratio_PDF", "", 1000, ratio_range_low, ratio_range_hgh);
  
  cout<<" -------------- Begin of Toy"<<endl;
  
  for(int i=1; i<=nToy; i++) {    
    if( i%(nToy/10)==0 ) cout<<TString::Format(" ---> processing Toy %5.2f", i*1./nToy)<<endl;
    
    double user_meas = 0;
    double user_pred = 0;
    
    for(int j=1; j<=size_meas; j++) {
      user_meas += map_pdf_func_meas_component[j]->GetRandom();      
    }
    
    for(int j=1; j<=size_pred; j++) {
      user_pred += map_pdf_func_pred_component[j]->GetRandom();      
    }

    double ratio = user_meas/user_pred;
    vec_ratio.push_back(ratio);

    h1_toy_ratio_PDF->Fill( ratio );
  }

  double value_sigma = 1-TMath::Prob(pow(nSigma,2), 1);
    
  int size_vec = vec_ratio.size();  
  h1_toy_ratio_PDF->Scale( 1./size_vec );
  
  sort( vec_ratio.begin(), vec_ratio.end() );
    
  int line = 0;
  double xlow_a = 0;
  double xlow_b = 0;
  for(int i=0; i<size_vec; i++) {
    line++;
    if( line*1./size_vec < (1-value_sigma)/2 ) {
      xlow_a = vec_ratio.at( i );
      xlow_b = vec_ratio.at( i+1 );
    }
    else {
      break;
    }
  }

  line = 0;
  double xhgh_a = 0;
  double xhgh_b = 0;
  for(int i=size_vec-1; i>=0; i--) {
    line++;
    if( line*1./size_vec < (1-value_sigma)/2 ) {
      xhgh_a = vec_ratio.at( i );
      xhgh_b = vec_ratio.at( i-1 );
    }
    else {
      break;
    }
  }
  
  toy_ratio_lower = ( xlow_a+xlow_b )/2;
  toy_ratio_upper = ( xhgh_a+xhgh_b )/2;

  //////////////////////////////////////////////
  
  if( flag_meas0 ) {
    //cout<<" checking ... "<<endl;
    toy_ratio_lower = 0;
    toy_ratio_upper = 0;

    line = 0;
    double xhgh_a = 0;
    double xhgh_b = 0;
    for(int i=size_vec-1; i>=0; i--) {
      line++;
      if( line*1./size_vec < (1-value_sigma) ) {
	xhgh_a = vec_ratio.at( i );
	xhgh_b = vec_ratio.at( i-1 );
      }
      else {
	break;
      }
    }
    
    toy_ratio_upper = ( xhgh_a+xhgh_b )/2;    
  }// flag_meas0

  
  cout<<" -------------- End of Toy"<<endl;
  
}
  
//////

void TRatio::Calculate_ratio_lower_upper(int nSigma, bool flag_gh)
{
  flag_ghh = flag_gh;
  
  double value_sigma = 1 - TMath::Prob( pow(nSigma, 2), 1);

  int line = 0;
  double val_typeA = 0;
  double val_typeB = 0;
  
  /////////////////////////////////////////////////////////////////////////////////////

  if(flag_gh) {
    line = 0;
    val_typeA = ratio_range_low;
    val_typeB = ratio_range_hgh;
    while( fabs(val_typeA-val_typeB)>analytic_result_precision ) {
      line++;
      double val_mid = (val_typeA+val_typeB)/2;
      double y_mid = func_gh_func_ratio_meas2pred->Integral(val_mid, ratio_range_hgh);
      if( y_mid > (1-value_sigma)/2 ) {
	val_typeA = val_mid;
      }
      else {
	val_typeB = val_mid;
      }
    
      if( line>30 ) {
	cerr<<" Error: cannot find ratio_sigma_upper"<<endl;
	break;
      }
    }
  
    ratio_upper = (val_typeA+val_typeB)/2;
  
    ////////////////////////
  
    line = 0;
    val_typeA = ratio_range_low;
    val_typeB = ratio_range_hgh;
    while( fabs(val_typeA-val_typeB)>analytic_result_precision ) {
      line++;
      double val_mid = (val_typeA+val_typeB)/2;
      double y_mid = func_gh_func_ratio_meas2pred->Integral(ratio_range_low, val_mid);
      if( y_mid > (1-value_sigma)/2 ) {
	val_typeB = val_mid;
      }
      else {
	val_typeA = val_mid;
      }
    
      if( line>30 ) {
	cerr<<" Error: cannot find ratio_sigma_lower"<<endl;
	break;
      }
    }
  
    ratio_lower = (val_typeA+val_typeB)/2;

    ///////////////////////////////////////////////

    if( flag_meas0 ) {
      //cout<<" ---> checking ..."<<endl;
      ratio_lower = 0;
      ratio_upper = 0;
    
      line = 0;
      val_typeA = ratio_range_low;
      val_typeB = ratio_range_hgh;
      while( fabs(val_typeA-val_typeB)>analytic_result_precision ) {
	line++;
	double val_mid = (val_typeA+val_typeB)/2;
	double y_mid = func_gh_func_ratio_meas2pred->Integral(val_mid, ratio_range_hgh);
	if( y_mid > (1-value_sigma) ) {
	  val_typeA = val_mid;
	}
	else {
	  val_typeB = val_mid;
	}
    
	if( line>30 ) {
	  cerr<<" Error: cannot find ratio_sigma_upper"<<endl;
	  break;
	}
      }
  
      ratio_upper = (val_typeA+val_typeB)/2;       
    
    }// flag_meas0    
  }
  
  /////////////////////////////////////////////////////////////////////////////////////

  if(!flag_gh) {
    line = 0;
    val_typeA = ratio_range_low;
    val_typeB = ratio_range_hgh;
    while( fabs(val_typeA-val_typeB)>analytic_result_precision ) {
      line++;
      double val_mid = (val_typeA+val_typeB)/2;
      double y_mid = func_ratio_meas2pred->Integral(val_mid, ratio_range_hgh);
      if( y_mid > (1-value_sigma)/2 ) {
	val_typeA = val_mid;
      }
      else {
	val_typeB = val_mid;
      }
    
      if( line>30 ) {
	cerr<<" Error: cannot find ratio_sigma_upper"<<endl;
	break;
      }
    }
  
    ratio_upper = (val_typeA+val_typeB)/2;
  
    ////////////////////////
  
    line = 0;
    val_typeA = ratio_range_low;
    val_typeB = ratio_range_hgh;
    while( fabs(val_typeA-val_typeB)>analytic_result_precision ) {
      line++;
      double val_mid = (val_typeA+val_typeB)/2;
      double y_mid = func_ratio_meas2pred->Integral(ratio_range_low, val_mid);
      if( y_mid > (1-value_sigma)/2 ) {
	val_typeB = val_mid;
      }
      else {
	val_typeA = val_mid;
      }
    
      if( line>30 ) {
	cerr<<" Error: cannot find ratio_sigma_lower"<<endl;
	break;
      }
    }
  
    ratio_lower = (val_typeA+val_typeB)/2;

    ///////////////////////////////////////////////

    if( flag_meas0 ) {
      //cout<<" ---> checking ..."<<endl;
      ratio_lower = 0;
      ratio_upper = 0;
    
      line = 0;
      val_typeA = ratio_range_low;
      val_typeB = ratio_range_hgh;
      while( fabs(val_typeA-val_typeB)>analytic_result_precision ) {
	line++;
	double val_mid = (val_typeA+val_typeB)/2;
	double y_mid = func_ratio_meas2pred->Integral(val_mid, ratio_range_hgh);
	if( y_mid > (1-value_sigma) ) {
	  val_typeA = val_mid;
	}
	else {
	  val_typeB = val_mid;
	}
    
	if( line>30 ) {
	  cerr<<" Error: cannot find ratio_sigma_upper"<<endl;
	  break;
	}
      }
  
      ratio_upper = (val_typeA+val_typeB)/2;       
    
    }// flag_meas0
  }
  
  ///////////////////////////////////////////////
  //
  // other special case: Prob(lower_range, ratio_M/P) < p-value required, then ...
  //
}

void TRatio::Func_ratio_meas2pred()
{
  TRatio::func_meas = func_summation_meas;
  TRatio::func_pred = func_summation_pred;
  
  TF1 *func_ratio_Y2X = new TF1("pdf_func_Y2X", TRatio::pdf_func_Y2X, ratio_range_low, ratio_range_hgh, 2);
  func_ratio_Y2X->SetParameters(data_range_low, data_range_hgh);
  func_ratio_Y2X->SetNpx(60000);
  func_ratio_meas2pred = func_ratio_Y2X;

  ////////////////////////////////////////////////////////////////////
  
  int N_digitize = gh_digitize;
  double step = (ratio_range_hgh - ratio_range_low)/N_digitize;
  
  gh_func_ratio_meas2pred = new TGraph();
  gh_func_ratio_meas2pred->SetName("gh_func_ratio_meas2pred");
    
  for(int idx=0; idx<=N_digitize; idx++) {
    double xval = ratio_range_low + idx*step;
    double yval = func_ratio_meas2pred->Eval( xval );
    gh_func_ratio_meas2pred->SetPoint( gh_func_ratio_meas2pred->GetN(), xval, yval );
  }

  func_gh_func_ratio_meas2pred = new TF1("func_gh_func_ratio_meas2pred",
					 [&](double *x, double *){ return gh_func_ratio_meas2pred->Eval(x[0]); },
					 ratio_range_low, ratio_range_hgh, 0);
  
  func_gh_func_ratio_meas2pred->SetNpx(N_digitize);
  
}

//////

void TRatio::Summation_meas_func()
{
  TString roostr = "";  

  if( size_meas==0 ) {
    cerr<<" ERROR: no meas"<<endl;
  }
  
  if( size_meas==1 ) {
    func_summation_meas = map_pdf_func_meas_component[1];
  }
  else {
    int line = 2;
    fconv_meas[line] = new TF1Convolution( map_pdf_func_meas_component[line-1], map_pdf_func_meas_component[line],
                                      data_range_low, data_range_hgh, true );
    fconv_meas[line]->SetNofPointsFFT(10000);
    roostr = TString::Format("func_conv_meas_%02d", line);
    func_conv_meas[line] = new TF1(roostr, fconv_meas[line], data_range_low, data_range_hgh, 0);
    func_conv_meas[line]->SetNpx(60000);
  
    if( size_meas>=3 ) {
      for( int idx=3; idx<=size_meas; idx++ ) {
        line++;
        fconv_meas[line] = new TF1Convolution( func_conv_meas[line-1], map_pdf_func_meas_component[line],
                                               data_range_low, data_range_hgh, true );
        fconv_meas[line]->SetNofPointsFFT(10000);
        roostr = TString::Format("func_conv_meas_%02d", line);
        func_conv_meas[line] = new TF1(roostr, fconv_meas[line], data_range_low, data_range_hgh, 0);
        func_conv_meas[line]->SetNpx(60000);
      }      
    }
    
    func_summation_meas = func_conv_meas[line];  
  }  
}

void TRatio::Summation_pred_func()
{
  TString roostr = "";  

  if( size_pred==0 ) {
    cerr<<" ERROR: no pred"<<endl;
  }
  
  if( size_pred==1 ) {
    func_summation_pred = map_pdf_func_pred_component[1];
  }
  else {
    int line = 2;
    fconv_pred[line] = new TF1Convolution( map_pdf_func_pred_component[line-1], map_pdf_func_pred_component[line],
                                      data_range_low, data_range_hgh, true );
    fconv_pred[line]->SetNofPointsFFT(10000);
    roostr = TString::Format("func_conv_pred_%02d", line);
    func_conv_pred[line] = new TF1(roostr, fconv_pred[line], data_range_low, data_range_hgh, 0);
    func_conv_pred[line]->SetNpx(60000);
  
    if( size_pred>=3 ) {
      for( int idx=3; idx<=size_pred; idx++ ) {
        line++;
        fconv_pred[line] = new TF1Convolution( func_conv_pred[line-1], map_pdf_func_pred_component[line],
                                               data_range_low, data_range_hgh, true );
        fconv_pred[line]->SetNofPointsFFT(10000);
        roostr = TString::Format("func_conv_pred_%02d", line);
        func_conv_pred[line] = new TF1(roostr, fconv_pred[line], data_range_low, data_range_hgh, 0);
        func_conv_pred[line]->SetNpx(60000);
      }      
    }
    
    func_summation_pred = func_conv_pred[line];  
  }  
}

//////
  
double TRatio::com_func_Y2X(double *x, double *par)
{
  double result = 0;
  
  double xcur = x[0];
  double Y2X = par[0];
  double range_low = par[1];
  double range_hgh = par[2];
  
  if( xcur<range_low || xcur>range_hgh ) {
    int flag = 0;
  }
  else {
    double eval_fX = func_pred->Eval( xcur );
    double eval_fY = func_meas->Eval( Y2X * xcur );
    double abs_X = fabs( xcur );
    result = abs_X * eval_fX * eval_fY;
  }
  
  return result;
}

double TRatio::pdf_func_Y2X(double *x, double *par)
{
  double result = 0;

  double Y2X = x[0];
  double range_low = par[0];
  double range_hgh = par[1];

  TF1 *roofunc_Y2X = new TF1("roofunc_Y2X", com_func_Y2X, range_low, range_hgh, 3);
  roofunc_Y2X->SetParameter(0, Y2X);
  roofunc_Y2X->SetParameter(1, range_low);
  roofunc_Y2X->SetParameter(2, range_hgh);
  roofunc_Y2X->SetNpx(10000);
  
  /// Exact solution (cannot reach tolerance because of roundoff error):
  /// gErrorIgnoreLevel = 6001;    
  double val_integration = roofunc_Y2X->Integral( range_low, range_hgh );    
  result = val_integration;

  delete roofunc_Y2X;    
  return result;
}


//////

void TRatio::Add_meas_component(double value, double weight)
{
  size_meas++;
  map_meas_pars[size_meas][0] = value;
  map_meas_pars[size_meas][1] = weight;

  TString roostr = TString::Format("map_pdf_func_meas_component_%02d", size_meas);
  map_pdf_func_meas_component[size_meas] = new TF1(roostr, pdf_func_sterling_poisson, data_range_low, data_range_hgh, 4);
  map_pdf_func_meas_component[size_meas]->SetParameters( value, weight,  data_range_low, data_range_hgh);
  map_pdf_func_meas_component[size_meas]->SetNpx(60000);

  sum_meas += value*weight;
}

void TRatio::Add_pred_component(double value, double weight)
{
  size_pred++;
  map_pred_pars[size_pred][0] = value;
  map_pred_pars[size_pred][1] = weight;

  TString roostr = TString::Format("map_pdf_func_pred_component_%02d", size_pred);
  map_pdf_func_pred_component[size_pred] = new TF1(roostr, pdf_func_sterling_poisson, data_range_low, data_range_hgh, 4);
  map_pdf_func_pred_component[size_pred]->SetParameters( value, weight,  data_range_low, data_range_hgh);
  map_pdf_func_pred_component[size_pred]->SetNpx(60000);

  sum_pred += value*weight;
}

//////

void TRatio::Set_data_range(double low, double hgh)
{
  data_range_low = low;
  data_range_hgh = hgh;
}

void TRatio::Set_ratio_range(double low, double hgh)
{
  ratio_range_low = low;
  ratio_range_hgh = hgh;
}

//////

double TRatio::pdf_func_sterling_poisson(double *x, double *par)
{
  double val_meas = par[0];
  double val_weight   = par[1];
  double xrange_low = par[2];
  double xrange_hgh = par[3];

  double xcur = x[0];
  double val_pred = xcur/val_weight;

  double result = 0;

  {
    if( xcur<xrange_low || xcur>xrange_hgh ) {
      int flag = 0;
    }
    else {
      if( val_meas>0 ) {

        /// 
        /// https://dlmf.nist.gov/5.11#i
        /// https://en.wikipedia.org/wiki/Stirling%27s_approximation    
        const int cn = 6;
        double ci[cn] = {1., 1./12, 1./288, -139./51840, -571./2488320, 163879./209018880};
        for(int i=0; i<cn; i++) {
          result += ci[i]/pow( val_meas, i );
        }
        result = exp( val_meas - val_pred + val_meas*log(val_pred/val_meas) ) / ( sqrt(2*TMath::Pi()*val_meas) * result  ) /val_weight;
        
        /// Jarrett's formula, DocDB 29109      
        // const int cn = 6;
        // double ci[cn] = {1., -1./12, 1./288, 139./51840, -571./2488320, -163879./209018880};
        // for(int i=0; i<cn; i++) {
        //   result += ci[i]/pow( val_meas, i+0.5 );
        // }
        // result = exp( val_meas - val_pred + val_meas*log(val_pred/val_meas) ) / ( sqrt(2*TMath::Pi()) ) * result /val_weight;
        
      }
      else if(val_meas==0) {// back to  Poisson distribution
        result = exp(-val_pred) / val_weight;
      }
      else {
        int flag = 0;
      }
    }
  }
  
  return result;
}

void TRatio::Clear()
{
  flag_clear++;
  
  size_meas = 0;
  size_pred = 0;

  map_pdf_func_meas_component.clear();
  map_pdf_func_pred_component.clear();

  for(int i=0; i<50; i++) {
    for(int j=0; j<2; j++) {
      map_meas_pars[i][j] = 0;
      map_pred_pars[i][j] = 0;
    }
  }
  
  data_range_low = 0;
  data_range_hgh = 0;
  ratio_range_low = 0;
  ratio_range_hgh = 0;
  analytic_result_precision = 1e-4;
  
  sum_meas = 0;
  sum_pred = 0;
  
  // func_summation_meas = NULL;
  // func_summation_pred = NULL;
  
  // fconv_meas.clear();
  // func_conv_meas.clear();  
  // fconv_pred.clear();
  // func_conv_pred.clear();

  // func_ratio_meas2pred = NULL;
  gh_digitize = 5000;
  // gh_func_ratio_meas2pred = NULL;
  // func_gh_func_ratio_meas2pred = NULL;
  flag_ghh = false;
  
  ratio_value = 0;
  ratio_lower = 0;
  ratio_upper = 0;
  
  toy_ratio_lower = 0;
  toy_ratio_upper = 0;

  // h1_toy_ratio_PDF = NULL;

  flag_meas0 = false;

  if(flag_clear==1) {
    func_summation_meas = NULL;
    func_summation_pred = NULL;
    
    func_ratio_meas2pred = NULL;
    gh_func_ratio_meas2pred = NULL;
    func_gh_func_ratio_meas2pred = NULL;
    
    h1_toy_ratio_PDF = NULL;
  }
  else {
    delete func_summation_meas;
    delete func_summation_pred;

    delete func_ratio_meas2pred;
    delete gh_func_ratio_meas2pred;
    delete func_gh_func_ratio_meas2pred;
    
    delete h1_toy_ratio_PDF;
  }
    
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// MAIN ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

void read_ratio_v02()
{
  gErrorIgnoreLevel = 6001;
  
  TString roostr = "";

  auto time1 = high_resolution_clock::now();

  /// declaration and definition
  TRatio *exampleA = new TRatio();

  /// clear and initialization
  exampleA->Clear();
  exampleA->Set_data_range(0, 1000);
  exampleA->Set_ratio_range(0, 2);

  /// set precision of CI by analytic method
  exampleA->Set_analytic_result_precision( 1e-4 );

  /// set the number of digitization on pdf of ratio by analytic mtethod
  exampleA->Set_analytic_gh_digitize( 5000 );
    
  /// measurement: event and weight
  exampleA->Add_meas_component(460, 1);
   
  /// prediction: event and weight
  exampleA->Add_pred_component(400, 0.25);
  exampleA->Add_pred_component(200, 1);
  exampleA->Add_pred_component(300, 0.5);
 
  /// necessary
  exampleA->Self_Check();
  
  /// print out inputs
  exampleA->Print_inputs();

  ////////////////////// Analytic method

  /// calculate the PDF of summation of meas, or pred
  exampleA->Summation_meas_func();  
  exampleA->Summation_pred_func();

  /// calculate the PDF of the ratio=meas/pred
  exampleA->Func_ratio_meas2pred();
  
  /// calculate the credible interval of the ratio
  exampleA->Calculate_ratio_lower_upper(1, true);// nSigma, digitize func_raito

  /// check time used
  auto time7 = high_resolution_clock::now();
  auto diff_1_7 = duration_cast<milliseconds>(time7 - time1);
  int cc_diff_1_7 = diff_1_7.count();
  cout<<endl<<TString::Format(" ---> Analytic method takes time (ms): %d", cc_diff_1_7)<<endl<<endl;
  
  ////////////////////// ToyMC method

  /// calculate the credible interval of the ratio by toy
  /// (nSigma, number of toys)
  exampleA->Toy_ResultsOfRatio(1, 1000000);
 
  /////////////////////////// self_check and print results
  ///////////////////////////

  cout<<endl<<TString::Format(" ---> Integration of analytic Ratio PDF (should be close to 1, else set range): %8.5f",
                              exampleA->Get_func_ratio_meas2pred()->Integral( exampleA->Get_ratio_range_low(), exampleA->Get_ratio_range_hgh() )
                              )<<endl<<endl;
 
  cout<<TString::Format(" ---> Ratio %8.4f, analytic lower/upper: %8.4f, %8.4f",
			      exampleA->Get_ratio(),
			      exampleA->Get_ratio_lower(),
			      exampleA->Get_ratio_upper()
			      )<<endl<<endl;
 
  cout<<TString::Format(" ---> Ratio %8.4f,    toyMC lower/upper: %8.4f, %8.4f",
			      exampleA->Get_ratio(),
			      exampleA->Get_toy_ratio_lower(),
			      exampleA->Get_toy_ratio_upper()
			      )<<endl<<endl;

  /////////////////////////// It will save time if don't draw figures
  ///////////////////////////
  
  roostr = "canv_data";
  TCanvas *canv_data = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_data, 0.15, 0.2,0.1,0.15);

  exampleA->Get_min_component()->Draw();

  for(int idx=1; idx<=exampleA->Get_size_pred(); idx++) {
    exampleA->Get_pred_component(idx)->Draw("same");
    exampleA->Get_pred_component(idx)->SetLineColor(kRed);
  }
  
  for(int idx=1; idx<=exampleA->Get_size_meas(); idx++) {
    exampleA->Get_meas_component(idx)->Draw("same");
    exampleA->Get_meas_component(idx)->SetLineColor(kBlue);
  }
  
  exampleA->Get_summation_pred_func()->Draw("same");    
  exampleA->Get_summation_pred_func()->SetLineColor(kGray+1);

  exampleA->Get_summation_meas_func()->Draw("same");    
  exampleA->Get_summation_meas_func()->SetLineColor(kBlack);
      
  ///////////////////////////
  
  roostr = "canv_analytical";
  TCanvas *canv_analytical = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_analytical, 0.15, 0.2,0.1,0.15);

  exampleA->Get_func_ratio_meas2pred()->Draw();
  
  ///////////////////////////
  
  roostr = "canv_toy";
  TCanvas *canv_toy = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_toy, 0.15, 0.2,0.1,0.15);

  exampleA->Get_toy_ratio_PDF()->Draw("hist");
  

}
