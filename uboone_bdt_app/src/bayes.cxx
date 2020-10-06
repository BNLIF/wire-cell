

#include "WCPLEEANA/bayes.h"
#include <cmath>

using namespace LEEana;

LEEana::Bayes::Bayes()
  : num_component(0)
  , mean(0)
  , acc_sigma2(0)
  , sigma(1)
  , f_conv(0)
  , f_conv_num(0)
  , g1(0)
{
}


LEEana::Bayes::~Bayes(){
  for (auto it = meas_pdf_vec.begin(); it != meas_pdf_vec.end(); it++){
    delete (*it);
  }
  for (auto it = f1_num_vec.begin(); it != f1_num_vec.end(); it++){
    delete (*it);
  }
  for (auto it = conv_vec.begin(); it != conv_vec.end(); it++){
    delete (*it);
  }
  // for (auto it = f_conv_vec.begin(); it!= f_conv_vec.end(); it++){
  //   delete (*it);
  // }
  for (auto it = f_test_vec.begin(); it!= f_test_vec.end(); it++){
    delete (*it);
  }
 
  if (f_conv_num != (TF1*)0) delete f_conv_num;
  if (f_conv != (TF1*)0) delete f_conv;
  if (g1 != (TGraph*)0) delete g1;
}

void LEEana::Bayes::add_meas_component(double meas, double sigma2, double weight, int flag){
  meas_vec.push_back(meas);
  sig2_vec.push_back(sigma2);

  if(meas == 0)
    weight_vec.push_back(weight);
  else
    weight_vec.push_back(1);

  mean += meas;
  acc_sigma2 += sigma2;

  double tmp_llimit = 0;
  double tmp_hlimit = 0;
  
  if (meas !=0){
    tmp_llimit = meas - 12*sqrt(sigma2) - 1;
    tmp_hlimit = meas + 12*sqrt(sigma2) + 1;
  }else{
    tmp_llimit = meas - 10;
    tmp_hlimit = meas + 10;
  }
  if (tmp_llimit <0) tmp_llimit = 0;


  // std::cout << tmp_llimit << " " << tmp_hlimit << " " << meas_vec.back() << " " << sig2_vec.back() << " " << weight_vec.back() << std::endl;
  
  TF1 *f1 = new TF1(Form("f_test_%d",num_component), Prop_Poisson_Pdf, tmp_llimit, tmp_hlimit, 5);
  f1->SetParameters(meas_vec.back(), sig2_vec.back(), weight_vec.back(), tmp_llimit, tmp_hlimit);
  f1->SetNpx(60000);
  f_test_vec.push_back(f1);
  
  // g1 = new TGraph();
  // g1->SetName(Form("g1_%d",num_component));
  // g1->SetPoint(g1->GetN(),tmp_llimit-0.2,0);
  // g1->SetPoint(g1->GetN(),tmp_llimit-0.1,0);
  // for (int idx = 0; idx <=5000;idx++){
  //   double xval = tmp_llimit + (tmp_hlimit - tmp_llimit)/5000.*(idx+0.5);
  //   double yval = f1->Eval(xval);
  //   g1->SetPoint(g1->GetN(),xval,yval);
  // }
  // g1->SetPoint(g1->GetN(),tmp_hlimit+0.1,0);
  // g1->SetPoint(g1->GetN(),tmp_hlimit+0.2,0);
   

  
  // TF1* f1_num = new TF1(Form("f1_num_%d",num_component),
  // 			    [&](double *x, double *){
  // 			      return g1->Eval(x[0]);
  // 			    },
  // 			    tmp_llimit, tmp_hlimit, 0);
  // f1_num->SetNpx(5000);

  // //  std::cout << f_conv_num->Eval(0.5) << " " << f1->Eval(0.5) << std::endl;
  
  // f1_num_vec.push_back(f1_num);
  
  //  num_component ++;
}


void LEEana::Bayes::do_convolution(){
  llimit = 0;
  if (mean != 0){
    hlimit = mean + 12*sqrt(acc_sigma2) + 1;
  }else{
    hlimit = mean + 10;
  }
  hlimit = std::max(hlimit, mean * 5 + 20);
  
  // std::cout << llimit << " " << hlimit << std::endl;
    
  for (size_t i=0; i!= meas_vec.size(); i++){
    double meas = meas_vec.at(i);
    double sigma2 = sig2_vec.at(i);
    double weight = weight_vec.at(i);

    // double tmp_llimit = 0;
    // double tmp_hlimit = 0;
  
    // if (meas !=0){
    //   tmp_llimit = meas - 12*sqrt(sigma2) - 1;
    //   tmp_hlimit = meas + 12*sqrt(sigma2) + 1;
    // }else{
    //   tmp_llimit = meas - 10;
    //   tmp_hlimit = meas + 10;
    // }
    // if (tmp_llimit <0) tmp_llimit = 0;

    TF1 *f1 = new TF1(Form("f_%d",num_component), Prop_Poisson_Pdf, llimit, hlimit, 5);
    f1->SetParameters(meas, sigma2, weight, llimit, hlimit);
    f1->SetNpx(60000);
    meas_pdf_vec.push_back(f1);
    
    if (i==0){
      f_conv = new TF1(Form("f_conv_%d",num_component), Prop_Poisson_Pdf, llimit, hlimit, 5);
      f_conv->SetParameters(meas, sigma2, weight, llimit, hlimit);
      f_conv->SetNpx(60000);
      f_conv_vec.push_back(f_conv);
    }else{
      TF1Convolution *conv = new TF1Convolution(f_conv, f1, llimit, hlimit, true);
      conv->SetNofPointsFFT(10000);
      conv_vec.push_back(conv);
      
      TF1 *new_f_conv = new TF1(Form("f_conv_%d",num_component), conv, llimit, hlimit, 0);
      new_f_conv->SetNpx(60000);
      f_conv_vec.push_back(new_f_conv);

      f_conv = new_f_conv;
    }
    
    num_component ++;
  }
    
  
    
  g1 = new TGraph();
  g1->SetPoint(g1->GetN(),llimit-0.2,0);
  g1->SetPoint(g1->GetN(),llimit-0.1,0);
  for (int idx = 0; idx <=5000;idx++){
    double xval = llimit + (hlimit - llimit)/5000.*(idx+0.5);
    double yval = f_conv->Eval(xval);
    g1->SetPoint(g1->GetN(),xval,yval);
  }
  g1->SetPoint(g1->GetN(),hlimit+0.1,0);
  g1->SetPoint(g1->GetN(),hlimit+0.2,0);
    


  f_conv_num = new TF1(Form("f_conv_num_%d",num_component),
		       [&](double *x, double *){
			 return g1->Eval(x[0]);
		       },
		       llimit, hlimit, 0);
  f_conv_num->SetNpx(5000);


}


std::pair<double, double> LEEana::Bayes::calculate_lower_upper(double nSigma){
  
  
  double value_sigma = 1 - TMath::Prob( pow(nSigma, 2), 1);
  gErrorIgnoreLevel = 6001;
 
  int line = 0;

  double ratio_range_low, ratio_range_high;
  f_conv_num->GetRange(ratio_range_low, ratio_range_high);


  double low_bound = ratio_range_low;
  double hig_bound = ratio_range_high;
  
  double val_typeA = 0;
  double val_typeB = 0;

  double low_integral = f_conv_num->Integral(ratio_range_low, mean);
  double full_integral = f_conv_num->Integral(ratio_range_low, ratio_range_high);

  // low limit ...
  double low_percentage = low_integral/full_integral;
  double analytic_result_precision = 1e-4;

  double level;

  //  std::cout << ratio_range_low << " " << mean << " " << ratio_range_high << " " << low_percentage << " " << value_sigma/2. << std::endl;
  
  if (low_percentage > value_sigma/2.){
  //  std::cout << ratio_range_low << " " << ratio_range_high << " " << low_integral/full_integral << std::endl;
    level = low_percentage - value_sigma/2.;

    line = 0;
    val_typeA = ratio_range_low;
    val_typeB = ratio_range_high;
    while( fabs(val_typeA-val_typeB)>analytic_result_precision ) {
      line++;
      double val_mid = (val_typeA+val_typeB)/2;
      double y_mid = f_conv_num->Integral(ratio_range_low,val_mid)/full_integral;
      if( y_mid > level ) {
	val_typeB = val_mid;
      } else {
	val_typeA = val_mid;
      }
    
      if( line>30 ) {
	std::cerr<<" Error: cannot find ratio_sigma_upper"<<std::endl;
	break;
      }
    }
    low_bound = (val_typeA+val_typeB)/2;

    level = 1 - level - value_sigma;
    
    line = 0;
    val_typeA = ratio_range_low;
    val_typeB = ratio_range_high;
    while( fabs(val_typeA-val_typeB)>analytic_result_precision ) {
      line++;
      double val_mid = (val_typeA+val_typeB)/2;
      double y_mid = f_conv_num->Integral(val_mid, ratio_range_high)/full_integral;
      if( y_mid > level ) {
	val_typeA = val_mid;
      } else {
	val_typeB = val_mid;
      }
    
      if( line>30 ) {
	std::cerr<<" Error: cannot find ratio_sigma_upper"<<std::endl;
	break;
      }
    }
  
    hig_bound = (val_typeA+val_typeB)/2;
    
  }else{

    level = 1-value_sigma;
    
    line = 0;
    val_typeA = ratio_range_low;
    val_typeB = ratio_range_high;
    while( fabs(val_typeA-val_typeB)>analytic_result_precision ) {
      line++;
      double val_mid = (val_typeA+val_typeB)/2;
      double y_mid = f_conv_num->Integral(val_mid, ratio_range_high)/full_integral;
      if( y_mid > level ) {
	val_typeA = val_mid;
      } else {
	val_typeB = val_mid;
      }

      
      
      if( line>30 ) {
	std::cerr<<" Error: cannot find ratio_sigma_upper"<<std::endl;
	break;
      }
    }
  
    hig_bound = (val_typeA+val_typeB)/2;
    
  }

  //  std::cout << low_bound <<  " " << hig_bound << std::endl;

  return std::make_pair(low_bound, hig_bound);
  
}

double LEEana::Bayes::get_covariance(){
  double sum = 0;
  double sum1 = 0;

  int npoints = f_conv_num->GetNpx();

  double max_val = f_conv_num->Eval(mean);
  
  for (int i=0;i!=10*npoints;i++){
    double x = llimit + (hlimit - llimit)/10./npoints *(i+0.5);
    double val = f_conv_num->Eval(x);
    if (num_component >1){
      if (x < mean - 6 * sqrt(acc_sigma2) || val < max_val*1e-3) continue;
    }
    double corr = 0;
    if (mean != 0)
      corr= 1./pow(x/mean, num_component-1);
    else
      corr= 1./pow(x, num_component-1);
    
    if (std::isinf(corr) || std::isnan(corr)) corr = 0;
    // corr = 1;
    
    sum  += val * pow(x-mean,2) *corr;
    sum1 += val * corr;
  }
  
  //  std::cout << std::endl;
  //std::cout << mean << " " << sum/sum1 << " " << acc_sigma2 << std::endl;
  
  return sum/sum1;
}

double LEEana::Bayes::get_covariance_mc(){

  double sum = 0;
  double sum1 = 0;

  //  for (Int_t i=0;i!=meas_pdf_vec.size();i++){
  //  std::cout << i << " " << meas_pdf_vec.at(i)->GetRandom() << std::endl;
  // }
  
  for (Int_t j=0;j!=100000;j++){
    double x = 0;
    for (int i=0;i!=f_test_vec.size();i++){
      x += f_test_vec.at(i)->GetRandom();
    }
    double corr;
    if (mean !=0)
      corr = 1./pow(x/mean,num_component-1);
    else
      corr = 1./pow(x,num_component-1);
    if (std::isnan(corr) || std::isinf(corr)) corr = 0;    
    
    sum += pow(x-mean,2) * corr;
    sum1 += corr;
  }
  
  double result = sum / sum1;

  //  std::cout << mean << " " << num_component << " " << result << " " << acc_sigma2 << std::endl;
  
  return result; 

}



double LEEana::Prop_Poisson_Pdf(double *x, double *par){
  double mu = x[0]; // expectation ...
  
  double meas = par[0]; // actual measurement
  double sigma2 = par[1]; // sigma2 of the measurement ...

  double weight = par[2]; // e.g. POT scaling ...

  double llimit = par[3];
  double hlimit = par[4];
  
  double eff_weight = 1;  // effective weight ...
  double eff_meas = meas; // effective measurement ...

  if (meas <0 ) return 0;
  
  if (meas > 0){
    eff_weight = sigma2/meas;
    eff_meas = meas/eff_weight;
  }
  

  //std::cout << mu << " " << meas << " " << sigma2 << " " << weight << std::endl;
  
  double result = 0;

  if (mu > llimit && mu < hlimit){
    if(mu >=0){
      if (eff_meas >0)
	result = std::exp(eff_meas - mu/eff_weight/weight + eff_meas * std::log(mu/eff_weight/eff_meas/weight));
      else
	result = std::exp(-mu/eff_weight/weight);
    }
    
    if (std::isnan(result) || result < 0 || std::isinf(result) ) result = 0;
  }
    
  return result;
}
