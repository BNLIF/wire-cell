#ifndef LEEANA_BAYES
#define LEEANA_BAYES

#include "TF1.h"
#include "TGraph.h"
#include "TF1Convolution.h"

namespace LEEana{
  double Prop_Poisson_Pdf(double *x, double *par);
  
  class Bayes{
  public:
    Bayes();
    ~Bayes();

    void add_meas_component(double meas, double sigma2, double weight, int flag = 1);
    void do_convolution();

    TF1* get_f_conv(){return f_conv;};
    TF1* get_f_conv_num(){return f_conv_num;};
    TF1 *get_f(int i){return meas_pdf_vec.at(i);};
    TGraph* get_graph(int i){return map_index_graph[i];};
    TF1 *get_f_num(int i){return f1_num_vec.at(i);};
    
    double get_mean(){return mean;};
    double get_covariance_mc();
    double get_covariance();
    

    std::pair<double, double> calculate_lower_upper(double nsigma);
    
  private:
    int num_component;
    
    std::vector<double> meas_vec;
    std::vector<double> sig2_vec;
    std::vector<double> weight_vec;

    std::vector<TF1*> meas_pdf_vec;
    std::map<int,TGraph*> map_index_graph;
    std::vector<TF1*> f1_num_vec;
    std::vector<TF1*> f_conv_vec;


    std::vector<TF1*> f_test_vec;
    
    std::vector<TF1Convolution*> conv_vec;
    
    TF1 *f_conv;
    
    TF1 *f_conv_num;
    TGraph *g1;

    double llimit; 
    double hlimit; 
    
    double mean;
    double acc_sigma2;
    double sigma;
    
  };
}

#endif
