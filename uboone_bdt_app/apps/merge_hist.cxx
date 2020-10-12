#include <iostream>

#include "WCPLEEANA/master_cov_matrix.h"
#include "WCPLEEANA/bayes.h"

#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h" 
#include "TAxis.h"
#include "TLegend.h"
#include "TMath.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"

using namespace std;
using namespace LEEana;

int main( int argc, char** argv )
{
  
  if (argc < 2){
    std::cout << "./merge_hist -r[#run, 0 for all] -e[1 for standard, 2 for Bayesian] -L[LEE strength]" << std::endl;
  }
  
  int run = 1; // run 1 ...
  int flag_err = 1;// 1 for standard, 2 for Bayesian ...
  float lee_strength = 0; // no LEE strength ...
  int flag_display = 0;
  
  for (Int_t i=1;i!=argc;i++){
    switch(argv[i][1]){
    case 'r':
      run = atoi(&argv[i][2]); // which run period
      break;
    case 'e':
      flag_err = atoi(&argv[i][2]); // error for plotting
      break;
    case 'l':
      lee_strength = atof(&argv[i][2]);
      break;
    case 'd':
      flag_display = atoi(&argv[i][2]);
      break;
    }
  }

  CovMatrix cov;

  // get data histograms ...
  // filetype, period, outfilename, external pot, fileno
  std::map<TString, std::tuple<int, int, TString, float, int, double> > map_inputfile_info = cov.get_map_inputfile_info();

  TFile *temp_file;
  TH1F *htemp;
  TTree *T;
  Double_t pot;
  std::map<TString, std::pair<TH1F*, double> > map_name_histogram;

  // data POT ...
  std::map<int, double> map_data_period_pot;
  //  std::vector<TH1F*> temp_histograms;
  
  // open all the histograms ...
  for (auto it = map_inputfile_info.begin(); it != map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    TString out_filename = std::get<2>(it->second);
    int file_no = std::get<4>(it->second);
    temp_file = new TFile(out_filename);
    T = (TTree*)temp_file->Get("T");
    T->SetBranchAddress("pot",&pot);
    T->GetEntry(0);

    if (filetype==5){
      map_data_period_pot[period] = pot;
    }
    
    // std::cout << pot << std::endl;
    
    std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > all_histo_infos;
    std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = cov.get_histograms(input_filename,0);
    std::copy(histo_infos.begin(), histo_infos.end(), std::back_inserter(all_histo_infos));
    std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos_err2 = cov.get_histograms(input_filename,1);
    std::copy(histo_infos_err2.begin(), histo_infos_err2.end(), std::back_inserter(all_histo_infos));
    std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos_cros = cov.get_histograms(input_filename,2);
    std::copy(histo_infos_cros.begin(), histo_infos_cros.end(), std::back_inserter(all_histo_infos));
    
    for (size_t i=0;i!=all_histo_infos.size();i++){
      htemp = (TH1F*)temp_file->Get(std::get<0>(all_histo_infos.at(i)));
      //      std::cout << std::get<0>(all_histo_infos.at(i)) << " " << htemp->GetSum() << std::endl;
      //      temp_histograms.push_back(htemp);
      map_name_histogram[std::get<0>(all_histo_infos.at(i))] = std::make_pair(htemp, pot);
    }
  }

  
 
  
  // create histograms for data, create histograms for predictions
  // obsch --> histograms (data, prediction, prediction_err2
  std::map<int, std::vector<TH1F*> > map_obsch_histos;
  // Bayesian error needed ...
  // obsch --> bin with overflow bin --> vector of all channels (merge certain channels) --> mean and err2 
  std::map<int, std::vector< std::vector< std::tuple<double, double, double, int, double> > > > map_obsch_bayes;
    
  for (auto it = map_inputfile_info.begin(); it != map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    TString out_filename = std::get<2>(it->second);
    int file_no = std::get<4>(it->second);
    if (filetype == 5){
      // name, nbin, lowlimit, highlimit, variable, channel cut, additional cut, weight
      std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = cov.get_histograms(input_filename,0);
      for (auto it1 = histo_infos.begin(); it1 != histo_infos.end(); it1++){
	int obsch = cov.get_obsch_name(std::get<5>(*it1));	
	htemp = map_name_histogram[std::get<0>(*it1)].first;

	std::vector<TH1F*> vec_histos;
	
	TH1F *hdata = (TH1F*)htemp->Clone(Form("data_%d",obsch));
	hdata->Reset();
	TH1F *hpred = (TH1F*)htemp->Clone(Form("pred_%d",obsch));
	hpred->Reset();
	TH1F *hpred_err2 = (TH1F*)htemp->Clone(Form("pred_err2_%d",obsch));
	hpred_err2->Reset();

	vec_histos.push_back(hdata);
	vec_histos.push_back(hpred);
	vec_histos.push_back(hpred_err2);

	// get histograms ...
	map_obsch_histos[obsch] = vec_histos;
	//map_obsch_bayes[obsch].resize(htemp->GetNbinsX()+1);
	// for (Int_t i=0;i!=htemp->GetNbinsX()+1;i++){
	//   std::vector< std::tuple<double, double, double> > temp;
	//   map_obsch_bayes[obsch].push_back(temp);
	// }
	
	
	//std::cout << std::get<5>(*it1) << " " << obsch << " " << htemp->GetSum() << std::endl;
      }
      
      break;
    }
  }
  
  // get data histograms ...
  cov.fill_data_histograms(run, map_obsch_histos, map_name_histogram);
  
  // get predictions and its uncertainties ...,
  cov.fill_pred_histograms(run, map_obsch_histos, map_obsch_bayes, map_name_histogram, lee_strength, map_data_period_pot);
  
  // get Bayesian errrors ...

  if (flag_err==2){
    std::cout << lee_strength << " " << run << std::endl;
    
    for (auto it = map_obsch_histos.begin(); it!= map_obsch_histos.end(); it++){
      TH1F *h1 = it->second.at(1);
      TH1F *h2 = it->second.at(2);
      int obsch = it->first;

      std::vector<std::vector< std::tuple<double, double, double, int, double> > >  bayes_inputs = map_obsch_bayes[obsch];

      //std::cout << obsch << " " << bayes_inputs.size() << " " << bayes_inputs.at(0).size() << " " << h1->GetNbinsX() << std::endl;

      //if (obsch !=1) continue;
      
      for (int i=0;i!=h1->GetNbinsX()+1;i++){
	Bayes bayes;
	//	if (i!=0) continue;
	//double temp = 0, temp1=0;
	for (auto it1 = bayes_inputs.begin(); it1!=bayes_inputs.end(); it1++){
	  bayes.add_meas_component(std::get<0>((*it1).at(i)), std::get<1>((*it1).at(i)), std::get<2>((*it1).at(i)));
	  // temp += std::get<0>((*it1).at(i));
	  // temp1 += std::get<1>((*it1).at(i));
	  //std::cout << i << " " << std::get<0>((*it1).at(i)) << " " << std::get<1>((*it1).at(i)) << " " << std::get<2>((*it1).at(i)) << " " << std::endl;
	}
	bayes.do_convolution();
	
	double cov = bayes.get_covariance();
	// double cov1 = bayes.get_covariance_mc();
	std::cout << obsch << " " << i << " "	  << h1->GetBinContent(i+1) << " " << cov  << " "  << h2->GetBinContent(i+1) << std::endl;
	// std::cout << temp << " " << temp1 << " " << h1->GetBinContent(i+1) << " " << h2->GetBinContent(i+1) << std::endl;
	h1->SetBinError(i+1,sqrt(cov));
      }
      // obsch --> bin with overflow bin --> vector of all channels (merge certain channels) --> mean and err2 
      //std::map<int, std::vector< std::vector< std::tuple<double, double, double> > > > map_obsch_bayes;
    }
  }

  if (flag_display == 1){
    // plotting ...
    TApplication theApp("theApp",&argc,argv);
    theApp.SetReturnFromRun(true);
    gStyle->SetOptStat(0);
    
    TCanvas c1("ToyMC","ToyMC",2000,800);
    c1.Divide(4,2);
    c1.Draw();
    
    if (flag_err==1){
      
      for (auto it = map_obsch_histos.begin(); it!= map_obsch_histos.end(); it++){
	TH1F *h1 = it->second.at(1);
	TH1F *h2 = it->second.at(2);
	for (int i=0;i!=h1->GetNbinsX()+1;i++){
	  h1->SetBinError(i+1,sqrt(h2->GetBinContent(i+1)));
	}
      }
      
      
      c1.cd(1);
      TGraphErrors *g10 = new TGraphErrors();
      TGraphErrors *g11 = new TGraphErrors();
      for (int i=0;i!=map_obsch_histos[1].at(0)->GetNbinsX()+1;i++){
	double x = map_obsch_histos[1].at(0)->GetBinCenter(i+1);
	double y = map_obsch_histos[1].at(0)->GetBinContent(i+1);
	double x_err = 0;
	double y_err = map_obsch_histos[1].at(0)->GetBinError(i+1);
	g10->SetPoint(i,x,y);
	g10->SetPointError(i,x_err,y_err);
	
	y = map_obsch_histos[1].at(1)->GetBinContent(i+1);
	y_err = map_obsch_histos[1].at(1)->GetBinError(i+1);
	
	g11->SetPoint(i,x,y);
	g11->SetPointError(i,x_err,y_err);
	
	//map_obsch_histos[1].at(0)->Draw();
	//map_obsch_histos[1].at(1)->Draw("same");
	//map_obsch_histos[1].at(1)->SetLineColor(2);
      }
      
      g10->Draw("A*");
      g10->SetMarkerStyle(20);
      g11->Draw("*same");
      g11->SetMarkerStyle(21);
      g11->SetMarkerColor(2);
      g11->SetLineColor(2);
      
      // for (Int_t i=0;i!=map_obsch_histos[1].at(1)->GetNbinsX()+1;i++){
      //   std::cout << i << " " << map_obsch_histos[1].at(1)->GetBinContent(i+1) << " " << map_obsch_histos[1].at(1)->GetBinError(i+1) << std::endl;
      // }
      
      c1.cd(2);
      // map_obsch_histos[3].at(0)->Draw();
      // map_obsch_histos[3].at(1)->Draw("same");
      // map_obsch_histos[3].at(1)->SetLineColor(2);
      
      // for (Int_t i=0;i!=map_obsch_histos[3].at(1)->GetNbinsX()+1;i++){
      //   std::cout << i << " " << map_obsch_histos[3].at(1)->GetBinContent(i+1) << " " << map_obsch_histos[3].at(1)->GetBinError(i+1) << std::endl;
    // }

    TGraphErrors *g30 = new TGraphErrors();
    TGraphErrors *g31 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[3].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[3].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[3].at(0)->GetBinContent(i+1);
      double x_err = 0;
      double y_err = map_obsch_histos[3].at(0)->GetBinError(i+1);
      g30->SetPoint(i,x,y);
      g30->SetPointError(i,x_err,y_err);

      y = map_obsch_histos[3].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[3].at(1)->GetBinError(i+1);

      g31->SetPoint(i,x,y);
      g31->SetPointError(i,x_err,y_err);
      
      //map_obsch_histos[3].at(0)->Draw();
      //map_obsch_histos[3].at(1)->Draw("same");
      //map_obsch_histos[3].at(1)->SetLineColor(2);
    }

    g30->Draw("A*");
    g30->SetMarkerStyle(20);
    g31->Draw("*same");
    g31->SetMarkerStyle(21);
    g31->SetMarkerColor(2);
    g31->SetLineColor(2);

    
    c1.cd(3);
    // map_obsch_histos[5].at(0)->Draw();
    // map_obsch_histos[5].at(1)->Draw("same");
    // map_obsch_histos[5].at(1)->SetLineColor(2);
    
    // for (Int_t i=0;i!=map_obsch_histos[5].at(1)->GetNbinsX()+1;i++){
    //   std::cout << i << " " << map_obsch_histos[5].at(1)->GetBinContent(i+1) << " " << map_obsch_histos[5].at(1)->GetBinError(i+1) << std::endl;
    // }

    TGraphErrors *g50 = new TGraphErrors();
    TGraphErrors *g51 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[5].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[5].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[5].at(0)->GetBinContent(i+1);
      double x_err = 0;
      double y_err = map_obsch_histos[5].at(0)->GetBinError(i+1);
      g50->SetPoint(i,x,y);
      g50->SetPointError(i,x_err,y_err);

      y = map_obsch_histos[5].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[5].at(1)->GetBinError(i+1);

      g51->SetPoint(i,x,y);
      g51->SetPointError(i,x_err,y_err);
      
      //map_obsch_histos[5].at(0)->Draw();
      //map_obsch_histos[5].at(1)->Draw("same");
      //map_obsch_histos[5].at(1)->SetLineColor(2);
    }

    g50->Draw("A*");
    g50->SetMarkerStyle(20);
    g51->Draw("*same");
    g51->SetMarkerStyle(21);
    g51->SetMarkerColor(2);
    g51->SetLineColor(2);

    
    c1.cd(5);
    // map_obsch_histos[2].at(0)->Draw();
    // map_obsch_histos[2].at(1)->Draw("same");
    // map_obsch_histos[2].at(1)->SetLineColor(2);

    TGraphErrors *g20 = new TGraphErrors();
    TGraphErrors *g21 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[2].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[2].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[2].at(0)->GetBinContent(i+1);
      double x_err = 0;
      double y_err = map_obsch_histos[2].at(0)->GetBinError(i+1);
      g20->SetPoint(i,x,y);
      g20->SetPointError(i,x_err,y_err);

      y = map_obsch_histos[2].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[2].at(1)->GetBinError(i+1);

      g21->SetPoint(i,x,y);
      g21->SetPointError(i,x_err,y_err);
      
      //map_obsch_histos[2].at(0)->Draw();
      //map_obsch_histos[2].at(1)->Draw("same");
      //map_obsch_histos[2].at(1)->SetLineColor(2);
    }

    g20->Draw("A*");
    g20->SetMarkerStyle(20);
    g21->Draw("*same");
    g21->SetMarkerStyle(21);
    g21->SetMarkerColor(2);
    g21->SetLineColor(2);
    

    
    c1.cd(6);
    // map_obsch_histos[4].at(0)->Draw();
    // map_obsch_histos[4].at(1)->Draw("same");
    // map_obsch_histos[4].at(1)->SetLineColor(2);

    TGraphErrors *g40 = new TGraphErrors();
    TGraphErrors *g41 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[4].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[4].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[4].at(0)->GetBinContent(i+1);
      double x_err = 0;
      double y_err = map_obsch_histos[4].at(0)->GetBinError(i+1);
      g40->SetPoint(i,x,y);
      g40->SetPointError(i,x_err,y_err);

      y = map_obsch_histos[4].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[4].at(1)->GetBinError(i+1);

      g41->SetPoint(i,x,y);
      g41->SetPointError(i,x_err,y_err);
      
      //map_obsch_histos[4].at(0)->Draw();
      //map_obsch_histos[4].at(1)->Draw("same");
      //map_obsch_histos[4].at(1)->SetLineColor(2);
    }

    g40->Draw("A*");
    g40->SetMarkerStyle(20);
    g41->Draw("*same");
    g41->SetMarkerStyle(21);
    g41->SetMarkerColor(2);
    g41->SetLineColor(2);
    
    c1.cd(7);
    // map_obsch_histos[6].at(0)->Draw();
    // map_obsch_histos[6].at(1)->Draw("same");
    // map_obsch_histos[6].at(1)->SetLineColor(2);

    TGraphErrors *g60 = new TGraphErrors();
    TGraphErrors *g61 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[6].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[6].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[6].at(0)->GetBinContent(i+1);
      double x_err = 0;
      double y_err = map_obsch_histos[6].at(0)->GetBinError(i+1);
      g60->SetPoint(i,x,y);
      g60->SetPointError(i,x_err,y_err);

      y = map_obsch_histos[6].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[6].at(1)->GetBinError(i+1);

      g61->SetPoint(i,x,y);
      g61->SetPointError(i,x_err,y_err);
      
      //map_obsch_histos[6].at(0)->Draw();
      //map_obsch_histos[6].at(1)->Draw("same");
      //map_obsch_histos[6].at(1)->SetLineColor(2);
    }

    g60->Draw("A*");
    g60->SetMarkerStyle(20);
    g61->Draw("*same");
    g61->SetMarkerStyle(21);
    g61->SetMarkerColor(2);
    g61->SetLineColor(2);

    
    c1.cd(4);
    // map_obsch_histos[7].at(0)->Draw();
    // map_obsch_histos[7].at(1)->Draw("same");
    // map_obsch_histos[7].at(1)->SetLineColor(2);

    TGraphErrors *g70 = new TGraphErrors();
    TGraphErrors *g71 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[7].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[7].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[7].at(0)->GetBinContent(i+1);
      double x_err = 0;
      double y_err = map_obsch_histos[7].at(0)->GetBinError(i+1);
      g70->SetPoint(i,x,y);
      g70->SetPointError(i,x_err,y_err);

      y = map_obsch_histos[7].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[7].at(1)->GetBinError(i+1);

      g71->SetPoint(i,x,y);
      g71->SetPointError(i,x_err,y_err);
      
      //map_obsch_histos[7].at(0)->Draw();
      //map_obsch_histos[7].at(1)->Draw("same");
      //map_obsch_histos[7].at(1)->SetLineColor(2);
    }

    g70->Draw("A*");
    g70->SetMarkerStyle(20);
    g71->Draw("*same");
    g71->SetMarkerStyle(21);
    g71->SetMarkerColor(2);
    g71->SetLineColor(2);
    
  }else if (flag_err == 2){

    

    
    c1.cd(1);
    TGraphAsymmErrors *g10 = new TGraphAsymmErrors();
    TGraphErrors *g11 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[1].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[1].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[1].at(0)->GetBinContent(i+1);
      double x_err = 0;
      auto result = cov.get_bayes_errors(y);
      double y_err = 0;
      g10->SetPoint(i,x,y);
      g10->SetPointError(i,0,0,result.first, result.second);
      y = map_obsch_histos[1].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[1].at(1)->GetBinError(i+1);
      g11->SetPoint(i,x,y);
      g11->SetPointError(i,x_err,y_err);
    }
    g10->Draw("A*");
    g10->SetMarkerStyle(20);
    g11->Draw("*same");
    g11->SetMarkerStyle(21);
    g11->SetMarkerColor(2);
    g11->SetLineColor(2);

    c1.cd(5);
    TGraphAsymmErrors *g20 = new TGraphAsymmErrors();
    TGraphErrors *g21 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[2].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[2].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[2].at(0)->GetBinContent(i+1);
      double x_err = 0;
      auto result = cov.get_bayes_errors(y);
      double y_err = 0;
      g20->SetPoint(i,x,y);
      g20->SetPointError(i,0,0,result.first, result.second);
      y = map_obsch_histos[2].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[2].at(1)->GetBinError(i+1);
      g21->SetPoint(i,x,y);
      g21->SetPointError(i,x_err,y_err);
    }
    g20->Draw("A*");
    g20->SetMarkerStyle(20);
    g21->Draw("*same");
    g21->SetMarkerStyle(21);
    g21->SetMarkerColor(2);
    g21->SetLineColor(2);

    c1.cd(2);
    TGraphAsymmErrors *g30 = new TGraphAsymmErrors();
    TGraphErrors *g31 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[3].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[3].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[3].at(0)->GetBinContent(i+1);
      double x_err = 0;
      auto result = cov.get_bayes_errors(y);
      double y_err = 0;
      g30->SetPoint(i,x,y);
      g30->SetPointError(i,0,0,result.first, result.second);
      y = map_obsch_histos[3].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[3].at(1)->GetBinError(i+1);
      g31->SetPoint(i,x,y);
      g31->SetPointError(i,x_err,y_err);
    }
    g30->Draw("A*");
    g30->SetMarkerStyle(20);
    g31->Draw("*same");
    g31->SetMarkerStyle(21);
    g31->SetMarkerColor(2);
    g31->SetLineColor(2);

    c1.cd(6);
    TGraphAsymmErrors *g40 = new TGraphAsymmErrors();
    TGraphErrors *g41 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[4].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[4].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[4].at(0)->GetBinContent(i+1);
      double x_err = 0;
      auto result = cov.get_bayes_errors(y);
      double y_err = 0;
      g40->SetPoint(i,x,y);
      g40->SetPointError(i,0,0,result.first, result.second);
      y = map_obsch_histos[4].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[4].at(1)->GetBinError(i+1);
      g41->SetPoint(i,x,y);
      g41->SetPointError(i,x_err,y_err);
    }
    g40->Draw("A*");
    g40->SetMarkerStyle(20);
    g41->Draw("*same");
    g41->SetMarkerStyle(21);
    g41->SetMarkerColor(2);
    g41->SetLineColor(2);

    c1.cd(3);
    TGraphAsymmErrors *g50 = new TGraphAsymmErrors();
    TGraphErrors *g51 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[5].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[5].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[5].at(0)->GetBinContent(i+1);
      double x_err = 0;
      auto result = cov.get_bayes_errors(y);
      double y_err = 0;
      g50->SetPoint(i,x,y);
      g50->SetPointError(i,0,0,result.first, result.second);
      y = map_obsch_histos[5].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[5].at(1)->GetBinError(i+1);
      g51->SetPoint(i,x,y);
      g51->SetPointError(i,x_err,y_err);
    }
    g50->Draw("A*");
    g50->SetMarkerStyle(20);
    g51->Draw("*same");
    g51->SetMarkerStyle(21);
    g51->SetMarkerColor(2);
    g51->SetLineColor(2);

    c1.cd(7);
    TGraphAsymmErrors *g60 = new TGraphAsymmErrors();
    TGraphErrors *g61 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[6].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[6].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[6].at(0)->GetBinContent(i+1);
      double x_err = 0;
      auto result = cov.get_bayes_errors(y);
      double y_err = 0;
      g60->SetPoint(i,x,y);
      g60->SetPointError(i,0,0,result.first, result.second);
      y = map_obsch_histos[6].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[6].at(1)->GetBinError(i+1);
      g61->SetPoint(i,x,y);
      g61->SetPointError(i,x_err,y_err);
    }
    g60->Draw("A*");
    g60->SetMarkerStyle(20);
    g61->Draw("*same");
    g61->SetMarkerStyle(21);
    g61->SetMarkerColor(2);
    g61->SetLineColor(2);

    c1.cd(4);
    TGraphAsymmErrors *g70 = new TGraphAsymmErrors();
    TGraphErrors *g71 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[7].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[7].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[7].at(0)->GetBinContent(i+1);
      double x_err = 0;
      auto result = cov.get_bayes_errors(y);
      double y_err = 0;
      g70->SetPoint(i,x,y);
      g70->SetPointError(i,0,0,result.first, result.second);
      y = map_obsch_histos[7].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[7].at(1)->GetBinError(i+1);
      g71->SetPoint(i,x,y);
      g71->SetPointError(i,x_err,y_err);
    }
    g70->Draw("A*");
    g70->SetMarkerStyle(20);
    g71->Draw("*same");
    g71->SetMarkerStyle(21);
    g71->SetMarkerColor(2);
    g71->SetLineColor(2);
    
  }
  theApp.Run();
  }


  
  
  // std::map<TString, std::pair<TString, int> > map_pred_histo_hist_err2_lee = cov.get_map_pred_histo_histo_err2_lee();
  // std::map<std::pair<TString, TString>, std::pair<TString, int> > map_pair_histo_histos_cross = cov.get_map_pair_hist_hists_cros();
  // std::map<int, std::set<std::set<std::pair<TString, int> > > > map_pred_obsch_histos = cov.get_map_pred_obsch_histos();

  if (flag_err == 1){
    // prediction ...
    TFile *file3 = new TFile("merge.root","RECREATE");
    
    TMatrixD* mat_collapse = cov.get_mat_collapse();
    mat_collapse->Write("mat_collapse");
    
    // additional covariance matrix ...
    TMatrixD* mat_add_cov = cov.get_add_cov_matrix();
    
    std::map<int, TH1F*> map_covch_histo;
    
    for (auto it = map_obsch_bayes.begin(); it != map_obsch_bayes.end(); it++){
      int obsch = it->first;
      std::vector<std::vector< std::tuple<double, double, double, int, double> > >  bayes_inputs = it->second;
      TH1F *htemp = (TH1F*)map_obsch_histos[obsch].at(0); // data histogram ...
      
      for (size_t i=0;i!=bayes_inputs.size(); i++){
	int covch = std::get<3>(bayes_inputs.at(i).front());
	// double add_sys = std::get<4>(bayes_inputs.at(i).front());
	int start_bin = cov.get_covch_startbin(covch);
	
	if (map_covch_histo.find(covch) == map_covch_histo.end()){
	  TH1F *hnew = (TH1F*)htemp->Clone(Form("histo_%d",covch));
	  hnew->Reset();
	  map_covch_histo[covch] = hnew;
      }
	TH1F *htemp1 = map_covch_histo[covch];
	for (size_t j=0;j!=bayes_inputs.at(i).size();j++){
	  htemp1->SetBinContent(j+1, htemp1->GetBinContent(j+1) + std::get<0>(bayes_inputs.at(i).at(j)));
	  htemp1->SetBinError(j+1, sqrt(pow(htemp1->GetBinError(j+1),2) + std::get<1>(bayes_inputs.at(i).at(j))));
	  
	  (*mat_add_cov)(start_bin + j, start_bin + j) += std::get<4>(bayes_inputs.at(i).at(j));
	  
	  //	if (std::get<4>(bayes_inputs.at(i).at(j)) > 0)
	  //  std::cout << obsch << " " << i << " " << start_bin + j << " " <<  std::get<4>(bayes_inputs.at(i).at(j)) << std::endl;
	}
	//      std::cout << obsch << " " << bayes_inputs.size() << " " << i << " " << covch << " " << start_bin << std::endl;
	//  break;
      }
    }
    
    
    for (auto it = map_covch_histo.begin(); it != map_covch_histo.end(); it++){
      it->second->SetDirectory(file3);
    }
    
    for (auto it = map_obsch_histos.begin(); it != map_obsch_histos.end(); it++){
      int obsch = it->first;
      TH1F *hdata = it->second.at(0);
      TH1F *hmc = it->second.at(1);
      
      hdata->SetName(Form("hdata_obsch_%d",obsch));
    hmc->SetName(Form("hmc_obsch_%d",obsch));
    hdata->SetDirectory(file3);
    hmc->SetDirectory(file3);
    }
    
    mat_add_cov->Write("cov_mat_add");
    
    file3->Write();
    file3->Close();
  }  

}
