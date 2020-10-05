#include <iostream>

#include "WCPLEEANA/master_cov_matrix.h"


#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
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
    }
  }

  CovMatrix cov;

  // get data histograms ...
  // filetype, period, outfilename, external pot, fileno
  std::map<TString, std::tuple<int, int, TString, float, int> > map_inputfile_info = cov.getp_map_inputfile_info();

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
  std::map<int, std::vector< std::vector< std::pair<double, double> > > > map_obsch_bayes;
    
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
	for (Int_t i=0;i!=htemp->GetNbinsX()+1;i++){
	  std::vector< std::pair<double, double> > temp;
	  map_obsch_bayes[obsch].push_back(temp);
	}
	
	
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

  for (auto it = map_obsch_histos.begin(); it!= map_obsch_histos.end(); it++){
    TH1F *h1 = it->second.at(1);
    TH1F *h2 = it->second.at(2);
    for (int i=0;i!=h1->GetNbinsX()+1;i++){
      h1->SetBinError(i+1,sqrt(h2->GetBinContent(i+1)));
    }
  }

  // plotting ...
  TApplication theApp("theApp",&argc,argv);
  theApp.SetReturnFromRun(true);
  gStyle->SetOptStat(0);

  TCanvas c1("ToyMC","ToyMC",2000,800);
  c1.Divide(4,2);
  c1.Draw();
  c1.cd(1);
  map_obsch_histos[1].at(0)->Draw();
  map_obsch_histos[1].at(1)->Draw("same");
  map_obsch_histos[1].at(1)->SetLineColor(2);

  // for (Int_t i=0;i!=map_obsch_histos[1].at(1)->GetNbinsX()+1;i++){
  //   std::cout << i << " " << map_obsch_histos[1].at(1)->GetBinContent(i+1) << " " << map_obsch_histos[1].at(1)->GetBinError(i+1) << std::endl;
  // }
  
  c1.cd(2);
  map_obsch_histos[3].at(0)->Draw();
  map_obsch_histos[3].at(1)->Draw("same");
  map_obsch_histos[3].at(1)->SetLineColor(2);

  // for (Int_t i=0;i!=map_obsch_histos[3].at(1)->GetNbinsX()+1;i++){
  //   std::cout << i << " " << map_obsch_histos[3].at(1)->GetBinContent(i+1) << " " << map_obsch_histos[3].at(1)->GetBinError(i+1) << std::endl;
  // }
  
  c1.cd(3);
  map_obsch_histos[5].at(0)->Draw();
  map_obsch_histos[5].at(1)->Draw("same");
  map_obsch_histos[5].at(1)->SetLineColor(2);

  // for (Int_t i=0;i!=map_obsch_histos[5].at(1)->GetNbinsX()+1;i++){
  //   std::cout << i << " " << map_obsch_histos[5].at(1)->GetBinContent(i+1) << " " << map_obsch_histos[5].at(1)->GetBinError(i+1) << std::endl;
  // }
  
  c1.cd(5);
  map_obsch_histos[2].at(0)->Draw();
  map_obsch_histos[2].at(1)->Draw("same");
  map_obsch_histos[2].at(1)->SetLineColor(2);
  
  c1.cd(6);
  map_obsch_histos[4].at(0)->Draw();
  map_obsch_histos[4].at(1)->Draw("same");
  map_obsch_histos[4].at(1)->SetLineColor(2);
  
  c1.cd(7);
  map_obsch_histos[6].at(0)->Draw();
  map_obsch_histos[6].at(1)->Draw("same");
  map_obsch_histos[6].at(1)->SetLineColor(2);
  
  c1.cd(4);
  map_obsch_histos[7].at(0)->Draw();
  map_obsch_histos[7].at(1)->Draw("same");
  map_obsch_histos[7].at(1)->SetLineColor(2);

  theApp.Run();


  
  // std::map<TString, std::pair<TString, int> > map_pred_histo_hist_err2_lee = cov.get_map_pred_histo_histo_err2_lee();
  // std::map<std::pair<TString, TString>, std::pair<TString, int> > map_pair_histo_histos_cross = cov.get_map_pair_hist_hists_cros();
  // std::map<int, std::set<std::set<std::pair<TString, int> > > > map_pred_obsch_histos = cov.get_map_pred_obsch_histos();

  // // prediction ...
 

  

}
