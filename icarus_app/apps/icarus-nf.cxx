#include <iostream>
#include <vector>
#include <set>
#include <algorithm>

#include "WCP2dToy/ExecMon.h"

#include "TFile.h"
#include "TH2F.h"
#include "TString.h"
#include "TMatrixD.h"
#include "TVirtualFFT.h"

#include <Eigen/IterativeLinearSolvers>

using namespace std;

int main(int argc, char* argv[])
{
  if (argc < 3) {
    cerr << "usage: icarus-coherent #filename #outfile" << endl;
  }else{

    WCP::ExecMon em("starting");
    TString filename = argv[1];
    TFile *file = new TFile(filename);
    
    TString outfile = argv[2];
    
    const int nticks = 4096;

    TH2F **hu = new TH2F*[4];
    hu[0] = (TH2F*)file->Get("hu_orig111"); // 1056 chs, from 13824 to 14880
    hu[1] = (TH2F*)file->Get("hu_orig120"); // 1056 chs, from 1152 to 2208
    hu[2] = (TH2F*)file->Get("hu_orig110"); // 1056 chs, from 0 to 1056
    hu[3] = (TH2F*)file->Get("hu_orig121"); // 1056 chs, from 14976 to 16032

    
    TFile *file3 = new TFile("processed_files/sum_all_28.root");
    // TFile *file3 = new TFile("processed_files/sum_0.root");
    TMatrixD *mat = (TMatrixD*)file3->Get("mat");
    //    TH1F *h_e = (TH1F*)file3->Get("h_e");
    //    TFile *file3 = new TFile("temp.root");
    //    TMatrixD *mat = (TMatrixD*)file3->Get("mat_u1");
    //    for (Int_t i=0;i!=4096;i++){
    //  for (Int_t j=i+1;j!=4096;j++){
    //	(*mat)(j,i) = (*mat)(i,j);
    //  }
    // }
    
    double max_val = 0;
    for (Int_t i=0;i!=4096;i++){
      if ((*mat)(i,i) > max_val) max_val = (*mat)(i,i);
    }
    std::cout <<"Diagonal Maximum: " << max_val << std::endl;
    double limit = 0.0002;
    double limit1 = 0.05;

    // double limit = 1.0;
    // double limit1 = 1.0;

    std::set<int> no_list;
    for (Int_t i=0;i!=4096;i++){
      if ((*mat)(i,i) <= max_val * limit || i >=2048){
	for (Int_t j=0;j!=4096;j++){
	  (*mat)(i,j) = 0;
	  (*mat)(j,i) = 0;
	}
	no_list.insert(i);
      }
    }
    
    std::cout << "Zero List Size: " << no_list.size() << " " << std::endl;
    
    for (Int_t i=0;i!=4096;i++){
      if ((*mat)(i,i)==0) continue;
      double val1 = sqrt((*mat)(i,i));
      for (Int_t j=i+1;j!=4096;j++){
	if ((*mat)(j,j) ==0) continue;
	double val2 = sqrt((*mat)(j,j));
	if (fabs((*mat)(i,j)) < limit1 * val1 * val2 ||
	    fabs((*mat)(j,i)) < limit1 * val1 * val2 ){
	  (*mat)(j,i) = 0;
	  (*mat)(i,j) = 0;
	}
      }
    }
    int count_non_zero =0;
    for (Int_t i=0;i!=4096;i++){
      for (Int_t j=0;j!=4096;j++){
	if ((*mat)(i,j) !=0) count_non_zero ++;
      }
    }
    std::cout << "Non-zero element: " << count_non_zero << std::endl;
    
        
    TH2F **hu1 = new TH2F*[4];
    for (Int_t i=0;i!=4;i++){
      hu1[i] = (TH2F*)hu[i]->Clone(Form("hu_%d",i));
      hu1[i]->Reset();
    }

    cout << em("open file") << std::endl;
    
    int n_u = 33;
    
    double lambda = 0.05;//0.05; // regularization strength ...
    //    double lambda1 = 0.05*400;//0.05; // regularization strength ...
    const int nsize = 4096;
    double lambda_n = nsize*32; //(nsize-no_list.size())*32; // noise part ...

    
    Eigen::VectorXd b(32*nsize+nsize); // measurement matrix
    Eigen::VectorXd b_2nd(32*nsize+nsize); // measurement matrix
    Eigen::VectorXd temp(nsize);
    Eigen::SparseMatrix<double> A(32*nsize+nsize,32*nsize+nsize); // cofficient ...
    Eigen::SparseMatrix<double> A_2nd(32*nsize+nsize,32*nsize+nsize); // cofficient ...
    Eigen::VectorXd result_init(32*nsize+nsize); // solution ...
    Eigen::VectorXd result(32*nsize+nsize); // solution ...

    Eigen::MatrixXd cov(nsize,nsize); // covariance matrix
    for (Int_t i=0;i!=nsize;i++){
      for (Int_t j=0;j!=nsize;j++){
	cov(i,j) = (*mat)(i,j);
      }
    }
    result_init.setZero();
    b.setZero();
    b_2nd.setZero();

    TH1F *ht = new TH1F("ht","ht",nticks,0,nticks);
    TH1 *h_real = 0;
    TH1 *h_imag = 0;

    TVirtualFFT::SetTransform(0);
    Int_t n = nticks;
    TVirtualFFT *ifft2 = TVirtualFFT::FFT(1,&n,"C2R M K");
    double temp_re[nticks],temp_im[nticks];
    TH1 *fb = 0;
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    std::vector<T> tripletList_2nd;
    tripletList.reserve(nsize*32*4);
    tripletList_2nd.reserve(count_non_zero * 35);
    
    cout << em("initiailization") << std::endl;
    
    for (Int_t i=28;i!=29;i++){

      Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
      Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver_2nd;
      
      for (Int_t j=0;j!=32;j++){
	for (Int_t k=0;k!=nticks;k++){
	  ht->SetBinContent(k+1,hu[2]->GetBinContent(i*32+j+1,k+1));	  
	}
	h_real = ht->FFT(0,"RE");
	h_imag = ht->FFT(0,"IM");

	// temp ...
	for (Int_t k=0;k!=nsize/2;k++){
	  temp(2*k)   = h_real->GetBinContent(k+1);
	  temp(2*k+1) = h_imag->GetBinContent(k+1);
	}
	Eigen::VectorXd temp1 = cov * temp;
	//Eigen::VectorXd temp1 = temp;

	// initialization ...
	for (Int_t k=0;k!=nsize;k++){
	  b(nsize*j+k) = temp(k);
	  b_2nd(nsize*j+k) = temp(k);
	  
	  result_init(nsize*j+k) = temp(k);
	  
	  tripletList.push_back(T(nsize*j+k,nsize*j+k,1.0 + lambda));
	  tripletList.push_back(T(nsize*j+k,nsize*32+k,1.0));

	  tripletList_2nd.push_back(T(nsize*j+k,nsize*j+k,1.0 + lambda));
	  tripletList_2nd.push_back(T(nsize*j+k,nsize*32+k,1.0));
	  
	  
	  b(32*nsize+k) += temp(k);
	  tripletList.push_back(T(nsize*32+k,nsize*j+k,1));
	  
	  if (no_list.find(k) == no_list.end()){
	    b_2nd(32*nsize+k) += temp1(k);
	    // last equation ...
	    for (Int_t q=0;q!=nsize;q++){
	      if (cov(k,q)!=0) {
		tripletList_2nd.push_back(T(nsize*32+k,nsize*j+q,cov(k,q)));
		//std::cout << 32 << " " << k << " " << j << " " << q << " " << cov(k,q) << std::endl;
	      }
	    }
	  }else{
	    b_2nd(32*nsize+k) += temp(k);
	    tripletList_2nd.push_back(T(nsize*32+k,nsize*j+k,1));
	  }
	}
	delete h_real;
	delete h_imag;
      } // loop over 32 channels ...
      
      // last piece ...
      for (Int_t k=0;k!=nsize;k++){
	tripletList.push_back(T(nsize*32+k,nsize*32+k, 32));
	
	if (no_list.find(k) == no_list.end()){
	  tripletList_2nd.push_back(T(nsize*32+k,nsize*32+k,lambda_n + 32*cov(k,k)));
	}else{
	  tripletList_2nd.push_back(T(nsize*32+k,nsize*32+k, 32));
	}
	
       	for (Int_t q=k+1;q!=nsize;q++){
	  if (no_list.find(q) == no_list.end()){
	    if (cov(q,k) != 0){
	      tripletList_2nd.push_back(T(nsize*32+q,nsize*32+k,32*cov(q,k)));
	      tripletList_2nd.push_back(T(nsize*32+k,nsize*32+q,32*cov(k,q)));
	    }
	  }
	}
      }
      A.setFromTriplets(tripletList.begin(), tripletList.end());
      A_2nd.setFromTriplets(tripletList_2nd.begin(), tripletList_2nd.end());

      
      
      //std::cout << tripletList.size() << " " << count_non_zero * 35 << std::endl;
      
      cout << em("loading content") << std::endl;
      
      solver.compute(A);
      result = solver.solveWithGuess(b,result_init);
      result_init = result;


      cout << em("solve first round") << std::endl;

      solver_2nd.setMaxIterations(200);
      solver_2nd.setTolerance(1e-3);
      
      solver_2nd.compute(A_2nd);
      result = solver_2nd.solveWithGuess(b_2nd,result_init);

      std::cout << "#interations:    " << solver_2nd.iterations() << std::endl;
      std::cout << "estimated error: " << solver_2nd.error() << std::endl;
      
      cout << em("solve 2nd round") << std::endl;


      
      if (std::isnan(solver_2nd.error()))  result = result_init;
      
      
      
      for (Int_t j=0;j!=32;j++){
	for (Int_t k=0;k!=nsize/2;k++){
	  temp_re[k] = result(nsize*j+2*k)/n;
	  temp_im[k] = result(nsize*j+2*k+1)/n;
	}
	ifft2->SetPointsComplex(temp_re,temp_im);
	ifft2->Transform();	
	fb = TH1::TransformHisto(ifft2,fb,"Re");
	for (Int_t k=0;k!=nsize;k++){
	  hu1[2]->SetBinContent(i*32+j+1,k+1, fb->GetBinContent(k+1));
	}
      }
      
    }

    
    
    
    TFile *file2 = new TFile(outfile,"RECREATE");
    for (Int_t i=0;i!=4;i++){
      hu[i]->SetDirectory(file2);
      hu1[i]->SetDirectory(file2);
    }
    //    mat->Write("mat");
    file2->Write();
    file2->Close();



    
  }
}
