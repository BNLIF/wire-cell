#include <iostream>
#include <vector>
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

    TMatrixD *mat = (TMatrixD*)file->Get("mat");
    
    TH2F **hu1 = new TH2F*[4];
    for (Int_t i=0;i!=4;i++){
      hu1[i] = (TH2F*)hu[i]->Clone(Form("hu_%d",i));
      hu1[i]->Reset();
    }

    cout << em("open file") << std::endl;
    
    int n_u = 33;
    
    float lambda = 0.05; // regularization strength ...
    float lambda_n = 0.;//32*32*1000; // noise part ...
    
    Eigen::VectorXf b(32*4096+4096); // measurement matrix
    Eigen::VectorXf temp(4096);
    Eigen::SparseMatrix<float> A(32*4096+4096,32*4096+4096); // cofficient ...
    Eigen::VectorXf result_init(32*4096+4096); // solution ...
    Eigen::VectorXf result(32*4096+4096); // solution ...

    Eigen::MatrixXf cov(4096,4096); // covariance matrix
    for (Int_t i=0;i!=4096;i++){
      for (Int_t j=0;j!=4096;j++){
	cov(i,j) = (*mat)(i,j);
      }
    }
    result_init.setZero();
    b.setZero();
    

    TH1F *ht = new TH1F("ht","ht",4096,0,4096);
    TH1 *h_real = 0;
    TH1 *h_imag = 0;

    TVirtualFFT::SetTransform(0);
    Int_t n = nticks;
    TVirtualFFT *ifft2 = TVirtualFFT::FFT(1,&n,"C2R M K");
    double temp_re[nticks],temp_im[nticks];
    TH1 *fb = 0;
    
    typedef Eigen::Triplet<float> T;
    std::vector<T> tripletList;
    tripletList.reserve(4096*4096*35);

    cout << em("initiailization") << std::endl;
    
    for (Int_t i=28;i!=29;i++){

      Eigen::BiCGSTAB<Eigen::SparseMatrix<float>> solver;
      for (Int_t j=0;j!=32;j++){
	for (Int_t k=0;k!=4096;k++){
	  ht->SetBinContent(k+1,hu[2]->GetBinContent(i*32+j+1,k+1));	  
	}
	h_real = ht->FFT(0,"RE");
	h_imag = ht->FFT(0,"IM");

	// temp ...
	for (Int_t k=0;k!=2048;k++){
	  temp(2*k)   = h_real->GetBinContent(k+1);
	  temp(2*k+1) = h_imag->GetBinContent(k+1);
	}
	//Eigen::VectorXf temp1 = cov * temp;
	Eigen::VectorXf temp1 = temp;

	// initialization ...
	for (Int_t k=0;k!=4096;k++){
	  b(4096*j+k) = temp(k);
	  result_init(4096*j+k) = temp(k);
	  b(32*4096+k) += temp1(k);

	  
	  // A.insert(4096*j+k,4096*j+k) = 1.0 + lambda;
	  // A.insert(4096*j+k,4096*32+k) = 1.0;
	  tripletList.push_back(T(4096*j+k,4096*j+k,1.0 + lambda));
	  tripletList.push_back(T(4096*j+k,4096*32+k,1.0));
	  
	  // last equation ...
	  for (Int_t q=0;q!=4096;q++){
	    //   A.insert(4096*32+k,4096*j+q) = cov(k,q);
	    //tripletList.push_back(T(4096*32+k,4096*j+q,cov(k,q)));
	  }
	  tripletList.push_back(T(4096*32+k,4096*j+k,1));
	}
	delete h_real;
	delete h_imag;
      } // loop over 32 channels ...
      
      // last piece ...
      for (Int_t k=0;k!=4096;k++){
       	for (Int_t q=0;q!=4096;q++){
       	  if (q==k){
      // 	    A.insert(4096*32+q,4096*32+k) = 1+cov(q,k);
	    //tripletList.push_back(T(4096*32+q,4096*32+k,lambda_n + 32*cov(q,k)));
	    tripletList.push_back(T(4096*32+q,4096*32+k,lambda_n + 32));
       	  }else{
      // 	    A.insert(4096*32+q,4096*32+k) = cov(q,k);
	    //	    tripletList.push_back(T(4096*32+q,4096*32+k,32*cov(q,k)));
       	  }
       	}
      }
      A.setFromTriplets(tripletList.begin(), tripletList.end());

      cout << em("loading content") << std::endl;
      
      solver.compute(A);

      cout << em("prepare solver") << std::endl;
     
      result = solver.solveWithGuess(b,result_init);

      cout << em("solve one") << std::endl;

      for (Int_t j=0;j!=32;j++){
	for (Int_t k=0;k!=n/2;k++){
	  temp_re[k] = result(4096*j+2*k)/n;
	  temp_im[k] = result(4096*j+2*k+1)/n;
	}
	ifft2->SetPointsComplex(temp_re,temp_im);
	ifft2->Transform();	
	fb = TH1::TransformHisto(ifft2,fb,"Re");
	for (Int_t k=0;k!=4096;k++){
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
