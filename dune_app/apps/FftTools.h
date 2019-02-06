// Functions for interpolation through FFT

void FftShiftSampling(std::vector<double> vin,
                      double smpLen, // maximum in x-axis, not bin number
                      double toffset,
                      std::vector<double>& vout){

  int inSmps = vin.size();
  std::vector<double> vRe(inSmps), vIm(inSmps);
  TVirtualFFT* fft_r2c = TVirtualFFT::FFT(1, &inSmps, "R2C ES K");
  fft_r2c->SetPoints(vin.data());  
  fft_r2c->Transform();
  fft_r2c->GetPointsComplex(vRe.data(), vIm.data());

  // FFT phase shift
  double f0 = 1./smpLen;
  for(int i=0; i<inSmps; i++){
    double fi = i * f0;
    double w = 2*TMath::Pi()*fi;
    TComplex z(vRe[i], vIm[i]);
    TComplex z1(0,w*toffset); // phase shift
    TComplex z2 = z * TComplex::Exp(z1);
    vRe[i] = z2.Re();
    vIm[i] = z2.Im();
  }

  std::vector<double> vReOut(inSmps);
  TVirtualFFT* fft_c2r = TVirtualFFT::FFT(1, &inSmps, "C2R M K");
  fft_c2r->SetPointsComplex(vRe.data(), vIm.data());
  fft_c2r->Transform();
  fft_c2r->GetPoints(vReOut.data());

  vout.clear();
  for(int i=0; i<inSmps; i++){
    vout.push_back(vReOut[i] / (1.0*inSmps));
  }

  delete fft_r2c;
  delete fft_c2r;
}

void FftUpSampling(std::vector<double> vin,
                  int inSmps1,
                  std::vector<double>& vout,
                  int outSmps1){

  int inSmps = inSmps1;
  int outSmps = outSmps1;
  vout.resize(outSmps1);

  if(outSmps<=inSmps){
  	for(int i=0; i<outSmps; i++){
  		if(i<inSmps) vout[i] = vin[i];
  		else vout[i] = vin[inSmps-1];
  	}
  	std::cout << "Warning: output sample size should be greater than the input." << std::endl;
  	return;
  }

  std::vector<double> vRe(inSmps), vIm(inSmps);
  TVirtualFFT* fft_r2c = TVirtualFFT::FFT(1, &inSmps, "R2C ES K");
  fft_r2c->SetPoints(vin.data());  
  fft_r2c->Transform();
  fft_r2c->GetPointsComplex(vRe.data(), vIm.data());

  // extend sample size
  vRe.resize(outSmps);
  vIm.resize(outSmps);

  // inverse FFT
  TVirtualFFT* fft_c2r = TVirtualFFT::FFT(1,&outSmps, "C2R M K");
  fft_c2r->SetPointsComplex(vRe.data(), vIm.data()); 
  fft_c2r->Transform();
  std::vector<double> vReOut(outSmps);
  fft_c2r->GetPoints(vReOut.data());

  for(int i=0; i<outSmps; i++){
    vout[i] = vReOut[i] / (1.0*inSmps);
  }

  delete fft_r2c;
  delete fft_c2r;
}
