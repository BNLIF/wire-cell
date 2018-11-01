void plot(){
  TFile* file = TFile::Open("mcs-tracks.root");
  TTree* T = (TTree*)file->Get("T");

  int N;
  std::vector<double>* x=0;
  std::vector<double>* y=0;
  std::vector<double>* z=0;
  T->SetBranchAddress("N", &N);
  T->SetBranchAddress("x", &x);
  T->SetBranchAddress("y", &y);
  T->SetBranchAddress("z", &z);
  
  TGraph2D* g1 = new TGraph2D();

  for(int I=0; I<50 /*tracks*/; I++){ 
    T->GetEntry(I);
    std::cout << N << " vertices"  << std::endl;

    for(int i=0; i<N; i++){
      g1->SetPoint(g1->GetN(), x->at(i), y->at(i), z->at(i));    
      //std::cout << x->at(i) << " " << y->at(i) << " " << z->at(i) << std::endl;
    }
  }
  g1->Draw("AP");
}
