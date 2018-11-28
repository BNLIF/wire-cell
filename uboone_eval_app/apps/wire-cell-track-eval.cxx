#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPolyLine3D.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TMath.h"

#include <iostream>
#include <map>

#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVector3.h"

#include "TStyle.h"
#include "TColor.h"
#include "TROOT.h"

using namespace std;

void set_plot_style()
{
  TStyle* myStyle = new TStyle("myStyle","My ROOT plot style");
  // plot style
  //myStyle->SetPalette(kInvertedDarkBodyRadiator); 
  //myStyle->SetPalette(kRainBow);
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  myStyle->SetNumberContours(NCont);
  myStyle->SetOptStat(0);
    
  myStyle->SetLabelFont(62,"xyz");
  myStyle->SetLabelSize(0.05,"xyz");
  myStyle->SetTitleFont(62, "xyz");
  myStyle->SetTitleSize(0.06,"xyz");
  myStyle->SetTitleOffset(1.3, "y");
  myStyle->SetTitleOffset(0.8, "x");
  myStyle->SetTitleOffset(1.5, "z");
    
  // only 5 in x to avoid label overlaps
  myStyle->SetNdivisions(505, "x");

  // set the margin sizes
  myStyle->SetPadTopMargin(0.05);
  myStyle->SetPadRightMargin(0.2); // increase for colz plots
  myStyle->SetPadBottomMargin(0.1);
  myStyle->SetPadLeftMargin(0.15);
  myStyle->SetFrameBorderMode(0);
  myStyle->SetCanvasBorderMode(0);
  myStyle->SetPadBorderMode(0);
  myStyle->SetCanvasColor(0);
  myStyle->SetPadColor(0);

  gROOT->SetStyle("myStyle");
  gROOT->ForceStyle();
}

int main(int argc, char* argv[])
{
  if(argc < 3){
    cout<<"Usage: wire-cell-track-eval input.root outputname.root [option: cluster_id]"<<endl;
    exit(1);
  }
  const char* inputroot = argv[1];
  const char* outputroot = argv[2];
  int cluster_check=-1; // specify cluster_id to check 
  if(argc==4) cluster_check=atoi(argv[3]);


  // 2D plot, for each wire plane's coordinate, theta_y: 0-pi/2, theta_xz (phi): 0-pi 
  TH2D* hntrack = new TH2D("hntrack","number of tracks",10,0,1,10,0,TMath::Pi());
  TH2D* hltrack = new TH2D("hltrack","length of tracks",10,0,1,10,0,TMath::Pi());
  TH2D* hntrack_true = new TH2D("hntrack_true","number of tracks",10,0,1,10,0,TMath::Pi());
  TH2D* hltrack_true = new TH2D("hltrack_true","length of tracks",10,0,1,10,0,TMath::Pi());

  // rotation matrix to convert collection plane (default detector) coordinate to inductions plane
  TMatrixD Muplane(3,3);
  TMatrixD Mvplane(3,3);
  double angle = 60./180*TMath::Pi();
  double cosine = TMath::Cos(angle);
  double sine = TMath::Sin(angle);
  Muplane(0,0) = 1.;
  Muplane(0,1) = 0;
  Muplane(0,2) = 0;
  Muplane(1,0) = 0;
  Muplane(1,1) = cosine;
  Muplane(1,2) = sine;
  Muplane(2,0) = 0;
  Muplane(2,1) = -sine;
  Muplane(2,2) = cosine;

  Mvplane(0,0) = 1.;
  Mvplane(0,1) = 0;
  Mvplane(0,2) = 0;
  Mvplane(1,0) = 0;
  Mvplane(1,1) = cosine;
  Mvplane(1,2) = -sine;
  Mvplane(2,0) = 0;
  Mvplane(2,1) = sine;
  Mvplane(2,2) = cosine;


  TH2D* hntrack_u = new TH2D("hntrack_u","number of tracks",10,0,1,10,0,TMath::Pi());
  TH2D* hltrack_u = new TH2D("hltrack_u","length of tracks",10,0,1,10,0,TMath::Pi());
  TH2D* hntrack_true_u = new TH2D("hntrack_true_u","number of tracks",10,0,1,10,0,TMath::Pi());
  TH2D* hltrack_true_u = new TH2D("hltrack_true_u","length of tracks",10,0,1,10,0,TMath::Pi());

  TH2D* hntrack_v = new TH2D("hntrack_v","number of tracks",10,0,1,10,0,TMath::Pi());
  TH2D* hltrack_v = new TH2D("hltrack_v","length of tracks",10,0,1,10,0,TMath::Pi());
  TH2D* hntrack_true_v = new TH2D("hntrack_true_v","number of tracks",10,0,1,10,0,TMath::Pi());
  TH2D* hltrack_true_v = new TH2D("hltrack_true_v","length of tracks",10,0,1,10,0,TMath::Pi());

    
  // check
  TH3D* hrec = new TH3D("hrec","",25,0,256,25,-115,117,100,0,1037);
  TPolyLine3D *leval1 = new TPolyLine3D(2);
  TPolyLine3D *leval2 = new TPolyLine3D(2);

  // read in the point cloud for each cluster_id
  TFile* f = new TFile(inputroot, "READ");
  TTree* clusters = (TTree*)f->Get("T_cluster");
  Double_t x=0;
  Double_t y=0;
  Double_t z=0;
  Int_t cluster_id=0;

  clusters->SetBranchAddress("x",&x);
  clusters->SetBranchAddress("y",&y);
  clusters->SetBranchAddress("z",&z);
  clusters->SetBranchAddress("cluster_id",&cluster_id);

  std::map<int, std::vector<double>> xpt;
  std::map<int, std::vector<double>> ypt;
  std::map<int, std::vector<double>> zpt;

  for(int i=0; i<clusters->GetEntries();i++)
    {
      clusters->GetEntry(i);
      auto it = xpt.find(cluster_id);
      if( it == xpt.end() ){
	std::vector<double> xvec;
	xvec.push_back(x);
	std::vector<double> yvec;
	yvec.push_back(y);
	std::vector<double> zvec;
	zvec.push_back(z);

	xpt[cluster_id] = xvec;
	ypt[cluster_id] = yvec;
	zpt[cluster_id] = zvec;
      }
      else{
	xpt[cluster_id].push_back(x);
	ypt[cluster_id].push_back(y);
	zpt[cluster_id].push_back(z);
      }
    }
  //cout<<clusters->GetEntries()<<endl;

  // PCA (Priciple Component Analysis): Matrix Decomposition --> Coordinate Rotation
  // n-dimentional points in nxn position covariance matrix;
  // diagonalization (rotation);
  // the biggest eigen value corresponds to the main axis --> primary direction of the point cloud
    
  // Step 1: calculate average position
  // Step 2: construct covariance matrix
  // merge the two steps into a single loop of the point cloud

  //cout<<xpt.size()<<endl;
  for(int i=0; i<xpt.size(); i++)
    { // each cluster
      double npts = xpt[i].size();
      double xx=0;
      double yy=0;
      double zz=0;
      double xy=0;
      double yz=0;
      double zx=0;
      double mx=0;
      double my=0;
      double mz=0;
      for(int j=0; j<npts; j++)
        { // each point
	  double x = xpt[i].at(j);
	  double y = ypt[i].at(j);
	  double z = zpt[i].at(j);
	  //check
	  if(i==cluster_check) hrec->Fill(x, y, z);

	  xx += x*x/npts;
	  yy += y*y/npts;
	  zz += z*z/npts;
	  xy += x*y/npts;
	  yz += y*z/npts;
	  zx += z*x/npts;
	  mx += x/npts;
	  my += y/npts;
	  mz += z/npts;
        } //each point   

      TMatrixDSym m(3);
      m(0,0) = xx - mx*mx; 
      m(0,1) = xy - mx*my;
      m(0,2) = zx - mz*mx;
      m(1,0) = xy - mx*my;
      m(1,1) = yy - my*my;
      m(1,2) = yz - my*mz;
      m(2,0) = zx - mz*mx;
      m(2,1) = yz - my*mz;
      m(2,2) = zz - mz*mz;

      TMatrixDSymEigen me(m);
      auto& eigenval = me.GetEigenValues();
      auto& eigenvec = me.GetEigenVectors();

      double maxeval =  eigenval(0);
      int maxevalaxis = 0;
      // check
      //cout<<"axis 0 eigenvalue: "<<maxeval<<endl;

      for(int k=1; k<3; k++)
        {
	  if(eigenval(k)>maxeval)
            {
	      maxevalaxis = k;
	      maxeval = eigenval(k);
            }
	  //cout<<"axis "<<i<<" eigenvalue: "<<eigenval(i)<<endl;
        }
        
      // PCA main direction
      double pcadir_x = eigenvec(0, maxevalaxis);
      double pcadir_y = eigenvec(1, maxevalaxis);
      double pcadir_z = eigenvec(2, maxevalaxis);

      // track direction
      TVector3 pca_dir(pcadir_x, pcadir_y, pcadir_z);
      TVector3 track_dir;
      track_dir = (1./pca_dir.Mag())*pca_dir;
        
      // theta_y: 0-90 degree, upgoing
      if(track_dir.Y()<0) track_dir *= -1.0;
        
      // check
      //if( (int)(track_dir.Mag()) != 1 ) cout<<"PCA Direction vector not normalized! [units?]"<<endl;

      // Loop over point cloud to find the edges (start, end) along PCA main direction
        
      TVector3 start(xpt[i].at(0), ypt[i].at(0), zpt[i].at(0));
      TVector3 end(start);
      double start_proj = track_dir.Dot(start);
      double end_proj = track_dir.Dot(end);


      for(int j=1; j<npts; j++)
        {//each point
	  TVector3 point(xpt[i].at(j), ypt[i].at(j), zpt[i].at(j));
	  double point_proj = track_dir.Dot(point);
	  if(point_proj < start_proj)
            {
	      start=point;
	      start_proj = point_proj;
            }
	  if(point_proj > end_proj)
            {
	      end=point;
	      end_proj = point_proj;
            } 
        }//each point
       
      // the start, end points are more precise than the PCA direction in case of a big blob (isochronous track)
      TVector3 recon_dir = end - start;
      if(recon_dir.Y()<0) recon_dir *= -1.0;
      double costheta_y = recon_dir.Y()/recon_dir.Mag();
      double phi = TMath::ACos(recon_dir.Z()/TMath::Sqrt(recon_dir.X()*recon_dir.X()+recon_dir.Z()*recon_dir.Z()));
        
      TVector3 recon_dir_u = Muplane*recon_dir;
      if(recon_dir_u.Y()<0) recon_dir_u *= -1.0;
      double costheta_y_u = recon_dir_u.Y()/recon_dir_u.Mag();
      double phi_u = TMath::ACos(recon_dir_u.Z()/TMath::Sqrt(recon_dir_u.X()*recon_dir_u.X()+recon_dir_u.Z()*recon_dir_u.Z()));
       
      TVector3 recon_dir_v = Mvplane*recon_dir;
      if(recon_dir_v.Y()<0) recon_dir_v *= -1.0;
      double costheta_y_v = recon_dir_v.Y()/recon_dir_v.Mag();
      double phi_v = TMath::ACos(recon_dir_v.Z()/TMath::Sqrt(recon_dir_v.X()*recon_dir_v.X()+recon_dir_v.Z()*recon_dir_v.Z()));
        
      TVector3 center = 0.5*(start+end);
      double length = recon_dir.Mag();
     
      // PCA track info (not used)
      TVector3 center2(mx, my, mz);
      TVector3 start2 = center2-0.5*length*track_dir;
      TVector3 end2 = center2+0.5*length*track_dir;

      // debug
      if(costheta_y >= 1 || costheta_y <=0 || phi >= TMath::Pi() || phi <= 0) {
	recon_dir.Print();
	cout<<"length: "<<length<<endl;
      }
        
      // check
      if(i==cluster_check || cluster_check==-1){
        //cout<<"Direction: "<<track_dir.X()<<" "<<track_dir.Y()<<" "<<track_dir.Z()<<endl;
        /* cout<<"Cluster_id:"<<i<<endl; */
        /* cout<<"Length: "<<length<<endl; */
        /* cout<<"Direction: "<<costheta_y<<" "<<phi<<endl; */
        /* cout<<"Center1: "<<center.X()<<" "<<center.Y()<<" "<<center.Z()<<endl; */
        /* cout<<"Center2: "<<mx<<" "<<my<<" "<<mz<<endl; */
        /* cout<<"Start1: "<<start.X()<<" "<<start.Y()<<" "<<start.Z()<<endl; */
        /* cout<<"End1: "<<end.X()<<" "<<end.Y()<<" "<<end.Z()<<endl; */
        /* cout<<"Start2: "<<start2.X()<<" "<<start2.Y()<<" "<<start2.Z()<<endl; */
        /* cout<<"End2: "<<end2.X()<<" "<<end2.Y()<<" "<<end2.Z()<<endl; */
        /* cout<<"Start Difference: "<<100*(start2.X()-start.X())/start2.X()<<" "<<100*(start2.Y()-start.Y())/start2.Y()<<" "<<100*(start2.Z()-start.Z())/start2.Z()<<endl; */
        /* cout<<"End Difference: "<<100*(end2.X()-end.X())/end2.X()<<" "<<100*(end2.Y()-end.Y())/end2.Y()<<" "<<100*(end2.Z()-end.Z())/end2.Z()<<endl; */

        leval1->SetPoint(0, start.X(), start.Y(), start.Z());
        leval1->SetPoint(1, end.X(), end.Y(), end.Z());
        leval2->SetPoint(0, start2.X(), start2.Y(), start2.Z());
        leval2->SetPoint(1, end2.X(), end2.Y(), end2.Z());


	
        // fill histogram
        if(length>100* 0.9){// 1m * 90% = 0.9m // PCA failure
	  hntrack->Fill(costheta_y, phi);
	  hltrack->Fill(costheta_y, phi, length/100.0);
	  hntrack_u->Fill(costheta_y_u, phi_u);
	  hltrack_u->Fill(costheta_y_u, phi_u, length/100.0);
	  hntrack_v->Fill(costheta_y_v, phi_v);
	  hltrack_v->Fill(costheta_y_v, phi_v, length/100.0);


        }
      }
    }// each cluster

  //average length
  /* for(int m=1; m<=hntrack->GetNbinsX(); m++) */
  /* { */
  /*     for(int n=1; n<=hntrack->GetNbinsY(); n++) */
  /*     { */
  /*         if(hntrack->GetBinContent(m,n)!=0) */
  /*         { */
  /*             hltrack->SetBinContent(m,n,hltrack->GetBinContent(m,n)/hntrack->GetBinContent(m,n)); */
  /*         } */
  /*         if(hntrack_u->GetBinContent(m,n)!=0) */
  /*         { */
  /*             hltrack_u->SetBinContent(m,n,hltrack_u->GetBinContent(m,n)/hntrack_u->GetBinContent(m,n)); */
  /*         } */
  /*         if(hntrack_v->GetBinContent(m,n)!=0) */
  /*         { */
  /*             hltrack_v->SetBinContent(m,n,hltrack_v->GetBinContent(m,n)/hntrack_v->GetBinContent(m,n)); */
  /*         } */
  /*     } */
  /* } */


  //check
  if(cluster_check!=-1){
    TCanvas* canv_check = new TCanvas("canv_check","",400,1600);
    hrec->Draw();
    leval1->Draw("same l");
    leval2->Draw("same l");
    leval1->SetLineColor(kRed);
    leval1->SetLineWidth(2);
    leval2->SetLineColor(kGreen);
    leval2->SetLineStyle(kDashed);
    leval2->SetLineWidth(2);
    canv_check->SaveAs("check.root");
  }

  TTree* tracks = (TTree*)f->Get("T_track");
  double x0,y0,z0;
  double x1,y1,z1;
  tracks->SetBranchAddress("x0",&x0);
  tracks->SetBranchAddress("y0",&y0);
  tracks->SetBranchAddress("z0",&z0);
  tracks->SetBranchAddress("x1",&x1);
  tracks->SetBranchAddress("y1",&y1);
  tracks->SetBranchAddress("z1",&z1);
   

  for(int nt=0; nt<tracks->GetEntries(); nt++)
    {
      tracks->GetEntry(nt);
      TVector3 start0(x0, y0, z0);
      TVector3 end1(x1, y1, z1);
      TVector3 dir = end1 - start0;
      if(dir.Y()<0) dir *= -1.0;
      double cosy = dir.Y()/dir.Mag();
      double phiz = TMath::ACos(dir.Z()/TMath::Sqrt(dir.X()*dir.X()+dir.Z()*dir.Z()));
      TVector3 dir_u = Muplane*dir;
      if(dir_u.Y()<0) dir_u *= -1.0;
      double cosy_u = dir_u.Y()/dir_u.Mag();
      double phiz_u = TMath::ACos(dir_u.Z()/TMath::Sqrt(dir_u.X()*dir_u.X()+dir_u.Z()*dir_u.Z()));
      TVector3 dir_v = Mvplane*dir;
      if(dir_v.Y()<0) dir_v *= -1.0;
      double cosy_v = dir_v.Y()/dir_v.Mag();
      double phiz_v = TMath::ACos(dir_v.Z()/TMath::Sqrt(dir_v.X()*dir_v.X()+dir_v.Z()*dir_v.Z()));


      hntrack_true->Fill(cosy, phiz);
      hltrack_true->Fill(cosy, phiz, dir.Mag()/100.0);
      hntrack_true_u->Fill(cosy_u, phiz_u);
      hltrack_true_u->Fill(cosy_u, phiz_u, dir_u.Mag()/100.0);
      hntrack_true_v->Fill(cosy_v, phiz_v);
      hltrack_true_v->Fill(cosy_v, phiz_v, dir_v.Mag()/100.0);
    }

  //average length
  /* for(int m=1; m<=hntrack->GetNbinsX(); m++) */
  /* { */
  /*     for(int n=1; n<=hntrack->GetNbinsY(); n++) */
  /*     { */
  /*         if(hntrack->GetBinContent(m,n)!=0) */
  /*         { */
  /*             hltrack->SetBinContent(m,n,hltrack->GetBinContent(m,n)/hntrack->GetBinContent(m,n)); */
  /*         } */
  /*         if(hntrack_u->GetBinContent(m,n)!=0) */
  /*         { */
  /*             hltrack_u->SetBinContent(m,n,hltrack_u->GetBinContent(m,n)/hntrack_u->GetBinContent(m,n)); */
  /*         } */
  /*         if(hntrack_v->GetBinContent(m,n)!=0) */
  /*         { */
  /*             hltrack_v->SetBinContent(m,n,hltrack_v->GetBinContent(m,n)/hntrack_v->GetBinContent(m,n)); */
  /*         } */
  /*     } */
  /* } */
    
  // Truth vs Recon
  TH2D* hntrack_comp = new TH2D("hntrack_comp","number of tracks",10,0,1,10,0,TMath::Pi());
  TH2D* hltrack_comp = new TH2D("hltrack_comp","length of tracks",10,0,1,10,0,TMath::Pi());
  TH2D* hntrack_comp_u = new TH2D("hntrack_comp_u","number of tracks",10,0,1,10,0,TMath::Pi());
  TH2D* hltrack_comp_u = new TH2D("hltrack_comp_u","length of tracks",10,0,1,10,0,TMath::Pi());
  TH2D* hntrack_comp_v = new TH2D("hntrack_comp_v","number of tracks",10,0,1,10,0,TMath::Pi());
  TH2D* hltrack_comp_v = new TH2D("hltrack_comp_v","length of tracks",10,0,1,10,0,TMath::Pi());
  for(int m=1; m<=hntrack_comp->GetNbinsX(); m++)
    {
      for(int n=1; n<=hntrack_comp->GetNbinsY(); n++)
        {
	  if(hntrack_true->GetBinContent(m,n)!=0){
	    hntrack_comp->SetBinContent(m,n, hntrack->GetBinContent(m,n)/hntrack_true->GetBinContent(m,n)-1.0 );
	    hltrack_comp->SetBinContent(m,n, hltrack->GetBinContent(m,n)/hltrack_true->GetBinContent(m,n)-1.0 );
	  }
	  if(hntrack_true_u->GetBinContent(m,n)!=0){
	    hntrack_comp_u->SetBinContent(m,n, hntrack_u->GetBinContent(m,n)/hntrack_true_u->GetBinContent(m,n)-1.0 );
	    hltrack_comp_u->SetBinContent(m,n, hltrack_u->GetBinContent(m,n)/hltrack_true_u->GetBinContent(m,n)-1.0 );
	  }
	  if(hntrack_true_v->GetBinContent(m,n)!=0){
	    hntrack_comp_v->SetBinContent(m,n, hntrack_v->GetBinContent(m,n)/hntrack_true_v->GetBinContent(m,n)-1.0 );
	    hltrack_comp_v->SetBinContent(m,n, hltrack_v->GetBinContent(m,n)/hltrack_true_v->GetBinContent(m,n)-1.0 );
	  }
        }
    }


  set_plot_style();
  TFile* output = new TFile(outputroot, "RECREATE");
  TCanvas* canv = new TCanvas("canv","",1200,1600);
  canv->Divide(3,2);
  canv->cd(1);
  hntrack->Draw("colz");
  canv->cd(2);
  hntrack_true->Draw("colz");
  canv->cd(3);
  hntrack_comp->Draw("colz");
  hntrack_comp->GetXaxis()->SetTitle("cos(theta_y)");
  hntrack_comp->GetYaxis()->SetTitle("phi_xz");
  canv->cd(4);
  hltrack->Draw("colz");
  canv->cd(5);
  hltrack_true->Draw("colz");
  canv->cd(6);
  hltrack_comp->Draw("colz");
  canv->Write();

  TCanvas* canv2 = new TCanvas("canv2","",1200,1600);
  canv2->Divide(3,2);
  canv2->cd(1);
  hntrack_comp_u->Draw("colz");
  canv2->cd(2);
  hntrack_comp_v->Draw("colz");
  canv2->cd(3);
  hntrack_comp->Draw("colz");
  hntrack_comp->GetXaxis()->SetTitle("cos(theta_y)");
  hntrack_comp->GetYaxis()->SetTitle("phi_xz");
  canv2->cd(4);
  hltrack_comp_u->Draw("colz");
  canv2->cd(5);
  hltrack_comp_v->Draw("colz");
  canv2->cd(6);
  hltrack_comp->Draw("colz");
  canv2->Write();

  hntrack->Write();
  hntrack_true->Write();
  hntrack_comp->Write();
  hntrack_u->Write();
  hntrack_true_u->Write();
  hntrack_comp_u->Write();
  hntrack_v->Write();
  hntrack_true_v->Write();
  hntrack_comp_v->Write();

  hltrack->Write();
  hltrack_true->Write();
  hltrack_comp->Write();
  hltrack_u->Write();
  hltrack_true_u->Write();
  hltrack_comp_u->Write();
  hltrack_v->Write();
  hltrack_true_v->Write();
  hltrack_comp_v->Write();


  return 0;
}
