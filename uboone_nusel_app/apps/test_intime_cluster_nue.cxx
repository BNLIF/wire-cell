#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TVector3.h"

#include <iostream>
#include <vector>

using namespace std;

bool boundary(TVector3 v)
{
    bool isBound1 = true;
    bool isBound2 = true;
    Double_t x = v.X();
    Double_t y = v.Y();
    Double_t z = v.Z();
    // X-Y plane: bottom - (80, -117) --- (256, -98), top - (100, 118) --- (256, 103)
    // 6 distances (sign) to boundary
    Double_t d1 = x;
    Double_t d2 = y-118;
    Double_t d3 = y+117;
    Double_t d4 = x-256;
    Double_t d5 = ((103-118)*x-(256-100)*y+256*118-103*100)/TMath::Sqrt((103-118)*(103-118)+(256-100)*(256-100));
    Double_t d6 = ((-98+117)*x-(256-80)*y+256*(-117)-(-98)*80)/TMath::Sqrt((-98+117)*(-98+117)+(256-80)*(256-80));
   
    //cout<<d1<<" "<<d2<<" "<<d3<<" "<<d4<<" "<<d5<<" "<<d6<<endl;
    Double_t fv = 5; // FC cut, cm
    if(d1>fv && d2<-fv && d3>fv && d4<-fv && d5>fv && d6<-fv) isBound1 = false; // contained 

    // X-Z plane: left - (0, 120) --- (11, 256), right - (1026, 256) --- (1037, 120)
    Double_t dd1 = z;
    Double_t dd2 = x-256;
    Double_t dd3 = x;
    Double_t dd4 = z - 1037;
    Double_t dd5 = ((120-256)*x-(1037-1026)*y+1037*256-120*1026)/TMath::Sqrt((120-256)*(120-256)+(1037-1026)*(1037-1026));
    Double_t dd6 = ((256-120)*x-(11-0)*y+11*120-256*0)/TMath::Sqrt((256-120)*(256-120)+(11-0)*(11-0));
    //cout<<dd1<<" "<<dd2<<" "<<dd3<<" "<<dd4<<" "<<dd5<<" "<<dd6<<endl;
    if(dd1>fv && dd2<-fv && dd3>fv && dd4<-fv && dd5>fv && dd6>fv) isBound2 = false; // contained 

    //cout<<"Bound1: "<<isBound1<<endl;
    //cout<<"Bound2: "<<isBound2<<endl;
    return isBound1||isBound2;
}

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
    if(argc < 2){
        std::cerr << "usage: apps /path/to/match.root" << std::endl;
        return 1;
    }
    TH1::AddDirectory(kFALSE);

    set_plot_style();

    TString fileinput = argv[1];
    TFile *file = new TFile(fileinput, "READ");

    TTree *proj = (TTree*)file->Get("T_proj");
    std::vector<int> *cluster_id = new std::vector<int>;
    std::vector<std::vector<int>> *channel_vec = new std::vector<std::vector<int>>;
    std::vector<std::vector<int>> *time_slice_vec = new std::vector<std::vector<int>>;
    std::vector<std::vector<int>> *charge_vec = new std::vector<std::vector<int>>;
    proj->SetBranchAddress("cluster_id", &cluster_id);
    proj->SetBranchAddress("channel", &channel_vec);
    proj->SetBranchAddress("time_slice", &time_slice_vec);
    proj->SetBranchAddress("charge", &charge_vec);

    // Trun info
    TTree* Trun = (TTree*)file->Get("Trun");
    int eventNo;
    int runNo;
    int subRunNo;
    unsigned int triggerBits;
    Trun->SetBranchAddress("eventNo", &eventNo);
    Trun->SetBranchAddress("runNo", &runNo);
    Trun->SetBranchAddress("subRunNo", &subRunNo);
    Trun->SetBranchAddress("triggerBits", &triggerBits);
    Trun->GetEntry(0);

    double lowerwindow = 0;
    double upperwindow = 0;
    if(triggerBits==2048) { lowerwindow = 3; upperwindow = 5; }// bnb
    if(triggerBits==512) { lowerwindow = 3.45; upperwindow = 5.45; } // extbnb

    // flash id and time
    TTree *flash = (TTree*)file->Get("T_flash");
    int flashes_id;
    double flash_time;
    flash->SetBranchAddress("flash_id", &flashes_id);
    flash->SetBranchAddress("time", &flash_time);
    //flash id - time
    std::map<int, double> flash_time_map;
    for(int i=0; i<flash->GetEntries(); i++)
    {
        flash->GetEntry(i);
        if(flash_time<upperwindow && flash_time>lowerwindow) flash_time_map[flashes_id] = flash_time;
    }


    // tpc-flash matching
    TTree *match = (TTree*)file->Get("T_match");
    int tpc_cluster_id;
    int flash_id;
    match->SetBranchAddress("tpc_cluster_id", &tpc_cluster_id);
    match->SetBranchAddress("flash_id", &flash_id);
    
    
    // filtering in-time flash
    std::map<int, int> intime_cluster_id; // cluster_id, flash_id

    // flash time filter
    for(int i=0; i<match->GetEntries(); i++)
    {
        match->GetEntry(i);
        if(flash_id != -1){
            auto it = flash_time_map.find(flash_id);
            if(it != flash_time_map.end()) // within (2, 6) us
            {
                intime_cluster_id[tpc_cluster_id] = flash_id;
            }
        }
    }


    // matched cluster id
    for(auto &it : intime_cluster_id)
    {
        //cout<<"Cluster: "<<it.first<<endl;
    }


    // clusters
    TTree *clusters = (TTree*)file->Get("T_cluster");
    TTree *intime_cluster;
    bool through = true; //any part of cluster touch the boundary within 5 cm
    bool thiscluster = false;
    bool matched = false;

    // read in
    TH2I* hucharge = new TH2I("hucharge","",2400, 0, 2400, 2400, 0, 2400);
    TH2I* hvcharge = new TH2I("hvcharge","",2400, 2400, 4800, 2400, 0, 2400);
    TH2I* hwcharge = new TH2I("hwcharge","",3456, 4800, 8256, 2400, 0, 2400);

    // charge scale
    TH1F* huscale = new TH1F("huscale","",1e4,0,1e5);
    TH1F* hvscale = new TH1F("hvscale","",1e4,0,1e5);
    TH1F* hwscale = new TH1F("hwscale","",1e4,0,1e5);
 
    // channel range
    int uchannel_min = 10000;
    int uchannel_max = 0;
    int vchannel_min = 10000;
    int vchannel_max = 0;
    int wchannel_min = 10000;
    int wchannel_max = 0;
    
    int uchannel0_min = 10000;
    int uchannel0_max = 0;
    int vchannel0_min = 10000;
    int vchannel0_max = 0;
    int wchannel0_min = 10000;
    int wchannel0_max = 0;

    // time range
    int utime_min = 3000;
    int utime_max = 0;
    int vtime_min = 3000;
    int vtime_max = 0;
    int wtime_min = 3000;
    int wtime_max = 0;
    
    int utime0_min = 3000;
    int utime0_max = 0;
    int vtime0_min = 3000;
    int vtime0_max = 0;
    int wtime0_min = 3000;
    int wtime0_max = 0;

    // 3D ranges
    double xmin=300;
    double xmax=0;
    double ymin=150;
    double ymax=-150;
    double zmin=1200;
    double zmax=0;

    double xmin0=300;
    double xmax0=0;
    double ymin0=150;
    double ymax0=-150;
    double zmin0=1200;
    double zmax0=0;
  
    //cout<<"READ: "<<proj->GetEntries()<<endl; // should be 1
    for(int entry=0; entry<proj->GetEntries(); entry++)
    {
        proj->GetEntry(entry);

        // # of charge hits
        //std::cout<<"Cluster: "<<cluster_id->size()<<" Channel: "<<channel_vec->size()<<" Time slices: "<<time_slice_vec->size()<<" Charge: "<<charge_vec->size()<<std::endl;
        
        for(int i=0; i<cluster_id->size(); i++)
        {
            //cout<<i<<endl;
            if(intime_cluster_id.find(cluster_id->at(i)) != intime_cluster_id.end())
            {
            std::vector<int> charge = charge_vec->at(i);
            std::vector<int> time = time_slice_vec->at(i);
            std::vector<int> channel = channel_vec->at(i);          

            cout<<"Cluster: "<<cluster_id->at(i)<<" Time: "<<flash_time_map[intime_cluster_id[cluster_id->at(i)]]<<endl;
            //cout<<"# of Channels: "<<channel.size()<<endl;
            //cout<<"# of time slices: "<<time.size()<<endl;
            //cout<<"# of charges : "<<charge.size()<<endl;

            if(charge.size()>1000){
		    matched = true;
            thiscluster=false;
 
            // on boundary
            TString String_clusterid;
            String_clusterid.Form("cluster_id==%d", cluster_id->at(i)); 
            gROOT->cd();
            intime_cluster = clusters->CopyTree(String_clusterid.Data());
            //cout<<intime_cluster->GetEntries()<<endl;
            double cx;
            double cy;
            double cz;
            intime_cluster->SetBranchAddress("x", &cx);
            intime_cluster->SetBranchAddress("y", &cy);
            intime_cluster->SetBranchAddress("z", &cz);

            for(int i=0; i<intime_cluster->GetEntries(); i++)
            {
                intime_cluster->GetEntry(i);
                TVector3 v0(cx, cy, cz);
                //if(boundary(v0)){
                if(0){
                    thiscluster=true;
                    cout<<"This cluter is on boundary!"<<endl;
                    break;
                }
                else{
                    xmin = xmin<cx?xmin:cx;
                    xmax = xmax>cx?xmax:cx;
                    ymin = ymin<cy?ymin:cy;
                    ymax = ymax>cy?ymax:cy;
                    zmin = zmin<cz?zmin:cz;
                    zmax = zmax>cz?zmax:cz;
                }
            }
            // update cluster range if not on boundary
            if(!thiscluster){
                xmin0 = xmin0<xmin?xmin0:xmin;
                xmax0 = xmax0>xmax?xmax0:xmax;
                ymin0 = ymin0<ymin?ymin0:ymin;
                ymax0 = ymax0>ymax?ymax0:ymax;
                zmin0 = zmin0<zmin?zmin0:zmin;
                zmax0 = zmax0>zmax?zmax0:zmax;

                for(int p=0; p<charge.size(); p++)
                {
                    if(channel.at(p)<2400) {
                        //cout<<"U plane: "<<channel.at(p)<<endl;
                        hucharge->SetBinContent(channel.at(p)+1, time.at(p), charge.at(p)+hucharge->GetBinContent(channel.at(p)+1, time.at(p)));
                        uchannel_min = uchannel_min<channel.at(p)?uchannel_min:channel.at(p);
                        uchannel_max = uchannel_max>channel.at(p)?uchannel_max:channel.at(p);
                        utime_min = utime_min<time.at(p)?utime_min:time.at(p);
                        utime_max = utime_max>time.at(p)?utime_max:time.at(p);
                    }
                    if(channel.at(p)>=2400 && channel.at(p)<4800) {
                        //cout<<"V plane: "<<channel.at(p)<<endl;
                        hvcharge->SetBinContent(channel.at(p)-2399, time.at(p), charge.at(p)+hvcharge->GetBinContent(channel.at(p)-2399, time.at(p)));
                        vchannel_min = vchannel_min<channel.at(p)?vchannel_min:channel.at(p);
                        vchannel_max = vchannel_max>channel.at(p)?vchannel_max:channel.at(p);
                        vtime_min = vtime_min<time.at(p)?vtime_min:time.at(p);
                        vtime_max = vtime_max>time.at(p)?vtime_max:time.at(p);
                    }
                    if(channel.at(p)>=4800) {
                        //cout<<"W plane: "<<channel.at(p)<<endl;
                        hwcharge->SetBinContent(channel.at(p)-4799, time.at(p), charge.at(p)+hwcharge->GetBinContent(channel.at(p)-4799, time.at(p)));
                        wchannel_min = wchannel_min<channel.at(p)?wchannel_min:channel.at(p);
                        wchannel_max = wchannel_max>channel.at(p)?wchannel_max:channel.at(p);
                        wtime_min = wtime_min<time.at(p)?wtime_min:time.at(p);
                        wtime_max = wtime_max>time.at(p)?wtime_max:time.at(p);
                    }
                }
                uchannel0_min = uchannel0_min<uchannel_min?uchannel0_min:uchannel_min;
                vchannel0_min = vchannel0_min<vchannel_min?vchannel0_min:vchannel_min;
                wchannel0_min = wchannel0_min<wchannel_min?wchannel0_min:wchannel_min;
                utime0_min = utime0_min<utime_min?utime0_min:utime_min;
                vtime0_min = vtime0_min<vtime_min?vtime0_min:vtime_min;
                wtime0_min = wtime0_min<wtime_min?wtime0_min:wtime_min; 
                uchannel0_max = uchannel0_max>uchannel_max?uchannel0_max:uchannel_max;
                vchannel0_max = vchannel0_max>vchannel_max?vchannel0_max:vchannel_max;
                wchannel0_max = wchannel0_max>wchannel_max?wchannel0_max:wchannel_max;
                utime0_max = utime0_max>utime_max?utime0_max:utime_max;
                vtime0_max = vtime0_max>vtime_max?vtime0_max:vtime_max;
                wtime0_max = wtime0_max>wtime_max?wtime0_max:wtime_max; 
            }

            // reset individual cluster range
            uchannel_min = 10000;
            uchannel_max = 0;
            vchannel_min = 10000;
            vchannel_max = 0;
            wchannel_min = 10000;
            wchannel_max = 0;
            utime_min = 3000;
            utime_max = 0;
            vtime_min = 3000;
            vtime_max = 0;
            wtime_min = 3000;
            wtime_max = 0;

            xmin = 300;
            xmax = 0;
            ymin = 150;
            ymax = -150;
            zmin = 1200;
            zmax = 0;
            
            through = through&&thiscluster; 
            } // cluster size
            else{
                cout<<"Extremely low visible energy event!"<<endl;
            }
            } // intime cluster
        }
    }
    if(through && matched){ 
        cout<<"All clusters at boundary!"<<endl;
    }

    if(!through && matched){
        cout<<"Run-subRun-eventNo: "<<runNo<<" "<<subRunNo<<" "<<eventNo<<endl;
        cout<<"One candidate: "<<endl;
        cout<<"U channel range: "<<uchannel0_min<<" - "<<uchannel0_max<<endl;
        cout<<"V channel range: "<<vchannel0_min<<" - "<<vchannel0_max<<endl;
        cout<<"W channel range: "<<wchannel0_min<<" - "<<wchannel0_max<<endl;
        cout<<"U time tick (2 us) range: "<<utime0_min<<" - "<<utime0_max<<endl;
        cout<<"V time tick (2 us) range: "<<vtime0_min<<" - "<<vtime0_max<<endl;
        cout<<"W time tick (2 us) range: "<<wtime0_min<<" - "<<wtime0_max<<endl;
        cout<<"X range: "<<xmin0<<" - "<<xmax0<<endl;
        cout<<"Y range: "<<ymin0<<" - "<<ymax0<<endl;
        cout<<"Z range: "<<zmin0<<" - "<<zmax0<<endl;

        // 2D charge plots
        Double_t x[2]={0.30, 0.99};
        Double_t y[2];

        for(int i=uchannel0_min+1; i<=uchannel0_max+1; i++)
        {
            for(int j=utime0_min; j<=utime0_max; j++)
            {
                int uq = hucharge->GetBinContent(i, j);
                if(uq>0) huscale->Fill(uq);
            }
        }
        for(int i=vchannel0_min+1; i<=vchannel0_max+1; i++)
        {
            for(int j=vtime0_min; j<=vtime0_max; j++)
            {
                int vq = hvcharge->GetBinContent(i-2399, j);
                if(vq>0) hvscale->Fill(vq);
            }
        }
        for(int i=wchannel0_min+1; i<=wchannel0_max+1; i++)
        {
            for(int j=wtime0_min; j<=wtime0_max; j++)
            {
                int wq = hwcharge->GetBinContent(i-4799, j);
                if(wq>0) hwscale->Fill(wq);
            }
        }

        TString fileoutput;
        fileoutput.Form("intime_chargeNue_%d_%d_%d", runNo, subRunNo, eventNo);
        //TFile* ff = new TFile(fileoutput+".root","RECREATE");
        TCanvas *c = new TCanvas("c","",1800, 500);
        c->UseCurrentStyle();
        c->Divide(3,1);

        c->cd(1);
        huscale->GetQuantiles(2, y, x);
        //cout<<"U charge range: "<<y[0]<<" "<<y[1]<<endl;
        //Zeroout(hucharge);
        hucharge->GetZaxis()->SetRangeUser(y[0], y[1]);
        hucharge->GetXaxis()->SetRangeUser((uchannel0_min-10)>0?(uchannel0_min-10):0, (uchannel0_max+10)<2400?(uchannel0_max+10):2400);
        hucharge->GetYaxis()->SetRangeUser((utime0_min-15)>0?(utime0_min-15):0, (utime0_max+15)<2400?(utime0_max+15):2400);
        hucharge->Draw("colz");
        hucharge->GetXaxis()->SetTitle("Wire");
        hucharge->GetYaxis()->SetTitle("Time tick [2 #mus]");
        hucharge->GetZaxis()->SetTitle("Charge [e-]");

        c->cd(2);
        hvscale->GetQuantiles(2, y, x);
        //cout<<"V charge range: "<<y[0]<<" "<<y[1]<<endl;
        //Zeroout(hvcharge);
        hvcharge->GetZaxis()->SetRangeUser(y[0], y[1]);
        hvcharge->GetXaxis()->SetRangeUser((vchannel0_min-10)>2400?(vchannel0_min-10):2400, (vchannel0_max+10)<4800?(vchannel0_max+10):4800);
        hvcharge->GetYaxis()->SetRangeUser((vtime0_min-15)>0?(vtime0_min-15):0, (vtime0_max+15)<2400?(vtime0_max+15):2400);
        hvcharge->Draw("colz");
        hvcharge->GetXaxis()->SetTitle("Wire");
        hvcharge->GetYaxis()->SetTitle("Time tick [2 #mus]");
        hvcharge->GetZaxis()->SetTitle("Charge [e-]");

        c->cd(3);
        hwscale->GetQuantiles(2, y, x);
        //cout<<"W charge range: "<<y[0]<<" "<<y[1]<<endl;
        //Zeroout(hwcharge);
        hwcharge->GetZaxis()->SetRangeUser(y[0], y[1]);
        hwcharge->GetXaxis()->SetRangeUser((wchannel0_min-10)>4800?(wchannel0_min-10):4800, (wchannel0_max+10)<8256?(wchannel0_max+10):8256);
        hwcharge->GetYaxis()->SetRangeUser((wtime0_min-15)>0?(wtime0_min-15):0, (wtime0_max+15)<2400?(wtime0_max+15):2400);
        hwcharge->Draw("colz");
        hwcharge->GetXaxis()->SetTitle("Wire");
        hwcharge->GetYaxis()->SetTitle("Time tick [2 #mus]");
        hwcharge->GetZaxis()->SetTitle("Charge [e-]");


        c->SaveAs(fileoutput+".png");

    }

    file->Close();
    
    return 0;
}
