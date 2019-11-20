#include <vector>
#include <iostream>
#include "math.h"
#include "TH1F.h"
#include "stdlib.h"
#include "TF1.h"
#include "TMath.h"
using namespace std;

namespace WCPDune{

  //remove pedestal, copied from DatauBooNEFrameDataSource.cxx
int get_baseline_rms(TH1F* h1, double& retBaseline, double& retRms){
   retBaseline=-1;
   retRms=-1;
   if(h1->GetEntries()==0){
     return -1;
   }

   float mean = h1->GetSum()/h1->GetNbinsX();
   float rms = 0;
   int nbin = h1->GetNbinsX();
   for (int j=0;j!=nbin;j++){
     rms += pow(h1->GetBinContent(j+1)-mean,2);
   }
   rms = sqrt(rms/h1->GetNbinsX());
   
   TH1F h2("h2","h2",100,mean-10*rms,mean+10*rms);
   for (int j=0;j!=nbin;j++){
     if (fabs(h1->GetBinContent(j+1)-mean)<6*rms)
       h2.Fill(h1->GetBinContent(j+1));
   }
   
   double par[3];
   // double xq = 0.5;
   // h2.GetQuantiles(1,&par[1],&xq);
   // xq = 0.5 + 0.34;
   // h2.GetQuantiles(1,&par[0],&xq);
   // xq = 0.5 - 0.34;
   // h2.GetQuantiels(1,&par[2],&xq);
   
   // par[2] = sqrt((pow(par[0]-par[1],2)+pow(par[2]-par[1],2))/2.);
   // par[0] = h2.GetSum();
   // f1->SetParameters(par);
   // h2.Fit(f1,"Q0","");
   // f1->GetParameters(par);
   
   retRms = rms;
   // a different way to calculate the mean ... 
   double xq = 0.5;
   //cout << "entries: " << h2.GetEntries() << " integral: " << h2.Integral() << endl;
   if(h2.GetEntries()==0){
     return -1;
   }
   h2.GetQuantiles(1,&par[1],&xq);
   
   //for (int j=0;j!=nbin;j++){
   //  h1->SetBinContent(j+1,h1->GetBinContent(j+1)-par[1]);
   //}
   retBaseline = par[1];
   return 0;
}

int LedgeIdentify(int channel, TH1F* h2, double baseline, double & LedgeStart, double & LedgeEnd){
	int ledge = 0;
        int UNIT = 5;    // rebin unit
        int CONTIN = 20; // length of the continuous region
        int JUMPS = 4;   // how many bins can accidental jump
        vector<int> averaged; // store the rebinned waveform
        int up = h2->GetNbinsX()/UNIT;
	int nticks =  h2->GetNbinsX();
	// rebin 
	for(int i=0;i<up;i++){ 
                int temp = 0;
                for(int j=0;j<UNIT;j++){
                   temp += h2->GetBinContent(i*UNIT+1+j);
                }
                averaged.push_back(temp);
        }
	// refine the selection cuts if there is a large signal
	if(h2->GetMaximum()-baseline>1000) { CONTIN = 16; JUMPS = 5; }
	// start judging
	int decreaseD = 0, tolerence=0;
        int start = 0, end = 0;
        for(int i=1;i<up-1;i++){
                if(averaged.at(i)<averaged.at(i-1)) {
                        if(decreaseD==0) start = i;
                        decreaseD +=1;
                }
		 else {
                        if(averaged.at(i+1)<averaged.at(i-1)&&tolerence<JUMPS&&decreaseD>0){ // we can ignore several jumps in the decreasing region
                                decreaseD+=2;
                                tolerence++;
                                i = i+1;
                        }
                        else{
                                if(decreaseD>CONTIN){
                                        ledge = 1;
                                        LedgeStart = start*UNIT;
                                        break;
                                }
                                else{
                                        decreaseD = 0;
                                        tolerence=0;
                                        start = 0;
                                        end = 0;
                                }
                        }
                }
        }
	// find the sharp start edge
	 if(ledge == 1&&LedgeStart>30){ 
                int edge = 0;
                int i = LedgeStart/UNIT-1;
                if(averaged.at(i)>averaged.at(i-1)&&averaged.at(i-1)>averaged.at(i-2)){ // find a edge
                        edge = 1;
                }
                if(edge == 0) ledge = 0; // if no edge, this is not ledge
                if((averaged.at(i)-averaged.at(i-2)<10*UNIT)&&(averaged.at(i)-averaged.at(i-3)<10*UNIT)) // slope cut
                        ledge = 0;
                if(averaged.at(LedgeStart/UNIT)-baseline*UNIT>200*UNIT) ledge = 0; // ledge is close to the baseline
        }
	// test the decay time
	if(ledge == 1&&LedgeStart>20){
                double height = 0;
                if(LedgeStart<5750) { // calculate the height of edge
                        double tempHeight = h2 ->GetBinContent(LedgeStart+1+200) +  h2 ->GetBinContent(LedgeStart+1+220) +  h2 ->GetBinContent(LedgeStart+1+180) +  h2 ->GetBinContent(LedgeStart+1+240);
                        height = h2 ->GetBinContent(LedgeStart+1) - tempHeight/4;
			height /= 0.7;
                }
                else height =  h2 ->GetBinContent(LedgeStart+1) -  h2 ->GetBinContent(6000);
                if(height<0) height = 80; // norminal value
                if(height>30&&LedgeStart<5900){ // test the decay with a relatively large height
                        double height50 = 0, height100 = 0;
                        height50 =  h2 ->GetBinContent(LedgeStart+51);
                        height100 =  h2 ->GetBinContent(LedgeStart+101);
                        double height50Pre =   h2 ->GetBinContent(LedgeStart+1)- height*(1-exp(-50/100.)); // minimum 100 ticks decay time
                        double height100Pre =   h2 ->GetBinContent(LedgeStart+1) - height*(1-exp(-100./100)); // minimum 100 ticks decay time
			// if the measured is much smaller than the predicted, this is not ledge
                        if(height50-height50Pre<-8) ledge = 0; 
                        if(height100-height100Pre<-8)  ledge = 0;
                }
        }

	// determine the end of ledge
	// case 1: find a jump of 10 ADC in the rebinned waveform
	// case 2: a continuous 20 ticks has an average close to baseline, and a RMS larger than 3 ADC
	// case 3: reaching the tick 6000 but neither 1 nor 2 occurs
	if(ledge ==1){
		LedgeEnd = 0;
		for(int i = LedgeStart/UNIT; i<up-1; i++){ // case 1
			if(averaged.at(i+1)-averaged.at(i)>50) { 
				LedgeEnd = i*UNIT+5;
				break;
			}
		}
		if(LedgeEnd == 0) { // not find a jump, case 2
			double tempA[20];
			for(int i = LedgeStart+80;i<nticks-20;i+=20){
				for(int j=i;j<20+i;j++)
					tempA[j-i] = h2->GetBinContent(j+1);
				if(TMath::Mean(20,tempA)-baseline<2&&TMath::RMS(20,tempA)>3){
					LedgeEnd = i;
					break;
				}
			}
		}
		if(LedgeEnd == 0) LedgeEnd = 6000;
	}
	// done, release the memory
        vector<int>(averaged).swap(averaged);
        return ledge;

}
int judgePlateau(int channel, TH1F* h2,double baseline, double & PlateauStart, double & PlateauStartEnd){
        int continueN = 0;
        int threshold = 200;
        int maximumF  = 50;
        int maxBin = h2->GetMaximumBin();
        for(int i=maxBin+10;i<5880&&i<maxBin+500;i++){
                int plateau = 1;
                int max = 0, min = 10000;
                for(int j=i;j<i+20;j++){
                        int binC = h2->GetBinContent(j+1);
                        if(binC<baseline+threshold||binC>h2->GetMaximum()-500) {
                                plateau = 0;
                                break;
                        }
                        if(binC>max) max = binC;
                        if(binC<min) min = binC;
                }
                if(plateau==1&&max-min<maximumF){ // plateau found
                        PlateauStart = i;
                        PlateauStartEnd = i+20;
                        for(int k = i+20; k<6000;k++){
                                if( h2->GetBinContent(k+1)<baseline+threshold){
                                        PlateauStartEnd = k-1;
                                        break;
                                }
                        }
                        return 1;
                }
        }
        return 0;
}

}//WCPDune
