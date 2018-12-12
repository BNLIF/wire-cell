// Created by Carlos (sarastce@mail.uc.edu), 12/12/2018
#include "TH1F.h"

namespace WireCellDune{


int StickyCodeIdent(TH1F* hsig){
// returns the number of sticky codes per channel

    int count(0), count1(0), count63(0); // mod64 = 0, 1, 63
    int nticks = hsig->GetNbinsX();
    for(int j=0; j<nticks; j++){
         int x = hsig->GetBinContent(j+1);
         x %= 64;
         if(x==0){
           count++;
         }
         else if(x==1){
           count1++;
         }
         else if(x==63){
           count63++;
         }
    }
    return count + count1 + count63;
}

}
