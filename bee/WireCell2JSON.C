#include "WCReader.cc"

void WireCell2JSON(TString filename, TString option = "truth", TString outfile="wc.json")
{
  //gROOT->Reset();
  //gROOT->ProcessLine(".x loadClasses.C" );
    gErrorIgnoreLevel=2001;
    WCReader r(filename, outfile);
    if (option == "mc") {
        r.DumpMC();
    }
    else {
        r.DumpSpacePoints(option);
    }
}
