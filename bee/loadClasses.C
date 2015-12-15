{
    TString include = ".include ";
    TString load = ".L ";

    TString prefix = ".";
    gROOT->ProcessLine( include + prefix );
    gROOT->ProcessLine( load + prefix + "/WCReader.cc+" );

}