#ifndef WCReader_H
#define WCReader_H

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include "TString.h"

class TTree;
class TFile;
class TDatabasePDG;

class WCReader {
public:
    TFile *rootFile;
    ofstream jsonFile;

    enum LIMITS {
        MAX_TRACKS = 30000,
    };

    int mc_Ntrack;  // number of tracks in MC
    int mc_id[MAX_TRACKS];  // track id; size == mc_Ntrack
    int mc_pdg[MAX_TRACKS];  // track particle pdg; size == mc_Ntrack
    int mc_process[MAX_TRACKS];  // track generation process code; size == mc_Ntrack
    int mc_mother[MAX_TRACKS];  // mother id of this track; size == mc_Ntrack
    float mc_startXYZT[MAX_TRACKS][4];  // start position of this track; size == mc_Ntrack
    float mc_endXYZT[MAX_TRACKS][4];  // end position of this track; size == mc_Ntrack
    float mc_startMomentum[MAX_TRACKS][4];  // start momentum of this track; size == mc_Ntrack
    float mc_endMomentum[MAX_TRACKS][4];  // end momentum of this track; size == mc_Ntrack
    std::vector<std::vector<int> > *mc_daughters;  // daughters id of this track; vector

    // int mc_isnu; // is neutrino interaction
    // int mc_nGeniePrimaries; // number of Genie primaries
    // int mc_nu_pdg; // pdg code of neutrino
    // int mc_nu_ccnc; // cc or nc
    // int mc_nu_mode; // mode: http://nusoft.fnal.gov/larsoft/doxsvn/html/MCNeutrino_8h_source.html
    // int mc_nu_intType; // interaction type
    // int mc_nu_target; // target interaction
    // int mc_hitnuc; // hit nucleon
    // int mc_hitquark; // hit quark
    // double mc_nu_Q2; // Q^2
    // double mc_nu_W; // W
    // double mc_nu_X; // X
    // double mc_nu_Y; // Y
    // double mc_nu_Pt; // Pt
    // double mc_nu_Theta; // angle relative to lepton
    // float mc_nu_pos[4];  // interaction position of nu
    // float mc_nu_mom[4];  // interaction momentum of nu

    // ----- derived ---
    std::map<int, int> trackIndex;
    std::vector<std::vector<int> > trackParents;
    std::vector<std::vector<int> > trackChildren;
    std::vector<std::vector<int> > trackSiblings;
    TDatabasePDG *dbPDG;

    //-------------------------------------
    WCReader();
    WCReader(const char* filename, const char* jsonFileName="wc.json");
    virtual ~WCReader();

    void DumpRunInfo();
    void DumpSpacePoints(TString option);
    void DumpMC();
    void DumpDeadArea();
    void DumpOp();

    void print_vector(ostream& out, vector<double>& v, TString desc, bool end=false);
    void print_vector_vector(ostream& out, vector<vector<double> >& vv, TString desc, bool end=false);

    bool IsPrimary(int i) { return mc_mother[i] == 0 ; }
    void ProcessTracks();
    double KE(float* momentum);  // KE
    TString PDGName(int pdg);
    double thresh_KE;
    bool KeepMC(int i);

    bool DumpMCJSON(int id, ostream& out);
    void DumpMCJSON(ostream& out = std::cout);


};

#endif