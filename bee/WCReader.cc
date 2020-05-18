#include "WCReader.h"

#include <iostream>
#include <vector>
#include <map>
#include <iomanip>

#include "TNamed.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TClonesArray.h"
#include "TH1F.h"

using namespace std;

WCReader::WCReader(){}

//----------------------------------------------------------------
WCReader::WCReader(const char* filename, const char* jsonFileName)
{
    dbPDG = new TDatabasePDG();
    mc_daughters = new std::vector<std::vector<int> >;  // daughters id of this track; vector
    thresh_KE = 1; // MeV
 
    rootFile = new TFile(filename);

    jsonFile.open(jsonFileName);
}

//----------------------------------------------------------------
WCReader::~WCReader()
{
    jsonFile.close();

    delete mc_daughters;

    rootFile->Close();
    delete rootFile;
}

//----------------------------------------------------------------
void WCReader::DumpRunInfo()
{
    int runNo=0, subRunNo=0, eventNo=0, detector=0;

    TTree *tr = (TTree*)rootFile->Get("Trun");
    if (tr) {
        tr->SetBranchAddress("runNo", &runNo);
        tr->SetBranchAddress("subRunNo", &subRunNo);
        tr->SetBranchAddress("eventNo", &eventNo);
        tr->SetBranchAddress("detector", &detector);
        tr->GetEntry(0);
    }

    jsonFile << fixed << setprecision(0);

    jsonFile << '"' << "runNo" << '"' << ":" << '"' << runNo << '"' << "," << endl;
    jsonFile << '"' << "subRunNo" << '"' << ":" << '"' << subRunNo << '"' << "," << endl;
    jsonFile << '"' << "eventNo" << '"' << ":" << '"' << eventNo << '"' << "," << endl;
    if (detector == 0) {
        jsonFile << '"' << "geom" << '"' << ":" << '"' << "uboone" << '"' << endl;
    }
    else if (detector == 1) {
        jsonFile << '"' << "geom" << '"' << ":" << '"' << "dune35t" << '"' << endl;
    }
    else if (detector == 2) {
        jsonFile << '"' << "geom" << '"' << ":" << '"' << "protodune" << '"' << endl;
    }
    else if (detector == 3) {
        jsonFile << '"' << "geom" << '"' << ":" << '"' << "dune10kt_workspace" << '"' << endl;
    }
}

void WCReader::DumpDeadArea()
{
    int nPoints = 0;
    const int MAX_POINTS = 10;
    double y[MAX_POINTS], z[MAX_POINTS];

    TTree *t = (TTree*)rootFile->Get("T_bad");
    if (t) {
        t->SetBranchAddress("bad_npoints", &nPoints);
        t->SetBranchAddress("bad_y", &y);
        t->SetBranchAddress("bad_z", &z);
    }

    jsonFile << "[";
    jsonFile << fixed << setprecision(1);
    int nEntries = t->GetEntries();
    for (int i=0; i<nEntries; i++) {
        t->GetEntry(i);
        if (i>0) jsonFile << "," << endl;
        jsonFile << "[";
        for (int j=0; j<nPoints; j++) {
            if (j>0) jsonFile << ",";
            jsonFile << "[" << y[j] << ", " << z[j] << "]";
        }
        jsonFile << "]";
    }
    jsonFile << "]" << endl;
}

void WCReader::DumpOp()
{
    double PE[32];
    double time;
    double total_PE;
    vector<double> *l1_fired_time = new std::vector<double>;
    vector<double> *l1_fired_pe = new std::vector<double>;

    TTree *t = (TTree*)rootFile->Get("T_flash");
    if (t) {
        t->SetBranchAddress("PE", &PE);
        t->SetBranchAddress("time", &time);
        t->SetBranchAddress("total_PE", &total_PE);
        t->SetBranchAddress("l1_fired_time",&l1_fired_time);
        t->SetBranchAddress("l1_fired_pe",&l1_fired_pe);
    }

    //For bee:
    vector<double> op_t;
    vector<double> op_peTotal;
    vector<vector<double> > op_pes;
    vector<vector<double> > op_l1_t;
    vector<vector<double> > op_l1_pe;

    int nEntries = t->GetEntries();
    int nFlash = nEntries;
    for (int i=0; i<nEntries; i++) {
        t->GetEntry(i);
        op_t.push_back(time);
        op_peTotal.push_back(total_PE);
        vector<double> tmp;
        op_pes.push_back(tmp);
        for (int j=0; j<32; j++) {
            op_pes[i].push_back(PE[j]);
        }
        op_l1_t.push_back(*l1_fired_time);
        op_l1_pe.push_back(*l1_fired_pe);
    }

    jsonFile << "{" << endl;

    jsonFile << fixed << setprecision(2);
    print_vector(jsonFile, op_t, "op_t");
    print_vector(jsonFile, op_peTotal, "op_peTotal");
    print_vector_vector(jsonFile, op_pes, "op_pes");
    print_vector_vector(jsonFile, op_l1_t, "op_l1_t");
    print_vector_vector(jsonFile, op_l1_pe, "op_l1_pe");


    // read the flash match tree
    int tpc_cluster_id;
    int flash_id;
    double pe_pred[32];
    t = (TTree*)rootFile->Get("T_match");
    if (t) {
        t->SetBranchAddress("tpc_cluster_id", &tpc_cluster_id);
        t->SetBranchAddress("flash_id", &flash_id);
        t->SetBranchAddress("pe_pred", &pe_pred);

        vector<vector<double> > op_cluster_ids;
        vector<vector<double> > op_pes_pred;
        vector<double> op_nomatching_cluster_ids;
        for (int i=0; i<nFlash; i++) { // each Flash has some properties
            // vector<double> tmp;
            op_cluster_ids.push_back(vector<double>());
            op_pes_pred.push_back(vector<double>());
            for (int j=0; j<32; j++) {
                op_pes_pred[i].push_back(0);
            }
        }

        nEntries = t->GetEntries();
        for (int i=0; i<nEntries; i++) {
            t->GetEntry(i);
            if (flash_id>=0 && flash_id<nFlash) {
                op_cluster_ids.at(flash_id).push_back(tpc_cluster_id);
                for (int j=0; j<32; j++) {
                    op_pes_pred.at(flash_id).at(j) += pe_pred[j];
                }
            }
            else if (flash_id == -1) {
                op_nomatching_cluster_ids.push_back(tpc_cluster_id);
            }

        }

        print_vector_vector(jsonFile, op_pes_pred, "op_pes_pred");
        jsonFile << fixed << setprecision(0);
        print_vector_vector(jsonFile, op_cluster_ids, "op_cluster_ids");
        print_vector(jsonFile, op_nomatching_cluster_ids, "op_nomatching_cluster_ids");
    }


    // always dump runinfo in the end
    DumpRunInfo();

    jsonFile << "}" << endl;


}

//----------------------------------------------------------------
void WCReader::DumpSpacePoints(TString option)
{
    double x=0, y=0, z=0, q=0, nq=1;
    int cluster_id=0, real_cluster_id=0;
    vector<double> vx, vy, vz, vq, vnq, vcluster_id, vreal_cluster_id;
    TTree * t = 0;

    if (option == "truth" || option == "true") {
        t = (TTree*)rootFile->Get("T_true");
    }
    else if (option == "rec_simple" || option == "simple") {
        t = (TTree*)rootFile->Get("T_rec");
    }
    else if (option == "rec_charge_blob" || option == "charge") {
        t = (TTree*)rootFile->Get("T_rec_charge");
    }
    else if (option == "rec_charge_cell" || option == "deblob") {
        t = (TTree*)rootFile->Get("T_rec_charge_blob");
    }
    else if (option == "cluster") {
        t = (TTree*)rootFile->Get("T_cluster");
    }
    else {
        cout << "WARNING: Wrong option: " << option << endl;
    }

    if (t) {
        t->SetBranchAddress("x", &x);
        t->SetBranchAddress("y", &y);
        t->SetBranchAddress("z", &z);

        // if (! (option.Contains("simple") ||  option.Contains("cluster")) ) {
        //     t->SetBranchAddress("q", &q);
        // }
        if (t->GetBranch("q")) {
            t->SetBranchAddress("q", &q);
        }
        if (option.Contains("charge")) {
            t->SetBranchAddress("nq", &nq);
        }
	bool flag_real_cluster_id = false;
        if (t->GetBranch("cluster_id")) {
            t->SetBranchAddress("cluster_id", &cluster_id);
            if (t->GetBranch("real_cluster_id")) {
	      flag_real_cluster_id = true;
                t->SetBranchAddress("real_cluster_id", &real_cluster_id);
            }
        }
        int nPoints = t->GetEntries();
        for (int i=0; i<nPoints; i++) {
	  t->GetEntry(i);
	  if (q<0) q = 1;
	  vx.push_back(x);
	  vy.push_back(y);
	  vz.push_back(z);
	  vq.push_back(q);
	  vnq.push_back(nq);
	  vcluster_id.push_back(cluster_id);
	  
	  if (flag_real_cluster_id)
	    vreal_cluster_id.push_back(real_cluster_id);
	  else
	    vreal_cluster_id.push_back(cluster_id);
        }
    }

    jsonFile << "{" << endl;

    jsonFile << fixed << setprecision(1);

    print_vector(jsonFile, vx, "x");
    print_vector(jsonFile, vy, "y");
    print_vector(jsonFile, vz, "z");

    jsonFile << fixed << setprecision(0);
    print_vector(jsonFile, vq, "q");
    print_vector(jsonFile, vnq, "nq");
    print_vector(jsonFile, vcluster_id, "cluster_id");
    print_vector(jsonFile, vreal_cluster_id, "real_cluster_id");


    jsonFile << '"' << "type" << '"' << ":" << '"' << option << '"' << "," << endl;

    // always dump runinfo in the end
    DumpRunInfo();

    jsonFile << "}" << endl;

}

//----------------------------------------------------------------
void WCReader::DumpVtx()
{
    double x=0, y=0, z=0;
    int type = -1;
    // 1: Steiner-inspired graph
    // 2: Steiner terminals
    // 3: end point (extreme point, connecting to one object)
    // 4: kink (connecting to two objects
    // 5: vertex (connecting to three or more objects)
    int flag_main = -1;
    int cluster_id=0;
    vector<int> *sub_cluster_ids = new std::vector<int>;

    vector<double> vx, vy, vz, vtype, vflag_main, vcluster_id;
    vector<vector<double> > vsub_cluster_ids;
    TTree * t = (TTree*)rootFile->Get("T_vtx");

    if (t) {
        t->SetBranchAddress("x", &x);
        t->SetBranchAddress("y", &y);
        t->SetBranchAddress("z", &z);
        t->SetBranchAddress("type", &type);
        t->SetBranchAddress("flag_main", &flag_main);
        t->SetBranchAddress("cluster_id", &cluster_id);
        t->SetBranchAddress("sub_cluster_ids", &sub_cluster_ids);
        
        int nPoints = t->GetEntries();
        for (int i=0; i<nPoints; i++) {
            t->GetEntry(i);
            vx.push_back(x);
            vy.push_back(y);
            vz.push_back(z);
            vtype.push_back(type);
            vflag_main.push_back(flag_main);
            vcluster_id.push_back(cluster_id);

            vsub_cluster_ids.push_back(vector<double>());
            int size = sub_cluster_ids->size();
            for (int j=0; j<size; j++) {
	      vsub_cluster_ids.back().push_back( sub_cluster_ids->at(j) );
            }
        }
    }

    jsonFile << "{" << endl;

    jsonFile << fixed << setprecision(1);

    print_vector(jsonFile, vx, "x");
    print_vector(jsonFile, vy, "y");
    print_vector(jsonFile, vz, "z");

    jsonFile << fixed << setprecision(0);
    print_vector(jsonFile, vtype, "type");
    print_vector(jsonFile, vflag_main, "flag_main");
    print_vector(jsonFile, vcluster_id, "cluster_id");
    print_vector_vector(jsonFile, vsub_cluster_ids, "sub_cluster_ids");

    // always dump runinfo in the end
    DumpRunInfo();

    jsonFile << "}" << endl;

}

//----------------------------------------------------------------
void WCReader::DumpMC()
{
    TTree *t = (TTree*)rootFile->Get("TMC");

    if (t) {
        t->SetBranchAddress("mc_Ntrack"       , &mc_Ntrack);
        t->SetBranchAddress("mc_id"           , &mc_id);
        t->SetBranchAddress("mc_pdg"          , &mc_pdg);
        t->SetBranchAddress("mc_process"      , &mc_process);
        t->SetBranchAddress("mc_mother"       , &mc_mother);
        t->SetBranchAddress("mc_daughters"    , &mc_daughters);
        t->SetBranchAddress("mc_startXYZT"    , &mc_startXYZT);
        t->SetBranchAddress("mc_endXYZT"      , &mc_endXYZT);
        t->SetBranchAddress("mc_startMomentum", &mc_startMomentum);
        t->SetBranchAddress("mc_endMomentum"  , &mc_endMomentum);
        t->GetEntry(0);
        ProcessTracks();
        DumpMCJSON(jsonFile);
    }
}

//----------------------------------------------------------------
void WCReader::ProcessTracks()
{
    // map track id to track index in the array
    for (int i=0; i<mc_Ntrack; i++) {
        trackIndex[mc_id[i]] = i;
    }

    // in trackParents, trackChildren, trackSiblings vectors, store track index (not track id)
    for (int i=0; i<mc_Ntrack; i++) {
        // currently, parent size == 1;
        // for primary particle, parent id = 0;
        vector<int> parents;
        if ( !IsPrimary(i) ) {
            parents.push_back(trackIndex[mc_mother[i]]);
        }
        trackParents.push_back(parents); // primary track will have 0 parents

        vector<int> children;
        int nChildren = (*mc_daughters).at(i).size();
        for (int j=0; j<nChildren; j++) {
            children.push_back(trackIndex[(*mc_daughters).at(i).at(j)]);
        }
        trackChildren.push_back(children);

    }

    // siblings
    for (int i=0; i<mc_Ntrack; i++) {
        vector<int> siblings;
        if ( IsPrimary(i) ) {
            for (int j=0; j<mc_Ntrack; j++) {
                if( IsPrimary(j) ) {
                    siblings.push_back(j);
                }
            }
        }
        else {
            // siblings are simply children of the mother
            int mother = trackIndex[mc_mother[i]];
            int nSiblings = trackChildren.at(mother).size();
            for (int j=0; j<nSiblings; j++) {
                siblings.push_back(trackChildren.at(mother).at(j));
            }
        }
        trackSiblings.push_back(siblings);
    }

}

//----------------------------------------------------------------
bool WCReader::DumpMCJSON(int id, ostream& out)
{
    int i = trackIndex[id];
    if (!KeepMC(i)) return false;

    int e = KE(mc_startMomentum[i])*1000;

    int nDaughter = (*mc_daughters).at(i).size();
    vector<int> saved_daughters;
    for (int j=0; j<nDaughter; j++) {
        int daughter_id = (*mc_daughters).at(i).at(j);
        // int e_daughter = KE(mc_startMomentum[ trackIndex[daughter_id] ])*1000;
        // if (e_daughter >= thresh_KE) {
        if ( KeepMC(trackIndex[daughter_id]) ) {
            saved_daughters.push_back(daughter_id);
        }
    }

    out << fixed << setprecision(1);
    out << "{";

    out << "\"id\":" << id << ",";
    out << "\"text\":" << "\"" << PDGName(mc_pdg[i]) << "  " << e << " MeV\",";
    out << "\"data\":{";
    out << "\"start\":[" << mc_startXYZT[i][0] << ", " <<  mc_startXYZT[i][1] << ", " << mc_startXYZT[i][2] << "],";
    out << "\"end\":[" << mc_endXYZT[i][0] << ", " <<  mc_endXYZT[i][1] << ", " << mc_endXYZT[i][2] << "]";
    out << "},";
    out << "\"children\":[";
    int nSavedDaughter = saved_daughters.size();
    if (nSavedDaughter == 0) {
        out << "],";
        out << "\"icon\":" << "\"jstree-file\"";
        out << "}";
        return true;
    }
    else {
        for (int j=0; j<nSavedDaughter; j++) {
            DumpMCJSON(saved_daughters.at(j), out);
            if (j!=nSavedDaughter-1) {
                out << ",";
            }
        }
        out << "]";
        out << "}";
        return true;
    }
}

//----------------------------------------------------------------
void WCReader::DumpMCJSON(ostream& out)
{
    out << "[";
    vector<int> primaries;
    for (int i=0; i<mc_Ntrack; i++) {
        if (IsPrimary(i)) {
            // int e = KE(mc_startMomentum[i])*1000;
            // if (e<thresh_KE) continue;
            if (KeepMC(i)) {
                primaries.push_back(i);
            }
        }
    }
    int size = primaries.size();
    // cout << size << endl;
    for (int i=0; i<size; i++) {
        if (DumpMCJSON(mc_id[primaries[i]], out) && i!=size-1) {
            out << ", ";
        }
    }

    out << "]";
}

//----------------------------------------------------------------
bool WCReader::KeepMC(int i)
{
    double e = KE(mc_startMomentum[i])*1000;
    double thresh_KE_em = 0.; // MeV  original 5
    double thresh_KE_np = 0; // MeV   original 10
    if (mc_pdg[i]==22 || mc_pdg[i]==11 || mc_pdg[i]==-11) {
        if (e>=thresh_KE_em) return true;
        else return false;
    }
    else if (mc_pdg[i]==2112 || mc_pdg[i]==2212 || mc_pdg[i]>1e9) {
        if (e>=thresh_KE_np) return true;
        else return false;
    }
    return true;
}

//----------------------------------------------------------------
double WCReader::KE(float* momentum)
{
    TLorentzVector particle(momentum);
    return particle.E()-particle.M();
}

//----------------------------------------------------------------
TString WCReader::PDGName(int pdg)
{
    TParticlePDG *p = dbPDG->GetParticle(pdg);
    if (p == 0) {
        if (pdg>1e9) {
            int z = (pdg - 1e9) / 10000;
            int a = (pdg - 1e9 - z*1e4) / 10;
            TString name;
            if (z == 18) name = "Ar";

            else if (z == 17) name = "Cl";
            else if (z == 19) name = "Ca";
            else if (z == 16) name = "S";
            else if (z == 15) name = "P";
            else if (z == 14) name = "Si";
            else if (z == 1) name = "H";
            else if (z == 2) name = "He";

            else return pdg;
            return Form("%s-%i", name.Data(), a);
        }
        return pdg;
    }
    else {
        return p->GetName();
    }
}

//----------------------------------------------------------------
void WCReader::print_vector(ostream& out, vector<double>& v, TString desc, bool end)
{
    int N = v.size();

    out << '"' << desc << '"' << ":[";
    for (int i=0; i<N; i++) {
        out << v[i];
        if (i!=N-1) {
            out << ",";
        }
    }
    out << "]";
    if (!end) out << ",";
    out << endl;
}

void WCReader::print_vector_vector(ostream& out, vector<vector<double> >& vv, TString desc, bool end)
{
    int NN = vv.size();

    out << '"' << desc << '"' << ":[";
    for (int i=0; i!=NN; i++) {
        out << "[";

        int N = vv[i].size();
        for (int j=0; j!=N; j++) {
            double value = vv[i][j];
            out << value;
            if (j!=N-1) {
                out << ",";
            }
        }

        out << "]";
        if (i!=NN-1) {
            out << ",";
        }

    }
    out << "]";
    if (!end) out << ",";
    out << endl;
}
