#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TH1F.h"

using namespace std;
//using namespace LEEana;

int main( int argc, char** argv )
{
  if (argc < 3) {
    std::cout << "Usage: prune_weight_trees #txtInputList #outputROOT #GenieFluxKnob" << std::endl;
    std::cout << "Optional GenieFluxKnob:" << std::endl;
    std::cout << " 1 UBGenieFluxSmallUni" << std::endl;
    std::cout << "     All_UBGenie" << std::endl;
    std::cout << "     AxFFCCQEshape_UBGenie" << std::endl;
    std::cout << "     DecayAngMEC_UBGenie" << std::endl;
    std::cout << "     NormCCCOH_UBGenie" << std::endl;
    std::cout << "     NormNCCOH_UBGenie" << std::endl;
    std::cout << "     RPA_CCQE_UBGenie" << std::endl;
    std::cout << "     RootinoFix_UBGenie" << std::endl;
    std::cout << "     ThetaDelta2NRad_UBGenie" << std::endl;
    std::cout << "     Theta_Delta2Npi_UBGenie" << std::endl;
    std::cout << "     TunedCentralValue_UBGenie" << std::endl;
    std::cout << "     VecFFCCQEshape_UBGenie" << std::endl;
    std::cout << "     XSecShape_CCMEC_UBGenie" << std::endl;
    std::cout << "     splines_general_Spline" << std::endl;
    std::cout << "     xsr_scc_Fa3_SCC" << std::endl;
    std::cout << "     xsr_scc_Fv3_SCC" << std::endl;
    std::cout << " 2 expskin_FluxUnisim" << std::endl;
    std::cout << " 3 horncurrent_FluxUnisim" << std::endl;
    std::cout << " 4 kminus_PrimaryHadronNormalization" << std::endl;
    std::cout << " 5 kplus_PrimaryHadronFeynmanScaling" << std::endl;
    std::cout << " 6 kzero_PrimaryHadronSanfordWang" << std::endl;
    std::cout << " 7 nucleoninexsec_FluxUnisim" << std::endl;
    std::cout << " 8 nucleonqexsec_FluxUnisim" << std::endl;
    std::cout << " 9 nucleontotxsec_FluxUnisim" << std::endl;
    std::cout << "10 piminus_PrimaryHadronSWCentralSplineVariation" << std::endl;
    std::cout << "11 pioninexsec_FluxUnisim" << std::endl;
    std::cout << "12 pionqexsec_FluxUnisim" << std::endl;
    std::cout << "13 piontotxsec_FluxUnisim" << std::endl;
    std::cout << "14 piplus_PrimaryHadronSWCentralSplineVariation" << std::endl;
    std::cout << "15 reinteractions_piminus_Geant4" << std::endl;
    std::cout << "16 reinteractions_piplus_Geant4" << std::endl;
    std::cout << "17 reinteractions_proton_Geant4" << std::endl;

    return -1;
  }

  string filelist(argv[1]);
  string outName(argv[2]);
  string SaveGenieFluxKnob(argv[3]);

  auto T_wgt = new TChain("wcpweights/T_wgt");

  std::ifstream in(filelist.c_str());
  std::string line;
  while (std::getline(in, line)) {
    T_wgt->Add(line.c_str());
  }
  // cout << T_wgt->GetEntries() << endl;

  Int_t run, subrun, event;   
  std::string* file_type = new std::string;
  Float_t weight_cv, weight_spline, weight_lee;
  auto mcweight = new std::map<std::string, std::vector<float>>;
  bool mcweight_filled;

  T_wgt->SetBranchAddress("run", &run);
  T_wgt->SetBranchAddress("subrun", &subrun);
  T_wgt->SetBranchAddress("event", &event);
  T_wgt->SetBranchAddress("file_type", &file_type);
  T_wgt->SetBranchAddress("weight_cv", &weight_cv);
  T_wgt->SetBranchAddress("weight_spline", &weight_spline);
  T_wgt->SetBranchAddress("weight_lee", &weight_lee);
  T_wgt->SetBranchAddress("mcweight", &mcweight);
  T_wgt->SetBranchAddress("mcweight_filled", &mcweight_filled);
  // T_wgt->Print();


  std::vector<std::string> UBGenieFluxSmallUni = {
    "All_UBGenie"
    , "AxFFCCQEshape_UBGenie"
    , "DecayAngMEC_UBGenie"
    , "NormCCCOH_UBGenie"
    , "NormNCCOH_UBGenie"
    , "RPA_CCQE_UBGenie"
    , "RootinoFix_UBGenie"
    , "ThetaDelta2NRad_UBGenie"
    , "Theta_Delta2Npi_UBGenie"
    , "TunedCentralValue_UBGenie"
    , "VecFFCCQEshape_UBGenie"
    , "XSecShape_CCMEC_UBGenie"
    , "splines_general_Spline"
    , "xsr_scc_Fa3_SCC"
    , "xsr_scc_Fv3_SCC"
  };

  std::vector<std::string> UBGenieFluxBigUni = {
    "expskin_FluxUnisim"
    , "horncurrent_FluxUnisim"
    , "kminus_PrimaryHadronNormalization"
    , "kplus_PrimaryHadronFeynmanScaling"
    , "kzero_PrimaryHadronSanfordWang"
    , "nucleoninexsec_FluxUnisim"
    , "nucleonqexsec_FluxUnisim"
    , "nucleontotxsec_FluxUnisim"
    , "piminus_PrimaryHadronSWCentralSplineVariation"
    , "pioninexsec_FluxUnisim"
    , "pionqexsec_FluxUnisim"
    , "piontotxsec_FluxUnisim"
    , "piplus_PrimaryHadronSWCentralSplineVariation"
    , "reinteractions_piminus_Geant4"
    , "reinteractions_piplus_Geant4"
    , "reinteractions_proton_Geant4"
  };

  auto ofile = TFile::Open(outName.c_str(), "recreate");
  // ofile->mkdir("wcpselection")->cd();


  // set branch
  TTree* UBTree = nullptr;
  std::map<std::string, std::vector<float> > weight_f1;
  std::vector<float> weight_f2;
  std::map<std::string, TH1F*> hwmap;

  if (SaveGenieFluxKnob == "UBGenieFluxSmallUni") {
    UBTree = new TTree("UBGenieFluxSmallUni","UBGenieFluxSmallUni");
    UBTree->Branch("run", &run, "run/I");
    UBTree->Branch("subrun", &subrun, "subrun/I");
    UBTree->Branch("event", &event, "event/I");
    UBTree->Branch("file_type", &file_type);
    UBTree->Branch("weight_cv", &weight_cv, "weight_cv/F");
    UBTree->Branch("weight_spline", &weight_spline, "weight_spline/F");
    UBTree->Branch("weight_lee", &weight_lee, "weight_lee/F");
    UBTree->Branch("mcweight_filled", &mcweight_filled);

    for (auto const& knob : UBGenieFluxSmallUni) {
      weight_f1.emplace(knob, std::vector<float>());
      UBTree->Branch(knob.c_str(), &(weight_f1.at(knob)) );
      hwmap.emplace(knob, new TH1F(("h_"+knob).c_str(), ("h_"+knob).c_str(), 1000, 0,2));
    }

  }
  else {
    int N = std::count(UBGenieFluxBigUni.begin(), UBGenieFluxBigUni.end(), SaveGenieFluxKnob);
    if (N!=1) {
      std::cerr << "Warning, cannot find the knob: " << SaveGenieFluxKnob << std::endl;
      return 0;
    }

    UBTree = new TTree(SaveGenieFluxKnob.c_str(), SaveGenieFluxKnob.c_str());
    UBTree->Branch("run", &run, "run/I");
    UBTree->Branch("subrun", &subrun, "subrun/I");
    UBTree->Branch("event", &event, "event/I");
    UBTree->Branch("file_type", &file_type);
    UBTree->Branch("weight_cv", &weight_cv, "weight_cv/F");
    UBTree->Branch("weight_spline", &weight_spline, "weight_spline/F");
    UBTree->Branch("weight_lee", &weight_lee, "weight_lee/F");
    UBTree->Branch("mcweight_filled", &mcweight_filled);
    UBTree->Branch(SaveGenieFluxKnob.c_str(), &weight_f2);
    hwmap.emplace(SaveGenieFluxKnob, new TH1F(("h_"+SaveGenieFluxKnob).c_str(), ("h_"+SaveGenieFluxKnob).c_str(), 1000, 0,2));
  }

  // loop through events
  for (size_t i=0; i<T_wgt->GetEntries(); i++) {
    T_wgt->GetEntry(i);
    // for (auto const& ele: *mcweight) {
    //   cout << ele.first << " " << ele.second.size() << endl;
    // }

    for (auto x: weight_f1) { x.second.clear(); }
    weight_f2.clear();
 
    if (mcweight_filled>0) {
      // cout << mcweight_filled << endl;
      if (SaveGenieFluxKnob=="UBGenieFluxSmallUni") {
        for(auto const& knob: UBGenieFluxSmallUni) {
          // weight_f1->push_back(mcweight->at(knob));
          weight_f1.at(knob) = mcweight->at(knob);

          for (auto const& w: weight_f1.at(knob)) {
            hwmap.at(knob)->Fill(w);
          }
        }
      }
      else {
         weight_f2 = mcweight->at(SaveGenieFluxKnob);

         for (auto const& w: weight_f2) {
           hwmap.at(SaveGenieFluxKnob)->Fill(w);
         }
      }
    }
 
    UBTree->Fill();
  } // end of T_wgt loop

  UBTree->Write();
  for (auto x: hwmap) {
    x.second->Write();
  }
  ofile->Close();
 
  return 0;
}
