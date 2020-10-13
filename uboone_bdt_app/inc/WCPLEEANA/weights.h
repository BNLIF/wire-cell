#ifndef UBOONE_LEE_WEIGHTS
#define UBOONE_LEE_WEIGHTS

#include <string>

namespace LEEana{
  struct WeightInfo{
    int run;
    int subrun;
    int event;
    std::string *file_type;
    float weight_cv;
    float weight_spline;
    float weight_lee;
    bool mcweight_filled;
    
    std::vector<float> *expskin_FluxUnisim;
    std::vector<float> *horncurrent_FluxUnisim;
    std::vector<float> *kminus_PrimaryHadronNormalization;
    std::vector<float> *kplus_PrimaryHadronFeynmanScaling;
    std::vector<float> *kzero_PrimaryHadronSanfordWang;
    std::vector<float> *nucleoninexsec_FluxUnisim;
    std::vector<float> *nucleonqexsec_FluxUnisim;
    std::vector<float> *nucleontotxsec_FluxUnisim;
    std::vector<float> *piminus_PrimaryHadronSWCentralSplineVariation;
    std::vector<float> *pioninexsec_FluxUnisim;
    std::vector<float> *pionqexsec_FluxUnisim;
    std::vector<float> *piontotxsec_FluxUnisim;
    std::vector<float> *piplus_PrimaryHadronSWCentralSplineVariation;

    std::vector<float> *All_UBGenie;
    std::vector<float> *AxFFCCQEshape_UBGenie;
    std::vector<float> *DecayAngMEC_UBGenie;
    std::vector<float> *NormCCCOH_UBGenie;
    std::vector<float> *NormNCCOH_UBGenie;
    std::vector<float> *RPA_CCQE_Reduced_UBGenie;
    std::vector<float> *RPA_CCQE_UBGenie;
    std::vector<float> *RootinoFix_UBGenie;
    std::vector<float> *ThetaDelta2NRad_UBGenie;
    std::vector<float> *Theta_Delta2Npi_UBGenie;
    std::vector<float> *TunedCentralValue_UBGenie;
    std::vector<float> *VecFFCCQEshape_UBGenie;
    std::vector<float> *XSecShape_CCMEC_UBGenie;
    std::vector<float> *splines_general_Spline;
    std::vector<float> *xsr_scc_Fa3_SCC;
    std::vector<float> *xsr_scc_Fv3_SCC;
    
  };

  void set_tree_address(TTree *T, WeightInfo& weight, TString option);
  void put_tree_address(TTree *T, WeightInfo& weight, TString option);
  int get_size(WeightInfo& weight, TString option);
}

int LEEana::get_size(WeightInfo& weight, TString option){
  if (weight.mcweight_filled){
    if (option == "expskin_FluxUnisim"){
      return weight.expskin_FluxUnisim->size();
  }else if (option == "horncurrent_FluxUnisim"){
      return weight.horncurrent_FluxUnisim->size();
  }else if (option == "kminus_PrimaryHadronNormalization"){
      return weight.kminus_PrimaryHadronNormalization->size();
  }else if (option == "kplus_PrimaryHadronFeynmanScaling"){
      return weight.kplus_PrimaryHadronFeynmanScaling->size();
  }else if (option == "kzero_PrimaryHadronSanfordWang"){
      return weight.kzero_PrimaryHadronSanfordWang->size();
  }else if (option == "nucleoninexsec_FluxUnisim"){
      return weight.nucleoninexsec_FluxUnisim->size();
  }else if (option == "nucleonqexsec_FluxUnisim"){
      return weight.nucleonqexsec_FluxUnisim->size();
  }else if (option == "nucleontotxsec_FluxUnisim"){
      return weight.nucleontotxsec_FluxUnisim->size();
  }else if (option == "piminus_PrimaryHadronSWCentralSplineVariation"){
      return weight.piminus_PrimaryHadronSWCentralSplineVariation->size();
  }else if (option == "pioninexsec_FluxUnisim"){
      return weight.pioninexsec_FluxUnisim->size();
  }else if (option == "pionqexsec_FluxUnisim"){
      return weight.pionqexsec_FluxUnisim->size();
  }else if (option == "piontotxsec_FluxUnisim"){
      return weight.piontotxsec_FluxUnisim->size();
  }else if (option == "piplus_PrimaryHadronSWCentralSplineVariation"){
      return weight.piplus_PrimaryHadronSWCentralSplineVariation->size();
  }else if (option == "UBGenieFluxSmallUni"){
      return weight.All_UBGenie->size();
    }
  }else{
    return 0;
  }
}

void LEEana::set_tree_address(TTree *T, WeightInfo& weight, TString option){
  T->SetBranchAddress("run",&weight.run);
  T->SetBranchAddress("subrun",&weight.subrun);
  T->SetBranchAddress("event",&weight.event);
  T->SetBranchAddress("file_type",&weight.file_type);
  T->SetBranchAddress("weight_cv",&weight.weight_cv);
  T->SetBranchAddress("weight_spline",&weight.weight_spline);
  T->SetBranchAddress("weight_lee",&weight.weight_lee);
  T->SetBranchAddress("mcweight_filled",&weight.mcweight_filled);

  if (option == "expskin_FluxUnisim"){
    T->SetBranchAddress("expskin_FluxUnisim",&weight.expskin_FluxUnisim);
  }else if (option == "horncurrent_FluxUnisim"){
    T->SetBranchAddress("horncurrent_FluxUnisim",&weight.horncurrent_FluxUnisim);
  }else if (option == "kminus_PrimaryHadronNormalization"){
    T->SetBranchAddress("kminus_PrimaryHadronNormalization",&weight.kminus_PrimaryHadronNormalization);
  }else if (option == "kplus_PrimaryHadronFeynmanScaling"){
    T->SetBranchAddress("kplus_PrimaryHadronFeynmanScaling",&weight.kplus_PrimaryHadronFeynmanScaling);
  }else if (option == "kzero_PrimaryHadronSanfordWang"){
    T->SetBranchAddress("kzero_PrimaryHadronSanfordWang",&weight.kzero_PrimaryHadronSanfordWang);
  }else if (option == "nucleoninexsec_FluxUnisim"){
    T->SetBranchAddress("nucleoninexsec_FluxUnisim",&weight.nucleoninexsec_FluxUnisim);
  }else if (option == "nucleonqexsec_FluxUnisim"){
    T->SetBranchAddress("nucleonqexsec_FluxUnisim",&weight.nucleonqexsec_FluxUnisim);
  }else if (option == "nucleontotxsec_FluxUnisim"){
    T->SetBranchAddress("nucleontotxsec_FluxUnisim",&weight.nucleontotxsec_FluxUnisim);
  }else if (option == "piminus_PrimaryHadronSWCentralSplineVariation"){
    T->SetBranchAddress("piminus_PrimaryHadronSWCentralSplineVariation",&weight.piminus_PrimaryHadronSWCentralSplineVariation);
  }else if (option == "pioninexsec_FluxUnisim"){
    T->SetBranchAddress("pioninexsec_FluxUnisim",&weight.pioninexsec_FluxUnisim);
  }else if (option == "pionqexsec_FluxUnisim"){
    T->SetBranchAddress("pionqexsec_FluxUnisim",&weight.pionqexsec_FluxUnisim);
  }else if (option == "piontotxsec_FluxUnisim"){
    T->SetBranchAddress("piontotxsec_FluxUnisim",&weight.piontotxsec_FluxUnisim);
  }else if (option == "piplus_PrimaryHadronSWCentralSplineVariation"){
    T->SetBranchAddress("piplus_PrimaryHadronSWCentralSplineVariation",&weight.piplus_PrimaryHadronSWCentralSplineVariation);
  }else if (option == "UBGenieFluxSmallUni"){
    
    T->SetBranchAddress("All_UBGenie",&weight.All_UBGenie);
    T->SetBranchAddress("AxFFCCQEshape_UBGenie",&weight.AxFFCCQEshape_UBGenie);
    T->SetBranchAddress("DecayAngMEC_UBGenie",&weight.DecayAngMEC_UBGenie);
    T->SetBranchAddress("NormCCCOH_UBGenie",&weight.NormCCCOH_UBGenie);
    T->SetBranchAddress("NormNCCOH_UBGenie",&weight.NormNCCOH_UBGenie);
    T->SetBranchAddress("RPA_CCQE_Reduced_UBGenie",&weight.RPA_CCQE_Reduced_UBGenie);
    T->SetBranchAddress("RPA_CCQE_UBGenie",&weight.RPA_CCQE_UBGenie);
    T->SetBranchAddress("RootinoFix_UBGenie",&weight.RootinoFix_UBGenie);
    T->SetBranchAddress("ThetaDelta2NRad_UBGenie",&weight.ThetaDelta2NRad_UBGenie);
    T->SetBranchAddress("Theta_Delta2Npi_UBGenie",&weight.Theta_Delta2Npi_UBGenie);
    T->SetBranchAddress("TunedCentralValue_UBGenie",&weight.TunedCentralValue_UBGenie);
    T->SetBranchAddress("VecFFCCQEshape_UBGenie",&weight.VecFFCCQEshape_UBGenie);
    T->SetBranchAddress("XSecShape_CCMEC_UBGenie",&weight.XSecShape_CCMEC_UBGenie);
    T->SetBranchAddress("splines_general_Spline",&weight.splines_general_Spline);
    T->SetBranchAddress("xsr_scc_Fa3_SCC",&weight.xsr_scc_Fa3_SCC);
    T->SetBranchAddress("xsr_scc_Fv3_SCC",&weight.xsr_scc_Fv3_SCC);
  }
}
void LEEana::put_tree_address(TTree *T, WeightInfo& weight, TString option){

  T->Branch("run",&weight.run,"run/I");
  T->Branch("subrun",&weight.subrun,"subrun/I");
  T->Branch("event",&weight.event,"event/I");
  T->Branch("file_type",&weight.file_type);
  T->Branch("weight_cv",&weight.weight_cv,"weight_cv/F");
  T->Branch("weight_spline",&weight.weight_spline,"weight_spline/F");
  T->Branch("weight_lee",&weight.weight_lee,"weight_lee/F");
  T->Branch("mcweight_filled",&weight.mcweight_filled,"mcweight_filled/O");

  if (option == "expskin_FluxUnisim"){
    T->Branch("expskin_FluxUnisim",&weight.expskin_FluxUnisim);
  }else if (option == "horncurrent_FluxUnisim"){
    T->Branch("horncurrent_FluxUnisim",&weight.horncurrent_FluxUnisim);
  }else if (option == "kminus_PrimaryHadronNormalization"){
    T->Branch("kminus_PrimaryHadronNormalization",&weight.kminus_PrimaryHadronNormalization);
  }else if (option == "kplus_PrimaryHadronFeynmanScaling"){
    T->Branch("kplus_PrimaryHadronFeynmanScaling",&weight.kplus_PrimaryHadronFeynmanScaling);
  }else if (option == "kzero_PrimaryHadronSanfordWang"){
    T->Branch("kzero_PrimaryHadronSanfordWang",&weight.kzero_PrimaryHadronSanfordWang);
  }else if (option == "nucleoninexsec_FluxUnisim"){
    T->Branch("nucleoninexsec_FluxUnisim",&weight.nucleoninexsec_FluxUnisim);
  }else if (option == "nucleonqexsec_FluxUnisim"){
    T->Branch("nucleonqexsec_FluxUnisim",&weight.nucleonqexsec_FluxUnisim);
  }else if (option == "nucleontotxsec_FluxUnisim"){
    T->Branch("nucleontotxsec_FluxUnisim",&weight.nucleontotxsec_FluxUnisim);
  }else if (option == "piminus_PrimaryHadronSWCentralSplineVariation"){
    T->Branch("piminus_PrimaryHadronSWCentralSplineVariation",&weight.piminus_PrimaryHadronSWCentralSplineVariation);
  }else if (option == "pioninexsec_FluxUnisim"){
    T->Branch("pioninexsec_FluxUnisim",&weight.pioninexsec_FluxUnisim);
  }else if (option == "pionqexsec_FluxUnisim"){
    T->Branch("pionqexsec_FluxUnisim",&weight.pionqexsec_FluxUnisim);
  }else if (option == "piontotxsec_FluxUnisim"){
    T->Branch("piontotxsec_FluxUnisim",&weight.piontotxsec_FluxUnisim);
  }else if (option == "piplus_PrimaryHadronSWCentralSplineVariation"){
    T->Branch("piplus_PrimaryHadronSWCentralSplineVariation",&weight.piplus_PrimaryHadronSWCentralSplineVariation);
  }else if (option == "UBGenieFluxSmallUni"){
    
    T->Branch("All_UBGenie",&weight.All_UBGenie);
    T->Branch("AxFFCCQEshape_UBGenie",&weight.AxFFCCQEshape_UBGenie);
    T->Branch("DecayAngMEC_UBGenie",&weight.DecayAngMEC_UBGenie);
    T->Branch("NormCCCOH_UBGenie",&weight.NormCCCOH_UBGenie);
    T->Branch("NormNCCOH_UBGenie",&weight.NormNCCOH_UBGenie);
    T->Branch("RPA_CCQE_Reduced_UBGenie",&weight.RPA_CCQE_Reduced_UBGenie);
    T->Branch("RPA_CCQE_UBGenie",&weight.RPA_CCQE_UBGenie);
    T->Branch("RootinoFix_UBGenie",&weight.RootinoFix_UBGenie);
    T->Branch("ThetaDelta2NRad_UBGenie",&weight.ThetaDelta2NRad_UBGenie);
    T->Branch("Theta_Delta2Npi_UBGenie",&weight.Theta_Delta2Npi_UBGenie);
    T->Branch("TunedCentralValue_UBGenie",&weight.TunedCentralValue_UBGenie);
    T->Branch("VecFFCCQEshape_UBGenie",&weight.VecFFCCQEshape_UBGenie);
    T->Branch("XSecShape_CCMEC_UBGenie",&weight.XSecShape_CCMEC_UBGenie);
    T->Branch("splines_general_Spline",&weight.splines_general_Spline);
    T->Branch("xsr_scc_Fa3_SCC",&weight.xsr_scc_Fa3_SCC);
    T->Branch("xsr_scc_Fv3_SCC",&weight.xsr_scc_Fv3_SCC);
  }
}



#endif
