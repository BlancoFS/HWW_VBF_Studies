#include "LatinoAnalysis/MultiDraw/interface/TTreeFunction.h"
#include "LatinoAnalysis/MultiDraw/interface/FunctionLibrary.h"
#include <vector>
#include "TVector2.h"
#include "TLorentzVector.h"
#include "Math/Vector4Dfwd.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"
#include <iostream>
#include <TMath.h>
#include <math.h>
#include "JHUGenMELA/MELA/interface/Mela.h"
#include "momemta/ConfigurationReader.h"
#include "momemta/MoMEMta.h"
#include "momemta/Types.h"

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/Rotation3D.h"
#include "Math/EulerAngles.h"
#include "Math/AxisAngle.h"
#include "Math/Quaternion.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/RotationZ.h"
#include "Math/RotationZYX.h"
#include "Math/LorentzRotation.h"
#include "Math/Boost.h"
#include "Math/BoostX.h"
#include "Math/BoostY.h"
#include "Math/BoostZ.h"
#include "Math/Transform3D.h"
#include "Math/Plane3D.h"
#include "Math/VectorUtil.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMath.h"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


using LorentzVectorM = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>;

using namespace momemta;


void normalizeInput_TMVA(LorentzVector& p4) {
  if (p4.M() > 0)
    return;

  // Increase the energy until M is positive                                                                                                                                 
  p4.SetE(p4.P());
  while (p4.M2() < 0) {
    double delta = p4.E() * 1e-5;
    p4.SetE(p4.E() + delta);
  };
}


float mlj_TMVA(float CleanJet_pt, float CleanJet_phi, float CleanJet_eta, float Lepton_pt, float Lepton_phi, float Lepton_eta){

  TLorentzVector j;
  TLorentzVector l;

  j.SetPtEtaPhiM(CleanJet_pt, CleanJet_eta, CleanJet_phi, 0.0);
  l.SetPtEtaPhiM(Lepton_pt, Lepton_eta, Lepton_phi, 0.0);

  double mass = (j + l).M();
  return mass;
}



float RecoMELA_TMVA(float nCleanJet, float nLepton, float PuppiMet_pt, float PuppiMet_phi, float Lepton_pt0, float Lepton_pt1, float Lepton_phi0, float Lepton_phi1, float Lepton_eta0, float Lepton_eta1,float CleanJet_pt0, float CleanJet_pt1, float CleanJet_phi0, float CleanJet_phi1, float CleanJet_eta0, float CleanJet_eta1, float Lepton_pdg0, float Lepton_pdg1, string name){


  Double_t LHCsqrts_= 13., mh_= 125.;
  TVar::VerbosityLevel verbosity_ = TVar::SILENT;

  static Mela* mela = new Mela(LHCsqrts_, mh_, verbosity_);

  TLorentzVector L1(0.,0.,0.,0.);
  TLorentzVector L2(0.,0.,0.,0.);
  TLorentzVector LL(0.,0.,0.,0.);
  TLorentzVector NuNu(0.,0.,0.,0.);
  TLorentzVector Higgs(0.,0.,0.,0.);
  TLorentzVector J1(0.,0.,0.,0.);
  TLorentzVector J2(0.,0.,0.,0.);

  if(nCleanJet >= 2 && nLepton > 1){

    if (Lepton_pdg0*Lepton_pdg1 != -11*13) return -9999.;

    L1.SetPtEtaPhiM(Lepton_pt0, Lepton_eta0, Lepton_phi0, 0.0);
    L2.SetPtEtaPhiM(Lepton_pt1, Lepton_eta1, Lepton_phi1, 0.0);

    J1.SetPtEtaPhiM(CleanJet_pt0, CleanJet_eta0, CleanJet_phi0, 0.0);
    J2.SetPtEtaPhiM(CleanJet_pt1, CleanJet_eta1, CleanJet_phi1, 0.0);

    LL = L1 + L2;

    double nunu_px = PuppiMet_pt*cos(PuppiMet_phi);
    double nunu_py = PuppiMet_pt*sin(PuppiMet_phi);
    double nunu_pz = LL.Pz();
    double nunu_m = 30.0; //Why 30? --> https://indico.cern.ch/event/850505/contributions/3593915/                                                                                                         

    double nunu_e = sqrt(nunu_px*nunu_px + nunu_py*nunu_py + nunu_pz*nunu_pz + nunu_m*nunu_m);
    NuNu.SetPxPyPzE(nunu_px, nunu_py, nunu_pz, nunu_e);
    Higgs = LL + NuNu;

    SimpleParticleCollection_t daughter;
    SimpleParticleCollection_t associated;
    SimpleParticleCollection_t mother;

    daughter.push_back(SimpleParticle_t(25, Higgs));

    associated.push_back(SimpleParticle_t(0,J1));
    associated.push_back(SimpleParticle_t(0,J2));

    if (Higgs.Pt() == 0 || Higgs.M()==0 || Lepton_pt0 < 10 || Lepton_pt1 < 10 || CleanJet_pt0 < 30 || CleanJet_pt1 < 30){
      return -9999;
    }

    mela->setCandidateDecayMode(TVar::CandidateDecay_WW); //Decay to WW                                                                                                                                    
    mela->setInputEvent(&daughter, &associated, 0, false);
    mela->setCurrentCandidateFromIndex(0);

    if (name=="VBF_QCD"){

      float RecoLevel_me_VBF_hsm = 0.;
      float RecoLevel_me_QCD_hsm = 0.;

      mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
      mela->computeProdP(RecoLevel_me_VBF_hsm, true);

      mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJQCD);
      mela->computeProdP(RecoLevel_me_QCD_hsm, true);

      float D_VBF_QCD = RecoLevel_me_VBF_hsm*RecoLevel_me_VBF_hsm / (RecoLevel_me_VBF_hsm*RecoLevel_me_VBF_hsm + RecoLevel_me_QCD_hsm*RecoLevel_me_QCD_hsm); 

      return D_VBF_QCD;

    }else if(name=="VBF_VH"){

      float RecoLevel_me_VBF_hsm = 0.;
      float RecoLevel_me_Had_WH_hsm = 0.;
      float RecoLevel_me_Had_ZH_hsm = 0.;
      float RecoLevel_me_VH_hsm = 0.;

      mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
      mela->computeProdP(RecoLevel_me_VBF_hsm, true);

      mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_ZH);
      mela->computeProdP(RecoLevel_me_Had_ZH_hsm, true);

      mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_WH);
      mela->computeProdP(RecoLevel_me_Had_WH_hsm, true);

      RecoLevel_me_VH_hsm = RecoLevel_me_Had_ZH_hsm + RecoLevel_me_Had_WH_hsm;

      float D_VBF_VH = (RecoLevel_me_VBF_hsm * RecoLevel_me_VBF_hsm) / (RecoLevel_me_VBF_hsm * RecoLevel_me_VBF_hsm  + 1e16 * (RecoLevel_me_VH_hsm * RecoLevel_me_VH_hsm));

      return D_VBF_VH;

    }else if(name=="VH_QCD"){

      float RecoLevel_me_QCD_hsm = 0.;
      float RecoLevel_me_Had_WH_hsm = 0.;
      float RecoLevel_me_Had_ZH_hsm = 0.;
      float RecoLevel_me_VH_hsm = 0.;

      mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJQCD);
      mela->computeProdP(RecoLevel_me_QCD_hsm, true);

      mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_ZH);
      mela->computeProdP(RecoLevel_me_Had_ZH_hsm, true);

      mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_WH);
      mela->computeProdP(RecoLevel_me_Had_WH_hsm, true);

      RecoLevel_me_VH_hsm = RecoLevel_me_Had_ZH_hsm + RecoLevel_me_Had_WH_hsm;

      float D_VH_QCD = (RecoLevel_me_VH_hsm * RecoLevel_me_VH_hsm) / ((RecoLevel_me_VH_hsm * RecoLevel_me_VH_hsm) + 1e-16 * (RecoLevel_me_QCD_hsm * RecoLevel_me_QCD_hsm));

      return D_VH_QCD;

    }else{

      return -9999.9;

    }


  }else{

    return -9999.;

  }

}


float RecoMoMEMta_TMVA(float nCleanJet, float nLepton, float PuppiMet_pt, float PuppiMet_phi, float Lepton_pt0, float Lepton_pt1, float Lepton_phi0, float Lepton_phi1, float Lepton_eta0, float Lepton_eta1,float CleanJet_pt0, float CleanJet_pt1, float CleanJet_phi0, float CleanJet_phi1, float CleanJet_eta0, float CleanJet_eta1, float Lepton_pdg0, float Lepton_pdg1){


  TLorentzVector L1(0.,0.,0.,0.);
  TLorentzVector L2(0.,0.,0.,0.);
  TLorentzVector LL(0.,0.,0.,0.);
  TLorentzVector NuNu(0.,0.,0.,0.);
  TLorentzVector Higgs(0.,0.,0.,0.);
  TLorentzVector J1(0.,0.,0.,0.);
  TLorentzVector J2(0.,0.,0.,0.);

  if(nCleanJet >= 2 && nLepton > 1){

    if (Lepton_pdg0*Lepton_pdg1 != -11*13) return -9999.;

    L1.SetPtEtaPhiM(Lepton_pt0, Lepton_eta0, Lepton_phi0, 0.0);
    L2.SetPtEtaPhiM(Lepton_pt1, Lepton_eta1, Lepton_phi1, 0.0);

    J1.SetPtEtaPhiM(CleanJet_pt0, CleanJet_eta0, CleanJet_phi0, 0.0);
    J2.SetPtEtaPhiM(CleanJet_pt1, CleanJet_eta1, CleanJet_phi1, 0.0);

    LL = L1 + L2;

    double nunu_px = PuppiMet_pt*cos(PuppiMet_phi);
    double nunu_py = PuppiMet_pt*sin(PuppiMet_phi);
    double nunu_pz = LL.Pz();
    double nunu_m = 30.0; //Why 30? --> https://indico.cern.ch/event/850505/contributions/3593915/                                                                                                         

    double nunu_e = sqrt(nunu_px*nunu_px + nunu_py*nunu_py + nunu_pz*nunu_pz + nunu_m*nunu_m);
    NuNu.SetPxPyPzE(nunu_px, nunu_py, nunu_pz, nunu_e);
    Higgs = LL + NuNu;


    momemta::Particle higgs { "higgs", LorentzVector(Higgs.Px(), Higgs.Py(), Higgs.Pz(), Higgs.E()), 25 }; // Higgs
    momemta::Particle Z { "Z", LorentzVector(Higgs.Px(), Higgs.Py(), Higgs.Pz(), Higgs.E()), 23 }; // Z, same 4 vector as Higgs                                                                            
    momemta::Particle jet1 { "jet1", LorentzVector(J1.Px(), J1.Py(), J1.Pz(), J1.E()), 1 };
    momemta::Particle jet2 { "jet2", LorentzVector(J2.Px(), J2.Py(), J2.Pz(), J2.E()), -1 };

    normalizeInput_TMVA(higgs.p4);
    normalizeInput_TMVA(Z.p4);
    normalizeInput_TMVA(jet1.p4);
    normalizeInput_TMVA(jet2.p4);


    logging::set_level(logging::level::off);

    // Higgs                                                                                                                                                                                               
    ConfigurationReader configuration_VBF("/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/qqH_hww_ME/higgs_jets.lua");
    MoMEMta weight_VBF(configuration_VBF.freeze());

    // DY                                                                                                                                                                                                 
    ConfigurationReader configuration_DY("/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/DY_ME/DY_ME.lua");
    MoMEMta weight_DY(configuration_DY.freeze());

    ParameterSet lua_parameters;
    lua_parameters.set("USE_TF", true);
    lua_parameters.set("USE_PERM", true);

    std::vector<std::pair<double, double>> weights_VBF = weight_VBF.computeWeights({higgs, jet1, jet2});
    std::vector<std::pair<double, double>> weights_DY = weight_DY.computeWeights({Z, jet1, jet2});

    double vbf = (double)weights_VBF.back().first;
    double dy = (double)weights_DY.back().first;

    return 150 * abs(vbf) / (150 * abs(vbf) + abs(dy));


  }else{
    
    return -9999.9;

  }
}



class TMVA_HWW: public multidraw::TTreeFunction {
public:
  TMVA_HWW(char const* name);

  char const* getName() const override { return "TMVA_HWW"; }
  TTreeFunction* clone() const override { return new TMVA_HWW(name_.c_str()); }

  unsigned getNdata() override { return 1; }
  double evaluate(unsigned) override;

protected:

  std::string name_;
  void bindTree_(multidraw::FunctionLibrary&) override;

  UIntValueReader*   nLepton{};
  UIntValueReader*   nCleanJet{};
  IntArrayReader*   Lepton_pdgId{};
  FloatArrayReader* Lepton_pt{};
  FloatArrayReader* Lepton_eta{};
  FloatArrayReader* Lepton_phi{};
  FloatArrayReader* CleanJet_pt{};
  FloatArrayReader* CleanJet_eta{};
  FloatArrayReader* CleanJet_phi{};
  FloatArrayReader* Jet_qgl{};
  FloatValueReader* metpt{};
  FloatValueReader* metphi{};
  FloatValueReader* mjj{};
  FloatValueReader* mll{};
  FloatValueReader* ptll{};
  FloatValueReader* detajj{};
  FloatValueReader* dphill{};
  FloatValueReader* dphijjmet{};
  FloatValueReader* dphilljj{};
  FloatValueReader* mti{};
  FloatValueReader* mtw1{};
  FloatValueReader* mtw2{};
  FloatValueReader* drll{};
  FloatValueReader* mth{};
  FloatValueReader* PuppiMET_pt{};
  FloatValueReader* PuppiMET_phi{};
  IntArrayReader* CleanJet_jetIdx{};
  FloatArrayReader* Jet_btagDeepFlavB{};
};


TMVA_HWW::TMVA_HWW(char const* name) :
  TTreeFunction()
{
  name_ = name;
}

double
TMVA_HWW::evaluate(unsigned)
{


  std::map<std::string,int> Use;
  Use["BDT"] = 1;
  Use["BDT1"] = 1;
  Use["BDTG4SK01F07"] = 1;
  Use["BDTG4SK01"] = 1;
  Use["BDTG4F07"] = 1;
  Use["BDT2"]  = 1;
  Use["BDTG4C3"] = 1;
  Use["BDTB"]   = 1;
  Use["BDTB2"]   = 1;
  Use["BDTB3"]  = 1;
  Use["BDTG4D3"]       = 1;
  Use["LikelihoodKDE"] = 1;
  Use["Likelihood"] = 1;

  unsigned njet = *nCleanJet->Get();

  float Ctot = log((abs(2 * Lepton_eta->At(0) - CleanJet_eta->At(0) - CleanJet_eta->At(1)) + abs(2 * Lepton_eta->At(1) - CleanJet_eta->At(0) - CleanJet_eta->At(1))) / *detajj->Get());
  float mlj_00 = mlj_TMVA(CleanJet_pt->At(0), CleanJet_phi->At(0), CleanJet_eta->At(0), Lepton_pt->At(0), Lepton_phi->At(0), Lepton_eta->At(0));
  float mlj_01 = mlj_TMVA(CleanJet_pt->At(0), CleanJet_phi->At(0), CleanJet_eta->At(0), Lepton_pt->At(1), Lepton_phi->At(1), Lepton_eta->At(1));
  float mlj_10 = mlj_TMVA(CleanJet_pt->At(1), CleanJet_phi->At(1), CleanJet_eta->At(1), Lepton_pt->At(0), Lepton_phi->At(0), Lepton_eta->At(0));
  float mlj_11 = mlj_TMVA(CleanJet_pt->At(1), CleanJet_phi->At(1), CleanJet_eta->At(1), Lepton_pt->At(1), Lepton_phi->At(1), Lepton_eta->At(1));

  float Jet_btagDeepFlavB_CleanJet_jetIdx_0_;
  float Jet_btagDeepFlavB_CleanJet_jetIdx_1_;
  if (njet==0){

  }else if (njet == 1){
    int jetIdx0 = CleanJet_jetIdx->At(0);
    Jet_btagDeepFlavB_CleanJet_jetIdx_0_ = jetIdx0 >= 0 ? Jet_btagDeepFlavB->At(jetIdx0) : -2;
  }
  else {
    int jetIdx0 = CleanJet_jetIdx->At(0);
    int jetIdx1 = CleanJet_jetIdx->At(1);
    Jet_btagDeepFlavB_CleanJet_jetIdx_0_ = jetIdx0 >= 0 ? Jet_btagDeepFlavB->At(jetIdx0) : -2;
    Jet_btagDeepFlavB_CleanJet_jetIdx_1_ = jetIdx1 >= 0 ? Jet_btagDeepFlavB->At(jetIdx1) : -2;
  }
  
  float D_VBF_QCD = RecoMELA_TMVA(*nCleanJet->Get(), *nLepton->Get(), *metpt->Get(), *metphi->Get(), Lepton_pt->At(0), Lepton_pt->At(1), Lepton_phi->At(0), Lepton_phi->At(1), Lepton_eta->At(0), Lepton_eta->At(1), CleanJet_pt->At(0), CleanJet_pt->At(1), CleanJet_phi->At(0), CleanJet_phi->At(1), CleanJet_eta->At(0), CleanJet_eta->At(1), Lepton_pdgId->At(0), Lepton_pdgId->At(1), "VBF_QCD");

  float D_VBF_VH = RecoMELA_TMVA(*nCleanJet->Get(), *nLepton->Get(), *metpt->Get(), *metphi->Get(), Lepton_pt->At(0), Lepton_pt->At(1), Lepton_phi->At(0), Lepton_phi->At(1), Lepton_eta->At(0), Lepton_eta->At(1), CleanJet_pt->At(0), CleanJet_pt->At(1), CleanJet_phi->At(0), CleanJet_phi->At(1), CleanJet_eta->At(0), CleanJet_eta->At(1), Lepton_pdgId->At(0), Lepton_pdgId->At(1), "VBF_VH");

  float D_QCD_VH = RecoMELA_TMVA(*nCleanJet->Get(), *nLepton->Get(), *metpt->Get(), *metphi->Get(), Lepton_pt->At(0), Lepton_pt->At(1), Lepton_phi->At(0), Lepton_phi->At(1), Lepton_eta->At(0), Lepton_eta->At(1), CleanJet_pt->At(0), CleanJet_pt->At(1), CleanJet_phi->At(0), CleanJet_phi->At(1), CleanJet_eta->At(0), CleanJet_eta->At(1), Lepton_pdgId->At(0), Lepton_pdgId->At(1), "VH_QCD");

  float D_VBF_DY = RecoMoMEMta_TMVA(*nCleanJet->Get(), *nLepton->Get(), *metpt->Get(), *metphi->Get(), Lepton_pt->At(0), Lepton_pt->At(1), Lepton_phi->At(0), Lepton_phi->At(1), Lepton_eta->At(0), Lepton_eta->At(1), CleanJet_pt->At(0), CleanJet_pt->At(1), CleanJet_phi->At(0), CleanJet_phi->At(1), CleanJet_eta->At(0), CleanJet_eta->At(1), Lepton_pdgId->At(0), Lepton_pdgId->At(1));


  float mjj_user{*mjj->Get()};
  float ctot_user = (float)Ctot;
  float detajj_user{*detajj->Get()};
  float drll_user{*drll->Get()};
  float jet1eta_user = CleanJet_eta->At(0);
  float jet2eta_user = CleanJet_eta->At(1);
  float puppimet_pt_user{*PuppiMET_pt->Get()};
  float puppimet_phi_user{*PuppiMET_phi->Get()};
  float mth_user{*mth->Get()};
  float ptll_user{*ptll->Get()};
  float mlj_00_user = (float)mlj_00;
  float mlj_01_user = (float)mlj_01;
  float mlj_10_user = (float)mlj_10;
  float mlj_11_user = (float)mlj_11;
  float mll_user{*mll->Get()};
  float btag_user = (float)Jet_btagDeepFlavB_CleanJet_jetIdx_0_;
  float btag_user_1 = (float)Jet_btagDeepFlavB_CleanJet_jetIdx_1_;
  float D_VBF_QCD_user = (float)D_VBF_QCD;
  float D_VBF_VH_user = (float)D_VBF_VH;
  float D_QCD_VH_user = (float)D_QCD_VH;
  float D_VBF_DY_user = (float)D_VBF_DY;
  
  TMVA::Reader *reader = new TMVA::Reader( "!Color:Silent" );

  reader->AddVariable("mjj", &mjj_user);
  reader->AddVariable("Ctot", &ctot_user);
  reader->AddVariable("detajj", &detajj_user);
  reader->AddVariable("drll", &drll_user);
  reader->AddVariable("jet1eta", &jet1eta_user);
  reader->AddVariable("jet2eta", &jet2eta_user);
  reader->AddVariable("PuppiMET_pt", &puppimet_pt_user);
  reader->AddVariable("PuppiMET_phi", &puppimet_phi_user);
  reader->AddVariable("mth", &mth_user);
  reader->AddVariable("ptll", &ptll_user);
  reader->AddVariable("mlj_00", &mlj_00_user);
  reader->AddVariable("mlj_01", &mlj_01_user);
  reader->AddVariable("mlj_10", &mlj_10_user);
  reader->AddVariable("mlj_11", &mlj_11_user);
  reader->AddVariable("mll", &mll_user);
  reader->AddVariable("btagDeepFlavB", &btag_user);
  reader->AddVariable("btagDeepFlavB_1", &btag_user_1);
  reader->AddVariable("D_VBF_QCD", &D_VBF_QCD_user);
  reader->AddVariable("D_VBF_VH", &D_VBF_VH_user);
  reader->AddVariable("D_QCD_VH", &D_QCD_VH_user);
  reader->AddVariable("D_VBF_DY", &D_VBF_DY_user);

  //TString dir    = "/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/WW/Full2016_v7/DNN/dataset/weights/";
  TString dir    = "/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/WW/Full2016_v7/DNN/dataset_for_test/weights/";
  TString prefix = "TMVAClassification";

  TString methodName = TString(name_) + TString(" method");
  TString weightfile = dir + prefix + TString("_") + TString(name_) + TString(".weights.xml");
  reader->BookMVA( methodName, weightfile );


  return reader->EvaluateMVA(methodName);


}

void
TMVA_HWW::bindTree_(multidraw::FunctionLibrary& _library)
{
  _library.bindBranch(nLepton, "nLepton");
  _library.bindBranch(nCleanJet, "nCleanJet");
  _library.bindBranch(Lepton_pdgId, "Lepton_pdgId");
  _library.bindBranch(Lepton_pt, "Lepton_pt");
  _library.bindBranch(Lepton_eta, "Lepton_eta");
  _library.bindBranch(Lepton_phi, "Lepton_phi");
  _library.bindBranch(CleanJet_pt, "CleanJet_pt");
  _library.bindBranch(CleanJet_eta, "CleanJet_eta");
  _library.bindBranch(CleanJet_phi, "CleanJet_phi");
  _library.bindBranch(metpt, "MET_pt");
  _library.bindBranch(metphi, "MET_phi");
  _library.bindBranch(Jet_qgl, "Jet_qgl");
  _library.bindBranch(mjj, "mjj");
  _library.bindBranch(mll, "mll");
  _library.bindBranch(ptll, "ptll");
  _library.bindBranch(detajj, "detajj");
  _library.bindBranch(dphill, "dphill");
  _library.bindBranch(dphijjmet, "dphijjmet");
  _library.bindBranch(dphilljj, "dphilljetjet");
  _library.bindBranch(PuppiMET_pt, "PuppiMET_pt");
  _library.bindBranch(PuppiMET_phi, "PuppiMET_phi");
  _library.bindBranch(mti, "mTi");
  _library.bindBranch(mth, "mth");
  _library.bindBranch(mtw1, "mtw2");
  _library.bindBranch(mtw2, "mtw1");
  _library.bindBranch(drll, "drll");
  _library.bindBranch(Jet_btagDeepFlavB,"Jet_btagDeepFlavB");
  _library.bindBranch(CleanJet_jetIdx,"CleanJet_jetIdx");
}
