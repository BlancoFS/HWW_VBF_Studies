#include "LatinoAnalysis/MultiDraw/interface/TTreeFunction.h"
#include "LatinoAnalysis/MultiDraw/interface/FunctionLibrary.h"
#include "TSystem.h"
#include "iostream"
#include "vector"
#include "TLorentzVector.h"
#include "TMath.h"
//#include "ZZMatrixElement/MELA/interface/Mela.h" //Use ZZMatrixElement or JHUGen(recommended)
#include "JHUGenMELA/MELA/interface/Mela.h"
#include "TSystem.h"
#include "TString.h"


float RecoMELA_VBF(float nCleanJet, float nLepton, float PuppiMet_pt, float PuppiMet_phi, float Lepton_pt0, float Lepton_pt1, float Lepton_phi0, float Lepton_phi1, float Lepton_eta0, float Lepton_eta1,float CleanJet_pt0, float CleanJet_pt1, float CleanJet_phi0, float CleanJet_phi1, float CleanJet_eta0, float CleanJet_eta1, float Lepton_pdg0, float Lepton_pdg1){

    
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

    if (Lepton_pdg0*Lepton_pdg1 != -11*13 || PuppiMet_pt==0) return -9999.;

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

    //mela->setCandidateDecayMode(TVar::CandidateDecay_WW); //Decay to WW
    mela->setCandidateDecayMode(TVar::CandidateDecay_Stable);
    mela->setInputEvent(&daughter, &associated, 0, false);
    mela->setCurrentCandidateFromIndex(0);

    float RecoLevel_me_VBF_hsm = 0.;

    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
    mela->computeProdP(RecoLevel_me_VBF_hsm, true);    

    float RecoLevel_me_QCD_hsm = 0.;

    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJQCD);
    mela->computeProdP(RecoLevel_me_QCD_hsm, true);

    mela->resetInputEvent(); 

    float vbf = RecoLevel_me_VBF_hsm;
    float ggh = RecoLevel_me_QCD_hsm;
    
    return (float)vbf*vbf/(vbf*vbf + ggh*ggh);
    
    
  }else{
    
    return -9999.;
    
  }

}

