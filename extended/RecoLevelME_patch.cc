#include "LatinoAnalysis/MultiDraw/interface/TTreeFunction.h"
#include "LatinoAnalysis/MultiDraw/interface/FunctionLibrary.h"
#include "TSystem.h"
#include "iostream"
#include "vector"
#include "TLorentzVector.h"
#include "TMath.h"
//#include "ZZMatrixElement/MELA/interface/Mela.h"
#include "JHUGenMELA/MELA/interface/Mela.h"
#include "TSystem.h"
#include <map>
#include "TString.h"

class RecoLevelME : public multidraw::TTreeFunction {
public:
  //Class Constructor 
  RecoLevelME(char const* name);
  //Class Destructor 
  ~RecoLevelME() {
  }
  //Functions from Multidraw namespace (TTreeFunction class)
  char const* getName() const override {return "RecoLevelME"; }
  TTreeFunction* clone() const override {return new RecoLevelME(name_.c_str());}
  unsigned getNdata() override {return 1; }
  //This function will return the required value
  double evaluate(unsigned) override;

protected:
  void bindTree_(multidraw::FunctionLibrary&) override;

  //name of the required ME
  std::string name_;

  //Needed variables to select the events
  UIntValueReader*  nCleanJet{};
  FloatArrayReader* CleanJet_pt{};
  FloatArrayReader* CleanJet_eta{};
  FloatArrayReader* CleanJet_phi{};
  //FloatArrayReader* CleanJet_mass{};

  IntArrayReader* Lepton_pdgId{};
  UIntValueReader*  nLepton{};
  FloatArrayReader* Lepton_pt{};
  FloatArrayReader* Lepton_eta{};
  FloatArrayReader* Lepton_phi{};

  FloatValueReader* MET_pt{};
  FloatValueReader* PuppiMET_pt{};
  FloatValueReader* PuppiMET_phi{};

  UIntValueReader*  nSubJet{};
  FloatArrayReader* SubJet_pt{};
  FloatArrayReader* SubJet_eta{};
  FloatArrayReader* SubJet_mass{};
  FloatArrayReader* SubJet_phi{};

private:

  Double_t LHCsqrts_= 13., mh_= 125.;
  TVar::VerbosityLevel verbosity_ = TVar::SILENT;
  
  static Mela* mela;

};
Mela* RecoLevelME :: mela = 0;

RecoLevelME::RecoLevelME(char const* name):
  TTreeFunction()
{
  name_ = name;
  if(mela == 0){
    mela = new Mela(LHCsqrts_, mh_, verbosity_);
  }
}

double
RecoLevelME::evaluate(unsigned)
{
  //Map to store the ME
  std::map<TString, float> MatrixElementsMap;

  //Initializing 4-vectors
  TLorentzVector L1(0.,0.,0.,0.);
  TLorentzVector L2(0.,0.,0.,0.);
  TLorentzVector LL(0.,0.,0.,0.);
  TLorentzVector NuNu(0.,0.,0.,0.);
  TLorentzVector Higgs(0.,0.,0.,0.);
  TLorentzVector J1(0.,0.,0.,0.);
  TLorentzVector J2(0.,0.,0.,0.);
  TLorentzVector J3(0.,0.,0.,0.);

  //Getting some values to select the events
  unsigned ncleanjet{*nCleanJet->Get()};
  unsigned nlep{*nLepton->Get()};
  unsigned nsubjet{*nSubJet->Get()};
  float Pmet_pt{*PuppiMET_pt->Get()};
  float Pmet_phi{*PuppiMET_phi->Get()};

  //Conditions to select the event
  if(ncleanjet>=2 && nlep>1){
    //STEP-1
    //4-vectors of the leptons

    int muons = 0;
    int electrons = 0;
    int lep1 = 0;
    int lep2 = 0;
    for (unsigned int ilep = 0; ilep<nlep; ilep++){
     if (abs(Lepton_pdgId->At(ilep)) == 13){
    	++muons;
    	if (muons == 1 && Lepton_pt->At(ilep) > 13){
    	  L1.SetPtEtaPhiM(Lepton_pt->At(ilep), Lepton_eta->At(ilep), Lepton_phi->At(ilep), 0.0);
	  lep1 = Lepton_pdgId->At(ilep);
    	}
      }
      if (abs(Lepton_pdgId->At(ilep)) == 11){
    	++electrons;
    	if (electrons == 1 && Lepton_pt->At(ilep) > 13){
    	  L2.SetPtEtaPhiM(Lepton_pt->At(ilep), Lepton_eta->At(ilep), Lepton_phi->At(ilep), 0.0);
	  lep2 = Lepton_pdgId->At(ilep);
    	}
     }
    }

    if (muons<1 || electrons<1){
      return -9999;
    }
    //L1.SetPtEtaPhiM(Lepton_pt->At(0), Lepton_eta->At(0), Lepton_phi->At(0), 0.0);
    //L2.SetPtEtaPhiM(Lepton_pt->At(1), Lepton_eta->At(1), Lepton_phi->At(1), 0.0);
    LL = L1 + L2;
    //Reconstructing Higgs 4 vector with MET
    double nunu_px = Pmet_pt*cos(Pmet_phi);
    double nunu_py = Pmet_pt*sin(Pmet_phi);
    double nunu_pz = LL.Pz();
    double nunu_m = 30.0; //Why 30? --> https://indico.cern.ch/event/850505/contributions/3593915/

    double nunu_e = sqrt(nunu_px*nunu_px + nunu_py*nunu_py + nunu_pz*nunu_pz + nunu_m*nunu_m);
    NuNu.SetPxPyPzE(nunu_px, nunu_py, nunu_pz, nunu_e);
    Higgs = LL + NuNu;
    //double hm = Higgs.M();

    int jetn = 0;
    bool use3jet = false;
    for (unsigned int ijet = 0; ijet<ncleanjet; ijet++){

      if (CleanJet_pt->At(ijet)>30){
	++jetn;
	if (jetn==1) J1.SetPtEtaPhiM(CleanJet_pt->At(0), CleanJet_eta->At(0), CleanJet_phi->At(0), 0.0);
	if (jetn==2) J2.SetPtEtaPhiM(CleanJet_pt->At(1), CleanJet_eta->At(1), CleanJet_phi->At(1), 0.0);
	//if (jetn==3) J3.SetPtEtaPhiM(CleanJet_pt->At(2), CleanJet_eta->At(2), CleanJet_phi->At(2), 0.0);
	//	if (jetn==3) use3jet = true;
      }

    }

    //if (jetn < 2) return -9999;


    //J1.SetPtEtaPhiM(CleanJet_pt->At(0), CleanJet_eta->At(0), CleanJet_phi->At(0), 0.0);
    //J2.SetPtEtaPhiM(CleanJet_pt->At(1), CleanJet_eta->At(1), CleanJet_phi->At(1), 0.0);

    SimpleParticleCollection_t daughter;
    SimpleParticleCollection_t associated;
    SimpleParticleCollection_t mother;

    daughter.push_back(SimpleParticle_t(25, Higgs)); 

    //daughter.push_back(SimpleParticle_t(13, L1));
    //daughter.push_back(SimpleParticle_t(11, L2));
    //daughter.emplace_back(12, NuNu);
    associated.push_back(SimpleParticle_t(0,J1));
    associated.push_back(SimpleParticle_t(0,J2));

    //if (use3jet){
      //associated.push_back(SimpleParticle_t(0,J3));
      //}


    if (Higgs.Pt() == 0 || Higgs.M()==0){
      return -9999;
    }

    //MELA MATRIX ELEMENTS CALCULATION (STEP-2)
    //mela->setCandidateDecayMode(TVar::CandidateDecay_Stable);
    mela->setCandidateDecayMode(TVar::CandidateDecay_WW);
    mela->setInputEvent(&daughter, &associated, 0, false);
    //mela->setInputEvent(&daughter_coll, &associated_coll, 0, 0);
    mela->setCurrentCandidateFromIndex(0);

    //->VBF Processes
    float RecoLevel_me_VBF_hsm = 0.;
    float RecoLevel_me_VBF_hm = 0.;
    float RecoLevel_me_VBF_hp = 0.;
    float RecoLevel_me_VBF_hl = 0.;
    float RecoLevel_me_VBF_HSMHiggs_WWZZ = 0.; //Interference term
    float RecoLevel_me_VBF_hMix = 0.;

    float Q2V1 = 0.;
    float Q2V2 = 0.;

    float costheta1 = 0.;
    float costheta2 = 0.;
    float costhetastar = 0.;
    
    float phi = 0.;
    float phi1 = 0.;

    float D_CP = 0.;
    float D_intH = 0.;
    float D_intL = 0.;

    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
    mela->computeVBFAngles(Q2V1, Q2V2, costheta1, costheta2, phi, costhetastar, phi1);
    mela->computeProdP(RecoLevel_me_VBF_hsm, true);
    MatrixElementsMap.insert({"Q2V1", Q2V1});
    MatrixElementsMap.insert({"Q2V2", Q2V2});
    MatrixElementsMap.insert({"costheta1", costheta1});
    MatrixElementsMap.insert({"costheta2", costheta2});
    MatrixElementsMap.insert({"costhetastar", costhetastar});
    MatrixElementsMap.insert({"phi", phi});
    MatrixElementsMap.insert({"phi1", phi1});
    MatrixElementsMap.insert({"RecoLevel_me_VBF_hsm", RecoLevel_me_VBF_hsm});
    
    mela->setProcess(TVar::H0minus, TVar::JHUGen, TVar::JJVBF);
    mela->computeProdP(RecoLevel_me_VBF_hm, true);
    MatrixElementsMap.insert({"RecoLevel_me_VBF_hm", RecoLevel_me_VBF_hm});

    mela->setProcess(TVar::H0hplus, TVar::JHUGen, TVar::JJVBF);
    mela->computeProdP(RecoLevel_me_VBF_hp, true);
    MatrixElementsMap.insert({"RecoLevel_me_VBF_hp", RecoLevel_me_VBF_hp});

    mela->setProcess(TVar::H0_g1prime2, TVar::JHUGen, TVar::JJVBF);
    mela->computeProdP(RecoLevel_me_VBF_hl, true);
    MatrixElementsMap.insert({"RecoLevel_me_VBF_hl", RecoLevel_me_VBF_hl});


    // CP mixing (maximum)
    //mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJVBF);
    //mela->selfDHzzcoupl[0][gHIGGS_VV_1][0] = 1; // a1
    //mela->selfDHzzcoupl[0][gHIGGS_VV_4][0] = 0.297979; // a3
    //mela->computeProdP(RecoLevel_me_VBF_hMix, false);


    //  ------- D_* anomalous couplings -------

    
    float ME_sm       = 0.;
    float ME_bsm      = 0.;
    float ME_bsm_mlt  = 0.;
    float KD_bsm      = 0.;
    float melaD0minus = 0.;
    float ME_int= 0.;
    float melaDCP= 0.;
    float melaD0hplus = 0.;
    float melaDint = 0.;
    float melaDL1 = 0.;
    float melaDL1int = 0.;
    float melaDL1Zg = 0.;
    float melaDL1Zgint = 0.;

    //get ME scalar
    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
    mela->computeProdP(ME_sm, false);

    //get ME bsm
    mela->setProcess(TVar::H0minus, TVar::JHUGen, TVar::JJVBF);
    mela->computeProdP(ME_bsm, false);
    ME_bsm_mlt = ME_bsm*pow(0.297979, 2);

    //compute D_BSM (eq.5 of HIG-17-011)
    KD_bsm = ME_sm / (ME_sm + ME_bsm);
    melaD0minus = ME_sm / (ME_sm + ME_bsm_mlt);

    MatrixElementsMap.insert({"KD_bsm", KD_bsm});
    MatrixElementsMap.insert({"melaD0minus", melaD0minus});

    mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJVBF);
    mela->selfDHzzcoupl[0][gHIGGS_VV_4][0]=1;
    mela->selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela->computeProdP(ME_int, false);

    //define D_CP
    melaDCP = (0.297979*(ME_int-(ME_sm + ME_bsm)))/(ME_sm + (pow(0.297979, 2)*ME_bsm));

    MatrixElementsMap.insert({"melaDCP", melaDCP});

 
    mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJVBF);
    mela->selfDHzzcoupl[0][gHIGGS_VV_2][0]=1;
    mela->computeProdP(ME_bsm, false);
    ME_bsm_mlt = ME_bsm*pow(0.271899, 2);
    melaD0hplus = ME_sm / (ME_sm + ME_bsm_mlt);

    MatrixElementsMap.insert({"melaD0hplus", melaD0hplus});

    mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJVBF);
    mela->selfDHzzcoupl[0][gHIGGS_VV_2][0]=1;
    mela->selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela->computeProdP(ME_int, false);
    melaDint = (0.271899*(ME_int-(ME_sm + ME_bsm)))/(ME_sm + (pow(0.271899, 2)*ME_bsm));

    MatrixElementsMap.insert({"melaDint", melaDint});

    mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJVBF);
    mela->selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=1;
    mela->computeProdP(ME_bsm, false);
    ME_bsm_mlt = ME_bsm*pow(2156.43, 2);
    melaDL1 = ME_sm / (ME_sm + ME_bsm_mlt);

    MatrixElementsMap.insert({"melaDL1", melaDL1});

    mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJVBF);
    mela->selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=1;
    mela->selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela->computeProdP(ME_int, false);
    melaDL1int = (2156.43*(ME_int-(ME_sm + ME_bsm)))/(ME_sm + (pow(2156.43, 2)*ME_bsm));

    MatrixElementsMap.insert({"melaDL1int", melaDL1int});

    mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJVBF);
    mela->selfDHzzcoupl[0][gHIGGS_ZA_1_PRIME2][0]=1;
    mela->computeProdP(ME_bsm, false);
    ME_bsm_mlt = ME_bsm*pow(4091.0, 2);
    melaDL1Zg = ME_sm / (ME_sm + ME_bsm_mlt);

    MatrixElementsMap.insert({"melaDL1Zg", melaDL1Zg});

    mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJVBF);
    mela->selfDHzzcoupl[0][gHIGGS_ZA_1_PRIME2][0]=1;
    mela->selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela->computeProdP(ME_int, false);
    melaDL1Zgint = (4091.0*(ME_int-(ME_sm + ME_bsm)))/(ME_sm + (pow(4091.0, 2)*ME_bsm));

    MatrixElementsMap.insert({"melaDL1Zgint", melaDL1Zgint});


    //->QCD Processes                                                                                                                                                                                      
    float RecoLevel_me_QCD_hsm = 0.;
    float RecoLevel_me_QCD_hm = 0.;
    float RecoLevel_me_QCD_hp = 0.;
    float RecoLevel_me_QCD_hl = 0.;

    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJQCD);
    mela->computeProdP(RecoLevel_me_QCD_hsm, true);
    MatrixElementsMap.insert({"RecoLevel_me_QCD_hsm", RecoLevel_me_QCD_hsm});

    mela->setProcess(TVar::H0minus, TVar::JHUGen, TVar::JJQCD);
    mela->computeProdP(RecoLevel_me_QCD_hm, true);
    MatrixElementsMap.insert({"RecoLevel_me_QCD_hm", RecoLevel_me_QCD_hm});

    mela->setProcess(TVar::H0hplus, TVar::JHUGen, TVar::JJQCD);
    mela->computeProdP(RecoLevel_me_QCD_hp, true);
    MatrixElementsMap.insert({"RecoLevel_me_QCD_hp", RecoLevel_me_QCD_hp});

    mela->setProcess(TVar::H0_g1prime2, TVar::JHUGen, TVar::JJQCD);
    mela->computeProdP(RecoLevel_me_QCD_hl, true);
    MatrixElementsMap.insert({"RecoLevel_me_QCD_hl", RecoLevel_me_QCD_hl});

    
    //->VH Process

    float RecoLevel_me_VH_hsm = 0.;
    float RecoLevel_me_Had_ZH_hsm = 0.;
    float RecoLevel_me_Had_WH_hsm = 0.;


    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_ZH);
    mela->computeProdP(RecoLevel_me_Had_ZH_hsm, true);
    MatrixElementsMap.insert({"RecoLevel_me_Had_ZH_hsm", RecoLevel_me_Had_ZH_hsm});
    
    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_WH);
    mela->computeProdP(RecoLevel_me_Had_WH_hsm, true);
    MatrixElementsMap.insert({"RecoLevel_me_Had_WH_hsm", RecoLevel_me_Had_WH_hsm});

    RecoLevel_me_VH_hsm = RecoLevel_me_Had_ZH_hsm + RecoLevel_me_Had_WH_hsm;

    MatrixElementsMap.insert({"RecoLevel_me_VH_hsm", RecoLevel_me_VH_hsm});

    // Finish

    mela->resetInputEvent(); 
    
    float required_matrixelement = MatrixElementsMap.find(name_)->second;

    return (double)required_matrixelement;
    
  }
  else return -9999;
}
void
RecoLevelME::bindTree_(multidraw::FunctionLibrary& _library)
{
  //CleanJets
  _library.bindBranch(nCleanJet, "nCleanJet");
  _library.bindBranch(CleanJet_pt, "CleanJet_pt");
  _library.bindBranch(CleanJet_eta, "CleanJet_eta");
  _library.bindBranch(CleanJet_phi, "CleanJet_phi");
  //_library.bindBranch(CleanJet_mass, "CleanJet_mass");
  //Leptons
  _library.bindBranch(Lepton_pdgId, "Lepton_pdgId");
  _library.bindBranch(nLepton, "nLepton");
  _library.bindBranch(Lepton_pt, "Lepton_pt");
  _library.bindBranch(Lepton_eta, "Lepton_eta");
  _library.bindBranch(Lepton_phi, "Lepton_phi");
  //MET
  _library.bindBranch(MET_pt, "MET_pt");
  _library.bindBranch(PuppiMET_pt, "PuppiMET_pt");
  _library.bindBranch(PuppiMET_phi, "PuppiMET_phi");
  //Subjets
  _library.bindBranch(nSubJet, "nSubJet");
  _library.bindBranch(SubJet_pt, "SubJet_pt");
  _library.bindBranch(SubJet_eta, "SubJet_eta");
  _library.bindBranch(SubJet_phi, "SubJet_phi");
  _library.bindBranch(SubJet_mass, "SubJet_mass");
}
