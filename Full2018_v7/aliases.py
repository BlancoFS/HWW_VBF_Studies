import os
import copy
import inspect

configurations = os.path.realpath(inspect.getfile(inspect.currentframe())) # this file
configurations = os.path.dirname(configurations) # Snowmass
configurations = os.path.dirname(configurations) # VBF
configurations = os.path.dirname(configurations) # Configurations

#aliases = {}

# imported from samples.py:
# samples, signals

mc = [skey for skey in samples if skey not in ('Fake', 'DATA', 'Dyemb')]
mc_emb = [skey for skey in samples if skey not in ('Fake', 'DATA')]

eleWP='mvaFall17V1Iso_WP90'
muWP = 'cut_Tight_HWWW_tthmva_80'

aliases['LepWPCut'] = {
    'expr': 'LepCut2l__ele_'+eleWP+'__mu_'+muWP,
    'samples': mc_emb + ['DATA']
}

aliases['gstarLow'] = {
    'expr': 'Gen_ZGstar_mass >0 && Gen_ZGstar_mass < 4',
    'samples': 'VgS'
}

aliases['gstarHigh'] = {
    'expr': 'Gen_ZGstar_mass <0 || Gen_ZGstar_mass > 4',
    'samples': 'VgS'
}

aliases['embedtotal'] = {
    'expr': 'embed_total_WP90V1',  # wrt. eleWP
    'samples': 'Dyemb'
}

# Fake leptons transfer factor 
aliases['fakeW'] = {
    'expr': 'fakeW2l_ele_'+eleWP+'_mu_'+muWP,
    'samples': ['Fake']
}
# And variations - already divided by central values in formulas !
aliases['fakeWEleUp'] = {
    'expr': 'fakeW2l_ele_'+eleWP+'_mu_'+muWP+'_EleUp',
    'samples': ['Fake']
}
aliases['fakeWEleDown'] = {
    'expr': 'fakeW2l_ele_'+eleWP+'_mu_'+muWP+'_EleDown',
    'samples': ['Fake']
}
aliases['fakeWMuUp'] = {
    'expr': 'fakeW2l_ele_'+eleWP+'_mu_'+muWP+'_MuUp',
    'samples': ['Fake']
}
aliases['fakeWMuDown'] = {
    'expr': 'fakeW2l_ele_'+eleWP+'_mu_'+muWP+'_MuDown',
    'samples': ['Fake']
}
aliases['fakeWStatEleUp'] = {
    'expr': 'fakeW2l_ele_'+eleWP+'_mu_'+muWP+'_statEleUp',
    'samples': ['Fake']
}
aliases['fakeWStatEleDown'] = {
    'expr': 'fakeW2l_ele_'+eleWP+'_mu_'+muWP+'_statEleDown',
    'samples': ['Fake']
}
aliases['fakeWStatMuUp'] = {
    'expr': 'fakeW2l_ele_'+eleWP+'_mu_'+muWP+'_statMuUp',
    'samples': ['Fake']
}
aliases['fakeWStatMuDown'] = {
    'expr': 'fakeW2l_ele_'+eleWP+'_mu_'+muWP+'_statMuDown',
    'samples': ['Fake']
}

# gen-matching to prompt only (GenLepMatch2l matches to *any* gen lepton)
aliases['PromptGenLepMatch2l'] = {
    'expr': 'Alt$(Lepton_promptgenmatched[0]*Lepton_promptgenmatched[1], 0)',
    'samples': mc
}

aliases['Top_pTrw'] = {
    'expr': '(topGenPt * antitopGenPt > 0.) * (TMath::Sqrt((0.103*TMath::Exp(-0.0118*topGenPt) - 0.000134*topGenPt + 0.973) * (0.103*TMath::Exp(-0.0118*antitopGenPt) - 0.000134*antitopGenPt + 0.973))) + (topGenPt * antitopGenPt <= 0.)',
    'samples': ['top']
}

aliases['nCleanGenJet'] = {
    'linesToAdd': ['.L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/Differential/ngenjet.cc+'],
    'class': 'CountGenJet',
    'samples': mc
}

##### DY Z pT reweighting
aliases['getGenZpt_OTF'] = {
    'linesToAdd':['.L %s/src/PlotsConfigurations/Configurations/patches/getGenZpt.cc+' % os.getenv('CMSSW_BASE')],
    'class': 'getGenZpt',
    'samples': ['DY']
}
handle = open('%s/src/PlotsConfigurations/Configurations/patches/DYrew30.py' % os.getenv('CMSSW_BASE'),'r')
exec(handle)
handle.close()
aliases['DY_NLO_pTllrw'] = {
    'expr': '('+DYrew['2018']['NLO'].replace('x', 'getGenZpt_OTF')+')*(nCleanGenJet == 0)+1.0*(nCleanGenJet > 0)',
    'samples': ['DY']
}
aliases['DY_LO_pTllrw'] = {
    'expr': '('+DYrew['2018']['LO'].replace('x', 'getGenZpt_OTF')+')*(nCleanGenJet == 0)+1.0*(nCleanGenJet > 0)',
    'samples': ['DY']
}

# Jet bins
# using Alt$(CleanJet_pt[n], 0) instead of Sum$(CleanJet_pt >= 30) because jet pt ordering is not strictly followed in JES-varied samples

# No jet with pt > 30 GeV
aliases['zeroJet'] = {
    'expr': 'Alt$(CleanJet_pt[0], 0) < 30.'
}

aliases['oneJet'] = {
    'expr': 'Alt$(CleanJet_pt[0], 0) > 30.'
}

aliases['multiJet'] = {
    'expr': 'Alt$(CleanJet_pt[1], 0) > 30.'
}

# B tagging

aliases['bVeto'] = {
    'expr': 'Sum$(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepFlavB[CleanJet_jetIdx] >  0.0494) == 0'
}

aliases['bReq'] = {
    'expr': 'Sum$(CleanJet_pt > 30. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepFlavB[CleanJet_jetIdx] >  0.0494) >= 1'
}

#aliases['topcr'] = {
#    'expr': 'mtw2>30 && mll>50 && ((zeroJet && !bVeto) || bReq)'
#}

aliases['topcr'] = {
    'expr': 'Sum$(CleanJet_pt>30) >= 2 && mpmet>20 && PuppiMET_pt>20 && mll>20 && ptll>30 && ((zeroJet && !bVeto) || bReq)'
}


btagSFSource = '%s/src/PhysicsTools/NanoAODTools/data/btagSF/DeepJet_102XSF_V1.csv' % os.getenv('CMSSW_BASE')

aliases['Jet_btagSF_deepflav_shape'] = {
    'linesToAdd': [
        'gSystem->Load("libCondFormatsBTauObjects.so");',
        'gSystem->Load("libCondToolsBTau.so");',
        'gSystem->AddIncludePath("-I%s/src");' % os.getenv('CMSSW_RELEASE_BASE'),
        '.L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/patches/btagsfpatch.cc+'
    ],
    'class': 'BtagSF',
    'args': (btagSFSource,'central','deepjet'),
    'samples': mc
}


aliases['bVetoSF'] = {
    'expr': 'TMath::Exp(Sum$(TMath::Log((CleanJet_pt>20 && abs(CleanJet_eta)<2.5)*Jet_btagSF_deepflav_shape[CleanJet_jetIdx]+1*(CleanJet_pt<20 || abs(CleanJet_eta)>2.5))))',
    'samples': mc
}

aliases['bReqSF'] = {
    'expr': 'TMath::Exp(Sum$(TMath::Log((CleanJet_pt>30 && abs(CleanJet_eta)<2.5)*Jet_btagSF_deepflav_shape[CleanJet_jetIdx]+1*(CleanJet_pt<30 || abs(CleanJet_eta)>2.5))))',
    'samples': mc
}

aliases['btagSF'] = {
    'expr': '(bVeto || (topcr && zeroJet))*bVetoSF + (topcr && !zeroJet)*bReqSF',
    'samples': mc
}
#hf
for shift in ['jes', 'lf', 'lfstats1', 'lfstats2', 'hfstats1', 'hfstats2', 'cferr1', 'cferr2']:
    aliases['Jet_btagSF_deepflav_shape_up_%s' % shift] = {
        'class': 'BtagSF',
        'args': (btagSFSource, 'up_' + shift,'deepjet'),
        'samples': mc
    }
    aliases['Jet_btagSF_deepflav_shape_down_%s' % shift] = {
        'class': 'BtagSF',
        'args': (btagSFSource, 'down_' + shift,'deepjet'),
        'samples': mc
    }
    
    for targ in ['bVeto', 'bReq']:
        alias = aliases['%sSF%sup' % (targ, shift)] = copy.deepcopy(aliases['%sSF' % targ])
        alias['expr'] = alias['expr'].replace('btagSF_deepflav_shape', 'btagSF_deepflav_shape_up_%s' % shift)
        
        alias = aliases['%sSF%sdown' % (targ, shift)] = copy.deepcopy(aliases['%sSF' % targ])
        alias['expr'] = alias['expr'].replace('btagSF_deepflav_shape', 'btagSF_deepflav_shape_down_%s' % shift)
        
    aliases['btagSF%sup' % shift] = {
        'expr': aliases['btagSF']['expr'].replace('SF', 'SF' + shift + 'up'),
        'samples': mc
    }
    
    aliases['btagSF%sdown' % shift] = {
        'expr': aliases['btagSF']['expr'].replace('SF', 'SF' + shift + 'down'),
        'samples': mc
    }
    


aliases['Jet_PUIDSF'] = {
  'expr' : 'TMath::Exp(Sum$((Jet_jetId>=2)*TMath::Log(Jet_PUIDSF_loose)))',
  'samples': mc
}

aliases['Jet_PUIDSF_up'] = {
  'expr' : 'TMath::Exp(Sum$((Jet_jetId>=2)*TMath::Log(Jet_PUIDSF_loose_up)))',
  'samples': mc
}

aliases['Jet_PUIDSF_down'] = {
  'expr' : 'TMath::Exp(Sum$((Jet_jetId>=2)*TMath::Log(Jet_PUIDSF_loose_down)))',
  'samples': mc
}

'''
# data/MC scale factors
aliases['SFweight'] = {
    'expr': ' * '.join(['SFweight2l', 'LepWPCut', 'LepSF2l__ele_' + eleWP + '__mu_' + muWP, 'btagSF', 'Jet_PUIDSF_loose']),
    'samples': mc
}
'''

# data/MC scale factors
aliases['SFweight'] = {
    'expr': ' * '.join(['SFweight2l', 'LepWPCut', 'LepSF2l__ele_' + eleWP + '__mu_' + muWP, 'btagSF', 'Jet_PUIDSF']),
    'samples': mc
}

# Muon ttHMVA SF needed for tau embedded samples
aliases['Muon_ttHMVA_SF'] = {
    'expr': '( (abs(Lepton_pdgId[0]) == 13)*(Lepton_tightMuon_cut_Tight_HWWW_tthmva_80_IdIsoSF[0]/Lepton_tightMuon_cut_Tight_HWWW_IdIsoSF[0])+(abs(Lepton_pdgId[0]) == 11) )*( (abs(Lepton_pdgId[1]) == 13)*(Lepton_tightMuon_cut_Tight_HWWW_tthmva_80_IdIsoSF[1]/Lepton_tightMuon_cut_Tight_HWWW_IdIsoSF[1])+ (abs(Lepton_pdgId[1]) == 11) )',
    'samples' : ['Dyemb']
}

# variations
aliases['SFweightEleUp'] = {
    'expr': 'LepSF2l__ele_'+eleWP+'__Up',
    'samples': mc_emb
}
aliases['SFweightEleDown'] = {
    'expr': 'LepSF2l__ele_'+eleWP+'__Do',
    'samples': mc_emb
}
aliases['SFweightMuUp'] = {
    'expr': 'LepSF2l__mu_'+muWP+'__Up',
    'samples': mc_emb
}
aliases['SFweightMuDown'] = {
    'expr': 'LepSF2l__mu_'+muWP+'__Do',
    'samples': mc_emb
}


aliases['Weight2MINLO'] = {
    'linesToAdd': ['.L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/Differential/weight2MINLO.cc+'],
    'class': 'Weight2MINLO',
    'args': '%s/src/LatinoAnalysis/Gardener/python/data/powheg2minlo/NNLOPS_reweight.root' % os.getenv('CMSSW_BASE'),
    'samples' : [skey for skey in samples if 'ggH_hww' in skey],
}

## GGHUncertaintyProducer wasn't run for GluGluHToWWTo2L2Nu_Powheg_M125 
thus = [
    'ggH_mu',
    'ggH_res',
    'ggH_mig01',
    'ggH_mig12',
    'ggH_VBF2j',
    'ggH_VBF3j',
    'ggH_pT60',
    'ggH_pT120',
    'ggH_qmtop'
]

for thu in thus:
    aliases[thu+'_2'] = {
        'linesToAdd': ['.L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/Differential/gghuncertainty.cc+'],
        'class': 'GGHUncertainty',
        'args': (thu,),
        'samples': ['ggH_hww']
    }


#### VARIABLES #####

aliases['mlljj'] = {
    'linesToAdd' : ['.L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/WW/Full2016_v7/extended/mlljj.cc+' ],
    'class' : 'mlljj'
}

aliases['mjjj'] = {
    'linesToAdd' : ['.L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/WW/Full2016_v7/extended/mjjj.cc+' ],
    'class' : 'mjjj'
}

aliases['jetdis'] = {
    'linesToAdd' : ['.L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/WW/Full2016_v7/extended/GetJetDis.cc+'],
    'class' : 'GetJetDis'
}



#### MELA FRAMEWORK #####

mes_reduced = ['RecoLevel_me_VBF_hsm', 'RecoLevel_me_QCD_hsm']

for me in mes_reduced:
    aliases[me]={
    'linesToAdd': [
        'gSystem->Load("%s/src/JHUGenMELA/MELA/data/%s/libmcfm_707.so","", kTRUE);'%(os.getenv('CMSSW_BASE'), os.getenv('SCRAM_ARCH')),
        'gSystem->Load("libJHUGenMELAMELA.so","", kTRUE);',
        '.L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/patches/RecoLevelME_patch.cc+'],
    'class': 'RecoLevelME',
    'args': (me,)
    }



aliases['D_VBF_QCD']     = { 'expr': '(RecoLevel_me_VBF_hsm**2)/(RecoLevel_me_VBF_hsm**2 + RecoLevel_me_QCD_hsm**2)' }
#aliases['D_VBF_VH']     = { 'expr': '(RecoLevel_me_VBF_hsm**2)/(RecoLevel_me_VBF_hsm**2 + 1e16 * (RecoLevel_me_VH_hsm**2))' }                                                                                                        
#aliases['D_VH_QCD']     = { 'expr': '(RecoLevel_me_VH_hsm**2)/(RecoLevel_me_VH_hsm**2 + 1e-16 * (RecoLevel_me_QCD_hsm**2))' }  


'''
###                                                                                                                                                                                                                                          
###                                                                                                                                                                                                                                          
### MATRIX ELEMENT METHOD WITH MoMEMta                                                                                                                                                                                                       
###                                                                                                                                                                                                                                          
###                                                                                                                                                                                                                                          


aliases['VBF_ME']={                                                                                                                                                                                                                          
    'linesToAdd': [                                                                                                                                                                                                                          
        'gSystem->Load("/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/lib/libmomemta.so.1.0.0","", kTRUE);',                                                                                                                             
        'gSystem->Load("/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/lib/libmomemta.so.1","", kTRUE);',                                                                                                                                 
        'gSystem->Load("/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/lib/libmomemta.so","", kTRUE);',                                                                                                                                   
        'gSystem->Load("/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/lib/libLHAPDF.so", "", kTRUE);',                                                                                                                                   
        '.L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/patches/RecoME_MoMEMta.cc+'],                                                                                                          
    'class': 'RecoME',                                                                                                                                                                                                                       
    'args': 'VBF'                                                                                                                                                                                                                            
}                                                                                                                                                                                                                                            

aliases['DY_ME']={                                                                                                                                                                                                                           
    'linesToAdd': [                                                                                                                                                                                                                          
        'gSystem->Load("/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/lib/libmomemta.so.1.0.0","", kTRUE);',                                                                                                                             
        'gSystem->Load("/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/lib/libmomemta.so.1","", kTRUE);',                                                                                                                                 
        'gSystem->Load("/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/lib/libmomemta.so","", kTRUE);',                                                                                                                                   
        'gSystem->Load("/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/lib/libLHAPDF.so", "", kTRUE);',                                                                                                                                   
        '.L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/patches/RecoME_MoMEMta.cc+'],                                                                                                 
    'class': 'RecoME',                                                                                                                                                                                                                       
    'args': 'DY'                                                                                                                                                                                                                             
}                                                                                                                                                                                                                                           
 
aliases['D_DY']     = { 'expr': '150 * abs(VBF_ME)/(150 * abs(VBF_ME) + abs(DY_ME))' } 
'''

######## 

aliases['BDTG4D3'] = {
    'linesToAdd': ['gSystem->Load("%s/src/JHUGenMELA/MELA/data/%s/libmcfm_707.so","", kTRUE);'%(os.getenv('CMSSW_BASE'), os.getenv('SCRAM_ARCH')),
                   'gSystem->Load("libJHUGenMELAMELA.so","", kTRUE);',
                   'gSystem->Load("/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/lib/libmomemta.so","", kTRUE);',
                   '.L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/src/PlotsConfigurations/Configurations/WW/Full2016_v7/DNN/TMVA_BDT_HWW.cc+' ],
    'class': 'TMVA_HWW',
    'args': 'BDTG4D3'
}
