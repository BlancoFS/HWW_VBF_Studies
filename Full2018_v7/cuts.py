#VBF cuts                                                                                                                                                                           

#-------------------------------------------------------------------------------                                                                                                    
# supercut                                                                                                                                                                          
#-------------------------------------------------------------------------------                                                                                                    
_tmp = [
     'Lepton_pdgId[0]*Lepton_pdgId[1] == -11*13',
     'Lepton_pt[0] > 25.',
     'Lepton_pt[1] > 15.',
     '(abs(Lepton_pdgId[1]) == 13 || Lepton_pt[1] > 13.)',
     '(nLepton >= 2 && Alt$(Lepton_pt[2], 0) < 10.)'
     ]

supercut = ' && '.join(_tmp)


def addcut(name, exprs):
    cuts[name] = ' && '.join(exprs)


#-------------------------------------------------------------------------------                                                                                                    
# VBF_Bkg No cuts                                                                                                                                                                   
#-------------------------------------------------------------------------------                                                                                                    
_tmp = [
     'Alt$(CleanJet_pt[1], 0) > 30.',
     'abs(CleanJet_eta[0]) < 4.5',
     'abs(CleanJet_eta[1]) < 4.5',
     'Sum$(CleanJet_pt>30) >= 2',
     #'Sum$(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.2217) == 0'
     ]

addcut('no_cuts', _tmp)


#-------------------------------------------------------------------------------                                                                                                                                                             
# VBF REGION                                                                                                                                                                                                                                 
#-------------------------------------------------------------------------------                                                                                                                                                             
_tmp = [
    'Alt$(CleanJet_pt[1], 0) > 30.',
    'abs(CleanJet_eta[0]) < 4.5',
    'abs(CleanJet_eta[1]) < 4.5',
    'Sum$(CleanJet_pt>30) >= 2',
    'Sum$(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepFlavB[CleanJet_jetIdx] > 0.0614) == 0',
    'mth > 60 && mth < 150',
    'PuppiMET_pt > 10',
    'mll > 12',
    'D_VBF_QCD > 0.1'
]

addcut('VBF', _tmp)

#-------------------------------------------------------------------------------                                                                                                                           
# VBF REGION                                                                                                                                                                                               
#-------------------------------------------------------------------------------                                                                                                                           
_tmp = [
    'Alt$(CleanJet_pt[1], 0) > 30.',
    'abs(CleanJet_eta[0]) < 4.5',
    'abs(CleanJet_eta[1]) < 4.5',
    'Sum$(CleanJet_pt>30) >= 2',
    'mth > 60',
    'D_VBF_QCD > 0.2'
]

addcut('VBF_min', _tmp)

#---------
# GGH
#----------

_tmp = [
    'Alt$(CleanJet_pt[1], 0) > 30.',
    'abs(CleanJet_eta[0]) < 4.5',
    'abs(CleanJet_eta[1]) < 4.5',
    'Sum$(CleanJet_pt>30) >= 2',
    'Sum$(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepFlavB[CleanJet_jetIdx] > 0.0614) == 0',
    'mth > 60 && mth < 150',
    'PuppiMET_pt > 10',
    'mll > 12',
    'D_VBF_QCD < 0.1'
]

addcut('GGH', _tmp)

#-------------------------------------------------------------------------------                                                                                                                           
# VBF REGION                                                                                                                                                                                               
#-------------------------------------------------------------------------------                                                                                                                           
_tmp = [
    'Alt$(CleanJet_pt[1], 0) > 30.',
    'abs(CleanJet_eta[0]) < 4.5',
    'abs(CleanJet_eta[1]) < 4.5',
    'Sum$(CleanJet_pt>30) >= 2',
    'Sum$(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepFlavB[CleanJet_jetIdx] > 0.0614) == 0',
    'mth > 60 && mth < 150',
    'PuppiMET_pt > 10',
    'mll > 12'
]

addcut('VBF_minimal_flav', _tmp)

#-------------------------------------------------------------------------------                                                                                                                                                             
# VBF REGION                                                                                                                                                                                                                                 
#-------------------------------------------------------------------------------                                                                                                                                                             
_tmp = [
    'Alt$(CleanJet_pt[1], 0) > 30.',
    'abs(CleanJet_eta[0]) < 4.5',
    'abs(CleanJet_eta[1]) < 4.5',
    'Sum$(CleanJet_pt>30) >= 2',
    'mth > 60 && mth < 150',
    'PuppiMET_pt > 10',
    'mll > 12'
]

addcut('VBF_minimal_nobtag', _tmp)

#-------------------------------------------------------------------------------                                                                                                                                                            
# VBF_Bkg DY                                                                                                                                                                                                                                
#-------------------------------------------------------------------------------                                                                                                                                                             
_tmp = [
     'Alt$(CleanJet_pt[1], 0) > 30.',
     'abs(CleanJet_eta[0]) < 4.5',
     'abs(CleanJet_eta[1]) < 4.5',
     'mth < 40',
     'Sum$(CleanJet_pt>20) >= 2',
     'Sum$(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepFlavB[CleanJet_jetIdx] > 0.0614) == 0',
     'mll < 80',
     ]

addcut('DY', _tmp)


#-------------------------------------------------------------------------------                                                                                                                                                             
# VBF_Bkg top                                                                                                                                                                                                                                
#-------------------------------------------------------------------------------                                                                                                                                                             
_tmp = [
     'Alt$(CleanJet_pt[1], 0) > 30.',
     'abs(CleanJet_eta[0]) < 4.5',
     'abs(CleanJet_eta[1]) < 4.5',
     'Sum$(CleanJet_pt>30) >= 2',
     'Sum$(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepFlavB[CleanJet_jetIdx] > 0.0614) >= 1',
     'mll > 80'
     ]

addcut('top', _tmp)

_tmp = [
     'Alt$(CleanJet_pt[1], 0) > 30.',
     'abs(CleanJet_eta[0]) < 4.5',
     'abs(CleanJet_eta[1]) < 4.5',
     'Sum$(CleanJet_pt>30) >= 2',
     'topcr'
     ]

addcut('top_cr', _tmp)



#-------------------------------------------------------------------------------                                                                                                                           
# WW                                                                                                                                                                                                       
#-------------------------------------------------------------------------------                                                                                                                           
_tmp = [
     'abs(CleanJet_eta[0]) < 4.7',
     'abs(CleanJet_eta[1]) < 4.7',
     'Sum$(CleanJet_pt>30) == 2',
     'bVeto',
     'mpmet > 20',
     'PuppiMET_pt > 20',
     'mll > 20',
     'ptll > 30',
     'BDTG4D3 < 0.5'
]

addcut('WW', _tmp)
