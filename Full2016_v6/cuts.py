 # VH2j cuts


#-------------------------------------------------------------------------------
# supercut
#-------------------------------------------------------------------------------
_tmp = [
     'Lepton_pdgId[0]*Lepton_pdgId[1] == -11*13',
     'Lepton_pt[0] > 25.',
     'Lepton_pt[1] > 13.',
     '(abs(Lepton_pdgId[1]) == 13 || Lepton_pt[1] > 13.)',
     '(nLepton >= 2 && Alt$(Lepton_pt[2], 0) < 10.)',
     ]

supercut = ' && '.join(_tmp)


def addcut(name, exprs):
    cuts[name] = ' && '.join(exprs)


#-------------------------------------------------------------------------------                                                                                                                           
# VBF_Higgs                                                                                                                                                                                                    
#-------------------------------------------------------------------------------                                                                                                                           
_tmp = [
     'Alt$(CleanJet_pt[1], 0) > 20.',
     'abs(CleanJet_eta[0]) < 2.5',
     'abs(CleanJet_eta[1]) < 2.5',
     'mth < 115',
     'drll < 2.',                                                                                                                                                                                 
     'Jet_qgl[0] < 0.6',
     'Jet_qgl[1] < 0.6',
     'bVeto'
     ]

addcut('Higgs', _tmp)


#-------------------------------------------------------------------------------
# Test
#-------------------------------------------------------------------------------
### _tmp = [
###      'VH2j_TMVAReader(Entry$) > 0.1',
###      ]
### 
### addcut('VH_2j_test', _tmp);