# variables

#variables['bdt'] = {      'name'  : 'VH2j_TMVAReader(Entry$)',# variable name
#                          'range' : (20, -1., 1.),            # variable range
#                          'xaxis' : 'BDT discriminant VH2j',  # x-axis name
#                          'fold'  : 3,                        # 0 = not fold (default), 1 = fold underflow bin, 2 = fold overflow bin, 3 = fold underflow and overflow
#                          'linesToAdd' : ['.L /afs/cern.ch/user/p/piedra/work/VH2jBDT/VH2j_TMVAReader.C+']}



variables['BDT_VBF'] = { 'name'  : 'BDTG4D3',
                         'range' : (20, -1., 1.),
                         'xaxis' : 'BDT VBF',
                         'fold'  : 3}

variables['BDT_VBF_fine'] = { 'name'  : 'BDTG4D3',
                              'range' : ([-1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0],),
                              'xaxis' : 'BDT VBF',
                              'fold'  : 3}

variables['BDT_VBF_2'] = { 'name'  : 'BDTG4D3',
                           'range' : ([-1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.925, 1.0],),
                           'xaxis' : 'BDT VBF',
                           'fold'  : 3}

#variables['BDT_VBF_015'] = { 'name'  : 'BDTG4D3',
#                             'range' : (8, 0.15, 1.),
#                             'xaxis' : 'BDT VBF',
#                             'fold'  : 0}

#variables['BDT_VBF_05'] = { 'name'  : 'BDTG4D3',
#                            'range' : (5, 0.5, 1.),
#                            'xaxis' : 'BDT VBF',
#                            'fold'  : 0}



'''
variables['BDT_VBF_bin'] = { 'name'  : 'BDTG4D3',
                             'range' : ([0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 1.0],),
                             'xaxis' : 'BDT VBF',
                             'fold'  : 0}


variables['BDT_VBF_bin2'] = { 'name'  : 'BDTG4D3',
                              'range' : ([0.65, 0.75, 0.85, 0.9, 0.95, 1.0],),
                              'xaxis' : 'BDT VBF',
                              'fold'  : 0}


variables['BDT_VBF_test'] = { 'name'  : 'BDTG4D3',
                              'range' : ([0.65, 0.7, 0.75, 0.825, 0.875, 0.925, 1.0],),
                              'xaxis' : 'BDT VBF',
                              'fold'  : 0}
'''

variables['D_VBF_QCD'] = { 'name'  : 'D_VBF_QCD',                                                                                                                                                                                            
                           'range' : (30, 0., 1.),                                                                                                                                                                                           
                           'xaxis' : 'D_VBF_QCD',                                                                                                                                                                                            
                           'fold'  : 3}

variables['D_VBF_VH'] = { 'name'  : 'D_VBF_VH',                                                                                                                                                                                              
                          'range' : (30, 0., 1.),                                                                                                                                                                                           
                          'xaxis' : 'D_VBF_VH',                                                                                                                                                                                             
                          'fold'  : 3}                                                                                                                                                                                                       

variables['D_VH_QCD'] = { 'name'  : 'D_VH_QCD',                                                                                                                                                                                              
                          'range' : (30, 0., 1.),                                                                                                                                                                                            
                          'xaxis' : 'D_VH_QCD',                                                                                                                                                                                              
                          'fold'  : 3}   


variables['D_VBF_DY'] = { 'name'  : 'D_DY',                                                                                                                                                                                                  
                          'range' : (30, 0., 1.),                                                                                                                                                                                            
                          'xaxis' : 'D_DY',                                                                                                                                                                                                  
                          'fold'  : 3}

'''

variables['events'] = {   'name'  : '1',
                          'range' : (1, 0, 2),
                          'xaxis' : 'events',
                          'fold'  : 3} 


#
# Centrality
#
variables['Ctot'] = {     'name': 'log((abs(2*Lepton_eta[0]-CleanJet_eta[0]-CleanJet_eta[1])+abs(2*Lepton_eta[1]-CleanJet_eta[0]-CleanJet_eta[1]))/detajj)',
                          'range' : (20, -4, 6),
                          #'range' : (10, -4, 2),
                          'xaxis' : 'Ctot',
                          'fold'  : 3}

variables['jet1_qgl'] = { 'name'  : 'Jet_qgl[0]',
                          'range' : (20, 0., 1.),
                          'xaxis' : 'Quark Gluon likelihood (1^{st} Jet)',
                          'fold'  : 3}

variables['jet2_qgl'] = { 'name'  : 'Jet_qgl[1]',
                          'range' : (20, 0., 1.),
                          'xaxis' : 'Quark Gluon likelihood (2^{nd} Jet)',
                          'fold'  : 3}

variables['detajj'] = {   'name'  : 'detajj',
                          'range' : (20, 0., 5.),
                          #'range' : (10, 3., 4.),
                          'xaxis' : '|#Delta#eta_{jj}|',
                          'fold'  : 3}

variables['detal1j1'] = { 'name'  : 'abs(Lepton_eta[0]-CleanJet_eta[0])',
                          #'range' : (20, 0., 4.),
                          'range' : (10, 0., 4.),
                          'xaxis' : '|#Delta#eta_{l1j1}|',
                          'fold'  : 3}

variables['detal2j2'] = { 'name'  : 'abs(Lepton_eta[1]-CleanJet_eta[1])',
                          'range' : (20, 0., 4.),
                          'xaxis' : '|#Delta#eta_{l2j2}|',
                          'fold'  : 3}
 
variables['detal1j2'] = { 'name'  : 'abs(Lepton_eta[0]-CleanJet_eta[1])',
                          'range' : (20, 0., 4.),
                          'xaxis' : '|#Delta#eta_{l1j2}|',
                          'fold'  : 3}

variables['detal2j1'] = { 'name'  : 'abs(Lepton_eta[1]-CleanJet_eta[0])',
                          'range' : (20, 0., 4.),
                          'xaxis' : '|#Delta#eta_{l2j1}|',
                          'fold'  : 3}
         
variables['detall'] = { 'name'  : 'abs(Lepton_eta[0]-Lepton_eta[1])',
                          'range' : (20, 0., 5.),
                          'xaxis' : '|#Delta#eta_{l1l2}|',
                          'fold'  : 3}

variables['detall_01'] = { 'name'  : 'abs(Lepton_eta[0]-Lepton_eta[1])',
                          'range' : (10, 0., 1.),
                          'xaxis' : '|#Delta#eta_{l1l2}|',
                          'fold'  : 0}

variables['drjj'] = {     'name'  : 'sqrt((CleanJet_eta[0]-CleanJet_eta[1])**2 + (CleanJet_phi[0]-CleanJet_phi[1])**2)',
                          'range' : (20, 0., 8.),
                          #'range' : (10, 0., 5.),
                          'xaxis' : '#DeltaR_{jj}',
                         #'linesToAdd': ['.L $CMSSW_BASE/src/PlotsConfigurations/Configurations/WW/Full2016_v6/extended/drjj.C+'], #if want to use a script
                          'fold'  : 3}
          
          

#variables['detaljmin'] = {'name'  : 'detaljmin(Lepton_eta[0], Lepton_eta[1], CleanJet_eta[0], CleanJet_eta[1])',
#                          'range' : (40, -4., 4.),
#                          'xaxis' : 'min #Delta#eta_{lj}',
#                          'fold'  : 0,
#                          'linesToAdd' : ['.L /afs/cern.ch/work/p/piedra/public/latinos/CMSSW_10_2_15_patch2/src/PlotsConfigurations/Configurations/VH2j/Full2016_v6/detaljmin.C+']}

variables['dphijjmet'] = {'name'  : 'abs(dphijjmet)',     
                          'range' : (20, 0., 3.2),
                          #'range' : (10, 0., 3.2),
                          'xaxis' : '#Delta#phi_{jjmet}',
                          'fold'  : 3}

variables['dphill'] = {   'name'  : 'abs(dphill)',     
                          'range' : (20, 0., 3.2),   
                          #'range' : (10, 0., 3.2),
                          'xaxis' : '#Delta#phi_{ll}',
                          'fold'  : 3}

variables['dphijj'] = {   'name'  : 'CleanJet_phi[0] - CleanJet_phi[1]',
                          'range' : (20, -3.2, 3.2),
                          #'range' : (10, 0., 3.2),                                                                                                                                                        
                          'xaxis' : '#Delta#phi_{jj}',
                          'fold'  : 3}

variables['dphill_01'] = {   'name'  : 'abs(dphill)',
                          'range' : (10, 0., 1.),
                          'xaxis' : '#Delta#phi_{ll}',
                          'fold'  : 3}

variables['dphillj'] = {  'name'  : 'abs(dphilljet)',     
                          'range' : (20, 0., 3.2),   
                          #'range' : (10, 0., 3.2),
                          'xaxis' : '#Delta#phi_{llj}',
                          'fold'  : 3}

variables['dphilljj'] = { 'name'  : 'abs(dphilljetjet)',     
                          'range' : (20, 0., 3.2),   
                          #'range' : (10, 0., 3.2),
                          'xaxis' : '#Delta#phi_{lljj}',
                          'fold'  : 3}

variables['drll'] = {     'name'  : 'drll',
                          'range' : (30, 0., 5.),
                          #'range' : (10, 0., 2.5),
                          'xaxis' : '#DeltaR_{ll}',
                          'fold'  : 3}

variables['drll_01'] = {  'name'  : 'drll',
                          'range' : (10, 0., 1.),                                                                                                                                                         
                          'xaxis' : '#DeltaR_{ll}',
                          'fold'  : 2}

variables['eta1'] = {     'name'  : 'Lepton_eta[0]',     
                          'range' : (40, -3.2, 3.2),   
                          #'range' : (15, -3.2, 3.2),
                          'xaxis' : '#eta 1st lepton',
                          'fold'  : 3}

variables['eta2'] = {     'name'  : 'Lepton_eta[1]',     
                          'range' : (40, -3.2, 3.2),
                          #'range' : (15, -3.2, 3.2),
                          'xaxis' : '#eta 2nd lepton',
                          'fold'  : 3}

variables['events'] = {   'name'  : '1',
                          'range' : (1, 0, 2),
                          'xaxis' : 'events',
                          'fold'  : 3}

variables['jeteta1'] = {  'name'  : 'CleanJet_eta[0]',
                          #'range' : (80, -4., 4.),
                          'range' : (15, -4., 4.),
                          'xaxis' : '#eta 1st jet',
                          'fold'  : 3}

variables['jeteta2'] = {  'name'  : 'CleanJet_eta[1]',
                          #'range' : (80, -4., 4.),
                          'range' : (15, -4., 4.),
                          'xaxis' : '#eta 2nd jet',
                          'fold'  : 3}

variables['jeteta3'] = {  'name'  : 'CleanJet_eta[2]',
                          #'range' : (80, -4., 4.),
                          'range' : (15, -4., 4.),
                          'xaxis' : '#eta 3rd jet',
                          'fold'  : 3}

variables['jetphi1'] = {  'name'  : 'CleanJet_phi[0]',
                          'range' : (40, -3.2, 3.2),
                          'xaxis' : '#phi 1st jet',
                          'fold'  : 3}

variables['jetphi2'] = {  'name'  : 'CleanJet_phi[1]',
                          'range' : (40, -3.2, 3.2),
                          'xaxis' : '#phi 2nd jet',
                          'fold'  : 3}

variables['jetphi3'] = {  'name'  : 'CleanJet_phi[2]',
                          'range' : (40, -3.2, 3.2),
                          'xaxis' : '#phi 3rd jet',
                          'fold'  : 3}

variables['GetjetDis'] = {  #'name'  : 'abs(CleanJet_eta[2]-(CleanJet_eta[0]+CleanJet_eta[1])/2)*(CleanJet_pt[2]>30)',
                          'name'  : 'abs(jetdis)',
                          'range' : (10, 0.0, 4.5),
                          'xaxis' : '|#eta^{*}| = |#eta_{3} - (#eta_{1} + #eta_{2})/2|',
                          'fold'  : 0}

variables['jetpt1'] = {   'name'  : 'CleanJet_pt[0]*(CleanJet_pt[0]>30)',     
                          #'range' : (40, 30., 190.),
                          'range' : (10, 30., 190.),
                          'xaxis' : 'p_{T} 1st jet',
                          'fold'  : 3}

variables['jetpt1_0j'] = {'name'  : 'CleanJet_pt[0]*(CleanJet_pt[0]<30)',     
                          'range' : (30, 0., 30.),   
                          'xaxis' : 'p_{T} 1st jet (p_{T} < 30 GeV) ',
                          'fold'  : 0}

variables['mjj'] = {      'name'  : 'mjj',
                          #'range' : (25, 200., 400.),
                          'range' : (30, 0., 4000.),
                          'xaxis' : 'm_{jj} [GeV]',
                          'fold'  : 3}

variables['mjj_CR'] = {      'name'  : 'mjj',
                          'range' : (20, 0., 400.),
                          'xaxis' : 'm_{jj} [GeV]',
                          'fold'  : 3}

variables['mjj_WW'] = {      'name'  : 'mjj',
                             'range' : (10, 0., 200.),
                             'xaxis' : 'm_{jj} [GeV]',
                             'fold'  : 3}

variables['mjj_GGH'] = {      'name'  : 'mjj',
                              'range' : (10, 0., 300.),
                              'xaxis' : 'm_{jj} [GeV]',
                              'fold'  : 3}

variables['mll'] = {      'name'  : 'mll',
                          'range' : (20, 10., 70.),
                          'xaxis' : 'm_{ll} [GeV]',
                          'fold'  : 3}

variables['mll_top'] = {      'name'  : 'mll',
                          'range' : (20, 80., 200.),
                          'xaxis' : 'm_{ll} [GeV]',
                          'fold'  : 3}

variables['mll10GeV'] = { 'name'  : 'mll',
                          #'range' : (20, 0., 200.),
                          'range' : (10, 0., 100.),
                          'xaxis' : 'm_{ll} [GeV]',
                          'fold'  : 3}

variables['mpmet'] = {    'name'  : 'mpmet',      
                          'range' : (50, 0., 150.),  
                          'xaxis' : 'min. (proj. tk. E_{T}^{miss}, proj. E_{T}^{miss}) [GeV]', 
                          'fold'  : 3}

variables['mth'] = {      'name'  : 'mth',
                          #'range' : (40, 0., 200.),
                          'range' : (10, 40., 125.),
                          'xaxis' : 'm_{T}^{H} [GeV]',
                          'fold'  : 3}

variables['mth_top'] = {      'name'  : 'mth',
                              #'range' : (40, 0., 200.),                                                                                                                                                   
                              'range' : (40, 0., 200.),
                              'xaxis' : 'm_{T}^{H} [GeV]',
                              'fold'  : 3}

variables['mth_DY'] = {      'name'  : 'mth',
                             'range' : (20, 0., 40.),
                             'xaxis' : 'm_{T}^{H} [GeV]',
                             'fold'  : 3}

variables['mth_WW'] = {      'name'  : 'mth',
                             'range' : (12, 50., 200.),
                             'xaxis' : 'm_{T}^{H} [GeV]',
                             'fold'  : 3}

variables['mth_GGH'] = {      'name'  : 'mth',
                             'range' : (10, 60., 145.),
                             'xaxis' : 'm_{T}^{H} [GeV]',
                             'fold'  : 3}

variables['njet'] = {     'name'  : 'Sum$(CleanJet_pt>30)',     
                          'range' : (5, 0, 5),   
                          'xaxis' : 'number of jets',
                          'fold'  : 3}

variables['bjet'] = {     'name'  : 'Sum$(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.2217)',     
                          'range' : (5, 0, 5),   
                          'xaxis' : 'number of b jets',
                          'fold'  : 2}

variables['nvtx'] = {     'name'  : 'PV_npvsGood',      
                          'range' : (50, 0, 50),  
                          'xaxis' : 'number of vertices', 
                          'fold'  : 3}

variables['phi1'] = {     'name'  : 'Lepton_phi[0]',
                          'range' : (40, -3.2, 3.2),
                          'xaxis' : '#phi 1st lepton',
                          'fold'  : 3}

variables['phi2'] = {     'name'  : 'Lepton_phi[1]',
                          'range' : (40, -3.2, 3.2),
                          'xaxis' : '#phi 2nd lepton',
                          'fold'  : 3}

variables['pfmet'] = {    'name'  : 'MET_pt',     
                          'range' : (50, 0., 150.),   
                          'xaxis' : 'PF MET [GeV]',
                          'fold'  : 3}

variables['pt1'] = {      'name'  : 'Lepton_pt[0]',     
                          #'range' : (40, 0., 200.),   
                          'range' : (10, 20., 100.),
                          'xaxis' : 'p_{T} 1st lepton [GeV]',
                          'fold'  : 3}

variables['pt2'] = {      'name'  : 'Lepton_pt[1]',     
                          #'range' : (40, 0., 200.),   
                          'range' : (10, 13., 100.),
                          'xaxis' : 'p_{T} 2nd lepton [GeV]',
                          'fold'  : 3}

variables['ptll'] = {     'name'  : 'ptll',
                          #'range' : (40, 30., 200.),
                          'range' : (10, 30., 200.),
                          'xaxis' : 'p_{T}^{ll} [GeV]',
                          'fold'  : 3}

variables['puppimet'] = { 'name'  : 'PuppiMET_pt',
                          'range' : (50, 0., 150.),
                          'xaxis' : 'puppi MET [GeV]',
                          'fold'  : 3}

variables['rawmet'] = {   'name'  : 'RawMET_pt',     
                          'range' : (50, 0., 150.),   
                          'xaxis' : 'raw MET [GeV]',
                          'fold'  : 3}

variables['TkMET'] = {    'name'  : 'TkMET_pt',     
                          'range' : (50, 0., 150.),   
                          'xaxis' : 'tracker MET [GeV]',
                          'fold'  : 3}

variables['ptlljjmet'] = { 'name'  : 'Lepton_pt[0] + Lepton_pt[1] + CleanJet_pt[0] + CleanJet_pt[1] + MET_pt',
                          'range' : (12, 100., 700.),
                          'xaxis' : 'p_{T} leptons+jets+MET',
                          'fold'  : 3}

variables['mlljj'] = {    'name'       : 'mlljj',
                          'range'      : (15, 0., 4000.),
                          'xaxis'      : 'm_{lljj}',
                          'fold'       : 2}

variables['mjjj'] = {    'name'       : 'mjjj',
                         'range'      : (10, 0., 3500.),
                         'xaxis'      : 'm_{jjj}',
                         'fold'       : 2}

variables['ptjjj'] = { 'name'  : '(CleanJet_pt[0] + CleanJet_pt[1] + CleanJet_pt[2])*(CleanJet_pt[2] > 30)',
                       'range' : (16, 100., 500.),
                       'xaxis' : 'p_{T} 3 jets',
                       'fold'  : 2}

variables['ptjj'] = { 'name'  : '(CleanJet_pt[0] + CleanJet_pt[1])',
                      'range' : (16, 0., 500.),
                      'xaxis' : 'p_{T} 2 jets',
                      'fold'  : 3}

'''

'''

variables['alphaj1'] = { 'name'  : 'alpha1',
                         'range' : (15, 0., 2.),
                         'xaxis' : '#alpha_{j1}',
                         'fold'  : 2}

variables['alphaj2'] = { 'name'  : 'alpha2',
                         'range' : (15, 0., 2.),
                         'xaxis' : '#alpha_{j2}',
                         'fold'  : 2}

variables['alphaj3'] = { 'name'  : 'alpha3',
                         'range' : (15, 0., 2.),
                         'xaxis' : '#alpha_{j3}',
                         'fold'  : 2}


variables['D_0-'] = { 'name'  : 'kd_vbf',
                         'range' : (20, 0., 1.),
                         'xaxis' : 'D_{0}^{-}',
                         'fold'  : 2}

variables['D_0-_2'] = { 'name'  : 'kd_vbf',
                        'range' : (20, 0.65, 1.),
                        'xaxis' : 'D_{0}^{-}',
                        'fold'  : 3}

variables['D_alt_VBF'] = { 'name'  : 'D_alt',
                           'range' : (20, 0., 1.),
                           'xaxis' : 'D_{alt}',
                           'fold'  : 0}

variables['D_alt_qcd'] = { 'name'  : 'D_alt_QCD',
                         'range' : (20, 0., 1.),
                         'xaxis' : 'D_{alt} QCD',
                         'fold'  : 2}

variables['MELA_fit'] = {      'name'  : 'D_alt',
                               'range' : ([0.8, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.],),
                               'xaxis' : 'D_{alt}',
                               'fold'  : 0}

variables['D_int_VBF'] = { 'name'  : 'D_int',
                           'range' : (20, -1., 1.),
                           'xaxis' : 'D_{int}',
                           'fold'  : 2}

variables['Q2V1'] = { 'name'  : 'sqrt(Q2V1)',
                      'range' : (25, 0., 500.),
                      'xaxis' : '#sqrt{q_{1}^{2 VBF}}',
                      'fold'  : 3}

variables['Q2V2'] = { 'name'  : 'sqrt(Q2V2)',
                      'range' : (25, 0., 500.),
                      'xaxis' : '#sqrt{q_{2}^{2 VBF}}',
                      'fold'  : 3}

variables['CosTheta1'] = { 'name'  : 'costheta1',
                           'range' : (20, -1., 1.),
                           'xaxis' : 'cos #theta_{1}^{VBF}',
                           'fold'  : 3}

variables['CosTheta2'] = { 'name'  : 'costheta2',
                           'range' : (20, -1., 1.),
                           'xaxis' : 'cos #theta_{2}^{VBF}',
                           'fold'  : 3}

variables['CosThetaStar'] = { 'name'  : 'costhetastar',
                           'range' : (20, -1., 1.),
                           'xaxis' : 'cos #theta_{*}^{VBF}',
                           'fold'  : 3}

variables['Phi'] = { 'name'  : 'phi',
                     'range' : (20, -3., 3.),
                     'xaxis' : '#Phi^{VBF}',
                     'fold'  : 3}

variables['Phi_1'] = { 'name'  : 'phi1',
                       'range' : (20, -3., 3.),
                       'xaxis' : '#Phi_{1}^{VBF}',
                       'fold'  : 3}


'''
