-- Load the library containing the matrix element
load_modules('/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/WW_leptonic_ME/build/libme_WW_leptonic_ME.so') -- path to your build folder

local lepton1 = declare_input("lepton1") -- Declare input particles, be awere of the Id of your leptons in the .cc file. They must agree.
local lepton2 = declare_input("lepton2")
local jet1 = declare_input("jet1")
local jet2 = declare_input("jet2")

-- Do not touch 
parameters = {
    energy = 13000.,
    top_mass = 173.,
    top_width = 1.491500,
    W_mass = 80.419002,
    W_width = 2.047600,
}

-- Do not touch at first
cuba = {
    relative_accuracy = 0.02,
    verbosity = 0,
    max_eval = 50000,
    n_start = 500,
    -- algorithm = "divonne"
}


BreitWignerGenerator.flatter_s13 = {
    ps_point = add_dimension(),
    mass = parameter('W_mass'),
    width = parameter('W_width')
}

BreitWignerGenerator.flatter_s24 = {
    ps_point = add_dimension(),
    mass = parameter('W_mass'),
    width = parameter('W_width')
}

BlockF.blockf = {
    p3 = lepton1.gen_p4,
    p4 = lepton2.gen_p4;

    s13 = 'flatter_s13::s',
    s24 = 'flatter_s24::s',
    q1 = add_dimension(),
    q2 = add_dimension()
}

Looper.looper = {
    solutions = 'blockf::solutions',
    path = Path('initial_state', 'WW', 'integrand')
}

--
-- Start of loop over solutions
--

    full_inputs = { lepton1.reco_p4, 'looper::particles/1',  lepton2.reco_p4, 'looper::particles/2'} -- keep order
    
    BuildInitialState.initial_state = {
        particles = full_inputs,
        do_transverse_boost = true
    }

    jacobians = { 'flatter_s13::jacobian', 'flatter_s24::jacobian', 'tf_p1::TF_times_jacobian', 'tf_p2::TF_times_jacobian', 'looper::jacobian' }

    MatrixElement.WW = {
      pdf = 'CT10nlo',
      pdf_scale = parameter('W_mass'),

      matrix_element = 'WW_leptonic_ME_sm_P1_Sigma_sm_gg_epvemumvmxuux',  -- Change name from your PROCESS/Subprocess folder. Name in the last part of .cc file
      matrix_element_parameters = {
          card = '/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/WW_leptonic_ME/Cards/param_card.dat' -- Change path to your Cards path
      },

      initialState = 'initial_state::partons',

      -- Order must agree the .cc file and full_inputs !!!!
      particles = {
        inputs = full_inputs,
        ids = {
          {
            pdg_id = -13,
            me_index = 1,
          },
          {
            pdg_id = 14,
            me_index = 2,
          },
          {
            pdg_id = 11,
            me_index = 3,
          },
          {
            pdg_id = -12,
            me_index = 4,
          },
        }
      },

      -- Other jacobians
      jacobians = jacobians
    }

    DoubleLooperSummer.integrand = { input = 'WW::output' }

integrand('integrand::sum')
