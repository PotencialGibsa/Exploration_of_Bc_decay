import FWCore.ParameterSet.Config as cms
from GeneratorInterface.EvtGenInterface.EvtGenSetting_cff import *
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *

generator = cms.EDFilter("Pythia8HadronizerFilter",
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    pythiaHepMCVerbosity = cms.untracked.bool(True),
    comEnergy = cms.double(13000.0),

    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        processParameters = cms.vstring( # put below any needed pythia parameter
            '541:m0 = 6.27447',
            '541:tau0 = 0.153',
            #'ProcessLevel:all = off',
        ),
        parameterSets = cms.vstring(
            'pythia8CommonSettings',
            'pythia8CP5Settings',
            'processParameters',
        ),
    ),

    ExternalDecays = cms.PSet(
        EvtGen130 = cms.untracked.PSet(
            decay_table            = cms.string('GeneratorInterface/EvtGenInterface/data/DECAY_2014_NOLONGLIFE.DEC'),
            particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt_2014.pdl'),
            list_forced_decays     = cms.vstring(
                'MyBc+',
                'MyBc-',
            ),
            operates_on_particles  = cms.vint32(541, -541),
            convertPythiaCodes     = cms.untracked.bool(False),
            user_decay_embedded    = cms.vstring(
'Particle   pi+           1.3957039e-01   0.0000000e+00',
'Particle   pi-           1.3957039e-01   0.0000000e+00',
'Particle   J/psi         3.0969000e+00   9.2600000e-05',
'Particle   K+            4.9367700e-01   0.0000000e+00',
'Particle   K-            4.9367700e-01   0.0000000e+00',
'Particle   rho+          7.7511000e-01   1.4910000e-01',
'Particle   rho-          7.7511000e-01   1.4910000e-01',
'Particle   B_c+          6.2744700e+00   0.0000000e+00',
'Particle   B_c-          6.2744700e+00   0.0000000e+00',
'Alias      MyBc+            B_c+',
'Alias      MyBc-            B_c-',
'ChargeConj MyBc-            MyBc+',
'Alias      Myrho+           rho+',
'Alias      Myrho-           rho-',
'ChargeConj Myrho-           Myrho+',
'Alias      MyJ/psi  J/psi',
'ChargeConj MyJ/psi  MyJ/psi',
'Decay MyJ/psi',
'1.000         mu+         mu-      PHOTOS VLL;',
'Enddecay',
'Decay Myrho+',
'1.000         pi+         pi0   PHOTOS PHSP;',
'Enddecay',
'CDecay Myrho-',

'Decay MyBc+',
'0.708      MyJ/psi     pi+     PHOTOS SVS;',
'0.056      MyJ/psi     K+      PHOTOS SVS;',
'0.236      MyJ/psi     Myrho+  PHOTOS PHSP;',
'Enddecay',
'CDecay MyBc-',
'End'   ),
        ),
        parameterSets = cms.vstring('EvtGen130'),
    ),

)

jpsifilter = cms.EDFilter("PythiaDauVFilter",
  verbose         = cms.untracked.int32(0),
  NumberDaughters = cms.untracked.int32(2),
  MotherID        = cms.untracked.int32(541),
  ParticleID      = cms.untracked.int32(443),
  DaughterIDs     = cms.untracked.vint32(13, -13),
  MinPt           = cms.untracked.vdouble(2.6, 2.6),
  MinEta          = cms.untracked.vdouble(-2.6, -2.6),
  MaxEta          = cms.untracked.vdouble( 2.6,  2.6)
)

ProductionFilterSequence = cms.Sequence(generator*jpsifilter)

