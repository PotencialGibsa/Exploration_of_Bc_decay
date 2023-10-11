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
'Particle   K+            4.9367700e-01   0.0000000e+00',
'Particle   K-            4.9367700e-01   0.0000000e+00', 
'Particle   J/psi         3.0969000e+00   9.2600000e-05',
'Particle   B_c+          6.2744700e+00   0.0000000e+00',
'Particle   B_c-          6.2744700e+00   0.0000000e+00',
'Particle   B+            5.2793400e+00   0.0000000e+00',
'Particle   B-            5.2793400e+00   0.0000000e+00',
'Particle   K*0           8.9555000e-01   4.7300000e-02',
'Particle   anti-K*0      8.9555000e-01   4.7300000e-02',
'Particle   B*+           5.3247000e+00   6.5819970e-06',
'Particle   B*-           5.3247000e+00   6.5819970e-06',
'Particle   B_s10         5.8287000e+00   5.0000000e-04',
'Particle   anti-B_s10    5.8287000e+00   5.0000000e-04',
'Particle   B_s2*0        5.8398600e+00   1.4900000e-03',
'Particle   anti-B_s2*0   5.8398600e+00   1.4900000e-03',


'Alias      MyK*0            K*0',
'Alias      My-anti-K*0      anti-K*0',
'ChargeConj My-anti-K*0      MyK*0',

'Alias      MyB+             B+',
'Alias      MyB-             B-',
'ChargeConj MyB-             MyB+',

'Alias      MyB*+            B*+',
'Alias      MyB*-            B*-',
'ChargeConj MyB*-            MyB*+',

'Alias      MyB_s10          B_s10',
'Alias      My-anti-B_s10    anti-B_s10',
'ChargeConj My-anti-B_s10    MyB_s10',

'Alias      MyB_s2*0         B_s2*0',
'Alias      My-anti-B_s2*0   anti-B_s2*0',
'ChargeConj My-anti-B_s2*0   MyB_s2*0',

'Alias      MyBc+            B_c+',
'Alias      MyBc-            B_c-',
'ChargeConj MyBc-            MyBc+',

'Alias      MyJ/psi  J/psi',
'ChargeConj MyJ/psi  MyJ/psi',

'Decay MyJ/psi',
'1.000         mu+         mu-      PHOTOS VLL;',
'Enddecay',

'Decay MyK*0',
'1.000 K- pi+       VSS;',
'Enddecay',
'CDecay My-anti-K*0',

'Decay MyB+',
'1.000 MyJ/psi             K+         SVS;',
'Enddecay',
'CDecay MyB-',

'Decay MyB*+',
'1.000        MyB+         gamma    VSP_PWAVE;',
'Enddecay',
'CDecay MyB*-',


'Decay MyB_s10',
'1.000        MyB*+      K-         VVS_PWAVE 0.0 0.0 0.0 0.0 1.0 0.0;',
'Enddecay',
'CDecay My-anti-B_s10',

'Decay MyB_s2*0', ## BR from BPH-16-003 results
'0.075         MyB*+     K-         TVS_PWAVE 0.0 0.0 1.0 0.0 0.0 0.0;', 
'0.925         MyB+      K-         TSS;',
'Enddecay',
'CDecay My-anti-B_s2*0',


'Decay MyBc+', #BR pulled out of a hat
'0.350         My-anti-K*0  MyB+        SVS;',
'0.200         My-anti-K*0  MyB*+       SVV_HELAMP 1.0 0.0 1.0 0.0 1.0 0.0;',
'0.130         MyB+         pi+   K-    PHSP;',
'0.130         MyB*+        pi+   K-    PHSP;',
'0.140         MyB_s2*0     pi+         STS;',
'0.050         MyB_s10      pi+         SVS;',
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
  MotherID        = cms.untracked.int32(521),
  ParticleID      = cms.untracked.int32(443),
  DaughterIDs     = cms.untracked.vint32(13, -13),
  MinPt           = cms.untracked.vdouble(2.6, 2.6),
  MinEta          = cms.untracked.vdouble(-2.6, -2.6),
  MaxEta          = cms.untracked.vdouble( 2.6,  2.6)
)

ProductionFilterSequence = cms.Sequence(generator*jpsifilter)

