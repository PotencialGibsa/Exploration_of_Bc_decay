import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList

process = cms.Process("Demo")
OUTNAME = 'BFinder'; OUTNAMEr = OUTNAME + '.root'

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100 # 900
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'/store/data/Run2018D/Charmonium/AOD/12Nov2019_UL2018-v1/40000/F9A4F086-B7CD-5F4F-98EE-DC65013F8339.root',
'/store/data/Run2017E/Charmonium/AOD/09Aug2019_UL2017-v1/240006/C1F5620A-3946-164A-88DF-E936DF0752F6.root',
'/store/data/Run2018B/Charmonium/AOD/12Nov2019_UL2018-v1/70000/FD464491-5BDD-E14F-A0E1-AE0700731F57.root',
'/store/data/Run2018C/Charmonium/AOD/12Nov2019_UL2018-v1/270000/F4916A7F-A35A-A744-A490-ABBE7C0A0ECD.root',
'/store/data/Run2018C/Charmonium/AOD/12Nov2019_UL2018-v1/270000/F3A21298-797B-624B-BFCC-A7FC1FC85D6D.root'
    )
)
# process.source.lumisToProcess = LumiList.LumiList(filename = 'python/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_MuonPhys.txt').getVLuminosityBlockRange()
##
### JSON FILE
## https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/PromptReco/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON_MuonPhys.txt


from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v28', '')## for 2016-2018 UL
# process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_Prompt_v11', '')## for 2018 D PROMPT RECO
# 94X_dataRun2_ReReco_EOY17_v2

## select generaltracks TRACKS TRACKS TRACKS TRACKS TRACKS
process.oniaSelectedTracks=cms.EDFilter("TrackSelector",
      src = cms.InputTag("generalTracks"),
      cut = cms.string('pt > 0.05 && abs(eta) <= 3.0'
                                '&& charge !=0'
                                '&& quality(\"loose\")'
                                )
)

#convert to PAT
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.CandidateSelectedTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
            src=cms.InputTag("oniaSelectedTracks"),
            particleType=cms.string('pi+')
)
from PhysicsTools.PatAlgos.producersLayer1.genericParticleProducer_cfi import patGenericParticles
process.patSelectedTracks = patGenericParticles.clone(src=cms.InputTag("CandidateSelectedTracks"))

#make muon with trigger matching embedded
### ==== Make PAT Muons ====
import PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi
process.oniaPATMuonsWithoutTrigger = PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi.patMuons.clone(
    muonSource = 'muons',
    # embed the tracks, so we don't have to carry them around
    embedTrack          = True,
    embedCombinedMuon   = True,
    embedStandAloneMuon = True,
    embedPFCandidate    = False,
    embedCaloMETMuonCorrs = cms.bool(False),
    embedTcMETMuonCorrs   = cms.bool(False),
    embedPfEcalEnergy     = cms.bool(False),
    # then switch off some features we don't need
    embedPickyMuon = False,
    embedTpfmsMuon = False,
    userIsolation = cms.PSet(),   # no extra isolation beyond what's in reco::Muon itself
    isoDeposits = cms.PSet(), # no heavy isodeposits
    addGenMatch = False,       # no mc: T&P doesn't take it from here anyway.
)
# Reset all these; the default in muonProducer_cfi is not empty, but wrong
process.oniaPATMuonsWithoutTrigger.userData.userInts.src    = []
process.oniaPATMuonsWithoutTrigger.userData.userFloats.src  = []
process.oniaPATMuonsWithoutTrigger.userData.userCands.src   = []
process.oniaPATMuonsWithoutTrigger.userData.userClasses.src = []
### ==== Unpack trigger, and match ====
from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import patTrigger as oniaPATTriggerTMP_i
process.oniaPATTriggerTMP = oniaPATTriggerTMP_i
process.oniaPATTriggerTMP.onlyStandAlone = True
process.oniaPATTrigger = cms.EDProducer("TriggerObjectFilterByCollection",
    src = cms.InputTag("oniaPATTriggerTMP"),
    collections = cms.vstring("hltL2MuonCandidates", "hltIterL3MuonCandidates", "hltHighPtTkMuonCands", "hltGlbTrkMuonCands")
)

### ==== Then perform a match for all HLT triggers of interest
process.PATmuonTriggerMatchHLT = cms.EDProducer( "PATTriggerMatcherDRDPtLessByR",
    src     = cms.InputTag( "oniaPATMuonsWithoutTrigger" ),
    matched = cms.InputTag( "oniaPATTrigger" ),
    matchedCuts = cms.string(""),
    maxDPtRel = cms.double( 0.5 ),
    maxDeltaR = cms.double( 0.5 ),
    resolveAmbiguities    = cms.bool( True ),
    resolveByMatchQuality = cms.bool( True ) #change with respect to previous tag
)

process.PATmuonMatchHLTL2   = process.PATmuonTriggerMatchHLT.clone(matchedCuts = cms.string('coll("hltL2MuonCandidates")'),
                                                   maxDeltaR = 0.3, maxDPtRel = 10.0)       #maxDeltaR Changed accordingly to Zoltan tuning. It was: 1.2
process.PATmuonMatchHLTL3   = process.PATmuonTriggerMatchHLT.clone(matchedCuts = cms.string('coll("hltIterL3MuonCandidates")'),
                                                   maxDeltaR = 0.1, maxDPtRel = 10.0)       #maxDeltaR Changed accordingly to Zoltan tuning. It was: 0.5
process.PATmuonMatchHLTL3T  = process.PATmuonTriggerMatchHLT.clone(matchedCuts = cms.string('coll("hltGlbTrkMuonCands")'),
                                                   maxDeltaR = 0.1, maxDPtRel = 10.0)       #maxDeltaR Changed accordingly to Zoltan tuning. It was: 0.5
process.PATmuonMatchHLTTkMu = process.PATmuonTriggerMatchHLT.clone(matchedCuts = cms.string('coll("hltHighPtTkMuonCands")'),
                                                   maxDeltaR = 0.1, maxDPtRel = 10.0)       #maxDeltaR Changed accordingly to Zoltan tuning. It was: 0.5

process.oniaPATTriggerMatchers1Mu = cms.Sequence(
      process.PATmuonMatchHLTL2 +
      process.PATmuonMatchHLTL3 +
      process.PATmuonMatchHLTL3T +
      process.PATmuonMatchHLTTkMu
)

oniaPATTriggerMatchers1MuInputTags = [
    cms.InputTag('PATmuonMatchHLTL2'),
    cms.InputTag('PATmuonMatchHLTL3'),
    cms.InputTag('PATmuonMatchHLTL3T'),
    cms.InputTag('PATmuonMatchHLTTkMu'),
]

## ==== Embed ====
process.oniaPATMuonsWithTrigger = cms.EDProducer( "PATTriggerMatchMuonEmbedder",
    src     = cms.InputTag(  "oniaPATMuonsWithoutTrigger" ),
    matches = cms.VInputTag()
)
process.oniaPATMuonsWithTrigger.matches += oniaPATTriggerMatchers1MuInputTags

## ==== Trigger Sequence ====
process.oniaPATTriggerMatching = cms.Sequence(
    process.oniaPATTriggerTMP * process.oniaPATTrigger *
    process.oniaPATTriggerMatchers1Mu *
    process.oniaPATMuonsWithTrigger
)

process.oniaPATMuonsWithTriggerSequence = cms.Sequence(
    process.oniaPATMuonsWithoutTrigger *
    process.oniaPATTriggerMatching
)

# process.load('HeavyFlavorAnalysis.Onia2MuMu.oniaPATMuonsWithTrigger_cff')
#select them, adjust as needed, currently soft muon selection with pT>3.5
process.oniaSelectedMuons = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('oniaPATMuonsWithTrigger'),
   cut = cms.string('muonID(\"TMOneStationTight\")' ),
   filter = cms.bool(True)
)


process.demo = cms.EDAnalyzer('Xb_frame',
    HLTriggerResults    = cms.InputTag("TriggerResults","","HLT"),
    beamSpotTag         = cms.InputTag("offlineBeamSpot"),
    VtxSample           = cms.InputTag('offlinePrimaryVerticesWithBS'),
    Trak                = cms.InputTag('patSelectedTracks'), #selectedPatTracks
    muons               = cms.InputTag('oniaSelectedMuons'), #oniaPATMuonsWithoutTrigger
    revtxtrks           = cms.InputTag('generalTracks'),
#      trackQualities     = cms.untracked.vstring('loose','tight','highPurity'),
    TriggerOSATag       = cms.untracked.InputTag("selectedPatTrigger", "", "PAT"),
    triggerObjects      = cms.InputTag("slimmedPatTrigger", "", "PAT"),
    fileName            = cms.untracked.string(OUTNAMEr),
    TriggerEvent = cms.InputTag("hltTriggerSummaryAOD"),
    TriggerFilters = cms.vstring(["hltJpsiTkVertexFilter"]),
#     TrackLabel = cms.InputTag('cleanPatTrackCands'),
#     trackRecoAlgorithm = cms.InputTag('generalTracks'),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(OUTNAMEr)
)


process.mySequence = cms.Sequence(
                   process.oniaPATMuonsWithTriggerSequence *
                   process.oniaSelectedMuons *
                   process.oniaSelectedTracks *
                   process.CandidateSelectedTracks *
                   process.patSelectedTracks *
                   process.demo
)

# process.dump=cms.EDAnalyzer('EventContentAnalyzer')
# process.p = cms.Path( process.dump * process.mySequence)
process.p = cms.Path( process.mySequence)
process.schedule = cms.Schedule(process.p)
#
# process.dump=cms.EDAnalyzer('EventContentAnalyzer')
# process.p = cms.Path(process.dump * process.demo)
