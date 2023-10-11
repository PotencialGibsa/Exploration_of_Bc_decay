import os
from CRABClient.UserUtilities import config
config = config()

task = '2015_MC_Bc_in_BKPI_WS_v6_trigger_matching_improved'

config.section_("General")
config.General.requestName = 'Bfinder_' + task
config.General.workArea = 'crab_projects_' + task
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'python/ULConfFile_cfg_RUN2_MC15.py'
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.inputDataset = '/MuOnia/Run2012B-22Jan2013-v1/AOD'
config.Data.inputDBS =	'global'
# config.Data.splitting = 'Automatic'
config.Data.splitting = 'FileBased'
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'CRAB3_Bfinder'
config.Data.outLFNDirBase = '/store/user/dshmygol/Bc_Bfinder/'+task

config.section_("Site")
# config.Site.storageSite = 'T2_RU_IHEP'
#config.Site.storageSite = 'T2_RU_JINR'
config.Site.storageSite = 'T3_CH_CERNBOX'

DS_names = [ '' , ## 5 items
"/BcToBuKPi_BuJPsiK_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/AODSIM",
"/BcToBuKPi_BuJPsiK_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8_ext1-v2/AODSIM"
#"/BcToJpsPi_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8_ext1-v2/AODSIM",
#"/BcToJpsPi_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/AODSIM"
]

if __name__ == '__main__':
    print 'multisubmit.\nunitsPerJob ~ 10 for maximum splitting\nHave you done scram b -j8?!\n'
    print 'STRONGLY SUGGEST: increase reportEvery to at least 800 before submission'
    print 'STRONGLY SUGGEST: remove frequent printouts (cout) from the code'
    import sys
    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException
    ####
    ####
    #units_per_job = 25  ## default
    units_per_job = 10 ## for recovery
    # units_per_job = 1400 ## for automatic
    ##
    ## ULbumm, B+mumu, B+->psiK+
    #
    n = 0
    if len(sys.argv) >= 2:
        ss = str(sys.argv[1])
	n = int(ss)
        #
    else:
        print '\nADD A NUMBER TO SET DATASET ... 0 - %i or write "all" to run on all 2018 data or "161718" to run on all datasets'%(len(DS_names) - 1)
        exit(0)

    def submit(cfg):
		try:
			print crabCommand('submit', config = cfg)
		except HTTPException, hte:
			print hte.headers
    #
    dset = DS_names[n]
    #
    #
    config.General.requestName = 'Bfinder_' + task + '_' + str(n)
    config.General.workArea = 'crab_projects'
    config.Data.inputDataset = dset
    print '\n', config.General.requestName
    print config.General.workArea
    print dset, '\n'
    #
    split_modifier = 1;
    # if ('8B-17' in dset) : split_modifier = 0.9
    #
    config.Data.unitsPerJob = int(units_per_job * split_modifier)
    #config.Data.lumiMask = lumi_mask
    submit(config)
