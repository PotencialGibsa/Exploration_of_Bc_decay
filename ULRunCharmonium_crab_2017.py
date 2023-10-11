import os
from CRABClient.UserUtilities import config
config = config()

task = '2017_Bc_in_BKPI_wrong_charge_trigger_matching'

config.section_("General")
config.General.requestName = 'Bfinder_' + task
config.General.workArea = 'crab_projects_' + task
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'python/ULConfFile_cfg_RUN2_Data.py'
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.inputDataset = '/MuOnia/Run2012B-22Jan2013-v1/AOD'
config.Data.inputDBS =	'global'
# config.Data.splitting = 'Automatic'
config.Data.splitting = 'LumiBased'
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'CRAB3_Bfinder'
config.Data.outLFNDirBase = '/store/user/dshmygol/Bc_Bfinder/'+task+'/2017/'

config.section_("Site")
# config.Site.storageSite = 'T2_RU_IHEP'
#config.Site.storageSite = 'T2_RU_JINR'
config.Site.storageSite = 'T3_CH_CERNBOX'

DS_names = [ '' , ## 5 items
'/Charmonium/Run2017B-09Aug2019_UL2017-v1/AOD', ## 1
'/Charmonium/Run2017C-09Aug2019_UL2017-v1/AOD', ## 2
'/Charmonium/Run2017D-09Aug2019_UL2017-v1/AOD', ## 3
'/Charmonium/Run2017E-09Aug2019_UL2017-v1/AOD', ## 4
'/Charmonium/Run2017F-09Aug2019_UL2017-v1/AOD', ## 5
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
    units_per_job = 14  ## default
    # units_per_job = 10 ## for recovery
    # units_per_job = 1400 ## for automatic
    #
    n = 0
    print(sys.argv, len(sys.argv))
    print(len(DS_names))
    if len(sys.argv) >= 2:
        ss = str(sys.argv[1])
        if ss == 'all':
            #
            for i in range(1, len(DS_names)):
                print '\nrunning\npython ULRunCharmonium_crab_2017.py %i'%i
                os.system('python ULRunCharmonium_crab_2017.py %i'%i)
            exit(0);
            #
        elif ss == '161718':
            #
            print '\nrunning\npython ULRunCharmonium_crab_2016.py all\npython ULRunCharmonium_crab_2017.py all\npython ULRunCharmonium_crab_2018.py all'
            os.system('python ULRunCharmonium_crab_2016.py all')
            os.system('python ULRunCharmonium_crab_2017.py all')
            os.system('python ULRunCharmonium_crab_2018.py all')
            exit(0);
            #
        elif ss == '1718':
            #
            print '\nrunning\npython ULRunCharmonium_crab_2017.py all\npython ULRunCharmonium_crab_2018.py all'
            os.system('python ULRunCharmonium_crab_2017.py all')
            os.system('python ULRunCharmonium_crab_2018.py all')
            exit(0);
            #
        else:
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
    # if dset == '/Charmonium/Run2018D-PromptReco-v2/AOD':
    #if n==3:
    #   config.Data.ignoreLocality = True
    #   config.Site.ignoreGlobalBlacklist = True
    #   config.Site.whitelist = ['T2_RU_*', 'T2_CH_CERN']
    #
    #
    config.General.requestName = 'Bfinder_' + task + '_' + ['','B','C','D','E','F'][n]
    #config.General.requestName = 'Bfinder_' + task + '_' + ['','B','F'][n]
    config.General.workArea = 'crab_projects'
    config.Data.inputDataset = dset
    print '\n', config.General.requestName
    print config.General.workArea
    print dset, '\n'
    #
    lumi_mask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_MuonJSON.txt'
    #lumi_mask = '/afs/cern.ch/work/m/msergeev/CMSSW_10_6_12/src/XbFrame/Xb_frame/crab_projects/crab_Bfinder_2017_Lxi4a_C/results/notFinishedLumis.json'
    #
    split_modifier = 1;
    # if ('8B-17' in dset) : split_modifier = 0.9
    #
    config.Data.unitsPerJob = int(units_per_job * split_modifier)
    #config.Data.lumiMask = lumi_mask
    submit(config)
