import os
from CRABClient.UserUtilities import config
config = config()

task = '2018_Bc_in_BKPI_wrong_charge_trigger_matching'

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
config.Data.outLFNDirBase = '/store/user/dshmygol/Bc_Bfinder/'+task+'/2018/'

config.section_("Site")
# config.Site.storageSite = 'T2_RU_IHEP'
#config.Site.storageSite = 'T2_RU_JINR'
config.Site.storageSite = 'T3_CH_CERNBOX'

DS_names = [ '' , ## 5 items
'/Charmonium/Run2018A-12Nov2019_UL2018_rsb-v1/AOD', ## 1
'/Charmonium/Run2018B-12Nov2019_UL2018-v1/AOD',
'/Charmonium/Run2018C-12Nov2019_UL2018_rsb_v2-v2/AOD',
'/Charmonium/Run2018D-12Nov2019_UL2018-v1/AOD',
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
    ##
    ## ULpl1, LLpl1 = Lb -> psi2S Lambda, JPP
    ##
    ## Uxi4a, Lxi4a = Xi- -> Lambda pi
    #
    ## ULbumm, B+mumu, B+->psiK+
    #
    ## ULBdxx , B0 -> [J p p p] K
    #
    n = 0
    if len(sys.argv) >= 2:
        ss = str(sys.argv[1])
        if ss == 'all':
            #
            for i in range(1, len(DS_names)):
                print '\nrunning\npython ULRunCharmonium_crab_2018.py %i'%i
                os.system('python ULRunCharmonium_crab_2018.py %i'%i)
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
    #
    #if n==3:
    #   config.Data.ignoreLocality = True
    #   config.Site.ignoreGlobalBlacklist = True
    #   config.Site.whitelist = ['T2_RU_*', 'T2_CH_CERN']
    #
    #
    config.General.requestName = 'Bfinder_' + task + '_' + ['','A','B','C','D'][n]
    #config.General.requestName = 'Bfinder_' + task + '_' + ['','A','','',''][n]
    config.General.workArea = 'crab_projects'
    config.Data.inputDataset = dset
    print '\n', config.General.requestName
    print config.General.workArea
    print dset, '\n'
    #
    lumi_mask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON_MuonPhys.txt'
    # lumi_mask = '/afs/cern.ch/work/s/spolikar/ULRUN2/CMSSW_10_6_12/src/XbFrame/Xb_frame/crab_projects/crab_Bfinder_2018_Lxi4a_C/results/notFinishedLumis.json'
    #
    split_modifier = 1;
    # if ('8B-17' in dset) : split_modifier = 0.9
    #
    config.Data.unitsPerJob = int(units_per_job * split_modifier)
    #config.Data.lumiMask = lumi_mask
    submit(config)
