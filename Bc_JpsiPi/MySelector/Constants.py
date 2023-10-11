#############################################################
PDG_MUON_MASS    =   0.1056583745
PDG_PION_MASS    =   0.13957061
PDG_PIOZ_MASS    =   0.1349770
PDG_KAON_MASS    =   0.493677
PDG_PROTON_MASS  =   0.9382720813
PDG_KSHORT_MASS  =   0.497611
PDG_KSHORT_DM    =   0.000013
PDG_KSHORT_TIME  =   0.8954 * 0.0000000001
PDG_KS_MASS      =   PDG_KSHORT_MASS
PDG_LAMBDA_MASS  =   1.115683
PDG_LAMBDA_DM    =   0.000006
PDG_LAMBDA_TIME  =   2.632 * 0.0000000001
PDG_SIGMA0_MASS  =   1.192642
PDG_XImunus_MASS =   1.32171
PDG_XImunus_DM   =   0.00007
PDG_XImunus_TIME =   1.639 * 0.0000000001
PDG_OMmunus_MASS =   1.67245
PDG_OMmunus_DM   =   0.00029
PDG_OMmunus_TIME =   0.821 * 0.0000000001
PDG_DPM_MASS     =   1.86965
PDG_DPM_DM       =   0.00005
PDG_DPM_TIME     =   1.040 * 0.000000000001
PDG_DZ_MASS      =   1.86483
PDG_DZ_DM        =   0.00005
PDG_DZ_TIME      =   0.4101 * 0.000000000001
PDG_DS_MASS      =   1.96834
PDG_DS_DM        =   0.00007
PDG_DS_TIME      =   0.504 * 0.000000000001
PDG_LAMCP_MASS   =   2.28646
PDG_LAMCP_DM     =   0.00031
PDG_LAMCP_TIME   =   2.00 * 0.0000000000001
PDG_XICZ_MASS    =   2.47087
PDG_XICZ_DM      =   0.00031
PDG_XICZ_TIME    =   1.12 * 0.0000000000001
PDG_XICP_MASS    =   2.46787
PDG_XICP_DM      =   0.00030
PDG_XICP_TIME    =   4.42 * 0.0000000000001
PDG_KSTARZ_MASS  =   0.89555
PDG_KSTARZ_GAMMA =   0.0473
PDG_KSTARP_MASS  =   0.89176
PDG_KSTARP_GAMMA =   0.0503
PDG_PHI_MASS     =   1.019461
PDG_PHI_GAMMA    =   0.004249
PDG_JPSI_MASS    =   3.096900
PDG_PSI2S_MASS   =   3.686097
PDG_X3872_MASS   =   3.87169
PDG_B_MASS      =   5.27934
PDG_BU_MASS      =   5.27932
PDG_BU_TIME      =   1.638 * 0.000000000001
PDG_B0_MASS      =   5.27963
PDG_B0_TIME      =   1.520 * 0.000000000001
PDG_BS_MASS      =   5.36689
PDG_BS_TIME      =   1.509 * 0.000000000001
PDG_BC_MASS      =   6.2749
PDG_BC_TIME      =   0.507 * 0.000000000001
PDG_LB_MASS      =   5.61960
PDG_LB_TIME      =   1.470 * 0.000000000001
PDG_XIBZ_MASS    =   5.7919
PDG_XIBZ_TIME    =   1.479 * 0.000000000001
PDG_XIBM_MASS    =   5.7970
PDG_XIBM_TIME    =   1.571 * 0.000000000001
PDG_OMBM_MASS    =   6.0461
PDG_OMBM_TIME    =   1.64 * 0.000000000001
PDG_C            =   29979245800. ### in cm/c
PDG_DSTR         =   2.01026
###  }}}
'''
md      = RooRealVar ( "md"     ,"M(D) [GeV]"               , PDG_DZ_MASS-0.03 , PDG_DZ_MASS+0.03  )
mds     = RooRealVar ( "mds"    ,"M(Dstar) [GeV]"           , 2.004, 2.019 )
#
dspt    = RooRealVar ( "dspt"   ,"dspt"                     , 3.0   , 33.0 )
dzpt    = RooRealVar ( "dzpt"   ,"dzpt"                     , 2.5   , 22.5  )
kspt    = RooRealVar ( "kspt"   ,"kspt"                     , 0.0   , 10.0 )
pipt    = RooRealVar ( "pipt"   ,"pipt"                     , 0.1   , 1.5 )
#
dset    = RooRealVar ( "dset"   ,"dset"                     , -2.6  , 2.6   )
dzet    = RooRealVar ( "dzet"   ,"dzet"                     , -2.6  , 2.6   )
kset    = RooRealVar ( "kset"   ,"kset"                     , -2.6  , 2.6   )
piet    = RooRealVar ( "piet"   ,"piet"                     , -2.6  , 2.6   )
#
dzds2   = RooRealVar ( "dzds2"  ,"DetSign"                  , 0.0   , 50    )
dzds3   = RooRealVar ( "dzds3"  ,"DetSign"                  , 0.0   , 80    )
#
dsvtxp  = RooRealVar ( "dsvtxp" ,"dsvtxp"                   , -0.1  , 1.1   )
dzvtxp  = RooRealVar ( "dzvtxp" ,"dzvtxp"                   , -0.1  , 1.1   )
#
ipsmin  = RooRealVar ( "ipsmin" ,"ipsmin"                   , 0.0   , 1000.0)
vari    = RooRealVar ( "vari"   ,"vari"                     , 0.0   , 5)
vdr     = RooRealVar ( "vdr"    ,"vdr"                      , 0.0   , 2)
vdz     = RooRealVar ( "vdz"    ,"vdz"                      , 0.0   , 5)
'''
